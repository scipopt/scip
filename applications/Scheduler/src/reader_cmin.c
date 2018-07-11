/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_cmin.c
 * @brief  cmin file reader
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h>
#endif
#include <ctype.h>

#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/cons_knapsack.h"

#include "cons_optcumulative.h"
#include "heur_optcumulative.h"
#include "reader_cmin.h"

#define READER_NAME             "cminreader"
#define READER_DESC             "file reader for cmin file format"
#define READER_EXTENSION        "cmin"

#define DEFAULT_FILENAME           "-"  /**< name of the file including best known solutions */
#define DEFAULT_DUALREDUCTION     TRUE  /**< add locks to avoid dual reductions */
#define DEFAULT_MIP              FALSE  /**< create a MIP formulation */
#define DEFAULT_INITIAL           TRUE  /**< should model constraints be in initial LP? */
#define DEFAULT_CIP               TRUE  /**< create a CIP formulation */
#define DEFAULT_RELAXATION           3  /**< which relaxation should be added to the maseter (0: none; 1: single; 2: edge-finding; 3: energetic-reasoning */

static const char delimchars[] = " \f\n\r\t\v";


/*
 * Local methods
 */

/** CP reading data */
struct CminInput
{
   SCIP_FILE*           file;
   char                 linebuf[SCIP_MAXSTRLEN];
   char*                token;
   int                  linenumber;
   int                  linepos;
   SCIP_Bool            haserror;
};
typedef struct CminInput CMININPUT;

/** issues an error message and marks the LP data to have errors */
static
void syntaxError(
   SCIP*                 scip,               /**< SCIP data structure */
   CMININPUT*            cmininput,          /**< CMIN reading data */
   const char*           msg                 /**< error message */
   )
{
   assert(cmininput != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error in line %d: %s ('%s')\n",
      cmininput->linenumber, msg, cmininput->token);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s\n", cmininput->linebuf);

   cmininput->haserror = TRUE;
}

/** gets the next line out of the file stream */
static
SCIP_Bool getNextLine(
   CMININPUT*              cmininput         /**< CMIN reading data */
   )
{
   char* pch;
   assert( cmininput->haserror == FALSE );

   /* clear the line */
   BMSclearMemoryArray(cmininput->linebuf, SCIP_MAXSTRLEN);

   if( SCIPfgets(cmininput->linebuf, (int) sizeof(cmininput->linebuf), cmininput->file) == NULL )
      return FALSE;

   cmininput->linenumber++;
   cmininput->linepos = 0;

   cmininput->linebuf[SCIP_MAXSTRLEN-1] = '\0';

   pch = strstr (cmininput->linebuf, "\n");
   if( pch != NULL )
      *pch = '\0';

   if (cmininput->linebuf[0] != '\0')
   {
      SCIPdebugMessage("input line %d: <%s>\n", cmininput->linenumber, cmininput->linebuf);
      return TRUE;
   }

   return FALSE;
}

/** returns whether the given character is a token delimiter */
static
SCIP_Bool isDelimChar(
   char                  c                   /**< input character */
   )
{
   return (c == '\0') || (strchr(delimchars, c) != NULL);
}

/** reads the next token from the input file into the token buffer; returns whether a token was read */
static
SCIP_Bool getNextToken(
   CMININPUT*            cmininput           /**< CMIN reading data */
   )
{
   char* buf;
   int tokenlen;

   assert(cmininput != NULL);

   /* skip delimiters */
   buf = cmininput->linebuf;
   while( isDelimChar(buf[cmininput->linepos]) )
   {
      if( buf[cmininput->linepos] == '\0' )
      {
         if( !getNextLine(cmininput) )
         {
            SCIPdebugMessage("(line %d) end of file\n", cmininput->linenumber);
            return FALSE;
         }
         assert(cmininput->linepos == 0);
      }
      else
         cmininput->linepos++;
   }
   assert(cmininput->linepos < SCIP_MAXSTRLEN);
   assert(!isDelimChar(buf[cmininput->linepos]));

   /* read value token */
   tokenlen = 0;
   while( isdigit(buf[cmininput->linepos]) )
   {
      assert(tokenlen < SCIP_MAXSTRLEN);
      assert(!isDelimChar(buf[cmininput->linepos]));
      cmininput->token[tokenlen] = buf[cmininput->linepos];
      tokenlen++;
      cmininput->linepos++;
   }

   assert(tokenlen < SCIP_MAXSTRLEN);
   cmininput->token[tokenlen] = '\0';

   SCIPdebugMessage("(line %d) read token: '%s'\n", cmininput->linenumber, cmininput->token);

   return TRUE;
}

/** method parses the best known solution for the total leftover out of the give file; if file does not exist or the
 *  problem is not listed the best known solution is set to -1 which means unknown
 */
static
SCIP_RETCODE findBestObjectiveValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            objval              /**< pointer to store best known solution */
   )
{
   SCIP_FILE* file;
   const char* probname;
   SCIP_Bool found;
   char* solufile;
   char buffer[SCIP_MAXSTRLEN];
   char* token;
   char* saveptr;

   assert(scip != NULL);
   assert(objval != NULL);

   /* get solution file  */
   SCIP_CALL( SCIPgetStringParam(scip, "reading/"READER_NAME"/filename", &solufile) );

   /* set objective value to unknown */
   (*objval) = SCIPinfinity(scip);

   if( *solufile == '-')
      return SCIP_OKAY;

   if( NULL == (file = SCIPfopen(solufile, "r")) )
   {
      SCIPwarningMessage(scip, "solution file <%s> does not exist\n", solufile);
      return SCIP_OKAY;
   }

   /* get problem name */
   probname = SCIPgetProbName(scip);
   found = FALSE;

   /* parse file line by line */
   while( SCIPfgets(buffer, (int) sizeof(buffer), file) != NULL )
   {
      char status[SCIP_MAXSTRLEN];

      SCIPdebugMessage("parse line <%s>\n", buffer);

      /* solution status */
      token = SCIPstrtok(buffer, " \n\t", &saveptr);
      (void)SCIPsnprintf(status, SCIP_MAXSTRLEN, "%s", token);

      /* problem name */
      token = SCIPstrtok(NULL, " \n\t", &saveptr);

      /* check if found the right line for the requested problem */
      if( strncmp(token, probname, strlen(token)+1) == 0 )
      {
         SCIPdebugMessage("found problem <%s> in solution file <%s>\n", probname, solufile);

         if( strncmp(status, "=opt=", 5) == 0 )
         {
            /* objective value */
            token = SCIPstrtok(NULL, " \n\t", &saveptr);

            (*objval) = atof(token);
            found = TRUE;

            SCIPdebugMessage("best known solution is <%g>\n", *objval);
         }

         break;
      }
   }

   SCIPfclose(file);

   if( !found )
      return SCIP_INVALIDCALL;

   return SCIP_OKAY;
}

#ifdef SCIP_DEBUG
/** display input data */
static
void displayInputData(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 durations,          /**< durations */
   int**                 demands,            /**< demands */
   int**                 costs,              /**< cost */
   int*                  capacities,         /**< machine capacities */
   int*                  releasedates,       /**< release dates */
   int*                  deadlinedates,      /**< deadline dates */
   int                   njobs,              /**< number of jobs */
   int                   nmachines           /**< number of machines */
   )
{
   int j;
   int m;

   for( j = 0; j < njobs; ++j )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Job %2d [%3d,%3d] ", j, releasedates[j], deadlinedates[j]);

      for( m = 0; m < nmachines; ++m )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "(%2d,%2d)" , durations[m][j], demands[m][j]);
      }
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "\n");
   }
}
#endif

/** initialize the sorted event point arrays */
static
void createSortedEventpoints(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  releasedates,       /**< release dates */
   int*                  deadlinedates,      /**< deadline dates */
   int*                  starttimes,         /**< array to store sorted start events */
   int*                  endtimes,           /**< array to store sorted end events */
   int*                  startindices,       /**< permutation with rspect to the start times */
   int*                  endindices,         /**< permutation with rspect to the end times */
   int                   njobs               /**< number of jobs */
   )
{
   int j;

   /* assign variables, start and endpoints to arrays */
   for ( j = 0; j < njobs; ++j )
   {
      starttimes[j] = releasedates[j];
      startindices[j] = j;

      endtimes[j] = deadlinedates[j];
      endindices[j] = j;
   }

   /* sort the arrays not-decreasing according to startsolvalues and endsolvalues (and sort the indices in the same way) */
   SCIPsortIntInt(starttimes, startindices, njobs);
   SCIPsortIntInt(endtimes, endindices, njobs);
}

/** computes the maximum energy for all variables which correspond to jobs which start between the given start time and
 *  end time
 *
 *  @return Maximum energy for the given time window
 */
static
SCIP_Longint computeMaxEnergy(
   int                   njobs,              /**< number of jobs */
   int*                  durations,          /**< durations */
   int*                  demands,            /**< demands */
   int*                  releasedates,       /**< release dates */
   int*                  deadlinedates,      /**< deadline dates */
   int                   starttime,          /**< start time */
   int                   endtime             /**< end time */
   )
{
   SCIP_Longint maxenergy;
   int v;

   assert(starttime < endtime);
   maxenergy = 0LL;

   for( v = 0; v < njobs; ++v )
   {
      /* collect jobs which run between the start and end time */
      if( deadlinedates[v] <= endtime && releasedates[v] >= starttime)
         maxenergy += (SCIP_Longint)(durations[v] * demands[v]); /*lint !e647*/
   }

   return maxenergy;
}

/** remove row which have a tightness which is smaller or equal to the given one
 *
 *  @return The number of remaining rows
 */
static
int removeRedundantRows(
   SCIP_Longint*         rowtightness,       /**< array containing the tightness for the previously selected rows */
   int*                  startidxs,          /**< array containing for each row the index for the start event */
   int                   nrows,              /**< current number of rows */
   SCIP_Longint          tightness           /**< tightness to use to detect redundant rows */
   )
{
   int keptrows;
   int j;

   keptrows = 0;

   for( j = 0; j < nrows; ++j )
   {
      rowtightness[keptrows] = rowtightness[j];
      startidxs[keptrows] = startidxs[j];

      /* only keep this row if the tightness is better as the (current) given one */
      if( rowtightness[j] > tightness )
         keptrows++;
   }

   return keptrows;
}

/** create interval relaxation for given sub-problem */
static
SCIP_RETCODE createIntervalRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   relaxation,         /**< a linear relaxation base on edge-finding idea or energetic-reasoning idea */
   int                   resource,           /**< resource id */
   SCIP_VAR**            vars,               /**< assignment variables */
   int*                  durations,          /**< durations */
   int*                  demands,            /**< demands */
   int                   capacity,           /**< machine capacity */
   int*                  releasedates,       /**< release dates */
   int*                  deadlinedates,      /**< deadline dates */
   int                   njobs               /**< number of jobs */
   )
{
   SCIP_Longint** rowtightness;
   int** startidxs;
   int* nrows;
   int* starttimes;
   int* endtimes;
   int* startindices;
   int* endindices;
   int starttime;
   int endtime;
   int i;
   int j;

   int naddedcons;

   SCIP_CALL( SCIPallocBufferArray(scip, &starttimes, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endtimes, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, njobs) );

   SCIP_CALL( SCIPallocBufferArray(scip, &nrows, njobs) );
   BMSclearMemoryArray(nrows, njobs);
   SCIP_CALL( SCIPallocBufferArray(scip, &rowtightness, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startidxs, njobs) );
   for( j = 0; j < njobs; ++j )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &rowtightness[j], njobs) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &startidxs[j], njobs) ); /*lint !e866*/
   }

   createSortedEventpoints(scip, releasedates, deadlinedates, starttimes, endtimes, startindices, endindices, njobs);

   starttime = -INT_MAX;

   /* check each startpoint of a job whether the capacity is kept or not */
   for( j = 0; j < njobs; ++j )
   {
      SCIP_Longint besttightness;

      assert(starttime <= starttimes[j]);

      /* if we hit the same start time again we skip the loop */
      if( starttime == starttimes[j])
         continue;

      starttime = starttimes[j];
      endtime = -INT_MAX;
      besttightness = 0LL;

      for( i = 0; i < njobs; ++i )
      {
         SCIP_Longint energy;
         SCIP_Longint maxenergy;
         SCIP_Longint tightness;

         assert(endtime <= endtimes[i]);

         /* if we hit the same end time again we skip the loop */
         if( endtime == endtimes[i] )
            continue;

         endtime = endtimes[i];

         /* skip all end times which are smaller than the start time */
         if( endtime <= starttime )
            continue;

         maxenergy = computeMaxEnergy(njobs, durations, demands, releasedates, deadlinedates, starttime, endtime);

         energy = (endtime - starttime) * capacity; /*lint !e647*/
         tightness = maxenergy - energy;

         /* check if the linear constraint is not trivially redundant */
         if( tightness > besttightness )
         {
            besttightness = tightness;

            nrows[i] = removeRedundantRows(rowtightness[i], startidxs[i], nrows[i], tightness);

            /* add row information */
            rowtightness[i][nrows[i]] = tightness;
            startidxs[i][nrows[i]] = j;
            nrows[i]++;
         }
      }
   }

   naddedcons = 0;

   for( j = njobs-1; j >= 0; --j )
   {
      for( i = 0; i < nrows[j]; ++i )
      {
         SCIP_CONS* cons;
         SCIP_Longint energy;
         char name[SCIP_MAXSTRLEN];
         int v;

         starttime = starttimes[startidxs[j][i]];
         endtime = endtimes[j];

         energy = (endtime - starttime) * capacity; /*lint !e647*/

         SCIPdebugMessage("create linear relaxation for time interval [%d,%d] <= %"SCIP_LONGINT_FORMAT" (tightness %"SCIP_LONGINT_FORMAT")\n",
            starttime, endtime, energy, rowtightness[j][i]);

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "relax_%d_start_%d_end_%d", resource, starttime, endtime);
         SCIP_CALL( SCIPcreateConsBasicKnapsack(scip, &cons, name, 0, NULL, NULL, energy) );

         for( v = 0; v < njobs; ++v)
         {
            int overlap;
            int duration;

            duration = durations[v];
            overlap = MIN(endtime - starttime, duration);
            assert(overlap > 0);
            overlap = MIN3(overlap, releasedates[v] + duration - starttime, endtime - deadlinedates[v] + duration);

            /* check for edge-finding idea */
            if( relaxation == 2 &&  releasedates[v] >= starttime && deadlinedates[v] <= endtime )
            {
               assert(duration == overlap);
               SCIP_CALL( SCIPaddCoefKnapsack(scip, cons, vars[v], (SCIP_Longint)(duration * demands[v])) ); /*lint !e647*/
            }
            else if( relaxation == 3 && overlap > 0 )
            {
               assert(overlap <= duration);
               SCIP_CALL( SCIPaddCoefKnapsack(scip, cons, vars[v], (SCIP_Longint)(overlap * demands[v])) ); /*lint !e647*/
            }
         }

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         naddedcons++;
      }
   }

   SCIPinfoMessage(scip, NULL, "add %d constraint as relaxation\n", naddedcons);

   /* free buffers */
   for( j = njobs-1; j >= 0; --j )
   {
      SCIPfreeBufferArray(scip, &startidxs[j]);
      SCIPfreeBufferArray(scip, &rowtightness[j]);
   }
   SCIPfreeBufferArray(scip, &startidxs);
   SCIPfreeBufferArray(scip, &rowtightness);
   SCIPfreeBufferArray(scip, &nrows);

   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &endtimes);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &starttimes);

   return SCIP_OKAY;
}

/** create MIP formulation and CP constraints */
static
SCIP_RETCODE createMipFormulation(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 durations,          /**< durations */
   int**                 demands,            /**< demands */
   int**                 costs,              /**< cost */
   int*                  capacities,         /**< machine capacities */
   int*                  releasedates,       /**< release dates */
   int*                  deadlinedates,      /**< deadline dates */
   int                   njobs,              /**< number of jobs */
   int                   nmachines           /**< number of machines */
   )
{
   SCIP_VAR*** vars;
   SCIP_CONS*** conss;
   SCIP_CONS* cons;

   char name[SCIP_MAXSTRLEN];

   int relaxation;
   int maxtime;
   int i;
   int j;
   int t;

   SCIP_CALL( SCIPgetIntParam(scip, "reading/"READER_NAME"/relaxation", &relaxation) );

   /* compute maximum time */
   maxtime = 0;
   for( j = 0; j < njobs; ++j )
   {
      maxtime = MAX(maxtime, deadlinedates[j]);
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &conss, nmachines) );

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nmachines) );

   /* create master problem  */
   for( i = 0; i < nmachines; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &vars[i], njobs) ); /*lint !e866*/

      for( j = 0; j < njobs; ++j )
      {
         SCIP_VAR* var;
         SCIP_Real objval;
         SCIP_Real ub;

         objval = (SCIP_Real)costs[i][j];

         /* construct variable name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d_machine_%d", j, i);

         if( releasedates[j] > deadlinedates[j] - durations[i][j] || demands[i][j] > capacities[i] )
            ub   = 0.0;
         else
            ub = 1.0;

         /* create binary variable */
         SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, ub, objval, SCIP_VARTYPE_BINARY) );
         vars[i][j] = var;

         /* add and release binary variable */
         SCIP_CALL( SCIPaddVar(scip, var) );
         SCIP_CALL( SCIPreleaseVar(scip, &var) );
      }
   }

   for( j = 0; j < njobs; ++j )
   {
      /* construct constraint name */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d_to_one_machine", j);
      SCIP_CALL( SCIPcreateConsBasicSetpart(scip, &cons, name, 0, NULL) );

      for( i = 0; i < nmachines; ++i )
      {
         SCIP_CALL( SCIPaddCoefSetppc(scip, cons, vars[i][j]) );
      }

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* create for each machines and time point a knapsack constraint */
   for( i = 0; i < nmachines; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &conss[i], maxtime) ); /*lint !e866*/

      for( t = 0; t < maxtime; ++t )
      {
         /* construct constraint name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "machines_%d_time_%d", i, t);

         SCIP_CALL( SCIPcreateConsBasicKnapsack(scip, &cons, name, 0, NULL, NULL, (SCIP_Longint)capacities[i]) );
         SCIP_CALL( SCIPaddCons(scip, cons) );

         conss[i][t] = cons;

         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   /* create for each job/machine/timepoint a binary variables */
   for( j = 0; j < njobs; ++j )
   {
      int est;

      est = releasedates[j];

      for( i = 0; i < nmachines; ++i )
      {
         int lst;

         lst = deadlinedates[j] - durations[i][j];

         /* construct constraint name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d_machine_%d", j, i);

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, vars[i][j], -1.0) );

         for( t = est; t <= lst; ++t )
         {
            SCIP_VAR* var;
            int h;

            /* construct variable name */
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d_machine_%d_time_%d", j, i, t);

            SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );

            SCIP_CALL( SCIPaddVar(scip, var) );

            /* add variable to linear constraint */
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, 1.0) );

            for( h = t ; h < MIN(t + durations[i][j], maxtime); ++h )
            {
               SCIP_CALL( SCIPaddCoefKnapsack(scip, conss[i][h], var, (SCIP_Longint)demands[i][j] ) );
            }

            SCIP_CALL( SCIPreleaseVar(scip, &var) );
         }

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   /* free constraint buffers */
   for( i = nmachines - 1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &conss[i]);
   }
   SCIPfreeBufferArray(scip, &conss);

   if( relaxation == 1 )
   {
      /* create for each optimal resource a singe constraint relaxation */
      for( i = 0; i < nmachines; ++i )
      {
         SCIP_Longint capacity;
         int est;
         int lct;

         est = INT_MAX;
         lct = INT_MIN;

         /* compute earliest start end latest completion time */
         for( j = 0; j < njobs; ++j )
         {
            if( demands[i][j] > 0 )
            {
               est = MIN(est, releasedates[j]);
               lct = MAX(lct, deadlinedates[j]);
            }
         }

         /* construct constraint name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "relax_%d", i);

         capacity = capacities[i] * (lct - est); /*lint !e647*/

         SCIP_CALL( SCIPcreateConsBasicKnapsack(scip, &cons, name, 0, NULL, NULL, capacity) );

         for( j = 0; j < njobs; ++j )
         {
            if( demands[i][j] > 0 )
            {
               SCIP_CALL( SCIPaddCoefKnapsack(scip, cons, vars[i][j], (SCIP_Longint)(durations[i][j] * demands[i][j]) ) ); /*lint !e647*/
            }
         }

         /* add the constraint to the problem and release it (locally) */
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }
   else if( relaxation >= 2 )
   {
      /* create for each optimal resource a singe constraint relaxation */
      for( i = 0; i < nmachines; ++i )
      {
         SCIP_CALL( createIntervalRelaxation(scip, relaxation, i, vars[i], durations[i], demands[i], capacities[i], releasedates, deadlinedates, njobs) );
      }
   }

   for( i = nmachines-1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &vars[i]);
   }
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** create MIP formulation */
static
SCIP_RETCODE createMipCpFormulation(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 durations,          /**< durations */
   int**                 demands,            /**< demands */
   int**                 costs,              /**< cost */
   int*                  capacities,         /**< machine capacities */
   int*                  releasedates,       /**< release dates */
   int*                  deadlinedates,      /**< deadline dates */
   int                   njobs,              /**< number of jobs */
   int                   nmachines           /**< number of machines */
   )
{
   SCIP_CONS*** conss;
   SCIP_CONS* cons;
   SCIP_VAR*** binvars;
   SCIP_VAR*** vars;
   SCIP_VAR* var;
   int* nnvars;
   int** localdemands;
   int** localdurations;

   char name[SCIP_MAXSTRLEN];

   int maxtime;
   int i;
   int j;
   int t;

   /* compute maximum time */
   maxtime = 0;
   for( j = 0; j < njobs; ++j )
   {
      maxtime = MAX(maxtime, deadlinedates[j]);
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &conss, nmachines) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, nmachines) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nmachines) );
   SCIP_CALL( SCIPallocBufferArray(scip, &localdemands, nmachines) );
   SCIP_CALL( SCIPallocBufferArray(scip, &localdurations, nmachines) );

   SCIP_CALL( SCIPallocBufferArray(scip, &nnvars, nmachines) );
   BMSclearMemoryArray(nnvars, nmachines);

   /* create for each machines and time point a knapsack constraint */
   for( i = 0; i < nmachines; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &conss[i], maxtime) ); /*lint !e866*/

      for( t = 0; t < maxtime; ++t )
      {
         /* construct constraint name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "machines_%d_time_%d", i, t);

         SCIP_CALL( SCIPcreateConsKnapsack(scip, &cons, name, 0, NULL, NULL, (SCIP_Longint)capacities[i],
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );

         conss[i][t] = cons;

         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &binvars[i], njobs) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &vars[i], njobs) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &localdemands[i], njobs) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &localdurations[i], njobs) ); /*lint !e866*/
   }

   for( j = 0; j < njobs; ++j )
   {
      SCIP_CONS* maccons;
      int est;

      est = releasedates[j];

      /* construct constraint name */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d", j);

      SCIP_CALL( SCIPcreateConsSetpart(scip, &cons, name, 0, NULL,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      for( i = 0; i < nmachines; ++i )
      {
         SCIP_CONS* lincons;
         int lst;
         int idx;

         lst = deadlinedates[j] - durations[i][j];

         if( est > lst )
            continue;

         idx = nnvars[i];
         nnvars[i]++;

         /* construct constraint name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d_machine_%d", j, i);

         SCIP_CALL( SCIPcreateConsSetpart(scip, &maccons, name, 0, NULL,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         /* construct variable name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d_choose_%d", j, i);

         SCIP_CALL( SCIPcreateVar(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
               TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, var) );
         binvars[i][idx] = var;
         SCIP_CALL( SCIPreleaseVar(scip, &var) );

         SCIP_CALL( SCIPgetNegatedVar(scip, binvars[i][idx], &var) );

         SCIP_CALL( SCIPaddCoefSetppc(scip, maccons, var) );

         /* add variable to set partitioning constraint */
         SCIP_CALL( SCIPaddCoefSetppc(scip, cons, binvars[i][idx]) );

         /* construct variable name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d_starts_%d", idx, i);

         SCIP_CALL( SCIPcreateVar(scip, &var, name, 0.0, (SCIP_Real)lst, 0.0, SCIP_VARTYPE_IMPLINT,
               TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

         SCIP_CALL( SCIPaddVar(scip, var) );
         vars[i][idx] = var;

         SCIP_CALL( SCIPreleaseVar(scip, &var) );

         /* construct constraint name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d_machine_%d", idx, i);

         SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, 0, NULL, NULL, 0.0, 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCoefLinear(scip, lincons, vars[i][idx], -1.0) );

         localdemands[i][idx] = demands[i][j];
         localdurations[i][idx] = durations[i][j];

         for( t = est; t <= lst; ++t )
         {
            int h;

            /* construct variable name */
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d_machine_%d_time_%d", j, i, t);

            SCIP_CALL( SCIPcreateVar(scip, &var, name, 0.0, 1.0, (SCIP_Real)costs[i][j], SCIP_VARTYPE_BINARY,
                  TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

            SCIP_CALL( SCIPaddVar(scip, var) );

            /* add variable to linear linking  constraint */
            SCIP_CALL( SCIPaddCoefLinear(scip, lincons, var, (SCIP_Real)t) );

            /* add variable to linear linking  constraint */
            SCIP_CALL( SCIPaddCoefSetppc(scip, maccons, var) );

            for( h = t ; h < MIN(t + durations[i][j], maxtime); ++h )
            {
               SCIP_CALL( SCIPaddCoefKnapsack(scip, conss[i][h], var, (SCIP_Longint)demands[i][j] ) );
            }

            SCIP_CALL( SCIPreleaseVar(scip, &var) );
         }

         SCIP_CALL( SCIPaddCons(scip, lincons) );
         SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

         SCIP_CALL( SCIPaddCons(scip, maccons) );
         SCIP_CALL( SCIPreleaseCons(scip, &maccons) );
      }


      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* create CP constraints */

   for( i = 0; i < nmachines; ++i )
   {
      /* construct constraint name */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "machine_%d", i);

      /* create machine choice constraint */
      SCIP_CALL( SCIPcreateConsOptcumulative(scip, &cons, name, nnvars[i], vars[i], binvars[i], localdurations[i], localdemands[i], capacities[i],
            FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }


   for( i = nmachines - 1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &localdurations[i]);
      SCIPfreeBufferArray(scip, &localdemands[i]);
      SCIPfreeBufferArray(scip, &vars[i]);
      SCIPfreeBufferArray(scip, &binvars[i]);
      SCIPfreeBufferArray(scip, &conss[i]);
   }

   SCIPfreeBufferArray(scip, &localdemands);
   SCIPfreeBufferArray(scip, &localdurations);
   SCIPfreeBufferArray(scip, &nnvars);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &binvars);
   SCIPfreeBufferArray(scip, &conss);

   return SCIP_OKAY;
}

/** create CIP formulation */
static
SCIP_RETCODE createCipFormulation(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 durations,          /**< durations */
   int**                 demands,            /**< demands */
   int**                 costs,              /**< cost */
   int*                  capacities,         /**< machine capacities */
   int*                  releasedates,       /**< release dates */
   int*                  deadlinedates,      /**< deadline dates */
   int                   njobs,              /**< number of jobs */
   int                   nmachines           /**< number of machines */
   )
{
   SCIP_CONS** conss;
   SCIP_CONS* cons;
   SCIP_VAR*** vars;
   SCIP_VAR*** binvars;
   SCIP_VAR* var;
   int* machines;

   char name[SCIP_MAXSTRLEN];

   SCIP_Real objval;
   SCIP_Bool dualreduction;

   int relaxation;
   int i;
   int j;

   SCIP_CALL( SCIPgetIntParam(scip, "reading/"READER_NAME"/relaxation", &relaxation) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/"READER_NAME"/dualreduction", &dualreduction) );

   SCIP_CALL( SCIPallocBufferArray(scip, &conss, njobs) );

   /* create for each job a set partitioning constraint */
   for( j = 0; j < njobs; ++j )
   {
      /* construct constraint name */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d", j);

      SCIP_CALL( SCIPcreateConsBasicSetpart(scip, &cons, name, 0, NULL) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      conss[j] = cons;
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, nmachines) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nmachines) );
   SCIP_CALL( SCIPallocBufferArray(scip, &machines, nmachines) );

   for( i = 0; i < nmachines; ++i )
   {
      int nvars;

      nvars = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &binvars[i], njobs) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &vars[i], njobs) ); /*lint !e866*/

      BMSclearMemoryArray(binvars[i], njobs); /*lint !e866*/
      BMSclearMemoryArray(vars[i], njobs); /*lint !e866*/

      for( j = 0; j < njobs; ++j )
      {
         /* check if job is scheduleable on that machine */
         if(releasedates[j] + durations[i][j] > deadlinedates[j] || demands[i][j] > capacities[i] )
            continue;

         /* construct variable name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d_choose_%d", j, i);

         SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, (SCIP_Real)costs[i][j], SCIP_VARTYPE_BINARY) );
         SCIP_CALL( SCIPaddVar(scip, var) );
         binvars[i][nvars] = var;
         SCIP_CALL( SCIPreleaseVar(scip, &var) );

         /* construct variable name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job_%d_starts_%d", j, i);

         SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, (SCIP_Real)releasedates[j], (SCIP_Real)(deadlinedates[j] - durations[i][j]), 0.0, SCIP_VARTYPE_INTEGER) );
         SCIP_CALL( SCIPaddVar(scip, var) );
         vars[i][nvars] = var;
         SCIP_CALL( SCIPreleaseVar(scip, &var) );

         if( !dualreduction )
         {
            SCIP_CALL( SCIPaddVarLocksType(scip, binvars[i][nvars], SCIP_LOCKTYPE_MODEL, 1, 1) );
            SCIP_CALL( SCIPaddVarLocksType(scip, vars[i][nvars], SCIP_LOCKTYPE_MODEL, 1, 1) );
         }

         /* add choice variable to set partitioning constraint */
         SCIP_CALL( SCIPaddCoefSetppc(scip, conss[j], binvars[i][nvars]) );

         demands[i][nvars] = demands[i][j];
         durations[i][nvars] = durations[i][j];

         nvars++;
      }

      machines[i] = nvars;

      if( nvars > 0 )
      {
         SCIP_Bool initial;
         SCIP_Bool separate;

         if( relaxation == 0 )
         {
            initial = FALSE;
            separate = FALSE;
         }
         else
         {
            initial = TRUE;
            separate = TRUE;
         }

         /* construct constraint name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "machine_%d", i);

         /* create machine choice constraint */
         SCIP_CALL( SCIPcreateConsOptcumulative(scip, &cons, name, nvars, vars[i], binvars[i], durations[i], demands[i], capacities[i],
               initial, separate, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   SCIP_CALL( SCIPinitHeurOptcumulative(scip, nmachines, njobs, machines, binvars, vars, durations, demands, capacities) );

   /* check for a given total leftover */
   SCIP_CALL( findBestObjectiveValue(scip, &objval) );

   SCIP_CALL( SCIPsetObjlimit(scip, objval) );

   /* free all buffers */
   for( i = 0; i < nmachines; ++i )
   {
      SCIPfreeBufferArray(scip, &vars[i]);
      SCIPfreeBufferArray(scip, &binvars[i]);
   }
   SCIPfreeBufferArray(scip, &machines );
   SCIPfreeBufferArray(scip, &vars );
   SCIPfreeBufferArray(scip, &binvars);
   SCIPfreeBufferArray(scip, &conss );

   return SCIP_OKAY;
}

/** read given file and create corresponding problem */
static
SCIP_RETCODE readFile(
   SCIP*                 scip,               /**< SCIP data structure */
   CMININPUT*            cmininput,          /**< CMIN reading data */
   const char*           filename            /**< file name */
   )
{
   int** durations;
   int** demands;
   int** costs;
   int* capacities;
   int* releasedates;
   int* deadlinedates;

   int njobs;
   int nmachines;

   int i;
   int j;

   /* read the number of jobs */
   if( !getNextToken(cmininput) )
   {
      syntaxError(scip, cmininput, "missing number if jobs\n");
      return SCIP_OKAY;
   }

   njobs = atoi(cmininput->token);

   /* read the number of machines */
   if( !getNextToken(cmininput) )
   {
      syntaxError(scip, cmininput, "missing number of machines");
      return SCIP_OKAY;
   }

   nmachines = atoi(cmininput->token);

   SCIPdebugMessage("njobs <%d> nmachines <%d>\n", njobs, nmachines);

   /* allocate all neccessary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &durations, nmachines) );
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, nmachines) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costs, nmachines) );
   for( i = 0; i < nmachines; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &durations[i], njobs) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &demands[i], njobs) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &costs[i], njobs) ); /*lint !e866*/
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &capacities, nmachines) );
   SCIP_CALL( SCIPallocBufferArray(scip, &releasedates, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &deadlinedates, njobs) );

   for( i = 0; i < nmachines && !cmininput->haserror; ++i )
   {
      for( j = 0; j < njobs; ++j )
      {
         /* get job duration */
         if( !getNextToken(cmininput) )
         {
            syntaxError(scip, cmininput, "missing job duration");
            break;
         }
         assert(cmininput->haserror == FALSE);

         durations[i][j] = atoi(cmininput->token);

         /* get demand  */
         if( !getNextToken(cmininput) )
         {
            syntaxError(scip, cmininput, "missing job demand\n");
            break;
         }
         assert(cmininput->haserror == FALSE);

         demands[i][j] = atoi(cmininput->token);

         /* get job cost */
         if( !getNextToken(cmininput) )
         {
            syntaxError(scip, cmininput, "missing job cost\n");
            break;
         }
         assert(cmininput->haserror == FALSE);

         costs[i][j] = atoi(cmininput->token);

         SCIPdebugMessage("job %2d: duration %d, demand %d, cost %d\n", j, durations[i][j], demands[i][j], costs[i][j]);
      }
   }

   /* parse the machine capacities */
   for( i = 0; i < nmachines && !cmininput->haserror; ++i)
   {
      /* get capacity */
      if( !getNextToken(cmininput) )
      {
         syntaxError(scip, cmininput, "missing machine capacity\n");
         break;
      }

      capacities[i] = atoi(cmininput->token);

      SCIPdebugMessage("capaciy of machine %d is %d\n", i, capacities[i]);
   }

   /* get release and deadline dates */
   for( j = 0; j < njobs && !cmininput->haserror; ++j )
   {

      /* get release date */
      if( !getNextToken(cmininput) )
      {
         syntaxError(scip, cmininput, "missing release date\n");
         break;
      }

      releasedates[j] = atoi(cmininput->token);

      /* get deadline data */
      if( !getNextToken(cmininput) )
      {
         syntaxError(scip, cmininput, "missing deadline date\n");
         break;
      }
      deadlinedates[j] = atoi(cmininput->token);

      SCIPdebugMessage("job %2d: [%d,%d]\n", j, releasedates[j], deadlinedates[j]);
   }

   if( !cmininput->haserror )
   {
      SCIP_Bool mip;
      SCIP_Bool cip;
      char* probname;
      char* str;

      SCIPdebug( displayInputData(scip, durations, demands, costs, capacities, releasedates, deadlinedates, njobs, nmachines) );

      SCIP_CALL( SCIPgetBoolParam(scip, "reading/"READER_NAME"/mip", &mip) );
      SCIP_CALL( SCIPgetBoolParam(scip, "reading/"READER_NAME"/cip", &cip) );

      /* construct problem name */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &str, filename, strlen(filename)+1) );

      SCIPsplitFilename(str, NULL, &probname, NULL, NULL);

      /* initialize problem data */
      SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

      SCIPfreeBufferArray(scip, &str);

      if( mip && cip )
      {
         SCIP_CALL( createMipCpFormulation(scip, durations, demands, costs, capacities, releasedates, deadlinedates, njobs, nmachines) );
      }
      else if( mip )
      {
         SCIP_CALL( createMipFormulation(scip, durations, demands, costs, capacities, releasedates, deadlinedates, njobs, nmachines) );
      }
      else if( cip )
      {
         SCIP_CALL( createCipFormulation(scip, durations, demands, costs, capacities, releasedates, deadlinedates, njobs, nmachines) );
      }
   }

   /* free all buffers in the reverse order since the buffers are handled by a stack */
   SCIPfreeBufferArray(scip, &deadlinedates);
   SCIPfreeBufferArray(scip, &releasedates);
   SCIPfreeBufferArray(scip, &capacities);
   for( i = nmachines - 1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &costs[i]);
      SCIPfreeBufferArray(scip, &demands[i]);
      SCIPfreeBufferArray(scip, &durations[i]);
   }
   SCIPfreeBufferArray(scip, &costs);
   SCIPfreeBufferArray(scip, &demands);
   SCIPfreeBufferArray(scip, &durations);

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCmin)
{  /*lint --e{715}*/

   CMININPUT cmininput;
   SCIP_FILE* file;

   /* try to open the file */
   if( NULL == (file = SCIPfopen(filename, "r")) )
   {
      perror(filename);
      return SCIP_NOFILE;
   }

   cmininput.linebuf[0] = '\0';
   cmininput.linenumber = 0;
   cmininput.linepos = 0;
   cmininput.haserror = FALSE;
   cmininput.file = file;

   SCIP_CALL( SCIPallocBufferArray(scip, &cmininput.token, SCIP_MAXSTRLEN) );
   cmininput.token[0] = '\0';

   SCIP_CALL( readFile(scip, &cmininput, filename) );

   SCIPfreeBufferArray(scip, &cmininput.token);

   /* close file */
   (void) SCIPfclose(file);

   if( cmininput.haserror )
      return SCIP_READERROR;

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the cp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCmin(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include cp reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadCmin) );

   SCIP_CALL( SCIPaddStringParam(scip,
         "reading/"READER_NAME"/filename", "name of the file including best known solutions",
         NULL, FALSE, DEFAULT_FILENAME, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/dualreduction", "add locks to avoid dual reductions?",
         NULL, FALSE, DEFAULT_DUALREDUCTION, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/mip", "create a MIP formulation?",
         NULL, FALSE, DEFAULT_MIP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/initial", "should model constraints be in initial LP?",
         NULL, TRUE, DEFAULT_INITIAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/cip", "create a CIP formulation?",
         NULL, FALSE, DEFAULT_CIP, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "reading/"READER_NAME"/relaxation", "which relaxation should be added to the maseter (0: none; 1: single; 2: edge-finding; 3: energetic-reasoning",
         NULL, FALSE, DEFAULT_RELAXATION, 0, 3, NULL, NULL) );

   return SCIP_OKAY;
}
