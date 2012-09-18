/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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
#define DEFAULT_DYNAMICCONSS     FALSE  /**< should model constraints be subject to aging? */
#define DEFAULT_DUALREDUCTION     TRUE  /**< add locks to avoid dual reductions */
#define DEFAULT_BRANCHPRIORITY       1  /**< branching priority for the binary choice variables */
#define DEFAULT_MIP              FALSE  /**< create a MIP formulation */
#define DEFAULT_INITIAL           TRUE  /**< should model constraints be in initial LP? */
#define DEFAULT_CIP               TRUE  /**< create a CIP formulation */

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

   if( SCIPfgets(cmininput->linebuf, sizeof(cmininput->linebuf), cmininput->file) == NULL )
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
   while( SCIPfgets(buffer, sizeof(buffer), file) != NULL )
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
   SCIP_CONS*** conss;
   SCIP_CONS* cons;
   SCIP_VAR* var;

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

   /* create for each machines and time point a knapsack constraint */
   for( i = 0; i < nmachines; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &conss[i], maxtime) );

      for( t = 0; t < maxtime; ++t )
      {
         /* construct constraint name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "machines(%d,%d)", i, t);

         SCIP_CALL( SCIPcreateConsKnapsack(scip, &cons, name, 0, NULL, NULL, (SCIP_Longint)capacities[i],
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
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

      /* construct constraint name */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job%d", j);

      SCIP_CALL( SCIPcreateConsSetpart(scip, &cons, name, 0, NULL,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      for( i = 0; i < nmachines; ++i )
      {
         int lst;

         lst = deadlinedates[j] - durations[i][j];

         for( t = est; t <= lst; ++t )
         {
            int h;

            /* construct variable name */
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%d,%d,%d", j, i, t);

            SCIP_CALL( SCIPcreateVar(scip, &var, name, 0.0, 1.0, (SCIP_Real)costs[i][j], SCIP_VARTYPE_BINARY,
                  TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

            SCIP_CALL( SCIPaddVar(scip, var) );

            /* add variable to set partitioning constraint */
            SCIP_CALL( SCIPaddCoefSetppc(scip, cons, var) );

            for( h = t ; h < MIN(t + durations[i][j], maxtime); ++h )
            {
               SCIP_CALL( SCIPaddCoefKnapsack(scip, conss[i][h], var, (SCIP_Longint)demands[i][j] ) );
            }

            SCIP_CALL( SCIPreleaseVar(scip, &var) );
         }
      }

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   for( i = nmachines - 1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &conss[i]);
   }
   SCIPfreeBufferArray(scip, &conss);

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
      SCIP_CALL( SCIPallocBufferArray(scip, &conss[i], maxtime) );

      for( t = 0; t < maxtime; ++t )
      {
         /* construct constraint name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "machines(%d,%d)", i, t);

         SCIP_CALL( SCIPcreateConsKnapsack(scip, &cons, name, 0, NULL, NULL, (SCIP_Longint)capacities[i],
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );

         conss[i][t] = cons;

         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &binvars[i], njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vars[i], njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &localdemands[i], njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &localdurations[i], njobs) );
   }

   for( j = 0; j < njobs; ++j )
   {
      SCIP_CONS* maccons;
      int est;

      est = releasedates[j];

      /* construct constraint name */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job%d", j);

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
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job(%d,%d)", j, i);

         SCIP_CALL( SCIPcreateConsSetpart(scip, &maccons, name, 0, NULL,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         /* construct variable name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job%d_choose%d", j, i);

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
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job%d_starts%d", idx, i);

         SCIP_CALL( SCIPcreateVar(scip, &var, name, 0.0, (SCIP_Real)lst, 0.0, SCIP_VARTYPE_IMPLINT,
               TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

         SCIP_CALL( SCIPaddVar(scip, var) );
         vars[i][idx] = var;

         SCIP_CALL( SCIPreleaseVar(scip, &var) );

         /* construct constraint name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job(%d,%d)", idx, i);

         SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, 0, NULL, NULL, 0.0, 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCoefLinear(scip, lincons, vars[i][idx], -1.0) );

         localdemands[i][idx] = demands[i][j];
         localdurations[i][idx] = durations[i][j];

         for( t = est; t <= lst; ++t )
         {
            int h;

            /* construct variable name */
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%d,%d,%d", j, i, t);

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
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "machine%d", i);

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
   SCIP_Bool dynamicconss;
   SCIP_VAR** vars;
   SCIP_VAR*** binvars;
   SCIP_VAR* var;

   int* localdemands;
   int* localdurations;
   int branchpriority;

   char name[SCIP_MAXSTRLEN];

   SCIP_Real objval;
   SCIP_Bool dualreduction;

   int i;
   int j;

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/"READER_NAME"/dynamicconss", &dynamicconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/"READER_NAME"/dualreduction", &dualreduction) );
   SCIP_CALL( SCIPgetIntParam(scip, "reading/"READER_NAME"/branchpriority", &branchpriority) );

   SCIP_CALL( SCIPallocBufferArray(scip, &conss, njobs) );

   /* create for each job a set partitioning constraint */
   for( j = 0; j < njobs; ++j )
   {
      /* construct constraint name */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job%d", j);

      SCIP_CALL( SCIPcreateConsSetpart(scip, &cons, name, 0, NULL,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      conss[j] = cons;
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, nmachines) );
   SCIP_CALL( SCIPallocBufferArray(scip, &localdemands, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &localdurations, njobs) );

   for( i = 0; i < nmachines; ++i )
   {
      int nvars;

      nvars = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &binvars[i], njobs) );

      BMSclearMemoryArray(binvars[i], njobs);

      for( j = 0; j < njobs; ++j )
      {
         /* check if job is scheduleable on that machine */
         if(releasedates[j] + durations[i][j] > deadlinedates[j] || demands[i][j] > capacities[i] )
            continue;

         /* construct variable name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job%d_choose%d", j, i);

         SCIP_CALL( SCIPcreateVar(scip, &var, name, 0.0, 1.0, (SCIP_Real)costs[i][j], SCIP_VARTYPE_BINARY,
               TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, var) );
         SCIP_CALL( SCIPchgVarBranchPriority(scip, var, branchpriority) );
         binvars[i][nvars] = var;
         SCIP_CALL( SCIPreleaseVar(scip, &var) );

         /* construct variable name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "job%d_starts%d", j, i);

         SCIP_CALL( SCIPcreateVar(scip, &var, name, (SCIP_Real)releasedates[j], (SCIP_Real)(deadlinedates[j] - durations[i][j]), 0.0, SCIP_VARTYPE_INTEGER,
               TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

         SCIP_CALL( SCIPaddVar(scip, var) );
         vars[nvars] = var;

         if( !dualreduction )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, binvars[i][nvars], 1, 1) );
            SCIP_CALL( SCIPaddVarLocks(scip, vars[nvars], 1, 1) );
         }

         SCIP_CALL( SCIPreleaseVar(scip, &var) );

         /* add choice variable to set partitioning constraint */
         SCIP_CALL( SCIPaddCoefSetppc(scip, conss[j], binvars[i][nvars]) );

         localdemands[nvars] = demands[i][j];
         localdurations[nvars] = durations[i][j];

         nvars++;
      }


      if( nvars > 0 )
      {
         /* construct constraint name */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "machine%d", i);

         /* create machine choice constraint */
         SCIP_CALL( SCIPcreateConsOptcumulative(scip, &cons, name, nvars, vars, binvars[i], localdurations, localdemands, capacities[i],
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   SCIP_CALL( SCIPinitHeurOptcumulative(scip, binvars, njobs, nmachines) );

   /* check for a given total leftover */
   SCIP_CALL( findBestObjectiveValue(scip, &objval) );

   SCIP_CALL( SCIPsetObjlimit(scip, objval) );

   /* free all buffers */
   SCIPfreeBufferArray(scip, &localdurations);
   SCIPfreeBufferArray(scip, &localdemands);
   for( i = 0; i < nmachines; ++i )
   {
      SCIPfreeBufferArray(scip, &binvars[i]);
   }
   SCIPfreeBufferArray(scip, &binvars);
   SCIPfreeBufferArray(scip, &vars );
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
      SCIP_CALL( SCIPallocBufferArray(scip, &durations[i], njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &demands[i], njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &costs[i], njobs) );
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
            syntaxError(scip, cmininput, "missing job duration\n");
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
      SCIP_CALL( SCIPcreateProb(scip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

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
         "reading/"READER_NAME"/dynamicconss", "should model constraints be subject to aging?",
         NULL, FALSE, DEFAULT_DYNAMICCONSS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/dualreduction", "add locks to avoid dual reductions?",
         NULL, FALSE, DEFAULT_DUALREDUCTION, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "reading/"READER_NAME"/branchpriority", "branching priority for the binary choice variables",
         NULL, FALSE, DEFAULT_BRANCHPRIORITY, -INT_MAX/4, INT_MAX/4, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/mip", "create a MIP formulation?",
         NULL, FALSE, DEFAULT_MIP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/initial", "should model constraints be in initial LP?",
         NULL, TRUE, DEFAULT_INITIAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/cip", "create a CIP formulation?",
         NULL, TRUE, DEFAULT_CIP, NULL, NULL) );

   return SCIP_OKAY;
}
