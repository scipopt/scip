/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_sm.c
 * @brief  scheduling problem file reader for RCPSP format
 * @author Michael Bastubbe
 * @author Stefan Heinz
 *
 * This reader is capabale of parsing resource-constrained project scheduling problem (RCPSP) instances. The <a
 * href="http://129.187.106.231/psplib/datasm.html">PSPlib</a> provides several instances set.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>


#include "heur_listscheduling.h"
#include "reader_sm.h"

#include "scip/cons_cumulative.h"
#include "scip/cons_linear.h"
#include "scip/cons_varbound.h"

#define READER_NAME             "smreader"
#define READER_DESC             "scheduling file reader for sm files (RCPSP format)"
#define READER_EXTENSION        "sm"


#define SM_MAX_LINELEN      65536     /**< size of the line buffer for reading or writing */

enum reading_states {
   ERROR = 0,
   NEXT,
   NJOBS,
   JOBS,
   NRESOURCES,
   RESOURCENAMES,
   RESOURCECAPACITIES,
   PRECEDENCES,
   END
};
typedef enum reading_states STATE;


/** data structure for resources constrained project scheduling problems */
struct SCIP_RcpspData
{
   SCIP_DIGRAPH*        precedencegraph;     /**< precedence graph of the jobs */
   const char**         jobnames;            /**< array of job names */
   const char**         resourcenames;       /**< array of resource names */
   int**                demands;             /**< resource demands matrix (job i needs demands[i][j] units of resource j) */
   int*                 durations;           /**< array of job durations */
   int*                 capacities;          /**< array of resource capacities */
   int                  njobs;               /**< number of jobs */
   int                  nresources;          /**< number of resources */
};
typedef struct SCIP_RcpspData SCIP_RCPSPDATA;


/*
 * Local methods
 */

#ifdef SCIP_DEBUG
/* print the resource constrained project scheduling data */
static
void outputRcpspData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RCPSPDATA*       rcpspdata           /**< pointer to resources constrained project scheduling data */
   )
{
   int i;
   int r;

   /* output jobs */
   SCIPinfoMessage(scip, NULL, "Number of jobs: %d\n", rcpspdata->njobs);
   for(i = 0; i < rcpspdata->njobs ; ++i )
   {
      SCIPinfoMessage(scip, NULL, "job: %-10s \n", rcpspdata->jobnames[i]);
      SCIPinfoMessage(scip, NULL, "  duration: %3d \n", rcpspdata->durations[i] );
      SCIPinfoMessage(scip, NULL, "  resource profile: ");

      for( r = 0; r < rcpspdata->nresources; ++r )
      {
         SCIPinfoMessage(scip, NULL, " %3d" , rcpspdata->demands[i][r] );
      }

      SCIPinfoMessage(scip, NULL, "\n");
   }
   SCIPinfoMessage(scip, NULL, "\n");

   /* output resource capacities */
   SCIPinfoMessage(scip, NULL, "Resource capacities: ");
   for( r = 0; r < rcpspdata->nresources; ++r )
   {
      if( r == 0 )
      {
         SCIPinfoMessage(scip, NULL, " %d " , rcpspdata->capacities[r] );
      }
      else
      {
         SCIPinfoMessage(scip, NULL, ", %d " , rcpspdata->capacities[r] );
      }
   }
   SCIPinfoMessage(scip, NULL, "\n");

   /* print precedence graph */
   SCIPinfoMessage(scip, NULL, "Precedences:\n");
   SCIPdigraphPrint(rcpspdata->precedencegraph, SCIPgetMessagehdlr(scip), NULL);
}
#endif

/** print error message */
static
void parseError(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   lineno,             /**< current line number of input file */
   const char*           msg,                /**< error message to display */
   const char*           erritem,            /**< token where the error occured, or NULL */
   STATE*                state               /**< pointer to current reading state */
   )
{
   assert(msg != NULL);
   assert(state != NULL);

   if( erritem != NULL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Line %d: %s <%s>\n", lineno, msg, erritem);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Line %d: %s\n", lineno, msg);
   }

   *state = ERROR;
}

/** check if we reached a section */
static
void checkForNewSection(
   char*                 linestr,            /**< current line */
   STATE*                state               /**< pointer to current reading state */
   )
{
   assert(linestr != NULL);
   assert(state != NULL);

   if( strncmp(linestr, "jobs", 4) == 0 )
      *state = NJOBS;
   else if( strncmp(linestr, "RESOURCES", 9) == 0 )
      *state = NRESOURCES;
   else if( strncmp(linestr, "PRECEDENCE", 4) == 0 )
      *state = PRECEDENCES;
   else if( strncmp(linestr, "REQUESTS", 4) == 0 )
      *state = JOBS;
   else if( strncmp(linestr, "RESOURCEAVAILABILITIES", 10) == 0 )
      *state = RESOURCENAMES;
}

/** parese number of resources */
static
SCIP_RETCODE getNResources(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   lineno,             /**< current line number of input file */
   char*                 linestr,            /**< current line */
   STATE*                state,              /**< pointer to current reading state */
   SCIP_RCPSPDATA*       rcpspdata           /**< pointer to resources constrained project scheduling data */
   )
{
   SCIP_Real nresources;
   char* endptr;
   char* number;

   assert(linestr != NULL);
   assert(state != NULL);

   if( strncmp(linestr, "RESOURCES", 4) == 0 )
      return SCIP_OKAY;

   /* truncate the line via ':' and ignore the first part */
   (void)SCIPstrtok(linestr, ":", &endptr);
   number = SCIPstrtok(NULL, ":", &endptr);

   if( !SCIPstrToRealValue(number, &nresources, &endptr) )
   {
      parseError(scip, lineno, "expexted number of resources", linestr, state);
      return SCIP_OKAY;
   }

   rcpspdata->nresources = (int)(nresources + 0.5);

   SCIP_CALL( SCIPallocBufferArray(scip, &rcpspdata->capacities, nresources) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rcpspdata->resourcenames, nresources) );

   *state = NEXT;

   return SCIP_OKAY;
}

/** parse number of jobs */
static
SCIP_RETCODE getNJobs(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   lineno,             /**< current line number of input file */
   char*                 linestr,            /**< current line */
   STATE*                state,              /**< pointer to current reading state */
   SCIP_RCPSPDATA*       rcpspdata           /**< pointer to resources constrained project scheduling data */
   )
{
   SCIP_Real njobs;
   char* endptr;
   char* number;

   assert(linestr != NULL);
   assert(state != NULL);

   /* truncate the line via ':' and ignore the first part */
   (void)SCIPstrtok(linestr, ":", &endptr);
   number = SCIPstrtok(NULL, ":", &endptr);

   if( !SCIPstrToRealValue(number, &njobs, &endptr) )
   {
      parseError(scip, lineno, "expexted number of jobs", linestr, state);
      return SCIP_OKAY;
   }

   rcpspdata->njobs = (int)(njobs + 0.5);

   SCIP_CALL( SCIPallocBufferArray(scip, &rcpspdata->jobnames, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rcpspdata->durations, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rcpspdata->demands, njobs) );

   *state = NEXT;

   return SCIP_OKAY;
}

/** pares resource capacities */
static
SCIP_RETCODE getResourcesNames(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   lineno,             /**< current line number of input file */
   char*                 linestr,            /**< current line */
   STATE*                state,              /**< pointer to current reading state */
   SCIP_RCPSPDATA*       rcpspdata           /**< pointer to resources constrained project scheduling data */
   )
{
   char* name;
   char* endptr;
   int r;

   assert(linestr != NULL);
   assert(state != NULL);

   if( strncmp(linestr, "RESOURCEAVAILABILITIES", 10) == 0 )
      return SCIP_OKAY;

   /* pares resource names */
   name = SCIPstrtok(linestr, "R", &endptr);
   r = 0;

   do
   {
      while(isspace(*name))
         name++;

      SCIP_CALL( SCIPduplicateBufferArray(scip, &rcpspdata->resourcenames[r], name, strlen(name) + 1) );
      r++;
   }
   while( (name = SCIPstrtok(NULL, "R", &endptr)) != NULL );

   *state = RESOURCECAPACITIES;

   return SCIP_OKAY;
}

/** parse resource capacities */
static
SCIP_RETCODE getResourcesCapacities(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   lineno,             /**< current line number of input file */
   char*                 linestr,            /**< current line */
   STATE*                state,              /**< pointer to current reading state */
   SCIP_RCPSPDATA*       rcpspdata           /**< pointer to resources constrained project scheduling data */
   )
{
   SCIP_Real value;
   int r;

   assert(linestr != NULL);
   assert(state != NULL);

   /* parse resources capacities */
   for( r = 0; r < rcpspdata->nresources; ++r )
   {
      if( SCIPstrToRealValue(linestr, &value, &linestr) )
         rcpspdata->capacities[r] = (int)(value + 0.5);
   }

   *state = END;

   return SCIP_OKAY;
}

/** parese job informations */
static
SCIP_RETCODE getJobs(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   lineno,             /**< current line number of input file */
   char*                 linestr,            /**< current line */
   STATE*                state,              /**< pointer to current reading state */
   SCIP_RCPSPDATA*       rcpspdata           /**< pointer to resources constrained project scheduling data */
   )
{
   char jobname[SCIP_MAXSTRLEN];
   int value;
   int jobid;
   int r;

   assert(linestr != NULL);
   assert(state != NULL);

   /* skip lines which are not of interest */
   if ( (!strncmp(linestr, "REQUESTS", 4) ) || ( !strncmp(linestr, "jobnr", 3) ) || ( !strncmp(linestr, "-", 1) )  )
   {
      *state = JOBS;
      return SCIP_OKAY;
   }

   /* parse job id */
   SCIPstrToIntValue(linestr, &value, &linestr);
   jobid = value - 1;

   /* construct job name */
   (void)SCIPsnprintf(jobname, SCIP_MAXSTRLEN, "%d" , jobid) ;

   /* copy job name */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &rcpspdata->jobnames[jobid], jobname, strlen(jobname) + 1) );

   /* skip next value */
   SCIPstrToIntValue(linestr, &value, &linestr);

   /* parse duration */
   SCIPstrToIntValue(linestr, &value, &linestr);
   rcpspdata->durations[jobid] = value;

   SCIP_CALL( SCIPallocBufferArray(scip, &rcpspdata->demands[jobid], rcpspdata->nresources) );

   /* parse demands */
   for( r = 0; r < rcpspdata->nresources; ++r )
   {
      SCIPstrToIntValue(linestr, &value, &linestr);
      rcpspdata->demands[jobid][r] = value;
   }

   /* check if we paresed the last job */
   if(jobid == rcpspdata->njobs - 1)
      *state = NEXT;

   return SCIP_OKAY;
}

/** get precedence relationship */
static
SCIP_RETCODE getPrecedence(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   lineno,             /**< current line number of input file */
   char*                 s,                  /**< current line */
   STATE*                state,              /**< pointer to current reading state */
   SCIP_RCPSPDATA*       rcpspdata           /**< pointer to resources constrained project scheduling data */
   )
{
   int nsuccessors;
   int value;
   int pred;
   int p;

   assert(s != NULL);
   assert(state != NULL);

   if( ( !strncmp(s, "PRECEDENCES", 3) ) ||    ( !strncmp(s, "jobnr", 4) ) )
   {
      *state = PRECEDENCES;
      return SCIP_OKAY;
   }

   /* create precedence graph if does not exist yet */
   if( rcpspdata->precedencegraph == NULL )
   {
      SCIP_CALL( SCIPdigraphCreate(&rcpspdata->precedencegraph, rcpspdata->njobs) );
   }

   /* parse predecessor */
   SCIPstrToIntValue(s, &value, &s);
   pred = value - 1;

   /* skip integer value */
   SCIPstrToIntValue(s, &value, &s);

   /* parse number of successors */
   SCIPstrToIntValue(s, &nsuccessors, &s);

   /* parse successors */
   for( p = 0; p < nsuccessors; ++p )
   {
      int succ;

      SCIPstrToIntValue(s, &value, &s);
      succ = value - 1;

      /* add precedence to digraph */
      SCIP_CALL( SCIPdigraphAddArc(rcpspdata->precedencegraph, pred, succ, (void*)(size_t)INT_MAX) );
   }

   if(pred == rcpspdata->njobs-1)
      *state = NEXT;

   return SCIP_OKAY;
}

/** compute trivial upper bound for makespan */
static
int computeUbmakespan(
   int*                  durations,          /**< array of durations */
   int                   njobs               /**< number og jobs */
   )
{
   int j;
   int ub;

   ub = 0;
   for( j = 0; j < njobs; ++j )
      ub += durations[j];

   return ub;
}

/** read file */
static
SCIP_RETCODE readFile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of input file */
   SCIP_RCPSPDATA*       rcpspdata           /**< pointer to resources constrained project scheduling data */
   )
{
   SCIP_FILE* fp;
   char       buf[SM_MAX_LINELEN];
   int        lineno = 0;
   char* s;
   STATE state = NEXT;

   assert(filename != NULL);

   if( NULL == (fp = SCIPfopen(filename, "r")) )
   {
      perror(filename);
      return SCIP_READERROR;
   }

   /* parse file line by line */
   while( state != END && state != ERROR && (NULL != SCIPfgets(buf, sizeof(buf), fp)) )
   {
      /* count line number */
      lineno++;

      if( NULL != (s = strpbrk(buf, "*\r\n")) )
         *s = '\0';
      else
      {
         parseError(scip, lineno, "line truncated", NULL, &state);
         break;
      }
      s = buf;

      /* remove white space */
      while(isspace(*s))
         s++;

      /* skip empty lines */
      if (*s == '\0')
         continue;

      if( state == NEXT )
      {
         checkForNewSection(s, &state);
      }

      SCIPdebugMessage("input line: <%s>\n", s);
      switch( state )
      {
      case ERROR:
         break;

      case NEXT:
         break;

      case NJOBS:
         SCIP_CALL( getNJobs(scip, lineno, s, &state, rcpspdata) );
         break;

      case JOBS:
         SCIP_CALL( getJobs(scip, lineno, s, &state, rcpspdata) );
         break;


      case NRESOURCES:
         SCIP_CALL( getNResources(scip, lineno, s, &state, rcpspdata) );
         break;

      case RESOURCENAMES:
         SCIP_CALL( getResourcesNames(scip, lineno, s, &state, rcpspdata) );
         break;

      case RESOURCECAPACITIES:
         SCIP_CALL( getResourcesCapacities(scip, lineno, s, &state, rcpspdata) );
         break;

      case PRECEDENCES:
         SCIP_CALL( getPrecedence(scip, lineno, s, &state, rcpspdata) );
         break;

      case END:
         parseError(scip, lineno, "additional characters after END", NULL, &state);
         break;

      default:
         SCIPerrorMessage("invalid reading state\n");
         SCIPABORT();
      }
   }
   SCIPfclose(fp);

   if( state != END && state != ERROR )
      parseError(scip, lineno, "unexpected EOF", NULL, &state);

   if( state == ERROR )
      return SCIP_READERROR;
   else
      return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopySm)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader handler */
   SCIP_CALL( SCIPincludeReaderSm(scip) );

   return SCIP_OKAY;
}

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeSm NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadSm)
{  /*lint --e{715}*/
   SCIP_RCPSPDATA rcpspdata;
   int j;

   /* initialize resources constrained project scheduling data */
   rcpspdata.precedencegraph = NULL;
   rcpspdata.jobnames = NULL;
   rcpspdata.durations = NULL;
   rcpspdata.demands = NULL;
   rcpspdata.capacities = NULL;
   rcpspdata.njobs = 0;
   rcpspdata.nresources = 0;

   /* read file */
   SCIP_CALL( readFile(scip, filename, &rcpspdata) );

   /* output rcpspdata to check it */
   SCIPdebug( outputRcpspData(scip, &rcpspdata) );

   /* create problem */
   SCIP_CALL( SCIPcreateSchedulingProblem(scip, filename, rcpspdata.jobnames, rcpspdata.resourcenames, rcpspdata.demands,
         rcpspdata.precedencegraph, rcpspdata.durations, rcpspdata.capacities, rcpspdata.njobs, rcpspdata.nresources) );

   (*result) = SCIP_SUCCESS;

   /* free buffer arrays */
   if( rcpspdata.njobs > 0 )
   {
      for( j = 0; j < rcpspdata.njobs; ++j )
      {
         SCIPfreeBufferArray(scip, &(rcpspdata.jobnames[j]));
         SCIPfreeBufferArray(scip, &(rcpspdata.demands[j]));
      }

      SCIPfreeBufferArray(scip, &rcpspdata.jobnames);
      SCIPfreeBufferArray(scip, &rcpspdata.durations);
      SCIPfreeBufferArray(scip, &rcpspdata.demands);
   }

   if( rcpspdata.nresources > 0 )
   {
      int r;

      for( r = 0; r < rcpspdata.nresources; ++r )
         SCIPfreeBufferArray(scip, &rcpspdata.resourcenames[r]);

      SCIPfreeBufferArray(scip, &rcpspdata.resourcenames);
      SCIPfreeBufferArray(scip, &rcpspdata.capacities);
   }

   if( rcpspdata.precedencegraph != NULL )
   {
      SCIPdigraphFree(&rcpspdata.precedencegraph);
   }

   return SCIP_OKAY;
}


/** problem writing method of reader */
#define readerWriteSm NULL


/*
 * reader specific interface methods
 */

/** includes the sch file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderSm(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create sch reader data */
   readerdata = NULL;

   /* include sch reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopySm) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadSm) );

   /* add reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/mipmodel", "create MIP model?",
         NULL, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates a cumulative scheduling problem */
SCIP_RETCODE SCIPcreateSchedulingProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           problemname,        /**< problem name */
   const char**          jobnames,           /**< job names, or NULL */
   const char**          resourcenames,      /**< resource names, or NULL */
   int**                 demands,            /**< demand matrix resource job demand */
   SCIP_DIGRAPH*         precedencegraph,    /**< direct graph to store the precedence conditions */
   int*                  durations,          /**< array to store the processing for each job */
   int*                  capacities,         /**< array to store the different capacities */
   int                   njobs,              /**< number of jobs to be parsed */
   int                   nresources          /**< number of capacities to be parsed */
   )
{
   SCIP_VAR** jobs;
   SCIP_VAR** vars;
   SCIP_VAR* var;

   SCIP_CONS* cons;

   char name[SCIP_MAXSTRLEN];

   int* consdurations;
   int* consdemands;

   int nvars;
   int ubmakespan;
   int i;
   int j;
   int r;

   assert( scip != NULL );
   assert( njobs >= 0 );

   SCIPdebugMessage( "start method SCIPcreateSchedulingSMProblem\n");

   /* create SCIP data structure */
   SCIP_CALL( SCIPcreateProb(scip, problemname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );


   /* compute a feasible upper bound on the makespan */
   ubmakespan = computeUbmakespan(durations, njobs);
   ubmakespan *= 100;

   /* allocate buffer for jobs and precedence constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &jobs, njobs) );

   /* create an activity constraint for each activity */
   for( j = 0; j < njobs - 1; ++j ) /* but not for last job which is the makespan (-1) */
   {
      /* construct variable name */
      if( jobnames != NULL )
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "start_%s", jobnames[j]);
      else
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "start_%d", j);

      /* create integer starting variable */
      SCIP_CALL( SCIPcreateVar(scip, &var, name, 0.0, (SCIP_Real)ubmakespan, 0.0, SCIP_VARTYPE_INTEGER,
            TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

      SCIP_CALL( SCIPaddVar(scip, var) );
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, var) );
      jobs[j] = var;
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   /* create makespan variable */
   SCIP_CALL( SCIPcreateVar(scip, &var, "makespan", 0.0, (SCIP_Real)ubmakespan, 1.0, SCIP_VARTYPE_INTEGER,
         TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPaddVar(scip, var) );
   SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, var) );

   jobs[njobs-1] = var;
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   /* precedence constraints */
   for( j = 0; j < njobs - 1; ++j )
   {
      SCIP_VAR* predvar;
      int nsuccessors;

      nsuccessors = SCIPdigraphGetNSuccessors(precedencegraph, j);

      predvar = jobs[j];
      assert(predvar != NULL);

      if( nsuccessors > 0 )
      {
         int* successors;
         void** distances;

         successors = SCIPdigraphGetSuccessors(precedencegraph, j);
         distances = SCIPdigraphGetSuccessorsDatas(precedencegraph, j);

         for( i = 0; i < nsuccessors; ++i )
         {
            SCIP_VAR* succvar;
            int distance;

            succvar = jobs[successors[i]];
            assert(succvar != NULL);

            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "precedences_(%d,%d)", j, successors[i]);

            distance = (int)(size_t)distances[i];

            if( distance == INT_MAX )
               distance = durations[j];

            SCIP_CALL( SCIPcreateConsVarbound(scip, &cons, name, predvar, succvar, -1.0,
                  -SCIPinfinity(scip), -distance,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
      }
      else
      {
         /* add precedence constraints for those jobs without successor */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "precedences_(%d,%d)", j, njobs);

         SCIP_CALL( SCIPcreateConsVarbound(scip, &cons, name, predvar, jobs[njobs-1], -1.0,
               -SCIPinfinity(scip), -durations[j],
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consdemands, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consdurations, njobs) );

   /* create resource constraints */
   for( r = 0; r < nresources; ++r )
   {
      nvars = 0;
      for( j = 0; j < njobs; ++j ) /* also makespan constraint! */
      {
         if( demands[j][r] > 0 )
         {
            vars[nvars] = jobs[j];
            consdemands[nvars] = demands[j][r];
            consdurations[nvars] = durations[j];
            nvars++;
         }
      }

      if( nvars > 0 )
      {
         /* construct constraint name */
         if( resourcenames != NULL )
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "R%s", resourcenames[r]);
         else
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "R%d", r);

         SCIP_CALL( SCIPcreateConsCumulative(scip, &cons, name,
               nvars, vars, consdurations, consdemands, capacities[r],
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   /* initialize the problem specific heuristic */
   SCIP_CALL( SCIPinitializeHeurListScheduling(scip, precedencegraph, jobs,
         durations, demands, capacities, njobs, nresources) );

   /* free buffer array */
   SCIPfreeBufferArray(scip, &consdurations);
   SCIPfreeBufferArray(scip, &consdemands);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &jobs);

   return SCIP_OKAY;
}
