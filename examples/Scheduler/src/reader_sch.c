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

/**@file   reader_sch.c
 * @brief  scheduling problem file reader for RCPSP/max format
 * @author Stefan Heinz
 *
 * This reader is capabale of parsing resource-constrained project scheduling problem with minimal and maximal time lags
 * (RCPSP/max) instances. The <a http://www.wior.uni-karlsruhe.de/LS_Neumann/Forschung/ProGenMax/rcpspmax.html">PSPlib</a>
 * provides several instances set.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#include "reader_sch.h"
#include "reader_sm.h"


#define READER_NAME             "schreader"
#define READER_DESC             "scheduling file reader for sch files (RCPSP/max format)"
#define READER_EXTENSION        "sch"


#define SCH_MAX_LINELEN      65536     /**< size of the line buffer for reading or writing */

/*
 * Local methods
 */

/** parse job id and check if only one job mode is present */
static
SCIP_RETCODE getJobId(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to search */
   int*                  job,                /**< pointer to store the parsed job id */
   char**                endptr              /**< pointer to store the final string position if successfully parsed */
   )
{
   int mode;

   /* get job id */
   SCIPstrToIntValue(str, job, endptr);

   /* get job mode */
   SCIPstrToIntValue(*endptr, &mode, endptr);

   if( mode != 1 )
   {
      SCIPwarningMessage(scip, "jobs with different modes are not supported\n");
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** parse job and capacities details */
static
SCIP_RETCODE parseDetails(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            file,               /**< file to parse */
   int*                  lineno,             /**< pointer to store line number of the file */
   int**                 demands,            /**< demand matrix resource job demand */
   SCIP_DIGRAPH*         precedencegraph,    /**< direct graph to store the precedence conditions */
   int*                  durations,          /**< array to store the processing for each job */
   int*                  capacities,         /**< array to store the different capacities */
   int                   njobs,              /**< number of jobs to be parsed */
   int                   nresources          /**< number of capacities to be parsed */
   )
{
   char buf[SCH_MAX_LINELEN];
   char* endptr;
   int j;

   /* get data for each job including a dummy job at the beginning and a dummy job at the end */
   for( j = 0; j < njobs; ++j )
   {
      if( NULL != SCIPfgets(buf, sizeof(buf), file) )
      {
         int* successors;
         int distance;
         int nsuccessors;
         int job;
         int s;

         /* get job id and check if only one mode is present */
         SCIP_CALL( getJobId(scip, buf, &job, &endptr) );

         SCIPdebugMessage("job %d -> ", j);

         /* get number of direct successors */
         SCIPstrToIntValue(endptr, &nsuccessors, &endptr);

         /* allocate buffer to temporarily collect the successors */
         SCIP_CALL( SCIPallocBufferArray(scip, &successors, nsuccessors) );

         /* parse successor job ids */
         for( s = 0; s < nsuccessors; ++s )
         {
            if( !SCIPstrToIntValue(endptr, &successors[s], &endptr) )
               return SCIP_READERROR;
         }

         /* parse distances between the job and its successor and add the arc with their data to the precedence graph */
         for( s = 0; s < nsuccessors; ++s )
         {
            char token[SCIP_MAXSTRLEN];
            char* tmpptr;

            SCIPstrCopySection(endptr, '[', ']', token, SCIP_MAXSTRLEN, &endptr);

            if( SCIPstrToIntValue(token, &distance, &tmpptr) )
            {
               SCIP_CALL( SCIPdigraphAddArc(precedencegraph, job, successors[s], (void*)(size_t)distance) );

               SCIPdebugPrintf(" %d[%d] ", successors[s], distance);
            }
            else
               return SCIP_READERROR;
         }

         SCIPdebugPrintf("\n");

         /* free the buffers */
         SCIPfreeBufferArray(scip, &successors);
      }
      else
         return SCIP_READERROR;

      (*lineno)++;
   }

   /* get data for each job including a dummy job at the beginning and a dummy job at the end */
   for( j = 0; j < njobs; ++j )
   {
      if( NULL != SCIPfgets(buf, sizeof(buf), file) )
      {
         int job;
         int r;

         /* get job id and check if only one mode is present */
         SCIP_CALL( getJobId(scip, buf, &job, &endptr) );

         /* get processing time */
         SCIPstrToIntValue(endptr, &durations[job], &endptr);

         SCIPdebugMessage("job %d has a processing times: %d\n", job, durations[job]);

         for( r = 0; r < nresources; ++r )
         {
            SCIPstrToIntValue(endptr, &demands[job][r], &endptr);
         }
      }
      else
         return SCIP_READERROR;

      (*lineno)++;
   }

   /* get resources capacities */
   if( nresources > 0 && NULL != SCIPfgets(buf, sizeof(buf), file) )
   {
      int r;

      SCIPdebugMessage("line %d %s", *lineno, buf);

      SCIPstrToIntValue(buf, &capacities[0], &endptr);

      SCIPdebugMessage("paresed capacities: <%d>", capacities[0]);

      for( r = 1; r < nresources; ++r )
      {
         SCIPstrToIntValue(endptr, &capacities[r], &endptr);

         SCIPdebugPrintf(", <%d>", capacities[r]);
      }

      SCIPdebugPrintf("\n");
   }
   else
      return SCIP_READERROR;

   (*lineno)++;

   return SCIP_OKAY;
}

/** read file and create problem */
static
SCIP_RETCODE readFile(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            file,               /**< file to pares */
   const char*           filename            /**< name of input file */
   )
{
   SCIP_RETCODE retcode;
   char buf[SCH_MAX_LINELEN];
   SCIP_DIGRAPH* precedencegraph;
   int** demands;
   int* durations;
   int* capacities;
   int lineno;
   int njobs;
   int nresources;
   int j;

   assert(scip != NULL);
   assert(file != NULL);
   assert(filename != NULL);

   lineno = 0;

   /* get number of jobs and resources */
   if( NULL != SCIPfgets(buf, sizeof(buf), file) )
   {
      char* endptr;
      int value;

      lineno++;

      SCIPstrToIntValue(buf, &value, &endptr);

      /* note that this format includes two dummy jobs */
      njobs = value + 2;

      /* get number of resources */
      SCIPstrToIntValue(endptr, &nresources, &endptr);
   }
   else
      return SCIP_READERROR;

   SCIP_CALL( SCIPallocBufferArray(scip, &capacities, nresources) );
   SCIP_CALL( SCIPallocBufferArray(scip, &durations, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, njobs) );

   for( j = 0; j < njobs; ++j )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &demands[j], nresources) );
      BMSclearMemoryArray(demands[j], nresources);
   }

   SCIP_CALL( SCIPdigraphCreate(&precedencegraph, njobs) );

   SCIPdebugMessage("problem has <%d> jobs and <%d> resources\n", njobs, nresources);

   retcode = parseDetails(scip, file, &lineno, demands, precedencegraph, durations, capacities, njobs, nresources);

   if( retcode == SCIP_OKAY )
   {
      SCIP_CALL( SCIPcreateSchedulingProblem(scip, filename, NULL, NULL, demands,
            precedencegraph, durations, capacities, njobs, nresources) );
   }

   /* free the precedence graph */
   SCIPdigraphFree(&precedencegraph);

   /* free buffer before evaluating the retcode */
   for( j = njobs - 1; j >= 0; --j )
   {
      SCIPfreeBufferArray(scip, &demands[j]);
   }
   SCIPfreeBufferArray(scip, &demands);
   SCIPfreeBufferArray(scip, &durations);
   SCIPfreeBufferArray(scip, &capacities);

   SCIP_CALL( retcode );

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopySch)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader handler */
   SCIP_CALL( SCIPincludeReaderSch(scip) );

   return SCIP_OKAY;
}

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeSch NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadSch)
{  /*lint --e{715}*/
   SCIP_FILE* file;
   SCIP_RETCODE retcode;

   if( NULL == (file = SCIPfopen(filename, "r")) )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* read file */
   retcode = readFile(scip, file, filename);

   /* close file */
   SCIPfclose(file);

   /* check retcode after the file was closed */
   SCIP_CALL( retcode );

   (*result) = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** problem writing method of reader */
#define readerWriteSch NULL


/*
 * reader specific interface methods
 */

/** includes the sch file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderSch(
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

   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopySch) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadSch) );

   /* add reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/mipmodel", "create MIP model?",
         NULL, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
