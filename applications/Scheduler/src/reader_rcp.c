/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_rcp.c
 * @brief  file reader for "pack" scheduling instances
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "reader_rcp.h"
#include "reader_sm.h"

/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "rcpreader"
#define READER_DESC             "reader for \"pack\" scheduling instances"
#define READER_EXTENSION        "rcp"

/**@} */


/**@name Local methods
 *
 * @{
 */

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
   char buf[SCIP_MAXSTRLEN];
   char* endptr;
   int j;

   /* get resources capacities */
   if( nresources > 0 && NULL != SCIPfgets(buf, (int) sizeof(buf), file) )
   {
      int r;

      SCIPdebugMessage("line %d %s", *lineno, buf);

      if( !SCIPstrToIntValue(buf, &capacities[0], &endptr) )
         return SCIP_READERROR;

      SCIPdebugMessage("paresed capacities: <%d>", capacities[0]);

      for( r = 1; r < nresources; ++r )
      {
         if( !SCIPstrToIntValue(endptr, &capacities[r], &endptr) )
            return SCIP_READERROR;

         SCIPdebugPrintf(", <%d>", capacities[r]);
      }

      SCIPdebugPrintf("\n");

      (*lineno)++;
   }
   else
      return SCIP_READERROR;

   /* get job details */
   for( j = 0; j < njobs; ++j )
   {
      if( NULL != SCIPfgets(buf, (int) sizeof(buf), file) )
      {
         int nsuccessors;
         int r;
         int s;

         /* get job duration */
         if( !SCIPstrToIntValue(buf, &durations[j], &endptr) )
            return SCIP_READERROR;

         SCIPdebugMessage("job %d: duration %d, demands (", j, durations[j]);

         /* parse resources demands */
         for( r = 0; r < nresources; ++r )
         {
            if( !SCIPstrToIntValue(endptr, &demands[j][r], &endptr) )
               return SCIP_READERROR;

            SCIPdebugPrintf(" %d ", demands[j][r]);
         }

         /* get number of successors */
         if( !SCIPstrToIntValue(endptr, &nsuccessors, &endptr) )
            return SCIP_READERROR;

         SCIPdebugPrintf("), successors %d:", nsuccessors);

         /* parse successor job ids */
         for( s = 0; s < nsuccessors; ++s )
         {
            int successor;

            if( !SCIPstrToIntValue(endptr, &successor, &endptr) )
               return SCIP_READERROR;

            /* add precedence to precedence graph */
            SCIP_CALL( SCIPdigraphAddArc(precedencegraph, j, successor-1, (void*)(size_t)INT_MAX) );

            SCIPdebugPrintf(" %d ", successor);
         }

         SCIPdebugPrintf("\n");
      }
      else
         return SCIP_READERROR;

      (*lineno)++;
   }

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
   char buf[SCIP_MAXSTRLEN];
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
   if( NULL != SCIPfgets(buf, (int) sizeof(buf), file) )
   {
      char* endptr;

      lineno++;

      /* get number of jobs */
      if( !SCIPstrToIntValue(buf, &njobs, &endptr) )
         return SCIP_READERROR;

      /* get number of resources */
      if( !SCIPstrToIntValue(endptr, &nresources, &endptr) )
         return SCIP_READERROR;
   }
   else
      return SCIP_READERROR;

   SCIP_CALL( SCIPallocBufferArray(scip, &capacities, nresources) );
   SCIP_CALL( SCIPallocBufferArray(scip, &durations, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, njobs) );

   for( j = 0; j < njobs; ++j )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &demands[j], nresources) ); /*lint !e866*/
      BMSclearMemoryArray(demands[j], nresources); /*lint !e866*/
   }

   SCIP_CALL( SCIPcreateDigraph(scip, &precedencegraph, njobs) );

   SCIPdebugMessage("problem has <%d> jobs and <%d> resources\n", njobs, nresources);

   retcode = parseDetails(scip, file, &lineno, demands, precedencegraph, durations, capacities, njobs, nresources);

   /* create problem */
   if( retcode == SCIP_OKAY )
   {
      SCIP_CALL( SCIPcreateSchedulingProblem(scip, filename, NULL, NULL, demands,
            precedencegraph, durations, capacities, njobs, nresources, TRUE) );
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

/**@} */

/**@name Callback methods of reader
 *
 * @{
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyRcp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader handler */
   SCIP_CALL( SCIPincludeReaderRcp(scip) );

   return SCIP_OKAY;
}

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeSch NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadRcp)
{  /*lint --e{715}*/
   SCIP_FILE* file;
   SCIP_RETCODE retcode;

   if( NULL == (file = SCIPfopen(filename, "r")) )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* read file and create problem */
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


/**@} */

/**@name Interface methods
 *
 * @{
 */

/*
 * reader specific interface methods
 */

/** includes the rcp file reader into SCIP */
SCIP_RETCODE SCIPincludeReaderRcp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* include sch reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyRcp, readerFreeSch, readerReadRcp, readerWriteSch, NULL) );

   return SCIP_OKAY;
}

/**@} */
