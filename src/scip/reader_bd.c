/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   reader_bd.c
 * @ingroup DEFPLUGINS_READER
 * @brief  Benders file reader - this is used to read multiple files that form the master and subproblems of a
 *         Benders' decomposition
 * @author Stephen J. Maher
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/benders_default.h"
#include "scip/def.h"
#include "scip/reader_bd.h"
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_reader.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_reader.h"


#define READER_NAME             "bendersreader"
#define READER_DESC             "file reader for multiple instance files for a Benders' decomposition"
#define READER_EXTENSION        "bd"

/** enum for the file types that are read by the BENDERS reader */
enum BendersFile_Section
{
   BENDERSFILE_SECTION_INIT = 0,
   BENDERSFILE_SECTION_MASTER = 1,
   BENDERSFILE_SECTION_SUB = 2
};
typedef enum BendersFile_Section BENDERSFILE_SECTION;


/** reads the given benders input file */
static
SCIP_RETCODE readBendersInputFile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_RETCODE retcode;
   SCIP_FILE* file;
   char tempfile[SCIP_MAXSTRLEN];
   char masterfile[SCIP_MAXSTRLEN];
   char** subproblemfiles = NULL;
   SCIP_Bool error;
   int nsubproblems = -1;
   int subprobcount = 0;
   int i;

   char* fromlastslash;
   char parent[SCIP_MAXSTRLEN];
   size_t parentlen;

   BENDERSFILE_SECTION section;

   assert(scip != NULL);
   assert(filename != NULL);

   retcode = SCIP_OKAY;

   /* finding the path to the file to that it can be prepended to the master and subproblem file names */
   fromlastslash = (char*) strrchr(filename, '/');

   if( fromlastslash == NULL )
      parentlen = 0;
   else
      parentlen = fromlastslash + 1 - filename;

   (void)SCIPstrncpy(parent, filename, (int)parentlen + 1);

   /* open input file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* read the file */
   error = FALSE;

   /* start parsing the file */
   section = BENDERSFILE_SECTION_INIT;
   while( !SCIPfeof(file) && !error )
   {
      char buffer[SCIP_MAXSTRLEN];
      int nread;

      /* get next line */
      if( SCIPfgets(buffer, (int) sizeof(buffer), file) == NULL )
         break;

      /* check if a new section begins */
      if( strncmp(buffer, "MASTER", 6) == 0 )
      {
         section = BENDERSFILE_SECTION_MASTER;
         continue;
      }
      else if( strncmp(buffer, "SUBPROBLEMS", 11) == 0 )
      {
         section = BENDERSFILE_SECTION_SUB;

         /* coverity[secure_coding] */
         nread = sscanf(buffer, "SUBPROBLEMS %1018d\n", &nsubproblems);
         if( nread < 1 )
         {
            error = TRUE;
            break;
         }

         /* if the number of subproblems is less than or equal 0, then this is an invalid entry */
         if( nsubproblems <= 0 )
         {
            SCIPwarningMessage(scip, "The value <%d> is not valid as the number of subproblems. Please fix.", nsubproblems);
            error = TRUE;
            break;
         }

         SCIPdebugMsg(scip, "There are %d subproblem files to be read\n", nsubproblems);

         /* allocating memory for the subproblem file names */
         SCIP_CALL_TERMINATE( retcode, SCIPallocBufferArray(scip, &subproblemfiles, nsubproblems), TERMINATE );

         continue;
      }

      /* base the parsing on the currently active section */
      switch (section)
      {
         case BENDERSFILE_SECTION_INIT:
            break;
         case BENDERSFILE_SECTION_MASTER:
            /* reading the master problem file name */
            nread = sscanf(buffer, "%1023s\n", tempfile);
            if( nread < 1 )
               error = TRUE;

            SCIPdebugMsg(scip, "The instance file for the master problem is <%s>\n", tempfile);

            /* prepending the parent path only if the master problem file is a relative path */
            if( tempfile[0] == '/' )
               (void) SCIPsnprintf(masterfile, SCIP_MAXSTRLEN, "%s", tempfile);
            else
               (void) SCIPsnprintf(masterfile, SCIP_MAXSTRLEN, "%s%s", parent, tempfile);

            break;
         case BENDERSFILE_SECTION_SUB:
            /* the number of subproblems must be known before reading the files */
            if( nsubproblems == -1 )
            {
               SCIPwarningMessage(scip, "The number of subproblems must be specified with the SUBPROBLEM marker\n");
               error = TRUE;
               break;
            }

            if( subprobcount >= nsubproblems )
            {
               SCIPerrorMessage("Error: There are more subproblem files listed than specified: "
                  "Specified %d subproblems.\n", nsubproblems);
               error = TRUE;
               break;
            }

            /* allocating the space for the file names */
            assert(subproblemfiles != NULL);
            SCIP_CALL_TERMINATE( retcode, SCIPallocBufferArray(scip, &(subproblemfiles[subprobcount]), SCIP_MAXSTRLEN), TERMINATE );

            /* read subproblem file name */
            /* coverity[secure_coding] */
            nread = sscanf(buffer, "%1023s\n", tempfile);
            if( nread < 1 )
               error = TRUE;

            SCIPdebugMsg(scip, "The instance file for the subproblem %d is <%s>\n", subprobcount, tempfile);

            /* prepending the parent path only if the subproblem file is a relative path */
            if( tempfile[0] == '/' )
               (void) SCIPsnprintf(subproblemfiles[subprobcount], SCIP_MAXSTRLEN, "%s", tempfile);
            else
               (void) SCIPsnprintf(subproblemfiles[subprobcount], SCIP_MAXSTRLEN, "%s%s", parent, tempfile);

            subprobcount++;
            break;

         default:
            break;
      }
   }

   /* close input file */
   SCIPfclose(file);

   /* compare nsubproblems and the subproblem file count. If there is a mismatch, then exit with an error */
   if( nsubproblems != subprobcount )
   {
      SCIPerrorMessage("Error: Number of subproblems specification is wrong: Specified %d subproblems, counted %d.\n",
         nsubproblems, subprobcount);
      error = TRUE;
   }

   /* create a decomposition and add it to the decomposition storage of SCIP */
   if( ! error )
   {
      SCIPinfoMessage(scip, NULL, "reading master problem file <%s>\n", masterfile);
      SCIPinfoMessage(scip, NULL, "============\n");

      /* reading the master problem instance file */
      SCIP_CALL( SCIPreadProb(scip, masterfile, NULL) );

      SCIPinfoMessage(scip, NULL, "\n\n");

      SCIPinfoMessage(scip, NULL, "reading the subproblem instances\n");
      SCIPinfoMessage(scip, NULL, "============\n");
      /* reading the subproblem instance files in the default Benders' decomposition */
      SCIP_CALL( SCIPcreateBendersDefaultFromFiles(scip, subproblemfiles, nsubproblems) );

      /* activating the Benders' constraint handler. The two-phase method is activated by default. If the user desires
       * not to use the two-phase method, then the setting in cons_benderslp must be explicitly changed.
       */
      SCIP_CALL( SCIPsetBoolParam(scip, "constraints/benderslp/active", TRUE) );
      SCIP_CALL( SCIPsetBoolParam(scip, "constraints/benders/active", TRUE) );

      /* changing settings that are required for Benders' decomposition
       *
       * TODO: review if all of these are necessary
       */
      SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
      SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxrounds", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "heuristics/trysol/freq", 1) );

      /* disabling aggregation since it can affect the mapping between the master and subproblem variables */
      SCIP_CALL( SCIPsetBoolParam(scip, "presolving/donotaggr", TRUE) );
      SCIP_CALL( SCIPsetBoolParam(scip, "presolving/donotmultaggr", TRUE) );

      SCIPinfoMessage(scip, NULL, "\n");
      SCIPinfoMessage(scip, NULL, "Original problem statistics\n");
      SCIPinfoMessage(scip, NULL, "============\n");
   }
   else
   {
      SCIPerrorMessage("Errors parsing the Benders' decomposition files in <%s>. The decomposition could not be created\n.", filename);
   }

TERMINATE:
   /* the subproblem files array only needs to be freed if there are subproblem files to read */
   if( nsubproblems > 0 )
   {
      for( i = subprobcount - 1; i >= 0; i-- )
         SCIPfreeBufferArray(scip, &(subproblemfiles[i]));

      SCIPfreeBufferArray(scip, &subproblemfiles);
   }
   assert(subproblemfiles == NULL);

   if( retcode != SCIP_OKAY )
   {
      SCIPfclose(file);
      return retcode;
   }

   if( error )
      return SCIP_READERROR;
   else
      return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyBenders)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);

   SCIP_STRINGEQ( SCIPreaderGetName(reader), READER_NAME, SCIP_INVALIDCALL );

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderBenders(scip) );

   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadBenders)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(result != NULL);

   SCIP_STRINGEQ( SCIPreaderGetName(reader), READER_NAME, SCIP_INVALIDCALL );

   *result = SCIP_DIDNOTRUN;

   if( filename == NULL )
   {
      SCIPerrorMessage("A file name must be specified for the Benders' instance reader.");

      return SCIP_READERROR;
   }

   SCIP_CALL( readBendersInputFile(scip, filename) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the benders file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderBenders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyBenders) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadBenders) );

   return SCIP_OKAY;
}
