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
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cmain.c
 * @brief  main file for AMPL interface
 * @author Stefan Vigerske
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "reader_nl.h"

static
SCIP_RETCODE run(
   const char*                nlfile,        /**< name of AMPL .nl file */
   const char*                setfile,       /**< SCIP settings file, or NULL to try default scip.set */
   SCIP_Bool                  interactive    /**< whether to start SCIP interactive shell instead of solving command */
)
{
   SCIP* scip;
   char buffer[SCIP_MAXSTRLEN];
   SCIP_Bool printstat;
   char* logfile = NULL;

   assert(nlfile != NULL);

   /* setup SCIP and print version information */
   SCIP_CALL( SCIPcreate(&scip) );

   /* we explicitly have to enable the use of a debug solution for this main SCIP instance */
   SCIPenableDebugSol(scip);

   SCIPprintVersion(scip, NULL);
   SCIPinfoMessage(scip, NULL, "\n");

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPincludeReaderNl(scip) );

   SCIP_CALL( SCIPaddBoolParam(scip, "display/statistics",
      "whether to print statistics on a solve",
      &printstat, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "display/logfile",
      "name of file to write SCIP log to (additionally to writing to stdout)",
      &logfile, FALSE, "", NULL, NULL) );
   assert(logfile != NULL);

   SCIPprintExternalCodes(scip, NULL);
   SCIPinfoMessage(scip, NULL, "\n");

   /* read setting file */
   if( setfile != NULL )
   {
      SCIP_CALL( SCIPreadParams(scip, setfile) );
   }

   /* set logfile, if not default empty string */
   if( *logfile )
      SCIPsetMessagehdlrLogfile(scip, logfile);

   /* setup commands to be executed by SCIP */

   SCIP_CALL( SCIPaddDialogInputLine(scip, "display param") );

   /* add .nl extension, if not given */
   (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "read %s%s", nlfile, (strlen(nlfile) < 3 || strcmp(nlfile+(strlen(nlfile)-3), ".nl") != 0) ? ".nl" : "");
   SCIP_CALL( SCIPaddDialogInputLine(scip, buffer) );

   if( !interactive )
   {
      /* SCIP_CALL( SCIPaddDialogInputLine(scip, "display problem") ); */

      SCIP_CALL( SCIPaddDialogInputLine(scip, "optimize") );

      SCIP_CALL( SCIPaddDialogInputLine(scip, "write amplsol") );

      if( printstat )
      {
         SCIP_CALL( SCIPaddDialogInputLine(scip, "display statistics") );
      }

      SCIP_CALL( SCIPaddDialogInputLine(scip, "quit") );
   }

   /* run SCIP */
   SCIP_CALL( SCIPstartInteraction(scip) );

   SCIP_CALL( SCIPfree(&scip) );

   return SCIP_OKAY;
}

/** main method starting SCIP */
int main(
   int                        argc,          /**< number of arguments from the shell */
   char**                     argv           /**< array of shell arguments */
   )
{
   SCIP_RETCODE retcode;
   SCIP_Bool interactive;
   const char* setfile;
   int i;

   assert(argc >= 1);
   if( argc == 1 )
   {
      SCIP* scip;

      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      SCIPprintVersion(scip, NULL);
      printf("\n");
      SCIP_CALL_ABORT( SCIPfree(&scip) );

      printf("Usage: %s <nl-file> { [-AMPL] | [<settings-file>] | -i }\n", argv[0]);
      printf("\n");
      printf("       -i starts the SCIP shell after reading settings and .nl file, instead of starting the SCIP solve\n");

      return 0;
   }

   /* check command line arguments after .nl file */
   interactive = FALSE;
   setfile = NULL;
   for( i = 2; i < argc; ++i )
   {
      if( strcmp(argv[i], "-AMPL") == 0 )
         continue;

      if( strcmp(argv[i], "-i") == 0 )
      {
         interactive = TRUE;
         continue;
      }

      setfile = argv[i];
   }

   if( setfile == NULL && SCIPfileExists("scip.set") )
      setfile = "scip.set";

   retcode = run(argv[1], setfile, interactive);

   /* evaluate retrun code of the SCIP process */
   if( retcode != SCIP_OKAY )
   {
      /* write error back trace */
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
