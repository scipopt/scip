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

/**@file   main.c
 * @brief  Main file for C compilation
 * @author Leon Eifler
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/scipshell.h"
#include "scip/message_default.h"
#include <string.h>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "cycplugins.h"
#include "probdata_cyc.h"
#include "reader_cyc.h"
#include "scip/debug.h"

#define COL_MAX_LINELEN 1024

/** Read the parameters from the command Line */
static
SCIP_RETCODE readParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< parameter file name */
   )
{
   if( SCIPfileExists(filename) )
   {
      /* read params from settingsfile */
      SCIPinfoMessage(scip, NULL, "reading user parameter file <%s>\n", filename);
      SCIP_CALL( SCIPreadParams(scip, filename) );
   }
   else
      SCIPinfoMessage(scip, NULL, "user parameter file <%s> not found - using default parameters\n", filename);

   return SCIP_OKAY;
}

/** execute the scip-program from the command-line */
static
SCIP_RETCODE fromCommandLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< input file name */
   const char*           soluname            /**< input file name */
   )
{
   SCIP_RETCODE retcode;
   /********************
    * Problem Creation *
    ********************/

   /** @note The message handler should be only fed line by line such the message has the chance to add string in front
    *        of each message
    */
   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "read problem <%s>\n", filename);
   SCIPinfoMessage(scip, NULL, "============\n");
   SCIPinfoMessage(scip, NULL, "\n");

   retcode = SCIPreadProb(scip, filename, NULL);

   switch( retcode )
   {
   case SCIP_NOFILE:
      SCIPinfoMessage(scip, NULL, "file <%s> not found\n", filename);
      return SCIP_OKAY;
   case SCIP_PLUGINNOTFOUND:
      SCIPinfoMessage(scip, NULL, "no reader for input file <%s> available\n", filename);
      return SCIP_OKAY;
   case SCIP_READERROR:
      SCIPinfoMessage(scip, NULL, "error reading file <%s>\n", filename);
      return SCIP_OKAY;
   default:
      SCIP_CALL( retcode );
   } /*lint !e788*/

   if( soluname != NULL )
   {
      retcode = SCIPreadProb(scip, soluname, NULL);

      switch( retcode )
      {
      case SCIP_NOFILE:
         SCIPinfoMessage(scip, NULL, "file <%s> not found\n", filename);
         return SCIP_OKAY;
      case SCIP_PLUGINNOTFOUND:
         SCIPinfoMessage(scip, NULL, "no reader for input file <%s> available\n", filename);
         return SCIP_OKAY;
      case SCIP_READERROR:
         SCIPinfoMessage(scip, NULL, "error reading file <%s>\n", filename);
         return SCIP_OKAY;
      default:
         SCIP_CALL( retcode );
      } /*lint !e788*/
   }

   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   SCIPinfoMessage(scip, NULL, "\nsolve problem\n");
   SCIPinfoMessage(scip, NULL, "=============\n\n");

   SCIP_CALL( SCIPsolve(scip) );

   /*******************
    * Solution Output *
    *******************/

   SCIPinfoMessage(scip, NULL, "\nprimal solution (transformed space):\n");
   SCIPinfoMessage(scip, NULL, "====================================\n\n");

   /**************
    * Statistics *
    **************/

   SCIPinfoMessage(scip, NULL, "\nStatistics\n");
   SCIPinfoMessage(scip, NULL, "==========\n\n");
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   return SCIP_OKAY;
}

/** process the arguments and set up the problem */
static
SCIP_RETCODE processArguments(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
   char* probname = NULL;
   char* soluname = NULL;
   char* settingsname = NULL;
   char* logname = NULL;
   char name_file[COL_MAX_LINELEN];

   SCIP_Bool quiet;
   SCIP_Bool paramerror;
   SCIP_Bool interactive;
   int i;

   /********************
    * Parse parameters *
    ********************/

   quiet = FALSE;
   paramerror = FALSE;
   interactive = (argc == 0);

   /*lint -e{850} read the arguments from commandLine */
   for( i = 1; i < argc; ++i )
   {
      if( strcmp(argv[i], "-l") == 0 )
      {
         i++;
         if( i < argc )
            logname = argv[i];
         else
         {
            printf("missing log filename after parameter '-l'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-q") == 0 )
         quiet = TRUE;
      else if( strcmp(argv[i], "-s") == 0 )
      {
         i++;
         if( i < argc )
            settingsname = argv[i];
         else
         {
            printf("missing settings filename after parameter '-s'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-f") == 0 )
      {
         i++;
         if( i < argc )
         {
            probname = argv[i];
            (void)SCIPsnprintf( name_file, SCIP_MAXSTRLEN, argv[i] );
         }
         else
         {
            printf("missing problem filename after parameter '-f'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-c") == 0 )
      {
         i++;
         if( i < argc )
         {
            SCIP_CALL( SCIPaddDialogInputLine(scip, argv[i]) );
            interactive = TRUE;
         }
         else
         {
            printf("missing command line after parameter '-c'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-x") == 0 )
      {
         i++;
         if( i < argc )
         {
            soluname = argv[i];
         }
         else
         {
            printf("missing solution filename after parameter '-x'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-b") == 0 )
      {
         i++;
         if( i < argc )
         {
            SCIP_FILE* file;

            file = SCIPfopen(argv[i], "r");
            if( file == NULL )
            {
               printf("cannot read command batch file <%s>\n", argv[i]);
               SCIPprintSysError(argv[i]);
               paramerror = TRUE;
            }
            else
            {
               while( !SCIPfeof(file) )
               {
                  char buffer[SCIP_MAXSTRLEN];

                  (void)SCIPfgets(buffer, (int) sizeof(buffer), file);
                  if( buffer[0] != '\0' )
                  {
                     SCIP_CALL( SCIPaddDialogInputLine(scip, buffer) );
                  }
               }
               SCIPfclose(file);
               interactive = TRUE;
            }
         }
         else
         {
            printf("missing command batch filename after parameter '-b'\n");
            paramerror = TRUE;
         }
      }
   }

   if( interactive && probname != NULL )
   {
      printf("cannot mix batch mode '-c' and '-b' with file mode '-f'\n");
      paramerror = TRUE;
   }

   if( !paramerror )
   {
      /***********************************
       * create log file message handler *
       ***********************************/

      if( quiet )
      {
         SCIPsetMessagehdlrQuiet(scip, quiet);
      }

      if( logname != NULL )
      {
         SCIPsetMessagehdlrLogfile(scip, logname);
      }

      /***********************************
       * Version and library information *
       ***********************************/

      SCIPprintVersion(scip, NULL);
      SCIPinfoMessage(scip, NULL, "\n");

      SCIPprintExternalCodes(scip, NULL);
      SCIPinfoMessage(scip, NULL, "\n");

      /*****************
       * Load settings *
       *****************/

      if( settingsname != NULL )
      {
         SCIP_CALL( readParams(scip, settingsname) );
      }
      else if( defaultsetname != NULL )
      {
         SCIP_CALL( readParams(scip, defaultsetname) );
      }
      /**************
       * Start SCIP *
       **************/

      if( probname != NULL )
      {
         /* run scip */
         SCIP_CALL( fromCommandLine(scip, probname, soluname) );

      }
      else
      {
         SCIPinfoMessage(scip, NULL, "\n");
         SCIP_CALL( SCIPstartInteraction(scip) );
      }
   }
   else
   {
      printf("\nsyntax: %s [-l <logfile>] [-q] [-s <settings>] [-f <problem>]\n"
         "  -l <logfile>  : copy output into log file\n"
         "  -q            : suppress screen messages\n"
         "  -s <settings> : load parameter settings (.set) file\n"
         "  -f <problem>  : load and solve problem file\n\n",
         argv[0]);
   }

   return SCIP_OKAY;
}

/** Set up the problem-structure and solve the clustering problem */
static
SCIP_RETCODE SCIPrunCyc(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;
   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

#ifdef WITH_DEBUG_SOLUTION
   SCIPdebugSolEnable(scip);
#endif
   /* include reader, problemdata*/
   SCIP_CALL( SCIPincludeCycPlugins(scip) );

   /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL( processArguments(scip, argc, argv, defaultsetname) );

   SCIPinfoMessage(scip, NULL, "\n");

   /* free scip */
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main method */
int
main(
   int                   argc,
   char**                argv
   )
{
   SCIP_RETCODE retcode;

   retcode = SCIPrunCyc(argc, argv, "scip.set");

   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
