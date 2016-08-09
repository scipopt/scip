
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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
#include "reader_spa.h"
#include "spaplugins.h"
#include "probdata_spa.h"
#include <string.h>
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
   SCIP_Bool outputorigsol = FALSE;
   SCIP_Real eps;
   char model;
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
   }


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
      }
   }

   /*******************
    * Problem Solving *
    *******************/
   model = SCIPspaGetModel(scip);
   if( model != 'p' )
   {
      SCIP_SOL* bestsol;
      /* solve problem */
      SCIP_CALL( SCIPgetRealParam(scip, "coherence_bound", &eps) );
      SCIPinfoMessage(scip, NULL, "\nsolve problem\n");
      SCIPinfoMessage(scip, NULL, "=============\n\n");

      SCIP_CALL( SCIPsolve(scip) );
      bestsol = SCIPgetBestSol(scip);
      /*******************
       * Solution Output *
       *******************/

      SCIP_CALL( SCIPgetBoolParam(scip, "misc/outputorigsol", &outputorigsol) );
      if ( outputorigsol )
      {
         SCIPinfoMessage(scip, NULL, "\nprimal solution (original space):\n");
         SCIPinfoMessage(scip, NULL, "=================================\n\n");

         if ( bestsol == NULL )
            SCIPinfoMessage(scip, NULL, "no solution available\n");
         else
         {
            SCIP_SOL* origsol;
            SCIP_CALL( SCIPcreateSolCopy(scip, &origsol, bestsol) );
            SCIP_CALL( SCIPretransformSol(scip, origsol) );
            SCIP_CALL( SCIPprintSol(scip, origsol, NULL, FALSE) );
            SCIPspaPrintSolutionValues(scip, origsol);
            SCIP_CALL( SCIPfreeSol(scip, &origsol) );
         }
      }
      else
      {
         SCIPinfoMessage(scip, NULL, "\nprimal solution (transformed space):\n");
         SCIPinfoMessage(scip, NULL, "====================================\n\n");
         SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );
      }

      /**************
       * Statistics *
       **************/

      SCIPinfoMessage(scip, NULL, "\nStatistics\n");
      SCIPinfoMessage(scip, NULL, "==========\n\n");
      SCIP_CALL( SCIPprintStatistics(scip, NULL) );
   }
   return SCIP_OKAY;
}

/** process the arguments and set up the problem */
static
SCIP_RETCODE processArguments(
   SCIP*                 scip,
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
)
{
   char* probname = NULL;
   char* soluname = NULL;
   char* file_ending;
   char filename[COL_MAX_LINELEN];
   char* settingsname = NULL;
   char* logname = NULL;
   char name_file[COL_MAX_LINELEN];
   char ending[COL_MAX_LINELEN] ;
   SCIP_Bool quiet;
   SCIP_Bool paramerror;
   SCIP_Bool interactive;
   int i;
   int ncluster;
   SCIP_Real eps=0;
   FILE* output;

   /********************
    * Parse parameters *
    ********************/

   quiet = FALSE;
   paramerror = (argc == 1);
   interactive = FALSE;

   /* read the arguments from commandLine */
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
         /*create the filename */
         snprintf(filename, 50, "solutions/");
         file_ending = strrchr(name_file, '/') + 1;
         strcat(filename, file_ending);

         /* if we do not just print the problem for poly-scip, print the solution to a file */
         if( SCIPspaGetModel(scip) != 'p' )
         {
            /* create the output-name */
            if( SCIPspaGetModel(scip) == 'e' )
            {
               SCIPgetRealParam( scip, "coherence_bound", &eps );
               snprintf(ending, 50, "_eps_%.2f", eps) ;
               strcat(filename, ending );
            }
            else if ( SCIPspaGetModel(scip) == 'w' )
            {
               SCIPgetRealParam( scip, "scale_coherence", &eps );
               snprintf(ending, 50, "_scale_%.4f", eps);
               strcat(filename, ending);
            }
            SCIP_CALL( SCIPgetIntParam(scip, "ncluster", &ncluster) );
            snprintf(ending, 50, "_%dcluster", ncluster);
            strcat(filename, ending);
            strcat(filename, ".sol");

            output = fopen(filename, "w+");
            SCIP_CALL( SCIPprintBestSol(scip, output, FALSE) );
            SCIP_CALL( SCIPprintStatistics(scip, output) );
            fclose(output);
         }
         else
         {
            SCIPinfoMessage(scip, NULL, "Do not solve the problem but write it to file for use in PolySCIP \n");
            SCIPpresolve(scip);
            snprintf(ending, 50, ".mps");
            strcat(filename, ending);
            SCIP_CALL( SCIPwriteOrigProblem(scip, filename, NULL, FALSE) );
         }
      }
      else
      {
         SCIPinfoMessage(scip, NULL, "\n");
         SCIPerrorMessage("Must specify .spa file to be read  \n");
         paramerror = TRUE;
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
SCIP_RETCODE SCIPrunSpa(
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

   /* add the problem-specifix parameters to scip */

   /* include reader, problemdata*/
   SCIP_CALL( SCIPincludeSpaPlugins(scip) );

   /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL( processArguments(scip, argc, argv, defaultsetname) );

   SCIPinfoMessage(scip, NULL, "\n");

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

   retcode = SCIPrunSpa(argc, argv, "scip.set");

   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
