/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scipshell.c
 * @brief  SCIP command line interface
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/scipshell.h"

/*
 * Message Handler
 */

/** message handler data */
struct SCIP_MessagehdlrData
{
   FILE*                 logfile;            /**< log file where to copy messages into */
   SCIP_Bool             quiet;              /**< should screen messages be suppressed? */
};

/** prints a message to the given file stream and writes the same messate to the log file */
static
void logMessage(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print message into */
   const char*           msg                 /**< message to print */
   )
{
   SCIP_MESSAGEHDLRDATA* messagehdlrdata;

   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);

   if( !messagehdlrdata->quiet || (file != stdout && file != stderr) )
   {
      fputs(msg, file);
      fflush(file);
   }
   if( messagehdlrdata->logfile != NULL && (file == stdout || file == stderr) )
   {
      fputs(msg, messagehdlrdata->logfile);
      fflush(messagehdlrdata->logfile);
   }
}

/** error message print method of message handler */
static
SCIP_DECL_MESSAGEERROR(messageErrorLog)
{
   logMessage(messagehdlr, file, msg);
}

/** warning message print method of message handler */
static
SCIP_DECL_MESSAGEWARNING(messageWarningLog)
{
   logMessage(messagehdlr, file, msg);
}

/** dialog message print method of message handler */
static
SCIP_DECL_MESSAGEDIALOG(messageDialogLog)
{
   logMessage(messagehdlr, file, msg);
}

/** info message print method of message handler */
static
SCIP_DECL_MESSAGEINFO(messageInfoLog)
{
   logMessage(messagehdlr, file, msg);
}



static
SCIP_RETCODE readParams(
   SCIP*                      scip,               /**< SCIP data structure */
   const char*                filename            /**< parameter file name */
   )
{
   if( SCIPfileExists(filename) )
   {
      SCIPinfoMessage(scip, NULL, "reading user parameter file <%s>\n", filename);
      SCIP_CALL( SCIPreadParams(scip, filename) );
   }
   else
      SCIPinfoMessage(scip, NULL, "user parameter file <%s> not found - using default parameters\n", filename);

   return SCIP_OKAY;
}

static
SCIP_RETCODE fromCommandLine(
   SCIP*                      scip,               /**< SCIP data structure */
   const char*                filename            /**< input file name */
   )
{
   /********************
    * Problem Creation *
    ********************/

   SCIPinfoMessage(scip, NULL, "\nread problem <%s>\n", filename);
   SCIPinfoMessage(scip, NULL, "============\n\n");
   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );


   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   SCIPinfoMessage(scip, NULL, "\nsolve problem\n");
   SCIPinfoMessage(scip, NULL, "=============\n\n");
   SCIP_CALL( SCIPsolve(scip) );

   SCIPinfoMessage(scip, NULL, "\nprimal solution:\n");
   SCIPinfoMessage(scip, NULL, "================\n\n");
   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );


   /**************
    * Statistics *
    **************/

   SCIPinfoMessage(scip, NULL, "\nStatistics\n");
   SCIPinfoMessage(scip, NULL, "==========\n\n");

   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   return SCIP_OKAY;
}

/** evaluates command line parameters and runs SCIP appropriately in the given SCIP instance */
SCIP_RETCODE SCIPprocessShellArguments(
   SCIP*                      scip,               /**< SCIP data structure */
   int                        argc,               /**< number of shell parameters */
   char**                     argv,               /**< array with shell parameters */
   const char*                defaultsetname      /**< name of default settings file */
   )
{
   char* probname = NULL;
   char* settingsname = NULL;
   char* logname = NULL;
   SCIP_Bool quiet;
   SCIP_Bool paramerror;
   SCIP_Bool interactive;
   int i;

   /********************
    * Parse parameters *
    ********************/
   
   quiet = FALSE;
   paramerror = FALSE;
   interactive = FALSE;
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
            probname = argv[i];
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

                  (void)SCIPfgets(buffer, sizeof(buffer), file);
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
      else
      {
         printf("invalid parameter <%s>\n", argv[i]);
         paramerror = TRUE;
      }
   }
   if( interactive && probname != NULL )
   {
      printf("cannot mix batch mode '-c' and '-b' with file mode '-f'\n");
      paramerror = TRUE;
   }

   if( !paramerror )
   {
      SCIP_MESSAGEHDLR* messagehdlr;
      SCIP_MESSAGEHDLRDATA* messagehdlrdata;
      SCIP_Bool error;

      /***********************************
       * create log file message handler *
       ***********************************/

      messagehdlr = NULL;
      messagehdlrdata = NULL;
      error = FALSE;
      if( logname != NULL || quiet )
      {
         SCIP_CALL( SCIPallocMemory(scip, &messagehdlrdata) );
         if( logname != NULL )
         {
            messagehdlrdata->logfile = fopen(logname, "a"); /* append to log file */
            if( messagehdlrdata->logfile == NULL )
            {
               SCIPerrorMessage("cannot open log file <%s> for writing\n", logname);
               error = TRUE;
            }
         }
         else
            messagehdlrdata->logfile = NULL;
         messagehdlrdata->quiet = quiet;
         SCIP_CALL( SCIPcreateMessagehdlr(&messagehdlr, FALSE, 
               messageErrorLog, messageWarningLog, messageDialogLog, messageInfoLog,
               messagehdlrdata) );
         SCIP_CALL( SCIPsetMessagehdlr(messagehdlr) );
      }

      if( !error )
      {
         /***********************************
          * Version and library information *
          ***********************************/
         
         SCIPprintVersion(NULL);
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
            SCIP_CALL( fromCommandLine(scip, probname) );
         }
         else
         {
            SCIPinfoMessage(scip, NULL, "\n");
            SCIP_CALL( SCIPstartInteraction(scip) );
         }

         /******************
          * Close log file *
          ******************/

         if( messagehdlrdata != NULL )
         {
            SCIP_CALL( SCIPsetDefaultMessagehdlr() );
            SCIP_CALL( SCIPfreeMessagehdlr(&messagehdlr) );
            if( messagehdlrdata->logfile != NULL )
               fclose(messagehdlrdata->logfile);
            SCIPfreeMemory(scip, &messagehdlrdata);
         }
      }
   }
   else
   {
      printf("\nsyntax: %s [-l <logfile>] [-q] [-s <settings>] [-f <problem>] [-b <batchfile>] [-c \"command\"]\n"
         "  -l <logfile>  : copy output into log file\n"
         "  -q            : suppress screen messages\n"
         "  -s <settings> : load parameter settings (.set) file\n"
         "  -f <problem>  : load and solve problem file\n"
         "  -b <batchfile>: load and execute dialog command batch file (can be used multiple times)\n"
         "  -c \"command\"  : execute single line of dialog commands (can be used multiple times)\n\n",
         argv[0]);
   }

   return SCIP_OKAY;
}

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
SCIP_RETCODE SCIPrunShell(
   int                        argc,               /**< number of shell parameters */
   char**                     argv,               /**< array with shell parameters */
   const char*                defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );


   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}
