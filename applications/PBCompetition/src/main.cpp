/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

#define PBCOMPETITION

#include <cstdio>
#include <cstring>
#include <ctime>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/scipshell.h"
#include "message_pbscip.h"
#include "event_bestsol.h"
#include "heur_indcoefdiving.h"
#include "heur_indrounding.h"
#include "heur_indoneopt.h"
#include "heur_indtwoopt.h"
#include "heur_smallcard.h"

#define SETOBJ FALSE
#define HEURISTICS_OFF FALSE

#define MAXMEMUSAGE 0.9
#define POSTTIME    3.0                           /**< time saved in the end to free every thing and display the solutions */


static void printUnsupportedSolution(SCIP *scip);

/** sets parameter for satisfiability problems */
static
SCIP_RETCODE loadSettingsNoObjPrePresolve(
      SCIP*                      scip                /**< SCIP data structure */
)
{
   /* TODO set parameters */
   return SCIP_OKAY;
}

/** sets parameter for pure satisfiability problems */
static
SCIP_RETCODE loadSettingsPureSat(
      SCIP*                      scip,               /**< SCIP data structure */
      SCIP_Bool                  hasnoobj            /**< is the problem a satisfiability problem */
)
{
   /* TODO set parameters */
   return SCIP_OKAY;
}

/** sets parameter for wbo instances */
static
SCIP_RETCODE loadSettingsWBO(
      SCIP*                      scip                /**< SCIP data structure */
)
{
   /* TODO set parameters */
   return SCIP_OKAY;
}


static
SCIP_RETCODE readParams(
      SCIP*                      scip,               /**< SCIP data structure */
      const char*                filename            /**< parameter file name */
)
{
   if( SCIPfileExists(filename) )
   {
      SCIPinfoMessage(scip, NULL, "reading user parameter file <%s>\n\n", filename);
      SCIP_CALL( SCIPreadParams(scip, filename) );
      return SCIP_OKAY;
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "user parameter file <%s> not found - using default parameters\n", filename);
   }

   return SCIP_OKAY;
}

#ifdef PBCOMPETITION
static
SCIP_RETCODE printSolution(
   SCIP*                      scip,               /**< SCIP data structure */
   SCIP_VAR**                 vars,               /**< SCIP problem variables */
   int                        nvars,              /**< number of problem variables */
   SCIP_Bool                  hasnoobj            /**< has the problem no original objective function? */
   )
{
   SCIP_MESSAGEHDLR* messagehdlr;
   SCIP_MESSAGEHDLRDATA* messagehdlrdata;
   SCIP_SOL* bestsol;
   SCIP_STATUS status;
   char* solutiontext;
   SCIP_Bool quiet;
   SCIP_Bool comment;
   int strlength;
   int charperline;
   int k;
   int v;
   int printlength;
   int n;

   n = 0;
   printlength = SCIP_MAXSTRLEN/8;

   messagehdlr = SCIPgetMessagehdlr(scip);
   assert(messagehdlr != NULL);
   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);

   /* store old quiet and comment value of message handler */
   quiet = SCIPmessagehdlrIsQuiet(messagehdlr);
   comment = messagehdlrdata->comment;

   /* turn on message handler */
   SCIPmessagehdlrSetQuiet(messagehdlr, FALSE);
   messagehdlrdata->comment = FALSE;

   status = SCIPgetStatus(scip);
   if( status == SCIP_STATUS_INFEASIBLE )
   {
      SCIPinfoMessage(scip, NULL, "s UNSATISFIABLE\n");
   }
   else
   {
      SCIP_Bool feasible;

      feasible = FALSE;

      bestsol = SCIPgetBestSol(scip);

      if( bestsol != NULL )
      {
         /* check if solution is feasible in original problem */
         SCIP_CALL( SCIPcheckSolOrig(scip, bestsol, &feasible, FALSE, FALSE) );

         /* if solution is not feasible in original problem; post UNKNOWN */
         if( !feasible )
         {
            SCIPinfoMessage(scip, NULL, "c internal error\n");
            SCIPinfoMessage(scip, NULL, "s UNKNOWN\n");
         }
         else
         {
            if( hasnoobj || status != SCIP_STATUS_OPTIMAL )
	    {
               SCIPinfoMessage(scip, NULL, "s SATISFIABLE\n");
	    }
            else
	    {
               SCIPinfoMessage(scip, NULL, "s OPTIMUM FOUND\n");
	    }

            charperline = printlength;
            k = 1;

            strlength = SCIP_MAXSTRLEN;
            SCIP_CALL( SCIPallocBufferArray(scip, &solutiontext, strlength) );
            n = sprintf(solutiontext, "v ");

            for( v = nvars - 1; v >= 0; --v )
            {

               assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_ORIGINAL);

               if( strstr(SCIPvarGetName(vars[v]), "andresultant") == NULL && strstr(SCIPvarGetName(vars[v]),
                     "indslack_") == NULL && strstr(SCIPvarGetName(vars[v]), "indicatorvar") == NULL )
               {
                  /* check if enough memory is left */
                  if( n * 1.2 > strlength )
                  {
                     strlength *= 2;
                     SCIP_CALL( SCIPreallocBufferArray(scip, &solutiontext, strlength) );
                  }

                  if( n > charperline )
                  {
                     n += sprintf(&solutiontext[n], "\nv ");
                     ++k;
                     charperline = k * printlength;
                  }

                  if( SCIPgetSolVal(scip, bestsol , vars[v]) > 0.5 )
                     n += sprintf(&solutiontext[n], "%s ", SCIPvarGetName(vars[v]));
                  else
                     n += sprintf(&solutiontext[n], "-%s ", SCIPvarGetName(vars[v]));
               }
            }

            /* if more than only "v " is in the output string, print it */
            if( n > 2 )
            {
               SCIPinfoMessage(scip, NULL, "%s", solutiontext);

               /* flush message buffer */
               SCIPinfoMessage(scip, NULL, "\n");
            }
            SCIPfreeBufferArray(scip, &solutiontext);
         }
      }
      else
      {
         SCIPinfoMessage(scip, NULL, "s UNKNOWN\n");
      }
   }

   /* reset old values of message handler data */
   SCIPmessagehdlrSetQuiet(messagehdlr, quiet);
   messagehdlrdata->comment = comment;

   return SCIP_OKAY;
}
#endif


/** run SCIP from command line */
static
SCIP_RETCODE fromCommandLine(
      SCIP*                      scip,               /**< SCIP data structure */
      const char*                filename,           /**< input file name */
      const char*                settingsfilename,   /**< settings file name */
      SCIP_Real                  timelimit,          /**< required time limit */
      clock_t                    starttime           /**< time the process started */

)
{
   SCIP_RETCODE retcode;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_Bool hasnoobj;
   SCIP_Bool hasindicator;
   SCIP_Bool puresat;
   int nconss;
   int nvars;
   int v;
   int h;
   int c;

   /********************
    * Problem Creation *
    ********************/

   /* set time limit */
   if( timelimit > 0.0 )
   {
      clock_t currtime;

      /* add reading time to solving time */
      SCIP_CALL( SCIPsetBoolParam(scip, "timing/reading", TRUE) );

      /* set time limit */
      currtime = clock();
      timelimit -= POSTTIME;   /* 3 seconds bonus, to finish solving and not get interrupted */
      timelimit -= (SCIP_Real) (currtime-starttime)/(SCIP_Real)CLOCKS_PER_SEC;
      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timelimit) );
   }

   /* use wall clock time */
   SCIP_CALL( SCIPsetIntParam(scip, "timing/clocktype", 2) );

   SCIPinfoMessage(scip, NULL, "reading problem <%s>\n", filename);

   retcode = SCIPreadProb(scip, filename, NULL);

   if( retcode == SCIP_PLUGINNOTFOUND )
   {
      /* try to read the problem as opb format */
      SCIP_RETCODE scipRetcode = SCIPreadProb(scip, filename, "opb");
      if( scipRetcode == SCIP_BIGINT)
      {
         printUnsupportedSolution(scip);
         return SCIP_OKAY;
      }
      SCIP_CALL(scipRetcode);
   }
   else
   {
      if( retcode == SCIP_BIGINT)
      {
         printUnsupportedSolution(scip);
         return SCIP_OKAY;
      }
      SCIP_CALL( retcode );
   }

   SCIPinfoMessage(scip, NULL, "problem read in %.2f\n", SCIPgetReadingTime(scip));

   /*******************
    * Problem Solving *
    *******************/

   vars = SCIPgetOrigVars(scip);
   nvars = SCIPgetNOrigVars(scip);
   assert(vars != NULL || nvars == 0);
   hasnoobj = TRUE;
   for( v = nvars - 1; v >= 0 && hasnoobj; --v )
   {
      assert(vars != NULL);

      if( !SCIPisZero(scip, SCIPvarGetObj(vars[v])) )
      {
         hasnoobj = FALSE;
         break;
      }
   }

   conss = SCIPgetOrigConss(scip);
   nconss = SCIPgetNOrigConss(scip);
   hasindicator = FALSE;
   for (c = 0; c < nconss; ++c)
   {
      conshdlr = SCIPconsGetHdlr(conss[c]);
      if ( strcmp(SCIPconshdlrGetName(conshdlr), "pseudoboolean") == 0 )
      {
         if ( SCIPgetIndVarPseudoboolean(scip, conss[c]) != NULL )
         {
            hasindicator = TRUE;
            break;
         }
      }
   }

#ifdef PBCOMPETITION
   /* create event handler for best solution found if an objective exists  */
   if( !hasnoobj )
   {
      SCIP_CALL( SCIPcreateEventHdlrBestsol(scip) );
   }
#endif

   if( hasnoobj )
   {
      SCIPinfoMessage(scip, NULL, "No objective function, only one solution is needed.\n");

      if( settingsfilename == NULL )
      {
         /* load settings for all SATUNSAT problems */
         SCIP_CALL( loadSettingsNoObjPrePresolve(scip) );
         SCIPinfoMessage(scip, NULL, "presolving settings loaded\n");
      }

      /*  insert objective function if no exists */
      if( SETOBJ ) /*lint !e774 !e506*/
      {
         for( v = nvars - 1; v >= 0; --v )
         {
            assert(vars != NULL);

            /* only change objective value of variable which wouldn't be fix by dual presolver */
            if( SCIPvarGetNLocksUp(vars[v]) > 0 )
            {
               SCIP_CALL( SCIPchgVarObj(scip, vars[v], 1.0) );
            }
         }
      }
   }

   if ( hasindicator )
   {
      /* load settings for WBO instances */
      SCIP_CALL( loadSettingsWBO(scip) );
   }
   /* start presolving */
   SCIP_CALL( SCIPpresolve(scip) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPwriteTransProblem(scip, "debug.lp", "lp", FALSE) );
   SCIP_CALL( SCIPwriteTransProblem(scip, "debug.cip", "cip", FALSE) );
#endif

   /* determine problem typ */
   puresat = FALSE;
   conshdlr = SCIPfindConshdlr(scip, "logicor");
   if( conshdlr != NULL && SCIPconshdlrGetNActiveConss(conshdlr) == SCIPgetNConss(scip) )
      puresat = TRUE;

#ifdef PBCOMPETITION
   // nonlinear = FALSE;
   conshdlr = SCIPfindConshdlr(scip, "and");
   // if( conshdlr != NULL && SCIPconshdlrGetNActiveConss(conshdlr) > 0 )
   //    nonlinear = TRUE;
#endif

   /* set setting for the branch-and-bound process */
   if( settingsfilename == NULL )
   {
      if( puresat )
         SCIP_CALL( loadSettingsPureSat(scip, hasnoobj) );
      else
         /* activate rapidlearning */
         SCIP_CALL( SCIPsetIntParam(scip, "separating/rapidlearning/freq", 0) );

      if( hasnoobj )
      {
         /* turn all heuristics off if desired */
         if( HEURISTICS_OFF ) /*lint !e774 !e506*/
         {
            SCIP_HEUR** heuristics;
            char parametername[SCIP_MAXSTRLEN];
            int nheuristics;

            /* turn off all heuristics */
            nheuristics = SCIPgetNHeurs(scip);
            heuristics = SCIPgetHeurs(scip);

            for( h = nheuristics - 1; h >= 0; --h )
            {
               sprintf(parametername, "heuristics/%s/freq", SCIPheurGetName(heuristics[h]));

               if( SCIPheurGetFreq(heuristics[h]) != -1 )
               {
                  SCIP_CALL( SCIPsetIntParam(scip, parametername, -1) );
               }
               assert( SCIPheurGetFreq(heuristics[h]) == -1 );
            }
         }
      }
   }

   /* in case no LP solver is available turn off the LP relaxation */
   if( strcmp(SCIPlpiGetSolverName(), "NONE") == 0 )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) );
      SCIP_CALL( SCIPsetBoolParam(scip, "constraints/indicator/sepaAlternativeLP", FALSE) );
   }

   /* write all non-default parameters in file */
   SCIPinfoMessage(scip, NULL, "- non default parameters ----------------------------------------------------------------------\n");
   SCIP_CALL( SCIPwriteParams( scip, NULL, TRUE, TRUE) );
   SCIPinfoMessage(scip, NULL, "-----------------------------------------------------------------------------------------------\n");

   /* start solving */
   SCIPinfoMessage(scip, NULL, "start solving\n");
   SCIP_CALL( SCIPsolve(scip) );

   /* printing solution */
#ifdef PBCOMPETITION
   SCIP_CALL( printSolution(scip, vars, nvars, hasnoobj));
#else
   SCIP_CALL( SCIPprintBestSol( scip, NULL, FALSE) );
#endif

   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   return SCIP_OKAY;
}

static void printUnsupportedSolution(SCIP *scip) {
   SCIP_MESSAGEHDLR* messagehdlr;
   SCIP_MESSAGEHDLRDATA* messagehdlrdata;

   messagehdlr = SCIPgetMessagehdlr(scip);
   assert(messagehdlr != NULL);
   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);

   /* turn on message handler */
   SCIPmessagehdlrSetQuiet(messagehdlr, FALSE);
   messagehdlrdata->comment = FALSE;
   SCIPinfoMessage(scip, NULL, "s UNSUPPORTED\n");
}



/** evaluates command line parameters and runs SCIP appropriately in the given SCIP instance */
static
SCIP_RETCODE processShellArguments(
      SCIP*                      scip,               /**< SCIP data structure */
      SCIP_MESSAGEHDLR*          defaultmessagehdlr, /**< standard message handler */
      int                        argc,               /**< number of shell parameters */
      char**                     argv,               /**< array with shell parameters */
      clock_t                    starttime,          /**< time the process started */
      const char*                defaultsetname      /**< name of default settings file */
)
{
   char* probname = NULL;
   char* settingsname = NULL;
   char* logname = NULL;
   SCIP_Bool quiet;
   SCIP_Bool paramerror;
   SCIP_Bool interactive;
   SCIP_Real timelimit;
   int i;

   /********************
    * Parse parameters *
    ********************/

   quiet = FALSE;
   paramerror = FALSE;
   timelimit = 0.0;

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
            SCIPerrorMessage("missing log filename after parameter '-l'\n");
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
            SCIPerrorMessage("missing settings filename after parameter '-s'\n");
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
            SCIPerrorMessage("missing problem filename after parameter '-f'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-t") == 0 )
      {
         i++;
         if( i < argc )
            /* get time limit */
            timelimit = atof(argv[i]);
         else
         {
            SCIPerrorMessage("missing amout of time '-t'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-m") == 0 )
      {
         i++;
         if( i < argc )
         {
            /* set time limit */
            SCIP_CALL( SCIPsetRealParam(scip, "limits/memory", atof(argv[i]) * MAXMEMUSAGE) );
         }
         else
         {
            SCIPerrorMessage("missing amout of memory '-m'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-d") == 0 )
      {
         i++;
         if( i < argc )
         {
            /* set display frequency limit */
            SCIP_CALL( SCIPsetIntParam(scip, "display/freq", atoi(argv[i])) );
         }
         else
         {
            SCIPerrorMessage("missing display frequency '-d'\n");
            paramerror = TRUE;
         }
      }
#ifndef PBCOMPETITION
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
            SCIPerrorMessage("missing command line after parameter '-c'\n");
            paramerror = TRUE;
         }
      }
#endif
      else
      {
         SCIPerrorMessage("invalid parameter <%s>\n", argv[i]);
         paramerror = TRUE;
      }
   }

   if( !paramerror )
   {
      /* change necessary flags of message handler */
      SCIPsetMessagehdlrQuiet(scip, quiet);
      SCIPsetMessagehdlrLogfile(scip, logname);

      if( !interactive )
      {
         /***********************
          * Version information *
          ***********************/
         if( probname == NULL )
         {
            SCIP_CALL( SCIPsetMessagehdlr(scip, defaultmessagehdlr) );
            SCIPmessageSetErrorPrintingDefault();
         }

         SCIPprintVersion(scip, NULL);
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
            SCIP_CALL( fromCommandLine(scip, probname, settingsname, timelimit, starttime) );
         }
         else
         {
            SCIP_CALL( SCIPstartInteraction(scip) );
         }
      }
      else
      {
         SCIPprintVersion(scip, NULL);
         SCIPinfoMessage(scip, NULL, "\n");

         SCIP_CALL( SCIPstartInteraction(scip) );
      }
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "syntax: %s [-l <logfile>] [-q] [-s <settings>] [-f <problem>] [-b <batchfile>]", argv[0]);
#ifndef PBCOMPETITION
      SCIPinfoMessage(scip, NULL, " [-c \"command\"]");
#endif
      SCIPinfoMessage(scip, NULL, " [-t <timelimit>] [-m <memlimit>] [-d <dispfreq>]\n");

      SCIPinfoMessage(scip, NULL, "  -l <logfile>   : copy output into log file\n");
      SCIPinfoMessage(scip, NULL, "  -q             : suppress screen messages\n");
      SCIPinfoMessage(scip, NULL, "  -s <settings>  : load parameter settings (.set) file\n");
      SCIPinfoMessage(scip, NULL, "  -f <problem>   : load and solve problem file\n");
      SCIPinfoMessage(scip, NULL, "  -d <batchfile> : load batchfile\n");
#ifndef PBCOMPETITION
      SCIPinfoMessage(scip, NULL, "  -c <command>   : execute command on command line\n");
#endif
      SCIPinfoMessage(scip, NULL, "  -t <timelimit> : sets the time limit\n");
      SCIPinfoMessage(scip, NULL, "  -m <memlimit>  : sets the memory limit\n");
      SCIPinfoMessage(scip, NULL, "  -d <dispfreq>  : sets display frequency\n");
   }

   return SCIP_OKAY;
}


/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runShell(
      int                        argc,               /**< number of shell parameters */
      char**                     argv,               /**< array with shell parameters */
      clock_t                    starttime,          /**< time the process started */
      const char*                defaultsetname      /**< name of default settings file */
)
{
   SCIP_MESSAGEHDLR* messagehdlr;
   SCIP_MESSAGEHDLR* defaultmessagehdlr;
   SCIP* scip;

   messagehdlr = NULL;
   scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* add heuristics for indicators */
   SCIP_CALL( SCIPincludeHeurIndcoefdiving(scip) );
   SCIP_CALL( SCIPincludeHeurIndrounding(scip) );
   SCIP_CALL( SCIPincludeHeurIndoneopt(scip) );
   SCIP_CALL( SCIPincludeHeurIndtwoopt(scip) );
   SCIP_CALL( SCIPincludeHeurSmallcard(scip) );

   defaultmessagehdlr = SCIPgetMessagehdlr(scip);
   /* capture default message handler */
   SCIPmessagehdlrCapture(defaultmessagehdlr);

   /* create own buffered message handler for PB output and overrides default scip message handler */
   SCIP_CALL( SCIPcreateMessagehdlrPbscip(&messagehdlr, TRUE, NULL, FALSE) );

   /* if we are running the PB competition we use this message handler to produce the required output */
   SCIP_CALL( SCIPsetMessagehdlr(scip, messagehdlr) );

   /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL( processShellArguments(scip, defaultmessagehdlr, argc, argv, starttime, defaultsetname) );

   /* release captured and own message handler */
   SCIP_CALL( SCIPmessagehdlrRelease(&messagehdlr) );
   SCIP_CALL( SCIPmessagehdlrRelease(&defaultmessagehdlr) );

   /********************
    * Deinitialization *
    ********************/
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}



int main(
   int                   argc,          /**< number of arguments */
   char**                argv           /**< string array with arguments */
   )
{
   SCIP_RETCODE retcode;

   clock_t starttime, endtime;

   starttime = clock();

   retcode = runShell(argc, argv, starttime, "scip.set");

   endtime = clock();

   printf("c Time complete: %g.\n",(double) (endtime-starttime)/(double)CLOCKS_PER_SEC);

   if( retcode != SCIP_OKAY )
      printf("s UNKNOWN \n");
  return 0;
}
