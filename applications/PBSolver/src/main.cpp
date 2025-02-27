/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   main.cpp
 * @brief  main file for the Pseudo-Boolean solver application
 * @author Alexander Hoen
 * @author Gioni Mexi
 *
 * @todo Add absolute feasibility tolerance parameter to enforce exact feasibility.
 * @todo Add separate integrality tolerance parameter to enforce exact integrality.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cstdio>
#include <cstring>
#include <ctime>
#include "scip/scipdefplugins.h"
#include "scip/scipshell.h"
#include "message_pb.h"
#include "event_bestsol.h"

#define SETOBJ         FALSE                 /**< insert objective function if no exists */
#define HEURISTICS_OFF FALSE                 /**< turn off heuristics */
#define MAXMEMUSAGE    0.9                   /**< maximal memory usage relative to given memory limit */
#define POSTTIME       3.0                   /**< time in seconds saved in the end to display solution and free everything */


/** sets parameters for pure satisfiability problems */
static
SCIP_RETCODE loadSettingsPureSat(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* TODO: set parameters */
   return SCIP_OKAY;
}

/** sets parameters for maximum satisfiability instances */
static
SCIP_RETCODE loadSettingsMaxSAT(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* TODO: set parameters */
   return SCIP_OKAY;
}

/** prints that the problem instance is unsupported */
static
SCIP_RETCODE printUnsupported(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_MESSAGEHDLR* messagehdlr = SCIPgetMessagehdlr(scip);
   assert(messagehdlr != NULL);
   SCIP_MESSAGEHDLRDATA* messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);
   SCIP_Bool quiet = SCIPmessagehdlrIsQuiet(messagehdlr);

   /* competition console log */
   if( messagehdlrdata->comment )
   {
      /* turn on message handler and remove comment declaration */
      SCIPmessagehdlrSetQuiet(messagehdlr, FALSE);
      messagehdlrdata->comment = FALSE;

      SCIPinfoMessage(scip, NULL, "s UNSUPPORTED\n");

      /* reset old values of message handler data */
      SCIPmessagehdlrSetQuiet(messagehdlr, quiet);
      messagehdlrdata->comment = TRUE;
   }

   return SCIP_OKAY;
}

/** prints the best primal solution */
static
SCIP_RETCODE printSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< SCIP problem variables */
   int                   nvars,              /**< number of problem variables */
   SCIP_Bool             hasobj              /**< does the problem have an objective function? */
   )
{
   SCIP_MESSAGEHDLR* messagehdlr = SCIPgetMessagehdlr(scip);
   assert(messagehdlr != NULL);
   SCIP_MESSAGEHDLRDATA* messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);
   SCIP_SOL* bestsol;
   SCIP_STATUS status;
   SCIP_Bool quiet = SCIPmessagehdlrIsQuiet(messagehdlr);
   SCIP_Bool feasible;
   char* solutiontext;
   int strlength;
   int printlength;
   int n;
   int v;

   /* enable message handler */
   SCIPmessagehdlrSetQuiet(messagehdlr, FALSE);

   /* competition console log */
   if( messagehdlrdata->comment )
   {
      /* disable comment prefix */
      messagehdlrdata->comment = FALSE;

      status = SCIPgetStatus(scip);
      if( status == SCIP_STATUS_INFEASIBLE )
         SCIPinfoMessage(scip, NULL, "s UNSATISFIABLE\n");
      else
      {
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
               if( !hasobj || status != SCIP_STATUS_OPTIMAL )
                  SCIPinfoMessage(scip, NULL, "s SATISFIABLE\n");
               else
                  SCIPinfoMessage(scip, NULL, "s OPTIMUM FOUND\n");

               strlength = 2 * SCIP_MAXSTRLEN;
               printlength = SCIP_MAXSTRLEN / 8;
               n = 0;
               v = 0;
               SCIP_CALL( SCIPallocBufferArray(scip, &solutiontext, strlength) );

               while( v < nvars )
               {
                  int printed = 0;

                  assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_ORIGINAL);

                  if( strstr(SCIPvarGetName(vars[v]), "andresultant_") == NULL
                     && strstr(SCIPvarGetName(vars[v]), "indslack_") == NULL
                     && strstr(SCIPvarGetName(vars[v]), "indicatorvar_") == NULL )
                  {
                     printed = SCIPsnprintf(solutiontext + n, strlength,
                           SCIPgetSolVal(scip, bestsol, vars[v]) > 0.5 ? " %s" : " -%s", SCIPvarGetName(vars[v]));
                     n += printed;
                     strlength -= printed;
                  }

                  if( n >= printlength || ( n >= 1 && v + 1 >= nvars ) )
                  {
                     if( strlength >= 1 )
                     {
                        strlength += n;
                        n = 0;
                        SCIPinfoMessage(scip, NULL, "v%s\n", solutiontext);
                     }
                     else
                     {
                        strlength += printed + SCIP_MAXSTRLEN;
                        n -= printed;
                        SCIP_CALL( SCIPreallocBufferArray(scip, &solutiontext, n + strlength) );

                        continue;
                     }
                  }
                  assert(strlength >= 1);

                  ++v;
               }

               SCIPfreeBufferArray(scip, &solutiontext);
            }
         }
         else
            SCIPinfoMessage(scip, NULL, "s UNKNOWN\n");
      }

      /* enable comment prefix */
      messagehdlrdata->comment = TRUE;
   }
   /* default console log */
   else
   {
      SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );
   }

   /* reset message handler */
   SCIPmessagehdlrSetQuiet(messagehdlr, quiet);

   return SCIP_OKAY;
}

/** run SCIP from command line */
static
SCIP_RETCODE fromCommandLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           settingsfilename,   /**< settings file name */
   const char*           problemfilename,    /**< problem file name */
   SCIP_Real             timelimit,          /**< required time limit */
   clock_t               startclock          /**< clock at which the process started */
   )
{
   SCIP_RETCODE retcode;
   SCIP_VAR** vars;
   SCIP_CONS** conshdlrconss;
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool hasobj;
   SCIP_Bool hasindicator;
   SCIP_Bool puresat;
   char* filenamecopy;
   char* extension;
   int nvars;
   int nconshdlrconss;
   int npuresatconss;
   int v;
   int c;

   /***********************
    * Settings Definition *
    ***********************/

   /* read parameter settings */
   if( settingsfilename != NULL )
   {
      SCIP_CALL( SCIPreadParams(scip, settingsfilename) );
   }

   /* define time limit */
   if( timelimit >= 0.0 )
   {
      /* get starting time and reserve finishing time */
      timelimit -= (SCIP_Real)(clock() - startclock) / (SCIP_Real)CLOCKS_PER_SEC + POSTTIME;

      /* stop immediately if time exceeded */
      if( timelimit < 0.0 )
         return SCIP_INVALIDCALL;

      /* set time limit */
      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timelimit) );
   }

   /* add reading time */
   SCIP_CALL( SCIPsetBoolParam(scip, "timing/reading", TRUE) );

   /* use wall clock */
   SCIP_CALL( SCIPsetIntParam(scip, "timing/clocktype", 2) );

   /********************
    * Problem Creation *
    ********************/

   /* read pseudoboolean problem */
   SCIPinfoMessage(scip, NULL, "reading problem <%s>\n", problemfilename);
   SCIP_CALL( SCIPduplicateBufferArray(scip, &filenamecopy, problemfilename, (int)strlen(problemfilename) + 1) );
   SCIPsplitFilename(filenamecopy, NULL, NULL, &extension, NULL);
   if( strcmp(extension, "opb") == 0 || strcmp(extension, "wbo") == 0 )
      retcode = SCIPreadProb(scip, problemfilename, extension);
   else
      retcode = SCIP_INVALIDDATA;
   SCIPfreeBufferArray(scip, &filenamecopy);

   /* declare unsupported problem */
   if( retcode == SCIP_INVALIDDATA )
   {
      SCIP_CALL( printUnsupported(scip) );
      return SCIP_OKAY;
   }

   /* detect unexpected error */
   SCIP_CALL( retcode );
   SCIPinfoMessage(scip, NULL, "problem read in %.3lf seconds\n", SCIPgetReadingTime(scip));

   /*******************
    * Problem Solving *
    *******************/

   vars = SCIPgetOrigVars(scip);
   nvars = SCIPgetNOrigVars(scip);
   assert(vars != NULL || nvars == 0);
   hasobj = FALSE;
   for( v = 0; v < nvars; ++v )
   {
      if( !SCIPisZero(scip, SCIPvarGetObj(vars[v])) )
      {
         hasobj = TRUE;
         break;
      }
   }

   /* detect soft constraints */
   conshdlr = SCIPfindConshdlr(scip, "pseudoboolean");
   conshdlrconss = SCIPconshdlrGetConss(conshdlr);
   nconshdlrconss = SCIPconshdlrGetNConss(conshdlr);
   assert(conshdlrconss != NULL || nconshdlrconss == 0);
   hasindicator = FALSE;
   for( c = 0; c < nconshdlrconss; ++c )
   {
      if( SCIPgetIndVarPseudoboolean(scip, conshdlrconss[c]) != NULL )
      {
         hasindicator = TRUE;
         break;
      }
   }

   /* create event handler for best solution found if an objective exists */
   if( hasobj )
   {
      if( (SCIPmessagehdlrGetData(SCIPgetMessagehdlr(scip)))->comment )
      {
         SCIP_CALL( SCIPcreateEventHdlrBestsol(scip) );
      }
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "problem without objective\n");
#if SETOBJ
      /* insert objective function if no exists */
      for( v = 0; v < nvars; ++v )
      {
         /* add objective coefficient if variable will not be fixed up by dual presolver */
         if( SCIPvarGetNLocksUp(vars[v]) >= 1 )
         {
            SCIP_CALL( SCIPchgVarObj(scip, vars[v], 1.0) );
         }
      }
#endif
   }

   if( hasindicator && settingsfilename == NULL )
   {
      /* load settings for maximum SAT */
      SCIP_CALL( loadSettingsMaxSAT(scip) );
   }

   /* start presolving */
   SCIP_CALL( SCIPpresolve(scip) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPwriteTransProblem(scip, "debug.cip", NULL, FALSE) );
#endif

   /* count number of constraints represented by or conditions */
   npuresatconss = 0;

   /* add logicor constraints */
   conshdlr = SCIPfindConshdlr(scip, "logicor");
   if( conshdlr != NULL )
      npuresatconss += SCIPconshdlrGetNCheckConss(conshdlr);

   /* add and constraints as representable by linear amount of logicor constraints */
   conshdlr = SCIPfindConshdlr(scip, "and");
   if( conshdlr != NULL )
      npuresatconss += SCIPconshdlrGetNCheckConss(conshdlr);

   /* add setppc constraints if all are covering */
   conshdlr = SCIPfindConshdlr(scip, "setppc");
   if( conshdlr != NULL )
   {
      nconshdlrconss = SCIPconshdlrGetNCheckConss(conshdlr);
      if( npuresatconss + nconshdlrconss == SCIPgetNCheckConss(scip) )
      {
         conshdlrconss = SCIPconshdlrGetCheckConss(conshdlr);
         assert(conshdlrconss != NULL || nconshdlrconss == 0);
         for( c = 0; c < nconshdlrconss; ++c )
         {
            if( SCIPgetTypeSetppc(scip, conshdlrconss[c]) != SCIP_SETPPCTYPE_COVERING )
               break;
         }
         if( c == nconshdlrconss )
            npuresatconss += nconshdlrconss;
      }
   }

   /* determine problem type */
   puresat = (npuresatconss == SCIPgetNCheckConss(scip));
   if( puresat )
      SCIPinfoMessage(scip, NULL, "problem pure SAT\n");

   /* set setting for the branch-and-bound process */
   if( settingsfilename == NULL )
   {
      if( puresat )
      {
         /* load settings for pure SAT */
         SCIP_CALL( loadSettingsPureSat(scip) );
      }
      else
      {
         /* activate rapidlearning */
         SCIP_CALL( SCIPsetIntParam(scip, "separating/rapidlearning/freq", 0) );
      }

#if HEURISTICS_OFF
      /* turn off heuristics */
      char parametername[SCIP_MAXSTRLEN];
      SCIP_HEUR** heuristics = SCIPgetHeurs(scip);
      int nheuristics = SCIPgetNHeurs(scip);
      int h;
      assert(heuristics != NULL);

      for( h = 0; h < nheuristics; ++h )
      {
         (void)SCIPsnprintf(parametername, SCIP_MAXSTRLEN, "heuristics/%s/freq", SCIPheurGetName(heuristics[h]));
         SCIP_CALL( SCIPsetIntParam(scip, parametername, -1) );
         assert(SCIPheurGetFreq(heuristics[h]) == -1);
      }
#endif
   }

   /* write non-default parameters to console */
   SCIPinfoMessage(scip, NULL, "\n- non default parameters ----------------------------------------------------------------------\n\n");
   SCIP_CALL( SCIPwriteParams( scip, NULL, TRUE, TRUE) );
   SCIPinfoMessage(scip, NULL, "-----------------------------------------------------------------------------------------------\n");

   /* start solving */
   SCIP_CALL( SCIPsolve(scip) );

   /* print resulting solution */
   SCIP_CALL( printSolution(scip, vars, nvars, hasobj) );

   /* print solving statistics */
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   return SCIP_OKAY;
}

/** evaluates command line parameters and runs SCIP appropriately in the given SCIP instance */
static
SCIP_RETCODE processShellArguments(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   clock_t               startclock,         /**< clock at which the process started */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
   SCIP_MESSAGEHDLR* messagehdlr;
   SCIP_MESSAGEHDLRDATA* messagehdlrdata;
   SCIP_Real timelimit = -1.0;
   SCIP_Bool quiet = FALSE;
   SCIP_Bool print = FALSE;
   SCIP_Bool paramerror = FALSE;
   SCIP_Bool interactive = FALSE;
   char* logname = NULL;
   char* settingsname = NULL;
   char* problemname = NULL;
   int i;

   /********************
    * Parse parameters *
    ********************/

   for( i = 1; i < argc; ++i )
   {
      /* check for a valid flag */
      if( argv[i][0] == '-' && strlen(argv[i]) == 2 )
      {
         switch( argv[i][1] )
         {
            /* set quiet flag */
            case 'q':
               quiet = TRUE;
               break;
            /* set print flag */
            case 'p':
               print = TRUE;
               break;
            /* get log filename */
            case 'l':
               if( ++i < argc )
                  logname = argv[i];
               else
               {
                  SCIPerrorMessage("missing log filename after parameter '-l'\n");
                  paramerror = TRUE;
               }
               break;
            /* get settings filename */
            case 's':
               if( ++i < argc )
                  settingsname = argv[i];
               else
               {
                  SCIPerrorMessage("missing settings filename after parameter '-s'\n");
                  paramerror = TRUE;
               }
               break;
            /* get problem filename */
            case 'f':
               if( ++i < argc )
                  problemname = argv[i];
               else
               {
                  SCIPerrorMessage("missing problem filename after parameter '-f'\n");
                  paramerror = TRUE;
               }
               break;
            /* set display frequency */
            case 'd':
               if( ++i < argc )
               {
                  SCIP_CALL( SCIPsetIntParam(scip, "display/freq", atoi(argv[i])) );
               }
               else
               {
                  SCIPerrorMessage("missing display frequency after parameter '-d'\n");
                  paramerror = TRUE;
               }
               break;
            /* get time limit */
            case 't':
               if( ++i < argc )
                  timelimit = atof(argv[i]);
               else
               {
                  SCIPerrorMessage("missing time limit after parameter '-t'\n");
                  paramerror = TRUE;
               }
               break;
            /* set memory limit */
            case 'm':
               if( ++i < argc )
               {
                  SCIP_CALL( SCIPsetRealParam(scip, "limits/memory", atof(argv[i]) * MAXMEMUSAGE) );
               }
               else
               {
                  SCIPerrorMessage("missing memory limit after parameter '-m'\n");
                  paramerror = TRUE;
               }
               break;
            /* command line input */
            case 'c':
               if( ++i < argc )
               {
                  SCIP_CALL( SCIPaddDialogInputLine(scip, argv[i]) );
                  interactive = TRUE;
               }
               else
               {
                  SCIPerrorMessage("missing command line after parameter '-c'\n");
                  paramerror = TRUE;
               }
               break;
            default:
               SCIPerrorMessage("invalid parameter '%s'\n", argv[i]);
               paramerror = TRUE;
               break;
         }
      }
      else
      {
         SCIPerrorMessage("invalid parameter '%s'\n", argv[i]);
         paramerror = TRUE;
      }
   }
   if( !paramerror )
   {
      /***********************
       * Version information *
       ***********************/

      /* set attributes of message handler */
      messagehdlr = SCIPgetMessagehdlr(scip);
      assert(messagehdlr != NULL);
      messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
      assert(messagehdlrdata != NULL);
      SCIPsetMessagehdlrQuiet(scip, quiet);
      SCIPsetMessagehdlrLogfile(scip, logname);
      messagehdlrdata->comment = !print;

      /* print version information */
      SCIPprintVersion(scip, NULL);
      SCIPinfoMessage(scip, NULL, "\n");
      SCIPprintExternalCodes(scip, NULL);
      SCIPinfoMessage(scip, NULL, "\n");

      /*****************
       * Load defaults *
       *****************/

      if( defaultsetname != NULL && SCIPfileExists(defaultsetname) )
      {
         SCIP_CALL( SCIPreadParams(scip, defaultsetname) );
      }

      if( !interactive && problemname != NULL )
      {
         /**************
          * Start SCIP *
          **************/

         SCIP_CALL( fromCommandLine(scip, settingsname, problemname, timelimit, startclock) );
      }
      else
      {
         SCIP_CALL( SCIPstartInteraction(scip) );
      }
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "syntax: %s [-q] [-p] [-l <logfile>] [-s <settings>] [-f <problem>] [-d <dispfreq>] [-t <timelimit>] [-m <memlimit>]\n", argv[0]);

      SCIPinfoMessage(scip, NULL, "   -q             : suppress screen messages\n");
      SCIPinfoMessage(scip, NULL, "   -p             : print default messages\n");
      SCIPinfoMessage(scip, NULL, "   -l <logfile>   : copy output into log file\n");
      SCIPinfoMessage(scip, NULL, "   -s <settings>  : load settings (.set) file\n");
      SCIPinfoMessage(scip, NULL, "   -f <problem>   : solve problem (.opb or .wbo) file\n");
      SCIPinfoMessage(scip, NULL, "   -d <dispfreq>  : log display frequency\n");
      SCIPinfoMessage(scip, NULL, "   -t <timelimit> : enforce time limit\n");
      SCIPinfoMessage(scip, NULL, "   -m <memlimit>  : enforce memory limit\n");
      SCIPinfoMessage(scip, NULL, "   -c <command>   : execute command line\n");
   }

   return SCIP_OKAY;
}

/** creates a SCIP instance with default plugins, evaluates command line parameters,
 *  runs SCIP appropriately, and frees the SCIP instance
 */
static
SCIP_RETCODE runShell(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   clock_t               startclock,         /**< clock at which the process started */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
   SCIP_MESSAGEHDLR* messagehdlr = NULL;
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create own buffered message handler for PB console log */
   SCIP_CALL( SCIPcreateMessagehdlrPbSolver(&messagehdlr, TRUE, NULL, FALSE) );

   /* set PB competition message handler */
   SCIP_CALL( SCIPsetMessagehdlr(scip, messagehdlr) );

   /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL( processShellArguments(scip, argc, argv, startclock, defaultsetname) );

   /* release captured and own message handler */
   SCIP_CALL( SCIPmessagehdlrRelease(&messagehdlr) );

   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );
   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main method starting the PBSolver */
int main(
   int                   argc,               /**< number of arguments */
   char**                argv                /**< string array with arguments */
   )
{
   SCIP_RETCODE retcode;
   clock_t startclock, endclock;

   startclock = clock();

   retcode = runShell(argc, argv, startclock, "scip.set");

   endclock = clock();

   if( retcode != SCIP_OKAY )
      printf("s UNKNOWN\n");

   printf("c Time complete (sec): %10.3lf\n", (double)(endclock - startclock) / (double)CLOCKS_PER_SEC);

  return 0;
}
