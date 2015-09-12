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

/**@file   benders.c
 * @brief  run Benders algorithm
 * @author Marc Pfetsch
 */

#include "benders.h"

/* default parameters */
#define DEFAULT_SOLVE_MASTER_APPROX     TRUE      /**< solve master problem approximately */
#define DEFAULT_COMPUTE_HEUR_SOL       FALSE      /**< compute heuristic solution for master problem */
#define DEFAULT_MASTER_GAP               0.1      /**< gap bound for approximately solving the master problem */
#define DEFAULT_REOPTIMIZATION          TRUE      /**< Use reoptimization to solve master problem? */

#define DEFAULT_MAXSTALLNODES            50L      /**< maximal number of stalling nodes in solving the separation problem */
#define DEFAULT_MASTERSTALLNODES       5000L      /**< stall nodes for the master problem */
#define DEFAULT_HEURSTALLNODES         1000L      /**< stall nodes for the heuristic */

/* other parameters */
#define MAXITERATIONS                   1000      /**< maximal number of iterations of main loop */



/** find minimum cardinality infeasible subsystem */
SCIP_RETCODE runBenders(
   SCIP*                 masterscip,         /**< master SCIP instance */
   BENDERS_CUTORACLE((*Oracle)),             /**< oracle for generation of a Benders cut */
   BENDERS_DATA*         data,               /**< user data for oracle */
   SCIP_Real             timelimit,          /**< time limit read from arguments */
   SCIP_Real             memlimit,           /**< memory limit read from arguments */
   int                   dispfreq,           /**< display frequency */
   SCIP_Bool             usereopt,           /**< Use reoptimization? */
   SCIP_Bool             solvemasterapprox,  /**< Solve master problem approximately? */
   SCIP_Longint          masterstallnodes,   /**< stall nodes for master problem if solvemasterapprox is true */
   SCIP_Real             mastergap,          /**< gap limit for master problem if solvemasterapprox is true */
   SCIP_VERBLEVEL        verblevel,          /**< verbosity level for output */
   SCIP_STATUS*          status              /**< status of optimization */
   )
{
   SCIP_CLOCK* totaltimeclock;
   SCIP_CLOCK* oracletimeclock;
   SCIP_CLOCK* mastertimeclock;
   SCIP_Bool masteroptimal = TRUE;
   const int maxIters = MAXITERATIONS;
   SCIP_VAR** mastervars;
   SCIP_Real* mastersolution;
   int nmastervars;
   int iter = 0;

   assert( status != NULL );
   *status = SCIP_STATUS_UNKNOWN;

   SCIP_CALL( SCIPtransformProb(masterscip) );
   SCIP_CALL( SCIPgetOrigVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, &mastersolution, nmastervars) );

   /* set output to console */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(masterscip, "display/verblevel", 5) );
#else
   SCIP_CALL( SCIPsetIntParam(masterscip, "display/verblevel", 0) );
#endif

   SCIP_CALL( SCIPsetRealParam(masterscip, "limits/memory", memlimit) );

   if ( dispfreq >= 0 )
      SCIP_CALL( SCIPsetIntParam(masterscip, "display/freq", dispfreq) );
   else
      SCIP_CALL( SCIPsetIntParam(masterscip, "display/freq", 1000) );

   /* possibly use reoptimization */
   if ( usereopt )
   {
      assert( SCIPgetNIntVars(masterscip) == 0 );
      if ( verblevel >= SCIP_VERBLEVEL_NORMAL )
         SCIPinfoMessage(masterscip, NULL, "\nUsing reoptimization.\n");
      SCIP_CALL( SCIPenableReoptimization(masterscip, TRUE) );
   }

   /* we use the time limit of the original SCIP version */
   SCIP_CALL( SCIPcreateClock(masterscip, &totaltimeclock) );
   SCIP_CALL( SCIPstartClock(masterscip, totaltimeclock) );

   SCIP_CALL( SCIPcreateClock(masterscip, &oracletimeclock) );
   SCIP_CALL( SCIPcreateClock(masterscip, &mastertimeclock) );

   /* output */
   if ( verblevel >= SCIP_VERBLEVEL_NORMAL )
   {
      if ( solvemasterapprox )
      {
         if ( ! SCIPisInfinity(masterscip, timelimit) )
            SCIPinfoMessage(masterscip, NULL, "\nApproximately solving master problem with time limit: %f ...\n", timelimit);
         else
            SCIPinfoMessage(masterscip, NULL, "\nApproximately solving master problem ...\n");
      }
      else
      {
         if ( ! SCIPisInfinity(masterscip, timelimit) )
            SCIPinfoMessage(masterscip, NULL, "\nOptimally solving master problem with time limit: %f ...\n", timelimit);
         else
            SCIPinfoMessage(masterscip, NULL, "\nOptimally solving master problem ...\n");
      }
   }

   /* print banner */
   if ( verblevel >= SCIP_VERBLEVEL_NORMAL )
      SCIPinfoMessage(masterscip, NULL, "  time |  iter | nconss| nvars |master|oracle| ncuts | dualbnd\n");

   /* iterate */
   do
   {
      SCIP_STATUS masterstatus = SCIP_STATUS_UNKNOWN;
      BENDERS_STATUS substatus = BENDERS_STATUS_UNKNOWN;
      SCIP_Bool success = FALSE;
      SCIP_Real currenttime;
      SCIP_Real subtimelimit;

      SCIP_SOL* mastersol = NULL;
      SCIP_Real mastersolobj;
      int ncuts;
      int v;

      SCIPdebugMessage("\n\nIteration %d:\n", iter);

      /* set current time limit */
      currenttime = SCIPgetClockTime(masterscip, totaltimeclock);
      if ( ! SCIPisInfinity(masterscip, timelimit) )
      {
         subtimelimit = timelimit - currenttime;
         if ( subtimelimit <= 0.1 )
         {
            SCIPdebugMessage("Time limit exceeded.\n");
            *status = SCIP_STATUS_TIMELIMIT;
            goto TERMINATE;
         }
      }
      else
         subtimelimit = SCIPinfinity(masterscip);
      SCIP_CALL( SCIPsetRealParam(masterscip, "limits/time", subtimelimit) );

      /* SCIP_CALL( SCIPprintOrigProblem(masterscip, NULL, "cip", FALSE) ); */

      /* set gap limit if we do not have to solve the master to optimality */
      if ( solvemasterapprox )
      {
         SCIP_CALL( SCIPsetLongintParam(masterscip, "limits/stallnodes", masterstallnodes) );
         SCIP_CALL( SCIPsetRealParam(masterscip, "limits/gap", mastergap) );
      }

      /* solve master problem */
      SCIP_CALL( SCIPstartClock(masterscip, mastertimeclock) );
      SCIP_CALL( SCIPsolve(masterscip) );

      /* possibly reset gap limit */
      if ( solvemasterapprox )
      {
         SCIP_CALL( SCIPsetLongintParam(masterscip, "limits/stallnodes", -1LL) );
         SCIP_CALL( SCIPsetRealParam(masterscip, "limits/gap", 0.0) );
      }
      SCIP_CALL( SCIPstopClock(masterscip, mastertimeclock) );

      masterstatus = SCIPgetStatus(masterscip);

      /* determine master problem solution status */
      masteroptimal = FALSE;
      switch ( masterstatus )
      {
      case SCIP_STATUS_OPTIMAL:
         masteroptimal = TRUE;
         break;

      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_UNBOUNDED:
         /* do nothing */
         break;

      case SCIP_STATUS_INFEASIBLE:
         SCIPinfoMessage(masterscip, NULL, "Master problem infeasible.\n");
         *status = SCIP_STATUS_INFEASIBLE;
         goto TERMINATE;

      case SCIP_STATUS_TIMELIMIT:
         *status = SCIP_STATUS_TIMELIMIT;
         SCIPinfoMessage(masterscip, NULL, "Time limit exceeded.\n");
         goto TERMINATE;

      case SCIP_STATUS_USERINTERRUPT:
         *status = SCIP_STATUS_USERINTERRUPT;
         SCIPinfoMessage(masterscip, NULL, "User interrupt.\n");
         return SCIP_OKAY;

      default:
         SCIPerrorMessage("Master problem returned with status %d. Exiting ...\n", masterstatus);
         return SCIP_ERROR;
      }

      mastersol = SCIPgetBestSol(masterscip);
      if ( mastersol == NULL )
      {
         SCIPerrorMessage("Benders master problem does not have a primal solution!\n");
         return SCIP_ERROR;
      }
      mastersolobj = SCIPgetSolOrigObj(masterscip, mastersol);

      /* copy solution */
      for (v = 0; v < nmastervars; ++v)
      {
         SCIP_Real val;

         val = SCIPgetSolVal(masterscip, mastersol, mastervars[v]);
         assert( SCIPisIntegral(masterscip, val) );
         mastersolution[v] = val;
      }

      /* compute current time limit */
      currenttime = SCIPgetClockTime(masterscip, totaltimeclock);
      if ( ! SCIPisInfinity(masterscip, timelimit) )
      {
         subtimelimit = timelimit - currenttime;
         if ( subtimelimit <= 0.1 )
         {
            SCIPdebugMessage("Time limit exceeded.\n");
            goto TERMINATE;
         }
         SCIPdebugMessage("Solving separation problem ... (time limit: %g)\n", subtimelimit);
      }
      else
      {
         subtimelimit = SCIPinfinity(masterscip);
         SCIPdebugMessage("Solving separation problem ...\n");
      }

      /* free solving data */
      SCIP_CALL( SCIPfreeTransform(masterscip) );

      /* check for Benders cuts */
      SCIP_CALL( SCIPstartClock(masterscip, oracletimeclock) );
      SCIP_CALL( Oracle(masterscip, nmastervars, mastervars, mastersolution, data, timelimit, &ncuts, &substatus) );
      SCIP_CALL( SCIPstopClock(masterscip, oracletimeclock) );

      if ( verblevel >= SCIP_VERBLEVEL_NORMAL )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, " ");
         SCIPdispTime(SCIPgetMessagehdlr(masterscip), NULL, SCIPgetClockTime(masterscip, totaltimeclock), 6);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispInt(SCIPgetMessagehdlr(masterscip), NULL, iter, 7);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispInt(SCIPgetMessagehdlr(masterscip), NULL, SCIPgetNConss(masterscip), 7);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispInt(SCIPgetMessagehdlr(masterscip), NULL, SCIPgetNVars(masterscip), 7);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispTime(SCIPgetMessagehdlr(masterscip), NULL, SCIPgetClockTime(masterscip, mastertimeclock), 6);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispTime(SCIPgetMessagehdlr(masterscip), NULL, SCIPgetClockTime(masterscip, oracletimeclock), 6);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispInt(SCIPgetMessagehdlr(masterscip), NULL, ncuts, 7);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|%13.6e\n", mastersolobj);
      }

      switch ( substatus )
      {
      case BENDERS_STATUS_ADDEDCUT:
         break;

      case BENDERS_STATUS_SUCESS:
         success = TRUE;
         break;

      case BENDERS_STATUS_TIMELIMIT:
         *status = SCIP_STATUS_TIMELIMIT;
         goto TERMINATE;

      default:
         SCIPerrorMessage("Subproblem returned with status %d. Exiting ...\n", substatus);
         return SCIP_ERROR;
      }

      if ( success )
         break;

      ++iter;
   }
   while ( iter < maxIters );

   SCIPdebugMessage("Solution process finished.\n");

   if ( iter >= maxIters )
   {
      *status = SCIP_STATUS_TOTALNODELIMIT;
      SCIPinfoMessage(masterscip, NULL, "Reached iteration limit.\n");
   }

   if ( masteroptimal )
   {
      assert( *status == SCIP_STATUS_UNKNOWN );
      *status = SCIP_STATUS_OPTIMAL;
   }

 TERMINATE:
   SCIPinfoMessage(masterscip, NULL, "\nTotal used time:\t %f\n", SCIPgetClockTime(masterscip, totaltimeclock));
   SCIPinfoMessage(masterscip, NULL, "Oracle time:\t\t %f\n", SCIPgetClockTime(masterscip, oracletimeclock));
   SCIPinfoMessage(masterscip, NULL, "Master problem time:\t %f\n", SCIPgetClockTime(masterscip, mastertimeclock));
   SCIPinfoMessage(masterscip, NULL, "Number of iterations:\t %d\n", iter);

   SCIPfreeBlockMemoryArray(masterscip, &mastersolution, nmastervars);

   SCIP_CALL( SCIPfreeClock(masterscip, &mastertimeclock) );
   SCIP_CALL( SCIPfreeClock(masterscip, &oracletimeclock) );
   SCIP_CALL( SCIPfreeClock(masterscip, &totaltimeclock) );

   return SCIP_OKAY;
}
