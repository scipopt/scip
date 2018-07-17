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

/**@file   benders.c
 * @brief  run Benders algorithm
 * @author Marc Pfetsch
 *
 * Run Benders algorithm using an oracle for solving the subproblems and solving the master problem to optimality.
 */

#include "benders.h"

/* other parameters */
#define MAXITERATIONS    10000               /**< maximal number of iterations of main loop */


/** output status */
static
SCIP_RETCODE printStatus(
   SCIP*                 masterscip,         /**< master problem SCIP instance */
   SCIP_STATUS           status              /**< solution status */
   )
{
   SCIPinfoMessage(masterscip, NULL, "SCIP Status        : ");
   switch ( status )
   {
   case SCIP_STATUS_UNKNOWN:
      SCIPinfoMessage(masterscip, NULL, "unknown");
      break;
   case SCIP_STATUS_USERINTERRUPT:
      SCIPinfoMessage(masterscip, NULL, "solving was interrupted [user interrupt]");
      break;
   case SCIP_STATUS_NODELIMIT:
      SCIPinfoMessage(masterscip, NULL, "solving was interrupted [node limit reached]");
      break;
   case SCIP_STATUS_TIMELIMIT:
      SCIPinfoMessage(masterscip, NULL, "solving was interrupted [time limit reached]");
      break;
   case SCIP_STATUS_MEMLIMIT:
      SCIPinfoMessage(masterscip, NULL, "solving was interrupted [memory limit reached]");
      break;
   case SCIP_STATUS_GAPLIMIT:
      SCIPinfoMessage(masterscip, NULL, "solving was interrupted [gap limit reached]");
      break;
   case SCIP_STATUS_OPTIMAL:
      SCIPinfoMessage(masterscip, NULL, "problem is solved [optimal solution found]");
      break;
   case SCIP_STATUS_INFEASIBLE:
      SCIPinfoMessage(masterscip, NULL, "problem is solved [infeasible]");
      break;
   case SCIP_STATUS_UNBOUNDED:
      SCIPinfoMessage(masterscip, NULL, "problem is solved [unbounded]");
      break;
   case SCIP_STATUS_TOTALNODELIMIT:
      SCIPinfoMessage(masterscip, NULL, "solving was interrupted [iteration limit reached]");
      break;
   case SCIP_STATUS_STALLNODELIMIT:
   case SCIP_STATUS_SOLLIMIT:
   case SCIP_STATUS_BESTSOLLIMIT:
   case SCIP_STATUS_RESTARTLIMIT:
   case SCIP_STATUS_INFORUNBD:
      SCIPerrorMessage("unexpected status code <%d>\n", status);
      return SCIP_INVALIDDATA;
   default:
      SCIPerrorMessage("invalid status code <%d>\n", status);
      return SCIP_INVALIDDATA;
   }
   SCIPinfoMessage(masterscip, NULL, "\n");

   return SCIP_OKAY;
}


/** output short statistics */
static
SCIP_RETCODE printShortStatistics(
   SCIP*                 masterscip,         /**< master problem SCIP instance */
   SCIP_STATUS           status,             /**< solution status */
   SCIP_CLOCK*           totaltimeclock,     /**< clock for total time */
   SCIP_Real             primalbound,        /**< primal bound */
   SCIP_Real             dualbound,          /**< dual bound */
   SCIP_Longint          ntotalnodes,        /**< total number of nodes */
   int                   niter               /**< number of iterations */
   )
{
   SCIP_Real gap = 1e20;

   if ( ! SCIPisInfinity(masterscip, primalbound) && ! SCIPisInfinity(masterscip, -dualbound) )
      gap = fabs(primalbound - dualbound)/(MAX3(fabs(primalbound), fabs(dualbound), 1.0));

   /* start output */
   SCIPinfoMessage(masterscip, NULL, "\n");

   SCIP_CALL( printStatus(masterscip, status) );
   SCIPinfoMessage(masterscip, NULL, "Solving Time (sec) : %.2f\n", SCIPgetClockTime(masterscip, totaltimeclock));
   SCIPinfoMessage(masterscip, NULL, "Solving Nodes      : %" SCIP_LONGINT_FORMAT " (total of %" SCIP_LONGINT_FORMAT " nodes in %d runs)\n",
      ntotalnodes, ntotalnodes, niter);
   SCIPinfoMessage(masterscip, NULL, "Primal Bound       : %+21.14e\n", primalbound);
   SCIPinfoMessage(masterscip, NULL, "Dual Bound         : %+21.14e\n", dualbound);
   if ( SCIPisInfinity(masterscip, gap) )
      SCIPinfoMessage(masterscip, NULL, "Gap                : infinite\n");
   else
      SCIPinfoMessage(masterscip, NULL, "Gap                : %.2f %%\n", 100.0 * gap);
   SCIPinfoMessage(masterscip, NULL, "\n");

   return SCIP_OKAY;
}


/** output long statistics */
static
SCIP_RETCODE printLongStatistics(
   SCIP*                 masterscip,         /**< master problem SCIP instance */
   SCIP_STATUS           status,             /**< solution status */
   SCIP_CLOCK*           totaltimeclock,     /**< clock for total time */
   SCIP_CLOCK*           oracletimeclock,    /**< clock for oracle */
   SCIP_CLOCK*           mastertimeclock,    /**< clock for master problem */
   SCIP_Real             primalbound,        /**< primal bound */
   SCIP_Real             dualbound,          /**< dual bound */
   SCIP_Longint          ntotalnodes,        /**< total number of nodes */
   SCIP_Longint          ntotalcuts,         /**< total number of cuts */
   int                   niter               /**< number of iterations */
   )
{
   SCIP_Real gap = 1e20;

   if ( ! SCIPisInfinity(masterscip, primalbound) && ! SCIPisInfinity(masterscip, -dualbound) )
      gap = fabs(primalbound - dualbound)/(MAX3(fabs(primalbound), fabs(dualbound), 1.0));

   /* start output */
   SCIPinfoMessage(masterscip, NULL, "\n");

   /* print main part of statistics */
   SCIP_CALL( printStatus(masterscip, status) );

   SCIPinfoMessage(masterscip, NULL, "Total Time         : %10.2f\n", SCIPgetClockTime(masterscip, totaltimeclock));
   SCIPinfoMessage(masterscip, NULL, "  solving          : %10.2f\n", SCIPgetClockTime(masterscip, totaltimeclock));
   SCIPinfoMessage(masterscip, NULL, "  master           : %10.2f (included in solving)\n", SCIPgetClockTime(masterscip, mastertimeclock));
   SCIPinfoMessage(masterscip, NULL, "  oracle           : %10.2f (included in solving)\n", SCIPgetClockTime(masterscip, oracletimeclock));

   SCIPinfoMessage(masterscip, NULL, "Original Problem   :\n");
   SCIPinfoMessage(masterscip, NULL, "  Problem name     : %s\n", SCIPgetProbName(masterscip));
   SCIPinfoMessage(masterscip, NULL, "  Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      SCIPgetNVars(masterscip), SCIPgetNOrigVars(masterscip), 0, 0, 0);
   SCIPinfoMessage(masterscip, NULL, "  Constraints      : %d initial, %d maximal\n", 1, SCIPgetNOrigConss(masterscip));
   SCIPinfoMessage(masterscip, NULL, "  Objective sense  : minimize\n");

   SCIPinfoMessage(masterscip, NULL, "Presolved Problem  :\n");
   SCIPinfoMessage(masterscip, NULL, "  Problem name     : %s\n", SCIPgetProbName(masterscip));
   SCIPinfoMessage(masterscip, NULL, "  Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      SCIPgetNVars(masterscip), SCIPgetNBinVars(masterscip), SCIPgetNIntVars(masterscip), SCIPgetNImplVars(masterscip), SCIPgetNContVars(masterscip));
   SCIPinfoMessage(masterscip, NULL, "  Constraints      : %d initial, %d maximal\n", SCIPgetNConss(masterscip), SCIPgetNOrigConss(masterscip));

   SCIPinfoMessage(masterscip, NULL, "Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children\n");
   SCIPinfoMessage(masterscip, NULL, "  %-17.17s: %10d %10d %10d %10d %10d %10d %10d %10d %10d %10d %10d %10d %10" SCIP_LONGINT_FORMAT " %10d\n",
      "benders", 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ntotalcuts, 0);

   SCIPinfoMessage(masterscip, NULL, "Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp    SB-Prop\n");
   SCIPinfoMessage(masterscip, NULL, "  %-17.17s: %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n", "benders",
      SCIPgetClockTime(masterscip, oracletimeclock), 0.0, SCIPgetClockTime(masterscip, oracletimeclock), 0.0, 0.0, 0.0, 0.0, 0.0);

   SCIPinfoMessage(masterscip, NULL, "B&B Tree           :\n");
   SCIPinfoMessage(masterscip, NULL, "  number of runs   : %10d\n", niter);
   SCIPinfoMessage(masterscip, NULL, "  nodes (total)    : %10" SCIP_LONGINT_FORMAT "\n", ntotalnodes);

   SCIPinfoMessage(masterscip, NULL, "Solution           :\n");
   SCIPinfoMessage(masterscip, NULL, "  Primal Bound     : %+21.14e\n", primalbound);
   SCIPinfoMessage(masterscip, NULL, "  Dual Bound       : %+21.14e\n", dualbound);
   if ( SCIPisInfinity(masterscip, gap) )
      SCIPinfoMessage(masterscip, NULL, "  Gap              :   infinite\n");
   else
      SCIPinfoMessage(masterscip, NULL, "  Gap              : %10.2f %%\n", 100.0 * gap);

#ifdef SCIP_OUTPUT
   SCIPinfoMessage(masterscip, NULL, "\nTotal used time:\t %f\n", SCIPgetClockTime(masterscip, totaltimeclock));
   SCIPinfoMessage(masterscip, NULL, "Oracle time:\t\t %f\n", SCIPgetClockTime(masterscip, oracletimeclock));
   SCIPinfoMessage(masterscip, NULL, "Master problem time:\t %f\n", SCIPgetClockTime(masterscip, mastertimeclock));
   SCIPinfoMessage(masterscip, NULL, "Number of iterations:\t %d\n", niter);
#endif

   return SCIP_OKAY;
}


/** run Benders algorithm using an oracle for the subproblems */
SCIP_RETCODE runBenders(
   SCIP*                 masterscip,         /**< master SCIP instance */
   BENDERS_CUTORACLE((*Oracle)),             /**< oracle for Benders subproblem */
   BENDERS_DATA*         data,               /**< user data for oracle */
   SCIP_Real             timelimit,          /**< time limit read from arguments */
   SCIP_Real             memlimit,           /**< memory limit read from arguments */
   int                   dispfreq,           /**< display frequency */
   SCIP_Bool             usereopt,           /**< Use reoptimization? */
   SCIP_Bool             solvemasterapprox,  /**< Solve master problem approximately? */
   SCIP_Longint          masterstallnodes,   /**< stall nodes for master problem if solvemasterapprox is true */
   SCIP_Real             mastergaplimit,     /**< gap limit for master problem if solvemasterapprox is true */
   SCIP_VERBLEVEL        verblevel,          /**< verbosity level for output */
   SCIP_STATUS*          status              /**< status of optimization */
   )
{  /*lint --e{788}*/
   SCIP_CLOCK* totaltimeclock;
   SCIP_CLOCK* oracletimeclock;
   SCIP_CLOCK* mastertimeclock;
   SCIP_Bool masteroptimal = TRUE;
   const int maxIters = MAXITERATIONS;
   SCIP_Longint ntotalnodes = 0LL;
   SCIP_Longint ntotalcuts = 0LL;
   SCIP_VAR** mastervars;
   SCIP_Real* mastersolution;
   SCIP_Real primalbound = 1e20;
   SCIP_Real dualbound = -1e20;
   SCIP_Real mastersolobj = 0.0;
   int nmastervars;
   int niter = 0;

   assert( status != NULL );
   *status = SCIP_STATUS_UNKNOWN;

   SCIP_CALL( SCIPgetOrigVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPallocClearBlockMemoryArray(masterscip, &mastersolution, nmastervars) );

   /* set output to console */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(masterscip, "display/verblevel", 5) );
#else
   SCIP_CALL( SCIPsetIntParam(masterscip, "display/verblevel", 0) );
#endif

   if ( ! SCIPisInfinity(masterscip, memlimit) )
   {
      SCIP_CALL( SCIPsetRealParam(masterscip, "limits/memory", memlimit) );
   }

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

   /* set up clocks */
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
            SCIPinfoMessage(masterscip, NULL, "\nApproximately solving master problem with time limit %.1f and gap limit %.2f%% ...\n", timelimit, 100.0 * mastergaplimit);
         else
            SCIPinfoMessage(masterscip, NULL, "\nApproximately solving master problem with gap limit %.2f%% ...\n", 100.0 * mastergaplimit);
      }
      else
      {
         if ( ! SCIPisInfinity(masterscip, timelimit) )
            SCIPinfoMessage(masterscip, NULL, "\nOptimally solving master problem with time limit: %.1f ...\n", timelimit);
         else
            SCIPinfoMessage(masterscip, NULL, "\nOptimally solving master problem ...\n");
      }
   }

   /* print banner */
   if ( verblevel >= SCIP_VERBLEVEL_NORMAL )
   {
      if ( solvemasterapprox )
         SCIPinfoMessage(masterscip, NULL, "  time | niter |nconss | nvars |master|totalnodes|oracle| ncuts |    dualbound |   gap\n");
      else
         SCIPinfoMessage(masterscip, NULL, "  time | niter |nconss | nvars |master|totalnodes|oracle| ncuts |    dualbound\n");
   }

   /* iterate */
   do
   {
      BENDERS_STATUS substatus = BENDERS_STATUS_UNKNOWN;
      SCIP_STATUS masterstatus = SCIP_STATUS_UNKNOWN;
      SCIP_Bool success = FALSE;
      SCIP_Real currenttime;
      SCIP_Real subtimelimit;
      SCIP_SOL* mastersol = NULL;
      SCIP_Real mastergap = 1e20;
      int ncuts = 0;
      int v;

      ++niter;
      SCIPdebugMessage("Iteration %d.\n", niter);

      /* --------- solve Benders subproblem */

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

      /* free solving data in order to add constraints */
      SCIP_CALL( SCIPfreeTransform(masterscip) );

      /* check for Benders cuts */
      SCIP_CALL( SCIPstartClock(masterscip, oracletimeclock) );
      SCIP_CALL( Oracle(masterscip, nmastervars, mastervars, mastersolution, data, timelimit, ntotalcuts, &ncuts, &substatus) );
      SCIP_CALL( SCIPstopClock(masterscip, oracletimeclock) );
      ntotalcuts += (SCIP_Longint) ncuts;

      switch ( substatus )
      {
      case BENDERS_STATUS_ADDEDCUT:
         break;

      case BENDERS_STATUS_SUCESS:
         success = TRUE;
         primalbound = mastersolobj;
         break;

      case BENDERS_STATUS_TIMELIMIT:
         *status = SCIP_STATUS_TIMELIMIT;
         goto TERMINATE;

      case BENDERS_STATUS_USERINTERRUPT:
         *status = SCIP_STATUS_USERINTERRUPT;
         goto TERMINATE;

      default:
         SCIPerrorMessage("Subproblem returned with status %d. Exiting ...\n", substatus);
         return SCIP_ERROR;
      }

      /* if not cuts could be found, the master problem is solved optimally */
      if ( success )
      {
         /* if last master problem was solved to optimality, we are done */
         if ( masteroptimal )
            break;
         /* otherwise, we have to resolve the master to optimality */
         solvemasterapprox = FALSE;

         if ( verblevel >= SCIP_VERBLEVEL_NORMAL )
            SCIPinfoMessage(masterscip, NULL, "Switching to optimally solving the master problem.\n");
      }

      /* --------- solve Benders master problem */

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
         SCIP_CALL( SCIPsetRealParam(masterscip, "limits/gap", mastergaplimit) );
      }

      /* solve master problem */
      SCIP_CALL( SCIPstartClock(masterscip, mastertimeclock) );
      SCIP_CALL( SCIPsolve(masterscip) );

      ntotalnodes += SCIPgetNTotalNodes(masterscip);

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
         if ( verblevel >= SCIP_VERBLEVEL_NORMAL )
            SCIPinfoMessage(masterscip, NULL, "Master problem infeasible.\n");
         *status = SCIP_STATUS_INFEASIBLE;
         goto TERMINATE;

      case SCIP_STATUS_TIMELIMIT:
         *status = SCIP_STATUS_TIMELIMIT;
         if ( verblevel >= SCIP_VERBLEVEL_NORMAL )
            SCIPinfoMessage(masterscip, NULL, "Time limit exceeded.\n");
         goto TERMINATE;

      case SCIP_STATUS_USERINTERRUPT:
         *status = SCIP_STATUS_USERINTERRUPT;
         if ( verblevel >= SCIP_VERBLEVEL_NORMAL )
            SCIPinfoMessage(masterscip, NULL, "User interrupt.\n");
         goto TERMINATE;

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
      mastergap = SCIPgetGap(masterscip);

      dualbound = MAX(dualbound, SCIPgetDualbound(masterscip));

      /* copy solution */
      for (v = 0; v < nmastervars; ++v)
      {
         SCIP_Real val;

         val = SCIPgetSolVal(masterscip, mastersol, mastervars[v]);
         assert( SCIPisIntegral(masterscip, val) );
         mastersolution[v] = val;
      }

      if ( verblevel >= SCIP_VERBLEVEL_NORMAL )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, " ");
         SCIPdispTime(SCIPgetMessagehdlr(masterscip), NULL, SCIPgetClockTime(masterscip, totaltimeclock), 6);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispInt(SCIPgetMessagehdlr(masterscip), NULL, niter, 7);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispInt(SCIPgetMessagehdlr(masterscip), NULL, SCIPgetNOrigConss(masterscip), 7);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispInt(SCIPgetMessagehdlr(masterscip), NULL, SCIPgetNOrigVars(masterscip), 7);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispTime(SCIPgetMessagehdlr(masterscip), NULL, SCIPgetClockTime(masterscip, mastertimeclock), 6);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispLongint(SCIPgetMessagehdlr(masterscip), NULL, ntotalnodes, 10);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispTime(SCIPgetMessagehdlr(masterscip), NULL, SCIPgetClockTime(masterscip, oracletimeclock), 6);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|");
         SCIPdispInt(SCIPgetMessagehdlr(masterscip), NULL, ncuts, 7);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "|%13.6e", mastersolobj);
         if ( solvemasterapprox )
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, " | %6.2f%%", 100.0 * mastergap);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(masterscip), NULL, "\n");
      }
   }
   while ( niter < maxIters );

   SCIPdebugMessage("Solution process finished.\n");

   if ( niter >= maxIters )
   {
      *status = SCIP_STATUS_TOTALNODELIMIT;
      if ( verblevel >= SCIP_VERBLEVEL_NORMAL )
         SCIPinfoMessage(masterscip, NULL, "Reached iteration limit.\n");
   }

   if ( masteroptimal )
   {
      assert( *status == SCIP_STATUS_UNKNOWN );
      *status = SCIP_STATUS_OPTIMAL;
   }

 TERMINATE:

   if ( verblevel >= SCIP_VERBLEVEL_NORMAL )
   {
      SCIP_CALL( printShortStatistics(masterscip, *status, totaltimeclock, primalbound, dualbound, ntotalnodes, niter) );
      SCIP_CALL( printLongStatistics(masterscip, *status, totaltimeclock, oracletimeclock, mastertimeclock, primalbound, dualbound, ntotalnodes, ntotalcuts, niter) );
   }

   SCIPfreeBlockMemoryArray(masterscip, &mastersolution, nmastervars);

   SCIP_CALL( SCIPfreeClock(masterscip, &mastertimeclock) );
   SCIP_CALL( SCIPfreeClock(masterscip, &oracletimeclock) );
   SCIP_CALL( SCIPfreeClock(masterscip, &totaltimeclock) );

   return SCIP_OKAY;
}
