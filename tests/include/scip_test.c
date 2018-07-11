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

#include "scip/nodesel_dfs.h"
#include "scip/cons_integral.h"

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecTest)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPinterruptSolve(scip) );
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTest)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPinterruptSolve(scip) );
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}

/* it can be called in SCIP_STAGE_PROBLEM and can get to
 *  SCIP_STAGE_TRANSFORMED
 *  SCIP_STAGE_PRESOLVING
 *  SCIP_STAGE_PRESOLVED
 *  SCIP_STAGE_SOLVING
 *  SCIP_STAGE_SOLVED
 *
 *  If stage == SCIP_STAGE_SOLVING and enableNLP is true, then SCIP will build its NLP
 */
SCIP_RETCODE TESTscipSetStage(SCIP* scip, SCIP_STAGE stage, SCIP_Bool enableNLP)
{
   /* no output nor warnings */
   /* @todo reset output to previous level! */
   SCIP_CALL( SCIPsetParam(scip, "display/verblevel", 0) );

   /* needs one node selector to call SCIPtransformProb, which is also
    * called by SCIP(pre)solve */
   SCIP_CALL( SCIPincludeNodeselDfs(scip) );
   /* we need to include integral conshdlr to supress a warning */
   SCIP_CALL( SCIPincludeConshdlrIntegral(scip) );

   switch( stage )
   {
      case SCIP_STAGE_TRANSFORMED:
         SCIP_CALL( SCIPtransformProb(scip) );
         break;

      case SCIP_STAGE_PRESOLVED:
         SCIP_CALL( SCIPpresolve(scip) );
         break;

      case SCIP_STAGE_SOLVED:
         SCIP_CALL( SCIPsolve(scip) );
         break;

      case SCIP_STAGE_PRESOLVING:
         SCIP_CALL( SCIPincludePresolBasic(scip, NULL, "presolTest",
                  "Presol to stop in PRESOLVING", 1, -1,
                  SCIP_PRESOLTIMING_ALWAYS, presolExecTest, NULL) );
         SCIP_CALL( SCIPpresolve(scip) );
         break;

      case SCIP_STAGE_SOLVING:
         SCIP_CALL( SCIPincludeHeurBasic(scip, NULL,
                  "heurTest", "heuristic to stop in SOLVING", '!', 1, 1, 0,
                  -1, SCIP_HEURTIMING_BEFORENODE, FALSE, heurExecTest, NULL) );

         /* enable NLP */
         if( enableNLP )
         {
            SCIP_CALL( SCIPpresolve(scip) );
            SCIPenableNLP(scip);
         }

         SCIP_CALL( SCIPsolve(scip) );
         break;

      default:
         fprintf(stderr, "Not implemented yet!\n");
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}
