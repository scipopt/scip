/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//#define SCIP_DEBUG
/**@file   branch_fullstrongprop.c
 * @brief  full strong LP branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_fullstrongprop.h"


#define BRANCHRULE_NAME          "fullstrongprop"
#define BRANCHRULE_DESC          "full strong branching with propagation"
#define BRANCHRULE_PRIORITY      0
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

#define DEFAULT_REEVALAGE        10LL        /**< number of intermediate LPs solved to trigger reevaluation of strong branching
                                              *   value for a variable that was already evaluated at the current node */


/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Longint          reevalage;          /**< number of intermediate LPs solved to trigger reevaluation of strong branching
                                              *   value for a variable that was already evaluated at the current node */
   int                   lastcand;           /**< last evaluated candidate of last branching rule execution */
};


/** gets strong branching information with propagation on column variable with fractional value */
static
SCIP_RETCODE getVarStrongbranchPropFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get strong branching values for */
   SCIP_Real             solval,             /**< ??????????????????? */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< pointer to store whether the downwards branch is infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< pointer to store whether the upwards branch is infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict,         /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   )
{
   SCIP_Bool cutoff;

   assert(lperror != NULL);

   if( downvalid != NULL )
      *downvalid = FALSE;
   if( upvalid != NULL )
      *upvalid = FALSE;
   if( downinf != NULL )
      *downinf = FALSE;
   if( upinf != NULL )
      *upinf = FALSE;
   if( downconflict != NULL )
      *downconflict = FALSE;
   if( upconflict != NULL )
      *upconflict = FALSE;

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      SCIPerrorMessage("cannot get strong branching information on non-COLUMN variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   assert(SCIPinProbing(scip));

   /* check if the solving process should be aborted */
   if( SCIPisStopped(scip) )
   {
      /* mark this as if the LP failed */
      *lperror = TRUE;
      return SCIP_OKAY;
   }

   assert(SCIPisFeasGE(scip, solval, SCIPvarGetLbLocal(var)));
   assert(SCIPisFeasLE(scip, solval, SCIPvarGetUbLocal(var)));

   SCIP_CALL( SCIPnewProbingNode(scip) );

   SCIP_CALL( SCIPchgVarUbProbing(scip, var, SCIPfloor(scip, solval)) );

   SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );

   if( cutoff )
   {
      if( downinf != NULL )
         *downinf = TRUE;

      *down = SCIPinfinity(scip);

      if( downvalid != NULL )
         *downvalid = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPsolveProbingLP(scip, -1, lperror) );

      SCIPdebugMessage("probing LP solution status: %d\n", SCIPgetLPSolstat(scip));

      if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         *down = SCIPgetLPObjval(scip);

         if( downvalid != NULL )
            *downvalid = TRUE;
      }
      else if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OBJLIMIT || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE )
      {
         if( downinf != NULL )
            *downinf = TRUE;

         *down = SCIPinfinity(scip);

         if( downvalid != NULL )
            *downvalid = TRUE;
      }
   }
   SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

   SCIP_CALL( SCIPnewProbingNode(scip) );

   if( SCIPisGT(scip, SCIPceil(scip, solval), SCIPvarGetLbLocal(var)) )
   {
      SCIP_CALL( SCIPchgVarLbProbing(scip, var, SCIPceil(scip, solval)) );
   }

   SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );

   if( cutoff )
   {
      if( upinf != NULL )
         *upinf = TRUE;

      *up = SCIPinfinity(scip);

      if( upvalid != NULL )
         *upvalid = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPsolveProbingLP(scip, -1, lperror) );

      SCIPdebugMessage("probing LP solution status: %d\n", SCIPgetLPSolstat(scip));

      if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         *up = SCIPgetLPObjval(scip);

         if( upvalid != NULL )
            *upvalid = TRUE;
      }
      else if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OBJLIMIT || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE )
      {
         if( upinf != NULL )
            *upinf = TRUE;

         *up = SCIPinfinity(scip);

         if( upvalid != NULL )
            *upvalid = TRUE;
      }
   }
   SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyFullstrongprop)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleFullstrongprop(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeFullstrongprop)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitFullstrongprop)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* initialize branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   branchruledata->lastcand = 0;

   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpFullstrongprop)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
#ifndef NDEBUG
   SCIP_Real cutoffbound;
#endif
   SCIP_Real lpobjval;
   SCIP_Real bestdown;
   SCIP_Real bestup;
   SCIP_Real bestscore;
   SCIP_Real provedbound;
   SCIP_Bool bestdownvalid;
   SCIP_Bool bestupvalid;
   SCIP_Bool allcolsinlp;
   SCIP_Bool exactsolve;
   int nlpcands;
   int npriolpcands;
   int bestcand;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of fullstrongprop branching\n");

   *result = SCIP_DIDNOTRUN;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get current LP objective bound of the local sub problem and global cutoff bound */
   lpobjval = SCIPgetLPObjval(scip);
#ifndef NDEBUG
   cutoffbound = SCIPgetCutoffbound(scip);
#endif

   /* check, if we want to solve the problem exactly, meaning that strong branching information is not useful
    * for cutting off sub problems and improving lower bounds of children
    */
   exactsolve = SCIPisExactSolve(scip);

   /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
   allcolsinlp = SCIPallColsInLP(scip);

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, &npriolpcands) );
   assert(nlpcands > 0);
   assert(npriolpcands > 0);

   /* if only one candidate exists, choose this one without applying strong branching */
   bestcand = 0;
   bestdown = lpobjval;
   bestup = lpobjval;
   bestdownvalid = TRUE;
   bestupvalid = TRUE;
   bestscore = -SCIPinfinity(scip);
   provedbound = lpobjval;
   if( nlpcands > 1 )
   {
      SCIP_Longint nodenum;
      SCIP_Real down;
      SCIP_Real up;
      SCIP_Real downgain;
      SCIP_Real upgain;
      SCIP_Real score;
      SCIP_Bool lperror;
      SCIP_Bool downvalid;
      SCIP_Bool upvalid;
      SCIP_Bool downinf;
      SCIP_Bool upinf;
      SCIP_Bool downconflict;
      SCIP_Bool upconflict;
      int nsbcalls;
      int i;
      int c;

      /* initialize strong branching */
      SCIP_CALL( SCIPstartProbing(scip) );

      /* get current node number */
      nodenum = SCIPgetNNodes(scip);

      /* search the full strong candidate
       * cycle through the candidates, starting with the position evaluated in the last run
       */
      nsbcalls = 0;
      for( i = 0, c = branchruledata->lastcand; i < nlpcands; ++i, ++c )
      {
         c = c % nlpcands;
         assert(lpcands[c] != NULL);

         {
            SCIPdebugMessage("applying strong branching on variable <%s> with solution %g\n",
               SCIPvarGetName(lpcands[c]), lpcandssol[c]);

            assert(SCIPisFeasGE(scip, lpcandssol[c], SCIPvarGetLbLocal(lpcands[c])));
            assert(SCIPisFeasLE(scip, lpcandssol[c], SCIPvarGetUbLocal(lpcands[c])));

            /* apply strong branching */
            SCIP_CALL( getVarStrongbranchPropFrac(scip, lpcands[c], lpcandssol[c], INT_MAX,
                  &down, &up, &downvalid, &upvalid, &downinf, &upinf, &downconflict, &upconflict, &lperror) );
            nsbcalls++;

            /* display node information line */
            if( SCIPgetDepth(scip) == 0 && nsbcalls % 100 == 0 )
            {
               SCIP_CALL( SCIPprintDisplayLine(scip, NULL, SCIP_VERBLEVEL_HIGH) );
            }

            /* check for an error in strong branching */
            if( lperror )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                  "(node %"SCIP_LONGINT_FORMAT") error in strong branching call for variable <%s> with solution %g\n",
                  SCIPgetNNodes(scip), SCIPvarGetName(lpcands[c]), lpcandssol[c]);
               break;
            }

            /* evaluate strong branching */
            down = MAX(down, lpobjval);
            up = MAX(up, lpobjval);
            downgain = down - lpobjval;
            upgain = up - lpobjval;
            assert(!allcolsinlp || exactsolve || !downvalid || downinf == SCIPisGE(scip, down, cutoffbound));
            assert(!allcolsinlp || exactsolve || !upvalid || upinf == SCIPisGE(scip, up, cutoffbound));
            assert(downinf || !downconflict);
            assert(upinf || !upconflict);

            /* check if there are infeasible roundings */
            if( downinf || upinf )
            {
               //assert(allcolsinlp);
               assert(!exactsolve);

               /* if for both infeasibilities, a conflict constraint was created, we don't need to fix the variable by
                * hand, but better wait for the next propagation round to fix them as an inference, and potentially
                * produce a cutoff that can be analyzed
                */
               if( allowaddcons && downinf == downconflict && upinf == upconflict )
               {
                  *result = SCIP_CONSADDED;
                  break; /* terminate initialization loop, because constraint was added */
               }
               else if( downinf && upinf )
               {
                  /* both roundings are infeasible -> node is infeasible */
                  *result = SCIP_CUTOFF;
                  SCIPdebugMessage(" -> variable <%s> is infeasible in both directions\n", SCIPvarGetName(lpcands[c]));
                  break; /* terminate initialization loop, because node is infeasible */
               }
               else if( downinf )
               {
                  SCIP_CALL( SCIPendProbing(scip) );

                  /* downwards rounding is infeasible -> change lower bound of variable to upward rounding */
                  if( SCIPisGT(scip, SCIPfeasCeil(scip, lpcandssol[c]), SCIPvarGetLbLocal(lpcands[c])) )
                  {
                     SCIP_CALL( SCIPchgVarLb(scip, lpcands[c], SCIPfeasCeil(scip, lpcandssol[c])) );
                  }
                  *result = SCIP_REDUCEDDOM;
                  SCIPdebugMessage(" -> variable <%s> is infeasible in downward branch\n", SCIPvarGetName(lpcands[c]));
                  break; /* terminate initialization loop, because LP was changed */
               }
               else
               {
                  assert(upinf);

                  SCIP_CALL( SCIPendProbing(scip) );

                  /* upwards rounding is infeasible -> change upper bound of variable to downward rounding */
                  if( SCIPisLT(scip, SCIPfeasFloor(scip, lpcandssol[c]), SCIPvarGetUbLocal(lpcands[c])) )
                  {
                     SCIP_CALL( SCIPchgVarUb(scip, lpcands[c], SCIPfeasFloor(scip, lpcandssol[c])) );
                  }
                  *result = SCIP_REDUCEDDOM;
                  SCIPdebugMessage(" -> variable <%s> is infeasible in upward branch\n", SCIPvarGetName(lpcands[c]));
                  break; /* terminate initialization loop, because LP was changed */
               }
            }
            else if( allcolsinlp && !exactsolve && downvalid && upvalid )
            {
               SCIP_Real minbound;

               /* the minimal lower bound of both children is a proved lower bound of the current subtree */
               minbound = MIN(down, up);
               provedbound = MAX(provedbound, minbound);
            }

            /* update pseudo cost values */
            assert(!downinf); /* otherwise, we would have terminated the initialization loop */
            assert(!upinf);
            SCIP_CALL( SCIPupdateVarPseudocost(scip, lpcands[c], 0.0-lpcandsfrac[c], downgain, 1.0) );
            SCIP_CALL( SCIPupdateVarPseudocost(scip, lpcands[c], 1.0-lpcandsfrac[c], upgain, 1.0) );
         }

         /* check for a better score, if we are within the maximum priority candidates */
         if( c < npriolpcands )
         {
            score = SCIPgetBranchScore(scip, lpcands[c], downgain, upgain);
            if( score > bestscore )
            {
               bestcand = c;
               bestdown = down;
               bestup = up;
               bestdownvalid = downvalid;
               bestupvalid = upvalid;
               bestscore = score;
            }
         }
         else
            score = 0.0;

         SCIPdebugMessage(" -> cand %d/%d (prio:%d) var <%s> (solval=%g, downgain=%g, upgain=%g, score=%g) -- best: <%s> (%g)\n",
            c, nlpcands, npriolpcands, SCIPvarGetName(lpcands[c]), lpcandssol[c], downgain, upgain, score,
            SCIPvarGetName(lpcands[bestcand]), bestscore);
      }

      /* end strong branching */
      if( SCIPinProbing(scip) )
      {
         SCIP_CALL( SCIPendProbing(scip) );
      }

      /* remember last evaluated candidate */
      branchruledata->lastcand = c;
   }

   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM && *result != SCIP_CONSADDED )
   {
      SCIP_NODE* downchild;
      SCIP_NODE* upchild;
      SCIP_VAR* var;

      assert(*result == SCIP_DIDNOTRUN);
      assert(0 <= bestcand && bestcand < nlpcands);
      assert(SCIPisLT(scip, provedbound, cutoffbound));

      var = lpcands[bestcand];

      /* perform the branching */
      SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%g, down=%g, up=%g, score=%g)\n",
         nlpcands, bestcand, SCIPvarGetName(var), lpcandssol[bestcand], bestdown, bestup, bestscore);
      SCIP_CALL( SCIPbranchVar(scip, var, &downchild, NULL, &upchild) );
      assert(downchild != NULL);
      assert(upchild != NULL);

      /* update the lower bounds in the children */
      if( allcolsinlp && !exactsolve )
      {
         SCIP_CALL( SCIPupdateNodeLowerbound(scip, downchild, bestdownvalid ? MAX(bestdown, provedbound) : provedbound) );
         SCIP_CALL( SCIPupdateNodeLowerbound(scip, upchild, bestupvalid ? MAX(bestup, provedbound) : provedbound) );
      }
      SCIPdebugMessage(" -> down child's lowerbound: %g\n", SCIPnodeGetLowerbound(downchild));
      SCIPdebugMessage(" -> up child's lowerbound: %g\n", SCIPnodeGetLowerbound(upchild));

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}


/*
 * branching specific interface methods
 */

/** creates the full strong LP branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleFullstrongprop(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create fullstrongprop branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
   branchruledata->lastcand = 0;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non-fundamental callbacks via specific setter functions*/
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyFullstrongprop) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeFullstrongprop) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitFullstrongprop) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpFullstrongprop) );

   /* fullstrongprop branching rule parameters */
   SCIP_CALL( SCIPaddLongintParam(scip,
         "branching/fullstrongprop/reevalage",
         "number of intermediate LPs solved to trigger reevaluation of strong branching value for a variable that was already evaluated at the current node",
         &branchruledata->reevalage, TRUE, DEFAULT_REEVALAGE, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
