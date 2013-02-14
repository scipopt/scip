/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_allfullstrong.c
 * @brief  all variables full strong LP branching rule
 * @author Tobias Achterberg
 *
 * The all variables full strong branching rule applies strong branching to every non-fixed variable
 * at the current node of the branch-and-bound search. The rule selects the candidate
 * which will cause the highest gain of the dual bound in the created sub-tree among all branching variables.
 *
 * For calculating the gain, a look-ahead is performed by solving the child node LPs which will result
 * from branching on a variable.
 *
 * For a more mathematical description and a comparison between the strong branching rule and other branching rules
 * in SCIP, we refer to
 *
 * @par
 * Tobias Achterberg@n
 * Constraint Integer Programming@n
 * PhD Thesis, Technische Universität Berlin, 2007@n
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_allfullstrong.h"


#define BRANCHRULE_NAME          "allfullstrong"
#define BRANCHRULE_DESC          "all variables full strong branching"
#define BRANCHRULE_PRIORITY      -1000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0


/** branching rule data */
struct SCIP_BranchruleData
{
   int                   lastcand;           /**< last evaluated candidate of last branching rule execution */
};


/** performs the all fullstrong branching */
static
SCIP_RETCODE branch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** pseudocands;
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
   int npseudocands;
   int npriopseudocands;
   int bestpseudocand;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

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

   /* get all non-fixed variables (not only the fractional ones) */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &pseudocands, &npseudocands, &npriopseudocands) );
   assert(npseudocands > 0);
   assert(npriopseudocands > 0);

   /* if only one candidate exists, choose this one without applying strong branching */
   bestpseudocand = 0;
   bestdown = lpobjval;
   bestup = lpobjval;
   bestdownvalid = TRUE;
   bestupvalid = TRUE;
   bestscore = -SCIPinfinity(scip);
   provedbound = lpobjval;
   if( npseudocands > 1 )
   {
      SCIP_Real solval;
      SCIP_Real down;
      SCIP_Real up;
      SCIP_Real downgain;
      SCIP_Real upgain;
      SCIP_Real score;
      SCIP_Bool integral;
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
      SCIP_CALL( SCIPstartStrongbranch(scip) );

      /* search the full strong candidate:
       * cycle through the candidates, starting with the position evaluated in the last run
       */
      nsbcalls = 0;
      for( i = 0, c = branchruledata->lastcand; i < npseudocands; ++i, ++c )
      {
         c = c % npseudocands;
         assert(pseudocands[c] != NULL);

         /* we can only apply strong branching on COLUMN variables that are in the current LP */
         if( !SCIPvarIsInLP(pseudocands[c]) )
            continue;

         solval = SCIPvarGetLPSol(pseudocands[c]);
         integral = SCIPisFeasIntegral(scip, solval);

         SCIPdebugMessage("applying strong branching on %s variable <%s>[%g,%g] with solution %g\n",
            integral ? "integral" : "fractional", SCIPvarGetName(pseudocands[c]), SCIPvarGetLbLocal(pseudocands[c]), 
            SCIPvarGetUbLocal(pseudocands[c]), solval);

         if( integral )
         {
            SCIP_CALL( SCIPgetVarStrongbranchInt(scip, pseudocands[c], INT_MAX, 
                  &down, &up, &downvalid, &upvalid, &downinf, &upinf, &downconflict, &upconflict, &lperror) );
         }
         else
         {
            SCIP_CALL( SCIPgetVarStrongbranchFrac(scip, pseudocands[c], INT_MAX, 
                  &down, &up, &downvalid, &upvalid, &downinf, &upinf, &downconflict, &upconflict, &lperror) );
         }
         nsbcalls++;

         /* display node information line in root node */
         if( SCIPgetDepth(scip) == 0 && nsbcalls % 100 == 0 )
         {
            SCIP_CALL( SCIPprintDisplayLine(scip, NULL, SCIP_VERBLEVEL_HIGH) );
         }

         /* check for an error in strong branching */
         if( lperror )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "(node %"SCIP_LONGINT_FORMAT") error in strong branching call for variable <%s> with solution %g\n", 
               SCIPgetNNodes(scip), SCIPvarGetName(pseudocands[c]), solval);
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
            assert(allcolsinlp);
            assert(!exactsolve);
            
            /* if for both infeasibilities, a conflict constraint was created, we don't need to fix the variable by hand,
             * but better wait for the next propagation round to fix them as an inference, and potentially produce a
             * cutoff that can be analyzed
             */
            if( allowaddcons && downinf == downconflict && upinf == upconflict )
            {
               *result = SCIP_CONSADDED;
               break; /* terminate initialization loop, because constraint was added */
            }
            else if( downinf && upinf )
            {
               if( integral )
               {
                  SCIP_Bool infeasible;
                  SCIP_Bool fixed;
                  
                  /* both bound changes are infeasible: variable can be fixed to its current value */
                  SCIP_CALL( SCIPfixVar(scip, pseudocands[c], solval, &infeasible, &fixed) );
                  assert(!infeasible);
                  assert(fixed);
                  *result = SCIP_REDUCEDDOM;
                  SCIPdebugMessage(" -> integral variable <%s> is infeasible in both directions\n",
                     SCIPvarGetName(pseudocands[c]));
                  break; /* terminate initialization loop, because LP was changed */
               }
               else
               {
                  /* both roundings are infeasible: the node is infeasible */
                  *result = SCIP_CUTOFF;
                  SCIPdebugMessage(" -> fractional variable <%s> is infeasible in both directions\n",
                     SCIPvarGetName(pseudocands[c]));
                  break; /* terminate initialization loop, because node is infeasible */
               }
            }
            else if( downinf )
            {
               SCIP_Real newlb;

               /* downwards rounding is infeasible -> change lower bound of variable to upward rounding */
               newlb = SCIPfeasCeil(scip, solval);
               if( SCIPvarGetLbLocal(pseudocands[c]) < newlb - 0.5 )
               {
                  SCIP_CALL( SCIPchgVarLb(scip, pseudocands[c], newlb) );
                  *result = SCIP_REDUCEDDOM;
                  SCIPdebugMessage(" -> variable <%s> is infeasible in downward branch\n", SCIPvarGetName(pseudocands[c]));
                  break; /* terminate initialization loop, because LP was changed */
               }
            }
            else
            {
               SCIP_Real newub;

               /* upwards rounding is infeasible -> change upper bound of variable to downward rounding */
               assert(upinf);
               newub = SCIPfeasFloor(scip, solval);
               if( SCIPvarGetUbLocal(pseudocands[c]) > newub + 0.5 )
               {
                  SCIP_CALL( SCIPchgVarUb(scip, pseudocands[c], newub) );
                  *result = SCIP_REDUCEDDOM;
                  SCIPdebugMessage(" -> variable <%s> is infeasible in upward branch\n", SCIPvarGetName(pseudocands[c]));
                  break; /* terminate initialization loop, because LP was changed */
               }
            }
         }
         else if( allcolsinlp && !exactsolve && downvalid && upvalid )
         {
            SCIP_Real minbound;
               
            /* the minimal lower bound of both children is a proved lower bound of the current subtree */
            minbound = MIN(down, up);
            provedbound = MAX(provedbound, minbound);
         }

         /* check for a better score, if we are within the maximum priority candidates */
         if( c < npriopseudocands )
         {
            if( integral )
            {
               SCIP_Real gains[3];

               gains[0] = downgain;
               gains[1] = 0.0;
               gains[2] = upgain;
               score = SCIPgetBranchScoreMultiple(scip, pseudocands[c], 3, gains);
            }
            else
               score = SCIPgetBranchScore(scip, pseudocands[c], downgain, upgain);

            if( score > bestscore )
            {
               bestpseudocand = c;
               bestdown = down;
               bestup = up;
               bestdownvalid = downvalid;
               bestupvalid = upvalid;
               bestscore = score;
            }
         }
         else
            score = 0.0;

         /* update pseudo cost values */
         if( !downinf )
         {
            SCIP_CALL( SCIPupdateVarPseudocost(scip, pseudocands[c],
                  solval-SCIPfeasCeil(scip, solval-1.0), downgain, 1.0) );
         }
         if( !upinf )
         {
            SCIP_CALL( SCIPupdateVarPseudocost(scip, pseudocands[c], 
                  solval-SCIPfeasFloor(scip, solval+1.0), upgain, 1.0) );
         }

         SCIPdebugMessage(" -> var <%s> (solval=%g, downgain=%g, upgain=%g, score=%g) -- best: <%s> (%g)\n",
            SCIPvarGetName(pseudocands[c]), solval, downgain, upgain, score,
            SCIPvarGetName(pseudocands[bestpseudocand]), bestscore);
      }

      /* remember last evaluated candidate */
      branchruledata->lastcand = c;

      /* end strong branching */
      SCIP_CALL( SCIPendStrongbranch(scip) );
   }

   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM && *result != SCIP_CONSADDED )
   {
      SCIP_NODE* downchild;
      SCIP_NODE* eqchild;
      SCIP_NODE* upchild;
      SCIP_VAR* var;

      assert(*result == SCIP_DIDNOTRUN);
      assert(0 <= bestpseudocand && bestpseudocand < npseudocands);
      assert(SCIPisLT(scip, provedbound, cutoffbound));

      var = pseudocands[bestpseudocand];

      /* perform the branching */
      SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s>[%g,%g] (solval=%g, down=%g, up=%g, score=%g)\n",
         npseudocands, bestpseudocand, SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetLPSol(var),
         bestdown, bestup, bestscore);
      SCIP_CALL( SCIPbranchVarVal(scip, var, SCIPvarGetLPSol(var), &downchild, &eqchild, &upchild) );

      /* update the lower bounds in the children */
      if( allcolsinlp && !exactsolve )
      {
         if( downchild != NULL )
         {
            SCIP_CALL( SCIPupdateNodeLowerbound(scip, downchild, bestdownvalid ? MAX(bestdown, provedbound) : provedbound) );
            SCIPdebugMessage(" -> down child's lowerbound: %g\n", SCIPnodeGetLowerbound(downchild));
         }
         if( eqchild != NULL )
         {
            SCIP_CALL( SCIPupdateNodeLowerbound(scip, eqchild, provedbound) );
            SCIPdebugMessage(" -> eq child's lowerbound:   %g\n", SCIPnodeGetLowerbound(eqchild));
         }
         if( upchild != NULL )
         {
            SCIP_CALL( SCIPupdateNodeLowerbound(scip, upchild, bestupvalid ? MAX(bestup, provedbound) : provedbound) );
            SCIPdebugMessage(" -> up child's lowerbound:   %g\n", SCIPnodeGetLowerbound(upchild));
         }
      }

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyAllfullstrong)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleAllfullstrong(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeAllfullstrong)
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
SCIP_DECL_BRANCHINIT(branchInitAllfullstrong)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* initialize branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   branchruledata->lastcand = 0;

   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpAllfullstrong)
{  /*lint --e{715}*/
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of allfullstrong branching\n");

   *result = SCIP_DIDNOTRUN;
   
   SCIP_CALL( branch(scip, branchrule, allowaddcons, result) );

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsAllfullstrong)
{  /*lint --e{715}*/
   assert(result != NULL);

   SCIPdebugMessage("Execps method of allfullstrong branching\n");

   *result = SCIP_DIDNOTRUN;

   if( SCIPhasCurrentNodeLP(scip) )
   {
      SCIP_CALL( branch(scip, branchrule, allowaddcons, result) );
   }

   return SCIP_OKAY;
}


/*
 * branching specific interface methods
 */

/** creates the all variables full strong LP branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleAllfullstrong(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create allfullstrong branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
   branchruledata->lastcand = 0;

   /* include allfullstrong branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non-fundamental callbacks via specific setter functions*/
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyAllfullstrong) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeAllfullstrong) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitAllfullstrong) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpAllfullstrong) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsAllfullstrong) );

   return SCIP_OKAY;
}
