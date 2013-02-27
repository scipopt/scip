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

/**@file   branch_pscost.c
 * @brief  pseudo costs branching rule
 * @author Tobias Achterberg
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_pscost.h"

#define BRANCHRULE_NAME          "pscost"
#define BRANCHRULE_DESC          "branching on pseudo cost values"
#define BRANCHRULE_PRIORITY      2000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

#define BRANCHRULE_STRATEGIES          "cdsu" /**< possible pseudo cost multiplication strategies for branching on external candidates */
#define BRANCHRULE_STRATEGY_DEFAULT       'u' /**< default pseudo cost multiplication strategy */
#define BRANCHRULE_SCOREMINWEIGHT_DEFAULT 0.8 /**< default weight for minimum of scores of a branching candidate */
#define BRANCHRULE_SCOREMAXWEIGHT_DEFAULT 1.3 /**< default weight for maximum of scores of a branching candidate */
#define BRANCHRULE_SCORESUMWEIGHT_DEFAULT 0.1 /**< default weight for sum of scores of a branching candidate */
#define BRANCHRULE_NCHILDREN_DEFAULT        2 /**< default number of children in n-ary branching */
#define BRANCHRULE_NARYMAXDEPTH_DEFAULT    -1 /**< default maximal depth where to do n-ary branching */
#define BRANCHRULE_NARYMINWIDTH_DEFAULT 0.001 /**< default minimal domain width in children when doing n-ary branching */
#define BRANCHRULE_NARYWIDTHFAC_DEFAULT   2.0 /**< default factor of domain width in n-ary branching */

#define WEIGHTEDSCORING(data, min, max, sum) \
   ((data)->scoreminweight * (min) + (data)->scoremaxweight * (max) + (data)->scoresumweight * (sum))

/** branching rule data */
struct SCIP_BranchruleData
{
   char                  strategy;           /**< strategy for computing score of external candidates */
   SCIP_Real             scoreminweight;     /**< weight for minimum of scores of a branching candidate */
   SCIP_Real             scoremaxweight;     /**< weight for maximum of scores of a branching candidate */
   SCIP_Real             scoresumweight;     /**< weight for sum of scores of a branching candidate */

   char                  updatestrategy;     /**< strategy used to update pseudo costs of continuous variables */

   int                   nchildren;          /**< targeted number of children in n-ary branching */
   int                   narymaxdepth;       /**< maximal depth where to do n-ary branching, -1 to turn off */
   SCIP_Real             naryminwidth;       /**< minimal domain width in children when doing n-ary branching, relative to global bounds */
   SCIP_Real             narywidthfactor;    /**< factor of domain width in n-ary branching */
};

/*
 * Local methods
 */

/** checks if a given branching candidate is better than a previous one and updates the best branching candidate accordingly */
static
SCIP_RETCODE updateBestCandidate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< branching rule data */
   SCIP_VAR**            bestvar,            /**< best branching candidate */
   SCIP_Real*            bestbrpoint,        /**< branching point for best branching candidate */
   SCIP_Real*            bestscore,          /**< score of best branching candidate */
   SCIP_VAR*             cand,               /**< branching candidate to consider */
   SCIP_Real             candscoremin,       /**< minimal score of branching candidate */
   SCIP_Real             candscoremax,       /**< maximal score of branching candidate */
   SCIP_Real             candscoresum,       /**< sum of scores of branching candidate */
   SCIP_Real             candsol             /**< proposed branching point of branching candidate */          
)
{
   SCIP_Real candbrpoint;
   SCIP_Real branchscore;

   SCIP_Real deltaminus;
   SCIP_Real deltaplus;

   SCIP_Real pscostdown;
   SCIP_Real pscostup;
   
   char strategy;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(bestvar != NULL);
   assert(bestbrpoint != NULL);
   assert(bestscore != NULL);
   assert(cand != NULL);

   /* a branching variable candidate should either be an active problem variable or a multi-aggregated variable */
   assert(SCIPvarIsActive(SCIPvarGetProbvar(cand)) ||
      SCIPvarGetStatus(SCIPvarGetProbvar(cand)) == SCIP_VARSTATUS_MULTAGGR);
   
   if( SCIPvarGetStatus(SCIPvarGetProbvar(cand)) == SCIP_VARSTATUS_MULTAGGR )
   {
      /* for a multi-aggregated variable, we call updateBestCandidate function recursively with all variables in the multi-aggregation */
      SCIP_VAR** multvars;
      int nmultvars;
      int i;
      SCIP_Bool success;
      SCIP_Real multvarlb;
      SCIP_Real multvarub;

      cand = SCIPvarGetProbvar(cand);
      multvars = SCIPvarGetMultaggrVars(cand);
      nmultvars = SCIPvarGetMultaggrNVars(cand);

      /* if we have a candidate branching point, then first register only aggregation variables
       * for which we can compute a corresponding branching point too (see also comments below)
       * if this fails, then register all (unfixed) aggregation variables, thereby forgetting about candsol
       */
      success = FALSE;
      if( candsol != SCIP_INVALID ) /*lint !e777*/
      {
         SCIP_Real* multscalars;
         SCIP_Real minact;
         SCIP_Real maxact;
         SCIP_Real aggrvarsol;
         SCIP_Real aggrvarsol1;
         SCIP_Real aggrvarsol2;

         multscalars = SCIPvarGetMultaggrScalars(cand);

         /* for computing the branching point, we need the current bounds of the multi-aggregated variable */
         minact = SCIPcomputeVarLbLocal(scip, cand);
         maxact = SCIPcomputeVarUbLocal(scip, cand);

         for( i = 0; i < nmultvars; ++i )
         {
            /* skip fixed variables */
            multvarlb = SCIPcomputeVarLbLocal(scip, multvars[i]);
            multvarub = SCIPcomputeVarUbLocal(scip, multvars[i]);
            if( SCIPisEQ(scip, multvarlb, multvarub) )
               continue;

            assert(multscalars != NULL);
            assert(multscalars[i] != 0.0);

            /* we cannot ensure that both the upper bound in the left node and the lower bound in the right node
             * will be candsol by a clever choice for the branching point of multvars[i],
             * but we can try to ensure that at least one of them will be at candsol
             */
            if( multscalars[i] > 0.0 )
            {
               /*    cand >= candsol
                * if multvars[i] >= (candsol - (maxact - multscalars[i] * ub(multvars[i]))) / multscalars[i]
                *                 = (candsol - maxact) / multscalars[i] + ub(multvars[i])
                */
               aggrvarsol1 = (candsol - maxact) / multscalars[i] + multvarub;

               /*     cand <= candsol
                * if multvars[i] <= (candsol - (minact - multscalar[i] * lb(multvars[i]))) / multscalars[i]
                *                 = (candsol - minact) / multscalars[i] + lb(multvars[i])
                */
               aggrvarsol2 = (candsol - minact) / multscalars[i] + multvarlb;
            }
            else
            {
               /*    cand >= candsol
                * if multvars[i] <= (candsol - (maxact - multscalars[i] * lb(multvars[i]))) / multscalars[i]
                *                 = (candsol - maxact) / multscalars[i] + lb(multvars[i])
                */
               aggrvarsol2 = (candsol - maxact) / multscalars[i] + multvarlb;

               /*    cand <= candsol
                * if multvars[i] >= (candsol - (minact - multscalar[i] * ub(multvars[i]))) / multscalars[i]
                *                 = (candsol - minact) / multscalars[i] + ub(multvars[i])
                */
               aggrvarsol1 = (candsol - minact) / multscalars[i] + multvarub;
            }

            /* by the above choice, aggrvarsol1 <= ub(multvars[i]) and aggrvarsol2 >= lb(multvars[i])
             * if aggrvarsol1 <= lb(multvars[i]) or aggrvarsol2 >= ub(multvars[i]), then choose the other one
             * if both are out of bounds, then give up
             * if both are inside bounds, then choose the one closer to 0.0 (someone has better idea???)
             */
            if( SCIPisFeasLE(scip, aggrvarsol1, multvarlb) )
            {
               if( SCIPisFeasGE(scip, aggrvarsol2, multvarub) )
                  continue;
               else
                  aggrvarsol = aggrvarsol2;
            }
            else
            {
               if( SCIPisFeasGE(scip, aggrvarsol2, multvarub) )
                  aggrvarsol = aggrvarsol1;
               else
                  aggrvarsol = REALABS(aggrvarsol1) < REALABS(aggrvarsol2) ? aggrvarsol1 : aggrvarsol2;
            }
            success = TRUE;

            SCIP_CALL( updateBestCandidate(scip, branchruledata, bestvar, bestbrpoint, bestscore,
                  multvars[i], candscoremin, candscoremax, candscoresum, aggrvarsol) );
         }
      }

      if( !success )
         for( i = 0; i < nmultvars; ++i )
         {
            /* skip fixed variables */
            multvarlb = SCIPcomputeVarLbLocal(scip, multvars[i]);
            multvarub = SCIPcomputeVarUbLocal(scip, multvars[i]);
            if( SCIPisEQ(scip, multvarlb, multvarub) )
               continue;

            SCIP_CALL( updateBestCandidate(scip, branchruledata, bestvar, bestbrpoint, bestscore,
               multvars[i], candscoremin, candscoremax, candscoresum, SCIP_INVALID) );
         }

      assert(*bestvar != NULL); /* if all variables were fixed, something is strange */
      
      return SCIP_OKAY;
   }
   
   /* select branching point for this variable */
   candbrpoint = SCIPgetBranchingPoint(scip, cand, candsol);
   assert(candbrpoint >= SCIPvarGetLbLocal(cand));
   assert(candbrpoint <= SCIPvarGetUbLocal(cand));
   assert(SCIPvarGetType(cand) == SCIP_VARTYPE_CONTINUOUS || !SCIPisIntegral(scip, candbrpoint));

   if( SCIPvarGetType(cand) == SCIP_VARTYPE_CONTINUOUS )
      strategy = (branchruledata->strategy == 'u' ? branchruledata->updatestrategy : branchruledata->strategy);
   else
      strategy = (branchruledata->strategy == 'u' ? 'l' : branchruledata->strategy);

   switch( strategy )
   {
   case 'l':
      if( SCIPisInfinity(scip,  SCIPgetSolVal(scip, NULL, cand)) || SCIPgetSolVal(scip, NULL, cand) <= SCIPadjustedVarUb(scip, cand, candbrpoint) )
         deltaminus = 0.0;
      else
         deltaminus = SCIPgetSolVal(scip, NULL, cand) - SCIPadjustedVarUb(scip, cand, candbrpoint);
      if( SCIPisInfinity(scip, -SCIPgetSolVal(scip, NULL, cand)) || SCIPgetSolVal(scip, NULL, cand) >= SCIPadjustedVarLb(scip, cand, candbrpoint) )
         deltaplus = 0.0;
      else
         deltaplus = SCIPadjustedVarLb(scip, cand, candbrpoint) - SCIPgetSolVal(scip, NULL, cand);
      break;

   case 'd':
      if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) )
         deltaminus = SCIPisInfinity(scip, candscoremax) ? SCIPinfinity(scip) : WEIGHTEDSCORING(branchruledata, candscoremin, candscoremax, candscoresum);
      else
         deltaminus = SCIPadjustedVarUb(scip, cand, candbrpoint) - SCIPvarGetLbLocal(cand);

      if( SCIPisInfinity(scip,  SCIPvarGetUbLocal(cand)) )
         deltaplus = SCIPisInfinity(scip, candscoremax) ? SCIPinfinity(scip) : WEIGHTEDSCORING(branchruledata, candscoremin, candscoremax, candscoresum);
      else
         deltaplus = SCIPvarGetUbLocal(cand) - SCIPadjustedVarLb(scip, cand, candbrpoint);
      break;
      
   case 's':
      if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) )
         deltaplus = SCIPisInfinity(scip, candscoremax) ? SCIPinfinity(scip) : WEIGHTEDSCORING(branchruledata, candscoremin, candscoremax, candscoresum);
      else
         deltaplus = SCIPadjustedVarUb(scip, cand, candbrpoint) - SCIPvarGetLbLocal(cand);

      if( SCIPisInfinity(scip,  SCIPvarGetUbLocal(cand)) )
         deltaminus = SCIPisInfinity(scip, candscoremax) ? SCIPinfinity(scip) : WEIGHTEDSCORING(branchruledata, candscoremin, candscoremax, candscoresum);
      else
         deltaminus = SCIPvarGetUbLocal(cand) - SCIPadjustedVarLb(scip, cand, candbrpoint);
      break;

   case 'v':
      deltaplus = SCIPisInfinity(scip, candscoremax) ? SCIPinfinity(scip) : WEIGHTEDSCORING(branchruledata, candscoremin, candscoremax, candscoresum);
      deltaminus = deltaplus;
      break;

   default :
      SCIPerrorMessage("branching strategy %c unknown\n", strategy);
      SCIPABORT();
      return SCIP_INVALIDDATA;  /*lint !e527*/
   }

   if( SCIPisInfinity(scip, deltaminus) || SCIPisInfinity(scip, deltaplus) )
   {
      branchscore = SCIPinfinity(scip);
   }
   else
   {
      pscostdown  = SCIPgetVarPseudocostVal(scip, cand, -deltaminus);
      pscostup    = SCIPgetVarPseudocostVal(scip, cand,  deltaplus);
      branchscore = SCIPgetBranchScore(scip, cand, pscostdown, pscostup);
      assert(!SCIPisNegative(scip, branchscore));
   }
   SCIPdebugMessage("branching score variable <%s>[%g,%g] = %g; score = %g; type=%d bestbrscore=%g\n",
      SCIPvarGetName(cand), SCIPvarGetLbLocal(cand), SCIPvarGetUbLocal(cand), branchscore, WEIGHTEDSCORING(branchruledata, candscoremin, candscoremax, candscoresum),
      SCIPvarGetType(cand), *bestscore);

   if( SCIPisInfinity(scip, branchscore) )
      branchscore = 0.9*SCIPinfinity(scip);
   
   if( SCIPisSumGT(scip, branchscore, *bestscore) )
   {
      (*bestscore)   = branchscore;
      (*bestvar)     = cand;
      (*bestbrpoint) = candbrpoint;
   }
   else if( SCIPisSumEQ(scip, branchscore, *bestscore)
      && !(SCIPisInfinity(scip, -SCIPvarGetLbLocal(*bestvar)) && SCIPisInfinity(scip, SCIPvarGetUbLocal(*bestvar))) )
   {
      /* if best candidate so far is not unbounded to both sides, maybe take new candidate */
      if( (SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand))     || SCIPisInfinity(scip, SCIPvarGetUbLocal(cand))) &&
          (SCIPisInfinity(scip, -SCIPvarGetLbLocal(*bestvar)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(*bestvar))) )
      { 
         /* if both variables are unbounded but one of them is bounded on one side, take the one with the larger bound on this side (hope that this avoids branching on always the same variable) */
         if( SCIPvarGetUbLocal(cand) > SCIPvarGetUbLocal(*bestvar) || SCIPvarGetLbLocal(cand) < SCIPvarGetLbLocal(*bestvar) )
         {
            (*bestvar)     = cand;
            (*bestbrpoint) = candbrpoint;
         }
      }
      else if( SCIPvarGetType(*bestvar) == SCIPvarGetType(cand) )
      { 
         /* if both have the same type, take the one with larger diameter */
         if( SCIPvarGetUbLocal(*bestvar) - SCIPvarGetLbLocal(*bestvar) < SCIPvarGetUbLocal(cand) - SCIPvarGetLbLocal(cand) )
         {
            (*bestvar)     = cand;
            (*bestbrpoint) = candbrpoint;
         }
      }
      else if( SCIPvarGetType(*bestvar) > SCIPvarGetType(cand) )
      { 
         /* take the one with better type ("more discrete") */
         (*bestvar)     = cand;
         (*bestbrpoint) = candbrpoint;
      }
   }

   return SCIP_OKAY;
}

/** selects the branching variable from given candidate array */
static
SCIP_RETCODE selectBranchVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_VAR**            cands,              /**< array of branching candidates */
   SCIP_Real*            candssol,           /**< array of candidate solution values */
   SCIP_Real*            candsscore,         /**< array of candidate scores */
   int                   ncands,             /**< the number of candidates */
   SCIP_VAR**            brvar,              /**< pointer to store the selected branching candidate or NULL if none */
   SCIP_Real*            brpoint             /**< pointer to store branching point of selected branching variable */
   )
{ /*lint --e{850}*/ 
   SCIP_BRANCHRULEDATA* branchruledata;

   SCIP_VAR* cand;
   SCIP_Real candsol;

   SCIP_Real bestbranchscore;

   SCIP_Real scoremin;
   SCIP_Real scoresum;
   SCIP_Real scoremax;

   SCIP_VAR** candssorted;
   int* candsorigidx;
   
   int i;
   int j;
   
   assert(brvar   != NULL);
   assert(brpoint != NULL);
   
   (*brvar)   = NULL;
   (*brpoint) = SCIP_INVALID;

   if( ncands == 0 )
      return SCIP_OKAY;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   
   /* sort branching candidates (in a copy), such that same variables are on consecutive positions */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &candssorted, cands, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candsorigidx, ncands) );
   for( i = 0; i < ncands; ++i )
      candsorigidx[i] = i;
   
   SCIPsortPtrInt((void**)candssorted, candsorigidx, SCIPvarComp, ncands);

   bestbranchscore = -1.0;

   for( i = 0; i < ncands; ++i )
   {
      cand = candssorted[i];

      /* there should be no fixed branching candidates */
      assert(!SCIPisEQ(scip, SCIPvarGetLbLocal(cand), SCIPvarGetUbLocal(cand)));

      /* compute min, sum, and max of all registered scores for this variables
       * set candsol to a valid value, if someone registered one */
      scoremin = candsscore[candsorigidx[i]];
      scoresum = scoremin;
      scoremax = scoremin;
      candsol  = candssol[candsorigidx[i]];
      for( j = i+1 ; j < ncands && SCIPvarCompare(candssorted[j], cand) == 0; ++j )
      {
         assert(candsscore[candsorigidx[j]] >= 0.0);
         scoresum += candsscore[candsorigidx[j]];
         if( candsscore[candsorigidx[j]] < scoremin )
            scoremin = candsscore[candsorigidx[j]];
         else if( candsscore[candsorigidx[j]] > scoremax )
            scoremax = candsscore[candsorigidx[j]];

         /* @todo if there are two valid externcandssol available for the same variable, should we take the one closer to the middle of the domain? */
         if( SCIPisInfinity(scip, REALABS(candsol)) )
            candsol = candssol[candsorigidx[j]];
      }
      /* set i to last occurrence of cand in candssorted (instead of first one as before), so in next round we look at another variable */
      i = j-1;
      assert(candssorted[i] == cand);
      
      /* check if new candidate is better than previous candidate (if any) */
      SCIP_CALL( updateBestCandidate(scip, branchruledata, brvar, brpoint, &bestbranchscore, cand, scoremin, scoremax, scoresum, candsol) );
      assert(*brvar != NULL);
   }
   
   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &candssorted);
   SCIPfreeBufferArray(scip, &candsorigidx);

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyPscost)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchrulePscost(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreePscost)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpPscost)
{  /*lint --e{715}*/
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real bestscore;
   SCIP_Real bestrootdiff;
   int nlpcands;
   int bestcand;
   int c;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of pscost branching\n");

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, NULL, &nlpcands) );
   assert(nlpcands > 0);

   bestcand = -1;
   bestscore = -SCIPinfinity(scip);
   bestrootdiff = 0.0;
   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_Real score;
      SCIP_Real rootsolval;
      SCIP_Real rootdiff;

      score = SCIPgetVarPseudocostScore(scip, lpcands[c], lpcandssol[c]);
      rootsolval = SCIPvarGetRootSol(lpcands[c]);
      rootdiff = REALABS(lpcandssol[c] - rootsolval);
      if( SCIPisSumGT(scip, score, bestscore) || (SCIPisSumEQ(scip, score, bestscore) && rootdiff > bestrootdiff) )
      {
         bestcand = c;
         bestscore = score;
         bestrootdiff = rootdiff;
      }
   }
   assert(0 <= bestcand && bestcand < nlpcands);
   assert(!SCIPisFeasIntegral(scip, lpcandssol[bestcand]));

   /* perform the branching */
   SCIPdebugMessage(" -> %d cands, selected cand %d: variable <%s> (solval=%g, score=%g)\n",
      nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandssol[bestcand], bestscore);

   /* perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, lpcands[bestcand], NULL, NULL, NULL) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** branching execution method for external candidates */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextPscost)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** externcands;
   SCIP_Real* externcandssol;
   SCIP_Real* externcandsscore;
   int nprioexterncands;
   SCIP_VAR* brvar;
   SCIP_Real brpoint;
   int nchildren;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPdebugMessage("Execext method of pscost branching\n");
   
   /* get branching candidates */
   SCIP_CALL( SCIPgetExternBranchCands(scip, &externcands, &externcandssol, &externcandsscore, NULL, &nprioexterncands, NULL, NULL, NULL) );
   assert(nprioexterncands > 0);
   
   /* get current update strategy for pseudo costs, if our multiplier rule is 'u' */
   if( branchruledata->strategy == 'u' )
   {
      SCIP_CALL( SCIPgetCharParam(scip, "branching/lpgainnormalize", &branchruledata->updatestrategy) );
   }

   /* select branching variable */
   SCIP_CALL( selectBranchVar(scip, branchrule, externcands, externcandssol, externcandsscore, nprioexterncands, &brvar, &brpoint) );
   
   if( brvar == NULL )
   {
      SCIPerrorMessage("branchExecextPscost failed to select a branching variable from %d candidates\n", nprioexterncands);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
  
   assert(SCIPvarIsActive(SCIPvarGetProbvar(brvar)));

   SCIPdebugMessage("branching on variable <%s>: new intervals: [%g, %g] and [%g, %g]\n",
      SCIPvarGetName(brvar), SCIPvarGetLbLocal(brvar), SCIPadjustedVarUb(scip, brvar, brpoint), SCIPadjustedVarLb(scip, brvar, brpoint), SCIPvarGetUbLocal(brvar));

   if( branchruledata->nchildren > 2 && SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) <= branchruledata->narymaxdepth )
   {
      /* do n-ary branching */
      SCIP_Real minwidth;

      minwidth = 0.0;
      if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(brvar)) && !SCIPisInfinity(scip, SCIPvarGetUbGlobal(brvar)) )
         minwidth = branchruledata->naryminwidth * (SCIPvarGetUbGlobal(brvar) - SCIPvarGetLbGlobal(brvar));

      SCIP_CALL( SCIPbranchVarValNary(scip, brvar, brpoint, branchruledata->nchildren, minwidth, branchruledata->narywidthfactor, &nchildren) );
   }
   else
   {
      /* do binary branching */
      SCIP_CALL( SCIPbranchVarValNary(scip, brvar, brpoint, 2, 0.0, 1.0, &nchildren) );
   }

   if( nchildren > 1 )
   {
      *result = SCIP_BRANCHED;
   }
   else
   {
      /* if there are no children, then variable should have been fixed by SCIPbranchVarVal */
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(brvar), SCIPvarGetUbLocal(brvar)));
      *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}


/*
 * branching specific interface methods
 */

/** creates the pseudo cost branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchrulePscost(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create pscost branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
   
   /* include allfullstrong branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non-fundamental callbacks via specific setter functions*/
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyPscost) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreePscost) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpPscost) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextPscost) );

   SCIP_CALL( SCIPaddCharParam(scip, "branching/"BRANCHRULE_NAME"/strategy",
         "strategy for utilizing pseudo-costs of external branching candidates (multiply as in pseudo costs 'u'pdate rule, or by 'd'omain reduction, or by domain reduction of 's'ibling, or by 'v'ariable score)",
         &branchruledata->strategy, FALSE, BRANCHRULE_STRATEGY_DEFAULT, BRANCHRULE_STRATEGIES, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/minscoreweight",
         "weight for minimum of scores of a branching candidate when building weighted sum of min/max/sum of scores",
         &branchruledata->scoreminweight, TRUE, BRANCHRULE_SCOREMINWEIGHT_DEFAULT, -SCIPinfinity(scip), SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/maxscoreweight",
         "weight for maximum of scores of a branching candidate when building weighted sum of min/max/sum of scores",
         &branchruledata->scoremaxweight, TRUE, BRANCHRULE_SCOREMAXWEIGHT_DEFAULT, -SCIPinfinity(scip), SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/sumscoreweight",
         "weight for sum of scores of a branching candidate when building weighted sum of min/max/sum of scores",
         &branchruledata->scoresumweight, TRUE, BRANCHRULE_SCORESUMWEIGHT_DEFAULT, -SCIPinfinity(scip), SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "branching/"BRANCHRULE_NAME"/nchildren",
         "number of children to create in n-ary branching",
         &branchruledata->nchildren, FALSE, BRANCHRULE_NCHILDREN_DEFAULT, 2, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "branching/"BRANCHRULE_NAME"/narymaxdepth",
         "maximal depth where to do n-ary branching, -1 to turn off",
         &branchruledata->narymaxdepth, FALSE, BRANCHRULE_NARYMAXDEPTH_DEFAULT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/naryminwidth",
         "minimal domain width in children when doing n-ary branching, relative to global bounds",
         &branchruledata->naryminwidth, FALSE, BRANCHRULE_NARYMINWIDTH_DEFAULT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/narywidthfactor",
         "factor of domain width in n-ary branching when creating nodes with increasing distance from branching value",
         &branchruledata->narywidthfactor, FALSE, BRANCHRULE_NARYWIDTHFAC_DEFAULT, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** selects a branching variable, due to pseudo cost, from the given candidate array and returns this variable together
 *  with a branching point */
SCIP_RETCODE SCIPselectBranchVarPscost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            branchcands,        /**< branching candidates */
   SCIP_Real*            branchcandssol,     /**< solution value for the branching candidates */
   SCIP_Real*            branchcandsscore,   /**< array of candidate scores */
   int                   nbranchcands,       /**< number of branching candidates */
   SCIP_VAR**            var,                /**< pointer to store the variable to branch on, or NULL if none */
   SCIP_Real*            brpoint             /**< pointer to store the branching point for the branching variable, will be fractional for a discrete variable */
   )
{
   SCIP_BRANCHRULE* branchrule;
   
   assert(scip != NULL);
   
   /* find branching rule */
   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);
   
   /* select branching variable */
   SCIP_CALL( selectBranchVar(scip, branchrule, branchcands, branchcandssol, branchcandsscore, nbranchcands, var, brpoint) );

   return SCIP_OKAY;
}
