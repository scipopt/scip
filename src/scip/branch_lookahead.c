/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*#define SCIP_DEBUG*/
/**@file   branch_lookahead.c
 * @brief  lookahead branching rule
 * @author Christoph Schubert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/branch_lookahead.h"

#include <assert.h>
#include <string.h>

#include "def.h"
#include "pub_branch.h"
#include "pub_message.h"
#include "pub_var.h"
#include "scip.h"
#include "type_branch.h"
#include "type_lp.h"
#include "type_result.h"
#include "type_retcode.h"
#include "type_scip.h"
#include "type_tree.h"
#include "type_var.h"

#define BRANCHRULE_NAME            "lookahead"
#define BRANCHRULE_DESC            "fullstrong branching with depth of 2" /* TODO CS: expand description */
#define BRANCHRULE_PRIORITY        536870911
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

/*
 * Data structures
 */

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Bool somerandomfield;
};

enum BoundStatus
{
   BOUNDSTATUS_NONE = 0,
   BOUNDSTATUS_UPPERBOUND,
   BOUNDSTATUS_LOWERBOUND,
   BOUNDSTATUS_BOTH
};
typedef enum BoundStatus BOUNDSTATUS;

typedef struct
{
   SCIP_Real             highestweight;
   SCIP_Real             sumofweights;
   int                   numberofweights;
} WeightData;

typedef struct
{
   int                   varindex;
   int                   ncutoffs;
   WeightData            upperbounddata;
   WeightData            lowerbounddata;
} ScoreData;

typedef struct
{
   SCIP_Real             objval;
   SCIP_Bool             cutoff;
   SCIP_Bool             lperror;
   SCIP_Bool             nobranch;
} BranchingResultData;

typedef struct
{
   SCIP_Real*            upperbounds;
   SCIP_Real*            lowerbounds;
   BOUNDSTATUS*          boundstatus;
   int*                  innervars;
   int                   ninnervars;
} InnerBoundData;

/*
 * Local methods
 */
static
SCIP_RETCODE initWeightData(
   WeightData*           weightdata
   )
{
   weightdata->highestweight = 0;
   weightdata->numberofweights = 0;
   weightdata->sumofweights = 0;
   return SCIP_OKAY;
}

static
SCIP_RETCODE initScoreData(
   ScoreData*            scoredata,
   int                   currentbranchvar
)
{
   scoredata->ncutoffs = 0;
   scoredata->varindex = currentbranchvar;
   SCIP_CALL( initWeightData(&scoredata->lowerbounddata) );
   SCIP_CALL( initWeightData(&scoredata->upperbounddata) );
   return SCIP_OKAY;
}

static
SCIP_RETCODE initBranchingResultData(
   SCIP*                 scip,
   BranchingResultData*  resultdata
)
{
   resultdata->objval = SCIPinfinity(scip);
   resultdata->cutoff = TRUE;
   resultdata->lperror = FALSE;
   resultdata->nobranch = FALSE;
   return SCIP_OKAY;
}

static
SCIP_RETCODE allocInnerBoundData(
   SCIP*                 scip,
   InnerBoundData**      innerbounddata
)
{
   int ntotalvars;

   ntotalvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBuffer(scip, innerbounddata) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*innerbounddata)->upperbounds, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*innerbounddata)->lowerbounds, ntotalvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &(*innerbounddata)->boundstatus, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*innerbounddata)->innervars, ntotalvars) );
   return SCIP_OKAY;
}

static
SCIP_RETCODE initInnerBoundData(
   InnerBoundData*       innerbounddata
)
{
   innerbounddata->ninnervars = 0;
   return SCIP_OKAY;
}

static
SCIP_RETCODE freeInnerBoundData(
   SCIP*                 scip,
   InnerBoundData**      innerbounddata
)
{
   SCIPfreeBufferArray(scip, &(*innerbounddata)->innervars);
   SCIPfreeCleanBufferArray(scip, &(*innerbounddata)->boundstatus);
   SCIPfreeBufferArray(scip, &(*innerbounddata)->lowerbounds);
   SCIPfreeBufferArray(scip, &(*innerbounddata)->upperbounds);
   SCIPfreeBuffer(scip, innerbounddata);
   return SCIP_OKAY;
}

/**
 * Executes the branching on the current probing node by adding a probing node with a new upper bound.
 */
static
SCIP_RETCODE executeBranchingOnUpperBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< variable to branch on */
   SCIP_Real             branchvarsolval,    /**< current (fractional) solution value of the variable */
   BranchingResultData*  resultdata
   )
{
   SCIP_Real oldupperbound;
   SCIP_Real oldlowerbound;
   SCIP_Real newupperbound;
   SCIP_LPSOLSTAT solstat;

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(!SCIPisFeasIntegral(scip, branchvarsolval));
   assert(resultdata != NULL);

   newupperbound = SCIPfeasFloor(scip, branchvarsolval);
   oldupperbound = SCIPvarGetUbLocal(branchvar);
   oldlowerbound = SCIPvarGetLbLocal(branchvar);

   SCIPdebugMessage("New upper bound: <%g>, old lower bound: <%g>, old upper bound: <%g>\n", newupperbound, oldlowerbound,
      oldupperbound);

   if( SCIPisFeasLT(scip, newupperbound, oldlowerbound) )
   {
      resultdata->lperror = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      if( SCIPisEQ(scip, oldupperbound, oldlowerbound) )
      {
         /* TODO: do smth with this info. */
         resultdata->nobranch = TRUE;
      }
      else if( SCIPisFeasLT(scip, newupperbound, oldupperbound) )
      {
         /* if the new upper bound is lesser than the old upper bound and also
          * greater than (or equal to) the old lower bound we set the new upper bound.
          * oldLowerBound <= newUpperBound < oldUpperBound */
         SCIP_CALL( SCIPchgVarUbProbing(scip, branchvar, newupperbound) );
      }

      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &resultdata->lperror, &resultdata->cutoff) );
      solstat = SCIPgetLPSolstat(scip);

      resultdata->lperror = resultdata->lperror || (solstat == SCIP_LPSOLSTAT_NOTSOLVED && resultdata->cutoff == FALSE) ||
            (solstat == SCIP_LPSOLSTAT_ITERLIMIT) || (solstat == SCIP_LPSOLSTAT_TIMELIMIT);
      assert(solstat != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

      if( !resultdata->lperror )
      {
         resultdata->objval = SCIPgetLPObjval(scip);
         resultdata->cutoff = resultdata->cutoff || SCIPisGE(scip, resultdata->objval, SCIPgetCutoffbound(scip));
         assert(((solstat != SCIP_LPSOLSTAT_INFEASIBLE) && (solstat != SCIP_LPSOLSTAT_OBJLIMIT)) || resultdata->cutoff);
      }
   }

   return SCIP_OKAY;
}

/**
 * Executes the branching on the current probing node by adding a probing node with a new lower bound.
 */
static
SCIP_RETCODE executeBranchingOnLowerBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< variable to branch on */
   SCIP_Real             branchvarsolval,    /**< current (fractional) solution value of the variable */
   BranchingResultData*  resultdata
   )
{
   SCIP_Real oldlowerbound;
   SCIP_Real oldupperbound;
   SCIP_Real newlowerbound;
   SCIP_LPSOLSTAT solstat;

   assert(scip != NULL );
   assert(branchvar != NULL );
   assert(!SCIPisFeasIntegral(scip, branchvarsolval));
   assert(resultdata != NULL );

   newlowerbound = SCIPfeasCeil(scip, branchvarsolval);
   oldlowerbound = SCIPvarGetLbLocal(branchvar);
   oldupperbound = SCIPvarGetUbLocal(branchvar);

   SCIPdebugMessage("New lower bound: <%g>, old lower bound: <%g>, old upper bound: <%g>\n", newlowerbound, oldlowerbound,
      oldupperbound);

   if( SCIPisFeasGT(scip, newlowerbound, oldupperbound) )
   {
      resultdata->cutoff = TRUE;
      resultdata->lperror = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      if( SCIPisEQ(scip, oldupperbound, oldlowerbound) )
      {
         /* TODO: do smth with this info. */
         resultdata->nobranch = TRUE;
      }
      else if( SCIPisFeasGT(scip, newlowerbound, oldlowerbound) )
      {
         /* if the new lower bound is greater than the old lower bound and also
          * lesser than (or equal to) the old upper bound we set the new lower bound.
          * oldLowerBound < newLowerBound <= oldUpperBound */
         SCIP_CALL( SCIPchgVarLbProbing(scip, branchvar, newlowerbound) );
      }

      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &resultdata->lperror, &resultdata->cutoff) );
      solstat = SCIPgetLPSolstat(scip);

      resultdata->lperror = resultdata->lperror || (solstat == SCIP_LPSOLSTAT_NOTSOLVED && resultdata->cutoff == FALSE) ||
            (solstat == SCIP_LPSOLSTAT_ITERLIMIT) || (solstat == SCIP_LPSOLSTAT_TIMELIMIT);
      assert(solstat != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

      if( !resultdata->lperror )
      {
         resultdata->objval = SCIPgetLPObjval(scip);
         resultdata->cutoff = resultdata->cutoff || SCIPisGE(scip, resultdata->objval, SCIPgetCutoffbound(scip));
         assert(((solstat != SCIP_LPSOLSTAT_INFEASIBLE) && (solstat != SCIP_LPSOLSTAT_OBJLIMIT)) || resultdata->cutoff);
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE calculateWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             upgain,
   SCIP_Real             downgain,
   SCIP_Real*            result
)
{
   SCIP_Real min;
   SCIP_Real max;
   SCIP_Real minweight = 4;
   SCIP_Real maxweight = 1;

   assert(scip != NULL);
   assert(result != NULL);

   min = MIN(downgain, upgain);
   max = MAX(upgain, downgain);

   *result = minweight * min + maxweight * max;

   SCIPdebugMessage("The calculated weight of <%g> and <%g> is <%g>.\n", upgain, downgain, *result);

   return SCIP_OKAY;
}

static
SCIP_RETCODE executeDeepBranchingOnVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lpobjval,           /**< objective value of the base lp */
   SCIP_VAR*             deepbranchvar,      /**< variable to branch up and down on */
   SCIP_Real             deepbranchvarsolval,/**< (fractional) solution value of the branching variable */
   SCIP_Bool*            fullcutoff,         /**< resulting decision whether this branch is cutoff */
   SCIP_Bool*            lperror,
   WeightData*           weightdata,         /**< container to be filled with the weight relevant data */
   int*                  ncutoffs,           /**< current (input) and resulting (output) number of cutoffs */
   InnerBoundData*       innerbounddata
)
{
   BranchingResultData* downresultdata;
   BranchingResultData* upresultdata;
   SCIP_Real downgain;
   SCIP_Real upgain;
   SCIP_Real currentweight;

   assert(scip != NULL);
   assert(deepbranchvar != NULL);
   assert(ncutoffs != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, &downresultdata) );
   SCIP_CALL( SCIPallocBuffer(scip, &upresultdata) );
   SCIP_CALL( initBranchingResultData(scip, downresultdata) );
   SCIP_CALL( initBranchingResultData(scip, upresultdata) );

   SCIPdebugMessage("Second level down branching on variable <%s>\n", SCIPvarGetName(deepbranchvar));
   SCIP_CALL( executeBranchingOnUpperBound(scip, deepbranchvar, deepbranchvarsolval, downresultdata) );

   if( downresultdata->lperror )
   {
      /* Something went wrong while solving the lp. Maybe exceeded the time-/iterlimit or we tried to add an upper
       * bound, which is lower than the current lower bound (may be the case if the lower bound was raised due to
       * propagation from other branches.) */
      *lperror = TRUE;
   }
   else
   {
      SCIPdebugMessage("Going back to layer 1.\n");
      /* go back one layer (we are currently in depth 2) */
      SCIP_CALL( SCIPbacktrackProbing(scip, 1) );

      SCIPdebugMessage("Second level up branching on variable <%s>\n", SCIPvarGetName(deepbranchvar));
      SCIP_CALL( executeBranchingOnLowerBound(scip, deepbranchvar, deepbranchvarsolval, upresultdata) );

      if( upresultdata->lperror )
      {
         /* Something went wrong while solving the lp. Maybe exceeded the time-/iterlimit or we tried to add an upper
          * bound, which is lower than the current lower bound (may be the case if the lower bound was raised due to
          * propagation from other branches.) */
         *lperror = TRUE;
      }
      else
      {
         SCIPdebugMessage("Going back to layer 1.\n");
         /* go back one layer (we are currently in depth 2) */
         SCIP_CALL( SCIPbacktrackProbing(scip, 1) );

         if( !downresultdata->cutoff && !upresultdata->cutoff )
         {
            downgain = downresultdata->objval - lpobjval;
            upgain = upresultdata->objval - lpobjval;

            SCIPdebugMessage("The difference between the objective values of the base lp and the upper bounded lp is <%g>\n",
               downgain);
            SCIPdebugMessage("The difference between the objective values of the base lp and the lower bounded lp is <%g>\n",
               upgain);

            assert(!SCIPisFeasNegative(scip, downgain));
            assert(!SCIPisFeasNegative(scip, upgain));

            SCIP_CALL( calculateWeight(scip, upgain, downgain, &currentweight) );

            weightdata->highestweight = MAX(weightdata->highestweight, currentweight);
            weightdata->sumofweights = weightdata->sumofweights + currentweight;
            weightdata->numberofweights++;

            SCIPdebugMessage("The sum of weights is <%g>.\n", weightdata->sumofweights);
            SCIPdebugMessage("The number of weights is <%i>.\n", weightdata->numberofweights);
            *fullcutoff = FALSE;
         }
         else if( downresultdata->cutoff && upresultdata->cutoff )
         {
            *fullcutoff = TRUE;
            *ncutoffs = *ncutoffs + 2;
         }
         else
         {
            *fullcutoff = FALSE;
            *ncutoffs = *ncutoffs + 1;
         }
      }
   }

   SCIPfreeBuffer(scip, &upresultdata);
   SCIPfreeBuffer(scip, &downresultdata);

   return SCIP_OKAY;
}

static
SCIP_RETCODE executeDeepBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lpobjval,           /**< objective value of the base lp */
   SCIP_Bool*            fullcutoff,         /**< resulting decision whether this branch is cutoff */
   SCIP_Bool*            lperror,
   WeightData*           weightdata,
   int*                  ncutoffs,
   InnerBoundData*       innerbounddata
)
{
   SCIP_VAR**  lpcands;
   SCIP_Real*  lpcandssol;
   int         nlpcands;
   int         j;

   assert(scip != NULL);
   assert(ncutoffs != NULL);

   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );

   SCIPdebugMessage("The deeper lp has <%i> variables with fractional value.\n", nlpcands);

   for( j = 0; j < nlpcands; j++ )
   {
      SCIP_VAR* deepbranchvar = lpcands[j];
      SCIP_Real deepbranchvarsolval = lpcandssol[j];

      SCIPdebugMessage("Start deeper branching on variable <%s> with solution value <%g>.\n",
         SCIPvarGetName(deepbranchvar), deepbranchvarsolval);

      SCIP_CALL( executeDeepBranchingOnVar(scip, lpobjval, deepbranchvar, deepbranchvarsolval,
         fullcutoff, lperror, weightdata, ncutoffs, innerbounddata) );

      if( *fullcutoff )
      {
         SCIPdebugMessage("The deeper lp on variable <%s> is cutoff, as both lps are cutoff.\n",
            SCIPvarGetName(deepbranchvar));
         break;
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE calculateAverageWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   WeightData            weightdata,         /**< calculation data for the average weight */
   SCIP_Real*            averageweight       /**< resulting average weight */
)
{
   assert(scip != NULL);
   assert(!SCIPisFeasNegative(scip, weightdata.sumofweights));
   assert(weightdata.numberofweights >= 0);
   assert(averageweight != NULL);

   if( weightdata.numberofweights > 0 )
   {
      *averageweight = (1 / weightdata.numberofweights) * weightdata.sumofweights;
   }
   else
   {
      *averageweight = 0;
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE calculateCurrentWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   ScoreData             scoredata,
   SCIP_Real*            highestweight,
   int*                  highestweightindex
)
{
   SCIP_Real averageweightupperbound = 0;
   SCIP_Real averageweightlowerbound = 0;
   SCIP_Real lambda;
   SCIP_Real totalweight;

   assert(scip != NULL);
   assert(!SCIPisFeasNegative(scip, scoredata.upperbounddata.highestweight));
   assert(!SCIPisFeasNegative(scip, scoredata.lowerbounddata.highestweight));
   assert(!SCIPisFeasNegative(scip, scoredata.ncutoffs));
   assert(highestweight != NULL);
   assert(highestweightindex != NULL);

   SCIP_CALL( calculateAverageWeight(scip, scoredata.upperbounddata, &averageweightupperbound) );
   SCIP_CALL( calculateAverageWeight(scip, scoredata.lowerbounddata, &averageweightlowerbound) );
   lambda = averageweightupperbound + averageweightlowerbound;

   assert(!SCIPisFeasNegative(scip, lambda));

   SCIPdebugMessage("The lambda value is <%g>.\n", lambda);

   totalweight = scoredata.lowerbounddata.highestweight + scoredata.upperbounddata.highestweight + scoredata.ncutoffs;
   if( SCIPisFeasGT(scip, totalweight, *highestweight) )
   {
      *highestweight = totalweight;
      *highestweightindex = scoredata.varindex;
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE handleNewBounds(
   SCIP*                 scip,
   BOUNDSTATUS*          boundstatus,
   SCIP_Real*            newlowerbounds,
   SCIP_Real*            newupperbounds,
   SCIP_RESULT*          result
)
{
   int i;
   int nprobvars;
   SCIP_VAR** probvars;

   nprobvars = SCIPgetNVars(scip);
   probvars = SCIPgetVars(scip);

   for( i = 0; i < nprobvars; i++ )
   {
      BOUNDSTATUS status;
      SCIP_VAR* branchvar;

      status = boundstatus[i];
      branchvar = probvars[i];
      if( status == BOUNDSTATUS_LOWERBOUND || status == BOUNDSTATUS_BOTH )
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         SCIP_Real newlowerbound = newlowerbounds[i];
         SCIP_CALL( SCIPtightenVarLb(scip, branchvar, newlowerbound, FALSE, &infeasible, &tightened) );
         if( infeasible )
         {
            *result = SCIP_DIDNOTFIND;
         }
         else if( *result != SCIP_DIDNOTFIND && tightened)
         {
            *result = SCIP_REDUCEDDOM;
         }
      }
      if( status == BOUNDSTATUS_UPPERBOUND || status == BOUNDSTATUS_BOTH )
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         SCIP_Real newupperbound = newupperbounds[i];
         SCIP_CALL( SCIPtightenVarUb(scip, branchvar, newupperbound, FALSE, &infeasible, &tightened) );

         if( infeasible )
         {
            *result = SCIP_DIDNOTFIND;
         }
         else if( *result != SCIP_DIDNOTFIND && tightened)
         {
            *result = SCIP_REDUCEDDOM;
         }
      }

      if( status != BOUNDSTATUS_NONE )
      {
         /* Reset the array s.t. it only contains zero values. */
         boundstatus[i] = BOUNDSTATUS_NONE;
      }
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE addNewUpperBound(
   SCIP*                 scip,
   SCIP_VAR*             branchvar,
   SCIP_Real             newupperbound,
   SCIP_Real*            newupperbounds,
   BOUNDSTATUS*          boundstatus
)
{
   int varindex;
   BOUNDSTATUS status;

   varindex = SCIPvarGetProbindex(branchvar);
   status = boundstatus[varindex];

   if( status == BOUNDSTATUS_UPPERBOUND || status == BOUNDSTATUS_BOTH )
   {
      /* we already have a valid new upper bound set, so we have to take the minimum of the previous and the new
       * "newupperbound". */
      SCIP_Real prevnewupperbound = newupperbounds[varindex];
      newupperbounds[varindex] = MIN(newupperbound, prevnewupperbound);
   }
   else
   {
      /* We either have no new bounds or only a lower bound for our var, so we can set the new upper bound directly. */
      newupperbounds[varindex] = newupperbound;
      if( status == BOUNDSTATUS_NONE )
      {
         boundstatus[varindex] = BOUNDSTATUS_UPPERBOUND;
      }
      else
      {
         boundstatus[varindex] = BOUNDSTATUS_BOTH;
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE addNewLowerBound(
   SCIP*                 scip,
   SCIP_VAR*             branchvar,
   SCIP_Real             newlowerbound,
   SCIP_Real*            newlowerbounds,
   BOUNDSTATUS*          boundstatus
)
{
   int varindex;
   BOUNDSTATUS status;

   varindex = SCIPvarGetProbindex(branchvar);
   status = boundstatus[varindex];

   if( status == BOUNDSTATUS_LOWERBOUND || status == BOUNDSTATUS_BOTH )
   {
      /* we already have a valid new lower bound set, so we have to take the minimum of the previous and the new
       * "newlowerbound". */
      SCIP_Real prevnewupperbound = newlowerbounds[varindex];
      newlowerbounds[varindex] = MAX(newlowerbound, prevnewupperbound);
   }
   else
   {
      /* We either have no new bounds or only a upper bound for our var, so we can set the new lower bound directly. */
      newlowerbounds[varindex] = newlowerbound;
      if( status == BOUNDSTATUS_NONE )
      {
         boundstatus[varindex] = BOUNDSTATUS_LOWERBOUND;
      }
      else
      {
         boundstatus[varindex] = BOUNDSTATUS_BOTH;
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE selectVarLookaheadBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            lpcands,            /**< array of fractional variables */
   SCIP_Real*            lpcandssol,         /**< array of fractional solution values */
   int                   nlpcands,           /**< number of fractional variables/solution values */
   int*                  bestcand,           /**< calculated index of the branching variable */
   SCIP_RESULT*          result              /**< pointer to store results of branching */){

   assert(scip != NULL);
   assert(lpcands != NULL);
   assert(lpcandssol != NULL);
   assert(bestcand != NULL);
   assert(result != NULL);

   if( nlpcands == 1)
   {
      /** if there is only one branching variable we can directly branch there */
      *bestcand = 0;
      return SCIP_OKAY;
   }

   if( SCIPgetDepthLimit(scip) <= (SCIPgetDepth(scip) + 2) )
   {
      SCIPdebugMessage("Cannot perform probing in selectVarLookaheadBranching, depth limit reached.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( nlpcands > 1 )
   {
      BranchingResultData* downbranchingresult;
      BranchingResultData* upbranchingresult;
      ScoreData* scoredata;

      int nglobalvars;

      SCIP_Real* newupperbounds;
      SCIP_Real* newlowerbounds;
      BOUNDSTATUS* boundstatus;

      SCIP_Real lpobjval;
      SCIP_Real highestscore = 0;
      int highestscoreindex = -1;
      int i;

      InnerBoundData* innerbounddata;

      nglobalvars = SCIPgetNVars(scip);

      SCIP_CALL( SCIPallocBuffer(scip, &downbranchingresult) );
      SCIP_CALL( SCIPallocBuffer(scip, &upbranchingresult) );
      SCIP_CALL( SCIPallocBuffer(scip, &scoredata) );

      SCIP_CALL( SCIPallocBufferArray(scip, &newupperbounds, nglobalvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newlowerbounds, nglobalvars) );
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &boundstatus, nglobalvars) );

      SCIP_CALL( allocInnerBoundData(scip, &innerbounddata) );

      SCIP_CALL( initBranchingResultData(scip, downbranchingresult) );
      SCIP_CALL( initBranchingResultData(scip, upbranchingresult) );

      lpobjval = SCIPgetLPObjval(scip);

      SCIPdebugMessage("The objective value of the base lp is <%g>.\n", lpobjval);

      SCIP_CALL( SCIPstartProbing(scip) );
      SCIPdebugMessage("Start Probing Mode\n");

      for( i = 0; i < nlpcands && !downbranchingresult->lperror && !upbranchingresult->lperror && !SCIPisStopped(scip); i++ )
      {
         SCIP_VAR* branchvar;
         SCIP_Real branchval;

         SCIP_CALL( initInnerBoundData(innerbounddata) );

         assert(lpcands[i] != NULL);

         branchvar = lpcands[i];
         branchval = lpcandssol[i];

         SCIPdebugMessage("Start branching on variable <%s>\n", SCIPvarGetName(branchvar));

         SCIP_CALL( initScoreData(scoredata, i) );
         scoredata->varindex = i;
         scoredata->ncutoffs = 0;

         SCIPdebugMessage("First level down branching on variable <%s>\n", SCIPvarGetName(branchvar));
         SCIP_CALL( executeBranchingOnUpperBound(scip, branchvar, branchval, downbranchingresult) );

         if( !downbranchingresult->lperror && !downbranchingresult->cutoff )
         {
            SCIP_CALL( executeDeepBranching(scip, lpobjval,
               &downbranchingresult->cutoff, &downbranchingresult->lperror, &scoredata->upperbounddata,
               &scoredata->ncutoffs, innerbounddata) );
         }
         if( downbranchingresult->lperror )
         {
            SCIPdebugMessage("There occurred an error while solving an lp of the upper bounded branch.\n");
            break;
         }

         SCIPdebugMessage("Going back to layer 0.\n");
         SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

         SCIPdebugMessage("First Level up branching on variable <%s>\n", SCIPvarGetName(branchvar));
         SCIP_CALL( executeBranchingOnLowerBound(scip, branchvar, branchval, upbranchingresult) );

         if( ! upbranchingresult->lperror && !upbranchingresult->cutoff )
         {
            SCIP_CALL( executeDeepBranching(scip, lpobjval,
               &upbranchingresult->cutoff, &upbranchingresult->lperror, &scoredata->lowerbounddata, &scoredata->ncutoffs,
               innerbounddata) );
         }
         if( upbranchingresult->lperror )
         {
            SCIPdebugMessage("There occurred an error while solving an lp of the lower bounded branch.\n");
            break;
         }

         SCIPdebugMessage("Going back to layer 0.\n");
         SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

         if( upbranchingresult->cutoff && downbranchingresult->cutoff )
         {
            *result = SCIP_CUTOFF;
            SCIPdebugMessage(" -> variable <%s> is infeasible in both directions\n", SCIPvarGetName(branchvar));
            break;
         }
         else if( upbranchingresult->cutoff )
         {
            SCIP_CALL( addNewUpperBound(scip, branchvar, branchval, newupperbounds, boundstatus) );
         }
         else if( downbranchingresult->cutoff )
         {
            SCIP_CALL( addNewLowerBound(scip, branchvar, branchval, newlowerbounds, boundstatus) );
         }
         else
         {
            SCIP_CALL( calculateCurrentWeight(scip, *scoredata,
               &highestscore, &highestscoreindex) );
         }
      }


      SCIPdebugMessage("End Probing Mode\n");
      SCIP_CALL( SCIPendProbing(scip) );

      if( downbranchingresult->lperror || upbranchingresult->lperror )
      {
         *result = SCIP_DIDNOTFIND;
      }
      else
      {
         SCIP_CALL( handleNewBounds(scip, boundstatus, newlowerbounds, newupperbounds, result) );
      }

      freeInnerBoundData(scip, &innerbounddata);

      SCIPfreeCleanBufferArray(scip, &boundstatus);
      SCIPfreeBufferArray(scip, &newlowerbounds);
      SCIPfreeBufferArray(scip, &newupperbounds);

      SCIPfreeBuffer(scip, &scoredata);
      SCIPfreeBuffer(scip, &upbranchingresult);
      SCIPfreeBuffer(scip, &downbranchingresult);

      if( highestscoreindex != -1 )
      {
         *bestcand = highestscoreindex;
      }

   }

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */


/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyLookahead)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   SCIP_CALL( SCIPincludeBranchruleLookahead(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeLookahead)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitLookahead)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitLookahead)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpLookahead)
{  /*lint --e{715}*/
   SCIP_VAR** tmplpcands;
   SCIP_VAR** lpcands;
   SCIP_Real* tmplpcandssol;
   SCIP_Real* lpcandssol;
   SCIP_Real* tmplpcandsfrac;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   int npriolpcands;
   int bestcand = -1;

   SCIPdebugMessage("Entering branchExeclpLookahead.\n");

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /*SCIPdebugMessage("Execlp method of lookahead branching\n");*/
   *result = SCIP_DIDNOTRUN;

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &tmplpcands, &tmplpcandssol, &tmplpcandsfrac, &nlpcands, &npriolpcands, NULL) );
   assert(nlpcands > 0);
   assert(npriolpcands > 0);

   /* copy LP banching candidates and solution values, because they will be updated w.r.t. the strong branching LP
    * solution
    */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcands, tmplpcands, nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcandssol, tmplpcandssol, nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcandsfrac, tmplpcandsfrac, nlpcands) );

   SCIPdebugMessage("The base lp has <%i> variables with fractional value.\n", nlpcands);

   SCIP_CALL( selectVarLookaheadBranching(scip, lpcands, lpcandssol, nlpcands, &bestcand, result) );

   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM && *result != SCIP_CONSADDED
      && 0 <= bestcand && bestcand < nlpcands )
   {
      SCIP_NODE* downchild = NULL;
      SCIP_NODE* upchild = NULL;
      SCIP_VAR* var;
      SCIP_Real val;

      assert(*result == SCIP_DIDNOTRUN);

      var = lpcands[bestcand];
      val = lpcandssol[bestcand];

      SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%g)\n",
         nlpcands, bestcand, SCIPvarGetName(var), val);
      SCIP_CALL( SCIPbranchVarVal(scip, var, val, &downchild, NULL, &upchild) );

      assert(downchild != NULL);
      assert(upchild != NULL);

      SCIPdebugMessage("Branched on variable <%s>\n", SCIPvarGetName(var));
      *result = SCIP_BRANCHED;
   }
   else
   {
      SCIPdebugMessage("Could not find any variable to branch\n");
   }

   SCIPfreeBufferArray(scip, &lpcandsfrac);
   SCIPfreeBufferArray(scip, &lpcandssol);
   SCIPfreeBufferArray(scip, &lpcands);

   SCIPdebugMessage("Exiting branchExeclpLookahead.\n");

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the lookahead branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleLookahead(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create lookahead branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
   /* TODO: (optional) create branching rule specific data here */

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyLookahead) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeLookahead) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitLookahead) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitLookahead) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpLookahead) );

   /* add lookahead branching rule parameters */

   return SCIP_OKAY;
}
