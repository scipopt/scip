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
#define SCIP_DEBUG
#define SCIP_STATISTICS
/**@file   branch_lookahead.c
 * @brief  lookahead branching rule
 * @author Christoph Schubert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_lookahead.h"
#include "scip/cons_setppc.h"
#include "scip/def.h"
#include "scip/pub_branch.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/type_branch.h"
#include "scip/type_lp.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"
#include "scip/type_tree.h"
#include "scip/type_var.h"

#define BRANCHRULE_NAME            "lookahead"
#define BRANCHRULE_DESC            "fullstrong branching over two levels"
#define BRANCHRULE_PRIORITY        536870911
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define DEFAULT_USEDIRECTDOMRED              TRUE
#define DEFAULT_USEIMPLIEDDOMRED             TRUE
#define DEFAULT_USEIMPLIEDBINARYCONSTRAINTS  TRUE

/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Bool             useimplieddomred;   /**< indicates whether the second level domain reduction data should be
                                              *   gathered and used */
   SCIP_Bool             usedirectdomred;    /**< indicates whether the first level domain reduction data should be
                                              *   gathered and used */
   SCIP_Bool             useimpliedbincons;  /**< indicates whether the data for the implied binary constraints should
                                              *   be gathered and used */

   SCIP_SOL*             prevbinsolution;    /**< the previous solution in the case that in the previous run only
                                              *   non-violating implied binary constraints were added */
   SCIP_VAR*             prevbinbranchvar;   /**< the previous branching decision in the case that in the previous run
                                              *   only non-violating implied binary constraints were added */
   SCIP_Real             prevbinbranchsol;   /**< the previous branching value in case that in the previous run only
                                              *   non-violating implied binary constraints were added */
#ifdef SCIP_STATISTICS
   int                   nfirstlvllps;       /**/
   int                   nsecondlvllps;      /**/
   int                   nfirstlvlcutoffs;   /**/
   int                   nsecondlvlcutoffs;  /**/
   int*                  nresults;
#endif
};

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
   WeightData*           upperbounddata;
   WeightData*           lowerbounddata;
} ScoreData;

typedef struct
{
   SCIP_Real             objval;
   SCIP_Bool             cutoff;
   SCIP_Bool             lperror;
   SCIP_Bool             nobranch;
} BranchingResultData;

/**
 * This enum is used to represent whether an upper bound, lower bound or both are set for a variable.
 */
typedef enum
{
   BOUNDSTATUS_NONE = 0,
   BOUNDSTATUS_UPPERBOUND,
   BOUNDSTATUS_LOWERBOUND,
   BOUNDSTATUS_BOTH
} BOUNDSTATUS;

/**
 * This struct collect the bounds, which can be used in the root problem. Therefore it contains two types of bounds:
 * - The bounds that occur when both up and down branches of a variable after a first level branch are cutoff. In this case
 *   the whole first level branch can be added as a restriction.
 * - The bounds that occur implicitly while branching on the second level. See SupposedBoundData for more information.
 */
typedef struct
{
   SCIP_Real*            upperbounds;        /**< The current upper bound for each active variable. Only contains meaningful
                                              *   data, if the corresponding boundstatus entry is BOUNDSTATUS_UPPERBOUND or
                                              *   BOUNDSTATUS_BOTH. */
   SCIP_Real*            lowerbounds;        /**< The current lower bound for each active variable. Only contains meaningful
                                              *   data, if the corresponding boundstatus entry is BOUNDSTATUS_LOWERBOUND or
                                              *   BOUNDSTATUS_BOTH. */
   BOUNDSTATUS*          boundstatus;        /**< The current boundstatus for each active variable. Depending on this value
                                              *   the corresponding upperbound and lowerbound values are meaningful.*/
   int*                  boundedvars;        /**< Contains the var indices that have entries in the other arrays. This array
                                              *   may be used to only iterate over the relevant variables. */
   int                   nboundedvars;       /**< The length of the boundedvars array. */
} ValidDomRedData;

/**
 * This struct collects the bounds, that are given implicitly on the second branching level.
 * Concrete: If a variable is regarded on both sides of the second level and is infeasible (in the same bound direction) on
 * both sides, the weaker bound can be applied.
 * Even more concrete: First level branching on variable x, second level branching on variable y (and may others). If the
 * constraint y <= 3 on the up branch of x and y <= 6 on the down branch of x are both infeasible, the y <= 3 bound can be
 * applied on the first level.
 */
typedef struct
{
   SCIP_Real*            upperbounds;        /**< The current upper bound for each active variable. Only contains meaningful
                                              *   data, if the corresponding boundstatus entry is BOUNDSTATUS_UPPERBOUND or
                                              *   BOUNDSTATUS_BOTH. */
   int*                  nupperboundupdates; /**< The number of times the corresponding upper bound was updated. */
   SCIP_Real*            lowerbounds;        /**< The current lower bound for each active variable. Only contains meaningful
                                              *   data, if the corresponding boundstatus entry is BOUNDSTATUS_LOWERBOUND or
                                              *   BOUNDSTATUS_BOTH. */
   int*                  nlowerboundupdates; /**< The number of times the corresponding lower bound was updated. */
   int*                  boundedvars;        /**< Contains the var indices that have entries in the other arrays. This array
                                              *   may be used to only iterate over the relevant variables. */
   int                   nboundedvars;       /**< The length of the boundedvars array. */
} SupposedDomRedData;

typedef struct
{
   SCIP_VAR**            eithervars;
   SCIP_VAR**            othervars;
   int                   nentries;
   int                   memsize;
} BinaryBoundData;

/*
 * Local methods for the data structures
 */
static
void initWeightData(
   WeightData*           weightdata
)
{
   weightdata->highestweight = 0;
   weightdata->numberofweights = 0;
   weightdata->sumofweights = 0;
}

static
void initScoreData(
   ScoreData*            scoredata,
   int                   currentbranchvar
)
{
   scoredata->ncutoffs = 0;
   scoredata->varindex = currentbranchvar;
   initWeightData(scoredata->lowerbounddata);
   initWeightData(scoredata->upperbounddata);
}

static
void initBranchingResultData(
   SCIP*                 scip,               /**< SCIP data structure */
   BranchingResultData*  resultdata
)
{
   resultdata->objval = SCIPinfinity(scip);
   resultdata->cutoff = TRUE;
   resultdata->lperror = FALSE;
   resultdata->nobranch = FALSE;
}

static
SCIP_RETCODE allocValidBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   ValidDomRedData**      validbounddata
)
{
   int ntotalvars;

   ntotalvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBuffer(scip, validbounddata) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*validbounddata)->upperbounds, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*validbounddata)->lowerbounds, ntotalvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &(*validbounddata)->boundstatus, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*validbounddata)->boundedvars, ntotalvars) );
   return SCIP_OKAY;
}

/**
 * Clears the given struct.
 * Assumptions:
 * - The boundstatus array was cleared, when the bounds were transferred to the valid bounds data structure.
 * - The upper-/lowerbounds arrays don't have to be reset, as these are only read in connection with the boundstatus array.
 * - The boundedvars array is only read in connection with the nboundedvars value, which will be set to 0.
 */
static
void initValidBoundData(
   ValidDomRedData*       validbounddata      /*< The struct that should get cleared.*/
)
{
   validbounddata->nboundedvars = 0;
}

static
void freeValidBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   ValidDomRedData**      validbounddata
)
{
   SCIPfreeBufferArray(scip, &(*validbounddata)->boundedvars);
   SCIPfreeCleanBufferArray(scip, &(*validbounddata)->boundstatus);
   SCIPfreeBufferArray(scip, &(*validbounddata)->lowerbounds);
   SCIPfreeBufferArray(scip, &(*validbounddata)->upperbounds);
   SCIPfreeBuffer(scip, validbounddata);
}

static
SCIP_RETCODE allocSupposedBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   SupposedDomRedData**   supposedbounddata
)
{
   int ntotalvars;

   ntotalvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBuffer(scip, supposedbounddata) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*supposedbounddata)->upperbounds, ntotalvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &(*supposedbounddata)->nupperboundupdates, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*supposedbounddata)->lowerbounds, ntotalvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &(*supposedbounddata)->nlowerboundupdates, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*supposedbounddata)->boundedvars, ntotalvars) );
   return SCIP_OKAY;
}

/**
 * Clears the given struct.
 * Assumptions:
 * - The boundstatus array was cleared, when the bounds were transferred to the valid bounds data structure.
 * - The upper-/lowerbounds arrays don't have to be reset, as these are only read in connection with the boundstatus array.
 * - The boundedvars array is only read in connection with the nboundedvars value, which will be set to 0.
 */
static
void initSupposedBoundData(
   SupposedDomRedData*    supposedbounddata   /*< The struct that should get cleared.*/
)
{
   supposedbounddata->nboundedvars = 0;
}

static
void freeSupposedBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   SupposedDomRedData**   supposedbounddata
)
{
   SCIPfreeBufferArray(scip, &(*supposedbounddata)->boundedvars);
   SCIPfreeCleanBufferArray(scip, &(*supposedbounddata)->nlowerboundupdates);
   SCIPfreeBufferArray(scip, &(*supposedbounddata)->lowerbounds);
   SCIPfreeCleanBufferArray(scip, &(*supposedbounddata)->nupperboundupdates);
   SCIPfreeBufferArray(scip, &(*supposedbounddata)->upperbounds);
   SCIPfreeBuffer(scip, supposedbounddata);
}

static
SCIP_RETCODE allocBinaryBoundData(
   SCIP*                 scip,
   BinaryBoundData**     binarybounddata,
   int                   nentries
)
{
   SCIP_CALL( SCIPallocBuffer(scip, binarybounddata) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*binarybounddata)->eithervars, nentries) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*binarybounddata)->othervars, nentries) );
   (*binarybounddata)->memsize = nentries;
   return SCIP_OKAY;
}

static
void initBinaryBoundData(
   BinaryBoundData*      binarybounddata
)
{
   binarybounddata->nentries = 0;
}

static
SCIP_RETCODE addBinaryBoundEntry(
   SCIP*                 scip,
   BinaryBoundData*      container,
   SCIP_VAR*             eithervar,
   SCIP_VAR*             othervar
)
{
   int emptyindex = container->nentries;

   if( emptyindex == container->memsize )
   {
      /* calculate new size, that can at least hold the old number of entries + 1 for the new entry */
      int newmemsize = SCIPcalcMemGrowSize(scip, emptyindex + 1);
      SCIP_CALL( SCIPreallocBufferArray(scip, &container->eithervars, newmemsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &container->othervars, newmemsize) );
      container->memsize = newmemsize;
   }

   container->eithervars[emptyindex] = eithervar;
   container->othervars[emptyindex] = othervar;
   container->nentries = emptyindex + 1;

   return SCIP_OKAY;
}

static
SCIP_Bool isBinaryBoundDataEmpty(
   BinaryBoundData*      container
)
{
   assert(container != NULL);

   return container->nentries == 0;
}

static
void freeBinaryBoundData(
   SCIP*                 scip,
   BinaryBoundData**     binarybounddata
)
{
   SCIPfreeBufferArray(scip, &(*binarybounddata)->othervars);
   SCIPfreeBufferArray(scip, &(*binarybounddata)->eithervars);
   SCIPfreeBuffer(scip, binarybounddata);
}

/*
 * Local methods for the logic
 */

/**
 * Executes the branching on the current probing node by adding a probing node with a new upper bound.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
static
SCIP_RETCODE executeDownBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< variable to branch on */
   SCIP_Real             branchvarsolval,    /**< Current (fractional) solution value of the variable. This value
                                              *   rounded down will be the upper bound of the new node. */
   BranchingResultData*  resultdata          /**< pointer to the result data which gets filled with the status */
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

   /* round the given value down, so that it can be used as the new upper bound */
   newupperbound = SCIPfeasFloor(scip, branchvarsolval);

   oldupperbound = SCIPvarGetUbLocal(branchvar);
   oldlowerbound = SCIPvarGetLbLocal(branchvar);
   SCIPdebugMessage("New upper bound: <%g>, old lower bound: <%g>, old upper bound: <%g>\n", newupperbound, oldlowerbound,
      oldupperbound);

   if( SCIPisFeasLT(scip, newupperbound, oldlowerbound) )
   {
      /* if lb > ub we can cutoff this node */
      resultdata->cutoff = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      if( SCIPisEQ(scip, oldupperbound, oldlowerbound) )
      {
         /*nothing has changed, so we cannot change any bound*/
         resultdata->nobranch = TRUE;
         /* TODO: do smth with this info. */
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
      assert(solstat != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

      /* for us an error occurred, if an error during the solving occurred, or the lp could not be solved but was not
       * cutoff, or if the iter or time limit was reached. */
      resultdata->lperror = resultdata->lperror || (solstat == SCIP_LPSOLSTAT_NOTSOLVED && resultdata->cutoff == FALSE)
         || (solstat == SCIP_LPSOLSTAT_ITERLIMIT) || (solstat == SCIP_LPSOLSTAT_TIMELIMIT);

      if( !resultdata->lperror )
      {
         /* if we have no error, we save the new objective value and the cutoff decision in the resultdata */
         resultdata->objval = SCIPgetLPObjval(scip);
         resultdata->cutoff = resultdata->cutoff || SCIPisGE(scip, resultdata->objval, SCIPgetCutoffbound(scip));
         assert(((solstat != SCIP_LPSOLSTAT_INFEASIBLE) && (solstat != SCIP_LPSOLSTAT_OBJLIMIT)) || resultdata->cutoff);
      }
   }

   return SCIP_OKAY;
}

/**
 * Executes the branching on the current probing node by adding a probing node with a new lower bound.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
static
SCIP_RETCODE executeUpBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< variable to branch on */
   SCIP_Real             branchvarsolval,    /**< Current (fractional) solution value of the variable. This value
                                              *   rounded up will be the lower bound of the new node. */
   BranchingResultData*  resultdata          /**< pointer to the result data which gets filled with the status */
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

   /* round the given value up, so that it can be used as the new lower bound */
   newlowerbound = SCIPfeasCeil(scip, branchvarsolval);

   oldlowerbound = SCIPvarGetLbLocal(branchvar);
   oldupperbound = SCIPvarGetUbLocal(branchvar);
   SCIPdebugMessage("New lower bound: <%g>, old lower bound: <%g>, old upper bound: <%g>\n", newlowerbound,
      oldlowerbound, oldupperbound);

   if( SCIPisFeasGT(scip, newlowerbound, oldupperbound) )
   {
      /* if lb > ub we can cutoff this node */
      resultdata->cutoff = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      if( SCIPisEQ(scip, oldupperbound, oldlowerbound) )
      {
         /*nothing has changed, so we cannot change any bound*/
         resultdata->nobranch = TRUE;
         /* TODO: do smth with this info. */
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
      assert(solstat != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

      /* for us an error occurred, if an error during the solving occurred, or the lp could not be solved but was not
       * cutoff, or if the iter or time limit was reached. */
      resultdata->lperror = resultdata->lperror || (solstat == SCIP_LPSOLSTAT_NOTSOLVED && resultdata->cutoff == FALSE)
         || (solstat == SCIP_LPSOLSTAT_ITERLIMIT) || (solstat == SCIP_LPSOLSTAT_TIMELIMIT);

      if( !resultdata->lperror )
      {
         /* if we have no error, we save the new objective value and the cutoff decision in the resultdata */
         resultdata->objval = SCIPgetLPObjval(scip);
         resultdata->cutoff = resultdata->cutoff || SCIPisGE(scip, resultdata->objval, SCIPgetCutoffbound(scip));
         assert(((solstat != SCIP_LPSOLSTAT_INFEASIBLE) && (solstat != SCIP_LPSOLSTAT_OBJLIMIT)) || resultdata->cutoff);
      }
   }

   return SCIP_OKAY;
}

/**
 * Handles the assignment of ne bounds (valid and supposed). Therefore the new bound is written directly over the old
 * bound. Analog the new bound status is written directly over the old one.
 *
 * @return TRUE, if no upper and lower bound for the given var were yet set; FALSE, else.
 */
static
SCIP_Bool addBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< The variable for which a bound should be added. */
   SCIP_Real             newbound,           /**< The value of the new bound. */
   SCIP_Bool             keepminbound,       /**< In case there is already a bound for the variable, this Indicates
                                              *   whether the min or the max value of the new and the old bound should
                                              *   be kept. */
   BOUNDSTATUS           boundtype,          /**< The type of the new bound. Must be BOUNDSTATUS_UPPERBOUND or
                                              *   BOUNDSTATUS_LOWERBOUND. */
   SCIP_Real*            oldbound,           /**< Pointer to the old bound. Depending on the oldboundstatus this may contain
                                              *   no meaningful data. Also gets the new bound set */
   BOUNDSTATUS*          oldboundstatus      /**< Pointer to the old boundstatus. Also gets the new status set*/
   )
{
   SCIP_Bool newboundadded = FALSE;

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(boundtype == BOUNDSTATUS_LOWERBOUND || boundtype == BOUNDSTATUS_UPPERBOUND);
   assert(oldbound != NULL);
   assert(oldboundstatus != NULL);

   if( *oldboundstatus == boundtype || *oldboundstatus == BOUNDSTATUS_BOTH )
   {
      /* we already have a valid bound with a fitting type set, so we can take min/max of this and the "newbound". */
      if (keepminbound)
      {
         *oldbound = MIN(newbound, *oldbound);
      }
      else
      {
         *oldbound = MAX(newbound, *oldbound);
      }
      SCIPdebugMessage("Updating an existent new bound. var <%s> type <%d> oldbound <%g> newbound <%g>.\n",
         SCIPvarGetName(branchvar), boundtype, *oldbound, newbound);
   }
   else
   {
      /* We either have no new bounds or only a bound with the other type for our var, so we can set the new bound directly. */
      SCIPdebugMessage("Adding new bound. var <%s> type <%d> newbound <%g>.\n", SCIPvarGetName(branchvar), boundtype,
         newbound);
      *oldbound = newbound;

      if( *oldboundstatus == BOUNDSTATUS_NONE )
      {
         *oldboundstatus = boundtype;
         newboundadded = TRUE;
      }
      else
      {
         *oldboundstatus = BOUNDSTATUS_BOTH;
      }
   }
   return newboundadded;
}

/**
 * Adds the given upper bound as a valid bound to the ValidDomRedData container.
 * A valid bound comes from a cutoff on the first level.
 */
static
void addValidUpperBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< the var to assign the new bound to */
   SCIP_Real             newupperbound,      /**< the new upper bound */
   ValidDomRedData*      validbounds         /**< the container to a add the bound to */
)
{
   int varindex;
   SCIP_Real* oldupperbound;
   BOUNDSTATUS* oldboundstatus;
   SCIP_Bool newboundadded;

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(validbounds != NULL);

   /* the container is indexed with the probindex */
   varindex = SCIPvarGetProbindex( branchvar );
   oldupperbound = &validbounds->upperbounds[varindex];
   oldboundstatus = &validbounds->boundstatus[varindex];

   /* Add the given bound to the container and keep the min of the new and the old bound.
    * We want to keep the min bound, as the minimum of both valid upper bounds is tighter and still valid. */
   newboundadded = addBound(scip, branchvar, newupperbound, TRUE, BOUNDSTATUS_UPPERBOUND, oldupperbound, oldboundstatus);

   if( newboundadded )
   {
      /* in case of a new entry add the varindex to the container */
      int nboundedvars = validbounds->nboundedvars;
      validbounds->boundedvars[nboundedvars] = varindex;
      validbounds->nboundedvars = nboundedvars + 1;
   }
}

/**
 * Adds the given lower bound as a valid bound to the ValidDomRedData container.
 * A valid bound comes from a cutoff on the first level.
 */
static
void addValidLowerBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< the var to assign the new bound to */
   SCIP_Real             newlowerbound,      /**< the new lower bound */
   ValidDomRedData*      validbounds         /**< the container to a add the bound to */
)
{
   int varindex;
   SCIP_Real* oldlowerbound;
   BOUNDSTATUS* oldboundstatus;
   SCIP_Bool newboundadded;

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(validbounds != NULL);

   /* the container is indexed with the probindex */
   varindex = SCIPvarGetProbindex( branchvar );
   oldlowerbound = &validbounds->lowerbounds[varindex];
   oldboundstatus = &validbounds->boundstatus[varindex];

   /* Add the given bound to the container and keep the max of the new and the old bound.
    * We want to keep the max bound, as the maximum of both valid lower bounds is tighter and still valid. */
   newboundadded = addBound(scip, branchvar, newlowerbound, FALSE, BOUNDSTATUS_LOWERBOUND, oldlowerbound, oldboundstatus);

   if( newboundadded )
   {
      /* in case of a new entry add the varindex to the container */
      int nboundedvars = validbounds->nboundedvars;
      validbounds->boundedvars[nboundedvars] = varindex;
      validbounds->nboundedvars = nboundedvars + 1;
   }
}

/**
 * Determines the status for a given number of upper and lower bound updates.
 * Used in context of the "supposed bound adding" as we only have the status indirectly through the number of updates of
 * the respective bound.
 */
static
BOUNDSTATUS getStatus(
   int                   nupperupdates,      /**< number of updates of the upper bound */
   int                   nlowerupdates       /**< number of updates of the lower bound */
)
{
   BOUNDSTATUS status;

   if( nupperupdates > 0 && nlowerupdates > 0 )
   {
      status = BOUNDSTATUS_BOTH;
   }
   else if( nupperupdates > 0 )
   {
      status = BOUNDSTATUS_UPPERBOUND;
   }
   else if( nlowerupdates > 0 )
   {
      status = BOUNDSTATUS_LOWERBOUND;
   }
   else
   {
      status = BOUNDSTATUS_NONE;
   }

   return status;
}

/**
 * Adds the given upper bound as a supposed bound to the SupposedDomRedData container.
 * A supposed upper bound is a cutoff of an up branch on the second level. We call it supposed, as both second level up
 * branches for the same first and second level variable have to be cutoff to count it as a valid upper bound.
 * Example: (first level var x, branched on 3.5; second level var y)
 * first level down: x <= 3 && y <= 6: (doesn't matter)
 *                   x <= 3 && y >= 7: cutoff
 * first level up:   x >= 4 && y <= 4: (doesn't matter)
 *                   x >= 4 && y >= 5: cutoff
 * In this case we have a new valid upper bound on the base level, namely y < 7 (= max{7,5}), as there is no way to get a
 * valid value for x if we chose y >= 7.
 */
static
void addSupposedUpperBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< the var to assign the new bound to */
   SCIP_Real             newupperbound,      /**< the (possibly) new upper bound */
   SupposedDomRedData*   supposedbounds      /**< the container to add the upper bound to */
)
{
   int varindex;
   int prevnupperupdates;
   int prevnlowerupdates;
   SCIP_Real* oldupperbound;
   SCIP_Bool newboundadded;
   BOUNDSTATUS oldboundstatus;

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(supposedbounds != NULL);

   /* the container is indexed with the probindex */
   varindex = SCIPvarGetProbindex( branchvar );
   oldupperbound = &supposedbounds->upperbounds[varindex];
   prevnupperupdates = supposedbounds->nupperboundupdates[varindex];
   prevnlowerupdates = supposedbounds->nlowerboundupdates[varindex];
   oldboundstatus = getStatus(prevnupperupdates, prevnlowerupdates);

   /* Add the given bound to the container and keep the max of the new and the old bound.
    * We want to keep the max bound, as the maximum of both second level upper bounds (where the supposed data comes from)
    * is the value we can transfer as a valid upper bound. */
   newboundadded = addBound(scip, branchvar, newupperbound, FALSE, BOUNDSTATUS_UPPERBOUND, oldupperbound, &oldboundstatus);

   /* increment the update counter */
   supposedbounds->nupperboundupdates[varindex] = prevnupperupdates + 1;
   assert(supposedbounds->nupperboundupdates[varindex] >= 1);
   assert(supposedbounds->nupperboundupdates[varindex] <= 2);

   if( newboundadded )
   {
      /* in case of a new entry add the varindex to the container */
      int nboundedvars = supposedbounds->nboundedvars;
      supposedbounds->boundedvars[nboundedvars] = varindex;
      supposedbounds->nboundedvars = nboundedvars + 1;
   }
}

/**
 * Adds the given lower bound as a supposed bound to the SupposedDomRedData container.
 * A supposed lower bound is a cutoff of a down branch on the second level. We call it supposed, as both second level down
 * branches for the same first and second level variable have to be cutoff to count it as a valid lower bound.
 * Example: (first level var x, branched on 3.5; second level var y)
 * first level down: x <= 3 && y <= 6: cutoff
 *                   x <= 3 && y >= 7: (doesn't matter)
 * first level up:   x >= 4 && y <= 4: cutoff
 *                   x >= 4 && y >= 5: (doesn'matter)
 * In this case we have a new valid lower bound on the base level, namely y > 4 (= min{4,6}), as there is no way to get a
 * valid value for x if we chose y <= 4.
 */
static
void addSupposedLowerBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< the var to assign the new bound to */
   SCIP_Real             newlowerbound,      /**< the (possibly) new lower bound */
   SupposedDomRedData*   supposedbounds      /**< the container to add the lower bound to */
)
{
   int varindex;
   int prevnupperupdates;
   int prevnlowerupdates;
   SCIP_Real* oldlowerbound;
   SCIP_Bool newboundadded;
   BOUNDSTATUS oldboundstatus;

   assert(scip != NULL);
   assert(branchvar != NULL);
   assert(supposedbounds != NULL);

   /* the container is indexed with the probindex */
   varindex = SCIPvarGetProbindex( branchvar );
   oldlowerbound = &supposedbounds->lowerbounds[varindex];
   prevnupperupdates = supposedbounds->nupperboundupdates[varindex];
   prevnlowerupdates = supposedbounds->nlowerboundupdates[varindex];
   oldboundstatus = getStatus(prevnupperupdates, prevnlowerupdates);

   /* Add the given bound to the container and keep the min of the new and the old bound.
    * We want to keep the min bound, as the minimum of both second level lower bounds (where the supposed data comes from)
    * is the value we can transfer as a valid lower bound. */
   newboundadded = addBound(scip, branchvar, newlowerbound, TRUE, BOUNDSTATUS_LOWERBOUND, oldlowerbound, &oldboundstatus);

   /* increment the update counter */
   supposedbounds->nlowerboundupdates[varindex] = prevnlowerupdates + 1;
   assert(supposedbounds->nlowerboundupdates[varindex] >= 1);
   assert(supposedbounds->nlowerboundupdates[varindex] <= 2);

   if( newboundadded )
   {
      /* in case of a new entry add the varindex to the container */
      int nboundedvars = supposedbounds->nboundedvars;
      supposedbounds->boundedvars[nboundedvars] = varindex;
      supposedbounds->nboundedvars = nboundedvars + 1;
   }
}

/**
 * Calculate the weight of a node for given up and down gains.
 *
 * @return the node weight for the given gains.
 */
static
SCIP_Real calculateNodeWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             upgain,             /**< the gain on the up branch */
   SCIP_Real             downgain            /**< the gain on the down branch */
)
{
   SCIP_Real result;
   SCIP_Real min;
   SCIP_Real max;
   /* the weights are chosen arbitrarily */
   SCIP_Real minweight = 4;
   SCIP_Real maxweight = 1;

   assert(scip != NULL);
   assert(SCIPisGE(scip, upgain, 0.0));
   assert(SCIPisGE(scip, downgain, 0.0));

   min = MIN(downgain, upgain);
   max = MAX(upgain, downgain);

   result = minweight * min + maxweight * max;

   SCIPdebugMessage("The calculated weight of <%g> and <%g> is <%g>.\n", upgain, downgain, result);

   return result;
}

/**
 * Create a name for the implied binary bounds.
 */
static
void createBinaryBoundConstraintName(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             eithervar,          /**< one variable for the constraint name */
   SCIP_VAR*             othervar,           /**< other variable for the constraint name */
   char*                 constraintname      /**< the char pointer to store the name in */

)
{
   const char* eithervarname;
   const char* othervarname;

   eithervarname = SCIPvarGetName(eithervar);
   othervarname = SCIPvarGetName(othervar);

   sprintf(constraintname, "lookahead_bin_%s_%s", eithervarname, othervarname);
}

/**
 * Checks whether the variables in the given solution violate a set packing constraint.
 * That means val[basevarforbound] + val[deepvarforbound] <= 1 is violated.
 *
 * @return TRUE, if the setpacking constraint is violated; FALSE, otherwise.
 */
static
SCIP_Bool isConstraintViolatedByBaseSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             baselpsol,          /**< the lp solution in the base node */
   SCIP_VAR*             basevarforbound,    /**< the first level branch var, ready for use in the grand child binary cons */
   SCIP_VAR*             deepvarforbound     /**< the second level branch var, ready for use in the grand child binary cons */
)
{
   SCIP_Real basesolval;
   SCIP_Real deepsolval;
   SCIP_Real sum;
   SCIP_Bool violated;

   assert(scip != NULL);
   assert(baselpsol != NULL);
   assert(basevarforbound != NULL);
   assert(deepvarforbound != NULL);
   assert(SCIPvarIsBinary(basevarforbound));
   assert(SCIPvarIsBinary(deepvarforbound));

   /* get values in the base lp solution of the given vars */
   basesolval = SCIPgetSolVal(scip, baselpsol, basevarforbound);
   deepsolval = SCIPgetSolVal(scip, baselpsol, deepvarforbound);
   sum = basesolval + deepsolval;

   /* check whether the constraint basesolval + deepsolval <= 1 is violated */
   violated = SCIPisGT(scip, sum, 1.0);
   SCIPdebugMessage("The given base lp values <%g> from <%s> (binary) and <%g> from <%s> (binary) have the sum <%g>.\n",
      basesolval, SCIPvarGetName(basevarforbound), deepsolval, SCIPvarGetName(deepvarforbound), sum);

   if( violated )
   {
      SCIPdebugMessage("Thus the implied binary constraint <%s> + <%s> <= 1 is violated by the base lp.\n",
         SCIPvarGetName(basevarforbound), SCIPvarGetName(deepvarforbound));
   }

   return violated;
}

/**
 * Adds the implied binary constraint either to the base problem or to the buffer for later use.
 * The constraint is added to the base problem, if it is violated by the baselpsol, that means, if
 * val[basevarforbound] + val[deepvarforbound] > 1.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
static
SCIP_RETCODE addGrandChildIntegerBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             baselpsol,          /**< the lp solution in the base node */
   SCIP_NODE*            basenode,           /**< the base node */
   SCIP_VAR*             basevarforbound,    /**< the first level branch var, ready for use in the grand child binary cons */
   SCIP_VAR*             deepvarforbound,    /**< the second level branch var, ready for use in the grand child binary cons */
   BinaryBoundData*      binarybounddata,    /**< pointer which gets filled with the data for implicit binary bounds */
   SCIP_Bool*            newconstadded       /**< pointer to store whether a constraint was added to the base node or not */
)
{
   assert(scip != NULL);
   assert(baselpsol != NULL);
   assert(basenode != NULL);
   assert(basevarforbound != NULL);
   assert(deepvarforbound != NULL);
   assert(binarybounddata != NULL);
   assert(newconstadded != NULL);

   if( isConstraintViolatedByBaseSolution(scip, baselpsol, basevarforbound, deepvarforbound) )
   {
      /* if the constraint (to be created) is violated by the base solution, add the constraint directly to the base node
       * and return from the branching rule. */
      SCIP_VAR* vars[2];
      SCIP_CONS* constraint;
      char constraintname[SCIP_MAXSTRLEN];
      vars[0] = basevarforbound;
      vars[1] = deepvarforbound;

      /* create a meaningful constraint name */
      createBinaryBoundConstraintName(scip, basevarforbound, deepvarforbound, constraintname);
      assert(constraintname != NULL);

      /* create setpacking constraint (x + y <= 1) */
      SCIP_CALL( SCIPcreateConsSetpack(scip, &constraint, constraintname, 2, vars, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE,
         FALSE, FALSE, FALSE, FALSE) );
      assert(constraint != NULL);

#ifdef PRINTNODECONS
      SCIPdebugMessage("Adding following constraint:\n");
      SCIP_CALL( SCIPprintCons(scip, constraint, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* add the constraint to the base node, so that it is contained in the next lp solution */
      SCIP_CALL( SCIPaddConsNode(scip, basenode, constraint, NULL) );
      /* release the constraint, as it is no longer needed */
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

      *newconstadded = TRUE;
   }
   else
   {
      /* add the constraint to the buffer and add them all later to the problem */
      SCIP_CALL( addBinaryBoundEntry(scip, binarybounddata, basevarforbound, deepvarforbound) );

      *newconstadded = FALSE;
   }

   return SCIP_OKAY;
}

/**
 * Executes the second level branching on a given variable.
 * Set the value of result to SCIP_CONSADDED, if a constraint was added to the base node.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
static
SCIP_RETCODE executeDeepBranchingOnVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< the lookahead branchruledata */
   SCIP_SOL*             baselpsol,          /**< the lp solution in the base node */
   SCIP_NODE*            basenode,           /**< the base node */
   SCIP_Real             lpobjval,           /**< objective value of the base lp */
   SCIP_VAR*             basevarforbound,    /**< the first level branch var, ready for use in the grand child binary cons.
                                              *   Can be NULL.*/
   SCIP_VAR*             deepbranchvar,      /**< variable to branch up and down on */
   SCIP_Real             deepbranchvarsolval,/**< (fractional) solution value of the branching variable */
   SCIP_Bool*            fullcutoff,         /**< resulting decision whether this branch is cutoff */
   SCIP_Bool*            lperror,            /**< Pointer that gets filled in case there occurred an error */
   WeightData*           weightdata,         /**< container to be filled with the weight relevant data */
   int*                  ncutoffs,           /**< current (input) and resulting (output) number of cutoffs */
   SupposedDomRedData*   supposedbounds,     /**< pointer which gets filled with the implicit domain reduction data from the
                                              *   second level */
   BinaryBoundData*      binarybounddata,    /**< pointer which gets filled with the data for implicit binary bounds */
   SCIP_RESULT*          result              /**< pointer to store results of branching */
)
{
   BranchingResultData* downresultdata;
   BranchingResultData* upresultdata;
   SCIP_Real currentweight;

   assert(scip != NULL);
   assert(deepbranchvar != NULL);
   assert(ncutoffs != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, &downresultdata) );
   SCIP_CALL( SCIPallocBuffer(scip, &upresultdata) );
   initBranchingResultData(scip, downresultdata);
   initBranchingResultData(scip, upresultdata);

   SCIPdebugMessage("Second level down branching on variable <%s>\n", SCIPvarGetName(deepbranchvar));
   /* execute the second level down branching */
   SCIP_CALL( executeDownBranching(scip, deepbranchvar, deepbranchvarsolval, downresultdata) );
#ifdef SCIP_STATISTICS
   branchruledata->nsecondlvllps++;
#endif

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
      /* execute the second level up branching */
      SCIP_CALL( executeUpBranching(scip, deepbranchvar, deepbranchvarsolval, upresultdata) );
#ifdef SCIP_STATISTICS
      branchruledata->nsecondlvllps++;
#endif

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
            /* if both branches are not cutoff we can calculate the weight for the current first level branching node */
            SCIP_Real downgain;
            SCIP_Real upgain;

            /* in SCIP we minimize, so the (non-negative) gain is difference between the new obj value and the base lp one */
            downgain = downresultdata->objval - lpobjval;
            upgain = upresultdata->objval - lpobjval;

            SCIPdebugMessage("The difference between the objective values of the base lp and the upper bounded lp is <%g>\n",
               downgain);
            SCIPdebugMessage("The difference between the objective values of the base lp and the lower bounded lp is <%g>\n",
               upgain);

            assert(!SCIPisFeasNegative(scip, downgain));
            assert(!SCIPisFeasNegative(scip, upgain));

            /* calculate the weight of both gains */
            currentweight = calculateNodeWeight(scip, upgain, downgain);

            /* add the new weight to the weight data */
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
#ifdef SCIP_STATISTICS
            branchruledata->nsecondlvlcutoffs = branchruledata->nsecondlvlcutoffs + 2;
#endif
         }
         else
         {
            *fullcutoff = FALSE;
            *ncutoffs = *ncutoffs + 1;
#ifdef SCIP_STATISTICS
            branchruledata->nsecondlvlcutoffs++;
#endif

            if( upresultdata->cutoff )
            {
               if( basevarforbound != NULL && branchruledata->useimpliedbincons && SCIPvarIsBinary(deepbranchvar) )
               {
                  /* if the up branching was cutoff and both branching variables are binary we can deduce a binary
                   * constraint */
                  SCIP_Bool consadded;
                  SCIP_CALL( addGrandChildIntegerBound(scip, baselpsol, basenode, basevarforbound, deepbranchvar,
                     binarybounddata, &consadded) );

                  if( consadded )
                  {
                     /* the deduced binary constraint is violated by the base lp solution and was added directly on the base
                      * node. So we set the result value to stop the branching and return to SCIP. */
                     *result = SCIP_CONSADDED;
                  }
               }

               if( branchruledata->useimplieddomred )
               {
                  /* we add the cutoff to the "supposed" buffer, so that it may be transferred to the "valid" buffer later on */
                  addSupposedUpperBound(scip, deepbranchvar, deepbranchvarsolval, supposedbounds);
               }
            }
            if( downresultdata->cutoff )
            {
               if( basevarforbound != NULL && branchruledata->useimpliedbincons && SCIPvarIsBinary(deepbranchvar) )
               {
                  /* if the down branching was cutoff and both branching variables are binary we can deduce a binary
                   * constraint */
                  SCIP_Bool consadded;
                  SCIP_VAR* deepvarforbound;

                  /* To create an implied binary bound, we need the negated var (1-y) for this up branching case.
                   * Description: The down branching condition for binary variables is y <= 0 <=> (1-y) <= 1. As this branch
                   * is cutoff we can build, together with the binary first level variable x and the branching constraint
                   * f(x) >= 1, the constraint f(x) + (1-y) <= 1.*/
                  SCIP_CALL( SCIPgetNegatedVar(scip, deepbranchvar, &deepvarforbound) );
                  SCIP_CALL( addGrandChildIntegerBound(scip, baselpsol, basenode, basevarforbound, deepvarforbound,
                     binarybounddata, &consadded) );

                  if( consadded )
                  {
                     /* the deduced binary constraint is violated by the base lp solution and was added directly on the base
                      * node. So we set the result value to stop the branching and return to SCIP. */
                     *result = SCIP_CONSADDED;
                  }
               }

               if( branchruledata->useimplieddomred )
               {
                  /* we add the cutoff to the "supposed" buffer, so that it may be transferred to the "valid" buffer later on */
                  addSupposedLowerBound(scip, deepbranchvar, deepbranchvarsolval, supposedbounds);
               }
            }
         }

      }
   }

   SCIPfreeBuffer(scip, &upresultdata);
   SCIPfreeBuffer(scip, &downresultdata);

   return SCIP_OKAY;
}

/**
 * Executes the second level branching over all lp candidates for one first level branching variable.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
static
SCIP_RETCODE executeDeepBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< the lookahead branchruledata */
   SCIP_SOL*             baselpsol,          /**< the lp solution in the base node */
   SCIP_NODE*            basenode,           /**< the base node */
   SCIP_Real             lpobjval,           /**< objective value of the base lp */
   SCIP_VAR*             basevarforbound,    /**< the first level branch var, ready for use in the grand child binary
                                              *   bounds. Can be NULL.*/
   SCIP_Bool*            fullcutoff,         /**< resulting decision whether this branch is cutoff */
   SCIP_Bool*            lperror,            /**< pointer which gets filled if an error occurrs during the lp solving */
   WeightData*           weightdata,         /**< pointer which gets filled with the relevant weight data */
   int*                  ncutoffs,           /**< pointer which gets filled with the number of second level cutoffs */
   SupposedDomRedData*   supposedbounds,     /**< pointer which gets filled with the implicit domain reduction data from the
                                              *   second level */
   BinaryBoundData*      binarybounddata,    /**< pointer which gets filled with the data for implicit binary bounds */
   SCIP_RESULT*          result              /**< pointer to store results of branching */
)
{
   SCIP_VAR**  lpcands;
   SCIP_Real*  lpcandssol;
   int         nlpcands;
   int         j;

   assert(scip != NULL);
   assert(ncutoffs != NULL);

   /* get the branching candidates for the current probing node */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );

   SCIPdebugMessage("The deeper lp has <%i> variables with fractional value.\n", nlpcands);

   for( j = 0; j < nlpcands && !*lperror && !SCIPisStopped(scip); j++ )
   {
      /* get the current variable and solution value */
      SCIP_VAR* deepbranchvar = lpcands[j];
      SCIP_Real deepbranchvarsolval = lpcandssol[j];

      SCIPdebugMessage("Start deeper branching on variable <%s> with solution value <%g>.\n",
         SCIPvarGetName(deepbranchvar), deepbranchvarsolval);

      /* execute the second level branching for a variable */
      SCIP_CALL( executeDeepBranchingOnVar(scip, branchruledata, baselpsol, basenode, lpobjval, basevarforbound, deepbranchvar, deepbranchvarsolval,
         fullcutoff, lperror, weightdata, ncutoffs, supposedbounds, binarybounddata, result) );

      if( *fullcutoff )
      {
         /* if both second level branches for a variable are cutoff, we stop the calculation, as the current first level
          * branch has to be cutoff */
         SCIPdebugMessage("The deeper lp on variable <%s> is cutoff, as both lps are cutoff.\n",
            SCIPvarGetName(deepbranchvar));
         break;
      }

      if( *result == SCIP_CONSADDED )
      {
         /* if we found a constraint, that violates the base lp we want to stop and start the branching rule again (handled
          * by SCIP) */
         SCIPdebugMessage("The deep branching is stopped, as an implied binary constraint was found, which is violated by the base LP.\n");
         break;
      }
   }

   return SCIP_OKAY;
}

/**
 * Calculates the average weight of the weight data given.
 */
static
SCIP_Real calculateAverageWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   WeightData*           weightdata          /**< calculation data for the average weight */
)
{
   SCIP_Real averageweight;

   assert(scip != NULL);
   assert(!SCIPisFeasNegative(scip, weightdata->sumofweights));
   assert(weightdata->numberofweights >= 0);

   if( weightdata->numberofweights > 0 )
   {
      averageweight = (1.0 / weightdata->numberofweights) * weightdata->sumofweights;
   }
   else
   {
      /* in case of no weights we have define the average as 0. */
      averageweight = 0;
   }
   return averageweight;
}

/**
 * Calculates the weight for the branching on a first level variable. All data needed for the weight calculation is stored
 * in the scoredata.
 */
static
void calculateCurrentWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   ScoreData*            scoredata,          /**< the information to calculate the new weight on */
   SCIP_Real*            highestweight,      /**< pointer to the current highest weight; gets updated with the new highest
                                              *   weight */
   int*                  highestweightindex  /**< pointer to the index with the highest weight; gets updated accordingly */
)
{
   SCIP_Real averageweightupperbound = 0;
   SCIP_Real averageweightlowerbound = 0;
   SCIP_Real lambda;
   SCIP_Real totalweight;

   assert(scip != NULL);
   assert(!SCIPisFeasNegative(scip, scoredata->upperbounddata->highestweight));
   assert(!SCIPisFeasNegative(scip, scoredata->lowerbounddata->highestweight));
   assert(!SCIPisFeasNegative(scip, (SCIP_Real)scoredata->ncutoffs));
   assert(highestweight != NULL);
   assert(highestweightindex != NULL);

   /* calculate the average weights of up and down branches */
   averageweightupperbound = calculateAverageWeight(scip, scoredata->upperbounddata);
   averageweightlowerbound = calculateAverageWeight(scip, scoredata->lowerbounddata);

   /* sum both averages to use it as a normalization for the number of cutoffs */
   lambda = averageweightupperbound + averageweightlowerbound;

   SCIPdebugMessage("The lambda value is <%g>.\n", lambda);
   assert(!SCIPisFeasNegative(scip, lambda));

   /* calculate the weight by adding up the max weights of up and down branches as well as the normalized number of cutoffs */
   totalweight = scoredata->lowerbounddata->highestweight + scoredata->upperbounddata->highestweight
      + lambda * scoredata->ncutoffs;

   if( SCIPisFeasGT(scip, totalweight, *highestweight) )
   {
      /* if the new weight is higher than the old one: replace it and update date index accordingly */
      *highestweight = totalweight;
      *highestweightindex = scoredata->varindex;
   }
}

/**
 * Adds the domain reductions found throughout the execution of the branching rule.
 * Domain reductions of a variable occur if:
 * - one branch on the first level is cutoff (directly or because both branches of a second level variable were cutoff)
 * - both second level branches in the same direction for the same first level variable are cutoff
 *
 * Sets the result pointer to SCIP_CUTOFF, if the reduction of at least one var led to an infeasible state.
 * Otherwise sets the result pointer to SCIP_REDUCEDDOM, if an actual domain reduction was executed.
 * Otherwise leaves the result as is.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 * @see ValidDomRedData
 * @see SupposedDomRedData
 */
static
SCIP_RETCODE addDomainReductions(
   SCIP*                 scip,               /**< SCIP data structure */
   ValidDomRedData*      validbounds,        /**< The struct containing all bounds that should be added. */
   SCIP_RESULT*          result              /**< Pointer to the result flag. Used to set the correct flags. */
)
{
   int i;
   int nboundedvars;
   SCIP_VAR** probvars;
   int nboundsadded = 0;

   assert(scip != NULL);
   assert(validbounds != NULL);
   assert(result != NULL);

   nboundedvars = validbounds->nboundedvars;
   probvars = SCIPgetVars(scip);

   /* Instead of iterating over all problem vars, we iterate only over those that have a new bound set. */
   SCIPdebugMessage("Trying to add domain reductions for <%d> variables.\n", nboundedvars);
   for( i = 0; i < nboundedvars && *result != SCIP_CUTOFF; i++ )
   {
      int probvarindex;
      BOUNDSTATUS status;
      SCIP_VAR* branchvar;

      /* The iteration index is only valid as an index for the boundedvars array. Every other array in validbounds is
       * indexed by the probindex contained in boundedvars. */
      probvarindex = validbounds->boundedvars[i];
      status = validbounds->boundstatus[probvarindex];
      branchvar = probvars[probvarindex];

      /* Handle the new lower bound (if present) */
      if( status == BOUNDSTATUS_LOWERBOUND || status == BOUNDSTATUS_BOTH )
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         SCIP_Real oldlowerbound;
         SCIP_Real proposedlowerbound;
         SCIP_Real newlowerbound;

         /* get the old and the new lower bound */
         oldlowerbound = SCIPvarGetLbLocal(branchvar);
         proposedlowerbound = validbounds->lowerbounds[probvarindex];

         /* add the new bound */
         SCIP_CALL( SCIPtightenVarLb(scip, branchvar, proposedlowerbound, FALSE, &infeasible, &tightened) );

         newlowerbound = SCIPvarGetLbLocal(branchvar);
         SCIPdebugMessage("Variable <%s>, old lower bound <%g>, proposed lower bound <%g>, new lower bound <%g>\n",
            SCIPvarGetName(branchvar), oldlowerbound, proposedlowerbound, newlowerbound);

         if( infeasible )
         {
            /* the domain reduction may result in an empty model (ub < lb) */
            *result = SCIP_CUTOFF;
            SCIPdebugMessage("The domain reduction of variable <%s> resulted in an empty model.\n",
               SCIPvarGetName(branchvar));
         }
         else if( tightened )
         {
            /* the lb is now strictly greater than before */
            *result = SCIP_REDUCEDDOM;
            SCIPdebugMessage("The lower bound of variable <%s> was successfully tightened.\n", SCIPvarGetName(branchvar));
            nboundsadded++;
         }
      }

      /* Handle the new upper bound (if present) */
      if( *result != SCIP_CUTOFF && (status == BOUNDSTATUS_UPPERBOUND || status == BOUNDSTATUS_BOTH) )
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         SCIP_Real oldupperbound;
         SCIP_Real proposedupperbound;
         SCIP_Real newupperbound;

         /* get the old and the new upper bound */
         oldupperbound = SCIPvarGetUbLocal(branchvar);
         proposedupperbound = validbounds->upperbounds[probvarindex];

         /* add the new bound */
         SCIP_CALL( SCIPtightenVarUb(scip, branchvar, proposedupperbound, FALSE, &infeasible, &tightened) );

         newupperbound = SCIPvarGetUbLocal(branchvar);
         SCIPdebugMessage("Variable <%s>, old upper bound <%g>, proposed upper bound <%g>, new upper bound <%g>\n",
            SCIPvarGetName(branchvar), oldupperbound, proposedupperbound, newupperbound);

         if( infeasible )
         {
            /* the domain reduction may result in an empty model (ub < lb) */
            *result = SCIP_CUTOFF;
            SCIPdebugMessage("The upper bound of variable <%s> could not be tightened.\n", SCIPvarGetName(branchvar));
         }
         else if( tightened )
         {
            /* the ub is now strictly smaller than before */
            *result = SCIP_REDUCEDDOM;
            SCIPdebugMessage("The upper bound of variable <%s> was successfully tightened.\n", SCIPvarGetName(branchvar));
            nboundsadded++;
         }
      }

      if( status != BOUNDSTATUS_NONE )
      {
         /* Reset the array s.t. it only contains zero values. Necessary for the CleanBufferArray usage. */
         validbounds->boundstatus[probvarindex] = BOUNDSTATUS_NONE;
      }
   }

   SCIPdebugMessage("Added <%d> real domain reductions to the problem.\n", nboundsadded);

   return SCIP_OKAY;
}
/**
 * Add the constraints found during the lookahead branching.
 * The implied binary bounds were found when two consecutive branching of binary variables were cutoff. Then these two
 * branching constraints can be combined into a single set packing constraint.
 */
static
SCIP_RETCODE handleImpliedBinaryBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            basenode,           /**< the base node to which the bounds should be added */
   BinaryBoundData*      binarybounddata     /**< the container with the bounds to be added */
)
{
   int i;
   int nentries;

   nentries = binarybounddata->nentries;

   SCIPdebugMessage("Adding <%i> implied binary bounds.\n", nentries);
   for( i = 0; i < nentries; i++ )
   {
      SCIP_VAR* eithervar;
      SCIP_VAR* othervar;
      SCIP_VAR* constraintvars[2];
      SCIP_CONS* constraint;
      char constraintname[SCIP_MAXSTRLEN];

      /* BinaryBoundData is sequentially filled, so we can access the arrays with the for variable */
      eithervar = binarybounddata->eithervars[i];
      othervar = binarybounddata->othervars[i];
      constraintvars[0] = eithervar;
      constraintvars[1] = othervar;

      /* create the constraint with a meaningful name */
      createBinaryBoundConstraintName(scip, eithervar, othervar, constraintname);
      SCIP_CALL( SCIPcreateConsSetpack(scip, &constraint, constraintname, 2, constraintvars, TRUE, TRUE, FALSE, FALSE,
         TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

#ifdef PRINTNODECONS
      SCIPdebugMessage("Adding following constraint:\n");
      SCIP_CALL( SCIPprintCons(scip, constraint, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* add the constraint to the given node */
      SCIP_CALL( SCIPaddConsNode(scip, basenode, constraint, NULL) );
      /* release the constraint, as it is no longer needed */
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );
   }
   return SCIP_OKAY;
}

/**
 * Transfers the valid bounds contained in the supposed bounds to the valid bounds.
 * Supposed bounds are found by buffering all second level cutoffs for one first level variable. If, for the same first
 * level variable, both upper or lower branchings for the same second level variable are cutoff, this second level
 * variable can already be cutoff in the base problem.
 * This method finds those "implicit" cutoffs and adds them to the struct which is later used to add the bounds to the
 * base problem.
 */
static
void transferBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   SupposedDomRedData*   supposedbounds,     /**< Bound data from the second level branches. Source for the transfer. */
   ValidDomRedData*      validbounds         /**< Bound data from the first level branches. Target for the transfer. */
)
{
   int i;
   int nofsecondlevelbounds = 0;
   SCIP_VAR** problemvars;

   /* the supposed bound data is indexed by the global prob index, so we can use the problemvars array directly */
   problemvars = SCIPgetVars(scip);

   SCIPdebugMessage("Transferring implicit bound data to the valid bound data.\n");
   SCIPdebugMessage("Number of entries <%d>\n", supposedbounds->nboundedvars);
   for(i = 0; i < supposedbounds->nboundedvars; i++ )
   {
      /* get all data from the supposedbounds */
      int boundedvarindex = supposedbounds->boundedvars[i];
      int nupperboundupdates = supposedbounds->nupperboundupdates[boundedvarindex];
      int nlowerboundupdates = supposedbounds->nlowerboundupdates[boundedvarindex];
      SCIP_VAR* boundedvar = problemvars[boundedvarindex];

      SCIPdebugMessage("Var: <%s>, nupperboundupdates: <%d>, nlowerboundupdates: <%d>\n", SCIPvarGetName(boundedvar),
         nupperboundupdates, nlowerboundupdates);

      /* add the supposed lower bounds only, if they were updated 2 times (once for each first level branch side) */
      if( nlowerboundupdates == 2 )
      {
         SCIP_Real lowerbound = supposedbounds->lowerbounds[boundedvarindex];

         SCIPdebugMessage("Adding second level lower bound for variable <%s>. Lower bound: <%g>\n",
            SCIPvarGetName(boundedvar), lowerbound);
         addValidLowerBound(scip, boundedvar, lowerbound, validbounds);
         nofsecondlevelbounds++;
      }

      /* add the supposed upper bounds only, if they were updated 2 times (once for each first level branch side) */
      if( nupperboundupdates == 2 )
      {
         SCIP_Real upperbound = supposedbounds->upperbounds[boundedvarindex];

         SCIPdebugMessage("Adding second level upper bound for variable <%s>. Upper bound: <%g>\n",
            SCIPvarGetName(boundedvar), upperbound);
         addValidUpperBound(scip, boundedvar, upperbound, validbounds);
         nofsecondlevelbounds++;
      }

      /* clean up afterwards, to reuse the same data structure for the next first level branching. */
      supposedbounds->nupperboundupdates[boundedvarindex] = 0;
      supposedbounds->nlowerboundupdates[boundedvarindex] = 0;
   }

   SCIPdebugMessage("Added <%d> bounds from the second level.\n", nofsecondlevelbounds);

}

/**
 * Selects a variable from a set of candidates by applying strong branching with a depth of 2.
 * If the branching generated additional bounds, like domain reductions from cutoffs, those are added and a suitable result
 * code is set.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
static
SCIP_RETCODE selectVarLookaheadBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< the lookahead branchruledata */
   SCIP_SOL*             baselpsol,          /**< the lp solution of the base node */
   SCIP_Bool*            depthtoosmall,      /**< pointer to store the result of the depth check */
   SCIP_VAR**            lpcands,            /**< array of fractional variables */
   SCIP_Real*            lpcandssol,         /**< array of fractional solution values */
   int                   nlpcands,           /**< number of fractional variables/solution values */
   int*                  bestcand,           /**< calculated index of the branching variable */
   SCIP_RESULT*          result              /**< pointer to store results of branching */
)
{
   assert(scip != NULL);
   assert(lpcands != NULL);
   assert(lpcandssol != NULL);
   assert(bestcand != NULL);
   assert(result != NULL);

   if( nlpcands == 1)
   {
      /* if there is only one branching variable we can directly branch there */
      *bestcand = 0;
      return SCIP_OKAY;
   }

   if( SCIPgetDepthLimit(scip) <= (SCIPgetDepth(scip) + 2) )
   {
      /* we need to branch at least 2 steps deep */
      SCIPdebugMessage("Cannot perform probing in selectVarLookaheadBranching, depth limit reached. Current:<%i>, Max:<%i>\n",
         SCIPgetDepthLimit(scip), SCIPgetDepth(scip));
      *depthtoosmall = TRUE;
      return SCIP_OKAY;
   }

   if( nlpcands > 1 )
   {
      /* declare all variables */
      BranchingResultData* downbranchingresult;
      BranchingResultData* upbranchingresult;
      ScoreData* scoredata;

      SCIP_NODE* basenode;
      SCIP_Real lpobjval;
      SCIP_Real highestscore = 0;
      int highestscoreindex = -1;
      int i;

      ValidDomRedData* validbounds = NULL;
      SupposedDomRedData* supposedbounds = NULL;
      BinaryBoundData* binarybounddata = NULL;

      /* allocate all structs */
      SCIP_CALL( SCIPallocBuffer(scip, &downbranchingresult) );
      SCIP_CALL( SCIPallocBuffer(scip, &upbranchingresult) );
      SCIP_CALL( SCIPallocBuffer(scip, &scoredata) );
      SCIP_CALL( SCIPallocBuffer(scip, &scoredata->lowerbounddata) );
      SCIP_CALL( SCIPallocBuffer(scip, &scoredata->upperbounddata) );

      if( branchruledata->usedirectdomred || branchruledata->useimplieddomred )
      {
         SCIP_CALL( allocValidBoundData(scip, &validbounds) );
         initValidBoundData(validbounds);
         if( branchruledata->useimplieddomred )
         {
            SCIP_CALL( allocSupposedBoundData(scip, &supposedbounds) );
         }
      }
      if( branchruledata->useimpliedbincons )
      {
         /* the initial number of entries in the struct is chosen arbitrarily */
         SCIP_CALL( allocBinaryBoundData(scip, &binarybounddata, (int)SCIPceil(scip, 0.1 * nlpcands)) );
         initBinaryBoundData(binarybounddata);
      }

      /* init all structs and variables */
      initBranchingResultData(scip, downbranchingresult);
      initBranchingResultData(scip, upbranchingresult);

      basenode = SCIPgetCurrentNode(scip);
      lpobjval = SCIPgetLPObjval(scip);

      SCIPdebugMessage("The objective value of the base lp is <%g>.\n", lpobjval);

      /* use the probing mode, so that we can execute our branching without influencing the original branching tree. */
      SCIP_CALL( SCIPstartProbing(scip) );
      SCIPdebugMessage("Start Probing Mode\n");

      for( i = 0; i < nlpcands && !downbranchingresult->lperror && !upbranchingresult->lperror && !SCIPisStopped(scip); i++ )
      {
         SCIP_VAR* branchvar;
         SCIP_Real branchval;

         assert(lpcands[i] != NULL);

         /* init theses structs for each var, as the contained data is read at the end of the loop */
         if( branchruledata->useimplieddomred )
         {
            initSupposedBoundData(supposedbounds);
         }
         initScoreData(scoredata, i);

         branchvar = lpcands[i];
         branchval = lpcandssol[i];

         SCIPdebugMessage("Start branching on variable <%s>\n", SCIPvarGetName(branchvar));

         SCIPdebugMessage("First level down branching on variable <%s>\n", SCIPvarGetName(branchvar));
         /* execute the down branching on first level for the variable "branchvar" */
         SCIP_CALL( executeDownBranching(scip, branchvar, branchval, downbranchingresult) );
#ifdef SCIP_STATISTICS
         branchruledata->nfirstlvllps++;
#endif

         /* if an error/cutoff occurred on the first level, we cannot branch deeper*/
         if( !downbranchingresult->lperror && !downbranchingresult->cutoff )
         {
            SCIP_VAR* basevarforbound = NULL;
            if( SCIPvarIsBinary(branchvar) )
            {
               /* To (possibly) create an implied binary bound later, we need the negated var (1-x) for this down branching
                * case.
                * Description: The down branching condition for binary variables is x <= 0. In case of a cutoff on the
                * second level after y <= 0 (or y >= 1), we can deduce a bound (1-x) + (1-y) (or just y) <= 1. */
               SCIP_CALL( SCIPgetNegatedVar(scip, branchvar, &basevarforbound) );
            }
            /* execute the branchings on the second level after the down branching on the first level */
            SCIP_CALL( executeDeepBranching(scip, branchruledata, baselpsol, basenode, lpobjval, basevarforbound,
               &downbranchingresult->cutoff, &downbranchingresult->lperror, scoredata->upperbounddata,
               &scoredata->ncutoffs, supposedbounds, binarybounddata, result) );
         }

         /* if an error occured during the down branch (and both following second level branchings), we want to stop
          * immediately */
         if( downbranchingresult->lperror )
         {
            *result = SCIP_DIDNOTFIND;
            SCIPdebugMessage("There occurred an error while solving an lp of the upper bounded branch.\n");
            break;
         }

         /* if a constraint was added to the problem we want to stop immediately (coming from implied binary bounds) */
         if( *result == SCIP_CONSADDED )
         {
            break;
         }

         /* reset the probing model */
         SCIPdebugMessage("Going back to layer 0.\n");
         SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

         SCIPdebugMessage("First Level up branching on variable <%s>\n", SCIPvarGetName(branchvar));
         /* execute the up branching on first level for the variable "branchvar" */
         SCIP_CALL( executeUpBranching(scip, branchvar, branchval, upbranchingresult) );
#ifdef SCIP_STATISTICS
         branchruledata->nfirstlvllps++;
#endif
         /* if an error/cutoff occurred on the first level, we cannot branch deeper*/
         if( !upbranchingresult->lperror && !upbranchingresult->cutoff )
         {
            SCIP_VAR* basevarforbound = NULL;
            if( SCIPvarIsBinary(branchvar) )
            {
               /* To (possibly) create an implied binary bound later, we need the var (x) for this up branching case.
                * Description: The down up condition for binary variables is x >= 1. In case of a cutoff on the second level
                * after y <= 0 (or y >= 1), we can deduce a bound x + (1-y) (or just y) <= 1. */
               basevarforbound = branchvar;
            }
            /* execute the branchings on the second level after the up branching on the first level */
            SCIP_CALL( executeDeepBranching(scip, branchruledata, baselpsol, basenode, lpobjval, basevarforbound,
               &upbranchingresult->cutoff, &upbranchingresult->lperror, scoredata->lowerbounddata, &scoredata->ncutoffs,
               supposedbounds, binarybounddata, result) );
         }
         /* if an error occured during the up branch (and both following second level branchings), we want to stop
          * immediately */
         if( upbranchingresult->lperror )
         {
            *result = SCIP_DIDNOTFIND;
            SCIPdebugMessage("There occurred an error while solving an lp of the lower bounded branch.\n");
            break;
         }

         /* if a constraint was added to the problem we want to stop immediately (coming from implied binary bounds) */
         if( *result == SCIP_CONSADDED )
         {
            break;
         }

         /* reset the probing model */
         SCIPdebugMessage("Going back to layer 0.\n");
         SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

         if( upbranchingresult->cutoff && downbranchingresult->cutoff )
         {
            /* if both first level branchings of one variable were cutoff, the whole base node can be cutoff */
            *result = SCIP_CUTOFF;
            SCIPdebugMessage(" -> variable <%s> is infeasible in both directions\n", SCIPvarGetName(branchvar));
#ifdef SCIP_STATISTICS
            branchruledata->nfirstlvlcutoffs = branchruledata->nfirstlvlcutoffs + 2;
#endif
            break;
         }
         else
         {
            /* If we reached this point, no errors occurred, not both branches were cutoff and no constraints were added
             * directly. That means we may transfer our gathered data for implied cutoffs to the concrete cutoffs. */
            if( branchruledata->useimplieddomred )
            {
               transferBoundData(scip, supposedbounds, validbounds);
            }

            if( upbranchingresult->cutoff )
            {
               /* if the up branching (on the lower bound) was cutoff, we can add this as a new upper bound for the var */
               if( branchruledata->usedirectdomred )
               {
                  addValidUpperBound(scip, branchvar, branchval, validbounds);
               }
#ifdef SCIP_STATISTICS
               branchruledata->nfirstlvlcutoffs++;
#endif
            }
            else if( downbranchingresult->cutoff )
            {
               /* if the down branching (on the upper bound) was cutoff, we can add this as a new lower bound for the var */
               if( branchruledata->usedirectdomred )
               {
                  addValidLowerBound(scip, branchvar, branchval, validbounds);
               }
#ifdef SCIP_STATISTICS
               branchruledata->nfirstlvlcutoffs++;
#endif
            }
            else
            {
               /* if neither of both branches was cutoff we can calculate the weight for the current variable */
               calculateCurrentWeight(scip, scoredata, &highestscore, &highestscoreindex);
            }
         }
      }

      SCIPdebugMessage("End Probing Mode\n");
      SCIP_CALL( SCIPendProbing(scip) );

      if( *result != SCIP_DIDNOTFIND && *result != SCIP_CUTOFF && *result != SCIP_CONSADDED
         && branchruledata->useimpliedbincons && !isBinaryBoundDataEmpty(binarybounddata) && highestscoreindex != -1 )
      {
         /* if we have no other result status set and found implied non violating binary bounds, we add those bounds and
          * save the branching variable together with its current value and the current solution. We do this, because we may
          * be called on the next iteration with the exact same (with the added bounds, but those didn't violate any rules,
          * so the solution will be the same). In this case we can save the execution time and return directly with the
          * already obtained branching decision. */
         SCIP_CALL( handleImpliedBinaryBounds(scip, basenode, binarybounddata) );
         SCIP_CALL( SCIPlinkLPSol(scip, branchruledata->prevbinsolution) );
         SCIP_CALL( SCIPunlinkSol(scip, branchruledata->prevbinsolution) );
         branchruledata->prevbinbranchvar = lpcands[highestscoreindex];
         branchruledata->prevbinbranchsol = lpcandssol[highestscoreindex];

         *result = SCIP_CONSADDED;
      }
      else if( *result != SCIP_DIDNOTFIND && *result != SCIP_CUTOFF && *result != SCIP_CONSADDED
         && (branchruledata->usedirectdomred || branchruledata->useimplieddomred) )
      {
         /* if we have no other result status set and found (potential) implied domain reductions, we add those here */
         SCIP_CALL( addDomainReductions(scip, validbounds, result) );
      }

      /* free the structs (in reverse order of allocation) */
      if( branchruledata->useimpliedbincons )
      {
         freeBinaryBoundData(scip, &binarybounddata);
      }
      if( branchruledata->usedirectdomred || branchruledata->useimplieddomred )
      {
         if( branchruledata->useimplieddomred )
         {
            freeSupposedBoundData(scip, &supposedbounds);
         }
         freeValidBoundData(scip, &validbounds);
      }

      SCIPfreeBuffer(scip, &scoredata->upperbounddata);
      SCIPfreeBuffer(scip, &scoredata->lowerbounddata);
      SCIPfreeBuffer(scip, &scoredata);
      SCIPfreeBuffer(scip, &upbranchingresult);
      SCIPfreeBuffer(scip, &downbranchingresult);

      /* save the highest score index in the result parameter */
      SCIPdebugMessage("The current best candidate is <%i>\n", *bestcand);
      if( highestscoreindex != -1 )
      {
         *bestcand = highestscoreindex;
         SCIPdebugMessage("The new best candidate is <%i>\n", *bestcand);
      }
   }

   return SCIP_OKAY;
}

/**
 * Executes the branching on a given variable with a given value.
 * If everything worked the result pointer is set to SCIP_BRANCHED.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 * */
static
SCIP_RETCODE branchOnVar(
   SCIP*                 scip                /**< SCIP data structure */,
   SCIP_VAR*             var,                /**< the variable to branch on */
   SCIP_Real             val,                /**< the value to branch on */
   SCIP_RESULT*          result              /**< the current result, will get overwritten with SCIP_BRANCHED */
)
{
   SCIP_NODE* downchild = NULL;
   SCIP_NODE* upchild = NULL;

   /*assert(*result == SCIP_DIDNOTRUN);*/

   SCIP_CALL( SCIPbranchVarVal(scip, var, val, &downchild, NULL, &upchild) );

   assert(downchild != NULL);
   assert(upchild != NULL);

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/**
 * We can use teh previous result, stored in the branchruledata, if the branchingvariable (as an indicator) is set and
 * the current lp solution is equal to the previous lp solution.
 *
 * @return \ref TRUE, if we can branch on the previous decision, \ref FALSE, else.
 */
static
SCIP_Bool isUsePreviousResult(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             currentsol,         /**< the current base lp solution */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< the branchruledata containing the prev branching decision */
)
{
   return branchruledata->prevbinbranchvar != NULL
      && SCIPareSolsEqual(scip, currentsol, branchruledata->prevbinsolution);
}

/**
 * Uses the results from the previous run saved in the branchruledata to branch.
 * This is the case, if in the previous run only non-violating constraints were added. In that case we can use the
 * branching decision we would have made then.
 * If everything worked the result pointer contains SCIP_BRANCHED.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
static
SCIP_RETCODE usePreviousResult(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< the branchruledata containing the prev branching decision */
   SCIP_Bool*            result              /**< the pointer to the branching result */
)
{
   SCIP_VAR* var;
   SCIP_Real val;

   SCIPdebugMessage("Branching based on previous solution.\n");

   /* extract the prev branching decision */
   var = branchruledata->prevbinbranchvar;
   val = branchruledata->prevbinbranchsol;

   /* execute the actual branching */
   SCIP_CALL( branchOnVar(scip, var, val, result) );

   SCIPdebugMessage("Result: Branched based on previous solution. Variable <%s>\n",
      SCIPvarGetName(branchruledata->prevbinbranchvar));

   /* reset the var pointer, as this is our indicator whether we should branch on prev data in the next call */
   branchruledata->prevbinbranchvar = NULL;

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

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPfreeMemoryArray(scip, &branchruledata->nresults);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitLookahead)
{  /*lint --e{715}*/
   int i;
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);

   for( i = 0; i < 18; i++)
   {
      branchruledata->nresults[i] = 0;
   }

   /* Create an empty solution. Gets filled in case of implied binary bounds. */
   SCIP_CALL( SCIPcreateSol(scip, &branchruledata->prevbinsolution, NULL) );
   branchruledata->prevbinbranchvar = NULL;

   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitLookahead)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   branchruledata = SCIPbranchruleGetData(branchrule);
   /* Free the solution that was used for implied binary bounds. */
   SCIP_CALL( SCIPfreeSol(scip, &branchruledata->prevbinsolution) );

   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpLookahead)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_SOL* baselpsol = NULL;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Entering branchExeclpLookahead.\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( branchruledata->useimpliedbincons )
   {
      SCIPdebugMessage("Creating an unlinked copy of the base lp solution.\n");
      /* create temporary solution */
      SCIP_CALL( SCIPcreateSol(scip, &baselpsol, NULL) );
      /* copy the current LP solution into our temporary solution */
      SCIP_CALL( SCIPlinkLPSol(scip, baselpsol) );
      /* unlink the solution, so that newly solved lps don't have any influence on our copy */
      SCIP_CALL( SCIPunlinkSol(scip, baselpsol) );
   }

   if( branchruledata->useimpliedbincons && isUsePreviousResult(scip, baselpsol, branchruledata) )
   {
      SCIP_CALL( usePreviousResult(scip, branchruledata, result) );
   }
   else
   {
      SCIP_VAR** tmplpcands;
      SCIP_VAR** lpcands;
      SCIP_Real* tmplpcandssol;
      SCIP_Real* lpcandssol;
      SCIP_Bool depthtoosmall = FALSE;
      int nlpcands;
      int bestcand = -1;

      /* get branching candidates and their solution values (integer variables with fractional value in the current LP) */
      SCIP_CALL( SCIPgetLPBranchCands(scip, &tmplpcands, &tmplpcandssol, NULL, &nlpcands, NULL, NULL) );
      assert(nlpcands > 0);

      /* copy LP banching candidates and solution values, because they will be updated w.r.t. the strong branching LP
       * solution during the second level branchings */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcands, tmplpcands, nlpcands) );
      SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcandssol, tmplpcandssol, nlpcands) );

      SCIPdebugMessage("The base lp has <%i> variables with fractional value.\n", nlpcands);

      /* execute the main logic */
      SCIP_CALL( selectVarLookaheadBranching(scip, branchruledata, baselpsol, &depthtoosmall, lpcands, lpcandssol, nlpcands,
         &bestcand, result) );

      SCIPdebugMessage("Result is <%i>\n", *result);

      if( *result != SCIP_CUTOFF /* a variable could not be branched in any direction or any of the calculated domain
                                  * reductions was infeasible */
         && *result != SCIP_REDUCEDDOM /* the domain of a variable was reduced by evaluating the calculated cutoffs */
         && *result != SCIP_CONSADDED /* implied binary constraints were already added */
         && *result != SCIP_DIDNOTFIND /* an lp error occurred on the way */
         && !depthtoosmall /* branching depth wasn't high enough */
         && (0 <= bestcand && bestcand < nlpcands) ) /* no valid candidate index could be found */
      {
         SCIP_VAR* var;
         SCIP_Real val;

         var = lpcands[bestcand];
         val = lpcandssol[bestcand];

         SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%g)\n",
            nlpcands, bestcand, SCIPvarGetName(var), val);

         SCIP_CALL( branchOnVar(scip, var, val, result) );

         SCIPdebugMessage("Result: Branched on variable <%s>\n", SCIPvarGetName(var));
      }
      else if( *result == SCIP_REDUCEDDOM)
      {
         SCIPdebugMessage("Result: Finished LookaheadBranching by reducing domains.\n");
      }
      else if( *result == SCIP_CUTOFF )
      {
         SCIPdebugMessage("Result: Finished LookaheadBranching by cutting of, as the current problem is infeasible.\n");
      }
      else if( *result == SCIP_CONSADDED )
      {
         SCIPdebugMessage("Result: Finished LookaheadBranching by adding constraints.\n");
      }
      else if( *result == SCIP_DIDNOTFIND )
      {
         SCIPdebugMessage("Result: An error occurred during the solving of one of the lps.\n");
      }
      else if( depthtoosmall )
      {
         SCIPdebugMessage("Result: The branching depth wasn't high enough for a 2 level branching.\n");
      }
      else
      {
         SCIPdebugMessage("Result: Could not find any variable to branch on.\n");
      }

      SCIPfreeBufferArray(scip, &lpcandssol);
      SCIPfreeBufferArray(scip, &lpcands);
   }

   if( branchruledata->useimpliedbincons )
   {
      SCIP_CALL( SCIPfreeSol(scip, &baselpsol) );
   }

   SCIPdebugMessage("Exiting branchExeclpLookahead.\n");

#ifdef SCIP_STATISTICS
   {
      int i;
      branchruledata->nresults[*result]++;
      for( i = 1; i < 18; i++ )
      {
         /* see type_result.h for the id <-> enum mapping */
         SCIPinfoMessage(scip, NULL, "Result with index <%i> was chosen <%i> times\n", i, branchruledata->nresults[i]);
      }
      SCIPinfoMessage(scip, NULL, "Solved <%i> lps on the first level and <%i> lps on the second level\n",
         branchruledata->nfirstlvllps, branchruledata->nsecondlvllps);
      SCIPinfoMessage(scip, NULL, "Cutoff <%i> branches on the first level and <%i> on the second level\n",
         branchruledata->nfirstlvlcutoffs, branchruledata->nsecondlvlcutoffs);
   }
#endif

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

#ifdef SCIP_STATISTICS
   branchruledata->nfirstlvlcutoffs = 0;
   branchruledata->nfirstlvllps = 0;
   branchruledata->nsecondlvlcutoffs = 0;
   branchruledata->nsecondlvllps = 0;
   /* 17 current number of possible result values and the index is 1 based, so 17 + 1 as array size with unused 0 element */
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nresults, 17 + 1) );
#endif

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
   SCIP_CALL( SCIPaddBoolParam(scip,
      "branching/lookahead/usedirectdomred",
      "should domain reductions found via cutoff on the first level be applied?",
      &branchruledata->usedirectdomred, TRUE, DEFAULT_USEDIRECTDOMRED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
      "branching/lookahead/useimplieddomred",
      "should domain reductions found via fitting cutoffs on the second level be applied?",
      &branchruledata->useimplieddomred, TRUE, DEFAULT_USEIMPLIEDDOMRED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
      "branching/lookahead/useimpliedbincons",
      "should implied binary constraints found via cutoffs on the second level be applied?",
      &branchruledata->useimpliedbincons, TRUE, DEFAULT_USEIMPLIEDBINARYCONSTRAINTS, NULL, NULL) );
   return SCIP_OKAY;
}
