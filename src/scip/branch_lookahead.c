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
/*
#define SCIP_DEBUG
#define SCIP_STATISTICS
*/
/**@file   branch_lookahead.c
 * @brief  lookahead branching rule
 * @author Christoph Schubert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_lookahead.h"
#include "scip/common_branch_lookahead.h"
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
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define DEFAULT_USEDIRECTDOMRED              TRUE
#define DEFAULT_USEIMPLIEDDOMRED             TRUE
#define DEFAULT_USEIMPLIEDBINARYCONSTRAINTS  TRUE
#define DEFAULT_ADDBINCONSROW                FALSE
#define DEFAULT_MAXNUMBERVIOLATEDCONS        10000

/*
 * Data structures
 */

/**
 * branching rule data
 */
struct SCIP_BranchruleData
{
   SCIP_Bool             useimplieddomred;   /**< indicates whether the second level domain reduction data should be
                                              *   gathered and used */
   SCIP_Bool             usedirectdomred;    /**< indicates whether the first level domain reduction data should be
                                              *   gathered and used */
   SCIP_Bool             useimpliedbincons;  /**< indicates whether the data for the implied binary constraints should
                                              *   be gathered and used */
   SCIP_Bool             addbinconsrow;      /**< Add the implied binary constraints as a row to the problem matrix */
   SCIP_SOL*             prevbinsolution;    /**< the previous solution in the case that in the previous run only
                                              *   non-violating implied binary constraints were added */
   SCIP_VAR*             prevbinbranchvar;   /**< the previous branching decision in the case that in the previous run
                                              *   only non-violating implied binary constraints were added */
   SCIP_Real             prevbinbranchsol;   /**< the previous branching value in case that in the previous run only
                                              *   non-violating implied binary constraints were added */
   int                   maxnviolatedcons;
   SCIP_Real             prevbestup;
   SCIP_Bool             prevbestupvalid;
   SCIP_Real             prevbestdown;
   SCIP_Bool             prevbestdownvalid;
   SCIP_Real             prevprovedbound;
   int                   restartindex;
#ifdef SCIP_STATISTICS
   int                   nbinconst;          /**< counter for the nubmer of "normal" binary constraints added */
   int                   nfirstlvllps;       /**< counter for the number of lps that were solved on the first level */
   int                   nsecondlvllps;      /**< counter for the number of lps that were solved on the second level */
   int                   nfirstlvlcutoffs;   /**< counter for the number of cutoffs that were made on the first level */
   int                   nsecondlvlcutoffs;  /**< counter for the number of cutoffs that were made on the second level */
   int                   nstoflvlcutoffs;    /**< counter for the number of first lvl cutoffs, that came from two second
                                              *   level cutoffs */
   int*                  nresults;           /**< Array of counters for each result state the lookahead branching finished.
                                              *   The first (0) entry is unused, as the result states are indexed 1-based
                                              *   and we use this index as our array index. */
#endif
};

typedef struct
{
   SCIP_Bool             addimpbinconst;     /**<  */
   SCIP_Bool             depthtoosmall;      /**<  */
   SCIP_Bool             lperror;            /**< Indicates whether the solving of a lp resulted in an error. */
   SCIP_Bool             firstlvlfullcutoff; /**<  */
   SCIP_Bool             domredcutoff;       /**<  */
   SCIP_Bool             domred;             /**<  */
   SCIP_Bool             propagationdomred;  /**<  */
   SCIP_Bool             limitreached;       /**<  */
   SCIP_Bool             maxnconsreached;    /**<  */
} STATUS;

/**
 * A container to hold the data needed to calculate the weight of branch after a first level node.
 */
typedef struct
{
   SCIP_Real             highestweight;      /**< the highest weight that occurred over all second level nodes */
   SCIP_Real             sumofweights;       /**< the sum of the weights of all second level nodes */
   int                   numberofweights;    /**< the number of all weights that could be calculated (in case of an
                                              *   infeasible second level node we cannot calculate a weight) */
} WEIGHTDATA;

/**
 * A container to hold the data needed to calculate the weight of a first level node.
 */
typedef struct
{
   int                   ncutoffs;           /**< counter for the the number of second level cutoffs */
   WEIGHTDATA*           upperbounddata;     /**< the WeightData of the down branch (branched on upper bound) */
   WEIGHTDATA*           lowerbounddata;     /**< the WeightData of the up branch (branched on lower bound) */
} SCOREDATA;

/**
 * A container to hold the result of a branching.
 */
typedef struct
{
   SCIP_Real             objval;             /**< The objective value of the solved lp. Only contains meaningful data, if
                                              *   cutoff == TRUE. */
   SCIP_Bool             cutoff;             /**< Indicates whether the node was infeasible and was cutoff. */
} BRANCHINGRESULTDATA;

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
} SUPPOSEDDOMREDDATA;

/**
 * A container to hold the implied binary bounds found, as these are added after the main loop of the branching rule. The
 * size of both arrays is changed dynamically, as we otherwise would need up to nvar*nvars number of entries to represent
 * all possible pairs. Therefore we hold the current max size of the arrays in "memsize" and reallocate the size if needed.
 * The arrays are consecutively filled, and therefore can be accessed 0-based up until the "nentries".
 */
typedef struct
{
   SCIP_VAR**            eithervars;         /**< An array containing one of the variables needed to create the implied
                                              *   binary constraint. */
   SCIP_VAR**            othervars;          /**< An array containing the other variable needed to create the implied binary
                                              *   constraint. */
   int                   nentries;           /**< The number of entries in both arrays. */
   int                   nviolatedentries;   /**< The number of entries that are violated by the base LP solution. */
   int                   memsize;            /**< The number of entries that currently fit in the arrays. Only for internal
                                              *   usage! */
} BINARYBOUNDDATA;


/*
 * Local methods for the data structures
 */

static
SCIP_RETCODE allocateStatus(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS**              status              /**< the status to be allocated */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, status) );

   (*status)->addimpbinconst = FALSE;
   (*status)->depthtoosmall = FALSE;
   (*status)->lperror = FALSE;
   (*status)->firstlvlfullcutoff = FALSE;
   (*status)->domred = FALSE;
   (*status)->domredcutoff = FALSE;
   (*status)->propagationdomred = FALSE;
   (*status)->limitreached = FALSE;
   (*status)->maxnconsreached = FALSE;

   return SCIP_OKAY;
}

static
void freeStatus(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS**              status              /**< the status to be freed */
   )
{
   SCIPfreeBuffer(scip, status);
}

/**
 * Initiates the WeightData struct.
 */
static
void initWeightData(
   WEIGHTDATA*           weightdata          /**< The struct to be initiated. */
   )
{
   weightdata->highestweight = 0;
   weightdata->numberofweights = 0;
   weightdata->sumofweights = 0;
}

/**
 * Initiates the ScoreData struct and the contained WeightData container.
 */
static
void initScoreData(
   SCOREDATA*            scoredata           /**< The struct to be initiated. */
   )
{
   scoredata->ncutoffs = 0;
   initWeightData(scoredata->lowerbounddata);
   initWeightData(scoredata->upperbounddata);
}

/**
 * Initiates the BranchingResultData struct.
 */
static
void initBranchingResultData(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGRESULTDATA*  resultdata          /**< The struct to be initiated. */
   )
{
   resultdata->objval = SCIPinfinity(scip);
   resultdata->cutoff = FALSE;
}

/**
 * Allocates buffer memory for the given SupposedDomRedData and the contained arrays.
 */
static
SCIP_RETCODE allocSupposedBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   SUPPOSEDDOMREDDATA**  supposedbounddata   /**< The struct to be allocated. */
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
 * Resets the given struct.
 * Assumptions:
 * - The boundstatus array was cleared, when the bounds were transferred to the valid bounds data structure.
 * - The upper-/lowerbounds arrays don't have to be reset, as these are only read in connection with the boundstatus array.
 * - The boundedvars array is only read in connection with the nboundedvars value, which will be set to 0.
 */
static
void initSupposedBoundData(
   SUPPOSEDDOMREDDATA*    supposedbounddata   /*< The struct that should get reset.*/
   )
{
   supposedbounddata->nboundedvars = 0;
}

static
void resetSupposedBoundData(
   SUPPOSEDDOMREDDATA*    supposedbounddata
   )
{
   int i;
   for( i = 0; i < supposedbounddata->nboundedvars; i++ )
   {
      int varindex = supposedbounddata->boundedvars[i];
      supposedbounddata->nupperboundupdates[varindex] = 0;
      supposedbounddata->nlowerboundupdates[varindex] = 0;
   }
   initSupposedBoundData(supposedbounddata);
}

/**
 * Frees the buffer memory of the given ValidDomRedData and the contained arrays.
 */
static
void freeSupposedBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   SUPPOSEDDOMREDDATA**  supposedbounddata   /**< The struct that should be freed. */
   )
{
   SCIPfreeBufferArray(scip, &(*supposedbounddata)->boundedvars);
   SCIPfreeCleanBufferArray(scip, &(*supposedbounddata)->nlowerboundupdates);
   SCIPfreeBufferArray(scip, &(*supposedbounddata)->lowerbounds);
   SCIPfreeCleanBufferArray(scip, &(*supposedbounddata)->nupperboundupdates);
   SCIPfreeBufferArray(scip, &(*supposedbounddata)->upperbounds);
   SCIPfreeBuffer(scip, supposedbounddata);
}

/**
 * Allocates buffer memory for the given SupposedDomRedData and the contained arrays. The "nentries" is just a starting
 * value. The size of the arrays will be reallocated if needed.
 */
static
SCIP_RETCODE allocBinaryBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYBOUNDDATA**     binarybounddata,    /**< The struct to be allocated. */
   int                   nentries            /**< The number of entries the contained arrays should support. */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, binarybounddata) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*binarybounddata)->eithervars, nentries) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*binarybounddata)->othervars, nentries) );
   (*binarybounddata)->memsize = nentries;
   (*binarybounddata)->nentries = 0;
   (*binarybounddata)->nviolatedentries = 0;
   return SCIP_OKAY;
}

static
SCIP_Bool isBinaryConstraintViolatedBySolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             eithervar,          /**< One of both vars for the constraint. */
   SCIP_VAR*             othervar,           /**< The other one of the vars for the constraint. */
   SCIP_SOL*             baselpsol           /**< the lp solution in the base node */
   )
{
   SCIP_Real eitherval = SCIPgetSolVal(scip, baselpsol, eithervar);
   SCIP_Real otherval = SCIPgetSolVal(scip, baselpsol, othervar);

   return SCIPisGT(scip, eitherval + otherval, 1);
}


/**
 * Adds the data for an implied binary constraint to the container.
 */
static
SCIP_RETCODE addBinaryBoundEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYBOUNDDATA*      container,          /**< The container for the implied binary constraints. */
   SCIP_VAR*             eithervar,          /**< One of both vars for the constraint. */
   SCIP_VAR*             othervar,           /**< The other one of the vars for the constraint. */
   SCIP_SOL*             baselpsol           /**< the lp solution in the base node */
   )
{
   int emptyindex;

   assert(scip != NULL);
   assert(container != NULL);
   assert(eithervar != NULL);
   assert(othervar != NULL);

   emptyindex = container->nentries;

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

   if( isBinaryConstraintViolatedBySolution(scip, eithervar, othervar, baselpsol) )
   {
      container->nviolatedentries++;
   }

   return SCIP_OKAY;
}

/**
 * Checks whether the given BinaryBoundData struct is empty.
 */
static
SCIP_Bool isBinaryBoundDataEmpty(
   BINARYBOUNDDATA*      container           /**< The container to be checked. */
   )
{
   assert(container != NULL);

   return container->nentries == 0;
}

/**
 * Frees the memory occupied by the BinaryBoundData and the contained arrays.
 */
static
void freeBinaryBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYBOUNDDATA**     binarybounddata     /**< The container to be freed. */
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
   BRANCHINGRESULTDATA*  resultdata,         /**< pointer to the result data which gets filled with the status */
   STATUS*               status
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
   SCIPdebugMessage("Proposed upper bound: <%g>, old bounds: [<%g>..<%g>], new bounds: [<%g>..<%g>]\n", newupperbound,
      oldlowerbound, oldupperbound, oldlowerbound, newupperbound);

   if( SCIPisFeasLT(scip, newupperbound, oldlowerbound) )
   {
      /* if lb > ub we can cutoff this node */
      resultdata->cutoff = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      if( SCIPisFeasLT(scip, newupperbound, oldupperbound) )
      {
         /* if the new upper bound is lesser than the old upper bound and also
          * greater than (or equal to) the old lower bound we set the new upper bound.
          * oldLowerBound <= newUpperBound < oldUpperBound */
         SCIP_CALL( SCIPchgVarUbProbing(scip, branchvar, newupperbound) );
      }

      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &status->lperror, &resultdata->cutoff) );
      solstat = SCIPgetLPSolstat(scip);
      assert(solstat != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

      /* for us an error occurred, if an error during the solving occurred, or the lp could not be solved but was not
       * cutoff, or if the iter or time limit was reached. */
      status->lperror = status->lperror || (solstat == SCIP_LPSOLSTAT_NOTSOLVED && resultdata->cutoff == FALSE);

      /* if we seem to have reached a {time, node, ...}-limit or the user cancelled the execution we want to stop further
       * calculations and instead return the current calculation state */
      status->limitreached = SCIPisStopped(scip);

      if( !status->limitreached && !status->lperror )
      {
         /* if we have no error, we save the new objective value and the cutoff decision in the resultdata */
         resultdata->objval = SCIPgetLPObjval(scip);
         resultdata->cutoff = resultdata->cutoff || SCIPisGE(scip, resultdata->objval, SCIPgetCutoffbound(scip));
         assert(solstat != SCIP_LPSOLSTAT_INFEASIBLE || resultdata->cutoff);
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
   BRANCHINGRESULTDATA*  resultdata,         /**< pointer to the result data which gets filled with the status */
   STATUS*               status
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
   SCIPdebugMessage("Proposed lower bound: <%g>, old bounds: [<%g>..<%g>], new bounds: [<%g>..<%g>]\n", newlowerbound,
      oldlowerbound, oldupperbound, newlowerbound, oldupperbound);

   if( SCIPisFeasGT(scip, newlowerbound, oldupperbound) )
   {
      /* if lb > ub we can cutoff this node */
      resultdata->cutoff = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      if( SCIPisFeasGT(scip, newlowerbound, oldlowerbound) )
      {
         /* if the new lower bound is greater than the old lower bound and also
          * lesser than (or equal to) the old upper bound we set the new lower bound.
          * oldLowerBound < newLowerBound <= oldUpperBound */
         SCIP_CALL( SCIPchgVarLbProbing(scip, branchvar, newlowerbound) );
      }

      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &status->lperror, &resultdata->cutoff) );
      solstat = SCIPgetLPSolstat(scip);
      assert(solstat != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

      /* for us an error occurred, if an error during the solving occurred, or the lp could not be solved but was not
       * cutoff, or if the iter or time limit was reached. */
      status->lperror = status->lperror || (solstat == SCIP_LPSOLSTAT_NOTSOLVED && resultdata->cutoff == FALSE);

      /* if we seem to have reached a {time, node, ...}-limit or the user cancelled the execution we want to stop further
       * calculations and instead return the current calculation state */
      status->limitreached = SCIPisStopped(scip);

      if( !status->limitreached && !status->lperror )
      {
         /* if we have no error, we save the new objective value and the cutoff decision in the resultdata */
         resultdata->objval = SCIPgetLPObjval(scip);
         resultdata->cutoff = resultdata->cutoff || SCIPisGE(scip, resultdata->objval, SCIPgetCutoffbound(scip));
         assert(solstat != SCIP_LPSOLSTAT_INFEASIBLE || resultdata->cutoff);
      }
   }

   return SCIP_OKAY;
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
   SUPPOSEDDOMREDDATA*   supposedbounds      /**< the container to add the upper bound to */
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
   SUPPOSEDDOMREDDATA*   supposedbounds      /**< the container to add the lower bound to */
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
 * Create a SetPacking constraint (x+y<=1) for use as an implied binary constraint.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
static
SCIP_RETCODE createImpliedBinaryConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< the lookahead branchruledata */
   SCIP_CONS**           constraint,         /**< Pointer to store the created constraint in */
   char*                 constraintname,     /**< name of the new constraint */
   SCIP_VAR**            consvars            /**< the variables that should be contained the constraint */
   )
{
   SCIP_Bool initial = branchruledata->addbinconsrow;
   SCIP_Bool separate = branchruledata->addbinconsrow;
   SCIP_Bool enforce = FALSE;
   SCIP_Bool check = FALSE;
   SCIP_Bool propagate = TRUE;
   SCIP_Bool local = TRUE;
   SCIP_Bool modifiable = FALSE;
   SCIP_Bool dynamic = FALSE;
   SCIP_Bool removable = FALSE;
   SCIP_Bool stickingatnode = FALSE;

   SCIP_CALL( SCIPcreateConsSetpack(scip, constraint, constraintname, 2, consvars, initial, separate, enforce,
         check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
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
   STATUS*               status,             /**< the status container */
   SCIP_SOL*             baselpsol,          /**< the lp solution in the base node */
   SCIP_NODE*            basenode,           /**< the base node */
   SCIP_Real             globallpobjval,     /**< objective value of the base lp */
   SCIP_Real             localpobjval,
   SCIP_VAR*             basevarforbound,    /**< the first level branch var, ready for use in the grand child binary cons.
                                              *   Can be NULL.*/
   SCIP_VAR*             deepbranchvar,      /**< variable to branch up and down on */
   SCIP_Real             deepbranchvarsolval,/**< (fractional) solution value of the branching variable */
   SCIP_Real             deepbranchvarfrac,
   SCIP_Real*            dualbound,
   SCIP_Bool*            fullcutoff,         /**< resulting decision whether this branch is cutoff */
   WEIGHTDATA*           weightdata,         /**< container to be filled with the weight relevant data */
   int*                  ncutoffs,           /**< current (input) and resulting (output) number of cutoffs */
   SUPPOSEDDOMREDDATA*   supposedbounds,     /**< pointer which gets filled with the implicit domain reduction data from the
                                              *   second level */
   BINARYBOUNDDATA*      binarybounddata     /**< pointer which gets filled with the data for implicit binary bounds */
   )
{
   BRANCHINGRESULTDATA* downresultdata;
   BRANCHINGRESULTDATA* upresultdata;
   SCIP_Real currentweight;

   assert(scip != NULL);
   assert(deepbranchvar != NULL);
   assert(ncutoffs != NULL);
   assert(dualbound != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, &downresultdata) );
   SCIP_CALL( SCIPallocBuffer(scip, &upresultdata) );
   initBranchingResultData(scip, downresultdata);
   initBranchingResultData(scip, upresultdata);

   SCIPdebugMessage("Second level down branching on variable <%s> and value <%g>\n", SCIPvarGetName(deepbranchvar),
      deepbranchvarsolval);
   /* execute the second level down branching */
   SCIP_CALL( executeDownBranching(scip, deepbranchvar, deepbranchvarsolval, downresultdata, status) );
   SCIPstatistic(
      branchruledata->nsecondlvllps++;
   )

   if( !status->limitreached && !status->lperror )
   {
      SCIPdebugMessage("Going back to layer 1.\n");
      /* go back one layer (we are currently in depth 2) */
      SCIP_CALL( SCIPbacktrackProbing(scip, 1) );

      SCIPdebugMessage("Second level up branching on variable <%s>\n", SCIPvarGetName(deepbranchvar));
      /* execute the second level up branching */
      SCIP_CALL( executeUpBranching(scip, deepbranchvar, deepbranchvarsolval, upresultdata, status) );
      SCIPstatistic(
         branchruledata->nsecondlvllps++;
      )

      if( !status->limitreached && !status->lperror )
      {
         SCIPdebugMessage("Going back to layer 1.\n");
         /* go back one layer (we are currently in depth 2) */
         SCIP_CALL( SCIPbacktrackProbing(scip, 1) );

         if( !downresultdata->cutoff && !upresultdata->cutoff )
         {
            /* if both branches are not cutoff we can calculate the weight for the current first level branching node */
            SCIP_Real globaldowngain;
            SCIP_Real globalupgain;
            SCIP_Real localdowngain;
            SCIP_Real localupgain;

            /* in SCIP we minimize, so the (non-negative) gain is difference between the new obj value and the base lp one */
            globaldowngain = downresultdata->objval - globallpobjval;
            globalupgain = upresultdata->objval - globallpobjval;
            localdowngain = downresultdata->objval - localpobjval;
            localupgain = upresultdata->objval - localpobjval;

            SCIPdebugMessage("The difference between the objective values of the base lp and the upper bounded lp is <%g>\n",
               globaldowngain);
            if( SCIPisNegative(scip, globaldowngain) )
            {
               SCIPdebugMessage("The difference is negative. To work on we overwrite it with 0.\n");
               globaldowngain = 0;
            }

            SCIPdebugMessage("The difference between the objective values of the base lp and the lower bounded lp is <%g>\n",
               globalupgain);
            if( SCIPisNegative(scip, globalupgain) )
            {
               SCIPdebugMessage("The difference is negative. To work on we overwrite it with 0.\n");
               globalupgain = 0;
            }
            assert(!SCIPisFeasNegative(scip, globaldowngain));
            assert(!SCIPisFeasNegative(scip, globalupgain));

            /* calculate the weight of both gains */
            currentweight = SCIPgetBranchScore(scip, NULL, globaldowngain, globalupgain);

            /* add the new weight to the weight data */
            weightdata->highestweight = MAX(weightdata->highestweight, currentweight);
            weightdata->sumofweights = weightdata->sumofweights + currentweight;
            weightdata->numberofweights++;

            SCIPdebugMessage("The sum of weights is <%g>.\n", weightdata->sumofweights);
            SCIPdebugMessage("The number of weights is <%i>.\n", weightdata->numberofweights);

            SCIP_CALL( SCIPupdateVarPseudocost(scip, deepbranchvar, 0.0-deepbranchvarfrac, localdowngain, 1.0) );
            SCIP_CALL( SCIPupdateVarPseudocost(scip, deepbranchvar, 1.0-deepbranchvarfrac, localupgain, 1.0) );

            *dualbound = MAX(*dualbound, MIN(downresultdata->objval, upresultdata->objval));
            *fullcutoff = FALSE;
         }
         else if( downresultdata->cutoff && upresultdata->cutoff )
         {
            *fullcutoff = TRUE;
            *ncutoffs = *ncutoffs + 2;
            SCIPstatistic(
               branchruledata->nsecondlvlcutoffs = branchruledata->nsecondlvlcutoffs + 2;
               branchruledata->nstoflvlcutoffs++;
            )
         }
         else
         {
            *fullcutoff = FALSE;
            *ncutoffs = *ncutoffs + 1;
            SCIPstatistic(
               branchruledata->nsecondlvlcutoffs++;
            )
            if( upresultdata->cutoff )
            {
               SCIP_Real localdowngain;

               localdowngain = downresultdata->objval - localpobjval;
               SCIP_CALL( SCIPupdateVarPseudocost(scip, deepbranchvar, 0.0-deepbranchvarfrac, localdowngain, 1.0) );

               *dualbound = MAX(*dualbound, downresultdata->objval);

               if( basevarforbound != NULL && branchruledata->useimpliedbincons && SCIPvarIsBinary(deepbranchvar) )
               {
                  /* if the up branching was cutoff and both branching variables are binary we can deduce a binary
                   * constraint. add the constraint to the buffer and add them all later to the problem */
                  SCIP_CALL( addBinaryBoundEntry(scip, binarybounddata, basevarforbound, deepbranchvar, baselpsol) );
                  if( binarybounddata->nviolatedentries >= branchruledata->maxnviolatedcons )
                  {
                     status->maxnconsreached = TRUE;
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
               SCIP_Real localupgain;

               localupgain = upresultdata->objval - localpobjval;
               SCIP_CALL( SCIPupdateVarPseudocost(scip, deepbranchvar, 1.0-deepbranchvarfrac, localupgain, 1.0) );

               *dualbound = MAX(*dualbound, upresultdata->objval);

               if( basevarforbound != NULL && branchruledata->useimpliedbincons && SCIPvarIsBinary(deepbranchvar) )
               {
                  /* if the down branching was cutoff and both branching variables are binary we can deduce a binary
                   * constraint */
                  SCIP_VAR* deepvarforbound;

                  /* To create an implied binary bound, we need the negated var (1-y) for this up branching case.
                   * Description: The down branching condition for binary variables is y <= 0 <=> (1-y) <= 1. As this branch
                   * is cutoff we can build, together with the binary first level variable x and the branching constraint
                   * f(x) >= 1, the constraint f(x) + (1-y) <= 1.*/
                  SCIP_CALL( SCIPgetNegatedVar(scip, deepbranchvar, &deepvarforbound) );
                  SCIP_CALL( addBinaryBoundEntry(scip, binarybounddata, basevarforbound, deepvarforbound, baselpsol) );
                  if( binarybounddata->nviolatedentries >= branchruledata->maxnviolatedcons )
                  {
                     status->maxnconsreached = TRUE;
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

static
SCIP_Bool isExecuteDeepBranchingLoop(
   STATUS*               status,
   SCIP_Bool             fullcutoff
   )
{
   SCIP_Bool result = !status->lperror && !fullcutoff && !status->maxnconsreached;;
   if( status->lperror )
   {
      /* an error occurred during one of the second level lps */
      SCIPdebugMessage("The deep branching is stopped, as an error occurred during one of the second level lps.\n");
   }
   if( fullcutoff )
   {
      /* if both second level branches for a variable are cutoff, we stop the calculation, as the current first level
       * branch has to be cutoff */
      SCIPdebugMessage("A deeper lp is cutoff, as both lps are cutoff.\n");
   }
   if( !status->maxnconsreached )
   {
      /*
       * if the maximum number of violated constraints to be added is reached we want to stop the current loop, add the
       * constraints and restart everything.
       */
      SCIPdebugMessage("The max number of violated constraints to be added is reached.\n");
   }
   return result;
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
   STATUS*               status,             /**< the status container */
   SCIP_SOL*             baselpsol,          /**< the lp solution in the base node */
   SCIP_NODE*            basenode,           /**< the base node */
   SCIP_Real             lpobjval,           /**< objective value of the base lp */
   SCIP_Real             locallpobjval,
   SCIP_VAR*             basevarforbound,    /**< the first level branch var, ready for use in the grand child binary
                                              *   bounds. Can be NULL.*/
   SCIP_Real*            dualbound,
   SCIP_Bool*            fullcutoff,         /**< resulting decision whether this branch is cutoff */
   WEIGHTDATA*           weightdata,         /**< pointer which gets filled with the relevant weight data */
   int*                  ncutoffs,           /**< pointer which gets filled with the number of second level cutoffs */
   SUPPOSEDDOMREDDATA*   supposedbounds,     /**< pointer which gets filled with the implicit domain reduction data from the
                                              *   second level */
   BINARYBOUNDDATA*      binarybounddata     /**< pointer which gets filled with the data for implicit binary bounds */
   )
{
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   int j;

   assert(scip != NULL);
   assert(ncutoffs != NULL);

   /* get the branching candidates for the current probing node */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL, NULL) );

   SCIPdebugMessage("The deeper lp has <%i> variables with fractional value.\n", nlpcands);

   for( j = 0; j < nlpcands && isExecuteDeepBranchingLoop(status, *fullcutoff) && !SCIPisStopped(scip); j++ )
   {
      /* get the current variable and solution value */
      SCIP_VAR* deepbranchvar = lpcands[j];
      SCIP_Real deepbranchvarsolval = lpcandssol[j];
      SCIP_Real deepbranchvarfrac = lpcandsfrac[j];

      SCIPdebugMessage("Start deeper branching on variable <%s> with solution value <%g> in [<%g>..<%g>].\n",
         SCIPvarGetName(deepbranchvar), deepbranchvarsolval, SCIPvarGetLbLocal(deepbranchvar),
         SCIPvarGetUbLocal(deepbranchvar));

      /* execute the second level branching for a variable */
      SCIP_CALL( executeDeepBranchingOnVar(scip, branchruledata, status, baselpsol, basenode, lpobjval, locallpobjval,
         basevarforbound, deepbranchvar, deepbranchvarsolval, deepbranchvarfrac, dualbound, fullcutoff, weightdata,
         ncutoffs, supposedbounds, binarybounddata) );
   }

   return SCIP_OKAY;
}

/**
 * Calculates the average weight of the weight data given.
 */
static
SCIP_Real calculateAverageWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   WEIGHTDATA*           weightdata          /**< calculation data for the average weight */
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
   SCOREDATA*            scoredata,          /**< the information to calculate the new weight on */
   SCIP_Real*            currentweight       /**< pointer to the current highest weight; gets updated with the new highest
                                              *   weight */
   )
{
   SCIP_Real averageweightupperbound = 0;
   SCIP_Real averageweightlowerbound = 0;
   SCIP_Real lambda;

   assert(scip != NULL);
   assert(!SCIPisFeasNegative(scip, scoredata->upperbounddata->highestweight));
   assert(!SCIPisFeasNegative(scip, scoredata->lowerbounddata->highestweight));
   assert(!SCIPisFeasNegative(scip, (SCIP_Real)scoredata->ncutoffs));
   assert(currentweight != NULL);

   /* calculate the average weights of up and down branches */
   averageweightupperbound = calculateAverageWeight(scip, scoredata->upperbounddata);
   averageweightlowerbound = calculateAverageWeight(scip, scoredata->lowerbounddata);

   /* sum both averages to use it as a normalization for the number of cutoffs */
   lambda = averageweightupperbound + averageweightlowerbound;

   SCIPdebugMessage("The lambda value is <%g>.\n", lambda);
   assert(!SCIPisFeasNegative(scip, lambda));

   /* calculate the weight by adding up the max weights of up and down branches as well as the normalized number of cutoffs */
   *currentweight = scoredata->lowerbounddata->highestweight + scoredata->upperbounddata->highestweight
      + lambda * scoredata->ncutoffs;
}

/**
 * Add the constraints found during the lookahead branching.
 * The implied binary bounds were found when two consecutive branching of binary variables were cutoff. Then these two
 * branching constraints can be combined into a single set packing constraint.
 */
static
SCIP_RETCODE handleImpliedBinaryBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< the lookahead branchruledata */
   SCIP_NODE*            basenode,           /**< the base node to which the bounds should be added */
   BINARYBOUNDDATA*      binarybounddata     /**< the container with the bounds to be added */
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
      SCIP_CALL( createImpliedBinaryConstraint(scip, branchruledata, &constraint, constraintname, constraintvars) );

#ifdef PRINTNODECONS
      SCIPdebugMessage("Adding following constraint:\n");
      SCIP_CALL( SCIPprintCons(scip, constraint, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* add the constraint to the given node */
      SCIP_CALL( SCIPaddConsNode(scip, basenode, constraint, NULL) );
      /* release the constraint, as it is no longer needed */
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

      SCIPstatistic(
         branchruledata->nbinconst++;
      )
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
   SUPPOSEDDOMREDDATA*   supposedbounds,     /**< Bound data from the second level branches. Source for the transfer. */
   VALIDDOMREDDATA*      validbounds         /**< Bound data from the first level branches. Target for the transfer. */
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

static
SCIP_Bool isExecuteFirstLevelBranching(
   STATUS*               status
   )
{
   return !status->lperror && !status->firstlvlfullcutoff && !status->limitreached && !status->maxnconsreached;
}

static
SCIP_Bool areBoundsChanged(
   SCIP*                 scip,
   SCIP_VAR*             var,
   SCIP_Real             lowerbound,
   SCIP_Real             upperbound
   )
{
   return SCIPvarGetLbLocal(var) != lowerbound || SCIPvarGetUbLocal(var) != upperbound;
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
   SCIP_VAR**            lpcands,            /**< array of fractional variables */
   SCIP_Real*            lpcandssol,         /**< array of fractional solution values */
   SCIP_Real*            lpcandsfrac,
   int                   nlpcands,           /**< number of fractional variables/solution values */
   int*                  start,
   STATUS*               status,             /**< a container to store the algo status in */
   int*                  bestcand,           /**< calculated index of the branching variable */
   SCIP_Real*            bestdown,
   SCIP_Bool*            bestdownvalid,
   SCIP_Real*            bestup,
   SCIP_Bool*            bestupvalid,
   SCIP_Real*            provedbound
   )
{
   assert(scip != NULL);
   assert(lpcands != NULL);
   assert(lpcandssol != NULL);
   assert(status != NULL);
   assert(bestcand != NULL);

   if( nlpcands == 1)
   {
      /* if there is only one branching variable we can directly branch there */
      *bestcand = 0;
   }
   else if( SCIPgetDepthLimit(scip) <= (SCIPgetDepth(scip) + 2) )
   {
      /* we need to branch at least 2 steps deep */
      SCIPdebugMessage("Cannot perform probing in selectVarLookaheadBranching, depth limit reached. Current:<%i>, Max:<%i>\n",
         SCIPgetDepthLimit(scip), SCIPgetDepth(scip));
      status->depthtoosmall = TRUE;
   }
   else if( nlpcands > 1 )
   {
      /* declare all variables */
      BRANCHINGRESULTDATA* downbranchingresult;
      BRANCHINGRESULTDATA* upbranchingresult;
      SCOREDATA* scoredata;

      SCIP_NODE* basenode;
      SCIP_Real lpobjval;
      SCIP_Real highestscore = 0;
      SCIP_Real highestscoreupperbound = SCIPinfinity(scip);
      SCIP_Real highestscorelowerbound = -SCIPinfinity(scip);
      int highestscoreindex = -1;
      int i;
      int c;

      VALIDDOMREDDATA* validbounds = NULL;
      SUPPOSEDDOMREDDATA* supposedbounds = NULL;
      BINARYBOUNDDATA* binarybounddata = NULL;

      /* allocate all structs */
      SCIP_CALL( SCIPallocBuffer(scip, &downbranchingresult) );
      SCIP_CALL( SCIPallocBuffer(scip, &upbranchingresult) );
      SCIP_CALL( SCIPallocBuffer(scip, &scoredata) );
      SCIP_CALL( SCIPallocBuffer(scip, &scoredata->lowerbounddata) );
      SCIP_CALL( SCIPallocBuffer(scip, &scoredata->upperbounddata) );

      if( branchruledata->usedirectdomred || branchruledata->useimplieddomred )
      {
         SCIP_CALL( allocValidBoundData(scip, &validbounds) );
         if( branchruledata->useimplieddomred )
         {
            SCIP_CALL( allocSupposedBoundData(scip, &supposedbounds) );
         }
      }
      if( branchruledata->useimpliedbincons )
      {
         /* the initial number of entries in the struct is chosen arbitrarily */
         SCIP_CALL( allocBinaryBoundData(scip, &binarybounddata, (int)SCIPceil(scip, 0.1 * nlpcands)) );
      }

      /* init all structs and variables */
      initBranchingResultData(scip, downbranchingresult);
      initBranchingResultData(scip, upbranchingresult);

      basenode = SCIPgetCurrentNode(scip);
      lpobjval = SCIPgetLPObjval(scip);
      *provedbound = lpobjval;

      SCIPdebugMessage("The objective value of the base lp is <%g>.\n", lpobjval);

      /* use the probing mode, so that we can execute our branching without influencing the original branching tree. */
      SCIP_CALL( SCIPstartProbing(scip) );
      SCIPenableVarHistory(scip);
      SCIPdebugMessage("Start Probing Mode\n");

      for( i = 0, c = *start; i < nlpcands && isExecuteFirstLevelBranching(status) && !SCIPisStopped(scip); i++, c++ )
      {
         SCIP_VAR* branchvar;
         SCIP_Real branchval;
         SCIP_Real branchfrac;
         SCIP_Real downdualbound;
         SCIP_Real updualbound;

         c = c % nlpcands;
         assert(lpcands[c] != NULL);

         /* init theses structs for each var, as the contained data is read at the end of the loop */
         if( branchruledata->useimplieddomred )
         {
            initSupposedBoundData(supposedbounds);
         }
         initScoreData(scoredata);

         branchvar = lpcands[c];
         branchval = lpcandssol[c];
         branchfrac = lpcandsfrac[c];

         SCIPdebugMessage("Start branching on variable <%s>\n", SCIPvarGetName(branchvar));

         if( isExecuteFirstLevelBranching(status) )
         {
            SCIPdebugMessage("First level down branching on variable <%s>\n", SCIPvarGetName(branchvar));
            /* execute the down branching on first level for the variable "branchvar" */
            SCIP_CALL( executeDownBranching(scip, branchvar, branchval, downbranchingresult, status) );
            SCIPstatistic(
               branchruledata->nfirstlvllps++;
            )

            /* if an error/cutoff occurred on the first level, we cannot branch deeper*/
            if( !status->limitreached && !status->lperror && !downbranchingresult->cutoff )
            {
               SCIP_VAR* basevarforbound = NULL;
               SCIP_Real locallpobjval;
               SCIP_Real gain;

               if( SCIPvarIsBinary(branchvar) )
               {
                  /* To (possibly) create an implied binary bound later, we need the negated var (1-x) for this down branching
                   * case.
                   * Description: The down branching condition for binary variables is x <= 0. In case of a cutoff on the
                   * second level after y <= 0 (or y >= 1), we can deduce a bound (1-x) + (1-y) (or just y) <= 1. */
                  SCIP_CALL( SCIPgetNegatedVar(scip, branchvar, &basevarforbound) );
               }

               downdualbound = downbranchingresult->objval;
               locallpobjval = downbranchingresult->objval;

               gain = MAX(0, locallpobjval - lpobjval);
               SCIP_CALL( SCIPupdateVarPseudocost(scip, branchvar, 0-branchfrac, gain, 1.0) );

               /* execute the branchings on the second level after the down branching on the first level */
               SCIP_CALL( executeDeepBranching(scip, branchruledata, status, baselpsol, basenode, lpobjval, locallpobjval,
                     basevarforbound, &downdualbound, &downbranchingresult->cutoff, scoredata->upperbounddata,
                     &scoredata->ncutoffs, supposedbounds, binarybounddata) );
            }

            /* reset the probing model */
            SCIPdebugMessage("Going back to layer 0.\n");
            SCIP_CALL( SCIPbacktrackProbing(scip, 0) );
         }

         if( isExecuteFirstLevelBranching(status) )
         {
            SCIPdebugMessage("First Level up branching on variable <%s>\n", SCIPvarGetName(branchvar));
            /* execute the up branching on first level for the variable "branchvar" */
            SCIP_CALL( executeUpBranching(scip, branchvar, branchval, upbranchingresult, status) );
            SCIPstatistic(
               branchruledata->nfirstlvllps++;
            )
            /* if an error/cutoff occurred on the first level, we cannot branch deeper*/
            if( !status->limitreached && !status->lperror && !upbranchingresult->cutoff )
            {
               SCIP_VAR* basevarforbound = NULL;
               SCIP_Real locallpobjval;
               SCIP_Real gain;

               if( SCIPvarIsBinary(branchvar) )
               {
                  /* To (possibly) create an implied binary bound later, we need the var (x) for this up branching case.
                   * Description: The down up condition for binary variables is x >= 1. In case of a cutoff on the second level
                   * after y <= 0 (or y >= 1), we can deduce a bound x + (1-y) (or just y) <= 1. */
                  basevarforbound = branchvar;
               }

               updualbound = upbranchingresult->objval;
               locallpobjval = upbranchingresult->objval;

               gain = MAX(0, locallpobjval - lpobjval);
               SCIP_CALL( SCIPupdateVarPseudocost(scip, branchvar, 1-branchfrac, gain, 1.0) );

               /* execute the branchings on the second level after the up branching on the first level */
               SCIP_CALL( executeDeepBranching(scip, branchruledata, status, baselpsol, basenode, lpobjval, locallpobjval,
                     basevarforbound, &updualbound, &upbranchingresult->cutoff, scoredata->lowerbounddata,
                     &scoredata->ncutoffs, supposedbounds, binarybounddata) );
            }

            /* reset the probing model */
            SCIPdebugMessage("Going back to layer 0.\n");
            SCIP_CALL( SCIPbacktrackProbing(scip, 0) );
         }

         if( upbranchingresult->cutoff && downbranchingresult->cutoff )
         {
            /* if both first level branchings of one variable were cutoff, the whole base node can be cutoff */
            status->firstlvlfullcutoff = TRUE;
            SCIPdebugMessage(" -> variable <%s> is infeasible in both directions\n", SCIPvarGetName(branchvar));
            SCIPstatistic(
               branchruledata->nfirstlvlcutoffs = branchruledata->nfirstlvlcutoffs + 2;
            )
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
               SCIPstatistic(
                  branchruledata->nfirstlvlcutoffs++;
               )

               *provedbound = MAX(*provedbound, downdualbound);
            }
            else if( downbranchingresult->cutoff )
            {
               /* if the down branching (on the upper bound) was cutoff, we can add this as a new lower bound for the var */
               if( branchruledata->usedirectdomred )
               {
                  addValidLowerBound(scip, branchvar, branchval, validbounds);
               }
               SCIPstatistic(
                  branchruledata->nfirstlvlcutoffs++;
               )

               *provedbound = MAX(*provedbound, updualbound);
            }
            else if( !status->limitreached && !status->maxnconsreached )
            {
               /* if neither of both branches was cutoff we can calculate the weight for the current variable */
               SCIP_Real currentScore;

               calculateCurrentWeight(scip, scoredata, &currentScore);

               if( SCIPisFeasGT(scip, currentScore, highestscore) )
               {
                  /* if the new weight is higher than the old one: replace it and update the index accordingly */
                  highestscore = currentScore;
                  highestscoreindex = c;
                  highestscorelowerbound = SCIPvarGetLbLocal(branchvar);
                  highestscoreupperbound = SCIPvarGetUbLocal(branchvar);
                  *bestdown = downdualbound;
                  *bestdownvalid = TRUE;
                  *bestup = updualbound;
                  *bestupvalid = TRUE;
               }

               *provedbound = MAX(*provedbound, MIN(updualbound, downdualbound));
            }
         }

         if( highestscoreindex != -1 && areBoundsChanged(scip, lpcands[highestscoreindex], highestscorelowerbound,
            highestscoreupperbound) )
         {
            /* in case the bounds of the current highest scored solution have changed due to domain propagation during the
             * lookahead branching we can/should not branch on this variable but instead return the SCIP_REDUCEDDOM result */
            status->propagationdomred = TRUE;
         }
      }

      SCIPdebugMessage("End Probing Mode\n");
      SCIP_CALL( SCIPendProbing(scip) );

      *start = c;

      if( !status->lperror && !status->depthtoosmall && !status->firstlvlfullcutoff
         && branchruledata->useimpliedbincons && !isBinaryBoundDataEmpty(binarybounddata) && highestscoreindex != -1 )
      {
         /* if we have no other result status set and found implied non violating binary bounds, we add those bounds and
          * save the branching variable together with its current value and the current solution. We do this, because we may
          * be called on the next iteration with the exact same (with the added bounds, but those didn't violate any rules,
          * so the solution will be the same). In this case we can save the execution time and return directly with the
          * already obtained branching decision. */
         SCIP_CALL( handleImpliedBinaryBounds(scip, branchruledata, basenode, binarybounddata) );
         status->addimpbinconst = TRUE;

         if( !status->maxnconsreached )
         {
            SCIP_CALL( SCIPlinkLPSol(scip, branchruledata->prevbinsolution) );
            SCIP_CALL( SCIPunlinkSol(scip, branchruledata->prevbinsolution) );
            branchruledata->prevbinbranchvar = lpcands[highestscoreindex];
            branchruledata->prevbinbranchsol = lpcandssol[highestscoreindex];
            branchruledata->prevbestup = *bestup;
            branchruledata->prevbestupvalid = *bestupvalid;
            branchruledata->prevbestdown = *bestdown;
            branchruledata->prevbestdownvalid = *bestdownvalid;
            branchruledata->prevprovedbound = *provedbound;
         }
      }

      if( !status->lperror && !status->depthtoosmall && !status->firstlvlfullcutoff
         && (branchruledata->usedirectdomred || branchruledata->useimplieddomred) )
      {
         /* if we have no other result status set and found (potential) implied domain reductions, we add those here */
         SCIP_CALL( addDomainReductions(scip, validbounds, &status->domredcutoff, &status->domred) );
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
            resetSupposedBoundData(supposedbounds);
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
   SCIPdebugMessage("Branching based on previous solution.\n");

   /* execute the actual branching */
   SCIP_CALL( branchOnVar(scip, branchruledata->prevbinbranchvar, branchruledata->prevbinbranchsol,
      branchruledata->prevbestdown,branchruledata->prevbestdownvalid, branchruledata->prevbestup,
      branchruledata->prevbestupvalid, branchruledata->prevprovedbound) );
   *result = SCIP_BRANCHED;

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

   SCIPstatistic(
         int i;

         for( i = 1; i < 18; i++ )
         {
            /* see type_result.h for the id <-> enum mapping */
            SCIPinfoMessage(scip, NULL, "Result <%s> was chosen <%i> times\n", getStatusString(i), branchruledata->nresults[i]);
         }
         SCIPinfoMessage(scip, NULL, "Solved <%i> lps on the first level and <%i> lps on the second level\n",
            branchruledata->nfirstlvllps, branchruledata->nsecondlvllps);
         SCIPinfoMessage(scip, NULL, "Cutoff <%i> branches on the first level and <%i> on the second level\n",
            branchruledata->nfirstlvlcutoffs, branchruledata->nsecondlvlcutoffs);
         SCIPinfoMessage(scip, NULL, "Cutoff <%i> branches on the first level based on a full cutoff on the second level\n",
            branchruledata->nstoflvlcutoffs);
         SCIPinfoMessage(scip, NULL, "Added <%i> binary constraints\n", branchruledata->nbinconst);
         SCIPfreeMemoryArray(scip, &branchruledata->nresults);
   )
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitLookahead)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);

   SCIPstatistic(
      {
         int i;
         for( i = 0; i < 18; i++)
         {
            branchruledata->nresults[i] = 0;
         }
      }
   )

   /* Create an empty solution. Gets filled in case of implied binary bounds. */
   SCIP_CALL( SCIPcreateSol(scip, &branchruledata->prevbinsolution, NULL) );
   branchruledata->prevbinbranchvar = NULL;
   branchruledata->restartindex = 0;

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

static
SCIP_RETCODE copyCurrentSolution(
   SCIP*                 scip,
   SCIP_SOL**            lpsol
   )
{
   /* create temporary solution */
   SCIP_CALL( SCIPcreateSol(scip, lpsol, NULL) );
   /* copy the current LP solution into our temporary solution */
   SCIP_CALL( SCIPlinkLPSol(scip, *lpsol) );
   /* unlink the solution, so that newly solved lps don't have any influence on our copy */
   SCIP_CALL( SCIPunlinkSol(scip, *lpsol) );

   return SCIP_OKAY;

}

/** get a copy of the fractional candidates we can branch on
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
static
SCIP_RETCODE getLPBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           lpcands,            /**< a pointer to store the variables */
   SCIP_Real**           lpcandssol,         /**< a pointer to store the solution values of the vars */
   SCIP_Real**           lpcandsfrac,
   int*                  nlpcands            /**< a pointer to store the number of candidates */
   )
{
   SCIP_VAR** tmplpcands;
   SCIP_Real* tmplpcandssol;
   SCIP_Real* tmplpcandsfrac;

   /* get branching candidates and their solution values (integer variables with fractional value in the current LP) */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &tmplpcands, &tmplpcandssol, &tmplpcandsfrac, nlpcands, NULL, NULL) );
   assert(*nlpcands > 0);

   /* copy LP banching candidates and solution values, because they will be updated w.r.t. the strong branching LP
    * solution during the second level branchings */
   SCIP_CALL( SCIPduplicateBufferArray(scip, lpcands, tmplpcands, *nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, lpcandssol, tmplpcandssol, *nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, lpcandsfrac, tmplpcandsfrac, *nlpcands) );

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpLookahead)
{  /*lint --e{715}*/
   /* TODO CS: handle the allowaddcons flag! */
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
      /* create a copy of the current lp solution to compare it with a previously  */
      SCIP_CALL( copyCurrentSolution(scip, &baselpsol) );
      SCIPdebugMessage("Created an unlinked copy of the base lp solution.\n");
   }

   if( branchruledata->useimpliedbincons && isUsePreviousResult(scip, baselpsol, branchruledata) )
   {
      /* in case we stopped the previous run without a branching decisions we have stored the decision and execute it
       * now */
      SCIP_CALL( usePreviousResult(scip, branchruledata, result) );
   }
   else
   {
      SCIP_VAR** lpcands;
      SCIP_Real* lpcandssol;
      SCIP_Real* lpcandsfrac;
      int nlpcands;
      int bestcand = 0;
      /* TODO: maybe add a struct for the resulting branching data below? */
      SCIP_Real bestdown;
      SCIP_Bool bestdownvalid = FALSE;
      SCIP_Real bestup;
      SCIP_Bool bestupvalid = FALSE;
      SCIP_Real provedbound;
      STATUS* status;

      /* get all fractional candidates we can branch on */
      SCIP_CALL( getLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands) );

      SCIPdebugMessage("The base lp has <%i> variables with fractional value.\n", nlpcands);

      /* creating a struct to store the algorithm status */
      SCIP_CALL( allocateStatus(scip, &status) );

      /* execute the main logic */
      SCIP_CALL( selectVarLookaheadBranching(scip, branchruledata, baselpsol, lpcands, lpcandssol, lpcandsfrac, nlpcands,
            &branchruledata->restartindex, status, &bestcand, &bestdown, &bestdownvalid, &bestup, &bestupvalid, &provedbound) );
      SCIPdebugMessage("Result is <%i>\n", *result);

      if( status->firstlvlfullcutoff || status->domredcutoff )
      {
         *result = SCIP_CUTOFF;
      }
      else if( status->addimpbinconst )
      {
         *result = SCIP_CONSADDED;
      }
      else if( status->domred || status->propagationdomred )
      {
         *result = SCIP_REDUCEDDOM;
      }
      else if( status->lperror )
      {
         *result = SCIP_DIDNOTFIND;
      }

      if( *result != SCIP_CUTOFF /* a variable could not be branched in any direction or any of the calculated domain
                                  * reductions was infeasible */
         && *result != SCIP_REDUCEDDOM /* the domain of a variable was reduced by evaluating the calculated cutoffs */
         && *result != SCIP_CONSADDED /* implied binary constraints were already added */
         && *result != SCIP_DIDNOTFIND /* an lp error occurred on the way */
         && !status->depthtoosmall /* branching depth wasn't high enough */
         && (0 <= bestcand && bestcand < nlpcands) ) /* no valid candidate index could be found */
      {
         SCIP_VAR* var;
         SCIP_Real val;

         var = lpcands[bestcand];
         val = lpcandssol[bestcand];

         SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%g)\n",
            nlpcands, bestcand, SCIPvarGetName(var), val);

         /* execute the branching as a result of the branching logic */
         SCIP_CALL( branchOnVar(scip, var, val, bestdown, bestdownvalid, bestup, bestupvalid, provedbound) );
         *result = SCIP_BRANCHED;
      }

      SCIPstatistic(
         if( *result == SCIP_BRANCHED )
         {
            SCIPdebugMessage("Result: Finished LookaheadBranching by branching.\n");
         }
         if( *result == SCIP_REDUCEDDOM )
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
         else if( status->depthtoosmall )
         {
            SCIPdebugMessage("Result: The branching depth wasn't high enough for a 2 level branching.\n");
         }
         else
         {
            SCIPdebugMessage("Result: Could not find any variable to branch on.\n");
         }
      )

      SCIPfreeBufferArray(scip, &lpcandsfrac);
      SCIPfreeBufferArray(scip, &lpcandssol);
      SCIPfreeBufferArray(scip, &lpcands);

      freeStatus(scip, &status);
   }

   if( branchruledata->useimpliedbincons )
   {
      SCIP_CALL( SCIPfreeSol(scip, &baselpsol) );
   }

   SCIPstatistic(
      branchruledata->nresults[*result]++;
   )

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

   SCIPstatistic(
      branchruledata->nfirstlvlcutoffs = 0;
      branchruledata->nfirstlvllps = 0;
      branchruledata->nsecondlvlcutoffs = 0;
      branchruledata->nsecondlvllps = 0;
      branchruledata->nbinconst = 0;
      branchruledata->nstoflvlcutoffs = 0;
      /* 17 current number of possible result values and the index is 1 based, so 17 + 1 as array size with unused 0 element */
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nresults, 17 + 1) );
   )

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
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/addbinconsrow",
         "should implied binary constraints be added as a row to the LP?",
         &branchruledata->addbinconsrow, TRUE, DEFAULT_ADDBINCONSROW, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead/maxnumberviolatedcons",
         "how many constraints that are violated by the base lp solution should be gathered until they are added?",
         &branchruledata->maxnviolatedcons, TRUE, DEFAULT_MAXNUMBERVIOLATEDCONS, -1, 100000000, NULL, NULL) );
   return SCIP_OKAY;
}
