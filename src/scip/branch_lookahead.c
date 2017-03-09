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
//#define PRINTNODECONS
/*
*/
#define SCIP_DEBUG
#define SCIP_STATISTIC

/**@file   branch_lookahead.c
 * @brief  lookahead branching rule
 * @author Christoph Schubert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_lookahead.h"
#include "scip/clock.h"
#include "scip/common_branch_lookahead.h"
#include "scip/cons_setppc.h"
#include "scip/cons_logicor.h"
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

#define DEFAULT_USEIMPLIEDBINARYCONSTRAINTS  TRUE
#define DEFAULT_ADDBINCONSROW                FALSE
#define DEFAULT_MAXNUMBERVIOLATEDCONS        10000
#define DEFAULT_REEVALAGE                    10LL
#define DEFAULT_STOPBRANCHING                TRUE
#define DEFAULT_PURGECONSONCUTOFF            TRUE
#define DEFAULT_FORCEBRANCHING               FALSE
#define DEFAULT_ONLYFULLSTRONG               FALSE
#define DEFAULT_RECURSIONDEPTH               2
#define DEFAULT_USEDOMAINREDUCTION           TRUE
#define DEFAULT_ADDNONVIOCONS                FALSE
#define DEFAULT_DOWNFIRST                    TRUE

/*
 * Data structures
 */


/**
 * A container to hold the result of a branching.
 */
typedef struct
{
   SCIP_Real             objval;             /**< The objective value of the solved lp. Only contains meaningful data, if
                                              *   cutoff == FALSE. */
   SCIP_Bool             cutoff;             /**< Indicates whether the node was infeasible and was cutoff. */
   SCIP_Bool             dualboundvalid;     /**< Us the value of the dual bound valid? That means, was the according LP
                                              *   or the sub problems solved to optimality? */
   SCIP_Real             dualbound;          /**< The best dual bound for this branching, may be changed by lower level
                                              *   branchings. */
} BRANCHINGRESULTDATA;

/**
 * The parameter that can be changed by the user/caller and alter the behaviour of the lookahead branching.
 */
typedef struct
{
   SCIP_Bool             usedomainreduction; /**< indicates whether the data for domain reductions should be gathered and
                                              *   used. */
   SCIP_Bool             usebincons;         /**< indicates whether the data for the implied binary constraints should
                                              *   be gathered and used */
   SCIP_Bool             addbinconsrow;      /**< Add the implied binary constraints as a row to the problem matrix */
   SCIP_Bool             stopbranching;      /**< indicates whether we should stop the first level branching after finding
                                              *   an infeasible first branch */
   SCIP_Bool             purgeconsoncutoff;  /**< Clear the list of already found constraints for a specific variable, if we
                                              *   found out, that the problem is infeasible. */
   SCIP_Longint          reevalage;          /**< The number of "normal" (not probing) lps that may have been solved before
                                              *   we stop using old data and start recalculating new first level data. */
   int                   maxnviolatedcons;   /**< The number of constraints we want to gather before restarting the run. Set
                                              *   to -1 for an unbounded list. */
   SCIP_Bool             forcebranching;     /**< Execute the lookahead logic even if only one branching candidate is given.
                                              *   May be used to calculate the score of a single candidate. */
   int                   recursiondepth;     /**< How deep should the recursion go? Default for Lookahead: 2 */
   SCIP_Bool             addnonviocons;      /**< Should constraints be added, that are not violated by the base LP? */
   SCIP_Bool             downfirst;          /**< Should the down branch be executed first? */
} CONFIGURATION;

#ifdef SCIP_STATISTIC
/**
 * The data used for some statistical analysis.
 */
typedef struct
{
   int                   ntotalresults;      /**< The total sum of the entries in nresults. */
   int*                  nresults;           /**< Array of counters for each result state the lookahead branching finished.
                                              *   The first (0) entry is unused, as the result states are indexed 1-based
                                              *   and we use this index as our array index. */
   int*                  nsinglecutoffs;     /**< The number of single cutoffs on a (probing) node per probingdepth. */
   int*                  nfullcutoffs;       /**< The number of double cutoffs on a (probing) node per probingdepth. */
   int*                  nlpssolved;         /**< The number of lps solved for a given probingdepth. */
   int*                  npropdomred;        /**< The number of domain reductions based on domain propagation per
                                              *   progingdepth. */
   int*                  nsinglecandidate;   /**< The number of times a single candidate was given to the recursion routine
                                              *   per probingdepth. */
   int*                  noldbranchused;     /**< The number of times old branching data is used (see the reevalage
                                              *   parameter in the CONFIGURATION struct) */
   int                   nbinconst;          /**< The number of binary constraints added to the base node. */
   int                   nbinconstvio;       /**< The number of binary constraints added to the base node, that are violated
                                              *   by the LP at that node. */
   int                   ndomred;            /**< The number of domain reductions added to the base node. */
   int                   ndomredvio;         /**< The number of domain reductions added to the base node, that are violated
                                              *   by the LP at that node. */
   int                   ndepthreached;      /**< The number of times the branching was aborted due to a too small depth. */
   int                   ndomredcons;        /**< The number of binary constraints ignored, as they would be dom reds. */
   int                   ncutoffproofnodes;  /**< The number of nodes needed to proof all found cutoffs. */
   int                   ndomredproofnodes;  /**< The number of nodes needed to proof all found domreds. */
} STATISTICS;

/**
 * Helper struct to store the statistical data needed in a single run.
 */
typedef struct
{
   int                   ncutoffproofnodes;  /**< The number of nodes needed to proof the current cutoff. */
   int                   ndomredproofnodes;  /**< The number of nodes needed to proof all current domreds. */
} LOCALSTATISTICS;
#endif

/**
 * The data that is preserved over multiple runs of the branching rule.
 */
typedef struct
{
   SCIP_SOL*             prevbinsolution;    /**< The previous solution for the case that in the previous run only
                                              *   non-violating implied binary constraints were added. */
   BRANCHINGDECISION*    prevdecision;       /**< The previous decision that gets used for the case that in the previous run
                                              *   only non-violating implied binary constraints were added.*/
   int                   restartindex;       /**< The index at which the iteration over the number of candidates starts. */
   SCIP_Longint*         lastbranchid;       /**< The node id at which the var was last branched on (for a given branching
                                              *   var). */
   SCIP_Longint*         lastbranchnlps;     /**< The number of (non-probing) LPs that where solved when the var was last
                                              *   branched on. */
   BRANCHINGRESULTDATA** lastbranchupres;    /**< The result of the last up branching for a given var. */
   BRANCHINGRESULTDATA** lastbranchdownres;  /**< The result of the last down branching for a given var. */
} PERSISTENTDATA;

/**
 * branching rule data
 */
struct SCIP_BranchruleData
{
   SCIP_Bool             isinitialized;      /**< indicates whether the fields in this struct are initialized */
   CONFIGURATION*        config;             /**< the parameter that influence the behaviour of the lookahead branching */
   PERSISTENTDATA*       persistent;         /**< the data that persists over multiple branching decisions */
#ifdef SCIP_STATISTIC
   STATISTICS*           statistics;         /**< statistical data container */
#endif
};

/**
 * information about the current status of the branching
 */
typedef struct
{
   SCIP_Bool             addbinconst;        /**< Was a binary constraint added? */
   SCIP_Bool             depthtoosmall;      /**< Was the remaining depth too small to branch on? */
   SCIP_Bool             lperror;            /**< Resulted a LP solving in an error? */
   SCIP_Bool             cutoff;             /**< Was the current node cutoff? */
   SCIP_Bool             domredcutoff;       /**< Was the current node cutoff due to domain reductions? */
   SCIP_Bool             domred;             /**< Were domain reductions added due to information obtained through
                                              *   branching? */
   SCIP_Bool             propagationdomred;  /**< Were domain reductions added due to domain propagation? */
   SCIP_Bool             limitreached;       /**< Was a limit (time, node, user, ...) reached? */
   SCIP_Bool             maxnconsreached;    /**< Was the max number of constraints(bin const and dom red) reached? */
} STATUS;

/**
 * all domain reductions found through cutoff of branches
 */
typedef struct
{
   SCIP_Real*            lowerbounds;        /**< The new lower bounds found for each variable in the problem. */
   SCIP_Bool*            lowerboundset;      /**< Indicates whether the lower bound may be added to the base node. */
   SCIP_Real*            upperbounds;        /**< The new upper bounds found for each variable in the problem. */
   SCIP_Bool*            upperboundset;      /**< Indicates whether the upper bound may be added to the base node. */
   SCIP_Bool*            baselpviolated;     /**< Indicates whether the base lp solution violates the new bounds of a var.*/
   int                   nviolatedvars;      /**< Tracks the number of vars that have a violated (by the base lp) new lower
                                              *   or upper bound. */
   int                   nchangedvars;       /**< Tracks the number of vars, that have a changed domain. (a change on both,
                                              *   upper and lower bound, counts as one.) */
#ifdef SCIP_STATISTIC
   int*                  lowerboundnproofs;  /**< The number of nodes needed to proof the lower bound for each variable. */
   int*                  upperboundnproofs;  /**< The number of nodes needed to proof the upper bound for each variable. */
#endif
} DOMAINREDUCTIONS;

/**
 * all constraints that were created and may be added to the base node
 */
typedef struct
{
   SCIP_CONS**           constraints;        /**< The array of constraints. Length is adjusted as needed. */
   int                   nconstraints;       /**< The number of entries in the array 'constraints'. */
   int                   memorysize;         /**< The number of entries that the array 'constraints' may hold before the
                                              *   array is reallocated. */
   int                   nviolatedcons;      /**< Tracks the number of constraints, that are violated by the base LP
                                              *   solution. */
} CONSTRAINTLIST;

/**
 * list of binary variables currently branched on
 * a down branching (x <= 0) is saved as the negated variable (1-x)
 * an up branching (x >= 1) is saved as the original variable (x)
 * these variables are used to build the binary constraint in case that a ('binary') brunch is cutoff
 */
typedef struct
{
   SCIP_VAR**            binaryvars;         /**< The binary variables currently branched on. */
   int                   nbinaryvars;        /**< The number of entries in 'nbinaryvars'. */
   int                   memorysize;         /**< The number of entries that the array 'binaryvars' may hold before the
                                              *   array is reallocated. */
} BINARYVARLIST;

/**
 * struct holding the relevant data for handling binary constraints
 */
typedef struct
{
   BINARYVARLIST*        binaryvars;         /**< The current binary vars, used to created the constraints. */
   CONSTRAINTLIST*       createdconstraints; /**< The created constraints. */
} BINCONSDATA;

/**
 * allocate the struct on the Buffer and initialize it with the default values
 */
static
SCIP_RETCODE allocDomainReductions(
   SCIP*                 scip,               /**< SCIP data structure */
   DOMAINREDUCTIONS**    domreds             /**< The struct that has to be allocated and initialized. */
   )
{
   int ntotalvars;
   int i;

   assert(scip != NULL);
   assert(domreds != NULL);

   /* The arrays saves the data for all variables in the problem via the ProbIndex. See SCIPvarGetProbindex() */
   ntotalvars = SCIPgetNVars(scip);

   /* Allocate the struct and the contained arrays. */
   SCIP_CALL( SCIPallocBuffer(scip, domreds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*domreds)->lowerbounds, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*domreds)->upperbounds, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*domreds)->lowerboundset, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*domreds)->upperboundset, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*domreds)->baselpviolated, ntotalvars) );
   SCIPstatistic(
      SCIP_CALL( SCIPallocBufferArray(scip, &(*domreds)->lowerboundnproofs, ntotalvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(*domreds)->upperboundnproofs, ntotalvars) );
   )

   /* Initialize the validity arrays to FALSE, such that the (undefined) starting values are not used. */
   for( i = 0; i < ntotalvars; i++ )
   {
      (*domreds)->lowerboundset[i] = FALSE;
      (*domreds)->upperboundset[i] = FALSE;
      (*domreds)->baselpviolated[i] = FALSE;
      SCIPstatistic(
         (*domreds)->lowerboundnproofs[i] = 0;
         (*domreds)->upperboundnproofs[i] = 0;
      )
   }

   /* At the start we have no domain reductions for any variable. */
   (*domreds)->nviolatedvars = 0;
   (*domreds)->nchangedvars = 0;

   return SCIP_OKAY;
}

/*
static
void printDomainReductions(
   SCIP*                 scip,
   DOMAINREDUCTIONS*     domreds
   )
{
   int ntotalvars;
   int i;
   SCIP_VAR** vars;

   assert(scip != NULL);
   assert(domreds != NULL);

   ntotalvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   for( i = 0; i < ntotalvars; i++)
   {
      if( domreds->lowerboundset[i] && domreds->upperboundset[i] )
      {
         SCIPinfoMessage(scip, NULL, "%5g <= %s <= %5g\n", domreds->lowerbounds[i], SCIPvarGetName(vars[i]), domreds->upperbounds[i]);
      }
      else if( domreds->lowerboundset[i] )
      {
         SCIPinfoMessage(scip, NULL, "%5g <= %s\n", domreds->lowerbounds[i], SCIPvarGetName(vars[i]));
      }
      else if( domreds->upperboundset[i] )
      {
         SCIPinfoMessage(scip, NULL, "         %s <= %5g\n", SCIPvarGetName(vars[i]), domreds->upperbounds[i]);
      }
   }
}
*/

static
void addLowerBoundProofNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< The variable the bound should be added for. */
   SCIP_Real             lowerbound,         /**< The new lower bound for the variable. */
   SCIP_SOL*             baselpsol,          /**< The LP solution of the base problem. Used to check whether the domain
                                              *   reduction is violated by it. */
   DOMAINREDUCTIONS*     domainreductions    /**< The struct the domain reduction should be added to. */
#ifdef SCIP_STATISTIC
   ,int                  nproofnodes         /**< The number of nodes needed to proof the new lower bound. */
#endif
)
{
   int varindex;
   SCIP_Real basesolutionval;
   SCIP_Real newlowerbound;
   SCIPstatistic( int newnproof; )

   assert(scip != NULL);
   assert(var != NULL);
   assert(baselpsol != NULL);
   assert(domainreductions != NULL);

   /* The arrays inside DOMAINREDUCTIONS are indexed via the problem index. */
   varindex = SCIPvarGetProbindex(var);

   /* We get the solution value to check whether the domain reduction is violated in the base LP */
   basesolutionval = SCIPgetSolVal(scip, baselpsol, var);

   /* If we have an old lower bound we take the stronger one, so the MAX is taken. Otherwise we use the new lower bound
    * directly. */
   if( domainreductions->lowerboundset[varindex] )
   {
      if( SCIPisGE(scip, domainreductions->lowerbounds[varindex], lowerbound) )
      {
         newlowerbound = domainreductions->lowerbounds[varindex];
         SCIPstatistic(
            if( SCIPisEQ(scip, domainreductions->lowerbounds[varindex], lowerbound) )
            {
               newnproof = MIN(domainreductions->lowerboundnproofs[varindex], nproofnodes);
            }
            else
            {
               newnproof = domainreductions->lowerboundnproofs[varindex];
            }
         )
      }
      else
      {
         newlowerbound = lowerbound;
         SCIPstatistic( newnproof = nproofnodes; )

         if( !domainreductions->upperboundset[varindex] )
         {
            domainreductions->nchangedvars++;
         }
      }
   }
   else
   {
      newlowerbound = lowerbound;
      SCIPstatistic( newnproof = nproofnodes; )
   }
   domainreductions->lowerbounds[varindex] = newlowerbound;
   domainreductions->lowerboundset[varindex] = TRUE;
   SCIPstatistic( domainreductions->lowerboundnproofs[varindex] = newnproof; )

   /* In case the new lower bound is greater than the base solution val and the base solution val is not violated by a
    * previously found bound, we increment the nviolatedvars counter and set the baselpviolated flag.  */
   if( SCIPisGT(scip, newlowerbound, basesolutionval) && !domainreductions->baselpviolated[varindex] )
   {
      domainreductions->baselpviolated[varindex] = TRUE;
      domainreductions->nviolatedvars++;
   }
}

/**
 * add a lower bound to the DOMAINREDUCTIONS struct
 */
static
void addLowerBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< The variable the bound should be added for. */
   SCIP_Real             lowerbound,         /**< The new lower bound for the variable. */
   SCIP_SOL*             baselpsol,          /**< The LP solution of the base problem. Used to check whether the domain
                                              *   reduction is violated by it. */
   DOMAINREDUCTIONS*     domainreductions    /**< The struct the domain reduction should be added to. */
   )
{
   /* We add the lower bound with number of proof nodes2, as this method is only called from the recursion directly. There it
    * is called in case that only one child node is cutoff. As proof nodes we count the cutoff node as well as the valid
    * node, as we need both to proof, that this domain reduction is valid. */
#ifdef SCIP_STATISTIC
   addLowerBoundProofNode(scip, var, lowerbound, baselpsol, domainreductions, 2);
#else
   addLowerBoundProofNode(scip, var, lowerbound, baselpsol, domainreductions);
#endif
}

/**
 * add an upper bound to the DOMAINREDUCTIONS struct
 */
static
void addUpperBoundProofNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< The variable the bound should be added for. */
   SCIP_Real             upperbound,         /**< The new upper bound for the variable. */
   SCIP_SOL*             baselpsol,          /**< The LP solution of the base problem. Used to check whether the domain
                                              *   reduction is violated by it. */
   DOMAINREDUCTIONS*     domainreductions    /**< The struct the domain reduction should be added to. */
#ifdef SCIP_STATISTIC
   ,int                  nproofnodes         /**< The number of nodes needed to proof the new lower bound. */
#endif
   )
{
   int varindex;
   SCIP_Real basesolutionval;
   SCIP_Real newupperbound;
   SCIPstatistic( int newnproof; )

   assert(scip != NULL);
   assert(var != NULL);
   assert(baselpsol != NULL);
   assert(domainreductions != NULL);

   /* The arrays inside DOMAINREDUCTIONS are indexed via the problem index. */
   varindex = SCIPvarGetProbindex(var);

   /* We get the solution value to check whether the domain reduction is violated in the base LP */
   basesolutionval = SCIPgetSolVal(scip, baselpsol, var);

   /* If we hava an old upper bound we take the stronger one, so the MIN is taken. Otherwise we use the new upper bound
    * directly. */
   if( domainreductions->upperboundset[varindex] )
   {
      if( SCIPisLE(scip, domainreductions->upperbounds[varindex], upperbound) )
      {
         newupperbound = domainreductions->upperbounds[varindex];

         SCIPstatistic(
            if( SCIPisEQ(scip, domainreductions->upperbounds[varindex], upperbound) )
            {
               newnproof = MIN(domainreductions->upperboundnproofs[varindex], nproofnodes);
            }
            else
            {
               newnproof = domainreductions->upperboundnproofs[varindex];
            }
         )
      }
      else
      {
         newupperbound = upperbound;
         SCIPstatistic( newnproof = nproofnodes; )
      }
   }
   else
   {
      newupperbound = upperbound;
      SCIPstatistic( newnproof = nproofnodes; )

      if( !domainreductions->lowerboundset[varindex] )
      {
         domainreductions->nchangedvars++;
      }
   }

   domainreductions->upperbounds[varindex] = newupperbound;
   domainreductions->upperboundset[varindex] = TRUE;
   SCIPstatistic( domainreductions->upperboundnproofs[varindex] = newnproof; )

   /* In case the new upper bound is lesser than the base solution val and the base solution val is not violated by a
    * previously found bound, we increment the nviolatedvars counter and set the baselpviolated flag.  */
   if( SCIPisLT(scip, newupperbound, basesolutionval) && !domainreductions->baselpviolated[varindex] )
   {
      domainreductions->baselpviolated[varindex] = TRUE;
      domainreductions->nviolatedvars++;
   }
}

static
void addUpperBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< The variable the bound should be added for. */
   SCIP_Real             upperbound,         /**< The new upper bound for the variable. */
   SCIP_SOL*             baselpsol,          /**< The LP solution of the base problem. Used to check whether the domain
                                              *   reduction is violated by it. */
   DOMAINREDUCTIONS*     domainreductions    /**< The struct the domain reduction should be added to. */
   )
{
   /* We add the upper bound with number of proof nodes 2, as this method is only called from the recursion directly. There it
    * is called in case that only one child node is cutoff. As proof nodes we count the cutoff node as well as the valid
    * node, as we need both to proof, that this domain reduction is valid. */
#ifdef SCIP_STATISTIC
   addUpperBoundProofNode(scip, var, upperbound, baselpsol, domainreductions, 2);
#else
   addUpperBoundProofNode(scip, var, upperbound, baselpsol, domainreductions);
#endif
}

/**
 * merges the domain reduction data from the two given branching childs data into the target parent data
 */
static
void applyDeeperDomainReductions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             baselpsol,          /**< The LP solution of the base problem. Used to check whether the domain
                                              *   reduction is violated by it. */
   DOMAINREDUCTIONS*     targetdomreds,      /**< The target that should be filled with the merged data. */
   DOMAINREDUCTIONS*     downdomreds,        /**< One of the source DOMAINREDUCTIONS. */
   DOMAINREDUCTIONS*     updomreds           /**< The other source DOMAINREDUCTIONS. */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(baselpsol != NULL);
   assert(targetdomreds != NULL);
   assert(downdomreds != NULL);
   assert(updomreds != NULL);

   /* as the bounds are tracked for all vars we have to iterate over all vars */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   assert(vars != NULL);
   assert(nvars > 0);

   for( i = 0; i < nvars; i++ )
   {
      assert(vars[i] != NULL);

      /* If not both child branches have a lower bound for a var, we cannot apply the lower bound to the parent */
      if( downdomreds->lowerboundset[i] && updomreds->lowerboundset[i] )
      {
         SCIP_Real newlowerbound;
         SCIPstatistic( int newnproofs; )

         /* If both child branches have a lower bound for a var, the MIN of both values represents a valid lower bound */
         newlowerbound = MIN(downdomreds->lowerbounds[i], updomreds->lowerbounds[i]);
         SCIPstatistic( newnproofs = downdomreds->lowerboundnproofs[i] + updomreds->lowerboundnproofs[i] + 2; )

         /* This MIN can now be added via the default add method */
#ifdef SCIP_STATISTIC
         addLowerBoundProofNode(scip, vars[i], newlowerbound, baselpsol, targetdomreds, newnproofs);
#else
         addLowerBoundProofNode(scip, vars[i], newlowerbound, baselpsol, targetdomreds);
#endif
      }

      /* If not both child branches have a lower bound for a var, we cannot apply the lower bound to the parent */
      if( downdomreds->upperboundset[i] && updomreds->upperboundset[i] )
      {
         SCIP_Real newupperbound;
         SCIPstatistic( int newnproofs; )

         /* If both child branches have an upper bound for a var, the MAX of both values represents a valid upper bound */
         newupperbound = MAX(downdomreds->upperbounds[i], updomreds->upperbounds[i]);
         SCIPstatistic( newnproofs = downdomreds->upperboundnproofs[i] + updomreds->upperboundnproofs[i] + 2; )

         /* This MAX can now be added via the default add method */
#ifdef SCIP_STATISTIC
         addUpperBoundProofNode(scip, vars[i], newupperbound, baselpsol, targetdomreds, newnproofs);
#else
         addUpperBoundProofNode(scip, vars[i], newupperbound, baselpsol, targetdomreds);
#endif
      }
   }
}

#ifdef SCIP_STATISTIC
static
void printStatistics(
   SCIP*                 scip,
   STATISTICS*           statistics,
   int                   recursiondepth
   );
#endif

/**
 * Applies the domain reductions the current node.
 */
static
SCIP_RETCODE applyDomainReductions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             baselpsol,          /**< The LP solution of the base problem. Used to check whether the domain
                                              *   reduction is violated by it. */
   DOMAINREDUCTIONS*     domreds,            /**< The domain reductions that should be applied the current node. */
   SCIP_Bool*            domredcutoff,       /**< Pointer to store whether a cutoff was found due to domain reductions */
   SCIP_Bool*            domred              /**< Pointer to store whether a strict cutoff was added. (Strict in the sense
                                              *   that the new domain is violated by the bese LP solution.) */
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics          /**< The statistics container. */
#endif
   )
{
   int i;
   SCIP_VAR** probvars;
   int nprobvars;
   int nboundsadded = 0;
   int nboundsaddedvio = 0;

   assert(scip != NULL);
   assert(baselpsol != NULL);
   assert(domreds != NULL);
   assert(domredcutoff != NULL);
   assert(domred != NULL);
   SCIPstatistic( assert(statistics != NULL) );

   /* initially we have no cutoff */
   *domredcutoff = FALSE;

   /* as the bounds are tracked for all vars we have to iterate over all vars */
   probvars = SCIPgetVars(scip);
   nprobvars = SCIPgetNVars(scip);

   assert(probvars != NULL);
   assert(nprobvars > 0);

   for( i = 0; i < nprobvars && !*domredcutoff; i++ )
   {
      SCIP_VAR* var;
      SCIP_Real baselpval;
      SCIP_Bool boundadded = FALSE;
      SCIP_Bool boundaddedvio = FALSE;

      var = probvars[i];

      assert(var != NULL);

      baselpval = SCIPgetSolVal(scip, baselpsol, var);

      if( !*domredcutoff && domreds->lowerboundset[i] )
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         SCIP_Real oldlowerbound;
         SCIP_Real proposedlowerbound;
         SCIP_Real newlowerbound;

         /* get the old and the new lower bound */
         oldlowerbound = SCIPvarGetLbLocal(var);
         proposedlowerbound = domreds->lowerbounds[i];

         /* add the new bound */
         SCIP_CALL( SCIPtightenVarLb(scip, var, proposedlowerbound, FALSE, &infeasible, &tightened) );

         newlowerbound = SCIPvarGetLbLocal(var);
         SCIPdebugMessage("Variable <%s>, old lower bound <%g>, proposed lower bound <%g>, new lower bound <%g>\n",
            SCIPvarGetName(var), oldlowerbound, proposedlowerbound, newlowerbound);

         if( infeasible )
         {
            /* the domain reduction may result in an empty model (ub < lb) */
            *domredcutoff = TRUE;
            SCIPdebugMessage("The domain reduction of variable <%s> resulted in an empty model.\n",
               SCIPvarGetName(var));
         }
         else if( tightened )
         {
            /* the lb is now strictly greater than before */
            SCIPdebugMessage("The lower bound of variable <%s> was successfully tightened.\n", SCIPvarGetName(var));
            boundadded = TRUE;
            *domred = TRUE;
            SCIPstatistic( statistics->ndomredproofnodes += domreds->lowerboundnproofs[i]; )

            if( SCIPisLT(scip, baselpval, newlowerbound) )
            {
               SCIPdebugMessage("The lower bound of variable <%s> is violated by the base lp value <%g>.\n",
                  SCIPvarGetName(var), baselpval);
               boundaddedvio = TRUE;
            }
         }
      }

      if( !*domredcutoff && domreds->upperboundset[i] )
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         SCIP_Real oldupperbound;
         SCIP_Real proposedupperbound;
         SCIP_Real newupperbound;

         /* get the old and the new upper bound */
         oldupperbound = SCIPvarGetUbLocal(var);
         proposedupperbound = domreds->upperbounds[i];

         /* add the new bound */
         SCIP_CALL( SCIPtightenVarUb(scip, var, proposedupperbound, FALSE, &infeasible, &tightened) );

         newupperbound = SCIPvarGetUbLocal(var);
         SCIPdebugMessage("Variable <%s>, old upper bound <%g>, proposed upper bound <%g>, new upper bound <%g>\n",
            SCIPvarGetName(var), oldupperbound, proposedupperbound, newupperbound);

         if( infeasible )
         {
            /* the domain reduction may result in an empty model (ub < lb) */
            *domredcutoff = TRUE;
            SCIPdebugMessage("The domain reduction of variable <%s> resulted in an empty model.\n",
               SCIPvarGetName(var));
         }
         else if( tightened )
         {
            /* the ub is now strictly smaller than before */
            SCIPdebugMessage("The upper bound of variable <%s> was successfully tightened.\n", SCIPvarGetName(var));
            boundadded = TRUE;
            *domred = TRUE;
            SCIPstatistic( statistics->ndomredproofnodes += domreds->upperboundnproofs[i]; )

            if( SCIPisGT(scip, baselpval, newupperbound) )
            {
               SCIPdebugMessage("The upper bound of variable <%s> is violated by the base lp value <%g>.\n",
                  SCIPvarGetName(var), baselpval);
               boundaddedvio = TRUE;
            }
         }
      }

      /* We increment the number of bounds added at most once per var */
      if( boundadded )
      {
         nboundsadded++;
      }

      /* We increment the number of bounds violated by the base lp at most once per var */
      if( boundaddedvio )
      {
         nboundsaddedvio++;
      }
   }

   SCIPstatistic(
      statistics->ndomred += nboundsadded;
      statistics->ndomredvio += nboundsaddedvio;
   )

   SCIPdebugMessage("Truly changed <%d> domains of the problem, <%d> of them are violated by the base lp.\n", nboundsadded,
      nboundsaddedvio);

   return SCIP_OKAY;
}

/**
 * frees the given DOMAINREDUCTIONS and all contained Arrays in the opposite order of allocation
 */
static
void freeDomainReductions(
   SCIP*                 scip,               /**< SCIP data structure */
   DOMAINREDUCTIONS**    domreds             /**< Pointer to the struct to be freed. */
   )
{
   assert(scip != NULL);
   assert(domreds != NULL);

   SCIPstatistic(
      SCIPfreeBufferArray(scip, &(*domreds)->upperboundnproofs);
      SCIPfreeBufferArray(scip, &(*domreds)->lowerboundnproofs);
   )
   SCIPfreeBufferArray(scip, &(*domreds)->baselpviolated);
   SCIPfreeBufferArray(scip, &(*domreds)->upperboundset);
   SCIPfreeBufferArray(scip, &(*domreds)->lowerboundset);
   SCIPfreeBufferArray(scip, &(*domreds)->upperbounds);
   SCIPfreeBufferArray(scip, &(*domreds)->lowerbounds);
   SCIPfreeBuffer(scip, domreds);
}

/**
 * Allocates and initializes the BINARYVARLIST struct.
 */
static
SCIP_RETCODE allocBinaryVarList(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYVARLIST**       list,               /**< Pointer to the list to be allocated and initialized. */
   int                   startsize           /**< The number of entries the list initially can hold. */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, list) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*list)->binaryvars, startsize) );

   /* We start with no entries and the (current) max length */
   (*list)->nbinaryvars = 0;
   (*list)->memorysize = startsize;

   return SCIP_OKAY;
}

/**
 * Appends a binary variable to the list, reallocating the list if necessary.
 */
static
SCIP_RETCODE appendToBinaryVarList(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYVARLIST*        list,               /**< The list to add the var to. */
   SCIP_VAR*             vartoadd            /**< The binary var to add to the list. */
   )
{
   assert(list != NULL);
   assert(vartoadd != NULL);
   assert(SCIPvarIsBinary(vartoadd));

   /* In case the list tries to hold more elements than it has space, reallocate  */
   if( list->memorysize == list->nbinaryvars )
   {
      /* resize the array, such that it can hold at least the new element */
      int newmemsize = SCIPcalcMemGrowSize(scip, list->memorysize + 1);
      SCIP_CALL( SCIPreallocBufferArray(scip, &list->binaryvars, newmemsize) );
      list->memorysize = newmemsize;
   }

   /* Set the new var at the first unused place, which is the length used as index */
   list->binaryvars[list->nbinaryvars] = vartoadd;
   list->nbinaryvars++;

   return SCIP_OKAY;
}

/**
 * Remove and return the last element from the list.
 */
static
SCIP_VAR* dropFromBinaryVarList(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYVARLIST*        list                /**< The list to remove the last element from. */
   )
{
   SCIP_VAR* lastelement;

   assert(list != NULL);
   assert(list->nbinaryvars > 0);
   assert(list->binaryvars[list->nbinaryvars-1] != NULL);

   /* get the last element and set the last pointer to NULL (maybe unnecessary, but feels cleaner) */
   lastelement = list->binaryvars[list->nbinaryvars-1];
   list->binaryvars[list->nbinaryvars-1] = NULL;

   assert(lastelement != NULL);

   /* decrement the number of entries */
   list->nbinaryvars--;

   return lastelement;
}

/**
 * Frees all resources allocated by a BINARYVARLIST in opposite order of allocation.
 */
static
void freeBinaryVarList(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYVARLIST**       list                /**< Pointer to the list to free */
   )
{
   assert(scip != NULL);
   assert(list != NULL);

   SCIPfreeBufferArray(scip, &(*list)->binaryvars);
   SCIPfreeBuffer(scip, list);
}

/**
 * Allocate and initialize the list holding the constraints.
 */
static
SCIP_RETCODE allocConstraintList(
   SCIP*                 scip,               /**< SCIP data structure */
   CONSTRAINTLIST**      conslist,           /**< Pointer to the list to be allocated and initialized. */
   int                   startsize           /**< The number of entries the list initially can hold. */
   )
{
   assert(scip != NULL);
   assert(conslist != NULL);
   assert(startsize > 0);

   SCIP_CALL( SCIPallocBuffer(scip, conslist) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*conslist)->constraints, startsize) );

   /* We start without any constraints */
   (*conslist)->nconstraints = 0;
   (*conslist)->memorysize = startsize;
   (*conslist)->nviolatedcons = 0;

   return SCIP_OKAY;
}

/**
 * Append an element to the end of the list of constraints.
 */
static
SCIP_RETCODE appendToConstraintList(
   SCIP*                 scip,               /**< SCIP data structure */
   CONSTRAINTLIST*       list,               /**< The list to add the element to. */
   SCIP_CONS*            constoadd           /**< The element to add to the list. */
   )
{
   assert(scip != NULL);
   assert(list != NULL);
   assert(constoadd != NULL);

   /* In case the list tries to hold more elements than it has space, reallocate  */
   if( list->memorysize == list->nconstraints )
   {
      /* resize the array, such that it can hold the new element */
      int newmemsize = SCIPcalcMemGrowSize(scip, list->memorysize + 1);
      SCIP_CALL( SCIPreallocBufferArray(scip, &list->constraints, newmemsize) );
      list->memorysize = newmemsize;
   }

   /* Set the new var at the first unused place, which is the length used as index */
   list->constraints[list->nconstraints] = constoadd;
   list->nconstraints++;

   return SCIP_OKAY;
}

/**
 * Free all resources of a constraint list in opposite order to the allocation.
 */
static
void freeConstraintList(
   SCIP*                 scip,               /**< SCIP data structure */
   CONSTRAINTLIST**      conslist            /**< Pointer to the list to be freed. */
   )
{
   assert(scip != NULL);
   assert(conslist != NULL);

   SCIPfreeBufferArray(scip, &(*conslist)->constraints);
   SCIPfreeBuffer(scip, conslist);
}

/**
 * Allocate and initialize the BINCONSDATA struct.
 */
static
SCIP_RETCODE allocBinConsData(
   SCIP*                 scip,               /**< SCIP data structure */
   BINCONSDATA**         consdata,           /**< Pointer to the struct to be allocated and initialized. */
   int                   maxdepth,           /**< The depth of the recursion as an upper bound of branch vars to hold. */
   int                   nstartcons          /**< The start size of the array containing the constraints. */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(maxdepth > 0);
   assert(nstartcons > 0);

   SCIP_CALL( SCIPallocBuffer(scip, consdata) );
   SCIP_CALL( allocBinaryVarList(scip, &(*consdata)->binaryvars, maxdepth) );
   SCIP_CALL( allocConstraintList(scip, &(*consdata)->createdconstraints, nstartcons) );

   return SCIP_OKAY;
}

/**
 * Free all resources in a BINCONSDATA in opposite order of allocation.
 */
static
void freeBinConsData(
   SCIP*                 scip,               /**< SCIP data structure */
   BINCONSDATA**         consdata            /**< Pointer to he struct to be freed. */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   freeConstraintList(scip, &(*consdata)->createdconstraints);
   freeBinaryVarList(scip, &(*consdata)->binaryvars);
   SCIPfreeBuffer(scip, consdata);
}

#ifdef SCIP_STATISTIC
static
void initStatistics(
   SCIP*                 scip,
   STATISTICS*           statistics,
   int                   recursiondepth
   )
{
   int i;

   statistics->ntotalresults = 0;
   statistics->nbinconst = 0;
   statistics->nbinconstvio = 0;
   statistics->ndomredvio = 0;
   statistics->ndepthreached = 0;
   statistics->ndomred = 0;
   statistics->ndomredcons = 0;
   statistics->ncutoffproofnodes = 0;
   statistics->ndomredproofnodes = 0;

   for( i = 0; i < 18; i++)
   {
      statistics->nresults[i] = 0;
   }

   for( i = 0; i < recursiondepth; i++ )
   {
      statistics->noldbranchused[i] = 0;
      statistics->nsinglecandidate[i] = 0;
      statistics->npropdomred[i] = 0;
      statistics->nfullcutoffs[i] = 0;
      statistics->nlpssolved[i] = 0;
      statistics->nsinglecutoffs[i] = 0;
   }
}
#endif

static
SCIP_RETCODE allocateStatus(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS**              status              /**< the status to be allocated */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, status) );

   (*status)->addbinconst = FALSE;
   (*status)->depthtoosmall = FALSE;
   (*status)->lperror = FALSE;
   (*status)->cutoff = FALSE;
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
 * Initiates the BranchingResultData struct.
 */
static
void initBranchingResultData(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGRESULTDATA*  resultdata          /**< The struct to be initiated. */
   )
{
   resultdata->objval = -SCIPinfinity(scip);
   resultdata->dualbound = -SCIPinfinity(scip);
   resultdata->cutoff = FALSE;
   resultdata->dualboundvalid = FALSE;
}

static
void copyBranchingResultData(
   BRANCHINGRESULTDATA*  sourcedata,
   BRANCHINGRESULTDATA*  targetdata
   )
{
   targetdata->cutoff = sourcedata->cutoff;
   targetdata->objval = sourcedata->objval;
}

/*
 * Local methods for the logic
 */

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

#ifdef SCIP_STATISTIC
static
SCIP_RETCODE allocateLocalStatistics(
   SCIP*                 scip,
   LOCALSTATISTICS**     localstats
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, localstats) );

   (*localstats)->ncutoffproofnodes = 0;
   (*localstats)->ndomredproofnodes = 0;

   return SCIP_OKAY;
}

static
void freeLocalStatistics(
   SCIP*                 scip,
   LOCALSTATISTICS**     localstats
   )
{
   SCIPfreeBuffer(scip, localstats);
}
#endif

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
   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "DownBranching: Var=<%s>, Proposed upper bound=<%g>, old bounds=[<%g>..<%g>], new bounds=[<%g>..<%g>]\n",
      SCIPvarGetName(branchvar), newupperbound, oldlowerbound, oldupperbound, oldlowerbound, newupperbound);

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
      /*SCIP_CALL( printCurrentSolution(scip) );*/

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
         resultdata->dualbound = SCIPgetLPObjval(scip);
         resultdata->dualboundvalid = TRUE;

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
   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "UpBranching: Var=<%s>, Proposed lower bound=<%g>, old bounds=[<%g>..<%g>], new bounds=[<%g>..<%g>]\n",
      SCIPvarGetName(branchvar), newlowerbound, oldlowerbound, oldupperbound, newlowerbound, oldupperbound);

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
      /*SCIP_CALL( printCurrentSolution(scip) );*/

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
         resultdata->dualbound = SCIPgetLPObjval(scip);
         resultdata->dualboundvalid = TRUE;
         resultdata->cutoff = resultdata->cutoff || SCIPisGE(scip, resultdata->objval, SCIPgetCutoffbound(scip));
         assert(solstat != SCIP_LPSOLSTAT_INFEASIBLE || resultdata->cutoff);
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE createBinaryConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,
   SCIP_CONS**           constraint,         /**< Pointer to store the created constraint in */
   char*                 constraintname,     /**< name of the new constraint */
   SCIP_VAR**            consvars,
   int                   nconsvars
   )
{
   SCIP_Bool initial = config->addbinconsrow;
   SCIP_Bool separate = config->addbinconsrow;
   SCIP_Bool enforce = FALSE;
   SCIP_Bool check = FALSE;
   SCIP_Bool propagate = TRUE;
   SCIP_Bool local = TRUE;
   SCIP_Bool modifiable = FALSE;
   SCIP_Bool dynamic = FALSE;
   SCIP_Bool removable = FALSE;
   SCIP_Bool stickingatnode = FALSE;

   SCIP_CALL( SCIPcreateConsLogicor(scip, constraint, constraintname, nconsvars, consvars, initial, separate, enforce,
         check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   return SCIP_OKAY;
}

/**
 * Create a name for the binary constraint.
 */
static
void createBinaryConstraintName(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            binaryvars,
   int                   nbinaryvars,
   char*                 constraintname      /**< the char pointer to store the name in */

   )
{
   int i;

   assert(scip != NULL);
   assert(binaryvars != NULL);
   assert(nbinaryvars > 0);
   assert(constraintname != NULL);

   sprintf(constraintname, "lookahead_bin_%s", SCIPvarGetName(binaryvars[0]));

   for( i = 1; i < nbinaryvars; i++ )
   {
      SCIP_VAR* var = binaryvars[i];
      const char* varname = SCIPvarGetName(var);

      sprintf(constraintname, "%s_%s", constraintname, varname);
   }
}

/**
 * Add the constraints found during the lookahead branching.
 * The implied binary bounds were found when two consecutive branching of binary variables were cutoff. Then these two
 * branching constraints can be combined into a single set packing constraint.
 */
static
SCIP_RETCODE addBinaryConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,
   BINCONSDATA*          binconsdata,
   SCIP_SOL*             baselpsol
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
#endif
   )
{
   /* If we only have one var for the constraint we can ignore it, as it is already added as a domain reduction. */
   if( binconsdata->binaryvars->nbinaryvars > 1 )
   {
      int i;
      SCIP_CONS* constraint;
      char constraintname[SCIP_MAXSTRLEN];
      SCIP_VAR** negatedvars;
      SCIP_Real lhssum = 0;

      SCIPdebugMessage("Adding binary constraint for <%i> vars.\n", binconsdata->binaryvars->nbinaryvars);

      SCIP_CALL( SCIPallocBufferArray(scip, &negatedvars, binconsdata->binaryvars->nbinaryvars) );
      for( i = 0; i < binconsdata->binaryvars->nbinaryvars; i++ )
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, binconsdata->binaryvars->binaryvars[i], &negatedvars[i]) );
         lhssum += SCIPgetSolVal(scip, baselpsol, negatedvars[i]);
      }

      if( config->addnonviocons || lhssum < 1 )
      {
         /* create a name for the new constraint */
         createBinaryConstraintName(scip, negatedvars, binconsdata->binaryvars->nbinaryvars, constraintname);
         /* create the constraint with the frehsly created name */
         SCIP_CALL( createBinaryConstraint(scip, config, &constraint, constraintname, negatedvars, binconsdata->binaryvars->nbinaryvars) );

#ifdef PRINTNODECONS
         SCIPinfoMessage(scip, NULL, "Created constraint:\n");
         SCIP_CALL( SCIPprintCons(scip, constraint, NULL) );
         SCIPinfoMessage(scip, NULL, "\n");
#endif

         SCIP_CALL( appendToConstraintList(scip, binconsdata->createdconstraints, constraint) );
      }

      /* the constraint we are building is a logic or: we have a list of binary variables that were
       * cutoff while we branched on with >= 1. So we have the constraint: x_1 + ... + x_n <= n-1.
       * Let y = (1-x), then we have an equivalent formulation: y_1 + ... + y_n >= 1. If the base lp
       * is violating this constraint we count this for our number of violated constraitns and bounds. */
      if( lhssum < 1 )
      {
         binconsdata->createdconstraints->nviolatedcons++;
      }

      SCIPfreeBufferArray(scip, &negatedvars);
   }
   else
   {
      SCIPstatistic(statistics->ndomredcons++;)
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE applyBinaryConstraints(
   SCIP*                 scip,
   SCIP_NODE*            basenode,
   CONSTRAINTLIST*       conslist,
   SCIP_Bool*            consadded
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
#endif
   )
{
   int i;

   SCIPdebugMessage("Adding <%i> binary constraints.\n", conslist->nconstraints);

   /* backwards to release the constraints in the reverse order of their creation. TODO: needed? */
   for( i = conslist->nconstraints-1; i >= 0; i-- )
   {
      SCIP_CONS* constraint = conslist->constraints[i];

#ifdef PRINTNODECONS
      SCIPprintCons(scip, constraint, NULL);
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* add the constraint to the given node */
      SCIP_CALL( SCIPaddConsNode(scip, basenode, constraint, NULL) );
      /* release the constraint, as it is no longer needed */
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

      *consadded = TRUE;
   }

   SCIPstatistic(
      statistics->nbinconst += conslist->nconstraints;
      statistics->nbinconstvio += conslist->nviolatedcons;
   );


   return SCIP_OKAY;
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

/** get a copy of the fractional candidates we can branch on
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
static
SCIP_RETCODE copyLPBranchCands(
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

   /* copy LP banching candidates and solution values, because they will be updated w.r.t. the strong branching LP
    * solution during the second level branchings */
   SCIP_CALL( SCIPduplicateBufferArray(scip, lpcands, tmplpcands, *nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, lpcandssol, tmplpcandssol, *nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, lpcandsfrac, tmplpcandsfrac, *nlpcands) );

   return SCIP_OKAY;
}

static
void freeLPBranchCands(
   SCIP*                 scip,
   SCIP_VAR***           lpcands,
   SCIP_Real**           lpcandssol,
   SCIP_Real**           lpcandsfrac
   )
{
   SCIPfreeBuffer(scip, lpcandsfrac);
   SCIPfreeBuffer(scip, lpcandssol);
   SCIPfreeBuffer(scip, lpcands);
}

static
SCIP_Bool isExecuteBranchingLoop(
   STATUS*               status
   )
{
   return !status->lperror && !status->cutoff && !status->limitreached && !status->maxnconsreached
         && !status->propagationdomred;
}

static
SCIP_Bool isUseOldBranching(
   SCIP*                 scip,
   CONFIGURATION*        config,
   PERSISTENTDATA*       persistent,
   SCIP_VAR*             branchvar
   )
{
   int varindex = SCIPvarGetProbindex(branchvar);

   SCIP_Longint currentnodeid = SCIPgetNNodes(scip);
   SCIP_Longint lastbranchingnodeid = persistent->lastbranchid[varindex];

   int reevalage = config->reevalage;
   SCIP_Longint currentnnodelps = SCIPgetNNodeLPs(scip);
   SCIP_Longint lastbranchingnnodelps = persistent->lastbranchnlps[varindex];

   return (currentnodeid == lastbranchingnodeid) && (currentnnodelps - lastbranchingnnodelps < reevalage);
}

static
void getOldBranching(
   SCIP*                 scip,
   PERSISTENTDATA*       persistent,
   SCIP_VAR*             branchvar,
   BRANCHINGRESULTDATA*  downbranchingresult,
   BRANCHINGRESULTDATA*  upbranchingresult
   )
{
   int varindex = SCIPvarGetProbindex(branchvar);

   copyBranchingResultData(persistent->lastbranchdownres[varindex], downbranchingresult);
   copyBranchingResultData(persistent->lastbranchupres[varindex], upbranchingresult);
}

static
void updateOldBranching(
   SCIP*                 scip,
   PERSISTENTDATA*       persistent,
   SCIP_VAR*             branchvar,
   BRANCHINGRESULTDATA*  downbranchingresult,
   BRANCHINGRESULTDATA*  upbranchingresult
   )
{
   int varindex = SCIPvarGetProbindex(branchvar);

   copyBranchingResultData(downbranchingresult, persistent->lastbranchdownres[varindex]);
   copyBranchingResultData(upbranchingresult, persistent->lastbranchupres[varindex]);

   persistent->lastbranchid[varindex] = SCIPgetNNodes(scip);
   persistent->lastbranchnlps[varindex] = SCIPgetNNodeLPs(scip);
}

static
SCIP_RETCODE selectVarRecursive(
   SCIP*                 scip,
   STATUS*               status,
   PERSISTENTDATA*       persistent,
   CONFIGURATION*        config,
   SCIP_SOL*             baselpsol,
   DOMAINREDUCTIONS*     domainreductions,
   BINCONSDATA*          binconsdata,
   SCIP_VAR**            lpcands,            /**< array of fractional variables */
   SCIP_Real*            lpcandssol,         /**< array of fractional solution values */
   SCIP_Real*            lpcandsfrac,
   int                   nlpcands,           /**< number of fractional variables/solution values */
   BRANCHINGDECISION*    decision,
   int                   recursiondepth
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
   ,LOCALSTATISTICS*     localstats
#endif
   );

static
SCIP_RETCODE executeDownBranchingRecursive(
   SCIP*                 scip,
   STATUS*               status,
   PERSISTENTDATA*       persistent,
   CONFIGURATION*        config,
   SCIP_SOL*             baselpsol,
   SCIP_VAR*             branchvar,
   SCIP_Real             branchval,
   SCIP_Real             branchvalfrac,
   SCIP_Real             localbaselpsolval,
   int                   probingdepth,
   int                   recursiondepth,
   DOMAINREDUCTIONS*     downdomainreductions,
   BINCONSDATA*          binconsdata,
   BRANCHINGRESULTDATA*  downbranchingresult,
   SCIP_Bool*            addeddomainreduction
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
   ,LOCALSTATISTICS*     localstats
#endif
   )
{
   if( config->usebincons && SCIPvarIsBinary(branchvar) )
   {
      /* In case that the branch variable is binary, add the negated var to the list.
       * This list is used to generate a set packing constraint for cutoff branches which were reached by only using
       * binary variables.
       * DownBranching on a binary variable x means: x <= 0
       * When this cutoff occurs we have that: x >= 1 <=> 1-x <= 0
       */
      SCIP_VAR* negbranchvar;

      SCIP_CALL( SCIPgetNegatedVar(scip, branchvar, &negbranchvar) );

      assert(negbranchvar != NULL);

      SCIP_CALL( appendToBinaryVarList(scip, binconsdata->binaryvars, negbranchvar) );
   }

   SCIPdebugMessage("Depth <%i>, Started down branching on var <%s> with var < <%g>\n", probingdepth, SCIPvarGetName(branchvar), branchval);

   SCIP_CALL( executeDownBranching(scip, branchvar, branchval, downbranchingresult, status) );
   SCIPstatistic( statistics->nlpssolved[probingdepth]++; )

   SCIPdebug(
      if( downbranchingresult->cutoff )
      {
         SCIPdebugMessage("Depth <%i>, The solved LP was infeasible and as such cutoff\n", probingdepth);
      }
   )

   if( !downbranchingresult->cutoff && !status->lperror && !status->limitreached && recursiondepth > 1 )
   {
      /* TODO: maybe reuse these variables in the upbranching case */
      SCIP_VAR** deeperlpcands;
      SCIP_Real* deeperlpcandssol;
      SCIP_Real* deeperlpcandsfrac;
      int deepernlpcands;
      SCIP_Real localdowngain;

      localdowngain = MAX(0, downbranchingresult->objval - localbaselpsolval);

      SCIP_CALL( SCIPupdateVarPseudocost(scip, branchvar, 0.0-branchvalfrac, localdowngain, 1.0) );

      SCIP_CALL( copyLPBranchCands(scip, &deeperlpcands, &deeperlpcandssol, &deeperlpcandsfrac, &deepernlpcands) );

      SCIPdebugMessage("Depth <%i>, Downbranching has <%i> candidates.\n", probingdepth, deepernlpcands);

      if( deepernlpcands > 0 )
      {
         BRANCHINGDECISION* deeperdecision;
         STATUS* deeperstatus;
         PERSISTENTDATA* deeperpersistent = NULL;
         SCIP_Real deeperlpobjval = SCIPgetLPObjval(scip);
         SCIPstatistic( LOCALSTATISTICS* deeperlocalstats; )

         SCIP_CALL( allocateStatus(scip, &deeperstatus) );

         SCIP_CALL( allocateBranchingDecision(scip, &deeperdecision, deeperlpobjval) );

         SCIPstatistic( SCIP_CALL( allocateLocalStatistics(scip, &deeperlocalstats) ); )

#ifdef SCIP_STATISTIC
         SCIP_CALL( selectVarRecursive(scip, deeperstatus, deeperpersistent, config, baselpsol, downdomainreductions, binconsdata,
               deeperlpcands, deeperlpcandssol, deeperlpcandsfrac, deepernlpcands, deeperdecision, recursiondepth - 1,
               statistics, deeperlocalstats) );
#else
         SCIP_CALL( selectVarRecursive(scip, deeperstatus, deeperpersistent, config, baselpsol, downdomainreductions, binconsdata,
               deeperlpcands, deeperlpcandssol, deeperlpcandsfrac, deepernlpcands, deeperdecision, recursiondepth - 1) );
#endif

         /* the proved dual bound of the deeper branching cannot be less than the current dual bound, as every deeper
          * node has more/tighter constraints and as such cannot be better than the base LP. */
         assert(SCIPisGE(scip, deeperdecision->provedbound, downbranchingresult->dualbound));
         downbranchingresult->dualbound = deeperdecision->provedbound;
         downbranchingresult->dualboundvalid = TRUE;

         SCIPstatistic(
            if( deeperlocalstats->ndomredproofnodes > 0 )
            {
               localstats->ndomredproofnodes += deeperlocalstats->ndomredproofnodes;
               *addeddomainreduction = TRUE;
            }
         )

         if( deeperstatus->cutoff )
         {
            /* TODO: set another flag? Rework in general
             * upbranchingresult->cutoff is TRUE, if the up child was directly infeasible (so here it is always
             * false, as we don't want to branch on an infeasible node)
             * deeperstatus->cutoff is TRUE, if any up/down child pair of the up child were cutoff
             * */
            downbranchingresult->cutoff = deeperstatus->cutoff;
            SCIPstatistic( localstats->ncutoffproofnodes += deeperlocalstats->ncutoffproofnodes; )

            SCIPdebugMessage("Depth <%i>, Both deeper children were cutoff, so the down branch is cutoff\n", probingdepth);
         }

         SCIPstatistic( freeLocalStatistics(scip, &deeperlocalstats); )
         freeStatus(scip, &deeperstatus);
         freeBranchingDecision(scip, &deeperdecision);
      }
      SCIPfreeBufferArray(scip, &deeperlpcandsfrac);
      SCIPfreeBufferArray(scip, &deeperlpcandssol);
      SCIPfreeBufferArray(scip, &deeperlpcands);
   }

   if( config->usebincons && downbranchingresult->cutoff && binconsdata->binaryvars->nbinaryvars == (probingdepth + 1) )
   {
#ifdef SCIP_STATISTIC
      addBinaryConstraint(scip, config, binconsdata, baselpsol, statistics);
#else
      addBinaryConstraint(scip, config, binconsdata, baselpsol);
#endif
   }

   if( config->usebincons && SCIPvarIsBinary(branchvar) )
   {
      SCIP_VAR* negbranchvar;

      SCIP_CALL( SCIPgetNegatedVar(scip, branchvar, &negbranchvar) );

      assert(negbranchvar != NULL);

#ifdef NDEBUG
      dropFromBinaryVarList(scip, binconsdata->binaryvars);
#else
      {
         SCIP_VAR* droppedelement;
         droppedelement = dropFromBinaryVarList(scip, binconsdata->binaryvars);
         assert(droppedelement == negbranchvar);
      }
#endif

   }

   /* reset the probing depth to undo the previous branching */
   SCIPbacktrackProbing(scip, probingdepth);

   return SCIP_OKAY;
}

static
SCIP_RETCODE executeUpBranchingRecursive(
   SCIP*                 scip,
   STATUS*               status,
   PERSISTENTDATA*       persistent,
   CONFIGURATION*        config,
   SCIP_SOL*             baselpsol,
   SCIP_VAR*             branchvar,
   SCIP_Real             branchval,
   SCIP_Real             branchvalfrac,
   SCIP_Real             localbaselpsolval,
   int                   probingdepth,
   int                   recursiondepth,
   DOMAINREDUCTIONS*     updomainreductions,
   BINCONSDATA*          binconsdata,
   BRANCHINGRESULTDATA*  upbranchingresult,
   SCIP_Bool*            addeddomainreduction
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
   ,LOCALSTATISTICS*     localstats
#endif
   )
{
   if( config->usebincons && SCIPvarIsBinary(branchvar) )
   {
      /* In case that the branch variable is binary, add the var to the list.
       * This list is used to generate a set packing constraint for cutoff branches which were reached by only using
       * binary variables.
       * UpBranching on a binary variable x means: x >= 1
       * When this cutoff occurs we have that: x <= 0
       */
      SCIP_CALL( appendToBinaryVarList(scip, binconsdata->binaryvars, branchvar) );
   }

   SCIPdebugMessage("Depth <%i>, Started up branching on var <%s> with var > <%g>\n", probingdepth, SCIPvarGetName(branchvar), branchval);

   SCIP_CALL( executeUpBranching(scip, branchvar, branchval, upbranchingresult, status) );
   SCIPstatistic( statistics->nlpssolved[probingdepth]++; )

   SCIPdebug(
      if( upbranchingresult->cutoff )
      {
         SCIPdebugMessage("Depth <%i>, The solved LP was infeasible and as such cutoff\n", probingdepth);
      }
   )

   if( !upbranchingresult->cutoff && !status->lperror && !status->limitreached && recursiondepth > 1 )
   {
      SCIP_VAR** deeperlpcands;
      SCIP_Real* deeperlpcandssol;
      SCIP_Real* deeperlpcandsfrac;
      int deepernlpcands;
      SCIP_Real localupgain;

      localupgain = MAX(0, upbranchingresult->objval - localbaselpsolval);

      SCIP_CALL( SCIPupdateVarPseudocost(scip, branchvar, 1.0-branchvalfrac, localupgain, 1.0) );

      SCIP_CALL( copyLPBranchCands(scip, &deeperlpcands, &deeperlpcandssol, &deeperlpcandsfrac, &deepernlpcands) );

      SCIPdebugMessage("Depth <%i>, Upbranching has <%i> candidates.\n", probingdepth, deepernlpcands);

      if( deepernlpcands > 0 )
      {
         BRANCHINGDECISION* deeperdecision;
         STATUS* deeperstatus;
         PERSISTENTDATA* deeperpersistent = NULL;
         SCIP_Real deeperlpobjval = SCIPgetLPObjval(scip);
         SCIPstatistic( LOCALSTATISTICS* deeperlocalstats; )

         SCIP_CALL( allocateStatus(scip, &deeperstatus) );

         SCIP_CALL( allocateBranchingDecision(scip, &deeperdecision, deeperlpobjval) );

         SCIPstatistic( SCIP_CALL( allocateLocalStatistics(scip, &deeperlocalstats) ); )

#ifdef SCIP_STATISTIC
         SCIP_CALL( selectVarRecursive(scip, deeperstatus, deeperpersistent, config, baselpsol, updomainreductions, binconsdata,
               deeperlpcands, deeperlpcandssol, deeperlpcandsfrac, deepernlpcands, deeperdecision, recursiondepth - 1,
               statistics, deeperlocalstats) );
#else
         SCIP_CALL( selectVarRecursive(scip, deeperstatus, deeperpersistent, config, baselpsol, updomainreductions, binconsdata,
               deeperlpcands, deeperlpcandssol, deeperlpcandsfrac, deepernlpcands, deeperdecision, recursiondepth - 1) );
#endif

         /* the proved dual bound of the deeper branching cannot be less than the current dual bound, as every deeper
          * node has more/tighter constraints and as such cannot be better than the base LP. */
         assert(SCIPisGE(scip, deeperdecision->provedbound, upbranchingresult->dualbound));
         upbranchingresult->dualbound = deeperdecision->provedbound;
         upbranchingresult->dualboundvalid = TRUE;

         SCIPstatistic(
            if( deeperlocalstats->ndomredproofnodes > 0 )
            {
               localstats->ndomredproofnodes += deeperlocalstats->ndomredproofnodes;
               *addeddomainreduction = TRUE;
            }
         )

         if( deeperstatus->cutoff )
         {
            /* TODO: set another flag? Rework in general
             * upbranchingresult->cutoff is TRUE, if the up child was directly infeasible (so here it is always
             * false, as we don't want to branch on an infeasible node)
             * deeperstatus->cutoff is TRUE, if any up/down child pair of the up child were cutoff
             * */
            upbranchingresult->cutoff = deeperstatus->cutoff;
            SCIPstatistic(
               localstats->ncutoffproofnodes += deeperlocalstats->ncutoffproofnodes;
            )

            SCIPdebugMessage("Depth <%i>, Both deeper children were cutoff, so the up branch is cutoff\n", probingdepth);
         }

         SCIPstatistic( freeLocalStatistics(scip, &deeperlocalstats); )
         freeStatus(scip, &deeperstatus);
         freeBranchingDecision(scip, &deeperdecision);
      }
      SCIPfreeBufferArray(scip, &deeperlpcandsfrac);
      SCIPfreeBufferArray(scip, &deeperlpcandssol);
      SCIPfreeBufferArray(scip, &deeperlpcands);
   }

   if( config->usebincons && upbranchingresult->cutoff && binconsdata->binaryvars->nbinaryvars == (probingdepth + 1) )
   {
#ifdef SCIP_STATISTIC
      addBinaryConstraint(scip, config, binconsdata, baselpsol, statistics);
#else
      addBinaryConstraint(scip, config, binconsdata, baselpsol);
#endif
   }

   if( config->usebincons && SCIPvarIsBinary(branchvar) )
   {


#ifdef NDEBUG
      dropFromBinaryVarList(scip, binconsdata->binaryvars);
#else
      {
         SCIP_VAR* droppedelement;
         droppedelement = dropFromBinaryVarList(scip, binconsdata->binaryvars);
         assert(droppedelement == branchvar);
      }
#endif
   }

   /* reset the probing depth to undo the previous branching */
   SCIPbacktrackProbing(scip, probingdepth);

   return SCIP_OKAY;
}

static
SCIP_RETCODE selectVarRecursive(
   SCIP*                 scip,
   STATUS*               status,
   PERSISTENTDATA*       persistent,
   CONFIGURATION*        config,
   SCIP_SOL*             baselpsol,
   DOMAINREDUCTIONS*     domainreductions,
   BINCONSDATA*          binconsdata,
   SCIP_VAR**            lpcands,            /**< array of fractional variables */
   SCIP_Real*            lpcandssol,         /**< array of fractional solution values */
   SCIP_Real*            lpcandsfrac,
   int                   nlpcands,           /**< number of fractional variables/solution values */
   BRANCHINGDECISION*    decision,
   int                   recursiondepth
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
   ,LOCALSTATISTICS*     localstats
#endif
   )
{
   int probingdepth;

   assert(scip != NULL);
   assert(status != NULL);
   assert(config != NULL);
   assert(!config->usedomainreduction || domainreductions != NULL);
   assert(!config->usebincons || binconsdata != NULL);
   assert(lpcands != NULL);
   assert(lpcandssol != NULL);
   assert(lpcandsfrac != NULL);
   assert(nlpcands > 0);
   assert(decision != NULL);
   assert(recursiondepth >= 1);
   SCIPstatistic( assert(statistics != NULL) );

   probingdepth = SCIPgetProbingDepth(scip);

   /* init default decision */
   decision->bestvar = lpcands[0];
   decision->bestval = lpcandssol[0];
   decision->bestdown = SCIPgetLPObjval(scip);
   decision->bestdownvalid = FALSE;
   decision->bestup = SCIPgetLPObjval(scip);
   decision->bestupvalid = FALSE;
   decision->provedbound = SCIPgetLPObjval(scip);

   if( !config->forcebranching && nlpcands == 1 )
   {
      SCIPdebugMessage("Only one candidate (<%s>) is given. This one is chosen without calculations.\n", SCIPvarGetName(lpcands[0]));
      SCIPstatistic( statistics->nsinglecandidate[probingdepth]++; )
   }
   else
   {

      if( SCIPgetDepthLimit(scip) <= (SCIPgetDepth(scip) + recursiondepth) )
      {
         /* we need at least 'recursiondepth' space for the branching */
         SCIPdebugMessage("Cannot perform probing in selectVarRecursive, depth limit reached. Current:<%i>, Max:<%i>\n",
            SCIPgetDepthLimit(scip), SCIPgetDepth(scip) + recursiondepth);
         status->depthtoosmall = TRUE;
         SCIPstatistic( statistics->ndepthreached++; )
      }
      else
      {
         int i;
         int c;
         SCIP_Real bestscore = -SCIPinfinity(scip);
         SCIP_Real bestscorelowerbound;
         SCIP_Real bestscoreupperbound;
         SCIP_Real localbaselpsolval = SCIPgetLPObjval(scip);

         bestscorelowerbound = SCIPvarGetLbLocal(decision->bestvar);
         bestscoreupperbound = SCIPvarGetUbLocal(decision->bestvar);

         SCIPdebugMessage("Started selectVarRecursive with depth <%i>/<%i> and <%i> candidates.\n", probingdepth+1, recursiondepth+probingdepth, nlpcands);

         for( i = 0, c = persistent ? persistent->restartindex : 0; i < nlpcands && isExecuteBranchingLoop(status) && !SCIPisStopped(scip); i++, c++ )
         {
            SCIP_Bool useoldbranching = FALSE;
            SCIP_Bool addeddomainreduction = FALSE;
            SCIP_VAR* branchvar;
            SCIP_Real branchval;
            SCIP_Real branchvalfrac;
            BRANCHINGRESULTDATA* downbranchingresult;
            BRANCHINGRESULTDATA* upbranchingresult;

            c = c % nlpcands;

            branchvar = lpcands[c];
            branchval = lpcandssol[c];
            branchvalfrac = lpcandsfrac[c];

            assert(branchvar != NULL);

            /* Reset the cutoffproofnodes, as the number of proof nodes from previous branching vars (which where not
             * cutoff, as we didn't break the loop) is not relevant for the min total sum of proof nodes. */
            SCIPstatistic( localstats->ncutoffproofnodes = 0; )

            SCIP_CALL( SCIPallocBuffer(scip, &downbranchingresult) );
            SCIP_CALL( SCIPallocBuffer(scip, &upbranchingresult) );

            initBranchingResultData(scip, downbranchingresult);
            initBranchingResultData(scip, upbranchingresult);

            if( persistent != NULL && isUseOldBranching(scip, config, persistent, branchvar) ){
               getOldBranching(scip, persistent, branchvar, downbranchingresult, upbranchingresult);
               useoldbranching = TRUE;
               SCIPstatistic( statistics->noldbranchused[probingdepth]++; )
            }
            else
            {
               DOMAINREDUCTIONS* downdomainreductions;
               DOMAINREDUCTIONS* updomainreductions;

               SCIPdebugMessage("Depth <%i>, Started branching on var <%s>\n", probingdepth, SCIPvarGetName(branchvar));

               if( config->usedomainreduction )
               {
                  SCIP_CALL( allocDomainReductions(scip, &downdomainreductions) );
                  SCIP_CALL( allocDomainReductions(scip, &updomainreductions) );
               }

               if( config->downfirst )
               {
#ifdef SCIP_STATISTIC
                  SCIP_CALL( executeDownBranchingRecursive(scip, status, persistent, config, baselpsol, branchvar, branchval,
                     branchvalfrac, localbaselpsolval, probingdepth, recursiondepth, downdomainreductions, binconsdata,
                     downbranchingresult, &addeddomainreduction, statistics, localstats) );

                  SCIP_CALL( executeUpBranchingRecursive(scip, status, persistent, config, baselpsol, branchvar, branchval,
                     branchvalfrac, localbaselpsolval, probingdepth, recursiondepth, updomainreductions, binconsdata,
                     upbranchingresult, &addeddomainreduction, statistics, localstats) );
#else
                  SCIP_CALL( executeDownBranchingRecursive(scip, status, persistent, config, baselpsol, branchvar, branchval,
                     branchvalfrac, localbaselpsolval, probingdepth, recursiondepth, downdomainreductions, binconsdata,
                     downbranchingresult, &addeddomainreduction) );

                  SCIP_CALL( executeUpBranchingRecursive(scip, status, persistent, config, baselpsol, branchvar, branchval,
                     branchvalfrac, localbaselpsolval, probingdepth, recursiondepth, updomainreductions, binconsdata,
                     upbranchingresult, &addeddomainreduction) );
#endif
               }
               else
               {
#ifdef SCIP_STATISTIC
                  SCIP_CALL( executeUpBranchingRecursive(scip, status, persistent, config, baselpsol, branchvar, branchval,
                     branchvalfrac, localbaselpsolval, probingdepth, recursiondepth, updomainreductions, binconsdata,
                     upbranchingresult, &addeddomainreduction, statistics, localstats) );

                  SCIP_CALL( executeDownBranchingRecursive(scip, status, persistent, config, baselpsol, branchvar, branchval,
                     branchvalfrac, localbaselpsolval, probingdepth, recursiondepth, downdomainreductions, binconsdata,
                     downbranchingresult, &addeddomainreduction, statistics, localstats) );
#else
                  SCIP_CALL( executeUpBranchingRecursive(scip, status, persistent, config, baselpsol, branchvar, branchval,
                     branchvalfrac, localbaselpsolval, probingdepth, recursiondepth, updomainreductions, binconsdata,
                     upbranchingresult, &addeddomainreduction) );

                  SCIP_CALL( executeDownBranchingRecursive(scip, status, persistent, config, baselpsol, branchvar, branchval,
                     branchvalfrac, localbaselpsolval, probingdepth, recursiondepth, downdomainreductions, binconsdata,
                     downbranchingresult, &addeddomainreduction) );
#endif
               }

               if( persistent != NULL )
               {
                  updateOldBranching(scip, persistent, branchvar, downbranchingresult, upbranchingresult);
               }

               if( config->usedomainreduction )
               {
                  applyDeeperDomainReductions(scip, baselpsol, domainreductions, downdomainreductions, updomainreductions);

                  freeDomainReductions(scip, &updomainreductions);
                  freeDomainReductions(scip, &downdomainreductions);
               }
            }

            {
               /* TODO: move this block to an own method when finished */

               if( upbranchingresult->cutoff && downbranchingresult->cutoff )
               {
                  SCIPdebugMessage("Depth <%i>, Both branches cutoff\n", probingdepth);

                  /* in a higher level this cutoff may be transferred as a domain reduction/valid bound */
                  status->cutoff = TRUE;
                  SCIPstatistic(
                     statistics->nfullcutoffs[probingdepth]++;
                     localstats->ncutoffproofnodes += 2;
                  )
               }
               else if( upbranchingresult->cutoff )
               {
                  SCIPdebugMessage("Depth <%i>, Up branch cutoff\n", probingdepth);

                  SCIPstatistic(
                     statistics->nsinglecutoffs[probingdepth]++;
                  )

                  if( config->usedomainreduction && !useoldbranching )
                  {
                     addUpperBound(scip, branchvar, branchval, baselpsol, domainreductions);
                     addeddomainreduction = TRUE;
                  }

                  if( downbranchingresult->dualboundvalid )
                  {
                     decision->provedbound = MAX(decision->provedbound, downbranchingresult->dualbound);
                  }
               }
               else if( downbranchingresult->cutoff )
               {
                  SCIPdebugMessage("Depth <%i>, Up branch cutoff\n", probingdepth);

                  SCIPstatistic(
                     statistics->nsinglecutoffs[probingdepth]++;
                  )

                  if( config->usedomainreduction && !useoldbranching )
                  {
                     addLowerBound(scip, branchvar, branchval, baselpsol, domainreductions);
                     addeddomainreduction = TRUE;
                  }
                  if( upbranchingresult->dualboundvalid )
                  {
                     decision->provedbound = MAX(decision->provedbound, upbranchingresult->dualbound);
                  }
               }
               else if( !status->limitreached )
               {
                  SCIP_Real downdualbound = downbranchingresult->dualbound;
                  SCIP_Real updualbound = upbranchingresult->dualbound;

                  /* TODO: currently scoring function is best dual bound. Maybe make it more generic? */
                  SCIP_Real score = MIN(downdualbound, updualbound);

                  SCIPdebugMessage("Depth <%i>, Neither branch is cutoff and no limit reached.\n", probingdepth);

                  if( SCIPisGT(scip, score, bestscore) )
                  {
                     bestscore = score;
                     bestscorelowerbound = SCIPvarGetLbLocal(decision->bestvar);
                     bestscoreupperbound = SCIPvarGetUbLocal(decision->bestvar);

                     decision->bestvar = branchvar;
                     decision->bestval = branchval;
                     decision->bestdown = downdualbound;
                     decision->bestdownvalid = TRUE;
                     decision->bestup = updualbound;
                     decision->bestupvalid = TRUE;
                  }

                  if( upbranchingresult->dualboundvalid && downbranchingresult->dualboundvalid )
                  {
                     decision->provedbound = MAX(decision->provedbound, MIN(updualbound, downdualbound));
                  }
                  else if( upbranchingresult->dualboundvalid )
                  {
                     decision->provedbound = MAX(decision->provedbound, updualbound);
                  }
                  else if( downbranchingresult->dualboundvalid )
                  {
                     decision->provedbound = MAX(decision->provedbound, downdualbound);
                  }

                  if( (config->maxnviolatedcons != -1) && (config->usebincons || config->usedomainreduction) && !useoldbranching )
                  {
                     int nimpliedbincons = 0;
                     int ndomreds = 0;

                     if( config->usebincons )
                     {
                        nimpliedbincons = binconsdata->createdconstraints->nviolatedcons;
                        SCIPdebugMessage("Depth <%i>, Found <%i> violating binary constraints.\n", probingdepth, nimpliedbincons);
                     }

                     if( config->usedomainreduction )
                     {
                        ndomreds = domainreductions->nviolatedvars;
                        SCIPdebugMessage("Depth <%i>, Found <%i> violating bound changes.\n", probingdepth, ndomreds);
                     }

                     if( nimpliedbincons + ndomreds >= config->maxnviolatedcons )
                     {
                        status->maxnconsreached = TRUE;
                     }
                  }
               }

               SCIPstatistic(
                  /* Increment the number of domredproofnodes by one, as we needed the current node as a proof node. */
                  if( addeddomainreduction )
                  {
                     localstats->ndomredproofnodes++;
                  }
               )
            }

            SCIPfreeBuffer(scip, &upbranchingresult);
            SCIPfreeBuffer(scip, &downbranchingresult);

            if( areBoundsChanged(scip, decision->bestvar, bestscorelowerbound, bestscoreupperbound) )
            {
               /* in case the bounds of the current highest scored solution have changed due to domain propagation during the
                * lookahead branching we can/should not branch on this variable but instead return the SCIP_REDUCEDDOM result */
               status->propagationdomred = TRUE;
               SCIPstatistic( statistics->npropdomred[probingdepth]++; )
            }
         }

         if( persistent != NULL )
         {
            persistent->restartindex = c;
         }
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE selectVarStart(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   STATUS*               status,
   BRANCHINGDECISION*    decision
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
   ,LOCALSTATISTICS*     localstats
#endif
   )
{
   CONFIGURATION* config;
   PERSISTENTDATA* persistent;
   int recursiondepth;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   DOMAINREDUCTIONS* domainreductions = NULL;
   BINCONSDATA* binconsdata = NULL;
   SCIP_SOL* baselpsol = NULL;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(status != NULL);
   assert(decision != NULL);
   SCIPstatistic( assert(statistics != NULL); )

   config = branchruledata->config;
   persistent = branchruledata->persistent;
   assert(config != NULL);
   assert(persistent != NULL);

   recursiondepth = config->recursiondepth;
   assert(recursiondepth > 0);

   SCIP_CALL( copyLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands) );

   if( config->usedomainreduction || config->usebincons )
   {
      SCIP_CALL( copyCurrentSolution(scip, &baselpsol) );
   }

   if( config->usedomainreduction )
   {
      SCIP_CALL( allocDomainReductions(scip, &domainreductions) );
   }

   if( config->usebincons )
   {
      SCIP_CALL( allocBinConsData(scip, &binconsdata, recursiondepth, SCIPceil(scip, 0.5*nlpcands)) );
   }

   SCIPstartProbing(scip);
   SCIPenableVarHistory(scip);

#ifdef SCIP_STATISTIC
   SCIP_CALL( selectVarRecursive(scip, status, persistent, config, baselpsol, domainreductions, binconsdata, lpcands, lpcandssol, lpcandsfrac, nlpcands, decision, recursiondepth, statistics, localstats) );
#else
   SCIP_CALL( selectVarRecursive(scip, status, persistent, config, baselpsol, domainreductions, binconsdata, lpcands, lpcandssol, lpcandsfrac, nlpcands, decision, recursiondepth) );
#endif

   SCIPendProbing(scip);

   if( (config->usebincons && binconsdata->createdconstraints->nconstraints > 0)
      || (config->usedomainreduction && domainreductions->nchangedvars > 0) )
   {
      SCIP_CALL( SCIPlinkLPSol(scip, persistent->prevbinsolution) );
      SCIP_CALL( SCIPunlinkSol(scip, persistent->prevbinsolution) );

      copyBranchingDecision(decision, persistent->prevdecision);
   }

   if( config->usebincons )
   {
      SCIP_NODE* basenode;

      assert(binconsdata->binaryvars->nbinaryvars == 0);

      basenode = SCIPgetCurrentNode(scip);

#ifdef SCIP_STATISTIC
      SCIP_CALL( applyBinaryConstraints(scip, basenode, binconsdata->createdconstraints, &status->addbinconst, statistics) );
#else
      SCIP_CALL( applyBinaryConstraints(scip, basenode, binconsdata->createdconstraints, &status->addbinconst) );
#endif
      freeBinConsData(scip, &binconsdata);

   }

   if( config->usedomainreduction )
   {
      if( !status->lperror && !status->depthtoosmall && !status->cutoff )
      {
#ifdef SCIP_STATISTIC
         SCIP_CALL( applyDomainReductions(scip, baselpsol, domainreductions, &status->domredcutoff, &status->domred, statistics) );
#else
         SCIP_CALL( applyDomainReductions(scip, baselpsol, domainreductions, &status->domredcutoff, &status->domred) );
#endif
      }
      freeDomainReductions(scip, &domainreductions);
   }

   if( config->usedomainreduction || config->usebincons )
   {
      SCIP_CALL( SCIPfreeSol(scip, &baselpsol) );
   }

   freeLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac);

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
   PERSISTENTDATA*       persistent
   )
{
   return isBranchingDecisionValid(scip, persistent->prevdecision)
      && SCIPareSolsEqual(scip, currentsol, persistent->prevbinsolution);
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
   BRANCHINGDECISION*    decision,
   SCIP_Bool*            result              /**< the pointer to the branching result */
   )
{
   SCIPdebugMessage("Branching based on previous solution.\n");

   /* execute the actual branching */
   SCIP_CALL( branchOnVar(scip, decision) );
   *result = SCIP_BRANCHED;

   SCIPdebugMessage("Result: Branched based on previous solution. Variable <%s>\n",
      SCIPvarGetName(decision->bestvar));

   /* reset the var pointer, as this is our indicator whether we should branch on prev data in the next call */
   decision->bestvar = NULL;

   return SCIP_OKAY;
}

static
SCIP_RETCODE initBranchruleData(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata
   )
{
   int nvars;
   int i;

   /* Create an empty solution. Gets filled in case of implied binary bounds. */
   SCIP_CALL( SCIPcreateSol(scip, &branchruledata->persistent->prevbinsolution, NULL) );

   /* The variables given by the SCIPgetVars() array are sorted with the binaries at first and the integer variables
    * directly afterwards. With the SCIPvarGetProbindex() method we can access the index of a given variable in the
    * SCIPgetVars() array and as such we can use it to access our arrays which should only contain binary and integer
    * variables.
    */
   nvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);

   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->persistent->lastbranchid, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->persistent->lastbranchnlps, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->persistent->lastbranchupres, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->persistent->lastbranchdownres, nvars) );

   branchruledata->persistent->prevdecision->bestvar = NULL;
   branchruledata->persistent->restartindex = 0;

   for( i = 0; i < nvars; i++ )
   {
      branchruledata->persistent->lastbranchid[i] = -1;
      branchruledata->persistent->lastbranchnlps[i] = 0;
      SCIP_CALL( SCIPallocMemory(scip, &branchruledata->persistent->lastbranchupres[i]) );
      SCIP_CALL( SCIPallocMemory(scip, &branchruledata->persistent->lastbranchdownres[i]) );
   }

   branchruledata->isinitialized = TRUE;

   return SCIP_OKAY;
}

#ifdef SCIP_STATISTIC
static
void printStatistics(
   SCIP*                 scip,
   STATISTICS*           statistics,
   int                   recursiondepth
   )
{
   if( statistics->ntotalresults > 0 )
   {
      int i;

      SCIPinfoMessage(scip, NULL, "Lookahead Branching was called <%i> times.\n", statistics->ntotalresults);
      for( i = 1; i < 18; i++ )
      {
         /* see type_result.h for the id <-> enum mapping */
         SCIPinfoMessage(scip, NULL, "Result <%s> was chosen <%i> times\n", getStatusString(i), statistics->nresults[i]);
      }

      for( i = 0; i < recursiondepth; i++ ) {
         SCIPinfoMessage(scip, NULL, "In depth <%i>, <%i> fullcutoffs and <%i> single cutoffs were found.\n", i, statistics->nfullcutoffs[i], statistics->nsinglecutoffs[i]);
         SCIPinfoMessage(scip, NULL, "In depth <%i>, <%i> LPs were solved.\n", i, statistics->nlpssolved[i]);
         SCIPinfoMessage(scip, NULL, "In depth <%i>, a decision was discarded <%i> times due to domain reduction because of propagation.\n", i, statistics->npropdomred[i]);
         SCIPinfoMessage(scip, NULL, "In depth <%i>, only one branching candidate was given <%i> times.\n", i, statistics->nsinglecandidate[i]);
         SCIPinfoMessage(scip, NULL, "In depth <%i>, old branching results were used in <%i> cases.\n", i, statistics->noldbranchused[i]);
      }

      SCIPinfoMessage(scip, NULL, "Depth limit was reached <%i> times.\n", statistics->ndepthreached);
      SCIPinfoMessage(scip, NULL, "Ignored <%i> binary constraints, that would be domain reductions.\n", statistics->ndomredcons);
      SCIPinfoMessage(scip, NULL, "Added <%i> binary constraints, of which <%i> where violated by the base LP.\n", statistics->nbinconst, statistics->nbinconstvio);
      SCIPinfoMessage(scip, NULL, "Reduced the domain of <%i> vars, <%i> of them where violated by the base LP.\n", statistics->ndomred, statistics->ndomredvio);
      SCIPinfoMessage(scip, NULL, "Needed <%i> additional nodes to proof the cutoffs of base nodes\n", statistics->ncutoffproofnodes);
      SCIPinfoMessage(scip, NULL, "Needed <%i> additional nodes to proof the domain reductions\n", statistics->ndomredproofnodes);
   }
}
#endif

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

   SCIPfreeMemory(scip, &branchruledata->persistent->prevdecision);
   SCIPfreeMemory(scip, &branchruledata->persistent);
   SCIPfreeMemory(scip, &branchruledata->config);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
#ifdef SCIP_STATISTIC
static
SCIP_DECL_BRANCHINIT(branchInitLookahead)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   int recursiondepth;

   branchruledata = SCIPbranchruleGetData(branchrule);
   recursiondepth = branchruledata->config->recursiondepth;

   SCIP_CALL( SCIPallocMemory(scip, &branchruledata->statistics) );
   /* 17 current number of possible result values and the index is 1 based, so 17 + 1 as array size with unused 0 element */
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nresults, 17 + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nsinglecutoffs, recursiondepth) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nfullcutoffs, recursiondepth) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nlpssolved, recursiondepth) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->npropdomred, recursiondepth) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nsinglecandidate, recursiondepth) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->noldbranchused, recursiondepth) );

   initStatistics(scip, branchruledata->statistics, recursiondepth);

   return SCIP_OKAY;
}
#endif

/** deinitialization method of branching rule (called before transformed problem is freed) */
#ifdef SCIP_STATISTIC
static
SCIP_DECL_BRANCHEXIT(branchExitLookahead)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   STATISTICS* statistics;
   int recursiondepth;

   branchruledata = SCIPbranchruleGetData(branchrule);
   statistics = branchruledata->statistics;
   recursiondepth = branchruledata->config->recursiondepth;

   printStatistics(scip, statistics, recursiondepth);

   SCIPfreeMemoryArray(scip, &statistics->noldbranchused);
   SCIPfreeMemoryArray(scip, &statistics->nsinglecandidate);
   SCIPfreeMemoryArray(scip, &statistics->npropdomred);
   SCIPfreeMemoryArray(scip, &statistics->nlpssolved);
   SCIPfreeMemoryArray(scip, &statistics->nfullcutoffs);
   SCIPfreeMemoryArray(scip, &statistics->nsinglecutoffs);
   SCIPfreeMemoryArray(scip, &statistics->nresults);
   SCIPfreeMemory(scip, &statistics);

   return SCIP_OKAY;
}
#endif

/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
static
SCIP_DECL_BRANCHEXITSOL(branchExitSolLookahead)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);

   if( branchruledata->isinitialized )
   {
      int nvars;
      int i;

      nvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);

      for( i = nvars-1; i >= 0; i--)
      {
         SCIPfreeMemory(scip, &branchruledata->persistent->lastbranchdownres[i]);
         SCIPfreeMemory(scip, &branchruledata->persistent->lastbranchupres[i]);
      }

      SCIPfreeMemoryArray(scip, &branchruledata->persistent->lastbranchdownres);
      SCIPfreeMemoryArray(scip, &branchruledata->persistent->lastbranchupres);
      SCIPfreeMemoryArray(scip, &branchruledata->persistent->lastbranchnlps);
      SCIPfreeMemoryArray(scip, &branchruledata->persistent->lastbranchid);

      /* Free the solution that was used for implied binary bounds. */
      SCIP_CALL( SCIPfreeSol(scip, &branchruledata->persistent->prevbinsolution) );

      branchruledata->isinitialized = FALSE;
   }


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

   if( !branchruledata->isinitialized )
   {
      initBranchruleData(scip, branchruledata);
   }

   if( branchruledata->config->usebincons || branchruledata->config->usedomainreduction )
   {
      /* create a copy of the current lp solution to compare it with a previously  */
      SCIP_CALL( copyCurrentSolution(scip, &baselpsol) );
      SCIPdebugMessage("Created an unlinked copy of the base lp solution.\n");
   }

   if( (branchruledata->config->usebincons || branchruledata->config->usedomainreduction)
      && isUsePreviousResult(scip, baselpsol, branchruledata->persistent) )
   {
      /* in case we stopped the previous run without a branching decisions we have stored the decision and execute it
       * now */
      SCIP_CALL( usePreviousResult(scip, branchruledata->persistent->prevdecision, result) );
   }
   else
   {
      SCIP_VAR** lpcands;
      SCIP_Real* lpcandssol;
      SCIP_Real* lpcandsfrac;
      int nlpcands;
      BRANCHINGDECISION* decision;
      STATUS* status;
      SCIP_Real lpobjval = SCIPgetLPObjval(scip);
      SCIPstatistic( LOCALSTATISTICS* localstats; )

      /* get all fractional candidates we can branch on */
      SCIP_CALL( copyLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands) );

      SCIPdebugMessage("The base lp has <%i> variables with fractional value.\n", nlpcands);

      /* create a struct to store the algorithm status */
      SCIP_CALL( allocateStatus(scip, &status) );
      /* create a struct to store the branching decision (in case there is one) */
      SCIP_CALL( allocateBranchingDecision(scip, &decision, lpobjval) );
      SCIPstatistic(
         /* create a struct to store the statistics needed for this single run */
         SCIP_CALL( allocateLocalStatistics(scip, &localstats) );
      )

      /* execute the main logic */
#ifdef SCIP_STATISTIC
      SCIP_CALL( selectVarStart(scip, branchruledata, status, decision, branchruledata->statistics, localstats) );
#else
      SCIP_CALL( selectVarStart(scip, branchruledata, status, decision) );
#endif

      if( status->cutoff || status->domredcutoff )
      {
         *result = SCIP_CUTOFF;
         SCIPstatistic(
            branchruledata->statistics->ncutoffproofnodes += localstats->ncutoffproofnodes;
         )
      }
      else if( status->addbinconst )
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

      SCIPdebugMessage("Result before branching is %s\n", getStatusString(*result));

      if( *result != SCIP_CUTOFF /* a variable could not be branched in any direction or any of the calculated domain
                                  * reductions was infeasible */
         && *result != SCIP_REDUCEDDOM /* the domain of a variable was reduced by evaluating the calculated cutoffs */
         && *result != SCIP_CONSADDED /* implied binary constraints were already added */
         && *result != SCIP_DIDNOTFIND /* an lp error occurred on the way */
         && !status->depthtoosmall /* branching depth wasn't high enough */
         && isBranchingDecisionValid(scip, decision)
         /*&& (0 <= bestcand && bestcand < nlpcands)*/ /* no valid candidate index could be found */
         )
      {
         SCIPinfoMessage(scip, NULL, " -> %d candidates, variable <%s> (solval=%g)\n",
            nlpcands, SCIPvarGetName(decision->bestvar), decision->bestval);

         /* execute the branching as a result of the branching logic */
         SCIP_CALL( branchOnVar(scip, decision) );

         *result = SCIP_BRANCHED;
      }

      SCIPdebugMessage("Result after branching is %s\n", getStatusString(*result));

      SCIPstatistic(
         if( *result == SCIP_BRANCHED )
         {
            SCIPdebugMessage("Result: Finished LookaheadBranching by branching.\n");
         }
         else if( *result == SCIP_REDUCEDDOM )
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
            SCIPABORT();
         }
      )

      SCIPstatistic( freeLocalStatistics(scip, &localstats); )
      freeBranchingDecision(scip, &decision);
      freeStatus(scip, &status);

      freeLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac);
   }

   if( branchruledata->config->usebincons )
   {
      SCIP_CALL( SCIPfreeSol(scip, &baselpsol) );
   }

   SCIPstatistic(
      branchruledata->statistics->ntotalresults++;
      branchruledata->statistics->nresults[*result]++;
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
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata->config) );

   /* needs to be allocated here, such that the previous decision can be filled and reset over multiple runs */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata->persistent) );
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata->persistent->prevdecision) );
   branchruledata->isinitialized = FALSE;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyLookahead) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeLookahead) );
   SCIPstatistic(SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitLookahead) ));
   SCIPstatistic(SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitLookahead) ));
   SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitSolLookahead) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpLookahead) );

   /* add lookahead branching rule parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/useimpliedbincons",
         "should implied binary constraints found via cutoffs on the second level be applied?",
         &branchruledata->config->usebincons, TRUE, DEFAULT_USEIMPLIEDBINARYCONSTRAINTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/addbinconsrow",
         "should implied binary constraints be added as a row to the LP?",
         &branchruledata->config->addbinconsrow, TRUE, DEFAULT_ADDBINCONSROW, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead/maxnumberviolatedcons",
         "how many constraints that are violated by the base lp solution should be gathered until they are added?",
         &branchruledata->config->maxnviolatedcons, TRUE, DEFAULT_MAXNUMBERVIOLATEDCONS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "branching/lookahead/reevalage",
         "number of intermediate LPs solved to trigger reevaluation of strong branching value for a variable that was already evaluated at the current node",
         &branchruledata->config->reevalage, TRUE, DEFAULT_REEVALAGE, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/stopbranching",
         "if the first level down branch is infeasible, should we stop evaluating the up branch?",
         &branchruledata->config->stopbranching, TRUE, DEFAULT_STOPBRANCHING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/purgeconsoncutoff",
         "if the first level down branch is infeasible, should we discard all constraints found for this variable?",
         &branchruledata->config->purgeconsoncutoff, TRUE, DEFAULT_PURGECONSONCUTOFF, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/forcebranching",
         "should lookahead branching be applied even if there is just a single candidate?",
         &branchruledata->config->forcebranching, TRUE, DEFAULT_FORCEBRANCHING, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead/recursiondepth",
         "In case of recursion, how deep should it go?",
         &branchruledata->config->recursiondepth, TRUE, DEFAULT_RECURSIONDEPTH, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/usedomainreduction",
         "should domain reductions found via cutoff be applied (only in recursion)?",
         &branchruledata->config->usedomainreduction, TRUE, DEFAULT_USEDOMAINREDUCTION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/addnonviocons",
         "should constraints be added, that are not violated by the base LP?",
         &branchruledata->config->addnonviocons, TRUE, DEFAULT_ADDNONVIOCONS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
            "branching/lookahead/downbranchfirst",
            "should the down branch be chosen first?",
            &branchruledata->config->downfirst, TRUE, DEFAULT_DOWNFIRST, NULL, NULL) );

   return SCIP_OKAY;
}
