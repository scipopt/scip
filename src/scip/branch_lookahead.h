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

/**@file   branch_lookahead.h
 * @ingroup BRANCHINGRULES
 * @brief  lookahead branching rule TODO CS: expand the description
 * @author Christoph Schubert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_LOOKAHEAD_H__
#define __SCIP_BRANCH_LOOKAHEAD_H__


#include "scip/scip.h"
#include "scip/common_branch_lookahead.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * The parameter that can be changed by the user/caller and alter the behaviour of the lookahead branching.
 */
struct Configuration
{
   SCIP_Bool             usedomainreduction; /**< indicates whether the data for domain reductions should be gathered and
                                              *   used. */
   SCIP_Bool             usebincons;         /**< indicates whether the data for the implied binary constraints should
                                              *   be gathered and used */
   SCIP_Bool             addbinconsrow;      /**< Add the implied binary constraints as a row to the problem matrix */
   SCIP_Bool             stopbranching;      /**< indicates whether we should stop the first level branching after finding
                                              *   an infeasible first branch */
   SCIP_Longint          reevalage;          /**< The number of "normal" (not probing) lps that may have been solved before
                                              *   we stop using old data and start recalculating new first level data. */
   int                   maxnviolatedcons;   /**< The number of constraints we want to gather before restarting the run. Set
                                              *   to -1 for an unbounded list. */
   SCIP_Bool             forcebranching;     /**< Execute the lookahead logic even if only one branching candidate is given.
                                              *   May be used to calculate the score of a single candidate. */
   int                   recursiondepth;     /**< How deep should the recursion go? Default for Lookahead: 2 */
   SCIP_Bool             addnonviocons;      /**< Should constraints be added, that are not violated by the base LP? */
   SCIP_Bool             downfirst;          /**< Should the down branch be executed first? */
   SCIP_Bool             abbreviated;        /**< Should the abbreviated version be used? */
   int                   maxncands;          /**< If abbreviated == TRUE, at most how many candidates should be handled? */
   SCIP_Bool             stopaftercutoff;    /**< Should a branching loop be stopped, if the cutoff of a variable is
                                              *   detected? */
};
typedef struct Configuration CONFIGURATION;

SCIP_RETCODE allocConfiguration(
   SCIP*                 scip,
   CONFIGURATION**       config
   );

void freeConfiguration(
   SCIP*                 scip,
   CONFIGURATION**       config
   );

/**
 * A container to hold the result of a branching.
 */
struct BranchingResultData
{
   SCIP_Real             objval;             /**< The objective value of the solved lp. Only contains meaningful data, if
                                              *   cutoff == FALSE. */
   SCIP_Bool             cutoff;             /**< Indicates whether the node was infeasible and was cutoff. */
   SCIP_Bool             dualboundvalid;     /**< Us the value of the dual bound valid? That means, was the according LP
                                              *   or the sub problems solved to optimality? */
   SCIP_Real             dualbound;          /**< The best dual bound for this branching, may be changed by lower level
                                              *   branchings. */
};
typedef struct BranchingResultData BRANCHINGRESULTDATA;

/**
 * The data that is preserved over multiple runs of the branching rule.
 */
struct PersistentData
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
};
typedef struct PersistentData PERSISTENTDATA;

/**
 * information about the current status of the branching
 */
struct Status
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
};
typedef struct Status STATUS;

SCIP_RETCODE allocateStatus(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS**              status              /**< the status to be allocated */
   );

void freeStatus(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS**              status              /**< the status to be freed */
   );

#ifdef SCIP_STATISTIC
/**
 * The data used for some statistical analysis.
 */
struct Statistics
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
};
typedef struct Statistics STATISTICS;

SCIP_RETCODE allocStatistics(
   SCIP*                 scip,
   STATISTICS**          statistics,
   int                   recursiondepth
   );

void printStatistics(
   SCIP*                 scip,
   STATISTICS*           statistics,
   int                   recursiondepth
   );

void freeStatistics(
   SCIP*                 scip,
   STATISTICS**          statistics
   );
/**
 * Helper struct to store the statistical data needed in a single run.
 */
struct LocalStatistics
{
   int                   ncutoffproofnodes;  /**< The number of nodes needed to proof the current cutoff. */
   int                   ndomredproofnodes;  /**< The number of nodes needed to proof all current domreds. */
};
typedef struct LocalStatistics LOCALSTATISTICS;

SCIP_RETCODE allocateLocalStatistics(
   SCIP*                 scip,
   LOCALSTATISTICS**     localstats
   );

void freeLocalStatistics(
   SCIP*                 scip,
   LOCALSTATISTICS**     localstats
   );

#endif

SCIP_RETCODE selectVarStart(
   SCIP*                 scip,
   CONFIGURATION*        config,
   PERSISTENTDATA*       persistent,
   STATUS*               status,
   BRANCHRULERESULT*     branchruleresult,
   SCIP_Real             lpobjval
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
   ,LOCALSTATISTICS*     localstats
#endif
   );



/** creates the Lookahead branching rule and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeBranchruleLookahead(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
