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
#define SCIP_DEBUG
*/
#define SCIP_STATISTIC

/**@file   branch_lookahead.c
 * @brief  lookahead branching rule
 * @author Christoph Schubert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "lpi/lpi.h"
#include "scip/branch_lookahead.h"
#include "scip/clock.h"
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

#define DEFAULT_USEBINARYCONSTRAINTS         TRUE
#define DEFAULT_USEDOMAINREDUCTION           TRUE
#define DEFAULT_ADDBINCONSROW                FALSE
#define DEFAULT_MAXNUMBERVIOLATEDCONS        10000
#define DEFAULT_REEVALAGE                    10LL
#define DEFAULT_STOPBRANCHING                TRUE
#define DEFAULT_FORCEBRANCHING               FALSE
#define DEFAULT_RECURSIONDEPTH               2
#define DEFAULT_ADDNONVIOCONS                FALSE
#define DEFAULT_DOWNFIRST                    TRUE
#define DEFAULT_ABBREVIATED                  FALSE
#define DEFAULT_MAXNCANDS                    4
#define DEFAULT_STOPAFTERCUTOFF              TRUE
#define DEFAULT_REUSEBASIS                   TRUE

/*
 * Data structures
 */


/** Holds the information needed for branching on a variable. */
typedef struct
{
   SCIP_VAR*             var;                /**< Variable to branch on. May be NULL.*/
   SCIP_Real             val;                /**< the fractional value of the variable to branch on */
   SCIP_Real             downdb;             /**< the highest dual bound found in the down branch */
   SCIP_Bool             downdbvalid;        /**< Indicator for the validity of the downdb value. Is FALSE, if no actual
                                              *   branching occurred or the value was determined by an LP not solved to
                                              *   optimality. */
   SCIP_Real             updb;               /**< the highest dual bound found in the up branch */
   SCIP_Bool             updbvalid;          /**< Indicator for the validity of the updb value. Is FALSE, if no actual
                                              *   branching occurred or the value was determined by an LP not solved to
                                              *   optimality. */
   SCIP_Real             proveddb;           /**< the highest dual bound found with optimally solved LPs */
} BRANCHINGDECISION;

/** Allocates a branching decision in the buffer and initiates it with default values. */
static
SCIP_RETCODE branchingDecisionAllocate(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGDECISION**   decision            /**< pointer to the decision to allocate and initiate */
   )
{
   assert(scip != NULL);
   assert(decision != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, decision) );
   (*decision)->var = NULL;
   (*decision)->downdb = -SCIPinfinity(scip);
   (*decision)->downdbvalid = FALSE;
   (*decision)->updb = -SCIPinfinity(scip);
   (*decision)->updbvalid = FALSE;
   (*decision)->proveddb = -SCIPinfinity(scip);

   return SCIP_OKAY;
}

/** Copies the data from the source to the target. */
static
void branchingDecisionCopy(
   BRANCHINGDECISION*    sourcedecision,     /**< the source branching decision */
   BRANCHINGDECISION*    targetdecision      /**< the target branching decision */
   )
{
   assert(sourcedecision != NULL);
   assert(targetdecision != NULL);

   targetdecision->var = sourcedecision->var;
   targetdecision->val = sourcedecision->val;
   targetdecision->downdb = sourcedecision->downdb;
   targetdecision->downdbvalid = sourcedecision->downdbvalid;
   targetdecision->updb = sourcedecision->updb;
   targetdecision->updbvalid = sourcedecision->updbvalid;
   targetdecision->proveddb = sourcedecision->proveddb;
}

/** Checks whether the given branching decision can be used to branch on. */
static
SCIP_Bool branchingDecisionIsValid(
   BRANCHINGDECISION*    decision            /**< the branching decision to check */
   )
{
   assert(decision != NULL);

   /* a branching decision is deemed valid, if the var point is not on the default NULL value (see the allocate method) */
   return decision->var != NULL;
}

/** Frees the allocated buffer memory of the branching decision. */
static
void branchingDecisionFree(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGDECISION**   decision            /**< pointer to the decision to be freed */
   )
{
   assert(scip != NULL);
   assert(decision != NULL);

   SCIPfreeBuffer(scip, decision);
}

/** A container to hold the result of a branching. */
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

/** Allocates a branching result in the buffer. */
static
SCIP_RETCODE branchingResultDataAllocate(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGRESULTDATA** resultdata          /**< pointer to the result to be allocated */
   )
{
   assert(scip != NULL);
   assert(resultdata != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, resultdata) );
   return SCIP_OKAY;
}

/** Initiates the branching result with default values. */
static
void branchingResultDataInit(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGRESULTDATA*  resultdata          /**< pointer to the result to be initialized */
   )
{
   assert(scip != NULL);
   assert(resultdata != NULL);

   resultdata->objval = -SCIPinfinity(scip);
   resultdata->dualbound = -SCIPinfinity(scip);
   resultdata->cutoff = FALSE;
   resultdata->dualboundvalid = FALSE;
}

/** Copies the data from the source to the target. */
static
void branchingResultDataCopy(
   BRANCHINGRESULTDATA*  sourcedata,         /**< the source branching result */
   BRANCHINGRESULTDATA*  targetdata          /**< the target branching result */
   )
{
   assert(sourcedata != NULL);
   assert(targetdata != NULL);

   targetdata->cutoff = sourcedata->cutoff;
   targetdata->objval = sourcedata->objval;
   targetdata->dualbound = sourcedata->dualbound;
   targetdata->dualboundvalid = sourcedata->dualboundvalid;
}

/** Frees the allocated buffer memory of the branching result. */
static
void branchingResultDataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGRESULTDATA** resultdata          /**< pointer to the result to be freed */
   )
{
   assert(scip != NULL);
   assert(resultdata != NULL);

   SCIPfreeBuffer(scip, resultdata);
}

/** The data that is preserved over multiple runs of the branching rule. */
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

/** The parameter that can be changed by the user/caller and alter the behaviour of the lookahead branching. */
typedef struct
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
   SCIP_Bool             reusebasis;         /**< If abbreviated == TRUE, should the solution lp-basis of the FSB run be
                                              *   used in the first abbreviated level?  */
} CONFIGURATION;

/** Allocates a configuration in the buffer and initiates it with the default values. */
static
SCIP_RETCODE configurationAllocate(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION**       config,             /**< pointer to the configuration to allocate in initialize */
   CONFIGURATION*        copysource          /**< Copy the settings from this config. May be NULL. */
   )
{
   assert(scip != NULL);
   assert(config != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, config) );

   /* if we have a source configuration, copy the data from there, otherwise take the global default values */
   if( copysource != NULL )
   {
      (*config)->addbinconsrow = copysource->addbinconsrow;
      (*config)->addnonviocons = copysource->addnonviocons;
      (*config)->downfirst = copysource->downfirst;
      (*config)->forcebranching = copysource->forcebranching;
      (*config)->maxnviolatedcons = copysource->maxnviolatedcons;
      (*config)->recursiondepth = copysource->recursiondepth;
      (*config)->reevalage = copysource->reevalage;
      (*config)->stopbranching = copysource->stopbranching;
      (*config)->usebincons = copysource->usebincons;
      (*config)->usedomainreduction = copysource->usedomainreduction;
      (*config)->abbreviated = copysource->abbreviated;
      (*config)->maxncands = copysource->maxncands;
      (*config)->stopaftercutoff = copysource->stopaftercutoff;
      (*config)->reusebasis = copysource->reusebasis;
   }
   else
   {
      (*config)->addbinconsrow = DEFAULT_ADDBINCONSROW;
      (*config)->addnonviocons = DEFAULT_ADDNONVIOCONS;
      (*config)->downfirst = DEFAULT_DOWNFIRST;
      (*config)->forcebranching = DEFAULT_FORCEBRANCHING;
      (*config)->maxnviolatedcons = DEFAULT_MAXNUMBERVIOLATEDCONS;
      (*config)->recursiondepth = DEFAULT_RECURSIONDEPTH;
      (*config)->reevalage = DEFAULT_REEVALAGE;
      (*config)->stopbranching = DEFAULT_STOPBRANCHING;
      (*config)->usebincons = DEFAULT_USEBINARYCONSTRAINTS;
      (*config)->usedomainreduction = DEFAULT_USEDOMAINREDUCTION;
      (*config)->abbreviated = DEFAULT_ABBREVIATED;
      (*config)->maxncands = DEFAULT_MAXNCANDS;
      (*config)->stopaftercutoff = DEFAULT_STOPAFTERCUTOFF;
      (*config)->reusebasis = DEFAULT_REUSEBASIS;
   }
   return SCIP_OKAY;
}

/** Frees the allocated buffer memory of the branching result. */
static
void configurationFree(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION**       config              /**< pointer to the configuration to free */
   )
{
   assert(scip != NULL);
   assert(config != NULL);

   SCIPfreeBuffer(scip, config);
}


#ifdef SCIP_STATISTIC
/** The data used for some statistical analysis. */
typedef struct
{
   int                   ntotalresults;      /**< The total sum of the entries in nresults. */
   int*                  nresults;           /**< Array of counters for each result state the lookahead branching finished.
                                              *   The first (0) entry is unused, as the result states are indexed 1-based
                                              *   and we use this index as our array index. */
   int*                  nsinglecutoffs;     /**< The number of single cutoffs on a (probing) node per probingdepth. */
   int*                  nfullcutoffs;       /**< The number of double cutoffs on a (probing) node per probingdepth. */
   int*                  nlpssolved;         /**< The number of all lps solved for a given probingdepth (incl. FSB). */
   int*                  nlpssolvedfsb;      /**< The number of lps solved by the initial FSB to get the FSB scores. */
   SCIP_Longint*         nlpiterations;      /**< The number of all lp iterations needed for a given probingdepth
                                              *   (incl. FSB). */
   SCIP_Longint*         nlpiterationsfsb;   /**< The number of lp iterations needed to get the FSB scores. */
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
   int                   recursiondepth;     /**< The recursiondepth of the LAB. Can be used to access the depth-dependent
                                              *   arrays contained in the statistics. */
} STATISTICS;

/* Initializes the statistics with the start values. */
static
void statisticsInit(
   SCIP*                 scip,               /**< SCIP data structure */
   STATISTICS*           statistics          /**< the statistics to be initialized */
   )
{
   int i;

   assert(scip != NULL);
   assert(statistics != NULL);
   assert(statistics->recursiondepth > 0);

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

   for( i = 0; i < statistics->recursiondepth; i++ )
   {
      statistics->noldbranchused[i] = 0;
      statistics->nsinglecandidate[i] = 0;
      statistics->npropdomred[i] = 0;
      statistics->nfullcutoffs[i] = 0;
      statistics->nlpssolved[i] = 0;
      statistics->nlpssolvedfsb[i] = 0;
      statistics->nlpiterations[i] = 0;
      statistics->nlpiterationsfsb[i] = 0;
      statistics->nsinglecutoffs[i] = 0;
   }
}

/** Allocates a static in the buffer and initiates it with the default values. */
static
SCIP_RETCODE statisticsAllocate(
   SCIP*                 scip,               /**< SCIP data structure */
   STATISTICS**          statistics,         /**< pointer to the statistics to be allocated */
   int                   recursiondepth      /**< the LAB recursion depth */
   )
{

   assert(scip != NULL);
   assert(statistics != NULL);
   assert(recursiondepth > 0);

   SCIP_CALL( SCIPallocBuffer(scip, statistics) );
   /* 17 current number of possible result values and the index is 1 based, so 17 + 1 as array size with unused 0 element */
   SCIP_CALL( SCIPallocBufferArray(scip, &(*statistics)->nresults, 17 + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*statistics)->nsinglecutoffs, recursiondepth) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*statistics)->nfullcutoffs, recursiondepth) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*statistics)->nlpssolved, recursiondepth) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*statistics)->nlpssolvedfsb, recursiondepth) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*statistics)->nlpiterations, recursiondepth) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*statistics)->nlpiterationsfsb, recursiondepth) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*statistics)->npropdomred, recursiondepth) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*statistics)->nsinglecandidate, recursiondepth) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*statistics)->noldbranchused, recursiondepth) );

   (*statistics)->recursiondepth = recursiondepth;

   statisticsInit(scip, *statistics);
   return SCIP_OKAY;
}

/* An array containing all human readable names of the possible results. */
static const char* names[18] = { "", "SCIP_DIDNOTRUN", "SCIP_DELAYED", "SCIP_DIDNOTFIND", "SCIP_FEASIBLE", "SCIP_INFEASIBLE",
   "SCIP_UNBOUNDED", "SCIP_CUTOFF", "SCIP_SEPARATED", "SCIP_NEWROUND", "SCIP_REDUCEDDOM", "SCIP_CONSADDED",
   "SCIP_CONSCHANGED", "SCIP_BRANCHED", "SCIP_SOLVELP", "SCIP_FOUNDSOL", "SCIP_SUSPENDED", "SCIP_SUCCESS" };

/** Returns a human readable name for the given result enum value. */
static
const char* getStatusString(
   SCIP_RESULT           result
)
{
   /* the result can be used as an array index, as it is internally handled as an int value. */
   return names[result];
}

/** Prints the content of the statistics to stdout.  */
static
void statisticsPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   STATISTICS*           statistics          /**< the statistics to print */
   )
{

   assert(scip != NULL);
   assert(statistics != NULL);
   assert(statistics->recursiondepth > 0);

   /* only print something, if we have any statistics */
   if( statistics->ntotalresults > 0 )
   {
      int i;

      SCIPinfoMessage(scip, NULL, "Lookahead Branching was called <%i> times.\n", statistics->ntotalresults);
      for( i = 1; i < 18; i++ )
      {
         /* see type_result.h for the id <-> enum mapping */
         SCIPinfoMessage(scip, NULL, "Result <%s> was chosen <%i> times\n", getStatusString(i), statistics->nresults[i]);
      }

      for( i = 0; i < statistics->recursiondepth; i++ )
      {
         SCIPinfoMessage(scip, NULL, "In depth <%i>, <%i> fullcutoffs and <%i> single cutoffs were found.\n",
            i, statistics->nfullcutoffs[i], statistics->nsinglecutoffs[i]);
         SCIPinfoMessage(scip, NULL, "In depth <%i>, <%i> LPs were solved, <%i> of them to calculate the FSB score.\n",
            i, statistics->nlpssolved[i], statistics->nlpssolvedfsb[i]);
         SCIPinfoMessage(scip, NULL, "In depth <%i>, <%"SCIP_LONGINT_FORMAT"> iterations were needed to solve the LPs, <%"SCIP_LONGINT_FORMAT"> of them to calculate the FSB score.\n",
            i, statistics->nlpiterations[i], statistics->nlpiterationsfsb[i]);
         SCIPinfoMessage(scip, NULL, "In depth <%i>, a decision was discarded <%i> times due to domain reduction because of propagation.\n",
            i, statistics->npropdomred[i]);
         SCIPinfoMessage(scip, NULL, "In depth <%i>, only one branching candidate was given <%i> times.\n",
            i, statistics->nsinglecandidate[i]);
         SCIPinfoMessage(scip, NULL, "In depth <%i>, old branching results were used in <%i> cases.\n",
            i, statistics->noldbranchused[i]);
      }

      SCIPinfoMessage(scip, NULL, "Depth limit was reached <%i> times.\n", statistics->ndepthreached);
      SCIPinfoMessage(scip, NULL, "Ignored <%i> binary constraints, that would be domain reductions.\n",
         statistics->ndomredcons);
      SCIPinfoMessage(scip, NULL, "Added <%i> binary constraints, of which <%i> where violated by the base LP.\n",
         statistics->nbinconst, statistics->nbinconstvio);
      SCIPinfoMessage(scip, NULL, "Reduced the domain of <%i> vars, <%i> of them where violated by the base LP.\n",
         statistics->ndomred, statistics->ndomredvio);
      SCIPinfoMessage(scip, NULL, "Needed <%i> additional nodes to proof the cutoffs of base nodes\n",
         statistics->ncutoffproofnodes);
      SCIPinfoMessage(scip, NULL, "Needed <%i> additional nodes to proof the domain reductions\n",
         statistics->ndomredproofnodes);
   }
}

/** Frees the allocated buffer memory of the statistics. */
static
void statisticsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   STATISTICS**          statistics          /**< pointer to the statistics to free */
   )
{

   assert(scip != NULL);
   assert(statistics != NULL);

   SCIPfreeBufferArray(scip, &(*statistics)->noldbranchused);
   SCIPfreeBufferArray(scip, &(*statistics)->nsinglecandidate);
   SCIPfreeBufferArray(scip, &(*statistics)->npropdomred);
   SCIPfreeBufferArray(scip, &(*statistics)->nlpiterationsfsb);
   SCIPfreeBufferArray(scip, &(*statistics)->nlpiterations);
   SCIPfreeBufferArray(scip, &(*statistics)->nlpssolvedfsb);
   SCIPfreeBufferArray(scip, &(*statistics)->nlpssolved);
   SCIPfreeBufferArray(scip, &(*statistics)->nfullcutoffs);
   SCIPfreeBufferArray(scip, &(*statistics)->nsinglecutoffs);
   SCIPfreeBufferArray(scip, &(*statistics)->nresults);
   SCIPfreeBuffer(scip, statistics);
}

/** Helper struct to store the statistical data needed in a single run. */
typedef struct
{
   int                   ncutoffproofnodes;  /**< The number of nodes needed to proof the current cutoff. */
   int                   ndomredproofnodes;  /**< The number of nodes needed to proof all current domreds. */
} LOCALSTATISTICS;

/** Allocates the local statistics in buffer memory and initializes it with default values. */
static
SCIP_RETCODE localStatisticsAllocate(
   SCIP*                 scip,               /**< SCIP data structure */
   LOCALSTATISTICS**     localstats          /**< pointer to the local statistics to allocate and initialize */
   )
{
   assert(scip != NULL);
   assert(localstats != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, localstats) );

   (*localstats)->ncutoffproofnodes = 0;
   (*localstats)->ndomredproofnodes = 0;

   return SCIP_OKAY;
}

/** Frees the allocated buffer memory of the local statistics. */
static
void localStatisticsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   LOCALSTATISTICS**     localstats          /**< pointer to the local statistics to be freed */
   )
{
   assert(scip != NULL);
   assert(localstats != NULL);

   SCIPfreeBuffer(scip, localstats);
}
#endif

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Bool             isinitialized;      /**< indicates whether the fields in this struct are initialized */
   CONFIGURATION*        config;             /**< the parameter that influence the behaviour of the lookahead branching */
   PERSISTENTDATA*       persistent;         /**< the data that persists over multiple branching decisions */
#ifdef SCIP_STATISTIC
   STATISTICS*           statistics;         /**< statistical data container */
#endif
};

/** all constraints that were created and may be added to the base node */
typedef struct
{
   SCIP_CONS**           constraints;        /**< The array of constraints. Length is adjusted as needed. */
   int                   nconstraints;       /**< The number of entries in the array 'constraints'. */
   int                   memorysize;         /**< The number of entries that the array 'constraints' may hold before the
                                              *   array is reallocated. */
   int                   nviolatedcons;      /**< Tracks the number of constraints, that are violated by the base LP
                                              *   solution. */
} CONSTRAINTLIST;

/** Allocate and initialize the list holding the constraints. */
static
SCIP_RETCODE constraintListAllocate(
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

/** Append an element to the end of the list of constraints. */
static
SCIP_RETCODE constraintListAppend(
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

/** Free all resources of a constraint list in opposite order to the allocation. */
static
void constraintListFree(
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

/** Allocates and initializes the BINARYVARLIST struct. */
static
SCIP_RETCODE binaryVarListAllocate(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYVARLIST**       list,               /**< Pointer to the list to be allocated and initialized. */
   int                   startsize           /**< The number of entries the list initially can hold. */
   )
{
   assert(scip != NULL);
   assert(list != NULL);
   assert(startsize > 0);

   SCIP_CALL( SCIPallocBuffer(scip, list) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*list)->binaryvars, startsize) );

   /* We start with no entries and the (current) max length */
   (*list)->nbinaryvars = 0;
   (*list)->memorysize = startsize;

   return SCIP_OKAY;
}

/** Appends a binary variable to the list, reallocating the list if necessary. */
static
SCIP_RETCODE binaryVarListAppend(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYVARLIST*        list,               /**< The list to add the var to. */
   SCIP_VAR*             vartoadd            /**< The binary var to add to the list. */
   )
{
   assert(scip != NULL);
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

/** Remove and return the last element from the list. */
static
SCIP_VAR* binaryVarListDrop(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYVARLIST*        list                /**< The list to remove the last element from. */
   )
{
   SCIP_VAR* lastelement;

   assert(scip != NULL);
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

/** Frees all resources allocated by a BINARYVARLIST in opposite order of allocation. */
static
void binaryVarListFree(
   SCIP*                 scip,               /**< SCIP data structure */
   BINARYVARLIST**       list                /**< Pointer to the list to free */
   )
{
   assert(scip != NULL);
   assert(list != NULL);

   SCIPfreeBufferArray(scip, &(*list)->binaryvars);
   SCIPfreeBuffer(scip, list);
}

/** struct holding the relevant data for handling binary constraints */
typedef struct
{
   BINARYVARLIST*        binaryvars;         /**< The current binary vars, used to created the constraints. */
   CONSTRAINTLIST*       createdconstraints; /**< The created constraints. */
} BINCONSDATA;

/** Allocate and initialize the BINCONSDATA struct. */
static
SCIP_RETCODE binConsDataAllocate(
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
   SCIP_CALL( binaryVarListAllocate(scip, &(*consdata)->binaryvars, maxdepth) );
   SCIP_CALL( constraintListAllocate(scip, &(*consdata)->createdconstraints, nstartcons) );

   return SCIP_OKAY;
}

/** Free all resources in a BINCONSDATA in opposite order of allocation. */
static
void binConsDataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   BINCONSDATA**         consdata            /**< Pointer to he struct to be freed. */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   constraintListFree(scip, &(*consdata)->createdconstraints);
   binaryVarListFree(scip, &(*consdata)->binaryvars);
   SCIPfreeBuffer(scip, consdata);
}

/** A struct holding information to speed up the solving time of the same problem. This is filled by the FSB scoring routine
 *  that is run to get the best candidates. This is read on the first level of the ALAB routine. */
typedef struct
{
   SCIP_LPISTATE*        lpistate;           /**< the basis information that may be set before another solve lp call */
   SCIP_LPINORMS*        lpinorms;           /**< the norms that may be set before another solve lp call */
   SCIP_Bool             primalfeas;         /**< indicates whether the solution was primal feasible */
   SCIP_Bool             dualfeas;           /**< indicates whether the solution was dual feasible */
} LPIMEMORY;

/** Allocates the lpi memory on the buffer and initializes it with default values. */
static
SCIP_RETCODE lpiMemoryAllocate(
   SCIP*                 scip,               /**< SCIP data structure */
   LPIMEMORY**           lpimemory           /**< the lpi memory to allocate and initialize */
   )
{
   assert(scip != NULL);
   assert(lpimemory != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, lpimemory) );

   (*lpimemory)->lpistate = NULL;
   (*lpimemory)->lpinorms = NULL;
   (*lpimemory)->primalfeas = FALSE;
   (*lpimemory)->dualfeas = FALSE;

   return SCIP_OKAY;
}

/** Copies the fields from the source to the target lpi memory. The content of the fields is not copied, that means
 *  afterwards the following holds true: source->lpistate == target->lpistate */
static
void lpiMemoryShallowCopy(
   LPIMEMORY*            source,             /**< the source struct for the copy procedure */
   LPIMEMORY*            target              /**< the target struct for the copy procedure */
   )
{
   assert(source != NULL);
   assert(target != NULL);

   target->lpinorms = source->lpinorms;
   target->lpistate = source->lpistate;
   target->primalfeas = source->primalfeas;
   target->dualfeas = source->dualfeas;
}

/** Checks that the lpi memory can be read into the lp solver. */
static
SCIP_Bool lpiMemoryIsReadable(
   LPIMEMORY*            memory,             /**< The lpi memory to check. May be NULL. */
   CONFIGURATION*        config              /**< the configuration containing the info whether we should use the data */
   )
{
   assert(config != NULL);
   assert(config->reusebasis || memory == NULL);

   /* the lpinorms may be NULL */
   return memory != NULL && memory->lpistate != NULL && config->reusebasis;
}

/** Checks that the lpi memory can be written to with previous data from the lp solver. */
static
SCIP_Bool lpiMemoryIsWritable(
   LPIMEMORY*            memory,             /**< The lpi memory to check. May be NULL. */
   CONFIGURATION*        config              /**< the configuration containing the info whether we should write the data */
   )
{
   assert(config->reusebasis || memory == NULL);

   /* the lpinorms may be NULL */
   return memory != NULL && config->reusebasis;
}

/** Frees the data contained in the given lpi memory. */
static
SCIP_RETCODE lpiMemoryClear(
   SCIP*                 scip,               /**< SCIP data structure */
   LPIMEMORY*            lpimemory           /**< the lpi memory to clear */
   )
{
   SCIP_LPI* lpi;
   BMS_BLKMEM* blkmem;

   assert(scip != NULL);
   assert(lpimemory != NULL);

   SCIP_CALL( SCIPgetLPI(scip, &lpi) );
   blkmem = SCIPblkmem(scip);

   if( lpimemory->lpistate != NULL )
   {
      SCIP_CALL( SCIPlpiFreeState(lpi, blkmem, &lpimemory->lpistate) );
      lpimemory->lpistate = NULL;
   }

   if( lpimemory->lpinorms != NULL )
   {
      SCIP_CALL( SCIPlpiFreeNorms(lpi, blkmem, &lpimemory->lpinorms) );
      lpimemory->lpinorms = NULL;
   }

   return SCIP_OKAY;
}

/** Frees the allocated buffer memory of the lpi memory without freeing the contained data. */
static
void lpiMemoryFree(
   SCIP*                 scip,               /**< SCIP data structure */
   LPIMEMORY**           lpimemory           /**< the lpi memory to free */
   )
{
   assert(scip != NULL);
   assert(lpimemory != NULL);

   SCIPfreeBuffer(scip, lpimemory);
}

/** A struct containing all information needed to branch on a variable. */
typedef struct
{
   SCIP_VAR*             branchvar;          /**< the variable to branch on */
   SCIP_Real             branchval;          /**< the fractional value to branch on */
   SCIP_Real             fracval;            /**< the fractional part of the value to branch on (val - floor(val)) */
   LPIMEMORY*            downlpimemory;      /**< the lpi memory containing the lp data from a previous down branch */
   LPIMEMORY*            uplpimemory;        /**< the lpi memory containing the lp data from a previous up branch */
} CANDIDATE;

/** Allocates the candidate on the buffer and initializes it with default values. */
static
SCIP_RETCODE candidateAllocate(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATE**           candidate           /**< the candidate to allocate and initialize */
   )
{
   assert(scip != NULL);
   assert(candidate != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, candidate) );
   SCIP_CALL( lpiMemoryAllocate(scip, &(*candidate)->downlpimemory) );
   SCIP_CALL( lpiMemoryAllocate(scip, &(*candidate)->uplpimemory) );

   (*candidate)->branchvar = NULL;

   return SCIP_OKAY;
}

/** Frees the allocated buffer memory of the candidate and clears the contained lpi memories. */
static
SCIP_RETCODE candidateFree(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATE**           candidate           /**< the candidate to free */
   )
{
   assert(scip != NULL);
   assert(candidate != NULL);

   /* if a candidate is freed, we no longer need the content of the lpi memory */
   SCIP_CALL( lpiMemoryClear(scip, (*candidate)->uplpimemory) );
   SCIP_CALL( lpiMemoryClear(scip, (*candidate)->downlpimemory) );

   lpiMemoryFree(scip, &(*candidate)->uplpimemory);
   lpiMemoryFree(scip, &(*candidate)->downlpimemory);
   SCIPfreeBuffer(scip, candidate);
   return SCIP_OKAY;
}

/** A struct acting as a fixed list of candidates */
typedef struct
{
   CANDIDATE**           candidates;         /**< the array of candidates */
   int                   ncandidates;        /**< the number of entries in candidates */
} CANDIDATELIST;

/** Allocates the candidate list on the buffer WITHOUT initializing the contained array of candidates. */
static
SCIP_RETCODE candidateListAllocate(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATELIST**       list                /**< the candidate list to allocate */
   )
{
   assert(scip != NULL);
   assert(list != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, list) );
   (*list)->candidates = NULL;
   (*list)->ncandidates = 0;
   return SCIP_OKAY;
}

/** Initializes the candidate list by allocating the contained array and the correct number of candidates. */
static
SCIP_RETCODE candidateListInit(
   SCIP*                 scip,               /**< SCIP data structure */
   CANDIDATELIST*        list,               /**< the candidate list to initialize */
   int                   ncandidates         /**< the number of candidates the list must hold */
   )
{
   int i;

   assert(scip != NULL);
   assert(list != NULL);
   assert(ncandidates > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &list->candidates, ncandidates) );
   list->ncandidates = ncandidates;

   for( i = 0; i < ncandidates; i++ )
   {
      SCIP_CALL( candidateAllocate(scip, &list->candidates[i]) );
   }

   return SCIP_OKAY;
}

/** Frees the allocated buffer memory of the candidate list and frees the contained candidates. */
static
SCIP_RETCODE candidateListFree(
   SCIP*                 scip,
   CANDIDATELIST**       list
   )
{
   int i;

   assert(scip != NULL);
   assert(list != NULL);

   for( i = (*list)->ncandidates-1; i >= 0; i-- )
   {
      SCIP_CALL( candidateFree(scip, &(*list)->candidates[i]) );
   }

   SCIPfreeBufferArray(scip, &(*list)->candidates);
   SCIPfreeBuffer(scip, list);
   return SCIP_OKAY;
}

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
 * allocate the struct on the Buffer and initialize it with the default values
 */
static
SCIP_RETCODE domainReductionsAllocate(
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

/**
 * frees the given DOMAINREDUCTIONS and all contained Arrays in the opposite order of allocation
 */
static
void domainReductionsFree(
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

/** Allocates the status on the buffer memory and initializes it with default values. */
static
SCIP_RETCODE statusAllocate(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS**              status              /**< the status to be allocated */
   )
{
   assert(scip != NULL);
   assert(status != NULL);

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

/** Frees the allocated buffer memory of the status. */
static
void statusFree(
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS**              status              /**< the status to be freed */
   )
{
   assert(scip != NULL);
   assert(status != NULL);
   SCIPfreeBuffer(scip, status);
}

/** A struct containing all data a LAB run has gathered. */
typedef struct
{
   BRANCHINGDECISION*    decision;           /**< the decision the LAB run made */
   SCIP_VAR**            candswithscore;     /**< The branching variables. Content may be NULL. */
   SCIP_Real*            candscores;         /**< the scores of the candidates */
   SCIP_Real*            candslpvalue;       /**< the branching value of the candidates */
   SCIP_Real*            candsvalfrac;       /**< the faractional value of the branching values */
   LPIMEMORY**           downlpimemories;    /**< the down lpi memory for each candidate */
   LPIMEMORY**           uplpimemories;      /**< the up lpi memory for each candidate */
   int                   ncandscores;        /**< the number of candidates */
} BRANCHRULERESULT;

/** Allocates the branchrule result without all the memory fields.  */
static
SCIP_RETCODE branchRuleResultAllocateReduced(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHRULERESULT**    branchruleresult    /**< the branchrule result to allocate */
)
{
   assert(scip != NULL);
   assert(branchruleresult != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, branchruleresult) );
   SCIP_CALL( branchingDecisionAllocate(scip, &(*branchruleresult)->decision) );
   (*branchruleresult)->candswithscore = NULL;
   (*branchruleresult)->candscores = NULL;
   (*branchruleresult)->candslpvalue = NULL;
   (*branchruleresult)->candsvalfrac = NULL;
   (*branchruleresult)->downlpimemories = NULL;
   (*branchruleresult)->uplpimemories = NULL;
   (*branchruleresult)->ncandscores = 0;

   return SCIP_OKAY;
}

/** Allocates the branchrule result with all the memory fields. */
static
SCIP_RETCODE branchRuleResultAllocateFull(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHRULERESULT**    branchruleresult,   /**< the branchrule result to allocate */
   int                   ncands              /**< the number of candidates that must be hold */
   )
{
   int i;

   assert(scip != NULL);
   assert(branchruleresult != NULL);
   assert(ncands > 0);

   SCIP_CALL( branchRuleResultAllocateReduced(scip, branchruleresult) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*branchruleresult)->downlpimemories, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*branchruleresult)->uplpimemories, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*branchruleresult)->candscores, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*branchruleresult)->candswithscore, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*branchruleresult)->candslpvalue, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*branchruleresult)->candsvalfrac, ncands) );
   (*branchruleresult)->ncandscores = ncands;

   for( i = 0; i < ncands; i++ )
   {
      (*branchruleresult)->candscores[i] = -1;
      (*branchruleresult)->candswithscore[i] = NULL;
      SCIP_CALL( lpiMemoryAllocate(scip, &(*branchruleresult)->downlpimemories[i]) );
      SCIP_CALL( lpiMemoryAllocate(scip, &(*branchruleresult)->uplpimemories[i]) );
   }

   return SCIP_OKAY;
}

/** Frees the buffer memory that was allocated for a reduced branchrule result. */
static
void branchRuleResultFreeReduced(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHRULERESULT**    branchruleresult    /**< the branchrule result to free */
   )
{
   assert(scip != NULL);
   assert(branchruleresult != NULL);

   branchingDecisionFree(scip, &(*branchruleresult)->decision);
   SCIPfreeBuffer(scip, branchruleresult);
}

/** Frees the buffer memory that was allocated for a full branchrule result. */
static
SCIP_RETCODE branchRuleResultFreeFull(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHRULERESULT**    branchruleresult    /**< the branchrule result to free */
   )
{
   int i;

   assert(scip != NULL);
   assert(branchruleresult != NULL);

   assert(&(*branchruleresult)->candslpvalue != NULL);
   assert(&(*branchruleresult)->candswithscore != NULL);
   assert(&(*branchruleresult)->candscores != NULL);

   for( i = (*branchruleresult)->ncandscores-1; i >= 0; i-- )
   {
      lpiMemoryFree(scip, &(*branchruleresult)->uplpimemories[i]);
      lpiMemoryFree(scip, &(*branchruleresult)->downlpimemories[i]);
   }

   SCIPfreeBufferArray(scip, &(*branchruleresult)->candsvalfrac);
   SCIPfreeBufferArray(scip, &(*branchruleresult)->candslpvalue);
   SCIPfreeBufferArray(scip, &(*branchruleresult)->candswithscore);
   SCIPfreeBufferArray(scip, &(*branchruleresult)->candscores);
   SCIPfreeBufferArray(scip, &(*branchruleresult)->uplpimemories);
   SCIPfreeBufferArray(scip, &(*branchruleresult)->downlpimemories);
   branchRuleResultFreeReduced(scip, branchruleresult);

   return SCIP_OKAY;
}

/** Container struct to keep the calculated score for each variable. */
typedef struct
{
   SCIP_Real*            scores;             /**< the scores for each problem variable */
} SCORECONTAINER;

static
SCIP_RETCODE scoreContainerAllocate(
   SCIP*                 scip,
   SCORECONTAINER**      scorecontainer
   )
{
   int ntotalvars;
   int i;

   assert(scip != NULL);
   assert(scorecontainer != NULL);

   /* The container saves the score for all variables in the problem via the ProbIndex. See SCIPvarGetProbindex() */
   ntotalvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBuffer(scip, scorecontainer) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*scorecontainer)->scores, ntotalvars) );

   for( i = 0; i < ntotalvars; i++ )
   {
      (*scorecontainer)->scores[i] = -1;
   }

   return SCIP_OKAY;
}

static
void scoreContainerSetScore(
   SCORECONTAINER*       scorecontainer,
   SCIP_VAR*             var,
   SCIP_Real             score
   )
{
   int probindex;

   assert(scorecontainer != NULL);
   assert(var != NULL);
   assert(score >= 0);

   probindex = SCIPvarGetProbindex(var);

   assert(scorecontainer->scores[probindex] == -1);

   scorecontainer->scores[probindex] = score;
}

static
void scoreContainerFree(
   SCIP*                 scip,
   SCORECONTAINER**      scorecontainer
   )
{
   SCIPfreeBufferArray(scip, &(*scorecontainer)->scores);
   SCIPfreeBuffer(scip, scorecontainer);
}


/*
 * Local methods for the logic
 */
static
SCIP_RETCODE selectVarRecursive(
   SCIP*                 scip,
   STATUS*               status,
   PERSISTENTDATA*       persistent,
   CONFIGURATION*        config,
   SCIP_SOL*             baselpsol,
   DOMAINREDUCTIONS*     domainreductions,
   BINCONSDATA*          binconsdata,
   CANDIDATELIST*        candidates,
   BRANCHRULERESULT*     branchruleresult,
   int                   recursiondepth,
   SCIP_Real             lpobjval
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
   ,LOCALSTATISTICS*     localstats
#endif
);

static
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

/** Adds the given lower bound to the container struct. */
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
   SCIPstatistic(assert(nproofnodes >= 0));

   /* The arrays inside DOMAINREDUCTIONS are indexed via the problem index. */
   varindex = SCIPvarGetProbindex(var);

   /* If we have an old lower bound we take the stronger one, so the MAX is taken. Otherwise we use the new lower bound
    * directly. */
   if( domainreductions->lowerboundset[varindex] )
   {
      if( SCIPisGE(scip, domainreductions->lowerbounds[varindex], lowerbound) )
      {
         /* the old lower bound is stronger (greater) than or equal to the given one, so we keep the older */
         newlowerbound = domainreductions->lowerbounds[varindex];
         SCIPstatistic(
            if( SCIPisEQ(scip, domainreductions->lowerbounds[varindex], lowerbound) )
            {
               /* if the given lower bound is equal to the old one we take the smaller number of proof nodes */
               newnproof = MIN(domainreductions->lowerboundnproofs[varindex], nproofnodes);
            }
            else
            {
               /* if the old bound is stronger (greater) we keep the old number of proof nodes */
               newnproof = domainreductions->lowerboundnproofs[varindex];
            }
         )
      }
      else
      {
         /* the new lower bound is stronger (greater) than the old one, so we update the bound and number of proof nodes */
         newlowerbound = lowerbound;
         SCIPstatistic( newnproof = nproofnodes; )
      }
   }
   else
   {
      /* we have no old lower bound, so we can directly take the new one together with the number of proof nodes */
      newlowerbound = lowerbound;
      SCIPstatistic( newnproof = nproofnodes; )

      if( !domainreductions->upperboundset[varindex] )
      {
         /* if we had no lower and upper bound set we now have a changed variable */
         domainreductions->nchangedvars++;
      }
   }

   /* set the calculated values */
   domainreductions->lowerbounds[varindex] = newlowerbound;
   domainreductions->lowerboundset[varindex] = TRUE;
   SCIPstatistic( domainreductions->lowerboundnproofs[varindex] = newnproof; )

   /* We get the solution value to check whether the domain reduction is violated in the base LP */
   basesolutionval = SCIPgetSolVal(scip, baselpsol, var);

   /* In case the new lower bound is greater than the base solution val and the base solution val is not violated by a
    * previously found bound, we increment the nviolatedvars counter and set the baselpviolated flag.  */
   if( SCIPisGT(scip, newlowerbound, basesolutionval) && !domainreductions->baselpviolated[varindex] )
   {
      domainreductions->baselpviolated[varindex] = TRUE;
      domainreductions->nviolatedvars++;
   }
}

/** Add a lower bound to the DOMAINREDUCTIONS struct. This is used as a wrapper to the 'addLowerBoundProofNode' method. */
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
   /* We add the lower bound with number of proof nodes 2, as this method is only called from the recursion directly. There it
    * is called in case that only one child node is cutoff. As proof nodes we count the cutoff node as well as the valid
    * node, as we need both to proof, that this domain reduction is valid. */
#ifdef SCIP_STATISTIC
   addLowerBoundProofNode(scip, var, lowerbound, baselpsol, domainreductions, 2);
#else
   addLowerBoundProofNode(scip, var, lowerbound, baselpsol, domainreductions);
#endif
}

/** Adds the given upper bound to the container struct. */
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
   SCIPstatistic(assert(nproofnodes >= 0));

   /* The arrays inside DOMAINREDUCTIONS are indexed via the problem index. */
   varindex = SCIPvarGetProbindex(var);

   /* If we have an old upper bound we take the stronger one, so the MIN is taken. Otherwise we use the new upper bound
    * directly. */
   if( domainreductions->upperboundset[varindex] )
   {
      if( SCIPisLE(scip, domainreductions->upperbounds[varindex], upperbound) )
      {
         /* the old upper bound is stronger (smaller) than or equal to the given one, so we keep the older */
         newupperbound = domainreductions->upperbounds[varindex];
         SCIPstatistic(
            if( SCIPisEQ(scip, domainreductions->upperbounds[varindex], upperbound) )
            {
               /* if the given upper bound is equal to the old one we take the smaller number of proof nodes */
               newnproof = MIN(domainreductions->upperboundnproofs[varindex], nproofnodes);
            }
            else
            {
               /* if the old bound is stronger (smaller) we keep the old number of proof nodes */
               newnproof = domainreductions->upperboundnproofs[varindex];
            }
         )
      }
      else
      {
         /* the new upper bound is stronger (smaller) than the old one, so we update the bound and number of proof nodes */
         newupperbound = upperbound;
         SCIPstatistic( newnproof = nproofnodes; )
      }
   }
   else
   {
      /* we have no old upper bound, so we can directly take the new one together with the number of proof nodes */
      newupperbound = upperbound;
      SCIPstatistic( newnproof = nproofnodes; )

      if( !domainreductions->lowerboundset[varindex] )
      {
         /* if we had no lower and upper bound set we now have a changed variable */
         domainreductions->nchangedvars++;
      }
   }

   /* set the calculated values */
   domainreductions->upperbounds[varindex] = newupperbound;
   domainreductions->upperboundset[varindex] = TRUE;
   SCIPstatistic( domainreductions->upperboundnproofs[varindex] = newnproof; )

   /* We get the solution value to check whether the domain reduction is violated in the base LP */
   basesolutionval = SCIPgetSolVal(scip, baselpsol, var);

   /* In case the new upper bound is lesser than the base solution val and the base solution val is not violated by a
    * previously found bound, we increment the nviolatedvars counter and set the baselpviolated flag.  */
   if( SCIPisLT(scip, newupperbound, basesolutionval) && !domainreductions->baselpviolated[varindex] )
   {
      domainreductions->baselpviolated[varindex] = TRUE;
      domainreductions->nviolatedvars++;
   }
}

/** Add a lower bound to the DOMAINREDUCTIONS struct. This is used as a wrapper to the 'addUpperBoundProofNode' method. */
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
#ifdef SCIP_STATISTIC
         int newnproofs;
#endif

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
#ifdef SCIP_STATISTIC
         int newnproofs;
#endif

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

/** Applies the domain reductions to the current node. */
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
         /* if we have a lower bound for the current var, apply it */
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

/** Copies the current LP solution into the given pointer. Needs to be freed after usage! */
static
SCIP_RETCODE copyCurrentSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            lpsol               /**< pointer to store the solution into */
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

/** Executes the branching on a given variable with a given value. */
static
SCIP_RETCODE branchOnVar(
   SCIP*                 scip                /**< SCIP data structure */,
   BRANCHINGDECISION*    decision            /**< the decision with all the needed data */
   )
{
   SCIP_VAR* bestvar = decision->var;
   SCIP_Real bestval = decision->val;

   SCIP_NODE* downchild = NULL;
   SCIP_NODE* upchild = NULL;

   SCIPdebugMessage("Effective branching on var <%s> with value <%g>. Old domain: [%g..%g].\n",
      SCIPvarGetName(bestvar), bestval, SCIPvarGetLbLocal(bestvar), SCIPvarGetUbLocal(bestvar));

   assert(!SCIPisIntegral(scip, bestval));

   /* branch on the given variable */
   SCIP_CALL( SCIPbranchVarVal(scip, bestvar, bestval, &downchild, NULL, &upchild) );

   assert(downchild != NULL);
   assert(upchild != NULL);

   /* update the lower bounds in the children; we must not do this if columns are missing in the LP
    * (e.g. because we are doing branch-and-price) or the problem should be solved exactly */
   if( SCIPallColsInLP(scip) && !SCIPisExactSolve(scip) )
   {
      SCIP_Real bestdown = decision->downdb;
      SCIP_Bool bestdownvalid = decision->downdbvalid;
      SCIP_Real bestup = decision->updb;
      SCIP_Bool bestupvalid = decision->updbvalid;
      SCIP_Real provedbound = decision->proveddb;

      /* update the lower bound for the LPs for further children of both created nodes */
      SCIP_CALL( SCIPupdateNodeLowerbound(scip, downchild, bestdownvalid ? MAX(bestdown, provedbound) : provedbound) );
      SCIP_CALL( SCIPupdateNodeLowerbound(scip, upchild, bestupvalid ? MAX(bestup, provedbound) : provedbound) );
   }
   SCIPdebugMessage(" -> down child's lowerbound: %g\n", SCIPnodeGetLowerbound(downchild));
   SCIPdebugMessage(" -> up child's lowerbound: %g\n", SCIPnodeGetLowerbound(upchild));

   return SCIP_OKAY;
}

/** Store the current lp solution in the lpi memory for further usage. */
static
SCIP_RETCODE storeInLPIMemory(
   SCIP*                 scip,               /**< SCIP data structure */
   LPIMEMORY*            lpimemory           /**< the lpi memory in which the data should be stored */
   )
{
   SCIP_LPI* lpi;
   BMS_BLKMEM* blkmem;

   assert(scip != NULL);
   assert(lpimemory != NULL);
   assert(lpimemory->lpistate == NULL);
   assert(lpimemory->lpinorms == NULL);

   SCIP_CALL( SCIPgetLPI(scip, &lpi) );
   blkmem = SCIPblkmem(scip);

   SCIP_CALL( SCIPlpiGetState(lpi, blkmem, &lpimemory->lpistate) );
   SCIP_CALL( SCIPlpiGetNorms(lpi, blkmem, &lpimemory->lpinorms) );
   lpimemory->primalfeas = SCIPlpiIsPrimalFeasible(lpi);
   lpimemory->dualfeas = SCIPlpiIsDualFeasible(lpi);

   assert(lpimemory->lpistate != NULL);
   /* the norms may be NULL */

   return SCIP_OKAY;
}

/** Sets the lp state and norms of the current node to the values stored in the lpi memory. */
static
SCIP_RETCODE restoreFromLPIMemory(
   SCIP*                 scip,               /**< SCIP data structure */
   LPIMEMORY*            lpimemory           /**< the lpi memory containing the stored data */
   )
{
   assert(scip != NULL);
   assert(lpimemory != NULL);
   assert(lpimemory->lpistate != NULL);

   /* as we solved the very same LP some time earlier and stored the state (the basis) and the norms we can now set those in
    * the LP solver, such that the solution does not (in best case) need any further calculation.
    * Some iterations may occur, as the conflict analysis may have added some constraints in the meantime. */
   SCIP_CALL( SCIPsetProbingLPState(scip, &(lpimemory->lpistate), &(lpimemory->lpinorms), lpimemory->primalfeas,
      lpimemory->dualfeas) );

   /* The state and norms will be freed later by the SCIP framework. Therefore they are set to NULL to enforce that I won't
    * free them on my own. */
   assert(lpimemory->lpistate == NULL);
   assert(lpimemory->lpinorms == NULL);

   return SCIP_OKAY;
}

/** Creates a new probing node with a new bound for the given candidate. */
static
SCIP_RETCODE executeBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,             /**< configuration containing flags changing the behavior */
   SCIP_Bool             downbranching,      /**< the branching direction */
   CANDIDATE*            candidate,          /**< the candidate to branch on */
   BRANCHINGRESULTDATA*  resultdata,         /**< pointer to the result data which gets filled with the status */
   STATUS*               status              /**< status will contain updated lperror and limit fields */
   )
{
   SCIP_Real oldupperbound;
   SCIP_Real oldlowerbound;
   SCIP_Real newbound;
   SCIP_LPSOLSTAT solstat;
   SCIP_VAR* branchvar;
   SCIP_Real branchval;

   assert(scip != NULL);
   assert(config != NULL);
   assert(candidate != NULL);
   assert(resultdata != NULL);
   assert(status != NULL);

   branchvar = candidate->branchvar;
   branchval = candidate->branchval;

   assert(branchvar != NULL);
   assert(!SCIPisFeasIntegral(scip, branchval));

   if( downbranching )
   {
      /* round the given value down, so that it can be used as the new upper bound */
      newbound = SCIPfeasFloor(scip, branchval);
   }
   else
   {
      /* round the given value up, so that it can be used as the new lower bound */
      newbound = SCIPfeasCeil(scip, branchval);
   }

   oldupperbound = SCIPvarGetUbLocal(branchvar);
   oldlowerbound = SCIPvarGetLbLocal(branchvar);

#ifdef SCIP_DEBUG
   if( downbranching )
   {
      SCIPdebugMessage("DownBranching: Var=<%s>, Proposed upper bound=<%g>, old bounds=[<%g>..<%g>], new bounds=[<%g>..<%g>]\n",
         SCIPvarGetName(branchvar), newbound, oldlowerbound, oldupperbound, oldlowerbound, newbound);
   }
   else
   {
      SCIPdebugMessage("UpBranching: Var=<%s>, Proposed lower bound=<%g>, old bounds=[<%g>..<%g>], new bounds=[<%g>..<%g>]\n",
         SCIPvarGetName(branchvar), newbound, oldlowerbound, oldupperbound, newbound, newbound);
   }
#endif

   if( (downbranching && SCIPisFeasLT(scip, newbound, oldlowerbound)) || (!downbranching && SCIPisFeasGT(scip, newbound, oldupperbound)) )
   {
      /* if lb > ub we can cutoff this node */
      resultdata->cutoff = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );

      if( downbranching )
      {
         /* down branching preparations */
         if( SCIPisFeasLT(scip, newbound, oldupperbound) )
         {
            /* if the new upper bound is lesser than the old upper bound and also
             * greater than (or equal to) the old lower bound we set the new upper bound.
             * oldLowerBound <= newUpperBound < oldUpperBound */
            SCIP_CALL( SCIPchgVarUbProbing(scip, branchvar, newbound) );
         }

         if( lpiMemoryIsReadable(candidate->downlpimemory, config) )
         {
            /* restore the stored LP data (e.g. the basis) from a previous run */
            SCIP_CALL( restoreFromLPIMemory(scip, candidate->downlpimemory) );
         }
      }
      else
      {
         /* up branching preparations */
         if( SCIPisFeasGT(scip, newbound, oldlowerbound) )
         {
            /* if the new lower bound is greater than the old lower bound and also
             * lesser than (or equal to) the old upper bound we set the new lower bound.
             * oldLowerBound < newLowerBound <= oldUpperBound */
            SCIP_CALL( SCIPchgVarLbProbing(scip, branchvar, newbound) );
         }

         if( lpiMemoryIsReadable(candidate->uplpimemory, config) )
         {
            /* restore the stored LP data (e.g. the basis) from a previous run */
            SCIP_CALL( restoreFromLPIMemory(scip, candidate->uplpimemory) );
         }
      }

      /* solve the prepared probing LP */
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
         resultdata->dualbound = SCIPgetLPObjval(scip);
         resultdata->dualboundvalid = TRUE;

         resultdata->cutoff = resultdata->cutoff || SCIPisGE(scip, resultdata->objval, SCIPgetCutoffbound(scip));
         assert(solstat != SCIP_LPSOLSTAT_INFEASIBLE || resultdata->cutoff);
      }
   }

   return SCIP_OKAY;
}

/** Executes the branching on the current probing node by adding a probing node with a new upper bound. */
static
SCIP_RETCODE executeDownBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,             /**< configuration containing flags changing the behavior */
   CANDIDATE*            candidate,          /**< the branching direction */
   BRANCHINGRESULTDATA*  resultdata,         /**< pointer to the result data which gets filled with the status */
   STATUS*               status              /**< status will contain updated lperror and limit fields */
   )
{
   SCIP_CALL( executeBranching(scip, config, TRUE, candidate, resultdata, status) );
   return SCIP_OKAY;
}

/** Executes the branching on the current probing node by adding a probing node with a new lower bound. */
static
SCIP_RETCODE executeUpBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,             /**< configuration containing flags changing the behavior */
   CANDIDATE*            candidate,          /**< the branching direction */
   BRANCHINGRESULTDATA*  resultdata,         /**< pointer to the result data which gets filled with the status */
   STATUS*               status              /**< status will contain updated lperror and limit fields */
   )
{
   SCIP_CALL( executeBranching(scip, config, FALSE, candidate, resultdata, status) );
   return SCIP_OKAY;
}

/** Creates a logic or constraint based on the given 'consvars'. This array has to consist of the negated
 * versions of the variables present on a cutoff "path" (path means all variables from the root directly
 * to the cutoff node).
 * Let x_1, ..., x_n be the variables on a path to a cutoff with the branchings x_i <= 1 for all i.
 * Summed up the constraints would look like x_1 + ... x_n <= n-1.
 * Let y_i = 1 - x_i. Then we have y_1 + ... + y_n >= 1 which is a logic or constraint.  */
static
SCIP_RETCODE createBinaryConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   CONFIGURATION*        config,             /**< configuration containing flags changing the behavior */
   SCIP_CONS**           constraint,         /**< Pointer to store the created constraint in */
   char*                 constraintname,     /**< name of the new constraint */
   SCIP_VAR**            consvars,           /**< array containing the negated binary vars */
   int                   nconsvars           /**< the number of elements in 'consvars' */
   )
{
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce = FALSE;
   SCIP_Bool check = FALSE;
   SCIP_Bool propagate = TRUE;
   SCIP_Bool local = TRUE;
   SCIP_Bool modifiable = FALSE;
   SCIP_Bool dynamic = FALSE;
   SCIP_Bool removable = FALSE;
   SCIP_Bool stickingatnode = FALSE;

   assert(scip != NULL);
   assert(config != NULL);
   assert(constraint != NULL);
   assert(constraintname != NULL);
   assert(consvars != NULL);
   assert(nconsvars > 0);

   initial = config->addbinconsrow;
   separate = config->addbinconsrow;

   /* creating a logic or constraint based on the list of vars in 'consvars'.
    * A logic or constraints looks like that: y_1 + ... + y_n >= 1. */
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
   SCIP_VAR**            binaryvars,         /**< the variables contained in the constraint */
   int                   nbinaryvars,        /**< the number of elements in 'binaryvars' */
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

      /* TODO: probably rework */
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
         SCIP_VAR* var = binconsdata->binaryvars->binaryvars[i];

         assert(SCIPvarIsBinary(var));

         SCIP_CALL( SCIPgetNegatedVar(scip, var, &negatedvars[i]) );
         lhssum += SCIPgetSolVal(scip, baselpsol, negatedvars[i]);
      }

      if( config->addnonviocons || lhssum < 1 )
      {
         /* create a name for the new constraint */
         createBinaryConstraintName(scip, negatedvars, binconsdata->binaryvars->nbinaryvars, constraintname);
         /* create the constraint with the freshly created name */
         SCIP_CALL( createBinaryConstraint(scip, config, &constraint, constraintname, negatedvars, binconsdata->binaryvars->nbinaryvars) );

#ifdef PRINTNODECONS
         SCIPinfoMessage(scip, NULL, "Created constraint:\n");
         SCIP_CALL( SCIPprintCons(scip, constraint, NULL) );
         SCIPinfoMessage(scip, NULL, "\n");
#endif

         SCIP_CALL( constraintListAppend(scip, binconsdata->createdconstraints, constraint) );
         /* the constraint we are building is a logic or: we have a list of binary variables that were
          * cutoff while we branched on with >= 1. So we have the constraint: x_1 + ... + x_n <= n-1.
          * Let y = (1-x), then we have an equivalent formulation: y_1 + ... + y_n >= 1. If the base lp
          * is violating this constraint we count this for our number of violated constraints and bounds. */
         if( lhssum < 1 )
         {
            binconsdata->createdconstraints->nviolatedcons++;
         }
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
   //SCIPdebugMessage("Bound change? Var: %s, SCIP LB: %g, My LB: %g, SCIP UB: %g, My UB: %g\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var), lowerbound, SCIPvarGetUbLocal(var), upperbound);
   return SCIPvarGetLbLocal(var) != lowerbound || SCIPvarGetUbLocal(var) != upperbound;
}

static
SCIP_Bool isExecuteBranchingLoop(
   STATUS*               status,
   CONFIGURATION*        config
   )
{
   return !status->lperror && (!config->stopaftercutoff || !status->cutoff) && !status->limitreached && !status->maxnconsreached
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

   branchingResultDataCopy(persistent->lastbranchdownres[varindex], downbranchingresult);
   branchingResultDataCopy(persistent->lastbranchupres[varindex], upbranchingresult);
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

   branchingResultDataCopy(downbranchingresult, persistent->lastbranchdownres[varindex]);
   branchingResultDataCopy(upbranchingresult, persistent->lastbranchupres[varindex]);

   persistent->lastbranchid[varindex] = SCIPgetNNodes(scip);
   persistent->lastbranchnlps[varindex] = SCIPgetNNodeLPs(scip);
}

#ifdef SCIP_STATISTIC
static
void mergeFSBStatistics(
   STATISTICS*           parentstatistics,
   STATISTICS*           fsbstatistics
   )
{
   int i;

   assert(parentstatistics->recursiondepth == fsbstatistics->recursiondepth);

   for( i = 0; i < parentstatistics->recursiondepth; i++ )
   {
      parentstatistics->nlpssolved[i] += fsbstatistics->nlpssolved[i];
      parentstatistics->nlpssolvedfsb[i] += fsbstatistics->nlpssolved[i];
      parentstatistics->nlpiterations[i] += fsbstatistics->nlpiterations[i];
      parentstatistics->nlpiterationsfsb[i] += fsbstatistics->nlpiterations[i];
   }
}
#endif

static
SCIP_RETCODE getFSBResult(
   SCIP*                 scip,
   CONFIGURATION*        parentconfig,
   SCIP_Real             lpobjval,
   BRANCHRULERESULT*     branchruleresult,
   SCORECONTAINER*       scorecontainer
#ifdef SCIP_STATISTIC
   ,STATISTICS*          parentstatistics
#endif
   )
{
   CONFIGURATION* config;
   STATUS* status;
   STATISTICS* statistics;
   LOCALSTATISTICS* localstats;

   SCIP_CALL( configurationAllocate(scip, &config, parentconfig) );
   SCIP_CALL( statusAllocate(scip, &status) );
   /* we need to allocate enough space for all possible depth, as there is currently a problem with setting the fsb stats.
    * E.G.: we want to gather statistics for the FSB run started on layer 2 (1-indexed). Then we start in probing depth 1
    * and solve the nodes in depth 2. */
   SCIP_CALL( statisticsAllocate(scip, &statistics, parentconfig->recursiondepth) );
   SCIP_CALL( localStatisticsAllocate(scip, &localstats) );

   /* We don't want any constraints to be added or domains to be reduced, as we are just interested in the score. */
   config->usebincons = FALSE;
   config->usedomainreduction = FALSE;
   /* Simple FSB is achieved by starting LAB with a max depth of 1 */
   config->recursiondepth = 1;
   /* We want the FSB score over all candidates, otherwise we get end in an endless loop. */
   config->abbreviated = FALSE;
   /* We want to get the scores for all variables, so we don't want to stop after a cutoff is found. */
   config->stopaftercutoff = FALSE;
   /* Even for one single candidate we want to get the corresponding values */
   config->forcebranching = TRUE;

   SCIP_CALL( selectVarStart(scip, config, NULL, status, branchruleresult, scorecontainer, lpobjval, statistics, localstats) );

   mergeFSBStatistics(parentstatistics, statistics);

   localStatisticsFree(scip, &localstats);
   statisticsFree(scip, &statistics);
   statusFree(scip, &status);
   configurationFree(scip, &config);

   return SCIP_OKAY;
}

/** Comparator used to order the candidates in a BRANCHRULERESULT by their FSB score. */
static
SCIP_DECL_SORTINDCOMP(branchRuleScoreComp)
{  /*lint --e{715}*/
   BRANCHRULERESULT* branchruleresult = (BRANCHRULERESULT*)dataptr;
   SCIP_Real score1;
   SCIP_Real score2;

   assert(branchruleresult != NULL);
   assert(0 <= ind1 && ind1 < branchruleresult->ncandscores);
   assert(0 <= ind2 && ind2 < branchruleresult->ncandscores);

   score1 = branchruleresult->candscores[ind1];
   score2 = branchruleresult->candscores[ind2];

   if( score1 == score2 )
   {
      return 0;
   }
   else if( score1 < score2 )
   {
      return -1;
   }
   else
   {
      return 1;
   }
}

static
SCIP_RETCODE getBestCandidates(
   SCIP*                 scip,
   CONFIGURATION*        config,
   CANDIDATELIST*        candidates,
   SCORECONTAINER*       scorecontainer
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
#endif
   )
{
   BRANCHRULERESULT* branchruleresult;
   SCIP_Real lpobjval;
   int nallcands;
   int nusedcands;
   int* permutation;
   int i;

   assert(scip != NULL);
   assert(candidates != NULL);

   SCIP_CALL( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, &nallcands, NULL, NULL) );
   assert(nallcands > 0);

   SCIPdebugMessage("All nlpcands: %i\n", nallcands);
   nusedcands = MIN(config->maxncands, nallcands);
   SCIPdebugMessage("Used nlpcands: %i\n", nusedcands);

   SCIP_CALL( candidateListInit(scip, candidates, nusedcands) );

   /* allocate the result struct to collect the FSB scores of all branching candidates */
   SCIP_CALL( branchRuleResultAllocateFull(scip, &branchruleresult, nallcands) );

   assert(branchruleresult->ncandscores == nallcands);

   lpobjval = SCIPgetLPObjval(scip);

   /* TODO: filter the "candidates" based on the presence of a score in the 'scorecontainer'.
    * Only those without a score need a new one. */
   {
      int i;
      SCIP_VAR** lpcands;
      SCIP_Real* lpcandssol;
      SCIP_Real* lpcandsfrac;
      int nlpcands;

      SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL, NULL) );

      for( i = 0; i < nlpcands; i++ )
      {
         SCIP_VAR* lpcand = lpcands[i];
         int probindex = SCIPvarGetProbindex(lpcand);
         int knownscore scorecontainer->scores[probindex];

         /* TODO: should I also store the lpimemory? */

         if( knownscore == -1 )
         {
            /* score is unknown and needs to be calculated */
         }
         else
         {
            /* score is known */
         }
      }
   }


   SCIPdebugMessage("Calculating the FSB result to get a score for all candidates.\n")
   /* Calculate all FSB scores and collect it in the result */;
   SCIP_CALL( getFSBResult(scip, config, lpobjval, branchruleresult, statistics) );

   /* allocate the permutation array that will contain the sorted indices (permutation[i] will contain the i-th index) */
   SCIP_CALL( SCIPallocBufferArray(scip, &permutation, branchruleresult->ncandscores) );

   SCIPdebugMessage("Sorting the candidates w.r.t. their FSB score.\n");
   /* sort the branching candidates according to their FSB score contained in the result struct */
   SCIPsortDown(permutation, branchRuleScoreComp, branchruleresult, branchruleresult->ncandscores);

   SCIPdebug(
      SCIPdebugMessage("All %i candidates, sorted by their FSB score:\n", branchruleresult->ncandscores);
      for( i = 0; i < branchruleresult->ncandscores; i++ )
      {
         int sortedindex = permutation[i];
         SCIP_VAR* var = branchruleresult->candswithscore[sortedindex];
         SCIP_Real score = branchruleresult->candscores[sortedindex];

         assert(var != NULL);

         SCIPdebugMessage(" Index %2i: Var %s Score %g\n", i, SCIPvarGetName(var), score);
      }
   )

   SCIPdebugMessage("Best %i candidates used for the lookahead branching:\n", candidates->ncandidates);
   for( i = 0; i < nusedcands; i++)
   {
      CANDIDATE* candidate = candidates->candidates[i];
      int sortedindex = permutation[i];
      int probindex;
      SCIP_Real score;

      assert(branchruleresult->candswithscore[sortedindex] != NULL);
      assert(candidate != NULL);

      /* set the branching information, such that we can re-use it in the lookahead branching */
      candidate->branchvar = branchruleresult->candswithscore[sortedindex];
      candidate->branchval = branchruleresult->candslpvalue[sortedindex];
      candidate->fracval = branchruleresult->candsvalfrac[sortedindex];

      if( config->reusebasis )
      {
         lpiMemoryShallowCopy(branchruleresult->downlpimemories[sortedindex], candidate->downlpimemory);
         lpiMemoryShallowCopy(branchruleresult->uplpimemories[sortedindex], candidate->uplpimemory);
      }

      SCIPdebugMessage(" Index %2i: Var %s Val %g\n", i, SCIPvarGetName(candidate->branchvar), candidate->branchval);
   }

   if( config->reusebasis )
   {
      for( i = candidates->ncandidates; i < nallcands; i++ )
      {
         int sortedindex = permutation[i];

         SCIP_CALL( lpiMemoryClear(scip, branchruleresult->downlpimemories[sortedindex]) );
         SCIP_CALL( lpiMemoryClear(scip, branchruleresult->uplpimemories[sortedindex]) );
      }
   }

   /* free the allocated resources */
   SCIPfreeBufferArray(scip, &permutation);
   SCIP_CALL( branchRuleResultFreeFull(scip, &branchruleresult) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE getAllCandidates(
   SCIP*                 scip,
   CANDIDATELIST*        candidates
   )
{
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   int i;

   assert(scip != NULL);
   assert(candidates != NULL);

   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL, NULL) );

   assert(lpcands != NULL);
   assert(lpcandssol != NULL);
   assert(lpcandsfrac != NULL);
   assert(nlpcands > 0);

   SCIP_CALL( candidateListInit(scip, candidates, nlpcands) );

   assert(nlpcands == candidates->ncandidates);

   for( i = 0; i < nlpcands; i++ )
   {
      CANDIDATE* candidate = candidates->candidates[i];

      candidate->branchvar = lpcands[i];
      candidate->branchval = lpcandssol[i];
      candidate->fracval = lpcandsfrac[i];
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE getCandidates(
   SCIP*                 scip,
   CONFIGURATION*        config,
   CANDIDATELIST*        candidates
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
#endif
   )
{
   /* we want to take the abbreviated list of candidates only on the first probing level */
   if( config->abbreviated && (!SCIPinProbing(scip) || SCIPgetProbingDepth(scip) == 0) )
   {
      /* call LAB with depth 1 to get the best (w.r.t. FSB score) candidates */
      SCIPdebugMessage("Getting the branching candidates by selecting the best %i candidates\n", candidates->ncandidates);
      SCIP_CALL( getBestCandidates(scip, config, candidates, statistics) );
   }
   else
   {
      /* get all candidates for the current node lp solution */
      SCIPdebugMessage("Getting the branching candidates by selecting all candidates.\n");
      SCIP_CALL( getAllCandidates(scip, candidates) );
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE executeDownBranchingRecursive(
   SCIP*                 scip,
   STATUS*               status,
   PERSISTENTDATA*       persistent,
   CONFIGURATION*        config,
   SCIP_SOL*             baselpsol,
   CANDIDATE*            candidate,
   SCIP_Real             localbaselpsolval,
   int                   probingdepth,
   int                   recursiondepth,
   DOMAINREDUCTIONS*     downdomainreductions,
   BINCONSDATA*          binconsdata,
   BRANCHINGRESULTDATA*  downbranchingresult,
   SCORECONTAINER*       scorecontainer,
   LPIMEMORY*            lpimemory,
   SCIP_Bool*            addeddomainreduction
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
   ,LOCALSTATISTICS*     localstats
#endif
   )
{
   SCIP_VAR* branchvar = candidate->branchvar;
   SCIP_Real branchval = candidate->branchval;
   SCIP_Real branchvalfrac = candidate->fracval;
   SCIP_Longint niterations;

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

      SCIP_CALL( binaryVarListAppend(scip, binconsdata->binaryvars, negbranchvar) );
   }

   SCIPdebugMessage("Depth <%i>, Started down branching on var <%s> with var < <%g>\n", probingdepth, SCIPvarGetName(branchvar), branchval);

   niterations = SCIPgetNLPIterations(scip);
   SCIP_CALL( executeDownBranching(scip, config, candidate, downbranchingresult, status) );
   SCIPstatistic(
      statistics->nlpssolved[probingdepth]++;
      statistics->nlpiterations[probingdepth] += SCIPgetNLPIterations(scip)-niterations;
   )

   SCIPdebug(
      if( downbranchingresult->cutoff )
      {
         SCIPdebugMessage("Depth <%i>, The solved LP was infeasible and as such cutoff\n", probingdepth);
      }
   )

   if( !downbranchingresult->cutoff && !status->lperror && !status->limitreached )
   {
      int deepernlpcands;
      SCIP_Real localdowngain;

      if( lpiMemoryIsWritable(lpimemory, config) )
      {
         SCIP_CALL( storeInLPIMemory(scip, lpimemory) );
         SCIPdebugMessage("Stored smth in LPIMemory after downbranching on var %s\n", SCIPvarGetName(branchvar));
      }

      localdowngain = MAX(0, downbranchingresult->objval - localbaselpsolval);

      SCIP_CALL( SCIPupdateVarPseudocost(scip, branchvar, 0.0-branchvalfrac, localdowngain, 1.0) );

      if( recursiondepth > 1 )
      {
         SCIP_CALL( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, &deepernlpcands, NULL, NULL) );

         SCIPdebugMessage("Depth <%i>, Downbranching has <%i> candidates.\n", probingdepth, deepernlpcands);

         if( deepernlpcands > 0 )
         {
            CANDIDATELIST* candidates;
            BRANCHRULERESULT* deeperbranchruleresult;
            STATUS* deeperstatus;
            PERSISTENTDATA* deeperpersistent = NULL;
            SCIP_Real deeperlpobjval = downbranchingresult->objval;
#ifdef SCIP_STATISTIC
            LOCALSTATISTICS* deeperlocalstats;
#endif

            SCIP_CALL( candidateListAllocate(scip, &candidates) );

#ifdef SCIP_STATISTIC
            SCIP_CALL( getCandidates(scip, config, candidates, statistics) );
#else
            SCIP_CALL( getCandidates(scip, config, candidates) );
#endif

            SCIP_CALL( statusAllocate(scip, &deeperstatus) );

            SCIP_CALL( branchRuleResultAllocateReduced(scip, &deeperbranchruleresult) );

#ifdef SCIP_STATISTIC
            SCIP_CALL( localStatisticsAllocate(scip, &deeperlocalstats) );
            SCIP_CALL( selectVarRecursive(scip, deeperstatus, deeperpersistent, config, baselpsol, downdomainreductions, binconsdata,
                  candidates, deeperbranchruleresult, scorecontainer, recursiondepth - 1, deeperlpobjval,
                  statistics, deeperlocalstats) );
#else
            SCIP_CALL( selectVarRecursive(scip, deeperstatus, deeperpersistent, config, baselpsol, downdomainreductions, binconsdata,
                  candidates, deeperbranchruleresult, scorecontainer, recursiondepth - 1, deeperlpobjval) );
#endif

            /* the proved dual bound of the deeper branching cannot be less than the current dual bound, as every deeper
             * node has more/tighter constraints and as such cannot be better than the base LP. */
            assert(SCIPisGE(scip, deeperbranchruleresult->decision->proveddb, downbranchingresult->dualbound));
            downbranchingresult->dualbound = deeperbranchruleresult->decision->proveddb;
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
               /* upbranchingresult->cutoff is TRUE, if the up child was directly infeasible (so here it is always
                * false, as we don't want to branch on an infeasible node)
                * deeperstatus->cutoff is TRUE, if any up/down child pair of the up child were cutoff
                * */
               downbranchingresult->cutoff = deeperstatus->cutoff;
               SCIPstatistic( localstats->ncutoffproofnodes += deeperlocalstats->ncutoffproofnodes; )

               SCIPdebugMessage("Depth <%i>, Both deeper children were cutoff, so the down branch is cutoff\n", probingdepth);
            }

#ifdef SCIP_STATISTIC
            localStatisticsFree(scip, &deeperlocalstats);
#endif
            statusFree(scip, &deeperstatus);
            branchRuleResultFreeReduced(scip, &deeperbranchruleresult);
            SCIP_CALL( candidateListFree(scip, &candidates) );
         }
      }
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
      binaryVarListDrop(scip, binconsdata->binaryvars);
#else
      {
         SCIP_VAR* droppedelement;
         droppedelement = binaryVarListDrop(scip, binconsdata->binaryvars);
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
   SCIP*                 scip,               /**< SCIP data structure */
   STATUS*               status,
   PERSISTENTDATA*       persistent,
   CONFIGURATION*        config,
   SCIP_SOL*             baselpsol,
   CANDIDATE*            candidate,
   SCIP_Real             localbaselpsolval,
   int                   probingdepth,
   int                   recursiondepth,
   DOMAINREDUCTIONS*     updomainreductions,
   BINCONSDATA*          binconsdata,
   BRANCHINGRESULTDATA*  upbranchingresult,
   SCORECONTAINER*       scorecontainer,
   LPIMEMORY*            lpimemory,          /**< The memory to store the lpi data in. May be NULL. */
   SCIP_Bool*            addeddomainreduction
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
   ,LOCALSTATISTICS*     localstats
#endif
   )
{
   SCIP_VAR* branchvar = candidate->branchvar;
   SCIP_Real branchval = candidate->branchval;
   SCIP_Real branchvalfrac = candidate->fracval;
   SCIP_Longint niterations;

   if( config->usebincons && SCIPvarIsBinary(branchvar) )
   {
      /* In case that the branch variable is binary, add the var to the list.
       * This list is used to generate a set packing constraint for cutoff branches which were reached by only using
       * binary variables.
       * UpBranching on a binary variable x means: x >= 1
       * When this cutoff occurs we have that: x <= 0
       */
      SCIP_CALL( binaryVarListAppend(scip, binconsdata->binaryvars, branchvar) );
   }

   SCIPdebugMessage("Depth <%i>, Started up branching on var <%s> with var > <%g>\n", probingdepth, SCIPvarGetName(branchvar), branchval);

   niterations = SCIPgetNLPIterations(scip);
   SCIP_CALL( executeUpBranching(scip, config, candidate, upbranchingresult, status) );
   SCIPstatistic(
      statistics->nlpssolved[probingdepth]++;
      statistics->nlpiterations[probingdepth] += SCIPgetNLPIterations(scip)-niterations;
   )

   SCIPdebug(
      if( upbranchingresult->cutoff )
      {
         SCIPdebugMessage("Depth <%i>, The solved LP was infeasible and as such cutoff\n", probingdepth);
      }
      else
      {
         SCIPdebugMessage("Depth <%i>, The solved LP feasible and has an objval <%g> or <%g>\n", probingdepth, upbranchingresult->objval, SCIPgetLPObjval(scip));
      }
   )

   if( !upbranchingresult->cutoff && !status->lperror && !status->limitreached )
   {
      int deepernlpcands;
      SCIP_Real localupgain;

      if( lpiMemoryIsWritable(lpimemory, config) )
      {
         SCIP_CALL( storeInLPIMemory(scip, lpimemory) );
         SCIPdebugMessage("Stored smth in LPIMemory after upbranching on var %s\n", SCIPvarGetName(branchvar));
      }

      localupgain = MAX(0, upbranchingresult->objval - localbaselpsolval);

      SCIP_CALL( SCIPupdateVarPseudocost(scip, branchvar, 1.0-branchvalfrac, localupgain, 1.0) );

      if( recursiondepth > 1 )
      {
         SCIP_CALL( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, &deepernlpcands, NULL, NULL) );

         SCIPdebugMessage("Depth <%i>, Upbranching has <%i> candidates.\n", probingdepth, deepernlpcands);

         if( deepernlpcands > 0 )
         {
            CANDIDATELIST* candidates;
            BRANCHRULERESULT* deeperbranchruleresult;
            STATUS* deeperstatus;
            PERSISTENTDATA* deeperpersistent = NULL;
            SCIP_Real deeperlpobjval = upbranchingresult->objval;
#ifdef SCIP_STATISTIC
            LOCALSTATISTICS* deeperlocalstats;
#endif

            SCIP_CALL( candidateListAllocate(scip, &candidates) );

#ifdef SCIP_STATISTIC
            SCIP_CALL( getCandidates(scip, config, candidates, statistics) );
#else
            SCIP_CALL( getCandidates(scip, config, candidates) );
#endif

            SCIPdebugMessage("Now the objval is <%g> or <%g>\n", upbranchingresult->objval, SCIPgetLPObjval(scip));

            SCIP_CALL( statusAllocate(scip, &deeperstatus) );

            SCIP_CALL( branchRuleResultAllocateReduced(scip, &deeperbranchruleresult) );

#ifdef SCIP_STATISTIC
            SCIP_CALL( localStatisticsAllocate(scip, &deeperlocalstats) );
            SCIP_CALL( selectVarRecursive(scip, deeperstatus, deeperpersistent, config, baselpsol, updomainreductions, binconsdata,
                  candidates, deeperbranchruleresult, scorecontainer, recursiondepth - 1, deeperlpobjval,
                  statistics, deeperlocalstats) );
#else
            SCIP_CALL( selectVarRecursive(scip, deeperstatus, deeperpersistent, config, baselpsol, updomainreductions, binconsdata,
                  candidates, deeperbranchruleresult, scorecontainer, recursiondepth - 1, deeperlpobjval) );
#endif

            /* The proved dual bound of the deeper branching cannot be less than the current dual bound, as every deeper
             * node has more/tighter constraints and as such cannot be better than the base LP. */
            assert(SCIPisGE(scip, deeperbranchruleresult->decision->proveddb, upbranchingresult->dualbound));
            upbranchingresult->dualbound = deeperbranchruleresult->decision->proveddb;
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
               /* upbranchingresult->cutoff is TRUE, if the up child was directly infeasible (so here it is always
                * false, as we don't want to branch on an infeasible node)
                * deeperstatus->cutoff is TRUE, if any up/down child pair of the up child were cutoff
                * */
               upbranchingresult->cutoff = deeperstatus->cutoff;
               SCIPstatistic(
                  localstats->ncutoffproofnodes += deeperlocalstats->ncutoffproofnodes;
               )

               SCIPdebugMessage("Depth <%i>, Both deeper children were cutoff, so the up branch is cutoff\n", probingdepth);
            }

#ifdef SCIP_STATISTIC
            localStatisticsFree(scip, &deeperlocalstats);
#endif
            statusFree(scip, &deeperstatus);
            branchRuleResultFreeReduced(scip, &deeperbranchruleresult);
            SCIP_CALL( candidateListFree(scip, &candidates) );
         }
      }
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
      binaryVarListDrop(scip, binconsdata->binaryvars);
#else
      {
         SCIP_VAR* droppedelement;
         droppedelement = binaryVarListDrop(scip, binconsdata->binaryvars);
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
   CANDIDATELIST*        candidates,
   BRANCHRULERESULT*     branchruleresult,
   SCORECONTAINER*       scorecontainer,
   int                   recursiondepth,
   SCIP_Real             lpobjval
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
   ,LOCALSTATISTICS*     localstats
#endif
   )
{
   BRANCHINGDECISION* decision;
   int nlpcands;
   int probingdepth;

   assert(scip != NULL);
   assert(status != NULL);
   assert(config != NULL);
   assert(!config->usedomainreduction || domainreductions != NULL);
   assert(!config->usebincons || binconsdata != NULL);
   assert(candidates != NULL);
   assert(candidates->ncandidates > 0);
   assert(branchruleresult != NULL);
   assert(branchruleresult->decision != NULL);
   assert(recursiondepth >= 1);
   SCIPstatistic( assert(statistics != NULL) );

   decision = branchruleresult->decision;
   nlpcands = candidates->ncandidates;
   probingdepth = SCIPgetProbingDepth(scip);

   /* init default decision */
   decision->var = candidates->candidates[0]->branchvar;
   decision->val = candidates->candidates[0]->branchval;
   decision->downdb = lpobjval;
   decision->downdbvalid = FALSE;
   decision->updb = lpobjval;
   decision->updbvalid = FALSE;
   decision->proveddb = lpobjval;

   if( !config->forcebranching && nlpcands == 1 )
   {
      SCIPdebugMessage("Only one candidate (<%s>) is given. This one is chosen without calculations.\n", SCIPvarGetName(decision->var));
      SCIPstatistic( statistics->nsinglecandidate[probingdepth]++; )
   }
   else
   {

      if( SCIP_MAXTREEDEPTH <= (SCIPgetDepth(scip) + recursiondepth) )
      {
         /* we need at least 'recursiondepth' space for the branching */
         SCIPdebugMessage("Cannot perform probing in selectVarRecursive, depth limit reached. Current:<%i>, Max:<%i>\n",
            SCIP_MAXTREEDEPTH, SCIPgetDepth(scip) + recursiondepth);
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
         SCIP_LPI* lpi;
         BRANCHINGRESULTDATA* downbranchingresult;
         BRANCHINGRESULTDATA* upbranchingresult;

         bestscorelowerbound = SCIPvarGetLbLocal(decision->var);
         bestscoreupperbound = SCIPvarGetUbLocal(decision->var);

         SCIP_CALL( branchingResultDataAllocate(scip, &downbranchingresult) );
         SCIP_CALL( branchingResultDataAllocate(scip, &upbranchingresult) );

         SCIP_CALL( SCIPgetLPI(scip, &lpi) );

         SCIPdebugMessage("Started selectVarRecursive with <%i> candidates.\n", nlpcands);

         for( i = 0, c = persistent ? persistent->restartindex : 0; i < nlpcands && isExecuteBranchingLoop(status, config) && !SCIPisStopped(scip); i++, c++ )
         {
            SCIP_Bool useoldbranching = FALSE;
            SCIP_Bool addeddomainreduction = FALSE;
            CANDIDATE* candidate;
            SCIP_VAR* branchvar;
            SCIP_Real branchval;
            SCIP_Real branchvalfrac;

            c = c % nlpcands;

            candidate = candidates->candidates[c];

            assert(candidate != NULL);

            branchvar = candidate->branchvar;
            branchval = candidate->branchval;
            branchvalfrac = candidate->fracval;

            assert(branchvar != NULL);

            /* Reset the cutoffproofnodes, as the number of proof nodes from previous branching vars (which where not
             * cutoff, as we didn't break the loop) is not relevant for the min total sum of proof nodes. */
            SCIPstatistic( localstats->ncutoffproofnodes = 0; )

            branchingResultDataInit(scip, downbranchingresult);
            branchingResultDataInit(scip, upbranchingresult);

            if( persistent != NULL && isUseOldBranching(scip, config, persistent, branchvar) ){
               getOldBranching(scip, persistent, branchvar, downbranchingresult, upbranchingresult);
               useoldbranching = TRUE;
               SCIPstatistic( statistics->noldbranchused[probingdepth]++; )

               SCIPdebugMessage("Depth <%i>, Used old branching results for the up and down branches of <%s>.\n", probingdepth, SCIPvarGetName(branchvar));
            }
            else
            {
               DOMAINREDUCTIONS* downdomainreductions;
               DOMAINREDUCTIONS* updomainreductions;

               LPIMEMORY* downlpimemory = NULL;
               LPIMEMORY* uplpimemory = NULL;

               SCIPdebugMessage("Depth <%i>, Started branching on var <%s>\n", probingdepth, SCIPvarGetName(branchvar));

               if( config->usedomainreduction )
               {
                  SCIP_CALL( domainReductionsAllocate(scip, &downdomainreductions) );
                  SCIP_CALL( domainReductionsAllocate(scip, &updomainreductions) );
               }

               if( branchruleresult->downlpimemories != NULL && branchruleresult->uplpimemories != NULL )
               {
                  downlpimemory = branchruleresult->downlpimemories[i];
                  uplpimemory = branchruleresult->uplpimemories[i];
               }

               if( config->downfirst )
               {
#ifdef SCIP_STATISTIC
                  SCIP_CALL( executeDownBranchingRecursive(scip, status, persistent, config, baselpsol, candidate,
                     localbaselpsolval, probingdepth, recursiondepth, downdomainreductions, binconsdata,
                     downbranchingresult, scorecontainer, downlpimemory, &addeddomainreduction, statistics, localstats) );

                  SCIP_CALL( executeUpBranchingRecursive(scip, status, persistent, config, baselpsol, candidate,
                     localbaselpsolval, probingdepth, recursiondepth, updomainreductions, binconsdata,
                     upbranchingresult, scorecontainer, uplpimemory, &addeddomainreduction, statistics, localstats) );
#else
                  SCIP_CALL( executeDownBranchingRecursive(scip, status, persistent, config, baselpsol, candidate,
                     localbaselpsolval, probingdepth, recursiondepth, downdomainreductions, binconsdata,
                     downbranchingresult, scorecontainer, downlpimemory, &addeddomainreduction) );

                  SCIP_CALL( executeUpBranchingRecursive(scip, status, persistent, config, baselpsol, candidate,
                     localbaselpsolval, probingdepth, recursiondepth, updomainreductions, binconsdata,
                     upbranchingresult, scorecontainer, uplpimemory, &addeddomainreduction) );
#endif
               }
               else
               {
#ifdef SCIP_STATISTIC
                  SCIP_CALL( executeUpBranchingRecursive(scip, status, persistent, config, baselpsol, candidate,
                     localbaselpsolval, probingdepth, recursiondepth, updomainreductions, binconsdata,
                     upbranchingresult, scorecontainer, uplpimemory, &addeddomainreduction, statistics, localstats) );

                  SCIP_CALL( executeDownBranchingRecursive(scip, status, persistent, config, baselpsol, candidate,
                     localbaselpsolval, probingdepth, recursiondepth, downdomainreductions, binconsdata,
                     downbranchingresult, scorecontainer, downlpimemory, &addeddomainreduction, statistics, localstats) );
#else
                  SCIP_CALL( executeUpBranchingRecursive(scip, status, persistent, config, baselpsol, candidate,
                     localbaselpsolval, probingdepth, recursiondepth, updomainreductions, binconsdata,
                     upbranchingresult, scorecontainer, uplpimemory, &addeddomainreduction) );

                  SCIP_CALL( executeDownBranchingRecursive(scip, status, persistent, config, baselpsol, candidate,
                     localbaselpsolval, probingdepth, recursiondepth, downdomainreductions, binconsdata,
                     downbranchingresult, scorecontainer, downlpimemory, &addeddomainreduction) );
#endif
               }

               if( persistent != NULL )
               {
                  updateOldBranching(scip, persistent, branchvar, downbranchingresult, upbranchingresult);
               }

               if( config->usedomainreduction )
               {
                  applyDeeperDomainReductions(scip, baselpsol, domainreductions, downdomainreductions, updomainreductions);

                  domainReductionsFree(scip, &updomainreductions);
                  domainReductionsFree(scip, &downdomainreductions);
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

                  if( scorecontainer != NULL )
                  {
                     scoreContainerSetScore(scorecontainer, branchvar, SCIPinfinity(scip));
                  }

                  if( branchruleresult->candscores != NULL )
                  {
                     branchruleresult->candscores[i] = SCIPinfinity(scip);
                     branchruleresult->candswithscore[i] = branchvar;
                     branchruleresult->candslpvalue[i] = branchval;
                     branchruleresult->candsvalfrac[i] = branchvalfrac;
                  }
               }
               else if( upbranchingresult->cutoff )
               {
                  /* TODO: rework scoring with one sided cutoff, see no cutoff case */
                  SCIP_Real score = downbranchingresult->dualboundvalid ? downbranchingresult->dualbound : downbranchingresult->objval;

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
                     decision->proveddb = MAX(decision->proveddb, downbranchingresult->dualbound);
                  }

                  if( scorecontainer != NULL )
                  {
                     scoreContainerSetScore(scorecontainer, branchvar, score)
                  }

                  if( branchruleresult->candscores != NULL )
                  {
                     branchruleresult->candscores[i] = score;
                     branchruleresult->candswithscore[i] = branchvar;
                     branchruleresult->candslpvalue[i] = branchval;
                     branchruleresult->candsvalfrac[i] = branchvalfrac;
                  }
               }
               else if( downbranchingresult->cutoff )
               {
                  /* TODO: rework scoring with one sided cutoff, see no cutoff case */
                  SCIP_Real score = upbranchingresult->dualboundvalid ? upbranchingresult->dualbound : upbranchingresult->objval;
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
                     decision->proveddb = MAX(decision->proveddb, upbranchingresult->dualbound);
                  }

                  if( scorecontainer != NULL )
                  {
                     scoreContainerSetScore(scorecontainer, branchvar, score)
                  }

                  if( branchruleresult->candscores != NULL )
                  {
                     branchruleresult->candscores[i] = score;
                     branchruleresult->candswithscore[i] = branchvar;
                     branchruleresult->candslpvalue[i] = branchval;
                     branchruleresult->candsvalfrac[i] = branchvalfrac;
                  }
               }
               else if( !status->limitreached )
               {
                  SCIP_Real downdualbound;
                  SCIP_Real downgain;
                  SCIP_Real updualbound;
                  SCIP_Real upgain;
                  SCIP_Real score;

                  SCIPdebugMessage("Depth <%i>, Neither branch is cutoff and no limit reached.\n", probingdepth);

                  /* TODO: currently scoring function is best dual bound. Maybe make it more generic? */
                  /* TODO: what if the dualbounds are not valid? Can this happen here? */
                  downdualbound = downbranchingresult->dualbound;
                  updualbound = upbranchingresult->dualbound;

                  downgain = MAX(0, downdualbound - lpobjval);
                  upgain = MAX(0, updualbound - lpobjval);

                  score = MIN(downgain, upgain);

                  if( scorecontainer != NULL )
                  {
                     scoreContainerSetScore(scorecontainer, branchvar, score);
                  }

                  if( branchruleresult->candscores != NULL )
                  {
                     branchruleresult->candswithscore[i] = branchvar;
                     branchruleresult->candscores[i] = score;
                     branchruleresult->candslpvalue[i] = branchval;
                     branchruleresult->candsvalfrac[i] = branchvalfrac;
                  }

                  if( SCIPisGT(scip, score, bestscore) )
                  {
                     SCIPdebugMessage("Depth <%i>, Old best var <%s> with LB <%g> and UB <%g>\n", probingdepth, SCIPvarGetName(decision->var), bestscorelowerbound, bestscoreupperbound);

                     bestscore = score;

                     decision->var = branchvar;
                     decision->val = branchval;
                     decision->downdb = downdualbound;
                     decision->downdbvalid = downbranchingresult->dualboundvalid;
                     decision->updb = updualbound;
                     decision->updbvalid = upbranchingresult->dualboundvalid;

                     bestscorelowerbound = SCIPvarGetLbLocal(decision->var);
                     bestscoreupperbound = SCIPvarGetUbLocal(decision->var);

                     SCIPdebugMessage("Depth <%i>, New best var <%s> with LB <%g> and UB <%g>\n", probingdepth, SCIPvarGetName(decision->var), bestscorelowerbound, bestscoreupperbound);
                  }

                  if( upbranchingresult->dualboundvalid && downbranchingresult->dualboundvalid )
                  {
                     decision->proveddb = MAX(decision->proveddb, MIN(updualbound, downdualbound));
                  }
                  else if( upbranchingresult->dualboundvalid )
                  {
                     decision->proveddb = MAX(decision->proveddb, updualbound);
                  }
                  else if( downbranchingresult->dualboundvalid )
                  {
                     decision->proveddb = MAX(decision->proveddb, downdualbound);
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

                     if( nimpliedbincons + ndomreds > config->maxnviolatedcons )
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

            if( areBoundsChanged(scip, decision->var, bestscorelowerbound, bestscoreupperbound) )
            {
               /* in case the bounds of the current highest scored solution have changed due to domain propagation during the
                * lookahead branching we can/should not branch on this variable but instead return the SCIP_REDUCEDDOM result */
               status->propagationdomred = TRUE;
               SCIPstatistic( statistics->npropdomred[probingdepth]++; )

               SCIPdebugMessage("Depth <%i>, Found a domain reduction via the domain propagation of SCIP.\n", probingdepth);
            }
         }

         branchingResultDataFree(scip, &upbranchingresult);
         branchingResultDataFree(scip, &downbranchingresult);

         if( persistent != NULL && !config->abbreviated )
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
   CONFIGURATION*        config,
   PERSISTENTDATA*       persistent,
   STATUS*               status,
   BRANCHRULERESULT*     branchruleresult,
   SCORECONTAINER*       scorecontainer,
   SCIP_Real             lpobjval
#ifdef SCIP_STATISTIC
   ,STATISTICS*          statistics
   ,LOCALSTATISTICS*     localstats
#endif
   )
{
   int recursiondepth;
   CANDIDATELIST* candidates;
   DOMAINREDUCTIONS* domainreductions = NULL;
   BINCONSDATA* binconsdata = NULL;
   SCIP_SOL* baselpsol = NULL;
   SCIP_Bool inprobing;

   assert(scip != NULL);
   assert(config != NULL);
   assert(status != NULL);
   assert(branchruleresult != NULL);
   assert(branchruleresult->decision != NULL);
   SCIPstatistic( assert(statistics != NULL); )

   inprobing = SCIPinProbing(scip);
   recursiondepth = config->recursiondepth;

   assert(recursiondepth > 0);

   SCIP_CALL( candidateListAllocate(scip, &candidates) );

#ifdef SCIP_STATISTIC
   SCIP_CALL( getCandidates(scip, config, candidates, statistics) );
#else
   SCIP_CALL( getCandidates(scip, config, candidates) );
#endif

   assert(candidates->ncandidates > 0);

   if( config->usedomainreduction || config->usebincons )
   {
      SCIP_CALL( copyCurrentSolution(scip, &baselpsol) );
   }

   if( config->usedomainreduction )
   {
      SCIP_CALL( domainReductionsAllocate(scip, &domainreductions) );
   }

   if( config->usebincons )
   {
      SCIP_CALL( binConsDataAllocate(scip, &binconsdata, recursiondepth, SCIPceil(scip, 0.5*candidates->ncandidates)) );
   }

   if( !inprobing )
   {
      SCIPdebugMessage("About to start probing.\n");
      SCIPstartProbing(scip);
      SCIPenableVarHistory(scip);
   }

#ifdef SCIP_STATISTIC
   SCIP_CALL( selectVarRecursive(scip, status, persistent, config, baselpsol, domainreductions, binconsdata, candidates,
         branchruleresult, scorecontainer, recursiondepth, lpobjval, statistics, localstats) );
#else
   SCIP_CALL( selectVarRecursive(scip, status, persistent, config, baselpsol, domainreductions, binconsdata, candidates,
         branchruleresult, scorecontainer, recursiondepth, lpobjval) );
#endif

   if( !inprobing )
   {
      SCIPendProbing(scip);
      SCIPdebugMessage("Ended probing.\n");
   }

   if( persistent != NULL && ((config->usebincons && binconsdata->createdconstraints->nconstraints > 0)
      || (config->usedomainreduction && domainreductions->nchangedvars > 0)) )
   {
      SCIP_CALL( SCIPlinkLPSol(scip, persistent->prevbinsolution) );
      SCIP_CALL( SCIPunlinkSol(scip, persistent->prevbinsolution) );

      branchingDecisionCopy(branchruleresult->decision, persistent->prevdecision);
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
      binConsDataFree(scip, &binconsdata);

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
      domainReductionsFree(scip, &domainreductions);
   }

   if( config->usedomainreduction || config->usebincons )
   {
      SCIP_CALL( SCIPfreeSol(scip, &baselpsol) );
   }

   SCIP_CALL( candidateListFree(scip, &candidates) );

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
   return branchingDecisionIsValid(persistent->prevdecision)
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
      SCIPvarGetName(decision->var));

   /* reset the var pointer, as this is our indicator whether we should branch on prev data in the next call */
   decision->var = NULL;

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

   branchruledata->persistent->prevdecision->var = NULL;
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
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nlpssolvedfsb, recursiondepth) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nlpiterations, recursiondepth) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nlpiterationsfsb, recursiondepth) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->npropdomred, recursiondepth) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->nsinglecandidate, recursiondepth) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->statistics->noldbranchused, recursiondepth) );

   branchruledata->statistics->recursiondepth = recursiondepth;

   statisticsInit(scip, branchruledata->statistics);

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

   branchruledata = SCIPbranchruleGetData(branchrule);
   statistics = branchruledata->statistics;

   statisticsPrint(scip, statistics);

   SCIPfreeMemoryArray(scip, &statistics->noldbranchused);
   SCIPfreeMemoryArray(scip, &statistics->nsinglecandidate);
   SCIPfreeMemoryArray(scip, &statistics->npropdomred);
   SCIPfreeMemoryArray(scip, &statistics->nlpiterationsfsb);
   SCIPfreeMemoryArray(scip, &statistics->nlpiterations);
   SCIPfreeMemoryArray(scip, &statistics->nlpssolvedfsb);
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
      int nlpcands;
      BRANCHRULERESULT* branchruleresult;
      SCORECONTAINER* scorecontainer = NULL;
      STATUS* status;
      SCIP_Real lpobjval;
#ifdef SCIP_STATISTIC
      LOCALSTATISTICS* localstats;
#endif

      /* get all fractional candidates we can branch on */
      SCIP_CALL( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, &nlpcands, NULL, NULL) );

      SCIPdebugMessage("The base lp has <%i> variables with fractional value.\n", nlpcands);

      /* create a struct to store the algorithm status */
      SCIP_CALL( statusAllocate(scip, &status) );
      /* create a struct to store the branching decision (in case there is one) */
      SCIP_CALL( branchRuleResultAllocateReduced(scip, &branchruleresult) );
      SCIP_CALL( scoreContainerAllocate(scip, &scorecontainer) );

      lpobjval = SCIPgetLPObjval(scip);

      /* execute the main logic */
#ifdef SCIP_STATISTIC
      /* create a struct to store the statistics needed for this single run */
      SCIP_CALL( localStatisticsAllocate(scip, &localstats) );
      SCIP_CALL( selectVarStart(scip, branchruledata->config, branchruledata->persistent, status, branchruleresult, scorecontainer, lpobjval, branchruledata->statistics, localstats) );
#else
      SCIP_CALL( selectVarStart(scip, branchruledata->config, branchruledata->persistent, status, branchruleresult, scorecontainer, lpobjval) );
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
         SCIPdebugMessage("domred: %s\n", status->domred ? "true" : "false");
         SCIPdebugMessage("propagationdomred: %s\n", status->propagationdomred ? "true" : "false");
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
         && branchingDecisionIsValid(branchruleresult->decision)
         /*&& (0 <= bestcand && bestcand < nlpcands)*/ /* no valid candidate index could be found */
         )
      {
         BRANCHINGDECISION* decision = branchruleresult->decision;

         SCIPdebugMessage(" -> %d candidates, variable <%s> (solval=%g)\n",
            nlpcands, SCIPvarGetName(decision->var), decision->val);

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

#ifdef SCIP_STATISTIC
      localStatisticsFree(scip, &localstats);
#endif
      scoreContainerFree(scip, &scorecontainer);
      branchRuleResultFreeReduced(scip, &branchruleresult);
      statusFree(scip, &status);
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
#ifdef SCIP_STATISTIC
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitLookahead) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitLookahead) );
#endif
   SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitSolLookahead) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpLookahead) );

   /* add lookahead branching rule parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/useimpliedbincons",
         "should implied binary constraints found via cutoffs on the second level be applied?",
         &branchruledata->config->usebincons, TRUE, DEFAULT_USEBINARYCONSTRAINTS, NULL, NULL) );
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
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/abbreviated",
         "should the abbreviated version of the LAB be used?",
         &branchruledata->config->abbreviated, TRUE, DEFAULT_ABBREVIATED, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead/maxncands",
         "If abbreviated, at most how many cands should be handled?",
         &branchruledata->config->maxncands, TRUE, DEFAULT_MAXNCANDS, 2, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/stopaftercutoff",
         "Should a branching loop be stopped, if a cutoff of one variable is detected?",
         &branchruledata->config->stopaftercutoff, TRUE, DEFAULT_STOPAFTERCUTOFF, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/lookahead/reusebasis",
         "If abbreviated, should an optimal FSB basis be reused in the first ALAB level?",
         &branchruledata->config->reusebasis, TRUE, DEFAULT_REUSEBASIS, NULL, NULL) );

   return SCIP_OKAY;
}
