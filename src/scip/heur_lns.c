/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_lns.c
 * @brief  "Large neighborhood search heuristic that orchestrates the popular neighborhoods Local Branching, RINS, RENS, DINS etc."
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <strings.h>
#include "scip/heur_lns.h"
#include "scipdefplugins.h"

#define HEUR_NAME             "lns"
#define HEUR_DESC             "Large neighborhood search heuristic that orchestrates the popular neighborhoods Local Branching, RINS, RENS, DINS etc."
#define HEUR_DISPCHAR         'L'
#define HEUR_PRIORITY         -1010000
#define HEUR_FREQ             20
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define NNEIGHBORHOODS 8

#define DEFAULT_SEED 113
/*
 * limit parameters for sub-SCIPs
 */
#define DEFAULT_NODESQUOT 0.1
#define DEFAULT_NODESOFFSET 500LL
#define DEFAULT_NSOLSLIM 3
#define DEFAULT_MINNODES 50LL
#define DEFAULT_MAXNODES 5000LL
#define LPLIMFAC 4.0

/*
 * parameters for the minimum improvement
 */
#define DEFAULT_MINIMPROVELOW 0.0001
#define DEFAULT_MINIMPROVEHIGH 0.1
#define MINIMPROVEFAC          1.5
#define DEFAULT_STARTMINIMPROVE 0.05
#define DEFAULT_ADJUSTMINIMPROVE TRUE

/*
 * bandit algorithm parameters
 */
#define DEFAULT_BESTSOLWEIGHT 3
#define DEFAULT_BANDITALGO 'e'    /**< the default bandit algorithm: (u)pper confidence bounds, (e)xp.3, epsilon (g)reedy */
#define DEFAULT_GAMMA 0.2         /**< default weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution for exp3 */
#define DEFAULT_BETA 0.0          /**< default gain offset between 0 and 1 at every observation for exp3 */
#define DEFAULT_GAINMEASURE      'b'/**< measure for the gain of a neighborhood? 'b'oolean, 'w'eighted boolean, 'g'ap based? */
#define DEFAULT_NORMBYEFFORT     TRUE /**< should the gain be normalized by the effort? */
#define GAINMEASURES "bwg"
#define DEFAULT_EPS        0.5       /**< probability for exploration in epsilon-greedy bandit algorithm */
#define DEFAULT_RESETWEIGHTS TRUE    /**< should the bandit algorithms be reset when a new problem is read? */
#define DEFAULT_SUBSCIPRANDSEEDS FALSE /**< should random seeds of sub-SCIPs be altered to increase diversification? */

/*
 * parameters to control variable fixing
 */
#define DEFAULT_USEREDCOST TRUE /**< should reduced cost scores be used for variable priorization? */
#define DEFAULT_USEDISTANCES TRUE  /**< should distances from fixed variables be used for variable priorization */
#define DEFAULT_DOMOREFIXINGS TRUE /**< should the LNS heuristic do more fixings by itself based on variable prioritization
                                      *  until the target fixing is reached? */
#define DEFAULT_ADJUSTFIXINGRATE TRUE   /**< should the heuristic adjust the target fixing rate based on the success? */
#define FIXINGRATE_DECAY 0.75  /**< geometric decay for fixing rate adjustments */
#define FIXINGRATE_STARTINC 0.2 /**< initial increment value for fixing rate */
#define DEFAULT_USESUBSCIPHEURS  FALSE   /**< should the heuristic activate other sub-SCIP heuristics during its search?  */
#define DEFAULT_GAINFILENAME "-" /**< file name to store all gains and the selection of the bandit */

/* individual neighborhood parameters */
#define MUTATIONSEED 121
#define CROSSOVERSEED 321

#define DEFAULT_MINFIXINGRATE_RENS 0.3
#define DEFAULT_MAXFIXINGRATE_RENS 0.7
#define DEFAULT_ACTIVE_RENS TRUE
#define DEFAULT_PRIORITY_RENS 1.0

#define DEFAULT_MINFIXINGRATE_RINS 0.2
#define DEFAULT_MAXFIXINGRATE_RINS 0.6
#define DEFAULT_ACTIVE_RINS TRUE
#define DEFAULT_PRIORITY_RINS 1.0

#define DEFAULT_MINFIXINGRATE_MUTATION 0.4
#define DEFAULT_MAXFIXINGRATE_MUTATION 0.9
#define DEFAULT_ACTIVE_MUTATION TRUE
#define DEFAULT_PRIORITY_MUTATION 1.0

#define DEFAULT_MINFIXINGRATE_GINS 0.3
#define DEFAULT_MAXFIXINGRATE_GINS 0.5
#define DEFAULT_ACTIVE_GINS TRUE
#define DEFAULT_PRIORITY_GINS 1.0

#define DEFAULT_MINFIXINGRATE_LOCALBRANCHING 0.0
#define DEFAULT_MAXFIXINGRATE_LOCALBRANCHING 0.9
#define DEFAULT_ACTIVE_LOCALBRANCHING TRUE
#define DEFAULT_PRIORITY_LOCALBRANCHING 1.0

#define DEFAULT_MINFIXINGRATE_PROXIMITY 0.0
#define DEFAULT_MAXFIXINGRATE_PROXIMITY 0.9
#define DEFAULT_ACTIVE_PROXIMITY TRUE
#define DEFAULT_PRIORITY_PROXIMITY 1.0

#define DEFAULT_MINFIXINGRATE_CROSSOVER 0.4
#define DEFAULT_MAXFIXINGRATE_CROSSOVER 0.9
#define DEFAULT_ACTIVE_CROSSOVER TRUE
#define DEFAULT_PRIORITY_CROSSOVER 1.0

#define DEFAULT_MINFIXINGRATE_ZEROOBJECTIVE 0.0
#define DEFAULT_MAXFIXINGRATE_ZEROOBJECTIVE 0.9
#define DEFAULT_ACTIVE_ZEROOBJECTIVE TRUE
#define DEFAULT_PRIORITY_ZEROOBJECTIVE 1.0

#define DEFAULT_MINFIXINGRATE_DINS 0.1
#define DEFAULT_MAXFIXINGRATE_DINS 0.5
#define DEFAULT_ACTIVE_DINS TRUE
#define DEFAULT_PRIORITY_DINS 1.0

#define DEFAULT_NSOLS_CROSSOVER 2 /**< parameter for the number of solutions that crossover should combine */
#define DEFAULT_NPOOLSOLS_DINS 5 /**< number of pool solutions where binary solution values must agree */

/* event handler properties */
#define EVENTHDLR_NAME         "Lns"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"
#define SCIP_EVENTTYPE_LNS (SCIP_EVENTTYPE_LPSOLVED | SCIP_EVENTTYPE_SOLFOUND | SCIP_EVENTTYPE_BESTSOLFOUND)

/*
 * Data structures
 */

/* additional neighborhood data structures */
typedef struct data_crossover DATA_CROSSOVER;
typedef struct data_mutation DATA_MUTATION;
typedef struct data_dins DATA_DINS;
typedef struct NH_FixingRate NH_FIXINGRATE;
typedef struct NH_Stats NH_STATS;
typedef struct Nh NH;
typedef struct VarPrio VARPRIO;

/** callback to let the neighborhood write its suggested variable fixings into a buffer data structure
 *
 * todo comments
 */
 #define DECL_VARFIXINGS(x) SCIP_RETCODE x ( \
   SCIP*                 scip,               \
   NH*                   neighborhood,       \
   SCIP_VAR**            varbuf,             \
   SCIP_Real*            valbuf,             \
   int*                  nfixings,           \
   SCIP_RESULT*          result              \
   )

/** callback for subproblem changes other than variable fixings
 *
 *  todo comments
 */
#define DECL_CHANGESUBSCIP(x) SCIP_RETCODE x (  \
   SCIP*                 sourcescip,         \
   SCIP*                 targetscip,         \
   SCIP_VAR**            subvars,            \
   int*                  ndomchgs,           \
   int*                  nchgobjs,           \
   int*                  naddedconss,        \
   SCIP_Bool*            success             \
   )

/** initialization callback for neighborhoods when a new problem is read
 *
 */
#define DECL_NHINIT(x) SCIP_RETCODE x ( \
   SCIP*                 scip,               \
   NH*                   neighborhood        \
   )

/** deinitialization callback for neighborhoods when exiting a problem
 *
 */
#define DECL_NHEXIT(x) SCIP_RETCODE x ( \
   SCIP*                 scip,               \
   NH*                   neighborhood        \
   )

/** deinitialization callback for neighborhoods before SCIP is freed
 *
 */
#define DECL_NHFREE(x) SCIP_RETCODE x ( \
   SCIP*                 scip,               \
   NH*                   neighborhood        \
   )

/** callback function for special sub-SCIP settings
 *
 * todo comments
 */
#define DECL_SETUPSUBSCIP(x) SCIP_RETCODE x (\
   SCIP*                 sourcescip,         \
   SCIP*                 targetscip          \
)

enum HistIndex
{
   HIDX_OPT = 0,
   HIDX_USR = 1,
   HIDX_NODELIM = 2,
   HIDX_STALLNODE = 3,
   HIDX_INFEAS = 4,
   HIDX_SOLLIM = 5,
   HIDX_OTHER = 6
};

typedef enum HistIndex HISTINDEX;

#define NHISTENTRIES 7

/** statistics for a neighborhood */
struct NH_Stats
{
   SCIP_CLOCK*           setupclock;
   SCIP_CLOCK*           submipclock;
   SCIP_Longint          usednodes;
   SCIP_Longint          lpiterations;
   SCIP_Real             oldupperbound;
   int                   nruns;
   int                   nrunsbestsol;
   int                   statushist[NHISTENTRIES];
   SCIP_Longint          nsolsfound;
   SCIP_Longint          nbestsolsfound;
   int                   presolrounds;
   int                   totalnbinfixings;
   int                   totalnintfixings;
   int                   totalnimplintfixings;
   int                   totalncontfixings;
   int                   totalnfixings;
};

/** fixing rate that can be automatically adjusted */
struct NH_FixingRate
{
   SCIP_Real             minfixingrate;
   SCIP_Real             targetfixingrate;
   SCIP_Real             increment;
   SCIP_Real             maxfixingrate;
};

/** neighborhood data structure with callbacks, statistics, fixingrate */
struct Nh
{
   char*                 name;               /**< the name of this neighborhood */
   NH_FIXINGRATE         fixingrate;         /**< fixing rate for this neighborhood */
   NH_STATS              stats;              /**< statistics for this neighborhood */
   DECL_VARFIXINGS       ((*varfixings));    /**< variable fixings callback for this neighborhood */
   DECL_CHANGESUBSCIP    ((*changesubscip)); /**< callback for subproblem changes other than variable fixings */
   DECL_SETUPSUBSCIP     ((*setupsubscip));  /**< callback for special sub-SCIP setup */
   DECL_NHINIT           ((*nhinit));        /**< initialization callback when a new problem is read */
   DECL_NHEXIT           ((*nhexit));        /**< deinitialization callback when exiting a problem */
   DECL_NHFREE           ((*nhfree));        /**< deinitialization callback before SCIP is freed */
   SCIP_Bool             active;             /**< is this neighborhood active or not? */
   SCIP_Real             priority;           /**< positive call priority to initialize bandit algorithms */
   union
   {
      DATA_MUTATION*     mutation;           /**< mutation data */
      DATA_CROSSOVER*    crossover;          /**< crossover data */
      DATA_DINS*         dins;               /**< dins data */
   }                     data;               /**< data object for neighborhood specific data */
};

/** mutation neighborhood data structure */
struct data_mutation
{
   SCIP_RANDNUMGEN*      rng;
};

struct data_crossover
{
   int                   nsols;              /**< parameter for the number of solutions that crossover should combine */
   SCIP_RANDNUMGEN*      rng;                /**< random number generator to draw from the solution pool */
};

struct data_dins
{
   int                   npoolsols;          /**< number of pool solutions where binary solution values must agree */
};

/**< specific data for the UCB bandit algorithm */
typedef struct BanditDataUcb BANDITDATAUCB;

/** primal heuristic data */
struct SCIP_HeurData
{
   NH**                  neighborhoods;      /**< array of neighborhoods with the best one at the first position */
   char                  banditalgo;         /**< the bandit algorithm: (u)pper confidence bounds, (e)xp.3, epsilon (g)reedy */
   char*                 gainfilename;       /**< file name to store all gains and the selection of the bandit */
   SCIP_BANDIT*          epsgreedynh;        /**< epsilon greedy selector for a neighborhood */
   FILE*                 gainfile;           /**< gain file pointer, or NULL */
   SCIP_BANDIT*          exp3;               /**< exp3 bandit algorithm */
   SCIP_Longint          nodesoffset;        /**< offset added to the nodes budget */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes in a single sub-SCIP */
   SCIP_Longint          targetnodes;        /**< targeted number of nodes to start a sub-SCIP */
   SCIP_Longint          minnodes;           /**< minimum number of nodes required to start a sub-SCIP */
   SCIP_Longint          usednodes;          /**< total number of nodes already spent in sub-SCIPs */
   SCIP_Real             nodesquot;          /**< fraction of nodes compared to the main SCIP for budget computation */
   SCIP_Real             startminimprove;    /**< initial factor by which LNS should at least improve the incumbent */
   SCIP_Real             minimprovelow;      /**< lower threshold for the minimal improvement over the incumbent */
   SCIP_Real             minimprovehigh;     /**< upper bound for the minimal improvement over the incumbent */
   SCIP_Real             minimprove;         /**< factor by which LNS should at least improve the incumbent */
   SCIP_Real             lplimfac;           /**< limit fraction of LPs per node to interrupt sub-SCIP */
   SCIP_Real             gamma;              /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution for exp3 */
   SCIP_Real             beta;               /**< gain offset between 0 and 1 at every observation for exp3 */
   SCIP_Real             eps;                /**< probability for exploration in epsilon-greedy bandit algorithm */
   int                   nneighborhoods;     /**< number of neighborhoods */
   int                   nactiveneighborhoods; /**< number of active neighborhoods */
   int                   ninitneighborhoods; /**< neighborhoods that were used at least one time */
   int                   nsolslim;           /**< limit on the number of improving solutions in a sub-SCIP call */
   int                   seed;               /**< initial random seed for bandit algorithms and random decisions by neighborhoods */
   int                   currneighborhood;   /**< index of currently selected neighborhood */
   int                   ndelayedcalls;      /**< the number of calls  */
   SCIP_Bool             useredcost;         /**< should reduced cost scores be used for variable prioritization? */
   SCIP_Bool             usedistances;       /**< should distances from fixed variables be used for variable prioritization */
   SCIP_Bool             domorefixings;      /**< should the LNS heuristic do more fixings by itself based on variable prioritization
                                               *  until the target fixing is reached? */
   SCIP_Bool             adjustfixingrate;   /**< should the heuristic adjust the target fixing rate based on the success? */
   SCIP_Bool             usesubscipheurs;    /**< should the heuristic activate other sub-SCIP heuristics during its search?  */
   SCIP_Bool             adjustminimprove;   /**< should the factor by which the minimum improvement is bound be dynamically updated? */
   SCIP_Bool             resetweights;       /**< should the bandit algorithms be reset when a new problem is read? */
   SCIP_Bool             subsciprandseeds;   /**< should random seeds of sub-SCIPs be altered to increase diversification? */
   SCIP_Bool             normbyeffort;       /**< should the gain be normalized by the effort? */
   char                  gainmeasure;        /**< measure for the gain of a neighborhood? 'b'oolean, 'w'eighted boolean,
                                               *  'e'ffort based? */
};

/** event handler data */
struct SCIP_EventData
{
   SCIP_VAR**            subvars;            /**< the variables of the subproblem */
   SCIP*                 sourcescip;         /**< original SCIP data structure */
   SCIP_HEUR*            heur;               /**< lns heuristic structure */
   SCIP_Longint          nodelimit;
   SCIP_Real             lplimfac;           /**< limit fraction of LPs per node to interrupt sub-SCIP */
   NH_STATS*             runstats;           /**< run statistics for the current neighborhood */
   SCIP_Bool             allgainsmode;       /**< true if solutions should only be checked for gain comparisons */
};

/** represents limits for the sub-SCIP solving process */
struct SolveLimits
{
   SCIP_Longint          nodelimit;          /**< maximum number of solving nodes for the sub-SCIP */
   SCIP_Real             memorylimit;        /**< memory limit for the sub-SCIP */
   SCIP_Real             timelimit;          /**< time limit for the sub-SCIP */
};

typedef struct SolveLimits SOLVELIMITS;

/** data structure that can be used for variable prioritization for additional fixings */
struct VarPrio
{
   SCIP*                 scip;               /**< SCIP data structure */
   SCIP_Real*            randscores;         /**< random scores for priorization */
   int*                  distances;          /**< breadth-first distances from already fixed variables */
   SCIP_Real*            redcostscores;      /**< reduced cost scores for fixing a variable to a reference value */
   unsigned int          useredcost:1;       /**< should reduced cost scores be used for variable prioritization? */
   unsigned int          usedistances:1;     /**< should distances from fixed variables be used for variable prioritization */
};

/*
 * Local methods
 */

/** Reset target fixing rate */
static
SCIP_RETCODE resetFixingRate(
   SCIP*                 scip,               /**< SCIP data structure */
   NH_FIXINGRATE*        fixingrate          /**< heuristic fixing rate */
   )
{
   assert(scip != NULL);
   assert(fixingrate != NULL);
   fixingrate->increment = FIXINGRATE_STARTINC;

   /* use the middle between the minimum and the maximum fixing rate */
   fixingrate->targetfixingrate = 0.5 * (fixingrate->minfixingrate + fixingrate->maxfixingrate);

   return SCIP_OKAY;
}

/** todo reset the currently active neighborhood */
static
void resetCurrentNeighborhood(
   SCIP_HEURDATA*        heurdata
   )
{
   assert(heurdata != NULL);
   heurdata->currneighborhood = -1;
   heurdata->ndelayedcalls = 0;
}

/** update increment for fixing rate */
static
void updateFixingRateIncrement(
   NH_FIXINGRATE*        fx                  /**< fixing rate */
   )
{
   fx->increment *= FIXINGRATE_DECAY;
   fx->increment = MAX(fx->increment, 0.01);
}


/** Increase fixing rate
 *
 *  decrease also the rate by which the target fixing rate is adjusted
 */
static
void increaseFixingRate(
   NH_FIXINGRATE*        fx                  /**< fixing rate */
   )
{
   fx->targetfixingrate += fx->increment;
   fx->targetfixingrate = MIN(fx->targetfixingrate, fx->maxfixingrate);
   updateFixingRateIncrement(fx);
}

/** Decrease fixing rate
 *
 *  decrease also the rate by which the target fixing rate is adjusted
 */
static
void decreaseFixingRate(
   NH_FIXINGRATE*        fx                  /**< fixing rate */
   )
{
   fx->targetfixingrate -= fx->increment;
   fx->targetfixingrate = MAX(fx->targetfixingrate, fx->minfixingrate);
   updateFixingRateIncrement(fx);
}

/** update fixing rate based on the results of the current run */
static
void updateFixingRate(
   SCIP*                 scip,               /**< SCIP data structure */
   NH*                   neighborhood,       /**< neighborhood */
   SCIP_STATUS           subscipstatus,      /**< status of the sub-SCIP run */
   NH_STATS*             runstats            /**< run statistics for this run */
   )
{
   NH_FIXINGRATE* fx;

   fx = &neighborhood->fixingrate;

   switch (subscipstatus) {
      case SCIP_STATUS_OPTIMAL:
      case SCIP_STATUS_INFEASIBLE:
      case SCIP_STATUS_INFORUNBD:
         decreaseFixingRate(fx);
         break;
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_USERINTERRUPT:
      case SCIP_STATUS_NODELIMIT:
         if( runstats->nbestsolsfound <= 0 )
            increaseFixingRate(fx);
         break;
      case SCIP_STATUS_SOLLIMIT:
         decreaseFixingRate(fx);
         break;
      default:
         break;
   }
}

/** increase target node limit */
static
void increaseTargetNodeLimit(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   heurdata->targetnodes *= 2;
   heurdata->targetnodes = MIN(heurdata->targetnodes, heurdata->maxnodes);
}

/** decrease target node limit */
static
void decreaseTargetNodeLimit(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   heurdata->targetnodes /= 2;
   heurdata->targetnodes = MAX(heurdata->targetnodes, heurdata->minnodes);
}

/** reset target node limit */
static
void resetTargetNodeLimit(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   heurdata->targetnodes = heurdata->minnodes;
}



/** update target node limit based on the current run results */
static
void updateTargetNodeLimit(
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_STATUS           subscipstatus       /**< status of the sub-SCIP run */
   )
{
   switch (subscipstatus) {
      case SCIP_STATUS_OPTIMAL:
      case SCIP_STATUS_INFEASIBLE:
      case SCIP_STATUS_INFORUNBD:
      case SCIP_STATUS_SOLLIMIT:
         decreaseTargetNodeLimit(heurdata);
         break;
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_NODELIMIT:
         increaseTargetNodeLimit(heurdata);
         break;
      case SCIP_STATUS_USERINTERRUPT:
         break;
      default:
         break;
   }
}

/** reset the minimum improvement for the sub-SCIPs */
static
void resetMinimumImprovement(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   assert(heurdata != NULL);
   heurdata->minimprove = heurdata->startminimprove;
}

/** increase minimum mprovement for the sub-SCIPs */
static
void increaseMinimumImprovement(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   assert(heurdata != NULL);

   heurdata->minimprove *= MINIMPROVEFAC;
   heurdata->minimprove = MIN(heurdata->minimprove, heurdata->minimprovehigh);
}

/** decrease the minimum improvement for the sub-SCIPs */
static
void decreaseMinimumImprovement(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   assert(heurdata != NULL);

   heurdata->minimprove /= MINIMPROVEFAC;
   SCIPdebugMessage("%.4f", heurdata->minimprovelow);
   heurdata->minimprove = MAX(heurdata->minimprove, heurdata->minimprovelow);
}

/** update the minimum improvement based on the status of the sub-SCIP */
static
void updateMinimumImprovement(
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_STATUS           subscipstatus,      /**< status of the sub-SCIP run */
   NH_STATS*             runstats            /**< run statistics for this run */
   )
{
   assert(heurdata != NULL);

   /* if the sub-SCIP status was infeasible, we rather want to make the sub-SCIP easier
    * with a smaller minimum improvement.
    *
    * If a solution limit was reached, we may, set it higher.
    */
   switch (subscipstatus) {
      case SCIP_STATUS_INFEASIBLE:
      case SCIP_STATUS_INFORUNBD:
         decreaseMinimumImprovement(heurdata);

         break;
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_BESTSOLLIMIT:
      case SCIP_STATUS_OPTIMAL:
         increaseMinimumImprovement(heurdata);
         break;
      case SCIP_STATUS_NODELIMIT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_USERINTERRUPT:
         if( runstats->nbestsolsfound <= 0 )
            decreaseMinimumImprovement(heurdata);
         break;
      default:
         break;
   }
}

/** Reset neighborhood statistics */
static
SCIP_RETCODE neighborhoodStatsReset(
   SCIP*                scip,                /**< SCIP data structure */
   NH_STATS*            stats                /**< neighborhood statistics */
   )
{
   assert(scip != NULL);
   assert(stats != NULL);

   stats->lpiterations = 0L;
   stats->nbestsolsfound = 0;
   stats->nruns = 0;
   stats->nrunsbestsol = 0;
   stats->nsolsfound = 0;
   stats->presolrounds = 0;
   stats->totalnbinfixings = 0;
   stats->totalncontfixings = 0;
   stats->totalnimplintfixings = 0;
   stats->totalnfixings = 0;
   stats->usednodes = 0L;

   BMSclearMemoryArray(stats->statushist, NHISTENTRIES);

   SCIP_CALL( SCIPresetClock(scip, stats->setupclock) );
   SCIP_CALL( SCIPresetClock(scip, stats->submipclock) );

   return SCIP_OKAY;
}

/** create a neighborhood of the specified name and include it into the LNS heuristic */
static
SCIP_RETCODE lnsIncludeNeighborhood(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the LNS heuristic */
   NH**                  neighborhood,       /**< neighborhood that should be created and included */
   const char*           name,               /**< name to distinguish this neighborhood */
   SCIP_Real             minfixingrate,      /**< default value for minfixingrate parameter of this neighborhood */
   SCIP_Real             maxfixingrate,      /**< default value for maxfixingrate parameter of this neighborhood */
   SCIP_Bool             active,             /**< default value for active parameter of this neighborhood */
   SCIP_Real             priority,           /**< positive call priority to initialize bandit algorithms */
   DECL_VARFIXINGS       ((*varfixings)),    /**< variable fixing callback for this neighborhood, or NULL */
   DECL_CHANGESUBSCIP    ((*changesubscip)), /**< subscip changes callback for this neighborhood, or NULL */
   DECL_SETUPSUBSCIP     ((*setupsubscip)),  /**< setup callback for this neighborhood, or NULL */
   DECL_NHINIT           ((*nhinit)),        /**< initialization callback for neighborhood, or NULL */
   DECL_NHEXIT           ((*nhexit)),        /**< deinitialization callback for neighborhood, or NULL */
   DECL_NHFREE           ((*nhfree))         /**< deinitialization callback before SCIP is freed, or NULL */
   )
{
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(neighborhood != NULL);
   assert(name != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, neighborhood) );
   assert(*neighborhood != NULL);

   SCIP_ALLOC( BMSduplicateMemoryArray(&(*neighborhood)->name, name, strlen(name)+1) );

   SCIP_CALL( SCIPcreateClock(scip, &(*neighborhood)->stats.setupclock) );
   SCIP_CALL( SCIPcreateClock(scip, &(*neighborhood)->stats.submipclock) );

   (*neighborhood)->changesubscip = changesubscip;
   (*neighborhood)->varfixings = varfixings;
   (*neighborhood)->setupsubscip = setupsubscip;
   (*neighborhood)->nhinit = nhinit;
   (*neighborhood)->nhexit = nhexit;
   (*neighborhood)->nhfree = nhfree;

   /* add parameters for this neighborhood */
   sprintf(paramname, "heuristics/lns/%s/minfixingrate", name);
   SCIP_CALL( SCIPaddRealParam(scip, paramname, "minimum fixing rate for this neighborhood",
         &(*neighborhood)->fixingrate.minfixingrate, TRUE, minfixingrate, 0.0, 1.0, NULL, NULL) );
   sprintf(paramname, "heuristics/lns/%s/maxfixingrate", name);
   SCIP_CALL( SCIPaddRealParam(scip, paramname, "maximum fixing rate for this neighborhood",
         &(*neighborhood)->fixingrate.maxfixingrate, TRUE, maxfixingrate, 0.0, 1.0, NULL, NULL) );
   sprintf(paramname, "heuristics/lns/%s/active", name);
   SCIP_CALL( SCIPaddBoolParam(scip, paramname, "is this neighborhood active?",
         &(*neighborhood)->active, TRUE, active, NULL, NULL) );
   sprintf(paramname, "heuristics/lns/%s/priority", name);
   SCIP_CALL( SCIPaddRealParam(scip, paramname, "positive call priority to initialize bandit algorithms",
         &(*neighborhood)->priority, TRUE, priority, 1e-2, 1.0, NULL, NULL) );

   heurdata->neighborhoods[heurdata->nneighborhoods++] = (*neighborhood);

   return SCIP_OKAY;
}

/** release all data and free a neighborhood */
static
SCIP_RETCODE lnsFreeNeighborhood(
   SCIP*                scip,               /**< SCIP data structure */
   NH**                 neighborhood        /**< pointer to neighborhood that should be freed */
   )
{
   NH* nhptr;
   assert(scip != NULL);
   assert(neighborhood != NULL);

   nhptr = *neighborhood;
   assert(nhptr != NULL);

   BMSfreeMemoryArray(&nhptr->name);

   /* release further, neighborhood specific data structures */
   if( nhptr->nhfree != NULL )
   {
      SCIP_CALL( nhptr->nhfree(scip, nhptr) );
   }

   SCIP_CALL( SCIPfreeClock(scip, &nhptr->stats.setupclock) );
   SCIP_CALL( SCIPfreeClock(scip, &nhptr->stats.submipclock) );

   SCIPfreeBlockMemory(scip, neighborhood);
   *neighborhood = NULL;

   return SCIP_OKAY;
}

/** initialize neighborhood specific data */
static
SCIP_RETCODE neighborhoodInit(
   SCIP*                 scip,               /**< SCIP data structure */
   NH*                   neighborhood        /**< neighborhood to initialize */
   )
{
   assert(scip != NULL);
   assert(neighborhood != NULL);

   if( neighborhood->nhinit != NULL )
   {
      SCIP_CALL( neighborhood->nhinit(scip, neighborhood) );
   }

   return SCIP_OKAY;
}

/** deinitialize neighborhood specific data */
static
SCIP_RETCODE neighborhoodExit(
   SCIP*                 scip,               /**< SCIP data structure */
   NH*                   neighborhood        /**< neighborhood to initialize */
   )
{
   assert(scip != NULL);
   assert(neighborhood != NULL);

   if( neighborhood->nhexit != NULL )
   {
      SCIP_CALL( neighborhood->nhexit(scip, neighborhood) );
   }

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE transferSolution(
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_EVENTDATA*       eventdata
   )
{
   SCIP*                 sourcescip;         /**< original SCIP data structure */
   SCIP_VAR**            subvars;            /**< the variables of the subproblem */
   SCIP_HEUR*            heur;               /**< lns heuristic structure */
   SCIP_SOL*             subsol;             /**< solution of the subproblem */
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_Bool  success;
   NH_STATS*  runstats;
   SCIP_SOL*  oldbestsol;

   assert(subscip != NULL);

   subsol = SCIPgetBestSol(subscip);
   assert(subsol != NULL);

   sourcescip = eventdata->sourcescip;
   subvars = eventdata->subvars;
   heur = eventdata->heur;
   runstats = eventdata->runstats;
   assert(sourcescip != NULL);
   assert(sourcescip != subscip);
   assert(heur != NULL);
   assert(subvars != NULL);
   assert(runstats != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(sourcescip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(sourcescip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(sourcescip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(sourcescip, newsol, nvars, vars, subsolvals) );

   oldbestsol = SCIPgetBestSol(sourcescip);
   /* try to add new solution to scip and free it immediately */

   if( eventdata->allgainsmode )
   {
      SCIP_CALL( SCIPcheckSol(sourcescip, newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );

      if( success )
      {
         runstats->nsolsfound++;
         if( SCIPgetSolTransObj(sourcescip, newsol) < SCIPgetCutoffbound(sourcescip) )
            runstats->nbestsolsfound++;
      }

      SCIP_CALL( SCIPfreeSol(sourcescip, &newsol) );
   }
   else
   {
      SCIP_CALL( SCIPtrySolFree(sourcescip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );

      if( success )
      {
         runstats->nsolsfound++;
         if( SCIPgetBestSol(sourcescip) != oldbestsol )
            runstats->nbestsolsfound++;
      }
   }

   SCIPfreeBufferArray(sourcescip, &subsolvals);

   return SCIP_OKAY;
}


/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * we interrupt the solution process
 */
static
SCIP_DECL_EVENTEXEC(eventExecLns)
{
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_LNS);
   assert(eventdata != NULL);

   /* treat the different atomic events */
   switch( SCIPeventGetType(event) )
   {
      case SCIP_EVENTTYPE_SOLFOUND:
      case SCIP_EVENTTYPE_BESTSOLFOUND:
         SCIP_CALL( transferSolution(scip, eventdata) );
         break;
      case SCIP_EVENTTYPE_LPSOLVED:
         /* interrupt solution process of sub-SCIP */
         if( SCIPgetNLPs(scip) > eventdata->lplimfac * eventdata->nodelimit )
         {
            SCIPdebugMsg(scip, "interrupt after  %" SCIP_LONGINT_FORMAT " LPs\n", SCIPgetNLPs(scip));
            SCIP_CALL( SCIPinterruptSolve(scip) );
         }
         break;
      default:
         break;
   }

   return SCIP_OKAY;
}

/** initialize neighborhood statistics before the next run */
static
void initRunStats(
   SCIP*                 scip,               /**< SCIP data structure */
   NH_STATS*             stats               /**< run statistics */
   )
{
   stats->nbestsolsfound = 0;
   stats->nsolsfound = 0;
   stats->lpiterations = 0L;
   stats->usednodes = 0L;
   stats->totalnfixings = 0;
   stats->oldupperbound = SCIPgetUpperbound(scip);
}

/** update run stats after the sub SCIP was solved */
static
void updateRunStats(
   SCIP*                 scip,               /**< SCIP data structure */
   NH_STATS*             stats,              /**< run statistics */
   SCIP*                 subscip             /**< sub-SCIP instance, or NULL */
   )
{
   /* treat an untransformed subscip as if none was created */
   if( subscip != NULL && !SCIPisTransformed(subscip) )
      subscip = NULL;

   stats->lpiterations = subscip != NULL ? SCIPgetNLPIterations(subscip) : 0L;
   stats->usednodes = subscip != NULL ? SCIPgetNNodes(subscip) : 0L;
}

/** get the histogram index for this status */
static
int getHistIndex(
   SCIP_STATUS           subscipstatus       /**< sub-SCIP status */
   )
{
   switch (subscipstatus) {
      case SCIP_STATUS_OPTIMAL:
         return HIDX_OPT;
      case SCIP_STATUS_INFEASIBLE:
         return HIDX_INFEAS;
      case SCIP_STATUS_NODELIMIT:
         return HIDX_NODELIM;
      case SCIP_STATUS_STALLNODELIMIT:
         return HIDX_STALLNODE;
      case SCIP_STATUS_SOLLIMIT:
         return HIDX_SOLLIM;
      case SCIP_STATUS_USERINTERRUPT:
         return HIDX_USR;
      default:
         return HIDX_OTHER;
   }
}

/** update the statistics of the neighborhood based on the sub-SCIP run */
static
void updateNeighborhoodStats(
   SCIP*                 scip,               /**< SCIP data structure */
   NH_STATS*             runstats,           /**< run statistics */
   NH*                   neighborhood,       /**< the selected neighborhood */
   SCIP_STATUS           subscipstatus       /**< sub-SCIP status */
   )
{
   NH_STATS* stats;
   stats = &neighborhood->stats;
   stats->lpiterations += runstats->lpiterations;
   stats->nbestsolsfound += runstats->nbestsolsfound;
   stats->nsolsfound += runstats->nsolsfound;

   if( runstats->nbestsolsfound > 0 )
      stats->nrunsbestsol += DEFAULT_BESTSOLWEIGHT;
   else if( runstats->nsolsfound > 0 )
      stats->nrunsbestsol++;

   stats->usednodes += runstats->usednodes;
   stats->nruns += 1;

   /* update the counter for the subscip status */
   ++stats->statushist[getHistIndex(subscipstatus)];
}


/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_HEURCOPY(heurCopyLns)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lns primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurCopyLns NULL
#endif

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolLns)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lns primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolLns NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolLns)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lns primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolLns NULL
#endif

static
SCIP_DECL_SORTINDCOMP(sortIndCompLns)
{
   VARPRIO* varprio;
   SCIP* scip;

   varprio = (VARPRIO*)dataptr;
   assert(varprio != NULL);

   scip = varprio->scip;
   assert(scip != NULL);

   if( ind1 == ind2 )
      return 0;

   /* priority is on distances, if enabled. The variable which is closer in a breadth-first search sense to
    * the already fixed variables has precedence */
   if( varprio->usedistances )
   {
      int dist1;
      int dist2;

      dist1 = varprio->distances[ind1];
      dist2 = varprio->distances[ind2];

      if( dist1 < 0 )
         dist1 = INT_MAX;

      if( dist2 < 0 )
         dist2 = INT_MAX;

      assert(varprio->distances != NULL);
      if( dist1 < dist2 )
         return -1;
      else if( dist1 > dist2 )
         return 1;
   }

   assert(!varprio->usedistances || varprio->distances[ind1] == varprio->distances[ind2]);

   /* if the indices tie considering reduced costs or distances are disabled -> use reduced cost information instead */
   if( varprio->useredcost )
   {
      assert(varprio->redcostscores != NULL);

      if( SCIPisLT(scip, varprio->redcostscores[ind1], varprio->redcostscores[ind2]) )
         return -1;
      else if( SCIPisGT(scip, varprio->redcostscores[ind1], varprio->redcostscores[ind2]) )
         return 1;
   }

   assert(!varprio->useredcost || SCIPisEQ(scip, varprio->redcostscores[ind1], varprio->redcostscores[ind2]));

   assert(varprio->randscores != NULL);

   if( varprio->randscores[ind1] < varprio->randscores[ind2] )
      return -1;
   else if( varprio->randscores[ind1] > varprio->randscores[ind2] )
      return 1;

   return ind1 - ind2;
}

/** Compute the reduced cost score for this variable in the reference solution */
static
SCIP_Real getVariableRedcostScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< the variable for which the score should be computed */
   SCIP_Real             refsolval           /**< solution value in reference solution */
   )
{
   SCIP_Real bestbound;
   SCIP_Real redcost;
   SCIP_Real score;
   assert(scip != NULL);
   assert(var != NULL);

   /* prefer column variables */
   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return SCIPinfinity(scip);

   redcost = SCIPvarGetBestRootRedcost(var);

   if( SCIPisNegative(scip, redcost) )
      bestbound = SCIPvarGetUbGlobal(var);
   else
      bestbound = SCIPvarGetLbGlobal(var);

   if( SCIPisInfinity(scip, REALABS(bestbound)) )
   {
      if( !SCIPisZero(scip, redcost) )
         return -SCIPinfinity(scip);
      else
         return 0.0;
   }

   score = redcost * (refsolval - bestbound);

   return score;
}

/** add variable and solution value to buffer data structure for variable fixings. The method checks if
 *  the value still lies within the variable bounds. The value stays unfixed otherwise.
 */
static
void add2variableBuffer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< (source) SCIP variable that should be added to the buffer */
   SCIP_Real             val,                /**< fixing value for this variable */
   SCIP_VAR**            varbuf,             /**< variable buffer to store variables that should be fixed */
   SCIP_Real*            valbuf,             /**< value buffer to store fixing values */
   int*                  nfixings,           /**< pointer to number of fixed buffer variables, will be increased by 1 */
   SCIP_Bool             integer             /**< is this an integer variable? */
   )
{
   assert(SCIPisFeasIntegral(scip, val) || ! SCIPvarIsIntegral(var));
   assert(*nfixings < SCIPgetNVars(scip));

   /* round the value to its nearest integer */
   if( integer )
      val = SCIPfloor(scip, val + 0.5);

   /* only add fixing if it is still valid within the global variable bounds. Invalidity
    * of this solution value may come from a dual reduction that was performed after the solution from which
    * this value originated was found
    */
   if( SCIPvarGetLbGlobal(var) <= val && val <= SCIPvarGetUbGlobal(var) )
   {
      varbuf[*nfixings] = var;
      valbuf[*nfixings] = val;
      ++(*nfixings);
   }
}

/** fix additional variables if the ones that the neighborhood found were not enough
 *
 *  todo use not always the best solution for the values, but a reference solution provided by the neighborhood itself
 *
 *  @note it may happen that the target fixing rate is not completely reached. This is the case if intermediate,
 *  dual reductions render the solution values of the incumbent solution (reference solution) infeasible for
 *  the current, global variable bounds.
 */
static
SCIP_RETCODE lnsFixMoreVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the LNS neighborhood */
   SCIP_VAR**            varbuf,
   SCIP_Real*            valbuf,
   int*                  nfixings,
   int                   ntargetfixings,
   SCIP_Bool*            success
   )
{
   VARPRIO varprio;
   SCIP_VAR** vars;
   SCIP_Real* redcostscores;
   SCIP_Real* solvals;
   SCIP_SOL* bestsol;
   SCIP_RANDNUMGEN* rng;
   int* distances;
   int* perm;
   SCIP_Real* randscores;
   int nbinvars;
   int nintvars;
   int nbinintvars;
   int nvars;
   int b;
   int nvarstoadd;

   assert(scip != NULL);
   assert(varbuf != NULL);
   assert(nfixings != NULL);
   assert(success != NULL);
   assert(heurdata != NULL);

   *success = FALSE;

   /* if the user parameter forbids more fixings, return immediately */
   if( ! heurdata->domorefixings )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   nbinintvars = nbinvars + nintvars;

   if( ntargetfixings >= nbinintvars )
      return SCIP_OKAY;

   bestsol = SCIPgetBestSol(scip);
   if( bestsol == NULL )
      return SCIP_OKAY;

   /* determine the number of required additional fixings */
   nvarstoadd = ntargetfixings - *nfixings;
   if( nvarstoadd == 0 )
      return SCIP_OKAY;

   varprio.usedistances = heurdata->usedistances;
   varprio.useredcost = heurdata->useredcost;
   varprio.scip = scip;
   rng = SCIPbanditGetRandnumgen(heurdata->epsgreedynh);
   assert(rng != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &randscores, nbinintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nbinintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &distances, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcostscores, nbinintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nbinintvars) );

   /* initialize variable graph distances from already fixed variables */
   if( *nfixings >= 1 )
   {
      SCIP_CALL( SCIPvariablegraphBreadthFirst(scip, NULL, varbuf, *nfixings, distances, INT_MAX, INT_MAX, ntargetfixings) );
   }
   else
   {
      /* initialize all equal distances to make them irrelevant */
      BMSclearMemoryArray(distances, nbinintvars);
   }

   /* todo filter already variables selected variables directly without penalty */
   varprio.randscores = randscores;
   varprio.distances = distances;
   varprio.redcostscores = redcostscores;

   SCIP_CALL( SCIPgetSolVals(scip, bestsol, nbinintvars, vars, solvals) );

   /* assign scores to every discrete variable of the problem */
   for( b = 0; b < nbinintvars; ++b )
   {
      SCIP_VAR* var = vars[b];
      randscores[b] = SCIPrandomGetReal(rng, 0.0, 1.0);
      perm[b] = b;
      redcostscores[b] = getVariableRedcostScore(scip, var, solvals[b]);
   }


   /* loop over already fixed variables and assign large penalties for their values */
   for( b = 0; b < *nfixings; ++b )
   {
      int probindex = SCIPvarGetProbindex(varbuf[b]);

      if( probindex < nbinintvars )
      {
         randscores[probindex] = 2.0;
         redcostscores[probindex] = SCIP_REAL_MAX;
      }
      distances[probindex] = -1;
   }

   /* use selection algorithm (order of the variables does not matter) for quickly completing the fixing */
   SCIPselectInd(perm, sortIndCompLns, &varprio, nvarstoadd, nbinintvars);

   /* loop over the first elements of the selection defined in permutation. They represent the best variables */
   for( b = 0; b < nvarstoadd; ++b )
   {
      int permindex = perm[b];
      assert(permindex >= 0);
      assert(permindex < nbinintvars);

      add2variableBuffer(scip, vars[permindex], solvals[permindex], varbuf, valbuf, nfixings, TRUE);
   }

   *success = TRUE;

   /* loop over perm array and pick the first variables as additional fixings */
   SCIPfreeBufferArray(scip, &solvals);
   SCIPfreeBufferArray(scip, &redcostscores);
   SCIPfreeBufferArray(scip, &distances);
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &randscores);

   return SCIP_OKAY;
}

/** call variable fixing callback for this neighborhood and orchestrate additional variable fixings, if necessary */
static
SCIP_RETCODE neighborhoodFixVariables(
  SCIP*                  scip,
  SCIP_HEURDATA*         heurdata,           /**< heuristic data of the LNS neighborhood */
  NH*                    neighborhood,
  SCIP_VAR**             varbuf,
  SCIP_Real*             valbuf,
  int*                   nfixings,
  SCIP_RESULT*           result              /**< pointer to store the result of the fixing operation */
  )
{
   int ntargetfixings;

   assert(scip != NULL);
   assert(neighborhood != NULL);
   assert(varbuf != NULL);
   assert(valbuf != NULL);
   assert(nfixings != NULL);
   assert(result != NULL);

   *nfixings = 0;

   *result = SCIP_DIDNOTRUN;

   if( neighborhood->varfixings != NULL )
   {
      SCIP_CALL( neighborhood->varfixings(scip, neighborhood, varbuf, valbuf, nfixings, result) );

      if( *result != SCIP_SUCCESS )
         return SCIP_OKAY;
   }

   assert(neighborhood->varfixings == NULL || *result != SCIP_DIDNOTRUN);

   ntargetfixings = (int)(neighborhood->fixingrate.targetfixingrate * (SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip)));
   SCIPdebugMsg(scip, "Neighborhood Fixings/Target: %d / %d\n",*nfixings, ntargetfixings);

   /** if too few fixings, use a strategy to select more variable fixings: randomized, LP graph, ReducedCost/Ps-Cost based, mix */
   if( (*result == SCIP_SUCCESS || *result == SCIP_DIDNOTRUN) && (*nfixings < ntargetfixings) )
   {
      SCIP_Bool success;
      SCIP_CALL( lnsFixMoreVariables(scip, heurdata, varbuf, valbuf, nfixings, ntargetfixings, &success) );

      if( success )
         *result = SCIP_SUCCESS;
      else if( *result == SCIP_SUCCESS )
         *result = SCIP_DIDNOTFIND;
      else
         *result = SCIP_DIDNOTRUN;

      SCIPdebugMsg(scip, "After additional fixings: %d / %d\n",*nfixings, ntargetfixings);
   }
   else
   {
      SCIPdebugMsg(scip, "No additional fixings performed\n");
   }

   return SCIP_OKAY;
}

/** change the sub-SCIP by restricting variable domains, changing objective coefficients, or adding constraints */
static
SCIP_RETCODE neighborhoodChangeSubscip(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   NH*                   neighborhood,       /**< neighborhood */
   SCIP_VAR**            targetvars,         /**< array of target SCIP variables aligned with source SCIP variables */
   int*                  ndomchgs,
   int*                  nchgobjs,
   int*                  naddedconss,
   SCIP_Bool*            success
   )
{
   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(neighborhood != NULL);
   assert(targetvars != NULL);
   assert(ndomchgs != NULL);
   assert(nchgobjs != NULL);
   assert(naddedconss != NULL);
   assert(success != NULL);

   *success = FALSE;
   *ndomchgs = 0;
   *nchgobjs = 0;
   *naddedconss = 0;

   /* call the change sub-SCIP callback of the neighborhood */
   if( neighborhood->changesubscip != NULL )
   {
      SCIP_CALL( neighborhood->changesubscip(sourcescip, targetscip, targetvars, ndomchgs, nchgobjs, naddedconss, success) );
   }
   else
   {
      *success = TRUE;
   }

   return SCIP_OKAY;
}

/** set sub-SCIP solving limits */
static
SCIP_RETCODE setLimits(
   SCIP*                 subscip,            /**< SCIP data structure */
   SOLVELIMITS*          solvelimits         /**< pointer to solving limits data structure */
   )
{
   assert(subscip != NULL);
   assert(solvelimits != NULL);

   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", solvelimits->nodelimit) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", solvelimits->nodelimit / 2));
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", solvelimits->timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", solvelimits->memorylimit) );

   return SCIP_OKAY;
}

/** determine limits for a sub-SCIP */
static
SCIP_RETCODE determineLimits(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< this heuristic */
   SOLVELIMITS*          solvelimits,        /**< pointer to solving limits data structure */
   SCIP_Bool*            runagain            /**< can we solve another sub-SCIP with these limits */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_Real initfactor;
   assert(scip != NULL);
   assert(heur != NULL);
   assert(solvelimits != NULL);
   assert(runagain != NULL);

   heurdata = SCIPheurGetData(heur);

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &solvelimits->timelimit) );
   if( !SCIPisInfinity(scip, solvelimits->timelimit) )
      solvelimits->timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &solvelimits->memorylimit) );

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( !SCIPisInfinity(scip, solvelimits->memorylimit) )
   {
      solvelimits->memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      solvelimits->memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( solvelimits->timelimit <= 0.0 || solvelimits->memorylimit <= 2.0*SCIPgetMemExternEstim(scip)/1048576.0 )
      *runagain = FALSE;

   /* calculate the search node limit of the heuristic  */
   solvelimits->nodelimit = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));
   solvelimits->nodelimit += heurdata->nodesoffset;
   solvelimits->nodelimit -= heurdata->usednodes;
   solvelimits->nodelimit -= 100 * SCIPheurGetNCalls(heur);

   /* use a smaller budget if not all neighborhoods have been initialized yet */
   assert(heurdata->ninitneighborhoods >= 0);
   initfactor = (heurdata->nactiveneighborhoods - heurdata->ninitneighborhoods + 1.0) / (heurdata->nactiveneighborhoods + 1.0);
   solvelimits->nodelimit = (SCIP_Longint)(solvelimits->nodelimit * initfactor);

   /* check whether we have enough nodes left to call subproblem solving */
   if( solvelimits->nodelimit < heurdata->targetnodes )
      *runagain = FALSE;

   return SCIP_OKAY;
}

/** select a neighborhood depending on the selected bandit algorithm */
static
SCIP_RETCODE selectNeighborhood(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the LNS neighborhood */
   int*                  neighborhoodidx     /**< pointer to store the selected neighborhood index */
   )
{
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(neighborhoodidx != NULL);

   *neighborhoodidx = -1;

   switch (heurdata->banditalgo) {
      case 'g':
         SCIP_CALL( SCIPselectBandit(scip, heurdata->epsgreedynh, neighborhoodidx) );
         break;
      case 'e':
         SCIP_CALL( SCIPselectBandit(scip, heurdata->exp3, neighborhoodidx) );

         break;
      case 'u':
         /*todo implement upper confidence bound selection */
         SCIPerrorMessage("Upper confidence bound selection not implemented yet");
         return SCIP_INVALIDCALL;
      default:
         SCIPerrorMessage("Unknown bandit algorithm '%c' selected\n", heurdata->banditalgo);
         return SCIP_INVALIDCALL;
   }

   assert(*neighborhoodidx >= 0);

   return SCIP_OKAY;
}

/** Calculate gain based on the selected gain measure */
static
SCIP_RETCODE getGain(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the LNS neighborhood */
   NH_STATS*             runstats,           /**< run statistics */
   SCIP_Real*            gainptr             /**< pointer to store the computed gain */
   )
{
   SCIP_Real gain = 0.0;
   assert(gainptr != NULL);

   /* compute the gain for this run based on the runstats */
   switch( heurdata->gainmeasure )
   {
   case 'b':
      if( runstats->nbestsolsfound > 0 )
         gain = 1.0;
      break;
   case 'w':
      if( runstats->nbestsolsfound > 0 )
         gain = 1.0;
      else if( runstats->nsolsfound > 0 )
         gain = 1.0 / DEFAULT_BESTSOLWEIGHT;
      break;
   case 'g':
      /* use the closed gap between the primal and dual bound as gain */
      if( runstats->nbestsolsfound > 0 )
      {
         if( SCIPisEQ(scip, SCIPgetUpperbound(scip), SCIPgetLowerbound(scip)) )
            gain = 1.0;
         else if( SCIPisInfinity(scip, runstats->oldupperbound) )
            gain = 1.0;
         else
         {
            gain = 1.0 - sqrt((SCIPgetUpperbound(scip) - SCIPgetLowerbound(scip)) / (runstats->oldupperbound - SCIPgetLowerbound(scip)));
         }
      }
      break;
   default:
      SCIPerrorMessage("Error, passed unknown character '%c' (only know \"%s\") for measuring gains\n", heurdata->gainmeasure, GAINMEASURES);
      SCIP_CALL( SCIP_INVALIDDATA );
      break;
   }

   if( heurdata->normbyeffort )
   {
      SCIP_Real effort;

      assert(runstats->usednodes >= 0);
      assert(runstats->totalnfixings >= 0);
      /* just add one node to avoid division by zero */
      effort = runstats->usednodes / (SCIP_Real)(heurdata->minnodes + 1.0);

      /* assume that every fixed variable linearly reduces the subproblem complexity */
      effort = (1.0 - (runstats->totalnfixings / (SCIP_Real)SCIPgetNVars(scip))) * effort;

      /* gain can be larger than 1.0 if a best solution was found within 0 nodes  */
      gain /= (effort + 1.0);
   }

   *gainptr = gain;
   return SCIP_OKAY;
}

/** update internal bandit algorithm statistics for future draws */
static
SCIP_RETCODE updateBanditAlgorithms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the LNS neighborhood */
   SCIP_Real             gain,               /**< measured gain */
   int                   neighborhoodidx     /**< the neighborhood that was chosen */
   )
{
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(neighborhoodidx >= 0);
   assert(neighborhoodidx < heurdata->nactiveneighborhoods);

   switch (heurdata->banditalgo) {
      case 'g':
         SCIP_CALL( SCIPupdateBandit(scip, heurdata->epsgreedynh, neighborhoodidx, gain) );
         break;
      case 'u':
         break;
      case 'e':
         SCIPdebugMsg(scip, "Rewarding Exp.3 algorithm action %d with gain %.2f\n", neighborhoodidx, gain);
         SCIP_CALL( SCIPupdateBandit(scip, heurdata->exp3, neighborhoodidx, gain) );
         break;
      default:
         break;
   }

   return SCIP_OKAY;
}

/** set up the sub-SCIP */
static
SCIP_RETCODE setupSubScip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_VAR**            subvars,            /**< array of sub-SCIP variables in the order of the main SCIP */
   SOLVELIMITS*          solvelimits,        /**< pointer to solving limits data structure */
   SCIP_HEUR*            heur,               /**< this heuristic */
   SCIP_Bool             objchgd             /**< did the objective change between the source and the target SCIP? */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_Real cutoff;
   SCIP_Real upperbound;

   heurdata = SCIPheurGetData(heur);

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console unless we are in debug mode */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* disable statistic timing inside sub SCIP */
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );

#ifdef LNS_SUBSCIPOUTPUT
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 1) );
   /* enable statistic timing inside sub SCIP */
      SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", heurdata->nsolslim) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
   if( ! heurdata->usesubscipheurs )
   {
      SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );
   }

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(subscip, "estimate") != NULL && !SCIPisParamFixed(subscip, "nodeselection/estimate/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* enable conflict analysis and restrict conflict pool */
   if( !SCIPisParamFixed(subscip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", TRUE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/maxstoresize") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "conflict/maxstoresize", 100) );
   }

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handlers; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no decutions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 10) );
   }

   /* add an objective cutoff */
   if( ! SCIPisInfinity(scip, SCIPgetUpperbound(scip)) )
   {
      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
      if( ! SCIPisInfinity(scip, -1.0 * SCIPgetLowerbound(scip)) )
      {
         cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip)
                            + heurdata->minimprove * SCIPgetLowerbound(scip);
      }
      else
      {
         if( SCIPgetUpperbound(scip) >= 0 )
            cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip);
         else
            cutoff = (1 + heurdata->minimprove) * SCIPgetUpperbound(scip);
      }
      cutoff = MIN(upperbound, cutoff);

      SCIPdebugMsg(scip, "Sub-SCIP cutoff: %15.9" SCIP_REAL_FORMAT " (%15.9" SCIP_REAL_FORMAT " in original space)\n",
         cutoff, SCIPretransformObj(scip, cutoff));

      /* if the objective changed between the source and the target SCIP, encode the cutoff as a constraint */
      if( !objchgd )
      {
         SCIP_CALL(SCIPsetObjlimit(subscip, cutoff));

         SCIPdebugMsg(scip, "Cutoff added as Objective Limit\n");
      }
      else
      {
         SCIP_CONS* objcons;
         int nvars;
         SCIP_VAR** vars;
         int i;

         vars = SCIPgetVars(scip);
         nvars = SCIPgetNVars(scip);

         SCIP_CALL( SCIPcreateConsLinear(subscip, &objcons, "objbound_of_origscip", 0, NULL, NULL, -SCIPinfinity(subscip), cutoff,
                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         for( i = 0; i < nvars; ++i)
         {
            if( !SCIPisFeasZero(subscip, SCIPvarGetObj(vars[i])) )
            {
               SCIP_CALL( SCIPaddCoefLinear(subscip, objcons, subvars[i], SCIPvarGetObj(vars[i])) );
            }
         }
         SCIP_CALL( SCIPaddCons(subscip, objcons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &objcons) );

         SCIPdebugMsg(scip, "Cutoff added as constraint\n");
      }
   }

   /* set solve limits for sub-SCIP */
   SCIP_CALL( setLimits(subscip, solvelimits) );

   /* change random seed of sub-SCIP */
   if( heurdata->subsciprandseeds )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "randomization/randomseedshift", (int)SCIPheurGetNCalls(heur)) );
   }

   SCIPdebugMsg(scip, "Solve Limits: %lld (%lld) nodes (stall nodes), %.1f sec., %d sols\n",
         solvelimits->nodelimit, solvelimits->nodelimit / 2, solvelimits->timelimit, heurdata->nsolslim);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLns)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** varbuf;
   SCIP_Real* valbuf;
   SCIP_VAR** vars;
   SCIP_VAR** subvars;
   NH_STATS runstats[NNEIGHBORHOODS];
   SCIP_STATUS subscipstatus[NNEIGHBORHOODS];
   SCIP* subscip = NULL;

   int nfixings;
   int nvars;
   int neighborhoodidx;
   int ntries;
   SCIP_Bool tryagain;
   NH* neighborhood;
   SOLVELIMITS solvelimits;
   SCIP_Bool success;
   SCIP_Bool run;
   SCIP_Bool allgainsmode;
   SCIP_Real gains[NNEIGHBORHOODS];
   int banditidx;
   int i;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   run = TRUE;
   /** check if budget allows a run of the next selected neighborhood */
   SCIP_CALL( determineLimits(scip, heur, &solvelimits, &run) );
   SCIPdebugMsg(scip, "Budget check: %" SCIP_LONGINT_FORMAT " %s\n", solvelimits.nodelimit, run ? "passed" : "must wait");

   if( !run )
      return SCIP_OKAY;

   allgainsmode = heurdata->gainfile != NULL;

   /* apply some other rules for a fair all gains mode; in normal execution mode, neighborhoods are iterated through */
   if( allgainsmode )
   {
      /* most neighborhoods require an incumbent solution */
      if( SCIPgetNSols(scip) < 2 )
      {
         SCIPdebugMsg(scip, "Not enough solutions for all gains mode\n");
         return SCIP_OKAY;
      }

      /* if the node is infeasible, or has no LP solution, which is required by some neighborhoods
       * if we are not in all gains mode, the neighborhoods delay themselves individually
       */
      if( nodeinfeasible || ! SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      {
         SCIPdebugMsg(scip, "Delay LNS heuristic until a feasible node with optimally solved LP relaxation\n");
         *result = SCIP_DELAYED;
         return SCIP_OKAY;
      }
   }

   /** use the neighborhood that requested a delay or select the next neighborhood to run based on the selected bandit algorithm */
   if( heurdata->currneighborhood >= 0 )
   {
      assert(!allgainsmode);
      banditidx = heurdata->currneighborhood;
      SCIPdebugMsg(scip, "Select delayed neighborhood %d (was delayed %d times)\n", banditidx, heurdata->ndelayedcalls);
   }
   else
   {
      SCIP_CALL( selectNeighborhood(scip, heurdata, &banditidx) );
      SCIPdebugMsg(scip, "Selected neighborhood %d with bandit algorithm\n", banditidx);
   }

   /* in all gains mode, we simply loop over all heuristics */
   if( ! allgainsmode )
      neighborhoodidx = banditidx;
   else
      neighborhoodidx = 0;

   assert(neighborhoodidx >= 0);
   assert(heurdata->nactiveneighborhoods > neighborhoodidx);

   /* allocate memory for variable fixings buffer */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varbuf, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &valbuf, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   /* initialize neighborhood statistics for a run */
   ntries = 1;
   do
   {
      SCIP_HASHMAP* varmapf;
      SCIP_EVENTHDLR* eventhdlr;
      SCIP_EVENTDATA eventdata;
      int ndomchgs;
      int nchgobjs;
      int naddedconss;
      int v;
      SCIP_RESULT fixresult;
      tryagain = FALSE;
      neighborhood = heurdata->neighborhoods[neighborhoodidx];
      SCIPdebugMsg(scip, "Running '%s' neighborhood %d\n", neighborhood->name, neighborhoodidx);

      initRunStats(scip, &runstats[neighborhoodidx]);
      gains[neighborhoodidx] = 0.0;

      subscipstatus[neighborhoodidx] = SCIP_STATUS_UNKNOWN;
      SCIP_CALL( SCIPstartClock(scip, neighborhood->stats.setupclock) );

      /** determine variable fixings and objective coefficients of this neighborhood */
      SCIP_CALL( neighborhoodFixVariables(scip, heurdata, neighborhood, varbuf, valbuf, &nfixings, &fixresult) );

      SCIPdebugMsg(scip, "Fix %d/%d variables\n", nfixings, nvars);

      /* Fixing was not successful, either because the fixing rate was not reached (and no additional variable
       * prioritization was used), or the neighborhood requested a delay, e.g., because no LP relaxation solution exists
       * at the current node
       *
       * The LNS heuristic keeps a delayed neighborhood active and delays itself.
       */
      if( fixresult != SCIP_SUCCESS )
      {
         SCIP_CALL( SCIPstopClock(scip, neighborhood->stats.setupclock) );

         /* to determine all gains, we cannot delay neighborhoods */
         if( allgainsmode )
         {
            if( ntries == heurdata->nactiveneighborhoods )
               break;

            neighborhoodidx = (neighborhoodidx + 1) % heurdata->nactiveneighborhoods;
            ntries++;
            tryagain = TRUE;

            continue;
         }


         /* delay the heuristic along with the selected neighborhood
          *
          * if the neighborhood has been delayed for too many consecutive calls, the delay is treated as a failure */
         if( fixresult == SCIP_DELAYED )
         {

            if( heurdata->ndelayedcalls > (SCIPheurGetFreq(heur) / 4 + 1) )
            {
               resetCurrentNeighborhood(heurdata);

               /* use SCIP_DIDNOTFIND to penalize the neighborhood with a bad reward */
               fixresult = SCIP_DIDNOTFIND;
            }
            else if( heurdata->currneighborhood == -1 )
            {
               heurdata->currneighborhood = neighborhoodidx;
               heurdata->ndelayedcalls = 1;
            }
            else
            {
               heurdata->ndelayedcalls++;
            }
         }

         if( fixresult == SCIP_DIDNOTRUN )
         {
            if( ntries < heurdata->nactiveneighborhoods )
            {
               neighborhoodidx = (neighborhoodidx + 1) % heurdata->nactiveneighborhoods;
               ntries++;
               tryagain = TRUE;

               SCIPdebugMsg(scip, "Neighborhood cannot run -> try next neighborhood %d\n", neighborhoodidx);
               continue;
            }
            else
               tryagain = FALSE;
               break;
         }

         assert(fixresult == SCIP_DIDNOTFIND);
         *result = fixresult;
         break;
      }

      *result = SCIP_DIDNOTFIND;

      neighborhood->stats.totalnfixings += nfixings;
      runstats[neighborhoodidx].totalnfixings = nfixings;

      SCIP_CALL( SCIPcreate(&subscip) );
      SCIP_CALL( SCIPhashmapCreate(&varmapf, SCIPblkmem(scip), nvars) );

      /** todo later: run global propagation for this set of fixings */
      SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapf, neighborhood->name, varbuf, valbuf, nfixings, FALSE, TRUE, &success, NULL) );

      /* store sub-SCIP variables in array for faster access */
      for( v = 0; v < nvars; ++v )
      {
         subvars[v] = (SCIP_VAR*)SCIPhashmapGetImage(varmapf, (void *)vars[v]);
         assert(subvars[v] != NULL);
      }

      SCIPhashmapFree(&varmapf);

      /** let the neighborhood add additional constraints, or restrict domains */
      SCIP_CALL( neighborhoodChangeSubscip(scip, subscip, neighborhood, subvars, &ndomchgs, &nchgobjs, &naddedconss, &success) );

      if( ! success )
      {
         SCIP_CALL( SCIPstopClock(scip, neighborhood->stats.setupclock) );

         if( !allgainsmode || ntries == heurdata->nactiveneighborhoods )
            break;

         neighborhoodidx = (neighborhoodidx + 1) % heurdata->nactiveneighborhoods;
         ntries++;
         tryagain = TRUE;

         SCIP_CALL( SCIPfree(&subscip) );

         continue;
      }

      /* set up sub-SCIP parameters */
      SCIP_CALL( setupSubScip(scip, subscip, subvars, &solvelimits, heur, nchgobjs > 0) );

      /* copy the necessary data into the event data to create new solutions */
      eventdata.nodelimit = solvelimits.nodelimit;
      eventdata.lplimfac = heurdata->lplimfac;
      eventdata.heur = heur;
      eventdata.sourcescip = scip;
      eventdata.subvars = subvars;
      eventdata.runstats = &runstats[neighborhoodidx];
      eventdata.allgainsmode = allgainsmode;

      /* include an event handler to transfer solutions into the main SCIP */
      SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecLns, NULL) );

      /*transform the problem before catching the events */
      SCIP_CALL( SCIPtransformProb(subscip) );
      SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LNS, eventhdlr, &eventdata, NULL) );

      SCIP_CALL( SCIPstopClock(scip, neighborhood->stats.setupclock) );

      /** todo alternatively: set up sub-SCIP and run presolving */
      /** todo was presolving successful enough regarding fixings? otherwise terminate */

      SCIP_CALL( SCIPstartClock(scip, neighborhood->stats.submipclock) );
      /* run sub-SCIP for the given budget, and collect statistics */
      SCIP_CALL_ABORT( SCIPsolve(subscip) );

      SCIP_CALL( SCIPstopClock(scip, neighborhood->stats.submipclock) );

      /* update statistics based on the sub-SCIP run results */
      updateRunStats(scip, &runstats[neighborhoodidx], subscip);
      subscipstatus[neighborhoodidx] = SCIPgetStatus(subscip);

      SCIP_CALL( getGain(scip, heurdata, &runstats[neighborhoodidx], &gains[neighborhoodidx]) );

      if( allgainsmode && ntries < heurdata->nactiveneighborhoods )
      {
         neighborhoodidx = (neighborhoodidx + 1) % heurdata->nactiveneighborhoods;
         ntries++;
         tryagain = TRUE;

         SCIP_CALL( SCIPfree(&subscip) );

         continue;
      }
   }
   while( tryagain && !SCIPisStopped(scip) );

   if( subscip != NULL )
   {
      SCIP_CALL( SCIPfree(&subscip) );
   }

   SCIPfreeBufferArray(scip, &subvars);
   SCIPfreeBufferArray(scip, &valbuf);
   SCIPfreeBufferArray(scip, &varbuf);

   /* update bandit index that may have changed unless we are in all gains mode */
   if( !allgainsmode )
      banditidx = neighborhoodidx;

   if( *result != SCIP_DELAYED )
   {
      /* decrease the number of neighborhoods that have not been initialized */
      if( neighborhood->stats.nruns == 0 )
         --heurdata->ninitneighborhoods;

      heurdata->usednodes += runstats[banditidx].usednodes;

      /** determine the success of this neighborhood, and update the target fixing rate for the next time */
      updateNeighborhoodStats(scip, &runstats[banditidx], heurdata->neighborhoods[banditidx], subscipstatus[banditidx]);

      /* adjust the fixing rate for this neighborhood
       * make no adjustments in all gains mode, because this only affects 1 of 8 heuristics
       */
      if( heurdata->adjustfixingrate && !allgainsmode )
         updateFixingRate(scip, heurdata->neighborhoods[banditidx], subscipstatus[banditidx], &runstats[banditidx]);

      /* similarly, update the minimum improvement for the LNS heuristic
       * make no adjustments in all gains mode
       */
      if( heurdata->adjustminimprove && !allgainsmode )
      {
         SCIPdebugMsg(scip, "Update Minimum Improvement: %.4f\n", heurdata->minimprove);
         updateMinimumImprovement(heurdata, subscipstatus[banditidx], &runstats[banditidx]);
         SCIPdebugMsg(scip, "--> %.4f\n", heurdata->minimprove);
      }

      /* update the target node limit based on the status of the selected algorithm */
      updateTargetNodeLimit(heurdata, subscipstatus[banditidx]);

      /* update the bandit algorithms by the measured gain */
      SCIP_CALL( updateBanditAlgorithms(scip, heurdata, gains[banditidx], banditidx) );

      resetCurrentNeighborhood(heurdata);
   }

   /* write single, measured gains and the bandit index to the gain file */
   if( allgainsmode )
   {
      for( i = 0; i < heurdata->nactiveneighborhoods; ++i )
      {
         fprintf(heurdata->gainfile, "%.4f,", gains[i]);
      }
      fprintf(heurdata->gainfile, "%d\n", banditidx);
   }

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */
static
DECL_VARFIXINGS(varFixingsRens)
{  /*lint --e{715}*/
   int nbinvars;
   int nintvars;
   SCIP_VAR** vars;
   int i;
   assert(scip != NULL);
   assert(varbuf != NULL);
   assert(nfixings != NULL);
   assert(valbuf != NULL);

   *result = SCIP_DELAYED;

   if( ! SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* get variable information */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   /* return if no binary or integer variables are present */
   if( nbinvars + nintvars == 0 )
      return SCIP_OKAY;

   /* loop over binary and integer variables; determine those that should be fixed in the sub-SCIP */
   for( i = 0; i < nbinvars + nintvars; ++i )
   {
      SCIP_VAR* var = vars[i];
      SCIP_Real lpsolval = SCIPgetSolVal(scip, NULL, var);
      assert((i < nbinvars && SCIPvarIsBinary(var)) || (i >= nbinvars && SCIPvarIsIntegral(var)));

      /* fix all binary and integer variables with integer LP solution value */
      if( SCIPisFeasIntegral(scip, lpsolval) )
         add2variableBuffer(scip, var, lpsolval, varbuf, valbuf, nfixings, TRUE);
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

static
DECL_CHANGESUBSCIP(changeSubscipRens)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nintvars;
   int nbinvars;
   int i;

   assert(SCIPhasCurrentNodeLP(sourcescip));
   assert(SCIPgetLPSolstat(sourcescip) == SCIP_LPSOLSTAT_OPTIMAL);

   /* get variable information */
   SCIP_CALL( SCIPgetVarsData(sourcescip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   /* restrict bounds of integer variables with fractional solution value */
   for( i = nbinvars; i < nbinvars + nintvars; ++i )
   {
      SCIP_VAR* var = vars[i];
      SCIP_Real lpsolval = SCIPgetSolVal(sourcescip, NULL, var);

      if( ! SCIPisFeasIntegral(sourcescip, lpsolval) )
      {
         SCIP_Real newlb = SCIPfloor(sourcescip, lpsolval);
         SCIP_Real newub = newlb + 1.0;

         /* only count this as a domain change if the new lower and upper bound are a further restriction */
         if( newlb > SCIPvarGetLbGlobal(subvars[i]) + 0.5 || newub < SCIPvarGetUbGlobal(subvars[i]) - 0.5 )
         {
            SCIP_CALL( SCIPchgVarLbGlobal(targetscip, subvars[i], newlb) );
            SCIP_CALL( SCIPchgVarUbGlobal(targetscip, subvars[i], newub) );
            (*ndomchgs)++;
         }
      }
   }

   *success = TRUE;

   return SCIP_OKAY;
}

/** collect fixings by matching solution values in a collection of solutions for all binary and integer variables,
 *  or for a custom set of variables
 */
static
SCIP_RETCODE fixMatchingSolutionValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sols,               /**< array of 2 or more solutions. It is okay for the array to contain one element
                                               *  equal to NULL to represent the current LP solution */
   int                   nsols,              /**< number of solutions in the array */
   SCIP_VAR**            vars,               /**< variable array for which solution values must agree */
   int                   nvars,              /**< number of variables, or -1 for all binary and integer variables */
   SCIP_VAR**            varbuf,             /**< buffer storage for variable fixings */
   SCIP_Real*            valbuf,             /**< buffer storage for fixing values */
   int*                  nfixings            /**< pointer to store the number of fixings */
   )
{
   int v;
   int nbinintvars;
   SCIP_SOL* firstsol;

   assert(scip != NULL);
   assert(sols != NULL);
   assert(nsols > 1);
   assert(varbuf != NULL);
   assert(valbuf != NULL);
   assert(nfixings != NULL);
   assert(*nfixings == 0);

   if( nvars == -1 || vars == NULL )
   {
      int nbinvars;
      int nintvars;
      SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );
      nbinintvars = nbinvars + nintvars;
      nvars = nbinintvars;
   }
   firstsol = sols[0];
   assert(nvars > 0);

   /* loop over integer and binary variables and check if their solution values match in all solutions */
   /* comment */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Real solval;
      SCIP_VAR* var;
      int s;

      var = vars[v];
      assert((v < SCIPgetNBinVars(scip) && SCIPvarIsBinary(var)) || (v >= SCIPgetNBinVars(scip) && SCIPvarIsIntegral(var)));
      solval = SCIPgetSolVal(scip, firstsol, var);

      /* determine if solution values match in all given solutions */
      for( s = 1; s < nsols; ++s )
      {
         SCIP_Real solval2 = SCIPgetSolVal(scip, sols[s], var);
         if( ! SCIPisFeasEQ(scip, solval, solval2) )
            break;
      }

      /* if we did not break early, all solutions agree on the solution value of this variable */
      if( s == nsols )
      {
         add2variableBuffer(scip, var, solval, varbuf, valbuf, nfixings, TRUE);
      }
   }

   return SCIP_OKAY;
}

static
DECL_VARFIXINGS(varFixingsRins)
{
   /*lint --e{715}*/
   int nbinvars;
   int nintvars;
   SCIP_VAR** vars;
   SCIP_SOL* incumbent;
   SCIP_SOL* sols[2];
   assert(scip != NULL);
   assert(varbuf != NULL);
   assert(nfixings != NULL);
   assert(valbuf != NULL);

   *result = SCIP_DELAYED;

   if( ! SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   incumbent = SCIPgetBestSol(scip);
   if( incumbent == NULL )
      return SCIP_OKAY;

   if( SCIPsolGetOrigin(incumbent) == SCIP_SOLORIGIN_ORIGINAL )
      return SCIP_OKAY;

   /* get variable information */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   /* return if no binary or integer variables are present */
   if( nbinvars + nintvars == 0 )
      return SCIP_OKAY;

   sols[0] = NULL;
   sols[1] = incumbent;

   SCIP_CALL( fixMatchingSolutionValues(scip, sols, 2, vars, nbinvars + nintvars, varbuf, valbuf, nfixings) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

static
DECL_NHINIT(nhInitCrossover)
{  /*lint --e{715}*/
   DATA_CROSSOVER* data;

   data = neighborhood->data.crossover;
   assert(data != NULL);

   if( data->rng != NULL )
      SCIPfreeRandom(scip, &data->rng);

   SCIP_CALL( SCIPcreateRandom(scip, &data->rng, CROSSOVERSEED + SCIPgetNVars(scip)) );

   return SCIP_OKAY;
}

static
DECL_NHEXIT(nhExitCrossover)
{  /*lint --e{715}*/
   DATA_CROSSOVER* data;
   data = neighborhood->data.crossover;

   assert(neighborhood != NULL);
   assert(data->rng != NULL);

   SCIPfreeRandom(scip, &data->rng);

   return SCIP_OKAY;
}

static
DECL_NHFREE(nhFreeCrossover)
{  /*lint --e{715}*/
   assert(neighborhood->data.crossover != NULL);
   SCIPfreeBlockMemory(scip, &neighborhood->data.crossover);

   return SCIP_OKAY;
}



static
DECL_VARFIXINGS(varFixingsCrossover)
{  /*lint --e{715}*/
   DATA_CROSSOVER* data;
   SCIP_RANDNUMGEN* rng;
   SCIP_SOL** sols;
   SCIP_SOL** scipsols;
   int nsols;
   int lastdraw;
   assert(scip != NULL);
   assert(varbuf != NULL);
   assert(nfixings != NULL);
   assert(valbuf != NULL);

   data = neighborhood->data.crossover;

   assert(data != NULL);
   nsols = data->nsols;

   *result = SCIP_DIDNOTRUN;

   /* return if the pool has not enough solutions */
   if( nsols > SCIPgetNSols(scip) )
      return SCIP_OKAY;

   rng = data->rng;
   lastdraw = SCIPgetNSols(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &sols, nsols) );
   scipsols = SCIPgetSols(scip);

   /* draw as many solutions from the pool as required by crossover, biased towards
    * better solutions; therefore, the sorting of the solutions by objective is implicitly used
    */
   while( nsols > 0 )
   {
      /* no need for randomization anymore, exactly nsols many solutions remain for the selection */
      if( lastdraw == nsols )
      {
         int s;

         /* fill the remaining slots 0,...,nsols - 1 by the solutions at the same places */
         for( s = 0; s < nsols; ++s )
            sols[s] = scipsols[s];

         nsols = 0;
      }
      else
      {
         int nextdraw;

         assert(nsols < lastdraw);

         /* draw from the lastdraw - nsols many solutions nsols - 1, ... lastdraw - 1 such that nsols many solution */
         nextdraw = SCIPrandomGetInt(rng, nsols - 1, lastdraw - 1);
         assert(nextdraw >= 0);

         sols[nsols - 1] = scipsols[nextdraw];
         nsols--;
         lastdraw = nextdraw;
      }
   }

   SCIP_CALL( fixMatchingSolutionValues(scip, sols, data->nsols, NULL, -1, varbuf, valbuf, nfixings) );

   *result = SCIP_SUCCESS;

   SCIPfreeBufferArray(scip, &sols);

   return SCIP_OKAY;
}

static
DECL_NHINIT(nhInitMutation)
{  /*lint --e{715}*/
   DATA_MUTATION* data;
   assert(scip != NULL);
   assert(neighborhood != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &neighborhood->data.mutation) );

   data = neighborhood->data.mutation;
   assert(data != NULL);

   SCIP_CALL( SCIPcreateRandom(scip, &data->rng, MUTATIONSEED) );

   return SCIP_OKAY;
}

static
DECL_NHEXIT(nhExitMutation)
{  /*lint --e{715}*/
   DATA_MUTATION* data;
   assert(scip != NULL);
   assert(neighborhood != NULL);
   data = neighborhood->data.mutation;
   assert(data != NULL);

   SCIPfreeRandom(scip, &data->rng);

   SCIPfreeBlockMemory(scip, &neighborhood->data.mutation);

   return SCIP_OKAY;
}

static
DECL_VARFIXINGS(varFixingsMutation)
{  /*lint --e{715}*/
   SCIP_RANDNUMGEN* rng;

   SCIP_VAR** vars;
   SCIP_VAR** varscpy;
   int i;
   int nvars;
   int nbinvars;
   int nintvars;
   int nbinintvars;
   int ntargetfixings;
   SCIP_SOL* incumbentsol;
   SCIP_Real targetfixingrate;

   assert(scip != NULL);
   assert(neighborhood != NULL);
   assert(neighborhood->data.mutation != NULL);
   assert(neighborhood->data.mutation->rng != NULL);
   rng = neighborhood->data.mutation->rng;

   *result = SCIP_DIDNOTRUN;

   /* get the problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   nbinintvars = nbinvars + nintvars;
   if( nbinintvars == 0 )
      return SCIP_OKAY;

   incumbentsol = SCIPgetBestSol(scip);
   if( incumbentsol == NULL )
      return SCIP_OKAY;

   targetfixingrate = neighborhood->fixingrate.targetfixingrate;
   ntargetfixings = (int)(targetfixingrate * nbinintvars) + 1;

   /* don't continue if number of discrete variables is too small to reach target fixing rate */
   if( nbinintvars <= ntargetfixings )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* copy variables into a buffer array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &varscpy, vars, nbinintvars) );

   /* partially perturb the array until the number of target fixings is reached */
   for( i = 0; *nfixings < ntargetfixings && i < nbinintvars; ++i )
   {
      int randint = SCIPrandomGetInt(rng, i, nbinintvars - 1);
      assert(randint < nbinintvars);

      if( randint > i )
      {
         SCIPswapPointers((void**)&varscpy[i], (void**)&varscpy[randint]);
      }
      /* copy the selected variables and their solution values into the buffer */
      add2variableBuffer(scip, varscpy[i], SCIPgetSolVal(scip, incumbentsol, varscpy[i]), varbuf, valbuf, nfixings, TRUE);
   }

   assert(i == nbinintvars || *nfixings == ntargetfixings);

   /* Not reaching the number of target fixings means that there is a significant fraction (at least 1 - targetfixingrate)
    * of variables for which the incumbent solution value does not lie within the global bounds anymore. This is a nonsuccess
    * for the neighborhood (additional fixings are not possible), which is okay because the incumbent solution is
    * significantly outdated
    */
   if( *nfixings == ntargetfixings )
      *result = SCIP_SUCCESS;

   /* free the buffer array */
   SCIPfreeBufferArray(scip, &varscpy);

   return SCIP_OKAY;
}

/** todo add local branching constraint */
static
SCIP_RETCODE addLocalBranchingConstraint(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_VAR**            subvars,            /**< array of sub SCIP variables in same order as source SCIP variables */
   int                   distance,           /**< right hand side of the local branching constraint */
   SCIP_Bool*            success,            /**< pointer to store of a local branching constraint has been successfully added */
   int*                  naddedconss         /**< pointer to increase the number of added constraints */
   )
{
   int nbinvars;
   int i;
   SCIP_SOL* referencesol;
   SCIP_CONS* localbranchcons;
   SCIP_VAR** vars;
   SCIP_Real* consvals;
   SCIP_Real rhs;

   assert(sourcescip != NULL);

   nbinvars = SCIPgetNBinVars(sourcescip);
   vars = SCIPgetVars(sourcescip);

   if( nbinvars <= 3 )
      return SCIP_OKAY;

   referencesol = SCIPgetBestSol(sourcescip);
   if( referencesol == NULL )
      return SCIP_OKAY;

   rhs = (SCIP_Real)distance;
   rhs = MAX(rhs, 2.0);


   SCIP_CALL( SCIPallocBufferArray(sourcescip, &consvals, nbinvars) );

   /* loop over binary variables and fill the local branching constraint */
   for( i = 0; i < nbinvars; ++i )
   {
      if( SCIPisEQ(sourcescip, SCIPgetSolVal(sourcescip, referencesol, vars[i]), 0.0) )
         consvals[i] = 1.0;
      else
      {
         consvals[i] = -1.0;
         rhs -= 1.0;
      }
   }

   /* create the local branching constraint in the target scip */
   SCIP_CALL( SCIPcreateConsBasicLinear(targetscip, &localbranchcons, "localbranch", nbinvars, subvars, consvals, -SCIPinfinity(sourcescip), rhs) );
   SCIP_CALL( SCIPaddCons(targetscip, localbranchcons) );
   SCIP_CALL( SCIPreleaseCons(targetscip, &localbranchcons) );

   *naddedconss = 1;
   *success = TRUE;

   SCIPfreeBufferArray(sourcescip, &consvals);

   return SCIP_OKAY;
}

static
DECL_CHANGESUBSCIP(changeSubscipLocalbranching)
{  /*lint --e{715}*/

   SCIP_CALL( addLocalBranchingConstraint(sourcescip, targetscip, subvars, (int)(0.2 * SCIPgetNBinVars(sourcescip)), success, naddedconss) );

   return SCIP_OKAY;
}

static
DECL_CHANGESUBSCIP(changeSubscipProximity)
{  /*lint --e{715}*/
   SCIP_SOL* referencesol;
   SCIP_VAR** vars;
   int nbinvars;
   int nintvars;
   int nvars;
   int i;

   SCIP_CALL( SCIPgetVarsData(sourcescip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   if( nbinvars == 0 )
      return SCIP_OKAY;

   referencesol = SCIPgetBestSol(sourcescip);
   if( referencesol == NULL )
      return SCIP_OKAY;

   /* loop over binary variables, set objective coefficients based on reference solution in a local branching fashion */
   for( i = 0; i < nbinvars; ++i )
   {
      SCIP_Real newobj;
      if( SCIPgetSolVal(sourcescip, referencesol, vars[i]) < 0.5 )
         newobj = -1.0;
      else
         newobj = 1.0;
      SCIP_CALL( SCIPchgVarObj(targetscip, subvars[i], newobj) );
   }

   /* loop over the remaining variables and change their objective coefficients to 0 */
   for( ; i < nvars; ++i )
   {
      SCIP_CALL( SCIPchgVarObj(targetscip, subvars[i], 0.0) );
   }

   *nchgobjs = nvars;
   *success = TRUE;

   return SCIP_OKAY;
}

static
DECL_CHANGESUBSCIP(changeSubscipZeroobjective)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nvars;
   int i;

   SCIP_CALL( SCIPgetVarsData(sourcescip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* loop over the variables and change their objective coefficients to 0 */
   for( i=0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPchgVarObj(targetscip, subvars[i], 0.0) );
   }

   *nchgobjs = nvars;
   *success = TRUE;

   return SCIP_OKAY;
}

/** compute tightened bounds for integer variables depending on how much the LP and the incumbent solution values differ */
static
void computeIntegerVariableBoundsDins(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP_VAR*             var,                /**< the variable for which bounds should be computed */
   SCIP_Real*            lbptr,              /**< pointer to store the lower bound in the DINS sub-SCIP */
   SCIP_Real*            ubptr               /**< pointer to store the upper bound in the DINS sub-SCIP */
   )
{
   SCIP_Real mipsol;
   SCIP_Real lpsol;

   SCIP_Real lbglobal;
   SCIP_Real ubglobal;
   SCIP_SOL* bestsol;

   /* get the bounds for each variable */
   lbglobal = SCIPvarGetLbGlobal(var);
   ubglobal = SCIPvarGetUbGlobal(var);

   assert(SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);
   /* get the current LP solution for each variable */
   lpsol = SCIPvarGetLPSol(var);

   /* get the current MIP solution for each variable */
   bestsol = SCIPgetBestSol(scip);
   mipsol = SCIPgetSolVal(scip, bestsol, var);


   /* if the solution values differ by 0.5 or more, the variable is rebounded, otherwise it is just copied */
   if( REALABS(lpsol - mipsol) >= 0.5 )
   {
      SCIP_Real range;

      *lbptr = lbglobal;
      *ubptr = ubglobal;

      /* create a equally sized range around lpsol for general integers: bounds are lpsol +- (mipsol-lpsol) */
      range = 2*lpsol-mipsol;

      if( mipsol >= lpsol )
      {
         range = SCIPfeasCeil(scip, range);
         *lbptr = MAX(*lbptr, range);

         /* when the bound new upper bound is equal to the current MIP solution, we set both bounds to the integral bound (without eps) */
         if( SCIPisFeasEQ(scip, mipsol, *lbptr) )
            *ubptr = *lbptr;
         else
            *ubptr = mipsol;
      }
      else
      {
         range = SCIPfeasFloor(scip, range);
         *ubptr = MIN(*ubptr, range);

         /* when the bound new upper bound is equal to the current MIP solution, we set both bounds to the integral bound (without eps) */
         if( SCIPisFeasEQ(scip, mipsol, *ubptr) )
            *lbptr = *ubptr;
         else
            *lbptr = mipsol;
      }

      /* the global domain of variables might have been reduced since incumbent was found: adjust lb and ub accordingly */
      *lbptr = MAX(*lbptr, lbglobal);
      *ubptr = MIN(*ubptr, ubglobal);
   }
   else
   {
      /* the global domain of variables might have been reduced since incumbent was found: adjust it accordingly */
      *lbptr = MAX(mipsol, lbglobal);
      *ubptr = MIN(mipsol, ubglobal);
   }
}

static
DECL_VARFIXINGS(varFixingsDins)
{
   DATA_DINS* data;
   SCIP_SOL* rootlpsol;
   SCIP_SOL** sols;
   int nsols;
   int nbinvars;
   int nintvars;
   SCIP_VAR** vars;
   int v;

   data = neighborhood->data.dins;
   assert(data != NULL);
   nsols = data->npoolsols;
   nsols = MIN(nsols, SCIPgetNSols(scip));

   *result = SCIP_DELAYED;

   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetBestSol(scip) == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   if( nbinvars + nintvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPcreateSol(scip, &rootlpsol, NULL) );

   /* save root solution LP values in solution */
   for( v = 0; v < nbinvars + nintvars; ++v )
   {
      SCIP_CALL( SCIPsetSolVal(scip, rootlpsol, vars[v], SCIPvarGetRootSol(vars[v])) );
   }

   /* add the node and the root LP solution */
   nsols += 2;

   SCIP_CALL( SCIPallocBufferArray(scip, &sols, nsols) );
   sols[0] = NULL; /* node LP solution */
   sols[1] = rootlpsol;

   /* copy the remaining MIP solutions after the LP solutions */
   BMScopyMemoryArray(&sols[2], SCIPgetSols(scip), nsols - 2);

   /* 1. Binary variables are fixed if their values agree in all the solutions */
   if( nbinvars > 0 )
   {
      SCIP_CALL( fixMatchingSolutionValues(scip, sols, nsols, vars, nbinvars, varbuf, valbuf, nfixings) );
   }

   /* 2. Integer variables are fixed if they have a very low distance between the incumbent and the root LP solution */
   for( v = nbinvars; v < nintvars; ++v )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      computeIntegerVariableBoundsDins(scip, vars[v], &lb, &ub);

      if( ub - lb < 0.5 )
      {
         assert(SCIPisFeasIntegral(scip, lb));
         add2variableBuffer(scip, vars[v], lb, varbuf, valbuf, nfixings, TRUE);

      }
   }

   *result = SCIP_SUCCESS;

   SCIPfreeBufferArray(scip, &sols);

   SCIPfreeSol(scip, &rootlpsol);

   return SCIP_OKAY;
}

static
DECL_CHANGESUBSCIP(changeSubscipDins)
{
   SCIP_VAR** vars;
   int nintvars;
   int nbinvars;
   int v;

   SCIP_CALL( SCIPgetVarsData(sourcescip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   /* 1. loop over integer variables and tighten the bounds */
   for( v = nbinvars; v < nintvars; ++v )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      computeIntegerVariableBoundsDins(sourcescip, vars[v], &lb, &ub);

      SCIP_CALL( SCIPchgVarLbGlobal(targetscip, subvars[v], lb) );
      SCIP_CALL( SCIPchgVarUbGlobal(targetscip, subvars[v], ub) );
      ++(*ndomchgs);
   }

   /* 2. add local branching constraint for binary variables */
   addLocalBranchingConstraint(sourcescip, targetscip, subvars, (int)(0.1 * SCIPgetNBinVars(sourcescip)), success, naddedconss );

   *success = TRUE;

   return SCIP_OKAY;
}

static
DECL_NHFREE(nhFreeDins)
{
   assert(neighborhood->data.dins != NULL);

   SCIPfreeBlockMemory(scip, &neighborhood->data.dins);

   return SCIP_OKAY;
}



/** include all neighborhoods */
static
SCIP_RETCODE includeNeighborhoods(
   SCIP*                scip,                /**< SCIP data structure */
   SCIP_HEURDATA*       heurdata             /**< heuristic data of the LNS heuristic */
   )
{
   NH* rens;
   NH* rins;
   NH* mutation;
   NH* localbranching;
   NH* crossover;
   NH* proximity;
   NH* zeroobjective;
   NH* dins;

   heurdata->nneighborhoods = 0;

   /* include the RENS neighborhood */
   SCIP_CALL( lnsIncludeNeighborhood(scip, heurdata, &rens, "rens",
         DEFAULT_MINFIXINGRATE_RENS, DEFAULT_MAXFIXINGRATE_RENS, DEFAULT_ACTIVE_RENS, DEFAULT_PRIORITY_RENS,
         varFixingsRens, changeSubscipRens, NULL, NULL, NULL, NULL) );

   /* include the RINS neighborhood */
   SCIP_CALL( lnsIncludeNeighborhood(scip, heurdata, &rins, "rins",
         DEFAULT_MINFIXINGRATE_RINS, DEFAULT_MAXFIXINGRATE_RINS, DEFAULT_ACTIVE_RINS, DEFAULT_PRIORITY_RINS,
         varFixingsRins, NULL, NULL, NULL, NULL, NULL) );

   /* include the mutation neighborhood */
   SCIP_CALL( lnsIncludeNeighborhood(scip, heurdata, &mutation, "mutation",
         DEFAULT_MINFIXINGRATE_MUTATION, DEFAULT_MAXFIXINGRATE_MUTATION, DEFAULT_ACTIVE_MUTATION, DEFAULT_PRIORITY_MUTATION,
         varFixingsMutation, NULL, NULL, nhInitMutation, nhExitMutation, NULL) );

   /* include the local branching neighborhood */
   SCIP_CALL( lnsIncludeNeighborhood(scip, heurdata, &localbranching, "localbranching",
         DEFAULT_MINFIXINGRATE_LOCALBRANCHING, DEFAULT_MAXFIXINGRATE_LOCALBRANCHING, DEFAULT_ACTIVE_LOCALBRANCHING, DEFAULT_PRIORITY_LOCALBRANCHING,
         NULL, changeSubscipLocalbranching, NULL, NULL, NULL, NULL) );

   /* include the crossover neighborhood */
   SCIP_CALL( lnsIncludeNeighborhood(scip, heurdata, &crossover, "crossover",
         DEFAULT_MINFIXINGRATE_CROSSOVER, DEFAULT_MAXFIXINGRATE_CROSSOVER, DEFAULT_ACTIVE_CROSSOVER, DEFAULT_PRIORITY_CROSSOVER,
         varFixingsCrossover, NULL, NULL,
         nhInitCrossover, nhExitCrossover, nhFreeCrossover) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &crossover->data.crossover) );
   crossover->data.crossover->rng = NULL;

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/lns/crossover/nsols", "the number of solutions that crossover should combine",
         &crossover->data.crossover->nsols, TRUE, DEFAULT_NSOLS_CROSSOVER, 2, 10, NULL, NULL) );

   /* include the Proximity neighborhood */
   SCIP_CALL( lnsIncludeNeighborhood(scip, heurdata, &proximity, "proximity",
         DEFAULT_MINFIXINGRATE_PROXIMITY, DEFAULT_MAXFIXINGRATE_PROXIMITY, DEFAULT_ACTIVE_PROXIMITY, DEFAULT_PRIORITY_PROXIMITY,
         NULL, changeSubscipProximity, NULL, NULL, NULL, NULL) );

   /* include the Zeroobjective neighborhood */
   SCIP_CALL( lnsIncludeNeighborhood(scip, heurdata, &zeroobjective, "zeroobjective",
         DEFAULT_MINFIXINGRATE_ZEROOBJECTIVE, DEFAULT_MAXFIXINGRATE_ZEROOBJECTIVE, DEFAULT_ACTIVE_ZEROOBJECTIVE, DEFAULT_PRIORITY_ZEROOBJECTIVE,
         NULL, changeSubscipZeroobjective,
         NULL, NULL, NULL, NULL) );

   /* include the DINS neighborhood */
   SCIP_CALL( lnsIncludeNeighborhood(scip, heurdata, &dins, "dins",
         DEFAULT_MINFIXINGRATE_DINS, DEFAULT_MAXFIXINGRATE_DINS, DEFAULT_ACTIVE_DINS, DEFAULT_PRIORITY_DINS,
         varFixingsDins, changeSubscipDins, NULL, NULL, NULL, nhFreeDins) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &dins->data.dins) );

   /* add DINS neighborhood parameters  */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/lns/dins/npoolsols",
         "number of pool solutions where binary solution values must agree",
         &dins->data.dins->npoolsols, TRUE, DEFAULT_NPOOLSOLS_DINS, 1, 100, NULL, NULL) );

   /* author bzfhende: TODO include the GINS neighborhood */

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitLns)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int i;
   SCIP_Real* priorities;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->epsgreedynh == NULL);
   heurdata->nactiveneighborhoods = heurdata->nneighborhoods;

   SCIP_CALL( SCIPallocBufferArray(scip, &priorities, heurdata->nactiveneighborhoods) );

   /* init neighborhoods for new problem by resetting their statistics and fixing rate */
   for( i = heurdata->nneighborhoods - 1; i >= 0; --i )
   {
      NH* neighborhood = heurdata->neighborhoods[i];

      SCIP_CALL( neighborhoodInit(scip, neighborhood) );

      SCIP_CALL( resetFixingRate(scip, &neighborhood->fixingrate) );

      SCIP_CALL( neighborhoodStatsReset(scip, &neighborhood->stats) );

      if( ! neighborhood->active )
      {
         if( heurdata->nactiveneighborhoods - 1 > i )
         {
            assert(heurdata->neighborhoods[heurdata->nactiveneighborhoods - 1]->active);
            SCIPswapPointers((void **)&heurdata->neighborhoods[i], (void **)&heurdata->neighborhoods[heurdata->nactiveneighborhoods - 1]);
         }
         heurdata->nactiveneighborhoods--;
      }
   }

   /* collect neighborhood priorities */
   for( i = 0; i < heurdata->nactiveneighborhoods; ++i )
      priorities[i] = heurdata->neighborhoods[i]->priority;


   /* create an exp3 bandit algorithm */
   if( heurdata->exp3 == NULL )
   {
      assert(heurdata->epsgreedynh == NULL);
      SCIP_CALL( SCIPcreateBanditExp3(scip, &heurdata->exp3,
            heurdata->nactiveneighborhoods, heurdata->gamma, heurdata->beta) );

      SCIP_CALL( SCIPcreateBanditEpsgreedy(scip, &heurdata->epsgreedynh, heurdata->eps, heurdata->nactiveneighborhoods) );
      SCIP_CALL( SCIPresetBandit(scip, heurdata->epsgreedynh, priorities, (unsigned int)(heurdata->seed + SCIPgetNVars(scip))) );
   }
   else if( heurdata->resetweights )
   {
      assert(heurdata->epsgreedynh != NULL);

      /* todo active neighborhoods might change between init calls, reset functionality must take this into account */
      SCIP_CALL( SCIPresetBandit(scip, heurdata->exp3, priorities, (unsigned int)(heurdata->seed + SCIPgetNVars(scip))) );
      SCIP_CALL( SCIPresetBandit(scip, heurdata->epsgreedynh, priorities, (unsigned int)(heurdata->seed + SCIPgetNVars(scip))) );
   }

   heurdata->usednodes = 0;
   heurdata->ninitneighborhoods = heurdata->nactiveneighborhoods;
   resetMinimumImprovement(heurdata);
   resetTargetNodeLimit(heurdata);
   resetCurrentNeighborhood(heurdata);

   SCIPfreeBufferArray(scip, &priorities);

   /* open gain file for reading */
   if( strncasecmp(heurdata->gainfilename, DEFAULT_GAINFILENAME, strlen(DEFAULT_GAINFILENAME) != 0 ) )
   {
      heurdata->gainfile = fopen(heurdata->gainfilename, "w");

      if( heurdata->gainfile == NULL )
      {
         SCIPerrorMessage("Error: Could not open gain file <%s>\n", heurdata->gainfilename);
         return SCIP_FILECREATEERROR;
      }

      SCIPdebugMsg(scip, "Writing gain information to <%s>\n", heurdata->gainfilename);
   }
   else
      heurdata->gainfile = NULL;

   return SCIP_OKAY;
}

/** todo print neighborhood statistics */
static
void printNeighborhoodStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   SCIP_Real* epsgreedyweights;
   int i;
   int j;
   HISTINDEX statusses[] = {HIDX_OPT, HIDX_INFEAS, HIDX_NODELIM, HIDX_STALLNODE, HIDX_SOLLIM, HIDX_USR, HIDX_OTHER};

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Neighborhoods      :%11s %11s %11s %11s %11s %11s %11s %11s %11s %4s %4s %4s %4s %4s %4s %4s\n",
            "Calls", "SetupTime", "SubmipTime", "SubmipNodes", "Sols", "Best", "Exp3", "EpsGreedy", "TgtFixRate",
            "Opt", "Inf", "Node", "Stal", "Sol", "Usr", "Othr");

   epsgreedyweights = SCIPgetWeightsEpsgreedy(heurdata->epsgreedynh);
   /* todo loop over neighborhoods and fill in statistics */
   for( i = 0; i < heurdata->nneighborhoods; ++i )
   {
      NH* neighborhood;
      SCIP_Real proba;
      neighborhood = heurdata->neighborhoods[i];
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "  %-17s:", neighborhood->name);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%11d", neighborhood->stats.nruns);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " %11.2f", SCIPgetClockTime(scip, neighborhood->stats.setupclock) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " %11.2f", SCIPgetClockTime(scip, neighborhood->stats.submipclock) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " %11" SCIP_LONGINT_FORMAT, neighborhood->stats.usednodes );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " %11d", neighborhood->stats.nsolsfound);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " %11d", neighborhood->stats.nbestsolsfound);

      if( neighborhood->active )
         proba = SCIPgetProbabilityExp3(heurdata->exp3, i);
      else
         proba = 0.0;
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " %11.5f", proba);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " %11.5f", epsgreedyweights[i]);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " %11.3f", neighborhood->fixingrate.targetfixingrate);

      /* loop over status histogram */
      for( j = 0; j < NHISTENTRIES; ++j )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " %4d", neighborhood->stats.statushist[statusses[j]]);

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "\n");
   }
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitLns)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   printNeighborhoodStatistics(scip, heurdata);


   /* free neighborhood specific data */
   for( i = 0; i < heurdata->nneighborhoods; ++i )
   {
      NH* neighborhood = heurdata->neighborhoods[i];

      SCIP_CALL( neighborhoodExit(scip, neighborhood) );
   }

   if( heurdata->gainfile != NULL )
   {
      fclose(heurdata->gainfile);
      heurdata->gainfile = NULL;
   }

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeLns)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* an exp3 is only initialized when a problem is read */
   if( heurdata->exp3 != NULL )
   {
      assert(heurdata->epsgreedynh != NULL);
      SCIPfreeBandit(scip, &heurdata->exp3);
      SCIP_CALL( SCIPfreeBandit(scip, &heurdata->epsgreedynh) );
   }

   /* free neighborhoods */
   for( i = 0; i < heurdata->nneighborhoods; ++i )
   {
      SCIP_CALL( lnsFreeNeighborhood(scip, &(heurdata->neighborhoods[i])) );
   }

   SCIPfreeBlockMemoryArray(scip, &heurdata->neighborhoods, NNEIGHBORHOODS);

   SCIPfreeBlockMemory(scip, &heurdata);

   return SCIP_OKAY;
}

/** creates the lns primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurLns(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create lns primal heuristic data */
   heurdata = NULL;

   heur = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* author bzfhende: TODO make this a user parameter? */
   heurdata->lplimfac = LPLIMFAC;
   heurdata->epsgreedynh = NULL;
   heurdata->exp3 = NULL;
   heurdata->gainfilename = NULL;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->neighborhoods, NNEIGHBORHOODS) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecLns, heurdata) );

   assert(heur != NULL);

   /* include all neighborhoods */
   SCIP_CALL( includeNeighborhoods(scip, heurdata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyLns) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeLns) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitLns) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitLns) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolLns) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolLns) );

   /* add lns primal heuristic parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes,  TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "offset added to the nodes budget",
         &heurdata->nodesoffset, FALSE, DEFAULT_NODESOFFSET, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start a sub-SCIP",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "fraction of nodes compared to the main SCIP for budget computation",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/startminimprove",
         "initial factor by which LNS should at least improve the incumbent",
         &heurdata->startminimprove, TRUE, DEFAULT_STARTMINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprovelow",
         "lower threshold for the minimal improvement over the incumbent",
         &heurdata->minimprovelow, TRUE, DEFAULT_MINIMPROVELOW, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprovehigh",
         "upper bound for the minimal improvement over the incumbent",
         &heurdata->minimprovehigh, TRUE, DEFAULT_MINIMPROVEHIGH, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nsolslim",
         "limit on the number of improving solutions in a sub-SCIP call",
         &heurdata->nsolslim, FALSE, DEFAULT_NSOLSLIM, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/banditalgo",
         "the bandit algorithm: (u)pper confidence bounds, (e)xp.3, epsilon (g)reedy",
         &heurdata->banditalgo, TRUE, DEFAULT_BANDITALGO, "ueg", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/gamma",
         "weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution for exp3",
         &heurdata->gamma, TRUE, DEFAULT_GAMMA, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/beta",
         "gain offset between 0 and 1 at every observation for exp3",
         &heurdata->beta, TRUE, DEFAULT_BETA, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/usedistances",
         "distances from fixed variables be used for variable prioritization",
         &heurdata->usedistances, TRUE, DEFAULT_USEDISTANCES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useredcost",
         "should reduced cost scores be used for variable prioritization?",
         &heurdata->useredcost, TRUE, DEFAULT_USEREDCOST, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/domorefixings",
         "should the LNS heuristic do more fixings by itself based on variable prioritization"
         "until the target fixing is reached?",
         &heurdata->domorefixings, TRUE, DEFAULT_DOMOREFIXINGS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/adjustfixingrate",
         "should the heuristic adjust the target fixing rate based on the success?",
         &heurdata->adjustfixingrate, TRUE, DEFAULT_ADJUSTFIXINGRATE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/usesubscipheurs",
         "should the heuristic activate other sub-SCIP heuristics during its search?",
         &heurdata->usesubscipheurs, TRUE, DEFAULT_USESUBSCIPHEURS, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/gainmeasure",
         "measure for the gain of a neighborhood? 'b'oolean, 'w'eighted boolean, 'g'ap based?",
         &heurdata->gainmeasure, TRUE, DEFAULT_GAINMEASURE, GAINMEASURES, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/seed",
         "initial random seed for bandit algorithms and random decisions by neighborhoods",
         &heurdata->seed, FALSE, DEFAULT_SEED, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/adjustminimprove",
         "should the factor by which the minimum improvement is bound be dynamically updated?",
         &heurdata->adjustminimprove, TRUE, DEFAULT_ADJUSTMINIMPROVE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/eps",
            "probability for exploration in epsilon-greedy bandit algorithm",
            &heurdata->eps, TRUE, DEFAULT_EPS, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/resetweights",
            "should the bandit algorithms be reset when a new problem is read?",
            &heurdata->resetweights, TRUE, DEFAULT_RESETWEIGHTS, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "heuristics/" HEUR_NAME "/gainfilename", "file name to store all gains and the selection of the bandit",
         &heurdata->gainfilename, TRUE, DEFAULT_GAINFILENAME, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/subsciprandseeds",
            "should random seeds of sub-SCIPs be altered to increase diversification?",
            &heurdata->subsciprandseeds, TRUE, DEFAULT_SUBSCIPRANDSEEDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/normbyeffort",
            "should the gain be normalized by the effort?",
            &heurdata->normbyeffort, TRUE, DEFAULT_NORMBYEFFORT, NULL, NULL) );

   return SCIP_OKAY;
}
