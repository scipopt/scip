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
 * @brief  lns primal heuristic
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/heur_lns.h"

#define HEUR_NAME             "lns"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         'L'
#define HEUR_PRIORITY         -1000000
#define HEUR_FREQ             20
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define NNEIGHBORHOODS 3
#define DEFAULT_NODESQUOT 0.05
#define DEFAULT_NODESOFFSET 500
#define DEFAULT_NSOLSLIM 3
#define DEFAULT_MINNODES 10LL
#define DEFAULT_MINIMPROVE 0.02
#define DEFAULT_MAXNODES 5000
#define LPLIMFAC 2.0
#define DEFAULT_INITSEED 113
#define MUTATIONSEED 121
#define DEFAULT_BESTSOLWEIGHT 3
#define DEFAULT_BANDITALGO 'e' /**< the default bandit algorithm: (u)pper confidence bounds, (e)xp.3, epsilon (g)reedy */
#define DEFAULT_GAMMA 0.0      /**< default weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution for exp3 */
#define DEFAULT_BETA 0.0       /**< default gain offset between 0 and 1 at every observation for exp3 */

#define DEFAULT_MINFIXINGRATE_RENS 0.3
#define DEFAULT_MAXFIXINGRATE_RENS 0.7
#define DEFAULT_MINFIXINGRATE_RINS 0.1
#define DEFAULT_MAXFIXINGRATE_RINS 0.5
#define DEFAULT_MINFIXINGRATE_MUTATION 0.4
#define DEFAULT_MAXFIXINGRATE_MUTATION 0.9
#define DEFAULT_MINFIXINGRATE_GINS 0.3
#define DEFAULT_MAXFIXINGRATE_GINS 0.5
#define DEFAULT_MINFIXINGRATE_LOCALBRANCHING 0.0
#define DEFAULT_MAXFIXINGRATE_LOCALBRANCHING 0.0
#define DEFAULT_MINFIXINGRATE_PROXIMITY 0.0
#define DEFAULT_MAXFIXINGRATE_PROXIMITY 0.0
#define DEFAULT_MINFIXINGRATE_CROSSOVER 0.4
#define DEFAULT_MAXFIXINGRATE_CROSSOVER 0.8

/* event handler properties */
#define EVENTHDLR_NAME         "Lns"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"
#define SCIP_EVENTTYPE_LNS (SCIP_EVENTTYPE_LPSOLVED | SCIP_EVENTTYPE_BESTSOLFOUND)

/*
 * Data structures
 */

/* additional crossover data structure */
typedef struct data_crossover DATA_CROSSOVER;
typedef struct data_mutation DATA_MUTATION;
typedef struct EpsGreedy EPSGREEDY;
typedef struct expThree EXPTHREE;
typedef struct NH_FixingRate NH_FIXINGRATE;
typedef struct NH_Stats NH_STATS;
typedef struct Nh NH;

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
   SCIP_Bool*            success             \
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

/** callback function for special sub-SCIP settings
 *
 * todo comments
 */
#define DECL_SETUPSUBSCIP(x) SCIP_RETCODE x (\
   SCIP*                 sourcescip,         \
   SCIP*                 targetscip          \
)

/** statistics for a neighborhood */
struct NH_Stats
{
   SCIP_Longint          usednodes;
   SCIP_Longint          lpiterations;
   SCIP_Real             totalgapclosed;
   int                   nruns;
   int                   nrunsbestsol;
   int                   nsolsfound;
   int                   nbestsolsfound;
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
   SCIP_Bool             active;             /**< is this neighborhood active or not? */
   union
   {
      DATA_MUTATION*     mutation;           /**< mutation data */
      DATA_CROSSOVER*    crossover;          /**< crossover data */
   }                     data;               /**< data object for neighborhood specific data */
};

/** mutation neighborhood data structure */
struct data_mutation
{
   SCIP_RANDNUMGEN*      rng;
};

/** current, unnormalized reward for item i */
#define DECL_EPSREWARD(x) SCIP_Real x ( \
   SCIP*                 scip,             \
   EPSGREEDY*            epsgreedy,        \
   int                   i                 \
   )

/** callback for the number of choices */
#define DECL_EPSNCHOICES(x) int x (        \
   SCIP*                 scip,             \
   EPSGREEDY*            epsgreedy         \
   )

/** structure that represents the adversarial bandit algorithm exp3 */
struct expThree
{
   int                   nactions;           /**< the number of actions to select from */
   int                   ndraws;             /**< the total number of draws for all arms */
   SCIP_Real*            weights;            /**< exponential weight for each arm */
   SCIP_Real*            cumulativegain;     /**< the cumulative gain for each arm */
   SCIP_Real             weightsum;          /**< the sum of all weights */
   SCIP_Real             gamma;              /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution */
   SCIP_Real             beta;               /**< gain offset between 0 and 1 at every observation */
   SCIP_RANDNUMGEN*      rng;                /**< random number generator */
};

struct EpsGreedy
{
   SCIP_Real             eps;                /**< epsilon parameter (between 0 and 1) to control epsilon greedy */
   SCIP_RANDNUMGEN*      rng;                /**< random number generator for randomized selection of routines  */
   DECL_EPSREWARD        ((*epsreward));     /**< reward callback for unnormalized reward of the i'th item */
   DECL_EPSNCHOICES      ((*epsnchoices));   /**< callback for the number of choices */
   void*                 userptr;            /**< user pointer for callback functions */
};

/** primal heuristic data */
struct SCIP_HeurData
{
   NH**                  neighborhoods;      /**< array of neighborhoods with the best one at the first position */
   char                  banditalgo;         /**< the bandit algorithm: (u)pper confidence bounds, (e)xp.3, epsilon (g)reedy */
   EPSGREEDY*            epsgreedynh;        /**< epsilon greedy selector for a neighborhood */
   EXPTHREE*             exp3;               /**< exp3 bandit algorithm */
   EPSGREEDY*            epsgreedyfilter;    /**< epsilon greedy selector for a filter strategy */
   SCIP_Longint          nodesoffset;        /**< offset added to the nodes budget */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes in a single sub-SCIP */
   SCIP_Longint          minnodes;           /**< minimum number of nodes required to start a sub-SCIP */
   SCIP_Longint          usednodes;          /**< total number of nodes already spent in sub-SCIPs */
   SCIP_Real             nodesquot;          /**< fraction of nodes compared to the main SCIP for budget computation */
   SCIP_Real             minimprove;         /**< factor by which LNS should at least improve the incumbent */
   SCIP_Real             lplimfac;           /**< limit fraction of LPs per node to interrupt sub-SCIP */
   SCIP_Real             gamma;              /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution for exp3 */
   SCIP_Real             beta;               /**< gain offset between 0 and 1 at every observation for exp3 */
   int                   nneighborhoods;     /**< number of neighborhoods */
   int                   nsolslim;           /**< limit on the number of improving solutions in a sub-SCIP call */
};

/** event handler data */
struct SCIP_EventData
{
   SCIP_VAR**            subvars;            /**< the variables of the subproblem */
   SCIP*                 sourcescip;         /**< original SCIP data structure */
   SCIP_HEUR*            heur;               /**< lns heuristic structure */
   SCIP_Longint          nodelimit;
   SCIP_Real             lplimfac;           /**< limit fraction of LPs per node to interrupt sub-SCIP */
};

/** represents limits for the sub-SCIP solving process */
struct SolveLimits
{
   SCIP_Longint          nodelimit;          /**< maximum number of solving nodes for the sub-SCIP */
   SCIP_Real             memorylimit;        /**< memory limit for the sub-SCIP */
   SCIP_Real             timelimit;          /**< time limit for the sub-SCIP */
};

typedef struct SolveLimits SOLVELIMITS;


/*
 * Local methods
 */

/** reset target fixing rate
 *
 *  todo use neighborhood specific values
 */
static
SCIP_RETCODE fixingRateReset(
   SCIP*                scip,                /**< SCIP data structure */
   NH_FIXINGRATE*       fixingrate           /**< heuristic fixing rate */
   )
{
   assert(scip != NULL);
   assert(fixingrate != NULL);

   /* use the middle between the minimum and the maximum fixing rate */
   fixingrate->targetfixingrate = 0.5 * (fixingrate->minfixingrate + fixingrate->maxfixingrate);

   return SCIP_OKAY;
}

/** todo reset neighborhood statistics */
static
void neighborhoodStatsReset(
   SCIP*                scip,                /**< SCIP data structure */
   NH_STATS*            stats                /**< neighborhood statistics */
   )
{
   assert(scip != NULL);
   assert(stats != NULL);

   stats->lpiterations = 0l;
   stats->nbestsolsfound = 0;
   stats->nruns = 0;
   stats->nrunsbestsol = 0;
   stats->nsolsfound = 0;
   stats->presolrounds = 0;
   stats->totalgapclosed = 0;
   stats->totalnbinfixings = 0;
   stats->totalncontfixings = 0;
   stats->totalnimplintfixings = 0;
   stats->totalnfixings = 0;
   stats->usednodes = 0l;
}

/** increase fixing rate */
static
void increaseFixingRate(
   NH_FIXINGRATE*        fx                  /**< fixing rate */
   )
{
   fx->targetfixingrate = 0.5 * (fx->targetfixingrate + fx->maxfixingrate);
}

/** todo decerase fixing rate */
static
void decreaseFixingRate(
   NH_FIXINGRATE*        fx                  /**< fixing rate */
   )
{
   fx->targetfixingrate = 0.5 * (fx->targetfixingrate + fx->minfixingrate);
}

/** todo update fixing rate */
static
SCIP_RETCODE updateFixingRate(
   SCIP*                scip,               /**< SCIP data structure */
   NH*                  neighborhood,       /**< neighborhood */
   SCIP_STATUS          subscipstatus,      /**< status of the sub-SCIP run */
   NH_STATS*            runstats            /**< run statistics for this run */
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
         increaseFixingRate(fx);
         break;
      default:
         break;
   }

   return SCIP_OKAY;
}

/** todo create a neighborhood of the specified name */
static
SCIP_RETCODE lnsIncludeNeighborhood(
   SCIP*                scip,                /**< SCIP data structure */
   SCIP_HEURDATA*       heurdata,            /**< heuristic data of the LNS heuristic */
   NH**                 neighborhood,        /**< neighborhood that should be created and included */
   const char*          name,                /**< name to distinguish this neighborhood */
   SCIP_Real            minfixingrate,       /**< default value for minfixingrate parameter of this neighborhood */
   SCIP_Real            maxfixingrate,       /**< default value for maxfixingrate parameter of this neighborhood */
   DECL_VARFIXINGS      ((*varfixings)),     /**< variable fixing callback for this neighborhood, or NULL */
   DECL_CHANGESUBSCIP   ((*changesubscip)),  /**< subscip changes callback for this neighborhood, or NULL */
   DECL_SETUPSUBSCIP    ((*setupsubscip)),   /**< setup callback for this neighborhood, or NULL */
   DECL_NHINIT          ((*nhinit)),         /**< initialization callback for neighborhood, or NULL */
   DECL_NHEXIT          ((*nhexit))          /**< deinitialization callback for neighborhood, or NULL */
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

   (*neighborhood)->changesubscip = changesubscip;
   (*neighborhood)->varfixings = varfixings;
   (*neighborhood)->setupsubscip = setupsubscip;
   (*neighborhood)->nhinit = nhinit;
   (*neighborhood)->nhexit = nhexit;

   sprintf(paramname, "heuristics/lns/%s/minfixingrate", name);
   /* author bzfhende: TODO add parameters for this neighborhood */
   SCIP_CALL( SCIPaddRealParam(scip, paramname, "minimum fixing rate for this neighborhood",
         &(*neighborhood)->fixingrate.minfixingrate, TRUE, minfixingrate, 0.0, 1.0, NULL, NULL) );
   sprintf(paramname, "heuristics/lns/%s/maxfixingrate", name);
      /* author bzfhende: TODO add parameters for this neighborhood */
   SCIP_CALL( SCIPaddRealParam(scip, paramname, "maximum fixing rate for this neighborhood",
         &(*neighborhood)->fixingrate.maxfixingrate, TRUE, maxfixingrate, 0.0, 1.0, NULL, NULL) );

   heurdata->neighborhoods[heurdata->nneighborhoods++] = (*neighborhood);

   return SCIP_OKAY;
}

/** release all data and free a neighborhood */
static
void lnsFreeNeighborhood(
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

   SCIPfreeBlockMemory(scip, neighborhood);
   *neighborhood = NULL;
}

/** todo initialize neighborhood specific data */
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

/** todo deinitialize neighborhood specific data */
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

/** create an epsilon greedy selector with the necessary callbacks */
static
SCIP_RETCODE epsGreedyCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   EPSGREEDY**           epsgreedy,          /**< pointer to store the epsilon greedy selector */
   unsigned int          initseed,           /**< initial seed for random number generation */
   void*                 userptr,            /**< user pointer that may be used inside the callbacks, or NULL, if not needed */
   DECL_EPSNCHOICES      ((*epsnchoices)),   /**< callback for the number of choices */
   DECL_EPSREWARD        ((*epsreward))      /**< reward callback for this epsilon greedy selector */
   )
{
   assert(scip != NULL);
   assert(epsgreedy != NULL);
   assert(epsnchoices != NULL);
   assert(epsreward != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, epsgreedy) );

   /* create random number generator for the selector */
   SCIP_CALL( SCIPrandomCreate(&(*epsgreedy)->rng, SCIPblkmem(scip), initseed) );

   (*epsgreedy)->epsnchoices = epsnchoices;
   (*epsgreedy)->epsreward = epsreward;
   (*epsgreedy)->userptr = userptr;

   return SCIP_OKAY;
}

/** free an epsilon greedy selector */
static
void epsGreedyFree(
   SCIP*                 scip,               /**< SCIP data structure */
   EPSGREEDY**           epsgreedy           /**< epsilon greedy selector pointer */
   )
{
   assert(scip != NULL);
   assert(epsgreedy != NULL);
   assert(*epsgreedy != NULL);

   SCIPrandomFree(&(*epsgreedy)->rng);

   SCIPfreeBlockMemory(scip, epsgreedy);
   *epsgreedy = NULL;
}

/** todo reset an epsilon greedy selector */
static
SCIP_RETCODE epsGreedyReset(
   SCIP*                scip                /**< SCIP data structure */
   )
{
   return SCIP_OKAY;
}

/** todo let the epsilon greedy selector choose its next move */
static
SCIP_RETCODE epsGreedySelect(
   SCIP*                scip,               /**< SCIP data structure */
   EPSGREEDY*           epsgreedy,          /**< epsilon greedy selector */
   int*                 i                   /**< pointer to store the selection (will be set to -1 if no choice is available) */
   )
{
   int nchoices;
   SCIP_Real rand;

   assert(i != NULL);
   assert(scip != NULL);
   assert(epsgreedy != NULL);

   nchoices = epsgreedy->epsnchoices(scip, epsgreedy);

   if( nchoices == 0 )
   {
      *i = -1;
      return SCIP_OKAY;
   }
   /** roll the dice to check if the best element should be picked, or an element at random */
   rand = SCIPrandomGetReal(epsgreedy->rng, 0.0, 1.0);

   /* todo make epsilon decrease over time */
   if( rand <= epsgreedy->eps )
   {
      /** pick the element with the largest reward */
      int j;
      SCIP_Real maxreward = epsgreedy->epsreward(scip, epsgreedy, 0);
      *i = 0;
      /* determine reward for every element */
      for( j = 1; j < nchoices; ++j )
      {
         SCIP_Real reward = epsgreedy->epsreward(scip, epsgreedy, j);

         if( maxreward < reward )
         {
            *i = j;
            maxreward = reward;
         }
      }
   }
   else
   {
      /* play one of the other arms at random */
      *i = SCIPrandomGetInt(epsgreedy->rng, 0, nchoices - 1);
   }

   return SCIP_OKAY;
}

static
DECL_EPSNCHOICES(epsNChoicesLns)
{
   SCIP_HEURDATA* heurdata = (SCIP_HEURDATA*)(epsgreedy->userptr);
   return heurdata->nneighborhoods;
}

static
DECL_EPSREWARD(epsRewardLns)
{
   SCIP_Real numerator;
   SCIP_Real denominator;
   SCIP_HEURDATA* heurdata = (SCIP_HEURDATA*)(epsgreedy->userptr);
   NH* neighborhood = heurdata->neighborhoods[i];
   assert(i < heurdata->nneighborhoods);

   numerator = neighborhood->stats.nrunsbestsol;
   denominator = MAX(1.0, neighborhood->stats.nruns * DEFAULT_BESTSOLWEIGHT);
   return numerator / denominator;
}


/** todo reset an exp3 bandit algorithm */
static
SCIP_RETCODE expThreeReset(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPTHREE*             exp3,               /**< the exp3 algorithm */
   unsigned int          randseed            /**< initial random seed */
   )
{
   int i;
   assert(exp3->nactions > 0);

   exp3->ndraws = 0;

   /* initialize weights */
   for( i = 0; i < exp3->nactions; ++i )
   {
      exp3->weights[i] = 1.0;
      exp3->cumulativegain[i] = 0.0;
   }
   exp3->weightsum = (SCIP_Real)exp3->nactions;


   /* author bzfhende: TODO reset random number generator */
   if( exp3->rng != NULL )
   {
      SCIPrandomFree(&exp3->rng);
   }

   SCIP_CALL( SCIPrandomCreate(&exp3->rng, SCIPblkmem(scip), SCIPinitializeRandomSeed(scip, randseed)) );

   return SCIP_OKAY;
}

/** create an exp3 bandit algorithm */
static
SCIP_RETCODE expThreeCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPTHREE**            exp3,               /**< pointer to store the exp3 algorithm */
   unsigned int          initseed,           /**< initial seed for the exp3 algorithm */
   int                   nactions,           /**< the number of actions */
   SCIP_Real             gamma,              /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution */
   SCIP_Real             beta                /**< gain offset between 0 and 1 at every observation */
   )
{
   assert(scip != NULL);
   assert(exp3 != NULL);
   assert(nactions > 0);

   SCIP_CALL( SCIPallocBlockMemory(scip, exp3) );

   (*exp3)->nactions = nactions;
   (*exp3)->rng = NULL;
   (*exp3)->gamma = gamma;
   (*exp3)->beta = beta;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*exp3)->weights, nactions) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*exp3)->cumulativegain, nactions) );

   SCIP_CALL( expThreeReset(scip, *exp3, initseed) );

   return SCIP_OKAY;
}

/** free an exp3 bandit algorithm */
static
void expThreeFree(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPTHREE**            exp3                /**< pointer to the exp3 algorithm */
   )
{
   assert(scip != NULL);
   assert(exp3 != NULL);
   assert(*exp3 != NULL);

   SCIPfreeBlockMemoryArray(scip, &(*exp3)->weights, (*exp3)->nactions);
   SCIPfreeBlockMemoryArray(scip, &(*exp3)->cumulativegain, (*exp3)->nactions);

   SCIPrandomFree(&(*exp3)->rng);

   SCIPfreeBlockMemory(scip, exp3);
}

/** draw the next action from the current probability distribution over the arms */
static
SCIP_RETCODE expThreeSelectAction(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPTHREE*             exp3,               /**< exp3 bandit algorithm */
   int*                  action              /**< the index of the selected action, between 0 and nactions - 1 */
   )
{
   SCIP_Real rand;
   SCIP_Real psum;
   SCIP_Real gammaoverk;
   SCIP_Real oneminusgamma;
   SCIP_Real* weights;
   SCIP_Real weightsum;
   int i;

   assert(scip != NULL);
   assert(exp3 != NULL);
   assert(action != NULL);

   /* draw a random number between 0 and 1 */
   rand = SCIPrandomGetReal(exp3->rng, 0, 1);

   /* initialize some local variables to speed up probability computations */
   oneminusgamma = 1 - exp3->gamma;
   gammaoverk = exp3->gamma / (SCIP_Real)exp3->nactions;
   weightsum = exp3->weightsum;
   weights = exp3->weights;
   psum = 0.0;

   /* loop over probability distribution until rand is reached
    * the loop terminates without looking at the last action,
    * which is then selected automatically if the target probability
    * is not reached earlier
    */
   for( i = 0; i < exp3->nactions - 1; ++i )
   {
      SCIP_Real prob;

      /* compute the probability for arm i as convex kombination of a uniform distribution and a weighted distribution */
      prob = oneminusgamma * weights[i] / weightsum + gammaoverk;
      psum += prob;

      /* break and select element if target probability is reached */
      if( rand <= psum )
         break;
   }

   *action = i;
   exp3->ndraws++;

   return SCIP_OKAY;
}

/** todo update the exp3 probability distribution after observing a gain */
static
SCIP_RETCODE expThreeUpdate(
   SCIP*                scip,               /**< SCIP data structure */
   EXPTHREE*            exp3,               /**< exp3 bandit algorithm */
   SCIP_Real            gain,               /**< the gain that has been observed for action i */
   int                  i                   /**< the last selected action, for which the gain has been observed */
   )
{
   SCIP_Real eta;
   SCIP_Real gainestim;
   SCIP_Real beta;
   SCIP_Real weightsum;
   SCIP_Real newweightsum;
   SCIP_Real* weights;
   SCIP_Real oneminusgamma;
   SCIP_Real gammaoverk;

   assert(scip != NULL);
   assert(exp3 != NULL);
   assert(i >= 0);
   assert(i < exp3->nactions);
   assert(exp3->ndraws > 0);

   /* the learning rate eta */
   eta = 1.0 / (SCIP_Real)exp3->nactions;

   beta = exp3->beta;
   oneminusgamma = 1 - exp3->gamma;
   gammaoverk = exp3->gamma * eta;
   weights = exp3->weights;
   weightsum = exp3->weightsum;
   newweightsum = weightsum;

   /* if beta is zero, only the observation for the current arm needs an update */
   if( SCIPisZero(scip, beta) )
   {
      SCIP_Real probai;
      probai = oneminusgamma * weights[i] / weightsum + gammaoverk;
      gainestim = gain / probai;
      newweightsum -= weights[i];
      weights[i] *= exp(eta * gainestim);
      newweightsum += weights[i];
   }
   else
   {
      int j;
      /* loop over all items and update their weights based on the influence of the beta parameter */
      for( j = 0; j < exp3->nactions; ++j )
      {
         SCIP_Real probaj;
         probaj = oneminusgamma * weights[j] / weightsum + gammaoverk;

         /* consider the gain only for the chosen arm i, use constant beta offset otherwise */
         if( j == i )
            gainestim = (gain + beta) / probaj;
         else
            gainestim = beta / probaj;

         newweightsum -= weights[j];
         weights[j] *= exp(eta * gainestim);
         newweightsum += weights[j];
      }
   }

   exp3->weightsum = newweightsum;

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

   assert(subscip != NULL);

   subsol = SCIPgetBestSol(subscip);
   assert(subsol != NULL);

   sourcescip = eventdata->sourcescip;
   subvars = eventdata->subvars;
   heur = eventdata->heur;
   assert(sourcescip != NULL);
   assert(sourcescip != subscip);
   assert(heur != NULL);
   assert(subvars != NULL);

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

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(sourcescip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );

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
   stats->nbestsolsfound = -SCIPgetNBestSolsFound(scip);
   stats->nsolsfound = -SCIPgetNSolsFound(scip);
   stats->lpiterations = 0l;
   stats->usednodes = 0l;
}

/** update run stats after the sub SCIP was solved */
static
void updateRunStats(
   SCIP*                 scip,               /**< SCIP data structure */
   NH_STATS*             stats,              /**< run statistics */
   SCIP*                 subscip             /**< sub-SCIP instance, or NULL */
   )
{
   stats->nbestsolsfound += SCIPgetNBestSolsFound(scip);
   stats->nsolsfound += SCIPgetNSolsFound(scip);
   stats->lpiterations = subscip != NULL ? SCIPgetNLPIterations(subscip) : 0l;
   stats->usednodes = subscip != NULL ? SCIPgetNNodes(subscip) : 0l;
}

/** update the statistics of the neighborhood based on the sub-SCIP run */
static
void updateNeighborhoodStats(
   SCIP*                 scip,               /**< SCIP data structure */
   NH_STATS*             runstats,           /**< run statistics */
   NH*                   neighborhood        /**< the selected neighborhood */
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

/** todo call variable fixing callback for this neighborhood */
static
SCIP_RETCODE neighborhoodFixVariables(
  SCIP*                  scip,
  NH*                    neighborhood,
  SCIP_VAR**             varbuf,
  SCIP_Real*             valbuf,
  int*                   nfixings,
  SCIP_Bool*             success
  )
{

   assert(scip != NULL);
   assert(neighborhood != NULL);
   assert(varbuf != NULL);
   assert(valbuf != NULL);
   assert(nfixings != NULL);
   assert(success != NULL);

   *nfixings = 0;

   *success = FALSE;
   if( neighborhood->varfixings != NULL )
   {
      SCIP_CALL( neighborhood->varfixings(scip, neighborhood, varbuf, valbuf, nfixings, success) );
   }
   else
   {
      *success = TRUE;
   }

   /** todo if too few fixings, use a strategy to select more variable fixings: randomized, LP graph, ReducedCost/Ps-Cost based, mix */

   if( *success && (*nfixings >= neighborhood->fixingrate.targetfixingrate * (SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip))) )
      *success = TRUE;

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

   /* calculate the maximal number of search nodes until heuristic is aborted */
   solvelimits->nodelimit = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));
   solvelimits->nodelimit += heurdata->nodesoffset;
   solvelimits->nodelimit -= heurdata->usednodes;
   solvelimits->nodelimit -= 100 * SCIPheurGetNCalls(heur);

   /* check whether we have enough nodes left to call subproblem solving */
   if( solvelimits->nodelimit < heurdata->minnodes )
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
         SCIP_CALL( epsGreedySelect(scip, heurdata->epsgreedynh, neighborhoodidx) );
         break;
      case 'e':
         SCIP_CALL( expThreeSelectAction(scip, heurdata->exp3, neighborhoodidx) );
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

/** update internal bandit algorithm statistics for future draws */
static
SCIP_RETCODE updateBanditAlgorithms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the LNS neighborhood */
   NH_STATS*             runstats,           /**< run statistics */
   int                   neighborhoodidx     /**< the neighborhood that was chosen */
   )
{
   SCIP_Real gain;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(runstats != NULL);
   assert(neighborhoodidx >= 0);
   assert(neighborhoodidx < heurdata->nneighborhoods);

   /* compute the gain for this run based on the runstats */
   if( runstats->nbestsolsfound > 0 )
      gain = 1.0;
   else if( runstats->nsolsfound > 0 )
      gain = 1.0 / DEFAULT_BESTSOLWEIGHT;
   else
      gain = 0.0;

   SCIP_CALL( expThreeUpdate(scip, heurdata->exp3, gain, neighborhoodidx) );

   return SCIP_OKAY;
}

/** set up the sub-SCIP */
static
SCIP_RETCODE setupSubScip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SOLVELIMITS*          solvelimits,        /**< pointer to solving limits data structure */
   SCIP_HEUR*            heur                /**< this heuristic */
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
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
#endif

   /* disable statistic timing inside sub SCIP */
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );

   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", heurdata->nsolslim) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

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
      SCIP_CALL(SCIPsetObjlimit(subscip, cutoff));
   }

   /* set solve limits for sub-SCIP */
   SCIP_CALL( setLimits(subscip, solvelimits) );

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
   NH_STATS runstats;
   SCIP_STATUS subscipstatus;

   int nfixings;
   int nvars;
   int neighborhoodidx;
   NH* neighborhood;
   SOLVELIMITS solvelimits;
   SCIP_Bool success;
   SCIP_Bool run;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DELAYED;

   run = TRUE;
   /** check if budget allows a run of the next selected neighborhood */
   SCIP_CALL( determineLimits(scip, heur, &solvelimits, &run) );
   SCIPdebugMsg(scip, "Budget check: %s\n", run ? "passed" : "must wait");

   if( !run )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /** todo select the next neighborhood to run based on the selected bandit algorithm */
   SCIP_CALL( selectNeighborhood(scip, heurdata, &neighborhoodidx) );
   assert(neighborhoodidx >= 0);
   assert(heurdata->nneighborhoods > neighborhoodidx);
   neighborhood = heurdata->neighborhoods[neighborhoodidx];
   SCIPdebugMsg(scip, "Selected '%s' neighborhood %d\n", neighborhood->name, neighborhoodidx);

   /** todo check if the selected neighborhood can run at this node */
   SCIP_CALL( SCIPallocBufferArray(scip, &varbuf, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &valbuf, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   /* initialize neighborhood statistics for a run */
   initRunStats(scip, &runstats);
   do
   {
      SCIP* subscip = NULL;
      SCIP_HASHMAP* varmapf;
      SCIP_EVENTHDLR* eventhdlr;
      SCIP_EVENTDATA eventdata;
      int v;

      subscipstatus = SCIP_STATUS_INFEASIBLE;
      /** determine variable fixings and objective coefficients of this neighborhood */
      SCIP_CALL( neighborhoodFixVariables(scip, neighborhood, varbuf, valbuf, &nfixings, &success) );

      SCIPdebugMsg(scip, "Fix %d/%d variables\n", nfixings, nvars);

      /* fixing was not successful */
      if( ! success )
      {
         updateRunStats(scip, &runstats, subscip);
         break;
      }


      *result = SCIP_DIDNOTFIND;

      neighborhood->stats.totalnfixings += nfixings;


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

      /** todo let the neighborhood add additional constraints, or restrict domains */
      SCIP_CALL( setupSubScip(scip, subscip, &solvelimits, heur) );

      /* copy the necessary data into the event data to create new solutions */
      eventdata.nodelimit = solvelimits.nodelimit;
      eventdata.lplimfac = heurdata->lplimfac;
      eventdata.heur = heur;
      eventdata.sourcescip = scip;
      eventdata.subvars = subvars;

      /* include an event handler to transfer solutions into the main SCIP */
      SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecLns, NULL) );

      /*transform the problem before catching the events */
      SCIP_CALL( SCIPtransformProb(subscip) );
      SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LNS, eventhdlr, &eventdata, NULL) );


      /** todo alternatively: set up sub-SCIP and run presolving */
      /** todo was presolving successful enough regarding fixings? otherwise terminate */

      /* run sub-SCIP for the given budget, and collect statistics */
      SCIP_CALL_ABORT( SCIPsolve(subscip) );

      /* update statistics based on the sub-SCIP run */
      updateRunStats(scip, &runstats, subscip);
      subscipstatus = SCIPgetStatus(subscip);
      heurdata->usednodes += SCIPgetNNodes(subscip);

      SCIP_CALL( SCIPfree(&subscip) );
   }
   while( FALSE );

   /** todo determine the success of this neighborhood, and the target fixing rate for the next time */
   updateNeighborhoodStats(scip, &runstats, neighborhood);
   updateFixingRate(scip, neighborhood, subscipstatus, &runstats);
   SCIP_CALL( updateBanditAlgorithms(scip, heurdata, &runstats, neighborhoodidx) );

   SCIPfreeBufferArray(scip, &subvars);
   SCIPfreeBufferArray(scip, &valbuf);
   SCIPfreeBufferArray(scip, &varbuf);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */
static
DECL_VARFIXINGS(varFixingsRens)
{
   int nbinvars;
   int nintvars;
   SCIP_VAR** vars;
   int i;
   assert(scip != NULL);
   assert(varbuf != NULL);
   assert(nfixings != NULL);
   assert(valbuf != NULL);

   if( ! SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

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
      if( SCIPisIntegral(scip, lpsolval) )
      {
         varbuf[*nfixings] = var;
         valbuf[*nfixings] = lpsolval;
         ++(*nfixings);
      }
   }

   *success = TRUE;

   return SCIP_OKAY;
}

static
DECL_CHANGESUBSCIP(changeSubscipRens)
{
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

      if( ! SCIPisIntegral(sourcescip, lpsolval) )
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

   return SCIP_OKAY;
}

static
DECL_VARFIXINGS(varFixingsRins)
{
   int nbinvars;
   int nintvars;
   SCIP_VAR** vars;
   int i;
   SCIP_SOL* incumbent;
   assert(scip != NULL);
   assert(varbuf != NULL);
   assert(nfixings != NULL);
   assert(valbuf != NULL);

   if( ! SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

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

   /* loop over binary and integer variables; determine those that should be fixed in the sub-SCIP */
   for( i = 0; i < nbinvars + nintvars; ++i )
   {
      SCIP_VAR* var = vars[i];
      SCIP_Real lpsolval = SCIPgetSolVal(scip, NULL, var);
      SCIP_Real incumbentsolval = SCIPgetSolVal(scip, incumbent, var);
      assert((i < nbinvars && SCIPvarIsBinary(var)) || (i >= nbinvars && SCIPvarIsIntegral(var)));

      /* fix all binary and integer variables with integer LP solution value */
      if( SCIPisEQ(scip, lpsolval, incumbentsolval) )
      {
         assert(SCIPisFeasIntegral(scip, incumbentsolval));
         incumbentsolval = SCIPfloor(scip, incumbentsolval + 0.5);
         varbuf[*nfixings] = var;
         valbuf[*nfixings] = incumbentsolval;
         ++(*nfixings);
      }
   }

   *success = TRUE;

   return SCIP_OKAY;
}

static
DECL_NHINIT(nhInitMutation)
{
   DATA_MUTATION* data;
   assert(scip != NULL);
   assert(neighborhood != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &neighborhood->data.mutation) );

   data = neighborhood->data.mutation;
   assert(data != NULL);

   SCIP_CALL( SCIPrandomCreate(&data->rng, SCIPblkmem(scip), SCIPinitializeRandomSeed(scip, MUTATIONSEED)) );

   return SCIP_OKAY;
}

static
DECL_NHEXIT(nhExitMutation)
{
   DATA_MUTATION* data;
   assert(scip != NULL);
   assert(neighborhood != NULL);
   data = neighborhood->data.mutation;
   assert(data != NULL);

   SCIPrandomFree(&data->rng);

   SCIPfreeBlockMemory(scip, &neighborhood->data.mutation);

   return SCIP_OKAY;
}

static
DECL_VARFIXINGS(varFixingsMutation)
{
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
   /* author bzfhende: TODO change target fixing rate to account for discrete variables specifically */
   if( nbinintvars <= ntargetfixings )
      return SCIP_OKAY;

   /* copy variables into a buffer array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &varscpy, vars, nbinvars + nintvars) );

   /* partially perturb the array until the number of target fixings are reached */
   for( i = 0; i < ntargetfixings; ++i )
   {
      int rand = SCIPrandomGetInt(rng, i, nbinintvars - 1);
      assert(rand < nbinintvars);

      if( rand > i )
      {
         SCIPswapPointers((void**)&varscpy[i], (void**)&varscpy[rand]);
      }
      /* copy the selected variables and their solution values into the buffer */
      varbuf[i] = varscpy[i];
      valbuf[i] = SCIPgetSolVal(scip, incumbentsol, varscpy[i]);
      assert(SCIPisIntegral(scip, valbuf[i]));
   }

   *nfixings = ntargetfixings;
   *success = TRUE;

   /* free the buffer array */
   SCIPfreeBufferArray(scip, &varscpy);

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

   heurdata->nneighborhoods = 0;

   /* include the RENS neighborhood */
   SCIP_CALL( lnsIncludeNeighborhood(scip, heurdata, &rens, "rens",
         DEFAULT_MINFIXINGRATE_RENS, DEFAULT_MAXFIXINGRATE_RENS, varFixingsRens, changeSubscipRens, NULL, NULL, NULL) );

   /* include the RINS neighborhood */
   SCIP_CALL( lnsIncludeNeighborhood(scip, heurdata, &rins, "rins",
         DEFAULT_MINFIXINGRATE_RINS, DEFAULT_MAXFIXINGRATE_RINS, varFixingsRins, NULL, NULL, NULL, NULL) );


   /* include the mutation neighborhood */
   SCIP_CALL( lnsIncludeNeighborhood(scip, heurdata, &mutation, "mutation",
         DEFAULT_MINFIXINGRATE_MUTATION, DEFAULT_MAXFIXINGRATE_MUTATION, varFixingsMutation, NULL, NULL, nhInitMutation, nhExitMutation) );


   /* TODO include the local branching neighborhood */

   /* author bzfhende: TODO include the crossover neighborhood */

   /* TODO include the DINS neighborhood */

   /* author bzfhende: TODO include the Zeroobjective neighborhood */

   /* author bzfhende: TODO include the Proximity neighborhood */

   /* author bzfhende: TODO include the GINS neighborhood/filter based on the LP relaxation data structure */

   /* author bzfhende: TODO include the Reduced-Cost neighborhood/filter */


   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitLns)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->epsgreedynh == NULL);


   /* init neighborhoods for new problem by resetting their statistics and fixing rate */
   for( i = 0; i < heurdata->nneighborhoods; ++i )
   {
      NH* neighborhood = heurdata->neighborhoods[i];

      SCIP_CALL( neighborhoodInit(scip, neighborhood) );

      SCIP_CALL( fixingRateReset(scip, &neighborhood->fixingrate) );

      neighborhoodStatsReset(scip, &neighborhood->stats);
   }

   /* create epsilon greedy bandit algorithm */
   SCIP_CALL( epsGreedyCreate(scip, &heurdata->epsgreedynh, DEFAULT_INITSEED, (void*)heurdata, epsNChoicesLns, epsRewardLns) );

   /* create an exp3 bandit algorithm */
   if( heurdata->exp3 == NULL )
   {
      SCIP_CALL( expThreeCreate(scip, &heurdata->exp3, DEFAULT_INITSEED, heurdata->nneighborhoods, heurdata->gamma, heurdata->beta) );
   }
   else
   {
      SCIP_CALL( expThreeReset(scip, heurdata->exp3, DEFAULT_INITSEED + SCIPgetNVars(scip)) );
   }

   heurdata->usednodes = 0;

   return SCIP_OKAY;
}

/** todo print neighborhood statistics */
static
void printNeighborhoodStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   int i;
   SCIPinfoMessage(scip, NULL, "Neighborhoods      :       Calls        Sols        Best        Exp3\n");

   /* todo loop over neighborhoods and fill in statistics */
   for( i = 0; i < heurdata->nneighborhoods; ++i )
   {
      NH* neighborhood;
      SCIP_Real proba;
      neighborhood = heurdata->neighborhoods[i];
      SCIPinfoMessage(scip, NULL, "  %-17s:", neighborhood->name);
      SCIPinfoMessage(scip, NULL, " %11d", neighborhood->stats.nruns);
      SCIPinfoMessage(scip, NULL, " %11d", neighborhood->stats.nsolsfound);
      SCIPinfoMessage(scip, NULL, " %11d", neighborhood->stats.nbestsolsfound);

      proba = (1 - heurdata->exp3->gamma) * heurdata->exp3->weights[i] / heurdata->exp3->weightsum + heurdata->exp3->gamma / heurdata->nneighborhoods;
      SCIPinfoMessage(scip, NULL, " %11.5f", proba);
      SCIPinfoMessage(scip, NULL, "\n");
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
   assert(heurdata->epsgreedynh != NULL);

   printNeighborhoodStatistics(scip, heurdata);

   epsGreedyFree(scip, &heurdata->epsgreedynh);

   /* free neighborhood specific data */
   for( i = 0; i < heurdata->nneighborhoods; ++i )
   {
      NH* neighborhood = heurdata->neighborhoods[i];

      SCIP_CALL( neighborhoodExit(scip, neighborhood) );
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
      expThreeFree(scip, &heurdata->exp3);

   /* free neighborhoods */
   for( i = 0; i < heurdata->nneighborhoods; ++i )
   {
      lnsFreeNeighborhood(scip, &(heurdata->neighborhoods[i]));
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
   heurdata->epsgreedyfilter = NULL;
   heurdata->exp3 = NULL;

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

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which LNS should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

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

   return SCIP_OKAY;
}
