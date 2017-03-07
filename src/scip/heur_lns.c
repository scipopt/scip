/* * * * *e * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*         ward         This file is part of the program and library             */
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

#include "scip/heur_lns.h"


#define HEUR_NAME             "lns"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         '?'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** callback to let the neighborhood write its suggested variable fixings into a buffer data structure
 *
 * todo comments
 */
 #define DECL_VARFIXINGS(x) SCIP_RETCODE x ( \
   SCIP*                 scip,               \
   SCIP_VAR**            varbuf,             \
   SCIP_Real*            valbuf,             \
   int*                  nfixings,           \
   SCIP_Bool*            success             \
   )

/** callback for changes subproblem changes other than variable fixings
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

/* additional crossover data structure */
typedef struct data_crossover DATA_CROSSOVER;

typedef struct EpsGreedy EPSGREEDY;
typedef struct NH_FixingRate NH_FIXINGRATE;
typedef struct NH_Stats NH_STATS;
typedef struct Nh NH;

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
   int                   nruns;
   int                   nsolsfound;
   int                   nbestsolsfound;
   SCIP_Longint          usednodes;
   SCIP_Longint          lpiterations;
   int                   presolrounds;
   int                   totalnbinfixings;
   int                   totalnintfixings;
   int                   totalnimplintfixings;
   int                   totalncontfixings;
   SCIP_Real             totalgapclosed;
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
   const char*           name;
   NH_FIXINGRATE*        fixingrate;
   NH_STATS              stats;
   DECL_VARFIXINGS       ((*varfixings));
   DECL_CHANGESUBSCIP    ((*changesubscip));
   DECL_SETUPSUBSCIP     ((*setupsubscip));
   SCIP_Bool             active;
};

/** current, unnormalized reward for item i */
#define DECL_EPSREWARD(x) SCIP_Real x ( \
   SCIP*                 scip,             \
   int                   i                 \
   )

/** todo callback for the number of choices */
#define DECL_EPSNCHOICES(x) int x (        \
   SCIP*                 scip              \
   )

struct EpsGreedy
{
   SCIP_Real             eps;                /**< epsilon parameter (between 0 and 1) to control epsilon greedy */
   SCIP_RANDNUMGEN*      rng;                /**< random number generator for randomized selection of routines  */
   DECL_EPSREWARD        ((*epsreward));     /**< reward callback for unnormalized reward of the i'th item */
   DECL_EPSNCHOICES      ((*epsnchoices));   /**< callback for the number of choices */
   void*                 userptr;            /**< user pointer for callback functions */
};


/** todo create fixing rate */
static
SCIP_RETCODE fixingRateCreate(
   SCIP*                scip                /**< SCIP data structure */
   )
{
   return SCIP_OKAY;
}

/** todo free fixing rate */
static
SCIP_RETCODE fixingRateFree (
   SCIP*                scip                /**< SCIP data structure */
   )
{
   return SCIP_OKAY;
}

/** todo reset fixing rate */
static
SCIP_RETCODE fixingRateReset (
   SCIP*                scip                /**< SCIP data structure */
   )
{
   return SCIP_OKAY;
}

/** todo update fixing rate */
static
SCIP_RETCODE fixingRateUpdate (
   SCIP*                scip                /**< SCIP data structure */
   )
{
   return SCIP_OKAY;
}

/** todo create a neighborhood of the specified name */
static
SCIP_RETCODE lnsIncludeNeighborhood(
   SCIP*                scip                /**< SCIP data structure */
   )
{
   return SCIP_OKAY;
}

/** todo free a neighborhood and release its data */
static
SCIP_RETCODE lnsFreeNeighborhood(
   SCIP*                scip                /**< SCIP data structure */
   )
{
   return SCIP_OKAY;
}


/** todo reset an epsilon greedy selector */
static
SCIP_RETCODE epsGreedyReset (
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

   nchoices = epsgreedy->epsnchoices(scip);

   if( nchoices == 0 )
   {
      *i = -1;
      return SCIP_OKAY;
   }
   /** roll the dice to check if the best element should be picked, or an element at random */
   rand = SCIPrandomGetReal(epsgreedy->rng, 0.0, 1.0);

   if( rand <= epsgreedy->eps )
   {
      /** pick the element with the largest reward */
      int j;
      SCIP_Real maxreward = epsgreedy->epsreward(scip, 0);
      *i = 0;
      /* determine reward for every element */
      for( j = 1; j < nchoices; ++j )
      {
         SCIP_Real reward = epsgreedy->epsreward(scip, j);

         if( maxreward < reward )
         {
            *i = j;
            maxreward = reward;
         }
      }
   }
   else
   {
      int j;
      SCIP_Real* rewards;
      SCIP_Real lastreward = 0.0;
      SCIP_Real selectedreward;

      SCIP_CALL( SCIPallocBufferArray(scip, &rewards, nchoices) );

      /* collect rewards and sum them up */
      for( j = 0; j < nchoices; ++j )
      {
         rewards[j] = epsgreedy->epsreward(scip, j) + lastreward;
         lastreward = rewards[j];
      }

      /* draw a new random number */
      rand = SCIPrandomGetReal(epsgreedy->rng, 0.0, 1.0);
      selectedreward = lastreward * rand;

      assert(lastreward >= selectedreward);
      /* loop over collected rewards and stop at the selected reward */
      for( j = 0; j < nchoices; ++j )
      {
         if( selectedreward <= rewards[j] )
         {
            *i = j;
            break;
         }
      }

      SCIPfreeBufferArray(scip, &rewards);

   }

   return SCIP_OKAY;
}


/** todo define neighborhood and callbacks */

/** primal heuristic data */
struct SCIP_HeurData
{
   NH**                  neighborhoods;      /**< array of neighborhoods with the best one at the first position */
   int                   nneighborhoods;     /**< number of neighborhoods */
   EPSGREEDY*            epsgreedynh;        /**< epsilon greedy selector for a neighborhood */
   EPSGREEDY*            epsgreedyfilter;    /**< epsilon greedy selector for a filter strategy */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


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

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_HEURFREE(heurFreeLns)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lns primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurFreeLns NULL
#endif


/** initialization method of primal heuristic (called after problem was transformed) */
#if 0
static
SCIP_DECL_HEURINIT(heurInitLns)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lns primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitLns NULL
#endif


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitLns)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lns primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitLns NULL
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


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLns)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lns primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/


   /** todo check if budget allows a run of the next selected neighborhood */




   /** todo check if the selected neighborhood can run at this node */


   /** todo determine variable fixings and objective coefficients of this neighborhood */

   /** todo if too few fixings, use a strategy to select more variable fixings: randomized, LP graph, ReducedCost/Ps-Cost based, mix */

   /** todo later: run global propagation for this set of fixings */

   /** todo let the neighborhood add additional constraints, or restrict domains */

   /** todo alternatively: set up sub-SCIP and run presolving */



   /** todo was presolving successful enough regarding fixings? otherwise terminate */

   /** todo run sub-SCIP for the given budget, and collect statistics */

   /** todo determine the success of this neighborhood, and the target fixing rate for the next time */

   /** todo update the remaining budget */

   /** todo select the next neighborhood to run based on epsilon greedy bandit mechanism */

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

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

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyLns, heurFreeLns, heurInitLns, heurExitLns, heurInitsolLns, heurExitsolLns, heurExecLns,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecLns, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyLns) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeLns) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitLns) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitLns) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolLns) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolLns) );
#endif

   /* add lns primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
