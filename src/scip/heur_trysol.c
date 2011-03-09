/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_trysol.c
 * @ingroup PRIMALHEURISTICS
 * @brief  primal heuristic that tries a given solution
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_trysol.h"


#define HEUR_NAME             "trysol"
#define HEUR_DESC             "try solution heuristic"
#define HEUR_DISPCHAR         'y'
#define HEUR_PRIORITY         -3000000     /* should process after all other heuristics */
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_DURINGLPLOOP




/*
 * Data structures
 */


/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*      sol;                /**< storing solution passed to heuristic (NULL if none) */
};




/*
 * Callback methods of primal heuristic
 */


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeTrysol)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );

   SCIPdebugMessage("free method of trysol primal heuristic.\n");

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeMemory(scip, &heurdata);

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolTrysol)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );

   SCIPdebugMessage("exitsol method of trysol primal heuristic.\n");

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free solution if one is still present */
   if ( heurdata->sol != NULL )
      SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTrysol)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Bool stored;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only run if solution present */
   if ( heurdata->sol == NULL )
      return SCIP_OKAY;

   SCIPdebugMessage("exec method of trysol primal heuristic.\n");
   *result = SCIP_DIDNOTFIND;

   /* try solution and free it - check everything, because we are not sure */
   SCIP_CALL( SCIPtrySolFree(scip, &heurdata->sol, TRUE, TRUE, TRUE, &stored) );
   assert( heurdata->sol == NULL );

   if ( stored )
      *result = SCIP_FOUNDSOL;

   return SCIP_OKAY;
}


#define heurInitTrysol NULL
#define heurExitTrysol NULL
#define heurInitsolTrysol NULL


/*
 * primal heuristic specific interface methods
 */

/** creates the trysol primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurTrySol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heurdata->sol = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, heurFreeTrysol, heurInitTrysol, heurExitTrysol,
         heurInitsolTrysol, heurExitsolTrysol, heurExecTrysol, heurdata) );

   return SCIP_OKAY;
}


/** pass solution to trysol heuristic */
SCIP_RETCODE SCIPheurPassSolTrySol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< trysol heuristic */
   SCIP_SOL*             sol                 /**< solution to be passed */
   )
{
   SCIP_HEURDATA* heurdata;

   assert( scip != NULL );
   assert( heur != NULL );
   assert( sol != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only store solution if we are not within our own SCIPtrySol() call */
   if ( heurdata->sol == NULL )
   {
      SCIP_CALL( SCIPcreateSolCopy(scip, &heurdata->sol, sol) );
      SCIPsolSetHeur(heurdata->sol, heur);
   }

   return SCIP_OKAY;
}
