/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_dps.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  dynamic partition search
 * @author Katrin Halbig
 *
 * The dynamic partition search (DPS) is a construction heuristic which additionally needs a
 * user decomposition.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/heur_dps.h"
#include "scip/pub_heur.h"
#include "scip/scip_heur.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"


#define HEUR_NAME             "dps"
#define HEUR_DESC             "primal heuristic for decomposable MIPs"
#define HEUR_DISPCHAR         '?' /** todo */
#define HEUR_PRIORITY         75000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyDps)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurDps(scip) );

   return SCIP_OKAY;
}


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeDps)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecDps)
{  /*lint --e{715}*/

   SCIPdebugMsg(scip, "initialize dps heuristic\n");

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the dps primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurDps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create dps primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   heur = NULL;

   /* include primal heuristic */

   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecDps, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyDps) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeDps) );

   /* add dps primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
