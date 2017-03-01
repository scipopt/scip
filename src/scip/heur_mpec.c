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

/**@file   heur_mpec.c
 * @brief  mpec primal heuristic
 * @author Felipe Serrano
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_mpec.h"
#include "nlpi/nlpi.h"


#define HEUR_NAME             "mpec"
#define HEUR_DESC             "regularization heuristic for convex and nonconvex MINLPs"
#define HEUR_DISPCHAR         'W'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_NLPI* nlpi;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_HASHMAP* var2idx;
};


/*
 * Local methods
 */


static
SCIP_RETCODE createNLP(
   SCIP* scip,
   SCIP_HEURDATA* heurdata
   )
{
   assert(heurdata != NULL);

   if( heurdata->nlpi != NULL )
      return SCIP_OKAY;

   /* @todo */

   return SCIP_OKAY;
}

static
SCIP_RETCODE freeNLP(
   SCIP* scip,
   SCIP_HEURDATA* heurdata
   )
{
   assert(heurdata != NULL);

   if( heurdata->nlpi == NULL )
      return SCIP_OKAY;

   assert(heurdata->nlpiprob != NULL);
   assert(heurdata->var2idx != NULL);

   SCIPhashmapFree(&heurdata->var2idx);
   SCIP_CALL( SCIPnlpiFreeProblem(heurdata->nlpi, &heurdata->nlpiprob) );

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyMpec)
{  /*lint --e{715}*/
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurMpec(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeMpec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolMpec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIP_CALL( createNLP(scip, heurdata) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolMpec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIP_CALL( freeNLP(scip, heurdata) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecMpec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetNIntVars(scip) >= 0 || SCIPgetNBinVars(scip) == 0
      || SCIPgetNNlpis(scip) == 0 || !SCIPisNLPConstructed(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the mpec primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurMpec(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata = NULL;
   SCIP_HEUR* heur = NULL;

   /* create mpec primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   BMSclearMemory(heurdata);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecMpec, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyMpec) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeMpec) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolMpec) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolMpec) );

   /* add mpec primal heuristic parameters */

   return SCIP_OKAY;
}
