/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_staircase.c
 * @brief  start heuristic for decompositions with staircase structure
 * @author Christine Tawfik
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/heur_staircase.h"
#include "scip/decomp.h"

#define HEUR_NAME             "staircase"
#define HEUR_DESC             "start heuristic for decompositions with staircase structure"
#define HEUR_DISPCHAR         'D'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

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

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_HEURCOPY(heurCopyStaircase)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurCopyStaircase NULL
#endif

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_HEURFREE(heurFreeStaircase)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurFreeStaircase NULL
#endif


/** initialization method of primal heuristic (called after problem was transformed) */
#if 0
static
SCIP_DECL_HEURINIT(heurInitStaircase)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitStaircase NULL
#endif


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitStaircase)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitStaircase NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolStaircase)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolStaircase NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolStaircase)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolStaircase NULL
#endif

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecStaircase)
{  /*lint --e{715}*/

   SCIP_DECOMPSTORE* decompstore;
   SCIP_DECOMP** decomps;
   int ndecomps;
   int d;

   /* query if there is a user decomposition (see methods in decomp.h) */
   decompstore = SCIPgetDecompstore(scip);
   assert(decompstore != NULL);

   decomps = SCIPdecompstoreGetDecomps(decompstore);
   ndecomps = SCIPdecompstoreGetNDecomps(decompstore);

   /* staircase heuristic requires a decomposition to work with */
   if( ndecomps == 0 )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Staircase heuristic works on %d decomposition%s\n", ndecomps, ndecomps == 1 ? "" : "s");



   /* author bzfhende
    *
    * TODO check if the decomposition has the appropriate staircase structure
    *    - no linking constraints
    *    - all linking variables link exactly 2 named blocks of constraints
    *    - the block connectivity graph is a path
    *
    *    if all the above assumptions are true, we should have an order of the blocks for the rolling horizon (if graph
    *    is a path, there are exactly two orders depending on which degree 1 node we start from)
    */

   /* todo comment */
   for( d = 0; d < ndecomps; ++d )
   {
      SCIP_DECOMP* decomp;
      SCIP_CONS** conss;
      SCIP_VAR** vars;
      int* conslabels;
      int* varlabels;
      int nlinkingconss = 0;
      int nlinkingvars = 0;
      int ncontlinkingvars = 0;
      int nconss;
      int nvars;
      int i;
      int v;

      nconss = SCIPgetNConss(scip);
      conss = SCIPgetConss(scip);
      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);

      decomp = decomps[d];

      /* does the decomposition have linking constraints? */
      SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varlabels, nvars) );

      SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);

      SCIPdecompGetVarsLabels(decomp, vars, varlabels, nvars);

      /* todo count number of linking constraints */
      for( i = 0; i < nconss; ++i )
      {
         if( conslabels[i] == SCIP_DECOMP_LINKCONS )
         {
            ++nlinkingconss;
            SCIPdebugMsg(scip, "Constraint %s is a linking constraint\n", SCIPconsGetName(conss[i]));

            SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );
         }
      }

      SCIPdebugMsg(scip, "Decomposition %d has %d linking constraints\n", d, nlinkingconss);

      /* todo count number of linking variables */
      for( v = 0; v < nvars; ++v )
      {
         if( varlabels[v] == SCIP_DECOMP_LINKVAR )
         {
            ++nlinkingvars;
            if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
               ++ncontlinkingvars;
         }
      }
      SCIPdebugMsg(scip, "Decomposition %d has %d linking variables (%d continuous)\n", d, nlinkingvars, ncontlinkingvars);

      SCIPfreeBufferArray(scip, &varlabels);
      SCIPfreeBufferArray(scip, &conslabels);

   }

   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
   /* author bzfhende
    *
    * TODO create or query a start solution
    *
    * - incumbent solution, if available
    * - zero solution (always available, perhaps infeasible)
    * - rounded LP relaxation solution (rounding based on locks or fractionality)
    * - partial solution (possible user input, see reader_mst.h)
    * - ...
    */

   /* author bzfhende
    *
    * TODO heuristic parameters
    *
    * - width of the horizon (previous blocks, subsequent blocks (with integer restrictions intact), subsequent relaxed blocks)
    *
    */

   /* author bzfhende
    *
    * TODO create a working solution by copying the start solution. This working solution is iteratively updated in the
    * main loop
    */


   /* author bzfhende
    *
    * TODO main loop: apply the rolling horizon algorithm.
    *
    * - in each iteration, formulate the corresponding sub-MIP in an auxiliary SCIP instance (have a look at the methods in scip_copy.c)
    * - solve this sub-MIP (have a look at heur_rins.c)
    * - update the working solution with the result from the sub-MIP
    */

   /* author bzfhende
    *
    * TODO if working solution is feasible, add it to SCIP
    */

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the staircase primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurStaircase(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create staircase primal heuristic data */
   heurdata = NULL;

   heur = NULL;

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyStaircase, heurFreeStaircase, heurInitStaircase, heurExitStaircase, heurInitsolStaircase, heurExitsolStaircase, heurExecStaircase,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecStaircase, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyStaircase) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeStaircase) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitStaircase) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitStaircase) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolStaircase) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolStaircase) );
#endif

   /* add staircase primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
