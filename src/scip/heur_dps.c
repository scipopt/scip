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
#include "scip/pub_dcmp.h"
#include "scip/pub_heur.h"
#include "scip/pub_misc.h"
#include "scip/scip_cons.h"
#include "scip/scip_dcmp.h"
#include "scip/scip_heur.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"


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
   int                   nblocks;            /**< number of blocks */
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
   SCIP_DECOMP** alldecomps;
   SCIP_DECOMP* decomp;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_VAR** sortedvars;
   SCIP_CONS** sortedconss;
   SCIP_HEURDATA* heurdata;
   int* sortedvarlabels;
   int* sortedconslabels;
   SCIP_Real memory; /* in MB */
   int ndecomps;
   int nvars;
   int nconss;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   /* -------------------------------------------------------------------- */
   SCIPdebugMsg(scip, "initialize dps heuristic\n");

   /* take the first transformed decomposition */
   SCIPgetDecomps(scip, &alldecomps, &ndecomps, FALSE);
   if( ndecomps == 0)
      return SCIP_OKAY;

   decomp = alldecomps[0];
   assert(decomp != NULL);
   SCIPdebugMsg(scip, "First transformed decomposition is selected\n");

   heurdata->nblocks = SCIPdecompGetNBlocks(decomp);
   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);

   /* if problem has no constraints, no variables or less than two blocks, return */
   if( nconss == 0 || nvars == 0 || heurdata->nblocks <= 1 )
   {
      SCIPdebugMsg(scip, "problem has no constraints, no variables or less than two blocks\n");
      return SCIP_OKAY;
   }

   /* estimate required memory for all blocks and terminate if not enough memory is available */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memory) );
   if( ((SCIPgetMemUsed(scip) + SCIPgetMemExternEstim(scip))/1048576.0) * (heurdata->nblocks/4.0 + 2) >= memory )
   {
      SCIPdebugMsg(scip, "The estimated memory usage for %d blocks is too large.\n", heurdata->nblocks);
      return SCIP_OKAY;
   }

   vars = SCIPgetVars(scip);
   conss = SCIPgetConss(scip);
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sortedvars, vars, nvars) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sortedconss, conss, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortedvarlabels, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortedconslabels, nconss) );

   /* get labels and sort in increasing order */
   SCIPdecompGetVarsLabels(decomp, sortedvars, sortedvarlabels, nvars);
   SCIPdecompGetConsLabels(decomp, sortedconss, sortedconslabels, nconss);
   SCIPsortIntPtr(sortedvarlabels, (void**)sortedvars, nvars);
   SCIPsortIntPtr(sortedconslabels, (void**)sortedconss, nconss);

   if( sortedvarlabels[0] == SCIP_DECOMP_LINKVAR ||
         sortedconslabels[0] != SCIP_DECOMP_LINKCONS ||
         heurdata->nblocks <= 1 )
   {
      SCIPdebugMsg(scip, "Problem has linking variables or no linking constraints or less than two blocks\n");
      goto TERMINATE;
   }

   /** ------------------------------------------------------------------------ */
   /** free memory */
TERMINATE:
   if( sortedconslabels != NULL )
   {
      SCIPfreeBufferArray(scip, &sortedconslabels);
   }

   if( sortedvarlabels != NULL )
   {
      SCIPfreeBufferArray(scip, &sortedvarlabels);
   }

   if( sortedconss != NULL )
   {
      SCIPfreeBufferArray(scip, &sortedconss);
   }

   if( sortedvars != NULL )
   {
      SCIPfreeBufferArray(scip, &sortedvars);
   }

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
