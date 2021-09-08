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
#include "scip/scipdefplugins.h"
#include "scip/scip_cons.h"
#include "scip/scip_dcmp.h"
#include "scip/scip_general.h"
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
   SCIP_CONS**           linkingconss;       /**< linking constraints */
   int                   nlinking;           /**< number of linking constraints */
   int                   nblocks;            /**< number of blocks */
};

struct Blockproblem
{
   SCIP*                 blockscip;          /**< SCIP data structure */
   SCIP_VAR**            blockvars;          /**< block variables
                                              *   blockvars[nblockvars] is first slack variable */
   SCIP_VAR**            slackvars;          /**< slack variables */
   SCIP_CONS**           linkingconss;       /**< linking constraints */
   int*                  linkingindices;     /**< indices of linking constraints in original problem */
   int                   nlinking;           /**< number of linking constraints */
   int                   nblockvars;         /**< number of block variables */
   int                   nslackvars;         /**< number of slack variables */
   SCIP_Real*            origobj;            /**< original objective coefficients */
};
typedef struct Blockproblem BLOCKPROBLEM;

/*
 * Local methods
 */

/** assigns linking variables to last block
 *
 * The labels are copied to newdecomp and the linking variables are assigned to the last block (i.e. highest block label).
 * Constraint labels and statistics are recomputed.
 */
static
SCIP_RETCODE assignLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          newdecomp,          /**< decomposition with assigned linking variables */
   SCIP_VAR**            sortedvars,         /**< sorted array of variables */
   SCIP_CONS**           sortedconss,        /**< sorted array of constraints */
   int*                  sortedvarlabels,    /**< sorted array of variable labels */
   int*                  sortedconslabels,   /**< sorted array of constraint labels */
   int                   nvars,              /**< number of variables */
   int                   nconss,             /**< number of constraints */
   int                   nlinkvars           /**< number of linking variables */
   )
{
   int newlabel;
   int maxgraphedge;
   int v;

   assert(scip != NULL);
   assert(newdecomp != NULL);
   assert(sortedvars != NULL);
   assert(sortedconss != NULL);
   assert(sortedvarlabels != NULL);
   assert(sortedconslabels != NULL);

   /* we do not need the block decomposition graph of the statistics */
   SCIP_CALL( SCIPgetIntParam(scip, "decomposition/maxgraphedge", &maxgraphedge) );
   if( !SCIPisParamFixed(scip, "decomposition/maxgraphedge") )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "decomposition/maxgraphedge", 0) );
   }

   /* copy the labels */
   SCIP_CALL( SCIPdecompSetVarsLabels(newdecomp, sortedvars, sortedvarlabels, nvars) );
   SCIP_CALL( SCIPdecompSetConsLabels(newdecomp, sortedconss, sortedconslabels, nconss) );

   /* assign linking variables */
   newlabel = sortedvarlabels[nvars - 1]; /* take always label of last block */
   assert(newlabel >= 0);
   for( v = 0; v < nlinkvars; v++ )
   {
      SCIP_CALL( SCIPdecompSetVarsLabels(newdecomp, &sortedvars[v], &newlabel, 1) );
   }
   SCIPdebugMsg(scip, "assigned %d linking variables\n", nlinkvars);

   /* recompute constraint labels and statistics */
   SCIP_CALL( SCIPcomputeDecompConsLabels(scip, newdecomp, sortedconss, nconss) );
   SCIP_CALL( SCIPcomputeDecompStats(scip, newdecomp, TRUE) );
   nlinkvars = SCIPdecompGetNBorderVars(newdecomp);

   /* get new labels and sort */
   SCIPdecompGetConsLabels(newdecomp, sortedconss, sortedconslabels, nconss);
   SCIPdecompGetVarsLabels(newdecomp, sortedvars, sortedvarlabels, nvars);
   SCIPsortIntPtr(sortedconslabels, (void**)sortedconss, nconss);
   SCIPsortIntPtr(sortedvarlabels, (void**)sortedvars, nvars);

   /* After assigning the linking variables, blocks can have zero constraints.
    * So the remaining variables are labeled as linking in SCIPcomputeDecompStats().
    * We assign this variables to the same label as above.
    */
   if( nlinkvars >= 1 )
   {
      assert(sortedvarlabels[0] == SCIP_DECOMP_LINKVAR);
      SCIPdebugMsg(scip, "assign again %d linking variables\n", nlinkvars);

      for( v = 0; v < nlinkvars; v++ )
      {
         SCIP_CALL( SCIPdecompSetVarsLabels(newdecomp, &sortedvars[v], &newlabel, 1) );
      }
      SCIP_CALL( SCIPcomputeDecompConsLabels(scip, newdecomp, sortedconss, nconss) );
      SCIP_CALL( SCIPcomputeDecompStats(scip, newdecomp, TRUE) );

      SCIPdecompGetConsLabels(newdecomp, sortedconss, sortedconslabels, nconss);
      SCIPdecompGetVarsLabels(newdecomp, sortedvars, sortedvarlabels, nvars);
      SCIPsortIntPtr(sortedconslabels, (void**)sortedconss, nconss);
      SCIPsortIntPtr(sortedvarlabels, (void**)sortedvars, nvars);
   }
   assert(sortedvarlabels[0] != SCIP_DECOMP_LINKVAR);

   SCIP_CALL( SCIPsetIntParam(scip, "decomposition/maxgraphedge", maxgraphedge) );

   return SCIP_OKAY;
}

/** creates a sub-SCIP and sets parameters */
static
SCIP_RETCODE createSubscip(
    SCIP*                scip,               /**< main SCIP data structure */
    SCIP**               subscip             /**< pointer to store created sub-SCIP */
   )
{
   assert(scip != NULL);
   assert(subscip != NULL);

   /* create a new SCIP instance */
   SCIP_CALL( SCIPcreate(subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(*subscip) );

   SCIP_CALL( SCIPcopyLimits(scip, *subscip) );

   /* avoid recursive calls */
   SCIP_CALL( SCIPsetSubscipsOff(*subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(*subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(*subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* disable expensive techniques */
   SCIP_CALL( SCIPsetIntParam(*subscip, "misc/usesymmetry", 0) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(*subscip, "misc/catchctrlc", FALSE) );

#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/freq", 100000000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(*subscip, "timing/statistictiming", FALSE) );
#endif

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(*subscip, "lp/checkdualfeas", FALSE) );

   return SCIP_OKAY;
}


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
   SCIP_DECOMP* assigneddecomp;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_VAR** sortedvars;
   SCIP_CONS** sortedconss;
   SCIP_HEURDATA* heurdata;
   BLOCKPROBLEM** blockproblem;
   int* sortedvarlabels;
   int* sortedconslabels;
   SCIP_Real memory; /* in MB */
   int ndecomps;
   int nvars;
   int nconss;
   int nblocks;
   int b;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   assigneddecomp = NULL;
   blockproblem = NULL;

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

   nblocks = SCIPdecompGetNBlocks(decomp);
   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);

   /* if problem has no constraints, no variables or less than two blocks, return */
   if( nconss == 0 || nvars == 0 || nblocks <= 1 )
   {
      SCIPdebugMsg(scip, "problem has no constraints, no variables or less than two blocks\n");
      return SCIP_OKAY;
   }

   /* estimate required memory for all blocks and terminate if not enough memory is available */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memory) );
   if( ((SCIPgetMemUsed(scip) + SCIPgetMemExternEstim(scip))/1048576.0) * (nblocks/4.0 + 2) >= memory )
   {
      SCIPdebugMsg(scip, "The estimated memory usage for %d blocks is too large.\n", nblocks);
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

   if( sortedvarlabels[0] == SCIP_DECOMP_LINKVAR )
   {
      /* create new decomposition; don't change the decompositions in the decompstore */
      SCIP_CALL( SCIPdecompCreate(&assigneddecomp, SCIPblkmem(scip), nblocks, SCIPdecompIsOriginal(decomp), SCIPdecompUseBendersLabels(decomp)) );

      SCIP_CALL( assignLinking(scip, assigneddecomp, sortedvars, sortedconss, sortedvarlabels, sortedconslabels, nvars, nconss, SCIPdecompGetNBorderVars(decomp)) );
      assert(SCIPdecompGetNBlocks(decomp) >= SCIPdecompGetNBlocks(assigneddecomp));
      decomp = assigneddecomp;

      /* number of blocks can get smaller */
      nblocks = SCIPdecompGetNBlocks(decomp);
   }

#ifdef SCIP_DEBUG
      char buffer[SCIP_MAXSTRLEN];
      SCIPdebugMsg(scip, "DPS used decomposition:\n%s\n", SCIPdecompPrintStats(decomp, buffer));
#endif

   if( sortedvarlabels[0] == SCIP_DECOMP_LINKVAR ||
         sortedconslabels[0] != SCIP_DECOMP_LINKCONS ||
         nblocks <= 1 )
   {
      SCIPdebugMsg(scip, "Problem has linking variables or no linking constraints or less than two blocks\n");
      goto TERMINATE;
   }

   /* initialize heurdata */
   heurdata->linkingconss = sortedconss;
   heurdata->nlinking = SCIPdecompGetNBorderConss(decomp);
   heurdata->nblocks = nblocks;

   /* allocate memory for blockproblems and initialize partially */
   SCIP_CALL( SCIPallocBufferArray(scip, &blockproblem, nblocks) );
   for( b = 0; b < nblocks; b++ )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &blockproblem[b]) );
      createSubscip(scip, &blockproblem[b]->blockscip);

      SCIP_CALL( SCIPallocBufferArray(scip, &blockproblem[b]->linkingconss, heurdata->nlinking) );
      SCIP_CALL( SCIPallocBufferArray(scip, &blockproblem[b]->linkingindices, heurdata->nlinking) );
      SCIP_CALL( SCIPallocBufferArray(scip, &blockproblem[b]->slackvars, heurdata->nlinking * 2) ); /* maximum two slacks per linking constraint */
      SCIP_CALL( SCIPallocBufferArray(scip, &blockproblem[b]->origobj, nvars) );
      blockproblem[b]->blockvars = NULL;
      blockproblem[b]->nblockvars = 0;
      blockproblem[b]->nlinking = 0;
      blockproblem[b]->nslackvars = 0;
   }

   /** ------------------------------------------------------------------------ */
   /** free memory */
TERMINATE:
   if( blockproblem != NULL )
   {
      for( b = nblocks - 1; b >= 0; b-- )
      {
         SCIPfreeBufferArray(scip, &(blockproblem[b])->origobj);
         SCIPfreeBufferArray(scip, &(blockproblem[b])->slackvars);
         SCIPfreeBufferArray(scip, &(blockproblem[b])->linkingindices);
         SCIPfreeBufferArray(scip, &(blockproblem[b])->linkingconss);
         SCIP_CALL( SCIPfree(&blockproblem[b]->blockscip) );
         SCIPfreeBlockMemory(scip, &blockproblem[b]);
      }
      SCIPfreeBufferArray(scip, &blockproblem);
   }

   if( assigneddecomp != NULL )
   {
      SCIPdecompFree(&assigneddecomp, SCIPblkmem(scip));
   }

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
