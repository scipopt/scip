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

/**@file   heur_redsize.c
 * @brief  primal heuristic that solves the problem with a sparser matrix as a submip
 * @author Leon Eifler
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_redsize.h"
#include "cycplugins.h"
#include "probdata_cyc.h"

#define HEUR_NAME             "redsize"
#define HEUR_DESC             "primal heuristic that solves the problem with a sparser matrix as a submip"
#define HEUR_DISPCHAR         'u'
#define HEUR_PRIORITY         536870911
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE           /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_REDUCTIONRATE 0.75           /**< default percentile of transition that gets deleted */

struct SCIP_HeurData
{
   SCIP_Real             reductionrate;      /**< percentile of transition that gets deleted */
};

/*
 * Local methods
 */

/** Add incomplete solution to main scip */
static
SCIP_RETCODE SCIPcycAddIncompleteSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure of subscip */
   SCIP_HEUR*            heur,               /**< pointer to heuristic */
   SCIP_SOL*             subsol,             /**< solution of subscip */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_VAR*** binvars;
   SCIP_Real** solclustering;
   SCIP_SOL* newsol;
   SCIP_Bool feasible;
   int nbins;
   int ncluster;
   int i;
   int t;

   nbins = SCIPcycGetNBins(scip);
   ncluster = SCIPcycGetNCluster(scip);
   binvars = SCIPcycGetBinvars(subscip);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solclustering, nbins) );

   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solclustering[i], ncluster) );

      for( t = 0; t < ncluster; ++t )
      {
         solclustering[i][t] = SCIPgetSolVal(subscip, subsol, binvars[i][t]);
      }
   }

   SCIP_CALL( assignVars(scip, newsol, solclustering, nbins, ncluster) );

   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, TRUE, &feasible) );

   if( feasible )
      *result = SCIP_FOUNDSOL;

   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeBlockMemoryArray(scip, &solclustering[i], ncluster);
   }

   SCIPfreeBlockMemoryArray(scip, &solclustering, nbins);

   return SCIP_OKAY;
}

/** set all the given percentile of nonzeros to zero */
static
SCIP_RETCODE SCIPreduceMatrixSize(
   SCIP*                scip,                /**< SCIP data structure */
   SCIP_Real**          matrix,              /**< the matrix */
   SCIP_Real            percentile,          /**< the percentile of entries to be deleted */
   SCIP_Real            scale,               /**< scaling between net flow and coherence */
   int                  size                 /**< the size of the matrix */
   )
{
   SCIP_Real* nonzeros;
   int*  idxnonzeros;
   int nnonzeros;
   int currentsize;
   int i;
   int j;
   int k;

   nnonzeros = 0;
   currentsize = size;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nonzeros, size) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idxnonzeros, size) );

   for( i = 0; i < size; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( !SCIPisZero(scip, matrix[i][j]) || !SCIPisZero(scip, matrix[j][i]) )
         {
            /* if we have non-zero entry, compute the impact and save it */
            nonzeros[nnonzeros] = MAX(scale * (matrix[i][j]+matrix[j][i]), matrix[i][j] - matrix[j][i]);
            idxnonzeros[nnonzeros] = i * size + j;

            nnonzeros++;

            /* realloc if necessary */
            if( currentsize < nnonzeros + 2 )
            {
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nonzeros, currentsize, currentsize + size) );
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &idxnonzeros, currentsize, currentsize + size) );
               currentsize += size;
            }
         }
      }
   }

   /* sort by the least impact */
   SCIPsortRealInt(nonzeros, idxnonzeros, nnonzeros);

   /* recompute the indizes and set them to 0 */
   for( i = 0; i < nnonzeros * percentile; ++i )
   {
      j = idxnonzeros[i] % size;
      k = idxnonzeros[i] / size;

      matrix[j][k] = 0.0;
      matrix[k][j] = 0.0;
   }

   SCIPfreeBlockMemoryArray(scip, &nonzeros, currentsize);
   SCIPfreeBlockMemoryArray(scip, &idxnonzeros, currentsize);

   return SCIP_OKAY;
}

/** main procedure of the heuristic, creates and solves a sub-SCIP */
static
SCIP_RETCODE SCIPapplyRedSize(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_RESULT*          result,             /**< result data structure */
   SCIP_Real             reductionrate,      /**< minimum percentage of integer variables that have to be fixed */
   SCIP_Longint          maxnodes            /**< maximum number of  nodes for the subproblem */
   )
{
   SCIP* subscip;
   SCIP_Real** cmatrix_orig;
   SCIP_Real** cmatrix_red;
   SCIP_SOL** subsols;

   SCIP_Real scale;

   int nbins;
   int ncluster;
   int nsubsols;
   int i;
   int j;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   assert(maxnodes >= 0);

   assert(0.0 <= reductionrate && reductionrate <= 1.0);

   nbins = SCIPcycGetNBins(scip);
   ncluster = SCIPcycGetNCluster(scip);
   cmatrix_orig = SCIPcycGetCmatrix(scip);
   scale = SCIPcycGetScale(scip);

   assert(nbins > 0 && ncluster > 0 && nbins >= ncluster);
   assert(cmatrix_orig != NULL);

   *result = SCIP_DIDNOTRUN;

   /* create subscip */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPincludeCycPlugins(subscip) );

#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", -1) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

   /* copy the original cmatrix */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cmatrix_red, nbins) );
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cmatrix_red[i], nbins) );
      for( j = 0; j < nbins; ++j )
      {
         cmatrix_red[i][j] = cmatrix_orig[i][j];
      }
   }

   /* delete entries from the copied cmatrix */
   SCIP_CALL( SCIPreduceMatrixSize(subscip, cmatrix_red, reductionrate, scale, nbins) );

   /* create probdata for the subscip */
   SCIP_CALL( SCIPcreateProbCyc(subscip, "subscip", nbins, ncluster, cmatrix_red) );

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* set nodelimit */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", maxnodes) );

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

   /* enable conflict analysis, disable analysis of boundexceeding LPs, and restrict conflict pool */
   if( !SCIPisParamFixed(subscip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", TRUE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useboundlp") )
   {
      SCIP_CALL( SCIPsetCharParam(subscip, "conflict/useboundlp", 'o') );
   }
   if( !SCIPisParamFixed(subscip, "conflict/maxstoresize") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "conflict/maxstoresize", 100) );
   }

   /* solve the subproblem */
   SCIP_CALL_ABORT( SCIPsolve(subscip) );

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   /* check, whether a solution was found; due to numerics, it might happen that not all solutions are feasible -> try
    * all solutions until one was accepted
    */
   nsubsols = SCIPgetNSols(subscip);
   subsols = SCIPgetSols(subscip);

   for( i = 0; i < nsubsols; ++i )
      SCIP_CALL( SCIPcycAddIncompleteSol(scip, subscip, heur, subsols[i], result) );

   SCIP_CALL( SCIPfree(&subscip) );

   /* free memory */
   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeBlockMemoryArray(scip, &cmatrix_red[i], nbins);
   }
   SCIPfreeBlockMemoryArray(scip, &cmatrix_red, nbins);

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRedsize)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurRedsize(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRedsize)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRedsize)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( SCIPapplyRedSize(scip, heur, result, heurdata->reductionrate, (SCIP_Longint) 1) );

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the oneopt primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRedsize(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create redsize primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
      HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
      HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRedsize, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRedsize) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRedsize) );

   /* add param */
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/reduction-rate",
         "percentile of transition probabilities that gets deleted in the submip",
         &heurdata->reductionrate, FALSE, DEFAULT_REDUCTIONRATE, 0.0, 1.0, NULL , NULL) );

   return SCIP_OKAY;
}
