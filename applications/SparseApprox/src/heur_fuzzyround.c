
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_fuzzyround.c
 * @brief  primal heuristic that constructs a feasible solution from the lp-relaxation. Round only on the bin-variables
 * and then reconstruct the rest of the variables accordingly.
 * @author Leon Eifler
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "probdata_spa.h"
#include "heur_fuzzyround.h"
#include "scip/cons_and.h"
/* @note If the heuristic runs in the root node, the timing is changed to (SCIP_HEURTIMING_DURINGLPLOOP |
 *       SCIP_HEURTIMING_BEFORENODE), see SCIP_DECL_HEURINITSOL callback.
 */

#define HEUR_NAME             "fuzzyround"
#define HEUR_DESC             "primal heuristic that constructs a feasible solution from the lp-relaxation"
#define HEUR_DISPCHAR         '&'
#define HEUR_PRIORITY         1000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE          /**< does the heuristic use a secondary SCIP instance? */


/*
 * Local methods
 */


/**
 * assign the variables in scip according to the found clusterassignment
 */

/**
 * assign the variables in scip according to the found clusterassignment
 */

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecFuzzyround)
{  /*lint --e{715}*/
   int i;
   int k;
   SCIP_Real maxlpval;
   int maxcluster;
   SCIP_VAR*** binvars;
   SCIP_SOL* sol;
   int nbins;
   int ncluster;
   SCIP_Real** clustering;
   int* binsincluster;
   SCIP_Bool feasible = FALSE;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);
   assert(nbins > 0);
   assert(ncluster > 0 && ncluster <= nbins);

   binvars = SCIPspaGetBinvars(scip);
   assert(binvars != NULL);

   /* allocate memory */

   SCIP_CALL( SCIPallocClearMemoryArray(scip, &clustering , nbins) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &binsincluster, ncluster) );
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &clustering[i], ncluster) );
   }
   /* for each bin, set the assignment with the highest lp-value to 1, the rest to 0 */
   for( i = 0; i < nbins; ++i ) {
      maxlpval = 0;
      maxcluster = -1;
      for (k = 0; k < ncluster; ++k)
      {
         assert( NULL != binvars[i][k]);
         if( SCIPisGT(scip, SCIPvarGetLPSol(binvars[i][k]), maxlpval) )
         {
            maxlpval = SCIPvarGetLPSol(binvars[i][k]);
            maxcluster = k;
            binsincluster[k]++;
         }
         else if( SCIPisEQ(scip, SCIPvarGetLPSol(binvars[i][k]), maxlpval) &&  maxcluster != -1 && binsincluster[maxcluster] > binsincluster[k] )
         {
            binsincluster[maxcluster]--;
            binsincluster[k]++;
            maxcluster = k;
         }
      }
      assert(maxcluster >= 0);
      clustering[i][maxcluster] = 1.0;
   }
   assert(isPartition(scip, clustering, nbins, ncluster));

   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
   assignVars(scip, sol, clustering, nbins, ncluster);
   SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, TRUE, TRUE, TRUE, TRUE, &feasible) );
   if( feasible )
      *result = SCIP_FOUNDSOL;
   else
      *result = SCIP_DIDNOTFIND;

   /** free allocated memory */

   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeMemoryArray(scip, &clustering[i]);
   }
   SCIPfreeMemoryArray(scip, &clustering);
   SCIPfreeMemoryArray(scip, &binsincluster);
   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the oneopt primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurFuzzyround(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_HEUR* heur;


   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
      HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
      HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecFuzzyround, NULL) );

   assert(heur != NULL);


   return SCIP_OKAY;
}
