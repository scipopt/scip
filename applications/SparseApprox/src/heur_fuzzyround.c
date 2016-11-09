
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
 * @brief  primal heuristic that constructs a feasible solution from the lp-relaxation
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
#define HEUR_FREQ             -1
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
static
SCIP_RETCODE assignVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< The SCIP solution */
   SCIP_Real**           clustering,         /**< The matrix with the clusterassignment */
   int                   nbins,              /**< The number of bins */
   int                   ncluster,           /**< The number of cluster */
   SCIP_Real**           cmatrix
)
{
   int i,j;
   int c;
   int c2;
   SCIP_VAR* var;
   SCIP_VAR*** binvars;
   SCIP_VAR****  edgevars;

   assert(nbins > 0 && ncluster > 0);

   binvars = SCIPspaGetBinvars(scip);
   edgevars = SCIPspaGetEdgevars(scip);

   for ( c = 0; c < ncluster; ++c )
   {
      /* set values of binary variables */
      for ( i = 0; i < nbins; ++i )
      {
         if( NULL != binvars[i][c] )
         {
            if( SCIPvarIsTransformed(binvars[i][c]) )
               var = binvars[i][c];
            else
               var = SCIPvarGetTransVar(binvars[i][c] );
            /* check if the clusterassignment ist feasible for the variable bounds. If not do not assign the variable */
            if( var != NULL && SCIPisLE(scip, SCIPvarGetLbGlobal(var), clustering[i][c]) && SCIPisGE(scip, SCIPvarGetUbGlobal(var), clustering[i][c]) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
               SCIP_CALL( SCIPsetSolVal( scip, sol, var, clustering[i][c]) );
            assert( SCIPisIntegral(scip, clustering[i][c]) );
         }
      }

      /* set the value for the edgevariables */
      for( i = 0; i < nbins; ++i )
      {
         for( j = 0; j < i; ++j )
         {
            if( NULL == edgevars[i][j][0] )
               continue;
            if( SCIPvarIsTransformed(edgevars[i][j][0]) )
               var = edgevars[i][j][0];
            else
               var = SCIPvarGetTransVar(edgevars[i][j][0]);
            if( var != NULL && SCIPisGE(scip, SCIPvarGetUbGlobal(var), clustering[j][c] * clustering[i][c]) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
            {
               if( SCIPisEQ(scip, 1.0, clustering[j][c] * clustering[i][c]) )
                  SCIP_CALL( SCIPsetSolVal( scip, sol, var, clustering[j][c] * clustering[i][c]  ) );
            }
            for( c2 = 0; c2 < ncluster; ++c2 )
            {
               if( NULL == edgevars[i][j][c] || c == c2 )
                  continue;
               if( c2 == c + 1 || ( c2 == 0 && c == ncluster -1) )
               {
                  if( SCIPvarIsTransformed(edgevars[i][j][1]) )
                     var = edgevars[i][j][1];
                  else
                     var = SCIPvarGetTransVar(edgevars[i][j][1]);
               }
               else if( c2 == c - 1 || ( c == 0 && c2 == ncluster -1) )
               {
                  if( SCIPvarIsTransformed(edgevars[j][i][1]) )
                     var = edgevars[j][i][1];
                  else
                     var = SCIPvarGetTransVar(edgevars[j][i][1]);
               }
               else
               {
                  if( SCIPvarIsTransformed(edgevars[i][j][2]) )
                     var = edgevars[i][j][2];
                  else
                     var = SCIPvarGetTransVar(edgevars[i][j][2]);
               }
               if( var != NULL && SCIPisGE(scip, SCIPvarGetUbGlobal(var), clustering[j][c2] * clustering[i][c]) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
               {
                  if( SCIPisEQ(scip, 1.0, clustering[j][c2] * clustering[i][c]) )
                     SCIP_CALL( SCIPsetSolVal( scip, sol, var, 1.0 ) );
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** checks if the assignment is finished, i.e. all columns have exactly one 1 and rest 0 values */
static
SCIP_Bool isPartition(
   SCIP_Real**                 clustering,  /**< The matrix containing the clustering */
   int                   nbins,              /**< The number of bins */
   int                   ncluster            /**< The number of clusters */
)
{
   int colsum;
   int i;
   int c;
   SCIP_Bool validassignment = TRUE;
   for( i = 0; i < nbins; ++i )
   {
      colsum = 0;
      for( c = 0; c < ncluster; ++c )
      {
         if( clustering[i][c] == -1 )
            validassignment = FALSE;
         colsum += clustering[i][c];
      }
      if( colsum != 1 )
         validassignment = FALSE;
   }
   return validassignment;
}
#endif

/**
 *  Initialize the q-matrix from a given (possibly incomplete) clusterassignment
 */
static
void computeIrrevMat(
   SCIP_Real**           clustering,         /**< The matrix containing the clusterassignment */
   SCIP_Real**           qmatrix,            /**< The matrix with the return-values, in each cell is the transition probability between two clusters */
   SCIP_Real**           cmatrix,            /**< The transition-matrix containg the probability-data */
   int                   nbins,              /**< The number of bins */
   int                   ncluster            /**< The number of possible clusters */
)
{
   int i,j,k,l;
   for( k = 0; k < ncluster; ++k )
   {
      for( l = 0; l < ncluster; ++l )
      {
         qmatrix[k][l] = 0;
         for( i = 0; i < nbins; ++i )
         {
            for( j = 0; j < nbins; ++j )
            {
               /* As -1 and 0 are both interpreted as 0, this check is necessary. Compute x_ik*x_jl*c_ij */
               qmatrix[k][l] += cmatrix[i][j] * clustering[i][k] * clustering[j][l];
            }
         }
      }
   }
}

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
   SCIP_Real** qmatrix;
   SCIP_Bool feasible;

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

   SCIPallocClearMemoryArray(scip, &clustering , nbins);
   SCIPallocClearMemoryArray(scip, &qmatrix, nbins);
   for( i = 0; i < nbins; ++i )
   {
      SCIPallocClearMemoryArray(scip, &clustering[i], ncluster);
      SCIPallocClearMemoryArray(scip, &qmatrix[i], ncluster);
   }
   /* for each bin, set the assignment with the highest lp-value to 1, the rest to 0 */
   for( i = 0; i < nbins; ++i )
   {
      maxlpval = 0;
      maxcluster = -1;
      for( k = 0; k < ncluster; ++k )
      {
         if( SCIPisGT(scip, SCIPvarGetLPSol(binvars[i][k]), maxlpval) )
         {
            maxlpval = SCIPvarGetLPSol(binvars[i][k]);
            maxcluster =  k;
         }
      }
      assert(maxcluster >= 0);
      clustering[i][maxcluster] = 1.0;
   }
   assert(isPartition(clustering, nbins, ncluster));

   computeIrrevMat(clustering, qmatrix, SCIPspaGetCmatrix(scip), nbins, ncluster);

   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
   assignVars(scip, sol, clustering, nbins, ncluster, qmatrix);
   SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, TRUE, TRUE, TRUE, &feasible) );
   if( feasible )
      *result = SCIP_FOUNDSOL;
   else
      *result = SCIP_DIDNOTFIND;

   /** free allocated memory */

   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeMemoryArray(scip, &clustering[i]);
      SCIPfreeMemoryArray(scip, &qmatrix[i]);
   }
   SCIPfreeMemoryArray(scip, &clustering);
   SCIPfreeMemoryArray(scip, &qmatrix);
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
