
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
#define HEUR_FREQ             25
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE          /**< does the heuristic use a secondary SCIP instance? */

/*
 * Data structures
 */


/*
 * Local methods
 */

/**  calculate the current epsI-value for a q-matrix */
static
SCIP_Real getIrrevBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           qmatrix,            /**< The irreversibility matrix*/
   int                   ncluster            /**< The number of cluster*/
)
{
   SCIP_Real epsI = SCIP_INVALID; /* this is the an upper bound to irreversibility */
   int i;
   int j;
   for( i = 0; i < ncluster; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         SCIP_Real temp;
         temp = REALABS(qmatrix[i][j] - qmatrix[j][i]);
         if( SCIPisPositive(scip, temp) && SCIPisLT(scip, temp, epsI) )
         {
            epsI = temp;
         }
      }
   }
   /* if we have no transitions at all then irreversibility should be set to 0 */
   if( epsI == SCIP_INVALID )
   {
      epsI = 0.0;
   }
   return epsI;
}

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
   SCIP_Real**           qmatrix             /**< The irreversibility matrix */
)
{
   int i,j;
   int c;
   int c2;
   SCIP_Real epsI;
   SCIP_Real q1;
   SCIP_Real q2;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_VAR* resultant;
   SCIP_VAR* targetvar;
   SCIP_VAR** indvars;
   SCIP_VAR** absvars;
   SCIP_VAR*** binvars;
   SCIP_VAR*****  edgevars;

   assert(nbins > 0 && ncluster > 0);

   indvars = SCIPspaGetIndvars(scip);
   absvars = SCIPspaGetAbsvars(scip);
   binvars = SCIPspaGetBinvars(scip);
   targetvar = SCIPspaGetTargetvar(scip);
   edgevars = SCIPspaGetEdgevars(scip);

   epsI = getIrrevBound(scip, qmatrix, ncluster);
   for ( c = 0; c < ncluster; ++c )
   {
      int isempty = 1;
      for( j = 0; j < nbins; ++j )
      {
         if( clustering[j][c] > 0 )
            isempty = 0;
      }
      /* set indicatorvar whether cluster is nonempty */

      if( indvars[c] != NULL && !isempty && SCIPisEQ(scip, SCIPvarGetUbGlobal(indvars[c]), 1.0) && SCIPvarGetStatus(indvars[c]) != SCIP_VARSTATUS_MULTAGGR )
         SCIP_CALL( SCIPsetSolVal(scip, sol, indvars[c], 1.0) );
      else if( indvars[c] != NULL && SCIPisZero(scip, SCIPvarGetLbGlobal(indvars[c])) && SCIPvarGetStatus(indvars[c]) != SCIP_VARSTATUS_MULTAGGR )
         SCIP_CALL( SCIPsetSolVal(scip, sol, indvars[c], 0.0) );

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
         for( j = 0; j < nbins; ++j )
         {
            if( j == i )
               continue;
            if( j < i )
            {
               if( NULL == edgevars[i][j][c][c] )
                  continue;
               if( SCIPvarIsTransformed(edgevars[i][j][c][c]) )
                  var = edgevars[i][j][c][c];
               else
                  var = SCIPvarGetTransVar(edgevars[i][j][c][c]);
               if( var != NULL && SCIPisGE(scip, SCIPvarGetUbGlobal(var), clustering[j][c] * clustering[i][c]) && SCIPisLE(scip, SCIPvarGetLbGlobal(var), clustering[j][c] * clustering[i][c]) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
                  SCIP_CALL( SCIPsetSolVal( scip, sol, var, clustering[j][c] * clustering[i][c]  ) );
            }
            for( c2 = 0; c2 < c; ++c2 )
            {
               if( NULL == edgevars[i][j][c][c2] )
                  continue;
               if( SCIPvarIsTransformed(edgevars[i][j][c][c2]) )
                  var = edgevars[i][j][c][c2];
               else
                  var = SCIPvarGetTransVar(edgevars[i][j][c][c]);
               if( var != NULL && SCIPisGE(scip, SCIPvarGetUbGlobal(var), clustering[j][c2] * clustering[i][c]) && SCIPisLE(scip, SCIPvarGetLbGlobal(var), clustering[j][c2] * clustering[i][c]) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
                  SCIP_CALL( SCIPsetSolVal( scip, sol, edgevars[i][j][c][c2], clustering[j][c2] * clustering[i][c] ) );
            }
         }
      }

      /* set variables of absolute value and q-variables */
      for ( i = 0; i <= c; ++i )
      {
         q1 = qmatrix[c][i];
         q2 = qmatrix[i][c];
         /* set the absolute-value vars. Always check if the found values are feasible for the bounds of the variables */
         if( absvars[i +ncluster * c] == NULL || absvars[c + ncluster * i] == NULL || SCIPvarGetStatus(absvars[i + ncluster * c]) == SCIP_VARSTATUS_MULTAGGR || SCIPvarGetStatus(absvars[c + ncluster * i]) == SCIP_VARSTATUS_MULTAGGR )
            continue;
         if( SCIPisGT(scip, q1, q2) )
         {
            if( SCIPisEQ(scip, SCIPvarGetUbGlobal(absvars[i + ncluster * c]), 1.0) )
               SCIP_CALL( SCIPsetSolVal(scip, sol, absvars[i + ncluster * c], 1.0) );
            if( SCIPisZero(scip, SCIPvarGetLbGlobal(absvars[c + ncluster * i])) )
               SCIP_CALL( SCIPsetSolVal(scip, sol, absvars[c + ncluster * i], 0.0) );
         }
         else if( SCIPisGT(scip, q2, q1) )
         {
            if( SCIPisEQ(scip, SCIPvarGetUbGlobal(absvars[c + ncluster * i]), 1.0) )
               SCIP_CALL( SCIPsetSolVal(scip, sol, absvars[c + ncluster * i], 1.0) );
            if( SCIPisZero(scip, SCIPvarGetLbGlobal(absvars[i + ncluster * c])) )
               SCIP_CALL( SCIPsetSolVal(scip, sol, absvars[i + ncluster * c], 0.0) );
         }
         else
         {
            SCIP_CALL( SCIPsetSolVal(scip, sol, absvars[c + ncluster * i], SCIPvarGetLbGlobal(absvars[c + ncluster * i]) ) );
            SCIP_CALL( SCIPsetSolVal(scip, sol, absvars[i + ncluster * c], SCIPvarGetLbGlobal(absvars[i + ncluster * c]) ) );
         }
      }
   }
   /* set the value of the target-function variable */
   if( SCIPisGT(scip, epsI, SCIPvarGetLbGlobal(targetvar)) && SCIPisLT(scip, epsI, SCIPvarGetUbGlobal(targetvar)) && SCIPvarGetStatus(targetvar) != SCIP_VARSTATUS_MULTAGGR )
      SCIP_CALL( SCIPsetSolVal( scip, sol, targetvar, epsI) );

   /* set the and variables that scip introduces in presoving */
   {
      int nconshdlrs = SCIPgetNConshdlrs(scip);
      SCIP_CONSHDLR** conshdlrs = SCIPgetConshdlrs(scip);
      for( i = 0; i < nconshdlrs ; ++i )
      {
         if ( 0 == strcmp( SCIPconshdlrGetName( conshdlrs[i] ) , "and" ) )
         {
            SCIP_CONS** cons = SCIPconshdlrGetConss(conshdlrs[i]);
            int ncons = SCIPconshdlrGetNConss(conshdlrs[i]);
            /* for each constraint, assign the value of the resultant */
            for( j = 0; j < ncons; ++j )
            {
               int k;
               int nvars = SCIPgetNVarsAnd(scip, cons[j]);
               SCIP_Real temp = 1.0;

               vars = SCIPgetVarsAnd(scip, cons[j]);
               resultant = SCIPgetResultantAnd(scip, cons[j]);
               assert( resultant != NULL );
               for( k = 0; k < nvars; ++k )
               {
                  temp = temp * SCIPgetSolVal(scip, sol, vars[k]);
                  assert(SCIPisIntegral(scip, temp));
               }
               if( SCIPvarGetStatus(resultant) != SCIP_VARSTATUS_MULTAGGR )
                  SCIP_CALL( SCIPsetSolVal(scip, sol, resultant, temp) );
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
