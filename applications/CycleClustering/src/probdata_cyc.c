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

/**@file   probdata_cyc.c
 * @brief  problem data for cycle clustering problem
 * @author Leon Eifler
 *
 * This file implements the problem data for the cycle clustering problem.
 *
 * The problem data contains original transition matrix, the scaling parameter that appears in the objective function,
 * and all variables that appear in the problem.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "probdata_cyc.h"

#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_quadratic.h"
#include "scip/var.h"
#include <assert.h>

struct SCIP_ProbData
{
   SCIP_VAR****          edgevars;           /**< variables for the edges (pairs of states) inside and between consecutive clusters */
   SCIP_DIGRAPH*         edgegraph;          /**< digraph which tells us which variables are actually created */
   SCIP_VAR***           binvars;            /**< variable-matrix belonging to the state(bin)-cluster assignment */
   SCIP_Real**           cmatrix;            /**< matrix to save the transition matrix */
   SCIP_Real             scale;              /**< the weight-factor for the coherence in the objective function */
   char                  model_variant;      /**< the model that is used. w for weighted objective, e for normal edge-representation */
   int                   nbins;              /**< the number of states */
   int                   ncluster;           /**< the number of clusters */
};

/** Check if the clustering has exactly one state in every cluster. */
SCIP_Bool isPartition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           solclustering,      /**< matrix with the clustering */
   int                   nbins,              /**< the number of bins */
   int                   ncluster            /**< the number of clusters */
   )
{
   int i;
   int j;

   /* check if the assignment violates paritioning, e.g. because we are in a subscip */
   for( i = 0; i < nbins; ++i )
   {
      SCIP_Real sum = 0.0;

      for( j = 0; j < ncluster; ++j )
      {
         if( !SCIPisIntegral(scip, solclustering[i][j]) )
            return FALSE;
         if( !SCIPisZero(scip, solclustering[i][j]) )
            sum += solclustering[i][j];
      }
      if( !SCIPisEQ(scip, sum, 1.0) )
         return FALSE;

   }

   return TRUE;
}

/** Assign the variables in scip according to the found clustering. */
SCIP_RETCODE assignVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the SCIP solution */
   SCIP_Real**           clustering,         /**< the matrix with the clusterassignment */
   int                   nbins,              /**< the number of bins */
   int                   ncluster            /**< the number of cluster */
   )
{
   SCIP_VAR* var;
   SCIP_VAR*** binvars;
   SCIP_VAR****  edgevars;
   int i;
   int j;
   int c;

   binvars = SCIPcycGetBinvars(scip);
   edgevars = SCIPcycGetEdgevars(scip);

   assert(nbins > 0 && ncluster > 0);
   assert(binvars != NULL);
   assert(edgevars != NULL);

   for ( c = 0; c < ncluster; ++c )
   {
      /* set values of state-variables */
      for ( i = 0; i < nbins; ++i )
      {
         if( NULL != binvars[i][c] )
         {
            if( SCIPvarIsTransformed(binvars[i][c]) )
               var = binvars[i][c];
            else
               var = SCIPvarGetTransVar(binvars[i][c] );

            /* check if the clusterassignment is feasible for the variable bounds. If not do not assign the variable */
            if( var != NULL && SCIPisLE(scip, SCIPvarGetLbGlobal(var), clustering[i][c]) &&
               SCIPisGE(scip, SCIPvarGetUbGlobal(var), clustering[i][c]) &&
               SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
            {
               SCIP_CALL( SCIPsetSolVal( scip, sol, var, clustering[i][c]) );
            }

            assert(SCIPisIntegral(scip, clustering[i][c]));
         }
      }

      /* set the value for the edgevariables for each pair of states */
      for( i = 0; i < nbins; ++i )
      {
         for( j = 0; j < i; ++j )
         {
            if( NULL == edgevars[i][j] || NULL == edgevars[j][i])
               continue;

            /* check if bins are in same cluster */
            if( SCIPisEQ(scip, 1.0, clustering[i][c] * clustering[j][c]) )
            {
               var = edgevars[i][j][0];
               if( var != NULL && SCIPisGE(scip, SCIPvarGetUbGlobal(var), clustering[j][c] * clustering[i][c])
                  && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
               {
                  SCIP_CALL( SCIPsetSolVal( scip, sol, var, 1.0 ) );
               }
            }

            /* check if bins are in consecutive clusters */
            else if( SCIPisEQ(scip, 1.0, clustering[i][c] * clustering[j][phi(c, ncluster)]) )
            {
               var = edgevars[i][j][1];
               if( var != NULL && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR &&
                  SCIPisGE(scip, SCIPvarGetUbGlobal(var), clustering[j][phi(c, ncluster)] * clustering[i][c]) )
               {
                  SCIP_CALL( SCIPsetSolVal( scip, sol, var, 1.0 ) );
               }
            }

            else if( SCIPisEQ(scip, 1.0, clustering[j][c] * clustering[i][phi(c, ncluster)]) )
            {
               var = edgevars[j][i][1];
               if( var != NULL && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR &&
                  SCIPisGE(scip, SCIPvarGetUbGlobal(var), clustering[j][c] * clustering[i][phi(c, ncluster)]) )
               {
                  SCIP_CALL( SCIPsetSolVal( scip, sol, var, 1.0 ) );
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** function that returns the successive cluster along the cycle */
int phi(
   int                   k,                  /**< the cluster */
   int                   ncluster            /**< the number of clusters*/
   )
{
   assert(k < ncluster && k >= 0);
   assert(ncluster > 0);

   return (k+1) % ncluster;
}

/** function that returns the predecessor-cluster along the cycle */
int phiinv(
   int                   k,                  /**< the cluster */
   int                   ncluster            /**< the number of clusters */
   )
{
   assert(k < ncluster && k >= 0);
   assert(ncluster > 0);

   if( k - 1 < 0 )
      return ncluster - 1;
   else
      return k - 1;
}

/** creates all the variables for the problem.  The constraints are added later, depending on the model that is used */
static
SCIP_RETCODE createVariables(
   SCIP*                 scip,               /**< SCIP data Structure */
   SCIP_PROBDATA*        probdata            /**< the problem data */
   )
{
   int i;
   int j;
   int c;
   int edgetype;
   int nbins = probdata->nbins;
   int ncluster = probdata->ncluster;
   char varname[SCIP_MAXSTRLEN];

   /* create variables for bins */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->binvars), nbins) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->binvars[i]), ncluster) ); /*lint !e866*/
   }

   for( i = 0; i < nbins; ++i )
   {
      for( c = 0; c < ncluster; ++c )
      {
         (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d_%d", i, c);
         SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->binvars[i][c], varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
         SCIP_CALL( SCIPaddVar(scip, probdata->binvars[i][c]) );
      }
   }

   /* create variables for the edges in each cluster combination. Index 0 are edges within cluster, 1 edges between
    * consequtive clusters and 2 edges between non-consequtive clusters
    */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->edgevars), nbins) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(probdata->edgevars[i]), nbins) ); /*lint !e866*/

      for( j = 0; j < nbins; ++j )
      {
         if( i == j || (SCIPisZero(scip, (probdata->cmatrix[i][j] - probdata->cmatrix[j][i]))
            && SCIPisZero(scip, (probdata->cmatrix[i][j] + probdata->cmatrix[j][i]) )) )
            continue;

         SCIP_CALL( SCIPdigraphAddArc(probdata->edgegraph, i, j, NULL) );

         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(probdata->edgevars[i][j]), 3) ); /*lint !e866*/

         for( edgetype = 0; edgetype < 3; ++edgetype )
         {
            if( edgetype == 0 && i < j )
               continue;

            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "y_%d_%d_%d", i, j, edgetype);
            SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->edgevars[i][j][edgetype], varname,
               0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
            SCIP_CALL( SCIPaddVar(scip, probdata->edgevars[i][j][edgetype]) );
         }
      }
   }

   return SCIP_OKAY;
}

/** create the problem without variable amount of clusters, use simpler non-facet-defining inequalities */
static
SCIP_RETCODE createProbSimplifiedTest(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata            /**< the problem data */
   )
{
   int i;
   int j;
   int c1;
   int c2;
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS* temp;
   int nbins = probdata->nbins;
   int ncluster = probdata->ncluster;
   SCIP_Real scale;

   SCIP_CALL( SCIPgetRealParam(scip, "cycleclustering/scale_coherence", &scale) );
   probdata->scale = scale;

   /*  create constraints */

   /* create the set-partitioning constraints of the bins */
   for( i = 0; i < nbins; ++i )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "setpart_%d", i+1 );
      SCIP_CALL( SCIPcreateConsSetpart(scip, &temp, consname, 0, NULL, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
         FALSE, FALSE, FALSE) );

      for ( c1 = 0; c1 < ncluster; ++c1 )
      {
         SCIP_CALL( SCIPaddCoefSetppc(scip, temp, probdata->binvars[i][c1]) );
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

   /* create constraints for the edge-cut variables */
   SCIPinfoMessage(scip, NULL, "Using edge-representation with simplified structure. Fixed amount of cluster. \n");
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( probdata->edgevars[i][j] == NULL )
            continue;

         /* these edges are not within a cluster */
         SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][0],
            (probdata->cmatrix[i][j] + probdata->cmatrix[j][i]) * scale) );

         /* these are the edges that are between consequtive clusters */
         SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][1], (probdata->cmatrix[i][j] - probdata->cmatrix[j][i])) );
         SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[j][i][1], (probdata->cmatrix[j][i] - probdata->cmatrix[i][j])) );

         /* create constraints that determine when the edge-variables have to be non-zero */
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            /* constraints for edges within clusters */
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d", i+1, j+1, c1+1 );
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], -1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], 1.0) );

            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );

            /* constraints for edges between clusters */
            for( c2 = 0; c2 < ncluster; ++c2 )
            {
               SCIP_VAR* var;

               if( c2 == c1 )
                  continue;

               if( c2 == c1 + 1 || ( c2 == 0 && c1 == ncluster -1) )
                  var = probdata->edgevars[i][j][1];
               else if( c2 == c1 - 1 || ( c1 == 0 && c2 == ncluster -1) )
                  var = probdata->edgevars[j][i][1];
               else
                  var = probdata->edgevars[i][j][2];

               /* if two bins are in a different cluster -> the corresponding edge must be cut */
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_inclusters_%d_%d", i+1, j+1, c1+1, c2+1 );
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

               SCIP_CALL( SCIPaddCoefLinear(scip, temp, var, -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c2], 1.0) );

               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );
            }
         }
      }
   }

   /*  only one cluster-pair at the same time can be active for an edge*/
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( NULL == probdata->edgevars[i][j]  )
            continue;

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "sumedge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, 1.0, 1.0 ) );

         for( c1 = 0; c1 < 3; ++c1 )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1], 1.0) );
         }
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][1], 1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }

   /* add constraint that ensures that each cluster is used */
   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cluster_%d_used", c1+1 );
      SCIP_CALL( SCIPcreateConsBasicLogicor(scip, &temp, consname, 0, NULL) );

      for( i = 0; i < nbins; ++i )
      {
         SCIP_CALL( SCIPaddCoefLogicor(scip, temp, probdata->binvars[i][c1]) );
      }

      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

   return SCIP_OKAY;
}

/** create the problem without variable amount of clusters, using three edge-variables for each pair of states.
 *  This is the tested default version.
 */
static
SCIP_RETCODE createProbSimplified(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata            /**< The problem data */
   )
{
   int i;
   int j;
   int c1;
   int nbins = probdata->nbins;
   int ncluster = probdata->ncluster;
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS* temp;
   SCIP_Real scale;

   SCIP_CALL( SCIPgetRealParam(scip, "cycleclustering/scale_coherence", &scale) );
   probdata->scale = scale;

   /* create constraints */

   /* create the set-partitioning constraints of the bins */
   for( i = 0; i < nbins; ++i )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "setpart_%d", i+1);
      SCIP_CALL( SCIPcreateConsSetpart(scip, &temp, consname, 0, NULL, TRUE, TRUE, TRUE, TRUE, TRUE,
         FALSE, FALSE, FALSE, FALSE, FALSE) );

      for ( c1 = 0; c1 < ncluster; ++c1 )
      {
         SCIP_CALL( SCIPaddCoefSetppc(scip, temp, probdata->binvars[i][c1]) );
      }

      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

   /* create constraints for the edge-variables */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
      "Using edge-representation with simplified structure. No variable amount of cluster. \n");

   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( NULL == probdata->edgevars[i][j] )
            continue;

         /* the general formulation is needed if there are more than 3 clusters. In the case of three clusters the
          * formulation is simplified
          */
         if( ncluster > 3 )
         {
            /* these edges are within a cluster */
            SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][0],
               (probdata->cmatrix[i][j] + probdata->cmatrix[j][i]) * scale) );

            /* these are the edges that are between consequtive clusters */
            SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][1],
               (probdata->cmatrix[i][j] - probdata->cmatrix[j][i])) );
            SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[j][i][1],
               (probdata->cmatrix[j][i] - probdata->cmatrix[i][j])) );

            /* create constraints that determine when the edge-variables have to be non-zero*/
            for( c1 = 0; c1 < ncluster; ++c1 )
            {
               /* constraints for edges within clusters and between clusters*/
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d", i, j, c1);
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][phiinv(c1, ncluster)], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][phi(c1, ncluster)], -1.0) );

               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_part2_%d", i, j, c1);
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][phi(c1, ncluster)], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][phiinv(c1, ncluster)], -1.0) );

               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d_%d", i, j, c1, phi(c1, ncluster));
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][1], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][phi(c1, ncluster)], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][phi(c1, ncluster)], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], -1.0) );

               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d_%d", i, j, c1, phiinv(c1, ncluster));
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][1], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][phiinv(c1,ncluster)], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][phiinv(c1, ncluster)], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], -1.0) );

               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );
            }
         }
         /* some variables become obsolete with three clusters */
         else
         {
            /* these are the edges that are between consequtive clusters */
            SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][1],
               (probdata->cmatrix[i][j] - probdata->cmatrix[j][i])
               - (probdata->cmatrix[i][j] + probdata->cmatrix[j][i]) * scale) );
            SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[j][i][1],
               (probdata->cmatrix[j][i] - probdata->cmatrix[i][j])
                  - (probdata->cmatrix[i][j] + probdata->cmatrix[j][i]) * scale) );

            SCIP_CALL( SCIPaddOrigObjoffset(scip, (probdata->cmatrix[i][j] + probdata->cmatrix[j][i]) * scale) );

            /* create constraints that determine when the edge-variables have to be non-zero*/
            for( c1 = 0; c1 < ncluster; ++c1 )
            {
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d", i, j, c1);
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][1], -1.0) );

               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][phi(c1, ncluster)], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][phiinv(c1, ncluster)], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][phiinv(c1, ncluster)], -1.0) );

               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d", j, i, c1);
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][1], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][phi(c1, ncluster)], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][phiinv(c1, ncluster)], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][phiinv(c1, ncluster)], -1.0) );

               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );
            }
         }
      }
   }

   /*  only one cluster-pair at the same time can be active for an edge*/
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( NULL == probdata->edgevars[i][j] || (SCIPisZero(scip, (probdata->cmatrix[i][j] - probdata->cmatrix[j][i]))
            && SCIPisZero(scip, (probdata->cmatrix[i][j] + probdata->cmatrix[j][i]) * scale) ) )
            continue;

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "sumedge_%d_%d",  i, j);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0) );

         for( c1 = 0; c1 < 2; ++c1 )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1], 1.0) );
         }

         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][1], 1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }

   /* add constraint that ensures that each cluster is used */
   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cluster_%d_used", c1);
      SCIP_CALL( SCIPcreateConsBasicLogicor(scip, &temp, consname, 0, NULL) );

      for( i = 0; i < nbins; ++i )
      {
         SCIP_CALL( SCIPaddCoefLogicor(scip, temp, probdata->binvars[i][c1]) );
      }

      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

   return SCIP_OKAY;
}

/** create the problem without variable amount of clusters, using quadratic formulations. This is inferior to the
 * simplified variant. Only useful for comparing relaxations.
 */
static
SCIP_RETCODE createProbQP(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata            /**< The problem data */
   )
{
   SCIP_VAR** edgevars;
   SCIP_CONS* temp;
   SCIP_Real scale;
   char varname[SCIP_MAXSTRLEN];
   char consname[SCIP_MAXSTRLEN];
   int i;
   int j;
   int c1;
   int c;
   int nbins = probdata->nbins;
   int ncluster = probdata->ncluster;

   SCIP_CALL( SCIPgetRealParam(scip, "cycleclustering/scale_coherence", &scale) );
   probdata->scale = scale;
   /* create variables for bins */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->binvars), nbins) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->binvars[i]), ncluster) ); /*lint !e866*/
   }

   for( i = 0; i < nbins; ++i )
   {
      for( c = 0; c < ncluster; ++c )
      {
         (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d_%d", i, c);
         SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->binvars[i][c], varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
         SCIP_CALL( SCIPaddVar(scip, probdata->binvars[i][c]) );
      }
   }

   /* create variables for the edges in each cluster combination. Index 0 are edges within cluster, 1 edges between
    * consequtive clusters and 2 edges between non-consequtive clusters
    */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(edgevars), (SCIP_Longint) ncluster * 2) );

   for( i = 0; i < 2 * ncluster; ++i )
   {
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "f_%d", i);

      SCIP_CALL( SCIPcreateVarBasic(scip, &edgevars[i], varname, -SCIPinfinity(scip), SCIPinfinity(scip),
         1.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, edgevars[i]) );
   }

   /* create variables for the edges in each cluster combination. Index 0 are edges within cluster, 1 edges between
    * consequtive clusters and 2 edges between non-consequtive clusters
    */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(probdata->edgevars), nbins) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(probdata->edgevars[i]), nbins) ); /*lint !e866*/

      for( j = 0; j < nbins; ++j )
         probdata->edgevars[i][j] = NULL;
   }

   /*
    *  create constraints
    */

   /* create the set-partitioning constraints of the bins */
   for( i = 0; i < nbins; ++i )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "setpart_%d", i+1);
      SCIP_CALL( SCIPcreateConsSetpart(scip, &temp, consname, 0, NULL, TRUE, TRUE, TRUE, TRUE, TRUE,
         FALSE, FALSE, FALSE, FALSE, FALSE) );

      for ( c1 = 0; c1 < ncluster; ++c1 )
      {
         SCIP_CALL( SCIPaddCoefSetppc(scip, temp, probdata->binvars[i][c1]) );
      }

      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

   /* add constraint that ensures that each cluster is used */
   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cluster_%d_used", c1);
      SCIP_CALL( SCIPcreateConsBasicLogicor(scip, &temp, consname, 0, NULL) );

      for( i = 0; i < nbins; ++i )
      {
         SCIP_CALL( SCIPaddCoefLogicor(scip, temp, probdata->binvars[i][c1]) );
      }

      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

   for( c = 0; c < ncluster; ++c)
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "irrev_%d", c);
      SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &temp, consname, 0, NULL, NULL, 0, NULL, NULL, NULL,
         -SCIPinfinity(scip), 0.0) );

      SCIP_CALL( SCIPaddLinearVarQuadratic(scip, temp, edgevars[c], 1.0) );

      for( i = 0; i < nbins; ++i )
      {
         for( j = 0; j < nbins; ++j )
         {
            if( i != j )
            {
               SCIP_CALL( SCIPaddBilinTermQuadratic(scip, temp, probdata->binvars[i][c],
                  probdata->binvars[j][phi(c,ncluster)], probdata->cmatrix[j][i] - probdata->cmatrix[i][j]) );
            }
         }
      }

      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

   for( c = 0; c < ncluster; ++c )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "coh_%d", c);
      SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &temp, consname, 0, NULL, NULL, 0, NULL, NULL, NULL,
         -SCIPinfinity(scip), 0.0) );

      SCIP_CALL( SCIPaddLinearVarQuadratic(scip, temp, edgevars[c+ncluster], 1.0) );

      for( i = 0; i < nbins; ++i)
      {
         for( j = 0; j < nbins; ++j )
         {
            if( i > j )
            {
               SCIP_CALL( SCIPaddBilinTermQuadratic(scip, temp, probdata->binvars[i][c], probdata->binvars[j][c],
                  -scale * (probdata->cmatrix[i][j] + probdata->cmatrix[j][i])) );
            }
         }
      }

      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

   for (i = 0; i < 2*ncluster; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &edgevars[i]) );
   }

   SCIPfreeBlockMemoryArray(scip, &edgevars, (SCIP_Longint) 2 * ncluster);

   return SCIP_OKAY;
}


/** create the problem with variable amount of clusters. Very large number of constraints not viable for large scale
 * problems.
 */
static
SCIP_RETCODE createProbOnlyEdge(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata            /**< The problem data */
   )
{
   SCIP_CONS* temp;
   SCIP_Real scale;
   SCIP_Real sign[3][3];
   char consname[SCIP_MAXSTRLEN];
   int nbins = probdata->nbins;
   int i;
   int j;
   int k;
   int l;

   SCIP_CALL( SCIPgetRealParam(scip, "cycleclustering/scale_coherence", &scale) );
   probdata->scale = scale;

   for( i = 0; i < 3; ++i )
   {
      for( j = 0; j < 3; ++j )
      {
         sign[i][j] = 1.0;
      }
      sign[i][i] = -1.0;
   }
   /*
    *  create constraints
    */

   /* create constraints for the edge-cut variables */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
      "Using edge-representation with simplified structure. No variable amount of cluster. \n");

   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < nbins; ++j )
      {
         /* set the objective weight for the edge-variables */

         /* these edges are not within a cluster */
         if( j < i )
         {
            SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][0], (probdata->cmatrix[i][j] + probdata->cmatrix[j][i]) * scale) );
            /* these are the edges that are between consequtive clusters */
            SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][1], (probdata->cmatrix[i][j] - probdata->cmatrix[j][i])) );
            SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[j][i][1], (probdata->cmatrix[j][i] - probdata->cmatrix[i][j])) );

            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d", i+1, j+1);
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][1], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][1], 1.0) );

            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );
         }

         for( k =  0; k < nbins; ++k )
         {
            if( i == k || i == j || k == j )
               continue;

            if( k < j && j < i )
            {
               for( l = 0; l < 3; l++ )
               {
                  (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "tri_%d_%d_%d", i, j, k);
                  SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], sign[l][0]) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][k][0], sign[l][1]) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][k][0], sign[l][2]) );

                  SCIP_CALL( SCIPaddCons(scip, temp) );
                  SCIP_CALL( SCIPreleaseCons(scip, &temp) );
               }
            }

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "tri_%d_%d_%d", i, j, k);
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[MAX(i,j)][MIN(i,j)][0], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][k][1], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][k][1], -1.0) );

            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "tri_%d_%d_%d", i, j, k);
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[MAX(i,j)][MIN(i,j)][0], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[k][j][1], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[k][i][1], -1.0) );

            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "tri_%d_%d_%d", i, j, k);
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[MAX(i,j)][MIN(i,j)][0], -1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][k][1], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][k][1], 1.0) );

            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "tri_%d_%d_%d", i, j, k);
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[MAX(i,j)][MIN(i,j)][0], -1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[k][j][1], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[k][i][1], 1.0) );

            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );
         }
      }
   }

   return SCIP_OKAY;
}

/** Scip callback to transform the problem */
static
SCIP_DECL_PROBTRANS(probtransCyc)
{
   int i;
   int j;
   int c;
   int edgetype;
   int nbins = sourcedata->nbins;
   int ncluster = sourcedata->ncluster;

   assert(scip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, targetdata) );

   (*targetdata)->nbins = sourcedata->nbins;
   (*targetdata)->ncluster = sourcedata->ncluster;
   (*targetdata)->model_variant = sourcedata->model_variant;
   (*targetdata)->scale = sourcedata->scale;

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*targetdata)->cmatrix), nbins) );
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*targetdata)->cmatrix[i]), nbins) ); /*lint !e866*/
   }
   /* copy the matrizes */
   for ( i = 0; i < nbins; ++i )
   {
      for ( j = 0; j < nbins; ++j )
      {
         (*targetdata)->cmatrix[i][j] = sourcedata->cmatrix[i][j];
      }
   }

   /* copy the variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*targetdata)->binvars), nbins) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->edgevars), nbins) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->edgevars[i]), nbins) ); /*lint !e866*/
      for( j = 0; j < nbins; ++j )
      {
         if( sourcedata->edgevars[i][j] == NULL || i == j)
            continue;
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->edgevars[i][j]), 3) ); /*lint !e866*/
         for( edgetype = 0; edgetype < 3; ++edgetype )
         {
            if( edgetype == 0 && i < j )
               continue;
            if( sourcedata->edgevars[i][j][edgetype] != NULL )
            {
               SCIP_CALL( SCIPtransformVar(scip, sourcedata->edgevars[i][j][edgetype],
                  &((*targetdata)->edgevars[i][j][edgetype])) );
            }
            else
               ((*targetdata)->edgevars[i][j][edgetype]) = NULL;
         }
      }
   }

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*targetdata)->binvars[i]), ncluster) ); /*lint !e866*/

      for( c = 0; c < ncluster; ++c )
      {
         if( sourcedata->binvars[i][c] != NULL )
         {
            SCIP_CALL( SCIPtransformVar(scip, sourcedata->binvars[i][c], &((*targetdata)->binvars[i][c])) );
         }
         else
            (*targetdata)->binvars[i][c] = NULL;
      }
   }

   SCIP_CALL( SCIPcopyDigraph(scip, &((*targetdata)->edgegraph), sourcedata->edgegraph) );

   return SCIP_OKAY;
}

/** delete-callback method of scip */
static
SCIP_DECL_PROBDELORIG(probdelorigCyc)
{
   int c;
   int edgetype;
   int i;
   int j;

   assert(probdata != NULL);
   assert(*probdata != NULL);

   /* release all the variables */

   /* binvars */
   for ( c = 0; c < (*probdata)->nbins; ++c )
   {
      for ( i = 0; i < (*probdata)->ncluster; ++i )
      {
         if( (*probdata)->binvars[c][i] != NULL )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &((*probdata)->binvars[c][i])) );
         }
      }
      SCIPfreeBlockMemoryArray( scip, &((*probdata)->binvars[c]), (*probdata)->ncluster);
   }

   /* cut-edge vars */
   for ( i = 0; i < (*probdata)->nbins; ++i )
   {
      for( j = 0; j < (*probdata)->nbins; ++j )
      {
         if( (*probdata)->edgevars[i][j] != NULL && j != i )
         {
            for ( edgetype = 0; edgetype < 3; ++edgetype )
            {
               if( edgetype == 0 && i < j )
                  continue;

               if( (*probdata)->edgevars[i][j][edgetype] != NULL )
               {
                  SCIP_CALL( SCIPreleaseVar( scip, &((*probdata)->edgevars[i][j][edgetype])) );
               }
            }

            SCIPfreeBlockMemoryArray(scip, &((*probdata)->edgevars[i][j]), 3);
         }
      }

      SCIPfreeBlockMemoryArray(scip, &((*probdata)->edgevars[i]), (*probdata)->nbins);
   }

   SCIPfreeBlockMemoryArray(scip, &((*probdata)->edgevars), (*probdata)->nbins);
   SCIPfreeBlockMemoryArray(scip, &((*probdata)->binvars), (*probdata)->nbins);

   SCIPdigraphFreeComponents((*probdata)->edgegraph);
   SCIPdigraphFree(&((*probdata)->edgegraph));

   for ( i = 0; i < (*probdata)->nbins; ++i )
   {
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->cmatrix[i]), (*probdata)->nbins);
   }
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->cmatrix, (*probdata)->nbins);

   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}

/** scip-callback to delete the transformed problem */
static
SCIP_DECL_PROBDELTRANS(probdeltransCyc)
{
   int c;
   int edgetype;
   int i;
   int j;

   assert(probdata != NULL);
   assert(*probdata != NULL);

   /* release all the variables */

   /* binvars */
   for ( i = 0; i < (*probdata)->nbins; ++i )
   {
      for ( c = 0; c < (*probdata)->ncluster; ++c )
      {
         if( (*probdata)->binvars[i][c] != NULL )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &((*probdata)->binvars[i][c])) );
         }
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->binvars[i]), (*probdata)->ncluster);
   }

   /* cut-edge vars */
   for ( i = 0; i < (*probdata)->nbins; ++i )
   {
      for( j = 0; j < (*probdata)->nbins; ++j )
      {
         if( (*probdata)->edgevars[i][j] != NULL && j != i )
         {
            for ( edgetype = 0; edgetype < 3; ++edgetype )
            {
               if( 0 == edgetype && j > i )
                  continue;

               if( (*probdata)->edgevars[i][j][edgetype] != NULL )
               {
                  SCIP_CALL( SCIPreleaseVar( scip, &((*probdata)->edgevars[i][j][edgetype])) );
               }
            }

            SCIPfreeBlockMemoryArray(scip, &((*probdata)->edgevars[i][j]), 3);
         }
      }

      SCIPfreeBlockMemoryArray(scip, &((*probdata)->edgevars[i]), (*probdata)->nbins);
   }

   SCIPfreeBlockMemoryArray(scip, &((*probdata)->edgevars), (*probdata)->nbins);
   SCIPfreeBlockMemoryArray(scip, &((*probdata)->binvars), (*probdata)->nbins);

   SCIPdigraphFreeComponents((*probdata)->edgegraph);
   SCIPdigraphFree(&((*probdata)->edgegraph));

   for ( i = 0; i < (*probdata)->nbins; ++i )
   {
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->cmatrix[i]), (*probdata)->nbins);
   }
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->cmatrix, (*probdata)->nbins);

   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}

/** callback method of scip for copying the probdata  */
static
SCIP_DECL_PROBCOPY(probcopyCyc)
{
   SCIP_Bool success;
   SCIP_VAR* var;
   int nbins;
   int ncluster;
   int edgetype;
   int i;
   int j;
   int c;

   assert(scip != NULL);
   assert(sourcescip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);

   /* set up data */
   SCIP_CALL( SCIPallocBlockMemory(scip, targetdata) );

   nbins = sourcedata->nbins;
   ncluster = sourcedata->ncluster;
   success = FALSE;

   (*targetdata)->nbins = nbins;
   (*targetdata)->ncluster = ncluster;
   (*targetdata)->model_variant = sourcedata->model_variant;
   (*targetdata)->scale = sourcedata->scale;

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*targetdata)->cmatrix), nbins) );
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*targetdata)->cmatrix[i]), nbins) ); /*lint !e866*/
   }
   /* copy the matrizes */
   for ( i = 0; i < nbins; ++i )
   {
      for ( j = 0; j < nbins; ++j )
      {
         (*targetdata)->cmatrix[i][j] = sourcedata->cmatrix[i][j];
      }
   }

   /* copy the variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*targetdata)->binvars), nbins) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->edgevars), nbins) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->edgevars[i]), nbins) ); /*lint !e866*/

      for( j = 0; j < nbins; ++j )
      {
         if( (sourcedata)->edgevars[i][j] == NULL || j == i )
            continue;

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*targetdata)->edgevars[i][j]), 3) );  /*lint !e866*/

         for( edgetype = 0; edgetype < 3; ++edgetype )
         {
            if( edgetype == 0 && j > i )
               continue;

            if( sourcedata->edgevars[i][j][edgetype] != NULL )
            {
               SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->edgevars[i][j][edgetype], &var) );

               if( SCIPvarIsActive(var) )
               {
                  SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->edgevars[i][j][edgetype]),
                     varmap, consmap, global, &success) );

                  assert(success);
                  assert((*targetdata)->edgevars[i][j][edgetype] != NULL);

                  SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->edgevars[i][j][edgetype]) );
               }
               else
                  (*targetdata)->edgevars[i][j][edgetype] = NULL;
            }
            else
               (*targetdata)->edgevars[i][j][edgetype] = NULL;
         }
      }
   }

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*targetdata)->binvars[i]), ncluster) ); /*lint !e866*/

      for( c = 0; c < ncluster; ++c )
      {
         if( sourcedata->binvars[i][c] != NULL )
         {
            SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->binvars[i][c], &var) );

            if( SCIPvarIsActive(var) )
            {
               SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->binvars[i][c]),
                  varmap, consmap, global, &success) );

               assert(success);
               assert((*targetdata)->binvars[i][c] != NULL);

               SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->binvars[i][c]) );
            }
            else
               (*targetdata)->binvars[i][c] = NULL;
         }
         else
            (*targetdata)->binvars[i][c] = NULL;
      }
   }

   SCIP_CALL( SCIPcopyDigraph(scip, &((*targetdata)->edgegraph), sourcedata->edgegraph) );

   assert(success);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/**
 * Create the probdata for an cyc-clustering problem
 */
SCIP_RETCODE SCIPcreateProbCyc(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */
   int                   nbins,              /**< number of bins */
   int                   ncluster,           /**< number of cluster */
   SCIP_Real**           cmatrix             /**< The transition matrix */
   )
{
   SCIP_PROBDATA* probdata = NULL;
   int i;
   int j;
   char model;

   assert(nbins > 0);  /* at least one node */
   assert(ncluster <= nbins);

   SCIP_CALL( SCIPcreateProbBasic(scip, name) );

   /* Set up the problem */
   SCIP_CALL( SCIPallocBlockMemory(scip, &probdata) );

   SCIP_CALL( SCIPcreateDigraph(scip, &(probdata->edgegraph), nbins) );

   /* allocate memory for the matrizes and create them from the edge-arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->cmatrix), nbins) );
   for ( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->cmatrix[i]), nbins) ); /*lint !e866*/
      for( j = 0; j < nbins; ++j )
      {
         probdata->cmatrix[i][j] = cmatrix[i][j];
      }
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Creating problem: %s \n", name);

   /* create variables */
   probdata->nbins=nbins;
   probdata->ncluster=ncluster;

   SCIP_CALL( SCIPgetCharParam(scip, "cycleclustering/model", &model) );

   /* create constraints depending on model selection */
   switch( model )
   {
   case 's':
      SCIP_CALL( createVariables(scip, probdata) );
      SCIP_CALL( createProbSimplified(scip, probdata) );
      break;
   case 't':
      SCIP_CALL( createVariables(scip, probdata) );
      SCIP_CALL( createProbSimplifiedTest(scip, probdata) );
      break;
   case 'e':
      SCIP_CALL( createVariables(scip, probdata) );
      SCIP_CALL( createProbOnlyEdge(scip, probdata) );
      break;
   case 'q':
      SCIP_CALL( createProbQP(scip, probdata) );
      break;
   default:
      SCIPABORT();
      break;
   }

   /** add callback methods to scip */
   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigCyc) );
   SCIP_CALL( SCIPsetProbCopy(scip, probcopyCyc) );
   SCIP_CALL( SCIPsetProbData(scip, probdata) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransCyc) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransCyc) );

   return SCIP_OKAY;
}

/**  Getter methods for the various parts of the probdata */


/** Returns the transition matrix*/
SCIP_Real** SCIPcycGetCmatrix(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);

   return probdata->cmatrix;
}

/** Returns the number of states */
int SCIPcycGetNBins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);

   return probdata->nbins;
}

/** Returns the number of clusters*/
int SCIPcycGetNCluster(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip!= NULL);

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);

   return probdata->ncluster;
}

/** Returns the state-variable-matrix*/
SCIP_VAR*** SCIPcycGetBinvars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip!= NULL);

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);
   assert(probdata->binvars != NULL);

   return probdata->binvars;
}

/** Returns the scaling parameter*/
SCIP_Real SCIPcycGetScale(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip!= NULL);

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);

   return probdata->scale;
}

/** Returns the edge variables */
SCIP_VAR**** SCIPcycGetEdgevars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip!= NULL);

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);
   assert(probdata->edgevars != NULL);

   return probdata->edgevars;
}

/** return one specific edge variable */
SCIP_VAR* getEdgevar(
   SCIP_VAR****          edgevars,           /**< edgevar data structure*/
   int                   state1,             /**< first state */
   int                   state2,             /**< second state */
   int                   direction           /**< direction, 0 = incluster, 1 = forward */
   )
{
   assert(edgevars != NULL);
   assert(edgevars[state1] != NULL);
   assert(edgevars[state1][state2] != NULL);
   assert(edgevars[state1][state2][direction] != NULL);

   return edgevars[state1][state2][direction];
}

/** check for an array of states, if all possible edge-combinations exist */
SCIP_Bool edgesExist(
   SCIP_VAR****          edgevars,           /**< edgevar data structure */
   int*                  states,             /**< state array */
   int                   nstates             /**< size of state array */
   )
{
   int i;
   int j;

   assert(edgevars != NULL);
   assert(states != NULL);

   for( i = 0; i < nstates; ++i )
   {
      assert(edgevars[states[i]] != NULL);

      for( j = 0; j < nstates; ++j )
      {
         if( j != i )
         {
            assert(edgevars[states[j]] != NULL);

            if( edgevars[states[i]][states[j]] == NULL )
               return FALSE;
         }
      }
   }

   return TRUE;
}

/** Returns the edge-graph */
SCIP_DIGRAPH* SCIPcycGetEdgeGraph(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip!= NULL);

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);
   assert(probdata->edgegraph != NULL);

   return probdata->edgegraph;
}


/** print the model-values like coherence in the clusters and transition-probabilities between clusters that are not
 *  evident from the scip-solution
 */
SCIP_RETCODE SCIPcycPrintSolutionValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< The solution containg the values */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_Real value;
   SCIP_Real objvalue = 0;
   SCIP_Real flow = 0;
   SCIP_Real coherence = 0;
   int c1;
   int c2;
   int c3;
   int i;
   int j;

   assert(scip!= NULL);

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\nDisplaying aggregated solution data: \n");

   for( c1  = 0; c1 < probdata->ncluster; ++c1 )
   {
      value = 0;

      for( i = 0; i < probdata->nbins; ++i )
      {
         for( j = 0; j < probdata->nbins; ++j )
         {
            if( j == i )
               continue;

            value+= probdata->cmatrix[i][j] * SCIPgetSolVal(scip, sol, probdata->binvars[i][c1])
               * SCIPgetSolVal(scip, sol, probdata->binvars[j][c1]);
         }
      }

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Coherence in cluster %d  :  %f \n", c1 + 1, value);
      coherence += value;
      objvalue += probdata->scale * value;
   }

   for( c1 = 0; c1 < probdata->ncluster; ++c1 )
   {
      for( c2 = 0; c2 < probdata->ncluster; ++c2 )
      {
         value = 0;

         for( i = 0; i < probdata->nbins; ++i )
         {
            for( j = 0; j < probdata->nbins; ++j )
            {
               value+= probdata->cmatrix[i][j] * SCIPgetSolVal(scip, sol, probdata->binvars[i][c1]) *
                  SCIPgetSolVal(scip, sol, probdata->binvars[j][c2]);
            }
         }
      }

      c3 = (c1 + 1) % probdata->ncluster;
      value = 0;

      for( i = 0; i < probdata->nbins; ++i )
      {
         for( j = 0; j < probdata->nbins; ++j )
         {
            value+= (probdata->cmatrix[i][j] - probdata->cmatrix[j][i])
               * SCIPgetSolVal(scip, sol, probdata->binvars[i][c1]) * SCIPgetSolVal(scip, sol, probdata->binvars[j][c3]);
         }
      }

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
         "  Net flow from %d to %d    :  %f \n", c1, phi(c1, probdata->ncluster), value);

      flow += value;
      objvalue += value;
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Total coherence         :  %f \n",  coherence);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Total net flow          :  %f \n",  flow);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Total objective value   :  %f \n",  objvalue);

   return SCIP_OKAY;
}
