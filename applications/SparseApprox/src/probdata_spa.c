/*
 * @file probdata_spa.c
 *
 *  Created on: Dec 1, 2015
 *      Author: bzfeifle
 */
#include "probdata_spa.h"
#include "scip/cons_orbitope.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_linear.h"
#include "scip/cons_indicator.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_logicor.h"
#include "scip/var.h"

#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>


struct SCIP_ProbData
{
   SCIP_VAR**            indicatorvar;     /**< variables indicating whether a cluster is empty */
   SCIP_VAR****          edgevars;         /**< variables for the edges cut by a partioning */
   SCIP_VAR***           binvars;          /**< variable-matrix belonging to the bin-cluster assignment */
   SCIP_Real**           cmatrix;          /**< matrix to save the transition matrix */
   SCIP_Real             scale;        /**< The weight for the coherence in the objective function */
   SCIP_Real             coherence;
   char                  model_variant;     /**< The model that is used. w for weighted objective, e for normal edge-representation */
   int                   nbins;
   int                   ncluster;
};

/** Creates all the variables for the problem. The constraints are added later, depending on the model that is used */
static
SCIP_RETCODE createVariables(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata            /**< The problem data */
)
{
   int i;
   int j;
   int c1;
   int nbins = probdata->nbins;
   int ncluster = probdata->ncluster;
   char varname[SCIP_MAXSTRLEN];


   /* create variables for bins */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->binvars), nbins) );
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->binvars[i]), ncluster) );
   }

   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < ncluster; ++j )
      {
         (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d_%d", i+1 , j+1 );
         SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->binvars[i][j], varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
         SCIP_CALL( SCIPaddVar(scip, probdata->binvars[i][j]) );
         SCIP_CALL( SCIPvarChgBranchPriority(probdata->binvars[i][j], 5) );
      }
   }
   SCIP_CALL( SCIPchgVarLbGlobal(scip, probdata->binvars[nbins - 1][0], 1.0) );


   /* create variables for the edges in each cluster combination. Index 0 are edges within cluster, 1 edges between consequtive clusters and 2 edges between non-consequtive clusters */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->edgevars), nbins) );
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->edgevars[i]), nbins) );
      for( j = 0; j < nbins; ++j )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->edgevars[i][j]), 3) );
         for( c1 = 0; c1 < 3; ++c1 )
         {
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "y_%d_%d_%d", i+1, j+1, c1+1 );
            SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->edgevars[i][j][c1], varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT) );
            SCIP_CALL( SCIPaddVar(scip, probdata->edgevars[i][j][c1]) );
         }
      }
   }
   /* add indicator variable that marks after which cluster the cycle is closed. Only used in model with variable amount of clusters*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->indicatorvar), ncluster) );
   for( i = 0; i < ncluster; ++i )
   {
      (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "ind_%d", i+1  );
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->indicatorvar[i], varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(scip, probdata->indicatorvar[i]) );
      SCIP_CALL( SCIPvarChgBranchPriority(probdata->indicatorvar[i], 10) );
   }

   return SCIP_OKAY;
}

/** create the problem without variable amount of clusters */
static
SCIP_RETCODE createProbSimplified(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata            /**< The problem data */
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

   SCIPgetRealParam(scip, "scale_coherence", &scale);
   probdata->scale = scale;

   /*
    *  create constraints
    */

   /* create the set-partitioning constraints of the bins */
   for( i = 0; i < nbins; ++i )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "setpart_%d", i+1 );
      SCIP_CALL( SCIPcreateConsSetpart(scip, &temp, consname, 0, NULL, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      for ( c1 = 0; c1 < ncluster; ++c1 )
      {
         SCIP_CALL( SCIPaddCoefSetppc(scip, temp, probdata->binvars[i][c1]) );
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

   /* create constraints for the edge-cut variables */
   SCIPinfoMessage(scip, NULL, "Using edge-representation with simplified structure. No variable amount of cluster. \n");
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( 0 == (probdata->cmatrix[i][j] - probdata->cmatrix[j][i]) )
            continue;
         /* set the objective weight for the edge-variables */

         /* these edges are not within a cluster */
         SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][0], (probdata->cmatrix[i][j] + probdata->cmatrix[j][i]) * scale ) );
         /* these are the edges that are between consequtive clusters */
         SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][1], (probdata->cmatrix[i][j] - probdata->cmatrix[j][i])) );
         SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[j][i][1], (probdata->cmatrix[j][i] - probdata->cmatrix[i][j])) );

         /* create constraints that determine when the edge-variables have to be non-zero */
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            /* constraints for edges within clusters */
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d", i+1, j+1, c1+1 );
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], -1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], 1.0) );
            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );

            /*(void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d", i+1, j+1, c1+1 );
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], -1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], 1.0) );
            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );

            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d", i+1, j+1, c1+1 );
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], -1.0) );
            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );*/

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

               /* if two bins are in a different cluster, then the corresponding edge must be cut*/
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_inclusters_%d_%d", i+1, j+1, c1+1, c2+1 );
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, var, -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c2], 1.0) );
               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );

               /*(void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_inclusters_%d_%d", i+1, j+1, c1+1, c2+1 );
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, var, 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c2], 1.0) );
               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_inclusters_%d_%d", i+1, j+1, c1+1, c2+1 );
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, var, 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c2], -1.0) );
               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );*/

            }
         }
      }
   }

   /*  only one cluster-pair at the same time can be active for an edge*/
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( NULL == probdata->edgevars[i][j][0] )
            continue;
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "sumedge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, 0.0, 1.0 ) );
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
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cluster_%d_used", c1+1 );
      SCIP_CALL( SCIPcreateConsBasicLogicor(scip, &temp, consname, 0, NULL) );
      /*      SCIP_CALL( SCIPaddCoefSetppc(scip, temp, probdata->indicatorvar[0]) );*/
      for( i = 0; i < nbins; ++i )
      {
         SCIP_CALL( SCIPaddCoefLogicor(scip, temp, probdata->binvars[i][c1]) );
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

   return SCIP_OKAY;
}

/** Create the model with the objective: irreversibility + weight * coherence with variable amount of clusters. Does not work with the current formulation */
static
SCIP_RETCODE createProbWeightedCoh(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata            /**< The problem data */
)
{
   int i;
   int j;
   int k;
   int c1;
   int c2;
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS* temp;
   int nbins = probdata->nbins;
   int ncluster = probdata->ncluster;
   SCIP_Real scale;

   SCIPgetRealParam(scip, "scale_coherence", &scale);
   probdata->scale = scale;

   /*
    *  create constraints
    */

   /* create the set-partitioning constraints of the bins */
   for( i = 0; i < nbins; ++i )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "setpart_%d", i+1 );
      SCIP_CALL( SCIPcreateConsSetpart(scip, &temp, consname, 0, NULL, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      for ( c1 = 0; c1 < ncluster; ++c1 )
      {
         SCIP_CALL( SCIPaddCoefSetppc(scip, temp, probdata->binvars[i][c1]) );
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

   /* create constraints for the edge-cut variables */
   SCIPinfoMessage(scip, NULL, "Using edge-representation with variable amount of clusters. \n");
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         /* these edges are not at all cut */
         SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][0], (probdata->cmatrix[i][j] + probdata->cmatrix[j][i]) * scale ) );
         /* these are the edges that are cut forward in the cycle */
         SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][1], (probdata->cmatrix[i][j] - probdata->cmatrix[j][i])) );
         /* these are the edges that are cut backward in the cycle */
         SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][2], (probdata->cmatrix[j][i] - probdata->cmatrix[i][j])) );

         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            /** in this variant the coherence enters the objective. Because of this we add the corresponding edgevars directly to the objective function */

            /* if two bins are in the same cluster, then the corresponding edge can not be cut */
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d", i+1, j+1, c1+1 );
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], -1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], 1.0) );
            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );


            /* set the objective weight for the edge-variables */
            for( c2 = 0; c2 < ncluster; ++c2 )
            {
               SCIP_VAR* var;
               if( c2 == c1 )
                  continue;
               if( c2 == c1 + 1 )
                  var = probdata->edgevars[i][j][1];
               else if( c2 == c1 - 1 )
                  var = probdata->edgevars[i][j][2];
               else
                  var = probdata->edgevars[i][j][3];

               if( c2 != 0 && c1 != 0)
               {
                  /* if two bins are in a different cluster, then the corresponding edge must be cut*/
                  (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_inclusters_%d_%d", i+1, j+1, c1+1, c2+1 );
                  SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, var, -1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c2], 1.0) );
                  SCIP_CALL( SCIPaddCons(scip, temp) );
                  SCIP_CALL( SCIPreleaseCons(scip, &temp) );
               }
               else
               {
                  int tmp = MAX(c1,c2);
                  /* if two bins are in a different cluster, then the corresponding edge must be cut*/
                  (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_inclusters_%d_%d", i+1, j+1, c1+1, c2+1 );
                  SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, var, -1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[tmp], -1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c2], 1.0) );
                  SCIP_CALL( SCIPaddCons(scip, temp) );
                  SCIP_CALL( SCIPreleaseCons(scip, &temp) );
               }

            }

            /* add constraints that are active when c1 is the return cluster */
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d_0", i+1, j+1, c1+1 );
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 2.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][1], -1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][0], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], 1.0) );
            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );

            /* add constraints that are active when c1 is the return cluster */
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_0_%d", i+1, j+1, c1+1 );
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 2.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][2], -1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][0], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], 1.0) );
            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );
         }
      }
   }

   /*  only one cluster-pair at the same time can be active for an edge*/
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( NULL == probdata->edgevars[i][j][0] )
            continue;
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "sumedge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, 1.0, 1.0 ) );
         for( c1 = 0; c1 < 4; ++c1 )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1], 1.0) );
         }

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }

   /* add constraint that ensures that each cluster is used */

   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cluster_%d_used", c1+1 );
      SCIP_CALL( SCIPcreateConsBasicLogicor(scip, &temp, consname, 0, NULL) );
      for( c2 = 0; c2 < c1; ++c2 )
      {
         SCIP_CALL( SCIPaddCoefLogicor(scip, temp, probdata->indicatorvar[c2]) );
      }
      for( i = 0; i < nbins; ++i )
      {
         SCIP_CALL( SCIPaddCoefLogicor(scip, temp, probdata->binvars[i][c1]) );
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }


   /* add constraint that nothing is used beyond return */
   for( i = 0; i < nbins; ++i )
   {
      for( k = 0; k < ncluster; ++k )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "return_%d_%d",  i+1, k+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0 ) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][k], 1.0) );
         for( c1 = 0; c1 < k; ++c1 )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], 1.0) );
         }
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }

   /* only one return point viable */
   (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "part_return");
   SCIP_CALL( SCIPcreateConsBasicSetpart(scip, &temp, consname, 0, NULL) );
   for( k = 0; k < ncluster; ++k )
   {
      SCIP_CALL( SCIPaddCoefSetppc(scip, temp, probdata->indicatorvar[k]) );
   }

   SCIP_CALL( SCIPaddCons(scip, temp) );
   SCIP_CALL( SCIPreleaseCons(scip, &temp) );

   return SCIP_OKAY;
}


/** create the problem without variable amount of clusters */
static
SCIP_RETCODE createProbExperimental(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata            /**< The problem data */
)
{
   int i;
   int j;
   int k;
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS* temp;
   int nbins = probdata->nbins;
   SCIP_Real scale;

   SCIPgetRealParam(scip, "scale_coherence", &scale);
   probdata->scale = scale;
   /*
    *  create constraints
    */

   /* create constraints for the edge-cut variables */
   SCIPinfoMessage(scip, NULL, "Using edge-representation with simplified structure. No variable amount of cluster. \n");
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < nbins; ++j )
      {
         /* these are the edges that are between consequtive clusters */
         if( i > j )
            SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][0], scale * (probdata->cmatrix[i][j] + probdata->cmatrix[j][i])) );
         SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][1], (probdata->cmatrix[i][j] - probdata->cmatrix[j][i])) );
         for( k = 0; k < nbins; ++k )
         {
            if( i > j && j > k )
            {
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_%d", i+1, j+1, k+1 );
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE ) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][k][0], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][k][0], -1.0) );
               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_%d", i+1, j+1, k+1 );
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE ) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][k][0], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][k][0], 1.0) );
               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_%d", i+1, j+1, k+1 );
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE ) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][k][0], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][k][0], 1.0) );
               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );
            }
            if( !(i > k && j < k) && !(i < k && j > k) )
            {
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_%d", i+1, j+1, k+1 );
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE ) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[MAX(i,j)][MIN(i,j)][0], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][k][1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][k][1], -1.0) );
               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );
            }
            if( !(i > k && i < j) && !(i < k && i > j))
            {
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_%d", i+1, j+1, k+1 );
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE ) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[MAX(j,k)][MIN(j,k)][0], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][k][1], 1.0) );
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
         if( NULL == probdata->edgevars[i][j][0] )
            continue;
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "sumedge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, 0.0, 1.0 ) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][1], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][1], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], 1.0) );


         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }

   return SCIP_OKAY;
}

/** create the problem without variable amount of clusters */
static
SCIP_RETCODE createProbExperimental2(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata            /**< The problem data */
)
{
   int i;
   int j;
   int l;
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS* temp;
   int nbins = probdata->nbins;
   probdata->ncluster = 2;
   for( i = 0; i < nbins; ++i )
   {
      SCIP_Real source = -(probdata->cmatrix[i][(i+1) % nbins] - probdata->cmatrix[(i+1) % nbins][i]);
      SCIP_Real sink = -(probdata->cmatrix[(i+1) % nbins][i] - probdata->cmatrix[i][(i+1) % nbins]);
      SCIP_CALL( SCIPchgVarObj(scip, probdata->binvars[i][0], source ) );
      SCIP_CALL( SCIPchgVarObj(scip, probdata->binvars[i][1], sink ) );
      for( j = 0; j < nbins; ++j )
      {
         source += MAX((probdata->cmatrix[i][j] - probdata->cmatrix[j][i]), 0);
         sink += MAX((probdata->cmatrix[j][i] - probdata->cmatrix[i][j]), 0);
      }

      if( source < 0 )
         SCIP_CALL( SCIPchgVarUb(scip, probdata->binvars[i][0], 0.0) );
      if( sink < 0 )
         SCIP_CALL( SCIPchgVarUb(scip, probdata->binvars[i][1], 0.0) );


      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "edge_%d",  i+1 );
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, 0.0, 1.0 ) );
      SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][0], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][1], 1.0) );
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );

      for( j = 0; j < nbins; ++j )
      {
         SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][0], (probdata->cmatrix[i][j] - probdata->cmatrix[j][i]) ) );
         /*if( j < i )
         {
            SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][1], 0.001 * (probdata->cmatrix[i][j] + probdata->cmatrix[j][i]) ) );
            SCIP_CALL( SCIPchgVarObj(scip, probdata->edgevars[i][j][2], 0.001 * (probdata->cmatrix[i][j] + probdata->cmatrix[j][i]) ) );
         }*/
      }
   }
   /*for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "edge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0 ) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][1], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][0], -1.0) );
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "edge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0 ) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][1], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][0], -1.0) );
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "edge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0 ) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][2], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][1], -1.0) );
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "edge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0 ) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][2], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][1], -1.0) );
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }*/

   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < nbins; ++j )
      {
         if( j == i )
            continue;
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "edge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &temp, consname, 0, NULL, NULL, 0, NULL, NULL, NULL, 0.0, 0.0 ) );
         SCIP_CALL( SCIPaddLinearVarQuadratic(scip, temp, probdata->edgevars[i][j][0], -1.0) );
         SCIP_CALL( SCIPaddBilinTermQuadratic(scip, temp, probdata->binvars[i][0], probdata->binvars[j][1], 1.0) );
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

        /* (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "edge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0 ) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][1], -1.0) );
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "edge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0 ) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], -1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][0], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][1], 1.0) );
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );*/
         for( l = 0; l < nbins; ++l )
         {
            if( j == l || l == i )
               continue;
            /*(void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "edge_%d_%d",  i+1, j+1 );
            SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 2.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][l][0], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][l][0], 1.0) );
            SCIP_CALL( SCIPaddCons(scip, temp) );
            SCIP_CALL( SCIPreleaseCons(scip, &temp) );*/
         }
      }
   }

   return SCIP_OKAY;
}

/** Scip callback to transform the problem */
static
SCIP_DECL_PROBTRANS(probtransSpa)
{
   int i;
   int j;
   int c1;
   int nbins = sourcedata->nbins;
   int ncluster = sourcedata->ncluster;
   assert(scip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, targetdata) );

   (*targetdata)->nbins = sourcedata->nbins;
   (*targetdata)->ncluster = sourcedata->ncluster;
   (*targetdata)->coherence = sourcedata->coherence;
   (*targetdata)->model_variant = sourcedata->model_variant;
   (*targetdata)->scale = sourcedata->scale;


   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->cmatrix), nbins) );
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->cmatrix[i]), nbins) );
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
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->binvars), nbins) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->edgevars), nbins) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->indicatorvar), ncluster) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->edgevars[i]), nbins) );
      for( j = 0; j < nbins; ++j )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->edgevars[i][j]), 3) );
         for( c1 = 0; c1 < 3; ++c1 )
         {
            if( sourcedata->edgevars[i][j][c1] != NULL )
            {
               SCIP_CALL( SCIPtransformVar(scip, sourcedata->edgevars[i][j][c1], &((*targetdata)->edgevars[i][j][c1])) );
            }
            else
               ((*targetdata)->edgevars[i][j][c1]) = NULL;
         }
      }
   }
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->binvars[i]), ncluster) );

      for( j = 0; j < ncluster; ++j )
      {
         if( sourcedata->binvars[i][j] != NULL )
         {
            SCIP_CALL( SCIPtransformVar(scip, sourcedata->binvars[i][j], &((*targetdata)->binvars[i][j])) );
         }
         else
            (*targetdata)->binvars[i][j] = NULL;
      }
   }

   for( i = 0; i < ncluster ; ++i )
   {
      if( sourcedata->indicatorvar[i] != NULL)
         SCIP_CALL( SCIPtransformVar(scip, sourcedata->indicatorvar[i], &((*targetdata)->indicatorvar[i])) );
      else
         (*targetdata)->indicatorvar[i] = NULL;
   }

   return SCIP_OKAY;
}

/** delete-callback method of scip */
static
SCIP_DECL_PROBDELORIG(probdelorigSpa)
{
   int c1;
   int i;
   int j;

   assert(probdata != NULL);
   assert(*probdata != NULL);

   /* release all the variables */

   /* binvars */
   for ( c1 = 0; c1 < (*probdata)->nbins; ++c1 )
   {
      for ( i = 0; i < (*probdata)->ncluster; ++i )
      {
         if( (*probdata)->binvars[c1][i] != NULL )
            SCIP_CALL( SCIPreleaseVar(scip, &((*probdata)->binvars[c1][i])) );
      }
      SCIPfreeMemoryArray( scip, &((*probdata)->binvars[c1]) );
   }
   /* non-emptiness vars */
   for ( c1 = 0; c1 < (*probdata)->ncluster; ++c1 )
   {
      if( (*probdata)->indicatorvar[c1] != NULL )
         SCIP_CALL( SCIPreleaseVar( scip, &((*probdata)->indicatorvar[c1])) );
   }

   /* cut-edge vars */

   for ( i = 0; i < (*probdata)->nbins; ++i )
   {
      for( j = 0; j < (*probdata)->nbins; ++j )
      {
         for ( c1 = 0; c1 < 3; ++c1 )
         {
            if( (*probdata)->edgevars[i][j][c1] != NULL )
               SCIP_CALL( SCIPreleaseVar( scip, &((*probdata)->edgevars[i][j][c1])) );
         }
         SCIPfreeMemoryArray( scip, &((*probdata)->edgevars[i][j]) );
      }
      SCIPfreeMemoryArray( scip, &((*probdata)->edgevars[i]) );
   }

   SCIPfreeMemoryArray( scip, &((*probdata)->edgevars) );
   SCIPfreeMemoryArray(scip, &((*probdata)->binvars));
   SCIPfreeMemoryArray(scip, &((*probdata)->indicatorvar));
   for ( i = 0; i < (*probdata)->nbins; ++i )
   {
      SCIPfreeMemoryArray(scip, &((*probdata)->cmatrix[i]));
   }
   SCIPfreeMemoryArray(scip, &(*probdata)->cmatrix);

   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


/** scip-callback to delete the transformed problem */
static
SCIP_DECL_PROBDELORIG(probdeltransSpa)
{
   int c1;
   int i;
   int j;

   assert(probdata != NULL);
   assert(*probdata != NULL);

   /* release all the variables */

   /* binvars */
   for ( c1 = 0; c1 < (*probdata)->nbins; ++c1 )
   {
      for ( i = 0; i < (*probdata)->ncluster; ++i )
      {
         if( (*probdata)->binvars[c1][i] != NULL )
            SCIP_CALL( SCIPreleaseVar(scip, &((*probdata)->binvars[c1][i])) );
      }
      SCIPfreeMemoryArray( scip, &((*probdata)->binvars[c1]) );
   }
   /* non-emptiness vars */
   for ( c1 = 0; c1 < (*probdata)->ncluster; ++c1 )
   {
      if( (*probdata)->indicatorvar[c1] != NULL )
         SCIP_CALL( SCIPreleaseVar( scip, &((*probdata)->indicatorvar[c1])) );
   }

   /* cut-edge vars */

   for ( i = 0; i < (*probdata)->nbins; ++i )
   {
      for( j = 0; j < (*probdata)->nbins; ++j )
      {
         for ( c1 = 0; c1 < 3; ++c1 )
         {
            if( (*probdata)->edgevars[i][j][c1] != NULL )
               SCIP_CALL( SCIPreleaseVar( scip, &((*probdata)->edgevars[i][j][c1])) );
         }
         SCIPfreeMemoryArray( scip, &((*probdata)->edgevars[i][j]) );
      }
      SCIPfreeMemoryArray( scip, &((*probdata)->edgevars[i]) );
   }

   SCIPfreeMemoryArray( scip, &((*probdata)->edgevars) );
   SCIPfreeMemoryArray(scip, &((*probdata)->binvars));
   SCIPfreeMemoryArray(scip, &((*probdata)->indicatorvar));

   for ( i = 0; i < (*probdata)->nbins; ++i )
   {
      SCIPfreeMemoryArray(scip, &((*probdata)->cmatrix[i]));
   }
   SCIPfreeMemoryArray(scip, &(*probdata)->cmatrix);

   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

/** callback method of scip for copying the probdata  */
static
SCIP_DECL_PROBCOPY(probcopySpa)
{
   int nbins;
   int ncluster;
   int i;
   int j;
   int c1;
   SCIP_Bool success;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(sourcescip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);

   /* set up data */
   SCIP_CALL( SCIPallocMemory(scip, targetdata) );

   nbins = sourcedata->nbins;
   ncluster = sourcedata->ncluster;

   (*targetdata)->nbins = nbins;
   (*targetdata)->ncluster = ncluster;
   (*targetdata)->coherence = sourcedata->coherence;
   (*targetdata)->model_variant = sourcedata->model_variant;
   (*targetdata)->scale = sourcedata->scale;


   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->cmatrix), nbins) );
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->cmatrix[i]), nbins) );
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
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->binvars), nbins) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->edgevars), nbins) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->indicatorvar), ncluster) );


   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->edgevars[i]), nbins) );
      for( j = 0; j < nbins; ++j )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->edgevars[i][j]), 3) );
         for( c1 = 0; c1 < 3; ++c1 )
         {
            if( sourcedata->edgevars[i][j][c1] != NULL )
            {
               SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->edgevars[i][j][c1], &var) );
               if( SCIPvarIsActive(var) )
               {
                  SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->edgevars[i][j][c1]), varmap, consmap, global, &success) );
                  assert(success);
                  assert((*targetdata)->edgevars[i][j][c1] != NULL);
                  SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->edgevars[i][j][c1]) );
               }
               else
                  (*targetdata)->edgevars[i][j][c1] = NULL;
            }
            else
               (*targetdata)->edgevars[i][j][c1] = NULL;
         }
      }
   }
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->binvars[i]), ncluster) );

      for( j = 0; j < ncluster; ++j )
      {
         if( sourcedata->binvars[i][j] != NULL )
         {
            SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->binvars[i][j], &var) );
            if( SCIPvarIsActive(var) )
            {
               SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->binvars[i][j]), varmap, consmap, global, &success) );
               assert(success);
               assert((*targetdata)->binvars[i][j] != NULL);

               SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->binvars[i][j]) );
            }
            else
               (*targetdata)->binvars[i][j] = NULL;
         }
         else
            (*targetdata)->binvars[i][j] = NULL;
      }
   }

   for( i = 0; i < ncluster; ++i )
   {
      {
         if( sourcedata->indicatorvar[i] != NULL )
         {
            SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->indicatorvar[i], &var) );

            if( SCIPvarIsActive(var) )
            {
               SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->indicatorvar[i], &var) );
               SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->indicatorvar[i]), varmap, consmap, global, &success) );
               assert(success);
               assert((*targetdata)->indicatorvar[i] != NULL);
               SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->indicatorvar[i]) );
            }
            else
               (*targetdata)->indicatorvar[i] = NULL;
         }
         else
            (*targetdata)->indicatorvar[i] = NULL;

      }
   }

   assert(success);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/**
 * Create the probdata for an spa-clustering problem
 */
SCIP_RETCODE SCIPcreateProbSpa(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */
   int                   nbins,              /**< number of bins */
   SCIP_Real**           cmatrix             /**< The transition matrix */
)
{
   SCIP_PROBDATA* probdata = NULL;
   int i;
   int j;
   SCIP_Real epsC;
   char model_variant;
   int ncluster;

   assert(nbins > 0);  /* at least one node */

   SCIP_CALL( SCIPcreateProbBasic(scip, name) );

   /* get the parameters for the coherence bound */
   SCIP_CALL( SCIPgetRealParam( scip, "coherence_bound", &epsC ) );

   /* get the model variant we wish to use */
   SCIP_CALL( SCIPgetCharParam(scip, "model_variant", &model_variant) );

   /* get the maximal amount of clusters */
   SCIP_CALL( SCIPgetIntParam(scip, "ncluster", &ncluster) );
   assert( ncluster <= nbins);

   /* Set up the problem */
   SCIP_CALL( SCIPallocMemory( scip, &probdata) );

   /* allocate memory for the matrizes and create them from the edge-arrays */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->cmatrix), nbins) );
   for ( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->cmatrix[i]), nbins) );
      for( j = 0; j < nbins; ++j )
      {
         probdata->cmatrix[i][j] = cmatrix[i][j];
      }
   }

   assert(epsC >= 0 && epsC <= 1.0);
   probdata->coherence = epsC;

   SCIPinfoMessage(scip, NULL, "Creating problem: %s \n", name);

   /* create variables */
   probdata->nbins=nbins;
   probdata->ncluster=ncluster;
   probdata->model_variant = model_variant;

   SCIP_CALL( createVariables(scip, probdata) );

   /* create constraints depending on model selection */

   switch (model_variant)
   {
   case 'w':
      SCIP_CALL( createProbWeightedCoh(scip, probdata) );
      break;
   case 's':
      SCIP_CALL( createProbSimplified(scip, probdata) );
      break;
   case 'e':
      SCIP_CALL( createProbExperimental(scip, probdata) );
      break;
   case 't':
      SCIP_CALL( createProbExperimental2(scip, probdata) );
      break;
   default:
      break;
   }

   /** add callback methods to scip */
   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigSpa) );
   SCIP_CALL( SCIPsetProbCopy(scip, probcopySpa) );
   SCIP_CALL( SCIPsetProbData(scip, probdata) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransSpa) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransSpa) );

   return SCIP_OKAY;
}



/**  Getter methods for the various parts of the probdata */
SCIP_Real** SCIPspaGetCmatrix(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->cmatrix;
}

int SCIPspaGetNrBins(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_PROBDATA* probdata;
   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->nbins;
}

int SCIPspaGetNrCluster(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_PROBDATA* probdata;
   assert(scip!= NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->ncluster;
}

SCIP_VAR*** SCIPspaGetBinvars(
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

SCIP_VAR** SCIPspaGetIndvars(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_PROBDATA* probdata;
   assert(scip!= NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   assert(probdata->indicatorvar != NULL);

   return probdata->indicatorvar;
}


SCIP_Real SCIPspaGetCoherence(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_PROBDATA* probdata;
   assert(scip!= NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->coherence;
}

char SCIPspaGetModel(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_PROBDATA* probdata;
   assert(scip!= NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->model_variant;
}

SCIP_Real SCIPspaGetScale(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_PROBDATA* probdata;
   assert(scip!= NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->scale;
}

SCIP_VAR**** SCIPspaGetEdgevars(
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

/** print the model-values like coherence in the clusters and transition-probabilities between clusters that are not evident from the scip-solution */
SCIP_RETCODE SCIPspaPrintSolutionValues(
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

   assert( scip!= NULL );
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   for( c1  = 0; c1 < probdata->ncluster; ++c1 )
   {
      value = 0;
      for( i = 0; i < probdata->nbins; ++i )
      {
         for( j = 0; j < probdata->nbins; ++j )
         {
            if( j == i )
               continue;
            value+= probdata->cmatrix[i][j] * SCIPgetSolVal(scip, sol, probdata->binvars[i][c1]) * SCIPgetSolVal(scip, sol, probdata->binvars[j][c1]);
         }
      }
      SCIPinfoMessage(scip, NULL, "Coherence in cluster %d is %f \n", c1 + 1, value);
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
               value+= probdata->cmatrix[i][j] * SCIPgetSolVal(scip, sol, probdata->binvars[i][c1]) * SCIPgetSolVal(scip, sol, probdata->binvars[j][c2]);
            }
         }
      }
      c3 = (c1 + 1) % probdata->ncluster;
      value = 0;
      for( i = 0; i < probdata->nbins; ++i )
      {
         for( j = 0; j < probdata->nbins; ++j )
         {
            value+= (probdata->cmatrix[i][j] - probdata->cmatrix[j][i]) * SCIPgetSolVal(scip, sol, probdata->binvars[i][c1]) * SCIPgetSolVal(scip, sol, probdata->binvars[j][c3]);
         }
      }
      SCIPinfoMessage(scip, NULL, "irrev_%d_%d is %f \n", c1 + 1, (c1 + 2) % probdata->ncluster, value);
      flow += value;
      objvalue += value;
   }
   SCIPinfoMessage(scip, NULL, "objvalue is %f \n",  objvalue);
   SCIPinfoMessage(scip, NULL, "Total coherence is %f \n",  coherence);
   SCIPinfoMessage(scip, NULL, "Total net flow is %f \n",  flow);
   return SCIP_OKAY;
}

