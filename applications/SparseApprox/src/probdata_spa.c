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
#include "scip/var.h"

#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>

#define BIG_M 1


struct SCIP_ProbData
{
   SCIP_VAR**            indicatorvar;     /**< variables indicating whether a cluster is empty */
   SCIP_VAR*****         edgevars;         /**< variables for the edges cut by a partioning */
   SCIP_VAR***           binvars;          /**< variable-matrix belonging to the bin-cluster assignment */
   SCIP_VAR*             targetvar;        /**< variable that enters the target-function */
   SCIP_VAR**            absvar;           /**< variables for linearizing the absolute value in the irreversibility constraints */
   SCIP_Real**           cmatrix;          /**< matrix to save the transition matrix */
   SCIP_Real             big_M;
   int                   nbins;
   int                   ncluster;
};

/**Create the transition matrix from a given array of edges and a singular distribution vector */
static
SCIP_RETCODE createCMatrix(
   SCIP*                 scip,               /**< SCIP data structure*/
   SCIP_Real**           cmatrix,            /**< The transition matrix*/
   SCIP_Real**           edges,              /**< Edge-Array*/
   SCIP_Real*            sd,                 /**< Stationary distribution vector*/
   int                   nedges,             /**< The number of edges*/
   int                   nbins               /**< The number of bins*/
)
{
   int i;
   int j;
   SCIP_Real* sum;
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &sum, nbins) );
   /* create the matrix and check if the given data is stochastic matrix*/
   for ( i = 0; i < nbins; ++i )
   {
      for ( j = 0; j < nbins; ++j )
      {
         cmatrix[i][j] = 0.0;
      }
   }

   for ( i = 0; i < nedges; ++i )
   {
      assert((int)edges[i][0] <= nbins && (int)edges[i][1]<= nbins);

      cmatrix[(int)edges[i][0] - 1][(int)edges[i][1] - 1] = sd[(int)edges[i][0] - 1] * edges[i][2];
      sum[(int)edges[i][0] - 1] += edges[i][2];
      /*if( edges[i][2] > 1.0 || edges[i][2] < 0 )
      {
         SCIPerrorMessage( "The given C-matrix is has an entry at %d, %d with illegal value %f \n", (int)edges[i][0], (int)edges[i][1] , edges[i][2] );
         return SCIP_READERROR;
      }*/
   }

   /*for ( i = 0; i < nbins; ++i )
   {
      if( SCIPisGT(scip, fabs(1.0 - sum[i]) , 10e-4) )
      {
         SCIPerrorMessage( "The given T-matrix is not stochastic, the row-sum in row %d is %f \n", i, sum[i] );
         return SCIP_READERROR;
      }
   }*/

   SCIPfreeMemoryArray(scip, &sum);
   return SCIP_OKAY;
}

/** returns the minimal value of the size x size - matrix that is non-zero. Only works on positive matrices */
SCIP_Real getMinNonZero(
   SCIP*                 scip,               /**< SCIP data structure*/
   SCIP_Real**           matrix,             /**< The matrix to be considered*/
   int                   size                /**< The size of the matrix*/
)
{
   int i;
   int j;
   SCIP_Real ret = SCIPinfinity(scip);
   for ( i = 0; i < size; ++i )
   {
      for ( j = 0; j < size; ++j )
      {
         if( matrix[i][j] > 0 && matrix[i][j] < ret )
         {
            ret = matrix[i][j];
         }
      }
   }
   assert(ret > 0);
   return ret;
}

/**  free memory allocated for a matrix of size nbins*nbins */
SCIP_RETCODE freeMatrix(
   SCIP_Real**           matrix,
   int                   nbins
)
{
   int i;
   for ( i = 0; i < nbins; ++i )
   {
      SCIPfreeMemoryArray(scip, &(matrix[i]));
   }
   SCIPfreeMemoryArray(scip, &matrix);
   return SCIP_OKAY;
}

static
SCIP_Real bigM(
   SCIP_Real**cmatrix,
   int nbins)
{
	SCIP_Real max = 0;
	int i,j;
	for( i = 0 ; i < nbins; ++i )
	{
	   for( j = 0; j < i; ++j )
	   {
            max += REALABS(cmatrix[i][j] - cmatrix[j][i]);
      }
	}
	return max;
}


static
SCIP_RETCODE createVariables(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata            /**< The problem data */
)
{
   int i;
   int j;
   int c1;
   int c2;
   int nbins = probdata->nbins;
   int ncluster = probdata->ncluster;
   char varname[SCIP_MAXSTRLEN];
   SCIP_Bool dummy;


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
   SCIP_CALL( SCIPfixVar(scip, probdata->binvars[0][0], 1.0, &dummy, &dummy  ) );
   SCIP_CALL( SCIPfixVar(scip, probdata->binvars[1][2], 0.0, &dummy, &dummy) );

   /* create variables for the edges in each cluster combination */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->edgevars), nbins) );
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->edgevars[i]), nbins) );
      for( j = 0; j < nbins; ++j )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->edgevars[i][j]), ncluster) );
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->edgevars[i][j][c1]), ncluster) );
            for( c2 = 0; c2 < ncluster; ++c2 )
            {
               (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "y_%d_%d_%d_%d", i+1, j+1, c1+1 , c2+1 );
               SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->edgevars[i][j][c1][c2], varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
               SCIP_CALL( SCIPaddVar(scip, probdata->edgevars[i][j][c1][c2]) );
            }
         }
      }
   }
   /* add indicator (1 if cluster is non-empty) and absolute value helper (for linearizing absolute value) variables */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->indicatorvar), ncluster) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->absvar), ncluster * ncluster) );
   for( i = 0; i < ncluster; ++i )
   {
      (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "ind_%d", i+1  );
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->indicatorvar[i], varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(scip, probdata->indicatorvar[i]) );
      for( j = 0; j < ncluster; ++j )
      {
         (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "abs_%d_%d", i+1 , j+1 );
         SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->absvar[i + j * ncluster], varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
         SCIP_CALL( SCIPvarChgBranchPriority(probdata->absvar[i + j * ncluster], 10) );
         SCIP_CALL( SCIPaddVar(scip, probdata->absvar[i + j * ncluster]) );
      }
   }
   /* add the lower bound on irreversibility variable that enters the target function */
   SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->targetvar, "epsI", 0.0, probdata->big_M, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, probdata->targetvar) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE createProbEdgeRep(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata,           /**< The problem data */
   SCIP_Real             epsC                /**< The coherence parameter */
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

   /*
    *  create constraints
    */

   /* create the set-partitioning constraints of the bins */

   for( i = 0; i < nbins; ++i )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "setpack_%d", i+1 );
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, 1.0, 1.0 ) );
      for ( c1 = 0; c1 < ncluster; ++c1 )
      {
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }
   /* create constraints for the edge-cut variables */

   SCIPinfoMessage(scip, NULL, "Using edge-representation \n");
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < nbins; ++j )
      {
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            {
               if( j < i )
               {
                  /* if two bins are in the same cluster, then the corresponding edge can not be cut */
                  (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d", i+1, j+1, c1+1 );
                  SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c1], -1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], 1.0) );
                  SCIP_CALL( SCIPaddCons(scip, temp) );
                  SCIP_CALL( SCIPreleaseCons(scip, &temp) );
               }
               for( c2 = 0; c2 < c1; ++c2 )
               {
                  if( j == i )
                     continue;
                  /* if two bins are in a different cluster, then the corresponding edge must be cut*/
                  (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d_%d", i+1, j+1, c1+1, c2+1 );
                  SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c2], -1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c2], 1.0) );
                  SCIP_CALL( SCIPaddCons(scip, temp) );
                  SCIP_CALL( SCIPreleaseCons(scip, &temp) );
               }
            }
         }
      }
   }

   /*    only one cluster-pair at the same time can be active for an edge*/
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "sumedge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, 0.0, 1.0 ) );
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            for( c2 = 0; c2 < c1; ++c2 )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c2], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][c1][c2], 1.0) );
            }
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c1], 1.0) );
         }
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }

   /* create constraints that manage the emptieness-indicator variables */
   for ( c1 = 0; c1 < ncluster; ++c1 )
   {
      for( i = 0; i < nbins; ++i )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "clopen_lb_%d_%d", c1 + 1, i+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], -1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip,temp,probdata->binvars[i][c1], 1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }

   /* create cohrerence constraints */
   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "coh_%d", c1 + 1 );

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], -epsC) );
      for( i = 0; i < nbins; ++i )
      {
         for ( j = 0; j < i; ++j )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c1], (probdata->cmatrix[i][j])) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c1], (probdata->cmatrix[j][i])) );
         }
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], (probdata->cmatrix[i][i])) );
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

   /* create irreversibility - constraints */
   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      for( c2 = 0; c2 < c1; ++c2 )
      {
         (void)SCIPsnprintf( consname, SCIP_MAXSTRLEN, "irrev1_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), probdata->big_M ) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->targetvar, 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c1 + ncluster * c2], probdata->big_M) );
         /* add all cut edges between cluster c1 and c2 */
         for( i = 0; i < nbins; ++i )
         {
            for( j = 0; j < i; ++j )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c2], (probdata->cmatrix[i][j]) - (probdata->cmatrix[j][i])) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][c1][c2], -(probdata->cmatrix[i][j]) + (probdata->cmatrix[j][i])) );
            }
         }
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );


         (void)SCIPsnprintf( consname, SCIP_MAXSTRLEN, "irrev2_%d_%d", c1 + 1, (c1 + 2) % ncluster );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), probdata->big_M) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->targetvar, 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c2 + ncluster * c1], probdata->big_M) );
         for( i = 0; i < nbins; ++i )
         {
            for( j = 0; j < i; ++j )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c2], (probdata->cmatrix[j][i]) - (probdata->cmatrix[i][j])) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][c1][c2], -(probdata->cmatrix[j][i]) + (probdata->cmatrix[i][j])) );
            }
         }
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }

   /* add constraints for the relation between the absvar and the indicator-variables */
   for ( c1 = 0; c1 < ncluster; ++c1 )
   {
      for( c2 = 0; c2 < c1; ++c2 )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "sum_d1_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c1 + ncluster * c2], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c2 + ncluster * c1], 1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "sum_d2_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c1 + ncluster * c2], -1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c2 + ncluster * c1], -1.0) );

         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c2], 1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }

   /* symmetry between clusters:
    *          c_{i+1}         <= c_{i}
    *    <==>  c_{i+1} - c_{i} <= 0            for all i = 1,...,|C|-1
    */
   for ( c1 = 0; c1 < ncluster-1; ++c1 )
   {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "sym_c%d_c%d", c1+1, c1+2 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], -1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1+1], 1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }

#if 0 /* todo: check that this constraints do exactly */
   /* add more constraints for better separation */
   for ( c1 = 0; c1 < ncluster; ++c1 )
   {
      for( c2 = 0; c2 < c1; ++c2 )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "d1_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c1 + ncluster * c2], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], -1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "d2_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c2 + ncluster * c1], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], -1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "d3_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c1 + ncluster * c2], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c2], -1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "d4_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c2 + ncluster * c1], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c2], -1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

      }
   }
#endif

   /* add constraint that ensures irreversibility between clusters is 0 if we have only one cluster
    * due to numerical reasons we use 1.0 instead of big_M; the constraint basically says:
    *
    *    epsilon <= sum_{k in Clusters} c_{k} - 1.0
    *
    * If more than one cluster is open, the constraint is trivially fulfilled, otherwise epsilon has to be zero.
    */
   (void)SCIPsnprintf( consname, SCIP_MAXSTRLEN, "nontrivial");
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), -1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->targetvar, 1.0));
   for ( c1 = 0; c1 < ncluster; ++c1 )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], -1.0) );
   }

   SCIP_CALL( SCIPaddCons(scip, temp) );
   SCIP_CALL( SCIPreleaseCons(scip, &temp) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE createProbEdgeRepRelax(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata,           /**< The problem data */
   SCIP_Real             epsC                /**< The coherence parameter */
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

   /*
    *  create constraints
    */

   /* create the set-partitioning constraints of the bins */

   for( i = 0; i < nbins; ++i )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "setpack_%d", i+1 );
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, 1.0, 1.0 ) );
      for ( c1 = 0; c1 < ncluster; ++c1 )
      {
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }
   /* create constraints for the edge-cut variables */

   SCIPinfoMessage(scip, NULL, "Using edge-representation \n");
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < nbins; ++j )
      {
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            {
               if( j < i )
               {
                  /* if two bins are in the same cluster, then the corresponding edge can not be cut */
                  (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d", i+1, j+1, c1+1 );
                  SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c1], -1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], 1.0) );
                  SCIP_CALL( SCIPaddCons(scip, temp) );
                  SCIP_CALL( SCIPreleaseCons(scip, &temp) );
               }
               for( c2 = 0; c2 < c1; ++c2 )
               {
                  if( j == i )
                     continue;
                  /* if two bins are in a different cluster, then the corresponding edge must be cut*/
                  (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d_%d", i+1, j+1, c1+1, c2+1 );
                  SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c2], -1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
                  SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c2], 1.0) );
                  SCIP_CALL( SCIPaddCons(scip, temp) );
                  SCIP_CALL( SCIPreleaseCons(scip, &temp) );
               }
            }
         }
      }
   }

   /*    only one cluster-pair at the same time can be active for an edge*/
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "sumedge_%d_%d",  i+1, j+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, 0.0, 1.0 ) );
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            for( c2 = 0; c2 < c1; ++c2 )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c2], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][c1][c2], 1.0) );
            }
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c1], 1.0) );
         }
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }

   /* create constraints that manage the emptieness-indicator variables */
   for ( c1 = 0; c1 < ncluster; ++c1 )
   {
      for( i = 0; i < nbins; ++i )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "clopen_lb_%d_%d", c1 + 1, i+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], -1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip,temp,probdata->binvars[i][c1], 1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }
   /* create cohrerence constraints */
   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "coh_%d", c1 + 1 );

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, epsC, SCIPinfinity(scip)) );
      for( i = 0; i < nbins; ++i )
      {
         for ( j = 0; j < i; ++j )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c1], (probdata->cmatrix[i][j])) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c1], (probdata->cmatrix[j][i])) );
         }
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], (probdata->cmatrix[i][i])) );
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }
   /* create irreversibility - constraints, relaxation with abs inside the sum */
   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      for( c2 = 0; c2 < c1; ++c2 )
      {
         (void)SCIPsnprintf( consname, SCIP_MAXSTRLEN, "irrev" );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0 ) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->targetvar, 1.0) );
         /* add all cut edges between cluster c1 and c2 */
         for( i = 0; i < nbins; ++i )
         {
            for( j = 0; j < i; ++j )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][c1][c2], -REALABS((probdata->cmatrix[i][j]) - (probdata->cmatrix[j][i]))) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[j][i][c1][c2], -REALABS(-(probdata->cmatrix[i][j]) + (probdata->cmatrix[j][i]))) );
            }
         }
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE createProbBinRep(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata,           /**< The problem data */
   SCIP_Real             epsC                /**< The coherence parameter */
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

   /*
    *  create constraints
    */

   /* create the set-partitioning constraints of the bins */

   for( i = 0; i < nbins; ++i )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "setpack_%d", i+1 );
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, 1.0, 1.0 ) );
      for ( c1 = 0; c1 < ncluster; ++c1 )
      {
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }
   /* create constraints for the edge-cut variables */


   /* create constraints that manage the emptieness-indicator variables */
   for ( c1 = 0; c1 < ncluster; ++c1 )
   {
      for( i = 0; i < nbins; ++i )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "clopen_lb_%d_%d", c1 + 1, i+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], -1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip,temp,probdata->binvars[i][c1], 1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }
   SCIPinfoMessage(scip, NULL, "Using bin-representation \n");
   /* create constraints that manage the emptieness-indicator variables */
   for ( c1 = 0; c1 < ncluster; ++c1 )
   {
      for( i = 0; i < nbins; ++i )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "clopen_lb_%d_%d", c1 + 1, i+1 );
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], -1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip,temp,probdata->binvars[i][c1], 1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }
   /*create constraints for the coherence in a cluster */
   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "coh_%d", c1 + 1 );

      SCIP_CALL( SCIPcreateConsQuadratic(scip, &temp, consname, 0, NULL, NULL, 0, NULL, NULL, NULL, 0.0, SCIPinfinity(scip), FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddLinearVarQuadratic(scip, temp, probdata->indicatorvar[c1], -epsC) );
      /* add all the bilinear terms to the constraint */
      for( i = 0; i < nbins; ++i )
      {
         for ( j = 0; j < nbins; ++j )
         {
            SCIP_CALL( SCIPaddBilinTermQuadratic(scip, temp, probdata->binvars[i][c1], probdata->binvars[j][c1], (probdata->cmatrix[i][j]) ) );
         }
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }
   /*create constraints for irreversibility. There are 2 constraints for each pair of clusters due to the linearization of the absolute value.*/
   for ( c1 = 0; c1 < ncluster; ++c1 )
   {
      for ( c2 = 0; c2 < c1; ++c2 )
      {
         (void)SCIPsnprintf( consname, SCIP_MAXSTRLEN, "irrev1_%d_%d", c1 + 1,c2 + 1 );
         /* use the absolute value for the irreversibility-measure */
         SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &temp, consname, 0, NULL, NULL, 0, NULL, NULL, NULL, -SCIPinfinity(scip), probdata->big_M ) );
         SCIP_CALL( SCIPaddLinearVarQuadratic(scip, temp, probdata->targetvar, 1.0) );
         SCIP_CALL( SCIPaddLinearVarQuadratic(scip, temp, probdata->absvar[c1 + ncluster * c2], probdata->big_M) );

         /* add all the bilinear terms to the constraint */
         for( i = 0; i < nbins; ++i )
         {
            for ( j = 0; j < nbins; ++j )
            {
               if( (probdata->cmatrix[i][j])>0 )
               {
                  SCIP_CALL( SCIPaddBilinTermQuadratic(scip, temp, probdata->binvars[i][c1], probdata->binvars[j][c2], (probdata->cmatrix[i][j])) );
                  SCIP_CALL( SCIPaddBilinTermQuadratic(scip, temp, probdata->binvars[i][c2], probdata->binvars[j][c1], -(probdata->cmatrix[i][j])) );
               }
            }
         }
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         /* add the second part of the constraint */
         (void)SCIPsnprintf( consname, SCIP_MAXSTRLEN, "irrev2_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &temp, consname, 0, NULL, NULL, 0, NULL, NULL, NULL, -SCIPinfinity(scip), probdata->big_M) );
         SCIP_CALL( SCIPaddLinearVarQuadratic(scip, temp, probdata->targetvar, 1.0) );
         SCIP_CALL( SCIPaddLinearVarQuadratic(scip, temp, probdata->absvar[c2 + ncluster * c1], probdata->big_M) );

         /* add all the bilinear terms to the constraint */
         for( i = 0; i < nbins; ++i )
         {
            for ( j = 0; j < nbins; ++j )
            {
               SCIP_CALL( SCIPaddBilinTermQuadratic(scip, temp, probdata->binvars[i][c1], probdata->binvars[j][c2], -(probdata->cmatrix[i][j])) );
               SCIP_CALL( SCIPaddBilinTermQuadratic(scip, temp, probdata->binvars[i][c2], probdata->binvars[j][c1], (probdata->cmatrix[i][j])) );
            }
         }
         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }

   /* add constraints for the relation between the absvar and the indicator-variables */
   for ( c1 = 0; c1 < ncluster; ++c1 )
   {
      for( c2 = 0; c2 < c1; ++c2 )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "sum_d1_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c1 + ncluster * c2], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c2 + ncluster * c1], 1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "sum_d2_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c1 + ncluster * c2], -1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c2 + ncluster * c1], -1.0) );

         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c2], 1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );
      }
   }

#if 0 /* todo: check that this constraints do exactly */
   /* add more constraints for better separation */
   for ( c1 = 0; c1 < ncluster; ++c1 )
   {
      for( c2 = 0; c2 < c1; ++c2 )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "d1_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c1 + ncluster * c2], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], -1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "d2_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c2 + ncluster * c1], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], -1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "d3_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c1 + ncluster * c2], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c2], -1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "d4_%d_%d", c1 + 1, c2 + 1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->absvar[c2 + ncluster * c1], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c2], -1.0) );

         SCIP_CALL( SCIPaddCons(scip, temp) );
         SCIP_CALL( SCIPreleaseCons(scip, &temp) );

      }
   }
#endif

   /* add constraint that ensures irreversibility between clusters is 0 if we have only one cluster */
   (void)SCIPsnprintf( consname, SCIP_MAXSTRLEN, "nontrivial");
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), -probdata->big_M) );
   SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->targetvar, 1.0));
   for ( c1 = 0; c1 < ncluster; ++c1 )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->indicatorvar[c1], -probdata->big_M) );
   }
   SCIP_CALL( SCIPaddCons(scip, temp) );
   SCIP_CALL( SCIPreleaseCons(scip, &temp) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE createProbMaxkCut(
   SCIP*                 scip,               /**< SCIP Data Structure */
   SCIP_PROBDATA*        probdata,           /**< The problem data */
   SCIP_Real             epsC                /**< The coherence parameter */
)
{
   int i;
   int j;
   int c1;
   SCIP_Real sum;
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS* temp;
   int nbins = probdata->nbins;
   int ncluster = probdata->ncluster;

   /* create the set-partitioning constraints of the bins */

   for( i = 0; i < nbins; ++i )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "setpack_%d", i+1 );
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, 1.0, 1.0 ) );
      for ( c1 = 0; c1 < ncluster; ++c1 )
      {
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }
   /* create constraints for the edge-cut variables */

   SCIPinfoMessage(scip, NULL, "Using edge-representation for max k-cut. \n");
   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      for( i = 0; i < nbins; ++i )
      {
         for( j = 0; j < nbins; ++j )
         {

            if( j < i )
            {
               /* if two bins are in the same cluster, then the corresponding edge can not be cut */
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d", i+1, j+1, c1+1 );
               SCIP_CALL( SCIPcreateConsLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0][0], -1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[j][c1], 1.0) );
               SCIP_CALL( SCIPaddCons(scip, temp) );
               SCIP_CALL( SCIPreleaseCons(scip, &temp) );
            }
         }
      }
   }

   /* create cohrerence constraints */
   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "coh_%d", c1 + 1 );

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, epsC * ncluster, SCIPinfinity(scip)) );
      for( i = 0; i < nbins; ++i )
      {
         for ( j = 0; j < i; ++j )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0][0], (probdata->cmatrix[i][j])) );
            SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0][0], (probdata->cmatrix[j][i])) );
         }
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->binvars[i][c1], (probdata->cmatrix[i][i])) );
      }
      SCIP_CALL( SCIPaddCons(scip, temp) );
      SCIP_CALL( SCIPreleaseCons(scip, &temp) );
   }
   /* create irreversibility - constraints, relaxation with abs inside the sum */

   sum = 0;
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         sum += REALABS((probdata->cmatrix[i][j]) - (probdata->cmatrix[j][i]));
      }
   }
   (void)SCIPsnprintf( consname, SCIP_MAXSTRLEN, "irrev" );
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &temp, consname, 0, NULL, NULL, -SCIPinfinity(scip), sum ) );
   SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->targetvar, 1.0) );

   /* add all cut edges between cluster c1 and c2 */
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         SCIP_CALL( SCIPaddCoefLinear(scip, temp, probdata->edgevars[i][j][0][0], REALABS((probdata->cmatrix[i][j]) - (probdata->cmatrix[j][i]))) );
      }
   }
   SCIP_CALL( SCIPaddCons(scip, temp) );
   SCIP_CALL( SCIPreleaseCons(scip, &temp) );

   return SCIP_OKAY;
}

static
SCIP_DECL_PROBTRANS(probtransSpa)
{
   int i;
   int j;
   int c1;
   int c2;
   int nbins = sourcedata->nbins;
   int ncluster = sourcedata->ncluster;
   assert(scip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, targetdata) );

   (*targetdata)->nbins = sourcedata->nbins;        /* copy length of array sets */
   (*targetdata)->ncluster = sourcedata->ncluster;            /* copy number of sets saved in array sets */

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
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->absvar), ncluster * ncluster) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->edgevars[i]), nbins) );
      for( j = 0; j < nbins; ++j )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->edgevars[i][j]), ncluster) );
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->edgevars[i][j][c1]), ncluster) );
            for( c2 = 0; c2 < ncluster; ++c2 )
            {
               if( sourcedata->edgevars[i][j][c1][c2] != NULL )
               {
                  SCIP_CALL( SCIPtransformVar(scip, sourcedata->edgevars[i][j][c1][c2], &((*targetdata)->edgevars[i][j][c1][c2])) );
               }
               else
                  ((*targetdata)->edgevars[i][j][c1][c2]) = NULL;
            }
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

   for( i = 0; i < ncluster * ncluster; ++i )
   {
      if( sourcedata->absvar[i] != NULL )
         SCIP_CALL( SCIPtransformVar(scip, sourcedata->absvar[i], &((*targetdata)->absvar[i])) );
      else
         (*targetdata)->absvar[i] = NULL;

      if( i < ncluster )
      {
         if( sourcedata->indicatorvar[i] != NULL)
            SCIP_CALL( SCIPtransformVar(scip, sourcedata->indicatorvar[i], &((*targetdata)->indicatorvar[i])) );
         else
            (*targetdata)->indicatorvar[i] = NULL;
      }
   }

   SCIP_CALL( SCIPtransformVar(scip, sourcedata->targetvar, &((*targetdata)->targetvar)) );

   return SCIP_OKAY;
}

/**
 * delete-callback method of scip
 */
static
SCIP_DECL_PROBDELORIG(probdelorigSpa)
{
   int c1;
   int c2;
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
   /* absolute value vars */
   for ( c1 = 0; c1 < ((*probdata)->ncluster) * ((*probdata)->ncluster); ++c1 )
   {
      if( (*probdata)->absvar[c1] != NULL )
         SCIP_CALL( SCIPreleaseVar( scip, &((*probdata)->absvar[c1])) );
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
         for ( c1 = 0; c1 < (*probdata)->ncluster; ++c1 )
         {
            for( c2 = 0; c2 < (*probdata)->ncluster; ++c2 )
            {
               if( (*probdata)->edgevars[i][j][c1][c2] != NULL )
                  SCIP_CALL( SCIPreleaseVar( scip, &((*probdata)->edgevars[i][j][c1][c2])) );
            }
            SCIPfreeMemoryArray( scip, &((*probdata)->edgevars[i][j][c1]) );
         }
         SCIPfreeMemoryArray( scip, &((*probdata)->edgevars[i][j]) );
      }
      SCIPfreeMemoryArray( scip, &((*probdata)->edgevars[i]) );
   }
   SCIP_CALL( SCIPreleaseVar( scip, &((*probdata)->targetvar)) );

   SCIPfreeMemoryArray( scip, &((*probdata)->edgevars) );
   SCIPfreeMemoryArray(scip, &((*probdata)->binvars));
   SCIPfreeMemoryArray(scip, &((*probdata)->indicatorvar));
   SCIPfreeMemoryArray(scip, &((*probdata)->absvar));

   freeMatrix((*probdata)->cmatrix, (*probdata)->nbins);

   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

static
SCIP_DECL_PROBDELORIG(probdeltransSpa)
{
   int c1;
   int c2;
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
   /* absolute value vars */
   for ( c1 = 0; c1 < ((*probdata)->ncluster) * ((*probdata)->ncluster); ++c1 )
   {
      if( (*probdata)->absvar[c1] != NULL )
         SCIP_CALL( SCIPreleaseVar( scip, &((*probdata)->absvar[c1])) );
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
         for ( c1 = 0; c1 < (*probdata)->ncluster; ++c1 )
         {
            for( c2 = 0; c2 < (*probdata)->ncluster; ++c2 )
            {
               if( (*probdata)->edgevars[i][j][c1][c2] != NULL )
                  SCIP_CALL( SCIPreleaseVar( scip, &((*probdata)->edgevars[i][j][c1][c2])) );
            }
            SCIPfreeMemoryArray( scip, &((*probdata)->edgevars[i][j][c1]) );
         }
         SCIPfreeMemoryArray( scip, &((*probdata)->edgevars[i][j]) );
      }
      SCIPfreeMemoryArray( scip, &((*probdata)->edgevars[i]) );
   }
   SCIP_CALL( SCIPreleaseVar( scip, &((*probdata)->targetvar)) );

   SCIPfreeMemoryArray( scip, &((*probdata)->edgevars) );
   SCIPfreeMemoryArray(scip, &((*probdata)->binvars));
   SCIPfreeMemoryArray(scip, &((*probdata)->indicatorvar));
   SCIPfreeMemoryArray(scip, &((*probdata)->absvar));

   freeMatrix((*probdata)->cmatrix, (*probdata)->nbins);

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
   int c2;
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
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->absvar), ncluster * ncluster) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->edgevars[i]), nbins) );
      for( j = 0; j < nbins; ++j )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->edgevars[i][j]), ncluster) );
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->edgevars[i][j][c1]), ncluster) );
            for( c2 = 0; c2 < ncluster; ++c2 )
            {
               if( sourcedata->edgevars[i][j][c1][c2] != NULL )
               {
                  SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->edgevars[i][j][c1][c2], &var) );
                  if( SCIPvarIsActive(var) )
                  {
                     SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->edgevars[i][j][c1][c2]), varmap, consmap, global, &success) );
                     assert(success);
                     assert((*targetdata)->edgevars[i][j][c1][c2] != NULL);
                     SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->edgevars[i][j][c1][c2]) );
                  }
                  else
                     (*targetdata)->edgevars[i][j][c1][c2] = NULL;
               }
               else
                  (*targetdata)->edgevars[i][j][c1][c2] = NULL;
            }
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

   for( i = 0; i < ncluster * ncluster; ++i )
   {
      if( sourcedata->absvar[i] != NULL )
      {
         SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->absvar[i], &var) );

         if( SCIPvarIsActive(var) )
         {
            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->absvar[i]), varmap, consmap, global, &success) );
            assert(success);
            assert((*targetdata)->absvar[i] != NULL);
            SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->absvar[i]) );
         }
         else
            (*targetdata)->absvar[i] = NULL;

      }
      else
         (*targetdata)->absvar[i] = NULL;

      if( i < ncluster )
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

   SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->targetvar, &var) );
   SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->targetvar), varmap, consmap, global, &success) );
   assert(success);
   assert((*targetdata)->targetvar != NULL);
   SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->targetvar) );

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
   int                   nedges,             /**< number of edges */
   int                   ncluster,           /**< number of possible cluster */
   SCIP_Real**           edges,              /**< array with start- and endpoints of the edges */
   SCIP_Real*            sd                  /**< array with the stationary distribution */
)
{
   SCIP_PROBDATA* probdata = NULL;
   int i;
   SCIP_Real epsC;
   char model_variant;

   assert(nbins > 0);  /* at least one node */
   assert(nedges >= 0); /* no negative number of edges */
   assert(ncluster > 0); /* at least one cluster */


   SCIP_CALL( SCIPcreateProbBasic(scip, name) );

   /* get the parameters for the coherence bound */
   SCIP_CALL( SCIPgetRealParam( scip, "coherence_bound", &epsC ) );

   SCIP_CALL( SCIPgetCharParam(scip, "model_variant", &model_variant) );
   /* Set up the problem */
   SCIP_CALL( SCIPallocMemory( scip, &probdata) );

   /* allocate memory for the matrizes and create them from the edge-arrays */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->cmatrix), nbins) );
   for ( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->cmatrix[i]), nbins) );
   }
   SCIP_CALL( createCMatrix( scip, probdata->cmatrix, edges, sd, nedges, nbins) );

   probdata->big_M = bigM(probdata->cmatrix, nbins);

   assert(epsC >= 0 && epsC <= 1.0);
   SCIPinfoMessage(scip, NULL, "Creating problem: %s \n", name);



   /* create variables */
   probdata->nbins=nbins;
   probdata->ncluster=ncluster;

   SCIP_CALL( createVariables(scip, probdata) );

   switch (model_variant)
   {
   case 'e':
      SCIP_CALL( createProbEdgeRep(scip, probdata, epsC) );
      break;
   case 'r':
      SCIP_CALL( createProbEdgeRepRelax(scip, probdata, epsC) );
      break;
   case 'b':
      SCIP_CALL( createProbBinRep(scip, probdata, epsC) );
      break;
   case 'c':
      SCIP_CALL( createProbMaxkCut(scip, probdata, epsC) );
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

SCIP_VAR** SCIPspaGetAbsvars(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_PROBDATA* probdata;
   assert(scip!= NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   assert(probdata->absvar != NULL);

   return probdata->absvar;
}

SCIP_VAR* SCIPspaGetTargetvar(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_PROBDATA* probdata;
   assert(scip!= NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   assert(probdata->targetvar != NULL);

   return probdata->targetvar;
}

SCIP_VAR***** SCIPspaGetEdgevars(
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

/**
 * print the model-values like coherence in the clusters and transition-probabilities between clusters that are not evident from the scip-solution
 */
SCIP_RETCODE SCIPspaPrintSolutionValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< The solution containg the values */
)
{
   SCIP_PROBDATA* probdata;
   SCIP_Real value;
   int c1;
   int c2;
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
            value+= probdata->cmatrix[i][j] * SCIPgetSolVal(scip, sol, probdata->binvars[i][c1]) * SCIPgetSolVal(scip, sol, probdata->binvars[j][c1]);
         }
      }
      SCIPinfoMessage(scip, NULL, "Coherence in cluster %d is %f \n", c1 + 1, value);
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
         SCIPinfoMessage(scip, NULL, "q_%d_%d is %f \n", c1 + 1, c2 + 1, value);
      }
   }
   return SCIP_OKAY;
}

