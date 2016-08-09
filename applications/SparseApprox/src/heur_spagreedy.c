/*
 * heur_spagreedy.c
 *
 *  Created on: Dec 7, 2015
 *      Author: bzfeifle
 */


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

/**@file   heur_spagreedy.c
 * @brief  Sparse Approximation primal heuristic
 * @author Leon Eifler
 */

#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "scip/misc.h"
#include "probdata_spa.h"
#include "heur_spagreedy.h"
#include "scip/cons_and.c"

#define HEUR_NAME             "spagreedy"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         'h'
#define HEUR_PRIORITY         536870911
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE                         /**< does the heuristic use a secondary SCIP instance? */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                  lasteffectrootdepth;/**< index of the last solution for which oneopt was performed */
   SCIP_Bool            local;
};


#ifndef NDEBUG
/** checks if the assignment is finished, i.e. all columns have exactly one 1 and rest 0 values */
static
SCIP_Bool isPartition(
   int**                 clusterassignment,  /**< The matrix containing the clustering */
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
         if( clusterassignment[c][i] == -1 )
            validassignment = FALSE;
         colsum += clusterassignment[c][i];
      }
      if( colsum != 1 )
         validassignment = FALSE;
   }
   return validassignment;
}
#endif

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
         temp = REALABS( qmatrix[i][j] - qmatrix[j][i] );
         if( SCIPisPositive(scip, temp) && SCIPisGT(scip, epsI, temp) )
            epsI = temp;
      }
   }
   /* if we have no transitions at all then irreversibility should be set to 0 */
   if( epsI == SCIP_INVALID )
   {
      epsI = 0.0;
   }
   return epsI;
}


/** Initialize the q-matrix from a given (possibly incomplete) clusterassignment */
static
void computeIrrevMat(
   int**                 clusterassignment,  /**< The matrix containing the (incomplete) clusterassignment */
   SCIP_Real**           qmatrix,            /**< The matrix with the return-values, in each cell is the irreversibility between two clusters */
   SCIP_Real**           cmatrix,            /**< The transition-matrix containg the probability-data */
   int                   nbins,              /**< The number of bins */
   int                   ncluster            /**< The number of possible clusters */
)
{
   int i;
   int j;
   int k;
   int l;
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
               if( clusterassignment[k][i] < 1 || clusterassignment[l][j] < 1 )
               {
                  continue;
               }
               qmatrix[k][l] += cmatrix[i][j];
            }
         }
      }
   }
}

/** Update the irreversibility matrix, after the clusterassignment[newcluster][newbin] was either set from 0 to 1 or from 1 to 0 */
static
void updateIrrevMat(
   int**                 clusterassignment,  /**< The matrix containing the (incomplete) clusterassignment */
   SCIP_Real**           qmatrix,            /**< The matrix with the return-values, in each cell is the irreversibility between two clusters */
   SCIP_Real**           cmatrix,            /**< The transition-matrix containg the probability-data */
   int                   newbin,             /**< The bin to be added to the assignment */
   int                   newcluster,         /**< The bluster in which the bin was changed */
   int                   nbins,              /**< The number of bins */
   int                   ncluster            /**< The number of clusters */
)
{
   int bin;
   int cluster;

   for( cluster = 0; cluster < ncluster; ++cluster )
   {
      for( bin = 0; bin < nbins; ++bin )
      {
         /* multiplier is 1 if clusterassignment is 1, and 0 if it is 0 (set to 0) or -1 (unassigned) */
         int temp = 0;
         if( clusterassignment[cluster][bin] == 1 )
            temp = 1;
         if( cluster != newcluster )
         {
            qmatrix[newcluster][cluster] +=  temp * cmatrix[newbin][bin];
            qmatrix[cluster][newcluster] +=  temp * cmatrix[bin][newbin];
         }
         else
         {
            if( bin == newbin )
               qmatrix[newcluster][newcluster] += cmatrix[newbin][bin];
            else
               qmatrix[newcluster][newcluster] += (cmatrix[newbin][bin] + cmatrix[bin][newbin]) * temp;
         }
      }
   }
}

static
void rotateSolution(
   int**                 clusterassignment,
   int*                  binsincluster,
   SCIP_Real**           qmatrix,
   SCIP_Real**           cmatrix,
   int                   nbins,
   int                   ncluster
)
{
   int i;
   SCIP_Real temp;
   temp = binsincluster[0];
   binsincluster[0] = binsincluster[1];
   binsincluster[1] = temp;
   if( qmatrix[0][1] - qmatrix[1][0] < 0 )
   {
      for( i = 0; i < nbins; ++i )
      {
         temp = clusterassignment[0][i];
         clusterassignment[0][i] = clusterassignment[1][i];
         clusterassignment[1][i] = temp;
      }
   }
   computeIrrevMat(clusterassignment, qmatrix, cmatrix, nbins, ncluster);
}

/**  assign the variables in scip according to the found clusterassignment */
static
SCIP_RETCODE assignVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< The SCIP solution */
   int**                 clusterassignment,  /**< The matrix with the clusterassignment */
   int*                  binsincluster,      /**< The Array with the number of bins in each cluster */
   int                   nbins,              /**< The number of bins */
   int                   ncluster,           /**< The number of cluster */
   SCIP_Real**           qmatrix             /**< The irreversibility matrix */
)
{
   int i,j;
   int c;
   int c2;
   SCIP_VAR* var;
   SCIP_VAR** indvars;
   SCIP_VAR*** binvars;
   SCIP_Real** matrixtest;
   SCIP_VAR***** edgevars;

   assert(nbins > 0 && ncluster > 0);

   indvars = SCIPspaGetIndvars(scip);
   binvars = SCIPspaGetBinvars(scip);
   edgevars = SCIPspaGetEdgevars(scip);
   SCIPallocClearMemoryArray(scip, &matrixtest, ncluster);
   for( i = 0; i < ncluster; ++i )
   {
      SCIPallocClearMemoryArray(scip, &matrixtest[i], ncluster);
   }

   for ( c = 0; c < ncluster; ++c )
   {
      /* set indicatorvar whether cluster is nonempty */
      if( NULL != indvars[c] && binsincluster[c] > 0 && SCIPisEQ(scip, SCIPvarGetUbGlobal(indvars[c]), 1.0) )
         SCIP_CALL( SCIPsetSolVal(scip, sol, indvars[c], 1.0) );
      else
      {
         if( NULL != indvars[c] && SCIPisZero(scip, SCIPvarGetLbGlobal(indvars[c])) )
            SCIP_CALL( SCIPsetSolVal(scip, sol, indvars[c], 0.0) );
      }
      /* set values of binary variables */
      for ( i = 0; i < nbins; ++i )
      {
         /* check if the clusterassignment ist feasible for the variable bounds. If not do not assign the variable */
         if( NULL != binvars[i][c] )
         {
            if( SCIPvarIsTransformed(binvars[i][c]) )
               var = binvars[i][c];
            else
               var = SCIPvarGetTransVar(binvars[i][c] );
            if( NULL != var && SCIPisLE(scip, SCIPvarGetLbGlobal(var), clusterassignment[c][i]) && SCIPisGE(scip, SCIPvarGetUbGlobal(var), clusterassignment[c][i]) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
               SCIP_CALL( SCIPsetSolVal( scip, sol, binvars[i][c], clusterassignment[c][i]) );
            assert( SCIPisIntegral(scip, clusterassignment[c][i]) );
         }
      }

      /* set the value for the edgevariables */
      for( i = 0; i < nbins; ++i )
      {
         for( j = 0; j < nbins; ++j )
         {
            if( i == j )
               continue;
            if( j < i )
            {
               if( NULL == edgevars[i][j][c][c] )
                  continue;
               if( SCIPvarIsTransformed(edgevars[i][j][c][c]) )
                  var = edgevars[i][j][c][c];
               else
                  var = SCIPvarGetTransVar(edgevars[i][j][c][c]);
               if( NULL != var && SCIPisGE(scip, SCIPvarGetUbGlobal(var), clusterassignment[c][j] * clusterassignment[c][i]) && SCIPisLE(scip, SCIPvarGetLbGlobal(var), clusterassignment[c][j] * clusterassignment[c][i]) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
                  SCIP_CALL( SCIPsetSolVal( scip, sol, var, clusterassignment[c][j] * clusterassignment[c][i]  ) );
            }
            c2 = (c + 1) % 3;
            if( NULL == edgevars[i][j][c][c2] )
               continue;
            if( SCIPvarIsTransformed(edgevars[i][j][c][c2]) )
               var = edgevars[i][j][c][c2];
            else
               var = SCIPvarGetTransVar(edgevars[i][j][c][c2]);
            if( NULL != var && SCIPisGE(scip, SCIPvarGetUbGlobal(var), clusterassignment[c2][j] * clusterassignment[c][i]) && SCIPisLE(scip, SCIPvarGetLbGlobal(var), clusterassignment[c2][j] * clusterassignment[c][i]) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
               SCIP_CALL( SCIPsetSolVal( scip, sol, var, clusterassignment[c2][j] * clusterassignment[c][i]  ) );
         }
      }
   }
   /* retransform the solution to original space, as the solution may be infeasible in transformed space due to presolving */
   SCIPretransformSol(scip, sol);
   /* free the allocated memory */
   for( i = 0; i < ncluster; ++i )
   {
      SCIPfreeMemoryArray(scip, &matrixtest[i]);
   }
   SCIPfreeMemoryArray(scip, &matrixtest);
   return SCIP_OKAY;
}

/** just search for the bin-pair with highest local irreversibility and asssign them to the first 2 clusters */
static
void assignFirstPair(
   SCIP*                 scip,               /**< Scip data structure */
   int**                 clusterassignment,  /**< The matrix with the Clusterassignment */
   SCIP_Real**           cmatrix,            /**< The transition matrix */
   SCIP_Real**           qmatrix,            /**< The irreversibility matrix */
   SCIP_Bool*            isassigned,         /**< TRUE, if the bin i was already assigned to a cluster*/
   int                   nbins,              /**< The number of bins*/
   int                   ncluster,           /**< The number of clusters*/
   int                   amountassigned,     /**< The total amount of bins already assigned*/
   int*                  binsincluster,      /**< The number of bins currently in a cluster*/
   SCIP_Real*            epsI                /**< The pairwise irreversibility bound*/
)
{
   int i,j;
   int c1;
   int maxind1, maxind2;
   SCIP_Real max = 0;

   maxind1 = 0;
   maxind2 = 0;
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( REALABS(cmatrix[i][j]-cmatrix[j][i]) > max )
         {
            max = REALABS(cmatrix[i][j]-cmatrix[j][i]);
            maxind1 = j;
            maxind2 = i;
         }
      }
   }
   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      clusterassignment[c1][maxind1] = 0;
      clusterassignment[c1][maxind2] = 0;

   }
   clusterassignment[0][maxind1] = 1;
   clusterassignment[1][maxind2] = 1;

   binsincluster[0]++;
   binsincluster[1]++;
   isassigned[maxind1] = TRUE;
   isassigned[maxind2] = TRUE;

   computeIrrevMat(clusterassignment, qmatrix, cmatrix, nbins, ncluster);
   *epsI = getIrrevBound(scip, qmatrix, ncluster);
}

/** Get the temporary irreversibility bound after newbin would be added to newcluster but don not change anything with the clustering */
static
SCIP_Real getTempIrrevBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           qmatrix,            /**< The irreversibility matrix */
   SCIP_Real**           cmatrix,            /**< The transition matrix */
   int**                 clusterassignment,  /**< The clusterassignment */
   int                   newbin,             /**< The bin that would be added to cluster */
   int                   newcluster,         /**< The cluster the bin would be added to */
   int                   nbins,              /**< The number of bins */
   int                   ncluster            /**< The number of cluster */
)
{
   int bin;
   int cluster;
   SCIP_Real epsI = SCIP_INVALID;
   SCIP_Real temp;
   int i;
   int j;

   for( cluster = 0; cluster < ncluster; ++cluster )
   {
      SCIP_Real q1 = qmatrix[newcluster][cluster];
      SCIP_Real q2 = qmatrix[cluster][newcluster];
      for( bin = 0; bin < nbins; ++bin )
      {

         /* multiplier is 1 if clusterassignment is 1, and 0 if it is 0 (set to 0) or -1 (unassigned) */
         temp = 0;
         if( clusterassignment[cluster][bin] == 1 )
            temp = 1;
         if( cluster != newcluster )
         {
            q1 +=  temp * cmatrix[newbin][bin];
            q2 +=  temp * cmatrix[bin][newbin];
         }
         else
         {
            if( bin == newbin )
               q1 += cmatrix[newbin][bin];
            else
               q2 += (cmatrix[newbin][bin] + cmatrix[bin][newbin]) * temp;
         }
      }

      temp = REALABS(q1-q2);
      if( SCIPisPositive(scip, temp) && SCIPisLT(scip, temp, epsI) )
         epsI = temp;
   }
   for( i = 0; i < ncluster; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( i == newcluster || j == newcluster )
            continue;
         temp = REALABS( qmatrix[i][j] - qmatrix[j][i] );
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

/* find and assign the next unassigned bin to an appropriate cluster */
static
SCIP_RETCODE assignNextBin(
   SCIP*                 scip,
   SCIP_Bool             localheur,
   int**                 clusterassignment,  /**< The matrix with the Clusterassignment */
   SCIP_Real**           cmatrix,            /**< The transition matrix */
   SCIP_Real**           qmatrix,            /**< The irreversibility matrix */
   SCIP_Bool*            isassigned,         /**< TRUE, if the bin i was already assigned to a cluster*/
   int                   nbins,              /**< The number of bins*/
   int                   ncluster,           /**< The number of cluster*/
   int                   amountassigned,     /**< The total amount of bins already assigned*/
   int*                  binsincluster,      /**< The number of bins currently in a cluster*/
   SCIP_Real*            epsI                /**< The pairwise irreversibility bound*/
)
{
   int i;
   int save;
   SCIP_Real* irrevbound;        /* gives the best achievable epsI-bound for each bin */
   SCIP_Bool** clusterispossible;/* Saves if a bin can be assigned to a cluster without violating feasibility */
   int* bestcluster;             /* saves the index of the cluster for which the best result was achieved */
   SCIP_Real tempirrev;
   /* allocate memory */
   SCIPallocClearMemoryArray(scip, &irrevbound, nbins);
   SCIPallocClearMemoryArray(scip, &bestcluster, nbins);
   SCIPallocClearMemoryArray(scip, &clusterispossible, ncluster);
   for( i = 0; i < ncluster; ++i )
   {
      SCIPallocClearMemoryArray(scip, &clusterispossible[i], nbins);
   }

   for( i = 0; i < nbins; ++i )
   {
      int c1;
      int c2;
      bestcluster[i] = 0;
      irrevbound[i] = 0;
      if( isassigned[i] )
         continue;

      /* check which clusters the bin can be assigned to */
      for( c1 = 0; c1 < ncluster; ++c1 )
      {
         if( 0 != clusterassignment[c1][i] )
            clusterispossible[c1][i] = TRUE;
         else
            clusterispossible[c1][i] = FALSE;
         /* if assignment to i would violate abs-var assignment then set clusterpossible to FALSE */

      }
      /* calculate the irrevbound for all possible clusterassignments */
      for( c2 = 0; c2 < ncluster; ++c2 )
      {
         save = clusterassignment[c2][i];
         if( !clusterispossible[c2][i] || clusterassignment[c2][i] == 0 )
            continue;
         clusterassignment[c2][i] = 1;

         /* save the best possible irrevbound for each bin */
         tempirrev = getTempIrrevBound(scip, qmatrix, cmatrix, clusterassignment, i, c2, nbins, ncluster);
         if( SCIPisGT(scip, tempirrev, irrevbound[i]) )
         {
            irrevbound[i] = tempirrev;
            bestcluster[i] = c2;
         }
         clusterassignment[c2][i] = save;

      }
      if( localheur && SCIPisGT(scip, irrevbound[i], *epsI) )
      {
         break;
      }
   }
   {
      SCIP_Real max = *epsI;
      int ind = -1;
      int c1;
      /* take the bin with the highest increase in irrev-bound */
      for( i = 0; i < nbins; ++i )
      {
         if( SCIPisLT(scip, max, irrevbound[i]) )
         {
            max = irrevbound[i];
            ind = i;
         }
      }
      /* If we can not increase the eps-bound then we assign the bin to a possible cluster with minimal coherence.
       * Take into account whether it is possible to assign the bin to a cluster due to fixed vars from scip.  */
      if( -1 == ind )
      {
         SCIP_Bool binfound = FALSE;
         int row = 0;
         int col = 0;
         int mincluster = 0;
         int clusterfailed = 0;
         SCIP_Bool* tried;
         SCIP_Real mincoherence;
         SCIPallocClearMemoryArray(scip, &tried, ncluster);
         while( !binfound )
         {
            mincoherence = 4;
            /*find the open cluster with minimal coherence */
            if( SCIPisZero(scip, qmatrix[ncluster - 1][ncluster - 1]) )
               qmatrix[ncluster - 1][ncluster - 1] = 0.0000001;
            for( row = 0; row < ncluster; ++row )
            {
               SCIP_Real qtemp;
               if( tried[row] )
                  continue;
               if( SCIPisZero(scip, qmatrix[row][row]) )
                  qtemp = 1;
               else
                  qtemp = qmatrix[row][row];
               if( SCIPisLT(scip, qtemp, mincoherence) )
               {
                  mincoherence = qtemp;
                  mincluster = row;
                  tried[row] = TRUE;
               }
            }
            col = 0;
            /* try to assign a bin to the mincluster. This is possible if there is a bin for which this cluster is possible and the bin is not assigned yet */
            while( (!clusterispossible[mincluster][col] || isassigned[col]) && col < nbins )
            {
               ++col;
            }
            /* if we can not find a bin for this cluster then try to find another one */
            if( nbins == col )
               clusterfailed++;
            else
            {
               /* if we found a bin and a cluster, then set this assignment */
               binfound = TRUE;
               ind = col;
               bestcluster[ind] = mincluster;
               assert( ind < nbins );
            }
            /* if we tried each cluster once and none was possible then we will assign a locally infeasible bin. */
            if( clusterfailed == ncluster )
               break;
         }

         /*check if there is no possible locally valid assignment. Then just assign the first non-assigned bin to the first cluster */
         if( clusterfailed == ncluster )
         {
            row = 0;
            col = 0;
            while(!binfound)
            {
               if( row == ncluster )
                  break;
               if( clusterassignment[row][col] != -1 )
               {
                  ++col;
                  if( col == nbins )
                  {
                     col = 0;
                     row++;
                  }
               }
               else
               {
                  ind = col;
                  bestcluster[ind] = row;
                  binfound = TRUE;
               }
            }
         }
         SCIPfreeMemoryArray(scip, &tried);
      }
      assert(!isassigned[ind] && ind > -1 && ind < nbins);
      /* assign this bin to the found cluster */
      isassigned[ind] = TRUE;
      for( c1 = 0; c1 < ncluster; ++c1 )
      {
         clusterassignment[c1][ind] = 0;
      }
      clusterassignment[bestcluster[ind]][ind] = 1;
      binsincluster[bestcluster[ind]]++;
      /* update the Irreversibility matrix */

      updateIrrevMat(clusterassignment, qmatrix, cmatrix, ind, bestcluster[ind], nbins, ncluster);

      *epsI = getIrrevBound(scip, qmatrix, ncluster);
   }

   /* free the allocated memory */
   for( i = 0; i < ncluster; ++i )
   {
      SCIPfreeMemoryArray(scip, &clusterispossible[i]);
   }
   SCIPfreeMemoryArray(scip, &clusterispossible);
   SCIPfreeMemoryArray(scip, &irrevbound);
   SCIPfreeMemoryArray(scip, &bestcluster);
   return SCIP_OKAY;
}

#if 0
/** permutes the clusterassignment-matrix such that it is an upper-triangle matrix. Use this method if the orbitope-constraint is used in the progam */
static
SCIP_RETCODE upperTriangle(
   SCIP*                 scip,               /**< SCIP data structure*/
   int**                 clusterassignment,  /**< The matrix containing the clusterassignment*/
   int*                  binsincluster,      /**< The number of bins in each cluster*/
   int                   nbins,              /**< The number of bins*/
   int                   ncluster            /**< The number of cluster*/
)
{
   int i = 0;
   int c = 0;
   int emptycluster = 0;
   int** temp;
   int* indexsort;
   SCIP_Bool* sorted;
   /* allocate memory for a temporary matrix to copy the permuataion out of */
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &sorted, ncluster) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &indexsort, ncluster) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &temp, ncluster) );
   for( i = 0; i < ncluster; ++i )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &temp[i], nbins) );
      indexsort[i] = i;
   }

   for( i = 0; i < ncluster; ++i )
   {
      /* get all empty clusters. Those will be the bottom rows of the matrix */
      if( 0 == binsincluster[i] )
      {
         emptycluster++;
         indexsort[ncluster - emptycluster] = i;
         sorted[i] = TRUE;
      }
   }
   i = 0;
   while( i < nbins )
   {
      int j = 0;
      while(clusterassignment[j][i] == 0)
      {
         ++j;
         assert(j < ncluster);
      }
      if( !sorted[j] )
      {
         indexsort[c] = j;
         sorted[j] = TRUE;
         ++c;
      }
      while(clusterassignment[j][i] == 1)
      {
         ++i;
         if( i == nbins )
            break;
      }
   }

   for( i = 0; i < ncluster; ++i )
   {
      if( indexsort[i] != -1 )
      {
         for( c = 0; c < nbins; ++c )
         {
            temp[i][c] = clusterassignment[indexsort[i]][c];
         }
      }

   }
   for( i = 0; i < ncluster; ++i )
   {
      binsincluster[i] = 0;
      for( c = 0; c < nbins; ++c )
      {
         clusterassignment[i][c] = temp[i][c];
         binsincluster[i] += clusterassignment[i][c];
      }
   }

   for( i = 0; i < ncluster; ++i )
   {
      SCIPfreeMemoryArray(scip, &temp[i]);
   }
   SCIPfreeMemoryArray(scip, &sorted);
   SCIPfreeMemoryArray(scip, &indexsort);
   SCIPfreeMemoryArray(scip, &temp);
   return SCIP_OKAY;
}
#endif


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopySpaGreedy)
{
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurSpaGreedy(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSpaGreedy)
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolSpaGreedy)
{
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* reset the timing mask to its default value */
   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitSpaGreedy)
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize last solution index */

   heurdata->lasteffectrootdepth = -1;

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSpaGreedy)
{
   SCIP_Real** cmatrix;       /* The transition matrixx */
   SCIP_Real** qmatrix;       /* The low-dimensional transition matrix between clusters */
   SCIP_VAR*** binvars;       /* SCIP variables */
   SCIP_VAR***** edgevars;
   int nbins;
   int ncluster;
   int** clusterassignment;   /* matrix for the assignment of the binary variables */
   int* binsincluster;        /* amount of bins in a given cluster */
   SCIP_Bool* isassigned;     /* TRUE if a bin has already bin assigned to a cluster */
   int i;
   int j;
   int k;
   int amountassigned = 0;    /* total amount of bins assigned */
   SCIP_SOL* sol;
   SCIP_Bool possible = TRUE;
   SCIP_Bool feasible = FALSE;
   SCIP_Real epsI = 0.0;
   SCIP_HEURDATA* heurdata;
   char model;

   *result = SCIP_DIDNOTRUN;

   /* for now: do not use heurisitc if weighted objective is used */
   model = SCIPspaGetModel(scip);

   heurdata = SCIPheurGetData(heur);
   if( SCIPgetEffectiveRootDepth(scip) == heurdata->lasteffectrootdepth )
      return SCIP_OKAY;

   heurdata->lasteffectrootdepth = SCIPgetEffectiveRootDepth(scip);
   /* get the problem data from scip */
   cmatrix = SCIPspaGetCmatrix(scip);
   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);
   binvars = SCIPspaGetBinvars(scip);
   edgevars = SCIPspaGetEdgevars(scip);
   assert( nbins > 0 && ncluster > 0 );


   /* allocate memory for the assignment */
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &clusterassignment, ncluster) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &binsincluster, ncluster) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &qmatrix, ncluster) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &isassigned, nbins) );

   for ( i = 0; i < ncluster; ++i )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &qmatrix[i], ncluster) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &clusterassignment[i], nbins) );
      for( j = 0; j < nbins; ++j )
      {
         /* unassigned is set to -1 so we can differentiate between unassigned and fixed in the branch and bound tree */
         clusterassignment[i][j] = -1;
      }
   }

   /* get the already fixed bin-variables from scip. An assignment of -1 one means unassigned. 0 is fixed to 0, 1 is fixed to 1 */
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < ncluster; ++j )
      {
         /* if the variable is not active or got deleted in presolving, find out the fixation over the constraints */
         if( NULL == binvars[i][j] )
         {
            for( k = 0; k < ncluster; ++k )
            {
               if( k != j )
               {
                  int h = 0;
                  while((binvars[h][k] == NULL || (SCIPvarGetStatus(binvars[h][k]) == SCIP_VARSTATUS_FIXED) || h == i) && h < nbins - 1)
                     h++;
                  if( h < nbins )
                  {
                     if( NULL != edgevars[i][h][j][k] && SCIPvarIsActive(edgevars[i][h][j][k]) )
                     {
                        clusterassignment[j][i] = 1;
                        binsincluster[j]++;
                        isassigned[i] = TRUE;
                        amountassigned++;
                     }
                     else
                        clusterassignment[j][i] = 0;
                  }
               }
               if( clusterassignment[j][i] != -1 )
                  break;
            }
         } else
         {
            /* if the bounds determine a fixed binary variable, then fix the variable in the clusterassignment */
            if( SCIPisEQ(scip, SCIPvarGetLbGlobal(binvars[i][j]), SCIPvarGetUbGlobal(binvars[i][j])) )
            {
               clusterassignment[j][i] = SCIPvarGetLbGlobal(binvars[i][j]);
               if( SCIPisEQ(scip, 1.0, clusterassignment[j][i]) )
               {
                  binsincluster[j]++;
                  isassigned[i] = TRUE;
                  amountassigned++;
               }
            }
         }
      }
   }
   /* check if the assignment violates paritioning, e.g. because we are in a subscip */
   for( i = 0; i < nbins; ++i )
   {
      int amountzeros = 0;
      int sum = 0;
      for( j = 0; j < ncluster; ++j )
      {
         if( 0 == clusterassignment[j][i] )
            amountzeros++;
         if( 1 == clusterassignment[j][i] )
            sum++;
      }
      if( ncluster == amountzeros || sum > 1 )
         possible = FALSE;
   }
   if( amountassigned < nbins && possible )
   {
      /* initialize the qmatrix and the lower irreversibility bound */
      computeIrrevMat(clusterassignment, qmatrix, cmatrix, nbins, ncluster);
      epsI = getIrrevBound(scip, qmatrix, ncluster);
      /* if no bins are assigned, then choose the first two bins manually */
      if( 0 == amountassigned )
      {
         assignFirstPair(scip, clusterassignment, cmatrix, qmatrix, isassigned, nbins, ncluster, amountassigned, binsincluster, &epsI);
         amountassigned = 2;
      }
      /* assign bins iteratively until all bins are assigned */
      while( amountassigned < nbins )
      {
         SCIP_CALL( assignNextBin(scip, heurdata->local, clusterassignment, cmatrix, qmatrix, isassigned, nbins, ncluster, amountassigned, binsincluster, &epsI ) );
         amountassigned++;
      }
      /* assert that the assignment is valid in the sense that it is a partition of the bins. Feasibility is not checked in this method */
      assert(isPartition(clusterassignment, nbins, ncluster));
      /* update the qmatrix */
      computeIrrevMat(clusterassignment, qmatrix, cmatrix, nbins, ncluster);

      /* if we use the model without absolute values, transform the found solution if necessary */
      if( model == 'w' || model == 's')
         rotateSolution(clusterassignment, binsincluster, qmatrix, cmatrix, nbins, ncluster);

      /* set the variables the problem to the found clustering and test feasibility */
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
      SCIP_CALL( assignVars( scip, sol, clusterassignment, binsincluster, nbins, ncluster, qmatrix) );
      SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, FALSE, FALSE, &feasible) );
   }
   if( feasible )
      *result = SCIP_FOUNDSOL;
   else
      *result = SCIP_DIDNOTFIND;

   /* free allocated memory */
   for ( i = 0; i < ncluster; ++i )
   {
      SCIPfreeMemoryArray(scip, &qmatrix[i]);
      SCIPfreeMemoryArray(scip, &clusterassignment[i]);
   }

   SCIPfreeMemoryArray(scip, &clusterassignment);
   SCIPfreeMemoryArray(scip, &binsincluster);
   SCIPfreeMemoryArray(scip, &qmatrix);
   SCIPfreeMemory(scip, &isassigned);

   return SCIP_OKAY;
}

/*
 * * primal heuristic specific interface methods
 */

/** creates the SpaGreedy - primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSpaGreedy(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create xyz primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */

   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
      HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
      HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecSpaGreedy, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopySpaGreedy) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeSpaGreedy) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolSpaGreedy) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitSpaGreedy) );

   SCIP_CALL( SCIPaddBoolParam(scip, "localheur", "If set to true, heuristic assigns bins as soon as any improvement is found", &heurdata->local, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
