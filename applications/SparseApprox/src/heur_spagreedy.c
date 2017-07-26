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
#include "scip/cons_and.h"

#define HEUR_NAME             "spagreedy"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         'h'
#define HEUR_PRIORITY         536870911
#define HEUR_FREQ             1
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

/**  calculate the current epsI-value for a q-matrix */
static
SCIP_Real getObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           qmatrix,            /**< The irreversibility matrix*/
   SCIP_Real             scale,
   int                   ncluster            /**< The number of cluster*/
)
{
   SCIP_Real objective = 0;
   int c,c2;

   for( c = 0; c < ncluster; ++c )
   {
      c2 = ( c + 1 ) % ncluster;
      objective += ( qmatrix[c][c2] - qmatrix[c2][c]);
      objective += scale * qmatrix[c][c];
   }
   /* if we have no transitions at all then irreversibility should be set to 0 */
   return objective;
}


/** Initialize the q-matrix from a given (possibly incomplete) clusterassignment */
static
void computeIrrevMat(
   SCIP_Real**           clusterassignment,  /**< The matrix containing the (incomplete) clusterassignment */
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
               if( clusterassignment[i][k] < 1 || clusterassignment[j][l] < 1 )
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
   SCIP_Real**           clusterassignment,  /**< The matrix containing the (incomplete) clusterassignment */
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
         if( clusterassignment[bin][cluster] == 1 )
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



/** Get the temporary objective value bound after newbin would be added to newcluster but dont not change anything with the clustering */
static
SCIP_Real getTempObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           qmatrix,            /**< The irreversibility matrix */
   SCIP_Real**           cmatrix,            /**< The transition matrix */
   SCIP_Real**           clusterassignment,  /**< The clusterassignment */
   int                   newbin,             /**< The bin that would be added to cluster */
   int                   newcluster,         /**< The cluster the bin would be added to */
   int                   nbins,              /**< The number of bins */
   int                   ncluster            /**< The number of cluster */
)
{
   SCIP_Real obj;
   SCIP_Real temp;
   int i;
   obj = getObjective(scip, qmatrix, SCIPspaGetScale(scip), ncluster);
   /* the coh in cluster changes as well as the flow to the next and the previous cluster */
   for( i = 0; i < nbins; ++i )
   {
      temp = (clusterassignment[i][phiinv(newcluster, ncluster)] < 1 ? 0 : 1);
      obj += (cmatrix[i][newbin] - cmatrix[newbin][i]) * temp;
      temp = (clusterassignment[i][phi(newcluster, ncluster)] < 1 ? 0 : 1);
      obj -= (cmatrix[i][newbin] - cmatrix[newbin][i]) * temp;
      temp = (clusterassignment[i][newcluster] < 1 ? 0 : 1);
      obj += (cmatrix[i][newbin] + cmatrix[newbin][i]) * temp;
   }

   return obj;
}

/* find and assign the next unassigned bin to an appropriate cluster */
static
SCIP_RETCODE assignNextBin(
   SCIP*                 scip,
   SCIP_Bool             localheur,
   SCIP_Real**           clusterassignment,  /**< The matrix with the Clusterassignment */
   SCIP_Real**           cmatrix,            /**< The transition matrix */
   SCIP_Real**           qmatrix,            /**< The irreversibility matrix */
   SCIP_Bool*            isassigned,         /**< TRUE, if the bin i was already assigned to a cluster*/
   int                   nbins,              /**< The number of bins*/
   int                   ncluster,           /**< The number of cluster*/
   int*                  amountassigned,     /**< The total amount of bins already assigned*/
   int*                  binsincluster,      /**< The number of bins currently in a cluster*/
   SCIP_Real*            epsI                /**< The pairwise irreversibility bound*/
)
{
   int i;
   int c;
   int c1;
   int save = -1;
   SCIP_Real* irrevbound;        /* gives the best achievable epsI-bound for each bin */
   SCIP_Bool** clusterispossible;/* Saves if a bin can be assigned to a cluster without violating feasibility */
   int* bestcluster;             /* saves the index of the cluster for which the best result was achieved */
   SCIP_Real tempirrev;
   /* allocate memory */
   SCIPallocClearMemoryArray(scip, &irrevbound, nbins);
   SCIPallocClearMemoryArray(scip, &bestcluster, nbins);
   SCIPallocClearMemoryArray(scip, &clusterispossible, nbins);
   for( i = 0; i < nbins; ++i )
   {
      SCIPallocClearMemoryArray(scip, &clusterispossible[i], ncluster);
   }

   /* make ceratin that each cluster is non-empty*/
   for( c = 0; c < ncluster; ++c )
   {
      tempirrev = 0;
      if( binsincluster[c] == 0 )
      {
         for( i = 0; i < nbins; ++i )
         {
            if( isassigned[i] )
               continue;
            irrevbound[i] = getTempObj(scip, qmatrix, cmatrix, clusterassignment, i, c, nbins, ncluster);
            if( irrevbound[i] > tempirrev )
            {
               save = i;
               tempirrev = irrevbound[i];
            }
            if( save == -1 )
               save = i;
         }
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            clusterassignment[save][c1] = 0;
         }
         clusterassignment[save][c] = 1;
         binsincluster[c]++;
         assert(binsincluster[c] == 1);
         isassigned[save] = TRUE;
         *amountassigned += 1;
         updateIrrevMat(clusterassignment, qmatrix, cmatrix, save, c, nbins, ncluster);
      }
   }

   for( i = 0; i < nbins; ++i )
   {
      int c2;
      bestcluster[i] = 0;
      irrevbound[i] = -SCIPinfinity(scip);
      if( isassigned[i] )
         continue;

      /* check which clusters the bin can be assigned to */
      for( c1 = 0; c1 < ncluster; ++c1 )
      {
         if( 0 != clusterassignment[i][c1] )
            clusterispossible[i][c1] = TRUE;
         else
            clusterispossible[i][c1] = FALSE;
         /* if assignment to i would violate abs-var assignment then set clusterpossible to FALSE */

      }
      /* calculate the irrevbound for all possible clusterassignments */
      for( c2 = 0; c2 < ncluster; ++c2 )
      {
         save = clusterassignment[i][c2];
         if( !clusterispossible[i][c2] || clusterassignment[i][c2] == 0 )
            continue;
         clusterassignment[i][c2] = 1;

         /* save the best possible irrevbound for each bin */
         tempirrev = getTempObj(scip, qmatrix, cmatrix, clusterassignment, i, c2, nbins, ncluster);
         if( SCIPisGT(scip, tempirrev, irrevbound[i]) )
         {
            irrevbound[i] = tempirrev;
            bestcluster[i] = c2;
         }
         clusterassignment[i][c2] = save;

      }
      if( localheur && SCIPisGT(scip, irrevbound[i], *epsI) )
      {
         break;
      }
   }
   {
      SCIP_Real max = -SCIPinfinity(scip);
      int ind = -1;
      /* take the bin with the highest increase in irrev-bound */
      for( i = 0; i < nbins; ++i )
      {
         if( SCIPisLT(scip, max, irrevbound[i]) )
         {
            max = irrevbound[i];
            ind = i;
         }
      }

      assert(!isassigned[ind] && ind > -1 && ind < nbins);
      /* assign this bin to the found cluster */
      isassigned[ind] = TRUE;
      for( c1 = 0; c1 < ncluster; ++c1 )
      {
         clusterassignment[ind][c1] = 0;
      }
      clusterassignment[ind][bestcluster[ind]] = 1;
      binsincluster[bestcluster[ind]]++;
      *amountassigned += 1;
      /* update the Irreversibility matrix */

      updateIrrevMat(clusterassignment, qmatrix, cmatrix, ind, bestcluster[ind], nbins, ncluster);

      *epsI = getObjective(scip, qmatrix, SCIPspaGetScale(scip), ncluster);
   }

   /* free the allocated memory */
   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeMemoryArray(scip, &clusterispossible[i]);
   }
   SCIPfreeMemoryArray(scip, &clusterispossible);
   SCIPfreeMemoryArray(scip, &irrevbound);
   SCIPfreeMemoryArray(scip, &bestcluster);
   return SCIP_OKAY;
}


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
   int nbins;
   int ncluster;
   SCIP_Real** clustering;   /* matrix for the assignment of the binary variables */
   int* binsincluster;        /* amount of bins in a given cluster */
   SCIP_Bool* isassigned;     /* TRUE if a bin has already bin assigned to a cluster */
   int i;
   int j;
   int c;
   int amountassigned;    /* total amount of bins assigned */
   SCIP_SOL* sol;
   SCIP_Bool possible = TRUE;
   SCIP_Bool feasible = FALSE;
   SCIP_Real epsI = 0.0;
   SCIP_HEURDATA* heurdata;

   *result = SCIP_DIDNOTRUN;
   amountassigned = 0;

   /* for now: do not use heurisitc if weighted objective is used */

   heurdata = SCIPheurGetData(heur);
   if( SCIPgetEffectiveRootDepth(scip) == heurdata->lasteffectrootdepth )
      return SCIP_OKAY;

   heurdata->lasteffectrootdepth = SCIPgetEffectiveRootDepth(scip);
   /* get the problem data from scip */
   cmatrix = SCIPspaGetCmatrix(scip);
   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);
   binvars = SCIPspaGetBinvars(scip);
   assert( nbins > 0 && ncluster > 0 );


   /* allocate memory for the assignment */
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &clustering, nbins) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &binsincluster, ncluster) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &qmatrix, ncluster) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &isassigned, nbins) );

   for ( i = 0; i < nbins; ++i )
   {
      if( i < ncluster )
         SCIP_CALL( SCIPallocClearMemoryArray(scip, &qmatrix[i], ncluster) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &clustering[i], ncluster) );
      for( j = 0; j < ncluster; ++j )
      {
         /* unassigned is set to -1 so we can differentiate between unassigned and fixed in the branch and bound tree */
         clustering[i][j] = -1;
      }
   }

   /* get the already fixed bin-variables from scip. An assignment of -1 one means unassigned. 0 is fixed to 0, 1 is fixed to 1 */
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < ncluster; ++j )
      {
         if( NULL == binvars[i][j] )
         {
            possible = FALSE;
            break;
         }
         /* if the bounds determine a fixed binary variable, then fix the variable in the clusterassignment */
         if( SCIPisEQ(scip, SCIPvarGetLbGlobal(binvars[i][j]), SCIPvarGetUbGlobal(binvars[i][j])) )
         {
            clustering[i][j] = SCIPvarGetLbGlobal(binvars[i][j]);
            if( SCIPisEQ(scip, 1.0, clustering[i][j]) )
            {
               binsincluster[j]++;
               isassigned[i] = TRUE;
               amountassigned += 1;
               for( c = 0; c < ncluster; ++c )
               {
                  if( clustering[i][c] == -1 )
                     clustering[i][c] = 0;
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
         if( 0 == clustering[i][j] )
            amountzeros++;
         if( 1 == clustering[i][j] )
            sum++;
      }
      if( ncluster == amountzeros || sum > 1 )
         possible = FALSE;
   }
   if( amountassigned < nbins && possible )
   {
      /* initialize the qmatrix and the lower irreversibility bound */
      computeIrrevMat(clustering, qmatrix, cmatrix, nbins, ncluster);
      epsI = getObjective(scip, qmatrix, SCIPspaGetScale(scip), ncluster);
      /* assign bins iteratively until all bins are assigned */
      while( amountassigned < nbins )
      {
         SCIP_CALL( assignNextBin(scip, heurdata->local, clustering, cmatrix, qmatrix, isassigned, nbins, ncluster, &amountassigned, binsincluster, &epsI ) );
      }
      /* assert that the assignment is valid in the sense that it is a partition of the bins. Feasibility is not checked in this method */
      assert(isPartition(scip,clustering, nbins, ncluster));
      /* update the qmatrix */
      computeIrrevMat(clustering, qmatrix, cmatrix, nbins, ncluster);

      /* if we use the model without absolute values, transform the found solution if necessary */
      /*if( model == 'w' || model == 's')
         rotateSolution(clustering, binsincluster, qmatrix, cmatrix, nbins, ncluster);*/

      /* set the variables the problem to the found clustering and test feasibility */
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
      SCIP_CALL( assignVars( scip, sol, clustering, nbins, ncluster) );
      SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, TRUE, TRUE, TRUE, TRUE, &feasible) );
      SCIPinfoMessage(scip, NULL, "\n" );
   }
   if( feasible )
      *result = SCIP_FOUNDSOL;
   else
      *result = SCIP_DIDNOTFIND;

   /* free allocated memory */
   for ( i = 0; i < nbins; ++i )
   {
      if( i < ncluster )
         SCIPfreeMemoryArray(scip, &qmatrix[i]);
      SCIPfreeMemoryArray(scip, &clustering[i]);
   }

   SCIPfreeMemoryArray(scip, &clustering);
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

   SCIP_CALL( SCIPaddBoolParam(scip, "localheur", "If set to true, heuristic assigns bins as soon as any improvement is found", &heurdata->local, FALSE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}
