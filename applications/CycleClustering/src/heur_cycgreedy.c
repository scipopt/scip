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

/**@file   heur_cycgreedy.c
 * @brief  Greedy primal heuristic. States are assigned to clusters iteratively. At each iteration all possible
 * assignments are computed and the one with the best change in objective value is selected.
 * @author Leon Eifler
 */

#include "heur_cycgreedy.h"

#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "scip/misc.h"
#include "probdata_cyc.h"
#include "scip/cons_and.h"

#define HEUR_NAME             "cycgreedy"
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
   int                   lasteffectrootdepth;/**< index of the last solution for which oneopt was performed */
   SCIP_Bool             local;              /**< the heuristic only computes assignments until any improvement is found */
};

/**  calculate the current objective value for a q-matrix */
static
SCIP_Real getObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           qmatrix,            /**< the irreversibility matrix*/
   SCIP_Real             scale,              /**< the scaling parameter in the objective function */
   int                   ncluster            /**< the number of cluster*/
   )
{
   SCIP_Real objective = 0.0;
   int c;
   int c2;

   for( c = 0; c < ncluster; ++c )
   {
      c2 = ( c + 1 ) % ncluster;
      objective += qmatrix[c][c2] - qmatrix[c2][c];
      objective += scale * qmatrix[c][c];
   }

   /* if we have no transitions at all then irreversibility should be set to 0 */
   return objective;
}

/** initialize the q-matrix from a given (possibly incomplete) clusterassignment */
static
void computeIrrevMat(
   SCIP_Real**           clusterassignment,  /**< the matrix containing the (incomplete) clusterassignment */
   SCIP_Real**           qmatrix,            /**< the returned matrix with the irreversibility between two clusters */
   SCIP_Real**           cmatrix,            /**< the transition-matrix containg the probability-data */
   int                   nbins,              /**< the number of bins */
   int                   ncluster            /**< the number of possible clusters */
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
               /* as -1 and 0 are both interpreted as 0, this check is necessary. Compute x_ik*x_jl*c_ij */
               if( clusterassignment[i][k] < 1 || clusterassignment[j][l] < 1 )
                  continue;

               qmatrix[k][l] += cmatrix[i][j];
            }
         }
      }
   }
}

/** update the irreversibility matrix, after the clusterassignment[newcluster][newbin] was either set
 *  from 0 to 1 or from 1 to 0
 */
static
void updateIrrevMat(
   SCIP_Real**           clusterassignment,  /**< the matrix containing the (incomplete) clusterassignment */
   SCIP_Real**           qmatrix,            /**< the returned matrix with the irreversibility between two clusters */
   SCIP_Real**           cmatrix,            /**< the transition-matrix containg the probability-data */
   int                   newbin,             /**< the bin to be added to the assignment */
   int                   newcluster,         /**< the bluster in which the bin was changed */
   int                   nbins,              /**< the number of bins */
   int                   ncluster            /**< the number of clusters */
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

/** get the temporary objective value bound after newbin would be added to newcluster
 *  but dont not change anything with the clustering
 */
static
SCIP_Real getTempObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           qmatrix,            /**< the irreversibility matrix */
   SCIP_Real**           cmatrix,            /**< the transition matrix */
   SCIP_Real**           clusterassignment,  /**< the clusterassignment */
   int                   newbin,             /**< the bin that would be added to cluster */
   int                   newcluster,         /**< the cluster the bin would be added to */
   int                   nbins,              /**< the number of bins */
   int                   ncluster            /**< the number of cluster */
   )
{
   SCIP_Real obj;
   SCIP_Real temp;
   int i;

   obj = getObjective(scip, qmatrix, SCIPcycGetScale(scip), ncluster);

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
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             localheur,          /**< should the heuristic only compute local optimal assignment */
   SCIP_Real**           clusterassignment,  /**< the matrix with the Clusterassignment */
   SCIP_Real**           cmatrix,            /**< the transition matrix */
   SCIP_Real**           qmatrix,            /**< the irreversibility matrix */
   SCIP_Bool*            isassigned,         /**< TRUE, if the bin i was already assigned to a cluster*/
   int                   nbins,              /**< the number of bins*/
   int                   ncluster,           /**< the number of cluster*/
   int*                  amountassigned,     /**< the total amount of bins already assigned*/
   int*                  binsincluster,      /**< the number of bins currently in a cluster*/
   SCIP_Real*            objective           /**< the objective */
   )
{
   SCIP_Real* binobjective;
   SCIP_Bool** clusterispossible;
   int* bestcluster;
   SCIP_Real tempobj;
   SCIP_Real max = -SCIPinfinity(scip);
   int i;
   int c;
   int c1;
   int c2;
   int save = -1;
   int ind = -1;

   /* allocate memory */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &binobjective, nbins) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &bestcluster, nbins) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &clusterispossible, nbins) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocClearBufferArray(scip, &clusterispossible[i], ncluster) ); /*lint !e866*/
   }

   /* make ceratin that each cluster is non-empty*/
   for( c = 0; c < ncluster; ++c )
   {
      tempobj = 0;

      if( binsincluster[c] == 0 )
      {
         for( i = 0; i < nbins; ++i )
         {
            /* if already assigned do nothing */
            if( isassigned[i] )
               continue;

            /* check if assigning this state is better than the previous best state */
            binobjective[i] = getTempObj(scip, qmatrix, cmatrix, clusterassignment, i, c, nbins, ncluster);

            if( binobjective[i] > tempobj )
            {
               save = i;
               tempobj = binobjective[i];
            }

            /* ensure that a state is assigned */
            if( save == -1 )
               save = i;
         }

         /* assign the found state to the cluster */
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            clusterassignment[save][c1] = 0;
         }

         clusterassignment[save][c] = 1;
         binsincluster[c]++;

         assert(binsincluster[c] == 1);

         isassigned[save] = TRUE;
         *amountassigned += 1;

         /* update the q-matrix */
         updateIrrevMat(clusterassignment, qmatrix, cmatrix, save, c, nbins, ncluster);
      }
   }

   /*phase 2: iteratively assign states such that at each iteration the highest objective improvement is achieved */
   for( i = 0; i < nbins; ++i )
   {
      bestcluster[i] = 0;
      binobjective[i] = -SCIPinfinity(scip);
   }

   for( i = 0; i < nbins; ++i )
   {
      if( isassigned[i] )
         continue;

      /* check which clusters the bin can be assigned to. -1 means unassigned, 0 means fixed to 0. */
      for( c1 = 0; c1 < ncluster; ++c1 )
      {
         /* if assignment to i would violate abs-var assignment then set clusterpossible to FALSE */
         if( 0 != clusterassignment[i][c1] )
            clusterispossible[i][c1] = TRUE;
         else
            clusterispossible[i][c1] = FALSE;
      }

      /* calculate the irrevbound for all possible clusterassignments */
      for( c2 = 0; c2 < ncluster; ++c2 )
      {
         if( !clusterispossible[i][c2] || clusterassignment[i][c2] == 0 )
            continue;

         /* temporarily assign i to c2 */
         save = (int) clusterassignment[i][c2];
         clusterassignment[i][c2] = 1;

         /* save the best possible irrevbound for each bin */
         tempobj = getTempObj(scip, qmatrix, cmatrix, clusterassignment, i, c2, nbins, ncluster);

         /* check if this is an improvement compared to the best known assignment */
         if( SCIPisGT(scip, tempobj, binobjective[i]) )
         {
            binobjective[i] = tempobj;
            bestcluster[i] = c2;
         }

         clusterassignment[i][c2] = save;
      }

      /* if localheur is true, then the heuristic assigns a state as soon as any improvement is found */
      if( localheur && SCIPisGT(scip, binobjective[i], *objective) )
         break;
   }

   /* take the bin with the highest increase in irrev-bound */
   for( i = 0; i < nbins; ++i )
   {
      if( SCIPisLT(scip, max, binobjective[i]) )
      {
         max = binobjective[i];
         ind = i;
      }
   }

   assert(!isassigned[ind] && ind > -1 && ind < nbins);

   /* assign this bin to the found cluster */
   for( c1 = 0; c1 < ncluster; ++c1 )
   {
      clusterassignment[ind][c1] = 0;
   }

   clusterassignment[ind][bestcluster[ind]] = 1;
   binsincluster[bestcluster[ind]]++;
   *amountassigned += 1;
   isassigned[ind] = TRUE;

   /* update the Irreversibility matrix */
   updateIrrevMat(clusterassignment, qmatrix, cmatrix, ind, bestcluster[ind], nbins, ncluster);
   *objective = getObjective(scip, qmatrix, SCIPcycGetScale(scip), ncluster);

   /* free the allocated memory */
   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeBufferArray(scip, &(clusterispossible[i]));
   }
   SCIPfreeBufferArray(scip, &clusterispossible);
   SCIPfreeBufferArray(scip, &bestcluster);
   SCIPfreeBufferArray(scip, &binobjective);

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyCycGreedy)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurCycGreedy(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeCycGreedy)
{  /*lint --e{715}*/
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
SCIP_DECL_HEUREXITSOL(heurExitsolCycGreedy)
{  /*lint --e{715}*/
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* reset the timing mask to its default value */
   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitCycGreedy)
{  /*lint --e{715}*/
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
SCIP_DECL_HEUREXEC(heurExecCycGreedy)
{  /*lint --e{715}*/
   SCIP_Real** cmatrix;                      /* the transition matrixx */
   SCIP_Real** qmatrix;                      /* the low-dimensional transition matrix between clusters */
   SCIP_VAR*** binvars;                      /* SCIP variables */
   SCIP_Real** clustering;                   /* matrix for the assignment of the binary variables */
   int* binsincluster;                       /* amount of bins in a given cluster */
   SCIP_Bool* isassigned;                    /* TRUE if a bin has already bin assigned to a cluster */
   SCIP_HEURDATA* heurdata;                  /* the heurdata */
   SCIP_SOL* sol;                            /* pointer to solution */
   SCIP_Bool possible = TRUE;                /* can the heuristic be run */
   SCIP_Bool feasible = FALSE;               /* is the solution feasible */
   SCIP_Real obj = 0.0;                      /* objective value */
   int amountassigned;                       /* total amount of bins assigned */
   int nbins;                                /* number of bins */
   int ncluster;                             /* number of cluster */
   int i;                                    /* running indices */
   int j;
   int c;

   *result = SCIP_DIDNOTRUN;
   amountassigned = 0;

   /* for now: do not use heurisitc if weighted objective is used */
   heurdata = SCIPheurGetData(heur);
   if( SCIPgetEffectiveRootDepth(scip) == heurdata->lasteffectrootdepth )
      return SCIP_OKAY;

   heurdata->lasteffectrootdepth = SCIPgetEffectiveRootDepth(scip);

   /* get the problem data from scip */
   cmatrix = SCIPcycGetCmatrix(scip);
   nbins = SCIPcycGetNBins(scip);
   ncluster = SCIPcycGetNCluster(scip);
   binvars = SCIPcycGetBinvars(scip);

   assert(nbins > 0 && ncluster > 0);

   /* allocate memory for the assignment */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &clustering, nbins) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &binsincluster, ncluster) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &qmatrix, ncluster) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &isassigned, nbins) );

   for ( i = 0; i < nbins; ++i )
   {
      if( i < ncluster )
      {
         SCIP_CALL( SCIPallocClearBufferArray(scip, &qmatrix[i], ncluster) ); /*lint !e866*/
      }

      SCIP_CALL( SCIPallocClearBufferArray(scip, &clustering[i], ncluster) ); /*lint !e866*/

      for( j = 0; j < ncluster; ++j )
      {
         /* unassigned is set to -1 so we can differentiate unassigned and fixed in the branch and bound tree */
         clustering[i][j] = -1;
      }
   }

   /* get the already fixed bin-variables from scip. An assignment of -1 one means unassigned.
    * 0 is fixed to 0, 1 is fixed to 1
    */
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
      obj = getObjective(scip, qmatrix, SCIPcycGetScale(scip), ncluster);

      /* assign bins iteratively until all bins are assigned */
      while( amountassigned < nbins )
      {
         SCIP_CALL( assignNextBin(scip, heurdata->local, clustering, cmatrix, qmatrix,
            isassigned, nbins, ncluster, &amountassigned, binsincluster, &obj ) );
      }

      /* assert that the assignment is valid in the sense that it is a partition of the bins.
       * Feasibility is not checked in this method
       */
      assert(isPartition(scip,clustering, nbins, ncluster));

      /* update the qmatrix */
      computeIrrevMat(clustering, qmatrix, cmatrix, nbins, ncluster);

      /* set the variables the problem to the found clustering and test feasibility */
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
      SCIP_CALL( assignVars( scip, sol, clustering, nbins, ncluster) );
      SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, TRUE, TRUE, TRUE, TRUE, &feasible) );
   }

   if( feasible )
      *result = SCIP_FOUNDSOL;
   else
      *result = SCIP_DIDNOTFIND;

   /* free allocated memory */
   for ( i = 0; i < nbins; ++i )
   {
      SCIPfreeBufferArray(scip, &clustering[i]);

      if( i < ncluster )
         SCIPfreeBufferArray(scip, &qmatrix[i]);
   }

   SCIPfreeBufferArray(scip, &isassigned);
   SCIPfreeBufferArray(scip, &qmatrix);
   SCIPfreeBufferArray(scip, &binsincluster);
   SCIPfreeBufferArray(scip, &clustering);

   return SCIP_OKAY;
}

/*
 * * primal heuristic specific interface methods
 */

/** creates the CycGreedy - primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurCycGreedy(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create greedy primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */

   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
      HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
      HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecCycGreedy, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyCycGreedy) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeCycGreedy) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolCycGreedy) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitCycGreedy) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "localheur", "If set to true, heuristic assigns bins as soon as any improvement is found",
      &heurdata->local, FALSE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}
