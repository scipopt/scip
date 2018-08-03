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

/**@file   heur_cyckerlin.c
 * @brief  improvement heuristic that exchanges binary variables between clusters.
 *  Similar to the famous kernighan/lin heuristic for graph partitioning
 * @author Leon Eifler
 */

/*---+---- 1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "heur_cyckerlin.h"

#include <assert.h>
#include <string.h>

#include "probdata_cyc.h"
#include "scip/pub_misc.h"

#define HEUR_NAME             "cyckerlin"
#define HEUR_DESC             "switch heuristic that tries to improve solution by trading states betweeen clusters"
#define HEUR_DISPCHAR         '@'
#define HEUR_PRIORITY         500
#define HEUR_FREQ             10
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define MAXPERMUTATIONS       5
#define DEFAULT_RANDSEED      177           /**< random seed */
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE          /**< does the heuristic use a secondary SCIP instance? */

struct SCIP_HeurData
{
   SCIP_SOL** candidates;
   int ncandidates;
   int candlength;
};

/*
 * Local methods
 */

/** external method that adds a solution to the list of candidate-solutions that should be improved */
SCIP_RETCODE addCandSolCyckerlin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< the given solution */
   )
{
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;

   heur = SCIPfindHeur(scip, "cyckerlin");

   if( heur == NULL )
      return SCIP_OKAY;

   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);
   assert(sol != NULL);

   /* realloc candidate array, if necessary */
   if( heurdata->candlength - 1 <= heurdata->ncandidates )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(heurdata->candidates), (SCIP_Longint) heurdata->candlength * 2) );
      heurdata->candlength *= 2;
   }

   heurdata->candidates[heurdata->ncandidates] = sol;
   heurdata->ncandidates++;

   return SCIP_OKAY;
}


/** get the bin-var assignment from scip and save it as a matrix */
static
SCIP_RETCODE getSolutionValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             bestsol,            /**< the solution */
   SCIP_Real**           solclustering,      /**< matrix to save the bin-vars*/
   SCIP_Bool**           binfixed,           /**< matrix to save if a bin is fixed in scip */
   int*                  clusterofbin,       /**< array containing the cluster of each bin */
   int*                  nbinsincluster      /**< number of bins in each cluster */
   )
{
   SCIP_VAR*** binvars;
   int nbins;
   int ncluster;
   int i;
   int c;

   assert(bestsol != NULL);

   nbins = SCIPcycGetNBins(scip);
   ncluster = SCIPcycGetNCluster(scip);
   binvars = SCIPcycGetBinvars(scip);

   assert(nbins > 0 && ncluster > 0 && nbins > ncluster);
   assert(binvars != NULL);

   /* get the bin-variable values from the solution */
   for( i = 0; i < nbins; ++i )
   {
      for( c = 0; c < ncluster; ++c )
      {
         binfixed[i][c] = FALSE;

         if( binvars[i][c] != NULL)
         {
            if( (SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(binvars[i][c]), SCIPvarGetLbGlobal(binvars[i][c]))) )
               binfixed[i][c] = TRUE;

            solclustering[i][c] = SCIPgetSolVal(scip, bestsol, binvars[i][c]);

            if( SCIPisFeasEQ(scip, solclustering[i][c], 1.0) )
            {
               clusterofbin[i] = c;
               nbinsincluster[c]++;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** Set a bin to a new cluster, update the qmatrix. */
static
void setBinToCluster(
   SCIP_Real**           solclustering,      /**< the matrix with the clustering of the bins */
   SCIP_Real**           cmatrix,            /**< the transition matrix*/
   SCIP_Real**           qmatrix,            /**< the matrix containing the transition probabilities between clusters*/
   int                   newbin,             /**< the bin to be changed*/
   int                   newcluster,         /**< the cluster where the bin is changed*/
   SCIP_Bool             setone,             /**< TRUE if the assignment is switched from 0 to 1, FALSE if 1 to 0*/
   int                   nbins,              /**< the number of bins*/
   int                   ncluster            /**< the number of clusters*/
   )
{
   int bin;
   int cluster;
   SCIP_Real sign = setone ? 1.0 : -1.0;

   for( cluster = 0; cluster < ncluster; ++cluster )
   {
      for( bin = 0; bin < nbins; ++bin )
      {
         if( cluster != newcluster )
         {
            qmatrix[newcluster][cluster] +=  sign * solclustering[bin][cluster] * cmatrix[newbin][bin];
            qmatrix[cluster][newcluster] +=  sign * solclustering[bin][cluster] * cmatrix[bin][newbin];
         }
         else
         {
            if( bin != newbin )
            {
               qmatrix[newcluster][newcluster] += sign * (cmatrix[newbin][bin] + cmatrix[bin][newbin])
                  * solclustering[bin][cluster];
            }
         }
      }
   }

   solclustering[newbin][newcluster] = (sign + 1.0) / 2.0;
}

/** initialize the q-matrix from a given (possibly incomplete) clusterassignment */
static
void computeIrrevMat(
   SCIP_Real**           clustering,         /**< the matrix containing the clusterassignment */
   SCIP_Real**           qmatrix,            /**< the returned matrix the transition probability between clusters */
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
               if( i == j )
                  continue;

               qmatrix[k][l] += cmatrix[i][j] * clustering[i][k] * clustering[j][l];
            }
         }
      }
   }
}

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
      objective += ( qmatrix[c][c2] - qmatrix[c2][c]);
      objective += scale * qmatrix[c][c];
   }

   /* if we have no transitions at all then irreversibility should be set to 0 */
   return objective;
}

/** exchange another bin to a different cluster. No bin may be changed twice */
static
SCIP_Bool switchNext(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           cmatrix,            /**< the transition matrix */
   SCIP_Real**           qmatrix,            /**< the irreversibility matrix */
   SCIP_Real**           clustering,         /**< the clusterassignement */
   SCIP_Bool**           binfixed,           /**< array containing information about fixedbins */
   SCIP_Bool*            binprocessed,       /**< has the bin already been switched? */
   int*                  clusterofbin,       /**< contains the cluster each bin is in at the moment */
   int*                  nbinsincluster,     /**< number of bins in each cluster */
   int*                  switchedbin,        /**< the bins swithced in each iteration */
   int*                  switchedcluster,    /**< the cluster to witch the bin was assigned in each iteration */
   SCIP_Real*            switchbound,        /**< the objective achieved in each iteration */
   SCIP_Real*            maxbound,           /**< the best objective value so far */
   int*                  bestlength,         /**< the amount of switches with the best objective value so far */
   int                   iteration           /**< which iteration are we in */
   )
{
   SCIP_Real irrevchg;
   SCIP_Real cohchg;
   SCIP_Real maxboundlocal;
   SCIP_Real scale;
   SCIP_Real oldobjective;
   int bin;
   int k;
   int i;
   int l;
   int indkmin;
   int indlmin;
   int maxbin;
   int maxcluster;
   int nbins = SCIPcycGetNBins(scip);
   int ncluster = SCIPcycGetNCluster(scip);

   scale = SCIPcycGetScale(scip);
   maxboundlocal = -SCIPinfinity(scip);
   oldobjective = getObjective(scip, qmatrix, scale, ncluster);
   maxbin = -1;
   maxcluster = -1;

   assert(isPartition(scip, clustering, nbins, ncluster));

   for( bin = 0; bin < nbins; ++bin )
   {
      if( binprocessed[bin] || nbinsincluster[clusterofbin[bin]] == 1 )
         continue;
      k = clusterofbin[bin];

      assert(SCIPisFeasEQ(scip, clustering[bin][k], 1.0));

      /* calculate the irreversibility and coherence after bin was moved from k to l */
      for( l = 0; l < ncluster; ++l )
      {
         if( binfixed[bin][k] || binfixed[bin][l] )
            continue;

         if( k == l )
            continue;

         if( k == 0 )
            indkmin = ncluster - 1;
         else
            indkmin = k - 1;

         if( l == 0 )
            indlmin = ncluster - 1;
         else
            indlmin = l - 1;

         irrevchg = 0;
         cohchg = 0;

         assert(SCIPisZero(scip, clustering[bin][l]));

         for( i = 0; i < nbins; ++i )
         {
            /*irreversibility change */
            irrevchg -= clustering[i][indkmin] * (cmatrix[i][bin] - cmatrix[bin][i]);
            irrevchg -= clustering[i][(k+1) % ncluster] * (cmatrix[bin][i] - cmatrix[i][bin]);

            clustering[bin][k] = 0;
            clustering[bin][l] = 1;

            irrevchg += clustering[i][indlmin] * (cmatrix[i][bin] - cmatrix[bin][i]);
            irrevchg += clustering[i][(l+1) % ncluster] * (cmatrix[bin][i] - cmatrix[i][bin]);

            clustering[bin][k] = 1;
            clustering[bin][l] = 0;

            /*coherence change */
            if( i != bin )
            {
               cohchg -= clustering[i][k] * (cmatrix[bin][i]+ cmatrix[i][bin]);
               cohchg += clustering[i][l] * (cmatrix[bin][i]+ cmatrix[i][bin]);
            }
         }

         if( oldobjective + irrevchg + scale * cohchg > maxboundlocal )
         {
            maxboundlocal = oldobjective + irrevchg + scale * cohchg;
            maxbin = bin;
            maxcluster = l;
         }
      }
   }

   if( maxbin == -1 )
      return FALSE;

   assert(maxbin >= 0 && maxcluster >= 0);
   assert(maxbin < nbins && maxcluster < ncluster);\

   /* assign the exchange and update all saving-structures */
   setBinToCluster(clustering, cmatrix, qmatrix, maxbin, clusterofbin[maxbin], FALSE, nbins, ncluster);
   setBinToCluster(clustering, cmatrix, qmatrix, maxbin, maxcluster, TRUE, nbins, ncluster);

   nbinsincluster[clusterofbin[maxbin]]--;
   nbinsincluster[maxcluster]++;

   clusterofbin[maxbin] = maxcluster;
   binprocessed[maxbin] = TRUE;

   switchedbin[iteration] = maxbin;
   switchedcluster[iteration] = maxcluster;
   switchbound[iteration] = getObjective(scip, qmatrix, scale, ncluster);

   if( switchbound[iteration] > *maxbound )
   {
      *maxbound = switchbound[iteration];
      *bestlength = iteration;
   }

   return TRUE;
}

/** Create a solution in scip from the clustering */
static
SCIP_RETCODE createSwitchSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic pointer */
   SCIP_Real**           cmatrix,            /**< the transition matrix */
   SCIP_Real**           qmatrix,            /**< the projected transition matrix using the clustering */
   SCIP_Bool**           binfixed,           /**< matrix that tells which bin-variables cannot be changed */
   SCIP_Real**           startclustering,    /**< the start-assignment */
   SCIP_RESULT*          result,             /**< result pointer */
   int                   nbins,              /**< the number of states */
   int                   ncluster            /**< the number of clusters */
   )
{
   SCIP_SOL* bestsol;
   SCIP_SOL* worksol;
   SCIP_Real** clustering;
   SCIP_Real** solclustering;
   SCIP_Bool* binprocessed;
   SCIP_Real max;
   SCIP_Real objective;
   SCIP_Real* switchbound;
   SCIP_Real maxbound;
   SCIP_Bool heurpossible = TRUE;
   SCIP_Bool feasible;
   int c;
   int i;
   int* nbinsincluster; /*lint !e771*/
   int* switchedbin;
   int* switchedcluster;
   int* clusterofbin;
   int bestlength;
   int nrswitches;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &binprocessed, nbins) );
   SCIP_CALL( SCIPallocBufferArray(scip, &switchbound, nbins) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nbinsincluster, ncluster) );
   SCIP_CALL( SCIPallocBufferArray(scip, &switchedbin, nbins) );
   SCIP_CALL( SCIPallocBufferArray(scip, &switchedcluster, nbins) );
   SCIP_CALL( SCIPallocBufferArray(scip, &clusterofbin, nbins) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &clustering, nbins) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &solclustering, nbins) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocClearBufferArray(scip, &clustering[i], ncluster) );    /*lint !e866*/
      SCIP_CALL( SCIPallocClearBufferArray(scip, &solclustering[i], ncluster) ); /*lint !e866*/
   }

   /* copy the solution so that we may change it and still keep the original one*/
   for( c = 0; c < ncluster; ++c )
   {
      nbinsincluster[c] = 0;
   }

   for( i = 0; i < nbins; ++i )
   {
      for( c = 0; c < ncluster; ++c )
      {
         solclustering[i][c] = startclustering[i][c];

         if( SCIPisFeasEQ(scip, startclustering[i][c], 1.0) )
         {
            clusterofbin[i] = c;
            nbinsincluster[c]++; /*lint !e771*/
         }
      }

      binprocessed[i] = FALSE;
   }

   bestsol = SCIPgetBestSol(scip);

   while( heurpossible )
   {
      /* we run the heuristic until we cannot find any more improvement */
      for( c = 0; c < ncluster; ++c )
      {
         nbinsincluster[c] = 0;
      }

      for( i = 0; i < nbins; ++i )
      {
         for( c = 0; c < ncluster; ++c )
         {
            clustering[i][c] = solclustering[i][c];

            if( SCIPisFeasEQ(scip, solclustering[i][c], 1.0) )
            {
               clusterofbin[i] = c;
               nbinsincluster[c]++;
            }
         }

         binprocessed[i] = FALSE;
      }

      bestlength = -1;
      nrswitches = nbins;

      /* initialize qmatrix */
      computeIrrevMat(solclustering, qmatrix, cmatrix, nbins, ncluster);
      maxbound = SCIPgetSolOrigObj(scip, bestsol);

      /* main part of the heuristic. States are switched until every state has been exchanged exactly once */
      for( i = 0; i < nrswitches; ++i )
      {
         if( !switchNext(scip, cmatrix, qmatrix, clustering, binfixed, binprocessed, clusterofbin,
            nbinsincluster, switchedbin, switchedcluster, switchbound, &maxbound, &bestlength, i) )
         {
            nrswitches = i;
            break;
         }
      }

      /* select the clustering with the best objective and reconstruct it from the start clustering */
      for( i = 0; i <= bestlength; ++i )
      {
         for( c = 0; c < ncluster; ++c )
         {
            solclustering[switchedbin[i]][c] = 0; /*lint !e771*/
         }

         solclustering[switchedbin[i]][switchedcluster[i]] = 1; /*lint !e771*/
         clusterofbin[switchedbin[i]] = switchedcluster[i];
      }

      computeIrrevMat(solclustering, qmatrix, cmatrix, nbins, ncluster);
      max = getObjective(scip, qmatrix, SCIPcycGetScale(scip), ncluster);
      objective = SCIPgetSolOrigObj(scip, bestsol);
      feasible = FALSE;

      /* if the solution is an improvement we add it to scip */
      if( max > objective )
      {
         SCIP_CALL( SCIPcreateSol(scip, &worksol, heur) );

         assert(isPartition(scip, solclustering, nbins, ncluster));

         SCIP_CALL( assignVars(scip, worksol, solclustering, nbins, ncluster) );
         SCIP_CALL( SCIPtrySolFree(scip, &worksol, FALSE, TRUE, TRUE, TRUE, TRUE, &feasible) );
      }

      if( feasible )
      {
         *result = SCIP_FOUNDSOL;
         objective = max;
      }
      else
      {
         *result = SCIP_DIDNOTFIND;
         heurpossible = FALSE;
      }
   }

   /* free memory */
   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeBufferArray(scip, &solclustering[i]);
      SCIPfreeBufferArray(scip, &clustering[i]);
   }

   SCIPfreeBufferArray(scip, &solclustering);
   SCIPfreeBufferArray(scip, &clustering);
   SCIPfreeBufferArray(scip, &clusterofbin);
   SCIPfreeBufferArray(scip, &switchedcluster);
   SCIPfreeBufferArray(scip, &switchedbin);
   SCIPfreeBufferArray(scip, &nbinsincluster);
   SCIPfreeBufferArray(scip, &switchbound);
   SCIPfreeBufferArray(scip, &binprocessed);

   return SCIP_OKAY;
}

/** method that randomly creates a different solution from a given solution. From each cluster, half the states are
 * randomly selected and added to the next cluster. */
static
SCIP_RETCODE permuteStartSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           startclustering,    /**< the solution to be permuted */
   SCIP_RANDNUMGEN*      rnd,                /**< a random number generator */
   int                   nbins,              /**< the number of states */
   int                   ncluster            /**< the number of clusters */
   )
{
   int i;
   int t;
   int c;
   int rndcluster;
   int pushed;
   int* binsincluster;
   int **bins;

   SCIP_CALL( SCIPallocBufferArray(scip, &binsincluster, ncluster) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bins, ncluster) );

   for( t = 0; t < ncluster; ++t )
   {
      binsincluster[t] = 0;

      for( i = 0; i < nbins; ++i )
      {
         if( SCIPisPositive(scip, startclustering[i][t]) )
            binsincluster[t]++;
      }

      SCIP_CALL( SCIPallocClearBufferArray(scip, &bins[t], binsincluster[t]) ); /*lint !e866*/

      c = 0;

      for( i = 0; i < nbins; ++i )
      {
         if( SCIPisPositive(scip, startclustering[i][t]) )
         {
            bins[t][c] = i;
            c++;
         }
      }
   }

   for( t = 0; t < ncluster; ++t )
   {
      pushed = 0;

      while(pushed < binsincluster[t] / 2) /*lint !e771*/
      {
         rndcluster = bins[t][SCIPrandomGetInt(rnd, 0, binsincluster[t] - 1)];

         if( rndcluster == nbins -1 )
            continue;
         if( SCIPisZero(scip, startclustering[rndcluster][t]) )
            continue;

         startclustering[rndcluster][t] = 0;
         startclustering[rndcluster][phi(t,ncluster)] = 1;
         pushed++;
      }

      SCIPfreeBufferArray(scip, &bins[t]);
   }

   SCIPfreeBufferArray(scip, &bins);
   SCIPfreeBufferArray(scip, &binsincluster);

   return SCIP_OKAY;
}

/** executes the exchange heuristic for a given solution */
static
SCIP_RETCODE runCyckerlin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic pointer */
   SCIP_SOL*             sol,                /**< given solution */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_Real** startclustering;
   SCIP_Bool** binfixed;
   SCIP_Real** cmatrix;
   SCIP_Real** qmatrix;
   SCIP_RANDNUMGEN* rnd;
   int* clusterofbin;
   int* nbinsincluster;
   int nbins;
   int ncluster;
   int i;
   int c;

   /* get problem variables */
   nbins = SCIPcycGetNBins(scip);
   ncluster = SCIPcycGetNCluster(scip);
   cmatrix = SCIPcycGetCmatrix(scip);
   SCIP_CALL( SCIPcreateRandom(scip, &rnd, DEFAULT_RANDSEED, TRUE) );

   assert(nbins >= 0);
   assert(ncluster >= 0);

   /* allocate Memory */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &startclustering, nbins) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &clusterofbin, nbins) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &binfixed, nbins) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &qmatrix, ncluster) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nbinsincluster, ncluster) );

   for( c = 0; c < ncluster; ++c )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &qmatrix[c], ncluster) ); /*lint !e866*/
   }
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocClearBufferArray(scip, &startclustering[i], ncluster) ); /*lint !e866*/
      SCIP_CALL( SCIPallocClearBufferArray(scip, &binfixed[i], ncluster) ); /*lint !e866*/
   }

   /* get the solution values from scip */
   SCIP_CALL( getSolutionValues(scip, sol, startclustering, binfixed, clusterofbin, nbinsincluster) );

   if( isPartition(scip, startclustering, nbins, ncluster) )
   {
      SCIP_CALL( createSwitchSolution(scip, heur, cmatrix, qmatrix, binfixed, startclustering, result, nbins, ncluster) );
      for( i = 0; i < MAXPERMUTATIONS; ++i )
      {
         SCIP_CALL( permuteStartSolution(scip, startclustering, rnd, nbins, ncluster) );

         assert(isPartition(scip, startclustering, nbins, ncluster));

         SCIP_CALL( createSwitchSolution(scip, heur, cmatrix, qmatrix, binfixed, startclustering, result, nbins, ncluster) );
      }
   }

   /* free all data-structures */
   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeBufferArray(scip, &binfixed[i]);
      SCIPfreeBufferArray(scip, &startclustering[i]);
   }

   for( c = 0; c < ncluster; ++c )
   {
      SCIPfreeBufferArray(scip, &qmatrix[c]);
   }

   SCIPfreeBufferArray(scip, &nbinsincluster);
   SCIPfreeBufferArray(scip, &qmatrix);
   SCIPfreeBufferArray(scip, &binfixed);
   SCIPfreeBufferArray(scip, &clusterofbin);
   SCIPfreeBufferArray(scip, &startclustering);

   SCIPfreeRandom(scip, &rnd);

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyCyckerlin)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurCycKerlin(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeCyckerlin)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   SCIPfreeMemoryArray(scip, &(heurdata->candidates));

   /* free heuristic data */
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolCyckerlin)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* reset the timing mask to its default value */
   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   heurdata->candlength = 0;

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitCyckerlin)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   heurdata->ncandidates = 0;
   heurdata->candlength = 10;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(heurdata->candidates), heurdata->candlength) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecCyckerlin)
{  /*lint --e{715}*/
   SCIP_Real objective;
   SCIP_HEURDATA* heurdata;
   int i;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   heurdata = SCIPheurGetData(heur);

   assert(NULL != heurdata);

   /* reset the timing mask to its default value (at the root node it could be different) */
   if( SCIPgetNNodes(scip) > 1 )
      SCIPheurSetTimingmask(heur, HEUR_TIMING);
   for( i = 0; i < heurdata->ncandidates; i++ )
   {
      objective = SCIPgetSolOrigObj(scip, heurdata->candidates[i]);

      if( !SCIPisZero(scip, objective) )
      {
         SCIP_CALL( runCyckerlin(scip, heur, heurdata->candidates[i], result) );
      }

      heurdata->candidates[i] = NULL;
   }

   heurdata->ncandidates = 0;
   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the oneopt primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurCycKerlin(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;

   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
      HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
      HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecCyckerlin, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyCyckerlin) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeCyckerlin) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolCyckerlin) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitCyckerlin) );

   return SCIP_OKAY;
}
