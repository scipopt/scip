
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

/**@file   heur_spakerlin.c
 * @brief  improvement heuristic that trades binary variables between clusters.
 * @author Leon Eifler
 */

/*---+---- 1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "probdata_spa.h"
#include "heur_spakerlin.h"
#include "scip/misc.h"

#define HEUR_NAME             "spakerlin"
#define HEUR_DESC             "switch heuristic that tries to improve solution by trading bins betweeen clusters, similar to the famous kernighan/lin heuristic"
#define HEUR_DISPCHAR         '@'
#define HEUR_PRIORITY         500
#define HEUR_FREQ             10
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define MAXROUNDS             5
#define MAXPERMUTATIONS       5
#define DEFAULT_RANDSEED       177           /**< random seed                                                                       */
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE          /**< does the heuristic use a secondary SCIP instance? */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   lastsolindex;       /**< index of the last solution for which oneopt was performed */
   int                   currentrounds;
};


/*
 * Local methods
 */


/** Get the bin-var assignment from scip and save it as a matrix */
static
SCIP_RETCODE getSolutionValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             bestsol,            /**< The solution */
   SCIP_Real**           solclustering,      /**< Matrix to save the bin-vars*/
   SCIP_Bool**           binfixed,           /**< Matrix to save if a bin is fixed in scip */
   int*                  clusterofbin,       /**< Array containing the cluster of each bin */
   int*                  nbinsincluster      /**< Number of bins in each cluster */
)
{
   int i;
   int c;
   int k;
   int j;
   SCIP_VAR*** binvars;
   SCIP_VAR**** edgevars;
   int nbins;
   int ncluster;
   char model;

   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);
   assert(nbins > 0 && ncluster > 0 && nbins > ncluster);
   binvars = SCIPspaGetBinvars(scip);
   edgevars = SCIPspaGetEdgevars(scip);
   assert(binvars != NULL && edgevars != NULL);
   SCIP_CALL( SCIPgetCharParam(scip, "model", &model) );

   if( model == 's' )
   {
      /* Get the bin-variable values from the solution */
      for( i = 0; i < nbins; ++i )
      {
         for( c = 0; c < ncluster; ++c )
         {
            binfixed[i][c] = FALSE;
            if( binvars[i][c] != NULL )
            {
               if( (SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(binvars[i][c]), SCIPvarGetLbGlobal(binvars[i][c]))) )
               {
                  solclustering[i][c] = SCIPgetSolVal(scip, bestsol, binvars[i][c]);
                  binfixed[i][c] = TRUE;
                  if( SCIPisFeasEQ(scip, solclustering[i][c], 1) )
                  {
                     clusterofbin[i] = c;
                     nbinsincluster[c]++;
                  }
               }else
               {
                  SCIP_Real solval = SCIPgetSolVal(scip, bestsol, binvars[i][c]);
                  assert(SCIPisFeasIntegral(scip, solval));
                  solclustering[i][c] = solval;
                  if( SCIPisFeasEQ(scip, solval, 1.0) )
                  {
                     clusterofbin[i] = c;
                     nbinsincluster[c]++;
                  }
               }
            }
            else
            {
               binfixed[i][c] = TRUE;
               for( k = 0; k < 3; ++k )
               {
                  if( k != c )
                  {
                     j = 0;
                     while((binvars[j][k] == NULL || SCIPvarGetStatus(binvars[j][k]) == SCIP_VARSTATUS_FIXED || j == i))
                     {
                        j++;
                        if ( j == nbins - 1 )
                           break;
                     }
                     if( j < nbins )
                     {
                        if( NULL != edgevars[i][j] && NULL != edgevars[i][j][k] && SCIPvarIsActive(edgevars[i][j][k]) )
                        {
                           solclustering[i][c] = 1;
                           clusterofbin[i] = c;
                           nbinsincluster[c]++;
                        }
                        else
                           solclustering[i][c] = 0;
                     }
                  }
               }
            }
         }
      }
   }
   else if ( model == 'e' )
   {
      SCIP_Bool processed[nbins];
      int nextbin;
      int currentbin;
      int currentcluster;
      int nprocessed;
      nprocessed = 0;
      currentbin = 0;
      for( i = 0; i < nbins; ++i )
      {
         processed[i] = FALSE;
      }
      currentcluster = 0;
      while( nprocessed < nbins && currentcluster < ncluster )
      {
         nextbin = -1;
         solclustering[currentbin][currentcluster] = 1;
         processed[currentbin] = TRUE;
         nprocessed++;
         for( i = currentbin + 1; i < nbins; ++i )
         {
            if( SCIPgetSolVal(scip, bestsol, edgevars[i][currentbin][0]) > 0 )
            {
               solclustering[i][currentcluster] = 1;
               processed[i] = TRUE;
               nprocessed++;
            }
            else if( !processed[i] && nextbin == -1 )
            {
               nextbin = i;
            }
         }
         currentcluster++;
         currentbin = nextbin;
      }
      for( i = 0; i < nbins; ++i )
      {
         if( !processed[i] )
            solclustering[i][0] = 1;
      }

   }
   else
      SCIPABORT();
   return SCIP_OKAY;
}


/** Set a bin to a new cluster, update the qmatrix. */
static
void setBinToCluster(
   SCIP_Real**           solclustering,      /**< The matrix with the clustering of the bins */
   SCIP_Real**           cmatrix,            /**< The transition matrix*/
   SCIP_Real**           qmatrix,            /**< The matrix containing the transition probabilities between clusters*/
   int                   newbin,             /**< The bin to be changed*/
   int                   newcluster,         /**< The cluster where the bin is changed*/
   SCIP_Bool             setone,             /**< TRUE if the assignment is switched from 0 to 1, FALSE if it is switched from 1 to 0*/
   int                   nbins,              /**< The number of bins*/
   int                   ncluster            /**< The number of clusters*/
)
{
   int bin;
   int cluster;

   int sign = setone ? 1 : -1;
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
               qmatrix[newcluster][newcluster] += sign * (cmatrix[newbin][bin] + cmatrix[bin][newbin]) * solclustering[bin][cluster];
         }
      }
   }
   solclustering[newbin][newcluster] = (sign + 1) / 2;
}

/**  Initialize the q-matrix from a given (possibly incomplete) clusterassignment */
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

/** exchange another bin to a different cluster. No bin may be changed twice */
static
SCIP_Bool switchNext(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           cmatrix,            /**< The transition matrix */
   SCIP_Real**           qmatrix,            /**< The irreversibility matrix */
   SCIP_Real**           clustering,         /**< The clusterassignement */
   SCIP_Bool**           binfixed,           /**< Array containing information about fixedbins */
   SCIP_Bool*            binprocessed,       /**< has the bin already been switched? */
   int*                  clusterofbin,       /**< Contains the cluster each bin is in at the moment */
   int*                  nbinsincluster,     /**< Number of bins in each cluster */
   int*                  switchedbin,        /**< The bins swithced in each iteration */
   int*                  switchedcluster,    /**< The cluster to witch the bin was assigned in each iteration */
   SCIP_Real*            switchbound,        /**< The objective achieved in each iteration */
   SCIP_Real*            maxbound,           /**< The best objective value so far */
   int*                  bestlength,         /**< The amount of switches with the best objective value so far */
   int                   iteration           /**< which iteration are we in */
)
{
   int bin;
   int k;
   int i;
   int l;
   int indkmin;
   int indlmin;
   SCIP_Real irrevchg;
   SCIP_Real cohchg;
   SCIP_Real maxboundlocal;
   SCIP_Real scale;
   SCIP_Real oldobjective;
   int maxbin;
   int maxcluster;
   int nbins = SCIPspaGetNrBins(scip);
   int ncluster = SCIPspaGetNrCluster(scip);


   scale = SCIPspaGetScale(scip);
   maxboundlocal = -SCIPinfinity(scip);
   oldobjective = getObjective(scip, qmatrix, scale, ncluster);
   maxbin = -1;
   maxcluster = -1;
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

/** */
static
SCIP_RETCODE createSwitchSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< Heuristic pointer */
   SCIP_Real**           cmatrix,            /**< The transition matrix */
   SCIP_Real**           qmatrix,            /**< The projected transition matrix using the clustering */
   SCIP_Bool**           binfixed,           /**< Matrix that tells which bin-variables cannot be changed */
   SCIP_Real**           startclustering,    /**< The start-assignment */
   SCIP_RESULT*          result,             /**< Result pointer */
   int                   nbins,              /**< The number of states */
   int                   ncluster            /**< The number of clusters */
)
{
   int c,i;
   int nbinsincluster[ncluster];
   SCIP_Real** clustering;
   SCIP_Real** solclustering;
   int clusterofbin[nbins];
   SCIP_Bool binprocessed[nbins];
   int switchedbin[nbins];
   int switchedcluster[nbins];
   SCIP_Real max;
   SCIP_Bool feasible;
   SCIP_Real objective;
   SCIP_Real switchbound[nbins];
   SCIP_Bool heurpossible = TRUE;
   int bestlength;
   int nrswitches;
   SCIP_Real maxbound;
   SCIP_SOL* bestsol;
   SCIP_SOL* worksol;

   /* allocate memory */
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &clustering, nbins) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &solclustering, nbins) );

   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &clustering[i], ncluster) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &solclustering[i], ncluster) );
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
            nbinsincluster[c]++;
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
         if( !switchNext(scip, cmatrix, qmatrix, clustering, binfixed, binprocessed, clusterofbin, nbinsincluster, switchedbin, switchedcluster, switchbound, &maxbound, &bestlength, i) )
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
            solclustering[switchedbin[i]][c] = 0;
         }
         solclustering[switchedbin[i]][switchedcluster[i]] = 1;
         clusterofbin[switchedbin[i]] = switchedcluster[i];
      }
      computeIrrevMat(solclustering, qmatrix, cmatrix, nbins, ncluster);
      max = getObjective(scip, qmatrix, SCIPspaGetScale(scip), ncluster);
      objective = SCIPgetSolOrigObj(scip, bestsol);
      feasible = FALSE;
      /* if the solution is an improvement we add it to scip */
      if( max > objective )
      {
         SCIP_CALL( SCIPcreateSol(scip, &worksol, heur) );
         assignVars(scip, worksol, solclustering, nbins, ncluster);
         SCIPtrySolFree(scip, &worksol, FALSE, TRUE, TRUE, TRUE, TRUE, &feasible);
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
   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeMemoryArray(scip, &clustering[i]);
      SCIPfreeMemoryArray(scip, &solclustering[i]);
   }
   SCIPfreeMemoryArray(scip, &clustering);
   SCIPfreeMemoryArray(scip, &solclustering);
   return SCIP_OKAY;
}

/** Method that randomly creates a different solution from a given solution. From each cluster, half the states are
 * randomly selected and added to the next cluster.*/
static
void permuteStartSolution(
   SCIP_Real**           startclustering,
   SCIP_RANDNUMGEN*      rand,
   int                   nbins,
   int                   ncluster
)
{
   int i,t;
   int c;
   int rnd;
   int pushed;
   int binsincluster[ncluster];
   int **bins;

   SCIPallocMemoryArray(scip, &bins, ncluster);

   for( t = 0; t < ncluster; ++t )
   {
      binsincluster[t] = 0;
      for( i = 0; i < nbins; ++i )
      {
         if( 0 < startclustering[i][t])
            binsincluster[t]++;
      }
      SCIPallocClearMemoryArray(scip, &bins[t], binsincluster[t]);
      c = 0;
      for( i = 0; i < nbins; ++i )
      {
         if( 0 < startclustering[i][t])
         {
            bins[t][c] = i;
            c++;
         }
      }
   }

   for( t = 0; t < ncluster; ++t )
   {
      pushed = 0;
      while(pushed < binsincluster[t] / 2)
      {
         rnd = bins[t][SCIPrandomGetInt(rand, 0, binsincluster[t] - 1)];
         if( rnd == nbins -1 )
            continue;
         if( startclustering[rnd][t] == 0 )
            continue;
         startclustering[rnd][t] = 0;
         startclustering[rnd][phi(t,ncluster)] = 1;
         pushed++;
      }
      SCIPfreeMemory(scip, &bins[t]);
   }
   SCIPfreeMemory(scip, &bins);
}
/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopySpakerlin)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurSpakerlin(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSpakerlin)
{   /*lint --e{715}*/
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
SCIP_DECL_HEUREXITSOL(heurExitsolSpakerlin)
{
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* reset the timing mask to its default value */
   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitSpakerlin)
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize last solution index */
   heurdata->lastsolindex = -1;
   heurdata->currentrounds = 0;

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSpakerlin)
{  /*lint --e{715}*/

   SCIP_SOL* bestsol;                        /* incumbent solution */
   SCIP_HEURDATA* heurdata;

   SCIP_Real** startclustering;                 /* the assignment given from the solution */
   SCIP_Bool** binfixed;                      /* The bins that are fixed from scip */
   SCIP_Real** cmatrix;
   SCIP_Real** qmatrix;
   SCIP_RANDNUMGEN* rand;

   int* clusterofbin;                         /* hold the cluster that each bin is in */
   int* nbinsincluster;                       /* The number of bins ins each cluster */

   SCIP_Real objective;                       /* The value of the objective function */

   int nbins;
   int ncluster;
   int c;
   int i;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* for now: do not use heurisitc if weighted objective is used */

   /* we only want to process each solution once */
   heurdata = SCIPheurGetData(heur);
   bestsol = SCIPgetBestSol(scip);

   if( bestsol == NULL || heurdata->currentrounds == MAXROUNDS )
      return SCIP_OKAY;

   if( heurdata->lastsolindex == SCIPsolGetIndex(bestsol) )
      heurdata->currentrounds++;
   else
      heurdata->currentrounds = 0;


   /* reset the timing mask to its default value (at the root node it could be different) */
   if( SCIPgetNNodes(scip) > 1 )
      SCIPheurSetTimingmask(heur, HEUR_TIMING);


   /* get problem variables */
   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);
   cmatrix = SCIPspaGetCmatrix(scip);
   SCIPrandomCreate(&rand, SCIPblkmem(scip),SCIPinitializeRandomSeed(scip, DEFAULT_RANDSEED));

   /* we do not want to run the heurtistic if there is no 'flow' between the clusters.
    * in case of a (ideally) full reversible problem there cannot be a better solution, in the other case, i.e., the
    * problem has irreversible parts, it seems the heuristic will not find solutions respecting the coherence conditions
    */
   objective = SCIPgetSolOrigObj(scip, bestsol);
   if( SCIPisZero(scip, objective) )
      return SCIP_OKAY;

   assert(nbins >= 0);
   assert(ncluster >= 0);

   /* allocate Memory */
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &startclustering, nbins) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &clusterofbin, nbins) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &binfixed, nbins) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &qmatrix, ncluster) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &nbinsincluster, ncluster) );

   for( c = 0; c < ncluster; ++c )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &qmatrix[c], ncluster) );
   }
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &startclustering[i], ncluster) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &binfixed[i], ncluster) );
   }

   /* get the solution values from scip */
   SCIP_CALL( getSolutionValues(scip, bestsol, startclustering, binfixed, clusterofbin, nbinsincluster) );

   if( isPartition(scip, startclustering, nbins, ncluster) )
   {
      SCIP_CALL( createSwitchSolution(scip, heur, cmatrix, qmatrix, binfixed, startclustering, result, nbins, ncluster) );
      for( i = 0; i < MAXPERMUTATIONS; ++i )
      {
         permuteStartSolution(startclustering, rand, nbins, ncluster);
         assert(isPartition(scip, startclustering, nbins, ncluster) );
         SCIP_CALL( createSwitchSolution(scip, heur, cmatrix, qmatrix, binfixed, startclustering, result, nbins, ncluster) );
      }
   }
   /* update the found solution, so that we do not run again immediatly */

   heurdata->lastsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));

   /* free all data-structures */

   /* free all data-structures */
   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeMemoryArray(scip, &startclustering[i]);
      SCIPfreeMemoryArray(scip, &binfixed[i]);
   }
   for( c = 0; c < ncluster; ++c )
   {
      SCIPfreeMemoryArray(scip, &qmatrix[c]);
   }

   SCIPfreeMemoryArray(scip, &qmatrix);
   SCIPfreeMemoryArray(scip, &startclustering);
   SCIPfreeMemoryArray(scip, &binfixed);
   SCIPfreeMemoryArray(scip, &clusterofbin);
   SCIPfreeMemoryArray(scip, &nbinsincluster);


   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the oneopt primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSpakerlin(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Oneopt primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
      HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
      HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecSpakerlin, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopySpakerlin) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeSpakerlin) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolSpakerlin) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitSpakerlin) );
   return SCIP_OKAY;
}
