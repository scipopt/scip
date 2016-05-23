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
 * @brief  improvement heuristic that trades spa-binary variables between clusters
 * @author Leon Eifler
 */

/*---+---- 1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "probdata_spa.h"
#include "heur_spakerlin.h"
#include "scip/cons_and.h"
/* @note If the heuristic runs in the root node, the timing is changed to (SCIP_HEURTIMING_DURINGLPLOOP |
 *       SCIP_HEURTIMING_BEFORENODE), see SCIP_DECL_HEURINITSOL callback.
 */

#define HEUR_NAME             "spakerlin"
#define HEUR_DESC             "switch heuristic that tries to improve solution by trading bins betweeen clusters, similar to the famous kernighan/lin heuristic"
#define HEUR_DISPCHAR         '@'
#define HEUR_PRIORITY         536870900
#define HEUR_FREQ             10
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE          /**< does the heuristic use a secondary SCIP instance? */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   lastsolindex;       /**< index of the last solution for which oneopt was performed */
};


/*
 * Local methods
 */

/**
 * This method changes the cluster of @p newbin and updates the qmatrix as well as the clustermatrix.
 */

static
void getSolutionValues(
   SCIP*                 scip,
   SCIP_SOL*             bestsol,
   SCIP_Real**           solclustering,
   SCIP_Bool**            binfixed,
   int*                  clusterofbin
)
{
   int i;
   int c;
   int k;
   int j;
   SCIP_VAR*** binvars;
   SCIP_VAR***** edgevars;
   int nbins;
   int ncluster;

   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);
   assert(nbins > 0 && ncluster > 0 && nbins > ncluster);
   binvars = SCIPspaGetBinvars(scip);
   edgevars = SCIPspaGetEdgevars(scip);
   assert(binvars != NULL && edgevars != NULL);

   /* Get the bin-variable values from the solution */
   for( i = 0; i < nbins; ++i )
   {
      for( c = 0; c < ncluster; ++c )
      {
         binfixed[i][c] = FALSE;
         if( binvars[i][c] != NULL )
         {
            if( (SCIPisEQ(scip, SCIPvarGetUbGlobal(binvars[i][c]), SCIPvarGetLbGlobal(binvars[i][c]))) )
            {
               solclustering[i][c] = SCIPgetSolVal(scip, bestsol, binvars[i][c]);
               binfixed[i][c] = TRUE;
               if( SCIPisEQ(scip, solclustering[i][c], 1) )
                  clusterofbin[i] = c;
            }else
            {
               SCIP_Real solval = SCIPgetSolVal(scip, bestsol, binvars[i][c]);
               assert(SCIPisIntegral(scip, solval));
               solclustering[i][c] = solval;
               if( SCIPisEQ(scip, solval, 1.0) )
                  clusterofbin[i] = c;
            }
         }
         else
         {
            binfixed[i][c] = TRUE;
            for( k = 0; k < ncluster; ++k )
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
                     if( NULL != edgevars[i][j][c][k] && SCIPvarIsActive(edgevars[i][j][c][k]) )
                     {
                        solclustering[i][c] = 1;
                        clusterofbin[i] = c;
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

static
SCIP_Bool isPartition(
   SCIP*                 scip,
   SCIP_Real**           solclustering,
   int                   nbins,
   int                   ncluster
)
{
   int i;
   int j;
   /* check if the assignment violates paritioning, e.g. because we are in a subscip */
   for( i = 0; i < nbins; ++i )
   {
      int sum = 0;
      for( j = 0; j < ncluster; ++j )
      {
         if( solclustering[i][j] )
            sum += solclustering[i][j];
      }
      if( !SCIPisEQ(scip, sum, 1.0) )
         return FALSE;
   }
   return TRUE;
}

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
            if( bin == newbin )
               qmatrix[newcluster][newcluster] += sign * cmatrix[newbin][bin];
            else
               qmatrix[newcluster][newcluster] += sign * (cmatrix[newbin][bin] + cmatrix[bin][newbin]) * solclustering[bin][cluster];
         }
      }
   }
   solclustering[newbin][newcluster] = (sign + 1) / 2;
}

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


/** exchange another bin to a different cluster. No bin may be changed twice */
static
SCIP_Bool switchNext(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           cmatrix,            /**< The transition matrix */
   SCIP_Real**           qmatrix,            /**< The irreversibility matrix */
   SCIP_Real**           clustering,         /**< The clusterassignement */
   SCIP_Bool**            binfixed,           /**< Array containing information about fixedbins */
   SCIP_Bool*            binprocessed,       /**< has the bin already been switched? */
   int*                  clusterofbin,       /**< Contains the cluster each bin is in at the moment */
   int*                  switchedbin,        /**< The bins swithced in each iteration */
   int*                  switchedcluster,    /**< The cluster to witch the bin was assigned in each iteration */
   SCIP_Real*            switchbound,        /**< The objective achieved in each iteration */
   SCIP_Real*            switchcoherence,    /**< Array to save the coherence after an exchange */
   SCIP_Real*            maxbound,           /**< The best objective value so far */
   SCIP_Real             coherence,          /**< The necessary cohernce. We do not want to produce infeasible solutions */
   int*                  bestlength,         /**< The amount of switches with the best objective value so far */
   int                   iteration           /**< which iteration are we in */
)
{
   int bin;
   int k;
   int i;
   int l;
   SCIP_Real qkl;
   SCIP_Real qlk;
   SCIP_Real maxboundlocal;
   SCIP_Real newcoherence;
   int maxbin;
   int maxcluster;
   int nbins = SCIPspaGetNrBins(scip);
   int ncluster = SCIPspaGetNrCluster(scip);

   maxboundlocal = -SCIPinfinity(scip);
   maxbin = -1;
   maxcluster = -1;
   for( bin = 0; bin < nbins; ++bin )
   {
      if( binprocessed[bin] )
         continue;
      k = clusterofbin[bin];
      assert(SCIPisEQ(scip, clustering[bin][k], 1.0));
      /* calculate the irreversibility after bin was moved from k to l */
      for( l = 0; l < ncluster; ++l )
      {
         qkl = qmatrix[k][l];
         qlk = qmatrix[l][k];
         if( binfixed[bin][k] || binfixed[bin][l] )
            continue;
         if( k == l )
            continue;
         assert(SCIPisZero(scip, clustering[bin][l]));
         for( i = 0; i < nbins; ++i )
         {
            qkl -= clustering[i][l] * cmatrix[bin][i];
            qlk -= clustering[i][l] * cmatrix[i][bin];

            clustering[bin][k] = 0;
            clustering[bin][l] = 1;

            qlk += clustering[i][k] * cmatrix[bin][i];
            qkl += clustering[i][k] * cmatrix[i][bin];

            clustering[bin][k] = 1;
            clustering[bin][l] = 0;

         }
         if( REALABS(qkl - qlk) > maxboundlocal )
         {
            maxboundlocal = REALABS(qkl - qlk);
            maxbin = bin;
            maxcluster = l;
         }
      }
   }
   if( maxbin == -1 )
      return FALSE;
   assert(maxbin >= 0 && maxcluster >= 0);
   assert(maxbin < nbins && maxcluster < ncluster);
   assert(maxboundlocal >= 0);
   /* assign the exchange and update all saving-structures */
   setBinToCluster(clustering, cmatrix, qmatrix, maxbin, clusterofbin[maxbin], FALSE, nbins, ncluster);
   setBinToCluster(clustering, cmatrix, qmatrix, maxbin, maxcluster, TRUE, nbins, ncluster);
   newcoherence = SCIPinfinity(scip);
   for( k = 0; k < ncluster; ++k )
   {
      if( qmatrix[k][k] < newcoherence )
         newcoherence = qmatrix[k][k];
   }
   switchcoherence[iteration] = newcoherence;
   clusterofbin[maxbin] = maxcluster;
   binprocessed[maxbin] = TRUE;
   switchedbin[iteration] = maxbin;
   switchedcluster[iteration] = maxcluster;
   switchbound[iteration] = getIrrevBound(scip, qmatrix, ncluster);
   if( switchbound[iteration] > *maxbound && switchcoherence[iteration] > coherence )
   {
      *maxbound = switchbound[iteration];
      *bestlength = iteration;
   }
   return TRUE;
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

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSpakerlin)
{  /*lint --e{715}*/

   SCIP_SOL* bestsol;                        /* incumbent solution */
   SCIP_SOL* worksol;                        /* working solution */
   SCIP_HEURDATA* heurdata;

   SCIP_Real** solclustering;                 /* the assignment given from the solution */
   SCIP_Real** clustering;                    /* the working cluster-assignment. We start with the one given by the solution */
   SCIP_Bool** binfixed;                      /* The bins that are fixed from scip */
   SCIP_Real** cmatrix;
   SCIP_Real** qmatrix;

   int* clusterofbin;                         /* hold the cluster that each bin is in */
   int* switchedbin;                          /* holds the bins that were exchanged in each iteration */
   int* switchedcluster;                      /* The clusters that the bins were switched to in each iteration */
   SCIP_Real* switchbound;                    /* The objective value after each exchange */
   SCIP_Real* switchcoherence;                /* The coherence after each exchange */
   SCIP_Bool* binprocessed;                   /* Has a bin been processed already ? */

   SCIP_Bool heurpossible = TRUE;             /* True if the heuristic can run */
   SCIP_Bool feasible;                        /* True if a found solution is feasible */
   SCIP_Real objective;                       /* The value of the objective function */
   SCIP_Real maxbound;
   SCIP_Real coherence;
   SCIP_Real max;

   int nbins;
   int ncluster;
   int c;
   int i;
   int nrswitches;
   int bestlength = -1;
   char model;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* for now: do not use heurisitc if weighted objective is used */
   model = SCIPspaGetModel(scip);
   if( model == 'w')
      return SCIP_OKAY;

   /* we only want to process each solution once */
   heurdata = SCIPheurGetData(heur);
   bestsol = SCIPgetBestSol(scip);

   if( bestsol == NULL || heurdata->lastsolindex == SCIPsolGetIndex(bestsol) )
      return SCIP_OKAY;


   /* reset the timing mask to its default value (at the root node it could be different) */
   if( SCIPgetNNodes(scip) > 1 )
      SCIPheurSetTimingmask(heur, HEUR_TIMING);


   /* get problem variables */
   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);
   cmatrix = SCIPspaGetCmatrix(scip);
   coherence = SCIPspaGetCoherence(scip);

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
   SCIPallocClearMemoryArray(scip, &solclustering, nbins);
   SCIPallocClearMemoryArray(scip, &clustering, nbins);
   SCIPallocClearMemoryArray(scip, &binfixed, nbins);
   SCIPallocClearMemoryArray(scip, &qmatrix, ncluster);
   SCIPallocClearMemoryArray(scip, &clusterofbin, nbins);

   nrswitches = nbins;

   SCIPallocClearMemoryArray(scip, &switchbound, nrswitches);
   SCIPallocClearMemoryArray(scip, &switchcoherence, nrswitches);
   SCIPallocClearMemoryArray(scip, &switchedbin, nrswitches);
   SCIPallocClearMemoryArray(scip, &switchedcluster, nrswitches);
   SCIPallocClearMemoryArray(scip, &binprocessed, nrswitches);


   for( c = 0; c < ncluster; ++c )
   {
      SCIPallocMemoryArray(scip, &qmatrix[c], ncluster);
   }
   for( i = 0; i < nbins; ++i )
   {
      SCIPallocClearMemoryArray(scip, &solclustering[i], ncluster);
      SCIPallocClearMemoryArray(scip, &clustering[i], ncluster);
      SCIPallocClearMemoryArray(scip, &binfixed[i], ncluster);
   }

   getSolutionValues(scip, bestsol, solclustering, binfixed, clusterofbin);

   if( !isPartition(scip, solclustering, nbins, ncluster) )
      heurpossible = FALSE;

   while( heurpossible )
   {
      /* copy the solution so that we may change it        */
      for( i = 0; i < nbins; ++i )
      {
         for( c = 0; c < ncluster; ++c )
         {
            clustering[i][c] = solclustering[i][c];
            if( SCIPisEQ(scip, solclustering[i][c], 1.0) )
               clusterofbin[i] = c;
         }
         binprocessed[i] = FALSE;
      }
      bestlength = -1;
      /* initialize qmatrix */
      computeIrrevMat(solclustering, qmatrix, cmatrix, nbins, ncluster);
      maxbound = getIrrevBound(scip, qmatrix, ncluster);
      for( i = 0; i < nrswitches; ++i )
      {
         if( !switchNext(scip, cmatrix, qmatrix, clustering, binfixed, binprocessed, clusterofbin, switchedbin, switchedcluster, switchbound, switchcoherence, &maxbound, coherence, &bestlength, i) )
         {
            nrswitches = i;
            break;
         }
         if( bestlength > -1 )
            assert( coherence <= switchcoherence[bestlength] );
      }
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
      max = getIrrevBound(scip, qmatrix, ncluster);
      feasible = FALSE;
      if( max > objective )
      {
         SCIP_CALL( SCIPcreateSol(scip, &worksol, heur) );
         assignVars(scip, worksol, solclustering, nbins, ncluster, qmatrix);
         SCIPtrySolFree(scip, &worksol, TRUE, TRUE, TRUE, TRUE, &feasible);
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
   /* update the found solution, so that we do not run again immediatly */

   heurdata->lastsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));

   /* free all data-structures */

   /* free all data-structures */
   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeMemoryArray(scip, &solclustering[i]);
      SCIPfreeMemoryArray(scip, &clustering[i]);
      SCIPfreeMemoryArray(scip, &binfixed[i]);
   }
   for( c = 0; c < ncluster; ++c )
   {
      SCIPfreeMemoryArray(scip, &qmatrix[c]);
   }

   SCIPfreeMemoryArray(scip, &qmatrix);
   SCIPfreeMemoryArray(scip, &switchbound);
   SCIPfreeMemoryArray(scip, &switchedbin);
   SCIPfreeMemoryArray(scip, &switchedcluster);
   SCIPfreeMemoryArray(scip, &switchcoherence);
   SCIPfreeMemoryArray(scip, &solclustering);
   SCIPfreeMemoryArray(scip, &clustering);
   SCIPfreeMemoryArray(scip, &binfixed);
   SCIPfreeMemoryArray(scip, &clusterofbin);
   SCIPfreeMemoryArray(scip, &binprocessed);

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
