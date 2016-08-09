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

/**@file   heur_spaswitch.c
 * @brief  improvement heuristic that trades spa-binary variables between clusters
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "probdata_spa.h"
#include "heur_spaswitch.h"
#include "scip/cons_and.h"

#define HEUR_NAME             "spaswitch"
#define HEUR_DESC             "switch heuristic that tries to improve solution by trading bins betweeen clusters"
#define HEUR_DISPCHAR         '!'
#define HEUR_PRIORITY         536870885
#define HEUR_FREQ             -1
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

static
SCIP_RETCODE updateSolAuxiliary(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< The SCIP solution */
   SCIP_Real**           clustering,         /**< The matrix with the clusterassignment */
   int                   nbins,              /**< The number of bins */
   int                   ncluster,           /**< The number of cluster */
   SCIP_Real**           qmatrix             /**< The irreversibility matrix */
)
{
   int c;
   int i;
   int q1;
   int q2;
   SCIP_Real epsI;

   SCIP_VAR** absvars;
   SCIP_VAR* targetvar;

   absvars = SCIPspaGetAbsvars(scip);
   targetvar = SCIPspaGetTargetvar(scip);
   epsI = getIrrevBound(scip, qmatrix, ncluster);
   for( c = 0; c < ncluster; ++c )
   {
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
   return SCIP_OKAY;
}

/**
 * assign the variables in scip according to the found clusterassignment.
 */
static
SCIP_RETCODE updateSolAssignment(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< The SCIP solution */
   SCIP_Real**           clustering,         /**< The matrix with the clusterassignment */
   int                   nbins,              /**< The number of bins */
   int                   ncluster,           /**< The number of cluster */
   SCIP_Real**           qmatrix,             /**< The irreversibility matrix */
   int                   bin,
   int                   oldcluster,
   int                   newcluster
)
{
   int i;
   int c;
   SCIP_VAR* var;
   SCIP_VAR*** binvars;
   SCIP_VAR*****  edgevars;

   assert(nbins > 0 && ncluster > 0);

   binvars = SCIPspaGetBinvars(scip);
   edgevars = SCIPspaGetEdgevars(scip);


   /** change the values of the binvars */

   if( SCIPvarIsTransformed(binvars[bin][oldcluster]) )
      var = binvars[bin][oldcluster];
   else
      var = SCIPvarGetTransVar(binvars[bin][oldcluster]);
   /* check if the clusterassignment ist feasible for the variable bounds. If not do not assign the variable */
   if( var != NULL && SCIPisLE(scip, SCIPvarGetLbGlobal(var), 0) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
      SCIP_CALL( SCIPsetSolVal( scip, sol, var, 0) );

   if( SCIPvarIsTransformed(binvars[bin][newcluster]) )
      var = binvars[bin][newcluster];
   else
      var = SCIPvarGetTransVar(binvars[bin][newcluster]);
   /* check if the clusterassignment ist feasible for the variable bounds. If not do not assign the variable */
   if( var != NULL && SCIPisGE(scip, SCIPvarGetUbGlobal(var), 1.0) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
      SCIP_CALL( SCIPsetSolVal( scip, sol, var, 1.0) );

   /* change the values of the edgevars */
   for( i = 0; i < nbins; ++i )
   {
      if( bin == i )
         continue;

      /** set the values for the old cluster */
      for( c = 0; c < ncluster; ++c )
      {
         if( NULL == edgevars[bin][i][1][0] )
            continue;
         if( c == oldcluster )
         {
            /* choose the right edgevariable w.r.t to the cluster/bin/i/oldcluster*/
            if( i < bin )
            {
               if( SCIPvarIsTransformed(edgevars[bin][i][oldcluster][oldcluster]) )
                  var = edgevars[bin][i][oldcluster][oldcluster];
               else
                  var = edgevars[bin][i][oldcluster][oldcluster];
            }
            else
            {
               if( SCIPvarIsTransformed(edgevars[i][bin][oldcluster][oldcluster]) )
                  var = edgevars[i][bin][oldcluster][oldcluster];
               else
                  var = edgevars[i][bin][oldcluster][oldcluster];
            }
         }
         else if( c < oldcluster )
         {
            if( SCIPvarIsTransformed(edgevars[bin][i][oldcluster][c]) )
               var = edgevars[bin][i][oldcluster][c];
            else
               var = edgevars[bin][i][oldcluster][c];
         }
         else
         {
            if( SCIPvarIsTransformed(edgevars[i][bin][c][oldcluster]) )
               var = edgevars[i][bin][c][oldcluster];
            else
               var = edgevars[i][bin][c][oldcluster];
         }
         /* set the solution value */
         if( var != NULL && SCIPisLE(scip, SCIPvarGetLbGlobal(var), 0) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
            SCIP_CALL( SCIPsetSolVal( scip, sol, var, 0.0  ) );
         /** set the values for the new cluster */
         if( c == newcluster )
         {
            /* choose the right edgevariable w.r.t to the cluster/bin/i/oldcluster*/
            if( i < bin )
            {
               if( SCIPvarIsTransformed(edgevars[bin][i][newcluster][newcluster]) )
                  var = edgevars[bin][i][newcluster][newcluster];
               else
                  var = edgevars[bin][i][newcluster][newcluster];
            }
            else
            {
               if( SCIPvarIsTransformed(edgevars[i][bin][newcluster][newcluster]) )
                  var = edgevars[i][bin][newcluster][newcluster];
               else
                  var = edgevars[i][bin][newcluster][newcluster];
            }
         }
         else if( c < newcluster )
         {
            if( SCIPvarIsTransformed(edgevars[bin][i][newcluster][c]) )
               var = edgevars[bin][i][newcluster][c];
            else
               var = edgevars[bin][i][newcluster][c];
         }
         else
         {
            if( SCIPvarIsTransformed(edgevars[i][bin][c][newcluster]) )
               var = edgevars[i][bin][c][newcluster];
            else
               var = edgevars[i][bin][c][newcluster];
         }
         if( var != NULL && SCIPisLE(scip, SCIPvarGetLbGlobal(var), 0) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
            SCIP_CALL( SCIPsetSolVal( scip, sol, var, clustering[i][c] ) );
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

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopySpaswitch)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurSpaswitch(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSpaswitch)
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
SCIP_DECL_HEUREXITSOL(heurExitsolSpaswitch)
{
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* reset the timing mask to its default value */
   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitSpaswitch)
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
SCIP_DECL_HEUREXEC(heurExecSpaswitch)
{  /*lint --e{715}*/

   SCIP_SOL* bestsol;                        /* incumbent solution */
   SCIP_SOL* worksol;                        /* working solution */
   SCIP_HEURDATA* heurdata;                  /* the heuristic data */
   SCIP_Bool feasible;                       /* True if a found solution is feasible */
   SCIP_Bool possible = TRUE;                /* True if it is possible to call the heuristic */

   SCIP_VAR** vars;                          /* SCIP variables */
   SCIP_VAR*** varmatrix;                    /* SCIP variables */
   SCIP_VAR***** edgevars;                   /* SCIP variables */

   SCIP_Real** solclustering;                /* the working cluster-assignment. We start with the one given by the solution */
   SCIP_Bool** binfixed;                     /* Is a bin fixed from scip */
   SCIP_Real** cmatrix;                      /* The transition matrix */
   SCIP_Real** qmatrix;                      /* The transition matrix projected down onto the current clustering */
   SCIP_Real* bincoherence;                  /* coherence influence of one bin on one cluster */
   SCIP_Real* mincoherence;                  /* minimal coherence influence for each cluster */
   SCIP_Real objective;                      /* value of the objective function */
   SCIP_Bool improvement = FALSE;            /* we switch bins until we can no longer find an imrpovement */
   SCIP_Real epsI;                           /* irreversibility */
   SCIP_Real qkl;
   SCIP_Real qlk;
   SCIP_Real neweps;
   int* clusterofbin;                        /* hold the cluster that each bin is in */
   int* minbin;
   int bin;
   int nbins;
   int ncluster;
   int c;
   int k,l;
   int i;
   int j;
   int bin1;
   int bin2;
   int nvars;
   int changecounter = 0;

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
   varmatrix = SCIPspaGetBinvars(scip);
   edgevars = SCIPspaGetEdgevars(scip);
   cmatrix = SCIPspaGetCmatrix(scip);

   /* we do not want to run the heurtistic if there is no 'flow' between the clusters.
    * in case of a (ideally) full reversible problem there cannot be a better solution, in the other case, i.e., the
    * problem has irreversible parts, it seems the heuristic will not find solutions respecting the coherence conditions
    */
   objective = SCIPgetSolOrigObj(scip, bestsol);
   if( SCIPisZero(scip, objective) )
      return SCIP_OKAY;

   /* create the working solution */
   SCIP_CALL( SCIPcreateSol(scip, &worksol, heur) );

   assert(nbins >= 0);
   assert(ncluster >= 0);
   assert(varmatrix != NULL);


   SCIP_CALL( SCIPcreateSol(scip, &worksol, heur) );
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real val;
      val = SCIPgetSolVal(scip, bestsol, vars[i]);
      SCIP_CALL( SCIPsetSolVal(scip, worksol, vars[i], val) );
   }


   /* allocate Memory */
   SCIPallocClearMemoryArray(scip, &bincoherence, nbins);
   SCIPallocClearMemoryArray(scip, &mincoherence, ncluster);
   SCIPallocClearMemoryArray(scip, &minbin, ncluster);
   SCIPallocClearMemoryArray(scip, &solclustering, nbins);
   SCIPallocClearMemoryArray(scip, &binfixed, nbins);
   SCIPallocClearMemoryArray(scip, &qmatrix, ncluster);
   SCIPallocClearMemoryArray(scip, &clusterofbin, nbins);

   /* Get the bin-variable values from the solution */
   for( i = 0; i < nbins; ++i )
   {
      SCIPallocClearMemoryArray(scip, &solclustering[i], ncluster);
      SCIPallocClearMemoryArray(scip, &binfixed[i], ncluster);

      for( c = 0; c < ncluster; ++c )
      {
         if( varmatrix[i][c] != NULL )
         {
            if( (SCIPisEQ(scip, SCIPvarGetUbGlobal(varmatrix[i][c]), SCIPvarGetLbGlobal(varmatrix[i][c]))) )
            {
               solclustering[i][c] = SCIPgetSolVal(scip, bestsol, varmatrix[i][c]);
               binfixed[i][c] = TRUE;

               if( SCIPisEQ(scip, solclustering[i][c], 1) )
                  clusterofbin[i] = c;
            }
            else
            {
               SCIP_Real solval = SCIPgetSolVal(scip, bestsol, varmatrix[i][c]);
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
                  while((varmatrix[j][k] == NULL || SCIPvarGetStatus(varmatrix[j][k]) == SCIP_VARSTATUS_FIXED || j == i))
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

   /* check if the assignment violates paritioning, e.g. because we are in a subscip */
   for( i = 0; i < nbins; ++i )
   {
      int amountzeros = 0;
      int sum = 0;
      for( j = 0; j < ncluster; ++j )
      {
         if( 0 == solclustering[i][j] )
            amountzeros++;
         if( 1 == solclustering[i][j] )
            sum++;
      }
      if( ncluster == amountzeros || sum > 1 )
         possible = FALSE;
   }

   /* initialize qmatrix */

   for( c = 0; c < ncluster; ++c )
   {
      SCIPallocMemoryArray(scip, &qmatrix[c], ncluster);
      mincoherence[c] = 1.0;
   }

   if( possible )
   {
      SCIP_Bool* checkedbins;

      /* allocate memory initialized to FALSE */
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &checkedbins, nbins) );

      computeIrrevMat(solclustering, qmatrix, cmatrix, nbins, ncluster);

      *result = SCIP_DIDNOTFIND;

      /* continue switching assignments until local optimum is reached */
      improvement = TRUE;
      while( improvement )
      {
         improvement = FALSE;

#ifndef NDEBUG
         for( i = 0; i < nbins; ++i )
            assert( SCIPisEQ(scip, solclustering[i][clusterofbin[i]], 1.0) );
#endif
         /* try to switch 1opt with any of the bins */

         for( bin = 0; bin < nbins; ++bin )
         {
            /* we do not want to switch a bin more than once */
            if( checkedbins[bin] )
               continue;

            for( k = 0; k < ncluster; ++k )
            {
               /* calculate the irreversibility after bin was moved from k to l */
               for( l = 0; l < ncluster; ++l )
               {
                  qkl = qmatrix[k][l];
                  qlk = qmatrix[l][k];

                  if( SCIPisZero(scip, solclustering[bin][k]) || binfixed[bin][k] || binfixed[bin][l] )
                     continue;

                  if( k == l )
                     continue;

                  /* do not open new clusters for now as the coherence would never be satisfied */
                  if( 1 == mincoherence[l] )
                     continue;

                  for( i = 0; i < nbins; ++i )
                  {
                     qkl -= solclustering[i][l] * cmatrix[bin][i];
                     qlk -= solclustering[i][l] * cmatrix[i][bin];

                     solclustering[bin][k] = 0;
                     solclustering[bin][l] = 1;

                     qlk += solclustering[i][k] * cmatrix[bin][i];
                     qkl += solclustering[i][k] * cmatrix[i][bin];

                     solclustering[bin][k] = 1;
                     solclustering[bin][l] = 0;

                  }
                  neweps = REALABS(qkl - qlk);

                  /* if we found an improvement, update the datastructures and try to add the solution to scip
                   * todo: add a paramter for minimal improvement
                   */
                  if( SCIPisGT(scip, neweps, 1.005*objective) && SCIPisEQ(scip, SCIPvarGetUbGlobal(varmatrix[bin][l]), 1) && SCIPisEQ(scip, SCIPvarGetLbGlobal(varmatrix[bin][k]), 0) )
                  {
                     setBinToCluster(solclustering, cmatrix, qmatrix, bin, k, FALSE, nbins, ncluster);
                     setBinToCluster(solclustering, cmatrix, qmatrix, bin, l, TRUE, nbins, ncluster);
                     clusterofbin[bin] = l;

                     SCIP_CALL( updateSolAssignment(scip, worksol, solclustering, nbins, ncluster, qmatrix, bin, k, l) );
                     SCIP_CALL( updateSolAuxiliary(scip, worksol, solclustering, nbins, ncluster, qmatrix) );

                     assert(isPartition(solclustering, nbins, ncluster));

                     epsI = getIrrevBound(scip, qmatrix, ncluster);
                     if( !SCIPisEQ(scip, neweps, epsI) )
                     {
                        SCIPdebugMessage("in 1opt: neweps: %f, epsi: %f \n", neweps, epsI);
                     }
                     SCIPtrySol(scip, worksol, FALSE, TRUE, TRUE, TRUE, &feasible);

                     if( feasible )
                     {
                        improvement = TRUE;
                        *result = SCIP_FOUNDSOL;
                        changecounter++;

                        /* mark the bin as switched */
                        checkedbins[bin] = TRUE;
                        objective = epsI;
                     }
                  }
               }
            }
         }

#ifndef NDEBUG
         for( i = 0; i < nbins; ++i )
            assert( SCIPisEQ(scip, solclustering[i][clusterofbin[i]], 1.0) );
#endif

         /* reset check array */
         for( bin = 0; bin < nbins; bin++ )
            checkedbins[bin] = FALSE;

         /* do a 2-opt i.e. trade any two bins between two clusters until local */
         for( bin1 = 0; bin1 < nbins; ++bin1 )
         {
            /* we do not want to switch a bin more than once */
            if( checkedbins[bin1] )
               continue;

            for( bin2 = 0; bin2 < nbins; ++bin2 )
            {
               /* we do not want to switch a bin more than once */
               if( checkedbins[bin2] )
                  continue;

               k = clusterofbin[bin1];
               l = clusterofbin[bin2];

               /* if the bins are the same or in the same cluster, do nothing */
               if( bin1 == bin2 || k == l || binfixed[bin1][k] || binfixed[bin1][l] || binfixed[bin2][k] || binfixed[bin2][l] )
                  continue;

               assert( SCIPisEQ(scip, solclustering[bin1][k], 1.0) && SCIPisEQ(scip, solclustering[bin2][l], 1.0) );

               qkl = qmatrix[k][l];
               qlk = qmatrix[l][k];

               /* calculate irreversibility after bin1 traded clusters with bin2 */
               for( i = 0; i < nbins; ++i )
               {
                  if( solclustering[i][l] == -1 || solclustering[i][k] == -1)
                     continue;
                  qkl -= solclustering[i][l] * cmatrix[bin1][i];
                  qlk -= solclustering[i][l] * cmatrix[i][bin1];

                  qkl -= solclustering[i][k] * cmatrix[i][bin2];
                  qlk -= solclustering[i][k] * cmatrix[bin2][i];

                  solclustering[bin1][k] = 0;
                  solclustering[bin1][l] = 1;
                  solclustering[bin2][l] = 0;
                  solclustering[bin2][k] = 1;

                  qlk += solclustering[i][k] * cmatrix[bin1][i];
                  qkl += solclustering[i][k] * cmatrix[i][bin1];

                  qlk += solclustering[i][l] * cmatrix[i][bin2];
                  qkl += solclustering[i][l] * cmatrix[bin2][i];

                  solclustering[bin1][k] = 1;
                  solclustering[bin1][l] = 0;
                  solclustering[bin2][l] = 1;
                  solclustering[bin2][k] = 0;
               }

               neweps = fabs(qkl - qlk);

               /* if we found an improvement, update the datastructures and try to add the solution to scip */
               if( SCIPisGT(scip, neweps, 1.005*objective) )
               {
                  setBinToCluster(solclustering, cmatrix, qmatrix, bin1, k, FALSE, nbins, ncluster);
                  setBinToCluster(solclustering, cmatrix, qmatrix, bin2, l, FALSE, nbins, ncluster);
                  setBinToCluster(solclustering, cmatrix, qmatrix, bin1, l, TRUE, nbins, ncluster);
                  setBinToCluster(solclustering, cmatrix, qmatrix, bin2, k, TRUE, nbins, ncluster);
                  clusterofbin[bin1] = l;
                  clusterofbin[bin2] = k;

                  SCIP_CALL( updateSolAssignment(scip, worksol, solclustering, nbins, ncluster, qmatrix, bin1, k, l) );
                  SCIP_CALL( updateSolAssignment(scip, worksol, solclustering, nbins, ncluster, qmatrix, bin2, l, k) );
                  SCIP_CALL( updateSolAuxiliary(scip, worksol, solclustering, nbins, ncluster, qmatrix) );

                  assert(isPartition(solclustering, nbins, ncluster));
                  epsI = getIrrevBound(scip, qmatrix, ncluster);
                  if( !SCIPisEQ(scip, neweps, epsI) )
                  {
                     SCIPdebugMessage("in 2opt: neweps: %f, epsi: %f \n", neweps, epsI);
                  }

                  SCIPtrySol(scip, worksol, FALSE, TRUE, TRUE, TRUE, &feasible);

                  if( feasible )
                  {
                     improvement = TRUE;
                     *result = SCIP_FOUNDSOL;
                     changecounter++;

                     /* mark the bin as switched */
                     checkedbins[bin1] = TRUE;
                     checkedbins[bin2] = TRUE;
                     objective = epsI;
                  }
               }
            }
         }
      }

      /* free memory */
      SCIPfreeMemoryArray(scip, &checkedbins);
   }

   /* update the found solution, so that we do not run again immediatly */

   heurdata->lastsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));

   /* free all data-structures */
   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeMemoryArray(scip, &solclustering[i]);
      SCIPfreeMemoryArray(scip, &binfixed[i]);
   }

   for( c = 0; c < ncluster; ++c )
   {
      SCIPfreeMemoryArray(scip, &qmatrix[c]);
   }

   SCIP_CALL( SCIPfreeSol(scip, &worksol) );

   SCIPfreeMemoryArray(scip, &qmatrix);
   SCIPfreeMemoryArray(scip, &bincoherence);
   SCIPfreeMemoryArray(scip, &mincoherence);
   SCIPfreeMemoryArray(scip, &minbin);
   SCIPfreeMemoryArray(scip, &solclustering);
   SCIPfreeMemoryArray(scip, &binfixed);
   SCIPfreeMemoryArray(scip, &clusterofbin);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the oneopt primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSpaswitch(
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
      HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecSpaswitch, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopySpaswitch) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeSpaswitch) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolSpaswitch) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitSpaswitch) );

   return SCIP_OKAY;
}
