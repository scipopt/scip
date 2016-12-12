
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
/**@file   sepa_edge.c
 * @brief  edge-separator. Separates triangle-inequalities in SparseApprox Problem
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "sepa_subtour.h"
#include "probdata_spa.h"
#include "scip/cons_linear.h"

#define SEPA_NAME              "subtour"
#define SEPA_DESC              "separator that elininates subtours of length smaller than |NCluster| in cycle-clusterign application"
#define SEPA_PRIORITY         536870911
#define SEPA_FREQ                     1
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE      /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE      /**< should separation method be delayed, if other separators found cuts? */
#define MAXROUNDS                     5

/** Consider three bins i,j,k. If there is a way from i to j and j and k are int the same cluster, then i and k can be also considered adjacent.
 *  This method adds new arcs to the capacity-graph based on this observation.
 */
static
SCIP_RETCODE contractAdjGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH*         capgraph,           /**< The directed graph containing all the edges */
   SCIP_Real***          adjmatrices,        /**< The adjacency matrices to save the distances for paths with fixed number of arcs */
   int**                 connectmatrix       /**< This matrix saves, after a contraction, which node was contracted over.*/
)
{
   int i;
   int j;
   int k;
   int nbins;
   int succ;
   SCIP_Real contractweight;
   SCIP_Bool** added;
   SCIP_VAR**** edgevars;

   nbins = SCIPdigraphGetNNodes(capgraph);
   edgevars = SCIPspaGetEdgevars(scip);

   assert(nbins > 0);
   assert(edgevars != NULL);

   SCIP_CALL( SCIPallocClearMemoryArray(scip, &added, nbins ) );
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &added[i], nbins ) );
   }

   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < SCIPdigraphGetNSuccessors(capgraph, i); ++j )
      {
         succ = SCIPdigraphGetSuccessors(capgraph, i)[j];
         if( added[i][succ] )
            continue;
         for( k = 0; k < nbins; ++k )
         {
            /* If the weight of y_{ij} + y^0_{jk} - 1 > y_{ik} then it is better to use the way from i over j to k. */
            if( edgevars[MAX(succ, k)][MIN(succ, k)] == NULL || i == j || i == k || j == k || added[i][k] )
               continue;
            contractweight = adjmatrices[0][i][succ] - (1 - SCIPvarGetLPSol(edgevars[MAX(succ, k)][MIN(succ, k)][0]));
            if( SCIPisLT(scip, adjmatrices[0][i][k], contractweight) )
            {
               if( SCIPisZero(scip, adjmatrices[0][i][k]) )
               {
                  SCIPdigraphAddArc(capgraph, i, k, NULL);
                  added[i][k] =  TRUE;
               }
               adjmatrices[0][i][k] = contractweight;
               connectmatrix[i][k] = succ;
            }
         }
      }
   }
#ifndef SCIP_DEBUG
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < nbins; ++j )
      {
         if( added[i][j] )
            assert( 0 < SCIPvarGetLPSol(edgevars[MAX(connectmatrix[i][j], j)][MIN(connectmatrix[i][j], j)][0]));
      }
   }
#endif

   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeMemoryArray(scip, &added[i]);
   }
   SCIPfreeMemoryArray(scip, &added);
   return SCIP_OKAY;
}

/** After finding a violation, construct and add all violated subtour cuts to scip */
static
SCIP_RETCODE addSubtourCuts(
   SCIP*                 scip,               /**< SCIP data structure. */
   SCIP_SEPA*            sepa,               /**< The subtour separator */
   SCIP_DIGRAPH*         capgraph,           /**< The directed edge-graph */
   SCIP_Real***          adjmatrices,        /**< The matrizes adjacency-matrix with the weight of all paths with 1,...,|Clutster| arcs */
   int**                 connectmatrix,      /**< The matrix that has the intermediate nodes of the contractions */
   int                   cyclelength,        /**< The length of the subtours to add */
   SCIP_Bool*            result              /**< Pointer to store the result of separation */
)
{
   SCIP_VAR**** edgevars;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_Bool* processed;
   SCIP_ROW* cut;
   int cycle[cyclelength];
   int nbins;
   int successor;
   int currnode;
   int contractions;
   int i;
   int k;
   int l;

   edgevars = SCIPspaGetEdgevars(scip);
   nbins = SCIPdigraphGetNNodes(capgraph);

   SCIP_CALL( SCIPallocClearMemoryArray(scip, &processed, nbins) );
   for( i = 0; i < nbins; ++i )
   {
      if( processed[i] )
         continue;
      if( SCIPisFeasGT(scip, adjmatrices[cyclelength - 1][i][i], cyclelength - 1) )
      {
         cycle[0] = i;
         processed[i] = TRUE;
         contractions = 0;
         for( k = 0; k < cyclelength - 1; ++k )
         {
            for( l = 0; l < SCIPdigraphGetNSuccessors(capgraph, cycle[k]); ++l )
            {
               currnode = cycle[k];
               successor = SCIPdigraphGetSuccessors(capgraph, currnode)[l];
               if( SCIPisEQ(scip, adjmatrices[0][currnode][successor] + adjmatrices[cyclelength - (k + 2)][successor][i], adjmatrices[cyclelength - (k + 1)][currnode][i]) )
               {
                  cycle[k + 1] = successor;
                  processed[successor] = TRUE;
                  if( connectmatrix[cycle[k]][cycle[k+1]] >= 0 )
                     contractions++;
               }
            }
         }
         if( connectmatrix[cycle[cyclelength - 1]][cycle[0]] >= 0)
            contractions++;
         (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "subtour_%d_legnth_%d_contracted_%d", i, cyclelength, contractions );

         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut,sepa, cutname, -SCIPinfinity(scip), cyclelength - 1 + contractions, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

         for( k = 0; k < cyclelength - 1; ++k )
         {

            if( connectmatrix[cycle[k]][cycle[k+1]] >= 0)
            {
               SCIP_CALL(SCIPaddVarToRow(scip, cut, edgevars[MAX(connectmatrix[cycle[k]][cycle[k+1]],cycle[k+1])][MIN(connectmatrix[cycle[k]][cycle[k+1]],cycle[k+1])][0], 1.0) );
               SCIP_CALL(SCIPaddVarToRow(scip, cut, edgevars[cycle[k]][connectmatrix[cycle[k]][cycle[k+1]]][1], 1.0) );
            }
            else
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[cycle[k]][cycle[k+1]][1], 1.0) );
         }
         if( connectmatrix[cycle[cyclelength - 1]][cycle[0]] >= 0)
         {
            SCIP_CALL(SCIPaddVarToRow(scip, cut, edgevars[cycle[cyclelength - 1]][connectmatrix[cycle[cyclelength - 1]][cycle[0]]][1], 1.0) );
            SCIP_CALL(SCIPaddVarToRow(scip, cut, edgevars[MAX(connectmatrix[cycle[cyclelength - 1]][cycle[0]],cycle[0])][MIN(connectmatrix[cycle[cyclelength - 1]][cycle[0]],cycle[0])][0], 1.0) );
         }
         else
            SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[cycle[cyclelength - 1]][cycle[0]][1], 1.0) );
         SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
         SCIPdebug( SCIPprintRow(scip, cut, NULL) );
         SCIPdebugMsg(scip, "Computed violation: %f \n", adjmatrices[cyclelength - 1][i][i] - cyclelength + 1);
         if( SCIPisCutEfficacious(scip, NULL, cut) )
         {
            SCIP_CALL( SCIPaddPoolCut(scip, cut) );
            *result = SCIP_SEPARATED;
         }
         SCIP_CALL( SCIPreleaseRow(scip, &cut) );

      }

   }
   SCIPfreeMemoryArray(scip, &processed);
   return SCIP_OKAY;
}

/** Adds all violated path cuts to scip. A path cut is violated if there exists a directed arc y_{ij} and a directed path from i to j of length k: 1 < k < |Clusters|. */
static
SCIP_RETCODE addPathCuts(
   SCIP*                 scip,               /**< SCIP data structure. */
   SCIP_SEPA*            sepa,               /**< The subtour separator */
   SCIP_DIGRAPH*         capgraph,           /**< The directed edge-graph */
   SCIP_Real***          adjmatrices,        /**< The matrizes adjacency-matrix with the weight of all paths with 1,...,|Clutster| arcs */
   int**                 connectmatrix,      /**< The matrix that has the intermediate nodes of the contractions */
   int                   pathlength,         /**< The length of the path to separate */
   SCIP_Bool*            result              /**< Pointer to store the result of separation */
)
{
   SCIP_VAR**** edgevars;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_Bool* processed;
   SCIP_ROW* cut;
   SCIP_Real pathcap;
   SCIP_Real edgecap;
   int path[pathlength + 1];
   int nbins;
   int successor;
   int currnode;
   int contractions;
   int i;
   int j;
   int k;
   int l;

   edgevars = SCIPspaGetEdgevars(scip);
   nbins = SCIPdigraphGetNNodes(capgraph);

   SCIP_CALL( SCIPallocClearMemoryArray(scip, &processed, nbins) );
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < SCIPdigraphGetNSuccessors(capgraph, i); ++j )
      {
         successor = SCIPdigraphGetSuccessors(capgraph, i)[j];
         pathcap = adjmatrices[pathlength - 1][i][j];
         edgecap = adjmatrices[0][i][j];
         if(SCIPisGT(scip, pathcap + edgecap, pathlength) )
         {
            path[0] = i;
            contractions = 0;
            for( k = 0; k < pathlength - 1; ++k )
            {
               for( l = 0; l < SCIPdigraphGetNSuccessors(capgraph, path[k]); ++l )
               {
                  currnode = path[k];
                  successor = SCIPdigraphGetSuccessors(capgraph, currnode)[l];
                  if( SCIPisEQ(scip, adjmatrices[0][currnode][successor] + adjmatrices[pathlength - (k + 2)][successor][j], adjmatrices[pathlength - (k + 1)][currnode][j]) )
                  {
                     path[k + 1] = successor;
                     processed[successor] = TRUE;
                     if( connectmatrix[path[k]][path[k+1]] >= 0 )
                        contractions++;
                     break;
                  }
               }
            }
            path[pathlength] = j;
            if( connectmatrix[path[pathlength - 1]][j] >= 0 )
               contractions++;
            if( connectmatrix[i][j] > 0 )
               contractions++;
            (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "pathcut_%d_legnth_%d_contractions_%d", i, pathlength, contractions );
            SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut,sepa, cutname, -SCIPinfinity(scip), pathlength + contractions, FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

            for( k = 0; k < pathlength; ++k )
            {
               if( connectmatrix[path[k]][path[k+1]] >= 0)
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[path[k]][connectmatrix[path[k]][path[k+1]]][1], 1.0) );
                  SCIP_CALL(SCIPaddVarToRow(scip, cut, edgevars[MAX(connectmatrix[path[k]][path[k+1]],path[k+1])][MIN(connectmatrix[path[k]][path[k+1]],path[k+1])][0], 1.0) );
               }
               else
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[path[k]][path[k+1]][1], 1.0) );
            }
            if( connectmatrix[i][j] >= 0)
            {
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][connectmatrix[i][j]][1], 1.0) );
               SCIP_CALL(SCIPaddVarToRow(scip, cut, edgevars[MAX(connectmatrix[i][j],j)][MIN(connectmatrix[i][j],j)][0], 1.0) );
            }
            else
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][j][1], 1.0) );
            SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
            if( SCIPisCutEfficacious(scip, NULL, cut) )
            {
               SCIP_CALL( SCIPaddPoolCut(scip, cut) );
               *result = SCIP_SEPARATED;
            }

            SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         }
      }

   }
   SCIPfreeMemoryArray(scip, &processed);
   return SCIP_OKAY;
}

/** Compute the next matrix with the weight off all the longest paths with fixed number of arcs. */
static
SCIP_Bool computeNextAdj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH*         capgraph,           /**< The directed edge-graph */
   SCIP_Real***          adjmatrices,        /**< The matrizes adjacency-matrix with the weight of all paths with 1,...,|Clutster| arcs */
   int                   cyclelength         /**< The current number of arcs in the paths */
)
{
   int i;
   int j;
   int l;
   int successor;
   int nnodes;
   SCIP_Bool foundviolation;

   foundviolation = FALSE;
   nnodes = SCIPdigraphGetNNodes(capgraph);
   for( i = 0; i < nnodes; ++i )
   {
      for( j = 0; j < nnodes; ++j )
      {
         for( l = 0; l < SCIPdigraphGetNSuccessors(capgraph, i); ++l )
         {
            successor = SCIPdigraphGetSuccessors(capgraph, i)[l];
            if( SCIPisFeasPositive(scip, adjmatrices[0][i][successor]) && SCIPisFeasPositive(scip, adjmatrices[cyclelength - 2][successor][j]) )
            {
               if( SCIPisGT(scip, adjmatrices[0][i][successor] +  adjmatrices[cyclelength - 2][successor][j], adjmatrices[cyclelength - 1][i][j]) )
                  adjmatrices[cyclelength - 1][i][j] = adjmatrices[cyclelength - 2][successor][j] + adjmatrices[0][i][successor];
            }
         }
      }
      /* Check if we have found a violated subtour constraint */
      if( SCIPisGT(scip, adjmatrices[cyclelength - 1][i][i], cyclelength - 1) )
         foundviolation = TRUE;
   }
   return foundviolation;
}


/** copy method for separator plugins (called when SCIP copies plugins) */

static
SCIP_DECL_SEPACOPY(sepaCopySubtour)
{   /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaSubtour(scip) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpSubtour)
{
   SCIP_VAR****    edgevars;
   SCIP_Real***    adjmatrices;
   int**     connectmatrix;
   SCIP_DIGRAPH*  capgraph;
   SCIP_Bool       violation;
   SCIP_Real  capacity;
   int nbins;
   int ncluster;
   int i;
   int j;
   int cyclelength;
   int rounds;

   rounds = SCIPsepaGetNCallsAtNode(sepa);
   SCIPdebugMsg(scip, "Calls at node: %d \n", rounds);
   if( (rounds >= MAXROUNDS) )
   {
    *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   edgevars = SCIPspaGetEdgevars(scip);
   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);

   assert(nbins > 0);
   assert(ncluster > 0);
   assert(NULL != edgevars);

   SCIP_CALL( SCIPallocClearMemoryArray(scip, &adjmatrices, ncluster) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &connectmatrix, nbins) );
   for( i = 0; i < nbins; ++i )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &connectmatrix[i], nbins) );
      for( j = 0; j < nbins; ++j )
      {
         connectmatrix[i][j] = -1;
      }
   }
   for( i = 0; i < ncluster; ++i )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &adjmatrices[i], nbins) );
      for( j = 0; j < nbins; ++j )
      {
         SCIP_CALL( SCIPallocClearMemoryArray(scip, &adjmatrices[i][j], nbins) );
      }
   }

   /*Create Digraph from the current LP-Solution */
   SCIP_CALL( SCIPdigraphCreate(&capgraph, nbins) );
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < nbins; ++j )
      {
         if( edgevars[i][j] != NULL && edgevars[i][j][1] != NULL )
         {
            capacity = SCIPvarGetLPSol(edgevars[i][j][1]);
            if( SCIPisFeasPositive(scip, capacity) )
            {
               SCIP_CALL( SCIPdigraphAddArc(capgraph, i, j, &capacity) );
               adjmatrices[0][i][j] = capacity;
            }
         }
      }
   }
   /** Contract undirected edges with sufficiently large lp-value */
   SCIP_CALL( contractAdjGraph(scip, capgraph, adjmatrices, connectmatrix) );

   cyclelength = 2;
   *result = SCIP_DIDNOTFIND;
   /* Loop until we have found a violation or until we have checked all possible violations */
   while( cyclelength < ncluster && *result==SCIP_DIDNOTFIND )
   {
      /* Compute the next adjacency matrix */
      violation = computeNextAdj(scip, capgraph, adjmatrices, cyclelength);
      /* If we found a violation separate it */
      if( violation )
         SCIP_CALL( addSubtourCuts(scip, sepa, capgraph, adjmatrices, connectmatrix, cyclelength, result) );
      SCIP_CALL( addPathCuts(scip, sepa, capgraph, adjmatrices, connectmatrix, cyclelength, result) );
      cyclelength++;
   }


   /* Free allocated memory */
   for( i = 0; i < nbins; ++i )
   {
      SCIPfreeMemoryArray(scip, &connectmatrix[i]);
   }
   for( i = 0; i < ncluster; ++i )
   {
      for( j = 0; j < nbins; ++j )
      {
         SCIPfreeMemoryArray(scip, &adjmatrices[i][j]);
      }
      SCIPfreeMemoryArray(scip, &adjmatrices[i]);
   }
   SCIPfreeMemoryArray(scip, &adjmatrices);
   SCIPfreeMemoryArray(scip, &connectmatrix);
   SCIPdigraphFreeComponents(capgraph);
   SCIPdigraphFree(&capgraph);
   return SCIP_OKAY;
}


/** creates the Subtour separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaSubtour(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_SEPA* sepa;


   /* include separator */

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
      SEPA_USESSUBSCIP, SEPA_DELAY,
      sepaExeclpSubtour, NULL,
      NULL) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopySubtour) );


   return SCIP_OKAY;
}
