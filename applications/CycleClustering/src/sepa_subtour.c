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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file   sepa_subtour.c
 * @brief  If there exists a transition forward along the cycle, then the state that the transition originates from can be reached only after
 * another ncluster - 1 transitions. Therefore cycles with a number of transitions smaller than that can be separated.
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "sepa_subtour.h"

#include "probdata_cyc.h"
#include "scip/cons_linear.h"
#include "scip/pub_misc.h"

#define SEPA_NAME              "subtour"
#define SEPA_DESC              "separator that elininates subtours of length smaller than |NCluster| in cycle-clusterign application"
#define SEPA_PRIORITY              1000
#define SEPA_FREQ                     1
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE      /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE      /**< should separation method be delayed, if other separators found cuts? */
#define MAXCUTS                     500
#define MAXROUNDS                    15

#ifdef SCIP_DEBUG
/** Print a cycle to the command line. For debugging purposes */
static
void printCycle(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  cycle,              /**< The cycle to be printed */
   int                   cyclelength,        /**< The length of the cycle */
   int                   nbins               /**< The number of states */
   )
{
   int i;

   SCIPinfoMessage(scip, NULL, "cycle_l%d_c: %d", cyclelength, cycle[0]);x
   for( i = 0; i < cyclelength; ++i )
   {
      SCIPinfoMessage(scip, NULL, " -> %d", cycle[i+1]);
   }
   SCIPinfoMessage(scip, NULL, "\n");
}
#endif

/** After finding a violation, construct and add all violated subtour cuts to scip */
static
SCIP_RETCODE addSubtourCuts(
   SCIP*                 scip,               /**< SCIP data structure. */
   SCIP_SEPA*            sepa,               /**< The subtour separator */
   SCIP_DIGRAPH*         capgraph,           /**< The directed edge-graph */
   SCIP_Real***          adjmatrices,        /**< The matrizes adjacency-matrix with the weight of all paths with 1,...,|Clutster| arcs */
   int                   cyclelength,        /**< The length of the subtours to add */
   SCIP_RESULT*          result,             /**< Pointer to store the result of separation */
   int*                  ncuts,              /**< Pointer to store number of cuts */
   int                   start               /**< The starting state */
   )
{
   SCIP_VAR**** edgevars;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_ROW* cut;
   int* cycle;
   SCIP_Bool* processed;
   SCIP_Bool doubleloop = FALSE;
   SCIP_Bool nullvars = FALSE;
   int nbins;
   int successor;
   int currnode;
   int k;
   int l;

   edgevars = SCIPcycGetEdgevars(scip);
   nbins = SCIPdigraphGetNNodes(capgraph);

   assert( SCIPisGT(scip, adjmatrices[cyclelength - 1][start][start], cyclelength - 1.0) );

   SCIP_CALL( SCIPallocClearMemoryArray(scip, &processed, nbins) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cycle, cyclelength + 1) );
   BMSclearMemoryArray(cycle, cyclelength + 1);

   cycle[0] = start;
   processed[start] = TRUE;

   /* iterate throguh all bins in the cycle */
   for( k = 0; k < cyclelength - 1; ++k )
   {
      currnode = cycle[k];
      if( currnode > nbins || currnode < 0 )
      {
         doubleloop = TRUE;
         break;
      }

      /* reconstruct the successor of the current bin from the adjacency matrices*/
      for( l = 0; l < SCIPdigraphGetNSuccessors(capgraph, cycle[k]); ++l )
      {
         successor = SCIPdigraphGetSuccessors(capgraph, currnode)[l];

         /* check if this successor of the current node is the one in the cycle. If so add it. */
         if( SCIPisFeasEQ(scip, adjmatrices[0][currnode][successor] + adjmatrices[cyclelength - (k + 2)][successor][start], adjmatrices[cyclelength - (k + 1)][currnode][start]) )
         {
            if( processed[successor] )
            {
               doubleloop = TRUE;
               break;
            }
            cycle[k + 1] = successor;
            processed[successor] = TRUE;
            break;
         }
      }
   }

   cycle[cyclelength] = start;

   if( !doubleloop )
   {
      /* construct the cut. Check for and add all necessary contractions */
      (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "subtour_%d_legnth_%d_contracted_%d", start, cyclelength );
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut,sepa, cutname, -SCIPinfinity(scip), cyclelength - 1.0, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

      for( k = 0; k < cyclelength; ++k )
      {
         if( NULL == edgevars[cycle[k]][cycle[k+1]] )
         {
            nullvars = TRUE;
            break;
         }

         SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[cycle[k]][cycle[k+1]][1], 1.0) );

         if( k > 0)
            SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[MAX(cycle[k],cycle[k+1])][MIN(cycle[k],cycle[k+1])][0], 1.0) );
      }

      SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
      SCIPdebug( printCycle(scip, cycle, cyclelength, nbins) );
      SCIPdebugMsg(scip, "Computed violation: %f \n", adjmatrices[cyclelength - 1][start][start] - cyclelength + 1);

      if( SCIPisCutEfficacious(scip, NULL, cut) && !nullvars )
      {
         SCIP_CALL( SCIPaddPoolCut(scip, cut) );
         *ncuts += 1;
         *result = SCIP_SEPARATED;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   }

   SCIPfreeBufferArray(scip, &cycle);
   SCIPfreeMemoryArray(scip, &processed);

   return SCIP_OKAY;
}

/** Detect if path inequalities are violated and if so, add them to scip */
static
SCIP_RETCODE addPathCuts(
   SCIP*                 scip,               /**< SCIP data structure. */
   SCIP_SEPA*            sepa,               /**< The subtour separator */
   SCIP_DIGRAPH*         capgraph,           /**< The directed edge-graph */
   SCIP_Real***          adjmatrices,        /**< The matrizes adjacency-matrix with the weight of all paths with 1,...,|Clutster| arcs */
   int                   pathlength,         /**< The length of the subtours to add */
   SCIP_RESULT*          result,             /**< Pointer to store the result of separation */
   int*                  ncuts,              /**< Pointer to store number of cuts */
   int                   start               /**< The start of the path */
   )
{
   SCIP_VAR**** edgevars;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_ROW* cut;
   int* path;
   SCIP_Bool nullvars = FALSE;
   SCIP_Real edgeweight;
   int currnode;
   int i;
   int j;
   int k;
   int nbins;
   int successor;

   edgevars = SCIPcycGetEdgevars(scip);
   nbins = SCIPdigraphGetNNodes(capgraph);

   SCIP_CALL( SCIPallocBufferArray(scip, &path, pathlength + 1) );

   path[0] =  start;

   for( j = 0; j < nbins; ++j )
   {
      if( j == start || edgevars[start][j] == NULL )
         continue;

      edgeweight = SCIPvarGetLPSol(edgevars[MAX(start,j)][MIN(start,j)][0]);

      if( edgeweight + adjmatrices[pathlength][start][j] > pathlength + 1 )
      {
         for( k = 0; k < pathlength; k++ )
         {
            currnode = path[k];

            for( i = 0; i < SCIPdigraphGetNSuccessors(capgraph, currnode); ++i )
            {
               successor = SCIPdigraphGetSuccessors(capgraph, currnode)[i];

               assert(successor < nbins && successor >= 0);

               if( SCIPisEQ(scip, adjmatrices[0][currnode][successor] + adjmatrices[pathlength - (k + 1)][successor][j], adjmatrices[pathlength - (k)][currnode][j]) )
               {
                  path[k + 1] = successor;
                  break;
               }
            }

            assert(path[k+1] >= 0 && path[k+1] < nbins);
         }

         path[pathlength + 1] = j;

         /* construct the cut. Check for and add all necessary contractions */
         (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "path_%d_legnth_%d_contracted_%d", start, pathlength );
         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut,sepa, cutname, -SCIPinfinity(scip), pathlength + 1.0, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

         for( k = 0; k <= pathlength; ++k )
         {
            if( NULL == edgevars[path[k]][path[k+1]] )
            {
               nullvars = TRUE;
               break;
            }
            SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[path[k]][path[k+1]][1], 1.0) );
            if( k > 0)
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[MAX(path[k],path[k+1])][MIN(path[k],path[k+1])][0], 1.0) );
         }

         SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[MAX(start,j)][MIN(start,j)][0], 1.0) );
         SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

         SCIPdebug( printCycle(scip, path, l, nbins) );

         if( SCIPisCutEfficacious(scip, NULL, cut) && !nullvars )
         {
            SCIP_CALL( SCIPaddPoolCut(scip, cut) );
            *ncuts += 1;
            *result = SCIP_SEPARATED;
         }

         SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         if( *ncuts >= MAXCUTS )
            goto TERMINATE;
      }
   }

TERMINATE:
   SCIPfreeBufferArray(scip, &path);

   return SCIP_OKAY;
}

/** Compute the next matrix with the weight off all the longest paths with fixed number of arcs. */
static
SCIP_Bool computeNextAdj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH*         capgraph,           /**< The directed edge-graph */
   SCIP_Real***          adjmatrices,        /**< The matrizes adjacency-matrix with the weight of all paths with 1,...,|Clutster| arcs */
   int                   cyclelength,        /**< The current number of arcs in the paths */
   int                   start               /**< The starting state */
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
   }

   /* Check if we have found a violated subtour constraint */
   if( SCIPisGT(scip, adjmatrices[cyclelength - 1][start][start], cyclelength - 1.0) )
      foundviolation = TRUE;

   return foundviolation;
}

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopySubtour)
{  /*lint --e{715}*/
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
{  /*lint --e{715}*/
   SCIP_VAR****    edgevars;
   SCIP_Real***    adjmatrices;
   SCIP_DIGRAPH*  capgraph;
   SCIP_Bool       violation;
   SCIP_Real  capacity;
   int ncuts;
   int nbins;
   int ncluster;
   int i;
   int j;
   int a;
   int k;
   int cyclelength;
   int rounds;

   rounds = SCIPsepaGetNCallsAtNode(sepa);
   ncluster = SCIPcycGetNCluster(scip);
   edgevars = SCIPcycGetEdgevars(scip);
   nbins = SCIPcycGetNBins(scip);
   ncuts = 0;

   if( rounds >= MAXROUNDS || ncluster < 4 )
   {
      *result =  SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   assert(nbins > 0);
   assert(ncluster > 0);
   assert(NULL != edgevars);

   /* allocate memory */
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &adjmatrices, ncluster) );
   for( i = 0; i < ncluster; ++i )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &adjmatrices[i], nbins) ); /*lint !e866*/
      for( j = 0; j < nbins; ++j )
      {
         SCIP_CALL( SCIPallocClearMemoryArray(scip, &adjmatrices[i][j], nbins) ); /*lint !e866*/
      }
   }

   /*Create Digraph from the current LP-Solution */
   for( a = 0; a < nbins; ++a )
   {
      SCIP_CALL( SCIPcreateDigraph(scip, &capgraph, nbins) );

      for( i = 0; i < nbins; ++i )
      {
         for( j = 0; j < nbins; ++j )
         {
            if( edgevars[i][j] != NULL && edgevars[i][j][1] != NULL )
            {
               capacity = SCIPvarGetLPSol(edgevars[i][j][1]);

               if ( i != a )
                  capacity += SCIPvarGetLPSol(edgevars[MAX(i,j)][MIN(i,j)][0]);

               if( SCIPisPositive(scip, capacity) )
               {
                  SCIP_CALL( SCIPdigraphAddArc(capgraph, i, j, &capacity) );
                  adjmatrices[0][i][j] = capacity;
               }
            }

            for( k = 1; k < ncluster; ++k )
            {
               adjmatrices[k][i][j] = 0;
            }
         }
      }

      /* a cyclelength of one does not make sense as there are no loops */
      cyclelength = 2;
      *result = SCIP_DIDNOTFIND;

      /* Iterate until we have found a sufficient number of cuts or until we have checked all possible violations */
      while( cyclelength < ncluster )
      {
         /* Compute the next adjacency matrix */
         violation = computeNextAdj(scip, capgraph, adjmatrices, cyclelength, a);

         /* If we found a violation separate it */
         if( violation && cyclelength != ncluster )
            SCIP_CALL( addSubtourCuts(scip, sepa, capgraph, adjmatrices, cyclelength, result, &ncuts, a) );
         if( cyclelength == ncluster - 1 )
            SCIP_CALL( addPathCuts(scip, sepa, capgraph, adjmatrices, cyclelength, result, &ncuts, a) );
         if( ncuts >= MAXCUTS )
            break;
         cyclelength++;
      }

      SCIPdigraphFreeComponents(capgraph);
      SCIPdigraphFree(&capgraph);
   }

   /* Free allocated memory */
   for( i = 0; i < ncluster; ++i )
   {
      for( j = 0; j < nbins; ++j )
      {
         SCIPfreeMemoryArray(scip, &adjmatrices[i][j]);
      }
      SCIPfreeMemoryArray(scip, &adjmatrices[i]);
   }
   SCIPfreeMemoryArray(scip, &adjmatrices);

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
