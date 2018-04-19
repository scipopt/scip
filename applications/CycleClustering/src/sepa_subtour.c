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
#define SEPA_FREQ                     5
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE      /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE      /**< should separation method be delayed, if other separators found cuts? */
#define MAXCUTS                    2000
#define MAXROUNDS                    15

#ifdef SCIP_DEBUG
/** Print a cycle to the command line. For debugging purposes */
static
void printCycle(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  cycle,              /**< The cycle to be printed */
   int                   cyclelength,        /**< The length of the cycle */
   int                   nstates               /**< The number of states */
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
   SCIP_Real***          adjacencymatrix,    /**< The matrizes adjacency-matrix with the weight of all paths with 1,...,|Clutster| arcs */
   SCIP_DIGRAPH*         adjacencygraph,     /**< The directed edge-graph */
   int**                 iscontracted,       /**< Information with intermediate contraction-nodes for contracted arcs */
   int                   cyclelength,        /**< The length of the subtours to add */
   SCIP_RESULT*          result,             /**< Pointer to store the result of separation */
   int*                  ncuts               /**< Pointer to store number of cuts */
   )
{
   SCIP_VAR**** edgevars;                    /* pointer for edge-variables from probdata */
   char cutname[SCIP_MAXSTRLEN];             /* cut-name */
   SCIP_ROW* cut;                            /* used to construct the cut */
   int** subtours;                           /* reconstruct the states in a subtour */
   int*  insubtour;                          /* save for each state if it is already in a different subtour */
   int nstates;                              /* store number of states from probdata */
   int currentnode;                          /* current node along the cycle */
   int successor;                            /* successor of current node */
   int intermediate;                         /* if arc is contracted, we use this for the intermediate node between current and successor */
   int anchor;                               /* state that we want to construct a violated subtour inequality from */
   int ncontractions;                        /* number of contracted arcs along the subtour */
   int liftabley;
   int liftablez;
   int c;
   int k;
   int l;
   SCIP_Bool isduplicate;                    /* is the subtour a duplicate of another one that has already been added */

   edgevars = SCIPcycGetEdgevars(scip);
   nstates = SCIPdigraphGetNNodes(adjacencygraph);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &insubtour, nstates) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &subtours, nstates) );

   for( k = 0; k < nstates; ++k )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &subtours[k], cyclelength + 1) ); /*lint !e866, !e776*/
      insubtour[k] = -1;
   }

   /* for each state, check if a subtour inequality is violated */
   for( anchor = 0; anchor < nstates; ++anchor )
   {
      /* while reconstructing the subtour, count the number of contractions */
      ncontractions = 0;
      /* a cycle inequality is violated if the following is true */
      if( SCIPisGT(scip, adjacencymatrix[cyclelength - 1][anchor][anchor], cyclelength - 1.0) )
      {
         subtours[anchor][0] = anchor;
         if( insubtour[anchor] == -1 )
            insubtour[anchor] = anchor;
         /*travere the cycle */
         for( k = 0; k < cyclelength -1; ++k )
         {
            currentnode = subtours[anchor][k];

            /* find the next state along the subtour */
            for( l = 0; l < SCIPdigraphGetNSuccessors(adjacencygraph, currentnode); l++ )
            {
               successor = SCIPdigraphGetSuccessors(adjacencygraph, currentnode)[l];
               /* check if this successor of the current node is the one in the cycle. If so add it. */
               if( SCIPisEQ(scip, adjacencymatrix[0][currentnode][successor] + adjacencymatrix[cyclelength - (k + 2)][successor][anchor],
                                  adjacencymatrix[cyclelength - (k + 1)][currentnode][anchor]) )
               {
                  subtours[anchor][k + 1] = successor;
                  insubtour[successor] = anchor;
                  if( iscontracted[currentnode][successor] != -1 )
                     ncontractions++;
                  break;
               }
            }
         }
         /* start and endnode are always the same in a cycle */
         subtours[anchor][cyclelength] = anchor;
         /* check last arc for a contraction */
         if( iscontracted[subtours[anchor][cyclelength - 1]][anchor] != -1 )
            ncontractions++;

         isduplicate = FALSE;
         /* if this anchor is already in another subtour, we check if the subtour is the same, since we don't want to add duplicates */
         if( insubtour[anchor] != anchor )
         {
            c = 0;
            isduplicate = TRUE;
            while( subtours[insubtour[anchor]][c] != anchor )
               c++;
            for( k = 0; k < cyclelength && isduplicate; ++k )
            {
               if( subtours[insubtour[anchor]][(k + c) % cyclelength] != subtours[anchor][k] )
                  isduplicate = FALSE;
            }
         }

         if( isduplicate )
            continue;

         liftabley = cyclelength - 1;
         liftablez = SCIPcycGetNCluster(scip) - cyclelength - 1;
         /* Now build the cut and add the subtour inequality */
         (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "subtour_%d_length_%d_contracted_%d", anchor, cyclelength, ncontractions );
         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut,sepa, cutname, -SCIPinfinity(scip), cyclelength + ncontractions - 1.0, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
         for( k = 0; k < cyclelength; ++k )
         {
            currentnode = subtours[anchor][k];
            successor = subtours[anchor][k+1];
            intermediate = iscontracted[currentnode][successor];

            if( intermediate != -1 )
            {
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[currentnode][intermediate][1], 1.0) );
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[MAX(intermediate, successor)][MIN(intermediate, successor)][0], 1.0) );
               if( liftabley > 0 && SCIPvarGetLPSol(edgevars[MAX(currentnode, intermediate)][MIN(currentnode, intermediate)][0]) > 0 )
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[MAX(currentnode, intermediate)][MIN(currentnode, intermediate)][0], 1.0) );
                  liftabley--;
               }
               if( liftablez > 0 && SCIPvarGetLPSol(edgevars[MAX(currentnode, intermediate)][MIN(currentnode, intermediate)][0]) > 0 )
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[intermediate][successor][1], 1.0) );
                  liftablez--;
               }
            }
            else
            {
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[currentnode][successor][1], 1.0) );
               if( liftabley > 0 && SCIPvarGetLPSol(edgevars[MAX(currentnode, successor)][MIN(currentnode, successor)][0]) > 0 )
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[MAX(currentnode, successor)][MIN(currentnode, successor)][0], 1.0) );
                  liftabley--;
               }
            }
         }
         SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
         SCIP_CALL( SCIPaddPoolCut(scip, cut) );

         /* Print for debugging purposes */
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, cut, NULL) ) );

         /* release data and increment cut counter */
         SCIP_CALL( SCIPreleaseRow(scip, &cut) );

         *result = SCIP_SEPARATED;
         (*ncuts)++;
      }
   }

   for( k = 0; k < nstates; ++k )
   {
      SCIPfreeBlockMemoryArray(scip, &(subtours[k]), cyclelength + 1);
   }
   SCIPfreeBlockMemoryArray(scip, &subtours, nstates);
   SCIPfreeBlockMemoryArray(scip, &insubtour, nstates);

   return SCIP_OKAY;
}

/** Detect if path inequalities are violated and if so, add them to scip */
static
SCIP_RETCODE addPathCuts(
   SCIP*                 scip,               /**< SCIP data structure. */
   SCIP_SEPA*            sepa,               /**< The subtour separator */
   SCIP_Real***          adjacencymatrix,    /**< The matrizes adjacency-matrix with the weight of all paths with 1,...,|Clutster| arcs */
   SCIP_DIGRAPH*         adjacencygraph,     /**< The directed edge-graph */
   int**                 iscontracted,       /**< Information with intermediate contraction-nodes for contracted arcs */
   int                   pathlength,         /**< The length of the subtours to add */
   SCIP_RESULT*          result,             /**< Pointer to store the result of separation */
   int*                  ncuts               /**< Pointer to store number of cuts */
   )
{
   SCIP_VAR**** edgevars;                    /* pointer for edge-variables from probdata */
   char cutname[SCIP_MAXSTRLEN];             /* cut-name */
   SCIP_ROW* cut;                            /* used to construct the cut */
   int* path;                                /* reconstruct the states in a violated path inequality */
   int nstates;                              /* store number of states from probdata */
   int currentnode;                          /* current node along the path */
   int successor;                            /* successor of current node */
   int intermediate;                         /* if arc is contracted, we use this for the intermediate node between current and successor */
   int start;                                /* starting-state that we want to construct a violated path inequality from */
   int end;                                  /* last state in path that we want to construct a violated path inequality from */
   int ncontractions;                        /* number of contracted arcs along the subtour */
   int k;
   int i;
   int j;
   int nz;
   int ny;

   edgevars = SCIPcycGetEdgevars(scip);
   nstates = SCIPdigraphGetNNodes(adjacencygraph);

   SCIP_CALL( SCIPallocMemoryArray(scip, &path, pathlength + 1) );

   for( start = 0; start < nstates; ++start )
   {
      path[0] =  start;
      for( j = 0; j < SCIPdigraphGetNSuccessors(adjacencygraph, start); ++j )
      {
         ncontractions = 0;

         end = SCIPdigraphGetSuccessors(adjacencygraph, start)[j];
         path[pathlength] = end;

         /* check if path-inequality is violated */
         if( SCIPisGT(scip, adjacencymatrix[pathlength - 1][start][end] + adjacencymatrix[0][start][end], (SCIP_Real) pathlength) )
         {
            /*reconstruct the path */
            for( k = 0; k < pathlength - 1; ++k )
            {
               currentnode = path[k];
               for( i = 0; i < SCIPdigraphGetNSuccessors(adjacencygraph, currentnode); ++i )
               {
                  successor = SCIPdigraphGetSuccessors(adjacencygraph, currentnode)[i];
                  if( SCIPisEQ(scip, adjacencymatrix[0][currentnode][successor] + adjacencymatrix[pathlength - (k + 2)][successor][end], adjacencymatrix[pathlength - (k + 1)][currentnode][end]) )
                  {
                     path[k + 1] = successor;
                     if( iscontracted[currentnode][successor] != -1 )
                        ncontractions++;
                     break;
                  }
               }
            }

            /* check the last arc along the path and the direct arc from start to end for contractions */
            if( iscontracted[path[pathlength - 1]][end] != -1 )
               ncontractions++;
            if( iscontracted[start][end] != -1 )
               ncontractions++;

            nz = pathlength;
            ny = 0;

            /*construct the corresponding inequality and add it to scip */
            (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "path_%d_%d_length_%d_contracted_%d", start, end, pathlength, ncontractions );
            SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut,sepa, cutname, -SCIPinfinity(scip), (SCIP_Real) pathlength + ncontractions, FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
            for( k = 0; k < pathlength; ++k )
            {
               currentnode = path[k];
               successor = path[k+1];
               intermediate = iscontracted[currentnode][successor];

               if( intermediate != -1 )
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[currentnode][intermediate][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[MAX(intermediate, successor)][MIN(intermediate, successor)][0], 1.0) );
                  if( nz < SCIPcycGetNCluster(scip) && SCIPisPositive(scip, SCIPvarGetLPSol(edgevars[intermediate][successor][1])) )
                  {
                     SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[intermediate][successor][1], 1.0) );
                     nz++;
                  }
                  if( ny < pathlength - 2 && SCIPisPositive(scip, SCIPvarGetLPSol(edgevars[MAX(currentnode, intermediate)][MIN(currentnode, intermediate)][0])) )
                  {
                     SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[MAX(currentnode, intermediate)][MIN(currentnode, intermediate)][0], 1.0) );
                     ny ++;
                  }
               }
               else
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[currentnode][successor][1], 1.0) );
                  if( ny < pathlength - 2 && SCIPisPositive(scip, SCIPvarGetLPSol(edgevars[MAX(currentnode, successor)][MIN(currentnode, successor)][0])) )
                  {
                     SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[MAX(currentnode, successor)][MIN(currentnode, successor)][0], 1.0) );
                     ny ++;
                  }
               }
            }

            /* add the direct arc from start to end */
            intermediate = iscontracted[start][end];
            if( iscontracted[start][end] != -1 )
            {
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[start][intermediate][1], 1.0) );
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[MAX(intermediate, end)][MIN(intermediate, end)][0], 1.0) );
            }
            else
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[start][end][1], 1.0) );

            SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

            /* print row if in debug mode */
            SCIPdebug( SCIPprintRow(scip, cut, NULL) );

            /* if an arc appears twice then the path inequality should not be used */
            if( SCIPgetRowMaxCoef(scip, cut) == 1 )
            {
               SCIP_CALL( SCIPaddPoolCut(scip, cut) );
               *result = SCIP_SEPARATED;
               (*ncuts)++;
            }
            SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         }
      }
   }

   SCIPfreeMemoryArray(scip, &path);

   return SCIP_OKAY;
}

/** Detect if path inequalities are violated and if so, add them to scip */
static
SCIP_RETCODE addTourCuts(
   SCIP*                 scip,               /**< SCIP data structure. */
   SCIP_SEPA*            sepa,               /**< The subtour separator */
   SCIP_Real***          adjacencymatrix,    /**< The matrizes adjacency-matrix with the weight of all paths with 1,...,|Clutster| arcs */
   SCIP_DIGRAPH*         adjacencygraph,     /**< The directed edge-graph */
   int**                 iscontracted,       /**< Information with intermediate contraction-nodes for contracted arcs */
   int                   tourlength,         /**< The length of the subtours to add */
   SCIP_RESULT*          result,             /**< Pointer to store the result of separation */
   int*                  ncuts               /**< Pointer to store number of cuts */
   )
{
   SCIP_VAR**** edgevars;                    /* pointer for edge-variables from probdata */
   char cutname[SCIP_MAXSTRLEN];             /* cut-name */
   SCIP_ROW* cut;                            /* used to construct the cut */
   int* tour;                                /* reconstruct the states in a violated tour inequality */
   int nstates;                              /* store number of states from probdata */
   int currentnode;                          /* current node along the tour */
   int successor;                            /* successor of current node */
   int intermediate;                         /* if arc is contracted, we use this for the intermediate node between current and successor */
   int start;                                /* starting-state that we want to construct a violated tour inequality from */
   int end;                                  /* last state in tour that we want to construct a violated tour inequality from */
   int ncontractions;                        /* number of contracted arcs along the subtour */
   int k;
   int i;
   int j;
   int nz;
   int ny;

   edgevars = SCIPcycGetEdgevars(scip);
   nstates = SCIPdigraphGetNNodes(adjacencygraph);

   SCIP_CALL( SCIPallocMemoryArray(scip, &tour, tourlength + 1) );

   for( start = 0; start < nstates; ++start )
   {
      tour[0] =  start;
      for( j = 0; j < SCIPdigraphGetNSuccessors(adjacencygraph, start); ++j )
      {
         ncontractions = 0;

         end = SCIPdigraphGetSuccessors(adjacencygraph, start)[j];
         tour[tourlength] = end;

         /* check if tour-inequality is violated */
         if( SCIPisGT(scip, adjacencymatrix[tourlength - 1][start][end] - adjacencymatrix[0][end][start], (SCIP_Real) tourlength - 1) )
         {
            /*reconstruct the tour */
            for( k = 0; k < tourlength - 1; ++k )
            {
               currentnode = tour[k];
               for( i = 0; i < SCIPdigraphGetNSuccessors(adjacencygraph, currentnode); ++i )
               {
                  successor = SCIPdigraphGetSuccessors(adjacencygraph, currentnode)[i];
                  if( SCIPisEQ(scip, adjacencymatrix[0][currentnode][successor] + adjacencymatrix[tourlength - (k + 2)][successor][end], adjacencymatrix[tourlength - (k + 1)][currentnode][end]) )
                  {
                     tour[k + 1] = successor;
                     if( iscontracted[currentnode][successor] != -1 )
                        ncontractions++;
                     break;
                  }
               }
            }

            /* check the last arc along the tour and the direct arc from start to end for contractions */
            if( iscontracted[tour[tourlength - 1]][end] != -1 )
               ncontractions++;
            if( iscontracted[end][start] != -1 )
               ncontractions++;

            nz = tourlength;
            ny = 0;

            /*construct the corresponding inequality and add it to scip */
            (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "tour_%d_%d_length_%d_contracted_%d", start, end, tourlength, ncontractions );
            SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut,sepa, cutname, -SCIPinfinity(scip), (SCIP_Real) tourlength + ncontractions - 1, FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
            for( k = 0; k < tourlength; ++k )
            {
               currentnode = tour[k];
               successor = tour[k+1];
               intermediate = iscontracted[currentnode][successor];

               if( intermediate != -1 )
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[currentnode][intermediate][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[MAX(intermediate, successor)][MIN(intermediate, successor)][0], 1.0) );
               }
               else
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[currentnode][successor][1], 1.0) );
               }
            }

            /* add the direct arc from start to end */
            intermediate = iscontracted[end][start];
            if( iscontracted[end][start] != -1 )
            {
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[end][intermediate][1], -1.0) );
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[MAX(intermediate, start)][MIN(intermediate, start)][0], 1.0) );
            }
            else
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[end][start][1], -1.0) );

            SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

            /* print row if in debug mode */
            SCIPdebug( SCIPprintRow(scip, cut, NULL) );

            /* if an arc appears twice then the tour inequality should not be used */
            if( SCIPgetRowMaxCoef(scip, cut) == 1 )
            {
               SCIP_CALL( SCIPaddPoolCut(scip, cut) );
               *result = SCIP_SEPARATED;
               (*ncuts)++;
            }
            SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         }
      }
   }

   SCIPfreeMemoryArray(scip, &tour);

   return SCIP_OKAY;
}



/** Compute the next matrix with the weight off all the longest paths with exactly narcs and store it in adjacencymatrix[narcs - 1].
 *  For this, simply compute 	\f{align*}{ d^{k}(currentnode,successor) = max_{l=1,\ldots,n} \{d^{k-1}(currentnode,l) + d^1(l,successor) \}	\f}.
 */
static
SCIP_Bool computeNextAdjacency
(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real***          adjacencymatrix,    /**< The max-distance matrices for all number of arcs less than narcs. Store next matrix here as well */
   SCIP_DIGRAPH*         adjacencygraph,     /**< The directed edge-graph */
   int                   narcs               /**< The current number of arcs in the paths */
)
{
   int currentnode;
   int intermediate;
   int successor;
   int l;
   int nnodes;
   SCIP_Bool foundviolation;

   foundviolation = FALSE;
   nnodes = SCIPdigraphGetNNodes(adjacencygraph);

   for( currentnode = 0; currentnode < nnodes; ++currentnode )
   {
      for( l = 0; l < SCIPdigraphGetNSuccessors(adjacencygraph, currentnode); ++l )
      {
         intermediate = SCIPdigraphGetSuccessors(adjacencygraph, currentnode)[l];
         for( successor = 0; successor < nnodes; ++successor )
         {
            if( SCIPisPositive(scip, adjacencymatrix[0][currentnode][intermediate]) && SCIPisPositive(scip, adjacencymatrix[narcs - 2][intermediate][successor]) )
            {
               if( SCIPisGT(scip, adjacencymatrix[0][currentnode][intermediate] +  adjacencymatrix[narcs - 2][intermediate][successor], adjacencymatrix[narcs - 1][currentnode][successor]) )
                  adjacencymatrix[narcs - 1][currentnode][successor] = adjacencymatrix[0][currentnode][intermediate] + adjacencymatrix[narcs - 2][intermediate][successor];
            }
         }
      }
   }

   /* Check if we have found a violated subtour constraint */
   for( currentnode = 0; currentnode < nnodes; ++currentnode )
   {
      if( SCIPisGT(scip, adjacencymatrix[narcs - 1][currentnode][currentnode], narcs - 1.0) )
         foundviolation = TRUE;
   }
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
   SCIP_VAR**** edgevars;                    /* the edgevars that we get from the probdata */
   SCIP_Real*** adjacencymatrix;             /* stores the length of max-distance paths between all nodes with fixed number arcs */
   SCIP_DIGRAPH* adjacencygraph;             /* stores all non-zero-length arcs w.r.t. current LP-solution */
   SCIP_DIGRAPH* edgegraph;                  /* pointer that we get from probdata, tells us which variables are not null */
   int** iscontracted;                       /* if it is better to take z_{ij} + y_{jk} rather than z_{ik} directly, this stores j at position (i,k) */
   SCIP_Bool violation;                      /* has a violated subtour ineq been found ? */
   int ncuts;                                /* how many cuts have been created */
   int nstates;                              /* number of states, get from probdata */
   int ncluster;                             /* number of clusters, get from probdata */
   int cyclelength;
   int rounds;
   int i;
   int j;
   int k;
   int state1;
   int state2;
   int state3;

   /* get problem information */
   rounds = SCIPsepaGetNCallsAtNode(sepa);
   ncluster = SCIPcycGetNCluster(scip);
   edgevars = SCIPcycGetEdgevars(scip);
   nstates = SCIPcycGetNBins(scip);
   edgegraph = SCIPcycGetEdgeGraph(scip);
   ncuts = 0;

   if( rounds >= MAXROUNDS )
   {
      *result =  SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   assert(nstates > 0);
   assert(ncluster > 0);
   assert(NULL != edgevars);
   assert(NULL != edgegraph);

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &adjacencymatrix, ncluster) );
   for( k = 0; k < ncluster; ++k )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &adjacencymatrix[k], nstates) ); /*lint !e866*/
      for( j = 0; j < nstates; ++j )
      {
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &adjacencymatrix[k][j], nstates) ); /*lint !e866*/
      }
   }

   /*Create Digraph from the current LP-Solution */
   SCIP_CALL( SCIPcreateDigraph(scip, &adjacencygraph, nstates) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &iscontracted, nstates) );


   /* get the values of the lp-solution */
   for( i = 0; i < nstates; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &iscontracted[i], nstates) );

      for( j = 0; j < nstates; ++j )
      {
         iscontracted[i][j] = -1;

         if( edgevars[i][j] != NULL && edgevars[i][j][1] != NULL )
            adjacencymatrix[0][i][j] = SCIPvarGetLPSol(edgevars[i][j][1]);
      }
   }

   /* contract the adjacency matrix if it is better to take z_{ij} + y_{jk} rather than z_{ik} directly, this stores j at position (i,k) */
   for( i = 0; i < nstates; ++i )
   {
      state1 = i;
      for( j = 0; j < SCIPdigraphGetNSuccessors(edgegraph, i); ++j )
      {
         state2 = SCIPdigraphGetSuccessors(edgegraph, i)[j];
         for( k = 0 ; k < SCIPdigraphGetNSuccessors(edgegraph, state2); ++k )
         {
            state3 = SCIPdigraphGetSuccessors(edgegraph, state2)[k];
            if( edgevars[state1][state2] == NULL || edgevars[state2][state3] == NULL || edgevars[state1][state3] == NULL )
               continue;
            if( SCIPisLT( scip, adjacencymatrix[0][state1][state3] ,SCIPvarGetLPSol(edgevars[state1][state2][1]) + SCIPvarGetLPSol(edgevars[MAX(state2, state3)][MIN(state2, state3)][0]) - 1) )
            {
               adjacencymatrix[0][state1][state3] = SCIPvarGetLPSol(edgevars[state1][state2][1]) + SCIPvarGetLPSol(edgevars[MAX(state2, state3)][MIN(state2, state3)][0]) - 1;
               iscontracted[state1][state3] = state2;
            }
         }
      }
   }

   /* save the contracted matrix as a digraph to be able to reuse it quicker */
   for( i = 0; i < nstates; ++i )
   {
      for( j = 0; j < nstates; ++j )
      {
         if( !SCIPisZero(scip, adjacencymatrix[0][i][j]) )
            SCIP_CALL( SCIPdigraphAddArc(adjacencygraph, i , j, NULL) );
      }
   }

   /* a cyclelength of one does not make sense as there are no loops */
   cyclelength = 2;
   *result = SCIP_DIDNOTFIND;

   /* Iterate until we have found a sufficient number of cuts or until we have checked all possible violations */
   while( cyclelength < ncluster )
   {
      /* Compute the next adjacency matrix */
      violation = computeNextAdjacency(scip, adjacencymatrix, adjacencygraph, cyclelength);

      /* If we found a violation separate it */
      if( violation )
         SCIP_CALL( addSubtourCuts(scip, sepa, adjacencymatrix, adjacencygraph, iscontracted, cyclelength, result, &ncuts) );

      /* check if any path-inequalities are violated and sepatare them */
      SCIP_CALL( addPathCuts(scip, sepa, adjacencymatrix, adjacencygraph, iscontracted, cyclelength, result, &ncuts) );

      if( cyclelength == ncluster - 1 )
         SCIP_CALL( addTourCuts(scip, sepa, adjacencymatrix, adjacencygraph, iscontracted, cyclelength, result, &ncuts) );

      /* stop if we added maximal number of cuts */
      if( ncuts >= MAXCUTS )
         break;

      cyclelength++;
   }

   SCIPdigraphFreeComponents(adjacencygraph);
   SCIPdigraphFree(&adjacencygraph);

   /* Free allocated memory */
   for( i = 0; i < nstates; ++i )
   {
      SCIPfreeBlockMemoryArray(scip, &iscontracted[i], nstates);
   }
   SCIPfreeBlockMemoryArray(scip, &iscontracted, nstates);

   for( i = 0; i < ncluster; ++i )
   {
      for( j = 0; j < nstates; ++j )
      {
         SCIPfreeBlockMemoryArray(scip, &adjacencymatrix[i][j], nstates);
      }
      SCIPfreeBlockMemoryArray(scip, &adjacencymatrix[i], nstates);
   }
   SCIPfreeBlockMemoryArray(scip, &adjacencymatrix, ncluster);

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
