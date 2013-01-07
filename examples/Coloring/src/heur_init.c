/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_init.c
 * @brief  initial primal heuristic for the vertex coloring problem
 * @author Gerald Gamrath
 *
 * This file implements a heuristic which computes a starting solution for the coloring problem. It
 * therefore computes maximal stable sets and creates one variable for each set, which is added to the
 * LP.
 *
 * The heuristic is called only one time: before solving the root node.
 *
 * It checks, whether a solution-file was read in and a starting solution already exists.  If this
 * is not the case, an initial possible coloring is computed by a greedy method.  After that, a
 * tabu-search is called, which tries to reduce the number of colors needed. The tabu-search algorithm
 * follows the description in
 *
 * "A Survey of Local Search Methods for Graph Coloring"@n
 * by P. Galinier and A. Hertz@n
 * Computers & Operations Research, 33 (2006)
 *
 * The tabu-search works as follows: given the graph and a number of colors it tries to color the
 * nodes of the graph with at most the given number of colors.  It starts with a random coloring. In
 * each iteration, it counts the number of violated edges, that is, edges for which both incident
 * nodes have the same color. It now switches one node to another color in each iteration, taking
 * the node and color, that cause the greatest reduction of the number of violated edges, or if no
 * such combination exists, the node and color that cause the smallest increase of that number.  The
 * former color of the node is forbidden for a couple of iterations in order to give the possibility
 * to leave a local minimum.
 *
 * As long as the tabu-search finds a solution with the given number of colors, this number is reduced
 * by 1 and the tabu-search is called another time. If no coloring was found after a given number
 * of iterations, the tabu-search is stopped and variables for all sets of the last feasible coloring
 * are created and added to the LP (after possible extension to maximal stable sets).
 *
 * The variables of these sets result in a feasible starting solution of the coloring problem.
 *
 * The tabu-search can be deactivated by setting the parameter <heuristics/initcol/usetabu> to
 * FALSE.  The number of iterations after which the tabu-search stops if no solution was yet found
 * can be changed by the param <heuristics/initcol/maxiter>. A great effect is also obtained by
 * changing the parameters <heuristics/initcol/tabubase> and <heuristics/initcol/tabugamma>, which
 * determine the number of iterations for which the former color of a node is forbidden; more
 * precisely, this number is \<tabubase\> + ncritical * \<tabugamma\>, where ncritical is the number
 * of nodes, which are incident to violated edges.  Finally, the level of output and the frequency of
 * status lines can be changed by <heuristics/initcol/output> and <heuristics/initcol/dispfreq>.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_init.h"
#include "pricer_coloring.h"
#include "probdata_coloring.h"
#include "reader_col.h"
#include "scip/cons_setppc.h"
#include "cons_storeGraph.h"
#include "tclique/tclique.h"

#define HEUR_NAME             "initcol"
#define HEUR_DESC             "initial primal heuristic for coloring"
#define HEUR_DISPCHAR         't'
#define HEUR_PRIORITY         1
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/* default values for parameters */
#define DEFAULT_USETABU    TRUE
#define DEFAULT_MAXITER    100000
#define DEFAULT_TABUBASE   50
#define DEFAULT_TABUGAMMA  0.9
#define DEFAULT_OUTPUT     1
#define DEFAULT_DISPFREQ   10000



/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Bool usetabu;   /* should the tabu search heuristic be used in order to improve the greedy-solution? */
   int maxiter;         /* maximal number of iterations to be performed in each tabu-run */
   int tabubase;        /* constant part of the tabu-duration */
   SCIP_Real tabugamma; /* factor for the linear part of the tabu-duration */
   int output;          /* verbosity level for the output of the tabu search, 0: no output, 1: normal, 2: high */
   int dispfreq;        /* frequency for displaying status information, only active with output verbosity level 2 */
};




/*
 * Local methods
 */



/** checks whether one of the nodes has no color respectively has color -1 in the given array */
static
SCIP_Bool hasUncoloredNode(
   int                   nnodes,             /**< the graph that should be colored */
   int*                  colors              /**< array of ints representing the colors */
   )
{
   int i;

   assert(colors != NULL);

   for( i = 0; i < nnodes; i++)
   {
      /* node not yet colored */
      if(colors[i] == -1)
      {
	return TRUE;
      }
   }
   return FALSE;
}


/** computes a stable set with a greedy-method and colors its nodes */
static
SCIP_RETCODE greedyStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH*        graph,              /**< pointer to graph data structure */
   int*                  colors,             /**< array of ints representing the different colors, -1 means uncolored */
   int                   nextcolor           /**< color in which the stable set will be colored */
   )
{
   SCIP_Bool indNode;
   int nnodes;
   int i;
   int j;
   int* degrees;
   int* sortednodes;
   int* values;
   int* stablesetnodes;
   int nstablesetnodes;

   assert(graph != NULL);
   assert(colors != NULL);

   /* get number of nodes */
   nnodes = tcliqueGetNNodes(graph);
   nstablesetnodes = 0;

   /* get the  degrees and weights for the nodes in the graph */
   degrees = tcliqueGetDegrees(graph);
   SCIP_CALL( SCIPallocBufferArray(scip, &stablesetnodes, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &values, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortednodes, nnodes) );

   /* set values to the nodes which are used for sorting them */
   /* value = degree of the node + number of nodes if the node is yet uncolored,
      therefore the yet colored nodes have lower values than the not yet colored nodes */
   for( i = 0; i < nnodes; i++ )
   {
      sortednodes[i] = i;
      values[i] = degrees[i] + ( colors[i] == -1 ? nnodes : 0);
   }

   /* sort the nodes w.r.t. the computed values */
   SCIPsortDownIntInt(values, sortednodes, nnodes);

   /* insert first node */
   stablesetnodes[0] = sortednodes[0];
   nstablesetnodes = 1;
   for( i = 1; i < nnodes; i++)
   {
      if( colors[sortednodes[i]] != -1 )
      {
         break;
      }
      indNode = TRUE;
      for( j = 0; j < nstablesetnodes; j++ )
      {
         if( tcliqueIsEdge(graph, sortednodes[i], stablesetnodes[j]) )
         {
            indNode = FALSE;
            break;
         }
      }
      if( indNode == TRUE )
      {
         stablesetnodes[nstablesetnodes] = sortednodes[i];
         nstablesetnodes++;
      }

   }
   for( i = 0; i < nstablesetnodes; i++ )
   {
      assert(colors[stablesetnodes[i]] == -1);
      colors[stablesetnodes[i]] = nextcolor;
   }
   SCIPfreeBufferArray(scip, &stablesetnodes);
   SCIPfreeBufferArray(scip, &sortednodes);
   SCIPfreeBufferArray(scip, &values);

   return SCIP_OKAY;
}


static
/** computes the initial coloring with a greedy method */
SCIP_RETCODE greedyInitialColoring(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH*        graph,              /**< pointer to graph data structure */
   int*                  colors,             /**< array of ints representing the different colors */
   int*                  ncolors             /**< number of colors needed */
   )
{
   int nnodes;
   int i;
   int color;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(colors != NULL);

   nnodes = COLORprobGetNNodes(scip);
   assert(nnodes > 0);

   for( i = 0; i < nnodes; i++ )
   {
      colors[i] = -1;
   }

   color = 0;
   /* create stable sets until all Nodes are covered */
   while( hasUncoloredNode(nnodes, colors) )
   {
      greedyStableSet(scip, graph, colors, color);
      color++;
   }
   *ncolors = color;

   return SCIP_OKAY;

}


#ifndef NDEBUG
/** computes the number of violated edges, that means the number of edges (i,j) where i and j have the same color */
static
int getNViolatedEdges(
   TCLIQUE_GRAPH*        graph,              /**< the graph */
   int*                  colors              /**< colors of the nodes */
   )
{
   int nnodes;
   int i;
   int* j;
   int cnt;

   assert(graph != NULL);
   assert(colors != NULL);

   /* get the number of nodes */
   nnodes = tcliqueGetNNodes(graph);
   cnt = 0;

   /* count the number of violated edges, only consider edges (i,j) with i > j since the graph is undirected bu */
   for( i = 0; i < nnodes; i++ )
   {
      for( j = tcliqueGetFirstAdjedge(graph,i);  j <= tcliqueGetLastAdjedge(graph,i) && *j < i; j++ )
      {
         if( colors[i] == colors[*j] )
            cnt++;
      }
   }
   return cnt;
}
#endif


/** runs tabu coloring heuristic, gets a graph and a number of colors
 *  and tries to color the graph with at most that many colors;
 *  starts with a random coloring and switches one node to another color in each iteration,
 *  forbidding the old color for a couple of iterations
 */
static
SCIP_Bool runTabuCol(
   TCLIQUE_GRAPH*        graph,              /**< the graph, that should be colored */
   int                   seed,               /**< seed for the first random coloring */
   int                   maxcolors,          /**< number of colors, which are allowed */
   int*                  colors,             /**< output: the computed coloring */
   SCIP_HEURDATA*        heurdata            /**< data of the heuristic */
   )
{
   int nnodes;
   int** tabu;
   int** adj;
   int obj;
   int bestobj;
   int i;
   int j;
   int node1;
   int node2;
   int color1;
   int color2;
   int* firstedge;
   int* lastedge;
   SCIP_Bool restrictive;
   int iter;
   int minnode;
   int mincolor;
   int minvalue;
   int ncritical;
   SCIP_Bool aspiration;
   int d;
   int oldcolor;

   assert(graph != NULL);
   assert(heurdata != NULL);

   if( heurdata->output >= 1 )
      printf("Running tabu coloring with maxcolors = %d...\n", maxcolors);

   // get size
   nnodes = tcliqueGetNNodes(graph);

   srand( seed );

   // init random coloring
   for( i = 0; i < nnodes; i++ )
   {
      int rnd = rand();
      colors[i] = rnd % maxcolors;
      assert( 0 <= colors[i] && colors[i] < maxcolors );
   }

   // init matrices
   SCIP_CALL( SCIPallocMemoryArray(scip, &tabu, nnodes) );   // stores iteration at which tabu node/color pair will expire to be tabu
   SCIP_CALL( SCIPallocMemoryArray(scip, &adj, nnodes) );    // stores number of adjacent nodes using specified color

   for( i = 0; i < nnodes; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tabu[i]), maxcolors) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(adj[i]), maxcolors) );
      for( j = 0; j < maxcolors; j++ )
      {
	 tabu[i][j] = 0;
         adj[i][j] = 0;
      }
   }

   // objective
   obj = 0;

   // init adj-matrix and objective
   for( node1 = 0; node1 < nnodes; node1++ )
   {
      color1 = colors[node1];
      firstedge = tcliqueGetFirstAdjedge(graph, node1);
      lastedge = tcliqueGetLastAdjedge(graph, node1);
      while(  firstedge <= lastedge )
      {
         node2 = *firstedge;
	 color2 = colors[node2];
	 assert( 0 <= color2 && color2 < maxcolors );
	 (adj[node1][color2])++;
	 if( color1 == color2 )
	    obj++;
         firstedge++;
      }
   }
   assert( obj % 2 == 0 );
   obj = obj / 2;
   assert( obj == getNViolatedEdges(graph, colors) );

   bestobj = obj;
   restrictive = FALSE;
   iter = 0;
   if( obj > 0 )
   {
      // perform predefined number of iterations
      for( iter = 1; iter <= heurdata->maxiter; iter++ )
      {
	 // find best 1-move among those with critical vertex
	 minnode = -1;
	 mincolor = -1;
	 minvalue = nnodes * nnodes;
	 ncritical = 0;
	 for( node1 = 0; node1 < nnodes; node1++ )
	 {
	    aspiration = FALSE;
	    color1 = colors[node1];
	    assert( 0 <= color1 && color1 < maxcolors );

	    // if node is critical (has incident violated edges)
	    if( adj[node1][color1] > 0 )
	    {
	       ncritical++;
	       // check all colors
	       for( j = 0; j < maxcolors; j++ )
	       {
		  // if color is new
		  if( j != color1 )
		  {
		     // change in the number of violated edges:
		     d = adj[node1][j] - adj[node1][color1];

		     // 'aspiration criterion': stop if we get feasible solution
		     if( obj + d == 0 )
		     {
                        if( heurdata->output >= 1 )
                           printf("   Feasible solution found after %d iterations!\n\n", iter);
			minnode = node1;
			mincolor = j;
			minvalue = d;
			aspiration = TRUE;
			break;
		     }

		     // if not tabu and better value
		     if( tabu[node1][j] < iter &&  d < minvalue )
		     {
			minnode = node1;
			mincolor = j;
			minvalue = d;
		     }
		  }
	       }
	    }
	    if( aspiration )
	       break;
	 }

	 // if no candidate could be found - tabu list is too restrictive: just skip current iteration
	 if( minnode == -1 )
	 {
	    restrictive = TRUE;
	    continue;
	 }
	 assert( minnode != -1 );
	 assert( mincolor >= 0 );

	 // perform changes
	 assert( colors[minnode] != mincolor );
	 oldcolor = colors[minnode];
	 colors[minnode] = mincolor;
	 obj += minvalue;
	 assert( obj == getNViolatedEdges(graph, colors) );
	 if( obj < bestobj )
	    bestobj = obj;

	 if( heurdata->output == 2 && (iter) % (heurdata->dispfreq) == 0 )
	 {
            printf("Iter: %d  obj: %d  critical: %d   node: %d  color: %d  delta: %d\n",
                                        iter, obj, ncritical, minnode, mincolor, minvalue);
	 }

	 // terminate if valid coloring has been found
	 if( obj == 0 )
	    break;

	 // update tabu list
	 assert( tabu[minnode][oldcolor] < iter );
	 tabu[minnode][oldcolor] = iter + (heurdata->tabubase) + (int) (((double) ncritical) * (heurdata->tabugamma));

	 // update adj matrix
	 for( firstedge = tcliqueGetFirstAdjedge(graph, minnode); firstedge <= tcliqueGetLastAdjedge(graph, minnode); firstedge++ )
	 {
	    (adj[*firstedge][mincolor])++;
	    (adj[*firstedge][oldcolor])--;
	 }
      }
   }
   if( heurdata->output == 2 )
   {
      printf("Best objective: %d\n ", bestobj);
      if( restrictive )
      {
         printf("\nTabu list is probably too restrictive.\n");
      }
      printf("\n");
   }
   if( heurdata->output >= 1 && bestobj != 0 )
   {
      printf("   No feasible solution found after %d iterations!\n\n", iter-1);
   }

   for( i = 0; i < nnodes; i++ )
   {
      SCIPfreeMemoryArray(scip, &(tabu[i]));
      SCIPfreeMemoryArray(scip, &(adj[i]));
   }
   SCIPfreeMemoryArray(scip, &tabu);
   SCIPfreeMemoryArray(scip, &adj);

   // check whether valid coloring has been found
   if( obj == 0 )
      return TRUE;
   return FALSE;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyInit)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeInit)
{
   SCIP_HEURDATA* heurdata;

   /* free heuristic rule data */
   heurdata = SCIPheurGetData(heur);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecInit)
{
   int i;
   int j;
   int k;
   int nnodes;
   SCIP_SOL* sol;
   SCIP_Bool stored;
   SCIP_Bool success;
   SCIP_Bool indnode;
   int* colors;
   int* bestcolors;
   int ncolors;
   int nstablesetnodes;
   int setnumber;
   SCIP_VAR* var;
   SCIP_CONS** constraints;
   TCLIQUE_GRAPH* graph;
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTFIND;
   nnodes = COLORprobGetNNodes(scip);
   graph = COLORprobGetGraph(scip);
   /* create stable sets if no solution was read */
   if( COLORprobGetNStableSets(scip) == 0 )
   {
      /* get memory for arrays */
      SCIP_CALL( SCIPallocMemoryArray(scip, &colors, nnodes) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &bestcolors, nnodes) );

      /* get the node-constraits */
      constraints = COLORprobGetConstraints(scip);
      assert(constraints != NULL);

      /* compute an initial coloring with a greedy method */
      SCIP_CALL( greedyInitialColoring(scip, graph, bestcolors, &ncolors) );

      if( heurdata->usetabu )
      {
         /* try to find better colorings with tabu search method */
         success = TRUE;
         while( success )
         {
            ncolors--;
            success = runTabuCol(graph, 0, ncolors, colors, heurdata);
            if( success )
            {
               for( i = 0; i < nnodes; i++ )
               {
                  bestcolors[i] = colors[i];
               }
            }

         }
      }

      /* create vars for the computed coloring */
      for( i = 0; i <= ncolors; i++ )
      {
         /* save nodes with color i in the array colors and the number of such nodes in nstablesetnodes */
         nstablesetnodes = 0;
         for( j = 0; j < nnodes; j++ )
         {
            if( bestcolors[j] == i )
            {
               colors[nstablesetnodes] = j;
               nstablesetnodes++;
            }
         }
         /* try to add more nodes to the stable set without violationg the stability */
         for( j = 0; j < nnodes; j++ )
         {
            indnode = TRUE;
            for( k = 0; k < nstablesetnodes; k++ )
            {
               if( j == colors[k] || tcliqueIsEdge(graph, j, colors[k]) )
               {
                  indnode = FALSE;
                  break;
               }
            }
            if( indnode == TRUE )
            {
               colors[nstablesetnodes] = j;
               nstablesetnodes++;
            }
         }
         /* create variable for the stable set and add it to SCIP */
         SCIPsortDownInt(colors, nstablesetnodes);
         SCIP_CALL( COLORprobAddNewStableSet(scip, colors, nstablesetnodes, &setnumber) );
         assert(setnumber != -1);

         /* create variable for the stable set and add it to SCIP */
         SCIP_CALL( SCIPcreateVar(scip, &var, NULL, 0, 1, 1, SCIP_VARTYPE_BINARY,
               TRUE, TRUE, NULL, NULL, NULL, NULL, (SCIP_VARDATA*)(size_t)setnumber) );

         SCIP_CALL( COLORprobAddVarForStableSet(scip, setnumber, var) );
         SCIP_CALL( SCIPaddVar(scip, var) );
         SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );

         for( j = 0; j < nstablesetnodes; j++ )
         {
            /* add variable to node constraints of nodes in the set */
            SCIP_CALL( SCIPaddCoefSetppc(scip, constraints[colors[j]], var) );
         }


      }

      SCIPfreeMemoryArray(scip, &colors);
      SCIPfreeMemoryArray(scip, &bestcolors);

   }
   /* create solution consisting of all yet created stable sets,
      that means all sets of the solution given by the solution file or created by the greedy and tabu search */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   assert(sol != NULL);
   for( i = 0; i < COLORprobGetNStableSets(scip); i++ )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, COLORprobGetVarForStableSet(scip, i), 1.0) );
   }
   SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, FALSE, FALSE, FALSE, &stored) );
   assert(stored);

   /* set maximal number of variables to be priced in each round */
   SCIPsetIntParam(scip, "pricers/coloring/maxvarsround", COLORprobGetNStableSets(scip)*COLORprobGetNNodes(scip)/50);

   *result = SCIP_FOUNDSOL;

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the init primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurInit(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create init primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   heur = NULL;
   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecInit, heurdata) );
   assert(heur != NULL);

   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyInit) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeInit) );

   /* add parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/initcol/usetabu",
         "should the tabu search heuristic be used in order to improve the greedy-solution?",
         &heurdata->usetabu, FALSE, DEFAULT_USETABU, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/initcol/maxiter",
         "maximal number of iterations to be performed in each tabu-run",
         &heurdata->maxiter, TRUE, DEFAULT_MAXITER, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/initcol/tabubase",
         "constant part of the tabu-duration",
         &heurdata->tabubase, TRUE, DEFAULT_TABUBASE, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/initcol/tabugamma",
         "factor for the linear part of the tabu-duration",
         &heurdata->tabugamma, TRUE, DEFAULT_TABUGAMMA, -100, 100, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/initcol/output",
         "verbosity level for the output of the tabu search, 0: no output, 1: normal, 2: high",
         &heurdata->output, FALSE, DEFAULT_OUTPUT, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/initcol/dispfreq",
         "frequency for displaying status information, only active with output verbosity level 2",
         &heurdata->dispfreq, TRUE, DEFAULT_DISPFREQ, 0, INT_MAX, NULL, NULL) );


   return SCIP_OKAY;
}
