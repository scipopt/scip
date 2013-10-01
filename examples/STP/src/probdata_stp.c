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

/**@file   probdata_stp.c
 * @brief  Problem data for stp problem
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Michael Winkler
 * @author Daniel Rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "probdata_stp.h"
#include <stdio.h>
#include "scip/scip.h"
#include "cons_stp.h"
#include "grph.h"
#include "portab.h"
#include "scip/cons_linear.h"
#include "scip/misc.h"
#include "scip/struct_misc.h"

#define   DEFAULT_COMPCENTRAL  1
#define   DEFAULT_EMITGRAPH    FALSE
#define   DEFAULT_REDUCTION    4

#define CENTER_OK    0
#define CENTER_DEG   1
#define CENTER_SUM   2
#define CENTER_MIN   3
#define CENTER_ALL   4

#define MODE_CUT    0
#define MODE_FLOW   1
#define MODE_PRICE  2



/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the graph of the steiner tree problem
 */
struct SCIP_ProbData
{
   int                   mode;               /**< solving mode selected by the user (Cut, Price, Flow) */
   SCIP_Bool             bigt;               /**< stores whether the 'T' model is being used (not relevant in the cut mode) */

   GRAPH*                graph;        	     /**< the graph */
   SCIP_CONS**           edgecons;           /**< array of constraints */
   SCIP_CONS**           pathcons;           /**< array of constraints */
   SCIP_VAR** 		 edgevars;	     /**< array of edge variables */
   SCIP_VAR**            flowvars;           /**< array of edge variables (needed only in the Flow mode) */
   int* 	         realterms;          /**< array of all terminals except the root */
   SCIP_Bool             emitgraph;          /**< emitgraph */
   int                   nedges;             /**< number of edges */
   int                   nterms;             /**< number of terminals */
   int                   realnterms;         /**< number of terminals except the root */
   int                   nlayers;            /**< number of layers */
   int                   nnodes;             /**< number of nodes */

   double*               xval;               /**< values of the edge variables */
   int                   nvars;              /**< number of variables */
   SCIP_Longint          lastlpiters;        /**< Branch and Cut */

};

/**@name Local methods
 *
 * @{
 */

/* what = CENTER_OK  : Do nothing
 *      = CENTER_DEG : find maximum degree
 *      = CENTER_SUM : find the minimum distance sum
 *      = CENTER_MIN : find the minimum largest distance
 *      = CENTER_ALL : find the minimum distance sum to all knots
 */
static
int central_terminal(
   GRAPH* g,
   int    what)
{
   PATH*   path;
   double* cost;
   int     i;
   int     k;
   int     center  = -1;
   int     degree  = 0;
   double  sum;
   double  max;
   double  minimum = FARAWAY;
   double  maximum = 0.0;
   double  oldval  = 0.0;

   assert(g         != NULL);
   assert(g->layers == 1);

   if (what == CENTER_OK)
      return g->source[0];

   /* Find knot with maximum degree.
    */
   if (what == CENTER_DEG)
   {
      degree = 0;

      for(i = 0; i < g->knots; i++)
      {
         if (Is_term(g->term[i]) && (g->grad[i] > degree))
         {
            degree = g->grad[i];
            center = i;
         }
      }
      assert(degree > 0);

      return center;
   }

   /* For the other medthods we need the shortest paths.
    */
   path = malloc((size_t)g->knots * sizeof(*path));
   cost = malloc((size_t)g->edges * sizeof(*cost));

   assert(path != NULL);
   assert(cost != NULL);

   for(i = 0; i < g->knots; i++)
      g->mark[i] = TRUE;

   for(i = 0; i < g->edges; i++)
      cost[i] = 1.0;

   for(i = 0; i < g->knots; i++)
   {
      if (!Is_term(g->term[i]))
         continue;

      graph_path_exec(g, FSP_MODE, i, cost, path);

      sum = 0.0;
      max = 0.0;

      for(k = 0; k < g->knots; k++)
      {
         assert((path[k].edge >= 0) || (k == i));
         assert((path[k].edge >= 0) || (path[k].dist == 0));

         if (Is_term(g->term[k]) || (what == CENTER_ALL))
         {
            sum += path[k].dist;

            if (path[k].dist > max)
               max = path[k].dist;
         }
      }

      if ((what == CENTER_SUM) || (what == CENTER_ALL))
      {
         if (sum < minimum)
         {
            minimum = sum;
            center  = i;
         }
         if (sum > maximum)
            maximum = sum;

         if (i == g->source[0])
            oldval = sum;
      }
      else
      {
         assert(what == CENTER_MIN);

         /* If the maximum distance to terminal ist shorter or if
          * it is of the same length but the degree of the knot is
          * higher, we change the center.
          */
         if (LT(max, minimum) || (EQ(max, minimum) && (g->grad[i] > degree)))
         {
            minimum = max;
            center  = i;
            degree  = g->grad[i];
         }
         if (max > maximum)
            maximum = max;

         if (i == g->source[0])
            oldval = max;
      }
   }
   assert(center >= 0);
   assert(Is_term(g->term[center]));

   free(cost);
   free(path);

   printf("Central Terminal is %d (min=%g, max=%g, old=%g)\n",
      center, minimum, maximum, oldval);

   return center;
}

/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   GRAPH*                graph               /**< graph */
   )
{
   assert(scip != NULL);
   assert(probdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, probdata) );

   (*probdata)->graph = graph;
   (*probdata)->xval = NULL;
   (*probdata)->lastlpiters = -1;

   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
   int e;
   int t;
   SCIPdebugPrintf ("probdataFREE \n");
   assert(scip != NULL);
   assert(probdata != NULL);

   /* release variables */
   for( e = 0; e < (*probdata)->nvars; ++e )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->edgevars[e]) );
   }
   SCIPfreeMemoryArrayNull(scip, &(*probdata)->edgevars);

   /* release path constraints */
   if( (*probdata)->mode == MODE_PRICE )
   {
      for( t = 0; t < (*probdata)->realnterms; ++t)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &((*probdata)->pathcons[t])) );
      }
   }
   else if( (*probdata)->mode == MODE_FLOW )
   {
      for( t = 0; t < (*probdata)->realnterms; ++t)
      {
	 /* release constraints and variables */
         for( e = 0; e < (*probdata)->nnodes - 1; ++e )
         {
            SCIP_CALL( SCIPreleaseCons(scip, &((*probdata)->pathcons[t * ((*probdata)->nnodes - 1 ) + e])) );
         }
         for( e = 0; e < (*probdata)->nedges; ++e )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->flowvars[t * (*probdata)->nedges + e]) );
	 }
         if( (*probdata)->bigt )
	    break;
      }
      SCIPfreeMemoryArrayNull(scip, &(*probdata)->flowvars);
   }

   /* release edge constraints (Price or Flow) */
   if( (*probdata)->mode != MODE_CUT)
   {
      for( t = 0; t < (*probdata)->realnterms; ++t)
      {
         for( e = 0; e < (*probdata)->nedges; ++e )
         {
            SCIP_CALL( SCIPreleaseCons(scip, &((*probdata)->edgecons[t * (*probdata)->nedges + e])) );
         }
	 if( (*probdata)->bigt )
	    break;

      }
      SCIPfreeMemoryArrayNull(scip, &((*probdata)->edgecons));
      SCIPfreeMemoryArrayNull(scip, &((*probdata)->pathcons));
   }

   SCIPfreeMemoryArrayNull(scip, &(*probdata)->xval);
   SCIPfreeMemoryArrayNull(scip, &(*probdata)->realterms);

   /* free probdata */
   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

/** print graph (in undirected form) in GML format */
static
SCIP_RETCODE probdataPrintGraph(
   GRAPH*                graph,              /**< Graph to be printed */
   const char*           filename,           /**< Name of the output file */
   SCIP_Bool*            edgemark            /**< Array of (undirected) edges to highlight */
   )
{
   char label[SCIP_MAXSTRLEN];
   FILE* file;
   int e;
   int n;
   int m;

   assert(graph != NULL);
   file = fopen((filename != NULL) ? filename : "graphX.gml", "w");

   for( e = 0; e < graph->edges; e += 2 )
   {
      assert(graph->tail[e] == graph->head[e + 1]);
      assert(graph->tail[e + 1] == graph->head[e]);
   }

   /* write GML format opening, undirected */
   SCIPgmlWriteOpening(file, FALSE);

   /* write all nodes, discriminate between root, terminals and the other nodes */
   e = 0;
   m = 0;
   for( n = 0; n < graph->knots; ++n )
   {
      if( n == graph->source[0] )
      {
         (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "Root");
	 SCIPgmlWriteNode(file, (unsigned int)n, label, "rectangle", "#666666", NULL);
	 m = 1;
      }
      else if( graph->term[n] == 0 )
      {
	 (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "Terminal %d", e + 1);
	 SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#ff0000", NULL);
	 e += 1;
      }
      else
      {
         (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "Node %d", n + 1 - e - m);
         SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#336699", NULL);
      }
   }

   /* write all edges (undirected) */
   for( e = 0; e < graph->edges; e += 2 )
   {
      (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%8.2f", graph->cost[e]);
      if( edgemark != NULL && edgemark[e / 2] == TRUE )
	 SCIPgmlWriteEdge(file, (unsigned int)graph->tail[e], (unsigned int)graph->head[e], label, "#ff0000");
      else
         SCIPgmlWriteEdge(file, (unsigned int)graph->tail[e], (unsigned int)graph->head[e], label, NULL);
   }

   /* write GML format closing */
   SCIPgmlWriteClosing(file);

   return SCIP_OKAY;
}

/** create constraints */
static
SCIP_RETCODE createConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   GRAPH* graph;
   char consname[SCIP_MAXSTRLEN];
   int t;
   int e;
   int k;
   int k2;
   int realnterms;
   int nedges;
   int nnodes;

   assert(scip != NULL);
   assert(probdata != NULL);

   SCIPdebugPrintf("createConstraints \n");
   graph = probdata->graph;
   nedges = probdata->nedges;
   nnodes = probdata->nnodes;
   realnterms = probdata->realnterms;

   /* create edge constraints (used by both Flow and Price) */
   if( !probdata->bigt )
   {
      /* create |T \ {root}|*|E| edge constraints (disaggregated) */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->edgecons), (realnterms) * (nedges)) );
      for( t = 0; t < realnterms; ++t )
      {
         for( e = 0; e < nedges; ++e )
         {
	    (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "EdgeConstraint%d_%d", t, e);
            SCIP_CALL( SCIPcreateConsLinear ( scip, &( probdata->edgecons[t * nedges + e] ), consname, 0, NULL, NULL,
	          -SCIPinfinity(scip), 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, probdata->mode == MODE_PRICE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, probdata->edgecons[t * nedges + e]) );
         }
      }
   }
   else
   {
      /* create |E| edge constraints ('T') */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->edgecons), nedges) );
      for( e = 0; e < nedges; ++e )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "EdgeConstraintT%d", e);
         SCIP_CALL( SCIPcreateConsLinear ( scip, &( probdata->edgecons[e] ), consname,
               0, NULL, NULL, -SCIPinfinity(scip), 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, probdata->mode == MODE_PRICE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, probdata->edgecons[e]) );
      }
   }

   /* Branch and Price mode */
   if( probdata->mode == MODE_PRICE )
   {
      /* create |T \ {root}| path constraints */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->pathcons), realnterms) );

      for( t = 0; t < realnterms; ++t )
      {
	 (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PathConstraint%d", t);
         SCIP_CALL( SCIPcreateConsLinear ( scip, &(probdata->pathcons[t]), consname, 0, NULL, NULL, 1.0, SCIPinfinity(scip),
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, probdata->pathcons[t]) );
      }
   }
   /* Flow mode */
   else if( probdata->mode == MODE_FLOW )
   {
      /* create path constraints */
      if( !probdata->bigt)
      {
         /* not in 'T' mode, so create |T \ {root} |*|V \ {root}| path constraints  */
         SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->pathcons), realnterms * (nnodes - 1)) );
         for( t = 0; t < realnterms; ++t )
         {
            k2 = 0;
            for( k = 0; k < nnodes; ++k )
            {
               /* if node k is not the root */
               if( k !=  graph->source[0])
               {
                  /* if node k is not the t-th terminal, set RHS = 0 */
                  if( k != probdata->realterms[t] )
                  {
                     (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PathConstraint%d_%d", t, k2);
                     SCIP_CALL( SCIPcreateConsLinear(scip, &( probdata->pathcons[t * (nnodes - 1) + k2] ), consname,
                           0, NULL, NULL, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
                     SCIP_CALL( SCIPaddCons(scip, probdata->pathcons[t * (nnodes - 1) + k2]) );
                  }

                  /* if node k is the t-th terminal, set RHS = 1 */
                  else
                  {
                     (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PathConstraint%d_%d", t, k2);
                     SCIP_CALL( SCIPcreateConsLinear(scip, &( probdata->pathcons[t * (nnodes - 1) + k2] ), consname,
                           0, NULL, NULL, 1.0, 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
                     SCIP_CALL( SCIPaddCons(scip, probdata->pathcons[t * (nnodes - 1) + k2]) );
                  }

                  k2 += 1;
               }
            }
         }
      }
      else
      {
         /* in 'T' mode, so create |V \ {root}| path constraints */
	 SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->pathcons), (nnodes - 1)) );
	 k2 = 0;
	 for( k = 0; k < nnodes; ++k )
         {
            /* if node k is not the root */
            if( k != graph->source[0])
	    {
	       /* if node k is not a terminal, set RHS = 0 */
	       if( graph->term[k] != 0 )
	       {
                  (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PathConstraintT%d", k2 + 1);
                  SCIP_CALL( SCIPcreateConsLinear(scip, &( probdata->pathcons[k2] ), consname,
                        0, NULL, NULL, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
                  SCIP_CALL( SCIPaddCons(scip, probdata->pathcons[k2]) );
	       }
	       /* if node k is a terminal, set RHS = 1 */
	       else
	       {
                  (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PathConstraintT%d", k2 + 1);
                  SCIP_CALL( SCIPcreateConsLinear(scip, &( probdata->pathcons[k2] ), consname,
                        0, NULL, NULL, 1.0, 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
                  SCIP_CALL( SCIPaddCons(scip, probdata->pathcons[k2]) );
	       }
	       k2 += 1;
	    }
         }
      }
   }

   return SCIP_OKAY;
}

/** create initial columns */
static
SCIP_RETCODE createVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Real             offset              /**< offset computed during the presolving */
   )
{
   GRAPH* graph;
   SCIP_VAR* offsetvar;
   char varname[SCIP_MAXSTRLEN];
   PATH* path;
   SCIP_VAR* var;
   double* edgecost;
   SCIP_Bool objint;
   int nnodes;
   int nedges;
   int e;
   int t;
   int k;
   int nvars;
   int nlayers;
   int k2;
   int tail;
   int root;
   int nflows;
   int realnterms;
   assert(scip != NULL);
   assert(probdata != NULL);

   graph = probdata->graph;
   SCIPdebugPrintf("createVariables \n");

   /* if the graph reduction solved the whole problem, NULL is returned */
   if( graph != NULL )
   {
      nedges = probdata->nedges;
      nvars = probdata->nvars;
      nlayers = probdata->nlayers;
      realnterms = probdata->realnterms;
      root = graph->source[0];
      nnodes = graph->knots;
      objint = SCIPisIntegral(scip, offset);

      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->xval, nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->edgevars, nvars) );

      /* Cut mode */
      if( probdata->mode == MODE_CUT )
      {
         for( k = 0; k < nlayers; ++k )
         {
            for( e = 0; e < nedges; ++e )
            {
               (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d_%d_%d", graph->tail[e], graph->head[e], k);
               SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->edgevars[k * nedges + e], varname, 0.0, 1.0, graph->cost[e], SCIP_VARTYPE_BINARY) );
               SCIP_CALL( SCIPaddVar(scip, probdata->edgevars[k * nedges + e]) );
            }
         }
      }
      /* Price or Flow mode */
      else
      {
         /* create and add the edge variables */
         for( e = 0; e < nedges; ++e )
         {
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "EdgeVar%d_%d", graph->tail[e], graph->head[e]);
            var = NULL;
            SCIP_CALL( SCIPcreateVarBasic(scip, &var, varname, 0.0, 1.0, graph->cost[e], SCIP_VARTYPE_BINARY) );
            objint = objint && SCIPisIntegral(scip, graph->cost[e]);
            SCIP_CALL( SCIPaddVar(scip, var) );
            SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );
            probdata->edgevars[e] = var;

            if ( !probdata->bigt )
	    {
               for( t = 0; t < realnterms; ++t )
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[t * nedges + e], var, -1.0) );
               }
	    }
	    else
	    {
	       SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[e], var, -realnterms) );
	    }
         }
      }
      /* Price mode */
      if( probdata->mode == MODE_PRICE )
      {

         /* the flow variables are not used in the Price mode */
	 probdata->flowvars = NULL;

         /* compute shortest paths to the root */
         SCIP_CALL( SCIPallocMemoryArray(scip, &path, graph->knots) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &edgecost, nedges) );
         for( e = 0; e < nedges; ++e )
            edgecost[e] = graph->cost[e];

	 for( e = 0; e < graph->knots; e++ )
            graph->mark[e] = 1;

         graph_path_exec(graph, FSP_MODE, root, edgecost, path);

         /* create and add initial path variables (one for each real terminal) */
         for( t = 0; t < realnterms; ++t )
         {
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "PathVar%d_0", t);
            var = NULL;
            SCIP_CALL( SCIPcreateVarBasic(scip, &var, varname, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, var) );
            SCIP_CALL( SCIPaddCoefLinear(scip, probdata->pathcons[t], var, 1.0) );
            tail = probdata->realterms[t];
            while( tail != root )
            {
	       if( !probdata->bigt )
                  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[t * nedges + path[tail].edge], var, 1.0) );
	       else
                  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[path[tail].edge], var, 1.0) );

               tail=graph->tail[path[tail].edge];
            }
         }

         SCIPfreeMemoryArray(scip, &edgecost);
         SCIPfreeMemoryArray(scip, &path);
      }
      /* Flow mode */
      else if( probdata->mode == MODE_FLOW )
      {
	 /* store the number of disparate flows (or commodities) in nflows */
	 if ( !probdata->bigt )
	    nflows = realnterms;
	 else
	    nflows = 1;

         SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->flowvars, nflows * nedges) );

         /* create and add the flow variables */
         for( e = 0; e < nedges; ++e )
         {
            for( t = 0; t < nflows; ++t )
            {
               var = NULL;
               (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "FlowVar%d.%d_%d", t, graph->tail[e], graph->head[e]);
               SCIP_CALL( SCIPcreateVarBasic(scip, &var, varname, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
               SCIP_CALL( SCIPaddVar(scip, var) );
               probdata->flowvars[t * nedges + e] = var;
               SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[t * nedges + e], probdata->flowvars[t * nedges + e], 1.0) );
            }
         }

         /* add the flow variables to the corresponding path constraints */
         for( t = 0; t < nflows; ++t )
         {
            k2 = 0;
            for( k = 0; k < nnodes; ++k )
            {
               if( k != root )
               {
                  e = graph->inpbeg[k];
                  while( e >= 0 )
                  {
                     SCIP_CALL( SCIPaddCoefLinear(scip, probdata->pathcons[t * (nnodes - 1)  + k2], probdata->flowvars[t * nedges + e], 1.0) );
                     e = graph->ieat[e];
                  }
                  e = graph->outbeg[k];
                  while( e >= 0 )
                  {
                     SCIP_CALL( SCIPaddCoefLinear(scip, probdata->pathcons[t * (nnodes - 1)  + k2 ], probdata->flowvars[t * nedges + e], -1.0) );
                     e = graph->oeat[e];
                  }
                  k2 += 1;
               }
            }
         }
      }

      /* if all edge costs and the offset are integral, tell SCIP the objective will be integral */
      if( objint )
         SCIP_CALL( SCIPsetObjIntegral(scip) );
   }
   else
   {
      probdata->edgevars = NULL;
      probdata->flowvars = NULL;
   }

   /* add offset */
   (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "OFFSET");
   SCIP_CALL( SCIPcreateVarBasic(scip, &offsetvar, varname, 1.0, 1.0, offset, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, offsetvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &offsetvar) );

   return SCIP_OKAY;
}


/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigStp)
{
   SCIPdebugPrintf("probdelorigStp \n");

   SCIPdebugMessage("free original problem data\n");

   if ( (*probdata)->mode == MODE_CUT )
      graph_mincut_exit();

   graph_path_exit();

   if( (*probdata)->graph != NULL )
      graph_free((*probdata)->graph);

   SCIP_CALL( probdataFree(scip, probdata) );
   SCIPdebugPrintf("probdelorigStpout \n");
   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransStp)
{
   int i;
   SCIPdebugPrintf("probtransStp \n");

   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->graph) );

   (*targetdata)->nlayers = sourcedata->nlayers;
   (*targetdata)->nedges = sourcedata->nedges;
   (*targetdata)->nnodes = sourcedata->nnodes;
   (*targetdata)->nterms = sourcedata->nterms;
   (*targetdata)->realnterms = sourcedata->realnterms;
   (*targetdata)->emitgraph = sourcedata->emitgraph;
   (*targetdata)->nvars = sourcedata->nvars;
   (*targetdata)->mode = sourcedata->mode;
   (*targetdata)->bigt = sourcedata->bigt;

   if( sourcedata->nedges > 0 )
   {
      /* transform variables */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->xval, sourcedata->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->edgevars, sourcedata->nvars) );
      SCIP_CALL( SCIPtransformVars(scip, sourcedata->nvars, sourcedata->edgevars, (*targetdata)->edgevars) );

      /* Cut mode */
      if( sourcedata->mode == MODE_CUT )
      {
         (*targetdata)->edgecons = NULL;
         (*targetdata)->pathcons = NULL;
      }
      /* Price or Flow mode */
      else
      {
         /* transform edge constraints */
	 if( sourcedata->bigt )
	 {
	    SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->edgecons, sourcedata->nedges) );
            SCIP_CALL( SCIPtransformConss(scip, sourcedata->nedges, sourcedata->edgecons, (*targetdata)->edgecons) );
	 }
	 else
	 {
            SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->edgecons, sourcedata->realnterms * sourcedata->nedges) );
            SCIP_CALL( SCIPtransformConss(scip, sourcedata->realnterms * sourcedata->nedges, sourcedata->edgecons, (*targetdata)->edgecons) );
	 }

	 /* transform constraints */
         if( sourcedata->mode == MODE_PRICE )
         {
	    SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->pathcons, sourcedata->realnterms) );
            SCIP_CALL( SCIPtransformConss(scip, sourcedata->realnterms, sourcedata->pathcons, (*targetdata)->pathcons) );
         }
         /* transform constraints and variables*/
         else if( sourcedata->mode == MODE_FLOW )
         {
	    if( sourcedata->bigt )
	    {
	       SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->flowvars, sourcedata->nedges) );
	       SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->pathcons, (sourcedata->nnodes - 1)) );
               SCIP_CALL( SCIPtransformConss(scip, (sourcedata->nnodes - 1), sourcedata->pathcons, (*targetdata)->pathcons) );
	       SCIP_CALL( SCIPtransformVars(scip, sourcedata->nedges, sourcedata->flowvars, (*targetdata)->flowvars) );
	    }
	    else
	    {
	       SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->flowvars, sourcedata->realnterms * sourcedata->nedges) );
               SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->pathcons,  sourcedata->realnterms * (sourcedata->nnodes - 1)) );
               SCIP_CALL( SCIPtransformConss(scip, sourcedata->realnterms * (sourcedata->nnodes - 1), sourcedata->pathcons, (*targetdata)->pathcons) );
	       SCIP_CALL( SCIPtransformVars(scip, sourcedata->nedges * sourcedata->realnterms, sourcedata->flowvars, (*targetdata)->flowvars) );
	    }
         }
      }

      /* transform array of (real) terminals */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->realterms, sourcedata->realnterms) );
      for( i = 0; i < sourcedata->realnterms; ++i )
      {
         (*targetdata)->realterms[i] = sourcedata->realterms[i];
      }
   }
   else
   {
      (*targetdata)->edgevars = NULL;
      (*targetdata)->xval = NULL;
      (*targetdata)->realterms = NULL;
      (*targetdata)->edgecons = NULL;
      (*targetdata)->pathcons = NULL;
      (*targetdata)->flowvars = NULL;
   }
   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransStp)
{
   SCIPdebugPrintf("probdeltransStp \n");

   SCIPdebugMessage("free transformed problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );
   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS* cons;
   SCIP_Real offset;
   PRESOL presolinfo;
   GRAPH* graph;
   SCIP_Bool print;
   int t;
   int k;
   int nedges;
   int nnodes;
   int realnterms;
   int compcentral;
   int reduction;
   char mode;

   assert(scip != NULL);

   presolinfo.fixed = 0;

   /* create graph */
   graph = graph_load(filename, &presolinfo);
   if( graph == NULL )
      return SCIP_READERROR;

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, graph) );

   /* get parameters */
   SCIP_CALL( SCIPgetCharParam(scip, "stp/mode", &mode) );
   SCIP_CALL( SCIPgetIntParam(scip, "stp/compcentral", &compcentral) );
   SCIP_CALL( SCIPgetIntParam(scip, "stp/reduction", &reduction) );
   SCIP_CALL( SCIPgetBoolParam(scip, "stp/emitgraph", &(probdata->emitgraph)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "stp/bigt", &(probdata->bigt)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "stp/printGraph", &print) );

   /* set solving mode */
   if( mode == 'f' )
      probdata->mode = MODE_FLOW;
   else if( mode == 'p' )
      probdata->mode = MODE_PRICE;
   else
      probdata->mode = MODE_CUT;


   /* create a problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, filename) );
   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigStp) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransStp) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransStp) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   if( probdata->mode == MODE_CUT )
      graph_mincut_init(graph);

   /* init shortest path algorithm (needed for reduction) */
   graph_path_init(graph);

   /* select a root node */
   graph->source[0] = central_terminal(graph, compcentral);

   /* print the graph */
   if( print )
   {
      SCIP_CALL( probdataPrintGraph(graph, "OriginalGraph.gml", NULL) );
   }

   /* presolving */
   offset = reduce(graph, reduction);
   probdata->graph = graph_pack(graph);
   graph = probdata->graph;

   /* if graph reduction solved the whole problem, NULL is returned */
   if( graph != NULL )
   {
      /* init shortest path algorithm (needed for creating path variables) */
      graph_path_exit();
      graph_path_init(graph);

      if( print )
      {
         SCIP_CALL( probdataPrintGraph(graph, "ReducedGraph.gml", NULL) );
      }
      nedges = graph->edges;
      nnodes = graph->knots;
      probdata->nnodes = nnodes;
      probdata->nedges = nedges;
      probdata->nterms = graph->terms;
      probdata->nlayers = graph->layers;
      probdata->nvars = probdata->nlayers * nedges;

      /* compute the real number of terminals (nterm-1 iff root is a terminal) */
      realnterms = graph->terms - 1 - graph->term[graph->source[0]];
      probdata->realnterms = realnterms;

      /* set up array of terminals (except for the root) */
      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->realterms, realnterms) );
      t = 0;
      for( k = 0; k < nnodes; ++k )
      {
	 if( graph->term[k] == 0 && k != graph->source[0] )
	 {
	    probdata->realterms[t] = k;
	    printf("realterms %d \n ", probdata->realterms[t]);
	    t += 1;
	 }
      }

      if( probdata->mode == MODE_CUT )
      {
	 /* create and add constraint for Branch and Cut */
         SCIP_CALL( SCIPcreateConsStp(scip, &cons, "stpcons", probdata->graph) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
      else
      {
         /* create and add constraints for Flow or Branch and Price */
         SCIP_CALL( createConstraints(scip, probdata) );
      }
   }
   /* graph reduction solved the whole problem, set vars to zero or NULL */
   else
   {
      probdata->pathcons = NULL;
      probdata->edgecons = NULL;
      probdata->nlayers = 0;
      probdata->nnodes = 0;
      probdata->nedges = 0;
      probdata->nterms = 0;
      probdata->realnterms = 0;
      probdata->nedges = 0;
      probdata->nvars = 0;
      probdata->realterms = NULL;
      probdata->xval = NULL;
   }

   /* create and add initial variables */
   SCIP_CALL( createVariables(scip, probdata, presolinfo.fixed + offset) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   return SCIP_OKAY;
}

/** returns the graph */
GRAPH* SCIPprobdataGetGraph(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->graph;
}


/** returns the number of variables */
int SCIPprobdataGetNVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->nvars;
}

/** returns the array with all variables */
SCIP_VAR** SCIPprobdataGetVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->edgevars;
}

/** returns the number of layers */
int SCIPprobdataGetNLayers(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->nlayers;
}

/** returns the number of edges */
int SCIPprobdataGetNEdges(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->nedges;
}

/** returns the number of terminals  */
int SCIPprobdataGetNTerms(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->nterms;
}

/** returns the number of terminals without the root node  */
int SCIPprobdataGetRNTerms(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->realnterms;
}

/** returns root  */
int SCIPprobdataGetRoot(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   GRAPH* graph;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = probdata->graph;
   assert(graph != NULL);

   return graph->source[0];
}


/** returns the variable for a given index */
SCIP_VAR* SCIPprobdataGetedgeVarByIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   idx
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->edgevars[idx];
}

/** returns the LP solution values */
double* SCIPprobdataGetXval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /*if( probdata->lastlpiters < SCIPgetNLPIterations(scip) )*/
   {
      SCIP_CALL_ABORT( SCIPgetSolVals(scip, sol, probdata->nvars, probdata->edgevars, probdata->xval) );

      /*probdata->lastlpiters = SCIPgetNLPIterations(scip);*/
   }

   return probdata->xval;
}


/** returns all edge constraints */
SCIP_CONS** SCIPprobdataGetEdgeConstraints(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->edgecons;
}

/** returns all path constraints */
SCIP_CONS** SCIPprobdataGetPathConstraints(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->pathcons;
}


/** returns the array with all variables */
int* SCIPprobdataGetRTerms(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->realterms;
}



/** returns the array with all edge variables */
SCIP_VAR** SCIPprobdataGetEdgeVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->edgevars;
}

/* returns if 'T' model is being used */
SCIP_Bool SCIPprobdataIsBigt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->bigt;
}

/** print (undirected) graph in GML format */
SCIP_RETCODE SCIPprobdataPrintGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of the output file */
   SCIP_SOL*             sol,                /**< solution to be printed; or NULL for LP solution */
   SCIP_Bool             printsol            /**< should solution be printed? */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_VAR** edgevars;
   SCIP_Bool* edgemark;
   int e;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   if( !printsol )
   {
      /* print the graph without highlighting a solution */
      SCIP_CALL( probdataPrintGraph( probdata->graph, filename, NULL) );
   }
   else
   {
      edgevars = probdata->edgevars;
      SCIP_CALL( SCIPallocBufferArray(scip, &edgemark, probdata->nedges / 2) );

      /* mark the edges used in the current solution */
      for( e = 0; e < probdata->graph->edges; e += 2 )
         if( !SCIPisZero(scip, SCIPgetSolVal(scip, sol, edgevars[e])) || !SCIPisZero(scip, SCIPgetSolVal(scip, sol, edgevars[e + 1])) )
            edgemark[e / 2] = TRUE;
         else
	    edgemark[e / 2] = FALSE;

      SCIP_CALL( probdataPrintGraph( probdata->graph, filename, edgemark) );
      SCIPfreeBufferArray(scip, &edgemark);
   }

   return SCIP_OKAY;
}

/** add new solution */
SCIP_RETCODE SCIPprobdataAddNewSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            nval,               /**< array [0..nvars], nval[v] = 1 if node v is in the solution, nval[v] = 0 if not */
   SCIP_SOL*             sol,                /**< the new solution */
   SCIP_HEUR*            heur,               /**< heuristic data */
   SCIP_Bool*            success             /**< denotes whether the new solution has been successfully added */
   )
{
   SCIP_PROBDATA* probdata;
   GRAPH* graph;
   SCIP_VAR** edgevars;
   SCIP_Real* edgecost;
   PATH* path;
   SCIP_VAR* var;
   SCIP_VAR** pathvars;
   char varname[SCIP_MAXSTRLEN];
   int e;
   int t;
   int nedges;
   int realnterms;
   int tail;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   edgevars = probdata->edgevars;
   graph = probdata->graph;

   assert(edgevars != NULL);
   assert(graph != NULL);

   /* create a new primal solution (initialized to zero) */
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   /* create path variables (Price mode) or set the flow vars (Flow mode) corresponding to the new solution */
   if( probdata->mode != MODE_CUT )
   {
      nedges = probdata->nedges;

      assert(nedges > 0);

      realnterms = probdata->realnterms;
      SCIP_CALL( SCIPallocMemoryArray(scip, &edgecost, nedges) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &path, graph->knots) );

      /* Flow mode */
      if ( probdata->mode == MODE_FLOW )
      {
         pathvars = NULL;
      }
      /* Price mode */
      else
      {
	 SCIP_CALL( SCIPallocMemoryArray(scip, &pathvars, realnterms) );
      }

      /* mark the tree generated by nvals */
      for( e = 0; e < nedges; e++ )
      {
         if( SCIPisEQ(scip, nval[e], 1.0) )
	    edgecost[e] = graph->cost[e]/nedges;
	 else
	    edgecost[e] = SCIPinfinity(scip);
      }

      for( e = 0; e < graph->knots; e++ )
         graph->mark[e] = 1; /* @todo mark only terminals? */
      graph_path_exec(graph, FSP_MODE, graph->source[0], edgecost, path);

      /* create and add path variables (Price mode) or set the flow variables (Flow mode) */
      for( t = 0; t < realnterms; ++t )
      {
         if( probdata->mode == MODE_PRICE )
         {
	    /* create a new path variable */
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "PathVar%d_X", t);
            var = NULL;
            SCIP_CALL( SCIPcreateVarBasic(scip, &var, varname, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, var) );
            SCIP_CALL( SCIPaddCoefLinear(scip, probdata->pathcons[t], var, 1.0) );
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, 1.0) );
	    pathvars[t] = var;
         }
         tail = probdata->realterms[t];
	 /* walk from terminal t to the root */
         while( tail != graph->source[0] )
         {
            if( !probdata->bigt )
	    {
               if( probdata->mode == MODE_PRICE )
               {
		  /* add the new path variable to the constraints corresponding to the current edge */
                  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[t * nedges + path[tail].edge], var, 1.0) );
               }
               else
               {
		  /* set the flow variable corresponding to the current edge */
		  SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->flowvars[t * nedges + path[tail].edge], 1.0) );
               }
	    }
            else
	    {
               if( probdata->mode == MODE_PRICE )
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[path[tail].edge], var, 1.0) );
               }
               else
               {
                  SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->flowvars[path[tail].edge], 1.0) );
               }
	    }
            tail = graph->tail[path[tail].edge];
         }
      }

      /* store new solution value */
      SCIP_CALL( SCIPsetSolVals(scip, sol, probdata->nvars, edgevars, nval) );

      /* try to add new solution to scip and free it immediately */
      SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, TRUE, TRUE, success) );

      /* free local arrays */
      SCIPfreeMemoryArrayNull(scip, &edgecost);
      SCIPfreeMemoryArrayNull(scip, &path);
      SCIPfreeMemoryArrayNull(scip, &pathvars);
   }
   /* Cut mode */
   else
   {
      /* store the new solution value */
      SCIP_CALL( SCIPsetSolVals(scip, sol, probdata->nvars, edgevars, nval) );

      /* try to add new solution to scip and free it immediately */
      SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, TRUE, TRUE, success) );
   }

   return SCIP_OKAY;
}
