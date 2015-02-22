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
   SCIP_CONS**           degcons;            /**< array of (node) degree constraints */
   SCIP_CONS**           edgecons;           /**< array of constraints */
   SCIP_CONS**           pathcons;           /**< array of constraints */
   SCIP_CONS*            hopcons;            /**< hop constraint */
   SCIP_CONS*            prizecons;          /**< prize constraint */
   SCIP_VAR** 		 edgevars;	     /**< array of edge variables */
   SCIP_VAR**            flowvars;           /**< array of edge variables (needed only in the Flow mode) */
   SCIP_VAR*             offsetvar;          /**< variable to model the objective offset */
   SCIP_Real             offset;             /**< offset of the problem, computed during the presolving */
   SCIP_Real*            xval;               /**< values of the edge variables */
   int* 	         realterms;          /**< array of all terminals except the root */
   SCIP_Bool             emitgraph;          /**< emitgraph */
   int                   nedges;             /**< number of edges */
   int                   nterms;             /**< number of terminals */
   int                   realnterms;         /**< number of terminals except the root */
   int                   nlayers;            /**< number of layers */
   int                   nnodes;             /**< number of nodes */
   int                   nvars;              /**< number of variables */
   int                   stp_type;           /**< STP type */
   SCIP_Longint          lastlpiters;        /**< Branch and Cut */
   SCIP_Bool             copy;               /**< is this the problem data of a copy/sub-MIP? */
   FILE*                 logfile;            /**< logfile for DIMACS challenge */
   FILE**                origlogfile;        /**< pointer to original problem data logfile pointer */

   /** for FiberSCIP **/
   SCIP_Bool             ug;                 /**< inidicat if this ug dual bound is set or not */
   int                   nSolvers;           /**< the number of solvers */
   SCIP_Real             ugDual;             /**< dual bound set by ug */
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
   if( graph != NULL )
      (*probdata)->stp_type = graph->stp_type;
   else
      (*probdata)->stp_type = STP_UNDIRECTED;
   (*probdata)->copy = FALSE;
   (*probdata)->logfile = NULL;
   (*probdata)->origlogfile = NULL;

   (*probdata)->ug = FALSE;
   (*probdata)->nSolvers =0;
   (*probdata)->ugDual = 0.0;

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
   SCIPdebugPrintf ("probdataFree \n");
   assert(scip != NULL);
   assert(probdata != NULL);

   /* release variables */
   for( e = 0; e < (*probdata)->nvars; ++e )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->edgevars[e]) );
   }
   SCIPfreeMemoryArrayNull(scip, &(*probdata)->edgevars);

   if( (*probdata)->offsetvar != NULL )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->offsetvar) );
   }

   /* Degree-Constrained STP? */
   if( (*probdata)->stp_type == STP_DEG_CONS )
   {
      assert((*probdata)->mode == MODE_CUT);

      /* release degree constraints */
      for( t = 0; t < (*probdata)->nnodes; ++t)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &((*probdata)->degcons[t])) );
      }
      SCIPfreeMemoryArrayNull(scip, &((*probdata)->degcons));
   }

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
         (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Root", n);
	 SCIPgmlWriteNode(file, (unsigned int)n, label, "rectangle", "#666666", NULL);
	 m = 1;
      }
      else if( graph->term[n] == 0 )
      {
	 (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Terminal %d", n, e + 1);
	 SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#ff0000", NULL);
	 e += 1;
      }
      else
      {
         (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Node %d", n, n + 1 - e - m);
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

/** create (edge-) HOP constraint (Cut Mode only) */
static
SCIP_RETCODE createHopConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )

{
   GRAPH* graph;
   SCIP_Real rhs;
   assert(scip != NULL);
   assert(probdata != NULL);

   SCIPdebugPrintf("createHopeConstraint \n");
   graph = probdata->graph;
   assert(graph != NULL);
   rhs = graph->hoplimit;
   printf("Hop limit: %f \n ", rhs);
   /* TODO: when presolving is enabled: set rhs = rhs - (number of fixed edges) */

   SCIP_CALL( SCIPcreateConsLinear ( scip, &(probdata->hopcons), "HopConstraint", 0, NULL, NULL,
         -SCIPinfinity(scip), rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPaddCons(scip, probdata->hopcons) );

   return SCIP_OKAY;
}

/** create (node-) degree constraints (Cut Mode only) */
static
SCIP_RETCODE createDegreeConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )

{
   GRAPH* graph;
   char consname[SCIP_MAXSTRLEN];
   int k;
   int nnodes;

   assert(scip != NULL);
   assert(probdata != NULL);

   SCIPdebugPrintf("createDegreeConstraints \n");
   graph = probdata->graph;
   nnodes = probdata->nnodes;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->degcons), nnodes ) );

   for( k = 0; k < nnodes; ++k )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "DegreeConstraint%d", k);
      SCIP_CALL( SCIPcreateConsLinear ( scip, &(probdata->degcons[k]), consname, 0, NULL, NULL,
            -SCIPinfinity(scip), graph->maxdeg[k], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(scip, probdata->degcons[k]) );
   }

   return SCIP_OKAY;
}




/** create Prize constraints (Cut Mode only) */
static
SCIP_RETCODE createPrizeConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )

{
   char consname[SCIP_MAXSTRLEN];


   assert(scip != NULL);
   assert(probdata != NULL);

   SCIPdebugPrintf("createPrizeConstraints \n");

#if 0
   SCIP_CALL( SCIPallocMemory(scip, &(probdata->prizecons)) );
#endif

   (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PrizeConstraint");
   SCIP_CALL( SCIPcreateConsLinear ( scip, &(probdata->prizecons), consname, 0, NULL, NULL,
         -SCIPinfinity(scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   // TODO SET TO 1 DISABLE WARNING OUTPUT FOR
   SCIP_CALL( SCIPaddCons(scip, probdata->prizecons) );


   return SCIP_OKAY;
}



/** create constraints (in Flow or Price Mode) */
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
      /* create |E| edge constraints (aggregated) */
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
      if( !probdata->bigt )
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
                  (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PathConstraintT%d", k2);
                  SCIP_CALL( SCIPcreateConsLinear(scip, &( probdata->pathcons[k2] ), consname,
                        0, NULL, NULL, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
                  SCIP_CALL( SCIPaddCons(scip, probdata->pathcons[k2]) );
	       }
	       /* if node k is a terminal, set RHS = 1 */
	       else
	       {
                  (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PathConstraintT%d", k2);
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
   char varname[SCIP_MAXSTRLEN];
   SCIP_VAR* var;
   SCIP_Real* edgecost;
   int tail;
   int k2;
   int e;
   int t;
   int k;

   assert(scip != NULL);
   assert(probdata != NULL);

   graph = probdata->graph;
   SCIPdebugPrintf("createVariables \n");

   /* if the graph reduction solved the whole problem, NULL is returned */
   if( graph != NULL )
   {
      int nedges = probdata->nedges;
      int nvars = probdata->nvars;
      int realnterms = probdata->realnterms;
      int root = graph->source[0];
      int nnodes = graph->knots;
      SCIP_Bool objint = SCIPisIntegral(scip, offset);

      assert(nedges = graph->edges);

      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->xval, nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->edgevars, nvars) );

      /* Cut mode */
      if( probdata->mode == MODE_CUT )
      {
	 assert(probdata->nlayers == 1);

         for( e = 0; e < nedges; ++e )
         {
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d_%d", graph->tail[e], graph->head[e]);
            SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->edgevars[e], varname, 0.0, 1.0, graph->cost[e], SCIP_VARTYPE_BINARY) );
            SCIP_CALL( SCIPaddVar(scip, probdata->edgevars[e]) );
            objint = objint && SCIPisIntegral(scip, graph->cost[e]);
         }

#if 0
         for( e = 0; e < nedges; e=e+2 )
         {
            SCIP_CONS* cons;
            cons = NULL;
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d_%d", graph->tail[e], graph->head[e]);
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, varname,
                  0, NULL, NULL, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->edgevars[e], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->edgevars[e+1], 1.0) );
         }
#endif
         /* Hop-Constrained STP */
         if( graph->stp_type == STP_HOP_CONS )
	 {
	    int hopfactor;
	    for( e = 0; e < nedges; ++e )
            {
               /* TODO: When presolving is used: MODIFY */
               hopfactor = 1;
               SCIP_CALL( SCIPaddCoefLinear(scip, probdata->hopcons, probdata->edgevars[e], hopfactor) );
	    }
	 }

         /* Degree-Constrained STP */
         if( graph->stp_type == STP_DEG_CONS )
	 {
	    for( k = 0; k < nnodes; ++k )
            {
	       for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
               {
		  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->degcons[k], probdata->edgevars[e], 1.0) );
		  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->degcons[k], probdata->edgevars[flipedge(e)], 1.0) );
	       }
	    }
	 }
	 /* PRIZECOLLECTING STP */
         if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT )
	 {
	    for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
            {
               if( !Is_term(graph->term[graph->head[e]]) )
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->prizecons, probdata->edgevars[e], 1.0) );
                  /*printf("add to cons %d %d \n", graph->tail[e], graph->head[e] );*/
		  /* variables are prefered to be branched on */
		  SCIP_CALL( SCIPchgVarBranchPriority( scip, probdata->edgevars[e], 1) );
               }
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
            SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) ); /* TODO c09 chg*/

            if( !probdata->bigt )
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
	    probdata->edgevars[e] = var;
         }
      }

      /* Price mode */
      if( probdata->mode == MODE_PRICE )
      {
	 PATH* path;

         /* the flow variables are not used in the Price mode */
	 probdata->flowvars = NULL;

         /* compute shortest paths to the root */
         SCIP_CALL( SCIPallocMemoryArray(scip, &path, graph->knots) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &edgecost, nedges) );
         for( e = 0; e < nedges; ++e )
	 {
            edgecost[e] = graph->cost[e];
	 }

	 for( e = 0; e < graph->knots; e++ )
	 {
            graph->mark[e] = 1;
	 }

         graph_path_exec(graph, FSP_MODE, root, edgecost, path);

         /* create and add initial path variables (one for each real terminal) */
         for( t = 0; t < realnterms; ++t )
         {
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "PathVar%d_0", t);
            var = NULL;
            SCIP_CALL( SCIPcreateVarBasic(scip, &var, varname, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, var) );
	    SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, probdata->pathcons[t], var, 1.0) );
            tail = probdata->realterms[t];
            while( tail != root )
            {
	       if( !probdata->bigt )
	       {
                  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[t * nedges + path[tail].edge], var, 1.0) );
	       }
	       else
	       {
                  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[path[tail].edge], var, 1.0) );
	       }

               tail = graph->tail[path[tail].edge];
            }
         }

         /* free local arrays */
         SCIPfreeMemoryArray(scip, &edgecost);
         SCIPfreeMemoryArray(scip, &path);
      }
      /* Flow mode */
      else if( probdata->mode == MODE_FLOW )
      {
	 /* store the number of disparate flows (commodities) in nflows */
	 int nflows;
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
                     SCIP_CALL( SCIPaddCoefLinear(scip, probdata->pathcons[t * (nnodes - 1)  + k2], probdata->flowvars[t * nedges + e], -1.0) );
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
   SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->offsetvar, varname, 1.0, 1.0, offset, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, probdata->offsetvar) );

   return SCIP_OKAY;
}


/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** copies user data of source SCIP for the target SCIP */
static
SCIP_DECL_PROBCOPY(probcopyStp)
{
   GRAPH* graphcopy;

   SCIPdebugPrintf("########################## probcopy ###########################\n");

   graphcopy = graph_copy(sourcedata->graph);

   graph_path_init(graphcopy);

   if( sourcedata->mode == MODE_CUT )
      graph_mincut_init(graphcopy);

   SCIP_CALL( probdataCreate(scip, targetdata, graphcopy) );

   (*targetdata)->mode = sourcedata->mode;
   (*targetdata)->bigt = sourcedata->bigt;
   (*targetdata)->nlayers = sourcedata->nlayers;
   (*targetdata)->nedges = sourcedata->nedges;
   (*targetdata)->nnodes = sourcedata->nnodes;
   (*targetdata)->nterms = sourcedata->nterms;
   (*targetdata)->realnterms = sourcedata->realnterms;
   (*targetdata)->emitgraph = sourcedata->emitgraph;
   (*targetdata)->nvars = sourcedata->nvars;
   (*targetdata)->copy = TRUE;
   (*targetdata)->offsetvar = NULL;

   if( sourcedata->offsetvar != NULL && SCIPvarIsActive(sourcedata->offsetvar) )
   {
      SCIP_Bool success;

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->offsetvar, &((*targetdata)->offsetvar), varmap, consmap, global, &success) );
      assert(success);

      SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->offsetvar) );
   }

   if( sourcedata->nedges > 0 )
   {
      SCIP_Bool success;
      int v;
      int c;
      int i;

      /* transform variables */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->xval, sourcedata->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->edgevars, sourcedata->nvars) );

      for( v = sourcedata->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->edgevars[v], &((*targetdata)->edgevars[v]), varmap, consmap, global, &success) );
         assert(success);

         SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->edgevars[v]) );
      }

      /* Cut mode */
      if( sourcedata->mode == MODE_CUT )
      {
         (*targetdata)->edgecons = NULL;
         (*targetdata)->pathcons = NULL;
	 if( sourcedata->stp_type == STP_DEG_CONS )
	 {
            SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->degcons, sourcedata->nnodes) );

            for( c = sourcedata->nnodes - 1; c >= 0; --c )
            {
               SCIP_CALL( SCIPgetConsCopy(sourcescip, scip, sourcedata->degcons[c], &((*targetdata)->degcons[c]),
                     SCIPconsGetHdlr(sourcedata->degcons[c]), varmap, consmap,
                     SCIPconsGetName(sourcedata->degcons[c]),
                     SCIPconsIsInitial(sourcedata->degcons[c]),
                     SCIPconsIsSeparated(sourcedata->degcons[c]),
                     SCIPconsIsEnforced(sourcedata->degcons[c]),
                     SCIPconsIsChecked(sourcedata->degcons[c]),
                     SCIPconsIsPropagated(sourcedata->degcons[c]),
                     SCIPconsIsLocal(sourcedata->degcons[c]),
                     SCIPconsIsModifiable(sourcedata->degcons[c]),
                     SCIPconsIsDynamic(sourcedata->degcons[c]),
                     SCIPconsIsRemovable(sourcedata->degcons[c]),
                     SCIPconsIsStickingAtNode(sourcedata->degcons[c]),
                     global, &success) );
               assert(success);

               SCIP_CALL( SCIPcaptureCons(scip, (*targetdata)->degcons[c]) );
            }
	 }
      }
      /* Price or Flow mode */
      else
      {
         /* transform edge constraints */
	 if( sourcedata->bigt )
	 {
	    SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->edgecons, sourcedata->nedges) );

            for( c = sourcedata->nedges - 1; c >= 0; --c )
            {
               SCIP_CALL( SCIPgetConsCopy(sourcescip, scip, sourcedata->edgecons[c], &((*targetdata)->edgecons[c]),
                     SCIPconsGetHdlr(sourcedata->edgecons[c]), varmap, consmap,
                     SCIPconsGetName(sourcedata->edgecons[c]),
                     SCIPconsIsInitial(sourcedata->edgecons[c]),
                     SCIPconsIsSeparated(sourcedata->edgecons[c]),
                     SCIPconsIsEnforced(sourcedata->edgecons[c]),
                     SCIPconsIsChecked(sourcedata->edgecons[c]),
                     SCIPconsIsPropagated(sourcedata->edgecons[c]),
                     SCIPconsIsLocal(sourcedata->edgecons[c]),
                     SCIPconsIsModifiable(sourcedata->edgecons[c]),
                     SCIPconsIsDynamic(sourcedata->edgecons[c]),
                     SCIPconsIsRemovable(sourcedata->edgecons[c]),
                     SCIPconsIsStickingAtNode(sourcedata->edgecons[c]),
                     global, &success) );
               assert(success);

               SCIP_CALL( SCIPcaptureCons(scip, (*targetdata)->edgecons[c]) );
            }
	 }
	 else
	 {
            SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->edgecons, sourcedata->realnterms * sourcedata->nedges) );
            for( c = sourcedata->realnterms * sourcedata->nedges - 1; c >= 0; --c )
            {
               SCIP_CALL( SCIPgetConsCopy(sourcescip, scip, sourcedata->edgecons[c], &((*targetdata)->edgecons[c]),
                     SCIPconsGetHdlr(sourcedata->edgecons[c]), varmap, consmap,
                     SCIPconsGetName(sourcedata->edgecons[c]),
                     SCIPconsIsInitial(sourcedata->edgecons[c]),
                     SCIPconsIsSeparated(sourcedata->edgecons[c]),
                     SCIPconsIsEnforced(sourcedata->edgecons[c]),
                     SCIPconsIsChecked(sourcedata->edgecons[c]),
                     SCIPconsIsPropagated(sourcedata->edgecons[c]),
                     SCIPconsIsLocal(sourcedata->edgecons[c]),
                     SCIPconsIsModifiable(sourcedata->edgecons[c]),
                     SCIPconsIsDynamic(sourcedata->edgecons[c]),
                     SCIPconsIsRemovable(sourcedata->edgecons[c]),
                     SCIPconsIsStickingAtNode(sourcedata->edgecons[c]),
                     global, &success) );
               assert(success);

               SCIP_CALL( SCIPcaptureCons(scip, (*targetdata)->edgecons[c]) );
            }
	 }

	 /* transform constraints */
         if( sourcedata->mode == MODE_PRICE )
         {
	    SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->pathcons, sourcedata->realnterms) );
            for( c = sourcedata->realnterms - 1; c >= 0; --c )
            {
               SCIP_CALL( SCIPgetConsCopy(sourcescip, scip, sourcedata->pathcons[c], &((*targetdata)->pathcons[c]),
                     SCIPconsGetHdlr(sourcedata->pathcons[c]), varmap, consmap,
                     SCIPconsGetName(sourcedata->pathcons[c]),
                     SCIPconsIsInitial(sourcedata->pathcons[c]),
                     SCIPconsIsSeparated(sourcedata->pathcons[c]),
                     SCIPconsIsEnforced(sourcedata->pathcons[c]),
                     SCIPconsIsChecked(sourcedata->pathcons[c]),
                     SCIPconsIsPropagated(sourcedata->pathcons[c]),
                     SCIPconsIsLocal(sourcedata->pathcons[c]),
                     SCIPconsIsModifiable(sourcedata->pathcons[c]),
                     SCIPconsIsDynamic(sourcedata->pathcons[c]),
                     SCIPconsIsRemovable(sourcedata->pathcons[c]),
                     SCIPconsIsStickingAtNode(sourcedata->pathcons[c]),
                     global, &success) );
               assert(success);

               SCIP_CALL( SCIPcaptureCons(scip, (*targetdata)->pathcons[c]) );
            }
         }
         /* transform constraints and variables */
         else if( sourcedata->mode == MODE_FLOW )
         {
	    if( sourcedata->bigt )
	    {
	       SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->flowvars, sourcedata->nedges) );
	       SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->pathcons, (sourcedata->nnodes - 1)) );
               for( c = sourcedata->nnodes - 2; c >= 0; --c )
               {
                  SCIP_CALL( SCIPgetConsCopy(sourcescip, scip, sourcedata->pathcons[c], &((*targetdata)->pathcons[c]),
                        SCIPconsGetHdlr(sourcedata->pathcons[c]), varmap, consmap,
                        SCIPconsGetName(sourcedata->pathcons[c]),
                        SCIPconsIsInitial(sourcedata->pathcons[c]),
                        SCIPconsIsSeparated(sourcedata->pathcons[c]),
                        SCIPconsIsEnforced(sourcedata->pathcons[c]),
                        SCIPconsIsChecked(sourcedata->pathcons[c]),
                        SCIPconsIsPropagated(sourcedata->pathcons[c]),
                        SCIPconsIsLocal(sourcedata->pathcons[c]),
                        SCIPconsIsModifiable(sourcedata->pathcons[c]),
                        SCIPconsIsDynamic(sourcedata->pathcons[c]),
                        SCIPconsIsRemovable(sourcedata->pathcons[c]),
                        SCIPconsIsStickingAtNode(sourcedata->pathcons[c]),
                        global, &success) );
                  assert(success);

                  SCIP_CALL( SCIPcaptureCons(scip, (*targetdata)->pathcons[c]) );
               }

               for( v = (*targetdata)->nedges - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->flowvars[v], &((*targetdata)->flowvars[v]), varmap, consmap, global, &success) );
                  assert(success);

                  SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->flowvars[v]) );
               }
	    }
	    else
	    {
	       SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->flowvars, sourcedata->realnterms * sourcedata->nedges) );
               SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->pathcons,  sourcedata->realnterms * (sourcedata->nnodes - 1)) );
               for( c = sourcedata->realnterms * (sourcedata->nnodes - 1) - 1; c >= 0; --c )
               {
                  SCIP_CALL( SCIPgetConsCopy(sourcescip, scip, sourcedata->pathcons[c], &((*targetdata)->pathcons[c]),
                        SCIPconsGetHdlr(sourcedata->pathcons[c]), varmap, consmap,
                        SCIPconsGetName(sourcedata->pathcons[c]),
                        SCIPconsIsInitial(sourcedata->pathcons[c]),
                        SCIPconsIsSeparated(sourcedata->pathcons[c]),
                        SCIPconsIsEnforced(sourcedata->pathcons[c]),
                        SCIPconsIsChecked(sourcedata->pathcons[c]),
                        SCIPconsIsPropagated(sourcedata->pathcons[c]),
                        SCIPconsIsLocal(sourcedata->pathcons[c]),
                        SCIPconsIsModifiable(sourcedata->pathcons[c]),
                        SCIPconsIsDynamic(sourcedata->pathcons[c]),
                        SCIPconsIsRemovable(sourcedata->pathcons[c]),
                        SCIPconsIsStickingAtNode(sourcedata->pathcons[c]),
                        global, &success) );
                  assert(success);

                  SCIP_CALL( SCIPcaptureCons(scip, (*targetdata)->pathcons[c]) );
               }

               for( v = sourcedata->nedges * sourcedata->realnterms - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->flowvars[v], &((*targetdata)->flowvars[v]), varmap, consmap, global, &success) );
                  assert(success);

                  SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->flowvars[v]) );
               }
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

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigStp)
{
   SCIPdebugPrintf("probdelorigStp \n");

   SCIPdebugMessage("free original problem data\n");

   if( (*probdata)->graph != NULL )
   {
      if ( (*probdata)->mode == MODE_CUT )
         graph_mincut_exit((*probdata)->graph);

      graph_path_exit((*probdata)->graph);

      graph_free((*probdata)->graph, TRUE);
   }

   /* free the (original) probdata */
   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransStp)
{
   SCIP_Real timelimit;

   SCIPdebugPrintf("probtransStp \n");

   /* adjust time limit to take into account reading time */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   timelimit -= SCIPgetReadingTime(scip);
   timelimit = MAX(0.0,timelimit);
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timelimit) );

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
   (*targetdata)->logfile = sourcedata->logfile;
   (*targetdata)->origlogfile = &(sourcedata->logfile);

   if( sourcedata->offsetvar != NULL )
   {
      SCIP_CALL( SCIPtransformVar(scip, sourcedata->offsetvar, &(*targetdata)->offsetvar) );
   }
   else
      (*targetdata)->offsetvar = NULL;

   if( sourcedata->nedges > 0 )
   {
      int i;

      /* transform variables */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->xval, sourcedata->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->edgevars, sourcedata->nvars) );
      SCIP_CALL( SCIPtransformVars(scip, sourcedata->nvars, sourcedata->edgevars, (*targetdata)->edgevars) );

      /* Cut mode */
      if( sourcedata->mode == MODE_CUT )
      {
         (*targetdata)->edgecons = NULL;
         (*targetdata)->pathcons = NULL;
	 if( sourcedata->stp_type == STP_DEG_CONS )
	 {
	    SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->degcons, sourcedata->nnodes) );
            SCIP_CALL( SCIPtransformConss(scip, sourcedata->nnodes, sourcedata->degcons, (*targetdata)->degcons) );
	 }
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



static
SCIP_DECL_PROBEXITSOL(probexitsolStp)
{

   SCIP_PROBDATA* probd;


   assert(scip != NULL);
   probd = SCIPgetProbData(scip);
#if 0
   GRAPH* graph;
   graph = probd->graph;
   assert(graph != NULL);


   SCIP_SOL* sol;
   SCIP_VAR** edgevars;
   GRAPH* graph;
   graph = probd->graph;
   if( graph != NULL && graph->stp_type == STP_GRID )
   {
      sol = SCIPgetBestSol(scip);

      /* print the coordinates of the best solution */
      if( sol != NULL &&  !(SCIPgetSubscipDepth(scip) > 0) )
      {
	 SCIP_QUEUE* queue;
         int**  coords;
         int*  ncoords;
         int*  nodecoords;
         int e;
         int i;
	 int* pnode;
	 int root;
	 int size = 0;
         int grid_dim;

	 coords = graph->grid_coordinates;
	 assert(coords != NULL);
	 ncoords = graph->grid_ncoords;
	 nodecoords = NULL;
	 root = graph->source[0];
	 assert(root >= 0);
	 grid_dim = graph->grid_dim;
	 assert(grid_dim > 1);
         edgevars = probd->edgevars;
	 assert(ncoords != NULL);

	 /* BFS until all terminals are reached */
	 SCIP_CALL( SCIPqueueCreate(&queue, size, 2) );
         for( e = 0; e < graph->edges; e++ )
            if( !SCIPisZero(scip, SCIPgetSolVal(scip, sol, edgevars[e])) )
               size++;
	 assert(size > 0);

	 SCIP_CALL( SCIPqueueInsert(queue, &root) );

	 printf("Coordinates of the best found solution: \n");
	 while( !SCIPqueueIsEmpty(queue) )
	 {
            pnode = (SCIPqueueRemove(queue));
            for( e = graph->outbeg[*pnode]; e != EAT_LAST; e = graph->oeat[e] )
            {
	       if( !SCIPisZero(scip, SCIPgetSolVal(scip, sol, edgevars[e])) )
	       {
	          graph_grid_coordinates(coords, &nodecoords, ncoords, graph->tail[e], grid_dim);
                  printf("(%d", nodecoords[0]);
                  for( i = 1; i < grid_dim; i++ )
                     printf(", %d", nodecoords[i]);
                  printf(") --> ");
                  graph_grid_coordinates(coords, &nodecoords, ncoords, graph->head[e], grid_dim);
                  printf("(%d", nodecoords[0]);
                  for( i = 1; i < grid_dim; i++ )
                     printf(", %d", nodecoords[i]);
                  printf(") \n");

		  SCIP_CALL( SCIPqueueInsert(queue, &(graph->head[e])) );
	       }
	    }
	 }

	 SCIPqueueFree(&queue);

         /*
           for( e = 0; e < graph->edges; e++ )
           if( !SCIPisZero(scip, SCIPgetSolVal(scip, sol, edgevars[e])) )
           {
           graph_grid_coordinates(coords, &nodecoords, ncoords, graph->tail[e], grid_dim);
           printf("(%d", nodecoords[0]);
           for( i = 1; i < grid_dim; i++ )
           printf(", %d", nodecoords[i]);
           printf(") --> ");
           graph_grid_coordinates(coords, &nodecoords, ncoords, graph->head[e], grid_dim);
           printf("(%d", nodecoords[0]);
           for( i = 1; i < grid_dim; i++ )
           printf(", %d", nodecoords[i]);
           printf(") \n");
           }
         */
         free(nodecoords);
      }
   }
#endif

   if( probd->logfile != NULL )
   {
      int success;
      SCIP_Real factor = 1.0;

      if( probd->stp_type ==  STP_MAX_NODE_WEIGHT )
         factor = -1.0;

      SCIPprobdataWriteLogLine(scip, "End\n");
      SCIPprobdataWriteLogLine(scip, "\n");
      SCIPprobdataWriteLogLine(scip, "SECTION Run\n");
      if( probd->ug )
      {
         SCIPprobdataWriteLogLine(scip, "Threads %d\n", probd->nSolvers);
         SCIPprobdataWriteLogLine(scip, "Time %.1f\n", SCIPgetTotalTime(scip));
         SCIPprobdataWriteLogLine(scip, "Dual %16.9f\n", factor * probd->ugDual);
      }
      else
      {
         SCIPprobdataWriteLogLine(scip, "Threads 1\n");
         SCIPprobdataWriteLogLine(scip, "Time %.1f\n", SCIPgetTotalTime(scip));
         SCIPprobdataWriteLogLine(scip, "Dual %16.9f\n", factor * SCIPgetDualbound(scip));
      }
      SCIPprobdataWriteLogLine(scip, "Primal %16.9f\n", factor * SCIPgetPrimalbound(scip));
      SCIPprobdataWriteLogLine(scip, "End\n");

      if( SCIPgetNSols(scip) > 0 )
      {
         SCIPprobdataWriteLogLine(scip, "\n");
         SCIPprobdataWriteLogLine(scip, "SECTION Finalsolution\n");

         SCIP_CALL( SCIPprobdataWriteSolution(scip, probd->logfile) );
         SCIPprobdataWriteLogLine(scip, "End\n");
      }

      success = fclose(probd->logfile);
      if( success != 0 )
      {
         SCIPerrorMessage("An error occurred while closing file <%s>\n", probd->logfile);
         return SCIP_FILECREATEERROR;
      }

      probd->logfile = NULL;
      *(probd->origlogfile) = NULL;
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
   int nedges;
   int nnodes;
   int realnterms;
   int compcentral;
   int reduction;
   char mode;
   char probtype[16];
   char printfs = FALSE;
   char* logfilename;
   char* tmpfilename;
   char* probname;
   assert(scip != NULL);

   presolinfo.fixed = 0;

   /* create graph */
   graph = graph_load(filename, &presolinfo);
   if( printfs )
      printf("load type :: %d \n\n", graph->stp_type);
   if( graph == NULL )
      return SCIP_READERROR;
   if( printfs )
      printf("fixed: %f \n\n", presolinfo.fixed );

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, graph) );

   /* get parameters */
   SCIP_CALL( SCIPgetCharParam(scip, "stp/mode", &mode) );
   SCIP_CALL( SCIPgetIntParam(scip, "stp/compcentral", &compcentral) );
   SCIP_CALL( SCIPgetIntParam(scip, "stp/reduction", &reduction) );
   SCIP_CALL( SCIPgetBoolParam(scip, "stp/emitgraph", &(probdata->emitgraph)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "stp/bigt", &(probdata->bigt)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "stp/printGraph", &print) );
   SCIP_CALL( SCIPgetStringParam(scip, "stp/logfile", &logfilename) );

   if( logfilename != NULL && logfilename[0] != '\0' )
   {
      probdata->logfile = fopen(logfilename, "w");

      if( probdata->logfile == NULL )
      {
         SCIPerrorMessage("cannot create file <%s> for writing\n", logfilename);
         SCIPprintSysError(logfilename);
         return SCIP_FILECREATEERROR;
      }
   }

   /* copy filename */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpfilename, filename, (int)strlen(filename)+1) );

   SCIPsplitFilename(tmpfilename, NULL, &probname, NULL, NULL);

   SCIPfreeBufferArray(scip, &tmpfilename);


   /* create a problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );
   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigStp) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransStp) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransStp) );
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolStp) );
   SCIP_CALL( SCIPsetProbCopy(scip, probcopyStp) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   SCIPprobdataWriteLogLine(scip, "SECTION Comment\n");
   SCIPprobdataWriteLogLine(scip, "Name %s\n", filename);

   switch( graph->stp_type )
   {
   case STP_UNDIRECTED:
      strcpy(probtype, "SPG");
      break;

   case STP_PRIZE_COLLECTING:
      strcpy(probtype, "PCSPG");
      break;

   case STP_ROOTED_PRIZE_COLLECTING:
      strcpy(probtype, "RPCST");
      break;

   case STP_NODE_WEIGHTS:
      strcpy(probtype, "NWSPG");
      break;

   case STP_DEG_CONS:
      strcpy(probtype, "DCST");
      break;

   case STP_GRID:
      strcpy(probtype, "RSMT");
      break;

   case STP_OBSTACLES_GRID:
      strcpy(probtype, "OARSMT");
      break;

   case STP_MAX_NODE_WEIGHT:
      strcpy(probtype, "MWCS");
      break;

   case STP_HOP_CONS:
      strcpy(probtype, "HCDST");
      break;

   default:
      strcpy(probtype, "UNKNOWN");
   }
   SCIPprobdataWriteLogLine(scip, "Problem %s\n", probtype);
   SCIPprobdataWriteLogLine(scip, "Program SCIP-Jack\n");
   SCIPprobdataWriteLogLine(scip, "Version 0.1\n");
   SCIPprobdataWriteLogLine(scip, "End\n");
   SCIPprobdataWriteLogLine(scip, "\n");
   SCIPprobdataWriteLogLine(scip, "SECTION Solutions\n");

   /* set solving mode */
   if( mode == 'f' )
   {
      assert(graph->stp_type == STP_UNDIRECTED);
      probdata->mode = MODE_FLOW;
   }
   else if( mode == 'p' )
   {
      assert(graph->stp_type == STP_UNDIRECTED);
      probdata->mode = MODE_PRICE;
   }
   else
      probdata->mode = MODE_CUT;

   assert(graph != NULL );
   /* init shortest path algorithm (needed for reduction) */

   graph_path_init(graph);

   /* select a root node */
   if( !(graph->stp_type == STP_DIRECTED) && compcentral != CENTER_DEG && graph->stp_type != STP_PRIZE_COLLECTING
      && graph->stp_type != STP_ROOTED_PRIZE_COLLECTING && graph->stp_type != STP_HOP_CONS )
      graph->source[0] = central_terminal(graph, compcentral);

   /* print the graph */
   if( print )
      SCIP_CALL( probdataPrintGraph(graph, "OriginalGraph.gml", NULL) );

   /* presolving */
   offset = reduce(graph, reduction, scip);

   graph_path_exit(graph);

   probdata->graph = graph_pack(graph);

   graph = probdata->graph;

   /* */
   if( graph != NULL )
   {
      int t;
      int k;

      /* init shortest path algorithm (needed for creating path variables) */
      graph_path_init(graph);

      if( probdata->mode == MODE_CUT )
         graph_mincut_init(graph);

      if( print )
         SCIP_CALL( probdataPrintGraph(graph, "ReducedGraph.gml", NULL) );

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
	    SCIPdebugMessage("realterms[%d] = %d \n ", t, probdata->realterms[t]);
	    t += 1;
	 }
      }

      if( probdata->mode == MODE_CUT )
      {
	 /* create and add constraint for Branch and Cut */
         SCIP_CALL( SCIPcreateConsStp(scip, &cons, "stpcons", probdata->graph) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

         /* if the problem is a HOP-Constrained-STP, an additional constraint is required */
	 if( graph->stp_type == STP_HOP_CONS )
	    SCIP_CALL( createHopConstraint(scip, probdata) );

	 /* if the problem is a Degree-Constrained-STP, additional constraints are required */
	 if( graph->stp_type == STP_DEG_CONS )
	    SCIP_CALL( createDegreeConstraints(scip, probdata) );

	 if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT)
            SCIP_CALL( createPrizeConstraints(scip, probdata) );
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

   /* add the new offset (offset) and the offset that had already been part of the problem before it was read in (presolinfo.fixed) */
   probdata->offset = presolinfo.fixed + offset;

   /* create and add initial variables */
   SCIP_CALL( createVariables(scip, probdata, probdata->offset ) );

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


/** returns offset of the problem */
SCIP_Real SCIPprobdataGetOffset(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->offset;
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
SCIP_Real* SCIPprobdataGetXval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_Real* vals;
   int e;
   int nedges;
   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   vals = probdata->xval;
   //assert(vals != NULL);
   nedges = probdata->nedges;
   assert(nedges >= 0);
   /*if( probdata->lastlpiters < SCIPgetNLPIterations(scip) )*/
   {
      SCIP_CALL_ABORT( SCIPgetSolVals(scip, sol, nedges, probdata->edgevars, vals) );
      /*probdata->lastlpiters = SCIPgetNLPIterations(scip);*/
   }

   for( e = 0; e < nedges; e++ )
      vals[e] = fmax(0.0, fmin(vals[e], 1.0));

   return vals;
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

      /* print the graph highlighting a solution */
      SCIP_CALL( probdataPrintGraph( probdata->graph, filename, edgemark) );

      SCIPfreeBufferArray(scip, &edgemark);
   }

   return SCIP_OKAY;
}

/** writes the best solution to a file */
SCIP_RETCODE SCIPprobdataWriteSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< file to write best solution to; or NULL, to write to stdout */
   )
{

   /*
     SCIP_PROBDATA* probd;
     probd = SCIPgetProbData(scip);
     graph = probd->graph;
   */
   SCIP_SOL* sol;
   SCIP_VAR** edgevars;
   GRAPH* graph;
   IDX** ancestors;
   IDX* curr;
   SCIP_PROBDATA* probdata;
   int  e;
   int  norgedges;
   int  norgnodes;
   int  nsolnodes;
   int  nsoledges;
   char* orgedges;
   char* orgnodes;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   graph = probdata->graph;//SCIPprobdataGetGraph(probdata);

   edgevars = probdata->edgevars;

   assert(graph != NULL);
   sol = SCIPgetBestSol(scip);
   nsolnodes = 0;
   nsoledges = 0;
   norgedges = graph->orgedges;
   norgnodes = graph->orgknots;
   assert(norgedges >= 0);
   assert(norgnodes >= 1);

   SCIP_CALL( SCIPallocBufferArray(scip, &orgedges, norgedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orgnodes, norgnodes) );

   for( e = 0; e < norgedges; e++ )
      orgedges[e] = FALSE;
   for( e = 0; e < norgnodes; e++ )
      orgnodes[e] = FALSE;
   ancestors = graph->ancestors;
   if( graph->stp_type == STP_UNDIRECTED || graph->stp_type == STP_DEG_CONS
      || graph->stp_type == STP_NODE_WEIGHTS || graph->stp_type == STP_HOP_CONS )
   {
      //printf("in: %d \n", norgnodes);
      curr = graph->fixedges;
      while( curr != NULL )
      {
	 //printf("index: %d max: %d \n", curr->index, norgedges);
	 if( orgedges[curr->index] == FALSE )
	 {
	    orgedges[curr->index] = TRUE;
	    nsoledges++;
	 }

	 //printf("indexorgtail: %d max: %d \n", graph->orgtail[curr->index], norgnodes);
	 if( orgnodes[graph->orgtail[curr->index]] == FALSE )
	 {
            orgnodes[graph->orgtail[curr->index]] = TRUE;
	    nsolnodes++;
	 }
         // printf("indexorghead: %d max: %d \n", graph->orghead[curr->index], norgnodes);
	 if( orgnodes[graph->orghead[curr->index]] == FALSE )
	 {
	    orgnodes[graph->orghead[curr->index]] = TRUE;
	    nsolnodes++;
	 }
         curr = curr->parent;
      }
      //printf("in2: %d \n", norgnodes);
      for( e = 0; e < graph->edges; e++ )
      {
	 if( !SCIPisZero(scip, SCIPgetSolVal(scip, sol, edgevars[e])) )
	 {
	    /* iterate through the list of ancestors */
            curr = ancestors[e];
            while( curr != NULL )
            {
               if( orgedges[curr->index] == FALSE )
               {
                  orgedges[curr->index] = TRUE;
                  nsoledges++;
               }
               if( orgnodes[graph->orgtail[curr->index]] == FALSE )
               {
                  orgnodes[graph->orgtail[curr->index]] = TRUE;
                  nsolnodes++;
               }
               if( orgnodes[graph->orghead[curr->index]] == FALSE )
               {
                  orgnodes[graph->orghead[curr->index]] = TRUE;
                  nsolnodes++;
               }
               curr = curr->parent;
            }
	 }
      }

      // printf("norgnodes: %d \n", norgnodes);
      //printf("norgedges: %d \n", norgedges);

      SCIPprobdataWriteLogLine(scip, "Vertices %d\n", nsolnodes);

      for( e = 0; e < norgnodes; e++ )
         if( orgnodes[e] == TRUE )
            SCIPinfoMessage(scip, file, "V %d\n", e + 1);

      SCIPprobdataWriteLogLine(scip, "Edges %d\n", nsoledges);
      if( graph->stp_type == STP_HOP_CONS )
      {
	 for( e = 0; e < norgedges; e++ )
	    if( orgedges[e] == TRUE )
               SCIPinfoMessage(scip, file, "E %d %d\n", graph->orgtail[e] + 1, graph->orghead[e] + 1);
      }
      else
      {
         for( e = 0; e < norgedges; e += 2 )
         {
            if( orgedges[e] == TRUE || orgedges[e + 1] == TRUE )
               SCIPinfoMessage(scip, file, "E %d %d\n", graph->orgtail[e] + 1, graph->orghead[e] + 1);
         }
      }

   }
   else if( graph->stp_type == STP_GRID )
   {
      int**  coords;
      int*  ncoords;
      int*  nodecoords;
      int*  nodenumber;
      int i;
      int nodecount;
      int grid_dim;
      char strdim[256];
      coords = graph->grid_coordinates;
      assert(coords != NULL);
      ncoords = graph->grid_ncoords;
      nodecoords = NULL;
      grid_dim = graph->grid_dim;
      assert(grid_dim > 1);
      assert(grid_dim < 256);
      assert(norgedges = graph->edges);
      assert(norgnodes = graph->knots);
      SCIP_CALL( SCIPallocBufferArray(scip, &nodenumber, norgnodes) );
      strcpy(strdim, "P");
      for( i = 1; i < grid_dim; i++ )
         strcat(strdim, "P");

      assert(ncoords != NULL);

      /* mark solution nodes and edges */
      for( e = 0; e < norgedges; e++ )
      {
         if( !SCIPisZero(scip, SCIPgetSolVal(scip, sol, edgevars[e])) )
         {
            nsoledges++;
            assert(orgedges[e] == FALSE);
            orgedges[e] = TRUE;
            if( orgnodes[graph->tail[e]] == FALSE )
            {
               orgnodes[graph->tail[e]] = TRUE;
               nsolnodes++;
            }
            if( orgnodes[graph->head[e]] == FALSE )
            {
               orgnodes[graph->head[e]] = TRUE;
               nsolnodes++;
            }
         }
      }

      SCIPprobdataWriteLogLine(scip, "Edges %d\n", nsoledges);
      SCIPprobdataWriteLogLine(scip, "Points %d\n", nsolnodes);
      nodecount = 0;
      for( e = 0; e < norgnodes; e++ )
      {
         if( orgnodes[e] == TRUE )
         {
	    nodenumber[e] = nodecount++;
            SCIPprobdataWriteLogLine(scip, "%s ", strdim);
            graph_grid_coordinates(coords, &nodecoords, ncoords, e, grid_dim);
            for( i = 0; i < grid_dim; i++ )
            {
               SCIPprobdataWriteLogLine(scip, "%d ", nodecoords[i]);
            }
            SCIPprobdataWriteLogLine(scip, "\n");
         }
      }
      assert(nodecount == nsolnodes);
      for( e = 0; e < norgedges; e += 2 )
      {
	 if( orgedges[e] == TRUE || orgedges[e + 1] == TRUE )
	    SCIPinfoMessage(scip, file, "E %d %d\n", nodenumber[graph->orgtail[e]] + 1, nodenumber[graph->orghead[e]] + 1);
      }
      SCIPfreeBufferArray(scip, &nodenumber);
   }
   else if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT
      || graph->stp_type == STP_ROOTED_PRIZE_COLLECTING )
   {
      int root;
      root = graph->source[0];
      assert(root >= 0);

      /* switch the terminal property (back to its original state), and mark the old terminals */
      /*
        for( k = 0; k < graph->knots; k++ )
        {
        if( Is_term(graph->term[k]) && k != root )
        {
        for( e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
        if( graph->tail[e] != root )
        break;
        assert(e != EAT_LAST);
        graph->term[graph->tail[e]] = 0;
        graph->term[k] = -2;
        }
        }

        printf("norgmodeledges: %d \n", graph->norgmodeledges);
        printf("norgmodelknots: %d \n", graph->norgmodelknots);
      */
      for( e = 0; e <= graph->edges; e++ )
      {
	 if( e == graph->edges || !SCIPisZero(scip, SCIPgetSolVal(scip, sol, edgevars[e])) )
	 {
	    /* iterate through the list of ancestors/fixed edged*/
	    if( e < graph->edges )
               curr = ancestors[e];
	    else
	       curr = graph->fixedges;
            while( curr != NULL )
            {
               if( e < graph->edges && graph->stp_type == STP_MAX_NODE_WEIGHT )
               {
                  if( !SCIPisZero(scip, SCIPgetSolVal(scip, sol, edgevars[flipedge(e)])) )
                  {
                     curr = curr->parent;
                     continue;
                  }
               }

               //	       assert(graph->head[curr->index] != root);
               //     if( curr->index] == root )
               //	  printf("rootedge: %d->%d \n",

	       if( curr->index < graph->norgmodeledges ) //if( graph->term[graph->head[curr->index]] != -2 )
	       {
                  // if( (graph->stp_type == STP_ROOTED_PRIZE_COLLECTING)? 1 : graph->tail[curr->index] != root )
                  //{
                  if( orgedges[curr->index] == FALSE )
                  {
                     orgedges[curr->index] = TRUE;
                     nsoledges++;
                  }
                  if( orgnodes[graph->orgtail[curr->index]] == FALSE )
                  {
                     orgnodes[graph->orgtail[curr->index]] = TRUE;
                     nsolnodes++;
                  }
                  if( orgnodes[graph->orghead[curr->index]] == FALSE )
                  {
                     orgnodes[graph->orghead[curr->index]] = TRUE;
                     nsolnodes++;
                  }
                  /*   }
                       else
                       {
                       assert(graph->tail[curr->index] == root);
                       if( orgnodes[graph->orghead[curr->index]] == FALSE )
                       {
                       orgnodes[graph->orghead[curr->index]] = TRUE;
                       nsolnodes++;
                       }
                       } */
	       }
	       else if( graph->orghead[curr->index] < graph->norgmodelknots )
	       {
                  if( orgnodes[graph->orghead[curr->index]] == FALSE )
                  {
		     if( graph->orghead[curr->index] == 2901 )
		     {
                        // printf("head: %d %d d\n", graph->orgtail[curr->index], graph->orghead[curr->index]);
		     }
                     orgnodes[graph->orghead[curr->index]] = TRUE;
                     nsolnodes++;
                  }
	       }
               curr = curr->parent;
            }
	 }
      }

      if( graph->stp_type == STP_ROOTED_PRIZE_COLLECTING && orgnodes[root] == FALSE )
      {
	 orgnodes[root] = TRUE;
	 assert(nsolnodes == 0);
	 nsolnodes = 1;
      }
      SCIPprobdataWriteLogLine(scip, "Vertices %d\n", nsolnodes);

      for( e = 0; e < norgnodes; e++ )
         if( orgnodes[e] == TRUE )
	    SCIPinfoMessage(scip, file, "V %d\n", e + 1);

      SCIPprobdataWriteLogLine(scip, "Edges %d\n", nsoledges);
      for( e = 0; e < norgedges; e += 2 )
	 if( orgedges[e] == TRUE || orgedges[e + 1] == TRUE )
	    SCIPinfoMessage(scip, file, "E %d %d\n", graph->orgtail[e] + 1, graph->orghead[e] + 1);
   }

   SCIPfreeBufferArray(scip, &orgnodes);
   SCIPfreeBufferArray(scip, &orgedges);

   return SCIP_OKAY;
}


/** writes a line to the log file */
void SCIPprobdataWriteLogLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   SCIP_PROBDATA* probdata;
   va_list ap;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   if( probdata->logfile == NULL )
      return;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintInfo(SCIPgetMessagehdlr(scip), probdata->logfile, formatstr, ap);
   va_end(ap);
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
   SCIP_VAR** edgevars;
   GRAPH* graph;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   edgevars = probdata->edgevars;

   graph = probdata->graph;
   assert(graph != NULL);
   assert(edgevars != NULL);

   /* create a new primal solution (initialized to zero) */
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   /* create path variables (Price mode) or set the flow vars (Flow mode) corresponding to the new solution */
   if( probdata->mode != MODE_CUT )
   {
      SCIP_Real* edgecost;
      SCIP_Real* flowvals;
      PATH* path;
      SCIP_VAR** pathvars;
      SCIP_VAR* var;
      char varname[SCIP_MAXSTRLEN];
      int realnterms = probdata->realnterms;
      int tail;
      int nedges = probdata->nedges;
      int e;
      int t;

      assert(nedges > 0);

      /* allocate memory for the values of the flow variables */
      SCIP_CALL( SCIPallocMemoryArray(scip, &flowvals, nedges * (probdata->bigt ? 1 : realnterms)) );
      BMSclearMemoryArray(flowvals, nedges * (probdata->bigt ? 1 : realnterms));

      /* allocate memory for the edgecost and the path array (both used for computing shortest paths) */
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
	    edgecost[e] = graph->cost[e] / nedges;
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
		  flowvals[t * nedges + path[tail].edge] = 1.0;
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
                  /* increment the flow variable corresponding to the current edge */
		  flowvals[path[tail].edge] += 1.0;
               }
	    }
            tail = graph->tail[path[tail].edge];
         }
      }

      /* store the new solution value */
      SCIP_CALL( SCIPsetSolVals(scip, sol, probdata->nvars, edgevars, nval) );
      if( probdata->mode == MODE_FLOW )
      {
	 SCIP_CALL( SCIPsetSolVals(scip, sol, nedges * (probdata->bigt ? 1 : realnterms) , probdata->flowvars, flowvals) );
      }

      /* try to add new solution to scip and free it immediately */
      SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, TRUE, TRUE, success) );

      /* free local arrays */
      SCIPfreeMemoryArrayNull(scip, &flowvals);
      SCIPfreeMemoryArrayNull(scip, &edgecost);
      SCIPfreeMemoryArrayNull(scip, &path);
      SCIPfreeMemoryArrayNull(scip, &pathvars);
   }
   /* Cut mode */
   else
   {
      SCIP_Bool feasible;
      int e;
      int nvars = probdata->nvars;
      int fails = 0;
      /* check whether the new solution is valid with respect to the original bounds */
#if 1
      // if( SCIPgetDepth(scip) != -1 )

      for( e = 0; e < nvars; e++ )
      {
         if( SCIPisGT(scip, nval[e], SCIPvarGetUbGlobal(edgevars[e])) ||  SCIPisGT(scip, SCIPvarGetLbGlobal(edgevars[e]), nval[e]) )
         {
            //  printf("XXXXXXXXXXXXXXXXx  \n solution violates orginal bounds (%d %f bounds: %f %f) \n  XXXXXXXXXXXX \n", e ,nval[e], SCIPvarGetLbGlobal(edgevars[e]), SCIPvarGetUbGlobal(edgevars[e]) );
            *success = FALSE;
            fails++;
            SCIP_CALL( SCIPfreeSol(scip, &sol) );
            return SCIP_OKAY;
         }
      }

#endif

      /* post-processing of solution for MWCS and PCSPG */
      if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT )
      {
         int k;

         for( k = 0; k < graph->knots; ++k )
         {
            /* is the kth node a terminal other than the root? */
            if( Is_term(graph->term[k]) && k != graph->source[0] )
            {
               int origterm;
               int edge1 = graph->inpbeg[k];
               int edge2 = graph->ieat[edge1];
               assert(graph->ieat[edge2] == EAT_LAST);

               if( !SCIPisZero(scip, graph->cost[edge2]) )
               {
                  int tmp = edge1;
                  edge1 = edge2;
                  edge2 = tmp;
               }
               assert(SCIPisZero(scip, graph->cost[edge2]));

               if( nval[edge2] > 0.5 )
               {
                  assert(nval[edge1] < 0.5);
                  continue;
               }
               assert(nval[edge1] > 0.5);
               assert(nval[edge2] < 0.5);

               origterm = graph->tail[edge2];

               for( e = graph->inpbeg[origterm]; e != EAT_LAST; e = graph->ieat[e] )
               {
                  if( nval[e] > 0.5 )
                  {
                     //printf("two edges to one terminal\n");
                     nval[edge1] = 0;
                     nval[edge2] = 1;
                     break;
                  }
               }
            }
         }
      }

#if 0
      if( fails > 0 )
      {
         int v;
         SCIP_Bool* edgemark;
         printf("solution violates orginal bounds ");


         SCIP_CALL( SCIPfreeSol(scip, &sol) );

         printf("nfails: %d       \n", fails);


         SCIP_CALL( SCIPallocBufferArray(scip, &edgemark, probdata->nedges) );
         for( v = 0; v < nvars; v++ )
            if(  SCIPisEQ(scip, SCIPvarGetUbGlobal(edgevars[v]), 0.0) )
            {
	       edgemark[v/2] = FALSE;

            }
            else
            {

	       edgemark[v/2] = TRUE;
            }
         SCIP_CALL( SCIPprobdataPrintGraph2(probdata->graph, "GraphSub.gml", edgemark) );
	 SCIPfreeBufferArray(scip, &edgemark);

         return SCIP_OKAY;
      }

#endif
      /* store the new solution value */
      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, edgevars, nval) );

      if( probdata->offsetvar != NULL && SCIPvarIsActive(probdata->offsetvar) )
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->offsetvar, 1.0) );
      }

      SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, TRUE, TRUE, &feasible) );

      /* printf("checked sol: feasible=%d\n", feasible); */

      /* check solution for feasibility in original problem space */
      if( !feasible )
      {
         SCIP_CALL( SCIPcheckSolOrig(scip, sol, &feasible, TRUE, TRUE) );

         /* printf("checked sol org: feasible=%d\n", feasible); */

         if( feasible )
         {
            SCIP_SOL* newsol;
            SCIP_VAR** origvars;
            SCIP_VAR* var;
            int norigvars;
            int v;

            SCIP_CALL( SCIPcreateOrigSol(scip, &newsol, SCIPsolGetHeur(sol)) );
            origvars = SCIPgetOrigVars(scip);
            norigvars = SCIPgetNOrigVars(scip);

            for( v = 0; v < norigvars; ++v )
            {
               var = origvars[v];
               SCIP_CALL( SCIPsetSolVal(scip, newsol, var, SCIPgetSolVal(scip, sol, var)) );
            }

            SCIP_CALL( SCIPfreeSol(scip, &sol) );
            sol = newsol;
         }
      }

      /* try to add new solution to scip and free it immediately */
      if( feasible )
      {
#ifndef NDEBUG
         SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, TRUE, TRUE, &feasible) );
         assert(feasible);
#endif

         SCIP_CALL( SCIPaddSolFree(scip, &sol, success) );
      }
      else
      {
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
	 *success = FALSE;
      }

   }

   return SCIP_OKAY;
}


/** print graph (in undirected form) in GML format */
SCIP_RETCODE SCIPprobdataPrintGraph2(
   const GRAPH*          graph,              /**< Graph to be printed */
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
         (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Root", n);
	 SCIPgmlWriteNode(file, (unsigned int)n, label, "rectangle", "#666666", NULL);
	 m = 1;
      }
      else if( graph->term[n] == 0 )
      {
	 (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Terminal %d", n, e + 1);
	 SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#ff0000", NULL);
	 e += 1;
      }
      else
      {
         (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Node %d", n, n + 1 - e - m);
         SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#336699", NULL);
      }
   }

   /* write all edges (undirected) */
   for( e = 0; e < graph->edges; e += 2 )
   {
      (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%8.2f", graph->cost[e]);
      if( edgemark != NULL && edgemark[e / 2] == TRUE )
	 SCIPgmlWriteEdge(file, (unsigned int)graph->tail[e], (unsigned int)graph->head[e], label, "#ff0000");
      //else
      // SCIPgmlWriteEdge(file, (unsigned int)graph->tail[e], (unsigned int)graph->head[e], label, NULL);
   }

   /* write GML format closing */
   SCIPgmlWriteClosing(file);

   return SCIP_OKAY;
}

/** returns problem type */
int SCIPprobdataGetType(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->stp_type;
}

/** writes end of log file */
SCIP_RETCODE SCIPprobdataWriteLogfileEnd(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   if( probdata->logfile != NULL )
   {
      int success;
      SCIP_Real factor = 1.0;

      if( probdata->stp_type ==  STP_MAX_NODE_WEIGHT )
         factor = -1.0;

      SCIPprobdataWriteLogLine(scip, "End\n");
      SCIPprobdataWriteLogLine(scip, "\n");
      SCIPprobdataWriteLogLine(scip, "SECTION Run\n");
      if( probdata->ug )
      {
         SCIPprobdataWriteLogLine(scip, "Threads %d\n", probdata->nSolvers);
         SCIPprobdataWriteLogLine(scip, "Time %.1f\n", SCIPgetTotalTime(scip));
         SCIPprobdataWriteLogLine(scip, "Dual %16.9f\n", factor * probdata->ugDual);
      }
      else
      {
         SCIPprobdataWriteLogLine(scip, "Threads 1\n");
         SCIPprobdataWriteLogLine(scip, "Time %.1f\n", SCIPgetTotalTime(scip));
         SCIPprobdataWriteLogLine(scip, "Dual %16.9f\n", factor * SCIPgetDualbound(scip));
      }
      SCIPprobdataWriteLogLine(scip, "Primal %16.9f\n", factor * SCIPgetPrimalbound(scip));
      SCIPprobdataWriteLogLine(scip, "End\n");

      if( SCIPgetNSols(scip) > 0 )
      {
         SCIPprobdataWriteLogLine(scip, "\n");
         SCIPprobdataWriteLogLine(scip, "SECTION Finalsolution\n");

         SCIP_CALL( SCIPprobdataWriteSolution(scip, probdata->logfile) );
         SCIPprobdataWriteLogLine(scip, "End\n");
      }

      success = fclose(probdata->logfile);
      if( success != 0 )
      {
         SCIPerrorMessage("An error occurred while closing file <%s>\n", probdata->logfile);
         return SCIP_FILECREATEERROR;
      }

      probdata->logfile = NULL;
   }


   return SCIP_OKAY;
}

/** writes end of log file */
void SCIPprobdataSetDualBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             dual                /**< dual bound */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   probdata->ug = TRUE;
   probdata->ugDual = dual;
}

/** writes end of log file */
void SCIPprobdataSetNSolvers(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nSolvers            /**< the number of solvers */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   probdata->nSolvers = nSolvers;
}
