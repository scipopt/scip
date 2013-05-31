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
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "probdata_stp.h"

#include "scip/scip.h"
#include "cons_stp.h"
#include "grph.h"
#include "portab.h"

#define   DEFAULT_COMPCENTRAL  1
#define   DEFAULT_EMITGRAPH    FALSE
#define   DEFAULT_REDUCTION    4

#define CENTER_OK    0
#define CENTER_DEG   1
#define CENTER_SUM   2
#define CENTER_MIN   3
#define CENTER_ALL   4


/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the graph of the steiner tree problem
 */
struct SCIP_ProbData
{
   GRAPH*                graph;        /**< the graph */
   int compcentral;
   int reduction;
   SCIP_Bool emitgraph;

   double* xval;
   SCIP_Longint lastlpiters;

   SCIP_VAR** vars;
   int nedges;
   int nlayers;
   int nvars;
};

#if 0
/* parameter list */
static PDSC jack_para_desc_tab[] =
{
   /* ACHTUNG: Tabelle muss alphabetisch sortiert sein !
    *
    *  Keyword        Type  Default       Outputformat
    */
   // {  "BACK_CUT",     DEF_BOL(0),            "Try Back-Cuts     : %s" },
   {  "COMP_CENTRAL", DEF_INT(1),            "Comp. Central Term: %d" },
   {  "EMIT_GRAPH",   DEF_BOL(0),            NULL                     },
   {  "REDUCTION",    DEF_INT(4),            "Reduction Level   : %d" },
   // {  "CREEP_FLOW",   DEF_BOL(1),            "Use Creep-Flow    : %s" },
   // {  "DISJUNCT_CUT", DEF_BOL(1),            "Only disjunct Cuts: %s" },
   // {  "FLOW_INIT",    DEF_BOL(0),            "Init all Flows    : %s" },
   // {  "FLOW_SEP",     DEF_BOL(1),            "Flow-Separator    : %s" },
   // {  "HEUR_INTER",   DEF_INT(10),           "Heuristic Interval: %d" },
   // {  "HEUR_RATIO",   DEF_DBL(0.8),          "Heuristic Ratio   : %g" },
   // {  "HOPS_SEP",     DEF_BOL(0),            "Hops Separator    : %s" },
   // {  "INCVIOL",      DEF_DBL(1.0),          "Min Violation Inc : %g" },
   // {  "MAX_HOPS",     DEF_INT(1000000),      "Max Hops          : %d" },
   // {  "MAX_SLACK",    DEF_DBL(0.0001),       "Maximum Slack     : %g" },
   // {  "MINPROG",      DEF_DBL(0.0001),       "Minimum Progress  : %g" },
   // {  "MINVIOL",      DEF_DBL(0.0001),       "Minimum Violation : %g" },
   // {  "MIN_WEIGHT",   DEF_DBL(0.0001),       "Minimum Weight    : %g" },
   // {  "MOST_VIOLATE", DEF_BOL(0),            NULL                     },
   // {  "MSTR_SEP",     DEF_BOL(0),            "MST Separator     : %s" },
   // {  "NESTED_CUT",   DEF_BOL(0),            "Try Nested-Cuts   : %s" },
   // {  "ONLY_BEST",    DEF_INT(0),            "Only best Inequals: %d" },
   // {  "PLAN_SEP",     DEF_BOL(0),            "Planar Separator  : %s" },
   // {  "RSTR_FACTOR",  DEF_DBL(0.6),          "Restructure Factor: %g" },
   // {  "RSTR_MIN",     DEF_INT(100),          "Restructure Min   : %d" },
   // {  "TAK_EVALRUNS", DEF_INT(10),           "Tak. Iter. Runs   : %d" },
   // {  "TAK_INITRUNS", DEF_INT(10),           "Tak. Init. Runs   : %d" },
   // {  "TAK_LEAFRUNS", DEF_INT(15),           "Tak. Leaf Runs    : %d" },
   // {  "TAK_ROOTRUNS", DEF_INT(50),           "Tak. Root Runs    : %d" },
   // {  "USE_RAYHEUR",  DEF_BOL(0),            "Use Ray. Heuristik: %s" },
   // {  "USE_TAKHEUR",  DEF_BOL(1),            "Use Tak. Heuristik: %s" },
   // {  "WAVE_SEP",     DEF_BOL(0),            "Do Wave Separation: %s" },
   // {  "WRITE_FIG",    DEF_INT(0),            "Write fig         : %d" },
};
#endif

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
   int v;

   assert(scip != NULL);
   assert(probdata != NULL);

   for( v = 0; v < (*probdata)->nvars; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->vars[v]) );
   }

   SCIPfreeMemoryArrayNull(scip, &(*probdata)->xval);
   SCIPfreeMemoryArrayNull(scip, &(*probdata)->vars);

   /* free probdata */
   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

/** create initial columns */
static
SCIP_RETCODE createVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Real             offset
   )
{
   GRAPH* graph;
   SCIP_VAR* offsetvar;
   char name[SCIP_MAXSTRLEN];
   int nlayers;
   int nedges;
   int nvars;
   int layer;
   int e;

   graph = probdata->graph;

   /* if graph reduction solved the whole problem, NULL is returned */
   if( graph != NULL )
   {
      nlayers = graph->layers;
      nedges = graph->edges;
      nvars = nlayers * nedges;
      probdata->nlayers = nlayers;
      probdata->nedges = nedges;
      probdata->nvars = nvars;

      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->vars, nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->xval, nvars) );

      for( layer = 0; layer < nlayers; ++layer )
      {
         for( e = 0; e < nedges; ++e )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d_%d_%d", graph->head[e], graph->tail[e], layer);

            SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->vars[layer * nedges + e], name, 0, 1, graph->cost[e], SCIP_VARTYPE_BINARY) );
            SCIP_CALL( SCIPaddVar(scip, probdata->vars[layer * nedges + e]) );
         }
      }
   }
   else
   {
      probdata->nlayers = 0;
      probdata->nedges = 0;
      probdata->nvars = 0;
      probdata->vars = NULL;
      probdata->xval = NULL;
   }

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "OFFSET");

   SCIP_CALL( SCIPcreateVarBasic(scip, &offsetvar, name, 1, 1, offset, SCIP_VARTYPE_CONTINUOUS) );
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
   SCIPdebugMessage("free original problem data\n");

   graph_mincut_exit();
   graph_path_exit();

   if( (*probdata)->graph != NULL )
      graph_free((*probdata)->graph);

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransStp)
{
   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->graph) );

   (*targetdata)->nlayers = sourcedata->nlayers;
   (*targetdata)->nedges = sourcedata->nedges;
   (*targetdata)->nvars = sourcedata->nvars;
   (*targetdata)->compcentral = sourcedata->compcentral;
   (*targetdata)->reduction = sourcedata->reduction;
   (*targetdata)->emitgraph = sourcedata->emitgraph;

   if( sourcedata->nvars > 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->vars, sourcedata->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->xval, sourcedata->nvars) );

      /* transform all variables */
      SCIP_CALL( SCIPtransformVars(scip, sourcedata->nvars, sourcedata->vars, (*targetdata)->vars) );
   }
   else
   {
      (*targetdata)->vars = NULL;
      (*targetdata)->xval = NULL;
   }

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransStp)
{
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

   assert(scip != NULL);
   presolinfo.fixed = 0;

   graph = graph_load(filename, &presolinfo);

   if( graph == NULL )
      return SCIP_READERROR;

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, graph) );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, filename) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigStp) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransStp) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransStp) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   SCIP_CALL( SCIPaddIntParam(scip, "stp/compcentral",
         "Comp. Central Term: 0 disable, 1 max. degree, 2 min. dist. sum to all terminals, 3 min. max. dist., 4 min. dist to all nodes",
         &probdata->compcentral, FALSE, DEFAULT_COMPCENTRAL, 0, 4, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "stp/reduction",
         "Reduction: 0 disable, 5 maximum",
         &probdata->reduction, FALSE, DEFAULT_REDUCTION, 0, 5, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "stp/emitgraph",
         "Emit graph",
         &probdata->emitgraph, FALSE, DEFAULT_EMITGRAPH, NULL, NULL) );

#if 0
   /* tell SCIP that the objective will be always integral */
   SCIP_CALL( SCIPsetObjIntegral(scip) );
#endif

   graph_path_init(graph);
   graph_mincut_init(graph);

   /* define root node */
   graph->source[0] = central_terminal(graph, probdata->compcentral);

   /* presolving */
   offset = reduce(graph, probdata->reduction);

   probdata->graph = graph_pack(graph);

   SCIP_CALL( createVariables(scip, probdata, presolinfo.fixed + offset) );

   /* if graph reduction solved the whole problem, NULL is returned */
   if( probdata->graph != NULL )
   {
      SCIP_CALL( SCIPcreateConsStp(scip, &cons, "stpcons", probdata->graph) );

      SCIP_CALL( SCIPaddCons(scip, cons) );

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

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

/** returns the array with all variables */
SCIP_VAR** SCIPprobdataGetVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->vars;
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

/** returns the variable for a given layer and edge */
SCIP_VAR* SCIPprobdataGetVarByLayerAndEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   layer,
   int                   edge
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->vars[layer * probdata->nedges + edge];
}

/** returns the variable for a given index */
SCIP_VAR* SCIPprobdataGetVarByIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   idx
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->vars[idx];
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

   if( sol == NULL && SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return NULL;

   /*if( probdata->lastlpiters < SCIPgetNLPIterations(scip) )*/
   {
      SCIP_CALL_ABORT( SCIPgetSolVals(scip, sol, probdata->nvars, probdata->vars, probdata->xval) );

      /*probdata->lastlpiters = SCIPgetNLPIterations(scip);*/
   }

   return probdata->xval;
}


/**@} */
