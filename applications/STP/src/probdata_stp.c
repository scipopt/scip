/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probdata_stp.c
 * @brief  Problem data for Steiner problems
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Michael Winkler
 * @author Daniel Rehfeldt
 *
 * This file implements the problem data for Steiner problems. For more details see \ref STP_PROBLEMDATA page.
 *
 * @page STP_PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins.
 *
 * The problem data contains the (preprocessed) graph, several constraints and further information.
 *
 * The function SCIPprobdataCreate(), which is called in the \ref reader_stp.c "reader plugin", parses the input file
 * reduces the graph, initializes the problem data structure and creates the problem in the SCIP environment.
 * See the body of the function SCIPprobdataCreate() for more details.
 *
 * A list of all interface methods can be found in probdata_stp.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "probdata_stp.h"
#include <stdio.h>
#include "cons_stp.h"
#include "dualascent.h"
#include "graph.h"
#include "reduce.h"
#include "stptest.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/misc.h"
#include "scip/struct_misc.h"
#ifdef WITH_UG
#include "branch_stp.h"
#endif

#define STP_SYM_PRIZE
#define STP_AGG_SYM

#define FLOWB       FALSE
#define USEOFFSETVAR TRUE

#define SYM_CONS_LIMIT 20000         /**< maximum number of symmetry inequalities for MWCSP and PCSPG */
#define CYC_CONS_LIMIT 10000         /**< maximum number of symmetry inequalities for PCSPG */

#define CUT_MAXNTERMINALS 500
#define CUT_MAXNEDGES     10000
#define CUT_MAXTOTNEDGES  50000

#ifdef WITH_UG
const char*
getBranchLinearConsName(const char* names, int i);

const char*
getBranchSetppcConsName(const char* names, int i);

extern
int getUgRank(void);
#endif


/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the graph of the Steiner tree problem
 */
struct SCIP_ProbData
{
   int                   mode;               /**< solving mode selected by the user (Cut, Price, Flow) */
   SCIP_Bool             bigt;               /**< stores whether the 'T' model is being used (not relevant in Cut mode) */
   GRAPH*                graph;              /**< the graph */
#ifdef WITH_UG
   GRAPH*                orggraph;          /**< copy of presolved graph needed for UG */
#endif
   SCIP_CONS**           degcons;            /**< array of (node) degree constraints */
   SCIP_CONS**           edgecons;           /**< array of constraints */
   SCIP_CONS**           pathcons;           /**< array of path constraints (needed only in the Price mode) */
#if FLOWB
   SCIP_CONS**           flowbcons;          /**< flow balance constraints */
#endif
   SCIP_CONS**           prizesymcons;       /**< prize-collecting symmetry constraints (to improve LP) */
   SCIP_CONS**           prizecyclecons;     /**< prize-collecting cycle constraints (to improve LP) */
   SCIP_CONS*            hopcons;            /**< hop constraint */
   SCIP_CONS*            budgetcons;         /**< budget constraint */
   SCIP_CONS*            prizecons;          /**< prize constraint */
   SCIP_VAR**            edgevars;           /**< array of edge variables */
   SCIP_VAR**            flowvars;           /**< array of edge variables (needed only in the Flow mode) */
#if USEOFFSETVAR
   SCIP_VAR*             offsetvar;          /**< variable to model the objective offset */
#endif
   SCIP_Real             offset;             /**< offset of the problem, computed during the presolving */
   SCIP_Real*            xval;               /**< values of the edge variables */
   SCIP_Longint          lastlpiters;        /**< Branch and Cut */
   int*                  pctermsorder;       /**< terminal order for PCSPG/MWCSP */
   int*                  realterms;          /**< array of all terminals except the root */
   int                   nedges;             /**< number of edges */
   int                   norgedges;          /**< number of original edges of the model */
   int                   nterms;             /**< number of terminals */
   int                   realnterms;         /**< number of terminals except the root */
   int                   nlayers;            /**< number of layers */
   int                   nnodes;             /**< number of nodes */
   int                   nnonterms;          /**< number of nodes */
   int                   nvars;              /**< number of variables */
   int                   stp_type;           /**< Steiner problem type */
   int                   minelims;           /**< minimal number of eliminations per reduction method */
   SCIP_Bool             emitgraph;          /**< emitgraph */
   SCIP_Bool             copy;               /**< is this the problem data of a copy/sub-MIP? */
   SCIP_Bool             usecyclecons;       /**< use 2-cycle inequalities for (R)PCSPG? */
   SCIP_Bool             usesymcons;         /**< use symmetry inequalities for PCSPG and MWCSP? */
   SCIP_Bool             graphHasVanished;   /**< has the graph vanished? */
   FILE*                 logfile;            /**< logfile for DIMACS challenge */
   FILE**                origlogfile;        /**< pointer to original problem data logfile pointer */
   FILE*                 intlogfile;         /**< logfile printing all intermediate solutions for DIMACS challenge */
   FILE**                origintlogfile;     /**< pointer to original problem data intlogfile pointer */

   /** for FiberSCIP **/
   SCIP_Bool             ug;                 /**< indicates if this ug dual bound is set or not */
   int                   nSolvers;           /**< the number of solvers */
   SCIP_Real             ugDual;             /**< dual bound set by ug */
};

/**@name Local methods
 *
 * @{
 */


static
void updateorgsol(
   GRAPH*                graph,              /**< graph data structure */
   IDX*                  curr,               /**< head of solution edge list */
   STP_Bool*             orgnodes,           /**< array to mark whether a node is part of the original solution */
   STP_Bool*             orgedges,           /**< array to mark whether an edge is part of the original solution */
   int*                  nsolnodes,          /**< pointer to store the number of nodes in the original solution */
   int*                  nsoledges           /**< pointer to store the number of edges in the original solution */
)
{
   int i;
   int e = 0;
   int k = 0;

   while( curr != NULL )
   {
      i = curr->index;

      if( orgedges[i] == FALSE )
      {
         orgedges[i] = TRUE;
         e++;
      }

      if( !(orgnodes[graph->orgtail[i]]) )
      {
         k++;
         orgnodes[graph->orgtail[i]] = TRUE;;
      }

      if( !(orgnodes[graph->orghead[i]]) )
      {
         k++;
         orgnodes[graph->orghead[i]] = TRUE;
      }

      curr = curr->parent;
   }

   (*nsolnodes) += k;
   (*nsoledges) += e;
}

static
void writeCommentSection(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const char*           filename            /**< problem file name */
)
{
   char probtype[16];

   SCIPprobdataWriteLogLine(scip, "SECTION Comment\n");
   SCIPprobdataWriteLogLine(scip, "Name %s\n", filename);

   switch( graph->stp_type )
   {
   case STP_SPG:
      strcpy(probtype, "SPG");
      break;
   case STP_SAP:
      strcpy(probtype, "SAP");
      break;

   case STP_PCSPG:
      strcpy(probtype, "PCSPG");
      break;

   case STP_RPCSPG:
      strcpy(probtype, "RPCST");
      break;

   case STP_NWSPG:
      strcpy(probtype, "NWSPG");
      break;

   case STP_DCSTP:
      strcpy(probtype, "DCST");
      break;

   case STP_RSMT:
      strcpy(probtype, "RSMT");
      break;

   case STP_OARSMT:
      strcpy(probtype, "OARSMT");
      break;

   case STP_MWCSP:
      strcpy(probtype, "MWCS");
      break;

   case STP_RMWCSP:
      strcpy(probtype, "RMWCS");
      break;

   case STP_BRMWCSP:
      strcpy(probtype, "BRMWCS");
      break;

   case STP_DHCSTP:
      strcpy(probtype, "HCDST");
      break;

   default:
      strcpy(probtype, "UNKNOWN");
   }
   SCIPprobdataWriteLogLine(scip, "Problem %s\n", probtype);
   SCIPprobdataWriteLogLine(scip, "Program SCIP-Jack\n");
   SCIPprobdataWriteLogLine(scip, "Version %s\n", VERSION_SCIPJACK);
   SCIPprobdataWriteLogLine(scip, "End\n");
   SCIPprobdataWriteLogLine(scip, "\n");
   SCIPprobdataWriteLogLine(scip, "SECTION Solutions\n");
}

/*
 * distinguishes a terminal as the root; with centertype
 *      = CENTER_OK  : Do nothing
 *      = CENTER_DEG : find maximum degree
 *      = CENTER_SUM : find the minimum distance sum
 *      = CENTER_MIN : find the minimum largest distance
 *      = CENTER_ALL : find the minimum distance sum to all knots
 */
static
SCIP_RETCODE findCentralTerminal(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  central_term,       /**< pointer to store the selected (terminal) vertex */
   int                   centertype          /**< type of root selection */
   )
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

   *central_term = g->source;

   if( centertype == STP_CENTER_OK || g->grad[g->source] == 0)
      return SCIP_OKAY;

   /* Find knot of maximum degree.
    */
   if( centertype == STP_CENTER_DEG )
   {
      degree = 0;

      for( i = 0; i < g->knots; i++ )
      {
         if( g->stp_type == STP_NWPTSPG && graph_nw_knotIsLeaf(g, i) )
            continue;

         if( Is_term(g->term[i]) && (g->grad[i] > degree) )
         {
            degree = g->grad[i];
            center = i;
         }
      }

      assert(degree > 0);
      *central_term = center;
      return SCIP_OKAY;
   }

   /* For the other methods we need the shortest paths */
   SCIP_CALL( SCIPallocBufferArray(scip, &path, g->knots) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, g->edges) );

   assert(path != NULL);
   assert(cost != NULL);

   for( i = 0; i < g->knots; i++ )
      g->mark[i] = TRUE;

   for( i = 0; i < g->edges; i++ )
      cost[i] = 1.0;

   for( i = 0; i < g->knots; i++ )
   {
      if (!Is_term(g->term[i]))
         continue;

      if( g->stp_type == STP_NWPTSPG && graph_nw_knotIsLeaf(g, i) )
         continue;

      graph_path_exec(scip, g, FSP_MODE, i, cost, path);

      sum = 0.0;
      max = 0.0;

      for( k = 0; k < g->knots; k++ )
      {
         assert((path[k].edge >= 0) || (k == i));
         assert((path[k].edge >= 0) || (path[k].dist == 0));

         if( Is_term(g->term[k]) || (centertype == STP_CENTER_ALL) )
         {
            sum += path[k].dist;

            if( path[k].dist > max )
               max = path[k].dist;
         }
      }

      if( (centertype == STP_CENTER_SUM) || (centertype == STP_CENTER_ALL) )
      {
         if( sum < minimum )
         {
            minimum = sum;
            center  = i;
         }
         if( sum > maximum )
            maximum = sum;

         if( i == g->source )
            oldval = sum;
      }
      else
      {
         assert(centertype == STP_CENTER_MIN);

         /* If the maximum distance to terminal ist shorter or if
          * it is of the same length but the degree of the knot is
          * higher, we change the center.
          */
         if( SCIPisLT(scip, max, minimum) || (SCIPisEQ(scip, max, minimum) && (g->grad[i] > degree)) )
         {
            minimum = max;
            center  = i;
            degree  = g->grad[i];
         }
         if( max > maximum )
            maximum = max;

         if( i == g->source )
            oldval = max;
      }
   }
   assert(center >= 0);
   assert(Is_term(g->term[center]));

   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &path);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Central Terminal is %d (min=%g, max=%g, old=%g)\n",
      center, minimum, maximum, oldval);

   *central_term = center;
   return SCIP_OKAY;
}

/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   GRAPH*                graph               /**< graph data structure */
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
      (*probdata)->stp_type = STP_SPG;
   (*probdata)->copy = FALSE;
   (*probdata)->logfile = NULL;
   (*probdata)->origlogfile = NULL;
   (*probdata)->intlogfile = NULL;
   (*probdata)->origintlogfile = NULL;
   (*probdata)->nedges = -1;
   (*probdata)->nnodes = -1;
   (*probdata)->nterms = -1;

   (*probdata)->ug = FALSE;
   (*probdata)->nSolvers = 0;
   (*probdata)->ugDual = 0.0;

   return SCIP_OKAY;
}


/** sets STP solving mode */
static
SCIP_RETCODE setStpSolvingMode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< pointer to problem data */
)
{
   char mode;

   SCIP_CALL( SCIPgetCharParam(scip, "stp/mode", &mode) );

   /* set STP solving mode */
   probdata->mode = STP_MODE_CUT;
   assert(mode != 'p' && "pricing mode currently not supported\n");

   if( mode == 'f' )
      probdata->mode = STP_MODE_FLOW;
   else
      assert(mode == 'c');

   return SCIP_OKAY;
}


/** presolves STP */
static
SCIP_RETCODE presolveStp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< pointer to problem data */
   SCIP_Real*            offset              /**< offset from STP reductions */
)
{
   GRAPH* packedgraph;
   GRAPH* graph = probdata->graph;
   SCIP_Real oldtimelimit;
   SCIP_Real presoltimelimit;
   int reduction;

   assert(graph);

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &oldtimelimit) );
   SCIP_CALL( SCIPgetRealParam(scip, "stp/pretimelimit", &presoltimelimit) );
   SCIP_CALL( SCIPgetIntParam(scip, "stp/reduction", &reduction) );

   if( presoltimelimit > -0.5 )
      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", presoltimelimit) );

   /* save original root */
   if( !graph_pc_isPcMw(graph) )
      graph->orgsource = graph->source;

   probdata->norgedges = graph->edges;

//#define UNIT_TEST_STP
#ifdef UNIT_TEST_STP
   SCIP_CALL( stptest_all(scip) );
#endif

   /* the actual presolving */
   SCIP_CALL( reduce(scip, graph, offset, reduction, probdata->minelims, TRUE) );

#ifdef STP_WRITE_RED_STATS
   graph_writeReductionStats(graph, SCIPgetProbName(scip), "~/redstats.txt");
   exit(1);
#endif

#ifdef WITH_UG
   SCIP_CALL( graph_pack(scip, graph, &packedgraph, offset, FALSE) );
#else
   SCIP_CALL( graph_pack(scip, graph, &packedgraph, offset, TRUE) );
#endif

   graph = packedgraph;
   assert(graph);

   probdata->stp_type = graph->stp_type;
   probdata->graph = graph;

   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", oldtimelimit) );

   return SCIP_OKAY;
}


/** frees the constraints of CUT model */
static
SCIP_RETCODE freeConstraintsCutModel(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< pointer to problem data */
   )
{
   assert(probdata->mode == STP_MODE_CUT);

   /* Degree-Constrained STP? */
   if( probdata->stp_type == STP_DCSTP )
   {
      /* release degree constraints */
      for( int t = 0; t < probdata->nnodes; ++t)
         SCIP_CALL( SCIPreleaseCons(scip, &(probdata->degcons[t])) );

      SCIPfreeMemoryArrayNull(scip, &(probdata->degcons));
   }

   if( probdata->graphHasVanished )
      return SCIP_OKAY;

   /* PC variant STP? */
   if( probdata->stp_type == STP_PCSPG || probdata->stp_type == STP_RPCSPG || probdata->stp_type == STP_MWCSP )
   {
#if FLOWB
      for( e = 0; e < probdata->nnonterms; e++ )
         SCIP_CALL( SCIPreleaseCons(scip, &(probdata->flowbcons[e])) );

      SCIPfreeMemoryArrayNull(scip, &(probdata->flowbcons));
#endif

      if( probdata->usecyclecons && probdata->stp_type == STP_PCSPG  )
      {
         /* release degree constraints */
         for( int e = 0; e < probdata->nedges; e++ )
            SCIP_CALL( SCIPreleaseCons(scip, &(probdata->prizecyclecons[e])) );

         SCIPfreeMemoryArrayNull(scip, &(probdata->prizecyclecons));
      }

      if( probdata->stp_type != STP_RPCSPG )
      {
         if( probdata->usesymcons )
         {
#ifdef STP_AGG_SYM
            const int e = probdata->realnterms - 1;

            for( int t = 0; t < e; ++t )
               SCIP_CALL( SCIPreleaseCons(scip, &(probdata->prizesymcons[t])) );
#else
            const int e = ((probdata->realnterms - 1) * probdata->realnterms) / 2;

            for( int t = 0; t < e; ++t )
               SCIP_CALL( SCIPreleaseCons(scip, &(probdata->prizesymcons[t])) );
#endif

            SCIPfreeMemoryArrayNull(scip, &(probdata->prizesymcons));
         }
         SCIP_CALL( SCIPreleaseCons(scip, &(probdata->prizecons)) );
      }
   }

   return SCIP_OKAY;
}


/** frees the constraints of non-CUT model */
static
SCIP_RETCODE freeConstraintsNonCutModel(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< pointer to problem data */
   )
{
   assert(probdata->mode != STP_MODE_CUT);

   /* release path constraints */
   if( probdata->mode == STP_MODE_PRICE )
   {
      for( int t = 0; t < probdata->realnterms; ++t )
         SCIP_CALL( SCIPreleaseCons(scip, &(probdata->pathcons[t])) );
   }
   else
   {
      assert(probdata->mode == STP_MODE_FLOW);

      for( int t = 0; t < probdata->realnterms; ++t )
      {
         /* release constraints and variables */
         for( int e = 0; e < probdata->nnodes - 1; ++e )
            SCIP_CALL( SCIPreleaseCons(scip, &(probdata->pathcons[t * (probdata->nnodes - 1 ) + e])) );

         for( int e = 0; e < probdata->nedges; ++e )
            SCIP_CALL( SCIPreleaseVar(scip, &probdata->flowvars[t * probdata->nedges + e]) );

         if( probdata->bigt )
            break;
      }

      SCIPfreeMemoryArrayNull(scip, &probdata->flowvars);
   }

   /* release edge constraints of both Price and Flow */
   for( int t = 0; t < probdata->realnterms; ++t )
   {
      for( int e = 0; e < probdata->nedges; ++e )
         SCIP_CALL( SCIPreleaseCons(scip, &(probdata->edgecons[t * probdata->nedges + e])));

      if( probdata->bigt )
         break;
   }

   SCIPfreeMemoryArrayNull(scip, &(probdata->edgecons));
   SCIPfreeMemoryArrayNull(scip, &(probdata->pathcons));

   return SCIP_OKAY;
}


/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
   SCIPdebugMessage("probdataFree \n");

   assert(scip != NULL);
   assert(probdata != NULL);

   for( int e = 0; e < (*probdata)->nvars; ++e )
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->edgevars[e]) );

   SCIPfreeMemoryArrayNull(scip, &(*probdata)->edgevars);

#if USEOFFSETVAR
   if( (*probdata)->offsetvar != NULL )
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->offsetvar) );
#endif

   if( (*probdata)->mode == STP_MODE_CUT )
   {
      SCIP_CALL( freeConstraintsCutModel(scip, *probdata) );
   }
   else
   {
      SCIP_CALL( freeConstraintsNonCutModel(scip, *probdata) );
   }

   if( (*probdata)->stp_type == STP_DHCSTP )
      SCIP_CALL( SCIPreleaseCons(scip, &((*probdata)->hopcons)) );

   if( (*probdata)->stp_type == STP_BRMWCSP )
      SCIP_CALL( SCIPreleaseCons(scip, &((*probdata)->budgetcons)) );

   SCIPfreeMemoryArrayNull(scip, &(*probdata)->xval);
   SCIPfreeMemoryArrayNull(scip, &(*probdata)->realterms);
   SCIPfreeMemoryArrayNull(scip, &(*probdata)->pctermsorder);

   /* free probdata */
   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


/** create (edge-) HOP constraint (cut mode only) */
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

   SCIPdebugMessage("createHopeConstraint \n");
   graph = probdata->graph;
   assert(graph != NULL);
   rhs = graph->hoplimit;
   SCIPdebugMessage("Hop limit: %f \n ", rhs);

   /* @todo with presolving contractions enabled one might set rhs = rhs - (number of fixed edges) */
   SCIP_CALL( SCIPcreateConsLinear ( scip, &(probdata->hopcons), "HopConstraint", 0, NULL, NULL,
         -SCIPinfinity(scip), rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPaddCons(scip, probdata->hopcons) );

   return SCIP_OKAY;
}


static
SCIP_RETCODE createBudgetConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   GRAPH* graph;
   assert(scip != NULL);
   assert(probdata != NULL);

   SCIPdebugMessage("create Budget constraint \n");
   graph = probdata->graph;
   assert(graph != NULL);
   assert(graph->budget > 0.0);

   SCIP_CALL( SCIPcreateConsLinear ( scip, &(probdata->budgetcons), "BudgetConstraint", 0, NULL, NULL,
         -SCIPinfinity(scip), graph->budget, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPaddCons(scip, probdata->budgetcons) );

   return SCIP_OKAY;
}

/** create (node-) degree constraints (cut mode only) */
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

   SCIPdebugMessage("createDegreeConstraints \n");
   graph = probdata->graph;
   nnodes = probdata->nnodes;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->degcons), nnodes ) );

   for( k = 0; k < nnodes; ++k )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "DegreeConstraint%d", k);
      SCIP_CALL( SCIPcreateConsLinear(scip, &(probdata->degcons[k]), consname, 0, NULL, NULL,
            -SCIPinfinity(scip), (SCIP_Real)graph->maxdeg[k], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(scip, probdata->degcons[k]) );
   }

   return SCIP_OKAY;
}


/** create Prize constraints (cut mode only) */
static
SCIP_RETCODE createPrizeConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )

{
   GRAPH* graph;
   int r;
   int ro2;
   int root;
   int nedges;
   int realnterms;
   char consname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(probdata != NULL);

   graph = probdata->graph;
   realnterms = probdata->realnterms;

   assert(graph != NULL);
   assert(graph->stp_type != STP_RPCSPG);

   root = graph->source;
   nedges = graph->edges;

   SCIPdebugMessage("createPrizeConstraints \n");

   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PrizeConstraint");
   SCIP_CALL( SCIPcreateConsLinear ( scip, &(probdata->prizecons), consname, 0, NULL, NULL, 1.0, 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
   SCIP_CALL( SCIPaddCons(scip, probdata->prizecons) );

   if( probdata->usesymcons )
   {
#ifdef STP_AGG_SYM
      ro2 = realnterms - 1;
#else
      ro2 = (realnterms * (realnterms - 1)) / 2;
#endif
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->prizesymcons), ro2) );

      for( r = 0; r < ro2; r++ )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PrizeSymConstraint%d", r);
#ifndef STP_SYM_PRIZE
         SCIP_CALL( SCIPcreateConsSetpack(scip, &(probdata->prizesymcons[r]), consname, 0, NULL, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ));
#else
         SCIP_CALL( SCIPcreateConsLinear ( scip, &(probdata->prizesymcons[r]), consname, 0, NULL, NULL,
                    -SCIPinfinity(scip), 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
#endif

         SCIP_CALL(SCIPaddCons(scip, probdata->prizesymcons[r]));
      }
   }

   if( probdata->usecyclecons && graph->stp_type == STP_PCSPG )
   {
      ro2 = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->prizecyclecons), nedges) );
      for( r = 0; r < nedges; r++ )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PrizeLPConstraint%d", ro2);
         if( root == graph->head[r] )
         {
            SCIP_CALL( SCIPcreateConsLinear ( scip, &(probdata->prizecyclecons[ro2]), consname, 0, NULL, NULL,
                  -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         }
         else
         {
            SCIP_CALL( SCIPcreateConsLinear ( scip, &(probdata->prizecyclecons[ro2]), consname, 0, NULL, NULL,
                  -SCIPinfinity(scip), 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         }
         SCIP_CALL( SCIPaddCons(scip, probdata->prizecyclecons[ro2++]) );
      }
   }

   return SCIP_OKAY;
}


#if FLOWB
/** create Flow Balance constraints (cut mode only) */
static
SCIP_RETCODE createFlowBalanceConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )

{
   GRAPH* graph;
   int r;
   int nnonterms;
   char consname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(probdata != NULL);

   graph = probdata->graph;

   assert(graph != NULL);

   if( graph->stp_type == STP_RPCSPG )
      nnonterms = graph->knots - graph->terms;
   else
      nnonterms = graph->knots - graph->terms;

   probdata->nnonterms = nnonterms;

   SCIPdebugMessage("createFlowBalanceConstraints \n");


   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->flowbcons), nnonterms) );
   for( r = 0; r < nnonterms; r++ )
   {
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "FlowBalanceConstraint%d", r);
      SCIP_CALL( SCIPcreateConsLinear ( scip, &(probdata->flowbcons[r]), consname, 0, NULL, NULL,
            -SCIPinfinity(scip), 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, probdata->flowbcons[r]) );
   }

   return SCIP_OKAY;
}

#endif

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
   assert(probdata->mode != STP_MODE_CUT);

   SCIPdebugMessage("createConstraints \n");
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
	          -SCIPinfinity(scip), 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, probdata->mode == STP_MODE_PRICE, FALSE, FALSE, FALSE) );
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
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, probdata->mode == STP_MODE_PRICE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, probdata->edgecons[e]) );
      }
   }

   /* Branch and Price mode */
   if( probdata->mode == STP_MODE_PRICE )
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
   else if( probdata->mode == STP_MODE_FLOW )
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
               if( k !=  graph->source)
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
            if( k != graph->source)
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
   SCIP_VAR* var;
   SCIP_Real* edgecost;
   int e;
   int t;
   int k;
   int k2;
   int tail;
   char varname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(probdata != NULL);

   t = 0;
   graph = probdata->graph;
   SCIPdebugMessage("createVariables \n");

   /* nontrivial problem? */
   if( !probdata->graphHasVanished )
   {
      int nedges = probdata->nedges;
      int nvars = probdata->nvars;
      int realnterms = probdata->realnterms;
      int root = graph->source;
      int nnodes = graph->knots;
      SCIP_Bool objint = SCIPisIntegral(scip, offset);

      assert(nedges == graph->edges);
      assert(realnterms >= 1);

      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->xval, nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->edgevars, nvars) );

      /* cut mode */
      if( probdata->mode == STP_MODE_CUT )
      {
         int a;
         assert(probdata->nlayers == 1);

         for( e = 0; e < nedges; ++e )
         {
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d_%d", graph->tail[e], graph->head[e]);
            SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->edgevars[e], varname, 0.0, 1.0, graph->cost[e], SCIP_VARTYPE_BINARY) );
            SCIP_CALL( SCIPaddVar(scip, probdata->edgevars[e]) );
            objint = objint && SCIPisIntegral(scip, graph->cost[e]);
         }

         /* Hop-Constrained STP */
         if( graph->stp_type == STP_DHCSTP )
         {
            SCIP_Real hopfactor;
            for( e = 0; e < nedges; ++e )
            {
               /* @note: When contractions are used in presolving: modify */
               hopfactor = 1.0;
               SCIP_CALL( SCIPaddCoefLinear(scip, probdata->hopcons, probdata->edgevars[e], hopfactor) );
            }
         }

         if( graph->stp_type == STP_BRMWCSP )
         {
            for( k = 0; k < nnodes; ++k )
            {
               const SCIP_Real kbudget = graph->costbudget[k];
               assert(kbudget >= 0.0);
               assert(graph->extended);

               if( graph_pc_knotIsDummyTerm(graph, k)  )
                  continue;

               if( SCIPisZero(scip, kbudget) )
                  continue;

               for( e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
                  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->budgetcons, probdata->edgevars[e], kbudget) );
            }
         }

         /* Degree-Constrained STP */
         if( graph->stp_type == STP_DCSTP )
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

#if FLOWB
            k2 = 0;

            for( k = 0; k < nnodes; k++ )
            {
               if( !Is_term(graph->term[k]) )
               {
                  /* incoming arcs */
                  for( a = graph->inpbeg[k]; a != EAT_LAST; a = graph->ieat[a] )
                  {
                     SCIP_CALL( SCIPaddCoefLinear(scip, probdata->flowbcons[k2], probdata->edgevars[a], 1.0) );
                  }
                  /* outgoing arcs */
                  for( a = graph->outbeg[k]; a != EAT_LAST; a = graph->oeat[a] )
                  {
                     SCIP_CALL( SCIPaddCoefLinear(scip, probdata->flowbcons[k2], probdata->edgevars[a], -1.0) );
                  }
                  k2++;
               }
            }
#endif

         if( graph_pc_isPcMw(graph) )
         {
            for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
            {
               const int head = graph->head[e];
               if( Is_term(graph->term[head]) && !graph_pc_knotIsFixedTerm(graph, head) )
               {
                  assert(graph->grad[head] == 2);
                  assert(graph->prize != NULL);

                  /* variables are preferred to be branched on */
                  SCIP_CALL(SCIPchgVarBranchPriority(scip, probdata->edgevars[e], 100 + (int)(100.0 * graph->prize[head])));
               }
            }
         }

         if( graph->stp_type == STP_PCSPG || graph->stp_type == STP_MWCSP )
         {
            int* pseudoterms = NULL;

            assert(graph->extended);
            assert(nnodes == probdata->nnodes);

            SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->pctermsorder), nnodes) );

            t = 0;
            for( int v = 0; v < nnodes; v++ )
            {
               if( Is_pseudoTerm(graph->term[v]) )
                  probdata->pctermsorder[v] = t++;
               else
                  probdata->pctermsorder[v] = nnodes;
            }
            assert(t == probdata->realnterms);

            if( probdata->usesymcons )
               SCIP_CALL( SCIPallocBufferArray(scip, &pseudoterms, probdata->realnterms) );

            if( probdata->usecyclecons && graph->stp_type == STP_PCSPG )
            {
               for( e = 0; e < nedges; e++ )
               {
                  const int head = graph->head[e];
                  SCIP_CALL(SCIPaddCoefLinear(scip, probdata->prizecyclecons[e], probdata->edgevars[e], 1.0));
                  SCIP_CALL(SCIPaddCoefLinear(scip, probdata->prizecyclecons[e], probdata->edgevars[flipedge(e)], 1.0));
                  if( root != head )
                  {
                     for( a = graph->inpbeg[head]; a != EAT_LAST; a = graph->ieat[a] )
                     {
                        SCIP_CALL(SCIPaddCoefLinear(scip, probdata->prizecyclecons[e], probdata->edgevars[a], -1.0));
                     }
                  }
               }
            }

            t = 0;
            k2 = 0;
            for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
            {
               const int head = graph->head[e];
               if( !Is_term(graph->term[head]) )
               {
                  SCIP_CALL(SCIPaddCoefLinear(scip, probdata->prizecons, probdata->edgevars[e], 1.0));

                  assert(graph->prize != NULL);
                  if( probdata->usesymcons )
                  {
                     assert(pseudoterms != NULL);

                     pseudoterms[t] = head;
#ifdef STP_AGG_SYM
#ifndef STP_SYM_PRIZE
                     for( int prev = 0; prev < t; prev++ )
                        SCIP_CALL( SCIPaddCoefSetppc(scip, probdata->prizesymcons[prev], probdata->edgevars[e]) );

                     /* skip the last term */
                     if( t == graph->terms - 2 )
                     {
                        assert(t == probdata->realnterms - 1);
                        continue;
                     }

                     for( a = graph->inpbeg[head]; a != EAT_LAST; a = graph->ieat[a] )
                        SCIP_CALL( SCIPaddCoefSetppc(scip, probdata->prizesymcons[t], probdata->edgevars[a]) );
#endif
#else
                     for( k = 0; k < t; k++ )
                     {
                        for( a = graph->inpbeg[pseudoterms[k]]; a != EAT_LAST; a = graph->ieat[a] )
                        {
                           SCIP_CALL(SCIPaddCoefSetppc(scip, probdata->prizesymcons[k2], probdata->edgevars[a]));
                        }

                        SCIP_CALL(SCIPaddCoefSetppc(scip, probdata->prizesymcons[k2], probdata->edgevars[e]));
                        k2++;
                     }
#endif
                     t++;
                  }
               }
            }

            if( probdata->usesymcons )
            {
#ifdef STP_SYM_PRIZE
               SCIP_Real* termprizes;
               SCIP_CALL( SCIPallocBufferArray(scip, &termprizes, probdata->realnterms) );

               assert(t == graph->terms - 1);
               assert(probdata->realnterms == t);

               for( int s = 0; s < probdata->realnterms; s++ )
               {
                  const int pterm = pseudoterms[s];
                  assert(Is_pseudoTerm(graph->term[pterm]));
                  assert(graph->prize[pterm] > 0.0);

                  termprizes[s] = graph->prize[pterm];
               }
               SCIPsortRealInt(termprizes, pseudoterms, probdata->realnterms);
#ifndef NDEBUG
               for( int s = 1; s < probdata->realnterms; s++ )
                  assert(graph->prize[pseudoterms[s - 1]] <= graph->prize[pseudoterms[s]]);
#endif

               for( int s = 0; s < probdata->realnterms; s++ )
               {
                  const int pterm = pseudoterms[s];
                  const int twin = graph_pc_getTwinTerm(graph, pterm);
                  const int twinrootedge = graph_pc_getRoot2PtermEdge(graph, twin);
                  const int prootedge = graph_pc_getRoot2PtermEdge(graph, pterm);

                  assert(graph->prize[pterm] == graph->cost[twinrootedge]);
                  assert(graph->cost[prootedge] == 0.0);

                  for( int prev = 0; prev < s; prev++ )
                     SCIP_CALL( SCIPaddCoefLinear(scip, probdata->prizesymcons[prev], probdata->edgevars[prootedge], 1.0) );

                  if( s == probdata->realnterms - 1 )
                     break;

                  SCIP_CALL(SCIPaddCoefLinear(scip, probdata->prizesymcons[s], probdata->edgevars[twinrootedge], -1.0));
#if 0
                  for( a = graph->inpbeg[pterm]; a != EAT_LAST; a = graph->ieat[a] )
                     SCIP_CALL( SCIPaddCoefSetppc(scip, probdata->prizesymcons[s], probdata->edgevars[a]) );
#endif
               }

               for( int p = 0; p < probdata->realnterms; p++ )
               {
                  const int index = pseudoterms[p];
                  assert(probdata->pctermsorder[index] >= 0);

                  probdata->pctermsorder[index] = p;
               }

               SCIPfreeBufferArray(scip, &termprizes);
#endif
               SCIPfreeBufferArray(scip, &pseudoterms);
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

            if( !probdata->bigt )
            {
               for( t = 0; t < realnterms; ++t )
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[t * nedges + e], var, -1.0) );
               }
            }
            else
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[e], var, (double ) -realnterms));
            }
            probdata->edgevars[e] = var;
         }
      }

      /* Price mode */
      if( probdata->mode == STP_MODE_PRICE )
      {
         PATH* path;

         /* the flow variables are not used in the Price mode */
         probdata->flowvars = NULL;

         /* compute shortest paths to the root */
         SCIP_CALL(SCIPallocBufferArray(scip, &path, graph->knots));
         SCIP_CALL(SCIPallocBufferArray(scip, &edgecost, nedges));
         for (e = 0; e < nedges; ++e)
            edgecost[e] = graph->cost[e];

         for (e = 0; e < graph->knots; e++)
            graph->mark[e] = 1;

         graph_path_exec(scip, graph, FSP_MODE, root, edgecost, path);

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
                  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[t * nedges + path[tail].edge], var, 1.0));
               }
               else
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, probdata->edgecons[path[tail].edge], var, 1.0));
               }

               tail = graph->tail[path[tail].edge];
            }
         }

         /* free local arrays */
         SCIPfreeBufferArray(scip, &edgecost);
         SCIPfreeBufferArray(scip, &path);
      }
      /* Flow mode */
      else if( probdata->mode == STP_MODE_FLOW )
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

#if USEOFFSETVAR
   /* add offset */
   (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "OFFSET");
   SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->offsetvar, varname, 1.0, 1.0, offset, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, probdata->offsetvar) );
#else
   SCIP_CALL( SCIPaddOrigObjoffset(scip, offset) );
#endif

   return SCIP_OKAY;
}


/** creates the actual MIP/IP model */
static
SCIP_RETCODE createModel(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata);

   if( !probdata->graphHasVanished )
   {
      GRAPH* const graph = probdata->graph;

      const int nnodes = graph_get_nNodes(graph);
      const int nedges = graph_get_nEdges(graph);
      /* compute the real number of terminals (basically nterms - 1) */
      const int realnterms = graph->terms - 1;
      const SCIP_Bool mw = (graph->stp_type == STP_MWCSP);
      const SCIP_Bool pc = (graph->stp_type == STP_PCSPG);

      assert(graph->terms > 1);

      SCIP_CALL(graph_path_init(scip, graph));
#ifdef WITH_UG
       SCIP_CALL( graph_path_init(scip, probdata->orggraph) );
 #endif

      if( probdata->mode == STP_MODE_CUT )
      {
         SCIP_CALL(graph_mincut_init(scip, graph));
      }

      probdata->nvars = nedges;
      probdata->nnodes = nnodes;
      probdata->nedges = nedges;
      probdata->nterms = graph->terms;
      probdata->nlayers = graph->layers;
      probdata->realnterms = realnterms;

      assert(Is_term(graph->term[graph->source]));

      SCIP_CALL(SCIPallocMemoryArray(scip, &probdata->realterms, realnterms));

      for( int k = 0, t = 0; k < nnodes; ++k )
      {
         if( graph->term[k] == 0 && k != graph->source )
         {
            probdata->realterms[t] = k;
            SCIPdebugMessage("realterms[%d] = %d \n ", t, probdata->realterms[t]);
            t++;
         }
      }

      if( probdata->mode == STP_MODE_CUT )
      {
         SCIP_CONS* cons;

         /* create and add constraint for Branch and Cut */
         SCIP_CALL(SCIPcreateConsStp(scip, &cons, "stpcons", probdata->graph));
         SCIP_CALL(SCIPaddCons(scip, cons));
         SCIP_CALL(SCIPreleaseCons(scip, &cons));

         /* if the problem is a HOP-Constrained-STP, an additional constraint is required */
         if( graph->stp_type == STP_DHCSTP )
         {
            SCIP_CALL(createHopConstraint(scip, probdata));
         }

         /* if the problem is a Degree-Constrained-STP, additional constraints are required */
         if( graph->stp_type == STP_DCSTP )
         {
            SCIP_CALL(createDegreeConstraints(scip, probdata));
         }

         if( graph->stp_type == STP_BRMWCSP )
         {
            SCIP_CALL(createBudgetConstraint(scip, probdata));
         }

         /* if the problem is a Prize-Collecting-STP or a Maximum Weight Connected Subgraph Problem, additional constraints are required */
         if( pc || mw )
         {
            SCIP_CALL(createPrizeConstraints(scip, probdata));
         }
#if FLOWB
          SCIP_CALL( createFlowBalanceConstraints(scip, probdata) );
 #endif
      }
      else
      {
         /* create and add constraints for Flow or Branch and Price */
         SCIP_CALL(createConstraints(scip, probdata));
      }
   }
   /* graph reduction solved the whole problem, set vars to zero or NULL */
   else
   {
      probdata->pathcons = NULL;
      probdata->edgecons = NULL;
      probdata->nlayers = 0;
      probdata->nnodes = 1;
      probdata->nedges = 0;
      probdata->nterms = 1;
      probdata->nnonterms = 0;
      probdata->realnterms = 0;
      probdata->nedges = 0;
      probdata->nvars = 0;
      probdata->realterms = NULL;
      probdata->xval = NULL;
   }
   probdata->pctermsorder = NULL;

   /* create and add initial variables */
   SCIP_CALL( createVariables(scip, probdata, probdata->offset) );

   return SCIP_OKAY;
}


/** creates a log file */
static
SCIP_RETCODE createLogfile(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   const char*           intlogfilename,     /**< user parameter */
   const char*           logfilename,        /**< user parameter */
   const char*           filename,           /**< file name */
   const char*           probname            /**< problem name */
   )
{

#ifdef WITH_UG
   if( getUgRank() != 0 )
   {
      logfilename = NULL;
      intlogfilename = NULL;
   }
#endif

   if( logfilename != NULL && logfilename[0] != '\0' )
   {
      char finalfilename[SCIP_MAXSTRLEN];

      if( strcmp("use_probname", logfilename) == 0 )
         (void) SCIPsnprintf(finalfilename, SCIP_MAXSTRLEN, "%s.stplog", probname);
      else
         (void) SCIPsnprintf(finalfilename, SCIP_MAXSTRLEN, "%s", logfilename);

      probdata->logfile = fopen(finalfilename, "w");

      if( probdata->logfile == NULL )
      {
         SCIPerrorMessage("cannot create file <%s> for writing\n", finalfilename);
         SCIPprintSysError(finalfilename);
         return SCIP_FILECREATEERROR;
      }
   }

   if( intlogfilename != NULL && intlogfilename[0] != '\0' )
   {
      char finalfilename[SCIP_MAXSTRLEN];

      if( strcmp("use_probname", intlogfilename) == 0 )
         (void) SCIPsnprintf(finalfilename, SCIP_MAXSTRLEN, "%s_int.stplog", probname);
      else
         (void) SCIPsnprintf(finalfilename, SCIP_MAXSTRLEN, "%s", intlogfilename);

      probdata->intlogfile = fopen(finalfilename, "w");

      if( probdata->intlogfile == NULL )
      {
         SCIPerrorMessage("cannot create file <%s> for writing\n", finalfilename);
         SCIPprintSysError(finalfilename);
         return SCIP_FILECREATEERROR;
      }
   }

   return SCIP_OKAY;
}

/** creates initial cuts for SCIP*/
static
SCIP_RETCODE createInitialCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   GRAPH* const graph = probdata->graph;
   SCIP_Real lpobjval;
   const SCIP_Bool mw = (graph->stp_type == STP_MWCSP);
   const SCIP_Bool pc = (graph->stp_type == STP_PCSPG);

   assert(STP_MODE_CUT == probdata->mode);

   if( pc || mw )
   {
      SCIP_CALL( SCIPStpDualAscentPcMw(scip, graph, NULL, &lpobjval, TRUE, TRUE, 1) );
   }
   else
   {
      if( graph->stp_type != STP_RPCSPG && graph->stp_type != STP_SPG && graph->stp_type != STP_RSMT && graph->stp_type != STP_OARSMT && graph->stp_type != STP_GSTP )
      {
         SCIP_CALL( SCIPStpDualAscent(scip, graph, NULL, &lpobjval, TRUE, FALSE, NULL, NULL, NULL, NULL, graph->source, FALSE, -1.0) );
      }
      else
      {
         SCIP_CALL( SCIPStpDualAscent(scip, graph, NULL, &lpobjval, TRUE, TRUE, NULL, NULL, NULL, NULL, graph->source, FALSE, -1.0) );
      }
   }

   return SCIP_OKAY;
}

/** set parameters */
static
SCIP_RETCODE setParams(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   symcons,            /**< user parameter */
   int                   cyclecons           /**< user parameter */
   )
{
   GRAPH* const graph = probdata->graph;
   const SCIP_Bool mw = (graph->stp_type == STP_MWCSP);
   const SCIP_Bool pc = (graph->stp_type == STP_PCSPG);

#ifndef WITH_UG
   if( (graph->edges > CUT_MAXTOTNEDGES) || ((graph->edges > CUT_MAXNEDGES) && (graph->terms > CUT_MAXNTERMINALS)) )
   {
      SCIP_CALL(SCIPsetIntParam(scip, "separating/aggregation/maxroundsroot", 3));
      SCIP_CALL(SCIPsetIntParam(scip, "separating/strongcg/maxroundsroot", 3));
      SCIP_CALL(SCIPsetIntParam(scip, "separating/gomory/maxroundsroot", 3));
      SCIP_CALL(SCIPsetIntParam(scip, "separating/zerohalf/maxroundsroot", 3));
   }
#endif

   SCIP_CALL(SCIPsetIntParam(scip, "separating/clique/freq", -1));

   if( pc )
   {
      SCIP_CALL(SCIPsetIntParam(scip, "constraints/logicor/presoltiming", 4));
   }

   if( graph->stp_type == STP_DHCSTP )
   {
      SCIP_CALL(SCIPsetIntParam(scip, "constraints/knapsack/propfreq", -1));
      SCIP_CALL(SCIPsetIntParam(scip, "constraints/logicor/propfreq", -1));
      SCIP_CALL(SCIPsetIntParam(scip, "separating/flowcover/freq", -1));
      SCIP_CALL(SCIPsetIntParam(scip, "heuristics/locks/freq", -1));
      SCIP_CALL(SCIPsetIntParam(scip, "heuristics/rounding/freq", -1));
   }

   probdata->usesymcons = FALSE;
   probdata->usecyclecons = FALSE;

   if( pc || mw )
   {
      if( cyclecons == STP_CONS_ALWAYS || (cyclecons == STP_CONS_AUTOMATIC && graph->edges <= CYC_CONS_LIMIT) )
         probdata->usecyclecons = TRUE;

      if( symcons == STP_CONS_ALWAYS || (symcons == STP_CONS_AUTOMATIC && graph->terms <= SYM_CONS_LIMIT) )
         probdata->usesymcons = TRUE;

      if( probdata->usesymcons  )
         SCIPdebugMessage("USE SYM CONS: %d \n", graph->terms);
      else
         SCIPdebugMessage("NO SYM CONS: \n");

      SCIP_CALL(SCIPsetIntParam(scip, "heuristics/trivial/freq", -1));
   }

   if( graph->knots > 1 )
   {
      /* we always to some basic presolving and thus assume that there is not problem with only one terminal
       * and additional nodes */
      assert(graph->terms > 1);

      probdata->graphHasVanished = FALSE;
   }
   else
   {
      assert(1 == graph->knots);
      assert(1 == graph->terms);
      assert(0 == graph->edges);

      probdata->graphHasVanished = TRUE;
   }

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

   SCIPdebugMessage("########################## probcopy ###########################\n");

   SCIP_CALL( graph_copy(scip, sourcedata->graph, &graphcopy) );
   SCIP_CALL( graph_path_init(scip, graphcopy) );

   if( sourcedata->mode == STP_MODE_CUT )
      SCIP_CALL( graph_mincut_init(scip, graphcopy) );

   SCIP_CALL( probdataCreate(scip, targetdata, graphcopy) );

#ifdef WITH_UG
   {
      GRAPH* orggraphcopy;
      SCIP_CALL( graph_copy(scip, sourcedata->orggraph, &orggraphcopy) );
      SCIP_CALL( graph_path_init(scip, orggraphcopy) );
      (*targetdata)->orggraph = orggraphcopy;
   }
#endif

   (*targetdata)->minelims = sourcedata->minelims;
   (*targetdata)->offset = sourcedata->offset;
   (*targetdata)->norgedges = sourcedata->norgedges;
   (*targetdata)->mode = sourcedata->mode;
   (*targetdata)->bigt = sourcedata->bigt;
   (*targetdata)->nlayers = sourcedata->nlayers;
   (*targetdata)->nedges = sourcedata->nedges;
   (*targetdata)->nnodes = sourcedata->nnodes;
   (*targetdata)->nterms = sourcedata->nterms;
   (*targetdata)->nnonterms = sourcedata->nnonterms;
   (*targetdata)->realnterms = sourcedata->realnterms;
   (*targetdata)->emitgraph = sourcedata->emitgraph;
   (*targetdata)->usesymcons = sourcedata->usesymcons;
   (*targetdata)->graphHasVanished = sourcedata->graphHasVanished;
   (*targetdata)->usecyclecons = sourcedata->usecyclecons;
   (*targetdata)->nvars = sourcedata->nvars;
   (*targetdata)->copy = TRUE;
#if USEOFFSETVAR
   (*targetdata)->offsetvar = NULL;

   if( sourcedata->offsetvar != NULL && SCIPvarIsActive(sourcedata->offsetvar) )
   {
      SCIP_Bool success;

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->offsetvar, &((*targetdata)->offsetvar), varmap, consmap, global, &success) );
      assert(success);

      SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->offsetvar) );
   }
#endif
   if( sourcedata->nedges > 0 )
   {
      SCIP_Bool success;
      int v;
      int c;

      /* transform variables */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->xval, sourcedata->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->edgevars, sourcedata->nvars) );

      for( v = sourcedata->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->edgevars[v], &((*targetdata)->edgevars[v]), varmap, consmap, global, &success) );
         assert(success);

         SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->edgevars[v]) );
      }

      /* cut mode */
      if( sourcedata->mode == STP_MODE_CUT )
      {
         (*targetdata)->edgecons = NULL;
         (*targetdata)->pathcons = NULL;
         if( sourcedata->stp_type == STP_DCSTP )
         {
            SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->degcons, sourcedata->nnodes));
            for( c = sourcedata->nnodes - 1; c >= 0; --c )
            {
               SCIP_CALL(
                     SCIPgetConsCopy(sourcescip, scip, sourcedata->degcons[c], &((*targetdata)->degcons[c]), SCIPconsGetHdlr(sourcedata->degcons[c]),
                           varmap, consmap, SCIPconsGetName(sourcedata->degcons[c]), SCIPconsIsInitial(sourcedata->degcons[c]),
                           SCIPconsIsSeparated(sourcedata->degcons[c]), SCIPconsIsEnforced(sourcedata->degcons[c]),
                           SCIPconsIsChecked(sourcedata->degcons[c]), SCIPconsIsPropagated(sourcedata->degcons[c]),
                           SCIPconsIsLocal(sourcedata->degcons[c]), SCIPconsIsModifiable(sourcedata->degcons[c]),
                           SCIPconsIsDynamic(sourcedata->degcons[c]), SCIPconsIsRemovable(sourcedata->degcons[c]),
                           SCIPconsIsStickingAtNode(sourcedata->degcons[c]), global, &success));
               assert(success);

               SCIP_CALL(SCIPcaptureCons(scip, (*targetdata)->degcons[c]));
            }
         }

         if( sourcedata->stp_type == STP_PCSPG || sourcedata->stp_type == STP_RPCSPG || sourcedata->stp_type == STP_MWCSP )
         {
#if FLOWB
            SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->flowbcons, sourcedata->nnonterms) );
            for( c = sourcedata->nnonterms - 1; c >= 0; --c )
            {
               SCIP_CALL( SCIPgetConsCopy(sourcescip, scip, sourcedata->flowbcons[c], &((*targetdata)->flowbcons[c]),
                           SCIPconsGetHdlr(sourcedata->flowbcons[c]), varmap, consmap,
                           SCIPconsGetName(sourcedata->flowbcons[c]),
                           SCIPconsIsInitial(sourcedata->flowbcons[c]),
                           SCIPconsIsSeparated(sourcedata->flowbcons[c]),
                           SCIPconsIsEnforced(sourcedata->flowbcons[c]),
                           SCIPconsIsChecked(sourcedata->flowbcons[c]),
                           SCIPconsIsPropagated(sourcedata->flowbcons[c]),
                           SCIPconsIsLocal(sourcedata->flowbcons[c]),
                           SCIPconsIsModifiable(sourcedata->flowbcons[c]),
                           SCIPconsIsDynamic(sourcedata->flowbcons[c]),
                           SCIPconsIsRemovable(sourcedata->flowbcons[c]),
                           SCIPconsIsStickingAtNode(sourcedata->flowbcons[c]),
                           global, &success) );
               assert(success);

               SCIP_CALL( SCIPcaptureCons(scip, (*targetdata)->flowbcons[c]) );
            }
#endif

            if( sourcedata->usecyclecons && sourcedata->stp_type == STP_PCSPG )
            {
               SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->prizecyclecons, sourcedata->nedges));
               for( c = sourcedata->nedges - 1; c >= 0; --c )
               {
                  SCIP_CALL(
                        SCIPgetConsCopy(sourcescip, scip, sourcedata->prizecyclecons[c], &((*targetdata)->prizecyclecons[c]),
                              SCIPconsGetHdlr(sourcedata->prizecyclecons[c]), varmap, consmap, SCIPconsGetName(sourcedata->prizecyclecons[c]),
                              SCIPconsIsInitial(sourcedata->prizecyclecons[c]), SCIPconsIsSeparated(sourcedata->prizecyclecons[c]),
                              SCIPconsIsEnforced(sourcedata->prizecyclecons[c]), SCIPconsIsChecked(sourcedata->prizecyclecons[c]),
                              SCIPconsIsPropagated(sourcedata->prizecyclecons[c]), SCIPconsIsLocal(sourcedata->prizecyclecons[c]),
                              SCIPconsIsModifiable(sourcedata->prizecyclecons[c]), SCIPconsIsDynamic(sourcedata->prizecyclecons[c]),
                              SCIPconsIsRemovable(sourcedata->prizecyclecons[c]), SCIPconsIsStickingAtNode(sourcedata->prizecyclecons[c]), global,
                              &success));
                  assert(success);

                  SCIP_CALL(SCIPcaptureCons(scip, (*targetdata)->prizecyclecons[c]));
               }
            }

            if( sourcedata->usesymcons && sourcedata->stp_type != STP_RPCSPG )
            {
#ifdef STP_AGG_SYM
               v = sourcedata->realnterms - 1;
#else
               v = ((sourcedata->realnterms - 1) * sourcedata->realnterms ) / 2;
#endif
               SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->prizesymcons, v));
               for( c = v - 1; c >= 0; --c )
               {
                  SCIP_CALL(
                        SCIPgetConsCopy(sourcescip, scip, sourcedata->prizesymcons[c], &((*targetdata)->prizesymcons[c]),
                              SCIPconsGetHdlr(sourcedata->prizesymcons[c]), varmap, consmap, SCIPconsGetName(sourcedata->prizesymcons[c]),
                              SCIPconsIsInitial(sourcedata->prizesymcons[c]), SCIPconsIsSeparated(sourcedata->prizesymcons[c]),
                              SCIPconsIsEnforced(sourcedata->prizesymcons[c]), SCIPconsIsChecked(sourcedata->prizesymcons[c]),
                              SCIPconsIsPropagated(sourcedata->prizesymcons[c]), SCIPconsIsLocal(sourcedata->prizesymcons[c]),
                              SCIPconsIsModifiable(sourcedata->prizesymcons[c]), SCIPconsIsDynamic(sourcedata->prizesymcons[c]),
                              SCIPconsIsRemovable(sourcedata->prizesymcons[c]), SCIPconsIsStickingAtNode(sourcedata->prizesymcons[c]), global,
                              &success));
                  assert(success);

                  SCIP_CALL(SCIPcaptureCons(scip, (*targetdata)->prizesymcons[c]));
               }
            }

            if( sourcedata->stp_type != STP_RPCSPG )
            {
               SCIP_CALL(
                     SCIPgetConsCopy(sourcescip, scip, sourcedata->prizecons, &((*targetdata)->prizecons), SCIPconsGetHdlr(sourcedata->prizecons),
                           varmap, consmap, SCIPconsGetName(sourcedata->prizecons), SCIPconsIsInitial(sourcedata->prizecons),
                           SCIPconsIsSeparated(sourcedata->prizecons), SCIPconsIsEnforced(sourcedata->prizecons),
                           SCIPconsIsChecked(sourcedata->prizecons), SCIPconsIsPropagated(sourcedata->prizecons),
                           SCIPconsIsLocal(sourcedata->prizecons), SCIPconsIsModifiable(sourcedata->prizecons),
                           SCIPconsIsDynamic(sourcedata->prizecons), SCIPconsIsRemovable(sourcedata->prizecons),
                           SCIPconsIsStickingAtNode(sourcedata->prizecons), global, &success));
               assert(success);

               SCIP_CALL(SCIPcaptureCons(scip, (*targetdata)->prizecons));
            }
         }

         if( sourcedata->stp_type == STP_DHCSTP )
         {
            SCIP_CALL(
                  SCIPgetConsCopy(sourcescip, scip, sourcedata->hopcons, &((*targetdata)->hopcons), SCIPconsGetHdlr(sourcedata->hopcons),
                        varmap, consmap, SCIPconsGetName(sourcedata->hopcons), SCIPconsIsInitial(sourcedata->hopcons),
                        SCIPconsIsSeparated(sourcedata->hopcons), SCIPconsIsEnforced(sourcedata->hopcons),
                        SCIPconsIsChecked(sourcedata->hopcons), SCIPconsIsPropagated(sourcedata->hopcons),
                        SCIPconsIsLocal(sourcedata->hopcons), SCIPconsIsModifiable(sourcedata->hopcons),
                        SCIPconsIsDynamic(sourcedata->hopcons), SCIPconsIsRemovable(sourcedata->hopcons),
                        SCIPconsIsStickingAtNode(sourcedata->hopcons), global, &success));
            assert(success);

            SCIP_CALL(SCIPcaptureCons(scip, (*targetdata)->hopcons));
         }

         if( sourcedata->stp_type == STP_BRMWCSP )
         {
            SCIP_CALL(
                  SCIPgetConsCopy(sourcescip, scip, sourcedata->budgetcons, &((*targetdata)->budgetcons), SCIPconsGetHdlr(sourcedata->budgetcons),
                        varmap, consmap, SCIPconsGetName(sourcedata->budgetcons), SCIPconsIsInitial(sourcedata->budgetcons),
                        SCIPconsIsSeparated(sourcedata->budgetcons), SCIPconsIsEnforced(sourcedata->budgetcons),
                        SCIPconsIsChecked(sourcedata->budgetcons), SCIPconsIsPropagated(sourcedata->budgetcons),
                        SCIPconsIsLocal(sourcedata->budgetcons), SCIPconsIsModifiable(sourcedata->budgetcons),
                        SCIPconsIsDynamic(sourcedata->budgetcons), SCIPconsIsRemovable(sourcedata->budgetcons),
                        SCIPconsIsStickingAtNode(sourcedata->budgetcons), global, &success));
            assert(success);

            SCIP_CALL(SCIPcaptureCons(scip, (*targetdata)->budgetcons));
         }

      }
      /* Price or Flow mode */
      else
      {
         /* transform edge constraints */
         if( sourcedata->bigt )
         {
            SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->edgecons, sourcedata->nedges));

            for( c = sourcedata->nedges - 1; c >= 0; --c )
            {
               SCIP_CALL(
                     SCIPgetConsCopy(sourcescip, scip, sourcedata->edgecons[c], &((*targetdata)->edgecons[c]),
                           SCIPconsGetHdlr(sourcedata->edgecons[c]), varmap, consmap, SCIPconsGetName(sourcedata->edgecons[c]),
                           SCIPconsIsInitial(sourcedata->edgecons[c]), SCIPconsIsSeparated(sourcedata->edgecons[c]),
                           SCIPconsIsEnforced(sourcedata->edgecons[c]), SCIPconsIsChecked(sourcedata->edgecons[c]),
                           SCIPconsIsPropagated(sourcedata->edgecons[c]), SCIPconsIsLocal(sourcedata->edgecons[c]),
                           SCIPconsIsModifiable(sourcedata->edgecons[c]), SCIPconsIsDynamic(sourcedata->edgecons[c]),
                           SCIPconsIsRemovable(sourcedata->edgecons[c]), SCIPconsIsStickingAtNode(sourcedata->edgecons[c]), global, &success));
               assert(success);

               SCIP_CALL(SCIPcaptureCons(scip, (*targetdata)->edgecons[c]));
            }
         }
         else
         {
            SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->edgecons, sourcedata->realnterms * sourcedata->nedges));
            for( c = sourcedata->realnterms * sourcedata->nedges - 1; c >= 0; --c )
            {
               SCIP_CALL(
                     SCIPgetConsCopy(sourcescip, scip, sourcedata->edgecons[c], &((*targetdata)->edgecons[c]),
                           SCIPconsGetHdlr(sourcedata->edgecons[c]), varmap, consmap, SCIPconsGetName(sourcedata->edgecons[c]),
                           SCIPconsIsInitial(sourcedata->edgecons[c]), SCIPconsIsSeparated(sourcedata->edgecons[c]),
                           SCIPconsIsEnforced(sourcedata->edgecons[c]), SCIPconsIsChecked(sourcedata->edgecons[c]),
                           SCIPconsIsPropagated(sourcedata->edgecons[c]), SCIPconsIsLocal(sourcedata->edgecons[c]),
                           SCIPconsIsModifiable(sourcedata->edgecons[c]), SCIPconsIsDynamic(sourcedata->edgecons[c]),
                           SCIPconsIsRemovable(sourcedata->edgecons[c]), SCIPconsIsStickingAtNode(sourcedata->edgecons[c]), global, &success));
               assert(success);

               SCIP_CALL(SCIPcaptureCons(scip, (*targetdata)->edgecons[c]));
            }
         }

         /* transform constraints */
         if( sourcedata->mode == STP_MODE_PRICE )
         {
            SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->pathcons, sourcedata->realnterms));
            for( c = sourcedata->realnterms - 1; c >= 0; --c )
            {
               SCIP_CALL(
                     SCIPgetConsCopy(sourcescip, scip, sourcedata->pathcons[c], &((*targetdata)->pathcons[c]),
                           SCIPconsGetHdlr(sourcedata->pathcons[c]), varmap, consmap, SCIPconsGetName(sourcedata->pathcons[c]),
                           SCIPconsIsInitial(sourcedata->pathcons[c]), SCIPconsIsSeparated(sourcedata->pathcons[c]),
                           SCIPconsIsEnforced(sourcedata->pathcons[c]), SCIPconsIsChecked(sourcedata->pathcons[c]),
                           SCIPconsIsPropagated(sourcedata->pathcons[c]), SCIPconsIsLocal(sourcedata->pathcons[c]),
                           SCIPconsIsModifiable(sourcedata->pathcons[c]), SCIPconsIsDynamic(sourcedata->pathcons[c]),
                           SCIPconsIsRemovable(sourcedata->pathcons[c]), SCIPconsIsStickingAtNode(sourcedata->pathcons[c]), global, &success));
               assert(success);

               SCIP_CALL(SCIPcaptureCons(scip, (*targetdata)->pathcons[c]));
            }
         }
         /* transform constraints and variables */
         else if( sourcedata->mode == STP_MODE_FLOW )
         {
            if( sourcedata->bigt )
            {
               SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->flowvars, sourcedata->nedges));
               SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->pathcons, (sourcedata->nnodes - 1)));
               for( c = sourcedata->nnodes - 2; c >= 0; --c )
               {
                  SCIP_CALL(
                        SCIPgetConsCopy(sourcescip, scip, sourcedata->pathcons[c], &((*targetdata)->pathcons[c]),
                              SCIPconsGetHdlr(sourcedata->pathcons[c]), varmap, consmap, SCIPconsGetName(sourcedata->pathcons[c]),
                              SCIPconsIsInitial(sourcedata->pathcons[c]), SCIPconsIsSeparated(sourcedata->pathcons[c]),
                              SCIPconsIsEnforced(sourcedata->pathcons[c]), SCIPconsIsChecked(sourcedata->pathcons[c]),
                              SCIPconsIsPropagated(sourcedata->pathcons[c]), SCIPconsIsLocal(sourcedata->pathcons[c]),
                              SCIPconsIsModifiable(sourcedata->pathcons[c]), SCIPconsIsDynamic(sourcedata->pathcons[c]),
                              SCIPconsIsRemovable(sourcedata->pathcons[c]), SCIPconsIsStickingAtNode(sourcedata->pathcons[c]), global, &success));
                  assert(success);

                  SCIP_CALL(SCIPcaptureCons(scip, (*targetdata)->pathcons[c]));
               }

               for( v = (*targetdata)->nedges - 1; v >= 0; --v )
               {
                  SCIP_CALL(
                        SCIPgetVarCopy(sourcescip, scip, sourcedata->flowvars[v], &((*targetdata)->flowvars[v]), varmap, consmap, global, &success));
                  assert(success);

                  SCIP_CALL(SCIPcaptureVar(scip, (*targetdata)->flowvars[v]));
               }
            }
            else
            {
               SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->flowvars, sourcedata->realnterms * sourcedata->nedges));
               SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->pathcons, sourcedata->realnterms * (sourcedata->nnodes - 1)));
               for( c = sourcedata->realnterms * (sourcedata->nnodes - 1) - 1; c >= 0; --c )
               {
                  SCIP_CALL(
                        SCIPgetConsCopy(sourcescip, scip, sourcedata->pathcons[c], &((*targetdata)->pathcons[c]),
                              SCIPconsGetHdlr(sourcedata->pathcons[c]), varmap, consmap, SCIPconsGetName(sourcedata->pathcons[c]),
                              SCIPconsIsInitial(sourcedata->pathcons[c]), SCIPconsIsSeparated(sourcedata->pathcons[c]),
                              SCIPconsIsEnforced(sourcedata->pathcons[c]), SCIPconsIsChecked(sourcedata->pathcons[c]),
                              SCIPconsIsPropagated(sourcedata->pathcons[c]), SCIPconsIsLocal(sourcedata->pathcons[c]),
                              SCIPconsIsModifiable(sourcedata->pathcons[c]), SCIPconsIsDynamic(sourcedata->pathcons[c]),
                              SCIPconsIsRemovable(sourcedata->pathcons[c]), SCIPconsIsStickingAtNode(sourcedata->pathcons[c]), global, &success));
                  assert(success);

                  SCIP_CALL(SCIPcaptureCons(scip, (*targetdata)->pathcons[c]));
               }

               for( v = sourcedata->nedges * sourcedata->realnterms - 1; v >= 0; --v )
               {
                  SCIP_CALL(
                        SCIPgetVarCopy(sourcescip, scip, sourcedata->flowvars[v], &((*targetdata)->flowvars[v]), varmap, consmap, global, &success));
                  assert(success);

                  SCIP_CALL(SCIPcaptureVar(scip, (*targetdata)->flowvars[v]));
               }
            }
         }
      }

      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->realterms, sourcedata->realnterms) );
      BMScopyMemoryArray((*targetdata)->realterms, sourcedata->realterms, sourcedata->realnterms);

      if( sourcedata->pctermsorder != NULL )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->pctermsorder, sourcedata->nnodes) );
         BMScopyMemoryArray((*targetdata)->pctermsorder, sourcedata->pctermsorder, sourcedata->nnodes);
      }
      else
         (*targetdata)->pctermsorder = NULL;
   }
   else
   {
      (*targetdata)->edgevars = NULL;
      (*targetdata)->xval = NULL;
      (*targetdata)->realterms = NULL;
      (*targetdata)->pctermsorder = NULL;
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
   SCIPdebugMessage("probdelorigStp \n");

   SCIPdebugMessage("free original problem data\n");

   if( !(*probdata)->graphHasVanished )
   {
      if ( (*probdata)->mode == STP_MODE_CUT )
         graph_mincut_exit(scip, (*probdata)->graph);

      graph_path_exit(scip, (*probdata)->graph);
   }

   graph_free(scip, &((*probdata)->graph), TRUE);

#ifdef WITH_UG
   if( (*probdata)->orggraph != NULL )
   {
      graph_path_exit(scip, (*probdata)->orggraph);
      graph_free(scip, &((*probdata)->orggraph), TRUE);
   }
#endif
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
   SCIP_Bool update;

   SCIPdebugMessage("probtransStp \n");

   SCIP_CALL( SCIPgetBoolParam(scip, "stp/countpresoltime", &update) );

   /* adjust time limit to take into account reading time */
   if( update )
   {
      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
      timelimit -= SCIPgetReadingTime(scip);
      timelimit = MAX(0.0,timelimit);
      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timelimit) );
   }

   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->graph) );
#ifdef WITH_UG
   (*targetdata)->orggraph = sourcedata->orggraph;
#endif
   (*targetdata)->nlayers = sourcedata->nlayers;
   (*targetdata)->nedges = sourcedata->nedges;
   (*targetdata)->nnodes = sourcedata->nnodes;
   (*targetdata)->nterms = sourcedata->nterms;
   (*targetdata)->realnterms = sourcedata->realnterms;
   (*targetdata)->nnonterms = sourcedata->nnonterms;
   (*targetdata)->emitgraph = sourcedata->emitgraph;
   (*targetdata)->usesymcons = sourcedata->usesymcons;
   (*targetdata)->graphHasVanished = sourcedata->graphHasVanished;
   (*targetdata)->usecyclecons = sourcedata->usecyclecons;
   (*targetdata)->nvars = sourcedata->nvars;
   (*targetdata)->mode = sourcedata->mode;
   (*targetdata)->offset = sourcedata->offset;
   (*targetdata)->norgedges = sourcedata->norgedges;
   (*targetdata)->minelims = sourcedata->minelims;
   (*targetdata)->bigt = sourcedata->bigt;
   (*targetdata)->logfile = sourcedata->logfile;
   (*targetdata)->origlogfile = &(sourcedata->logfile);
   (*targetdata)->intlogfile = sourcedata->intlogfile;
   (*targetdata)->origintlogfile = &(sourcedata->intlogfile);
#if USEOFFSETVAR
   if( sourcedata->offsetvar != NULL )
   {
      SCIP_CALL( SCIPtransformVar(scip, sourcedata->offsetvar, &(*targetdata)->offsetvar) );
   }
   else
      (*targetdata)->offsetvar = NULL;
#endif
   if( sourcedata->nedges > 0 )
   {
      int i;

      /* transform variables */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->xval, sourcedata->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->edgevars, sourcedata->nvars) );
      SCIP_CALL( SCIPtransformVars(scip, sourcedata->nvars, sourcedata->edgevars, (*targetdata)->edgevars) );

      /* cut mode */
      if( sourcedata->mode == STP_MODE_CUT )
      {
         (*targetdata)->edgecons = NULL;
         (*targetdata)->pathcons = NULL;

	      if( sourcedata->stp_type == STP_DCSTP )
	      {
	         SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->degcons, sourcedata->nnodes) );
            SCIP_CALL( SCIPtransformConss(scip, sourcedata->nnodes, sourcedata->degcons, (*targetdata)->degcons) );
	      }

#if FLOWB
            SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->flowbcons, sourcedata->nnonterms) );
            SCIP_CALL( SCIPtransformConss(scip, sourcedata->nnonterms, sourcedata->flowbcons, (*targetdata)->flowbcons) );
#endif

         if( sourcedata->stp_type == STP_PCSPG || sourcedata->stp_type == STP_MWCSP )
         {
            if( sourcedata->usecyclecons && sourcedata->stp_type == STP_PCSPG )
            {
               SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->prizecyclecons, sourcedata->nedges) );
               SCIP_CALL( SCIPtransformConss(scip, sourcedata->nedges, sourcedata->prizecyclecons, (*targetdata)->prizecyclecons) );
            }

            if( sourcedata->usesymcons )
            {
#ifdef STP_AGG_SYM
               i = sourcedata->realnterms - 1;
#else
               i = ((sourcedata->realnterms - 1) * (sourcedata->realnterms)) / 2;
#endif
               SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->prizesymcons, i) );
               SCIP_CALL( SCIPtransformConss(scip, i, sourcedata->prizesymcons, (*targetdata)->prizesymcons) );
            }
            SCIP_CALL( SCIPtransformCons(scip, sourcedata->prizecons, &(*targetdata)->prizecons) );
         }

         if( sourcedata->stp_type == STP_DHCSTP )
            SCIP_CALL( SCIPtransformCons(scip, sourcedata->hopcons, &(*targetdata)->hopcons) );

         if( sourcedata->stp_type == STP_BRMWCSP )
            SCIP_CALL( SCIPtransformCons(scip, sourcedata->budgetcons, &(*targetdata)->budgetcons) );
      }
      /* Price or Flow mode */
      else
      {
         /* transform edge constraints */
         if( sourcedata->bigt )
         {
            SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->edgecons, sourcedata->nedges));
            SCIP_CALL(SCIPtransformConss(scip, sourcedata->nedges, sourcedata->edgecons, (*targetdata)->edgecons));
         }
         else
         {
            SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->edgecons, sourcedata->realnterms * sourcedata->nedges));
            SCIP_CALL(SCIPtransformConss(scip, sourcedata->realnterms * sourcedata->nedges, sourcedata->edgecons, (*targetdata)->edgecons));
         }

         /* transform constraints */
         if( sourcedata->mode == STP_MODE_PRICE )
         {
            SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->pathcons, sourcedata->realnterms));
            SCIP_CALL(SCIPtransformConss(scip, sourcedata->realnterms, sourcedata->pathcons, (*targetdata)->pathcons));
         }
         /* transform constraints and variables*/
         else if( sourcedata->mode == STP_MODE_FLOW )
         {
            if( sourcedata->bigt )
            {
               SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->flowvars, sourcedata->nedges));
               SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->pathcons, (sourcedata->nnodes - 1)));
               SCIP_CALL(SCIPtransformConss(scip, (sourcedata->nnodes - 1), sourcedata->pathcons, (*targetdata)->pathcons));
               SCIP_CALL(SCIPtransformVars(scip, sourcedata->nedges, sourcedata->flowvars, (*targetdata)->flowvars));
            }
            else
            {
               SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->flowvars, sourcedata->realnterms * sourcedata->nedges));
               SCIP_CALL(SCIPallocMemoryArray(scip, &(*targetdata)->pathcons, sourcedata->realnterms * (sourcedata->nnodes - 1)));
               SCIP_CALL(SCIPtransformConss(scip, sourcedata->realnterms * (sourcedata->nnodes - 1), sourcedata->pathcons, (*targetdata)->pathcons));
               SCIP_CALL(SCIPtransformVars(scip, sourcedata->nedges * sourcedata->realnterms, sourcedata->flowvars, (*targetdata)->flowvars));
            }
         }
      }
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->realterms, sourcedata->realnterms) );
      BMScopyMemoryArray((*targetdata)->realterms, sourcedata->realterms, sourcedata->realnterms);

      if( sourcedata->pctermsorder != NULL )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->pctermsorder, sourcedata->nnodes) );
         BMScopyMemoryArray((*targetdata)->pctermsorder, sourcedata->pctermsorder, sourcedata->nnodes);
      }
      else
         (*targetdata)->pctermsorder = NULL;
   }
   else
   {
      (*targetdata)->edgevars = NULL;
      (*targetdata)->xval = NULL;
      (*targetdata)->realterms = NULL;
      (*targetdata)->pctermsorder = NULL;
      (*targetdata)->edgecons = NULL;
      (*targetdata)->pathcons = NULL;
      (*targetdata)->flowvars = NULL;
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_PROBINITSOL(probinitsolStp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(probdata != NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_PROBEXITSOL(probexitsolStp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(probdata != NULL);

   if( probdata->logfile != NULL )
   {
      int success;
      SCIP_Real factor = 1.0;

#ifdef WITH_UG
      if( getUgRank() != 0 )
         return SCIP_ERROR;
#endif

      if( probdata->stp_type == STP_MWCSP || probdata->stp_type == STP_RMWCSP || probdata->stp_type == STP_BRMWCSP )
         factor = -1.0;

      SCIPprobdataWriteLogLine(scip, "End\n");
      SCIPprobdataWriteLogLine(scip, "\n");
      SCIPprobdataWriteLogLine(scip, "SECTION Run\n");
      if( probdata->ug )
      {
         SCIPprobdataWriteLogLine(scip, "Threads %d\n", probdata->nSolvers);
         SCIPprobdataWriteLogLine(scip, "Time %.1f\n", SCIPgetTotalTime(scip));
      }
      else
      {
         SCIPprobdataWriteLogLine(scip, "Threads 1\n");
         SCIPprobdataWriteLogLine(scip, "Time %.1f\n", SCIPgetTotalTime(scip));
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
      *(probdata->origlogfile) = NULL;
   }

   if( probdata->intlogfile != NULL )
   {
      const int success = fclose(probdata->intlogfile);

#ifdef WITH_UG
      if( getUgRank() != 0 )
         return SCIP_ERROR;
#endif

      if( success != 0 )
      {
         SCIPerrorMessage("An error occurred while closing file <%s>\n", probdata->intlogfile);
         return SCIP_FILECREATEERROR;
      }

      probdata->intlogfile = NULL;
      *(probdata->origintlogfile) = NULL;
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
SCIP_RETCODE SCIPprobdataSetDefaultParams(
   SCIP*                 scip                /**< SCIP data structure */
)
{

   SCIP_CALL( SCIPsetSubscipsOff(scip, TRUE) );

   /* set STP-specific default parameters */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/freq", 1) );
   SCIP_CALL( SCIPsetIntParam(scip, "limits/maxsol", 400) );
   SCIP_CALL( SCIPsetIntParam(scip, "lp/rowagelimit", 30) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxroundsroot", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxrounds", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxstallroundsroot", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxcutsroot", 100000) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxcuts", 1000) );   // todo tune

   SCIP_CALL( SCIPsetRealParam(scip, "separating/minefficacyroot", 0.01) ); // todo tune
   SCIP_CALL( SCIPsetRealParam(scip, "separating/minorthoroot", 0.4) ); // todo tune > 0.4
   SCIP_CALL( SCIPsetRealParam(scip, "separating/minortho", 0.4) ); // todo tune > 0.4 best soplex: 0.8
   SCIP_CALL( SCIPsetRealParam(scip, "separating/objparalfac", 0.01) ); // todo tune < 0.1

   SCIP_CALL( SCIPsetRealParam(scip, "separating/intsupportfac", 0.0) );
   SCIP_CALL( SCIPsetIntParam(scip, "branching/relpscost/maxproprounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/alns/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/coefdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/feaspump/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/fracdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/farkasdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/guideddiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/linesearchdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/nlpdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/objpscostdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/pscostdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/randrounding/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rootsoldiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/shiftandpropagate/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/shifting/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/subnlp/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/undercover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/veclendiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/zirounding/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/oneopt/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/locks/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/probing/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/pseudoobj/timingmask", 5) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/redcost/freq", -1) );
   SCIP_CALL( SCIPsetRealParam(scip, "branching/relpscost/maxreliable", 1.0) );

   SCIP_CALL( SCIPsetBoolParam(scip, "timing/reading", TRUE) );

   // todo test properly; normal dfs?
   SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/restartdfs/stdpriority", 400000) );

   return SCIP_OKAY;
}

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_Real offset = 0.0;
   PRESOL presolinfo;
   GRAPH* graph;
   SCIP_Bool printGraph;
   int symcons;
   int cyclecons;
   int usedacuts;
   int compcentral;
   char* intlogfilename;
   char* logfilename;
   char* probname;
   char tmpfilename[SCIP_MAXSTRLEN];

   assert(scip != NULL);

   presolinfo.fixed = 0;

   /* read graph from file */
   SCIP_CALL( graph_load(scip, &graph, filename, &presolinfo) );

   SCIPdebugMessage("load type :: %d \n\n", graph->stp_type);
   SCIPdebugMessage("fixed: %f \n\n", presolinfo.fixed );

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, graph) );

   /* get parameters */
   SCIP_CALL( SCIPgetIntParam(scip, "stp/compcentral", &compcentral) );
   SCIP_CALL( SCIPgetIntParam(scip, "stp/usesymcons", &(symcons)) );
   SCIP_CALL( SCIPgetIntParam(scip, "stp/usecyclecons", &(cyclecons)) );
   SCIP_CALL( SCIPgetIntParam(scip, "stp/usedacuts", &(usedacuts)) );
   SCIP_CALL( SCIPgetIntParam(scip, "stp/minelims", &(probdata->minelims)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "stp/emitgraph", &(probdata->emitgraph)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "stp/bigt", &(probdata->bigt)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "stp/printGraph", &printGraph) );
   SCIP_CALL( SCIPgetStringParam(scip, "stp/logfile", &logfilename) );
   SCIP_CALL( SCIPgetStringParam(scip, "stp/intlogfile", &intlogfilename) );

   /* copy filename */
   (void) SCIPsnprintf(tmpfilename, SCIP_MAXSTRLEN, "%s", filename);
   SCIPsplitFilename(tmpfilename, NULL, &probname, NULL, NULL);

   SCIP_CALL( createLogfile(scip, probdata, intlogfilename, logfilename, filename, probname) );

   /* create a problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );
   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigStp) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransStp) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransStp) );
   SCIP_CALL( SCIPsetProbInitsol(scip, probinitsolStp) );
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolStp) );
   SCIP_CALL( SCIPsetProbCopy(scip, probcopyStp) );

   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   setStpSolvingMode(scip, probdata);

   if( printGraph )
   {
      SCIP_CALL( graph_writeGml(graph, "OriginalGraph.gml", NULL) );
   }

   probdata->graph = graph;

   /* reduce the graph (and do some house-holding) */
   presolveStp(scip, probdata, &offset);

   graph = probdata->graph;

   writeCommentSection(scip, graph, filename);

   if( printGraph )
   {
      SCIP_CALL( graph_writeGml(graph, "ReducedGraph.gml", NULL) );
   }

   /* select a root node */
   if( graph->stp_type == STP_SPG || graph->stp_type == STP_NWPTSPG || graph->stp_type == STP_NWSPG )
   {
      SCIP_CALL( findCentralTerminal(scip, graph, &(graph->source), compcentral) );
   }

#ifdef WITH_UG
   {
      GRAPH* newgraph;
      graph_copy(scip, graph, &(newgraph));
      probdata->orggraph = graph;
      graph = newgraph;
   }
#endif

   SCIP_CALL( setParams(scip, probdata, symcons, cyclecons) );

   /* setting the offset to the fixed value given in the input file plus the fixings
    * given by the reduction techniques */
   probdata->offset = presolinfo.fixed + offset;

   SCIP_CALL( createModel(scip, probdata) );

#ifdef PRINT_PRESOL
   (void)SCIPsnprintf(presolvefilename, SCIP_MAXSTRLEN, "presol/%s-presolve.stp", probname);
   SCIP_CALL( SCIPwriteOrigProblem(scip, presolvefilename, NULL, FALSE) );
#endif

   /* todo test */
   assert(usedacuts == STP_CONS_ALWAYS || usedacuts == STP_CONS_AUTOMATIC);

   if( probdata->mode == STP_MODE_CUT && usedacuts != STP_CONS_NEVER )
   {
      assert(usedacuts == STP_CONS_ALWAYS || usedacuts == STP_CONS_AUTOMATIC);

      SCIP_CALL( createInitialCuts(scip, probdata) );
   }

   return SCIP_OKAY;
}

/** sets the probdata graph */
void SCIPprobdataSetGraph(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   GRAPH*                graph               /**< graph data structure */
   )
{
   assert(probdata != NULL);

   probdata->graph = graph;
}


/** returns the graph */
GRAPH* SCIPprobdataGetGraph(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->graph;
}

/** returns the graph */
GRAPH* SCIPprobdataGetGraph2(
   SCIP*                 scip            /**< problem data */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->graph;
}

/** sets the offset */
void SCIPprobdataSetOffset(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Real             offset              /**< the offset value */
   )
{
   assert(probdata != NULL);

   probdata->offset = offset;
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

/** returns the number of nodes */
int SCIPprobdataGetNNodes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->nnodes;
}

/** returns the number of terminals */
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

/** returns the number of terminals without the root node */
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

/** returns root */
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

   return graph->source;
}

/** returns numer of original edges */
int SCIPprobdataGetNorgEdges(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->norgedges;
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
   int                   idx                 /**< index of the edge */
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
   SCIP_SOL*             sol                 /**< solution to get values from */
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

   nedges = probdata->nedges;
   assert(nedges >= 0);

   SCIP_CALL_ABORT( SCIPgetSolVals(scip, sol, nedges, probdata->edgevars, vals) );

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

/** returns array */
int* SCIPprobdataGetPctermsorder(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->pctermsorder;
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
      SCIP_CALL( graph_writeGml( probdata->graph, filename, NULL) );
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
      SCIP_CALL( graph_writeGml( probdata->graph, filename, edgemark) );

      SCIPfreeBufferArray(scip, &edgemark);
   }

   return SCIP_OKAY;
}


/** writes the best solution to the intermediate solution file */
SCIP_RETCODE SCIPprobdataWriteIntermediateSolution(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);

   if( probdata->intlogfile != NULL )
   {
      SCIP_CALL( SCIPprobdataWriteSolution(scip, probdata->intlogfile) );
   }

   return SCIP_OKAY;
}

/** writes the best solution to a file */
SCIP_RETCODE SCIPprobdataWriteSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< file to write best solution to; or NULL, to write to stdout */
   )
{
   SCIP_SOL* sol;
   SCIP_VAR** edgevars;
   GRAPH* graph;
   IDX** ancestors;
   IDX* curr;
   SCIP_PROBDATA* probdata;
   int  e;
   int  k;
   int  norgedges;
   int  norgnodes;
   int  nsolnodes;
   int  nsoledges;
   STP_Bool* orgedges;
   STP_Bool* orgnodes;

   assert(scip != NULL);

#ifdef WITH_UG
   if( getUgRank() != 0 )
      return SCIP_OKAY;
#endif

   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);

   edgevars = probdata->edgevars;
#ifdef WITH_UG
   graph = probdata->orggraph;
#else
   graph = probdata->graph;
#endif

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
   for( k = 0; k < norgnodes; k++ )
      orgnodes[k] = FALSE;

   ancestors = graph->ancestors;

   if( graph->stp_type == STP_SPG || graph->stp_type == STP_SAP ||graph->stp_type == STP_DCSTP
      || graph->stp_type == STP_NWSPG || graph->stp_type == STP_DHCSTP || graph->stp_type == STP_GSTP
      || graph->stp_type == STP_RSMT || graph->stp_type == STP_NWPTSPG  )
   {
      GRAPH* solgraph;
      SCIP_QUEUE* queue;
      int* nodechild;
      int* edgeancestor;
      int* pnode;
      int tail;
      int head;

      /* iterate through the list of fixed edges */
      updateorgsol(graph, graph_get_fixedges(graph), orgnodes, orgedges, &nsolnodes, &nsoledges);

      for( e = 0; e < graph->edges; e++ )
         if( !SCIPisZero(scip, SCIPgetSolVal(scip, sol, edgevars[e])) )
            /* iterate through the list of ancestors */
            updateorgsol(graph, ancestors[e], orgnodes, orgedges, &nsolnodes, &nsoledges);

#if 1
      SCIP_CALL( SCIPallocBufferArray(scip, &nodechild, norgnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &edgeancestor, 2 * nsoledges) );

      /* initialize new graph */
      SCIP_CALL( graph_init(scip, &solgraph, nsolnodes, 2 * nsoledges, 1) );

      /* add vertices to new graph */
      for( k = 0; k < norgnodes; k++ )
      {
         if( orgnodes[k] )
         {
            nodechild[k] = solgraph->knots;
            graph_knot_add(solgraph, -1);
         }
         else
         {
            nodechild[k] = -1;
         }
      }

      assert(nsolnodes == solgraph->knots);
      assert(nodechild[graph->orgsource] >= 0);

      /* set root of new graph */
      solgraph->source = nodechild[graph->orgsource];

      /* add edges to new graph */
      for( e = 0; e < norgedges; e++ )
      {
         if( !orgedges[e] )
            continue;

         tail = nodechild[graph->orgtail[e]];
         head = nodechild[graph->orghead[e]];

         assert(tail >= 0);
         assert(head >= 0);

         edgeancestor[solgraph->edges] = e;
         edgeancestor[solgraph->edges + 1] = flipedge(e);
         graph_edge_add(scip, solgraph, tail, head, 1.0, 1.0);
      }

      for( e = 0; e < norgedges; e++ )
         orgedges[e] = FALSE;

      for( k = 0; k < norgnodes; k++ )
         orgnodes[k] = FALSE;

      SCIP_CALL( SCIPqueueCreate(&queue, nsolnodes, 1.1) );
      SCIP_CALL( SCIPqueueInsert(queue, &(solgraph->source)) );
      solgraph->mark[solgraph->source] = FALSE;

      nsolnodes = 1;

      /* BFS from root */
      while( !SCIPqueueIsEmpty(queue) )
      {
         pnode = (SCIPqueueRemove(queue));
         k = *pnode;

         /* traverse outgoing arcs */
         for( e = solgraph->outbeg[k]; e != EAT_LAST; e = solgraph->oeat[e] )
         {
            head = solgraph->head[e];

            /* vertex not labeled yet? */
            if( solgraph->mark[head] )
            {
               orgedges[edgeancestor[e]] = TRUE;
               solgraph->mark[head] = FALSE;
               nsolnodes++;
               SCIP_CALL(SCIPqueueInsert(queue, &(solgraph->head[e])));
            }
         }
      }

      nsoledges = nsolnodes - 1;

      SCIPqueueFree(&queue);

      graph_free(scip, &solgraph, TRUE);

      for( e = 0; e < norgedges; e++ )
      {
         if( orgedges[e] )
         {
            orgnodes[graph->orgtail[e]] = TRUE;
            orgnodes[graph->orghead[e]] = TRUE;
         }
      }

      SCIPfreeBufferArray(scip, &edgeancestor);
      SCIPfreeBufferArray(scip, &nodechild);
#endif
      if( graph->stp_type == STP_GSTP )
      {
         norgnodes -= graph->terms;
         nsolnodes -= graph->terms;
         nsoledges -= graph->terms;
         assert(nsolnodes >= 0);
         assert(nsoledges >= 1);
      }

      SCIPprobdataWriteLogLine(scip, "Vertices %d\n", nsolnodes);

      for( k = 0; k < norgnodes; k++ )
         if( orgnodes[k] == TRUE )
            SCIPinfoMessage(scip, file, "V %d\n", k + 1);

      SCIPprobdataWriteLogLine(scip, "Edges %d\n", nsoledges);

      if( graph->stp_type == STP_DHCSTP )
      {
         for( e = 0; e < norgedges; e++ )
            if( orgedges[e] == TRUE )
               SCIPinfoMessage(scip, file, "E %d %d\n", graph->orgtail[e] + 1, graph->orghead[e] + 1);
      }
      else
      {
         for( e = 0; e < norgedges; e += 2 )
         {
            if( graph->stp_type == STP_GSTP && (Is_term(graph->term[graph->orgtail[e]]) || Is_term(graph->term[graph->orghead[e]])) )
               continue;
            if( orgedges[e] == TRUE || orgedges[e + 1] == TRUE )
               SCIPinfoMessage(scip, file, "E %d %d\n", graph->orgtail[e] + 1, graph->orghead[e] + 1);
         }
      }
   }
   else if( graph->stp_type == STP_RSMT )
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
      for( k = 0; k < norgnodes; k++ )
      {
         if( orgnodes[k] == TRUE )
         {
            nodenumber[k] = nodecount++;
            SCIPprobdataWriteLogLine(scip, "%s ", strdim);
            SCIP_CALL( graph_grid_coordinates(scip, coords, &nodecoords, ncoords, k, grid_dim) );
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
   else if( graph_pc_isPcMw(graph))
   {
      int* solnodequeue;
      assert(graph->source >= 0 && nsolnodes == 0);

      SCIP_CALL( SCIPallocBufferArray(scip, &solnodequeue, norgnodes) );

      /* cover RPCSPG/RMWCSP with single node solution */
      if( graph_pc_isRootedPcMw(graph) && graph->orgsource != -1 )
      {
         SCIPdebugMessage("graph->orgsource %d \n", graph->orgsource);
         solnodequeue[nsolnodes++] = graph->orgsource;
         orgnodes[graph->orgsource] = TRUE;
      }

      for( e = 0; e <= graph->edges; e++ )
      {
         if( e == graph->edges || !SCIPisZero(scip, SCIPgetSolVal(scip, sol, edgevars[e])) )
         {
            /* iterate through the list of ancestors/fixed edges */
            if( e < graph->edges )
               curr = ancestors[e];
            else
               curr = graph_get_fixedges(graph);

            while (curr != NULL)
            {
               const int ancestoredge = curr->index;
               if( e < graph->edges && (graph->stp_type == STP_MWCSP || graph->stp_type == STP_RMWCSP || graph->stp_type == STP_BRMWCSP) )
               {
                  if( !SCIPisZero(scip, SCIPgetSolVal(scip, sol, edgevars[flipedge(e)])) )
                  {
                     curr = curr->parent;
                     continue;
                  }
               }

               if( ancestoredge < graph->norgmodeledges )
               {
                  if( !orgedges[ancestoredge] )
                  {
                     orgedges[ancestoredge] = TRUE;
                     nsoledges++;
                  }
                  if( !orgnodes[graph->orgtail[ancestoredge]] )
                  {
                     orgnodes[graph->orgtail[ancestoredge]] = TRUE;
                     solnodequeue[nsolnodes++] = graph->orgtail[ancestoredge];
                  }
                  if( !orgnodes[graph->orghead[ancestoredge]] )
                  {
                     orgnodes[graph->orghead[ancestoredge]] = TRUE;
                     solnodequeue[nsolnodes++] = graph->orghead[ancestoredge];
                  }
               }
               else if( graph->orghead[ancestoredge] < graph->norgmodelknots )
               {
                  if( !orgnodes[graph->orghead[ancestoredge]])
                  {
                     orgnodes[graph->orghead[ancestoredge]] = TRUE;
                     solnodequeue[nsolnodes++] = graph->orghead[ancestoredge];
                  }
               }
               curr = curr->parent;
            }
         }
      }

      SCIP_CALL( graph_sol_markPcancestors(scip, graph->pcancestors, graph->orgtail, graph->orghead, norgnodes,
            orgnodes, orgedges, solnodequeue, &nsolnodes, &nsoledges ) );

      SCIPprobdataWriteLogLine(scip, "Vertices %d\n", nsolnodes);

      for( k = 0; k < norgnodes; k++ )
         if( orgnodes[k] == TRUE )

	   SCIPinfoMessage(scip, file, "V %d\n", k + 1);

      SCIPprobdataWriteLogLine(scip, "Edges %d\n", nsoledges);

      for( e = 0; e < norgedges; e += 2 )
         if( orgedges[e] == TRUE || orgedges[e + 1] == TRUE )
            SCIPinfoMessage(scip, file, "E %d %d\n", graph->orgtail[e] + 1, graph->orghead[e] + 1);

      SCIPfreeBufferArray(scip, &solnodequeue);
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

#ifdef WITH_UG
   if( getUgRank() != 0 )
      return;
#endif

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   if( probdata->logfile != NULL )
   {
      va_start(ap, formatstr); /*lint !e826*/
      SCIPmessageVFPrintInfo(SCIPgetMessagehdlr(scip), probdata->logfile, formatstr, ap);
      va_end(ap);
   }
   if( probdata->intlogfile != NULL )
   {
      va_start(ap, formatstr); /*lint !e826*/
      SCIPmessageVFPrintInfo(SCIPgetMessagehdlr(scip), probdata->intlogfile, formatstr, ap);
      va_end(ap);
   }
}

/** add new solution */
SCIP_RETCODE SCIPprobdataAddNewSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            nval,               /**< array [0..nvars], nval[v] = 1 if node v is in the solution, nval[v] = 0 if not */
   SCIP_HEUR*            heur,               /**< heuristic data */
   SCIP_Bool*            success             /**< denotes whether the new solution has been successfully added */
   )
{
   SCIP_SOL* sol = NULL;
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
   if( probdata->mode != STP_MODE_CUT )
   {
      SCIP_Real* edgecost;
      SCIP_Real* flowvals;
      PATH* path;
      SCIP_VAR** pathvars;
      SCIP_VAR* var = NULL;
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
      if ( probdata->mode == STP_MODE_FLOW )
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
         graph->mark[e] = 1;
      graph_path_exec(scip, graph, FSP_MODE, graph->source, edgecost, path);

      /* create and add path variables (Price mode) or set the flow variables (Flow mode) */
      for( t = 0; t < realnterms; ++t )
      {
         if( probdata->mode == STP_MODE_PRICE )
         {
	    /* create a new path variable */
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "PathVar%d_X", t);
            var = NULL;
            SCIP_CALL( SCIPcreateVarBasic(scip, &var, varname, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, var) );
            SCIP_CALL( SCIPaddCoefLinear(scip, probdata->pathcons[t], var, 1.0) );
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, 1.0) );
            assert(var != NULL);
            assert(pathvars != NULL);
            pathvars[t] = var;
         }
         tail = probdata->realterms[t];

         /* walk from terminal t to the root */
         while( tail != graph->source )
         {
            if( !probdata->bigt )
            {
               if( probdata->mode == STP_MODE_PRICE )
               {
                  /* add the new path variable to the constraints corresponding to the current edge */
                  assert(var != NULL);
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
               if( probdata->mode == STP_MODE_PRICE )
               {
                  assert(var != NULL);
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
         if( probdata->mode == STP_MODE_PRICE )
         {
            assert(var != NULL);
            SCIP_CALL( SCIPreleaseVar(scip, &var) );
         }
      }

      /* store the new solution value */
      SCIP_CALL( SCIPsetSolVals(scip, sol, probdata->nvars, edgevars, nval) );
      if( probdata->mode == STP_MODE_FLOW )
      {
         SCIP_CALL( SCIPsetSolVals(scip, sol, nedges * (probdata->bigt ? 1 : realnterms) , probdata->flowvars, flowvals) );
      }

      /* try to add new solution to scip and free it immediately */
      SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, FALSE, TRUE, TRUE, TRUE, success) );

      /* free local arrays */
      SCIPfreeMemoryArrayNull(scip, &flowvals);
      SCIPfreeMemoryArrayNull(scip, &edgecost);
      SCIPfreeMemoryArrayNull(scip, &path);
      SCIPfreeMemoryArrayNull(scip, &pathvars);
   }
   /* cut mode */
   else
   {
      SCIP_Bool feasible;
      int e;
      int nvars = probdata->nvars;

      /* check whether the new solution is valid with respect to the original bounds */
      for( e = 0; e < nvars; e++ )
      {
         if( SCIPisGT(scip, nval[e], SCIPvarGetUbGlobal(edgevars[e])) ||  SCIPisGT(scip, SCIPvarGetLbGlobal(edgevars[e]), nval[e]) )
         {
	         SCIPdebugMessage("solution not valid wrt to global bounds: %d %d nval %f ub: %f \n",
	               graph->tail[e],  graph->head[e], nval[e], SCIPvarGetUbGlobal(edgevars[e]) );
            *success = FALSE;
            SCIP_CALL( SCIPfreeSol(scip, &sol) );
            return SCIP_OKAY;
         }
      }

      /* post-processing of solution for MWCS and PCSPG */
      if( graph->stp_type == STP_PCSPG || graph->stp_type == STP_MWCSP )
      {
         int k;

         for( k = 0; k < graph->knots; ++k )
         {
            /* is the kth node a terminal other than the root? */
            if( Is_term(graph->term[k]) && k != graph->source )
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
                     nval[edge1] = 0;
                     nval[edge2] = 1;
                     break;
                  }
               }
            }
         }
      }

      /* store the new solution value */
      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, edgevars, nval) );
#if USEOFFSETVAR
      if( probdata->offsetvar != NULL && SCIPvarIsActive(probdata->offsetvar) )
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->offsetvar, 1.0) );
      }
#endif
      SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, FALSE, TRUE, TRUE, TRUE, &feasible) );

      SCIPdebugMessage("checked sol: feasible=%u\n", feasible);

      /* check solution for feasibility in original problem space */
      if( !feasible )
      {
         SCIP_CALL( SCIPcheckSolOrig(scip, sol, &feasible, TRUE, TRUE) );
         SCIPdebugMessage("checked sol org: feasible=%u\n", feasible);

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
         SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, FALSE, TRUE, TRUE, TRUE, &feasible) );
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

      if( probdata->stp_type == STP_MWCSP || probdata->stp_type == STP_RMWCSP || probdata->stp_type == STP_BRMWCSP )
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

/** branching information from UG */
void initReceivedSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const int             lLinearConsNames,   /**< number of linear constraints */
   const char*           linearConsNames,    /**< linear constraints string */
   const int             lSetppcConsNames,   /**< number of setppc constraints */
   const char*           setppcConsNames     /**< number of setppc constraints */
   )
{

#ifdef WITH_UG
   SCIP_PROBDATA* probdata;
   GRAPH* graph;
   int* nodestate;
   int nnodes;

   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPprobdataSetDefaultParams(scip) );

   probdata = SCIPgetProbData(scip);

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);
   assert(probdata->orggraph != NULL);

   graph_mincut_exit(scip, graph);
   graph_path_exit(scip, graph);
   assert(graph->ancestors == NULL);

   graph_free(scip, &graph, TRUE);
   graph_copy(scip, probdata->orggraph, &graph);

   probdata->graph = graph;

   assert(graph != NULL && probdata->mode == STP_MODE_CUT);

   nnodes = graph->knots;
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &nodestate, nnodes) );
   SCIPStpBranchruleInitNodeState(graph, nodestate);

   for( int i = 0; i < lLinearConsNames; i++ )
   {
      const char* consname = getBranchLinearConsName(linearConsNames, i);
      SCIPdebugMessage("add lin cons %s \n", consname);
      if( consname != NULL)
         SCIP_CALL_ABORT( STPStpBranchruleParseConsname(nodestate, consname) );
   }

   for( int i = 0; i < lSetppcConsNames; i++ )
   {
      const char* consname = getBranchSetppcConsName(setppcConsNames, i);
      SCIPdebugMessage("add ppc cons %s \n", consname);
      if( consname != NULL)
         SCIP_CALL_ABORT( STPStpBranchruleParseConsname(nodestate, consname) );
   }

   for( int k = 0; k < nnodes; k++ )
   {
      if( nodestate[k] == BRANCH_STP_VERTEX_TERM && !Is_term(graph->term[k]) )
      {
         SCIPdebugMessage("UG make term: %d \n", k);

         if( graph_pc_isPcMw(graph) )
         {
            if( Is_pseudoTerm(graph->term[k]) )
            {
               graph_pc_enforcePterm(graph, k);
            }
            else if( graph_pc_isRootedPcMw(graph) )
            {
               graph->prize[k] = FARAWAY;
               graph_knot_chg(graph, k, 0);
            }
         }
         else
         {
            graph_knot_chg(graph, k, 0);
         }
      }
      else if( nodestate[k] == BRANCH_STP_VERTEX_KILLED )
      {
         assert(!Is_term(graph->term[k]));

         for( int e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
         {
            if( graph->cost[e] < BLOCKED )
               graph->cost[e] = BLOCKED;

            if( graph->cost[flipedge(e)] < BLOCKED )
               graph->cost[flipedge(e)] = BLOCKED;
         }

         if( Is_pseudoTerm(graph->term[k]) )
         {
            const int twinterm = graph_pc_getTwinTerm(graph, k);
            const int root2twin = graph_pc_getRoot2PtermEdge(graph, twinterm);
            const int k2twin = graph->term2edge[k];

            assert(graph_pc_isPcMw(graph));
            assert(k2twin >= 0 && graph->tail[k2twin] == k && graph->head[k2twin] == twinterm);
            assert(root2twin >= 0 && graph->head[root2twin] == twinterm);
            assert(SCIPisEQ(scip, graph->cost[root2twin], graph->prize[k]));

            graph->cost[k2twin] = 0.0;

            if( !graph_pc_isRootedPcMw(graph) )
            {
               const int root2k = graph_pc_getRoot2PtermEdge(graph, twinterm);
               assert(graph->tail[root2k] == graph->source && graph->head[root2k] == k);

               graph->cost[root2k] = 0.0;
            }
         }
      }
   } /* for k = 0 : nnnodes */

   SCIPfreeBufferArray(scip, &nodestate);

   SCIP_CALL_ABORT( graph_path_init(scip, graph) );
   SCIP_CALL_ABORT( graph_mincut_init(scip, graph) );

#else
   assert(0 && "only call me when using UG");
#endif
}
