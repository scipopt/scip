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

/**@file   cons_stp.c
 * @brief  Constraint handler for Steiner problems
 * @author Gerald Gamrath
 * @author Daniel Rehfeldt
 * @author Michael Winkler
 *
 * This file checks solutions for feasibility and separates violated model constraints. For more details see \ref STP_CONS page.
 *
 * @page STP_CONS Separating violated constraints
 *
 * In this file a constraint handler checking solutions for feasibility and separating violated model constraints is implemented.
 * The separation problem for the cut inequalities described in \ref STP_PROBLEM can be solved by a max-flow algorithm in
 * polynomial time.  Regarding the variable values of a given LP solution as capacities on the edges, one can check for each
 * \f$ t \in T \setminus \{r\} \f$, with \f$ r \f$ being the root, whether the minimal \f$ (r, t) \f$-cut is less than one. In this case,
 * a violated cut inequality has been found, otherwise none exists. In order to calculate such a minimal cut an adaptation of Hao
 * and Orlin's preflow-push algorithm (see A Faster Algorithm for Finding the Minimum Cut in a Directed Graph) is used. Furthermore, the file implements a dual ascent heuristic, based on a concept described
 * in "A dual ascent approach for Steiner tree problems on a directed graph." by R. Wong.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "cons_stp.h"
#include "probdata_stp.h"
#include "graph.h"
#include "portab.h"
#include "branch_stp.h"
#include "prop_stp.h"

#include "scip/scip.h"
#include "scip/misc.h"
#include "scip/cons_linear.h"
#include <time.h>
#if 0
#ifdef WITH_UG
#define ADDCUTSTOPOOL 1
#else
#define ADDCUTSTOPOOL 0
#endif
#endif

#define ADDCUTSTOPOOL 0

#define Q_NULL     -1         /* NULL element of queue/list */

/**@name Constraint handler properties
 *
 * @{
 */

#define CONSHDLR_NAME          "stp"
#define CONSHDLR_DESC          "steiner tree constraint handler"
#define CONSHDLR_SEPAPRIORITY   9999999 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  9999999 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             0 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ            1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define DEFAULT_MAXROUNDS            20 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS     INT_MAX /**< maximal number of cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT INT_MAX /**< maximal number of cuts separated per separation round in the root node */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP

#define PCIMPLICATIONS_ALLOC_FACTOR 4
#define DEFAULT_BACKCUT        FALSE /**< Try Back-Cuts FALSE*/
#define DEFAULT_CREEPFLOW      TRUE  /**< Use Creep-Flow */
#define DEFAULT_DISJUNCTCUT    FALSE /**< Only disjunct Cuts FALSE */
#define DEFAULT_NESTEDCUT      FALSE /**< Try Nested-Cuts FALSE*/
#define DEFAULT_FLOWSEP        TRUE  /**< Try Flow-Cuts */
#define DEFAULT_INFLOWSEP       TRUE  /**< Try in-flow Cuts */
#define DEFAULT_INFLOWTERMSEP   TRUE  /**< Try terminal in-flow Cuts */
#define DEFAULT_OUTFLOWSEP      TRUE
#define DEFAULT_BALANCEFLOWSEP  TRUE


/* *
#define FLOW_FACTOR     100000
#define CREEP_VALUE     1         this is the original value todo check what is better
*/

#define FLOW_FACTOR     1000000
#define CREEP_VALUE     10



/**@} */

/*
 * Data structures
 */

/** @brief Constraint data for  \ref cons_stp.c "Stp" constraints */
struct SCIP_ConsData
{
   GRAPH*                graph;              /**< graph data structure */
};

/** @brief Constraint handler data for \ref cons_stp.c "Stp" constraint handler */
struct SCIP_ConshdlrData
{
   int*                  pcimplstart;        /**< start for each proper potential terminal */
   int*                  pcimplverts;        /**< all vertices */
   int                   pcimplnppterms;     /**< number of poper potential terminals used */
   SCIP_Bool             backcut;            /**< should backcuts be applied? */
   SCIP_Bool             creepflow;          /**< should creepflow cuts be applied? */
   SCIP_Bool             disjunctcut;        /**< should disjunction cuts be applied? */
   SCIP_Bool             nestedcut;          /**< should nested cuts be applied? */
   SCIP_Bool             flowsep;            /**< should flow separation be applied? */
   SCIP_Bool             inflowsep;          /**< should unit in-flow separation be applied? */
   SCIP_Bool             intermflowsep;      /**< should unit terminal in-flow separation be applied? */
   SCIP_Bool             outflowsep;         /**< should single-edge out-flow separation be applied? */
   SCIP_Bool             balanceflowsep;     /**< should flow-balance separation be applied? */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of cuts separated per separation round in the root node */
};


/**@name Local methods
 *
 * @{
 */


/** initialize (R)PC implications */
static
SCIP_RETCODE init_pcmwimplications(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
)
{
   const int nnodes = g->knots;
   const int nppterms = graph_pc_nProperPotentialTerms(g);
   int* start;
   int* verts;
   int* termmark;
   int* visitlist;
   SCIP_Real* dist;
   STP_Bool* visited;
   int nspares;
   int termscount;
   int nimplications;
   const int maxnimplications = PCIMPLICATIONS_ALLOC_FACTOR * g->edges;
   const int slotsize = ((nppterms == 0) ? 0 : maxnimplications / nppterms);

   assert(g != NULL && conshdlrdata != NULL);
   assert(graph_pc_isPcMw(g) && !g->extended);
   assert(!conshdlrdata->pcimplstart);
   assert(!conshdlrdata->pcimplverts);
   assert(0 == conshdlrdata->pcimplnppterms);

   assert(slotsize >= 1 && slotsize * nppterms <= maxnimplications);

   SCIP_CALL(SCIPallocBufferArray(scip, &dist, nnodes));
   SCIP_CALL(SCIPallocBufferArray(scip, &visited, nnodes));
   SCIP_CALL(SCIPallocBufferArray(scip, &visitlist, nnodes));
   SCIP_CALL(SCIPallocBufferArray(scip, &termmark, nnodes));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(conshdlrdata->pcimplstart), nppterms + 1));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(conshdlrdata->pcimplverts), maxnimplications));

   start = conshdlrdata->pcimplstart;
   verts = conshdlrdata->pcimplverts;
   conshdlrdata->pcimplnppterms = nppterms;

   for( int i = 0; i < nnodes; i++ )
   {
      visited[i] = FALSE;
      g->path_state[i] = UNKNOWN;
      dist[i] = FARAWAY;
   }

   graph_pc_termMarkProper(g, termmark);

   start[0] = 0;
   nspares = 0;
   termscount = 0;
   nimplications = 0;

   /* main loop: initialize implication lists */
   for( int i = 0; i < nnodes; i++ )
   {
      int nvisits;
      int nadded;

      if( !Is_term(g->term[i]) || graph_pc_knotIsFixedTerm(g, i) || graph_pc_termIsNonLeafTerm(g, i) )
         continue;

      assert(i != g->source);
      assert(g->path_heap && g->path_state);

      (void) graph_sdWalksConnected(scip, g, termmark, g->cost, NULL, i, 1000, dist, visitlist, &nvisits, visited, TRUE);

      assert(nvisits >= 1 && visitlist[0] == i);
      assert(nspares >= 0);

      for( int j = 1; j < MIN(nvisits, slotsize + nspares + 1); j++ )
      {
         const int vert = visitlist[j];
         assert(nimplications < maxnimplications);

         verts[nimplications++] = vert;
      }

      nadded = nimplications - start[termscount];
      assert(nadded >= 0);

      if( nadded > slotsize )
         nspares -= nadded - slotsize;
      else
         nspares += slotsize - nadded;

      assert(termscount < nppterms);
      start[++termscount] = nimplications;
   }
   assert(termscount == nppterms);

#ifndef WITH_UG
   printf("number of implications %d \n", nimplications);
#endif

   SCIPfreeBufferArray(scip, &termmark);
   SCIPfreeBufferArray(scip, &visitlist);
   SCIPfreeBufferArray(scip, &visited);
   SCIPfreeBufferArray(scip, &dist);

   return SCIP_OKAY;
}

/** returns inconing flow for given node */
static
SCIP_Real get_inflow(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      xval,               /**< edge values */
   int                   vert                /**< the vertex */
)
{
   double insum = 0.0;

   for( int e = g->inpbeg[vert]; e != EAT_LAST; e = g->ieat[e] )
      insum += xval[e];

   return insum;
}

/** add a cut */
static
SCIP_RETCODE cut_add(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      xval,               /**< edge values */
   int*                  capa,               /**< edges capacities (scaled) */
   const int             updatecapa,         /**< update capacities? */
   int*                  ncuts,              /**< pointer to store number of cuts */
   SCIP_Bool             local,              /**< is the cut local? */
   SCIP_Bool*            success             /**< pointer to store whether add cut be added */
   )
{
   SCIP_ROW* row;
   SCIP_VAR** vars = SCIPprobdataGetVars(scip);
   SCIP_Real sum = 0.0;
   SCIP_Bool inccapa = FALSE;
   unsigned int i;
   int* gmark = g->mark;
   int* ghead = g->head;
   int* gtail = g->tail;
   unsigned int nedges = (unsigned int) g->edges;

   assert(g->knots > 0);

   (*success) = FALSE;

   assert(g != NULL);
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "twocut", 1.0, SCIPinfinity(scip), local, FALSE, TRUE) );

   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   assert(gmark[g->source]);

   for( i = 0; i < nedges; i++ )
   {
      if( !gmark[ghead[i]] && gmark[gtail[i]] )
      {
         if( updatecapa )
         {
            if( capa[i] < FLOW_FACTOR )
               inccapa = TRUE;

            capa[i] = FLOW_FACTOR;

            if( !inccapa )
            {
               SCIP_CALL(SCIPflushRowExtensions(scip, row));
               SCIP_CALL(SCIPreleaseRow(scip, &row));
               return SCIP_OKAY;
            }
         }

         if( xval != NULL )
         {
            sum += xval[i];

            if( SCIPisFeasGE(scip, sum, 1.0) )
            {
               SCIP_CALL(SCIPflushRowExtensions(scip, row));
               SCIP_CALL(SCIPreleaseRow(scip, &row));
               return SCIP_OKAY;
            }
         }
         SCIP_CALL(SCIPaddVarToRow(scip, row, vars[i], 1.0));
      }
   }

   assert(sum < 1.0);

   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   /* checks whether cut is sufficiently violated */
   if( SCIPisCutEfficacious(scip, NULL, row) )
   {
      SCIP_Bool infeasible;

      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

      SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

#if ADDCUTSTOPOOL
      /* if at root node, add cut to pool */
      if( !infeasible )
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif
      (*ncuts)++;
      (*success) = TRUE;
   }

   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}


static
int graph_next_term(
   const GRAPH*          g,                  /**< graph data structure */
   int                   terms,              /**< number of terminals */
   int*                  term,               /**< terminal array */
   const int*            w,                  /**< awake level */
   const SCIP_Bool       firstrun            /**< first run?  */
   )
{
   int i;
   int k;
   int t;
   int wmax;
   int mindist = g->knots + 1;

   assert(term != NULL);

   if( firstrun ) // todo randomize?
   {
      assert(w[term[terms - 1]] == 0);
      return term[terms - 1];
   }

   k = -1;

   for( i = 0; (i < terms); i++ )
   {
      assert(w[term[i]] >= 0);

      if( w[term[i]] == 0 )
      {
         assert(g->mincut_dist[term[i]] < g->knots + 1);

         if( g->mincut_dist[term[i]] < mindist )
         {
            k = i;
            mindist = g->mincut_dist[term[i]];
         }
      }
   }

   if( k == -1 )
   {
      wmax = 0;

      for( i = 0; (i < terms); i++ )
      {
         if( w[term[i]] > wmax )
         {
            k = i;
            wmax = w[term[i]];
            mindist = g->mincut_dist[term[i]];
         }
         else if( w[term[i]] == wmax && g->mincut_dist[term[i]] < mindist )
         {
            assert(wmax != 0);

            k = i;
            mindist = g->mincut_dist[term[i]];
         }
      }
   }

   assert(k >= 0);
   assert(k < terms);

   t       = term[k];
   term[k] = term[terms - 1];

   return t;
}

static
void set_capacity(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Bool       creep_flow,         /**< creep flow? */
   const int             flip,               /**< reverse the flow? */
   int*                  capa,               /**< edges capacities (scaled) */
   const SCIP_Real*      xval                /**< edge values */
   )
{
   int k;
   int krev;
   int nedges = g->edges;

   assert(g    != NULL);
   assert(xval != NULL);

   for( k = 0; k < nedges; k += 2 )
   {
      krev = k + 1;
      if( !flip )
      {
         capa[k]     = (int)(xval[k    ]
            * FLOW_FACTOR + 0.5);
         capa[krev] = (int)(xval[krev]
            * FLOW_FACTOR + 0.5);
      }
      else
      {
         capa[k]     = (int)(xval[krev]
            * FLOW_FACTOR + 0.5);
         capa[krev] = (int)(xval[k    ]
            * FLOW_FACTOR + 0.5);
      }

      if( creep_flow )
      {
         capa[k] += CREEP_VALUE;
         capa[krev] += CREEP_VALUE;
      }
   }
}

/** separate PCSPG/MWCS implications */
static
SCIP_RETCODE sep_implicationsPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   int                   maxcuts,            /**< maximal number of cuts */
   int*                  ncuts               /**< pointer to store number of cuts */
   )
{
   GRAPH* g = SCIPprobdataGetGraph2(scip);
   SCIP_Real* nodeinflow;
   SCIP_VAR** vars = SCIPprobdataGetVars(scip);
   SCIP_ROW* row = NULL;
   const SCIP_Real* xval = SCIPprobdataGetXval(scip, NULL);
   int* verts;
   int* start;
   const int nnodes = g->knots;
   int cutscount;
   int ptermcount;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conshdlrdata != NULL);
   assert(g != NULL);
   assert(xval != NULL);

   /* nothing to separate? */
   if( graph_pc_nNonFixedTerms(g) == 0 )
      return SCIP_OKAY;

   /* initialize? */
   if( conshdlrdata->pcimplstart == NULL )
   {
      assert(conshdlrdata->pcimplverts == NULL);
      graph_pc_2org(scip, g);
      SCIP_CALL( init_pcmwimplications(scip, g, conshdlrdata) );
      graph_pc_2trans(scip, g);
   }

   verts = conshdlrdata->pcimplverts;
   start = conshdlrdata->pcimplstart;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodeinflow, nnodes) );

   /* initialize node sums */
   for( int i = 0; i < nnodes; i++ )
      nodeinflow[i] = get_inflow(g, xval, i);

   cutscount = 0;
   ptermcount = 0;

   assert(g->extended);

   /* main separation loop */
   for( int i = 0; i < nnodes; i++ )
   {
      int maxnode;
      SCIP_Real maxflow;
      const SCIP_Real inflow = nodeinflow[i];

      if( !Is_pseudoTerm(g->term[i]) )
         continue;

      ptermcount++;

      if( SCIPisFeasGE(scip, inflow, 1.0) )
         continue;

      maxnode = -1;
      maxflow = 0.0;
      for( int j = start[ptermcount - 1]; j < start[ptermcount]; j++ )
      {
         const int vert = verts[j];
         if( SCIPisFeasGT(scip, nodeinflow[vert], inflow) && nodeinflow[vert] > maxflow )
         {
            maxnode = vert;
            maxflow = nodeinflow[vert];
         }
      }

      /* separate? */
      if( maxnode >= 0 )
      {
         SCIP_Bool infeasible;


#if 0
         SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "pcimplicate", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE));
         SCIP_CALL(SCIPcacheRowExtensions(scip, row));

         for( int e = g->inpbeg[maxnode]; e != EAT_LAST; e = g->ieat[e] )
            SCIP_CALL(SCIPaddVarToRow(scip, row, vars[e], 1.0));

         for( int e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
            SCIP_CALL(SCIPaddVarToRow(scip, row, vars[e], -1.0));
#else
         {
            const int twinterm = graph_pc_getTwinTerm(g, i);
            const int rootedge = graph_pc_getRoot2PtermEdge(g, twinterm);
            assert(rootedge >= 0);

            SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "pcimplicate", -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE));
            SCIP_CALL(SCIPcacheRowExtensions(scip, row));

            for( int e = g->inpbeg[maxnode]; e != EAT_LAST; e = g->ieat[e] )
               SCIP_CALL(SCIPaddVarToRow(scip, row, vars[e], 1.0));

            SCIP_CALL(SCIPaddVarToRow(scip, row, vars[rootedge], 1.0));
         }
#endif

         SCIP_CALL(SCIPflushRowExtensions(scip, row));

         SCIP_CALL(SCIPaddRow(scip, row, FALSE, &infeasible));

#if ADDCUTSTOPOOL
         /* add cut to pool */
         if( !infeasible )
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

         SCIP_CALL(SCIPreleaseRow(scip, &row));

         if( *ncuts + cutscount++ >= maxcuts )
            break;
      }
   }
   assert((*ncuts + cutscount > maxcuts) || ptermcount == graph_pc_nProperPotentialTerms(g));
   assert((*ncuts + cutscount > maxcuts) || ptermcount == conshdlrdata->pcimplnppterms);


   *ncuts += cutscount;
   SCIPdebugMessage("PcImplication Separator: %d Inequalities added\n", cutscount);

   SCIPfreeBufferArray(scip, &nodeinflow);

   return SCIP_OKAY;
}


#if 0
/** separate degree-2 cuts */
static
SCIP_RETCODE sep_deg2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   maxcuts,            /**< maximal number of cuts */
   int*                  ncuts               /**< pointer to store number of cuts */
   )
{
   GRAPH* g;
   SCIP_VAR** vars;
   SCIP_ROW* row = NULL;
   SCIP_Real* xval;
   int cutscount = 0;
   int nnodes;
   const SCIP_Bool* deg2bounded = SCIPStpPropGet2BoundedArr(scip);

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conshdlrdata != NULL);
   assert(deg2bounded != NULL);

   vars = SCIPprobdataGetVars(scip);
   g = consdata->graph;
   assert(g != NULL);

   xval = SCIPprobdataGetXval(scip, NULL);
   assert(xval != NULL);

   nnodes = g->knots;

   for( int i = 0; i < nnodes; i++ )
   {
      double inoutsum;

      if( Is_term(g->term[i]) )
         continue;

      if( !deg2bounded[i] )
         continue;

      inoutsum = 0.0;

      for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         inoutsum += xval[e] + xval[flipedge_Uint(e)];
         assert(flipedge_Uint(e) == (unsigned) flipedge(e));
      }

      if( SCIPisFeasGT(scip, inoutsum, 2.0) )
      {
         SCIP_Bool infeasible;

         SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "deg2", -SCIPinfinity(scip), 2.0, FALSE, FALSE, TRUE));

         SCIP_CALL(SCIPcacheRowExtensions(scip, row));

         for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         {
            SCIP_CALL(SCIPaddVarToRow(scip, row, vars[e], 1.0));
            SCIP_CALL(SCIPaddVarToRow(scip, row, vars[flipedge_Uint(e)], 1.0));
            assert(flipedge_Uint(e) == (unsigned) flipedge(e));
         }

         SCIP_CALL(SCIPflushRowExtensions(scip, row));

         SCIP_CALL(SCIPaddRow(scip, row, FALSE, &infeasible));

#if ADDCUTSTOPOOL
         /* add cut to pool */
         if( !infeasible )
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

         cutscount++;

         SCIP_CALL(SCIPreleaseRow(scip, &row));

         if( *ncuts + cutscount >= maxcuts )
            break;
      }
   }

   printf("Deg2 Separator: %d Inequalities added\n", cutscount);
   *ncuts += cutscount;

   return SCIP_OKAY;
}
#endif



/** separate in-flow cuts:
 *  input of a non-terminal vertex has to be <= 1.0 */
static
SCIP_RETCODE sep_flowIn(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      xval,               /**< LP-solution values */
   int                   vertex,             /**< vertex */
   SCIP_VAR**            vars,               /**< variables */
   int*                  cutcount            /**< counts cuts */
)
{
   SCIP_Real sum = 0.0;

   assert(xval && cutcount && vars);
   assert(!Is_term(g->term[vertex]));

   for( int k = g->inpbeg[vertex]; k != EAT_LAST; k = g->ieat[k] )
      sum += xval[k];

   if( SCIPisFeasGT(scip, sum, 1.0) )
   {
      SCIP_ROW* row = NULL;
      SCIP_Bool infeasible;

      SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "inflow", -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE));
      SCIP_CALL(SCIPcacheRowExtensions(scip, row));

      for( int k = g->inpbeg[vertex]; k != EAT_LAST; k = g->ieat[k] )
         SCIP_CALL(SCIPaddVarToRow(scip, row, vars[k], 1.0));

      SCIP_CALL(SCIPflushRowExtensions(scip, row));
      SCIP_CALL(SCIPaddRow(scip, row, FALSE, &infeasible));

#if ADDCUTSTOPOOL
      /* if at root node, add cut to pool */
      if( !infeasible )
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

      (*cutcount)++;

      SCIP_CALL(SCIPreleaseRow(scip, &row));
   }

   return SCIP_OKAY;
}


/** separate terminal in-flow cuts
 *  at terminal input sum == 1
 *  basically a cut (starcut) */
static
SCIP_RETCODE sep_flowTermIn(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      xval,               /**< LP-solution values */
   int                   vertex,             /**< vertex */
   SCIP_VAR**            vars,               /**< variables */
   int*                  cutcount            /**< counts cuts */
)
{
   SCIP_Real sum = 0.0;

   assert(Is_term(g->term[vertex]));

   for( int k = g->inpbeg[vertex]; k != EAT_LAST; k = g->ieat[k] )
      sum += xval[k];

   if( !SCIPisFeasEQ(scip, sum, 1.0) )
   {
      SCIP_ROW* row = NULL;
      SCIP_Bool infeasible;

      SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "term", 1.0, 1.0, FALSE, FALSE, TRUE));

      SCIP_CALL(SCIPcacheRowExtensions(scip, row));

      for( int k = g->inpbeg[vertex]; k != EAT_LAST; k = g->ieat[k] )
         SCIP_CALL(SCIPaddVarToRow(scip, row, vars[k], 1.0));

      SCIP_CALL(SCIPflushRowExtensions(scip, row));
      SCIP_CALL(SCIPaddRow(scip, row, FALSE, &infeasible));

#if ADDCUTSTOPOOL
      /* add cut to pool */
      if( !infeasible )
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

      (*cutcount)++;

      SCIP_CALL(SCIPreleaseRow(scip, &row));
   }

   return SCIP_OKAY;
}


/** separate flow-balance constraints
 *  incoming flow <= outgoing flow */
static
SCIP_RETCODE sep_flowBalance(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      xval,               /**< LP-solution values */
   int                   vertex,             /**< vertex */
   SCIP_VAR**            vars,               /**< variables */
   int*                  cutcount            /**< counts cuts */
)
{
   SCIP_ROW* row = NULL;
   SCIP_Real sum = 0.0;

   assert(!Is_term(g->term[vertex]));

   for( int k = g->inpbeg[vertex]; k != EAT_LAST; k = g->ieat[k] )
      sum -= xval[k];

   for( int k = g->outbeg[vertex]; k != EAT_LAST; k = g->oeat[k] )
      sum += xval[k];

   if( SCIPisFeasNegative(scip, sum) )
   {
      SCIP_Bool infeasible;

      SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "flowbalance", 0.0, (g->terms == 2) ? 0.0 : SCIPinfinity(scip), FALSE, FALSE, TRUE));
      SCIP_CALL(SCIPcacheRowExtensions(scip, row));

      for( int k = g->inpbeg[vertex]; k != EAT_LAST; k = g->ieat[k] )
         SCIP_CALL(SCIPaddVarToRow(scip, row, vars[k], -1.0));

      for( int k = g->outbeg[vertex]; k != EAT_LAST; k = g->oeat[k] )
         SCIP_CALL(SCIPaddVarToRow(scip, row, vars[k], 1.0));

      SCIP_CALL(SCIPflushRowExtensions(scip, row));
      SCIP_CALL(SCIPaddRow(scip, row, FALSE, &infeasible));

#if ADDCUTSTOPOOL
      /* if at root node, add cut to pool */
      if( !infeasible )
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

      (*cutcount)++;

      SCIP_CALL(SCIPreleaseRow(scip, &row));
   }

   return SCIP_OKAY;
}


/** separate
 * the value of each outgoing edge needs to be smaller than the sum of the in-going edges */
static
SCIP_RETCODE sep_flowEdgeOut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      xval,               /**< LP-solution values */
   int                   vertex,             /**< vertex */
   SCIP_VAR**            vars,               /**< variables */
   int*                  cutcount            /**< counts cuts */
)
{
   const int i = vertex;

   for( int ijedge = g->outbeg[i]; ijedge != EAT_LAST; ijedge = g->oeat[ijedge] )
   {
      const int j = g->head[ijedge];
      SCIP_Real sum = -xval[ijedge];

      for( int e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
         if( g->tail[e] != j )
            sum += xval[e];

      if( SCIPisFeasNegative(scip, sum) )
      {
         SCIP_Bool infeasible;
         SCIP_ROW* row = NULL;

         SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "edgeflow", 0.0, SCIPinfinity(scip), FALSE, FALSE, TRUE));
         SCIP_CALL(SCIPcacheRowExtensions(scip, row));

         SCIP_CALL(SCIPaddVarToRow(scip, row, vars[ijedge], -1.0));

         for( int e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
         {
            if( g->tail[e] != j )
               SCIP_CALL(SCIPaddVarToRow(scip, row, vars[e], 1.0));
         }

         SCIP_CALL(SCIPflushRowExtensions(scip, row));
         SCIP_CALL(SCIPaddRow(scip, row, FALSE, &infeasible));

#if ADDCUTSTOPOOL
         /* add cut to pool */
         if( !infeasible )
            SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

         (*cutcount)++;
         SCIP_CALL(SCIPreleaseRow(scip, &row));
      }
   }

   return SCIP_OKAY;
}


/** separate flow-cuts */
static
SCIP_RETCODE sep_flow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   maxcuts,            /**< maximal number of cuts */
   int*                  ncuts               /**< pointer to store number of cuts */
   )
{
   const GRAPH* g = consdata->graph;
   SCIP_VAR** vars;
   SCIP_Real* xval;
   int count = 0;
   const SCIP_Bool flowsep = conshdlrdata->flowsep;
   const SCIP_Bool inflowsep = conshdlrdata->inflowsep;
   const SCIP_Bool intermflowsep = conshdlrdata->intermflowsep;
   const SCIP_Bool outflowsep = conshdlrdata->outflowsep;
   const SCIP_Bool balanceflowsep = conshdlrdata->balanceflowsep;

   assert(scip && conshdlr && g);

   vars = SCIPprobdataGetVars(scip);
   xval = SCIPprobdataGetXval(scip, NULL);
   assert(xval);

   for( int i = 0; i < g->knots; i++ )
   {
      if( i == g->source )
         continue;

      if( intermflowsep && Is_term(g->term[i]) )
      {
         SCIP_CALL( sep_flowTermIn(scip, conshdlr, g, xval, i, vars, &count) );

         if( *ncuts + count >= maxcuts )
            break;
      }

      /* flow cuts disabled? */
      if( !flowsep )
         continue;

      if( outflowsep )
      {
         SCIP_CALL( sep_flowEdgeOut(scip, conshdlr, g, xval, i, vars, &count) );

         if( *ncuts + count >= maxcuts )
            break;
      }

      /* from here on consider only non terminals */
      if( Is_term(g->term[i]) )
         continue;

      if( inflowsep )
      {
         SCIP_CALL( sep_flowIn(scip, conshdlr, g, xval, i, vars, &count) );

         if( *ncuts + count >= maxcuts )
            break;
      }

      if( balanceflowsep )
      {
         SCIP_CALL( sep_flowBalance(scip, conshdlr, g, xval, i, vars, &count) );

         if( *ncuts + count >= maxcuts )
            break;
      }
   }

   SCIPdebugMessage("Flow Separator: %d Inequalities added\n", count);

   *ncuts += count;

   return SCIP_OKAY;
}

/** separate 2-cuts */
static
SCIP_RETCODE sep_2cut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   const int*            termorg,            /**< original terminals or NULL */
   int                   maxcuts,            /**< maximal number of cuts */
   int*                  ncuts               /**< pointer to store number of cuts */
   )
{
#if 0
   const SCIP_Bool nested_cut   = conshdlrdata->nestedcut;
   const SCIP_Bool back_cut     = conshdlrdata->backcut;
   const SCIP_Bool creep_flow   = conshdlrdata->creepflow;
   const SCIP_Bool disjunct_cut = conshdlrdata->disjunctcut;
#endif
   /* we do not longer support any other flow as they slow everything down and are of little use anyway todo remove user parameter */
   const SCIP_Bool flowsep      = conshdlrdata->flowsep;
   const SCIP_Bool nested_cut   = FALSE;
   const SCIP_Bool creep_flow   = TRUE;
   const SCIP_Bool disjunct_cut = FALSE;
   const SCIP_Bool intree = (SCIPgetDepth(scip) > 0);

   SCIP_VAR** vars;
   GRAPH*  g;
   SCIP_Real* xval;
   int*    w;
   int*    capa;
   int*    term;
   int*    start;
   int*    excess;
   int*    rootcut;
   int*    edgearr;
   int*    headarr;
   int*    residual;
   int*    edgecurr;
   int*    headactive;
   int*    edgeflipped;
   int*    headinactive;
   int     i;
   int     k;
   int     e;
   int     root;
   int     head;
   int     count;
   int     terms;
   int     nedges;
   int     nnodes;
   int     newnedges;
   int     rootcutsize;
   SCIP_Bool rerun;
   SCIP_Bool addedcut;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conshdlrdata != NULL);

   g = consdata->graph;
   assert(g != NULL);

   root = g->source;
   excess = g->mincut_e;
   nedges = g->edges;
   nnodes = g->knots;
   addedcut = FALSE;
   residual = g->mincut_r;
   edgecurr = g->mincut_numb;
   headactive = g->mincut_head;
   headinactive = g->mincut_head_inact;

   assert(residual != NULL);
   assert(edgecurr != NULL);
   assert(headactive != NULL);
   assert(headinactive != NULL);

   xval = SCIPprobdataGetXval(scip, NULL);
   assert(xval != NULL);

   assert(creep_flow == TRUE);
   assert(nested_cut == FALSE);
   assert(disjunct_cut == FALSE);

   /* for 2-terminal nets no cuts are necessary if flows are given */
   if( flowsep && (g->terms == 2) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &capa, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &w, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &term, g->terms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearr, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &headarr, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgeflipped, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &start, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rootcut, nnodes + 1) );

#ifdef STP_MAXFLOW_TIME
   clock_t startt, endt;
   double cpu_time_used;
   startt = clock();
#endif

   vars = SCIPprobdataGetVars(scip);
   assert(vars != NULL);

   assert(nedges >= nnodes);

   for( k = 0; k < nnodes; k++ )
   {
      w[k] = 0;
      excess[k] = 0;
   }

   for( e = 0; e < nedges; e += 2 )
   {
      const int erev = e + 1;

      if( intree && SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[erev]) < 0.5 )
      {
         capa[e] = 0;
         capa[erev] = 0;
         residual[e] = 0;
         residual[erev] = 0;

         headarr[e] = 1;
         headarr[erev] = 1;
      }
      else
      {
         capa[e]     = (int)(xval[e] * FLOW_FACTOR + 0.5) + CREEP_VALUE;
         capa[erev]  = (int)(xval[erev] * FLOW_FACTOR + 0.5) + CREEP_VALUE;
         residual[e] = capa[e];
         residual[erev] = capa[erev];

         headarr[e] = SCIPisFeasLT(scip, xval[e], 1.0) ? 1 : 0;
         headarr[erev] = SCIPisFeasLT(scip, xval[erev], 1.0) ? 1 : 0;
      }
      edgearr[e] = -1;
      edgearr[erev] = -1;
   }

   /*
    * bfs along 0 edges from the root
    * */

   w[root] = 1;
   rootcutsize = 0;
   rootcut[rootcutsize++] = root;

   /* bfs loop */
   for( i = 0; i < rootcutsize; i++ )
   {
      assert(rootcutsize <= nnodes);

      k = rootcut[i];

      assert(k < nnodes);

      /* traverse outgoing arcs */
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         head = g->head[e];

         /* head not been added to root cut yet? */
         if( w[head] == 0 )
         {
            if( headarr[e] == 0 )
            {
               w[head] = 1;
               rootcut[rootcutsize++] = head;
            }
            else
            {
               /* push as much as possible out of perpetually dormant nodes (possibly to other dormant nodes) */
               assert(w[head] == 0);
#ifndef NDEBUG
               residual[e] = 0;
#endif
               excess[head] += capa[e];
            }
         }
      }
   }

   i = 0;
   terms = 0;

   /* fill auxiliary adjacent vertex/edges arrays and get useable terms */
   for( k = 0; k < nnodes; k++ )
   {
      headactive[k] = Q_NULL;
      headinactive[k] = Q_NULL;

      start[k] = i;

      /* non-dormant node? */
      if( w[k] == 0 )
      {
         edgecurr[k] = i;
         for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            if( w[g->head[e]] == 0 && capa[e] != 0 )
            {
               edgearr[e] = i;
               residual[i] = capa[e];
               headarr[i++] = g->head[e];
            }
         }

         /* unreachable node? */
         if( edgecurr[k] == i )
         {
            w[k] = 1;
         }
         else if( Is_term(g->term[k]) )
         {
            term[terms++] = k;
         }
      }
      else
      {
         edgecurr[k] = -1;
      }
   }

   newnedges = i;
   start[nnodes] = i;

   /* initialize edgeflipped */
   for( e = nedges - 1; e >= 0; e-- )
   {
      if( edgearr[e] >= 0 )
      {
         i = edgearr[e];
         edgeflipped[i] = edgearr[flipedge(e)];
      }
   }

   SCIPdebugMessage("Cut Pretest: %d eliminations\n", g->terms - terms - 1);

#ifdef STP_MAXFLOW_TIME
   endt = clock();
   cpu_time_used = ((double) (endt - startt)) / CLOCKS_PER_SEC;
   startt = clock();
#endif

   count = 0;
   rerun = FALSE;

   while( terms > 0 )
   {
      if( ((unsigned) terms) % 32 == 0 && SCIPisStopped(scip) )
         break;

      /* look for reachable terminal */
      i = graph_next_term(g, terms, term, w, !rerun);

      terms--;

      assert(g->term[i]       == 0);
      assert(g->source != i);

      if( nested_cut && !disjunct_cut )
         set_capacity(g, creep_flow, 0, capa, xval);

      do
      {
#ifdef STP_MAXFLOW_WRITE
         /* write flow problem in extended dimacs format */
         FILE *fptr;

         fptr = fopen("flow", "w+");
         assert(fptr != NULL);

         fprintf(fptr, "p max %d %d \n", nnodes, nedges);
         fprintf(fptr, "n %d s \n", g->source + 1);
         fprintf(fptr, "n %d t \n", i + 1);

         for( k = 0; k < nnodes; k++ )
         {
            for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
            {
               fprintf(fptr, "a %d %d %d \n", k + 1, g->head[e] + 1, capa[e]);
            }
         }

         fprintf(fptr, "x\n");

         fclose(fptr);
#endif
         // declare cuts on branched-on (artificial) terminals as local
         const SCIP_Bool localcut = (termorg != NULL && termorg[i] != g->term[i]);

         /* non-trivial cut? */
         if( w[i] != 1 )
         {
            graph_mincut_exec(g, root, i, nnodes, newnedges, rootcutsize, rootcut, capa, w, start, edgeflipped, headarr, rerun);

            /* cut */
            for( k = nnodes - 1; k >= 0; k-- )
               g->mark[k] = (w[k] != 0);

            assert(g->mark[root]);
         }
         else
         {
            SCIP_Real flowsum = 0.0;

            assert(rerun);

            for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
               flowsum += xval[e];

            if( SCIPisFeasGE(scip, flowsum, 1.0) )
               continue;

            for( k = nnodes - 1; k >= 0; k-- )
               g->mark[k] = TRUE;

            g->mark[i] = FALSE;
         }

         rerun = TRUE;

         SCIP_CALL( cut_add(scip, conshdlr, g, xval, capa, nested_cut || disjunct_cut, ncuts, localcut, &addedcut) );
         if( addedcut )
         {
            count++;

            if( *ncuts >= maxcuts )
               goto TERMINATE;
         }
         else
            break;
      }
      while( nested_cut );               /* Nested Cut is CONSTANT ! */
   } /* while terms > 0 */


#ifdef STP_MAXFLOW_TIME
   endt = clock();
   cpu_time_used = ((double) (endt - startt)) / CLOCKS_PER_SEC;
#endif

#if 0
      /*
       * back cuts currently not supported
       *  */
      /* back cuts enabled? */
      if( back_cut )
      {
         for( k = 0; k < nnodes; k++ )
            w[k] = 0;

         if( !nested_cut || disjunct_cut )
            set_capacity(g, creep_flow, 1, capa, xval);

         terms = tsave;

         while( terms > 0 )
         {
            /* look for reachable terminal */
            i = graph_next_term(g, terms, term, w, TRUE);

            terms--;

            assert(g->term[i]       == 0);
            assert(g->source != i);

            if( nested_cut && !disjunct_cut )
               set_capacity(g, creep_flow, 1, capa, xval);

            rerun = FALSE;

            do
            {
               graph_mincut_exec(g, i, g->source, nedges, capa, w, start, edgearr, headarr, rerun);

               rerun = TRUE;

               for( k = 0; k < nnodes; k++ )
               {
                  g->mark[k] = (w[k] != 0) ? 0 : 1; // todo not the other way around??
                  w[k] = 0;
               }

               SCIP_CALL( cut_add(scip, conshdlr, g, xval, capa, nested_cut || disjunct_cut, ncuts, &addedcut) );
               if( addedcut )
               {
                  count++;

                  if( *ncuts >= maxcuts )
                     goto TERMINATE;
               }
               else
                  break;
#if 0
               if (nested_cut || disjunct_cut)
                  for(k = p->beg[p->rcnt - 1]; k < p->nzcnt; k++)
                     capa[p->ind[k] % nedges
                        + (((p->ind[k] % nedges) % 2)
                           ? -1 : 1)] = FLOW_FACTOR;
#endif
            }
            while( nested_cut );                /* Nested Cut is CONSTANT todo why not only one round? seems to make no sense whatsoever */

            rerun = FALSE;
         }
      }
#endif
 TERMINATE:
   SCIPfreeBufferArray(scip, &rootcut);
   SCIPfreeBufferArray(scip, &start);
   SCIPfreeBufferArray(scip, &edgeflipped);
   SCIPfreeBufferArray(scip, &headarr);
   SCIPfreeBufferArray(scip, &edgearr);

   SCIPfreeBufferArray(scip, &term);
   SCIPfreeBufferArray(scip, &w);

   SCIPfreeBufferArray(scip, &capa);

   SCIPdebugMessage("2-cut Separator: %d Inequalities added\n", count);

   return SCIP_OKAY;
}



/**@} */


/**@name Callback methods
 *
 * @{
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyStp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrStp(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeStp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolStp)
{  /*lint --e{715}*/
#ifdef WITH_UG
   SCIPStpConshdlrSetGraph(scip, SCIPprobdataGetGraph2(scip));
#endif
   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXITSOL(consExitsolStp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemoryArrayNull(scip, &(conshdlrdata->pcimplstart));
   SCIPfreeMemoryArrayNull(scip, &(conshdlrdata->pcimplverts));

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteStp)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransStp)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create constraint data for target constraint */
   SCIP_CALL( SCIPallocBlockMemory(scip, &targetdata) );

   targetdata->graph = sourcedata->graph;

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpStp)
{  /*lint --e{715}*/
#if 0
   SCIP_PROBDATA* probdata;
   GRAPH* graph;

   SCIP_Real lpobjval;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   SCIP_CALL( SCIPdualAscentPcStp(scip, graph, NULL, &lpobjval, TRUE, 1) );
#endif

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpStp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   GRAPH* g;
   int* termorg = NULL;
   int* nodestatenew = NULL;
   int maxcuts;
   int ncuts = 0;
   const SCIP_Bool atrootnode = (SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0);
   SCIP_Bool chgterms;
#ifndef NDEBUG
   int nterms;
#endif

   *result = SCIP_DIDNOTRUN;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   maxcuts = atrootnode ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts;

   assert(nconss == 1);
   consdata = SCIPconsGetData(conss[0]);

   assert(consdata != NULL);

   g = consdata->graph;
   assert(g != NULL);

#ifndef NDEBUG
   nterms = g->terms;
#endif

   chgterms = (!atrootnode && (graph_typeIsSpgLike(g) || graph_pc_isPcMw(g)));

   SCIP_CALL( sep_flow(scip, conshdlr, conshdlrdata, consdata, maxcuts, &ncuts) );

   if( graph_pc_isPcMw(g) && g->stp_type != STP_BRMWCSP )
      SCIP_CALL( sep_implicationsPcMw(scip, conshdlr, conshdlrdata, maxcuts, &ncuts) );

   /* change graph according to branch-and-bound terminal changes  */
   if( chgterms )
   {
      SCIP_Bool conflict = FALSE;
      const int nnodes = g->knots;

      SCIP_CALL(SCIPallocBufferArray(scip, &nodestatenew, nnodes));
      SCIP_CALL(SCIPallocBufferArray(scip, &termorg, nnodes));
      BMScopyMemoryArray(termorg, g->term, nnodes);

      SCIPStpBranchruleInitNodeState(g, nodestatenew);
      SCIP_CALL( SCIPStpBranchruleGetVertexChgs(scip, nodestatenew, &conflict) );

      assert(!conflict);

      for( int k = 0; k < nnodes; k++ )
      {
         if( nodestatenew[k] == BRANCH_STP_VERTEX_TERM && !Is_term(g->term[k]) )
            graph_knot_chg(g, k, STP_TERM);
      }
   }

   SCIP_CALL( sep_2cut(scip, conshdlr, conshdlrdata, consdata, termorg, maxcuts, &ncuts) );

   if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   /* restore graph */
   if( chgterms )
   {
      const int nnodes = g->knots;

      for( int k = 0; k < nnodes; k++ )
      {
         if( g->term[k] != termorg[k] )
            graph_knot_chg(g, k, termorg[k]);
      }
   }

#ifndef NDEBUG
   assert(g->terms == nterms);
#endif

   SCIPfreeBufferArrayNull(scip, &termorg);
   SCIPfreeBufferArrayNull(scip, &nodestatenew);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpStp)
{  /*lint --e{715}*/
   SCIP_Bool feasible;
   SCIP_CONSDATA* consdata;

   for( int i = 0; i < nconss; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);

      SCIP_CALL( SCIPStpValidateSol(scip, consdata->graph, SCIPprobdataGetXval(scip, NULL), FALSE, &feasible) );

      if( !feasible )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }

   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsStp)
{  /*lint --e{715}*/
   SCIP_Bool feasible;

   assert(nconss == 1);

   for( int i = 0; i < nconss; i++ )
   {
      const SCIP_CONSDATA* consdata = SCIPconsGetData(conss[i]);

      SCIP_CALL( SCIPStpValidateSol(scip, consdata->graph, SCIPprobdataGetXval(scip, NULL), FALSE, &feasible) );

      if( !feasible )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckStp)
{ /*lint --e{715}*/
   const GRAPH* g = SCIPprobdataGetGraph2(scip);
   SCIP_Bool feasible;

   assert(g != NULL);

   SCIP_CALL(SCIPStpValidateSol(scip, g, SCIPprobdataGetXval(scip, sol), FALSE, &feasible));

   if( !feasible )
   {
      *result = SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropStp)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   GRAPH* graph;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   /* for degree constrained model, check whether problem is infeasible */
   if( graph->stp_type == STP_DCSTP )
   {
      int k;
      int nnodes;
      int degsum;
      int* maxdegs;

      nnodes = graph->knots;
      maxdegs = graph->maxdeg;

      assert(maxdegs != NULL);

      degsum = 0;
      for( k = 0; k < nnodes; k++ )
      {
         if( Is_term(graph->term[k]) )
         {
            assert(maxdegs[k] > 0);
            degsum += maxdegs[k] - 1;
         }
         else
         {
            assert(maxdegs[k] >= 0);
            degsum += MAX(maxdegs[k] - 2, 0);
         }
      }

      if( degsum < graph->terms - 2 )
         *result = SCIP_CUTOFF;
      else
	 *result = SCIP_DIDNOTFIND;
   }
   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockStp)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   vars = SCIPprobdataGetVars(scip);
   nvars = SCIPprobdataGetNVars(scip);

   for( v = 0; v < nvars; ++v )
      SCIP_CALL( SCIPaddVarLocksType(scip, vars[v], SCIP_LOCKTYPE_MODEL, 1, 1) );

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyStp)
{  /*lint --e{715}*/
   const char* consname;
   SCIP_PROBDATA* probdata;
   GRAPH* graph;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   consname = SCIPconsGetName(sourcecons);

   /* creates and captures a and constraint */
   SCIP_CALL( SCIPcreateConsStp(scip, cons, consname, graph) );

   *valid = TRUE;

   return SCIP_OKAY;
}


/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the handler for stp constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrStp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create stp constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );

   conshdlr = NULL;
   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpStp, consEnfopsStp, consCheckStp, consLockStp,
         conshdlrdata) );
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyStp, consCopyStp) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteStp) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransStp) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolStp) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolStp) );

   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpStp) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropStp, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpStp, NULL, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeStp) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/stp/backcut", "Try Back-Cuts",
         &conshdlrdata->backcut, TRUE, DEFAULT_BACKCUT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/stp/creepflow", "Use Creep-Flow",
         &conshdlrdata->creepflow, TRUE, DEFAULT_CREEPFLOW, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/stp/disjunctcut", "Only disjunct Cuts",
         &conshdlrdata->disjunctcut, TRUE, DEFAULT_DISJUNCTCUT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/stp/nestedcut", "Try Nested-Cuts",
         &conshdlrdata->nestedcut, TRUE, DEFAULT_NESTEDCUT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/stp/flowsep", "Try Flow-Cuts",
         &conshdlrdata->flowsep, TRUE, DEFAULT_FLOWSEP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/stp/inflowsep", "Try Unit Inflow-Cuts",
         &conshdlrdata->inflowsep, TRUE, DEFAULT_INFLOWSEP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/stp/intermflowsep", "Try terminal Unit Inflow-Cuts",
         &conshdlrdata->intermflowsep, TRUE, DEFAULT_INFLOWTERMSEP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/stp/outflowsep", "Try single edge Outflow-Cuts",
         &conshdlrdata->outflowsep, TRUE, DEFAULT_OUTFLOWSEP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/stp/balanceflowsep", "Try Flow-balance Cuts",
         &conshdlrdata->balanceflowsep, TRUE, DEFAULT_BALANCEFLOWSEP, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         &conshdlrdata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxroundsroot",
         "maximal number of separation rounds per node in the root node (-1: unlimited)",
         &conshdlrdata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxsepacuts",
         "maximal number of cuts separated per separation round",
         &conshdlrdata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxsepacutsroot",
         "maximal number of cuts separated per separation round in the root node",
         &conshdlrdata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );

   conshdlrdata->pcimplstart = NULL;
   conshdlrdata->pcimplverts = NULL;
   conshdlrdata->pcimplnppterms = 0;

   return SCIP_OKAY;
}

/** creates and captures a stp constraint */
SCIP_RETCODE SCIPcreateConsStp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   GRAPH*                graph               /**< graph data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the stp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("stp constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->graph = graph;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, TRUE, TRUE, TRUE, TRUE,
         FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** add cut corresponding to contraction */
SCIP_RETCODE SCIPStpAddContractionCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             edge,               /**< edge */
   SCIP_VAR*             revedge,            /**< reversed edge */
   SCIP_Bool             localcut            /**< add local cut? */
)
{
   SCIP_ROW* row = NULL;
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool infeasible;

   if( SCIPvarGetLbLocal(edge) > 0.5 || SCIPvarGetUbLocal(edge) < 0.5 || SCIPvarGetLbLocal(revedge) > 0.5 || SCIPvarGetUbLocal(revedge) < 0.5 )
   {
      printf("cannot add contraction cut \n");
      return SCIP_OKAY;
   }

   conshdlr = SCIPfindConshdlr(scip, "stp");
   assert(conshdlr != NULL);
   assert(SCIPconshdlrGetNConss(conshdlr) > 0);

   SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "contraction", 1.0, SCIPinfinity(scip), localcut, FALSE, TRUE));
   SCIP_CALL(SCIPcacheRowExtensions(scip, row));

   SCIP_CALL(SCIPaddVarToRow(scip, row, edge, 1.0));
   SCIP_CALL(SCIPaddVarToRow(scip, row, revedge, 1.0));

   SCIP_CALL(SCIPflushRowExtensions(scip, row));

   SCIP_CALL(SCIPaddRow(scip, row, FALSE, &infeasible));

#if ADDCUTSTOPOOL
   /* add cut to pool */
   if( !infeasible )
   SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

   SCIP_CALL(SCIPreleaseRow(scip, &row));

   return SCIP_OKAY;
}

/** sets graph */
void SCIPStpConshdlrSetGraph(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLR* conshdlr;

   conshdlr = SCIPfindConshdlr(scip, "stp");
   assert(conshdlr != NULL);
   assert(SCIPconshdlrGetNConss(conshdlr) > 0);

   consdata = SCIPconsGetData(SCIPconshdlrGetConss(conshdlr)[0]);

   assert(consdata != NULL);

   consdata->graph = SCIPprobdataGetGraph2(scip);
   assert(consdata->graph != NULL);
}

/** returns implications start array */
int* SCIPStpGetPcImplStarts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, "stp");
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->pcimplstart;
}

/** returns number implications starts */
int SCIPStpGetPcImplNstarts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, "stp");
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->pcimplnppterms;
}


/** returns implications vertices array */
int* SCIPStpGetPcImplVerts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, "stp");
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->pcimplverts;
}


/**@} */
