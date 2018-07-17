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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_stp.c
 * @brief  Constraint handler for Steiner problems
 * @author Gerald Gamrath
 * @author Daniel Rehfeldt
 * @author Michael Winkler
 *
 * This file checks solutions for feasibility and separates violated model constraints. For more details see \ref CONS page.
 *
 * @page CONS Separating violated constraints
 *
 * In this file a constraint handler checking solutions for feasibility and separating violated model constraints is implemented.
 * The separation problem for the cut inequalities described in \ref PROBLEM can be solved by a max-flow algorithm in
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
#include "grph.h"
#include "heur_prune.h"
#include "heur_ascendprune.h"
#include "portab.h"
#include "branch_stp.h"

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
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  9999999 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             0 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ            1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define DEFAULT_MAXROUNDS            10 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS     INT_MAX /**< maximal number of cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT INT_MAX /**< maximal number of cuts separated per separation round in the root node */


#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP

#define DEFAULT_BACKCUT        FALSE /**< Try Back-Cuts FALSE*/
#define DEFAULT_CREEPFLOW      TRUE  /**< Use Creep-Flow */
#define DEFAULT_DISJUNCTCUT    FALSE /**< Only disjunct Cuts FALSE */
#define DEFAULT_NESTEDCUT      FALSE /**< Try Nested-Cuts FALSE*/
#define DEFAULT_FLOWSEP        TRUE  /**< Try Flow-Cuts */

#define DEFAULT_DAMAXDEVIATION 0.25  /**< max deviation for dual ascent */
#define DA_MAXDEVIATION_LOWER 0.01  /**< lower bound for max deviation for dual ascent */
#define DA_MAXDEVIATION_UPPER 0.9  /**< upper bound for max deviation for dual ascent */
#define DA_EPS (5.0 * 1e-7)

/* *
#define FLOW_FACTOR     100000
#define CREEP_VALUE     1         this is the original value todo check what is better
*/

#define FLOW_FACTOR     1000000
#define CREEP_VALUE     10

/* do depth-first search */
#define DFS


#ifdef BITFIELDSARRAY
#define ARRLENGTH 32
#define SetBit(Arr, pos)     ( Arr[(pos/ARRLENGTH)] |= (1 << (pos%ARRLENGTH)) )
#define CleanBit(Arr, pos)   ( Arr[(pos/ARRLENGTH)] &= ~(1 << (pos%ARRLENGTH)) )
#define BitTrue(Arr, pos)    ( Arr[(pos/ARRLENGTH)] & (1 << (pos%ARRLENGTH)) )
#endif


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
   SCIP_Bool             backcut;            /**< should backcuts be applied? */
   SCIP_Bool             creepflow;          /**< should creepflow cuts be applied? */
   SCIP_Bool             disjunctcut;        /**< should disjunction cuts be applied? */
   SCIP_Bool             nestedcut;          /**< should nested cuts be applied? */
   SCIP_Bool             flowsep;            /**< should flow separation be applied? */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of cuts separated per separation round in the root node */
};


/**@name Local methods
 *
 * @{
 */

/** returns whether node realtail is active or leads to active node other than dfsbase */
static
SCIP_Bool is_active(
   const int*            active,             /**< active nodes array */
   int                   realtail,           /**< vertex to start from */
   int                   dfsbase             /**< DFS source vertex */
   )
{
   int curr;

   for( curr = active[realtail]; curr != 0 && curr != dfsbase + 1; curr = active[curr - 1] )
   {
      assert(curr >= 0);
   }

   return (curr == 0);
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
      if( gmark[gtail[i]] && !gmark[ghead[i]] ) // todo bool array?
      {
         if( updatecapa )
         {
            if( capa[i] < FLOW_FACTOR )
               inccapa = TRUE;

            SCIPdebugMessage("set capa[%d] from %6d to %6d\n", i, capa[i], FLOW_FACTOR);
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

/** separate */
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
   GRAPH*  g;
   SCIP_VAR** vars;
   SCIP_ROW* row = NULL;
   SCIP_Real* xval;
   SCIP_Real sum;
   int    i;
   int    k;
   int    j;
   int    ind;
   int    layer;
   int    count = 0;
   unsigned int    flowsep;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conshdlrdata != NULL);

   vars = SCIPprobdataGetVars(scip);
   flowsep = conshdlrdata->flowsep;

   /* get the graph */
   g = consdata->graph;
   assert(g != NULL);

   xval = SCIPprobdataGetXval(scip, NULL);
   assert(xval != NULL);

   for(i = 0; i < g->knots; i++)
   {
      for(layer = 0; layer < g->layers; layer++)
      {
         /* continue at root */
         if( i == g->source )
            continue;

         /* at terminal: input sum == 1
          * basically a cut (starcut))
          */
         if( g->term[i] == layer )
         {
            sum = 0.0;

            for( k = g->inpbeg[i]; k != EAT_LAST; k = g->ieat[k] )
            {
               ind  = layer * g->edges + k;
               sum += (xval != NULL) ? xval[ind] : 0.0;
            }

            if( !SCIPisFeasEQ(scip, sum, 1.0) )
            {
               SCIP_Bool infeasible;

               SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "term", 1.0,
                     1.0, FALSE, FALSE, TRUE) );

               SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

               for(k = g->inpbeg[i]; k != EAT_LAST; k = g->ieat[k])
               {
                  ind  = layer * g->edges + k;

                  SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], 1.0) );
               }

               SCIP_CALL( SCIPflushRowExtensions(scip, row) );

               SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

#if ADDCUTSTOPOOL
               /* add cut to pool */
               if( !infeasible )
                  SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

               count++;

               SCIP_CALL( SCIPreleaseRow(scip, &row) );

               if( *ncuts + count >= maxcuts )
                  goto TERMINATE;
            }
         }

         /* flow cuts disabled? */
         if( !flowsep )
            continue;

         /* the value of each outgoing edge needs to be smaller than the sum of the ingoing edges */
         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
         {
            ind = layer * g->edges + j;
            sum = (xval != NULL) ? -xval[ind] : -1.0;

            for( k = g->inpbeg[i]; k != EAT_LAST; k = g->ieat[k] )
            {
               ind  = layer * g->edges + k;
               sum += (xval != NULL) ? xval[ind] : 0.0;
            }
            if( SCIPisFeasNegative(scip, sum) )
            {
               SCIP_Bool infeasible;

               SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "flow", 0.0, SCIPinfinity(scip),
                     FALSE, FALSE, TRUE) );

               SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

               ind = layer * g->edges + j;

               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], -1.0) );

               for( k = g->inpbeg[i]; k != EAT_LAST; k = g->ieat[k] )
               {
                  ind  = layer * g->edges + k;

                  SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], 1.0) );
               }

               SCIP_CALL( SCIPflushRowExtensions(scip, row) );

               SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

#if ADDCUTSTOPOOL
               /* add cut to pool */
               if( !infeasible )
                  SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

               count++;

               SCIP_CALL( SCIPreleaseRow(scip, &row) );

               if( *ncuts + count >= maxcuts )
                  goto TERMINATE;
            }
         }

         /* consider only non terminals */
         if( g->term[i] == layer )
            continue;

         /* input of a vertex has to be <= 1.0 */
         sum   = 0.0;

         for( k = g->inpbeg[i]; k != EAT_LAST; k = g->ieat[k] )
         {
            ind  = layer * g->edges + k;
            sum += (xval != NULL) ? xval[ind] : 1.0;
         }
         if( SCIPisFeasGT(scip, sum, 1.0) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "infl", -SCIPinfinity(scip),
                  1.0, FALSE, FALSE, TRUE) );

            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

            for( k = g->inpbeg[i]; k != EAT_LAST; k = g->ieat[k] )
            {
               ind  = layer * g->edges + k;

               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], 1.0) );
            }

            SCIP_CALL( SCIPflushRowExtensions(scip, row) );

            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

#if ADDCUTSTOPOOL
            /* if at root node, add cut to pool */
            if( !infeasible )
               SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

            count++;

            SCIP_CALL( SCIPreleaseRow(scip, &row) );

            if( *ncuts + count >= maxcuts )
               goto TERMINATE;
         }

         /* incoming flow <= outgoing flow */
         sum   = 0.0;

         for( k = g->inpbeg[i]; k != EAT_LAST; k = g->ieat[k] )
         {
            ind = layer * g->edges + k;
            sum -= (xval != NULL) ? xval[ind] : 1.0;
         }
         for( k = g->outbeg[i]; k != EAT_LAST; k = g->oeat[k] )
         {
            ind = layer * g->edges + k;
            sum += (xval != NULL) ? xval[ind] : 0.0;
         }
         if( SCIPisFeasNegative(scip, sum) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "bala", 0.0,
                  (g->terms == 2) ? 0.0 : SCIPinfinity(scip), FALSE, FALSE, TRUE) );

            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

            for( k = g->inpbeg[i]; k != EAT_LAST; k = g->ieat[k] )
            {
               ind = layer * g->edges + k;

               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], -1.0) );
            }
            for( k = g->outbeg[i]; k != EAT_LAST; k = g->oeat[k] )
            {
               ind = layer * g->edges + k;

               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], 1.0) );
            }

            SCIP_CALL( SCIPflushRowExtensions(scip, row) );

            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

#if ADDCUTSTOPOOL
            /* if at root node, add cut to pool */
            if( !infeasible )
               SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

            count++;

            SCIP_CALL( SCIPreleaseRow(scip, &row) );

            if( *ncuts + count >= maxcuts )
               goto TERMINATE;
         }
      }
   }

 TERMINATE:
   SCIPdebugMessage("In/Out Separator: %d Inequalities added\n", count);

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
   int     newnnodes;
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

#if 0
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
#if 1 /* for debug */
               residual[e] = 0;
#endif
               excess[head] += capa[e];
            }
         }
      }
   }

   i = 0;
   terms = 0;
   newnnodes = 0;

   /* fill auxiliary adjacent vertex/edges arrays and get useable terms */
   for( k = 0; k < nnodes; k++ )
   {
      headactive[k] = Q_NULL;
      headinactive[k] = Q_NULL;

      start[k] = i;

      /* non-dormant node? */
      if( w[k] == 0 )
      {
         newnnodes++;
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
            newnnodes--;
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

#if 0
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
#if 0
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


#if 0
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


static
SCIP_RETCODE dualascent_init(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   const int* RESTRICT   start,              /**< CSR start array [0,...,nnodes] */
   const int* RESTRICT   edgearr,            /**< CSR ancestor edge array */
   int                   root,               /**< the root */
   SCIP_Bool             is_pseudoroot,      /**< is the root a pseudo root? */
   int                   ncsredges,          /**< number of CSR edges */
   int*                  gmark,              /**< array for marking nodes */
   int* RESTRICT         active,             /**< active vertices mark */
   SCIP_PQUEUE*          pqueue,             /**< priority queue */
   GNODE**               gnodearr,           /**< array containing terminal nodes*/
   SCIP_Real* RESTRICT   rescap,             /**< residual capacity */
   SCIP_Real*            dualobj,            /**< dual objective */
   int*                  augmentingcomponent /**< augmenting component */
)
{
   const int nnodes = g->knots;
   *dualobj = 0.0;
   *augmentingcomponent = -1;

   for( int i = 0; i < ncsredges; i++ )
      rescap[i] = g->cost[edgearr[i]];

   /* mark terminals as active, add all except root to pqueue */
   for( int i = 0, k = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) )
      {
         active[i] = 0;
         assert(g->grad[i] > 0);
         if( i != root )
         {
            SCIP_Real warmstart = FALSE;
            gnodearr[k]->number = i;
            gnodearr[k]->dist = g->grad[i];

            /* for variants with dummy terminals */
            if( g->grad[i] == 2 )
            {
               int a;

               for( a = g->inpbeg[i]; a != EAT_LAST; a = g->ieat[a] )
                  if( g->cost[a] == 0.0 )
                     break;

               if( a != EAT_LAST )
               {
                  const int tail = g->tail[a];
                  gnodearr[k]->dist += g->grad[tail] - 1;

                  if( is_pseudoroot )
                  {
                     SCIP_Bool zeroedge = FALSE;
                     for( a = g->inpbeg[tail]; a != EAT_LAST; a = g->ieat[a] )
                        if( g->cost[a] == 0.0 )
                        {
                           zeroedge = TRUE;
                           gnodearr[k]->dist += g->grad[g->tail[a]] - 1;
                        }

                     /* warmstart possible? */
                     if( !zeroedge )
                     {
                        int j;
                        int end;
                        int prizearc;
                        SCIP_Real prize;

                        if( rescap[start[i]] == 0.0 )
                           prizearc = start[i] + 1;
                        else
                           prizearc = start[i];

                        prize = rescap[prizearc];
                        assert(prize > 0.0);

                        for( j = start[tail], end = start[tail + 1]; j != end; j++ )
                           if( rescap[j] < prize )
                              break;

                        if( j == end )
                        {
                           warmstart = TRUE;
                           *dualobj += prize;
                           rescap[prizearc] = 0.0;
                           for( j = start[tail], end = start[tail + 1]; j != end; j++ )
                              rescap[j] -= prize;
                        }
                     }
                  }
               }

               assert(gnodearr[k]->dist > 0);
            }
            if( !warmstart )
               SCIP_CALL(SCIPpqueueInsert(pqueue, gnodearr[k]));
            else if( *augmentingcomponent == -1 )
            {
               SCIP_CALL(SCIPpqueueInsert(pqueue, gnodearr[k]));
               *augmentingcomponent = i;
            }
            k++;
         }
      }
      else
      {
         active[i] = -1;
      }
   }

   for( int i = 0; i < nnodes + 1; i++ )
      gmark[i] = FALSE;

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

#if 1
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
#endif
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

   SCIP_CALL( sep_flow(scip, conshdlr, conshdlrdata, consdata, maxcuts, &ncuts) );

   /* change graph according to branch-and-bound terminal changes  */
   if( !atrootnode && g->stp_type == STP_SPG )
   {
      const int nnodes = g->knots;

      SCIP_CALL(SCIPallocBufferArray(scip, &nodestatenew, nnodes));
      SCIP_CALL(SCIPallocBufferArray(scip, &termorg, nnodes));
      BMScopyMemoryArray(termorg, g->term, nnodes);

      SCIPStpBranchruleInitNodeState(g, nodestatenew);
      SCIP_CALL( SCIPStpBranchruleApplyVertexChgs(scip, nodestatenew, NULL) );

      for( int k = 0; k < nnodes; k++ )
         if( nodestatenew[k] == BRANCH_STP_VERTEX_TERM && !Is_term(g->term[k]) )
            graph_knot_chg(g, k, 0);
   }

   SCIP_CALL( sep_2cut(scip, conshdlr, conshdlrdata, consdata, termorg, maxcuts, &ncuts) );

   if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   /* restore graph */
   if( !atrootnode && g->stp_type == STP_SPG )
   {
      for( int k = 0; k < g->knots; k++ )
         if( g->term[k] != termorg[k] )
            graph_knot_chg(g, k, termorg[k]);
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
   int i;

   for( i = 0; i < nconss; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);

      SCIP_CALL( SCIPStpValidateSol(scip, consdata->graph, SCIPprobdataGetXval(scip, NULL), &feasible) );

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

      SCIP_CALL( SCIPStpValidateSol(scip, consdata->graph, SCIPprobdataGetXval(scip, NULL), &feasible) );

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

   SCIP_CALL(SCIPStpValidateSol(scip, g, SCIPprobdataGetXval(scip, sol), &feasible));

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
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         &conshdlrdata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/maxroundsroot",
         "maximal number of separation rounds per node in the root node (-1: unlimited)",
         &conshdlrdata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/maxsepacuts",
         "maximal number of cuts separated per separation round",
         &conshdlrdata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/maxsepacutsroot",
         "maximal number of cuts separated per separation round in the root node",
         &conshdlrdata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );


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

/** sets graph */
void SCIPStpConshdlrSetGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g                   /**< graph data structure */
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

/* dual ascent heuristic */
SCIP_RETCODE SCIPStpDualAscent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real* RESTRICT   redcost,            /**< array to store reduced costs or NULL */
   SCIP_Real* RESTRICT   nodearrreal,        /**< real vertices array for internal computations or NULL */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Bool             addcuts,            /**< should dual ascent add Steiner cuts? */
   SCIP_Bool             ascendandprune,     /**< should the ascent-and-prune heuristic be executed? */
   GNODE**               gnodearrterms,      /**< gnode terminals array for internal computations or NULL */
   const int*            result,             /**< solution array or NULL */
   int* RESTRICT         edgearrint,         /**< int edges array for internal computations or NULL */
   int* RESTRICT         nodearrint,         /**< int vertices array for internal computations or NULL */
   int                   root,               /**< the root */
   SCIP_Bool             is_pseudoroot,      /**< is the root a pseudo root? */
   SCIP_Real             damaxdeviation,     /**< maximum deviation for dual-ascent ( -1.0 for default) */
   STP_Bool* RESTRICT    nodearrchar         /**< STP_Bool vertices array for internal computations or NULL */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_PQUEUE* pqueue;
   SCIP_VAR** vars;
   SCIP_Real* RESTRICT rescap;
   GNODE** gnodearr;
   int* RESTRICT edgearr;
   int* RESTRICT tailarr;
   int* RESTRICT start;
   int* RESTRICT stackarr;
   int* RESTRICT cutverts;
   int* RESTRICT unsatarcs;
   int* RESTRICT unsattails;
   int* RESTRICT gmark;
   int* RESTRICT active;
   SCIP_Real dualobj;
   SCIP_Real currscore;
   const SCIP_Real maxdeviation = (damaxdeviation > 0.0) ? damaxdeviation : DEFAULT_DAMAXDEVIATION;
   const int nnodes = g->knots;
   const int nterms = g->terms;
   const int nedges = g->edges;
   int ncsredges;
   int norgcutverts;
   int stacklength;
   int augmentingcomponent;
   const SCIP_Bool addconss = (SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE);

   /* should currently not  be activated */
   assert(addconss || !addcuts);
   assert(g != NULL);
   assert(scip != NULL);
   assert(objval != NULL);
   assert(Is_term(g->term[root]));
   assert(maxdeviation >= DA_MAXDEVIATION_LOWER && maxdeviation <= DA_MAXDEVIATION_UPPER);
   assert(damaxdeviation == -1.0 || damaxdeviation > 0.0);

   if( nnodes == 1 )
      return SCIP_OKAY;

   if( addcuts )
   {
      vars = SCIPprobdataGetVars(scip);
      assert(vars != NULL);

      if( !addconss )
      {
         conshdlr = SCIPfindConshdlr(scip, "stp");
         assert(conshdlr != NULL);
      }
   }
   else
   {
      vars = NULL;
   }

   /* if specified root is not a terminal, take default root */
   if( !Is_term(g->term[root]) )
      root = g->source;

#ifdef BITFIELDSARRAY
   u_int32_t* bitarr;
   SCIP_CALL( SCIPallocBufferArray(scip, &bitarr, nedges / ARRLENGTH + 1) );
#endif

   stacklength = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &unsattails, nedges) );

   if( redcost == NULL )
      SCIP_CALL( SCIPallocBufferArray(scip, &rescap, nedges) );
   else
      rescap = redcost;

   if( nodearrint == NULL )
      SCIP_CALL( SCIPallocBufferArray(scip, &cutverts, nnodes) );
   else
      cutverts = nodearrint;

   if( edgearrint == NULL )
      SCIP_CALL( SCIPallocBufferArray(scip, &unsatarcs, nedges) );
   else
      unsatarcs = edgearrint;

   if( gnodearrterms == NULL )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
      for( int i = 0; i < nterms - 1; i++ )
         SCIP_CALL( SCIPallocBlockMemory(scip, &gnodearr[i]) ); /*lint !e866*/
   }
   else
   {
      gnodearr = gnodearrterms;
   }

   SCIP_CALL( SCIPpqueueCreate(&pqueue, nterms, 2.0, GNODECmpByDist) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &active, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &edgearr, nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &tailarr, nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &start, nnodes + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &gmark, nnodes + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &stackarr, nnodes) );

   /* fill auxiliary adjacent vertex/edges arrays */
   graph_get_csr(g, edgearr, tailarr, start, &ncsredges);

   /* initialize priority queue and res. capacity */
   SCIP_CALL( dualascent_init(scip, g, start, edgearr, root, is_pseudoroot, ncsredges, gmark, active, pqueue,
         gnodearr, rescap, &dualobj, &augmentingcomponent) );

   /* mark whether an arc is satisfied (has capacity 0) */
   for( int i = 0; i < ncsredges; i++ )
   {
#ifdef BITFIELDSARRAY
      if( SCIPisZero(scip, rescap[i]) )
         SetBit(bitarr, i);
      else
         CleanBit(bitarr, i);
#else
      if( rescap[i] == 0.0 )
      {
         if( active[tailarr[i] - 1] == 0 )
            tailarr[i] = 0;
         else
            tailarr[i] *= -1;
      }
#endif
   }

   norgcutverts = 0;

   /* (main) dual ascent loop */
   while( SCIPpqueueNElems(pqueue) > 0 && !SCIPisStopped(scip) )
   {
      /* get active vertex of minimum score */
      GNODE* const gnodeact = (GNODE*) SCIPpqueueRemove(pqueue);
      const SCIP_Real prio1 = gnodeact->dist;
      const SCIP_Real prio2 = (SCIPpqueueNElems(pqueue) > 0) ? ((GNODE*) SCIPpqueueFirst(pqueue))->dist : FARAWAY;
      const int v = gnodeact->number;
      SCIP_Real degsum = g->grad[v];
      int ncutverts = 0;
      int nunsatarcs = 0;

      SCIP_Bool firstrun = TRUE;

      SCIPdebugMessage("DA: START WITH v %d prio1 %f prio2 %f \n", v, prio1, prio2);

      /* perform augmentation as long as priority of root component does not exceed max deviation */
      for( ; ; )
      {
         assert(stacklength == 0);

         /* 1. step: BFS from v (or connected component) on saturated, incoming arcs */

         if( firstrun )
         {
            firstrun = FALSE;
            gmark[v + 1] = TRUE;
            cutverts[ncutverts++] = v;
            assert(stacklength < nnodes);
            stackarr[stacklength++] = v;
         }
         /* not in first processing of root component: */
         else
         {
            for( int i = norgcutverts; i < ncutverts; i++ )
            {
               const int s = cutverts[i];

               assert(gmark[s + 1]);
               assert(active[s] != 0);
               assert(stacklength < nnodes);

               stackarr[stacklength++] = s;
            }
         }
#ifdef DFS
         while( stacklength )
         {
            const int node = stackarr[--stacklength];
#else
         for( int n = 0; n < stacklength; n++ )
         {
            int end;

            assert(n < nnodes);
            node = stackarr[n];
#endif

            /* traverse incoming arcs */
            for( int i = start[node], end = start[node + 1]; i != end; i++ )
            {
               int tail = tailarr[i];

               /* zero reduced-cost arc? */
               if( tail <= 0 )
               {
                  tail *= -1;
                  if( !gmark[tail] )
                  {
                     /* if an active vertex has been hit (other than v), break */
                     if( 0 == tail )
                     {
                        const int realtail = g->tail[edgearr[i]];

                        /* v should not be processed */
                        if( realtail == v )
                           continue;

                        /* is realtail active or does realtail lead to an active vertex other than v? */
                        if( is_active(active, realtail, v) )
                        {
                           active[v] = realtail + 1;
                           stacklength = 0;
                           goto ENDOFLOOP;
                        }

                        tail = realtail + 1;

                        /* have we processed tail already? */
                        if( gmark[tail] )
                           continue;
                     }

                     assert(tail > 0);

                     gmark[tail] = TRUE;
                     tail--;
                     cutverts[ncutverts++] = tail;
                     degsum += g->grad[tail];

                     assert(stacklength < nnodes);
                     stackarr[stacklength++] = tail;
                  } /* marked */
               } /* zero reduced-cost arc */
               else if( !gmark[tail] )
               {
                  unsattails[nunsatarcs] = tail;
                  unsatarcs[nunsatarcs++] = i;
               }
            }
         }
#ifndef DFS
         stacklength = 0;
#endif
         currscore = degsum - (ncutverts - 1);

         /* guiding solution provided? */
         if( result != NULL )
         {
            int nsolarcs = 0;
            for( int i = 0; i < nunsatarcs; i++ )
            {
               const int a = unsatarcs[i];

               assert(tailarr[a] > 0);

               if( !(gmark[tailarr[a]]) )
               {
                  if( result[edgearr[a]] == CONNECT )
                     nsolarcs++;
               }
            }

            assert(nsolarcs > 0);
            assert(currscore <= nedges);

            if( nsolarcs > 1 )
              currscore += (SCIP_Real) ((nsolarcs - 1) * (g->knots * 2.0));
         }
         else
         {
            assert(SCIPisGE(scip, currscore, prio1));
         }

         SCIPdebugMessage("DA: deviation %f \n", (currscore - prio1) / prio1);
         SCIPdebugMessage("DA: currscore %f prio1 %f prio2 %f \n", currscore, prio1, prio2);

         /* augmentation criteria met? */
         if( ((currscore - prio1) / prio1) <= maxdeviation || currscore <= prio2 )
         {
            SCIP_CONS* cons = NULL;
            SCIP_ROW* row = NULL;

            int shift = 0;
            SCIP_Real min = FARAWAY;
            SCIP_Bool isactive = FALSE;

            /* 2. step: get minimum residual capacity among cut-arcs */

            /* adjust array of unsatisfied arcs */

            for( int i = 0; i < nunsatarcs; i++ )
            {
               const int tail = unsattails[i];

               if( gmark[tail] )
               {
                  shift++;
               }
               else
               {
                  const int a = unsatarcs[i];

                  assert(tailarr[a] > 0);
                  assert(rescap[a] > 0);

                  if( rescap[a] < min )
                     min = rescap[a];
                  if( shift )
                  {
                     unsattails[i - shift] = tail;
                     unsatarcs[i - shift] = a;
                  }
               }
            }

            assert(SCIPisLT(scip, min, FARAWAY));
            nunsatarcs -= shift;

            norgcutverts = ncutverts;

            /* 3. step: perform augmentation */

            /* create constraints/cuts ? */
            if( addcuts )
            {
               if( addconss )
               {
                  SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "da", 0, NULL, NULL,
                        1.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
               }
               else
               {
                  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "da", 1.0,
                        SCIPinfinity(scip), FALSE, FALSE, TRUE) );

                  SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
               }
            }

            shift = 0;

            /* update (dual) objective */
            dualobj += min;

            for( int i = 0; i < nunsatarcs; i++ )
            {
               const int a = unsatarcs[i];
               assert(a >= 0);

               if( addcuts )
               {
                  assert(vars != NULL);

                  if( addconss )
                     SCIP_CALL( SCIPaddCoefLinear(scip, cons, vars[edgearr[a]], 1.0) );
                  else
                     SCIP_CALL( SCIPaddVarToRow(scip, row, vars[edgearr[a]], 1.0) );
               }
               rescap[a] -= min;

               assert(SCIPisGE(scip, rescap[a], 0.0));

               if( rescap[a] <= DA_EPS )
               {
                  int tail = unsattails[i];

                  rescap[a] = 0.0;

                  assert(tail > 0);
                  assert(tailarr[a] > 0);

                  tailarr[a] *= -1;

                  if( active[tail - 1] >= 0 && is_active(active, tail - 1, v) )
                  {
                     assert(tail - 1 != v);
                     tailarr[a] = 0;
                     if( !isactive )
                     {
                        isactive = TRUE;
                        active[v] = tail;
                     }
                  }


                  if( !(gmark[tail])  )
                  {
                     assert(tail != 0);

                     gmark[tail] = TRUE;
                     tail--;
                     degsum += g->grad[tail];
                     cutverts[ncutverts++] = tail;
                  }

                  shift++;
               }
               else if( shift )
               {
                  unsattails[i - shift] = unsattails[i];
                  unsatarcs[i - shift] = a;
               }
            }

            if( addcuts )
            {
               if( addconss )
               {
                  SCIP_CALL( SCIPaddCons(scip, cons) );
                  SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               }
               else
               {
                  SCIP_Bool infeasible;

                  SCIP_CALL( SCIPflushRowExtensions(scip, row) );
                  SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
                  SCIP_CALL( SCIPreleaseRow(scip, &row) );

                  assert(!infeasible);
               }
            }

            if( isactive )
            {
               stacklength = 0;
               goto ENDOFLOOP;
            }
            nunsatarcs -= shift;

            continue;
         }
         else
         {
            SCIP_Bool insert = TRUE;

            if( is_pseudoroot )
            {
               int i = start[v];
               const int end = start[v + 1];

               assert(end - i == 2);

               for( ; i != end; i++ )
                  if( rescap[i] != 0.0 )
                     break;

               if( i == end )
               {
                  if( augmentingcomponent == -1 )
                     augmentingcomponent = v;

                  if( augmentingcomponent != v )
                     insert = FALSE;
               }
            }

            if( insert )
            {
               /* reinsert active vertex */
               gnodeact->dist = currscore;
               SCIP_CALL( SCIPpqueueInsert(pqueue, gnodeact) );
            }
         }

         ENDOFLOOP:

         for( int i = 0; i < ncutverts; i++ )
            gmark[cutverts[i] + 1] = FALSE;

         for( int i = 0; i < nnodes + 1; i++ )
         {
            assert(!gmark[i]);
         }

         break;
      } /* augmentation loop */
   } /* dual ascent loop */

   SCIPdebugMessage("DA: dualglobal: %f \n", dualobj);
   *objval = dualobj;

   for( int i = ncsredges; i < nedges; i++ )
   {
      edgearr[i] = i;
      rescap[i] = g->cost[i];
   }

   /* re-extend rescap array */
   for( int i = 0; i < ncsredges; i++ )
   {
      if( edgearr[i] != i  )
      {
         SCIP_Real bufferedval = rescap[i];
         int a = i;

         rescap[i] = g->cost[i];
         while( edgearr[a] != a )
         {
            const int shift = edgearr[a];
            const SCIP_Real min = rescap[shift];

            rescap[shift] = bufferedval;
            bufferedval = min;
            edgearr[a] = a;
            a = shift;
         }
      }
   }

#ifdef BITFIELDSARRAY
   SCIPfreeBufferArray(scip, &bitarr);
#endif

   SCIPfreeMemoryArray(scip, &stackarr);
   SCIPfreeMemoryArray(scip, &gmark);
   SCIPfreeMemoryArray(scip, &start);
   SCIPfreeMemoryArray(scip, &tailarr);
   SCIPfreeMemoryArray(scip, &edgearr);
   SCIPfreeMemoryArray(scip, &active);

   SCIPpqueueFree(&pqueue);

   if( gnodearrterms == NULL )
   {
      for( int i = nterms - 2; i >= 0; i-- )
         SCIPfreeBlockMemory(scip, &gnodearr[i]);
      SCIPfreeBufferArray(scip, &gnodearr);
   }

   /* call Ascend-And-Prune? */
   if( ascendandprune )
   {
       SCIP_Bool success;
       STP_Bool* RESTRICT mynodearrchar = NULL;

       if( nodearrchar == NULL )
          SCIP_CALL( SCIPallocBufferArray(scip, &mynodearrchar, nnodes) );
       else
          mynodearrchar = nodearrchar;

       SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, g, rescap, unsatarcs, cutverts, root, mynodearrchar, &success, TRUE) );

       if( nodearrchar == NULL )
          SCIPfreeBufferArray(scip, &mynodearrchar);
   }

   if( edgearrint == NULL )
      SCIPfreeBufferArray(scip, &unsatarcs);

   if( nodearrint == NULL )
      SCIPfreeBufferArray(scip, &cutverts);

   if( redcost == NULL )
      SCIPfreeBufferArray(scip, &rescap);

   SCIPfreeBufferArray(scip, &unsattails);

   return SCIP_OKAY;
}

/** dual ascent heuristic for PCSPG and MWCSP */
SCIP_RETCODE SCIPStpDualAscentPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            redcost,            /**< array to store reduced costs or NULL */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Bool             addcuts,            /**< should dual-ascent add Steiner cuts? */
   SCIP_Bool             ascendandprune,     /**< perform ascend-and-prune and add solution? */
   int                   nruns               /**< number of dual ascent runs */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_PQUEUE* pqueue;
   SCIP_VAR** vars;
   GRAPH* transgraph;
   SCIP_Real min;
   SCIP_Real prio1;
   SCIP_Real offset;
   SCIP_Real dualobj;
   SCIP_Real currscore;
   SCIP_Real maxdeviation;
   SCIP_Real* rescap;
   GNODE* gnodeact;
   GNODE** gnodearr;
   int s;
   int i;
   int k;
   int v;
   int a;
   int tail;
   int pnode;
   int shift;
   int root;
   int nnodes;
   int nterms;
   int nedges;
   int degsum;
   int ncutverts;
   int pseudoroot;
   int nunsatarcs;
   int stacklength;
   int norgcutverts;
   int* cutverts;
   int* stackarr;
   STP_Bool* origedge;
   int* unsatarcs;
   STP_Bool firstrun;
   STP_Bool* sat;
   STP_Bool* active;
   const SCIP_Bool addconss = (SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE);

   /* should currently not  be activated */
   assert(addconss || !addcuts);

   assert(g != NULL);
   assert(scip != NULL);
   assert(nruns >= 0);
   assert(objval != NULL);

   if( g->knots == 1 )
      return SCIP_OKAY;

   if( addcuts )
   {
      vars = SCIPprobdataGetVars(scip);
      assert(vars != NULL);
      if( !addconss )
      {
         conshdlr = SCIPfindConshdlr(scip, "stp");
         assert(conshdlr != NULL);
      }
   }
   else
   {
      vars = NULL;
   }

   root = g->source;
   degsum = 0;
   offset = 0.0;
   dualobj = 0.0;

   ncutverts = 0;
   norgcutverts = 0;
   maxdeviation = DEFAULT_DAMAXDEVIATION;

   SCIP_CALL( graph_pc_getSap(scip, g, &transgraph, &offset) );

   nnodes = transgraph->knots;
   nedges = transgraph->edges;
   nterms = transgraph->terms;
   pseudoroot = nnodes - 1;

   if( redcost == NULL )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &rescap, nedges) );
   }
   else
   {
      rescap = redcost;
   }

   stacklength = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, &stackarr, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sat, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &active, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutverts, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &unsatarcs, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &origedge, nedges) );

   for( i = 0; i < nedges; i++ )
      if( !Is_term(transgraph->term[transgraph->tail[i]]) && transgraph->head[i] == pseudoroot )
         origedge[i] = FALSE;
      else if( transgraph->tail[i] == pseudoroot && !Is_term(transgraph->term[transgraph->head[i]])  )
         origedge[i] = FALSE;
      else
         origedge[i] = TRUE;

   for( i = 0; i < nterms - 1; i++ )
   {
      SCIP_CALL( SCIPallocBuffer(scip, &gnodearr[i]) ); /*lint !e866*/
   }

   SCIP_CALL( SCIPpqueueCreate( &pqueue, nnodes, 2.0, GNODECmpByDist) );

   k = 0;
   /* mark terminals as active, add all except root to pqueue */
   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(transgraph->term[i]) )
      {
         active[i] = TRUE;
         assert(transgraph->grad[i] > 0);
         if( i != root  )
         {
            gnodearr[k]->number = i;
            gnodearr[k]->dist = transgraph->grad[i];

            for( a = transgraph->inpbeg[i]; a != EAT_LAST; a = transgraph->ieat[a] )
               if( SCIPisEQ(scip, transgraph->cost[a], 0.0) )
                  break;

            if( a != EAT_LAST )
               gnodearr[k]->dist += transgraph->grad[transgraph->tail[a]] - 1;

            assert(gnodearr[k]->dist > 0);

            SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[k++]) );
         }
      }
      else
      {
         active[i] = FALSE;
      }
      transgraph->mark[i] = FALSE;
   }

   for( i = 0; i < nedges; i++ )
   {
      rescap[i] = transgraph->cost[i];
      if( SCIPisZero(scip, rescap[i]) )
         sat[i] = TRUE;
      else
         sat[i] = FALSE;
   }

   /* dual ascent loop */
   while( SCIPpqueueNElems(pqueue) > 0 && !SCIPisStopped(scip) )
   {
      /* get active vertex of minimum score */
      gnodeact = (GNODE*) SCIPpqueueRemove(pqueue);

      v = gnodeact->number;
      prio1 = gnodeact->dist;

      firstrun = TRUE;
      nunsatarcs = 0;

      /* perform augmentation as long as ... */
      for( ; ; )
      {
         assert(stacklength == 0);
         /* 1. step: BFS from v (or connected component) on saturated, incoming arcs */

         if( firstrun )
         {
            degsum = transgraph->grad[v];
            ncutverts = 0;
            firstrun = FALSE;
            nunsatarcs = 0;
            transgraph->mark[v] = TRUE;
            cutverts[ncutverts++] = v;
            stackarr[stacklength++] = v;
         }
         /* not in first processing of root component: */
         else
         {
            for( i = norgcutverts; i < ncutverts; i++ )
            {
               s = cutverts[i];
               assert(transgraph->mark[s]);
               if( active[s] )
               {
                  active[v] = FALSE;
                  stacklength = 0;
                  goto ENDOFLOOP;
               }

               stackarr[stacklength++] = s;
            }
         }

         while( stacklength )
         {
            pnode = stackarr[--stacklength];

            /* traverse incoming arcs */
            for( a = transgraph->inpbeg[pnode]; a != EAT_LAST; a = transgraph->ieat[a] )
            {
               tail = transgraph->tail[a];
               if( sat[a] )
               {
                  if( !transgraph->mark[tail] )
                  {
                     /* if an active vertex has been hit, break */
                     if( active[tail] )
                     {
                        active[v] = FALSE;
                        stacklength = 0;
                        goto ENDOFLOOP;
                     }

                     degsum += transgraph->grad[tail];
                     transgraph->mark[tail] = TRUE;
                     cutverts[ncutverts++] = tail;
                     stackarr[stacklength++] = tail;
                  }
               }
               else if( !transgraph->mark[tail] )
               {
                  unsatarcs[nunsatarcs++] = a;
               }
            }
         }

         currscore = degsum - (ncutverts - 1);

         assert(SCIPisGE(scip, currscore, prio1));

         /* augmentation criteria met? */
         if( SCIPisLE(scip, (currscore - prio1) / prio1, maxdeviation) || (SCIPpqueueNElems(pqueue) == 0) )
         {
            SCIP_Bool in = FALSE;
            SCIP_ROW* row;
            SCIP_CONS* cons = NULL;

            /* 2. pass: get minimum residual capacity among cut-arcs */

            /* adjust array of unsatisfied arcs */
            min = FARAWAY;
            shift = 0;

            for( i = 0; i < nunsatarcs; i++ )
            {
               a = unsatarcs[i];
               if( transgraph->mark[transgraph->tail[a]] )
               {
                  shift++;
               }
               else
               {

                  assert(!sat[a]);
                  if( SCIPisLT(scip, rescap[a], min) )
                     min = rescap[a];
                  if( shift != 0 )
                     unsatarcs[i - shift] = a;
               }
            }

            assert(SCIPisLT(scip, min, FARAWAY));
            nunsatarcs -= shift;

            if( nunsatarcs > 0)
               assert(!transgraph->mark[transgraph->tail[unsatarcs[nunsatarcs-1]]]);

            norgcutverts = ncutverts;


            /* 3. pass: perform augmentation */


            /* create constraint/row */

            if( addcuts )
            {
               if( addconss )
               {
                  SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "da", 0, NULL, NULL,
                        1.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );
               }
               else
               {
                  SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "da", 1.0, SCIPinfinity(scip), FALSE, FALSE, TRUE));
                  SCIP_CALL(SCIPcacheRowExtensions(scip, row));
               }
            }

            dualobj += min;
            for( i = 0; i < nunsatarcs; i++ )
            {
               a = unsatarcs[i];
               if( a == -1 )
                  continue;

               if( addcuts && origedge[a] )
               {
                  assert(vars != NULL);
                  assert(cons != NULL);

                  if( g->tail[a] == root && g->head[a] == v )
                     in = TRUE;

                  if( addconss )
                     SCIP_CALL( SCIPaddCoefLinear(scip, cons, vars[a], 1.0) );
                  else
                     SCIP_CALL( SCIPaddVarToRow(scip, row, vars[a], 1.0) );
               }
               rescap[a] -= min;

               assert(SCIPisGE(scip, rescap[a], 0.0));

               if( SCIPisEQ(scip, rescap[a], 0.0) )
               {
                  sat[a] = TRUE;
                  if( !(transgraph->mark[transgraph->tail[a]]) )
                  {
                     tail = transgraph->tail[a];
                     degsum += transgraph->grad[tail];
                     transgraph->mark[tail] = TRUE;
                     cutverts[ncutverts++] = tail;
                  }
               }
            }

            if( addcuts )
            {
               assert(vars != NULL);

               if( !in )
               {
                  for( i = g->outbeg[root]; i != EAT_LAST; i = g->oeat[i] )
                     if( g->head[i] == v )
                     {
                        if( addconss )
                           SCIP_CALL( SCIPaddCoefLinear(scip, cons, vars[i], 1.0) );
                        else
                           SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i], 1.0) );
                     }
               }

               if( addconss )
               {
                  SCIP_CALL( SCIPaddCons(scip, cons) );
                  SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               }
               else
               {
                  SCIP_Bool infeasible;
                  assert(row != NULL);

                  SCIP_CALL( SCIPflushRowExtensions(scip, row) );
                  SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
                  SCIP_CALL( SCIPreleaseRow(scip, &row) );

                  assert(!infeasible);
               }
            }

            continue;
         }
         else
         {
            /* reinsert active vertex */
            gnodeact->dist = currscore;
            SCIP_CALL( SCIPpqueueInsert(pqueue, gnodeact) );
         }

      ENDOFLOOP:

         for( i = 0; i < ncutverts; i++ )
            transgraph->mark[cutverts[i]] = FALSE;

         break;
      } /* augmentation loop */
   } /* dual ascent loop */


   *objval = dualobj + offset;
   SCIPdebugMessage("DA: dualglobal: %f \n", *objval + SCIPprobdataGetOffset(scip));

   /* call dual Ascend-And-Prune? */
   if( ascendandprune )
   {
      SCIP_Bool success;
      SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, g, rescap, unsatarcs, cutverts, -1, active, &success, TRUE));
   }

   /* free memory */
   SCIPpqueueFree(&pqueue);

   for( i = nterms - 2; i >= 0; i-- )
      SCIPfreeBuffer(scip, &gnodearr[i]);

   SCIPfreeBufferArray(scip, &origedge);
   SCIPfreeBufferArray(scip, &unsatarcs);
   SCIPfreeBufferArray(scip, &cutverts);
   SCIPfreeBufferArray(scip, &gnodearr);
   SCIPfreeBufferArray(scip, &active);
   SCIPfreeBufferArray(scip, &sat);
   SCIPfreeBufferArray(scip, &stackarr);

   if( redcost == NULL )
      SCIPfreeBufferArray(scip, &rescap);

   graph_free(scip, &transgraph, TRUE);

   return SCIP_OKAY;

}

/**@} */
