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
 * This file implements a constraint handler checking solutions for feasibility and separating violated model constraints as
 * described in: "Solving the Steiner Tree Problem in Graphs to Optimality" by T. Koch and A. Martin.
 * The separation problem for the cut inequalities described in \ref PROBLEM can be solved by a max flow algorithm in
 * polynomial time.  Regarding the variable values of a given LP solution as capacities on the edges, one can check for each
 * \f$t ∈ T \ {r}\f$, with \f$r\f$ being the root, whether the minimal \f$(r, t)\f$-cut is less than one. In this case,
 * a violated cut inequality has been found, otherwise none exists. In order to calculate such a minimal cut, Hao
 * and Orlin's preflow-push algorithm is used.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "cons_stp.h"
#include "probdata_stp.h"

#include "scip/scip.h"
#include "grph.h"
#include "portab.h"

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

#define DEFAULT_MAXROUNDS             5 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS     INT_MAX /**< maximal number of cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT INT_MAX /**< maximal number of cuts separated per separation round in the root node */


#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP

#define DEFAULT_BACKCUT       FALSE /**< Try Back-Cuts */
#define DEFAULT_CREEPFLOW      TRUE /**< Use Creep-Flow */
#define DEFAULT_DISJUNCTCUT   FALSE /**< Only disjunct Cuts */
#define DEFAULT_NESTEDCUT     FALSE /**< Try Nested-Cuts */
#define DEFAULT_FLOWSEP        TRUE /**< Try Flow-Cuts */

#define FLOW_FACTOR     100000
#define CREEP_VALUE     1
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
   SCIP_Bool             disjunctcut;        /**< should disjunktion cuts be applied? */
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


/** add a cut */
static
SCIP_RETCODE cut_add(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   const GRAPH*          g,                  /**< graph data structure */
   const int             layer,              /**< current layer, set to zero usually */
   const SCIP_Real*      xval,               /**< edge values */
   int*                  capa,               /**< edges capacities (scaled) */
   const int             updatecapa,         /**< update capacities? */
   int*                  ncuts,              /**< pointer to store number of cuts */
   SCIP_Bool*            success             /**< pointer to store whether add cut be added */
   )
{
   SCIP_ROW* row;
   SCIP_VAR** vars = SCIPprobdataGetVars(scip);
   SCIP_Real sum = 0.0;
   SCIP_Bool inccapa = FALSE;
   int i;
   int ind;
   (*success) = FALSE;

   assert(scip != NULL);
   assert(g         != NULL);
   assert((layer >= 0) && (layer < g->layers));

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "2cut", 1.0, SCIPinfinity(scip), FALSE, FALSE, TRUE) );

   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for( i = 0; i < g->edges; i++ )
   {
      if( (g->mark[g->source[layer]] == g->mark[g->tail[i]])
         && (g->mark[g->tail[i]] != g->mark[g->head[i]]) )
      {
         ind = layer * g->edges + i;

         if( updatecapa )
         {
            if( capa[i] < FLOW_FACTOR )
               inccapa = TRUE;

            SCIPdebugMessage("set capa[%d] from %6d to %6d\n", i, capa[i], FLOW_FACTOR);
            capa[i] = FLOW_FACTOR;

            if( !inccapa )
            {
               SCIP_CALL( SCIPflushRowExtensions(scip, row) );
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
               return SCIP_OKAY;
            }
         }

         if( xval != NULL )
         {
            sum += xval[ind];

            if( SCIPisFeasGE(scip, sum, 1.0) )
            {
               SCIP_CALL( SCIPflushRowExtensions(scip, row) );
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
               return SCIP_OKAY;
            }
         }
         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], 1.0) );
      }
   }
   assert(sum < 1.0);

   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   /* checks, if cut is sufficiently violated */
   if( SCIPisCutEfficacious(scip, NULL, row) )
   {
      SCIP_Bool infeasible;

      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

      SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &infeasible) );
      (*ncuts)++;
      (*success) = TRUE;
   }

   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}

static
int graph_next_term(
   int           terms,
   int*          term,
   const int*    w
   )
{
   int wmax = 0;
   int i;
   int k = -1;
   int t;

   assert(term != NULL);

   for( i = 0; (i < terms); i++ )
   {
      if( w[term[i]] == 0 )
      {
         k = i;
         break;
      }
      if( w[term[i]] > wmax )
      {
         k    = i;
         wmax = w[term[i]];
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
   const int             layer,              /**< current layer, usually set to zero  */
   const SCIP_Bool       creep_flow,         /**< creep flow? */
   const int             flip,               /**< reverse the flow? */
   int*                  capa,               /**< edges capacities (scaled) */
   const SCIP_Real*      xval                /**< edge values */
   )
{
   int k;

   assert(g    != NULL);
   assert(xval != NULL);

   for( k = 0; k < g->edges; k += 2 )
   {
      if( !flip )
      {
         capa[k]     = (int)(xval[layer * g->edges + k    ]
            * FLOW_FACTOR + 0.5);
         capa[k + 1] = (int)(xval[layer * g->edges + k + 1]
            * FLOW_FACTOR + 0.5);
      }
      else
      {
         capa[k]     = (int)(xval[layer * g->edges + k + 1]
            * FLOW_FACTOR + 0.5);
         capa[k + 1] = (int)(xval[layer * g->edges + k    ]
            * FLOW_FACTOR + 0.5);
      }

      if( creep_flow && (capa[k] == 0) && (capa[k + 1] == 0) )
      {
         capa[k]     = CREEP_VALUE;
         capa[k + 1] = CREEP_VALUE;
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
         if( i == g->source[layer] )
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

               SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &infeasible) );
               count++;

               SCIP_CALL( SCIPreleaseRow(scip, &row) );

               if( *ncuts + count >= maxcuts )
                  goto TERMINATE;
            }
         }
         /* no flows ? */
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

               SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &infeasible) );
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

            SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &infeasible) );
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
                  (g->locals[layer] == 2) ? 0.0 : SCIPinfinity(scip), FALSE, FALSE, TRUE) );

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

            SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &infeasible) );
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
   int                   maxcuts,            /**< maximal number of cuts */
   int*                  ncuts               /**< pointer to store number of cuts */
   )
{
   const SCIP_Bool nested_cut   = conshdlrdata->nestedcut;
   const SCIP_Bool back_cut     = conshdlrdata->backcut;
   const SCIP_Bool creep_flow   = conshdlrdata->creepflow;
   const SCIP_Bool disjunct_cut = conshdlrdata->disjunctcut;
   const SCIP_Bool flowsep      = conshdlrdata->flowsep;

   GRAPH*  g;
   SCIP_Real* xval;
   SCIP_Real* cost;
   PATH*   path;
   int*    w;
   int*    capa;
   int*    term;
   int     terms = 0;
   int     tsave;
   int     i;
   int     k;
   int     layer;
   int     count = 0;
   int     rerun = FALSE;
   int     nedges;
   int     nnodes;
   SCIP_Bool addedcut;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conshdlrdata != NULL);

   g = consdata->graph;
   assert(g != NULL);

   nedges = g->edges;
   nnodes = g->knots;
   addedcut = FALSE;

   xval = SCIPprobdataGetXval(scip, NULL);
   assert(xval != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &capa, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &w, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &term, g->terms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );

   for( layer = 0; layer < g->layers; layer++ )
   {
      /* For 2-terminal nets no cuts are necessary if flows are given */
      if( flowsep && (g->locals[layer] == 2) )
         continue;

      for( i = 0; i < nedges; i++ )
         cost[i] = SCIPisFeasLT(scip, xval[layer * nedges + i], 1.0) ? 1.0 : 0.0;

      for( i = 0; i < nnodes; i++ )
      {
	 w[i] = 0;
         g->mark[i] = TRUE;
      }

      graph_path_exec(scip, g, FSP_MODE, g->source[layer], cost, path);

      /* search all terminals not connected to the root by the LP solution */
      for( i = 0, count = 0; i < nnodes; i++ )
      {
         if( (g->term[i] == layer) && (i != g->source[layer]) )
         {
            if( SCIPisPositive(scip, path[i].dist) )
               term[terms++] = i;
            else
               count++;
         }
      }
      SCIPdebugMessage("Cut Pretest: %d eliminations\n", count);

      count = 0;
      tsave = terms;

      /* from source to terminal */
      if( !nested_cut || disjunct_cut )
         set_capacity(g, layer, creep_flow, 0, capa, xval);

      while( terms > 0 )
      {
         if( SCIPisStopped(scip) && terms % 100 == 0 )
            break;

         /* look for reachable terminal */
         i = graph_next_term(terms, term, w);

         terms--;

         assert(g->term[i]       == layer);
         assert(g->source[layer] != i);

         if( nested_cut && !disjunct_cut )
            set_capacity(g, layer, creep_flow, 0, capa, xval);

         do
         {
            graph_mincut_exec(g, g->source[layer], i, capa, w, rerun);

            rerun = TRUE;

            /* cut */
            for( k = 0; k < nnodes; k++ )
               g->mark[k] = (w[k] != 0);

	    SCIP_CALL( cut_add(scip, conshdlr, g, layer, xval, capa, nested_cut || disjunct_cut, ncuts, &addedcut) );
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
      }

      /* back cuts enabled? */
      if( back_cut )
      {
         if( !nested_cut || disjunct_cut )
            set_capacity(g, layer, creep_flow, 1, capa, xval);

         terms = tsave;

         while( terms > 0 )
         {
            /* look for reachable terminal */
            i = graph_next_term(terms, term, w);

            terms--;

            assert(g->term[i]       == layer);
            assert(g->source[layer] != i);

            if( nested_cut && !disjunct_cut )
               set_capacity(g, layer, creep_flow, 1, capa, xval);

            rerun = FALSE;

            do
            {
               graph_mincut_exec(g, i, g->source[layer], capa, w, rerun);

               rerun = TRUE;

               for( k = 0; k < nnodes; k++ )
                  g->mark[k] = (w[k] != 0) ? 1 : 0;

	       SCIP_CALL( cut_add(scip, conshdlr, g, layer, xval, capa, nested_cut || disjunct_cut, ncuts, &addedcut) );
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
            while( nested_cut );                /* Nested Cut is CONSTANT ! */

            rerun = FALSE;
         }
      }
   }

 TERMINATE:
   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &term);
   SCIPfreeBufferArray(scip, &w);
   SCIPfreeBufferArray(scip, &cost);
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

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpStp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int maxcuts;
   int ncuts = 0;
   int i;

   *result = SCIP_DIDNOTRUN;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   maxcuts = SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0 ?
      conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts;

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[i]);

      SCIP_CALL( sep_flow(scip, conshdlr, conshdlrdata, consdata, maxcuts, &ncuts) );

      SCIP_CALL( sep_2cut(scip, conshdlr, conshdlrdata, consdata, maxcuts, &ncuts) );
   }

   if( ncuts > 0 )
      *result = SCIP_SEPARATED;

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

      SCIP_CALL( SCIPvalidateStpSol(scip, consdata->graph, SCIPprobdataGetXval(scip, NULL), &feasible) );

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
   SCIP_CONSDATA* consdata;
   int i;

   for( i = 0; i < nconss; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);

      SCIP_CALL( SCIPvalidateStpSol(scip, consdata->graph, SCIPprobdataGetXval(scip, NULL), &feasible) );

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
{  /*lint --e{715}*/
   SCIP_Bool feasible;
   SCIP_CONSDATA* consdata;
   int i;

   for( i = 0; i < nconss; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);

      SCIP_CALL( SCIPvalidateStpSol(scip, consdata->graph, SCIPprobdataGetXval(scip, sol), &feasible) );

      if( !feasible )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
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
   if( graph->stp_type == STP_DEG_CONS )
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
      SCIP_CALL( SCIPaddVarLocks(scip, vars[v], 1, 1) );

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

/**@} */
