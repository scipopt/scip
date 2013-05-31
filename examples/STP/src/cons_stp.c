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

/**@file   cons_stp.c
 * @brief  Constraint handler stores the local branching decision data
 * @author Timo Berthold
 * @author Stefan Heinz
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
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ            1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS         0 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP

#define DEFAULT_BACKCUT       FALSE /**< Try Back-Cuts */
#define DEFAULT_CREEPFLOW      TRUE /**< Use Creep-Flow */
#define DEFAULT_DISJUNCTCUT   FALSE /**< Only disjunct Cuts */
#define DEFAULT_NESTEDCUT     FALSE /**< Try Nested-Cuts */
#define DEFAULT_FLOWSEP        TRUE /**< Try Flow-Cuts */

#define FLOW_FACTOR     100000
#define CREEP_VALUE     1

#define C_T_CAPA     0
#define C_T_2CUT     1
#define C_T_TERM     2
#define C_T_INFL     3
#define C_T_BALA     4
#define C_T_FLOW     5
#define C_T_WAVE     6
#define C_T_MSTR     7
#define C_T_PLAN     8
#define C_T_HOPS     9
#define C_T_LAST    10  /* Ist immer der Letzte */


/**@} */

/*
 * Data structures
 */

/** @brief Constraint data for  \ref cons_stp.c "Stp" constraints */
struct SCIP_ConsData
{
   GRAPH*                graph;
};

/** @brief Constraint handler data for \ref cons_stp.c "Stp" constraint handler */
struct SCIP_ConshdlrData
{
   SCIP_Bool backcut;
   SCIP_Bool creepflow;
   SCIP_Bool disjunctcut;
   SCIP_Bool nestedcut;
   SCIP_Bool flowsep;
};


/**@name Local methods
 *
 * @{
 */


static int cut_add(
   SCIP*         scip,
   SCIP_CONSHDLR* conshdlr,
   const GRAPH*  g,
   const int     layer,
   const int     type,
   const double* xval,
   int* ncuts
   )
{
   int     i;
   int     ind;
   int ret = 0;
   double  sum   = 0.0;
   SCIP_ROW* row;
   SCIP_VAR** vars = SCIPprobdataGetVars(scip);

   assert(scip != NULL);
   assert(g         != NULL);
   assert(layer     >= 0);

   assert((layer >= 0) && (layer < g->layers));

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "2cut", 1.0, SCIPinfinity(scip), FALSE, FALSE, TRUE) );

   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for(i = 0; i < g->edges; i++)
   {
      if ((g->mark[g->source[layer]] == g->mark[g->tail[i]])
       && (g->mark[g->tail[i]] != g->mark[g->head[i]]))
      {
         ind = layer * g->edges + i;

         if (xval != NULL)
         {
            sum += xval[ind];

            if( SCIPisFeasGE(scip, sum, 1.0) )
            {
               SCIP_CALL( SCIPflushRowExtensions(scip, row) );
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
               return 0;
            }
         }
         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], 1.0) );
      }
   }
   assert(sum < 1.0);

   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   /* checks, if cut is violated enough */
   if( SCIPisCutEfficacious(scip, NULL, row) )
   {
      SCIP_Bool infeasible;

      SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &infeasible) );
      (*ncuts)++;
      ret = 1;
   }

   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return ret;
}

static int graph_next_term(
   int           terms,
   int*          term,
   const int*    w)
{
   int wmax = 1;
   int i;
   int k    = -1;
   int t;

   assert(term != NULL);

   for(i = 0; (i < terms); i++)
   {
      if (w[term[i]] == 0)
      {
         k = i;
         break;
      }
      if (w[term[i]] > wmax)
      {
         k    = i;
         wmax = w[term[i]];
      }
   }
   assert(k >= 0);
   assert(k < terms);

   t       = term[k];
   term[k] = term[terms - 1];

   return(t);
}

static void set_capacity(
   const GRAPH*  g,
   const int     layer,
   const int     creep_flow,
   const int     flip,
   int*          capa,
   const double* xval)
{
   int k;

   assert(g    != NULL);
   assert(xval != NULL);

   for(k = 0; k < g->edges; k += 2)
   {
      if (!flip)
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

      if (creep_flow && (capa[k] == 0) && (capa[k + 1] == 0))
      {
         capa[k]     = CREEP_VALUE;
         capa[k + 1] = CREEP_VALUE;
      }
   }
}

static
int sep_flow(
   SCIP* scip,
   SCIP_CONSHDLR* conshdlr,
   SCIP_CONSHDLRDATA* conshdlrdata,
   SCIP_CONSDATA* consdata,
   int* ncuts
   )
{
   GRAPH*  g;
   double* xval;
   int    i;
   int    k;
   int    j;
   int    ind;
   int    layer;
   int    count = 0;
   double sum;
   SCIP_ROW* row = NULL;
   SCIP_VAR** vars;
   int flowsep;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conshdlrdata != NULL);

   vars = SCIPprobdataGetVars(scip);
   flowsep = conshdlrdata->flowsep;


   g = consdata->graph;
   assert(g != NULL);

   xval = SCIPprobdataGetXval(scip, NULL);
   assert(xval != NULL);

   for(i = 0; i < g->knots; i++)
   {
      for(layer = 0; layer < g->layers; layer++)
      {
         /* Fuer die Quelle koennen wir nichts tun.
          */
         if (i == g->source[layer])
            continue;

         /* Bei Terminal: Summe Input == 1
          * (Ist eigendlich ein Schnitt (Starcut))
          */
         if (g->term[i] == layer)
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "term", 1.0,
                  1.0, FALSE, FALSE, TRUE) );

            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

            sum = 0.0;

            for(k = g->inpbeg[i]; k != EAT_LAST; k = g->ieat[k])
            {
               ind  = layer * g->edges + k;
               sum += (xval != NULL) ? xval[ind] : 0.0;

               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], 1.0) );
            }
            if (fabs(sum - 1.0) > EPSILON)
            {
               SCIP_Bool infeasible;

               SCIP_CALL( SCIPflushRowExtensions(scip, row) );

               SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &infeasible) );
               count++;
            }

            SCIP_CALL( SCIPreleaseRow(scip, &row) );

         }
         /* Etwa keine Fluesse ?
          */
         if (!flowsep)
            continue;

         /* Jeder der Rausgeht, muss kleiner der Summe derer
          * die reingehen sein.
          */
         for(j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j])
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "flow", 0.0, SCIPinfinity(scip),
                  FALSE, FALSE, TRUE) );

            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

            ind = layer * g->edges + j;
            sum = (xval != NULL) ? -xval[ind] : -1.0;

            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], -1.0) );

            for(k = g->inpbeg[i]; k != EAT_LAST; k = g->ieat[k])
            {
               ind  = layer * g->edges + k;
               sum += (xval != NULL) ? xval[ind] : 0.0;

               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], 1.0) );
            }
            if (sum + EPSILON < 0.0)
            {
               SCIP_Bool infeasible;

               SCIP_CALL( SCIPflushRowExtensions(scip, row) );

               SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &infeasible) );
               count++;
            }
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
         }

         /* Ab hier nur noch nicht Terminals
          */
         if (g->term[i] == layer)
            continue;

         /* Der Input eines Knotens kann nicht groesser als 1.0 sein
          */
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "infl", -SCIPinfinity(scip),
               1.0, FALSE, FALSE, TRUE) );

         SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

         sum   = 0.0;

         for(k = g->inpbeg[i]; k != EAT_LAST; k = g->ieat[k])
         {
            ind  = layer * g->edges + k;
            sum += (xval != NULL) ? xval[ind] : 1.0;

            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], 1.0) );
         }
         if (sum - EPSILON > 1.0)
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPflushRowExtensions(scip, row) );

            SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &infeasible) );
            count++;
         }

         SCIP_CALL( SCIPreleaseRow(scip, &row) );

         /* Was reingeht, muss mindestens auch rausgehen,
          * bei 2-termialen Netzen muss es gleich sein.
          */
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "bala", 0.0,
               (g->locals[layer] == 2) ? 0.0 : SCIPinfinity(scip), FALSE, FALSE, TRUE) );

         SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
         sum   = 0.0;

         for(k = g->inpbeg[i]; k != EAT_LAST; k = g->ieat[k])
         {
            ind  = layer * g->edges + k;
            sum -= (xval != NULL) ? xval[ind] : 1.0;

            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], -1.0) );
         }
         for(k = g->outbeg[i]; k != EAT_LAST; k = g->oeat[k])
         {
            ind  = layer * g->edges + k;
            sum += (xval != NULL) ? xval[ind] : 0.0;

            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[ind], 1.0) );
         }
         if (sum + EPSILON < 0.0)
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPflushRowExtensions(scip, row) );

            SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &infeasible) );
            count++;
         }
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }
   }
   printf("In/Out Separator: %d Inequalities added\n", count);

   return(count);
}


static
int sep_2cut(
   SCIP* scip,
   SCIP_CONSHDLR* conshdlr,
   SCIP_CONSHDLRDATA* conshdlrdata,
   SCIP_CONSDATA* consdata,
   int* ncuts
   )
{
   const int nested_cut   = conshdlrdata->nestedcut;
   const int back_cut     = conshdlrdata->backcut;
   const int creep_flow   = conshdlrdata->creepflow;
   const int disjunct_cut = conshdlrdata->disjunctcut;
   const int flowsep      = conshdlrdata->flowsep;

   GRAPH*  g;
   double* xval;
   PATH*   path;
   int*    term;
   int     terms = 0;
   int     tsave;
   int*    capa;
   double* cost;
   int*    w;
   int     i;
   int     k;
   int     layer;
   int     count = 0;
   int     rerun = FALSE;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conshdlrdata != NULL);

   g = consdata->graph;
   assert(g != NULL);

   xval = SCIPprobdataGetXval(scip, NULL);
   assert(xval != NULL);

   capa  = malloc((size_t)g->edges * sizeof(int));
   cost  = malloc((size_t)g->edges * sizeof(double));
   w     = calloc((size_t)g->knots, sizeof(int));
   term  = malloc((size_t)g->terms * sizeof(int));
   path  = malloc((size_t)g->knots * sizeof(PATH));

   assert(capa != NULL);
   assert(cost != NULL);
   assert(w    != NULL);
   assert(term != NULL);
   assert(path != NULL);

   for(layer = 0; layer < g->layers; layer++)
   {
      /* Fuer 2-terminale Netze brauchen wir keine Schnitte wenn die Flows
       * da sind.
       */
      if (flowsep && (g->locals[layer] == 2))
         continue;

      for(i = 0; i < g->edges; i++)
         cost[i] = (xval[layer * g->edges + i] < (1.0 - SCIPepsilon(scip))) ? 1.0 : 0.0;

      for(i = 0; i < g->knots; i++)
         g->mark[i] = TRUE;

      graph_path_exec(g, FSP_MODE, g->source[layer], cost, path);

      for(i = 0, count = 0; i < g->knots; i++)
      {
         if ((g->term[i] == layer) && (i != g->source[layer]))
         {
            if (SCIPisPositive(scip, path[i].dist))
               term[terms++] = i;
            else
               count++;
         }
      }
      printf("Cut Pretest: %d eliminations\n", count);

      count = 0;
      tsave = terms;

      /* from source to terminal
       */
      if (!nested_cut || disjunct_cut)
         set_capacity(g, layer, creep_flow, 0, capa, xval);

      while(terms > 0)
      {
         /* Wir suchen ein Terminal zu dem wir gehen koennen
          */
         i = graph_next_term(terms, term, w);

         terms--;

         assert(g->term[i]       == layer);
         assert(g->source[layer] != i);

         if (nested_cut && !disjunct_cut)
            set_capacity(g, layer, creep_flow, 0, capa, xval);

         do
         {
            graph_mincut_exec(g, g->source[layer], i, capa, w, rerun);

            rerun = TRUE;

            /* Die Welt wird zweigeteilt, daher "Schnitt"
             */
            for(k = 0; k < g->knots; k++)
               g->mark[k] = (w[k] != 0);

            if (cut_add(scip, conshdlr, g, layer, C_T_2CUT, xval, ncuts))
               count++;
            else
               break;
#if 0
            if (nested_cut || disjunct_cut)
               for(k = p->beg[p->rcnt - 1]; k < p->nzcnt; k++)
                  capa[p->ind[k] % g->edges] = FLOW_FACTOR;
#endif
         }
         while(nested_cut);               /* Nested Cut ist KONSTANT ! */
      }

      /* Auch den Rueckweg versuchen ?
       */
      if (back_cut)
      {
         if (!nested_cut || disjunct_cut)
            set_capacity(g, layer, creep_flow, 1, capa, xval);

         terms = tsave;

         while(terms > 0)
         {
            /* Wir suchen ein Terminal zu dem wir gehen koennen
             */
            i = graph_next_term(terms, term, w);

            terms--;

            assert(g->term[i]       == layer);
            assert(g->source[layer] != i);

            if (nested_cut && !disjunct_cut)
               set_capacity(g, layer, creep_flow, 1, capa, xval);

            rerun = FALSE;

            do
            {
               graph_mincut_exec(g, i, g->source[layer], capa, w, rerun);

               rerun = TRUE;

               for(k = 0; k < g->knots; k++)
                  g->mark[k] = (w[k] != 0) ? 1 : 0;

               if (cut_add(scip, conshdlr, g, layer, C_T_2CUT, xval, ncuts))
                  count++;
               else
                  break;
#if 0
               if (nested_cut || disjunct_cut)
                  for(k = p->beg[p->rcnt - 1]; k < p->nzcnt; k++)
                     capa[p->ind[k] % g->edges
                        + (((p->ind[k] % g->edges) % 2)
                           ? -1 : 1)] = FLOW_FACTOR;
#endif
            }
            while(nested_cut);                /* Nested Cut ist KONSTANT ! */

            rerun = FALSE;
         }
      }
   }
   free(path);
   free(term);
   free(w);
   free(cost);
   free(capa);

   printf("2-cut Separator: %d Inequalities added\n", count);

   return(count);
}



/**@} */


/**@name Callback methods
 *
 * @{
 */

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
   int ncuts = 0;
   int i;

   *result = SCIP_DIDNOTRUN;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[i]);

      sep_flow(scip, conshdlr, conshdlrdata, consdata, &ncuts);

      sep_2cut(scip, conshdlr, conshdlrdata, consdata, &ncuts);
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

      feasible = validate(consdata->graph, SCIPprobdataGetXval(scip, NULL));

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

      feasible = validate(consdata->graph, SCIPprobdataGetXval(scip, NULL));

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

      feasible = validate(consdata->graph, SCIPprobdataGetXval(scip, sol));

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
   {
      SCIP_CALL( SCIPaddVarLocks(scip, vars[v], 1, 1) );
   }

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
   SCIPallocMemory(scip, &conshdlrdata);
   /* TODO: (optional) create constraint handler specific data here */

   conshdlr = NULL;
   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpStp, consEnfopsStp, consCheckStp, consLockStp,
         conshdlrdata) );
   assert(conshdlr != NULL);

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


   return SCIP_OKAY;
}

/** creates and captures a stp constraint */
SCIP_RETCODE SCIPcreateConsStp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   GRAPH*                graph
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

   SCIPdebugMessage("created constraint: ");


   return SCIP_OKAY;
}

/**@} */
