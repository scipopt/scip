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

/**@file   mincut.c
 * @brief  Minimum cut routines for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file encompasses minimum cut routines for Steiner tree problems.
 *
 */

#include "mincut.h"
#include "probdata_stp.h"
#include "portab.h"

#define Q_NULL     -1         /* NULL element of queue/list */
#define ADDCUTSTOPOOL FALSE

/* *
#define FLOW_FACTOR     100000
#define CREEP_VALUE     1         this is the original value todo check what is better
*/

#define FLOW_FACTOR     1000000
#define CREEP_VALUE     10




/*
 * Local methods
 */

//#define STP_MAXFLOW_WRITE

#ifdef STP_MAXFLOW_WRITE
/* writes flow problem in extended dimacs format */
static
void writeFlowProb(
   const GRAPH*          g,                  /**< graph data structure */
   const int*            capa,               /**< edge capacity */
   int                   sinkterm            /**< sink */
   )
{
   FILE *fptr;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);

   fptr = fopen("flow", "w+");
   assert(fptr != NULL);

   fprintf(fptr, "p max %d %d \n", nnodes, nedges);
   fprintf(fptr, "n %d s \n", g->source + 1);
   fprintf(fptr, "n %d t \n", sinkterm + 1);

   for( int k = 0; k < nnodes; k++ )
   {
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         fprintf(fptr, "a %d %d %d \n", k + 1, g->head[e] + 1, capa[e]);
      }
   }

   fprintf(fptr, "x\n");

   fclose(fptr);
}
#endif



/** add a cut */
static
SCIP_RETCODE lpcutAdd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   const GRAPH*          g,                  /**< graph data structure */
   const int*            nodes_inRootComp,   /**< for each node: in root component? */
   const STP_Bool*       edges_isRemoved,    /**< for each edge: removed? */
   const SCIP_Real*      xvals,              /**< edge values */
   int*                  capa,               /**< edges capacities (scaled) */
   const int             updatecapa,         /**< update capacities? */
   SCIP_Bool             local,              /**< is the cut local? */
   SCIP_Bool*            success             /**< pointer to store whether add cut be added */
   )
{
   SCIP_ROW* row;
   SCIP_VAR** vars = SCIPprobdataGetVars(scip);
   SCIP_Real sum = 0.0;
   const int nedges = graph_get_nEdges(g);
   const int* const gtail = g->tail;
   const int* const ghead = g->head;

   assert(g);
   assert(g->knots > 0);
   assert(xvals);
   assert(!updatecapa);

   (*success) = FALSE;

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "twocut", 1.0, SCIPinfinity(scip), local, FALSE, TRUE) );

   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   assert(nodes_inRootComp[g->source]);

   for( int i = 0; i < nedges; i++ )
   {
      if( !nodes_inRootComp[ghead[i]] && nodes_inRootComp[gtail[i]] )
      {
#ifdef STP_USE_ADVANCED_FLOW
         if( updatecapa )
         {
            const SCIP_Bool inccapa = (capa[e] < FLOW_FACTOR);

            capa[e] = FLOW_FACTOR;

            if( !inccapa )
            {
               SCIP_CALL(SCIPflushRowExtensions(scip, row));
               SCIP_CALL(SCIPreleaseRow(scip, &row));
               return SCIP_OKAY;
            }
         }
#endif

         if( edges_isRemoved[i] )
         {
            assert(EQ(xvals[i], 0.0));
            continue;
         }

         sum += xvals[i];

         if( SCIPisFeasGE(scip, sum, 1.0) )
         {
            SCIP_CALL( SCIPflushRowExtensions(scip, row) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
          //  printf("Discard cut! %.18f >= 1.0 \n", sum);

            return SCIP_OKAY;
         }

         SCIP_CALL(SCIPaddVarToRow(scip, row, vars[i], 1.0));
      }
   }

   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   /* checks whether cut is sufficiently violated */
   if( SCIPisCutEfficacious(scip, NULL, row) )
   {
      SCIP_Bool infeasible;

      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

      SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
      assert(!infeasible);

#if ADDCUTSTOPOOL
      /* if at root node, add cut to pool */
      if( !infeasible )
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif
      (*success) = TRUE;
   }

   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}


/** sets (integer) capacity for edges*/
static
void lpcutSetEdgeCapacity(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      xval,               /**< edge values */
   SCIP_Bool             creep_flow,         /**< creep flow? */
   SCIP_Bool             flip,               /**< reverse the flow? */
   int* RESTRICT         capa                /**< edges capacities (scaled) */
   )
{
   int k;
   int krev;
   int nedges = g->edges;

   assert(0 && "currently not supported!");
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


/** returns next promising terminal for computing cut */
static
int getNextSinkTerm(
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
      int todo; // probably got to randomize this!
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


/*
 * Interface methods
 */


/** separates Steiner cuts for LP */
SCIP_RETCODE mincut_separateLp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   const int*            termorg,            /**< original terminals or NULL */
   GRAPH*                g,                  /**< graph data structure */
   int                   maxncuts,           /**< maximal number of cuts */
   int*                  ncuts               /**< pointer to store number of cuts */
   )
{
   /* we do not longer support any other flow, since they slow everything down and are of
    * little use anyway todo remove user parameter */
#ifdef SCIP_DISABLED
   const SCIP_Bool nested_cut   = conshdlrdata->nestedcut;
   const SCIP_Bool back_cut     = conshdlrdata->backcut;
   const SCIP_Bool creep_flow   = conshdlrdata->creepflow;
   const SCIP_Bool disjunct_cut = conshdlrdata->disjunctcut;
#endif
   const SCIP_Bool nested_cut   = FALSE;
   const SCIP_Bool creep_flow   = TRUE;
   const SCIP_Bool disjunct_cut = FALSE;
   const SCIP_Bool intree = (SCIPgetDepth(scip) > 0);

   SCIP_VAR** vars = SCIPprobdataGetVars(scip);
   STP_Bool* edges_isRemoved;
   const SCIP_Real* xval = SCIPprobdataGetXval(scip, NULL);
   int* w;
   int* capa;
   int* terms;
   int* csr_start;
   int* excess;
   int* rootcut;
   int* edgearr;
   int* csr_headarr;
   int* residual;
   int* edgecurr;
   int* csr_edgeflipped;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   const int root = g->source;
   int cutcount;
   int ntermcands;
   int csr_nedges;
   int rootcutsize;
   SCIP_Bool rerun;

   assert(conshdlr != NULL);
   assert(xval && vars);

   excess = g->mincut_e;
   residual = g->mincut_r;
   edgecurr = g->mincut_numb;
   assert(excess && residual && edgecurr);
   assert(creep_flow == TRUE);
   assert(nested_cut == FALSE);
   assert(disjunct_cut == FALSE);

   SCIP_CALL( SCIPallocBufferArray(scip, &capa, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &w, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &terms, g->terms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearr, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &csr_headarr, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &csr_edgeflipped, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &csr_start, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rootcut, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edges_isRemoved, nedges) );

   /* sets excess, g->mincut_head,  g->mincut_head_inact */
   graph_mincut_setDefaultVals(g);

#ifdef STP_MAXFLOW_TIME
   clock_t startt, endt;
   double cpu_time_used;
   startt = clock();
#endif

   for( int i = 0; i < nedges; i++ )
   {
      edges_isRemoved[i] = (SCIPvarGetUbGlobal(vars[i]) < 0.5);
   }

   for( int k = 0; k < nnodes; k++ )
   {
      w[k] = 0;
   }


   // method set capa from LP todo

   for( int e = 0; e < nedges; e += 2 )
   {
      const int erev = e + 1;
      SCIP_Bool eIsRemoved;

      if( intree )
         eIsRemoved = SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[erev]) < 0.5;
      else
         eIsRemoved = edges_isRemoved[e] && edges_isRemoved[erev];

      assert(eIsRemoved || !edges_isRemoved[e] || !edges_isRemoved[erev]);

      if( eIsRemoved )
      {
         int todo; // what is the point of setting residual here???
         capa[e] = 0;
         capa[erev] = 0;
         residual[e] = 0;
         residual[erev] = 0;

         csr_headarr[e] = 1;
         csr_headarr[erev] = 1;
      }
      else
      {
         capa[e]     = (int)(xval[e] * FLOW_FACTOR + 0.5) + CREEP_VALUE;
         capa[erev]  = (int)(xval[erev] * FLOW_FACTOR + 0.5) + CREEP_VALUE;
         residual[e] = capa[e];
         residual[erev] = capa[erev];

         csr_headarr[e] = SCIPisFeasGE(scip, xval[e], 1.0) ? 0 : 1;
         csr_headarr[erev] = SCIPisFeasGE(scip, xval[erev], 1.0) ? 0 : 1;
      }
   }

   /*
    * BFS along 0 edges from the root
    * */

   w[root] = 1;
   rootcutsize = 0;
   rootcut[rootcutsize++] = root;

   /* BFS loop */
   for( int i = 0; i < rootcutsize; i++ )
   {
      const int k = rootcut[i];

      assert(rootcutsize <= nnodes);
      assert(k < nnodes);

      /* traverse outgoing arcs */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];

         /* head not been added to root cut yet? */
         if( w[head] == 0 )
         {
            if( csr_headarr[e] == 0 )
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

   for( int e = 0; e < nedges; e++ )
      edgearr[e] = -1;

   csr_nedges = 0;
   ntermcands = 0;

   /* fill auxiliary adjacent vertex/edges arrays and get useable terms */
   for( int k = 0; k < nnodes; k++ )
   {
      csr_start[k] = csr_nedges;

      /* non-dormant node? */
      if( w[k] == 0 )
      {
         edgecurr[k] = csr_nedges;
         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];
            if( w[head] == 0 && capa[e] != 0 )
            {
               edgearr[e] = csr_nedges;
               residual[csr_nedges] = capa[e];
               csr_headarr[csr_nedges++] = head;
            }
         }

         /* unreachable node? */
         if( edgecurr[k] == csr_nedges )
         {
            w[k] = 1;
         }
         else if( Is_term(g->term[k]) )
         {
            terms[ntermcands++] = k;
         }
      }
      else
      {
         edgecurr[k] = -1;
      }
   }

   csr_start[nnodes] = csr_nedges;

   /* initialize edgeflipped */
   for( int e = 0; e < nedges; e++ )
   {
      if( edgearr[e] >= 0 )
      {
         const int csr_pos = edgearr[e];
         csr_edgeflipped[csr_pos] = edgearr[flipedge(e)];
      }
   }

   SCIPdebugMessage("Cut Pretest: %d eliminations\n", g->terms - ntermcands - 1);

#ifdef STP_MAXFLOW_TIME
   endt = clock();
   cpu_time_used = ((double) (endt - startt)) / CLOCKS_PER_SEC;
   startt = clock();
#endif

   cutcount = 0;
   rerun = FALSE;

   while( ntermcands > 0 )
   {
      int sinkterm;
      SCIP_Bool addedcut = FALSE;

      if( ((unsigned) ntermcands) % 8 == 0 && SCIPisStopped(scip) )
         break;

      /* look for non-reachable terminal */
      sinkterm = getNextSinkTerm(g, ntermcands, terms, w, !rerun);
      ntermcands--;

      assert(Is_term(g->term[sinkterm]) && root != sinkterm);

      if( nested_cut && !disjunct_cut )
         lpcutSetEdgeCapacity(g, xval, creep_flow, 0, capa);

      do
      {
         /* declare cuts on branched-on (artificial) terminals as local */
         const SCIP_Bool localcut = (termorg != NULL && termorg[sinkterm] != g->term[sinkterm]);

#ifdef STP_MAXFLOW_WRITE
         writeFlowProb(g, capa, sinkterm);
#endif

         /* non-trivial cut? */
         if( w[sinkterm] != 1 )
         {
            graph_mincut_exec(g, root, sinkterm, nnodes, csr_nedges, rootcutsize, rootcut, capa, w, csr_start, csr_edgeflipped, csr_headarr, rerun);

            assert(w[root] != 0);

            SCIP_CALL( lpcutAdd(scip, conshdlr, g, w, edges_isRemoved, xval, capa, nested_cut || disjunct_cut, localcut, &addedcut) );
         }
         else
         {
            SCIP_Real flowsum = 0.0;

            assert(rerun);

            for( int e = g->inpbeg[sinkterm]; e != EAT_LAST; e = g->ieat[e] )
               flowsum += xval[e];

            if( SCIPisFeasGE(scip, flowsum, 1.0) )
               continue;

            for( int k = 0; k < nnodes; k++ )
               g->mark[k] = TRUE;

            g->mark[sinkterm] = FALSE;

            SCIP_CALL( lpcutAdd(scip, conshdlr, g, g->mark, edges_isRemoved, xval, capa, nested_cut || disjunct_cut, localcut, &addedcut) );
         }

         rerun = TRUE;

         if( addedcut )
         {
            cutcount++;
            (*ncuts)++;

            if( *ncuts >= maxncuts )
               goto TERMINATE;
         }
         else
         {
            break;
         }
      }
      while( nested_cut );               /* Nested Cut is CONSTANT ! */
   } /* while terms > 0 */


#ifdef STP_MAXFLOW_TIME
   endt = clock();
   cpu_time_used = ((double) (endt - startt)) / CLOCKS_PER_SEC;
#endif

#ifdef SCIP_DISABLED
      /*
       * back cuts currently not supported
       *  */
      /* back cuts enabled? */
      if( back_cut )
      {
         for( k = 0; k < nnodes; k++ )
            w[k] = 0;

         if( !nested_cut || disjunct_cut )
            lpcutSetEdgeCapacity(g, creep_flow, 1, capa, xval);

         ntermcands = tsave;

         while( ntermcands > 0 )
         {
            /* look for reachable terminal */
            i = getNextSinkTerm(g, ntermcands, terms, w, TRUE);

            ntermcands--;

            assert(g->terms[i]       == 0);
            assert(g->source != i);

            if( nested_cut && !disjunct_cut )
               lpcutSetEdgeCapacity(g, creep_flow, 1, capa, xval);

            rerun = FALSE;

            do
            {
               graph_mincut_exec(g, i, g->source, nedges, capa, w, csr_start, edgearr, csr_headarr, rerun);

               rerun = TRUE;

               for( k = 0; k < nnodes; k++ )
               {
                  g->mark[k] = (w[k] != 0) ? 0 : 1; // todo not the other way around??
                  w[k] = 0;
               }

               SCIP_CALL( lpcutAdd(scip, conshdlr, g, xval, capa, nested_cut || disjunct_cut, ncuts, &addedcut) );
               if( addedcut )
               {
                  cutcount++;

                  if( *ncuts >= maxncuts )
                     goto TERMINATE;
               }
               else
                  break;
#ifdef SCIP_DISABLED
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
   SCIPfreeBufferArray(scip, &edges_isRemoved);
   SCIPfreeBufferArray(scip, &rootcut);
   SCIPfreeBufferArray(scip, &csr_start);
   SCIPfreeBufferArray(scip, &csr_edgeflipped);
   SCIPfreeBufferArray(scip, &csr_headarr);
   SCIPfreeBufferArray(scip, &edgearr);

   SCIPfreeBufferArray(scip, &terms);
   SCIPfreeBufferArray(scip, &w);

   SCIPfreeBufferArray(scip, &capa);

   SCIPdebugMessage("2-cut Separator: %d Inequalities added\n", cutcount);

   return SCIP_OKAY;
}
