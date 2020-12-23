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

/** minimum cut helper */
typedef struct minimum_cut_helper
{
   const SCIP_Real*      xval;
   int*                  nodes_wakeState;
   int*                  edges_capa;
   int*                  terms;
   int*                  csr_start;
   int*                  rootcut;
   int*                  csr_edgeDefaultToCsr;
   int*                  csr_headarr;
   int*                  csr_edgeflipped;
   STP_Bool*             edges_isRemoved;    /**< only used for LP cuts */
   int                   ntermcands;
   int                   rootcutsize;
   int                   csr_nedges;
   SCIP_Bool             isLpcut;            /**< cut for LP? */
} MINCUT;



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


/** initializes */
static
SCIP_RETCODE mincutInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   MINCUT**              mincut              /**< minimum cut */
)
{
   MINCUT* mcut;
   int* nodes_wakeState;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);

   assert(scip);
   assert(nnodes > 0 && nedges > 0);

   SCIP_CALL( SCIPallocMemory(scip, mincut) );
   mcut = *mincut;

   SCIP_CALL( SCIPallocBufferArray(scip, &(mcut->edges_capa), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mcut->nodes_wakeState), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mcut->terms), g->terms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mcut->csr_edgeDefaultToCsr), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mcut->csr_headarr), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mcut->csr_edgeflipped), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mcut->csr_start), nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mcut->rootcut), nnodes + 1) );

   mcut->xval = NULL;
   mcut->edges_isRemoved = NULL;
   mcut->isLpcut = FALSE;
   mcut->ntermcands = -1;
   mcut->rootcutsize = -1;
   mcut->csr_nedges = -1;
   nodes_wakeState = mcut->nodes_wakeState;

   for( int k = 0; k < nnodes; k++ )
   {
      nodes_wakeState[k] = 0;
   }


   return SCIP_OKAY;
}


/** initializes for LP */
static
SCIP_RETCODE mincutInitForLp(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   MINCUT*               mincut              /**< minimum cut */
)
{
   SCIP_VAR** vars = SCIPprobdataGetVars(scip);
   STP_Bool* RESTRICT edges_isRemoved;
   const int nedges = graph_get_nEdges(g);

   assert(mincut && scip);
   assert(!mincut->isLpcut);
   assert(!mincut->edges_isRemoved);

   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->edges_isRemoved), nedges) );
   mincut->isLpcut = TRUE;
   edges_isRemoved = mincut->edges_isRemoved;

   for( int i = 0; i < nedges; i++ )
   {
      edges_isRemoved[i] = (SCIPvarGetUbGlobal(vars[i]) < 0.5);
   }

   mincut->xval = SCIPprobdataGetXval(scip, NULL);
   assert(mincut->xval);

   return SCIP_OKAY;
}


/** initializes for LP */
static
SCIP_RETCODE mincutPrepareForLp(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   MINCUT*               mincut              /**< minimum cut */
)
{
   SCIP_VAR** vars = SCIPprobdataGetVars(scip);
   int* RESTRICT nodes_wakeState = mincut->nodes_wakeState;
   int* RESTRICT edges_capa = mincut->edges_capa;
   int* RESTRICT terms = mincut->terms;
   int* RESTRICT csr_start = mincut->csr_start;
   int* RESTRICT rootcut = mincut->rootcut;
   int* RESTRICT csr_edgeDefaultToCsr = mincut->csr_edgeDefaultToCsr;
   int* RESTRICT csr_headarr = mincut->csr_headarr;
   int* RESTRICT csr_edgeflipped = mincut->csr_edgeflipped;
   const STP_Bool* const edges_isRemoved = mincut->edges_isRemoved;
   const SCIP_Real* const xval = mincut->xval;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   const int root = g->source;
   int rootcutsize;
   int csr_nedges;
   int ntermcands;
   const SCIP_Bool intree = (SCIPgetDepth(scip) > 0);
   int* RESTRICT excess = g->mincut_e;
   int* RESTRICT residual = g->mincut_r;
   int* RESTRICT edgecurr = g->mincut_numb;

   assert(residual && edgecurr);
   assert(xval && vars && edges_isRemoved);
   assert(mincut->isLpcut);

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
         edges_capa[e] = 0;
         edges_capa[erev] = 0;
         residual[e] = 0;
         residual[erev] = 0;

         csr_headarr[e] = 1;
         csr_headarr[erev] = 1;
      }
      else
      {
         edges_capa[e]     = (int)(xval[e] * FLOW_FACTOR + 0.5) + CREEP_VALUE;
         edges_capa[erev]  = (int)(xval[erev] * FLOW_FACTOR + 0.5) + CREEP_VALUE;
         residual[e] = edges_capa[e];
         residual[erev] = edges_capa[erev];

         /* NOTE: here we misuse csr_headarr to mark the free edges for the BFS */
         csr_headarr[e] = SCIPisFeasGE(scip, xval[e], 1.0) ? 0 : 1;
         csr_headarr[erev] = SCIPisFeasGE(scip, xval[erev], 1.0) ? 0 : 1;
      }
   }

   /*
    * BFS along 0 edges from the root
    * */


   nodes_wakeState[root] = 1;
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
         if( nodes_wakeState[head] == 0 )
         {
            if( csr_headarr[e] == 0 )
            {
               nodes_wakeState[head] = 1;
               rootcut[rootcutsize++] = head;
            }
            else
            {
               /* push as much as possible out of perpetually dormant nodes (possibly to other dormant nodes) */
               assert(nodes_wakeState[head] == 0);
#ifndef NDEBUG
               residual[e] = 0;
#endif
               excess[head] += edges_capa[e];
            }
         }
      }
   }

   for( int e = 0; e < nedges; e++ )
      csr_edgeDefaultToCsr[e] = -1;

   csr_nedges = 0;
   ntermcands = 0;

   /* fill auxiliary adjacent vertex/edges arrays and get useable terms */
   for( int k = 0; k < nnodes; k++ )
   {
      csr_start[k] = csr_nedges;

      /* non-dormant node? */
      if( nodes_wakeState[k] == 0 )
      {
         edgecurr[k] = csr_nedges;
         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];
            if( nodes_wakeState[head] == 0 && edges_capa[e] != 0 )
            {
               csr_edgeDefaultToCsr[e] = csr_nedges;
               residual[csr_nedges] = edges_capa[e];
               csr_headarr[csr_nedges++] = head;
            }
         }

         /* unreachable node? */
         if( edgecurr[k] == csr_nedges )
         {
            nodes_wakeState[k] = 1;
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

   mincut->ntermcands = ntermcands;
   mincut->rootcutsize = rootcutsize;
   mincut->csr_nedges = csr_nedges;

   /* initialize edgeflipped */
   for( int e = 0; e < nedges; e++ )
   {
      if( csr_edgeDefaultToCsr[e] >= 0 )
      {
         const int csr_pos = csr_edgeDefaultToCsr[e];
         csr_edgeflipped[csr_pos] = csr_edgeDefaultToCsr[flipedge(e)];
      }
   }

   SCIPdebugMessage("Cut Pretest: %d eliminations\n", g->terms - ntermcands - 1);

   return SCIP_OKAY;
}




/** returns next promising terminal for computing cut */
static
int mincutGetNextSinkTerm(
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Bool             firstrun,            /**< first run?  */
   MINCUT*               mincut              /**< minimum cut */
   )
{
   int* RESTRICT termcands = mincut->terms;
   const int* const w = mincut->nodes_wakeState;
   int k;
   int t;
   int mindist = g->knots + 1;
   const int ntermcands = mincut->ntermcands;

   assert(termcands != NULL);
   assert(ntermcands > 0);

   if( firstrun ) // todo randomize?
   {
      int todo; // probably got to randomize this! just swap with a random element
      assert(w[termcands[ntermcands - 1]] == 0);
      return termcands[ntermcands - 1];
   }

   k = -1;

   for( int i = 0; i < ntermcands; i++ )
   {
      assert(w[termcands[i]] >= 0);

      if( w[termcands[i]] == 0 )
      {
         assert(g->mincut_dist[termcands[i]] < g->knots + 1);

         if( g->mincut_dist[termcands[i]] < mindist )
         {
            k = i;
            mindist = g->mincut_dist[termcands[i]];
         }
      }
   }

   if( k == -1 )
   {
      int wmax = 0;

      for( int i = 0; i < ntermcands; i++ )
      {
         if( w[termcands[i]] > wmax )
         {
            k = i;
            wmax = w[termcands[i]];
            mindist = g->mincut_dist[termcands[i]];
         }
         else if( w[termcands[i]] == wmax && g->mincut_dist[termcands[i]] < mindist )
         {
            assert(wmax != 0);

            k = i;
            mindist = g->mincut_dist[termcands[i]];
         }
      }
   }

   assert(k >= 0);
   assert(k < ntermcands);

   t = termcands[k];
   termcands[k] = termcands[ntermcands - 1];

   return t;
}


/** executes */
static
void mincutExec(
   const GRAPH*          g,                  /**< the graph */
   int                   sinkterm,
   SCIP_Bool             wasRerun,
   MINCUT*               mincut              /**< minimum cut */
)
{
   const int root = g->source;
   const int nnodes = graph_get_nNodes(g);

   graph_mincut_exec(g, root, sinkterm, nnodes, mincut->csr_nedges, mincut->rootcutsize,
         mincut->rootcut, mincut->edges_capa, mincut->nodes_wakeState,
         mincut->csr_start, mincut->csr_edgeflipped, mincut->csr_headarr, wasRerun);
}


/** frees */
static
void mincutFree(
   SCIP*                 scip,               /**< SCIP data structure */
   MINCUT**              mincut              /**< minimum cut */
   )
{
   MINCUT* mcut;

   assert(scip && mincut);
   assert(*mincut);
   mcut = *mincut;

   SCIPfreeBufferArrayNull(scip, &(mcut->edges_isRemoved));
   SCIPfreeBufferArray(scip, &(mcut->rootcut));
   SCIPfreeBufferArray(scip, &(mcut->csr_start));
   SCIPfreeBufferArray(scip, &(mcut->csr_edgeflipped));
   SCIPfreeBufferArray(scip, &(mcut->csr_headarr));
   SCIPfreeBufferArray(scip, &(mcut->csr_edgeDefaultToCsr));
   SCIPfreeBufferArray(scip, &(mcut->terms));
   SCIPfreeBufferArray(scip, &(mcut->nodes_wakeState));
   SCIPfreeBufferArray(scip, &(mcut->edges_capa));

   SCIPfreeMemory(scip, mincut);
}


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
   MINCUT* mincut;
   const SCIP_Real* xval;
   int* nodes_wakeState;
   STP_Bool* edges_isRemoved;
   const int nnodes = graph_get_nNodes(g);
   int cutcount;
   SCIP_Bool wasRerun;

#ifdef STP_MAXFLOW_TIME
   clock_t startt, endt;
   double cpu_time_used;
   startt = clock();
#endif

   assert(conshdlr && ncuts);
   assert(creep_flow == TRUE);
   assert(nested_cut == FALSE);
   assert(disjunct_cut == FALSE);

   SCIP_CALL( mincutInit(scip, g, &mincut) );
   SCIP_CALL( mincutInitForLp(scip, g, mincut) );

   /* sets excess, g->mincut_head,  g->mincut_head_inact */
   graph_mincut_setDefaultVals(g);

#ifdef STP_MAXFLOW_TIME
   endt = clock();
   cpu_time_used = ((double) (endt - startt)) / CLOCKS_PER_SEC;
   startt = clock();
#endif

   SCIP_CALL( mincutPrepareForLp(scip, g, mincut) );

   xval = mincut->xval;
   nodes_wakeState = mincut->nodes_wakeState;
   edges_isRemoved = mincut->edges_isRemoved;
   cutcount = 0;
   wasRerun = FALSE;

   assert(xval && nodes_wakeState && edges_isRemoved);

   while( mincut->ntermcands > 0 )
   {
      int sinkterm;
      SCIP_Bool addedcut = FALSE;

      if( ((unsigned) mincut->ntermcands) % 8 == 0 && SCIPisStopped(scip) )
         break;

      /* look for non-reachable terminal */
      sinkterm = mincutGetNextSinkTerm(g, !wasRerun, mincut);
      mincut->ntermcands--;

      assert(Is_term(g->term[sinkterm]) && root != sinkterm);

      if( nested_cut && !disjunct_cut )
         lpcutSetEdgeCapacity(g, xval, creep_flow, 0, mincut->edges_capa);

      do
      {
         /* declare cuts on branched-on (artificial) terminals as local */
         const SCIP_Bool localcut = (termorg != NULL && termorg[sinkterm] != g->term[sinkterm]);

#ifdef STP_MAXFLOW_WRITE
         writeFlowProb(g, edges_capa, sinkterm);
#endif

         /* non-trivial cut? */
         if( nodes_wakeState[sinkterm] != 1 )
         {
            mincutExec(g, sinkterm, wasRerun, mincut);

            assert(nodes_wakeState[root] != 0);

            SCIP_CALL( lpcutAdd(scip, conshdlr, g, nodes_wakeState, edges_isRemoved, xval,
                  mincut->edges_capa, nested_cut || disjunct_cut, localcut, &addedcut) );
         }
         else
         {
            SCIP_Real flowsum = 0.0;
            assert(wasRerun);

            for( int e = g->inpbeg[sinkterm]; e != EAT_LAST; e = g->ieat[e] )
               flowsum += xval[e];

            if( SCIPisFeasGE(scip, flowsum, 1.0) )
               continue;

            for( int k = 0; k < nnodes; k++ )
               g->mark[k] = TRUE;

            g->mark[sinkterm] = FALSE;

            SCIP_CALL( lpcutAdd(scip, conshdlr, g, g->mark, edges_isRemoved, xval,
                  mincut->edges_capa, nested_cut || disjunct_cut, localcut, &addedcut) );
         }

         wasRerun = TRUE;

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
            nodes_wakeState[k] = 0;

         if( !nested_cut || disjunct_cut )
            lpcutSetEdgeCapacity(g, creep_flow, 1, edges_capa, xval);

         ntermcands = tsave;

         while( ntermcands > 0 )
         {
            /* look for reachable terminal */
            i = mincutGetNextSinkTerm(g, ntermcands, terms, nodes_wakeState, TRUE);

            ntermcands--;

            assert(g->terms[i]       == 0);
            assert(g->source != i);

            if( nested_cut && !disjunct_cut )
               lpcutSetEdgeCapacity(g, creep_flow, 1, edges_capa, xval);

            wasRerun = FALSE;

            do
            {
               graph_mincut_exec(g, i, g->source, nedges, edges_capa, nodes_wakeState, csr_start, csr_edgeDefaultToCsr, csr_headarr, wasRerun);

               wasRerun = TRUE;

               for( k = 0; k < nnodes; k++ )
               {
                  g->mark[k] = (nodes_wakeState[k] != 0) ? 0 : 1; // todo not the other way around??
                  nodes_wakeState[k] = 0;
               }

               SCIP_CALL( lpcutAdd(scip, conshdlr, g, xval, edges_capa, nested_cut || disjunct_cut, ncuts, &addedcut) );
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
                     edges_capa[p->ind[k] % nedges
                        + (((p->ind[k] % nedges) % 2)
                           ? -1 : 1)] = FLOW_FACTOR;
#endif
            }
            while( nested_cut );                /* Nested Cut is CONSTANT todo why not only one round? seems to make no sense whatsoever */

            wasRerun = FALSE;
         }
      }
#endif
 TERMINATE:
   mincutFree(scip, &mincut);

   SCIPdebugMessage("2-cut Separator: %d Inequalities added\n", cutcount);

   return SCIP_OKAY;
}
