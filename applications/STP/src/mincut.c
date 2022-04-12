/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
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

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


//#define SCIP_DEBUG


#include "mincut.h"
#include "probdata_stp.h"
#include "misc_stp.h"
#include "stpvector.h"
#include "reduce.h"
#include "portab.h"

#define Q_NULL     -1         /* NULL element of queue/list */
#define ADDCUTSTOPOOL FALSE
#define TERMSEPA_SPARSE_MAXRATIO 4
#define TERMSEPA_MAXCUTSIZE_DEFAULT 5
#define TERMSEPA_MAXNCUTS_DEFAULT   100
#define TERMSEPA_NROOTCANDS 3

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
   int*                  termsepa_termToCopy;
   int*                  terms_minsepasize;     /**< size of smallest separator for given terminal (non-defined for Steiner nodes) */
   int*                  terms_mincompsize;     /**< size of smallest component for given terminal (non-defined for Steiner nodes) */
   STP_Bool*             edges_isRemoved;    /**< only used for LP cuts */
   int                   ntermcands;
   int                   rootcutsize;
   int                   csr_nedges;
   int                   root;
   int                   termsepa_nnodes;
   int                   termsepa_nedges;
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator or NULL */
   SCIP_Bool             isLpcut;            /**< cut for LP? */
} MINCUT;



/** single separator info */
typedef struct terminal_separator
{
   int                   sinkterm;
   int                   nsinknodes;
} TSEPA;


/** storage */
struct terminal_separator_storage
{
   int*                  nsepas;             /**< of size maxsepasize + 1 */
   int*                  currsepa_n;         /**< of size maxsepasize + 1 */
   TSEPA*                sepas;
   int*                  sepaterms_csr;
   int*                  sepastarts_csr;
   int                   nsepaterms_csr;
   int                   nsepas_all;
   int                   root;
   int                   maxnsepas;          /**< maximum number of separators to compute */
   int                   maxsepasize;        /**< maximum size of separator */
};

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



#ifdef SCIP_DEBUG
static inline
void debugPrintCutNodes(
   const GRAPH*          g,                  /**< the graph */
   const MINCUT*         mincut              /**< minimum cut */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int* const nodes_wakeState = mincut->nodes_wakeState;

   printf("root component nodes: \n");

   for( int i = 0; i < nnodes; i++ )
   {
      if( nodes_wakeState[i] )
         graph_knot_printInfo(g, i);
   }

   printf("remaining nodes: \n");

   for( int i = 0; i < nnodes; i++ )
   {
      if( !nodes_wakeState[i] )
         graph_knot_printInfo(g, i);
   }
}

static inline
void debugPrintCutEdges(
   const GRAPH*          g,                  /**< the graph */
   const MINCUT*         mincut              /**< minimum cut */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int* const nodes_wakeState = mincut->nodes_wakeState;

   printf("cut edges: \n");

   for( int i = 0; i < nnodes; i++ )
   {
      if( !nodes_wakeState[i] )
         continue;

      for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];

         if( nodes_wakeState[head] == 0 )
         {
            graph_edge_printInfo(g, e);
         }
      }
   }
}


static inline
void debugPrintCsrCutEdges(
   const GRAPH*          g,                  /**< the graph */
   const MINCUT*         mincut              /**< minimum cut */
)
{
   const int* const edges_capa = mincut->edges_capa;
   const int* const nodes_wakeState = mincut->nodes_wakeState;
   const int* const csr_start = mincut->csr_start;
   const int* const csr_headarr = mincut->csr_headarr;
   const int nnodes_extended = mincut->termsepa_nnodes;

   assert(!mincut->isLpcut);

   printf("CSR cut edges: \n");

   for( int i = 0; i < nnodes_extended; i++ )
   {
      const int start = csr_start[i];
      const int end = csr_start[i + 1];

      if( !nodes_wakeState[i] )
         continue;

      for( int j = start; j != end; j++ )
      {
         const int head = csr_headarr[j];

         if( nodes_wakeState[head] == 0 )
         {
            printf("%d->%d, c=%d \n", i, head, edges_capa[j]);
         }
      }
   }
}



/** prints extended graph */
static inline
void debugPrintCsr(
   const GRAPH*          g,                  /**< the graph */
   const MINCUT*         mincut              /**< minimum cut */
)
{
   const int* const residual = g->mincut_r;
   const int* const csr_start = mincut->csr_start;
   const int* const csr_headarr = mincut->csr_headarr;
   const int* const csr_edgeflipped = mincut->csr_edgeflipped;
   const int nnodes_extended = mincut->termsepa_nnodes;

   assert(!mincut->isLpcut);

   for( int i = 0; i < nnodes_extended; i++ )
   {
      const int start = csr_start[i];
      const int end = csr_start[i + 1];

      for( int j = start; j != end; j++ )
      {
         const int head = csr_headarr[j];

         printf("edge %d: %d->%d, r=%d (fliped=%d) \n", j, i, head, residual[j], csr_edgeflipped[j]);
      }
   }
}
#endif


#ifndef NDEBUG
/** valid flip-edges? */
static inline
SCIP_Bool csrFlipedgesAreValid(
   const GRAPH*          g,                  /**< the graph */
   const MINCUT*         mincut              /**< minimum cut */
)
{
   const int* const residual = g->mincut_r;
   const int* const csr_start = mincut->csr_start;
   const int* const csr_headarr = mincut->csr_headarr;
   const int* const csr_edgeflipped = mincut->csr_edgeflipped;
   const int nnodes_extended = mincut->termsepa_nnodes;

   assert(!mincut->isLpcut);

   for( int i = 0; i < nnodes_extended; i++ )
   {
      const int start = csr_start[i];
      const int end = csr_start[i + 1];

      for( int j = start; j != end; j++ )
      {
         assert(csr_edgeflipped[j] >= 0);

         if( csr_headarr[csr_edgeflipped[j]] != i )
         {
            return FALSE;
         }

         assert(residual[j] > 0 || residual[csr_edgeflipped[j]] > 0);
      }
   }

   return TRUE;
}


/** checks cut */
static
SCIP_Bool termsepaCutIsCorrect(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   int                   ncutterms,          /**< number of cut terminals */
   const int*            cutterms,           /**< stores cut terminals */
   int                   sinkterm,           /**< sink terminal */
   const MINCUT*         mincut              /**< minimum cut */
)
{
   SCIP_Bool isCorrect = TRUE;
   const int* const nodes_wakeState = mincut->nodes_wakeState;
   STP_Vectype(int) stack = NULL;
   SCIP_Bool* nodes_isBlocked;
   SCIP_Bool* nodes_isVisited;
   const int nnodes = graph_get_nNodes(g);
   int nvisted = 1;
   int nsinknodes = 0;

   assert(ncutterms > 0);
   assert(nodes_wakeState[sinkterm] == 0);

   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &nodes_isBlocked, nnodes) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &nodes_isVisited, nnodes) );

   for( int i = 0; i < nnodes; i++ )
   {
      nodes_isBlocked[i] = FALSE;
      nodes_isVisited[i] = FALSE;
   }

   for( int i = 0; i < ncutterms; i++ )
   {
      assert(graph_knot_isInRange(g, cutterms[i]));
      assert(!nodes_isBlocked[cutterms[i]]);

      nodes_isBlocked[cutterms[i]] = TRUE;
   }

   assert(!nodes_isBlocked[mincut->root]);

   for( int i = 0; i < nnodes; i++ )
   {
      if( nodes_wakeState[i] == 0 )
      {
         nsinknodes++;
      }
   }

   /* DFS from sink terminal */
   StpVecReserve(scip, stack, nnodes);
   StpVecPushBack(scip, stack, sinkterm);
   nodes_isVisited[sinkterm] = TRUE;

   while( StpVecGetSize(stack)  )
   {
      const int k = stack[StpVecGetSize(stack) - 1];
      StpVecPopBack(stack);

      assert(nodes_isVisited[k]);

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         if( !nodes_isVisited[head] && !nodes_isBlocked[head] )
         {
            nvisted++;
            nodes_isVisited[head] = TRUE;
            StpVecPushBack(scip, stack, head);
         }
      }
   }

   StpVecFree(scip, stack);

   if( nvisted > nsinknodes )
   {
      SCIPerrorMessage("nvisted=%d nsinknodes=%d\n", nvisted, nsinknodes);
      isCorrect = FALSE;
   }

   if( nodes_isVisited[mincut->root] )
   {
      SCIPerrorMessage("root is in cut! \n");
      isCorrect = FALSE;
   }


   SCIPfreeMemoryArray(scip, &nodes_isVisited);
   SCIPfreeMemoryArray(scip, &nodes_isBlocked);


   return isCorrect;
}
#endif


/** gets infinity */
static
int termsepaGetCapaInf(
   const GRAPH*          g,                  /**< the graph */
   const MINCUT*         mincut              /**< minimum cut */
)
{
   return g->terms;
}


/** updates; takes root candidate if greater than current */
static
void updateTerminalSource(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   rootcand,           /**< new root candidate */
   const GRAPH*          g,                  /**< the graph */
   int*                  root,               /**< incumbent  (IN/OUT) */
   int*                  rootcompsize        /**< size of component (IN/OUT) */
)
{
   STP_Vectype(int) queue = NULL;
   SCIP_Bool* nodes_isVisited;
   const int nnodes = graph_get_nNodes(g);
   int compsize = 0;

   assert(graph_knot_isInRange(g, rootcand) && Is_term(g->term[rootcand]));

   SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &nodes_isVisited, nnodes) );
   StpVecReserve(scip, queue, nnodes);
   StpVecPushBack(scip, queue, rootcand);
   nodes_isVisited[rootcand] = TRUE;

   /* BFS loop stopping at roots */
   for( int i = 0; i < StpVecGetSize(queue); i++ )
   {
      const int node = queue[i];
      assert(nodes_isVisited[node]);

      compsize++;

      if( Is_term(g->term[node]) && node != rootcand )
         continue;

      for( int e = g->outbeg[node]; e >= 0; e = g->oeat[e] )
      {
         const int head = g->head[e];
         if( !nodes_isVisited[head] )
         {
            nodes_isVisited[head] = TRUE;
            StpVecPushBack(scip, queue, head);
         }
      }
   }

   for( int i = 0; i < StpVecGetSize(queue); i++ )
   {
      const int node = queue[i];
      assert(nodes_isVisited[node]);
      nodes_isVisited[node] = FALSE;
   }

   StpVecFree(scip, queue);

#ifndef NDEBUG
   for( int i = 0; i < nnodes; i++ )
   {
      assert(nodes_isVisited[i] == FALSE);
   }
#endif

   if( compsize > *rootcompsize )
   {
      *rootcompsize = compsize;
      *root = rootcand;
   }

   SCIPfreeCleanBufferArray(scip, &nodes_isVisited);
}


/** finds a root terminal */
static
int termsepaFindTerminalSource(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   const MINCUT*         mincut              /**< minimum cut */
)
{
   int source = -1;

   /* NOTE: hack to allow for stable unit tests  */
   if( !mincut->randnumgen )
   {
      source = g->source;
   }
   else
   {
      int* terms;
      int sourcecompsize = 0;
      const int ntries = MIN(TERMSEPA_NROOTCANDS, g->terms);

      SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &terms, g->terms) );
      graph_getTerms(g, terms);

      for( int i = 0; i < ntries; i++ )
      {
         const int pos = SCIPrandomGetInt(mincut->randnumgen, 0, g->terms - 1);
         const int cand = terms[pos];

         updateTerminalSource(scip, cand, g, &source, &sourcecompsize );
      }

      assert(sourcecompsize > 0);

      SCIPfreeMemoryArray(scip, &terms);
   }

   assert(source >= 0);
   assert(Is_term(g->term[source]));

   return source;
}


/** collect cut nodes */
static
void termsepaCollectCutNodes(
   const GRAPH*          g,                  /**< the graph */
   const TERMSEPAS*      termsepas,          /**< terminal separator storage */
   const MINCUT*         mincut,             /**< minimum cut */
   int                   sinkterm,           /**< sink terminal */
   int*                  cutterms,           /**< terminals */
   int*                  ncutterms,          /**< number of terminals */
   SCIP_Bool*            cutIsGood           /**< is cut nodes */
)
{
   const int* const edges_capa = mincut->edges_capa;
   const int* const nodes_wakeState = mincut->nodes_wakeState;
   const int* const csr_start = mincut->csr_start;
   const int* const csr_headarr = mincut->csr_headarr;
   const int nnodes_extended = mincut->termsepa_nnodes;
   const int capa_inf = termsepaGetCapaInf(g, mincut);
   const int maxsepasize = termsepas->maxsepasize;
   SCIP_Bool isGood = TRUE;
   int n = 0;

   for( int i = 0; i < nnodes_extended && isGood; i++ )
   {
      const int start = csr_start[i];
      const int end = csr_start[i + 1];

      if( !nodes_wakeState[i] )
         continue;

      for( int j = start; j != end; j++ )
      {
         const int head = csr_headarr[j];

         if( nodes_wakeState[head] == 0 && edges_capa[j] > 0 )
         {
            assert(n <= maxsepasize);
            assert(edges_capa[j] == 1 || edges_capa[j] == capa_inf);

            if( edges_capa[j] == capa_inf || n == maxsepasize )
            {
               isGood = FALSE;
               break;
            }

            cutterms[n++] = i;
            assert(Is_term(g->term[i]));
         }
      }
   }

   *cutIsGood = isGood;
   *ncutterms = n;
}



/** helper; traverses and does additional optional stuff */
static
SCIP_RETCODE termsepaTraverseSinkComp(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   SCIP_Bool             removeTerms,        /**< remove terminals reachable from sink via non-terminal paths? */
   SCIP_Bool             updateVisitedTerms, /**< do update? */
   int                   ncutterms,          /**< size of the separator */
   const int*            cutterms,           /**< the separator nodes (all terminals) */
   int                   sinkterm,           /**< sink terminal of current cut */
   MINCUT*               mincut,             /**< minimum cut */
   int*                  ncompnodes          /**< number of nodes in component (OUT) */
)
{
   int* RESTRICT termcands = mincut->terms;
   STP_Vectype(int) queue = NULL;
   SCIP_Bool* nodes_isVisited;
   const int nnodes = graph_get_nNodes(g);

   assert(ncutterms > 0);
   assert(mincut->nodes_wakeState[sinkterm] == 0);

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &nodes_isVisited, nnodes) );

   for( int i = 0; i < ncutterms; i++ )
   {
      assert(graph_knot_isInRange(g, cutterms[i]));
      assert(!nodes_isVisited[cutterms[i]]);

      nodes_isVisited[cutterms[i]] = TRUE;
   }

   assert(!nodes_isVisited[mincut->root]);

   /* BFS from sink terminal */
   StpVecReserve(scip, queue, nnodes);
   StpVecPushBack(scip, queue, sinkterm);
   nodes_isVisited[sinkterm] = TRUE;

   /* BFS loop */
   for( int i = 0; i < StpVecGetSize(queue); i++ )
   {
      const int node = queue[i];
      assert(nodes_isVisited[node]);

      if( removeTerms && Is_term(g->term[node]) && node != sinkterm )
      {
         const int ntermcands = mincut->ntermcands;

         assert(mincut->nodes_wakeState[node] == 0);

         /* todo more efficiently */
         for( int t = 0; t < ntermcands; t++ )
         {
            if( node == termcands[t] )
            {
              // printf("removing terminal %d \n", node);

               SWAP_INTS(termcands[mincut->ntermcands - 1], termcands[t]);
               mincut->ntermcands--;
               break;
            }
         }

         continue;
      }

      for( int e = g->outbeg[node]; e >= 0; e = g->oeat[e] )
      {
         const int head = g->head[e];
         if( !nodes_isVisited[head] )
         {
            nodes_isVisited[head] = TRUE;
            StpVecPushBack(scip, queue, head);
         }
      }
   }

   if( ncompnodes )
      *ncompnodes = StpVecGetSize(queue);

   for( int i = 0; i < StpVecGetSize(queue); i++ )
   {
      const int node = queue[i];
      assert(nodes_isVisited[node]);
      nodes_isVisited[node] = FALSE;

      if( updateVisitedTerms && Is_term(g->term[node]) && node != sinkterm )
      {
         assert(mincut->terms_minsepasize && mincut->terms_mincompsize);

         if( mincut->terms_minsepasize[node] > ncutterms )
            mincut->terms_minsepasize[node] = ncutterms;

         if( mincut->terms_mincompsize[node] > *ncompnodes )
            mincut->terms_mincompsize[node] = *ncompnodes;
      }
   }

   for( int i = 0; i < ncutterms; i++ )
   {
      assert(nodes_isVisited[cutterms[i]]);
      nodes_isVisited[cutterms[i]] = FALSE;
   }

   StpVecFree(scip, queue);

#ifndef NDEBUG
   for( int i = 0; i < nnodes; i++ )
   {
      assert(nodes_isVisited[i] == FALSE);
   }
#endif

   SCIPfreeCleanBufferArray(scip, &nodes_isVisited);

   return SCIP_OKAY;
}


/** removes terminals in current cut */
static
SCIP_RETCODE termsepaRemoveCutTerminals(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   int                   ncutterms,          /**< size of the separator */
   const int*            cutterms,           /**< the separator nodes (all terminals) */
   int                   sinkterm,           /**< sink terminal of current cut */
   MINCUT*               mincut              /**< minimum cut */
)
{
   const SCIP_Bool removeTerms = TRUE;
   const SCIP_Bool updateVisitedTerms = FALSE;

   SCIP_CALL( termsepaTraverseSinkComp(scip, g, removeTerms, updateVisitedTerms, ncutterms, cutterms, sinkterm, mincut, NULL) );

   return SCIP_OKAY;
}


/** gets number of nodes */
static
SCIP_RETCODE termsepaGetCompNnodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   int                   ncutterms,          /**< size of the separator */
   const int*            cutterms,           /**< the separator nodes (all terminals) */
   int                   sinkterm,           /**< sink terminal of current cut */
   MINCUT*               mincut,             /**< minimum cut */
   int*                  ncompnodes          /**< number of nodes in component */
)
{
   const SCIP_Bool removeTerms = FALSE;
   const SCIP_Bool updateVisitedTerms = TRUE;

   SCIP_CALL( termsepaTraverseSinkComp(scip, g, removeTerms, updateVisitedTerms, ncutterms, cutterms, sinkterm, mincut, ncompnodes) );

   return SCIP_OKAY;
}


/** stores cut
 *  NOTE: this methods is call once the cut vertices are already stored in the CSR array */
static
SCIP_RETCODE termsepaStoreCutFinalize(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   int                   sinkterm,           /**< sink terminal */
   MINCUT*               mincut,             /**< minimum cut */
   int                   ncutterms,          /**< number of cut terminals */
   const int*            cutterms,           /**<  cut terminals */
   TERMSEPAS*            termsepas,          /**< terminal separator storage */
   SCIP_Bool*            success             /**< success? */
)
{
   TSEPA* sepas = termsepas->sepas;
   int nsinknodes = 0;
   const int nsepas_all = termsepas->nsepas_all;

   assert(0 <= nsepas_all && nsepas_all < termsepas->maxnsepas);
   assert(0 <= ncutterms && ncutterms <= termsepas->maxsepasize);
   assert(termsepas->nsepaterms_csr + ncutterms <= termsepas->maxnsepas * termsepas->maxsepasize);

   *success = TRUE;

   SCIP_CALL( termsepaGetCompNnodes(scip, g, ncutterms, cutterms, sinkterm, mincut, &nsinknodes) );

   /* NOTE: if the sink is already contained in a separated component and does not make anything smaller, we don't want to take it */
   if( ncutterms >= mincut->terms_minsepasize[sinkterm] && nsinknodes >= mincut->terms_mincompsize[sinkterm] )
   {
      SCIPdebugMessage("discarding already covered component because of SEPA size: %d>=%d, COMP size: %d>=%d \n", ncutterms
           , mincut->terms_minsepasize[sinkterm], nsinknodes, mincut->terms_mincompsize[sinkterm]);

      *success = FALSE;
      return SCIP_OKAY;
   }

   sepas[nsepas_all].sinkterm = sinkterm;
   sepas[nsepas_all].nsinknodes = nsinknodes;

#ifndef NDEBUG
   {
      const int* const nodes_wakeState = mincut->nodes_wakeState;
      int nsinknodes_dbg = 0;
      const int nnodes = graph_get_nNodes(g);
      for( int i = 0; i < nnodes; i++ )
      {
         if( nodes_wakeState[i] == 0 )
            nsinknodes_dbg++;
      }

      assert(nsinknodes_dbg >= nsinknodes);
   }
#endif

   termsepas->sepastarts_csr[nsepas_all + 1] = termsepas->sepastarts_csr[nsepas_all] + ncutterms;
   termsepas->nsepaterms_csr += ncutterms;
   termsepas->nsepas_all++;
   termsepas->nsepas[ncutterms]++;

   assert(termsepas->sepastarts_csr[nsepas_all + 1] == termsepas->nsepaterms_csr);

#ifdef SCIP_DEBUG
   if( ncutterms <= TERMSEPA_MAXCUTSIZE )
   {
      printf("terminal cut of size %d for sink %d \n", ncutterms, sinkterm);
      printf("nsinknodes=%d  \n", nsinknodes);

      for( int i = 0; i < ncutterms; i++ )
      {
         printf("%d \n", cutterms[i]);
      }
   }
#endif

   return SCIP_OKAY;
}


/** tries to add cut */
static
SCIP_RETCODE termsepaStoreCutTry(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   int                   sinkterm,           /**< sink terminal */
   MINCUT*               mincut,             /**< minimum cut */
   TERMSEPAS*            termsepas           /**< terminal separator storage */
)
{
   int* cutterms;
   int ncutterms;
   SCIP_Bool isGoodCut;

   assert(!mincut->isLpcut);
   assert(graph_knot_isInRange(g, sinkterm) && Is_term(g->term[sinkterm]));
   assert(mincut->nodes_wakeState[sinkterm] == 0);

   cutterms = &(termsepas->sepaterms_csr[termsepas->nsepaterms_csr]);
   termsepaCollectCutNodes(g, termsepas, mincut, sinkterm, cutterms, &ncutterms, &isGoodCut);

   if( isGoodCut )
   {
      assert(termsepaCutIsCorrect(scip, g, ncutterms, cutterms, sinkterm, mincut));

      SCIP_CALL( termsepaStoreCutFinalize(scip, g, sinkterm, mincut, ncutterms, cutterms, termsepas, &isGoodCut) );

      if( isGoodCut )
      {
         SCIP_CALL( termsepaRemoveCutTerminals(scip, g, ncutterms, cutterms, sinkterm, mincut) );
      }
   }

   return SCIP_OKAY;
}

/** initializes */
static
SCIP_RETCODE termsepaCsrInit(
   const GRAPH*          g,                  /**< the graph */
   MINCUT*               mincut              /**< minimum cut */
)
{
   int* RESTRICT csr_edgeDefaultToCsr = mincut->csr_edgeDefaultToCsr;
   const int nedges = graph_get_nEdges(g);

   assert(nedges >= 2);

   for( int e = 0; e < nedges; e++ )
      csr_edgeDefaultToCsr[e] = -1;

   return SCIP_OKAY;
}


/** creates copies of potential separator terminals */
static
void termsepaCsrAddTermCopies(
   const GRAPH*          g,                  /**< the graph */
   MINCUT*               mincut              /**< minimum cut */
)
{
   int* RESTRICT excess = g->mincut_e;
   int* RESTRICT nodes_termToCopy = mincut->termsepa_termToCopy;
   int* RESTRICT terms = mincut->terms;
   const int* const nodes_wakeState = mincut->nodes_wakeState;
   const int root = mincut->root;
   const int nnodes_org = graph_get_nNodes(g);
   int nsepaterms = 0;
   int ntermcands = 0;

   assert(nodes_wakeState && nodes_termToCopy);

   for( int k = 0; k < nnodes_org; k++ )
   {
      SCIP_Bool isSeparator;
      assert(nodes_termToCopy[k] == -1);

      if( !Is_term(g->term[k]) || k == root )
      {
         continue;
      }

      isSeparator = FALSE;

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];

         /* head not in root cut? */
         if( nodes_wakeState[head] == 0 )
         {
            isSeparator = TRUE;
            break;
         }
      }

      if( nodes_wakeState[k] == 0 )
      {
         SCIPdebugMessage("add terminal candidate %d \n", k);
         terms[ntermcands++] = k;
      }

      if( isSeparator )
      {
         nodes_termToCopy[k] = nnodes_org + nsepaterms;
         nsepaterms++;

         assert(excess[nodes_termToCopy[k]] == 0);

         if( nodes_wakeState[k] != 0 )
         {
            /* NOTE: 1 is the default value for separator edges...so we basically push everything out of k */
            excess[nodes_termToCopy[k]] = 1;
         }

         SCIPdebugMessage("adding separator terminal %d (copyindex=%d) \n", k, nodes_termToCopy[k]);

      }
   }

   // todo check on large test set ... does not seem to help so far
#ifdef SCIP_DISABLED
   if( mincut->randnumgen  )
   {
      SCIPrandomPermuteIntArray(mincut->randnumgen, terms, 0, ntermcands);
      printf("PERMUTE \n\n \n");
   }
#endif

   assert(nnodes_org + nsepaterms <= mincut->termsepa_nnodes);

   mincut->termsepa_nnodes = nnodes_org + nsepaterms;
   mincut->ntermcands = ntermcands;
}


/** adds CSR edges */
static
void termsepaCsrAddEdges(
   const GRAPH*          g,                  /**< the graph */
   MINCUT*               mincut              /**< minimum cut */
)
{
   int* RESTRICT edges_capa = mincut->edges_capa;
   int* RESTRICT residual = g->mincut_r;
   int* RESTRICT edgecurr = g->mincut_numb;
   int* RESTRICT csr_edgeDefaultToCsr = mincut->csr_edgeDefaultToCsr;
   int* RESTRICT csr_start = mincut->csr_start;
   int* RESTRICT csr_headarr = mincut->csr_headarr;
   int* RESTRICT nodes_wakeState = mincut->nodes_wakeState;
   const int* const nodes_termToCopy = mincut->termsepa_termToCopy;
   const int capa_infinity = termsepaGetCapaInf(g, mincut);
   const int capa_one = 1;
   const int nnodes_org = graph_get_nNodes(g);
   int csr_nedges = 0;

   assert(residual && edgecurr && edges_capa);

   /* add edges from original nodes */
   for( int k = 0; k < nnodes_org; k++ )
   {
      const SCIP_Bool kIsSepaTerm = nodes_termToCopy[k] >= 0;
      edgecurr[k] = -1;
      csr_start[k] = csr_nedges;

      /* add edge from terminal to copy if existent */
      if( kIsSepaTerm )
      {
         const int kCopy = nodes_termToCopy[k];

         assert(Is_term(g->term[k]) && k != mincut->root);
         assert(nnodes_org <= kCopy && kCopy < mincut->termsepa_nnodes);

         residual[csr_nedges] = capa_one;
         csr_headarr[csr_nedges++] = kCopy;
      }

      /* non-dormant node or potential separator? */
      if( kIsSepaTerm || nodes_wakeState[k] == 0 )
      {
         assert(k != mincut->root);
         edgecurr[k] = csr_nedges;

         /* NOTE: we go two times over the incident edges so that the edges to
          * the terminal copies are at the start */
         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];

            /* is head a separator terminal? */
            if( nodes_termToCopy[head] >= 0 )
            {
               const int headCopy = nodes_termToCopy[head];

               assert(Is_term(g->term[head]) && head != mincut->root);
               assert(nnodes_org <= headCopy && headCopy < mincut->termsepa_nnodes);

               residual[csr_nedges] = 0;
               csr_headarr[csr_nedges++] = headCopy;
            }
         }

         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];
            const SCIP_Bool headIsSepaTerm = nodes_termToCopy[head] >= 0;

            if( kIsSepaTerm && headIsSepaTerm )
               continue;

            if( headIsSepaTerm || nodes_wakeState[head] == 0 )
            {
               csr_edgeDefaultToCsr[e] = csr_nedges;
               residual[csr_nedges] = kIsSepaTerm ? 0 : capa_infinity;
               csr_headarr[csr_nedges++] = head;
            }
         }

         /* unreachable node? */
         if( edgecurr[k] == csr_nedges )
         {
            assert(g->grad[k] == 0);
            nodes_wakeState[k] = 1;
         }
      }
   }

   /* add edges from copy terminals */
   for( int k = 0; k < nnodes_org; k++ )
   {
      const SCIP_Bool kIsSepaTerm = nodes_termToCopy[k] >= 0;

      if( kIsSepaTerm )
      {
         const int kCopy = nodes_termToCopy[k];

         assert(nodes_wakeState[kCopy] == 0);
         assert(Is_term(g->term[k]) && k != mincut->root);
         assert(nnodes_org <= kCopy && kCopy < mincut->termsepa_nnodes);

         edgecurr[kCopy] = csr_nedges;
         csr_start[kCopy] = csr_nedges;

         /* edge from copy to k */
         residual[csr_nedges] = 0;
         csr_headarr[csr_nedges++] = k;

         // should not be necessary todo remove
#ifdef SCIP_DISABLED
         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];
            /* is head a separator terminal? */
            if( nodes_termToCopy[head] >= 0 )
            {
               const int headCopy = nodes_termToCopy[head];

               assert(Is_term(g->term[head]) && head != root);
               assert(nnodes_org <= headCopy && headCopy < mincut->termsepa_nnodes);

               residual[csr_nedges] = 0;
               csr_headarr[csr_nedges++] = headCopy;
            }
         }
#endif

         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];
            const SCIP_Bool headIsSepaTerm = nodes_termToCopy[head] >= 0;

            if( headIsSepaTerm || nodes_wakeState[head] == 0 )
            {
               residual[csr_nedges] = capa_infinity;
               csr_headarr[csr_nedges++] = head;
            }
         }
      }
   }

   assert(csr_nedges <= mincut->termsepa_nedges);

   mincut->termsepa_nedges = csr_nedges;
   csr_start[mincut->termsepa_nnodes] = csr_nedges;

   BMScopyMemoryArray(edges_capa, residual, csr_nedges);

   mincut->csr_nedges = csr_nedges;
}


/** adds CSR reverse edges */
static
void termsepaCsrAddReverseEdges(
   const GRAPH*          g,                  /**< the graph */
   MINCUT*               mincut              /**< minimum cut */
)
{
   int* RESTRICT csr_edgeflipped = mincut->csr_edgeflipped;
   const int* const csr_edgeDefaultToCsr = mincut->csr_edgeDefaultToCsr;
   const int* const csr_start = mincut->csr_start;
   const int* const csr_headarr = mincut->csr_headarr;
   const int nnodes_org = graph_get_nNodes(g);
   const int nedges_org = graph_get_nEdges(g);
   const int nnodes_extended = mincut->termsepa_nnodes;

   /* go over all terminal copies */
   for( int copyterm = nnodes_org; copyterm < nnodes_extended; copyterm++ )
   {
      const int copyedges_start = csr_start[copyterm];
      const int copyedges_end = csr_start[copyterm + 1];

      assert(copyedges_start < copyedges_end);

      for( int copyedge = copyedges_start; copyedge != copyedges_end; copyedge++ )
      {
         int antiedge;
         const int copyneighbor = csr_headarr[copyedge];
         const int neighboredges_start = csr_start[copyneighbor];
         const int neighboredges_end = csr_start[copyneighbor + 1];

         assert(neighboredges_start < neighboredges_end);

         /* NOTE: copy-edges are at the beginning, so should be half-way efficient */
         for( antiedge = neighboredges_start; antiedge != neighboredges_end; antiedge++ )
         {
            if( csr_headarr[antiedge] == copyterm )
            {
               break;
            }
         }

         assert(antiedge < neighboredges_end);
         assert(csr_edgeflipped[antiedge] == -1 || csr_edgeflipped[antiedge] == copyedge);
         assert(csr_edgeflipped[copyedge] == -1 || csr_edgeflipped[copyedge] == antiedge);

         csr_edgeflipped[antiedge] = copyedge;
         csr_edgeflipped[copyedge] = antiedge;
      }
   }

   for( int e = 0; e < nedges_org; e++ )
   {
      if( csr_edgeDefaultToCsr[e] >= 0 )
      {
         const int csr_pos = csr_edgeDefaultToCsr[e];
         assert(csr_edgeflipped[csr_pos] == -1);
         csr_edgeflipped[csr_pos] = csr_edgeDefaultToCsr[flipedge(e)];
      }
   }

   assert(csrFlipedgesAreValid(g, mincut));

}


/** gets maximum number of nodes for extended terminal separation graph  */
static
int termsepaCsrGetMaxNnodes(
   const GRAPH*          g                   /**< the graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int nterms = graph_get_nTerms(g);

   assert(nterms >= 2);

   return nnodes + nterms - 1;
}


/** gets maximum number of edges for extended terminal separation graph */
static
int termsepaCsrGetMaxNedges(
   int                   root,               /**< root of enlarged graph */
   const GRAPH*          g                   /**< the graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   int nedges_enlarged = nedges;

   assert(graph_knot_isInRange(g, root));
   assert(Is_term(g->term[root]));
   assert(nedges >= 2);

   for( int i = 0; i < nnodes; i++ )
   {
      if( !Is_term(g->term[i]) || i == root )
         continue;

      /* NOTE: the + 2 is for the two arcs between the terminal and its copy */
      nedges_enlarged += 2 * g->grad[i] + 2;
   }

   return nedges_enlarged;
}


/** builds and marks root component (nodes reachable from root via non-terminal paths) */
static
void termsepaBuildRootcomp(
   const GRAPH*          g,                  /**< the graph */
   MINCUT*               mincut              /**< minimum cut */
)
{
   int* RESTRICT nodes_wakeState = mincut->nodes_wakeState;
   int* RESTRICT rootcut = mincut->rootcut;
   const int root = mincut->root;
   int rootcutsize = 0;

   nodes_wakeState[root] = 1;
   rootcut[rootcutsize++] = root;

   /* BFS loop, ignoring terminals */
   for( int i = 0; i < rootcutsize; i++ )
   {
      const int k = rootcut[i];

      assert(rootcutsize <= g->knots);
      assert(k < g->knots);

      if( Is_term(g->term[k]) && k != root )
         continue;

      /* traverse outgoing arcs */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];

         /* head not been added to root cut yet? */
         if( nodes_wakeState[head] == 0 )
         {
            nodes_wakeState[head] = 1;
            SCIPdebugMessage("add to root cut: %d \n", head);

            rootcut[rootcutsize++] = head;
         }
      }
   }

   mincut->rootcutsize = rootcutsize;
}


/** builds CSR representation of enlarged graph; also build terminal candidates */
static
SCIP_RETCODE termsepaBuildCsr(
   const GRAPH*          g,                  /**< the graph */
   MINCUT*               mincut              /**< minimum cut */
)
{
   SCIP_CALL( termsepaCsrInit(g, mincut) );

   /* create copies for potential separator terminals */
   termsepaCsrAddTermCopies(g, mincut);

   termsepaCsrAddEdges(g, mincut);

   termsepaCsrAddReverseEdges(g, mincut);

#ifdef SCIP_DEBUG
   //debugPrintCsr(g, mincut);
#endif

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
   int* nodes_wakeState;
   STP_Bool* RESTRICT edges_isRemoved;
   const int nedges = graph_get_nEdges(g);
   const int nnodes = graph_get_nNodes(g);

   assert(mincut && scip);
   assert(mincut->isLpcut);
   assert(!mincut->edges_isRemoved);

   mincut->xval = SCIPprobdataGetXval(scip, NULL);
   assert(mincut->xval);

   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->edges_capa), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->nodes_wakeState), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->terms), g->terms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->csr_edgeDefaultToCsr), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->csr_headarr), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->csr_edgeflipped), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->csr_start), nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->rootcut), nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->edges_isRemoved), nedges) );

   edges_isRemoved = mincut->edges_isRemoved;
   for( int i = 0; i < nedges; i++ )
   {
      edges_isRemoved[i] = (SCIPvarGetUbGlobal(vars[i]) < 0.5);
   }

   nodes_wakeState = mincut->nodes_wakeState;
   for( int k = 0; k < nnodes; k++ )
   {
      nodes_wakeState[k] = 0;
   }

#ifndef NDEBUG
   for( int i = 0; i < nedges; i++ )
   {
      mincut->csr_edgeflipped[i] = -1;
   }
#endif

   return SCIP_OKAY;
}


/** initializes for terminal separation */
static
SCIP_RETCODE mincutInitForTermSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   MINCUT*               mincut              /**< minimum cut */
)
{
   int* RESTRICT nodes_termToCopy;
   int* RESTRICT nodes_wakeState;
   int* RESTRICT terms_sepasize;
   int* RESTRICT terms_mincompsize;
   int nnodes_enlarged;
   int nedges_enlarged;
   const int nedges = graph_get_nEdges(g);
   const int nnodes = graph_get_nNodes(g);

   assert(mincut && scip);
   assert(!mincut->isLpcut);
   assert(!mincut->edges_isRemoved);

   mincut->root = termsepaFindTerminalSource(scip, g, mincut);
   nnodes_enlarged = termsepaCsrGetMaxNnodes(g);
   nedges_enlarged = termsepaCsrGetMaxNedges(mincut->root, g);

   SCIPdebugMessage("selected source %d \n", mincut->root);

   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->edges_capa), nedges_enlarged) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->nodes_wakeState), nnodes_enlarged) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->terms), g->terms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->csr_edgeDefaultToCsr), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->csr_headarr), nedges_enlarged) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->csr_edgeflipped), nedges_enlarged) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->csr_start), nnodes_enlarged + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->rootcut), nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->termsepa_termToCopy), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->terms_minsepasize), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->terms_mincompsize), nnodes) );

   nodes_termToCopy = mincut->termsepa_termToCopy;
   terms_sepasize = mincut->terms_minsepasize;
   terms_mincompsize = mincut->terms_mincompsize;
   for( int i = 0; i < nnodes; i++ )
   {
      nodes_termToCopy[i] = -1;
      terms_sepasize[i] = nnodes;
      terms_mincompsize[i] = nnodes;
   }

   mincut->termsepa_nnodes = nnodes_enlarged;
   mincut->termsepa_nedges = nedges_enlarged;

   SCIP_CALL( graph_mincut_reInit(scip, nnodes_enlarged, nedges_enlarged, g) );

   nodes_wakeState = mincut->nodes_wakeState;
   for( int k = 0; k < nnodes_enlarged; k++ )
   {
      nodes_wakeState[k] = 0;
   }

#ifndef NDEBUG
   for( int i = 0; i < nedges_enlarged; i++ )
   {
      mincut->csr_edgeflipped[i] = -1;
      mincut->edges_capa[i] = -1;
   }
#endif

   return SCIP_OKAY;
}


/** initializes */
static
SCIP_RETCODE mincutInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator or NULL */
   SCIP_Bool             isLpcut,            /**< for LP cut? */
   GRAPH*                g,                  /**< the graph */
   MINCUT**              mincut              /**< minimum cut */
)
{
   MINCUT* mcut;

   assert(scip);
   assert(graph_get_nNodes(g) > 0 && graph_get_nEdges(g) > 0);

   SCIP_CALL( SCIPallocMemory(scip, mincut) );
   mcut = *mincut;

   mcut->xval = NULL;
   mcut->edges_isRemoved = NULL;
   mcut->isLpcut = FALSE;
   mcut->ntermcands = -1;
   mcut->rootcutsize = -1;
   mcut->csr_nedges = -1;
   mcut->root = g->source;
   mcut->termsepa_termToCopy = NULL;
   mcut->terms_minsepasize = NULL;
   mcut->terms_mincompsize = NULL;
   mcut->termsepa_nnodes = -1;
   mcut->termsepa_nedges = -1;
   mcut->randnumgen = randnumgen;
  // mcut->termsepa_inEdges = NULL;
  // mcut->termsepa_inNeighbors = NULL;
  // mcut->termsepa_inEdgesStart = NULL;

   if( isLpcut )
   {
      mcut->isLpcut = TRUE;
      SCIP_CALL( mincutInitForLp(scip, g, mcut) );
   }
   else
   {
      mcut->isLpcut = FALSE;
      SCIP_CALL( mincutInitForTermSepa(scip, g, mcut) );
   }


   return SCIP_OKAY;
}


/** prepares for LP */
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

#ifdef SCIP_DISABLED
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
         edges_capa[e] = 0;
         edges_capa[erev] = 0;
      //   residual[e] = 0;
      //   residual[erev] = 0;

         csr_headarr[e] = 1;
         csr_headarr[erev] = 1;
      }
      else
      {
         edges_capa[e]     = (int)(xval[e] * FLOW_FACTOR + 0.5) + CREEP_VALUE;
         edges_capa[erev]  = (int)(xval[erev] * FLOW_FACTOR + 0.5) + CREEP_VALUE;
     //    residual[e] = edges_capa[e];
     //    residual[erev] = edges_capa[erev];

         /* NOTE: here we misuse csr_headarr to mark the free edges for the BFS */
         csr_headarr[e] = SCIPisFeasGE(scip, xval[e], 1.0) ? 0 : 1;
         csr_headarr[erev] = SCIPisFeasGE(scip, xval[erev], 1.0) ? 0 : 1;
      }
   }
#endif

   for( int e = 0; e < nedges; e++ )
   {
      SCIP_Bool eIsRemoved = edges_isRemoved[e];

      if( intree && !eIsRemoved )
      {
         // todo remove...does not really work together with non-local cuts!
         eIsRemoved = (SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[flipedge(e)]) < 0.5);
      }

      if( eIsRemoved )
      {
         edges_capa[e] = 0;
         csr_headarr[e] = 1;
      }
      else
      {
         edges_capa[e] = (int)(xval[e] * FLOW_FACTOR + 0.5) + CREEP_VALUE;

         /* NOTE: here we misuse csr_headarr to mark the free edges for the BFS */
         csr_headarr[e] = SCIPisFeasGE(scip, xval[e], 1.0) ? 0 : 1;
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
            if( nodes_wakeState[head] == 0 )
            {
               if( edges_capa[e] == 0 && edges_capa[flipedge(e)] == 0 )
                  continue;

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



/** prepares for terminal separator comptation */
static
SCIP_RETCODE mincutPrepareForTermSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   MINCUT*               mincut              /**< minimum cut */
)
{
   assert(!mincut->isLpcut);

   termsepaBuildRootcomp(g, mincut);

   /* build enlarged graph (and also terminal cut candidates) */
   SCIP_CALL( termsepaBuildCsr(g, mincut) );

   return SCIP_OKAY;
}


/** returns next promising terminal for computing cut */
static
int mincutGetNextSinkTerm(
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Bool             firstrun,           /**< first run?  */
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

   if( firstrun )
   {
      if( mincut->randnumgen && ntermcands > 1 )
      {
         const int pos = SCIPrandomGetInt(mincut->randnumgen, 0, ntermcands - 1);
         assert(0 <= pos && pos <= ntermcands - 1);

         SWAP_INTS(termcands[ntermcands - 1], termcands[pos]);
      }

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
   int                   sinkterm,           /**< sink terminal */
   SCIP_Bool             wasRerun,           /**< not the first run? */
   MINCUT*               mincut              /**< minimum cut */
)
{
   const int* capa;
   int nnodes;
   const int root = mincut->root;

   if( mincut->isLpcut )
   {
      capa = mincut->edges_capa;
      nnodes = graph_get_nNodes(g);
      assert(g->source == root);
   }
   else
   {
      /* NOTE: in this case mincut->edges_capa refers to the extended CSR graph,
       * and cannot be handled by graph_mincut_exec....anyway only needed for debug checks*/
      capa = NULL;
      nnodes = mincut->termsepa_nnodes;
   }

   graph_mincut_exec(g, root, sinkterm, nnodes, mincut->csr_nedges, mincut->rootcutsize,
         mincut->rootcut, capa, mincut->nodes_wakeState,
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

   SCIPfreeBufferArrayNull(scip, &(mcut->terms_mincompsize));
   SCIPfreeBufferArrayNull(scip, &(mcut->terms_minsepasize));
   SCIPfreeBufferArrayNull(scip, &(mcut->termsepa_termToCopy));
   SCIPfreeBufferArrayNull(scip, &(mcut->edges_isRemoved));
   SCIPfreeBufferArray(scip, &(mcut->rootcut));
   SCIPfreeBufferArray(scip, &(mcut->csr_start));
   SCIPfreeBufferArray(scip, &(mcut->csr_edgeflipped));
   SCIPfreeBufferArray(scip, &(mcut->csr_headarr));
   SCIPfreeBufferArray(scip, &(mcut->csr_edgeDefaultToCsr));
   SCIPfreeBufferArray(scip, &(mcut->terms));
   SCIPfreeBufferArray(scip, &(mcut->nodes_wakeState));
   SCIPfreeBufferArrayNull(scip, &(mcut->edges_capa));

   SCIPfreeMemory(scip, mincut);
}


#ifdef SCIP_DISABLED
/** checks */
static
SCIP_RETCODE lpcutRemoveReachableTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   int                   sinkterm,
   MINCUT*               mincut              /**< minimum cut */
)
{
   SCIP_Bool* nodes_isVisited;
   STP_Vectype(int) queue = NULL;
   const int* const nodes_wakeState = mincut->nodes_wakeState;
   const SCIP_Real* const xval = mincut->xval;
   const int nnodes = graph_get_nNodes(g);

   assert(xval);
   assert(mincut->isLpcut);

   SCIP_CALL( SCIPallocBufferArray(scip, &nodes_isVisited, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      nodes_isVisited[i] = FALSE;

   printf("check for removables from %d \n", sinkterm);

   /* BFS from sink terminal */
   StpVecReserve(scip, queue, 8);
   StpVecPushBack(scip, queue, sinkterm);
   nodes_isVisited[sinkterm] = TRUE;

   /* BFS loop */
   for( int i = 0; i < StpVecGetSize(queue); i++ )
   {
      const int node = queue[i];

      assert(StpVecGetSize(queue) <= nnodes);
      assert(graph_knot_isInRange(g, node));

      /* traverse outgoing arcs */
      for( int e = g->outbeg[node]; e >= 0; e = g->oeat[e] )
      {
         const int head = g->head[e];

         if( nodes_wakeState[head] == 0 && !nodes_isVisited[head] )
         {
            if( SCIPisFeasGE(scip, xval[e], 1.0) )
            {
               nodes_isVisited[head] = TRUE;
               StpVecPushBack(scip, queue, head);

               if( Is_term(g->term[head]) )
               {
                  printf("found terminal %d \n", head);

                  const int ntermcands = mincut->ntermcands;

                  assert(mincut->nodes_wakeState[node] == 0);

                  /* todo more efficiently */
                  for( int t = 0; t < ntermcands; t++ )
                  {
                     if( head == mincut->terms[t] )
                     {
                        printf("removing terminal %d \n", head);

                        SWAP_INTS(mincut->terms[mincut->ntermcands - 1], mincut->terms[t]);
                        mincut->ntermcands--;
                        break;
                     }
                  }

               }
            }
         }
      }
   }

   StpVecFree(scip, queue);

   SCIPfreeBufferArray(scip, &nodes_isVisited);

   return SCIP_OKAY;
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

   SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "twocut", 1.0, SCIPinfinity(scip), local, FALSE, TRUE) );
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

         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i], 1.0) );
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
         capa[k] = (int) (xval[k] * FLOW_FACTOR + 0.5);
         capa[krev] = (int) (xval[krev] * FLOW_FACTOR + 0.5);
      }
      else
      {
         capa[k] = (int) (xval[krev] * FLOW_FACTOR + 0.5);
         capa[krev] = (int) (xval[k] * FLOW_FACTOR + 0.5);
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


/** initializes */
SCIP_RETCODE mincut_termsepasInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   int                   maxnsepas,          /**< maximum number of separators to compute */
   int                   maxsepasize,        /**< maximum size of separator */
   TERMSEPAS**           termsepas           /**< to initialize */
)
{
   TERMSEPAS* tsepas;
   assert(scip && g);
   assert(maxnsepas >= 1);
   assert(maxsepasize >= 2);

   SCIP_CALL( SCIPallocMemory(scip, termsepas) );
   tsepas = *termsepas;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(tsepas->nsepas), maxsepasize + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(tsepas->currsepa_n), maxsepasize + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(tsepas->sepas), maxnsepas) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(tsepas->sepastarts_csr), maxnsepas + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(tsepas->sepaterms_csr), maxnsepas * maxsepasize + 1) );
   tsepas->maxnsepas = maxnsepas;
   tsepas->maxsepasize = maxsepasize;
   tsepas->nsepas_all = 0;
   tsepas->nsepaterms_csr = 0;
   tsepas->sepastarts_csr[0] = 0;
   tsepas->root = -1;

   for( int i = 0; i <= maxsepasize; i++ )
   {
      tsepas->nsepas[i] = 0;
      tsepas->currsepa_n[i] = -1;
   }

   return SCIP_OKAY;
}


/** frees */
void mincut_termsepasFree(
   SCIP*                 scip,               /**< SCIP */
   TERMSEPAS**           termsepas           /**< to free */
)
{
   TERMSEPAS* tsepas;
   assert(scip && termsepas);

   tsepas = *termsepas;
   assert(tsepas);

   SCIPfreeMemoryArray(scip, &(tsepas->sepaterms_csr));
   SCIPfreeMemoryArray(scip, &(tsepas->sepastarts_csr));
   SCIPfreeMemoryArray(scip, &(tsepas->sepas));
   SCIPfreeMemoryArray(scip, &(tsepas->currsepa_n));
   SCIPfreeMemoryArray(scip, &(tsepas->nsepas));

   SCIPfreeMemory(scip, termsepas);
}


/** returns number of all separators */
int mincut_termsepasGetNall(
   const TERMSEPAS*      termsepas           /**< terminal separators */
)
{
   assert(termsepas);
   assert(termsepas->nsepas_all >= 0);

   return termsepas->nsepas_all;
}


/** returns number of separators per given size */
int mincut_termsepasGetN(
   const TERMSEPAS*      termsepas,          /**< terminal separators */
   int                   sepasize            /**< size */
)
{
   assert(termsepas);
   assert(sepasize >= 1 && sepasize <= termsepas->maxsepasize);
   assert(termsepas->nsepas[sepasize] >= 0);
   assert(termsepas->nsepas[sepasize] <= termsepas->maxnsepas);

   return termsepas->nsepas[sepasize];
}


/** Returns first separator of given size. Returns NULL if none is available. */
const int* mincut_termsepasGetFirst(
   int                   sepasize,           /**< size */
   TERMSEPAS*            termsepas,          /**< terminal separators */
   int*                  sinkterm,           /**< the sink */
   int*                  nsinknodes          /**< number of sink-side nodes */
)
{
   assert(termsepas && sinkterm && nsinknodes);
   assert(sepasize >= 1 && sepasize <= termsepas->maxsepasize);

   termsepas->currsepa_n[sepasize] = -1;

   return mincut_termsepasGetNext(sepasize, termsepas, sinkterm, nsinknodes);
}


/** Returns next separator of given size. Returns NULL if none is available. */
const int* mincut_termsepasGetNext(
   int                   sepasize,           /**< size */
   TERMSEPAS*            termsepas,          /**< terminal separators */
   int*                  sinkterm,           /**< the sink */
   int*                  nsinknodes          /**< number of sink-side nodes */
)
{
   assert(termsepas && sinkterm && nsinknodes);
   assert(sepasize >= 1 && sepasize <= termsepas->maxsepasize);
   assert(termsepas->currsepa_n[sepasize] >= -1);

   *sinkterm = -1;
   *nsinknodes = -1;

   if( termsepas->currsepa_n[sepasize] >= termsepas->nsepas_all )
   {
      assert(termsepas->currsepa_n[sepasize] == termsepas->nsepas_all);
   }
   else
   {
      const int* const starts = termsepas->sepastarts_csr;
      const int* const terms = termsepas->sepaterms_csr;
      const int startpos = termsepas->currsepa_n[sepasize] + 1;
      int s;

      assert(0 <= startpos && startpos <= termsepas->nsepas_all);

      for( s = startpos; s < termsepas->nsepas_all; s++ )
      {
         const int size = starts[s + 1] - starts[s];
         assert(sepasize >= 1 && sepasize <= termsepas->maxsepasize);

         if( size == sepasize )
            break;
      }

      termsepas->currsepa_n[sepasize] = s;

      if( s < termsepas->nsepas_all )
      {
         TSEPA tsepa = termsepas->sepas[s];
         *sinkterm = tsepa.sinkterm;
         *nsinknodes = tsepa.nsinknodes;

         return &(terms[starts[s]]);
      }
   }

   return NULL;
}


/** gets terminal separator source */
int mincut_termsepasGetSource(
   const TERMSEPAS*      termsepas           /**< terminal separators */
)
{
   assert(termsepas);

   return termsepas->root;
}


/** is it promising to look for terminal separators? */
SCIP_Bool mincut_findTerminalSeparatorsIsPromising(
   const GRAPH*          g                   /**< graph data structure */
   )
{
   int nnodes;
   int nedges;

   assert(g);

   if( g->terms < 4 )
   {
      return FALSE;
   }

   graph_get_nVET(g, &nnodes, &nedges, NULL);

   if( nedges == 0 )
   {
      return FALSE;
   }

   assert(nnodes > 0 && nedges >= 2);
   assert(nedges % 2 == 0);

   nedges /= 2;

   return ((nedges / nnodes) <= TERMSEPA_SPARSE_MAXRATIO);
}


/** searches for (small) terminal separators */
SCIP_RETCODE mincut_findTerminalSeparators(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator or NULL */
   GRAPH*                g,                  /**< graph data structure */
   TERMSEPAS*            termsepas           /**< terminal separator storage */
   )
{
   MINCUT* mincut;
   int* nodes_wakeState;
   SCIP_Bool wasRerun;

   assert(scip && g && termsepas);
   assert(mincut_termsepasGetNall(termsepas) == 0);

   if( g->terms < 3 )
   {
      SCIPdebugMessage("not enough terminals (%d), aborting terminal separator check \n", g->terms);
      return SCIP_OKAY;
   }

   SCIP_CALL( reduce_unconnected(scip, g) );

   SCIP_CALL( mincutInit(scip, randnumgen, FALSE, g, &mincut) );
   termsepas->root = mincut->root;

   /* sets excess, g->mincut_head, g->mincut_head_inact */
   graph_mincut_setDefaultVals(g);

   SCIP_CALL( mincutPrepareForTermSepa(scip, g, mincut) );
   nodes_wakeState = mincut->nodes_wakeState;
   wasRerun = FALSE;

   assert(nodes_wakeState);
   assert(termsepas->nsepas_all == 0);

#ifdef SCIP_DEBUG
   graph_printInfoReduced(g);
   SCIPdebugMessage("number of terminal separation candidates: %d \n",  mincut->ntermcands );
#endif

   // todo probably want to bound the maximum number of iterations!
   while( mincut->ntermcands > 0 && termsepas->nsepas_all < termsepas->maxnsepas )
   {
      int sinkterm;

      if( ((unsigned) mincut->ntermcands) % 32 == 0 && SCIPisStopped(scip) )
         break;

      /* look for non-reachable terminal */
      sinkterm = mincutGetNextSinkTerm(g, !wasRerun, mincut);
      mincut->ntermcands--;

      SCIPdebugMessage("computing cut for sink terminal %d \n", sinkterm);

      assert(Is_term(g->term[sinkterm]) && mincut->root != sinkterm);

      /* non-trivial cut? */
      if( nodes_wakeState[sinkterm] != 1 )
      {
         mincutExec(g, sinkterm, wasRerun, mincut);
         assert(nodes_wakeState[mincut->root] != 0);

         SCIP_CALL( termsepaStoreCutTry(scip, g, sinkterm, mincut, termsepas) );
#ifdef SCIP_DEBUG
         debugPrintCsrCutEdges(g, mincut);
#endif
      }
      else
      {
         SCIPdebugMessage("cut is trivial \n");

         assert(wasRerun);
      }

      wasRerun = TRUE;
   }

   SCIPdebugMessage("number of separators: %d \n", termsepas->nsepas_all);

   mincutFree(scip, &mincut);

   return SCIP_OKAY;
}


/** separates Steiner cuts for LP */
SCIP_RETCODE mincut_separateLp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator  or NULL   */
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

   SCIP_CALL( mincutInit(scip, randnumgen, TRUE, g, &mincut) );

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

      if( ((unsigned) mincut->ntermcands) % 32 == 0 && SCIPisStopped(scip) )
         break;

      /* look for non-reachable terminal */
      sinkterm = mincutGetNextSinkTerm(g, !wasRerun, mincut);
      mincut->ntermcands--;

      assert(Is_term(g->term[sinkterm]) && g->source != sinkterm);

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
            assert(nodes_wakeState[g->source] != 0);

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
