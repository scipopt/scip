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

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


//#define SCIP_DEBUG


#include "mincut.h"
#include "probdata_stp.h"
#include "misc_stp.h"
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
   int*                  termsepa_termToCopy;
   STP_Bool*             edges_isRemoved;    /**< only used for LP cuts */
   int                   ntermcands;
   int                   rootcutsize;
   int                   csr_nedges;
   int                   root;
   int                   termsepa_nnodes;
   int                   termsepa_nedges;
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
            printf("%d->%d \n", i, head);
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
      }
   }

   return TRUE;
}
#endif

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
         /* NOTE: 1 is the default value for separator edges...so we basically push everything out of k */
         excess[nodes_termToCopy[k]] = 1;

         SCIPdebugMessage("adding separator terminal %d (copyindex=%d) \n", k, nodes_termToCopy[k]);

      }
   }

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
   int* RESTRICT residual = g->mincut_r;
   int* RESTRICT edgecurr = g->mincut_numb;
   int* RESTRICT csr_edgeDefaultToCsr = mincut->csr_edgeDefaultToCsr;
   int* RESTRICT csr_start = mincut->csr_start;
   int* RESTRICT csr_headarr = mincut->csr_headarr;
   int* RESTRICT nodes_wakeState = mincut->nodes_wakeState;
   const int* const nodes_termToCopy = mincut->termsepa_termToCopy;
   const int capa_infinity = g->terms;
   const int capa_one = 1;
   const int root = mincut->root;
   const int nnodes_org = graph_get_nNodes(g);
   int csr_nedges = 0;

   assert(residual && edgecurr);

   /* add edges from original nodes */
   for( int k = 0; k < nnodes_org; k++ )
   {
      const SCIP_Bool kIsRealTerm = (Is_term(g->term[k]) && k != root);
      edgecurr[k] = -1;
      csr_start[k] = csr_nedges;

      /* add edge from terminal to copy if existent */
      if( nodes_termToCopy[k] >= 0 )
      {
         const int kCopy = nodes_termToCopy[k];

         assert(Is_term(g->term[k]) && k != root);
         assert(nnodes_org <= kCopy && kCopy < mincut->termsepa_nnodes);

         residual[csr_nedges] = capa_one;
         csr_headarr[csr_nedges++] = kCopy;
      }

      /* non-dormant node? */
      if( nodes_wakeState[k] == 0 )
      {
         assert(k != root);
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

               assert(Is_term(g->term[head]) && head != root);
               assert(nnodes_org <= headCopy && headCopy < mincut->termsepa_nnodes);

               residual[csr_nedges] = 0;
               csr_headarr[csr_nedges++] = headCopy;
            }
         }

         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];
            if( nodes_wakeState[head] == 0 )
            {
               csr_edgeDefaultToCsr[e] = csr_nedges;
               residual[csr_nedges] = kIsRealTerm ? 0 : capa_infinity;
               csr_headarr[csr_nedges++] = head;
            }
         }

         /* unreachable node? */
         if( edgecurr[k] == csr_nedges )
         {
            assert(0 && "should not happen");
         //   nodes_wakeState[k] = 1;
         }
      }
   }

   /* add edges from copy terminals */
   for( int k = 0; k < nnodes_org; k++ )
   {
      if( nodes_termToCopy[k] >= 0 )
      {
         const int kCopy = nodes_termToCopy[k];

         assert(nodes_wakeState[kCopy] == 0);
         assert(Is_term(g->term[k]) && k != root);
         assert(nnodes_org <= kCopy && kCopy < mincut->termsepa_nnodes);

         edgecurr[kCopy] = csr_nedges;
         csr_start[kCopy] = csr_nedges;

         /* edge from copy to k */
         residual[csr_nedges] = 0;
         csr_headarr[csr_nedges++] = k;

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

         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];
            if( nodes_wakeState[head] == 0 )
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
         assert(csr_edgeflipped[antiedge] == -1);
         assert(csr_edgeflipped[copyedge] == -1);

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
int termsepaGetMaxNnodes(
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
int termsepaGetMaxNedges(
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
   debugPrintCsr(g, mincut);
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
   const int nnodes_enlarged = termsepaGetMaxNnodes(g);
   const int nedges_enlarged = termsepaGetMaxNedges(mincut->root, g);
   const int nedges = graph_get_nEdges(g);
   const int nnodes = graph_get_nNodes(g);

   assert(mincut && scip);
   assert(!mincut->isLpcut);
   assert(!mincut->edges_isRemoved);

   mincut->edges_capa = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->nodes_wakeState), nnodes_enlarged) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->terms), g->terms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->csr_edgeDefaultToCsr), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->csr_headarr), nedges_enlarged) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->csr_edgeflipped), nedges_enlarged) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->csr_start), nnodes_enlarged + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->rootcut), nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(mincut->termsepa_termToCopy), nnodes) );

   nodes_termToCopy = mincut->termsepa_termToCopy;
   for( int i = 0; i < nnodes; i++ )
   {
      nodes_termToCopy[i] = -1;
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
   }
#endif

   return SCIP_OKAY;
}


/** initializes */
static
SCIP_RETCODE mincutInit(
   SCIP*                 scip,               /**< SCIP data structure */
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
   mcut->termsepa_nnodes = -1;
   mcut->termsepa_nedges = -1;
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
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator  or NULL   */
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
      if( randnumgen && ntermcands > 1 )
      {
         const int pos = SCIPrandomGetInt(randnumgen, 0, ntermcands - 1);
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
   int                   sinkterm,
   SCIP_Bool             wasRerun,
   MINCUT*               mincut              /**< minimum cut */
)
{
   int nnodes;
   const int root = mincut->root;

   if( mincut->isLpcut )
   {
      nnodes = graph_get_nNodes(g);
   }
   else
   {
      nnodes = mincut->termsepa_nnodes;
   }

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


/** searches for (small) terminal separators */
SCIP_RETCODE mincut_findTerminalSeparators(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
   )
{
   MINCUT* mincut;
   int* nodes_wakeState;
   SCIP_Bool wasRerun;

   SCIP_CALL( mincutInit(scip, FALSE, g, &mincut) );

   /* sets excess, g->mincut_head,  g->mincut_head_inact */
   graph_mincut_setDefaultVals(g);

   SCIP_CALL( mincutPrepareForTermSepa(scip, g, mincut) );
   nodes_wakeState = mincut->nodes_wakeState;
   wasRerun = FALSE;

   assert(nodes_wakeState);

   printf("ntermcands=%d \n",  mincut->ntermcands );

   while( mincut->ntermcands > 0 )
   {
      int sinkterm;

      if( ((unsigned) mincut->ntermcands) % 32 == 0 && SCIPisStopped(scip) )
         break;

      /* look for non-reachable terminal */
      sinkterm = mincutGetNextSinkTerm(g, !wasRerun, NULL, mincut);
      mincut->ntermcands--;

      printf("computing cut for sink terminal %d \n", sinkterm);

      assert(Is_term(g->term[sinkterm]) && g->source != sinkterm);

       /* non-trivial cut? */
       if( nodes_wakeState[sinkterm] != 1 )
       {
          mincutExec(g, sinkterm, wasRerun, mincut);
          assert(nodes_wakeState[g->source] != 0);
#ifdef SCIP_DEBUG
          debugPrintCsrCutEdges(g, mincut);
#endif
       }
       else
       {
          printf("cut is trivial \n");

          assert(wasRerun);
       }

       wasRerun = TRUE;
   }

   mincutFree(scip, &mincut);

   return SCIP_OKAY;
}


/*
 * Interface methods
 */


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

   SCIP_CALL( mincutInit(scip, TRUE, g, &mincut) );

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
      sinkterm = mincutGetNextSinkTerm(g, !wasRerun, randnumgen, mincut);
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
