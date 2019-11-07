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

/**@file   heur_local.c
 * @brief  Improvement heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements several local heuristics, including vertex insertion, key-path exchange and key-vertex elimination,
 * ("Fast Local Search for Steiner Trees in Graphs" by Uchoa and Werneck). Other heuristics are for PCSTP and MWCSP.
 *
 * A list of all interface methods can be found in heur_local.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_local.h"
#include "heur_tm.h"
#include "probdata_stp.h"
#include "cons_stp.h"


/* @note if heuristic is running in root node timing is changed there to (SCIP_HEURTIMING_DURINGLPLOOP |
 *       SCIP_HEURTIMING_BEFORENODE), see SCIP_DECL_HEURINITSOL callback
 */

#define HEUR_NAME             "local"
#define HEUR_DESC             "improvement heuristic for STP"
#define HEUR_DISPCHAR         '-'
#define HEUR_PRIORITY         100
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)

#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_DURINGROOT    TRUE
#define DEFAULT_MAXFREQLOC    FALSE
#define DEFAULT_MAXNBESTSOLS  30
#define DEFAULT_NBESTSOLS     15
#define DEFAULT_MINNBESTSOLS  10
#define LOCAL_MAXRESTARTS  6

#define GREEDY_MAXRESTARTS  3  /**< Max number of restarts for greedy PC/MW heuristic if improving solution has been found. */
#define GREEDY_EXTENSIONS_MW 6   /**< Number of extensions for greedy MW heuristic. MUST BE HIGHER THAN GREEDY_EXTENSIONS */
#define GREEDY_EXTENSIONS    5  /**< Number of extensions for greedy PC heuristic. */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   nfails;             /**< number of fails */
   int                   maxnsols;           /**< maximal number of best solutions to improve */
   int                   nbestsols;          /**< number of best solutions to improve */
   int*                  lastsolindices;     /**< indices of a number of best solutions already tried */
   SCIP_Bool             maxfreq;            /**< should the heuristic be called with maximum frequency? */
   SCIP_Bool             duringroot;         /**< should the heuristic be called during the root node? */
};


/*
 * Local methods
 */

/** recursive methode for a DFS ordering of graph 'g' */
static
void dfsorder(
   const GRAPH*          graph,
   const int*            edges,
   const int*            node,
   int*                  counter,
   int*                  dfst
   )
{
   int oedge = graph->outbeg[*node];

   while( oedge >= 0 )
   {
      if( edges[oedge] >= 0 )
      {
         dfsorder(graph, edges, &(graph->head[oedge]), counter, dfst);
      }
      oedge = graph->oeat[oedge];
   }

   dfst[*counter] = *node;
   (*counter)++;
}


static
SCIP_Real getNewPrizeNode(
   const GRAPH*          graph,
   const STP_Bool*       steinertree,
   const int*            graphmark,
   int                   node,
   STP_Bool*             prizemark,
   int*                  prizemarklist,
   int*                  prizemarkcount
   )
{
   SCIP_Real prizesum = 0.0;
   assert(graph_pc_isPcMw(graph));

   if( graphmark[node] && !steinertree[node] && Is_pterm(graph->term[node]) && !prizemark[node] )
   {
      prizesum += graph->prize[node];
      prizemark[node] = TRUE;
      prizemarklist[(*prizemarkcount)++] = node;
   }

   return prizesum;
}

static
SCIP_Real getNewPrize(
   const GRAPH*          graph,
   const STP_Bool*       steinertree,
   const int*            graphmark,
   int                   edge,
   STP_Bool*             prizemark,
   int*                  prizemarklist,
   int*                  prizemarkcount
   )
{
   SCIP_Real prizesum = 0.0;

   if( graph_pc_isPcMw(graph) )
   {
      const int mhead = graph->head[edge];
      const int mtail = graph->tail[edge];

      prizesum += getNewPrizeNode(graph, steinertree, graphmark, mhead, prizemark, prizemarklist, prizemarkcount);
      prizesum += getNewPrizeNode(graph, steinertree, graphmark, mtail, prizemark, prizemarklist, prizemarkcount);
   }

   return prizesum;
}


/** computes lowest common ancestors for all pairs {vbase(v), vbase(u)} such that {u,w} is a boundary edge,
 * first call should be with u := root */
static
SCIP_RETCODE lca(
   SCIP*                 scip,
   const GRAPH*          graph,
   int                   u,
   UF*                   uf,
   STP_Bool*             nodesmark,
   int*                  steineredges,
   IDX**                 lcalists,
   PHNODE**              boundpaths,
   int*                  heapsize,
   int*                  vbase
   )
{
   int* uboundpaths; /* boundary-paths (each one represented by its boundary edge) having node 'u' as an endpoint */
   int ancestor;
   int v;
   int i;
   int oedge; /* outgoing edge */
   IDX* curr;
   uf->parent[u] = u;

   for( oedge = graph->outbeg[u]; oedge != EAT_LAST; oedge = graph->oeat[oedge] )
   {
      v = graph->head[oedge];
      if( steineredges[oedge] == 0 )
      {
         SCIP_CALL( lca(scip, graph, v, uf, nodesmark, steineredges, lcalists, boundpaths, heapsize, vbase) );
         SCIPStpunionfindUnion(uf, u, v, FALSE);
         uf->parent[SCIPStpunionfindFind(uf, u)] = u;
      }
   }
   nodesmark[u] = TRUE;

   /* iterate through all boundary-paths having one endpoint in the voronoi region of node 'u' */
   SCIP_CALL( SCIPpairheapBuffarr(scip, boundpaths[u], heapsize[u], &uboundpaths) );
   for( i = 0; i < heapsize[u]; i++ )
   {
      oedge = uboundpaths[i];
      v = vbase[graph->head[oedge]];
      if( nodesmark[v] )
      {
         ancestor = uf->parent[SCIPStpunionfindFind(uf, v)];

         /* if the ancestor of 'u' and 'v' is one of the two, the boundary-edge is already in boundpaths[u] */
         if( ancestor != u && ancestor != v)
         {
            SCIP_CALL( SCIPallocBlockMemory(scip, &curr) );
            curr->index = oedge;
            curr->parent = lcalists[ancestor];
            lcalists[ancestor] = curr;
         }
      }
   }

   /* free the boundary-paths array */
   SCIPfreeBufferArray(scip, &uboundpaths);

   return SCIP_OKAY;
}

/** submethod for local extend */
static
SCIP_RETCODE addToCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const PATH*           path,               /**< shortest data structure array */
   int                   i,                  /**< node */
   int                   greedyextensions,   /**< greedyextensions */
   int*                  nextensions,        /**< nextensions */
   GNODE*                candidates,         /**< candidates */
   SCIP_PQUEUE*          pqueue              /**< pqueue */
   )
{

   assert(!graph_pc_knotIsFixedTerm(graph, i));

   if( *nextensions < greedyextensions )
   {
      candidates[*nextensions].dist = graph->prize[i] - path[i].dist;
      candidates[*nextensions].number = i;

      SCIP_CALL( SCIPpqueueInsert(pqueue, &(candidates[(*nextensions)++])) );
   }
   else
   {
      /* get candidate vertex of minimum value */
      GNODE* min = (GNODE*) SCIPpqueueFirst(pqueue);
      if( SCIPisLT(scip, min->dist, graph->prize[i] - path[i].dist) )
      {
         min = (GNODE*) SCIPpqueueRemove(pqueue);
         min->dist = graph->prize[i] - path[i].dist;
         min->number = i;
         SCIP_CALL( SCIPpqueueInsert(pqueue, min) );
      }
   }

   return SCIP_OKAY;
}

/** checks whether node is crucial, i.e. a terminal or a vertex with degree at least 3 (w.r.t. the steinertree) */
static
STP_Bool nodeIsCrucial(
   const GRAPH* graph,
   int* steineredges,
   int node
   )
{
   assert(graph != NULL);
   assert(steineredges != NULL);

   if( graph->term[node] == -1 )
   {
      int counter = 0;
      int e = graph->outbeg[node];
      while( e >= 0 )
      {
         /* check if the adjacent node is in the ST */
         if( steineredges[e] > -1 || steineredges[flipedge(e)] > -1 )
         {
            counter++;
         }
         e = graph->oeat[e];
      }

      if( counter < 3 )
      {
         return FALSE;
      }
   }

   return TRUE;
}


/** perform local vertex insertion heuristic on given Steiner tree */
static
SCIP_Bool solIsTrivialPcMw(
   const GRAPH*          graph,              /**< graph data structure */
   const int*            solEdges            /**< Steiner tree edges */
)
{
   int e;
   const int root = graph->source;

   assert(graph_pc_isPcMw(graph));

   for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
      if( !Is_term(graph->term[graph->head[e]]) && solEdges[e] )
         break;

   if( e == EAT_LAST )
   {
      SCIPdebugMessage("trivial solution given \n");
      return TRUE;
   }

   return FALSE;
}


/** perform local vertex insertion heuristic on given Steiner tree */
static
void markSolTreeNodes(
   const GRAPH*          graph,              /**< graph data structure */
   const int*            solEdges,           /**< Steiner tree edges */
   NODE*                 linkcutNodes,       /**< Steiner tree nodes */
   STP_Bool*             solNodes            /**< Steiner tree nodes */
   )
{
   const int nnodes = graph->knots;
   const int nedges = graph->edges;

   for( int i = 0; i < nnodes; i++ )
   {
      solNodes[i] = FALSE;
      SCIPlinkcuttreeInit(&linkcutNodes[i]);
   }

   /* create a link-cut tree representing the current Steiner tree */
   for( int e = 0; e < nedges; e++ )
      if( solEdges[e] == CONNECT )
         SCIPlinkcuttreeLink(&linkcutNodes[graph->head[e]], &linkcutNodes[graph->tail[e]], flipedge(e));

   /* mark current Steiner tree nodes */
   for( int e = 0; e < nedges; e++ )
   {
      if( solEdges[e] == CONNECT )
      {
         solNodes[graph->tail[e]] = TRUE;
         solNodes[graph->head[e]] = TRUE;
      }
   }
}


/** perform local vertex insertion heuristic on given Steiner tree */
static
SCIP_RETCODE localVertexInsertion(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   STP_Bool*             solNodes,           /**< Steiner tree nodes */
   NODE*                 linkcutNodes,       /**< Steiner tree nodes */
   int*                  solEdges            /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   )
{
   int* insert = NULL;
   int* adds = NULL;
   int* cuts = NULL;
   int* cuts2 = NULL;
   int* solDegree = NULL;
   int i = 0;
   int newnode = 0;
   int newnverts = 0;
   const int nnodes = graph->knots;
   const int nedges = graph->edges;
   const int root = graph->source;
   const STP_Bool pc = graph_pc_isPc(graph);
   const STP_Bool mw = (graph->stp_type == STP_MWCSP);
   const STP_Bool mwpc = graph_pc_isPcMw(graph);
#ifndef NDEBUG
   const SCIP_Real initialobj = graph_sol_getObj(graph->cost, solEdges, 0.0, nedges);
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &insert, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &adds, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cuts, nnodes) );

   if( mw )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &cuts2, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &solDegree, nnodes) );

      BMSclearMemoryArray(solDegree, nnodes);

      for( int e = 0; e < nedges; e++ )
      {
         if( solEdges[e] == CONNECT )
         {
            solDegree[graph->tail[e]]++;
            solDegree[graph->head[e]]++;
         }
      }
   }

   for( ;; )
   {
      SCIP_Real diff;

      /* if vertex i is not in the current ST and has at least two adjacent nodes, it might be added */
      if( !solNodes[i] && graph->grad[i] > 1 && (!mwpc || !Is_term(graph->term[i])) )
      {
         NODE* v;
         int counter;
         int lastnodeidx;
         int insertcount = 0;

         /* if an outgoing edge of vertex i points to the current ST, SCIPlinkcuttreeLink the edge to a list */
         for( int oedge = graph->outbeg[i]; oedge != EAT_LAST; oedge = graph->oeat[oedge])
            if( solNodes[graph->head[oedge]] && (!mwpc || !Is_term(graph->term[graph->head[oedge]])) )
               insert[insertcount++] = oedge;

         /* if there are less than two edges connecting node i and the current tree, continue */
         if( insertcount <= 1 )
            goto ENDOFLOOP;

         if( mw )
            SCIPlinkcuttreeInit(&linkcutNodes[i]);

         /* the node to insert */
         v = &linkcutNodes[i];

         SCIPlinkcuttreeLink(v, &linkcutNodes[graph->head[insert[0]]], insert[0]);

         lastnodeidx = graph->head[insert[0]];

         if( mw )
         {
            assert(!SCIPisPositive(scip, graph->prize[i]));

            diff = -1.0;
            assert(solDegree != NULL);
            solDegree[lastnodeidx]++;
         }
         else
            diff = graph->cost[v->edge];

         counter = 0;

         /* try to add edges between new vertex and tree */
         for( int k = 1; k < insertcount; k++ )
         {
            NODE* firstnode;
            int firstnodidx;
            SCIPlinkcuttreeEvert(v);

            /* next vertex in the current Steiner tree adjacent to vertex i resp. v (the one being scrutinized for possible insertion) */
            firstnodidx = graph->head[insert[k]];
            firstnode = &linkcutNodes[firstnodidx];

            if( mw )
            {
               NODE* chainfirst;
               NODE* chainlast;
               SCIP_Real minweight;

               assert(solDegree != NULL);

               minweight = SCIPlinkcuttreeFindMinChain(scip, graph->prize, graph->head, solDegree, firstnode, &chainfirst, &chainlast);

               if( SCIPisLT(scip, minweight, graph->prize[i]) )
               {
                  assert(chainfirst != NULL && chainlast != NULL);
                  for( NODE* mynode = chainfirst; mynode != chainlast; mynode = mynode->parent )
                  {
                     int mynodeidx = graph->head[mynode->edge];
                     solNodes[mynodeidx] = FALSE;
                     solDegree[mynodeidx] = 0;
                  }

                  SCIPlinkcuttreeCut(chainfirst);
                  SCIPlinkcuttreeCut(chainlast);

                  SCIPlinkcuttreeLink(v, firstnode, insert[k]);
                  solDegree[graph->head[insert[k]]]++;

                  diff = graph->prize[i] - minweight;
                  break;
               }
            }
            else
            {
               /* if there is an edge with cost greater than that of the current edge... */
               NODE* max = SCIPlinkcuttreeFindMax(scip, graph->cost, firstnode);
               if( SCIPisGT(scip, graph->cost[max->edge], graph->cost[insert[k]]) )
               {
                  diff += graph->cost[insert[k]];
                  diff -= graph->cost[max->edge];
                  cuts[counter] = max->edge;
                  SCIPlinkcuttreeCut(max);
                  SCIPlinkcuttreeLink(v, firstnode, insert[k]);
                  assert(v->edge == insert[k]);
                  adds[counter++] = v->edge;
               }
            }
         }

         if( pc && Is_pterm(graph->term[i]) )
            diff -= graph->prize[i];

         /* if the new tree is more expensive than the old one, restore the latter */
         if( mw )
         {
            if( SCIPisLT(scip, diff, 0.0) )
            {
               assert(solDegree != NULL);

               SCIPlinkcuttreeEvert(v);
               solDegree[lastnodeidx]--;
               SCIPlinkcuttreeCut(&linkcutNodes[graph->head[insert[0]]]);
            }
            else
            {
               solNodes[i] = TRUE;
               newnverts++;
            }
         }
         else
         {
            if( !SCIPisNegative(scip, diff) )
            {
               SCIPlinkcuttreeEvert(v);
               for( int k = counter - 1; k >= 0; k-- )
               {
                  SCIPlinkcuttreeCut(&linkcutNodes[graph->head[adds[k]]]);
                  SCIPlinkcuttreeEvert(&linkcutNodes[graph->tail[cuts[k]]]);
                  SCIPlinkcuttreeLink(&linkcutNodes[graph->tail[cuts[k]]], &linkcutNodes[graph->head[cuts[k]]], cuts[k]);
               }

               /* finally, cut the edge added first (if it had been cut during the insertion process, it would have been restored above) */
               SCIPlinkcuttreeEvert(v);
               SCIPlinkcuttreeCut(&linkcutNodes[graph->head[insert[0]]]);
            }
            else
            {
               SCIPlinkcuttreeEvert(&linkcutNodes[root]);
               adds[counter] = insert[0];
               newnode = i;
               solNodes[i] = TRUE;
               newnverts++;
               SCIPdebugMessage("ADDED VERTEX \n");
            }
         }
      }

      ENDOFLOOP:

      if( i < nnodes - 1 )
         i++;
      else
         i = 0;

      if( newnode == i )
         break;
   }

   /* free buffer memory */
   if( mw )
   {
      SCIPfreeBufferArray(scip, &solDegree);
      SCIPfreeBufferArray(scip, &cuts2);
   }
   SCIPfreeBufferArray(scip, &cuts);
   SCIPfreeBufferArray(scip, &adds);
   SCIPfreeBufferArray(scip, &insert);

   for( int e = 0; e < nedges; e++ )
      solEdges[e] = UNKNOWN;

   if( newnverts > 0  )
   {
      if( mwpc )
         SCIP_CALL( SCIPStpHeurTMPrunePc(scip, graph, graph->cost, solEdges, solNodes) );
      else
         SCIP_CALL( SCIPStpHeurTMPrune(scip, graph, graph->cost, 0, solEdges, solNodes) );

      for( i = 0; i < nnodes; i++ )
         SCIPlinkcuttreeInit(&linkcutNodes[i]);

      /* create a link-cut tree representing the current Steiner tree */
      for( int e = 0; e < nedges; e++ )
      {
         if( solEdges[e] == CONNECT )
         {
            assert(solNodes[graph->tail[e]]);
            assert(solNodes[graph->head[e]]);
            SCIPlinkcuttreeLink(&linkcutNodes[graph->head[e]], &linkcutNodes[graph->tail[e]], flipedge(e));
         }
      }
      SCIPlinkcuttreeEvert(&linkcutNodes[root]);
   }
   else
   {
      SCIPlinkcuttreeEvert(&linkcutNodes[root]);
      for( i = 0; i < nnodes; i++ )
      {
         if( solNodes[i] && linkcutNodes[i].edge != -1 )
            solEdges[flipedge(linkcutNodes[i].edge)] = 0;
      }
   }

#ifndef NDEBUG
   {
      const SCIP_Real newobj = graph_sol_getObj(graph->cost, solEdges, 0.0, nedges);
      SCIPdebugMessage("vertex inclusion obj before/after: %f/%f \n", initialobj, newobj);
      assert(SCIPisLE(scip, newobj, initialobj));
   }
#endif

   return SCIP_OKAY;
}


/** perform local vertex insertion heuristic on given Steiner tree */
static
SCIP_RETCODE localKeyVertexHeuristics(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   STP_Bool*             solNodes,           /**< Steiner tree nodes */
   NODE*                 linkcutNodes,       /**< Steiner tree nodes */
   int*                  solEdges            /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   )
{
   IDX* blists_curr;
   IDX** blists_start;  /* array [1,..,nnodes],
                         * if node i is in the current ST, blists_start[i] points to a linked list of all nodes having i as their base */
   PATH* mst;           /* minimal spanning tree structure */
   PATH* vnoi;
   GRAPH* supergraph;
   IDX** lvledges_start;  /* horizontal edges */
   IDX* lvledges_curr;
   PHNODE** boundpaths;
   UF uf;  /* union-find*/
   SCIP_Real* memdist;
   SCIP_Real kpcost;
   SCIP_Real mstcost;
   SCIP_Real edgecost;
   SCIP_Real kpathcost;
   int* vbase;
   int* state;
   int* kpedges;
   int* kpnodes;
   int* dfstree;
   int* newedges;
   int* memvbase;
   int* heapsize;
   int* boundedges;
   int* meminedges;
   int* supernodes;
   int* supernodesid;
   int* prizemarklist = NULL;
   int* const graphmark = graph->mark;
   int e;
   int i;
   int k;
   int l;
   int node;
   int edge;
   int count;
   int nruns;
   int newedge;
   int oldedge;
   int adjnode;
   int nkpedges;
   int nstnodes;
   int nkpnodes;
   int crucnode;   /* current crucial node*/
   int nresnodes;
   int kptailnode;  /* tail node of the current keypath*/
   int localmoves = 2;
   int newpathend = -1;
   int nsupernodes;
   int nboundedges;
   int rootpathstart;
   int prizemarkcount;
   const int probtype = graph->stp_type;
   const int root = graph->source;
   const int nnodes = graph->knots;
   const int nedges = graph->edges;
   STP_Bool* pinned;
   STP_Bool* scanned;
   STP_Bool* nodesmark;
   STP_Bool* prizemark = NULL;
   const STP_Bool mwpc = graph_pc_isPcMw(graph);
   SCIP_Bool solimproved = FALSE;

#ifndef NDEBUG
   const SCIP_Real initialobj = graph_sol_getObj(graph->cost, solEdges, 0.0, graph->edges);
   SCIP_Real objimprovement = 0.0;
#endif

   /* allocate memory */

   /* memory needed for both Key-Path Elimination and Exchange */
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );


   /* only needed for Key-Path Elimination */
   SCIP_CALL( SCIPallocBufferArray(scip, &newedges, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lvledges_start, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundedges, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &supernodesid, nnodes) );

   /* only needed for Key-Path Exchange */

   /* memory needed for both Key-Path Elimination and Exchange */
   if( mwpc )
   {
      SCIP_CALL(SCIPallocBufferArray(scip, &prizemark, nnodes));
      SCIP_CALL(SCIPallocBufferArray(scip, &prizemarklist, nnodes));
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &scanned, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heapsize, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blists_start, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &memvbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &memdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &meminedges, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundpaths, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pinned, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dfstree, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodesmark, nnodes) );


   for( k = 0; k < nnodes; k++ )
      graphmark[k] = (graph->grad[k] > 0);

   graphmark[root] = TRUE;

   if( mwpc )
      for( k = 0; k < nnodes; k++ )
         prizemark[k] = FALSE;

   /* initialize data structures */
   SCIP_CALL( SCIPStpunionfindInit(scip, &uf, nnodes) );

   /* main loop */
   for( nruns = 0; nruns < LOCAL_MAXRESTARTS && localmoves > 0; nruns++ )
   {
      localmoves = 0;

      BMSclearMemoryArray(blists_start, nnodes);

      /* find a DFS order of the ST nodes */
      nstnodes = 0;
      dfsorder(graph, solEdges, &(root), &nstnodes, dfstree);

      /* compute a voronoi diagram with the ST nodes as bases */
      graph_voronoi(scip, graph, graph->cost, graph->cost, solNodes, vbase, vnoi);

      state = graph->path_state;

      /* initialize data structures  */
      for( k = 0; k < nnodes; k++ )
      {
         assert(state[k] == CONNECT || !graphmark[k]);

         pinned[k] = FALSE;
         scanned[k] = FALSE;
         nodesmark[k] = FALSE;

         /* initialize pairing heaps */
         heapsize[k] = 0;
         boundpaths[k] = NULL;

         lvledges_start[k] = NULL;

         if( !graphmark[k] )
            continue;

         /* link all nodes to their (respective) voronoi base */
         SCIP_CALL( SCIPallocBlockMemory(scip, &blists_curr) );
         blists_curr->index = k;
         blists_curr->parent = blists_start[vbase[k]];
         blists_start[vbase[k]] = blists_curr;
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &supernodes, nstnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &kpnodes, nstnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &kpedges, nstnodes) );

      if( mwpc )
      {
         assert(graph->extended);

         for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
         {
            k = graph->head[e];
            if( Is_term(graph->term[k]) )
            {
               int pterm;

               if( !graph_pc_knotIsFixedTerm(graph, k) )
               {
                  graphmark[k] = FALSE;
                  pterm = graph->head[graph->term2edge[k]];
                  assert(Is_pterm(graph->term[pterm]));

                  pinned[pterm] = TRUE;
               }
            }
         }

         if( !graph_pc_isRootedPcMw(graph) )
            graphmark[root] = FALSE;
      }

      /* for each node, store all of its outgoing boundary-edges in a (respective) heap*/
      for( e = 0; e < nedges; e += 2 )
      {
         if( graph->oeat[e] == EAT_FREE )
            continue;

         node = graph->tail[e];
         adjnode = graph->head[e];
         newedges[e] = UNKNOWN;
         newedges[e + 1] = UNKNOWN;

         /* is edge 'e' a boundary-edge? */
         if( vbase[node] != vbase[adjnode] && graphmark[node] && graphmark[adjnode] )
         {
            edgecost = vnoi[node].dist + graph->cost[e] + vnoi[adjnode].dist;

            /* add the boundary-edge 'e' and its reversed to the corresponding heaps */
            SCIP_CALL( SCIPpairheapInsert(scip, &boundpaths[vbase[node]], e, edgecost, &(heapsize[vbase[node]])) );
            SCIP_CALL( SCIPpairheapInsert(scip, &boundpaths[vbase[adjnode]], flipedge(e), edgecost, &(heapsize[vbase[adjnode]])) );
         }
      }

      /* find LCAs for all edges */
      SCIP_CALL( lca(scip, graph, root, &uf, nodesmark, solEdges, lvledges_start, boundpaths, heapsize, vbase) );

      /* henceforth, the union-find structure will be used on the ST */
      SCIPStpunionfindClear(scip, &uf, nnodes);

      /* henceforth, nodesmark will be used to mark the current supervertices (except for the one representing the root-component) */
      for( i = 0; dfstree[i] != root; i++ )
         nodesmark[dfstree[i]] = FALSE;
      nodesmark[dfstree[i]] = FALSE;

      for( k = 0; k < nnodes; k++ )
         assert(!nodesmark[k]);

      /* main loop visiting all nodes of the current ST in post-order */
      for( i = 0; dfstree[i] != root; i++ )
      {
         crucnode = dfstree[i];
         scanned[crucnode] = TRUE;

         SCIPdebugMessage("iteration %d (crucial node: %d) \n", i, crucnode);

         /*  has the node been temporarily removed from the ST? */
         if( !graphmark[crucnode] )
            continue;

         /* is node 'crucnode' a removable crucial node? (i.e. not pinned or a terminal) */
         if( !pinned[crucnode] && !Is_term(graph->term[crucnode]) && nodeIsCrucial(graph, solEdges, crucnode) )
         {
            for( k = 0; k < nnodes; k++ )
               assert(state[k] == CONNECT || !graphmark[k]);

            /* find all (unique) key-paths starting in node 'crucnode' */
            k = UNKNOWN;
            kpcost = 0.0;
            nkpnodes = 0;
            nkpedges = 0;
            nsupernodes = 0;
            for( edge = graph->outbeg[crucnode]; edge != EAT_LAST; edge = graph->oeat[edge] )
            {
               /* check whether the outgoing edge is in the ST */
               if( (solEdges[edge] > -1 && solNodes[graph->head[edge]]) || (solEdges[flipedge(edge)] > -1 && solNodes[graph->tail[edge]]) )
               {
                  kpcost += graph->cost[edge];

                  /* check whether the current edge leads to the ST root*/
                  if( solEdges[flipedge(edge)] > -1 )
                  {
                     k = flipedge(edge);
                     kpedges[nkpedges++] = k;
                     assert( edge == linkcutNodes[crucnode].edge );
                  }
                  else
                  {
                     kpedges[nkpedges++] = edge;
                     adjnode = graph->head[edge];
                     e = edge;

                     /* move along the key-path until its end (i.e. a crucial or pinned node) is reached */
                     while( !pinned[adjnode] && !nodeIsCrucial(graph, solEdges, adjnode) && solNodes[adjnode] )
                     {
                        /* update the union-find data structure */
                        SCIPStpunionfindUnion(&uf, crucnode, adjnode, FALSE);

                        kpnodes[nkpnodes++] = adjnode;

                        for( e = graph->outbeg[adjnode]; e != EAT_LAST; e = graph->oeat[e] )
                        {
                           if( solEdges[e] > -1 )
                           {
                              kpcost += graph->cost[e];
                              kpedges[nkpedges++] = e;
                              break;
                           }
                        }

                        /* assert that each leaf of the ST is a terminal */


                        if( e == EAT_LAST )
                        {
                           localmoves = 0;

                           goto TERMINATE;
                        }
                        assert(e != EAT_LAST);
                        adjnode = graph->head[e];
                     }
                     /* does the last node on the path belong to a removed component? */
                     if( !solNodes[adjnode] )
                     {
                        kpcost -= graph->cost[e];
                        nkpedges--;
                        adjnode = graph->tail[e];
                        if( adjnode != crucnode )
                        {
                           supernodes[nsupernodes++] = adjnode;
                           nodesmark[adjnode] = TRUE;
                        }
                     }
                     else
                     {
                        supernodes[nsupernodes++] = adjnode;
                        nodesmark[adjnode] = TRUE;
                     }
                  }
               }
            }

            /* traverse the key-path leading to the root-component */
            rootpathstart = nkpnodes;
            if( k != -1 )
            {
               /* begin with the edge starting in the root-component of node 'crucnode' */
               e = k;
               adjnode = graph->tail[e];
               while( !pinned[adjnode] && !nodeIsCrucial(graph, solEdges, adjnode) && solNodes[adjnode] )
               {
                  /* update the union-find data structure */
                  kpnodes[nkpnodes++] = adjnode;

                  for( e = graph->inpbeg[adjnode]; e != EAT_LAST; e = graph->ieat[e] )
                  {
                     if( solEdges[e] > -1 )
                     {
                        assert(solNodes[graph->tail[e]]);
                        kpcost += graph->cost[e];
                        kpedges[nkpedges++] = e;
                        break;
                     }
                  }

                  assert( e != EAT_LAST );
                  adjnode = graph->tail[e];
               }
               supernodes[nsupernodes++] = adjnode;
            }

            /* the last of the key-path nodes to be stored is the current key-node */
            kpnodes[nkpnodes++] = crucnode;

            /* number of reset nodes */
            nresnodes = 0;

            /* reset all nodes (referred to as 'C' henceforth) whose bases are internal nodes of the current key-paths */
            for( k = 0; k < nkpnodes; k++ )
            {
               /* reset all nodes having the current (internal) keypath node as their voronoi base */
               blists_curr = blists_start[kpnodes[k]];
               while( blists_curr != NULL )
               {
                  node = blists_curr->index;
                  assert(graphmark[node]);

                  /* store all relevant data */
                  memvbase[nresnodes] = vbase[node];
                  memdist[nresnodes] =  vnoi[node].dist;
                  meminedges[nresnodes] = vnoi[node].edge;
                  nresnodes++;

                  /* reset data */
                  vbase[node] = UNKNOWN;
                  vnoi[node].dist = FARAWAY;
                  vnoi[node].edge = UNKNOWN;
                  state[node] = UNKNOWN;
                  blists_curr = blists_curr->parent;
               }
            }

            /* add vertical boundary-paths between the child components and the root-component (wrt node 'crucnode') */
            nboundedges = 0;
            for( k = 0; k < nsupernodes - 1; k++ )
            {
               l = supernodes[k];
               edge = UNKNOWN;
               while( boundpaths[l] != NULL )
               {
                  SCIP_CALL( SCIPpairheapDeletemin(scip, &edge, &edgecost, &boundpaths[l], &heapsize[l]) );

                  node = (vbase[graph->head[edge]] == UNKNOWN)? UNKNOWN : SCIPStpunionfindFind(&uf, vbase[graph->head[edge]]);

                  /* check whether edge 'edge' represents a boundary-path having an endpoint in the kth-component and in the root-component respectively */
                  if( node != UNKNOWN && !nodesmark[node] && graphmark[node] )
                  {
                     boundedges[nboundedges++] = edge;
                     SCIP_CALL( SCIPpairheapInsert(scip, &boundpaths[l], edge, edgecost, &heapsize[l]) );
                     break;
                  }
               }
            }

            /* add horizontal boundary-paths (between the  child-components) */
            lvledges_curr = lvledges_start[crucnode];
            while( lvledges_curr != NULL )
            {
               edge = lvledges_curr->index;
               k = vbase[graph->tail[edge]];
               l = vbase[graph->head[edge]];
               node = (l == UNKNOWN)? UNKNOWN : SCIPStpunionfindFind(&uf, l);
               adjnode = (k == UNKNOWN)? UNKNOWN : SCIPStpunionfindFind(&uf, k);

               /* check whether the current boundary-path connects two child components */
               if( node != UNKNOWN && nodesmark[node] && adjnode != UNKNOWN && nodesmark[adjnode] )
               {
                  assert(graphmark[node]);
                  assert(graphmark[adjnode]);
                  boundedges[nboundedges++] = edge;
               }
               lvledges_curr = lvledges_curr->parent;
            }

            /* try to connect the nodes of C (directly) to COMP(C), as a preprocessing for graph_voronoiRepair */
            count = 0;
            for( k = 0; k < nkpnodes; k++ )
            {
               blists_curr = blists_start[kpnodes[k]];
               assert( blists_curr != NULL );
               while( blists_curr != NULL )
               {
                  node = blists_curr->index;

                  /* iterate through all outgoing edges of 'node' */
                  for( edge = graph->inpbeg[node]; edge != EAT_LAST; edge = graph->ieat[edge] )
                  {
                     adjnode = graph->tail[edge];

                     /* check whether the adjacent node is not in C and allows a better voronoi assignment of the current node */
                     if( state[adjnode] == CONNECT && SCIPisGT(scip, vnoi[node].dist, vnoi[adjnode].dist + graph->cost[edge])
                        && graphmark[vbase[adjnode]] && graphmark[adjnode] )
                     {
                        vnoi[node].dist = vnoi[adjnode].dist + graph->cost[edge];
                        vbase[node] = vbase[adjnode];
                        vnoi[node].edge = edge;
                     }
                  }
                  if( vbase[node] != UNKNOWN )
                  {
                     heap_add(graph->path_heap, state, &count, node, vnoi);
                  }
                  blists_curr = blists_curr->parent;
               }
            }

            /* if there are no key-path nodes, something has gone wrong */
            assert(nkpnodes != 0);

            graph_voronoiRepairMult(scip, graph, graph->cost, &count, vbase, boundedges, &nboundedges, nodesmark, &uf, vnoi);

            /* create a supergraph, having the endpoints of the key-paths incident to the current crucial node as (super-) vertices */
            SCIP_CALL( graph_init(scip, &supergraph, nsupernodes, nboundedges * 2, 1) );
            supergraph->stp_type = STP_SPG;

            /* add vertices to the supergraph */
            for( k = 0; k < nsupernodes; k++ )
            {
               supernodesid[supernodes[k]] = k;
               graph_knot_add(supergraph, graph->term[supernodes[k]]);
            }

            /* the (super-) vertex representing the current root-component of the ST */
            k = supernodes[nsupernodes - 1];

            /* add edges to the supergraph */
            for( l = 0; l < nboundedges; l++ )
            {
               edge = boundedges[l];
               node = SCIPStpunionfindFind(&uf, vbase[graph->tail[edge]]);
               adjnode = SCIPStpunionfindFind(&uf, vbase[graph->head[edge]]);

               /* if node 'node' or 'adjnode' belongs to the root-component, take the (temporary) root-component identifier instead */
               node = ((nodesmark[node])? node : k);
               adjnode = ((nodesmark[adjnode])? adjnode : k);

               /* compute the cost of the boundary-path pertaining to the boundary-edge 'edge' */
               edgecost = vnoi[graph->tail[edge]].dist + graph->cost[edge] + vnoi[graph->head[edge]].dist;
               graph_edge_add(scip, supergraph, supernodesid[node], supernodesid[adjnode], edgecost, edgecost);
            }

            /* compute a MST on the supergraph */
            SCIP_CALL( SCIPallocBufferArray(scip, &mst, nsupernodes) );
            SCIP_CALL( graph_path_init(scip, supergraph) );
            graph_path_exec(scip, supergraph, MST_MODE, nsupernodes - 1, supergraph->cost, mst);

            /* compute the cost of the MST */
            mstcost = 0.0;
            prizemarkcount = 0;

            /* compute the cost of the MST */
            for( l = 0; l < nsupernodes - 1; l++ )
            {
               /* compute the edge in the original graph corresponding to the current MST edge */
               if( mst[l].edge % 2  == 0 )
                  edge = boundedges[mst[l].edge / 2 ];
               else
                  edge = flipedge(boundedges[mst[l].edge / 2 ]);

               mstcost += graph->cost[edge];
               mstcost -= getNewPrize(graph, solNodes, graphmark, edge, prizemark, prizemarklist, &prizemarkcount);

               assert( newedges[edge] != crucnode && newedges[flipedge(edge)] != crucnode );

               /* mark the edge (in the original graph) as visited */
               newedges[edge] = crucnode;

               /* traverse along the boundary-path belonging to the boundary-edge 'edge' */
               for( node = graph->tail[edge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
               {
                  e = vnoi[node].edge;

                  /* if edge 'e' and its reversed have not been visited yet */
                  if( newedges[e] != crucnode && newedges[flipedge(e)] != crucnode )
                  {
                     newedges[e] = crucnode;
                     mstcost += graph->cost[e];
                     mstcost -= getNewPrize(graph, solNodes, graphmark, e, prizemark, prizemarklist, &prizemarkcount);
                  }
               }
               for( node = graph->head[edge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
               {
                  e = flipedge(vnoi[node].edge);

                  /* if edge 'e' and its reversed have not been visited yet */
                  if( newedges[vnoi[node].edge] != crucnode && newedges[e] != crucnode )
                  {
                     newedges[e] = crucnode;
                     mstcost += graph->cost[e];
                     mstcost -= getNewPrize(graph, solNodes, graphmark, e, prizemark, prizemarklist, &prizemarkcount);
                  }
               }
            }

            for( int pi = 0; pi < prizemarkcount; pi++ )
               prizemark[prizemarklist[pi]] = FALSE;

            if( SCIPisLT(scip, mstcost, kpcost) )
            {
               localmoves++;
               solimproved = TRUE;

               SCIPdebugMessage("found improving solution in KEY VERTEX ELIMINATION (round: %d) \n ", nruns);

#ifndef NDEBUG
               assert((kpcost - mstcost) >= 0.0);
               objimprovement += (kpcost - mstcost);
#endif

               /* unmark the original edges spanning the supergraph */
               for( e = 0; e < nkpedges; e++ )
               {
                  assert(solEdges[kpedges[e]] != -1);
                  solEdges[kpedges[e]] = -1;
               }

               /* mark all ST nodes except for those belonging to the root-component as forbidden */
               for( k = rootpathstart; k < nkpnodes; k++ )
               {
                  graphmark[kpnodes[k]] = FALSE;
                  solNodes[kpnodes[k]] = FALSE;
               }

               for( k = 0; k < i; k++ )
               {
                  node = SCIPStpunionfindFind(&uf, dfstree[k]);
                  if( nodesmark[node] || node == crucnode )
                  {
                     graphmark[dfstree[k]] = FALSE;
                     solNodes[dfstree[k]] = FALSE;
                  }
               }

               /* add the new edges reconnecting the (super-) components */
               for( l = 0; l < nsupernodes - 1; l++ )
               {
                  if( mst[l].edge % 2  == 0 )
                     edge = boundedges[mst[l].edge / 2 ];
                  else
                     edge = flipedge(boundedges[mst[l].edge / 2 ]);

                  /* change the orientation within the target-component if necessary */
                  if( !nodesmark[vbase[graph->head[edge]]] )
                  {
                     node = vbase[graph->head[edge]];
                     k = SCIPStpunionfindFind(&uf, node);
                     assert(nodesmark[k]);
                     while( node != k )
                     {
                        /* the ST edge pointing towards the root */
                        e = linkcutNodes[node].edge;

                        assert(solEdges[e] == -1 && solEdges[flipedge(e)] != -1 );
                        solEdges[e] = CONNECT;
                        solEdges[flipedge(e)] = UNKNOWN;
                        node = graph->head[e];
                     }
                  }

                  /* is the vbase of the current boundary-edge tail in the root-component? */
                  if( !nodesmark[SCIPStpunionfindFind(&uf, vbase[graph->tail[edge]])] )
                  {

                     solEdges[edge] = CONNECT;

                     for( node = graph->tail[edge], adjnode = graph->head[edge]; node != vbase[node]; adjnode = node, node = graph->tail[vnoi[node].edge] )
                     {
                        graphmark[node] = FALSE;

                        if( solEdges[flipedge(vnoi[node].edge)] == CONNECT )
                           solEdges[flipedge(vnoi[node].edge)] = UNKNOWN;

                        solEdges[vnoi[node].edge] = CONNECT;
                     }

                     assert(!nodesmark[node] && vbase[node] == node);
                     assert( graphmark[node] == TRUE );

                     /* is the pinned node its own component identifier? */
                     if( !Is_term(graph->term[node]) && scanned[node] && !pinned[node] && SCIPStpunionfindFind(&uf, node) == node )
                     {
                        graphmark[graph->head[edge]] = FALSE;
                        oldedge = edge;

                        for( edge = graph->outbeg[node]; edge != EAT_LAST; edge = graph->oeat[edge] )
                        {
                           adjnode = graph->head[edge];
                           /* check whether edge 'edge' leads to an ancestor of terminal 'node' */
                           if( solEdges[edge] == CONNECT && graphmark[adjnode] && solNodes[adjnode]  && SCIPStpunionfindFind(&uf, adjnode) != node )
                           {

                              assert(scanned[adjnode]);
                              /* meld the heaps */
                              SCIPpairheapMeldheaps(scip, &boundpaths[node], &boundpaths[adjnode], &heapsize[node], &heapsize[adjnode]);

                              /* update the union-find data structure */
                              SCIPStpunionfindUnion(&uf, node, adjnode, FALSE);

                              /* move along the key-path until its end (i.e. until a crucial node is reached) */
                              while( !nodeIsCrucial(graph, solEdges, adjnode) && !pinned[adjnode] )
                              {
                                 for( e = graph->outbeg[adjnode]; e != EAT_LAST; e = graph->oeat[e] )
                                 {
                                    if( solEdges[e] != -1 )
                                       break;
                                 }

                                 /* assert that each leaf of the ST is a terminal */
                                 assert( e != EAT_LAST );
                                 adjnode = graph->head[e];
                                 if( !solNodes[adjnode]  )
                                    break;
                                 assert(scanned[adjnode]);
                                 assert(SCIPStpunionfindFind(&uf, adjnode) != node);

                                 /* update the union-find data structure */
                                 SCIPStpunionfindUnion(&uf, node, adjnode, FALSE);

                                 /* meld the heaps */
                                 SCIPpairheapMeldheaps(scip, &boundpaths[node], &boundpaths[adjnode], &heapsize[node], &heapsize[adjnode]);
                              }
                           }
                        }
                        edge = oldedge;
                     }

                     /* mark the start node (lying in the root-component of the ST) of the current boundary-path as pinned,
                      * so that it may not be removed later on */
                     pinned[node] = TRUE;

                     for( node = graph->head[edge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                     {
                        graphmark[node] = FALSE;
                        if( solEdges[vnoi[node].edge] == CONNECT )
                           solEdges[vnoi[node].edge] = -1;

                        solEdges[flipedge(vnoi[node].edge)] = CONNECT;
                     }
                  }
                  else
                  {

                     solEdges[edge] = CONNECT;

                     for( node = graph->tail[edge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                     {
                        graphmark[node] = FALSE;
                        if( solEdges[vnoi[node].edge] != CONNECT && solEdges[flipedge(vnoi[node].edge)] != CONNECT )
                        {
                           solEdges[vnoi[node].edge] = CONNECT;
                        }

                     }

                     for( node = graph->head[edge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                     {
                        graphmark[node] = FALSE;

                        solEdges[flipedge(vnoi[node].edge)] = CONNECT;
                        solEdges[vnoi[node].edge] = UNKNOWN;
                     }
                  }
               }

               for( k = 0; k < nkpnodes; k++ )
               {
                  assert(graphmark[kpnodes[k]] == FALSE);
                  assert(solNodes[kpnodes[k]] == FALSE);
               }
               assert(!graphmark[crucnode]);
            }
            else /* improving solution found */
            {
               /* no improving solution has been found during the move */

               /* meld the heap pertaining to 'crucnode' and all heaps pertaining to descendant key-paths of node 'crucnode' */
               for( k = 0; k < rootpathstart; k++ )
               {
                  SCIPpairheapMeldheaps(scip, &boundpaths[crucnode], &boundpaths[kpnodes[k]], &heapsize[crucnode], &heapsize[kpnodes[k]]);
               }
               for( k = 0; k < nsupernodes - 1; k++ )
               {
                  SCIPpairheapMeldheaps(scip, &boundpaths[crucnode], &boundpaths[supernodes[k]], &heapsize[crucnode], &heapsize[supernodes[k]]);

                  /* update the union-find data structure */
                  SCIPStpunionfindUnion(&uf, crucnode, supernodes[k], FALSE);
               }
            }

            /* free the supergraph and the MST data structure */
            graph_path_exit(scip, supergraph);
            graph_free(scip, &supergraph, TRUE);
            SCIPfreeBufferArray(scip, &mst);

            /* unmark the descendant supervertices */
            for( k = 0; k < nsupernodes - 1; k++ )
            {
               nodesmark[supernodes[k]] = FALSE;
            }

            /* debug test; to be deleted later on */
            for( k = 0; k < nnodes; k++ )
            {
               assert( !nodesmark[k] );
            }

            /* restore the original voronoi diagram */
            l = 0;
            for( k = 0; k < nkpnodes; k++ )
            {
               /* restore data of all nodes having the current (internal) key-path node as their voronoi base */
               blists_curr = blists_start[kpnodes[k]];
               while( blists_curr != NULL )
               {
                  node = blists_curr->index;
                  vbase[node] = memvbase[l];
                  vnoi[node].dist = memdist[l];
                  vnoi[node].edge = meminedges[l];
                  l++;
                  blists_curr = blists_curr->parent;
               }
            }

            assert(l == nresnodes);
         }

         /** Key-Path Exchange */
         if( probtype != STP_MWCSP )
         {
            /* if the node has just been eliminated, skip Key-Path Exchange */
            if( !graphmark[crucnode] )
               continue;

            /* is crucnode a crucial or pinned vertex? */
            if( (!nodeIsCrucial(graph, solEdges, crucnode) && !pinned[crucnode]) )
               continue;

            if( Is_term(graph->term[crucnode]) || pinned[crucnode] )
            {
               for( edge = graph->outbeg[crucnode]; edge != EAT_LAST; edge = graph->oeat[edge] )
               {
                  adjnode = graph->head[edge];
                  /* check whether edge 'edge' leads to an ancestor of terminal 'crucnode' */
                  if( solEdges[edge] == CONNECT && solNodes[adjnode] && graphmark[adjnode] )
                  {
                     assert( SCIPStpunionfindFind(&uf, adjnode) != crucnode);
                     assert(scanned[adjnode]);
                     /* meld the heaps */
                     SCIPpairheapMeldheaps(scip, &boundpaths[crucnode], &boundpaths[adjnode], &heapsize[crucnode], &heapsize[adjnode]);

                     /* update the union-find data structure */
                     SCIPStpunionfindUnion(&uf, crucnode, adjnode, FALSE);

                     /* move along the key-path until its end (i.e. until a crucial node is reached) */
                     while( !nodeIsCrucial(graph, solEdges, adjnode) && !pinned[adjnode] )
                     {
                        for( e = graph->outbeg[adjnode]; e != EAT_LAST; e = graph->oeat[e] )
                        {
                           if( solEdges[e] != -1 )
                              break;
                        }

                        /* assert that each leaf of the ST is a terminal */
                        assert( e != EAT_LAST );
                        adjnode = graph->head[e];
                        if( !solNodes[adjnode] || !graphmark[adjnode] )
                           break;
                        assert(scanned[adjnode]);
                        assert(SCIPStpunionfindFind(&uf, adjnode) != crucnode);

                        /* update the union-find data structure */
                        SCIPStpunionfindUnion(&uf, crucnode, adjnode, FALSE);

                        /* meld the heaps */
                        SCIPpairheapMeldheaps(scip, &boundpaths[crucnode], &boundpaths[adjnode], &heapsize[crucnode], &heapsize[adjnode]);
                     }
                  }
               }

            }

            /* counts the internal nodes of the keypath */
            nkpnodes = 0;

            for( k = 0; k < nnodes; k++ )
               assert(state[k] == CONNECT || !graphmark[k]);

            /* find the (unique) key-path containing the parent of the current crucial node 'crucnode' */
            kptailnode = graph->head[linkcutNodes[crucnode].edge];
            kpathcost = graph->cost[linkcutNodes[crucnode].edge];

            while( !nodeIsCrucial(graph, solEdges, kptailnode) && !pinned[kptailnode] )
            {
               kpathcost += graph->cost[linkcutNodes[kptailnode].edge];

               kpnodes[nkpnodes++] = kptailnode;
               kptailnode = graph->head[linkcutNodes[kptailnode].edge];
            }

            /* counts the reset nodes during voronoi repair */
            nresnodes = 0;

            /* reset all nodes (henceforth referred to as 'C') whose bases are internal nodes of the current keypath */
            for( k = 0; k < nkpnodes; k++ )
            {
               /* reset all nodes having the current (internal) keypath node as their voronoi base */
               blists_curr = blists_start[kpnodes[k]];
               while( blists_curr != NULL )
               {
                  node = blists_curr->index;
                  memvbase[nresnodes] = vbase[node];
                  memdist[nresnodes] =  vnoi[node].dist;
                  meminedges[nresnodes] = vnoi[node].edge;
                  nresnodes++;
                  vbase[node] = UNKNOWN;
                  vnoi[node].dist = FARAWAY;
                  vnoi[node].edge = UNKNOWN;
                  state[node] = UNKNOWN;
                  blists_curr = blists_curr->parent;
               }
            }

            edgecost = UNKNOWN;
            e = UNKNOWN;
            while( boundpaths[crucnode] != NULL )
            {
               SCIP_CALL( SCIPpairheapDeletemin(scip, &e, &edgecost, &boundpaths[crucnode], &(heapsize[crucnode])) );
               assert( e != UNKNOWN );
               k = vbase[graph->tail[e]];
               l = vbase[graph->head[e]];

               assert(graphmark[k]);
               node = (l == UNKNOWN || !graphmark[l] )? UNKNOWN : SCIPStpunionfindFind(&uf, l);

               /* does the boundary-path end in the root component? */
               if( node != UNKNOWN && node != crucnode && graphmark[l] )
               {
                  SCIP_CALL( SCIPpairheapInsert(scip, &boundpaths[crucnode], e, edgecost, &(heapsize[crucnode])) );
                  break;
               }
            }

            if( boundpaths[crucnode] == NULL )
            {
               oldedge = UNKNOWN;
            }
            else
            {
               oldedge = e;
            }

            /* counts the nodes connected during the following 'preprocessing' */
            count = 0;

            /* try to connect the nodes of C (directly) to COMP(C), as a preprocessing for voronoi-repair */
            for( k = 0; k < nkpnodes; k++ )
            {
               blists_curr = blists_start[kpnodes[k]];
               assert( blists_curr != NULL );
               while( blists_curr != NULL )
               {
                  node = blists_curr->index;

                  /* iterate through all outgoing edges of 'node' */
                  for( edge = graph->inpbeg[node]; edge != EAT_LAST; edge = graph->ieat[edge] )
                  {
                     adjnode = graph->tail[edge];

                     /* check whether the adjacent node is not in C and allows a better voronoi assignment of the current node */
                     if( state[adjnode] == CONNECT && SCIPisGT(scip, vnoi[node].dist, vnoi[adjnode].dist + graph->cost[edge])
                        && graphmark[vbase[adjnode]] && graphmark[adjnode] )
                     {
                        vnoi[node].dist = vnoi[adjnode].dist + graph->cost[edge];
                        vbase[node] = vbase[adjnode];
                        vnoi[node].edge = edge;
                     }
                  }
                  if( vbase[node] != UNKNOWN )
                  {
                     heap_add(graph->path_heap, state, &count, node, vnoi);
                  }
                  blists_curr = blists_curr->parent;
               }
            }
            if( nkpnodes > 0 )
               assert(count > 0);
            newedge = UNKNOWN;

            /* if there is no key path, nothing has to be repaired */
            if( nkpnodes > 0 )
               graph_voronoiRepair(scip, graph, graph->cost, &count, vbase, vnoi, &newedge, crucnode, &uf);
            else
               newedge = linkcutNodes[crucnode].edge;

            if( oldedge != UNKNOWN && newedge != UNKNOWN && SCIPisLT(scip, edgecost,
                  vnoi[graph->tail[newedge]].dist + graph->cost[newedge] + vnoi[graph->head[newedge]].dist) )
               newedge = oldedge;

            if( oldedge != UNKNOWN && newedge == UNKNOWN )
               newedge = oldedge;

            assert( newedge != UNKNOWN );

            edgecost = vnoi[graph->tail[newedge]].dist + graph->cost[newedge] + vnoi[graph->head[newedge]].dist;

            if( mwpc )
            {
               prizemarkcount = 0;
               edgecost -= getNewPrize(graph, solNodes, graphmark, newedge, prizemark, prizemarklist, &prizemarkcount);

               for( node = graph->tail[newedge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                  edgecost -= getNewPrize(graph, solNodes, graphmark, vnoi[node].edge, prizemark, prizemarklist, &prizemarkcount);

               for( node = graph->head[newedge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                  edgecost -= getNewPrize(graph, solNodes, graphmark, vnoi[node].edge, prizemark, prizemarklist, &prizemarkcount);

               for( int pi = 0; pi < prizemarkcount; pi++ )
                  prizemark[prizemarklist[pi]] = FALSE;
            }

            if( SCIPisLT(scip, edgecost, kpathcost) )
            {
               solimproved = TRUE;
               node = SCIPStpunionfindFind(&uf, vbase[graph->head[newedge]]);

               SCIPdebugMessage( "ADDING NEW KEY PATH (%f )\n", edgecost - kpathcost );

#ifndef NDEBUG
               assert((kpathcost - edgecost) >= 0.0);
               objimprovement += (kpathcost - edgecost);
#endif

               localmoves++;

               /* remove old keypath */
               assert(  solEdges[flipedge(linkcutNodes[crucnode].edge)] != UNKNOWN );
               solEdges[flipedge(linkcutNodes[crucnode].edge)] = UNKNOWN;
               solNodes[crucnode] = FALSE;
               graphmark[crucnode] = FALSE;

               for( k = 0; k < nkpnodes; k++ )
               {
                  assert(  solEdges[flipedge(linkcutNodes[kpnodes[k]].edge)] != UNKNOWN );
                  solEdges[flipedge(linkcutNodes[kpnodes[k]].edge)] = UNKNOWN;
                  solNodes[kpnodes[k]] = FALSE;
                  graphmark[kpnodes[k]] = FALSE;
               }
               assert(graphmark[kptailnode]);

               if( node == crucnode )
                  newedge = flipedge(newedge);

               for( node = graph->tail[newedge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
               {
                  graphmark[node] = FALSE;

                  solEdges[flipedge(vnoi[node].edge)] = CONNECT;
                  solEdges[vnoi[node].edge] = UNKNOWN;
               }

               for( node = graph->head[newedge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
               {
                  graphmark[node] = FALSE;

                  solEdges[vnoi[node].edge] = CONNECT;
               }

               solEdges[flipedge(newedge)] = CONNECT;

               newpathend = vbase[graph->tail[newedge]];
               assert(node == vbase[graph->head[newedge]] );

               /* flip all edges on the ST path between the endnode of the new key-path and the current crucial node */
               k = newpathend;
               assert(SCIPStpunionfindFind(&uf, newpathend) == crucnode);

               while( k != crucnode )
               {
                  assert(graphmark[k]);
                  assert( solEdges[flipedge(linkcutNodes[k].edge)] != -1);
                  solEdges[flipedge(linkcutNodes[k].edge)] = UNKNOWN;

                  solEdges[linkcutNodes[k].edge] = CONNECT;

                  k = graph->head[linkcutNodes[k].edge];
               }

               for( k = 0; k < i; k++ )
               {
                  if( crucnode == SCIPStpunionfindFind(&uf, dfstree[k]) )
                  {
                     graphmark[dfstree[k]] = FALSE;
                     solNodes[dfstree[k]] = FALSE;
                  }
               }

               /* update union find */
               if( !Is_term(graph->term[node]) && scanned[node] && !pinned[node] && SCIPStpunionfindFind(&uf, node) == node )
               {
                  for( edge = graph->outbeg[node]; edge != EAT_LAST; edge = graph->oeat[edge] )
                  {
                     adjnode = graph->head[edge];
                     /* check whether edge 'edge' leads to an ancestor of terminal 'node' */
                     if( solEdges[edge] == CONNECT && solNodes[adjnode]  && graphmark[adjnode] && SCIPStpunionfindFind(&uf, adjnode) != node )
                     {
                        assert(scanned[adjnode]);
                        /* meld the heaps */
                        SCIPpairheapMeldheaps(scip, &boundpaths[node], &boundpaths[adjnode], &heapsize[node], &heapsize[adjnode]);

                        /* update the union-find data structure */
                        SCIPStpunionfindUnion(&uf, node, adjnode, FALSE);

                        /* move along the key-path until its end (i.e. until a crucial node is reached) */
                        while( !nodeIsCrucial(graph, solEdges, adjnode) && !pinned[adjnode] )
                        {
                           for( e = graph->outbeg[adjnode]; e != EAT_LAST; e = graph->oeat[e] )
                           {
                              if( solEdges[e] != -1 )
                                 break;
                           }

                           /* assert that each leaf of the ST is a terminal */
                           assert( e != EAT_LAST );
                           adjnode = graph->head[e];
                           if( !solNodes[adjnode]  )
                              break;
                           assert(scanned[adjnode]);
                           assert(SCIPStpunionfindFind(&uf, adjnode) != node);

                           /* update the union-find data structure */
                           SCIPStpunionfindUnion(&uf, node, adjnode, FALSE);

                           /* meld the heaps */
                           SCIPpairheapMeldheaps(scip, &boundpaths[node], &boundpaths[adjnode], &heapsize[node], &heapsize[adjnode]);
                        }
                     }
                  }

               }
               pinned[node] = TRUE;
            }

            /* restore the original voronoi digram */
            l = 0;
            for( k = 0; k < nkpnodes; k++ )
            {
               /* reset all nodes having the current (internal) keypath node as their voronoi base */
               blists_curr = blists_start[kpnodes[k]];
               while( blists_curr != NULL )
               {
                  node = blists_curr->index;
                  vbase[node] = memvbase[l];
                  vnoi[node].dist = memdist[l];
                  vnoi[node].edge = meminedges[l];
                  l++;
                  blists_curr = blists_curr->parent;
               }
            }
            assert(l == nresnodes);
         }
      }


      /**********************************************************/

   TERMINATE:

      SCIPStpunionfindClear(scip, &uf, nnodes);

      /* free data structures */
      SCIPfreeBufferArray(scip, &kpedges);
      SCIPfreeBufferArray(scip, &kpnodes);
      SCIPfreeBufferArray(scip, &supernodes);

      for( k = nnodes - 1; k >= 0; k-- )
      {
         if( boundpaths[k] != NULL )
            SCIPpairheapFree(scip, &boundpaths[k]);

         blists_curr = blists_start[k];
         lvledges_curr = lvledges_start[k];
         while( lvledges_curr != NULL )
         {
            lvledges_start[k] = lvledges_curr->parent;
            SCIPfreeBlockMemory(scip, &lvledges_curr);
            lvledges_curr = lvledges_start[k];
         }

         while( blists_curr != NULL )
         {
            blists_start[k] = blists_curr->parent;
            SCIPfreeBlockMemory(scip, &blists_curr);
            blists_curr = blists_start[k];
         }
      }

      /* has there been a move during this run? */
      if( localmoves > 0 )
      {
         for( i = 0; i < nnodes; i++ )
         {
            solNodes[i] = FALSE;
            graphmark[i] = (graph->grad[i] > 0);
            SCIPlinkcuttreeInit(&linkcutNodes[i]);
         }

         graphmark[root] = TRUE;

         /* create a link-cut tree representing the current Steiner tree */
         for( e = 0; e < nedges; e++ )
         {
            assert(graph->head[e] == graph->tail[flipedge(e)]);

            /* if edge e is in the tree, so are its incident vertices */
            if( solEdges[e] != -1 )
            {
               solNodes[graph->tail[e]] = TRUE;
               solNodes[graph->head[e]] = TRUE;
               SCIPlinkcuttreeLink(&linkcutNodes[graph->head[e]], &linkcutNodes[graph->tail[e]], flipedge(e));
            }
         }
         assert( linkcutNodes[root].edge == -1 );
         linkcutNodes[root].edge = -1;
      }
   }

   /* free data structures */
   SCIPStpunionfindFreeMembers(scip, &uf);
   SCIPfreeBufferArray(scip, &nodesmark);
   SCIPfreeBufferArray(scip, &dfstree);
   SCIPfreeBufferArray(scip, &pinned);
   SCIPfreeBufferArray(scip, &boundpaths);
   SCIPfreeBufferArray(scip, &meminedges);
   SCIPfreeBufferArray(scip, &memdist);
   SCIPfreeBufferArray(scip, &memvbase);
   SCIPfreeBufferArray(scip, &blists_start);
   SCIPfreeBufferArray(scip, &heapsize);
   SCIPfreeBufferArray(scip, &scanned);
   SCIPfreeBufferArrayNull(scip, &prizemarklist);
   SCIPfreeBufferArrayNull(scip, &prizemark);
   SCIPfreeBufferArray(scip, &supernodesid);
   SCIPfreeBufferArray(scip, &boundedges);
   SCIPfreeBufferArray(scip, &lvledges_start);
   SCIPfreeBufferArray(scip, &newedges);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &vnoi);
   /******/

   if( solimproved )
   {
      SCIP_CALL( SCIPStpHeurTMpruneEdgeSol(scip, graph, solEdges) );
   }

#ifndef NDEBUG
   {
      const SCIP_Real newobj = graph_sol_getObj(graph->cost, solEdges, 0.0, nedges);
      SCIPdebugMessage("key vertex heuristic obj before/after: %f/%f \n", initialobj, newobj);
      assert(SCIPisLE(scip, newobj + objimprovement, initialobj));
   }
#endif

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyLocal)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPStpIncludeHeurLocal(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeLocal)
{   /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolLocal)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);

   heurdata->nfails = 1;
   heurdata->nbestsols = DEFAULT_NBESTSOLS;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(heurdata->lastsolindices), heurdata->maxnsols) );

   for( int i = 0; i < heurdata->maxnsols; i++ )
      heurdata->lastsolindices[i] = -1;

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolLocal)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->lastsolindices != NULL);
   SCIPfreeMemoryArray(scip, &(heurdata->lastsolindices));

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLocal)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_PROBDATA* probdata;
   GRAPH* graph;                             /* graph structure */
   SCIP_SOL* newsol;                         /* new solution */
   SCIP_SOL* impsol;                         /* new improved solution */
   SCIP_SOL** sols;                          /* solutions */
   SCIP_VAR** vars;                          /* SCIP variables */
   SCIP_Real pobj;
   SCIP_Real* nval;
   SCIP_Real* xval;
   int v;
   int min;
   int root;
   int nvars;
   int nsols;                                /* number of all solutions found so far */
   int nedges;
   int* results;
   int* lastsolindices;
   SCIP_Bool feasible;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   lastsolindices = heurdata->lastsolindices;
   assert(lastsolindices != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   *result = SCIP_DIDNOTRUN;

   /* the local heuristics may not work correctly for several problem variants*/
   if( graph->stp_type != STP_SPG && graph->stp_type != STP_RSMT && graph->stp_type != STP_OARSMT &&
      graph->stp_type != STP_PCSPG && graph->stp_type != STP_RPCSPG && graph->stp_type != STP_GSTP
      && graph->stp_type != STP_MWCSP )
      return SCIP_OKAY;

   /* don't run local in a Subscip */
   if( SCIPgetSubscipDepth(scip) > 0 )
      return SCIP_OKAY;

   /* no solution available? */
   if( SCIPgetBestSol(scip) == NULL )
      return SCIP_OKAY;

   root = graph->source;
   sols = SCIPgetSols(scip);
   nsols = SCIPgetNSols(scip);
   nedges = graph->edges;

   assert(heurdata->maxnsols >= 0);

   min = MIN(heurdata->maxnsols, nsols);

   /* only process each solution once */
   for( v = 0; v < min; v++ )
   {
      if( SCIPsolGetIndex(sols[v]) != lastsolindices[v] )
      {
         /* shift all solution indices right of the new solution index */
         for( int i = min - 1; i >= v + 1; i-- )
            lastsolindices[i] = lastsolindices[i - 1];
         break;
      }
   }

   /* no new solution available? */
   if( v == min )
      return SCIP_OKAY;

   newsol = sols[v];
   lastsolindices[v] = SCIPsolGetIndex(newsol);

   /* solution not good enough? */
   if( (v > heurdata->nbestsols && !(heurdata->maxfreq)) && graph->stp_type != STP_MWCSP )
      return SCIP_OKAY;

   /* has the new solution been found by this very heuristic? */
   if( SCIPsolGetHeur(newsol) == heur )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   vars = SCIPprobdataGetVars(scip);
   nvars = SCIPprobdataGetNVars(scip);
   xval = SCIPprobdataGetXval(scip, newsol);

   if( vars == NULL )
      return SCIP_OKAY;

   assert(vars != NULL);
   assert(xval != NULL);

   /* for PC/MW: test whether solution is trivial */
   if( graph->stp_type == STP_PCSPG || graph->stp_type == STP_MWCSP )
   {
      int e;
      assert(graph->extended);

      for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
         if( !Is_term(graph->term[graph->head[e]]) && SCIPisEQ(scip, xval[e], 1.0) )
            break;
      if( e == EAT_LAST )
         return SCIP_OKAY;
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &results, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nval, nvars) );

   /* set solution array */
   for( int e = 0; e < nedges; e++ )
   {
      if( SCIPisEQ(scip, xval[e], 1.0) )
         results[e] = CONNECT;
      else
         results[e] = UNKNOWN;
   }

   if( !graph_sol_valid(scip, graph, results) )
   {
      SCIPfreeBufferArray(scip, &nval);
      SCIPfreeBufferArray(scip, &results);
      return SCIP_OKAY;
   }

   /* pruning necessary? */
   if( SCIPsolGetHeur(newsol) == NULL ||
      !(strcmp(SCIPheurGetName(SCIPsolGetHeur(newsol)), "rec") == 0 ||
         strcmp(SCIPheurGetName(SCIPsolGetHeur(newsol)), "TM") == 0) )
   {
      const int nnodes = graph->knots;
      STP_Bool* steinertree;
      SCIP_CALL( SCIPallocBufferArray(scip, &steinertree, nnodes) );
      assert(graph_sol_valid(scip, graph, results));

      graph_sol_setVertexFromEdge(graph, results, steinertree);

      for( int e = 0; e < nedges; e++ )
         results[e] = UNKNOWN;

      if( graph_pc_isPcMw(graph) )
         SCIP_CALL( SCIPStpHeurTMPrunePc(scip, graph, graph->cost, results, steinertree) );
      else
         SCIP_CALL( SCIPStpHeurTMPrune(scip, graph, graph->cost, 0, results, steinertree) );

      SCIPfreeBufferArray(scip, &steinertree);
   }

   /* execute local heuristics */
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, results) );

#if 0
   if( graph_pc_isPcMw(graph) )
      SCIP_CALL( SCIPStpHeurLocalExtendPcMwImp(scip, graph, results) );
#endif

   /* can we connect the network */
   for( v = 0; v < nvars; v++ )
      nval[v] = (results[v % nedges] == (v / nedges)) ? 1.0 : 0.0;

   SCIP_CALL( SCIPStpValidateSol(scip, graph, nval, &feasible) );

   /* solution feasible? */
   if( feasible )
   {
      assert(nedges == nvars);

      pobj = 0.0;

      for( v = 0; v < nedges; v++ )
         pobj += graph->cost[v] * nval[v];

      /* has solution been improved? */
      if( SCIPisGT(scip, SCIPgetSolOrigObj(scip, newsol) - SCIPprobdataGetOffset(scip), pobj) )
      {
         SCIP_SOL* bestsol;
         SCIP_Bool success;

         bestsol = sols[0];
         impsol = NULL;
         SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, impsol, heur, &success) );

         if( success )
         {
            *result = SCIP_FOUNDSOL;

            if( heurdata->nbestsols < heurdata->maxnsols && SCIPisGT(scip, SCIPgetSolOrigObj(scip, bestsol) - SCIPprobdataGetOffset(scip), pobj) )
            {
               heurdata->nfails = 0;
               heurdata->nbestsols++;
            }
            SCIPdebugMessage("success in local: old: %f new: %f \n", (SCIPgetSolOrigObj(scip, bestsol) - SCIPprobdataGetOffset(scip)), pobj);
         }
      }
   }

   if( *result != SCIP_FOUNDSOL )
   {
      heurdata->nfails++;
      if( heurdata->nbestsols > DEFAULT_MINNBESTSOLS && heurdata->nfails > 1 && graph->stp_type != STP_MWCSP )
         heurdata->nbestsols--;

      SCIPdebugMessage("fail! %d \n", heurdata->nbestsols);
   }

   SCIPfreeBufferArray(scip, &nval);
   SCIPfreeBufferArray(scip, &results);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */



/** perform local heuristics on a given Steiner tree todo delete cost parameter */
SCIP_RETCODE SCIPStpHeurLocalRun(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   int*                  solEdges            /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   )
{
   NODE* linkcutNodes;
   const int root = graph->source;
   const int nnodes = graph->knots;
   const int probtype = graph->stp_type;
   STP_Bool* solNodes;
   const STP_Bool mw = (probtype == STP_MWCSP);
   const STP_Bool mwpc = graph_pc_isPcMw(graph);
#ifndef NDEBUG
   const SCIP_Real initialobj = graph_sol_getObj(graph->cost, solEdges, 0.0, graph->edges);
#endif

   assert(graph && solEdges);
   assert(graph_valid(scip, graph));

   if( graph->grad[root] == 0 || graph->terms == 1 )
      return SCIP_OKAY;

   if( mwpc )
   {
      assert(graph->extended);

      if( solIsTrivialPcMw(graph, solEdges) )
         return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &linkcutNodes, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solNodes, nnodes) );

   if( mwpc )
      SCIP_CALL( SCIPStpHeurLocalExtendPcMw(scip, graph, graph->cost, solEdges, solNodes) );

   markSolTreeNodes(graph, solEdges, linkcutNodes, solNodes);

   assert(linkcutNodes[root].edge == -1);

   /* Call first local heuristic? */
   if( probtype == STP_SPG || probtype == STP_RSMT || probtype == STP_OARSMT || probtype == STP_GSTP || (mwpc) )
   {
      SCIP_CALL( localVertexInsertion(scip, graph, solNodes, linkcutNodes, solEdges) );
   }

   assert(graph_sol_valid(scip, graph, solEdges));

   /* run Key-Vertex Elimination & Key-Path Exchange heuristics? */
   if( !mw )
   {
      SCIP_CALL( localKeyVertexHeuristics(scip, graph, solNodes, linkcutNodes, solEdges) );
   }

#ifndef NDEBUG
   {
      const SCIP_Real newobj = graph_sol_getObj(graph->cost, solEdges, 0.0, graph->edges);
      assert(SCIPisLE(scip, newobj, initialobj));
      assert(graph_sol_valid(scip, graph, solEdges));
   }
#endif

   SCIPfreeBufferArray(scip, &solNodes);
   SCIPfreeBufferArray(scip, &linkcutNodes);

   return SCIP_OKAY;
}

/** Implication based local heuristic for (R)PC and MW */
SCIP_RETCODE SCIPStpHeurLocalExtendPcMwImp(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int*                  result              /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
)
{
   const int* starts = SCIPStpGetPcImplStarts(scip);
   const int* verts = SCIPStpGetPcImplVerts(scip);

   assert(graph_pc_isPcMw(graph));

   if( starts != NULL )
   {
      const int nnodes = graph->knots;
      STP_Bool* stvertex;
      int nfound = 0;
      int ptermcount = 0;

      assert(graph->extended);
      assert(verts != NULL);

      SCIPallocBufferArray(scip, &stvertex, nnodes);

      graph_sol_setVertexFromEdge(graph, result, stvertex);

      for( int i = 0; i < nnodes; i++ )
      {
         if( !Is_pterm(graph->term[i]) )
            continue;

         assert(!graph_pc_knotIsFixedTerm(graph, i));

         ptermcount++;

         if( stvertex[i] )
            continue;

         for( int j = starts[ptermcount - 1]; j < starts[ptermcount]; j++ )
         {
            const int vert = verts[j];
            if( stvertex[vert] )
            {
               /* now connect the vertex */

               graph_knot_printInfo(graph, i);
               nfound++;
               break;
            }
         }
      }

      assert(ptermcount == graph_pc_nPotentialTerms(graph));

      if( nfound > 0 )
      {
         printf("nfound: %d \n\n\n", nfound);
         /* todo prune! */
         //return SCIP_ERROR;
      }
      else
         printf("none %d \n", 0);

      SCIPfreeBufferArray(scip, &stvertex);
   }
   return SCIP_OKAY;
}

/** Greedy Extension local heuristic for (R)PC and MW */
SCIP_RETCODE SCIPStpHeurLocalExtendPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge cost array */
   int*                  stedge,             /**< initialized array to indicate whether an edge is part of the Steiner tree */
   STP_Bool*             stvertex            /**< uninitialized array to indicate whether a vertex is part of the Steiner tree */
   )
{
   GNODE candidates[MAX(GREEDY_EXTENSIONS, GREEDY_EXTENSIONS_MW)];
   int candidatesup[MAX(GREEDY_EXTENSIONS, GREEDY_EXTENSIONS_MW)];

   PATH* path;
   PATH* orgpath;
   SCIP_PQUEUE* pqueue;
   SCIP_Real bestsolval;

   int nextensions;
   const int greedyextensions = (graph->stp_type == STP_MWCSP) ? GREEDY_EXTENSIONS_MW : GREEDY_EXTENSIONS;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;
   const int root = graph->source;
   STP_Bool* stvertextmp;
   SCIP_Bool extensions = FALSE;

#ifndef NDEBUG
   const SCIP_Real initialobj = graph_sol_getObj(graph->cost, stedge, 0.0, nedges);
#endif

#ifdef NEW
   SCIP_Real* costbiased;
   SCIP_Real* prizebiased;
   SCIP_CALL( SCIPallocBufferArray(scip, &costbiased, graph->edges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &prizebiased, graph->knots) );
#endif

   assert(scip != NULL);
   assert(graph != NULL);
   assert(stedge != NULL);
   assert(cost != NULL);
   assert(stvertex != NULL);
   assert(graph->extended);

   graph_pc_2transcheck(graph);
   SCIP_CALL( SCIPallocBufferArray(scip, &stvertextmp, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orgpath, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );

   /* initialize solution vertex array with FALSE */
   BMSclearMemoryArray(stvertex, nnodes);

   stvertex[root] = TRUE;

   for( int j = 0; j < nnodes; j++ )
      path[j].edge = UNKNOWN;

   for( int e = 0; e < nedges; e++ )
      if( stedge[e] == CONNECT )
      {
         path[graph->head[e]].edge = e;
         stvertex[graph->head[e]] = TRUE;
      }

   for( int e = 0; e < nedges; e++ )
      if( stedge[e] == CONNECT )
         assert(stvertex[graph->tail[e]]);

#ifdef NEW
   graph_pc_getBiased(scip, graph, TRUE, costbiased, prizebiased);
   graph_path_st_pcmw_extendBiased(scip, graph, costbiased, prizebiased, path, stvertex, &extensions);
#else
   graph_path_st_pcmw_extend(scip, graph, cost, FALSE, path, stvertex, &extensions);
#endif

   BMScopyMemoryArray(orgpath, path, nnodes);

   /*** compute solution value and save greedyextensions many best unconnected nodes  ***/

   SCIP_CALL( SCIPpqueueCreate(&pqueue, greedyextensions, 2.0, GNODECmpByDist) );

   assert(orgpath[root].edge == UNKNOWN);

   bestsolval = 0.0;
   nextensions = 0;
   for( int i = 0; i < nnodes; i++ )
   {
      if( graph->grad[i] == 0 || root == i )
         continue;

      if( Is_term(graph->term[i]) && !graph_pc_knotIsFixedTerm(graph, i) )
         continue;

      if( stvertex[i] )
      {
         assert(orgpath[i].edge >= 0);

         bestsolval += graph->cost[orgpath[i].edge];

         if( Is_pterm(graph->term[i]) )
            bestsolval -= graph->prize[i];
      }
      else if( orgpath[i].edge != UNKNOWN && Is_pterm(graph->term[i]) )
      {
         SCIP_CALL( addToCandidates(scip, graph, path, i, greedyextensions, &nextensions, candidates, pqueue) );
      }
   }

   for( int restartcount = 0; restartcount < GREEDY_MAXRESTARTS && !graph_pc_isRootedPcMw(graph); restartcount++ )
   {
      int l = 0;
      SCIP_Bool extensionstmp = FALSE;
      int extcount = nextensions;

      /* write extension candidates into array, from max to min */
      while( SCIPpqueueNElems(pqueue) > 0 )
      {
         GNODE* min = (GNODE*) SCIPpqueueRemove(pqueue);
         assert(extcount > 0);
         candidatesup[--extcount] = min->number;
      }
      assert(extcount == 0);

      /* iteratively insert new subpaths and try to improve solution */
      for( ; l < nextensions; l++ )
      {
         const int extensioncand = candidatesup[l];
         if( !stvertex[extensioncand] )
         {
            SCIP_Real newsolval = 0.0;
            int k = extensioncand;

            BMScopyMemoryArray(stvertextmp, stvertex, nnodes);
            BMScopyMemoryArray(path, orgpath, nnodes);

            /* add new extension */
            while( !stvertextmp[k] )
            {
               stvertextmp[k] = TRUE;
               assert(orgpath[k].edge != UNKNOWN);
               k = graph->tail[orgpath[k].edge];
               assert(k != extensioncand);
            }
#ifdef NEW
            assert(graph_sol_valid(scip, graph, stedge));
            graph_path_st_pcmw_extendBiased(scip, graph, costbiased, prizebiased, path, stvertextmp, &extensionstmp);

#else
            graph_path_st_pcmw_extend(scip, graph, cost, TRUE, path, stvertextmp, &extensionstmp);
#endif

            for( int j = 0; j < nnodes; j++ )
            {
               if( graph->grad[j] == 0 || root == j )
                  continue;

               if( Is_term(graph->term[j]) && !graph_pc_knotIsFixedTerm(graph, j) )
                  continue;

               if( stvertextmp[j] )
               {
                  assert(path[j].edge >= 0);

                  newsolval += graph->cost[path[j].edge];

                  if( Is_pterm(graph->term[j]) )
                     newsolval -= graph->prize[j];
               }
            }

            /* new solution value better than old one? */
            if( SCIPisLT(scip, newsolval, bestsolval) )
            {
               extensions = TRUE;
               bestsolval = newsolval;
               BMScopyMemoryArray(stvertex, stvertextmp, nnodes);
               BMScopyMemoryArray(orgpath, path, nnodes);

               /* save greedyextensions many best unconnected nodes  */
               nextensions = 0;

               for( int j = 0; j < nnodes; j++ )
                  if( !stvertex[j] && Is_pterm(graph->term[j]) && path[j].edge != UNKNOWN )
                     SCIP_CALL( addToCandidates(scip, graph, path, j, greedyextensions, &nextensions, candidates, pqueue) );

               break;
            } /* if new solution value better than old one? */
         } /* if !stvertex[i] */
      } /* for l < nextension */

      /* no more extensions performed? */
      if( l == nextensions )
         break;
   } /* main loop */

   /* have vertices been added? */
   if( extensions )
   {
      for( int e = 0; e < nedges; e++ )
         stedge[e] = UNKNOWN;
      SCIP_CALL( SCIPStpHeurTMPrunePc(scip, graph, graph->cost, stedge, stvertex) );
   }

   SCIPpqueueFree(&pqueue);
   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &orgpath);
   SCIPfreeBufferArray(scip, &stvertextmp);

#ifdef NEW
   SCIPfreeBufferArray(scip, &prizebiased);
   SCIPfreeBufferArray(scip, &costbiased);
#endif

#ifndef NDEBUG
   assert(SCIPisLE(scip, graph_sol_getObj(graph->cost, stedge, 0.0, nedges), initialobj));
#endif

   return SCIP_OKAY;
}

/** Greedy Extension local heuristic for (R)PC and MW */
SCIP_RETCODE SCIPStpHeurLocalExtendPcMwOut(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   int*                  stedge,             /**< initialized array to indicate whether an edge is part of the Steiner tree */
   STP_Bool*             stvertex            /**< uninitialized array to indicate whether a vertex is part of the Steiner tree */
   )
{
   int candidates[GREEDY_EXTENSIONS];
   int ncandidates;
   DHEAP* dheap;
   STP_Bool* stvertextmp;
   SCIP_Real* dist;
   int* pred;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;
   SCIP_Bool extensions = FALSE;
   int maxnode;
   const SCIP_Bool isexended = graph->extended;

#ifndef NDEBUG
   const SCIP_Real initialobj = graph_sol_getObj(graph->cost, stedge, 0.0, nedges);
#endif

   assert(scip && graph && stedge && stvertex);

   graph_pc_2orgcheck(graph);

   graph_sol_setVertexFromEdge(graph, stedge, stvertex);

   /* compute candidates for extension */

   maxnode = -1;
   ncandidates = 0;

   for( int k = 0; k < nnodes; k++ )
      if( graph->mark[k] && !stvertex[k] && Is_term(graph->term[k]) && !graph_pc_termIsNonLeaf(graph, k) )
      {
         assert(graph->mark[k]);

         if( maxnode == -1 || graph->prize[k] > graph->prize[maxnode] )
            maxnode = k;
      }

   if( maxnode != -1 )
   {
      SCIP_RANDNUMGEN* randnumgen;
      int shift;

      SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );

      SCIP_CALL( SCIPallocBufferArray(scip, &dist, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pred, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &stvertextmp, nnodes) );

      graph_heap_create(scip, nnodes, NULL, NULL, &dheap);
      graph_init_csr(scip, graph);

      shift = SCIPrandomGetInt(randnumgen, 0, nnodes - 1);
      ncandidates = 1;
      candidates[0] = maxnode;

      for( int k = 0; k < nnodes && ncandidates < GREEDY_EXTENSIONS; k++ )
      {
         const int node = (k + shift) % nnodes;
         if( graph->mark[k] && !stvertex[node] && Is_term(graph->term[node])
            && !graph_pc_termIsNonLeaf(graph, node) && node != maxnode )
         {
            assert(graph->mark[node]);
            candidates[ncandidates++] = node;
         }
      }

      SCIPfreeRandom(scip, &randnumgen);
   }

   /* main loop */
   for( int k = 0; k < ncandidates; k++ )
   {
      const int cand = candidates[k];
      SCIP_Bool success = FALSE;

      if( stvertex[cand] )
      {
         assert(k > 0);
         continue;
      }

      graph_path_st_pcmw_extendOut(scip, graph, cand, stvertex, dist, pred, stvertextmp, dheap, &success);

      if( success )
         extensions = TRUE;
   }

   /* have vertices been added? */
   if( extensions )
   {
      graph_pc_2trans(graph);

      for( int e = 0; e < nedges; e++ )
         stedge[e] = UNKNOWN;
      SCIP_CALL( SCIPStpHeurTMPrunePc(scip, graph, graph->cost, stedge, stvertex) );
   }

   if( maxnode != -1 )
   {
      graph_heap_free(scip, TRUE, TRUE, &dheap);
      graph_free_csr(scip, graph);

      SCIPfreeBufferArray(scip, &stvertextmp);
      SCIPfreeBufferArray(scip, &pred);
      SCIPfreeBufferArray(scip, &dist);
   }

#ifndef NDEBUG
   assert(SCIPisLE(scip, graph_sol_getObj(graph->cost, stedge, 0.0, nedges), initialobj));
#endif

   if( isexended && !graph->extended )
      graph_pc_2trans(graph);

   if( !isexended && graph->extended )
      graph_pc_2org(graph);

   return SCIP_OKAY;
}


/** creates the local primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPStpIncludeHeurLocal(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Local primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecLocal, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyLocal) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeLocal) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolLocal) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolLocal) );

   /* add local primal heuristic parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "stp/duringroot",
         "should the heuristic be called during the root node?",
         &heurdata->duringroot, FALSE, DEFAULT_DURINGROOT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/maxfreq",
         "should the heuristic be executed at maximum frequeny?",
         &heurdata->maxfreq, FALSE, DEFAULT_MAXFREQLOC, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxnsols",
         "maximum number of best solutions to improve",
         &heurdata->maxnsols, FALSE, DEFAULT_MAXNBESTSOLS, 1, 50, NULL, NULL) );

   return SCIP_OKAY;
}
