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
#define DEFAULT_MAXNBESTSOLS  25
#define DEFAULT_NBESTSOLS     10
#define DEFAULT_MINNBESTSOLS  6
#define LOCAL_MAXRESTARTS  5

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

/** perform local heuristics on a given Steiner tree todo delete cost parameter */
SCIP_RETCODE SCIPStpHeurLocalRun(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< arc cost array todo delete this parameter and use graph->cost */
   int*                  best_result         /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   )
{
   NODE* nodes;
   PATH* vnoi;
   int* vbase;
   int* const graphmark = graph->mark;
   int e;
   int i;
   int k;
   const int root = graph->source;
   const int nnodes = graph->knots;
   const int nedges = graph->edges;
   const int probtype = graph->stp_type;
   int newnverts;
   STP_Bool* steinertree;
   const STP_Bool pc = ((probtype == STP_PCSPG) || (probtype == STP_RPCSPG));
   const STP_Bool mw = (probtype == STP_MWCSP);
   const STP_Bool mwpc =  (pc || mw);
#ifndef NDEBUG
   const SCIP_Real initialobj = graph_sol_getObj(graph->cost, best_result, 0.0, nedges);
#endif

   assert(graph != NULL);
   assert(cost != NULL);
   assert(best_result != NULL);
   assert(graph_valid(graph));

   newnverts = 0;

   if( graph->grad[root] == 0 || graph->terms == 1 )
      return SCIP_OKAY;

   /* for PC variants test whether solution is trivial */
   if( mwpc )
   {
      for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
         if( !Is_term(graph->term[graph->head[e]]) && best_result[e] )
            break;

      if( e == EAT_LAST )
      {
         SCIPdebugMessage("Local heuristic: return trivial \n");
         return SCIP_OKAY;
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &steinertree, nnodes) );

   if( mwpc )
   {
      SCIP_Bool dummy;
      SCIP_CALL( SCIPStpHeurLocalExtendPcMw(scip, graph, cost, vnoi, best_result, steinertree, &dummy) );
   }

   for( i = 0; i < nnodes; i++ )
   {
      steinertree[i] = FALSE;
      SCIPlinkcuttreeInit(&nodes[i]);
   }

   /* create a link-cut tree representing the current Steiner tree */
   for( e = 0; e < nedges; e++ )
   {
      assert(graph->head[e] == graph->tail[flipedge(e)]);

      /* if edge e is in the tree, so are its incident vertices */
      if( best_result[e] == CONNECT )
      {
         steinertree[graph->tail[e]] = TRUE;
         steinertree[graph->head[e]] = TRUE;
         SCIPlinkcuttreeLink(&nodes[graph->head[e]], &nodes[graph->tail[e]], flipedge(e));
      }
   }

   assert( nodes[root].edge == -1 );
   nodes[root].edge = -1;

   /* VERTEX INSERTION */
   if( probtype == STP_SPG || probtype == STP_RSMT || probtype == STP_OARSMT || probtype == STP_GSTP || (mwpc) )
   {
      int newnode;
      int* insert;
      int* adds;
      int* cuts;
      int* cuts2;
      int* stdeg;

      SCIP_CALL( SCIPallocBufferArray(scip, &insert, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &adds, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &cuts, nnodes) );

      if( mw )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &cuts2, nnodes) );
         SCIP_CALL( SCIPallocBufferArray(scip, &stdeg, nnodes) );

         BMSclearMemoryArray(stdeg, nnodes);

         for( e = 0; e < nedges; e++ )
            if( best_result[e] == CONNECT )
            {
               stdeg[graph->tail[e]]++;
               stdeg[graph->head[e]]++;
            }
      }
      else
      {
         cuts2 = NULL;
         stdeg = NULL;
      }

      i = 0;
      newnode = 0;

      for( ;; )
      {
         SCIP_Real diff;

         /* if vertex i is not in the current ST and has at least two adjacent nodes, it might be added */
         if( !steinertree[i] && graph->grad[i] > 1 && (!mwpc || !Is_term(graph->term[i])) )
         {
            NODE* v;
            int counter;
            int lastnodeidx;
            int insertcount = 0;

            /* if an outgoing edge of vertex i points to the current ST, SCIPlinkcuttreeLink the edge to a list */
            for( int oedge = graph->outbeg[i]; oedge != EAT_LAST; oedge = graph->oeat[oedge])
               if( steinertree[graph->head[oedge]] && (!mwpc || !Is_term(graph->term[graph->head[oedge]])) )
                  insert[insertcount++] = oedge;

            /* if there are less than two edges connecting node i and the current tree, continue */
            if( insertcount <= 1 )
               goto ENDOFLOOP;

            if( mw )
               SCIPlinkcuttreeInit(&nodes[i]);

            /* the node to insert */
            v = &nodes[i];

            SCIPlinkcuttreeLink(v, &nodes[graph->head[insert[0]]], insert[0]);

            lastnodeidx = graph->head[insert[0]];

            if( mw )
            {
               assert(!SCIPisPositive(scip, graph->prize[i]));

               diff = -1.0;
               assert(stdeg != NULL);
               stdeg[lastnodeidx]++;
            }
            else
               diff = graph->cost[v->edge];

            counter = 0;

            /* try to add edges between new vertex and tree */
            for( k = 1; k < insertcount; k++ )
            {
               NODE* firstnode;
               int firstnodidx;
               SCIPlinkcuttreeEvert(v);

               /* next vertex in the current Steiner tree adjacent to vertex i resp. v (the one being scrutinized for possible insertion) */
               firstnodidx = graph->head[insert[k]];
               firstnode = &nodes[firstnodidx];

               if( mw )
               {
                  NODE* chainfirst;
                  NODE* chainlast;
                  SCIP_Real minweight;

                  assert(stdeg != NULL);

                  minweight = SCIPlinkcuttreeFindMinChain(scip, graph->prize, graph->head, stdeg, firstnode, &chainfirst, &chainlast);

                  if( SCIPisLT(scip, minweight, graph->prize[i]) )
                  {
                     assert(chainfirst != NULL && chainlast != NULL);
                     for( NODE* mynode = chainfirst; mynode != chainlast; mynode = mynode->parent )
                     {
                        int mynodeidx = graph->head[mynode->edge];
                        steinertree[mynodeidx] = FALSE;
                        stdeg[mynodeidx] = 0;
                     }

                     SCIPlinkcuttreeCut(chainfirst);
                     SCIPlinkcuttreeCut(chainlast);

                     SCIPlinkcuttreeLink(v, firstnode, insert[k]);
                     stdeg[graph->head[insert[k]]]++;

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
                  assert(stdeg != NULL);

                  SCIPlinkcuttreeEvert(v);
                  stdeg[lastnodeidx]--;
                  SCIPlinkcuttreeCut(&nodes[graph->head[insert[0]]]);
               }
               else
               {
                  steinertree[i] = TRUE;
                  newnverts++;
               }
            }
            else
            {
               if( !SCIPisNegative(scip, diff) )
               {
                  SCIPlinkcuttreeEvert(v);
                  for (k = counter - 1; k >= 0; k--)
                  {
                     SCIPlinkcuttreeCut(&nodes[graph->head[adds[k]]]);
                     SCIPlinkcuttreeEvert(&nodes[graph->tail[cuts[k]]]);
                     SCIPlinkcuttreeLink(&nodes[graph->tail[cuts[k]]],
                           &nodes[graph->head[cuts[k]]], cuts[k]);
                  }

                  /* finally, cut the edge added first (if it had been cut during the insertion process, it would have been restored above) */
                  SCIPlinkcuttreeEvert(v);
                  SCIPlinkcuttreeCut(&nodes[graph->head[insert[0]]]);
               }
               else
               {
                  SCIPlinkcuttreeEvert(&nodes[root]);
                  adds[counter] = insert[0];
                  newnode = i;
                  steinertree[i] = TRUE;
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
         SCIPfreeBufferArray(scip, &stdeg);
         SCIPfreeBufferArray(scip, &cuts2);
      }
      SCIPfreeBufferArray(scip, &cuts);
      SCIPfreeBufferArray(scip, &adds);
      SCIPfreeBufferArray(scip, &insert);

      for( e = 0; e < nedges; e++ )
         best_result[e] = UNKNOWN;

      if( newnverts > 0  )
      {
         if( mwpc )
            SCIP_CALL( SCIPStpHeurTMPrunePc(scip, graph, graph->cost, best_result, steinertree) );
         else
            SCIP_CALL( SCIPStpHeurTMPrune(scip, graph, graph->cost, 0, best_result, steinertree) );

         for( i = 0; i < nnodes; i++ )
            SCIPlinkcuttreeInit(&nodes[i]);

         /* create a link-cut tree representing the current Steiner tree */
         for( e = 0; e < nedges; e++ )
         {
            if( best_result[e] == CONNECT )
            {
               assert(steinertree[graph->tail[e]]);
               assert(steinertree[graph->head[e]]);
               SCIPlinkcuttreeLink(&nodes[graph->head[e]], &nodes[graph->tail[e]], flipedge(e));
            }
         }
         SCIPlinkcuttreeEvert(&nodes[root]);
      }
      else
      {
         SCIPlinkcuttreeEvert(&nodes[root]);
         for( i = 0; i < nnodes; i++ )
         {
            if( steinertree[i] && nodes[i].edge != -1 )
               best_result[flipedge(nodes[i].edge)] = 0;
         }
      }

#ifdef SCIP_DEBUG
      {
         SCIP_Real obj = 0.0;
         for( e = 0; e < nedges; e++ )
            obj += (best_result[e] > -1) ? graph->cost[e] : 0.0;
         printf("ObjAfterVertexInsertion=%.12e\n", obj);

      }
#endif
   }

   assert(graph_sol_valid(scip, graph, best_result));

   /* Key-Vertex Elimination & Key-Path Exchange */
   if( !mw )
   {
      IDX* blists_curr;
      IDX** blists_start;  /* array [1,..,nnodes],
                            * if node i is in the current ST, blists_start[i] points to a linked list of all nodes having i as their base */
      PATH* mst;           /* minimal spanning tree structure */
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
      STP_Bool* pinned;
      STP_Bool* scanned;
      STP_Bool* nodesmark;

#ifdef SCIP_DEBUG
      SCIP_Real obj = 0.0;
      for( e = 0; e < nedges; e++ )
         obj += (best_result[e] > -1) ? graph->cost[e] : 0.0;
      printf(" ObjBEFKEYVertexELimination=%.12e\n", obj);
#endif

      for( k = 0; k < nnodes; k++ )
         graphmark[k] = (graph->grad[k] > 0);

      graphmark[root] = TRUE;

      /* allocate memory */

      /* only needed for Key-Path Elimination */
      SCIP_CALL( SCIPallocBufferArray(scip, &newedges, nedges) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lvledges_start, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &boundedges, nedges) );
      SCIP_CALL( SCIPallocBufferArray(scip, &supernodesid, nnodes) );

      /* only needed for Key-Path Exchange */

      /* memory needed for both Key-Path Elimination and Exchange */
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

      /* initialize data structures */
      SCIP_CALL( SCIPStpunionfindInit(scip, &uf, nnodes) );

      for( nruns = 0; nruns < LOCAL_MAXRESTARTS && localmoves > 0; nruns++ )
      {
         localmoves = 0;

         BMSclearMemoryArray(blists_start, nnodes);

         /* find a DFS order of the ST nodes */
         nstnodes = 0;
         dfsorder(graph, best_result, &(root), &nstnodes, dfstree);
         assert(root == graph->source);

         /* compute a voronoi diagram with the ST nodes as bases */
         graph_voronoi(scip, graph, graph->cost, graph->cost, steinertree, vbase, vnoi);

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
            for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
            {
               k = graph->head[e];
               if( Is_term(graph->term[k]) )
               {
                  graphmark[k] = FALSE;
                  for( l = graph->outbeg[k]; l != EAT_LAST; l = graph->oeat[l] )
                  {
                     if( !Is_term(graph->term[graph->head[l]]) )
                     {
                        assert(graph->head[l] != root);
                        pinned[graph->head[l]] = TRUE;
                     }
                  }
               }
            }
            if( probtype != STP_RPCSPG )
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
         SCIP_CALL( lca(scip, graph, root, &uf, nodesmark, best_result, lvledges_start, boundpaths, heapsize, vbase) );

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
            if( !pinned[crucnode] && !Is_term(graph->term[crucnode]) && nodeIsCrucial(graph, best_result, crucnode) )
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
                  if( (best_result[edge] > -1 && steinertree[graph->head[edge]]) || (best_result[flipedge(edge)] > -1 && steinertree[graph->tail[edge]]) )
                  {
                     kpcost += graph->cost[edge];

                     /* check whether the current edge leads to the ST root*/
                     if( best_result[flipedge(edge)] > -1 )
                     {
                        k = flipedge(edge);
                        kpedges[nkpedges++] = k;
                        assert( edge == nodes[crucnode].edge );
                     }
                     else
                     {
                        kpedges[nkpedges++] = edge;
                        adjnode = graph->head[edge];
                        e = edge;

                        /* move along the key-path until its end (i.e. a crucial or pinned node) is reached */
                        while( !pinned[adjnode] && !nodeIsCrucial(graph, best_result, adjnode) && steinertree[adjnode] )
                        {
                           /* update the union-find data structure */
                           SCIPStpunionfindUnion(&uf, crucnode, adjnode, FALSE);

                           kpnodes[nkpnodes++] = adjnode;

                           for( e = graph->outbeg[adjnode]; e != EAT_LAST; e = graph->oeat[e] )
                           {
                              if( best_result[e] > -1 )
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
                        if( !steinertree[adjnode] )
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
                  while( !pinned[adjnode] && !nodeIsCrucial(graph, best_result, adjnode) && steinertree[adjnode] )
                  {
                     /* update the union-find data structure */
                     kpnodes[nkpnodes++] = adjnode;

                     for( e = graph->inpbeg[adjnode]; e != EAT_LAST; e = graph->ieat[e] )
                     {
                        if( best_result[e] > -1 )
                        {
                           assert(steinertree[graph->tail[e]]);
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
                     assert( (vbase[graph->tail[edge]] == UNKNOWN)? UNKNOWN : SCIPStpunionfindFind(&uf, vbase[graph->tail[edge]]) == l );

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

               /* compute the cost of the MST */
               for( l = 0; l < nsupernodes - 1; l++ )
               {
                  /* compute the edge in the original graph corresponding to the current MST edge */
                  if( mst[l].edge % 2  == 0 )
                     edge = boundedges[mst[l].edge / 2 ];
                  else
                     edge = flipedge(boundedges[mst[l].edge / 2 ]);

                  mstcost += graph->cost[edge];
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
                     }
                  }
               }

               if( SCIPisLT(scip, mstcost, kpcost) )
               {
                  localmoves++;
                  SCIPdebugMessage("found improving solution in KEY VERTEX ELIMINATION (round: %d) \n ", nruns);

                  /* unmark the original edges spanning the supergraph */
                  for( e = 0; e < nkpedges; e++ )
                  {
                     assert(best_result[kpedges[e]] != -1);
                     best_result[kpedges[e]] = -1;
                  }

                  /* mark all ST nodes except for those belonging to the root-component as forbidden */
                  for( k = rootpathstart; k < nkpnodes; k++ )
                  {
                     graphmark[kpnodes[k]] = FALSE;
                     steinertree[kpnodes[k]] = FALSE;
                  }

                  for( k = 0; k < i; k++ )
                  {
                     node = SCIPStpunionfindFind(&uf, dfstree[k]);
                     if( nodesmark[node] || node == crucnode )
                     {
                        graphmark[dfstree[k]] = FALSE;
                        steinertree[dfstree[k]] = FALSE;
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
                           e = nodes[node].edge;

                           assert(best_result[e] == -1 && best_result[flipedge(e)] != -1 );
                           best_result[e] = CONNECT;
                           best_result[flipedge(e)] = UNKNOWN;
                           node = graph->head[e];
                        }
                     }

                     /* is the vbase of the current boundary-edge tail in the root-component? */
                     if( !nodesmark[SCIPStpunionfindFind(&uf, vbase[graph->tail[edge]])] )
                     {

                        best_result[edge] = CONNECT;

                        for( node = graph->tail[edge], adjnode = graph->head[edge]; node != vbase[node]; adjnode = node, node = graph->tail[vnoi[node].edge] )
                        {
                           graphmark[node] = FALSE;

                           if( best_result[flipedge(vnoi[node].edge)] == CONNECT )
                              best_result[flipedge(vnoi[node].edge)] = UNKNOWN;

                           best_result[vnoi[node].edge] = CONNECT;
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
                              if( best_result[edge] == CONNECT && graphmark[adjnode] && steinertree[adjnode]  && SCIPStpunionfindFind(&uf, adjnode) != node )
                              {

                                 assert(scanned[adjnode]);
                                 /* meld the heaps */
                                 SCIPpairheapMeldheaps(scip, &boundpaths[node], &boundpaths[adjnode], &heapsize[node], &heapsize[adjnode]);

                                 /* update the union-find data structure */
                                 SCIPStpunionfindUnion(&uf, node, adjnode, FALSE);

                                 /* move along the key-path until its end (i.e. until a crucial node is reached) */
                                 while( !nodeIsCrucial(graph, best_result, adjnode) && !pinned[adjnode] )
                                 {
                                    for( e = graph->outbeg[adjnode]; e != EAT_LAST; e = graph->oeat[e] )
                                    {
                                       if( best_result[e] != -1 )
                                          break;
                                    }

                                    /* assert that each leaf of the ST is a terminal */
                                    assert( e != EAT_LAST );
                                    adjnode = graph->head[e];
                                    if( !steinertree[adjnode]  )
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
                           if( best_result[vnoi[node].edge] == CONNECT )
                              best_result[vnoi[node].edge] = -1;

                           best_result[flipedge(vnoi[node].edge)] = CONNECT;

                        }
                     }
                     else
                     {

                        best_result[edge] = CONNECT;

                        for( node = graph->tail[edge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                        {
                           graphmark[node] = FALSE;
                           if( best_result[vnoi[node].edge] != CONNECT && best_result[flipedge(vnoi[node].edge)] != CONNECT )
                              best_result[vnoi[node].edge] = CONNECT;

                        }

                        for( node = graph->head[edge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                        {
                           graphmark[node] = FALSE;

                           best_result[flipedge(vnoi[node].edge)] = CONNECT;
                           best_result[vnoi[node].edge] = UNKNOWN;
                        }
                     }
                  }

                  for( k = 0; k < nkpnodes; k++ )
                  {
                     assert(graphmark[kpnodes[k]] == FALSE);
                     assert(steinertree[kpnodes[k]] == FALSE);
                  }
                  assert(!graphmark[crucnode]);
               }
               else
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
               if( (!nodeIsCrucial(graph, best_result, crucnode) && !pinned[crucnode]) )
                  continue;

               if( Is_term(graph->term[crucnode]) || pinned[crucnode] )
               {
                  for( edge = graph->outbeg[crucnode]; edge != EAT_LAST; edge = graph->oeat[edge] )
                  {
                     adjnode = graph->head[edge];
                     /* check whether edge 'edge' leads to an ancestor of terminal 'crucnode' */
                     if( best_result[edge] == CONNECT && steinertree[adjnode] && graphmark[adjnode] )
                     {
                        assert( SCIPStpunionfindFind(&uf, adjnode) != crucnode);
                        assert(scanned[adjnode]);
                        /* meld the heaps */
                        SCIPpairheapMeldheaps(scip, &boundpaths[crucnode], &boundpaths[adjnode], &heapsize[crucnode], &heapsize[adjnode]);

                        /* update the union-find data structure */
                        SCIPStpunionfindUnion(&uf, crucnode, adjnode, FALSE);

                        /* move along the key-path until its end (i.e. until a crucial node is reached) */
                        while( !nodeIsCrucial(graph, best_result, adjnode) && !pinned[adjnode] )
                        {
                           for( e = graph->outbeg[adjnode]; e != EAT_LAST; e = graph->oeat[e] )
                           {
                              if( best_result[e] != -1 )
                                 break;
                           }

                           /* assert that each leaf of the ST is a terminal */
                           assert( e != EAT_LAST );
                           adjnode = graph->head[e];
                           if( !steinertree[adjnode] || !graphmark[adjnode] )
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
               kptailnode = graph->head[nodes[crucnode].edge];
               kpathcost = graph->cost[nodes[crucnode].edge];

               while( !nodeIsCrucial(graph, best_result, kptailnode) && !pinned[kptailnode] )
               {
                  kpathcost += graph->cost[nodes[kptailnode].edge];

                  kpnodes[nkpnodes++] = kptailnode;
                  kptailnode = graph->head[nodes[kptailnode].edge];
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
                  newedge = nodes[crucnode].edge;

               if( oldedge != UNKNOWN && newedge != UNKNOWN && SCIPisLT(scip, edgecost,
                     vnoi[graph->tail[newedge]].dist + graph->cost[newedge] + vnoi[graph->head[newedge]].dist) )
                  newedge = oldedge;
               if( oldedge != UNKNOWN && newedge == UNKNOWN )
                  newedge = oldedge;

               assert( newedge != UNKNOWN );
               edgecost = vnoi[graph->tail[newedge]].dist + graph->cost[newedge] + vnoi[graph->head[newedge]].dist;
               if( SCIPisLT(scip, edgecost, kpathcost) )
               {
                  node = SCIPStpunionfindFind(&uf, vbase[graph->head[newedge]]);

                  SCIPdebugMessage( "ADDING NEW KEY PATH (%f )\n", edgecost - kpathcost );

                  localmoves++;

                  /* remove old keypath */
                  assert(  best_result[flipedge(nodes[crucnode].edge)] != UNKNOWN );
                  best_result[flipedge(nodes[crucnode].edge)] = UNKNOWN;
                  steinertree[crucnode] = FALSE;
                  graphmark[crucnode] = FALSE;

                  for( k = 0; k < nkpnodes; k++ )
                  {
                     assert(  best_result[flipedge(nodes[kpnodes[k]].edge)] != UNKNOWN );
                     best_result[flipedge(nodes[kpnodes[k]].edge)] = UNKNOWN;
                     steinertree[kpnodes[k]] = FALSE;
                     graphmark[kpnodes[k]] = FALSE;
                  }
                  assert(graphmark[kptailnode]);

                  if( node == crucnode )
                     newedge = flipedge(newedge);

                  for( node = graph->tail[newedge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                  {
                     graphmark[node] = FALSE;

                     best_result[flipedge(vnoi[node].edge)] = CONNECT;
                     best_result[vnoi[node].edge] = UNKNOWN;
                  }

                  for( node = graph->head[newedge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                  {
                     graphmark[node] = FALSE;

                     best_result[vnoi[node].edge] = CONNECT;
                  }

                  best_result[flipedge(newedge)] = CONNECT;

                  newpathend = vbase[graph->tail[newedge]];
                  assert(node == vbase[graph->head[newedge]] );

                  /* flip all edges on the ST path between the endnode of the new key-path and the current crucial node */
                  k = newpathend;
                  assert(SCIPStpunionfindFind(&uf, newpathend) == crucnode);

                  while( k != crucnode )
                  {
                     assert(graphmark[k]);
                     assert( best_result[flipedge(nodes[k].edge)] != -1);
                     best_result[flipedge(nodes[k].edge)] = UNKNOWN;

                     best_result[nodes[k].edge] = CONNECT;

                     k = graph->head[nodes[k].edge];
                  }

                  for( k = 0; k < i; k++ )
                  {
                     if( crucnode == SCIPStpunionfindFind(&uf, dfstree[k]) )
                     {
                        graphmark[dfstree[k]] = FALSE;
                        steinertree[dfstree[k]] = FALSE;
                     }
                  }

                  /* update union find */
                  if( !Is_term(graph->term[node]) && scanned[node] && !pinned[node] && SCIPStpunionfindFind(&uf, node) == node )
                  {
                     for( edge = graph->outbeg[node]; edge != EAT_LAST; edge = graph->oeat[edge] )
                     {
                        adjnode = graph->head[edge];
                        /* check whether edge 'edge' leads to an ancestor of terminal 'node' */
                        if( best_result[edge] == CONNECT && steinertree[adjnode]  && graphmark[adjnode] && SCIPStpunionfindFind(&uf, adjnode) != node )
                        {
                           assert(scanned[adjnode]);
                           /* meld the heaps */
                           SCIPpairheapMeldheaps(scip, &boundpaths[node], &boundpaths[adjnode], &heapsize[node], &heapsize[adjnode]);

                           /* update the union-find data structure */
                           SCIPStpunionfindUnion(&uf, node, adjnode, FALSE);

                           /* move along the key-path until its end (i.e. until a crucial node is reached) */
                           while( !nodeIsCrucial(graph, best_result, adjnode) && !pinned[adjnode] )
                           {
                              for( e = graph->outbeg[adjnode]; e != EAT_LAST; e = graph->oeat[e] )
                              {
                                 if( best_result[e] != -1 )
                                    break;
                              }

                              /* assert that each leaf of the ST is a terminal */
                              assert( e != EAT_LAST );
                              adjnode = graph->head[e];
                              if( !steinertree[adjnode]  )
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
               steinertree[i] = FALSE;
               graphmark[i] = (graph->grad[i] > 0);
               SCIPlinkcuttreeInit(&nodes[i]);
            }

            graphmark[root] = TRUE;

            /* create a link-cut tree representing the current Steiner tree */
            for( e = 0; e < nedges; e++ )
            {
               assert(graph->head[e] == graph->tail[flipedge(e)]);

               /* if edge e is in the tree, so are its incident vertices */
               if( best_result[e] != -1 )
               {
                  steinertree[graph->tail[e]] = TRUE;
                  steinertree[graph->head[e]] = TRUE;
                  SCIPlinkcuttreeLink(&nodes[graph->head[e]], &nodes[graph->tail[e]], flipedge(e));
               }
            }
            assert( nodes[root].edge == -1 );
            nodes[root].edge = -1;
         }
      }

      /* free data structures */
      SCIPStpunionfindFree(scip, &uf);
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
      SCIPfreeBufferArray(scip, &supernodesid);
      SCIPfreeBufferArray(scip, &boundedges);
      SCIPfreeBufferArray(scip, &lvledges_start);
      SCIPfreeBufferArray(scip, &newedges);
      /******/
   }

#ifdef SCIP_DEBUG
   {
      SCIP_Real obj = 0.0;
      for( e = 0; e < nedges; e++ )
         obj += (best_result[e] > -1) ? graph->cost[e] : 0.0;

      printf(" ObjAfterHeurLocal=%.12e\n", obj);
   }
#endif

#ifndef NDEBUG
   assert(SCIPisLE(scip, graph_sol_getObj(graph->cost, best_result, 0.0, nedges), initialobj));
#endif

   SCIPfreeBufferArray(scip, &steinertree);
   SCIPfreeBufferArray(scip, &nodes);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &vnoi);

   return SCIP_OKAY;
}

/** Greedy Extension local heuristic for (R)PC and MW */
SCIP_RETCODE SCIPStpHeurLocalExtendPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge cost array*/
   PATH*                 path,               /**< shortest data structure array */
   int*                  stedge,             /**< initialized array to indicate whether an edge is part of the Steiner tree */
   STP_Bool*             stvertex,           /**< uninitialized array to indicate whether a vertex is part of the Steiner tree */
   SCIP_Bool*            extensions          /**< pointer to store whether extensions have been made */
   )
{
   GNODE candidates[MAX(GREEDY_EXTENSIONS, GREEDY_EXTENSIONS_MW)];
   int candidatesup[MAX(GREEDY_EXTENSIONS, GREEDY_EXTENSIONS_MW)];

   PATH* orgpath;
   SCIP_PQUEUE* pqueue;
   SCIP_Real bestsolval;
   int i;
   int root;
   int nedges;
   int nnodes;
   int nextensions;
   int greedyextensions;
   STP_Bool* stvertextmp;

   assert(scip != NULL);
   assert(path != NULL);
   assert(graph != NULL);
   assert(stedge != NULL);
   assert(cost != NULL);
   assert(stvertex != NULL);

   root = graph->source;
   nnodes = graph->knots;
   nedges = graph->edges;

   graph_pc_2transcheck(graph);
   SCIP_CALL( SCIPallocBufferArray(scip, &stvertextmp, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orgpath, nnodes) );

   /* initialize solution vertex array with FALSE */
   BMSclearMemoryArray(stvertex, nnodes);

   stvertex[graph->source] = TRUE;
   path[graph->source].edge = UNKNOWN;

   for( int e = 0; e < nedges; e++ )
      if( stedge[e] == CONNECT )
      {
         i = graph->head[e];
         path[i].edge = e;
         stvertex[i] = TRUE;
      }

   graph_path_st_pcmw_extend(scip, graph, cost, path, stvertex, extensions);

   BMScopyMemoryArray(orgpath, path, nnodes);

   if( graph->stp_type == STP_MWCSP )
      greedyextensions = GREEDY_EXTENSIONS_MW;
   else
      greedyextensions = GREEDY_EXTENSIONS;

   /*** compute solution value and save greedyextensions many best unconnected nodes  ***/

   SCIP_CALL( SCIPpqueueCreate(&pqueue, greedyextensions, 2.0, GNODECmpByDist) );

   bestsolval = 0.0;
   nextensions = 0;
   for( i = 0; i < nnodes; i++ )
   {
      if( graph->grad[i] == 0 )
         continue;

      if( Is_term(graph->term[i]) )
      {
         if( i != root )
         {
            int e;
            SCIP_Bool connect = FALSE;
            SCIP_Real ecost = -1.0;
            for( e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
            {
               if( root == graph->tail[e] )
                  ecost = graph->cost[e];
               else if( stvertex[graph->tail[e]] )
                  connect = TRUE;
            }

            if( !connect )
               bestsolval += ecost;
         }
      }
      else if( stvertex[i] )
      {
         bestsolval += graph->cost[orgpath[i].edge];
      }
      else if( orgpath[i].edge != UNKNOWN && Is_pterm(graph->term[i]) )
      {
         if( nextensions < greedyextensions )
         {
            candidates[nextensions].dist = graph->prize[i] - path[i].dist;
            candidates[nextensions].number = i;

            SCIP_CALL( SCIPpqueueInsert(pqueue, &(candidates[nextensions++])) );
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
      }
   }

   for( int restartcount = 0; restartcount < GREEDY_MAXRESTARTS;  restartcount++ )
   {
      int l = 0;
      SCIP_Bool extensionstmp = FALSE;

      i = nextensions;

      /* write extension candidates into array, from max to min */
      while( SCIPpqueueNElems(pqueue) > 0 )
      {
         GNODE* min = (GNODE*) SCIPpqueueRemove(pqueue);
         assert(i > 0);
         candidatesup[--i] = min->number;
      }
      assert(i == 0);

      /* iteratively insert new subpaths and try to improve solution */
      for( ; l < nextensions; l++ )
      {
         i = candidatesup[l];
         if( !stvertex[i] )
         {
            SCIP_Real newsolval = 0.0;
            int k = i;

            BMScopyMemoryArray(stvertextmp, stvertex, nnodes);
            BMScopyMemoryArray(path, orgpath, nnodes);

            /* add new extension */
            while( !stvertextmp[k] )
            {
               stvertextmp[k] = TRUE;
               assert(orgpath[k].edge != UNKNOWN);
               k = graph->tail[orgpath[k].edge];
            }

            graph_path_st_pcmw_extend(scip, graph, cost, path, stvertextmp, &extensionstmp);

            for( int j = 0; j < nnodes; j++ )
            {
               if( graph->grad[j] == 0 )
                  continue;

               if( Is_term(graph->term[j]) )
               {
                  if( j != root )
                  {
                     int e;
                     SCIP_Bool connect = FALSE;
                     SCIP_Real ecost = -1.0;
                     for( e = graph->inpbeg[j]; e != EAT_LAST; e = graph->ieat[e] )
                     {
                        if( root == graph->tail[e] )
                           ecost = graph->cost[e];
                        else if( stvertextmp[graph->tail[e]] )
                           connect = TRUE;
                     }

                     if( !connect )
                        newsolval += ecost;
                  }
               }
               else if( stvertextmp[j] )
               {
                  newsolval += graph->cost[path[j].edge];
               }
            }

            /* new solution value better than old one? */
            if( SCIPisLT(scip, newsolval, bestsolval) )
            {
               *extensions = TRUE;

               bestsolval = newsolval;
               BMScopyMemoryArray(stvertex, stvertextmp, nnodes);

               /* save greedyextensions many best unconnected nodes  */
               nextensions = 0;

               for( int j = 0; j < nnodes; j++ )
               {
                  orgpath[j].edge = path[j].edge;
                  orgpath[j].dist = path[j].dist;
                  if( !stvertex[j] && Is_pterm(graph->term[j]) && path[j].edge != UNKNOWN )
                  {
                     if( nextensions < greedyextensions )
                     {
                        candidates[nextensions].dist = graph->prize[j] - path[j].dist;
                        candidates[nextensions].number = j;

                        SCIP_CALL( SCIPpqueueInsert(pqueue, &(candidates[nextensions++])) );
                     }
                     else
                     {
                        /* get candidate vertex of minimum value */
                        GNODE* min = (GNODE*) SCIPpqueueFirst(pqueue);
                        if( SCIPisLT(scip, min->dist, graph->prize[j] - path[j].dist) )
                        {
                           min = (GNODE*) SCIPpqueueRemove(pqueue);
                           min->dist = graph->prize[j] - path[j].dist;
                           min->number = j;
                           SCIP_CALL( SCIPpqueueInsert(pqueue, min) );
                        }
                     }
                  }
               }

               break;
            } /* if new solution value better than old one? */
         } /* if !stvertex[i] */
      } /* for l < nextension */
      /* no more extensions performed? */
      if( l == nextensions )
         break;
   } /* main loop */

   /* have vertices been added? */
   if( *extensions )
   {
      SCIPdebugMessage("SCIPStpHeurLocalExtendPcMw found extensions \n");
      for( int e = 0; e < nedges; e++ )
         stedge[e] = UNKNOWN;
      SCIP_CALL( SCIPStpHeurTMPrunePc(scip, graph, cost, stedge, stvertex) );
   }

   SCIPpqueueFree(&pqueue);
   SCIPfreeBufferArray(scip, &orgpath);
   SCIPfreeBufferArray(scip, &stvertextmp);

#ifdef SCIP_DEBUG
   {
      SCIP_Real t = 0.0;
      for( int e = 0; e < nedges; e++ )
         if( stedge[e] == CONNECT )
            t += graph->cost[e];

      SCIPdebugMessage("SCIPStpHeurLocalExtendPcMw: exit real cost %f \n", t);
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
   int i;
   int e;
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
         for( i = min - 1; i >= v + 1; i-- )
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

   /* for PC variants: test whether solution is trivial */
   if( graph->stp_type == STP_PCSPG || graph->stp_type == STP_RPCSPG || graph->stp_type == STP_MWCSP )
   {
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
   for( e = 0; e < nedges; e++ )
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

      for( v = nnodes - 1; v >= 0; --v )
         steinertree[v] = FALSE;

      for( e = nedges - 1; e >= 0; --e )
      {
         if( results[e] == CONNECT )
         {
            steinertree[graph->tail[e]] = TRUE;
            steinertree[graph->head[e]] = TRUE;
         }
         results[e] = UNKNOWN;
      }

      if( graph->stp_type == STP_PCSPG || graph->stp_type == STP_RPCSPG || graph->stp_type == STP_MWCSP )
         SCIP_CALL( SCIPStpHeurTMPrunePc(scip, graph, graph->cost, results, steinertree) );
      else
         SCIP_CALL( SCIPStpHeurTMPrune(scip, graph, graph->cost, 0, results, steinertree) );

      SCIPfreeBufferArray(scip, &steinertree);
   }

   /* execute local heuristics */
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, graph->cost, results) );

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
