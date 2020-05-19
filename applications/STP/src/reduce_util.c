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

/**@file   reduce_util.c
 * @brief  utility methods for Steiner tree reductions
 * @author Daniel Rehfeldt
 *
 * This file implements utility methods for Steiner tree problem reduction techniques.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG
#include "reduce.h"
#include "portab.h"

#define STP_TPATHS_NTERMBASES 4


/** storage for edge on complete graph */
typedef struct complete_edge
{
   int                   tail;              /**< tail vertex */
   int                   head;              /**< head vertex */
   SCIP_Real             cost;              /**< edge cost */
} CEDGE;


/** lightweight minimum spanning tree structure that allows to add vertices to given MST on complete graph (in CSR format) */
struct dynamic_complete_minimum_spanning_tree
{
   CEDGE*                edgestore;         /**< storage for edges (of size maxnnodes) */
   SCIP_Real*            adjcost_buffer;    /**< distances buffer (of size maxnnodes) */
   SCIP_Bool*            nodemark;          /**< array for marking nodes (of size maxnnodes) */
   int                   maxnnodes;         /**< maximum number of nodes that can be handled */
};


/** Steiner nodes to terminal paths
 * NOTE: all arrays are of size STP_TPATHS_NTERMBASES * nnodes */
struct nodes_to_terminal_paths
{
   PATH*                 termpaths;         /**< path data (leading to first, second, ... terminal) */
   int*                  termbases;         /**< terminals to each non terminal */
   int*                  state;             /**< array to mark the state of each node during calculation */
   int                   nnodes;            /**< number of nodes of underlying graph */
};


/** see reduce.h */
struct node_one_hop_star
{
   int*                  edgeId;            /**< IDs for each adjacent edge of current node (of size maxNodeDegree) */
   int*                  edgesSelected;     /**< list of currently selected edges (of size maxNodeDegree) */
   int*                  edgesSelectedPos;  /**< list of position of currently selected edges w.r.t. edgeId (of size maxNodeDegree) */
   int                   nodeDegree;        /**< degree of current node */
   int                   starDegree;        /**< degree of current star */
   int                   maxNodeDegree;     /**< maximum allowed node degree */
   int                   starcenter;        /**< node for which the star is created */
   SCIP_Bool             allStarsChecked;   /**< have all stars been checked? */
};


/** see reduce.h */
struct special_distance_graph
{
   GRAPH*                distgraph;          /**< (complete) distance graph */
   PATH*                 sdmst;              /**< MST on sdgraph */
   SCIP_Real*            mstcosts;           /**< maximum MST edge costs in descending order */
   int*                  nodes_id;           /**< number of each node in original graph */
   STP_Bool*             halfedge_isInMst;   /**< signifies whether edge of original graph is part of MST
                                                  NOTE: operates on edges / 2! */
   SCIP_Real             mstmaxcost;         /**< maximum edge cost */
   int                   nnodesorg;          /**< number of nodes of original graph */
   int                   nedgesorg;          /**< number of edges of original graph */
   SCIP_Bool             mstcostsReady;      /**< are the mstcosts already available? */
};


/** allocates TPATHS data */
static
SCIP_RETCODE tpathsAlloc(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   TPATHS**              tpaths              /**< the terminal paths */
)
{

   TPATHS* tp;
   const int nnodes = graph_get_nNodes(g);
   assert(nnodes >= 1);

   SCIP_CALL( SCIPallocMemory(scip, tpaths) );

   tp = *tpaths;
   tp->nnodes = nnodes;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(tp->termpaths), nnodes * STP_TPATHS_NTERMBASES) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(tp->termbases), nnodes * STP_TPATHS_NTERMBASES) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(tp->state), nnodes * STP_TPATHS_NTERMBASES) );

   return SCIP_OKAY;
}


/** allocates TPATHS data */
static
SCIP_RETCODE tpathsBuild(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph */
   TPATHS*               tpaths              /**< the terminal paths */
)
{
   int* heap;
   const int nnodes = graph_get_nNodes(g);

   assert(nnodes >= 1);
   assert(STP_TPATHS_NTERMBASES == 4);

   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1)  );
   graph_get4nextTerms(scip, g, g->cost, g->cost, tpaths->termpaths, tpaths->termbases, heap, tpaths->state);
   SCIPfreeBufferArray(scip, &heap);

   return SCIP_OKAY;
}

/** recursive method for adding node to MST */
static
void dcmstInsert(
   const CSR*            org_mst,            /**< the base MST */
   const SCIP_Real       adjcosts[],         /**< (undirected) adjacency costs for new node */
   int                   root,               /**< the current root */
   CEDGE                 new_mst[],          /**< new MST */
   SCIP_Bool             new_nodemarked[],   /**< array */
   CEDGE*                max_path_edge,      /**< pointer to maximum edge on path to new node */
   int*                  new_nedges          /**< pointer to current number of edges */
)
{
   CEDGE root2new = { .tail = root, .head = org_mst->nnodes, .cost = adjcosts[root] };
   const int* const org_start = org_mst->start;
   const int* const org_head = org_mst->head;
   const SCIP_Real* const org_cost = org_mst->cost;

   assert(new_nodemarked[root]);

   /* visit all neighbors or root in the original MST */
   for( int i = org_start[root]; i != org_start[root + 1]; ++i )
   {
      const int w = org_head[i];

      /* node not visited yet? */
      if( !new_nodemarked[w] )
      {
         const SCIP_Real costroot2w = org_cost[i];

         new_nodemarked[w] = TRUE;
         dcmstInsert(org_mst, adjcosts, w, new_mst, new_nodemarked, max_path_edge, new_nedges);

         assert(max_path_edge->tail >= 0);
         assert(*new_nedges >= 0 && *new_nedges < org_mst->nnodes);

         if( max_path_edge->cost < costroot2w )
         {
            new_mst[(*new_nedges)++] = *max_path_edge;

            if( costroot2w < root2new.cost )
            {
               root2new.tail = root;
               root2new.head = w;
               root2new.cost = costroot2w;
            }
         }
         else
         {
            const int nedges = (*new_nedges);

            new_mst[nedges].tail = root;
            new_mst[nedges].head = w;
            new_mst[nedges].cost = costroot2w;

            (*new_nedges)++;

            if( max_path_edge->cost < root2new.cost )
            {
               root2new = *max_path_edge;
            }
         }
      }
   }

   *max_path_edge = root2new;
}


/** add node to MST */
static inline
void dcmstAddNode(
   const CSR*            mst_in,             /**< source */
   const SCIP_Real*      adjcosts,           /**< (undirected) adjacency costs for new node */
   DCMST*                dmst                /**< underlying structure */
)
{
   CEDGE max_path_edge = { .tail = -1, .head = -1, .cost = -1.0 };
   CEDGE* const edgestore = dmst->edgestore;
   SCIP_Bool* const nodemark = dmst->nodemark;
   int nedges_new = 0;
   const int nnodes_in = mst_in->nnodes;

   assert(nnodes_in >= 1);

   nodemark[0] = TRUE;

   for( int i = 1; i < nnodes_in; ++i )
      nodemark[i] = FALSE;

   dcmstInsert(mst_in, adjcosts, 0, edgestore, nodemark, &max_path_edge, &nedges_new);

   assert(nedges_new == nnodes_in - 1);

   edgestore[nedges_new] = max_path_edge;
}


/** transforms edge-store to CSR  */
static inline
void dcmstGetCSRfromStore(
   const DCMST*          dmst,               /**< underlying structure */
   CSR*                  mst_out             /**< target */
)
{
   const CEDGE* const edgestore = dmst->edgestore;
   int* const mst_start = mst_out->start;
   int* const mst_head = mst_out->head;
   SCIP_Real* const mst_cost = mst_out->cost;
   const int mst_nnodes = mst_out->nnodes;

   /* undirected edges */
   const int mst_nedges = mst_nnodes - 1;

   assert(mst_nnodes <= dmst->maxnnodes);
   assert(2 * mst_nedges == mst_out->nedges_max);

   BMSclearMemoryArray(mst_start, mst_nnodes + 1);

   for( int i = 0; i < mst_nedges; ++i )
   {
      const int v1 = edgestore[i].tail;
      const int v2 = edgestore[i].head;

      assert(v1 >= 0 && v1 < mst_nnodes);
      assert(v2 >= 0 && v2 < mst_nnodes);

      mst_start[v1]++;
      mst_start[v2]++;
   }

   assert(mst_start[mst_nnodes] == 0);

   for( int i = 1; i <= mst_nnodes; ++i )
   {
      mst_start[i] += mst_start[i - 1];
   }

   assert(mst_start[mst_nnodes] == mst_out->nedges_max);

   for( int i = 0; i < mst_nedges; ++i )
   {
      const int v1 = edgestore[i].tail;
      const int v2 = edgestore[i].head;
      const SCIP_Real cost = edgestore[i].cost;

      assert(mst_start[v1] >= 1);
      assert(mst_start[v2] >= 1);

      mst_head[--mst_start[v1]] = v2;
      mst_cost[mst_start[v1]] = cost;

      mst_head[--mst_start[v2]] = v1;
      mst_cost[mst_start[v2]] = cost;
   }
}


/** gets weight of MST from DCMST store */
static inline
SCIP_Real dcmstGetWeightFromStore(
   SCIP*                 scip,               /**< SCIP */
   int                   mst_nedges,         /**< number of edges */
   const DCMST*          dmst                /**< underlying structure */
)
{
   const CEDGE* const edgestore = dmst->edgestore;
   SCIP_Real weight = 0.0;

   assert(mst_nedges < dmst->maxnnodes);

   for( int i = 0; i < mst_nedges; ++i )
   {
#ifndef NDEBUG
      const int v1 = edgestore[i].tail;
      const int v2 = edgestore[i].head;

      assert(v1 >= 0 && v1 < dmst->maxnnodes);
      assert(v2 >= 0 && v2 < dmst->maxnnodes);
      assert(GE(dmst->edgestore[i].cost, 0.0));
      assert(LE(dmst->edgestore[i].cost, FARAWAY));
#endif

      weight += edgestore[i].cost;
   }

   return weight;
}


/** sets star position array to initial setting for current star degree */
static inline
void starSelectedPositionsReset(
   STAR*                 star                /**< the star */
)
{
   int* const edgesSelectedPos =  star->edgesSelectedPos;
   const int starDegree = star->starDegree;

   for( int i = 0; i < starDegree; i++ )
   {
      edgesSelectedPos[i] = i;
   }
}


/** fills array star->edgesSelected by using the current positions */
static inline
void starSelectedEdgesUpdate(
   STAR*                 star                /**< the star (in/out) */
)
{
   int* const edgesSelected = star->edgesSelected;
   const int* const edgesSelectedPos = star->edgesSelectedPos;
   const int* const edgeId = star->edgeId;
   const int starDegree = star->starDegree;

   assert(starDegree >= 3);

   for( int i = 0; i < starDegree; i++ )
   {
      const int pos = edgesSelectedPos[i];
      edgesSelected[i] = edgeId[pos];
   }
}


/** moves to next star */
static inline
void starSelectedPositionsSetNext(
   STAR*                 star                /**< the star (in/out) */
)
{
   int pos;
   const int nodeDegree = star->nodeDegree;
   const int starDegree = star->starDegree;
   int* const edgesSelectedPos = star->edgesSelectedPos;

   assert(3 <= starDegree && starDegree <= nodeDegree);

   /* all current positions are stored in edgesSelectedPos[0,...,starDegree-1] */

   /* check for each position, bottom-up, whether it can be increased without it hitting the border */
   for( pos = starDegree - 1; pos >= 0; pos-- )
   {
      const SCIP_Bool isLastPos = (pos == (starDegree - 1));
      const int border = isLastPos ? nodeDegree : edgesSelectedPos[pos + 1];

      /* still space? */
      if( edgesSelectedPos[pos] < border - 1 )
      {
         break;
      }
   }

   if( pos >= 0 )
   {
      edgesSelectedPos[pos]++;

      /* adapt all following positions */
      for( int i = pos + 1; i < starDegree; i++ )
      {
         edgesSelectedPos[i] = edgesSelectedPos[i - 1] + 1;
      }
   }
   else
   {
      assert(pos == -1);
      star->starDegree--;
      starSelectedPositionsReset(star);
   }
}


/** have we just move past the last star? */
static inline
SCIP_Bool starIsDeg2(
   const STAR*           star                /**< the star (in/out) */
)
{
   if( star->starDegree <= 2 )
   {
      assert(star->starDegree == 2);
      return TRUE;
   }

   return FALSE;
}


/** allocates memory */
static
SCIP_RETCODE sdgraphAlloc(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   const int nterms = g->terms;
   SDGRAPH* g_sd;

   SCIP_CALL( SCIPallocMemory(scip, sdgraph) );
   g_sd = *sdgraph;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->sdmst), nterms) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->mstcosts), nterms) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->nodes_id), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->halfedge_isInMst), nedges / 2) );

   g_sd->nnodesorg = nnodes;
   g_sd->nedgesorg = nedges;
   g_sd->mstcostsReady = FALSE;

   return SCIP_OKAY;
}


/** builds distance graph */
static
SCIP_RETCODE sdgraphBuildDistgraph(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH*              g_sd,               /**< the SD graph */
   VNOI**                vnoi,               /**< Voronoi */
   int**                 distedge2org        /**< array of size nedges / 2 */
)
{
   GRAPH* distgraph;
   int* RESTRICT distnodes_id = g_sd->nodes_id;
   SCIP_Real* RESTRICT nodes_vdist;
   int* RESTRICT nodes_vbase;
   int* RESTRICT edgeorg;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   const int nterms = g->terms;
   const SCIP_Longint terms2 = (nterms - 1) * nterms;
   int maxnedges;
   int nnodes_distgraph;

   assert(nterms >= 1);

   SCIP_CALL( SCIPallocBufferArray(scip, distedge2org, nedges / 2) );
   edgeorg = *distedge2org;
    for( int e = 0; e < nedges / 2; e++ )
       edgeorg[e] = UNKNOWN;

   SCIP_CALL( graph_vnoiInit(scip, g, TRUE, vnoi) );
   graph_vnoiCompute(scip, g, *vnoi);
   nodes_vbase = (*vnoi)->nodes_base;
   nodes_vdist = (*vnoi)->nodes_dist;

   if( nedges >= terms2 )
   {
      assert(terms2 <= INT_MAX);
      maxnedges = terms2;
   }
   else
   {
      maxnedges = nedges;
   }

   SCIP_CALL( graph_init(scip, &(g_sd->distgraph), nterms, maxnedges, 1) );
   distgraph = g_sd->distgraph;
   nnodes_distgraph = 0;

   /* add the nodes */
   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) )
      {
         graph_knot_add(distgraph, STP_TERM_NONE);
         distnodes_id[k] = nnodes_distgraph++;
      }
      else
      {
         distnodes_id[k] = UNKNOWN;
      }
   }

   assert(distgraph->knots == nnodes_distgraph);
   assert(distgraph->knots == nterms);

   graph_knot_chg(distgraph, 0, STP_TERM);
   distgraph->source = 0;

   /* add the edges */
   for( int tail = 0; tail < nnodes; tail++ )
   {
      const int vbase_tail = nodes_vbase[tail];

      for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         const int vbase_head = nodes_vbase[head];

         assert(tail == g->tail[e]);

         if( vbase_tail != vbase_head )
         {
            const SCIP_Real distance = g->cost[e] + nodes_vdist[tail] + nodes_vdist[head];
            int ne;

            assert(Is_term(g->term[vbase_tail]));
            assert(Is_term(g->term[vbase_head]));
            assert(distnodes_id[vbase_tail] >= 0);
            assert(distnodes_id[vbase_head] >= 0);

            /* find the corresponding edge in the distance graph */
            for( ne = distgraph->outbeg[distnodes_id[vbase_tail]]; ne != EAT_LAST; ne = distgraph->oeat[ne] )
            {
               if( distgraph->head[ne] == distnodes_id[vbase_head] )
                  break;
            }

            /* edge exists already? */
            if( ne != EAT_LAST )
            {
               assert(ne >= 0);
               assert(distgraph->head[ne] == distnodes_id[vbase_head]);
               assert(distgraph->tail[ne] == distnodes_id[vbase_tail]);

               if( distgraph->cost[ne] > distance )
               {
                  distgraph->cost[ne]            = distance;
                  distgraph->cost[Edge_anti(ne)] = distance;
                  edgeorg[ne / 2] = e;
                  assert(ne <= maxnedges);
               }
            }
            else
            {
               edgeorg[distgraph->edges / 2] = e;
               graph_edge_add(scip, distgraph, distnodes_id[vbase_tail], distnodes_id[vbase_head], distance, distance);
               assert(distgraph->edges <= maxnedges);
            }
         }
      }
   }

   assert(graph_valid(scip, distgraph));

   return SCIP_OKAY;
}


/** builds MST costs (ordered) for distance graph */
static
void sdgraphBuildMstcosts(
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   const GRAPH* const distgraph = g_sd->distgraph;
   const PATH* const mst = g_sd->sdmst;
   SCIP_Real* RESTRICT mstcosts = g_sd->mstcosts;
   const int nnodes_distgraph = graph_get_nNodes(distgraph);

   assert(mst[0].edge == UNKNOWN);
   assert(nnodes_distgraph >= 1);

   for( int k = 1; k < nnodes_distgraph; k++ )
   {
      const int e = mst[k].edge;
      assert(e >= 0);
      assert(GE(distgraph->cost[e], 0.0));

      mstcosts[k - 1] = distgraph->cost[e];
   }

   SCIPsortDownReal(mstcosts, nnodes_distgraph - 1);

   /* debug sentinel */
   if( nnodes_distgraph > 1 )
   {
      mstcosts[nnodes_distgraph - 1] = -FARAWAY;
   }
   else
   {
      mstcosts[0] = 0.0;
   }
}


/** builds MST on distance graph */
static
SCIP_RETCODE sdgraphBuildMst(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   const VNOI*           vnoi,               /**< Voronoi */
   const int*            distedge2org,       /**< array of size nedges / 2 */
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   PATH* RESTRICT mst = g_sd->sdmst;
   GRAPH* distgraph = g_sd->distgraph;
   STP_Bool* RESTRICT orgedges_isInMst = g_sd->halfedge_isInMst;
   const int* const nodes_vbase = vnoi->nodes_base;
   const int* const nodes_vpred = vnoi->nodes_pred;
   SCIP_Real maxcost;
   const int nedges = graph_get_nEdges(g);
   const int nnodes_distgraph = graph_get_nNodes(distgraph);

   for( int k = 0; k < nnodes_distgraph; k++ )
      distgraph->mark[k] = TRUE;

   for( int e = 0; e < nedges / 2; e++ )
      orgedges_isInMst[e] = FALSE;

   SCIP_CALL( graph_path_init(scip, distgraph) );
   graph_path_exec(scip, distgraph, MST_MODE, distgraph->source, distgraph->cost, mst);
   graph_path_exit(scip, distgraph);

   assert(mst[0].edge == -1);

   maxcost = 0.0;
   for( int k = 1; k < nnodes_distgraph; k++ )
   {
      const int e = mst[k].edge;
      const int ne = distedge2org[e / 2];
      const SCIP_Real cost = distgraph->cost[e];

      assert(e >= 0);

      if( cost > maxcost )
         maxcost = cost;

      orgedges_isInMst[ne / 2] = TRUE;

      for( int v = g->head[ne]; v != nodes_vbase[v]; v = g->tail[nodes_vpred[v]] )
      {
         orgedges_isInMst[nodes_vpred[v] / 2] = TRUE;
      }

      for( int v = g->tail[ne]; v != nodes_vbase[v]; v = g->tail[nodes_vpred[v]] )
      {
         orgedges_isInMst[nodes_vpred[v]/ 2] = TRUE;
      }

      assert(e != EAT_LAST);
   }

   g_sd->mstmaxcost = maxcost;

   return SCIP_OKAY;
}

/** finalizes distance graph */
static
void sdgraphFinalize(
   SCIP*                 scip,               /**< SCIP data structure */
   VNOI**                vnoi,               /**> Voronoi */
   int**                 edgeorg             /**< array of size nedges / 2 */
)
{
   graph_vnoiFree(scip, vnoi);
   SCIPfreeBufferArray(scip, edgeorg);
}


/** apply pseudo eliminations provided */
SCIP_RETCODE reduce_applyPseudoDeletions(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const SCIP_Bool*      pseudoDelNodes,     /**< node with pseudo deletable nodes */
   GRAPH*                graph,              /**< graph data structure */
   SCIP_Real*            offsetp,            /**< offset pointer (for PC) */
   int*                  nelims              /**< number of eliminations */
)
{
   int adjvert[STP_DELPSEUDO_MAXGRAD];
   SCIP_Real cutoffs[STP_DELPSEUDO_MAXNEDGES];
   SCIP_Real cutoffsrev[STP_DELPSEUDO_MAXNEDGES];
   const PATH* nodeTo3TermsPaths = redcostdata->nodeTo3TermsPaths;
   const SCIP_Real* rootToNodeDist = redcostdata->rootToNodeDist;
   const SCIP_Real* redcost = redcostdata->redEdgeCost;
   const SCIP_Real cutoffbound = redcostdata->cutoff;
   int* nodetouchcount;
   const int nnodes = graph_get_nNodes(graph);
   SCIP_Bool success;
   const SCIP_Bool isPc = graph_pc_isPc(graph);
   const SCIP_Bool isExtendedOrg = graph->extended;

   assert(GE(cutoffbound, 0.0));
   assert(nodeTo3TermsPaths && rootToNodeDist && redcost);

   SCIP_CALL( SCIPallocBufferArray(scip, &nodetouchcount, nnodes) );

   for( int k = 0; k < nnodes; k++ )
      nodetouchcount[k] = 0;

   if( isPc )
      graph_pc_2orgcheck(scip, graph);

   *nelims = 0;

   for( int degree = 3; degree <= STP_DELPSEUDO_MAXGRAD; degree++ )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         SCIP_Real prize = -1.0;
         int edgecount = 0;
         SCIP_Bool rpc3term = FALSE;

         if( !pseudoDelNodes[k] || nodetouchcount[k] > 0 )
            continue;

         if( isPc && degree == 3 && graph_pc_knotIsNonLeafTerm(graph, k) && graph->grad[k] == 3 )
         {
            rpc3term = TRUE;
            prize = graph->prize[k];
         }
         else if( (degree != graph->grad[k] || Is_anyTerm(graph->term[k])) )
         {
            continue;
         }

         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
            nodetouchcount[graph->head[e]]++;

         if( rpc3term )
         {
            for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
            {
               const int head = graph->head[e];
               if( !graph_pc_knotIsDummyTerm(graph, head) )
                  adjvert[edgecount++] = head;
            }
         }
         else
         {
            for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
               adjvert[edgecount++] = graph->head[e];
         }

         assert(edgecount == degree);

         edgecount = 0;
         for( int i = 0; i < degree - 1; i++ )
         {
            const int vert = adjvert[i];
            for( int i2 = i + 1; i2 < degree; i2++ )
            {
               const int vert2 = adjvert[i2];

               assert(edgecount < STP_DELPSEUDO_MAXNEDGES);

               cutoffs[edgecount] = cutoffbound - (rootToNodeDist[vert] + nodeTo3TermsPaths[vert2].dist);
               cutoffsrev[edgecount] = cutoffbound - (rootToNodeDist[vert2] + nodeTo3TermsPaths[vert].dist);

               edgecount++;
            }
         }

         assert(edgecount > 0);

#ifdef SCIP_DEBUG
         SCIPdebugMessage("try pseudo-deletion of ");
         graph_knot_printInfo(graph, k);
#endif

         /* now try to eliminate */
         SCIP_CALL( graph_knot_delPseudo(scip, graph, redcost, cutoffs, cutoffsrev, k, &success) );

         if( success )
         {
            (*nelims)++;
            graph->mark[k] = FALSE;

            SCIPdebugMessage("deletion successful! \n");

            if( rpc3term )
            {
               assert(isPc);
               assert(offsetp);
               assert(GE(prize, 0.0));

               *offsetp += prize;
            }
         }
         else
         {
            for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
            {
               nodetouchcount[graph->head[e]]--;
               assert(nodetouchcount[graph->head[e]] >= 0);
            }
         }
      }
   }

   assert(graph_valid(scip, graph));

   if( isPc && isExtendedOrg != graph->extended )
      graph_pc_2trans(scip, graph);


   SCIPfreeBufferArray(scip, &nodetouchcount);

   return SCIP_OKAY;
}

/** initializes dynamic MST structure */
SCIP_RETCODE reduce_dcmstInit(
   SCIP*                 scip,               /**< SCIP */
   int                   maxnnodes,          /**< maximum number of nodes that can be handled */
   DCMST**               dcmst               /**< to be initialized */
)
{
   DCMST* mst;

   assert(scip && dcmst);
   assert(maxnnodes >= 1);

   SCIP_CALL( SCIPallocMemory(scip, dcmst) );

   mst = *dcmst;

   mst->maxnnodes = maxnnodes;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mst->edgestore), maxnnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mst->adjcost_buffer), maxnnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mst->nodemark), maxnnodes) );


   return SCIP_OKAY;
}


/** adds node to CSR "mst_in" and saves result in "mst_out" */
void reduce_dcmstAddNode(
   SCIP*                 scip,               /**< SCIP */
   const CSR*            mst_in,             /**< source */
   const SCIP_Real*      adjcosts,           /**< (undirected) adjacency costs for new node */
   DCMST*                dmst,               /**< underlying structure */
   CSR*                  mst_out             /**< target */
)
{
   assert(mst_in && adjcosts && dmst && mst_out);

   assert(reduce_dcmstMstIsValid(scip, mst_in));

   assert(mst_out->nnodes == mst_in->nnodes + 1);
   assert(mst_out->nedges_max == mst_in->nedges_max + 2);
   assert(mst_in->nnodes < dmst->maxnnodes);

   dcmstAddNode(mst_in, adjcosts, dmst);

   dcmstGetCSRfromStore(dmst, mst_out);

   assert(mst_out->nnodes == mst_in->nnodes + 1);
   assert(reduce_dcmstMstIsValid(scip, mst_out));
}


/** Adds node to CSR "mst".
 *  NOTE: There needs to be enough space in CSR arrays for one more node! */
void reduce_dcmstAddNodeInplace(
   SCIP*                 scip,               /**< SCIP */
   const SCIP_Real*      adjcosts,           /**< (undirected) adjacency costs for new node */
   DCMST*                dmst,               /**< underlying structure */
   CSR*                  mst                 /**< source/target */
)
{
   assert(mst && adjcosts && dmst);

   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(mst->nnodes < dmst->maxnnodes);

   dcmstAddNode(mst, adjcosts, dmst);

   mst->nnodes += 1;
   mst->nedges_max += 2;

   dcmstGetCSRfromStore(dmst, mst);

   assert(reduce_dcmstMstIsValid(scip, mst));
}

/** computes MST on 0 node */
void reduce_dcmstGet0NodeMst(
   SCIP*                 scip,               /**< SCIP */
   CSR*                  mst                 /**< MST */
)
{
   int* const start = mst->start;

   assert(mst->nnodes == 0);
   assert(mst->nedges_max == 0);

   start[0] = 0;

   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(EQ(0.0, reduce_dcmstGetWeight(scip, mst)));
}


/** computes MST on 1 node */
void reduce_dcmstGet1NodeMst(
   SCIP*                 scip,               /**< SCIP */
   CSR*                  mst                 /**< MST */
)
{
   int* const start = mst->start;

   assert(mst->nnodes == 1);
   assert(mst->nedges_max == 0);

   start[0] = 0;
   start[1] = 0;

   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(EQ(0.0, reduce_dcmstGetWeight(scip, mst)));
}


/** computes MST on 2 nodes */
void reduce_dcmstGet2NodeMst(
   SCIP*                 scip,               /**< SCIP */
   SCIP_Real             edgecost,           /**< edge cost */
   CSR*                  mst                 /**< MST */
)
{
   SCIP_Real* const cost = mst->cost;
   int* const start = mst->start;
   int* const head = mst->head;

   assert(edgecost > 0.0);
   assert(mst->nnodes == 2);
   assert(mst->nedges_max == 2);

   start[0] = 0;
   start[1] = 1;
   start[2] = 2;

   head[0] = 1;
   head[1] = 0;

   cost[0] = edgecost;
   cost[1] = edgecost;

   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(EQ(edgecost, reduce_dcmstGetWeight(scip, mst)));
}


/** computes MST on 3 nodes */
void reduce_dcmstGet3NodeMst(
   SCIP*                 scip,               /**< SCIP */
   SCIP_Real             edgecost01,         /**< edge cost */
   SCIP_Real             edgecost02,         /**< edge cost */
   SCIP_Real             edgecost12,         /**< edge cost */
   CSR*                  mst                 /**< MST */
)
{
   assert(0 && "implement me");
}


/** gets weight of MST extended along given vertex */
SCIP_Real reduce_dcmstGetExtWeight(
   SCIP*                 scip,               /**< SCIP */
   const CSR*            mst,                /**< MST for which to compute extended weight */
   const SCIP_Real*      adjcosts,           /**< (undirected) adjacency costs for new node */
   DCMST*                dmst                /**< underlying structure */
)
{
   SCIP_Real weight;
#ifndef NDEBUG
   /* since we have a tree, |E_{ext}| = |E| + 1 = |V| */
   const int nedges_ext = mst->nnodes;
#endif

   assert(scip && adjcosts && dmst);
   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(mst->nedges_max % 2 == 0);
   assert((mst->nedges_max / 2) + 1 == nedges_ext);

   dcmstAddNode(mst, adjcosts, dmst);

   weight = dcmstGetWeightFromStore(scip, mst->nnodes, dmst);

   if( GT(weight, FARAWAY) )
      weight = FARAWAY;

   assert(GE(weight, 0.0));

   return weight;
}


/** gets weight of MST */
SCIP_Real reduce_dcmstGetWeight(
   SCIP*                 scip,               /**< SCIP */
   const CSR*            mst_in              /**< source */
)
{
   SCIP_Real weight = 0.0;
   const int nedges = mst_in->nedges_max;
   const SCIP_Real* cost = mst_in->cost;

   assert(scip);
   assert(reduce_dcmstMstIsValid(scip, mst_in));

   for( int i = 0; i < nedges; i++ )
   {
      assert(cost[i] >= 0.0);

      weight += cost[i];
   }

   weight /= 2.0;

   if( GT(weight, FARAWAY) )
      weight = FARAWAY;

   assert(GE(weight, 0.0));

   return weight;
}


/** returns maximum number of nodes */
int reduce_dcmstGetMaxnnodes(
   const DCMST*          dmst                /**< underlying structure */
)
{
   assert(dmst);

   return dmst->maxnnodes;
}


/** Returns buffer of size 'reduce_dcmstGetMaxnnodes'.
  * NOTE: buffer is never used within any other function, apart from allocation and freeing.
  * NOTE: in debug mode the array is initialized to -1.0 */
SCIP_Real* reduce_dcmstGetAdjcostBuffer(
   const DCMST*          dmst                /**< underlying structure */
)
{
   assert(dmst);
   assert(dmst->adjcost_buffer);

#ifndef NDEBUG
   for( int i = 0; i < dmst->maxnnodes; i++ )
      dmst->adjcost_buffer[i] = -1.0;
#endif

   return dmst->adjcost_buffer;
}


/** frees dynamic MST structure */
void reduce_dcmstFree(
   SCIP*                 scip,               /**< SCIP */
   DCMST**               dcmst               /**< to be initialized */
)
{
   assert(scip && dcmst);

   SCIPfreeMemoryArray(scip, &((*dcmst)->nodemark));
   SCIPfreeMemoryArray(scip, &((*dcmst)->adjcost_buffer));
   SCIPfreeMemoryArray(scip, &((*dcmst)->edgestore));

   SCIPfreeMemory(scip, dcmst);
}


/** is the CSR a valid MST on any underlying graph (with number of nodes and edges of the CSR)? */
SCIP_Bool reduce_dcmstMstIsValid(
   SCIP*                 scip,               /**< SCIP */
   const CSR*            cmst                /**< the MST candidate */
)
{
   SCIP_Bool* visited;
   const int* const head_csr = cmst->head;
   const int nnodes = cmst->nnodes;
   SCIP_Bool isValid = TRUE;

#ifndef NDEBUG
   const int* const start_csr = cmst->start;
#endif

   if( nnodes == 0 )
   {
      assert(cmst->nedges_max == 0);
      assert(start_csr[0] == 0);

      return TRUE;
   }

   assert(nnodes >= 1);
   assert(cmst->nedges_max % 2 == 0);
   assert(start_csr[0] == 0);

   if( !graph_csr_isValid(cmst, FALSE) )
   {
      SCIPdebugMessage("CSR is broken! \n");
      return FALSE;
   }

   if( cmst->nnodes != (cmst->nedges_max / 2) + 1 )
   {
      SCIPdebugMessage("wrong nodes/edges ratio \n");
      return FALSE;
   }

   if( nnodes == 1 )
   {
      return TRUE;
   }

   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &visited, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      visited[i] = FALSE;

   for( int i = 0; i < cmst->nedges_max; i++ )
   {
      const int head = head_csr[i];

      assert(head >= 0 && head < nnodes);

      visited[head] = TRUE;
   }

   for( int i = 0; i < nnodes; i++ )
   {
      if( !visited[i] )
      {
         SCIPdebugMessage("mst does not contain node %d \n", i);

         isValid = FALSE;
         break;
      }
   }

   SCIPfreeMemoryArray(scip, &visited);

   return isValid;
}


/** initializes STAR structure */
SCIP_RETCODE reduce_starInit(
   SCIP*                 scip,               /**< SCIP */
   int                   maxdegree,          /**< maximum node degree that can be handled */
   STAR**                star                /**< the star */
)
{
   STAR* s;

   assert(scip && star);
   assert(maxdegree >= 3);

   SCIP_CALL( SCIPallocMemory(scip, star) );

   s = *star;

   s->starcenter = -1;
   s->nodeDegree = -1;
   s->starDegree = -1;
   s->maxNodeDegree = maxdegree;
   s->allStarsChecked = FALSE;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(s->edgeId), maxdegree) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(s->edgesSelected), maxdegree) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(s->edgesSelectedPos), maxdegree) );

   return SCIP_OKAY;
}


/** frees STAR structure */
void reduce_starFree(
   SCIP*                 scip,               /**< SCIP */
   STAR**                star                /**< the star */
)
{
   STAR* s;
   assert(scip && star);

   s = *star;
   assert(s);

   SCIPfreeMemoryArray(scip, &(s->edgesSelectedPos));
   SCIPfreeMemoryArray(scip, &(s->edgesSelected));
   SCIPfreeMemoryArray(scip, &(s->edgeId));

   SCIPfreeMemory(scip, star);
}


/** resets star data structure with new node data */
void reduce_starReset(
   const GRAPH*          g,                  /**< graph */
   int                   node,               /**< the node (degree <= STP_DELPSEUDO_MAXGRAD) */
   STAR*                 star                /**< the star */
)
{
   assert(g && star);
   assert(0 <= node && node < g->knots);
   assert(g->grad[node] <= star->maxNodeDegree);

   star->nodeDegree = g->grad[node];
   star->starDegree = star->nodeDegree;
   star->allStarsChecked = FALSE;
   star->starcenter = node;

   for( int e = g->outbeg[node], i = 0; e != EAT_LAST; e = g->oeat[e], i++ )
   {
      assert(i < star->nodeDegree);
      star->edgeId[i] = e;
   }

   /* initially, select the entire star */
   starSelectedPositionsReset(star);
}


/** gets center */
int reduce_starGetCenter(
   const STAR*           star                /**< the star (in/out) */
)
{
   assert(star);

   return star->starcenter;
}

/** gets next star */
const int* reduce_starGetNext(
   STAR*                 star,               /**< the star (in/out) */
   int*                  nedges              /**< number of edges of next star (out) */
)
{
   assert(star && nedges);
   assert(!reduce_starAllAreChecked(star));

   *nedges = star->starDegree;
   starSelectedEdgesUpdate(star);
   starSelectedPositionsSetNext(star);

   /* just finished? */
   if( starIsDeg2(star) )
      star->allStarsChecked = TRUE;

   assert(3 <= *nedges && *nedges <= star->maxNodeDegree);

   return star->edgesSelected;
}


/** gets ruled out edges after termination */
const int* reduce_starGetRuledOutEdges(
   STAR*                 star,               /**< the star */
   int*                  nedges              /**< number of edges of next star (out) */
)
{
   assert(star);
   assert(reduce_starAllAreChecked(star));

   // todo fill later
   // todo also add a method to abort early if no ruled-out edges can be found anymore

   return NULL;
}


/** sets current star to ruled-out */
void reduce_starSetRuledOut(
   STAR*                 star                /**< the star */
)
{
   assert(star);
   assert(!reduce_starAllAreChecked(star));


   // todo fill later for edge rule out
}


/** sets current star to failed */
void reduce_starSetFailed(
   STAR*                 star                /**< the star */
)
{
   assert(star);
   assert(!reduce_starAllAreChecked(star));

   // todo fill later for edge rule out
}


/** have all stars been checked? */
SCIP_Bool reduce_starAllAreChecked(
   const STAR*           star                /**< the star */
)
{
   assert(star);

   return star->allStarsChecked;
}


/** builds reduced costs data structure and returns it.
 * NOTE: memory needs to be freed again! */
SCIP_RETCODE reduce_redcostdataInit(
   SCIP*                 scip,               /**< SCIP */
   int                   nnodes,             /**< number of nodes */
   int                   nedges,             /**< number of edges */
   SCIP_Real             cutoff,             /**< reduced cost cutoff value or -1.0 if not used */
   int                   redCostRoot,        /**< graph root for reduced cost calculation */
   REDCOST*              redcostdata         /**< data to initialize */
)
{
   SCIP_Real* redEdgeCost;
   SCIP_Real* rootToNodeDist;
   PATH* nodeTo3TermsPaths;
   int* nodeTo3TermsBases;

   assert(scip);
   assert(nnodes >= 0);
   assert(nedges >= 0);
   assert(nedges % 2 == 0);
   assert(redCostRoot >= 0);
   assert(GE(cutoff, 0.0) || EQ(cutoff, 0.0));
   assert(redcostdata);

   SCIP_CALL( SCIPallocMemoryArray(scip, &redEdgeCost, nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &rootToNodeDist, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodeTo3TermsPaths, 3 * nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodeTo3TermsBases, 3 * nnodes) );

   redcostdata->redEdgeCost = redEdgeCost;
   redcostdata->rootToNodeDist = rootToNodeDist;
   redcostdata->nodeTo3TermsPaths = nodeTo3TermsPaths;
   redcostdata->nodeTo3TermsBases = nodeTo3TermsBases;
   redcostdata->cutoff = cutoff;
   redcostdata->redCostRoot = redCostRoot;

#ifndef NDEBUG
   redcostdata->nnodes = nnodes;
   redcostdata->nedges = nedges;
#endif

   return SCIP_OKAY;
}


/** frees member arrays */
void reduce_redcostdataFreeMembers(
   SCIP*                 scip,               /**< SCIP */
   REDCOST*              redcostdata         /**< data */
)
{
   assert(scip && redcostdata);

   SCIPfreeMemoryArray(scip, &(redcostdata->nodeTo3TermsBases));
   SCIPfreeMemoryArray(scip, &(redcostdata->nodeTo3TermsPaths));
   SCIPfreeMemoryArray(scip, &(redcostdata->rootToNodeDist));
   SCIPfreeMemoryArray(scip, &(redcostdata->redEdgeCost));
}


/** initializes SD graph */
SCIP_RETCODE reduce_sdgraphInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   VNOI* vnoi;
   int* edgeorg;
   assert(scip && g && sdgraph);

   SCIP_CALL( sdgraphAlloc(scip, g, sdgraph) );

   // todo already build biased!
   SCIP_CALL( sdgraphBuildDistgraph(scip, g, *sdgraph, &vnoi, &edgeorg) );
   SCIP_CALL( sdgraphBuildMst(scip, g, vnoi, edgeorg, *sdgraph) );

   sdgraphFinalize(scip, &vnoi, &edgeorg);

   return SCIP_OKAY;
}


/** return maximum MST edge cost */
SCIP_Real reduce_sdgraphGetMaxCost(
   const SDGRAPH*        sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(GE(sdgraph->mstmaxcost, 0.0));

   return sdgraph->mstmaxcost;
}


/** initializes SD graph */
const STP_Bool* reduce_sdgraphGetMstHalfMark(
   const SDGRAPH*        sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(sdgraph->halfedge_isInMst);

   return sdgraph->halfedge_isInMst;
}


/** MST costs in descending order available? */
SCIP_Bool reduce_sdgraphHasOrderedMstCosts(
   const SDGRAPH*         sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(sdgraph->mstcosts);

   return sdgraph->mstcostsReady;
}


/** initializes all MST costs in descending order */
void reduce_sdgraphInitOrderedMstCosts(
   SDGRAPH*              sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(sdgraph->mstcosts);

   if( !sdgraph->mstcostsReady )
   {
      sdgraphBuildMstcosts(sdgraph);
      sdgraph->mstcostsReady = TRUE;
   }

   assert(reduce_sdgraphHasOrderedMstCosts(sdgraph));
   assert(sdgraph->mstcosts[0] == sdgraph->mstmaxcost);
}


/** returns all MST costs in descending order */
const SCIP_Real* reduce_sdgraphGetOrderedMstCosts(
   const SDGRAPH*         sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(reduce_sdgraphHasOrderedMstCosts(sdgraph));
   assert(sdgraph->mstcosts[0] == sdgraph->mstmaxcost);

   return sdgraph->mstcosts;
}


/** frees SD graph */
void reduce_sdgraphFree(
   SCIP*                 scip,               /**< SCIP */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   SDGRAPH* g_sd;
   assert(scip && sdgraph);

   g_sd = *sdgraph;
   assert(g_sd);

   SCIPfreeMemoryArray(scip, &(g_sd->sdmst));
   SCIPfreeMemoryArray(scip, &(g_sd->mstcosts));
   SCIPfreeMemoryArray(scip, &(g_sd->nodes_id));
   SCIPfreeMemoryArray(scip, &(g_sd->halfedge_isInMst));

   graph_free(scip, &(g_sd->distgraph), TRUE);

   SCIPfreeMemory(scip, sdgraph);
}


/** initializes TPATHS structure */
SCIP_RETCODE reduce_tpathsInit(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph NOTE: will make the graph, thus not const :(
                                                  terrible design */
   TPATHS**              tpaths              /**< the terminal paths */
)
{

   assert(scip);
   assert(graph_isMarked(g));
   assert(STP_TPATHS_NTERMBASES == 4);

   SCIP_CALL( tpathsAlloc(scip, g, tpaths) );
   SCIP_CALL( tpathsBuild(scip, g, *tpaths) );

   return SCIP_OKAY;
}


/** frees TPATHS structure */
void reduce_tpathsFree(
   SCIP*                 scip,               /**< SCIP */
   TPATHS**              tpaths              /**< the terminal paths */
)
{
   TPATHS* tp;
   assert(scip && tpaths);

   tp = *tpaths;
   assert(tp);

   SCIPfreeMemoryArray(scip, &(tp->state));
   SCIPfreeMemoryArray(scip, &(tp->termbases));
   SCIPfreeMemoryArray(scip, &(tp->termpaths));

   SCIPfreeMemory(scip, tpaths);
}
