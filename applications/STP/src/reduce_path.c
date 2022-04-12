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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reduce_path.c
 * @brief  Path deletion reduction test for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements an improved version of the so called "path substitution test" by Polzin and Vahdati Daneshmand.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG

#include <assert.h>
#include "reduce.h"
#include "extreduce.h"

#define SP_MAXNDISTS 10000
#define SP_MAXLENGTH 11
#define SP_MAXDEGREE 8
#define SP_MAXDEGREE_FAST 5
#define SP_MAXDEGREE_PC 10
#define SP_MAXNSTARTS 10000
#define VNODES_UNSET -1
#define SP_MAXNPULLS 4


/** path replacement */
typedef struct path_replacement
{
   EXTPERMA*             extperma;           /**< extension data or NULL */
   SDPROFIT*             sdprofit;           /**< SD profit */
   STP_Vectype(SCIP_Real) firstneighborcosts; /**< edge costs from tail of path to non-path neighbors */
   STP_Vectype(int)      firstneighbors;     /**< first neighbors (only used in case SDs are being used) */
   STP_Vectype(int)      currneighbors;      /**< current neighbors */
   STP_Vectype(int)      pathedges;          /**< edges of path */
   STP_Vectype(int)      visitednodes;       /**< visited nodes */
   SCIP_Real* RESTRICT   sp_dists;           /**< distances to neighbors of path start node */
   int* RESTRICT         sp_starts;          /**< CSR like starts for each node in the path, pointing to sp_dists */
   int* RESTRICT         nodes_index;        /**< maps each node to index in 0,1,2,..., or to VNODES_UNSET, VNODES_INPATH  */
   SCIP_Bool* RESTRICT   nodes_isTerm;       /**< terminal node? */
   SCIP_Bool* RESTRICT   nodeindices_isPath; /**< is a node index in the path? */
   SCIP_Bool* RESTRICT   firstneighbors_isKilled; /**< for each first neighbor index: is ruled-out? */
   SCIP_Real             pathcost;           /**< cost of path */
   int                   nfirstneighbors;    /**< number of neighbors of first path node */
   int                   pathtail;           /**< first node of path */
   int                   failneighbor;       /**< temporary */
   int                   nnodes;             /**< number of nodes */
   int                   maxdegree;
   SCIP_Bool             probIsPc;           /**< prize-collecting problem? */
   SCIP_Bool             useSd;              /**< use special distances? */
} PR;


/*
 * Local methods
 */

/** gets profit for given node */
inline static
SCIP_Real getSdProfit(
   const SDPROFIT*      sdprofit,           /**< the SD profit */
   const int*           nodes_index,        /**< index array */
   int                  node,               /**< node to get profit for */
   int                  nonsource           /**< node that should not be a source */
)
{
   const int source1 = sdprofit->nodes_biassource[node];

   assert(GE(sdprofit->nodes_bias[node], 0.0));
   assert(LE(sdprofit->nodes_bias[node], FARAWAY));
   assert(GE(sdprofit->nodes_bias2[node], 0.0));
   assert(LE(sdprofit->nodes_bias2[node], FARAWAY));
   assert(GE(sdprofit->nodes_bias[node], sdprofit->nodes_bias2[node]));

   if( source1 == node )
      return sdprofit->nodes_bias[node];

   if( source1 != nonsource && nodes_index[source1] < 0 )
   {
      return sdprofit->nodes_bias[node];
   }
   else
   {
      const int source2 = sdprofit->nodes_biassource2[node];

      if( source2 != nonsource && nodes_index[source2] < 0 )
      {
         return sdprofit->nodes_bias2[node];
      }
   }

   return 0.0;
}


/** gets head of path */
static inline
int pathGetHead(
   const GRAPH*          g,                  /**< graph data structure */
   const PR*             pr                  /**< path replacement data structure */
   )
{
   const int npathedges = StpVecGetSize(pr->pathedges);
   return g->head[pr->pathedges[npathedges - 1]];
}


/** deletes given edge */
static inline
void deleteEdge(
   SCIP*                 scip,               /**< SCIP */
   int                   edge,               /**< edge to delete */
   GRAPH*                g,                  /**< graph data structure */
   PR*                   pr                  /**< path replacement data structure */
   )
{
   if( pr->useSd )
   {
      EXTPERMA* const extperma = pr->extperma;
      assert(extperma && extperma->distdata_default);

      if( extperma->solIsValid )
      {
         const int* const result = extperma->result;
         if( result[edge] == CONNECT || result[flipedge(edge)] == CONNECT )
            extperma->solIsValid = FALSE;
      }

      extreduce_edgeRemove(scip, edge, g, extperma->distdata_default, extperma);
   }
   else
   {
      graph_edge_delFull(scip, g, edge, TRUE);
   }
}


/** removes nodes of degree one and unmarks */
static
void deletenodesDeg1(
   SCIP*                 scip,               /**< SCIP */
   EXTPERMA*             extperma,           /**< extension data */
   GRAPH*                g                   /**< graph data structure (in/out) */
)
{
   const int nnodes = graph_get_nNodes(g);
   SCIP_Bool rerun = TRUE;

   while( rerun )
   {
      rerun = FALSE;
      for( int i = 0; i < nnodes; ++i )
      {
         if( g->grad[i] == 1 && !Is_term(g->term[i]) )
         {
            const int sibling = g->head[g->outbeg[i]];

            extreduce_edgeRemove(scip, g->outbeg[i], g, extperma->distdata_default, extperma);
            assert(g->grad[i] == 0);

            g->mark[i] = FALSE;

            if( Is_term(g->term[sibling]) )
               continue;

            if( g->grad[sibling] == 0 )
               g->mark[sibling] = FALSE;
            else if( g->grad[sibling] == 1 )
               rerun = TRUE;
         }
      }
   }
}

/** tries to rule out neighbors of path head */
static inline
void ruleOutFromHead(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   PR*                   pr,                 /**< path replacement data structure */
   SCIP_Bool*            needFullRuleOut     /**< pointer to store whether full rule-out is needed */
   )
{
   const SCIP_Real* const sp_dists = pr->sp_dists;
   const SCIP_Real pathcost = pr->pathcost;
   const int ncurrneighbors = StpVecGetSize(pr->currneighbors);
   const int nfirst = pr->nfirstneighbors;
   int failneighbor = VNODES_UNSET;
   const SCIP_Bool useSd = pr->useSd;
   DISTDATA* const distdata = useSd ? pr->extperma->distdata_default : NULL;

   assert(FALSE == *needFullRuleOut);
   SCIPdebugMessage("starting head rule-out with pathcost=%f \n", pr->pathcost);

   for( int i = 0; i < ncurrneighbors; i++ )
   {
      int j;
      const int neighbor = pr->currneighbors[i];
      const int sp_start = pr->sp_starts[pr->nodes_index[neighbor]];
      int nfails = 0;

      /* NOTE: we need to be able to rule out all but one of the path-tail neighbors */

      assert(pr->nodes_index[neighbor] >= 0);

      for( j = 0; j < nfirst; j++ )
      {
         if( pr->firstneighbors_isKilled[j] )
            continue;

         SCIPdebugMessage("dist for neighbor=%d, idx=%d: %f \n", neighbor, j, sp_dists[sp_start + j]);
         if( GT(sp_dists[sp_start + j], pathcost) )
         {
            if( useSd )
            {
               const int fneighbor = pr->firstneighbors[j];
               SCIP_Real sd = extreduce_distDataGetSdDoubleEq(scip, g, pathcost, fneighbor, neighbor, distdata);

               if( LT(sd, pathcost) )
                  continue;

               if( EQ(sd, pathcost) )
               {
                  const int startedge = pr->pathedges[0];
                  sd = extreduce_distDataGetSdDoubleForbiddenSingle(scip, g, startedge, neighbor, fneighbor, distdata);
                  if( LE(sd, pathcost) )
                     continue;
               }
            }

            /* second fail? */
            if( nfails++ > 0 )
               break;
         }
      }

      if( nfails > 1 )
      {
         SCIPdebugMessage("neighbor %d not head-ruled-out \n", neighbor);

         /* second time that a current neighbor could not be ruled-out? */
         if( failneighbor != VNODES_UNSET )
         {
            pr->failneighbor = VNODES_UNSET;
            *needFullRuleOut = TRUE;
            return;
         }
         failneighbor = neighbor;
      }
   }

   pr->failneighbor = failneighbor;
}


/** tries to rule out neighbors of path tail */
static inline
void ruleOutFromTailSingle(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   PR*                   pr,                 /**< path replacement data structure */
   SCIP_Bool*            isRuledOut          /**< ruled-out? */
)
{
   const SCIP_Real* const sp_dists = pr->sp_dists;
   const int* const sp_starts = pr->sp_starts;
   const int pathhead = pathGetHead(g, pr);
   const int pathhead_idx = pr->nodes_index[pathhead];
   const SCIP_Bool useSd = pr->useSd;
   DISTDATA* const distdata = useSd ? pr->extperma->distdata_default : NULL;

   assert(!*isRuledOut);

   SCIPdebugMessage("try tail rule-out with pathcost=%f \n", pr->pathcost);

   for( int i = 0; i < pr->nfirstneighbors; i++ )
   {
      if( !pr->firstneighbors_isKilled[i] )
      {
         const SCIP_Real extpathcost = pr->pathcost + pr->firstneighborcosts[i];
         const SCIP_Real altpathcost = sp_dists[sp_starts[pathhead_idx] + i];

         SCIPdebugMessage("extpathcost for idx=%d: %f \n", i, extpathcost);
         SCIPdebugMessage("althpathost for idx=%d: %f \n", i, altpathcost);

         if( GT(altpathcost, extpathcost) )
         {
            if( useSd )
            {
               const int fneighbor = pr->firstneighbors[i];
               SCIP_Real sd = extreduce_distDataGetSdDoubleEq(scip, g, extpathcost, fneighbor, pathhead, distdata);

               if( LT(sd, extpathcost) )
                  continue;

               if( EQ(sd, extpathcost) )
               {
                  const int startedge = pr->pathedges[0];
                  sd = extreduce_distDataGetSdDoubleForbiddenSingle(scip, g, startedge, pathhead, fneighbor, distdata);
                  if( LE(sd, extpathcost) )
                     continue;
               }
            }

            SCIPdebugMessage("no rule-out for initial neighbor idx=%d \n", i);
            return;
         }
      }
   }

   *isRuledOut = TRUE;
}


/** tries to rule out neighbors of path tail */
static inline
void ruleOutFromTailCombs(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   PR*                   pr,                 /**< path replacement data structure */
   SCIP_Bool*            isRuledOut          /**< ruled-out? */
)
{
   const SCIP_Real* const sp_dists = pr->sp_dists;
   const int* const sp_starts = pr->sp_starts;
   const int pathhead = pathGetHead(g, pr);
   const int pathhead_idx = pr->nodes_index[pathhead];
   int nfails = 0;
   const SCIP_Bool useSd = pr->useSd;
   DISTDATA* const distdata = useSd ? pr->extperma->distdata_default : NULL;

   SCIPdebugMessage("try full tail rule-out with pathcost=%f \n", pr->pathcost);

   for( int i = 0; i < pr->nfirstneighbors; i++ )
   {
      const SCIP_Real altpathcost = sp_dists[sp_starts[pathhead_idx] + i];

      SCIPdebugMessage("althpathost for idx=%d: %f \n", i, altpathcost);

      if( GT(altpathcost, pr->pathcost) )
      {
         if( useSd )
         {
            const int fneighbor = pr->firstneighbors[i];
            SCIP_Real sd = extreduce_distDataGetSdDoubleEq(scip, g, pr->pathcost, fneighbor, pathhead, distdata);

            if( LT(sd, pr->pathcost) )
            {
               pr->firstneighbors_isKilled[i] = TRUE;
               continue;
            }

            if( EQ(sd, pr->pathcost) )
            {
               const int startedge = pr->pathedges[0];
               sd = extreduce_distDataGetSdDoubleForbiddenSingle(scip, g, startedge, pathhead, fneighbor, distdata);
               if( LE(sd, pr->pathcost) )
               {
                  pr->firstneighbors_isKilled[i] = TRUE;
                  continue;
               }
            }
         }

         SCIPdebugMessage("no full rule-out for initial neighbor idx=%d \n", i);
         nfails++;
      }
      else
      {
         pr->firstneighbors_isKilled[i] = TRUE;
      }
   }

   if( nfails > 1 )
   {
      SCIPdebugMessage("...no full rule-out \n");
      assert(!*isRuledOut);
   }
   else
   {
      *isRuledOut = TRUE;
   }
}


/** adds new path node */
static inline
void addPathNode(
   SCIP*                 scip,               /**< SCIP */
   int                   node,               /**< node to add */
   PR*                   pr                  /**< path replacement data structure */
   )
{
   assert(pr->nodes_index[node] == VNODES_UNSET);

   pr->nodes_index[node] = StpVecGetSize(pr->visitednodes);
   StpVecPushBack(scip, pr->visitednodes, node);
   pr->nodeindices_isPath[pr->nodes_index[node]] = TRUE;
}

/** adds new NON-path node */
static inline
void addNonPathNode(
   SCIP*                 scip,               /**< SCIP */
   int                   node,               /**< node to add */
   PR*                   pr                  /**< path replacement data structure */
   )
{
   assert(pr->nodes_index[node] == VNODES_UNSET);

   pr->nodes_index[node] = StpVecGetSize(pr->visitednodes);
   StpVecPushBack(scip, pr->visitednodes, node);
   pr->nodeindices_isPath[pr->nodes_index[node]] = FALSE;
}


/** adds new path head edge to failed node */
static inline
void addPathHeadEdge(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   PR*                   pr                  /**< path replacement data structure */
   )
{
   const DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const dcsr_range = dcsr->range;
   const int* const dcsr_heads = dcsr->head;
   const int pathhead = pathGetHead(g, pr);
   const int extnode = pr->failneighbor;
   int j;

   assert(pr->nodes_index[extnode] >= 0);
   assert(!pr->nodeindices_isPath[pr->nodes_index[extnode]]);
   assert(extnode >= 0);

   SCIPdebugMessage("adding new path edge to node %d \n", extnode);

   for( j = dcsr_range[pathhead].start; j != dcsr_range[pathhead].end; j++ )
   {
      if( dcsr_heads[j] == extnode )
         break;
   }
   assert(j != dcsr_range[pathhead].end);
   assert(EQ(dcsr->cost[j], g->cost[dcsr->edgeid[j]]));

   StpVecPushBack(scip, pr->pathedges, dcsr->edgeid[j]);
   pr->pathcost += dcsr->cost[j];
   pr->nodeindices_isPath[pr->nodes_index[extnode]] = TRUE;

   if( pr->probIsPc )
   {
      assert(LE(g->prize[pathhead], dcsr->cost[j]));
      pr->pathcost -= g->prize[pathhead];
   }

   assert(pathGetHead(g, pr) == extnode);
}


/** adds first two path nodes */
static inline
void addInitialPathNodes(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   int                   startedge_tail,     /**< tail of start edge */
   int                   startdege_head,     /**< head of start edge */
   PR*                   pr                  /**< path replacement data structure */
   )
{
   const DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const dcsr_range = dcsr->range;
   const int* const dcsr_heads = dcsr->head;
   SCIP_Real* RESTRICT sp_dists = pr->sp_dists;
   int* RESTRICT sp_starts = pr->sp_starts;
   int starts_final = pr->nfirstneighbors;

   /* add tail node */
   addPathNode(scip, startedge_tail, pr);
   sp_starts[starts_final + 1] = sp_starts[starts_final] + pr->nfirstneighbors;
   for( int i = sp_starts[starts_final], j = 0; i != sp_starts[starts_final + 1]; i++, j++ )
      sp_dists[i] = pr->firstneighborcosts[j];

   /* add head node*/
   addPathNode(scip, startdege_head, pr);
   starts_final++;
   sp_starts[starts_final + 1] = sp_starts[starts_final] + pr->nfirstneighbors;
   for( int i = sp_starts[starts_final]; i != sp_starts[starts_final + 1]; i++ )
      sp_dists[i] = FARAWAY;

   /* update distances by using single edges */
   for( int j = dcsr_range[startdege_head].start; j != dcsr_range[startdege_head].end; j++ )
   {
      const int head = dcsr_heads[j];
      const int head_index = pr->nodes_index[head];

      /* unvisited or path tail? */
      if( head_index < 0 || head == startedge_tail )
         continue;

      assert(0 <= head_index && head_index < pr->nfirstneighbors);
      assert(EQ(sp_dists[sp_starts[starts_final] + head_index], FARAWAY));
      SCIPdebugMessage("setting distance from first path node(%d) for initial neighbor (idx=%d) to %f \n",
            startdege_head, head_index, dcsr->cost[j]);
      sp_dists[sp_starts[starts_final] + head_index] = dcsr->cost[j];
   }
}


/** collects neighbors */
static inline
void pathneighborsCollect(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   PR*                   pr                  /**< path replacement data structure */
   )
{
   const DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const dcsr_range = dcsr->range;
   SCIP_Real* RESTRICT sp_dists = pr->sp_dists;
   const int* const dcsr_heads = dcsr->head;
   int* RESTRICT nodes_index = pr->nodes_index;
   int* RESTRICT sp_starts = pr->sp_starts;
   const int basenode = pathGetHead(g, pr);
   const int nfirstneighbors = pr->nfirstneighbors;

   SCIPdebugMessage("extending path to node %d \n", basenode);

   StpVecClear(pr->currneighbors);

   for( int i = dcsr_range[basenode].start; i != dcsr_range[basenode].end; i++ )
   {
      const int head = dcsr_heads[i];
      int head_index = nodes_index[head];

      if( head_index >= 0 && pr->nodeindices_isPath[head_index] )
         continue;

      assert(head_index < StpVecGetSize(pr->visitednodes));
      SCIPdebugMessage("adding neighbor %d  head_index=%d\n", head, head_index);

      StpVecPushBack(scip, pr->currneighbors, head);

      if( head_index == VNODES_UNSET )
      {
         head_index = StpVecGetSize(pr->visitednodes);
         SCIPdebugMessage("mapping new neighbor %d->%d \n", head, head_index);
         nodes_index[head] = head_index;
         StpVecPushBack(scip, pr->visitednodes, head);
         pr->nodeindices_isPath[head_index] = FALSE;
         assert(head_index + 1 < SP_MAXNSTARTS);

         sp_starts[head_index + 1] = sp_starts[head_index] + nfirstneighbors;

         for( int j = sp_starts[head_index]; j != sp_starts[head_index + 1]; j++ )
            sp_dists[j] = FARAWAY;
      }
   }
}


/** updates distances from neighbors */
static inline
void pathneighborsUpdateDistances(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   PR*                   pr                  /**< path replacement data structure */
   )
{
   STP_Vectype(int) currneighbors;
   const SDPROFIT* const sdprofit = pr->sdprofit;
   const DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const dcsr_range = dcsr->range;
   const SCIP_Real* const dcsr_costs = dcsr->cost;
   SCIP_Real* RESTRICT sp_dists = pr->sp_dists;
   const int* const dcsr_heads = dcsr->head;
   const int* const nodes_index = pr->nodes_index;
   const int* const sp_starts = pr->sp_starts;
   const int nfirstneighbors = pr->nfirstneighbors;
   const int nneighbors = StpVecGetSize(pr->currneighbors);
   const int pathtail = pr->pathtail;
   const int pathhead = pathGetHead(g, pr);
   const int nloops = MIN(SP_MAXNPULLS, nneighbors);

   /* NOTE: to keep the following loop easy */
   StpVecPushBack(scip, pr->currneighbors, pathhead);
   currneighbors = pr->currneighbors;

   for( int loop = 0; loop < nloops; loop++ )
   {
      SCIP_Bool hasUpdates = FALSE;
      /* compute distances to each of the initial neighbors */
      for( int iter = 0; iter <= nneighbors; iter++ )
      {
         const int node = currneighbors[iter];
         const int node_start = sp_starts[nodes_index[node]];

         assert(nodes_index[node] >= 0);

         for( int i = dcsr_range[node].start; i != dcsr_range[node].end; i++ )
         {
            const int head = dcsr_heads[i];

            if( nodes_index[head] != VNODES_UNSET )
            {
               SCIP_Real head_profit = 0.0;
               const int head_start = sp_starts[nodes_index[head]];

               if( head == pathtail && node == pathhead )
                  continue;

               if( nodes_index[head] > nfirstneighbors )
               {
                  assert(head != pathtail);
                  head_profit = getSdProfit(sdprofit, nodes_index, head, node);
               }

               for( int k = 0; k < nfirstneighbors; k++ )
               {
                  SCIP_Real newdist = dcsr_costs[i] + sp_dists[head_start + k];
                  SCIP_Real profitBias = MIN(head_profit, dcsr_costs[i]);

                  profitBias = MIN(profitBias, sp_dists[head_start + k]);
                  //printf("%f \n", profitBias);

                  newdist -= profitBias;

                  if( LT(newdist, sp_dists[node_start + k]) )
                  {
                     SCIPdebugMessage("(l=%d) updating distance for node %d to orgindex %d to %f \n",
                           loop, node, k, newdist);
                     sp_dists[node_start + k] = newdist;
                     hasUpdates = TRUE;
                  }
               }
            }
         }
      }

      if( !hasUpdates )
         break;
   }

   StpVecPopBack(currneighbors);
}


/** preprocess for doing path extension later on */
static inline
void pathExendPrepare(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   startedge,          /**< first edge */
   PR*                   pr                  /**< path replacement data structure */
   )
{
   const int basetail = g->tail[startedge];
   const int basehead = g->head[startedge];
   const DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const dcsr_range = dcsr->range;
   const int* const dcsr_heads = dcsr->head;
   const SCIP_Real* const dcsr_costs = dcsr->cost;
   SCIP_Real* RESTRICT sp_dists = pr->sp_dists;
   int* RESTRICT sp_starts = pr->sp_starts;
   const SCIP_Bool isPc = pr->probIsPc;

   assert(0 == StpVecGetSize(pr->pathedges));
   assert(0 == StpVecGetSize(pr->visitednodes));
   assert(0 == StpVecGetSize(pr->currneighbors));
   assert(0 == StpVecGetSize(pr->firstneighborcosts));
   assert(0 == StpVecGetSize(pr->firstneighbors));

   SCIPdebugMessage("---Checking edge %d->%d \n\n", basetail, basehead);

   pr->pathcost = g->cost[startedge];
   pr->pathtail = basetail;
   StpVecPushBack(scip, pr->pathedges, startedge);

   for( int i = dcsr_range[basetail].start; i != dcsr_range[basetail].end; i++ )
   {
      const int head = dcsr_heads[i];
      assert(pr->nodes_index[head] == VNODES_UNSET);
      assert(!isPc || !graph_pc_knotIsDummyTerm(g, head));

      if( head == basehead )
         continue;

      SCIPdebugMessage("mapping first neighbor %d->%d \n", head, StpVecGetSize(pr->visitednodes));

      addNonPathNode(scip, head, pr);
      StpVecPushBack(scip, pr->currneighbors, head);
      StpVecPushBack(scip, pr->firstneighborcosts, dcsr_costs[i]);

      if( pr->useSd )
         StpVecPushBack(scip, pr->firstneighbors, head);
   }

   pr->nfirstneighbors = StpVecGetSize(pr->currneighbors);
   SCIPdebugMessage("Having %d initial neighbors \n", pr->nfirstneighbors);

   sp_starts[0] = 0;
   for( int i = 0; i < pr->nfirstneighbors; i++ )
   {
      const int neighbor = pr->currneighbors[i];
      const SCIP_Real basecost = pr->firstneighborcosts[i];
      sp_starts[i + 1] = sp_starts[i] + pr->nfirstneighbors;

      /* set 2-edge distances via path tail */
      if( isPc )
      {
         for( int j = sp_starts[i], k = 0; j != sp_starts[i + 1]; j++, k++ )
         {
            assert(LE(g->prize[basetail], MIN(basecost, pr->firstneighborcosts[k])));
            sp_dists[j] = basecost + pr->firstneighborcosts[k] - g->prize[basetail];
         }

         assert(EQ(sp_dists[sp_starts[i] + i], 2.0 * basecost - g->prize[basetail]));
      }
      else
      {
         for( int j = sp_starts[i], k = 0; j != sp_starts[i + 1]; j++, k++ )
            sp_dists[j] = basecost + pr->firstneighborcosts[k];

         assert(EQ(sp_dists[sp_starts[i] + i], 2.0 * basecost));
      }

#ifdef SCIP_DEBUG
      for( int j = sp_starts[i]; j != sp_starts[i + 1]; j++ )
         printf("%d->%d dist=%f \n", i, j - sp_starts[i], sp_dists[j]);
#endif

      /* set self-distance */
      sp_dists[sp_starts[i] + i] = 0.0;

      /* update distances by using single edges */
      for( int j = dcsr_range[neighbor].start; j != dcsr_range[neighbor].end; j++ )
      {
         const int head = dcsr_heads[j];
         const int head_index = pr->nodes_index[head];

         /* unvisited or path tail? */
         if( head_index < 0 || head == basetail )
            continue;

         assert(0 <= head_index && head_index < pr->nfirstneighbors);

#ifdef SCIP_DEBUG
         if( LT(dcsr_costs[j], sp_dists[sp_starts[i] + head_index]) )
            printf("updating %d->%d distold=%f distnew=%f \n", i, head_index, sp_dists[sp_starts[i] + head_index], dcsr_costs[j]);
#endif
         if( LT(dcsr_costs[j], sp_dists[sp_starts[i] + head_index]) )
            sp_dists[sp_starts[i] + head_index] = dcsr_costs[j];
      }
   }

   for( int i = 0; i < pr->nfirstneighbors; i++ )
      pr->firstneighbors_isKilled[i] = FALSE;

   addInitialPathNodes(scip, g, basetail, basehead, pr);

   /* NOTE: firstneighborcosts are only used for a subpath, thus we can already subtract the prize here */
   if( isPc )
   {
      for( int i = 0; i < pr->nfirstneighbors; i++ )
      {
         pr->firstneighborcosts[i] -= g->prize[basetail];
         assert(GE(pr->firstneighborcosts[i], 0.0));
      }
   }
}


/** tries to eliminate edges of path starting from edge */
static inline
void pathExend(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   PR*                   pr,                 /**< path replacement data structure */
   SCIP_Bool*            isExendible,        /**< extendible? */
   SCIP_Bool*            isRedundant         /**< redundant? */
   )
{
   const int npathedges = StpVecGetSize(pr->pathedges);
   const int pathhead = pathGetHead(g, pr);
   SCIP_Bool needCombRuleOut = FALSE;
   SCIP_Bool isSingleRuledOut = FALSE;
   SCIP_Bool isCombRuledOut = FALSE;

   assert(*isExendible);
   assert(!(*isRedundant));

   if( npathedges > SP_MAXLENGTH || g->grad[pathhead] > pr->maxdegree || pr->nodes_isTerm[pathhead] )
   {
      *isExendible = FALSE;
      return;
   }

   pathneighborsCollect(scip, g, pr);
   pathneighborsUpdateDistances(scip, g, pr);

   if( StpVecGetSize(pr->currneighbors) == 0 )
   {
      // todo really true?
      *isRedundant = TRUE;
      return;
   }

   ruleOutFromHead(scip, g, pr, &needCombRuleOut);

   /* NOTE: if we have exactly one neighbor, and a failed neighbor, they are the same and we have to extend */
   if( StpVecGetSize(pr->currneighbors) == 1 && pr->failneighbor != VNODES_UNSET )
   {
      addPathHeadEdge(scip, g, pr);
      return;
   }

   /* we try combination rule-out anyway! */
   ruleOutFromTailCombs(scip, g, pr, &isCombRuledOut);

   if( needCombRuleOut && !isCombRuledOut )
   {
      *isExendible = FALSE;
      return;
   }

   ruleOutFromTailSingle(scip, g, pr, &isSingleRuledOut);

   if( !isSingleRuledOut )
   {
      *isExendible = FALSE;
      return;
   }

   if( isCombRuledOut )
   {
      *isRedundant = TRUE;
      return;
   }

   if( pr->failneighbor == VNODES_UNSET )
   {
      *isRedundant = TRUE;
      return;
   }

   /* if neither of the above cases holds, we extend along the failed neighbor */
   addPathHeadEdge(scip, g, pr);
}


/** initializes */
static
SCIP_RETCODE prInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data or NULL */
   PR**                  pathreplace         /**< to initialize */
   )
{
   PR* pr;
   const int nnodes = graph_get_nNodes(g);

   SCIP_CALL( SCIPallocMemory(scip, pathreplace) );
   pr = *pathreplace;

   pr->useSd = (extperma != NULL);
   pr->probIsPc = graph_pc_isPc(g);
   pr->maxdegree = pr->probIsPc ? SP_MAXDEGREE_PC : SP_MAXDEGREE;
   if( extperma && extperma->mode == extred_fast )
      pr->maxdegree = SP_MAXDEGREE_FAST;

   SCIP_CALL( SCIPallocMemoryArray(scip, &pr->nodes_isTerm, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &pr->nodeindices_isPath, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &pr->nodes_index, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &pr->sp_starts, SP_MAXNSTARTS) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &pr->sp_dists, SP_MAXNDISTS) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &pr->firstneighbors_isKilled, pr->maxdegree + 1) );

   SCIP_CALL( reduce_sdprofitInit(scip, g, &pr->sdprofit) );

   for( int i = 0; i < nnodes; i++ )
      pr->nodes_index[i] = VNODES_UNSET;

   graph_getIsTermArray(g, pr->nodes_isTerm);

   pr->extperma = extperma;
   pr->currneighbors = NULL;
   pr->firstneighborcosts = NULL;
   pr->pathedges = NULL;
   pr->firstneighbors = NULL;
   pr->visitednodes = NULL;
   pr->nfirstneighbors = -1;
   pr->nnodes = nnodes;
   StpVecReserve(scip, pr->currneighbors, pr->maxdegree + 1);
   StpVecReserve(scip, pr->firstneighborcosts, pr->maxdegree + 1);
   StpVecReserve(scip, pr->firstneighbors, pr->maxdegree + 1);
   StpVecReserve(scip, pr->pathedges, SP_MAXLENGTH);
   StpVecReserve(scip, pr->visitednodes, SP_MAXLENGTH);

   return SCIP_OKAY;
}


/** cleans temporary data */
static
void prClean(
   PR*                   pr                 /**< to be cleaned */
   )
{
   const int nvisited = StpVecGetSize(pr->visitednodes);

   for( int i = 0; i < nvisited; i++ )
   {
      const int node = pr->visitednodes[i];
      assert(node >= 0 && node < pr->nnodes && pr->nodes_index[node] != VNODES_UNSET);
      pr->nodes_index[node] = VNODES_UNSET;
   }

   StpVecClear(pr->firstneighbors);
   StpVecClear(pr->firstneighborcosts);
   StpVecClear(pr->currneighbors);
   StpVecClear(pr->pathedges);
   StpVecClear(pr->visitednodes);

#ifndef NDEBUG
   for( int i = 0; i < pr->nnodes; i++ )
      assert(pr->nodes_index[i] == VNODES_UNSET);

   for( int i = 0; i < SP_MAXNSTARTS; i++ )
      pr->sp_starts[i] = -1;
#endif
}


/** frees */
static
void prFree(
   SCIP*                 scip,               /**< SCIP data structure */
   PR**                  pathreplace         /**< to be freed */
   )
{
   PR* pr = *pathreplace;

   StpVecFree(scip, pr->pathedges);
   StpVecFree(scip, pr->currneighbors);
   StpVecFree(scip, pr->visitednodes);
   StpVecFree(scip, pr->firstneighborcosts);
   StpVecFree(scip, pr->firstneighbors);

   SCIPfreeMemoryArray(scip, &pr->firstneighbors_isKilled);
   SCIPfreeMemoryArray(scip, &pr->sp_dists);
   SCIPfreeMemoryArray(scip, &pr->sp_starts);
   SCIPfreeMemoryArray(scip, &pr->nodes_index);
   SCIPfreeMemoryArray(scip, &pr->nodeindices_isPath);
   SCIPfreeMemoryArray(scip, &pr->nodes_isTerm);
   reduce_sdprofitFree(scip, &pr->sdprofit);

   SCIPfreeMemory(scip, pathreplace);
}


/** is execution of path replacement reduction method promising? */
static
SCIP_Bool prIsPromising(
   const GRAPH*          g                   /**< graph structure */
   )
{
   // todo
   return TRUE;
}


/** tries to eliminate edges of path starting from edge */
static inline
SCIP_RETCODE processPath(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   startedge,          /**< edge to start from (head) */
   PR*                   pr,                 /**< path replacement data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  nelims              /**< pointer to number of reductions */
   )
{
   SCIP_Bool pathIsExtendable = TRUE;
   SCIP_Bool pathIsRedundant = FALSE;
   const int tail = g->tail[startedge];
   const int head = g->head[startedge];
   const int maxdeg = pr->maxdegree;
   assert(StpVecGetSize(pr->pathedges) == 0);

   if( g->grad[tail] > maxdeg || g->grad[head] > maxdeg || pr->nodes_isTerm[tail] || pr->nodes_isTerm[head] )
      return SCIP_OKAY;

   if( g->grad[tail] <= 1 )
      return SCIP_OKAY;

   pathExendPrepare(scip, g, startedge, pr);

   while( pathIsExtendable )
   {
      assert(!pathIsRedundant);
      pathExend(scip, g, pr, &pathIsExtendable, &pathIsRedundant);

      if( pathIsRedundant )
      {
         SCIPdebugMessage("deleting edge %d-%d \n", g->tail[startedge], g->head[startedge]);

         deleteEdge(scip, startedge, g, pr);

         (*nelims)++;
         break;
      }
   }

   prClean(pr);

   return SCIP_OKAY;
}


/** executes reduction method */
static
SCIP_RETCODE pathreplaceExec(
   SCIP*                 scip,               /**< SCIP data structure */
   PR*                   pr,                 /**< path replacement data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  nelims              /**< pointer to number of reductions */
   )
{
   const int nedges = graph_get_nEdges(g);
   const SCIP_Bool isPc = pr->probIsPc;
   SD* const sddata = pr->useSd ? pr->extperma->distdata_default->sdistdata : NULL;
   // todo have edge stack?

   for( int e = 0; e < nedges; e++ )
   {
      /* deleted edge? */
      if( g->oeat[e] == EAT_FREE )
         continue;

      /* dummy edge? */
      if( isPc && graph_pc_edgeIsExtended(g, e) )
         continue;

      if( sddata && reduce_sdgraphEdgeIsInMst(sddata->sdgraph, e) )
         continue;

      SCIP_CALL( processPath(scip, e, pr, g, nelims) );

      // todo if edgestack: check whether edge was killed, otherwise go different!
      // afterwards deactivate edge!
   }


   return SCIP_OKAY;
}

/*
 * Interface methods
 */


/** paths elimination while using special distances
 * NOTE: SD-repair needs to be set-up! */
SCIP_RETCODE reduce_pathreplaceExt(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   int*                  nelims              /**< pointer to number of reductions */
   )
{
   PR* pathreplace;
   SCIP_Bool hasDcsr = (g->dcsr_storage != NULL);

   assert(scip && g && extperma);

   if( !prIsPromising(g) )
   {
      return SCIP_OKAY;
   }

   if( !hasDcsr )
      SCIP_CALL( graph_init_dcsr(scip, g) );

   SCIP_CALL( prInit(scip, g, extperma, &pathreplace) );
   SCIP_CALL( pathreplaceExec(scip, pathreplace, g, nelims) );
   prFree(scip, &pathreplace);

   if( !hasDcsr )
      graph_free_dcsr(scip, g);

   if( *nelims > 0 )
   {
      deletenodesDeg1(scip, extperma, g);
   }

   assert(graph_valid(scip, g));

   return SCIP_OKAY;
}


/** paths elimination */
SCIP_RETCODE reduce_pathreplace(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  nelims              /**< pointer to number of reductions */
   )
{
   PR* pathreplace;

   if( !prIsPromising(g) )
   {
      return SCIP_OKAY;
   }

   SCIP_CALL( graph_init_dcsr(scip, g) );
   SCIP_CALL( prInit(scip, g, NULL, &pathreplace) );

   SCIP_CALL( pathreplaceExec(scip, pathreplace, g, nelims) );

   prFree(scip, &pathreplace);
   graph_free_dcsr(scip, g);

   if( *nelims > 0 )
   {
      reduce_nodesDeg1(scip, g);
   }

   assert(graph_valid(scip, g));

   return SCIP_OKAY;
}
