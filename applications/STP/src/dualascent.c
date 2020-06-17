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

/**@file   cons_stp.c
 * @brief  Dual-ascent dual heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file includes dual-ascent for classic Steiner tree and some variants.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include "scip/cons_linear.h"
#include "dualascent.h"
#include "probdata_stp.h"
#include "graph.h"
#include "heur_ascendprune.h"
#include "scip/scip.h"
#include "scip/misc.h"


#define DEFAULT_DAMAXDEVIATION  0.25  /**< max deviation for dual ascent */
#define DA_MAXDEVIATION_LOWER   0.01  /**< lower bound for max deviation for dual ascent */
#define DA_MAXDEVIATION_UPPER   0.9   /**< upper bound for max deviation for dual ascent */
#define DA_EPS                  (5e-7)

/* do depth-first search */
#define DFS

#ifndef RESTRICT
#define RESTRICT restrict
#endif

#ifdef BITFIELDSARRAY
#define ARRLENGTH 32
#define SetBit(Arr, pos)     ( Arr[(pos/ARRLENGTH)] |= (1 << (pos%ARRLENGTH)) )
#define CleanBit(Arr, pos)   ( Arr[(pos/ARRLENGTH)] &= ~(1 << (pos%ARRLENGTH)) )
#define BitTrue(Arr, pos)    ( Arr[(pos/ARRLENGTH)] & (1 << (pos%ARRLENGTH)) )
#endif



/** internal data for path based dual-ascent */
typedef struct dual_ascent_paths
{
   DIJK*                 dijklimited;        /**< Dijkstra data */
   SCIP_Real*            costs_reversed;     /**< reversed edge costs */
   SCIP_Real             distlimit;          /**< limit */
   int                   startnode;          /**< node */
} DAPATHS;




/**@} */


/**@name Local methods
 *
 * @{
 */

#ifndef NDEBUG
/** can all terminal be reached? */
static
SCIP_Bool allTermsReachable(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   const SCIP_Real*      redcost             /**< array to store reduced costs */
   )
{
   int* RESTRICT queue;
   STP_Bool* RESTRICT scanned;
   int qsize;
   const int nnodes = graph_get_nNodes(g);
   const int root = g->source;
   int termscount;

   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &queue, nnodes ) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &scanned, nnodes) );

   BMSclearMemoryArray(scanned, nnodes);

   termscount = 1;
   qsize = 0;
   scanned[root] = TRUE;
   queue[qsize++] = root;

   /* DFS */
   while( qsize > 0 )
   {
      const int k = queue[--qsize];
      scanned[k] = TRUE;

      for( int a = g->outbeg[k]; a != EAT_LAST; a = g->oeat[a] )
      {
         const int head = g->head[a];

         if( SCIPisZero(scip, redcost[a]) )
         {
            /* vertex not visited yet? */
            if( !scanned[head] )
            {
               scanned[head] = TRUE;
               queue[qsize++] = head;

               if( Is_term(g->term[head]) )
                  termscount++;
            }
         }
      }
   }

   SCIPfreeMemoryArray(scip, &scanned);
   SCIPfreeMemoryArray(scip, &queue);

  // printf("%d vs %d \n", termscount, g->terms);

   return (termscount == g->terms);
}



#endif


/** sets shortest path parameters: start node and distance limit */
static
void dapathsSetRunParams(
   const GRAPH*          transgraph,         /**< transformed SAP graph */
   DAPATHS*              dapaths             /**< to be initialized */
   )
{
   int start = -1;
   SCIP_Real maxprize = 0.0;
   const int nnodes = graph_get_nNodes(transgraph);
   const int root = transgraph->source;

   for( int i = 0; i < nnodes; i++ )
   {
      if( !Is_term(transgraph->term[i]) )
         continue;

      if( i == root )
         continue;

      assert(transgraph->grad[i] == 2);

      for( int e = transgraph->inpbeg[i]; e != EAT_LAST; e = transgraph->ieat[e] )
      {
         if( GT(transgraph->cost[e], maxprize) )
         {
            start = i;
            maxprize = transgraph->cost[e];
         }
      }
   }

   assert(graph_knot_isInRange(transgraph, start ));
   assert(GT(maxprize, 0.0));

   dapaths->startnode = start;

   assert(transgraph->outbeg[root] >= 0);
   dapaths->distlimit = transgraph->cost[transgraph->outbeg[root]];

#ifndef NDEBUG
   assert(GT(dapaths->distlimit, 0.0));
   assert(LT(dapaths->distlimit, FARAWAY));

   for( int e = transgraph->outbeg[root]; e != EAT_LAST; e = transgraph->oeat[e] )
   {
      assert(EQ(transgraph->cost[e], dapaths->distlimit));
   }
#endif

   //printf("maxprize=%f \n", maxprize);
   //printf("distlimit=%f \n", dapaths->distlimit);
}


/** initializes */
static
SCIP_RETCODE dapathsInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          transgraph,         /**< transformed SAP graph */
   DAPATHS*              dapaths             /**< to be initialized */
   )
{
   const int nedges = graph_get_nEdges(transgraph);
   SCIP_CALL( graph_dijkLimited_init(scip, transgraph, &(dapaths->dijklimited)) );

   SCIP_CALL( SCIPallocBufferArray(scip, &(dapaths->costs_reversed), nedges) );
   graph_get_edgeRevCosts(transgraph, dapaths->costs_reversed);

   dapathsSetRunParams(transgraph, dapaths);

   return SCIP_OKAY;
}


/** runs */
static
void dapathsRunShortestPaths(
   const GRAPH*          transgraph,         /**< transformed SAP graph */
   DAPATHS*              dapaths             /**< to be initialized */
   )
{
   const SCIP_Real* revcosts = dapaths->costs_reversed;
   const int start = dapaths->startnode;
   const SCIP_Real distlimit = dapaths->distlimit;
   DIJK* dijklimited = dapaths->dijklimited;

   assert(GT(distlimit, 0.0));

   graph_pathLimitedExec(transgraph, revcosts, start, distlimit, dijklimited);
}


/** computes reduced costs */
static
void dapathsComputeRedCosts(
   const GRAPH*          transgraph,         /**< transformed SAP graph */
   const DAPATHS*        dapaths,            /**< to be initialized */
   SCIP_Real* RESTRICT   redcost,            /**< array to store reduced costs */
   SCIP_Real*            objval             /**< pointer to store (dual) objective value */
)
{
   const int nnodes = graph_get_nNodes(transgraph);
   const int nedges = graph_get_nEdges(transgraph);
   const DIJK* dijklimited = dapaths->dijklimited;
   const SCIP_Real* const node_distance = dijklimited->node_distance;
   const SCIP_Real distlimit = dapaths->distlimit;

   BMScopyMemoryArray(redcost, transgraph->cost, nedges);

   for( int i = 0; i < nnodes; ++i )
   {
      const SCIP_Real dist_i = MIN(node_distance[i], distlimit);
      for( int e = transgraph->outbeg[i]; e != EAT_LAST; e = transgraph->oeat[e] )
      {
         const int j = transgraph->head[e];
         const SCIP_Real dist_j = node_distance[j];

         if( LT(dist_j, distlimit) )
         {
            const SCIP_Real offset = MAX(0.0, dist_i - dist_j);
            assert(LE(offset, redcost[e]));

            redcost[e] -= offset;
         }
      }
   }

   *objval = distlimit;
}


/** frees */
static
void dapathsFreeMembers(
   SCIP*                 scip,               /**< SCIP data structure */
   DAPATHS*              dapaths             /**< to be initialized */
   )
{
   SCIPfreeBufferArray(scip, &(dapaths->costs_reversed));
   graph_dijkLimited_free(scip, &(dapaths->dijklimited));
}


/** returns whether node realtail is active or leads to active node other than dfsbase */
static inline
SCIP_Bool is_active(
   const int*            active,             /**< active nodes array */
   int                   realtail,           /**< vertex to start from */
   int                   dfsbase             /**< DFS source vertex */
   )
{
   int curr;

   for( curr = active[realtail]; curr != 0 && curr != dfsbase + 1; curr = active[curr - 1] )
   {
      assert(curr >= 0);
   }

   return (curr == 0);
}


/**@} */

/**@name Interface methods
 *
 * @{
 */

/** initializes */
static
SCIP_RETCODE initDualAscent(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   const int* RESTRICT   start,              /**< CSR start array [0,...,nnodes] */
   const int* RESTRICT   edgearr,            /**< CSR ancestor edge array */
   int                   root,               /**< the root */
   SCIP_Bool             is_pseudoroot,      /**< is the root a pseudo root? */
   int                   ncsredges,          /**< number of CSR edges */
   int*                  gmark,              /**< array for marking nodes */
   int* RESTRICT         active,             /**< active vertices mark */
   SCIP_PQUEUE*          pqueue,             /**< priority queue */
   GNODE*                gnodearr,           /**< array containing terminal nodes*/
   SCIP_Real* RESTRICT   rescap,             /**< residual capacity */
   SCIP_Real*            dualobj,            /**< dual objective */
   int*                  augmentingcomponent /**< augmenting component */
)
{
   const int nnodes = g->knots;

   *dualobj = 0.0;
   *augmentingcomponent = -1;

   for( int i = 0; i < ncsredges; i++ )
      rescap[i] = g->cost[edgearr[i]];

   /* mark terminals as active, add all except root to pqueue */
   for( int i = 0, termcount = 0; i < nnodes; i++ )
   {
      if( !Is_term(g->term[i]) )
      {
         active[i] = -1;
         continue;
      }

      active[i] = 0;
      assert(g->grad[i] > 0);

      if( i != root )
      {
         assert(termcount < g->terms - 1);
         assert(gnodearr);

         gnodearr[termcount].number = i;
         gnodearr[termcount].dist = g->grad[i];

         /* for variants with dummy terminals */
         if( g->grad[i] == 2 )
         {
            int a;

            for( a = g->inpbeg[i]; a != EAT_LAST; a = g->ieat[a] )
               if( SCIPisZero(scip, g->cost[a]) )
                  break;

            if( a != EAT_LAST )
            {
               const int tail = g->tail[a];
               gnodearr[termcount].dist += g->grad[tail] - 1;

               if( is_pseudoroot )
               {
                  for( a = g->inpbeg[tail]; a != EAT_LAST; a = g->ieat[a] )
                  {
                     if( SCIPisZero(scip, g->cost[a]) )
                     {
                        gnodearr[termcount].dist += g->grad[g->tail[a]] - 1;
                     }
                  }
               }
            }

            assert(gnodearr[termcount].dist > 0);
         }

         SCIP_CALL(SCIPpqueueInsert(pqueue, &(gnodearr[termcount])));

         if( *augmentingcomponent == -1 )
            *augmentingcomponent = i;

         termcount++;
      }
   }

   assert(*augmentingcomponent >= 0);

   for( int i = 0; i < nnodes + 1; i++ )
      gmark[i] = FALSE;

   return SCIP_OKAY;
}

/** dual ascent heuristic */
SCIP_RETCODE dualascent_exec(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real* RESTRICT   redcost,            /**< array to store reduced costs or NULL */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Bool             addcuts,            /**< should dual ascent add Steiner cuts? */
   SCIP_Bool             ascendandprune,     /**< should the ascent-and-prune heuristic be executed? */
   const int*            result,             /**< solution array or NULL */
   int                   root,               /**< the root */
   SCIP_Bool             is_pseudoroot,      /**< is the root a pseudo root? */
   SCIP_Real             damaxdeviation      /**< maximum deviation for dual-ascent ( -1.0 for default) */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_PQUEUE* pqueue;
   SCIP_VAR** vars;
   SCIP_Real* RESTRICT rescap;
   GNODE* gnodearr;
   int* RESTRICT edgearr;
   int* RESTRICT tailarr;
   int* RESTRICT start;
   int* RESTRICT stackarr;
   int* RESTRICT cutverts;
   int* RESTRICT unsatarcs;
   int* RESTRICT unsattails;
   int* RESTRICT gmark;
   int* RESTRICT active;
   SCIP_Real dualobj;
   SCIP_Real currscore;
   const SCIP_Real maxdeviation = (damaxdeviation > 0.0) ? damaxdeviation : DEFAULT_DAMAXDEVIATION;
   const int nnodes = g->knots;
   const int nterms = g->terms;
   const int nedges = g->edges;
   int ncsredges;
   int norgcutverts;
   int stacklength;
   int augmentingcomponent;
   const SCIP_Bool addconss = (SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE);

   /* should currently not  be activated */
   assert(addconss || !addcuts);
   assert(g != NULL);
   assert(scip != NULL);
   assert(objval != NULL);
   assert(Is_term(g->term[root]));
   assert(maxdeviation >= DA_MAXDEVIATION_LOWER && maxdeviation <= DA_MAXDEVIATION_UPPER);
   assert(damaxdeviation == -1.0 || damaxdeviation > 0.0);

   if( nnodes == 1 )
      return SCIP_OKAY;

   if( addcuts )
   {
      vars = SCIPprobdataGetVars(scip);
      assert(vars != NULL);

      if( !addconss )
      {
         conshdlr = SCIPfindConshdlr(scip, "stp");
         assert(conshdlr != NULL);
      }
   }
   else
   {
      vars = NULL;
   }

   /* if specified root is not a terminal, take default root */
   if( root < 0 || !Is_term(g->term[root]) )
      root = g->source;

#ifdef BITFIELDSARRAY
   u_int32_t* bitarr;
   SCIP_CALL( SCIPallocBufferArray(scip, &bitarr, nedges / ARRLENGTH + 1) );
#endif

   stacklength = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &unsattails, nedges) );

   if( redcost == NULL )
      SCIP_CALL( SCIPallocBufferArray(scip, &rescap, nedges) );
   else
      rescap = redcost;

   SCIP_CALL( SCIPallocBufferArray(scip, &cutverts, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &unsatarcs, nedges) );

   if( nterms > 1 )
      SCIP_CALL( SCIPallocMemoryArray(scip, &gnodearr, nterms - 1) );
   else
      gnodearr = NULL;

   SCIP_CALL( SCIPpqueueCreate(&pqueue, nterms, 2.0, GNODECmpByDist) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &active, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &edgearr, nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &tailarr, nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &start, nnodes + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &gmark, nnodes + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &stackarr, nnodes) );

   /* fill auxiliary adjacent vertex/edges arrays */
   graph_get_csr(g, edgearr, tailarr, start, &ncsredges);

   /* initialize priority queue and res. capacity */
   SCIP_CALL( initDualAscent(scip, g, start, edgearr, root, is_pseudoroot, ncsredges, gmark, active, pqueue,
         gnodearr, rescap, &dualobj, &augmentingcomponent) );

   /* mark whether an arc is satisfied (has capacity 0) */
   for( int i = 0; i < ncsredges; i++ )
   {
#ifdef BITFIELDSARRAY
      if( SCIPisZero(scip, rescap[i]) )
         SetBit(bitarr, i);
      else
         CleanBit(bitarr, i);
#else
      if( rescap[i] == 0.0 )
      {
         if( active[tailarr[i] - 1] == 0 )
            tailarr[i] = 0;
         else
            tailarr[i] *= -1;
      }
#endif
   }

   norgcutverts = 0;

   /* (main) dual ascent loop */
   while( SCIPpqueueNElems(pqueue) > 0 && !SCIPisStopped(scip) )
   {
      /* get active vertex of minimum score */
      GNODE* const gnodeact = (GNODE*) SCIPpqueueRemove(pqueue);
      const SCIP_Real prio1 = gnodeact->dist;
      const SCIP_Real prio2 = (SCIPpqueueNElems(pqueue) > 0) ? ((GNODE*) SCIPpqueueFirst(pqueue))->dist : FARAWAY;
      const int v = gnodeact->number;
      SCIP_Real degsum = g->grad[v];
      int ncutverts = 0;
      int nunsatarcs = 0;

      SCIP_Bool firstrun = TRUE;

      SCIPdebugMessage("DA: START WITH v %d prio1 %f prio2 %f \n", v, prio1, prio2);

      /* perform augmentation as long as priority of root component does not exceed max deviation */
      for( ; ; )
      {
         assert(stacklength == 0);

         /* 1. step: BFS from v (or connected component) on saturated, incoming arcs */

         if( firstrun )
         {
            firstrun = FALSE;
            gmark[v + 1] = TRUE;
            cutverts[ncutverts++] = v;
            assert(stacklength < nnodes);
            stackarr[stacklength++] = v;
         }
         /* not in first processing of root component: */
         else
         {
            for( int i = norgcutverts; i < ncutverts; i++ )
            {
               const int s = cutverts[i];

               assert(gmark[s + 1]);
               assert(active[s] != 0);
               assert(stacklength < nnodes);

               stackarr[stacklength++] = s;
            }
         }
#ifdef DFS
         while( stacklength )
         {
            const int node = stackarr[--stacklength];
#else
         for( int n = 0; n < stacklength; n++ )
         {
            int end;

            assert(n < nnodes);
            node = stackarr[n];
#endif

            /* traverse incoming arcs */
            for( int i = start[node], end = start[node + 1]; i != end; i++ )
            {
               int tail = tailarr[i];

               /* zero reduced-cost arc? */
               if( tail <= 0 )
               {
                  tail *= -1;
                  if( !gmark[tail] )
                  {
                     /* if an active vertex has been hit (other than v), break */
                     if( 0 == tail )
                     {
                        const int realtail = g->tail[edgearr[i]];

                        /* v should not be processed */
                        if( realtail == v )
                           continue;

                        /* is realtail active or does realtail lead to an active vertex other than v? */
                        if( is_active(active, realtail, v) )
                        {
                           active[v] = realtail + 1;
                           stacklength = 0;
                           goto ENDOFLOOP;
                        }

                        tail = realtail + 1;

                        /* have we processed tail already? */
                        if( gmark[tail] )
                           continue;
                     }

                     assert(tail > 0);

                     gmark[tail] = TRUE;
                     tail--;
                     cutverts[ncutverts++] = tail;
                     degsum += g->grad[tail];

                     assert(stacklength < nnodes);
                     stackarr[stacklength++] = tail;
                  } /* marked */
               } /* zero reduced-cost arc */
               else if( !gmark[tail] )
               {
                  unsattails[nunsatarcs] = tail;
                  unsatarcs[nunsatarcs++] = i;
               }
            }
         }
#ifndef DFS
         stacklength = 0;
#endif
         currscore = degsum - (ncutverts - 1);

         /* guiding solution provided? */
         if( result != NULL )
         {
            int nsolarcs = 0;
            for( int i = 0; i < nunsatarcs; i++ )
            {
               const int a = unsatarcs[i];

               assert(tailarr[a] > 0);

               if( !(gmark[tailarr[a]]) )
               {
                  if( result[edgearr[a]] == CONNECT )
                     nsolarcs++;
               }
            }

            assert(nsolarcs > 0);
            assert(currscore <= nedges);

            if( nsolarcs > 1 )
              currscore += (SCIP_Real) ((nsolarcs - 1) * (g->knots * 2.0));
         }
         else
         {
            assert(SCIPisGE(scip, currscore, prio1));
         }

         SCIPdebugMessage("DA: deviation %f \n", (currscore - prio1) / prio1);
         SCIPdebugMessage("DA: currscore %f prio1 %f prio2 %f \n", currscore, prio1, prio2);

         /* augmentation criteria met? */
         if( ((currscore - prio1) / prio1) <= maxdeviation || currscore <= prio2 )
         {
            SCIP_CONS* cons = NULL;
            SCIP_ROW* row = NULL;

            int shift = 0;
            SCIP_Real min = FARAWAY;
            SCIP_Bool isactive = FALSE;

            /* 2. step: get minimum residual capacity among cut-arcs */

            /* adjust array of unsatisfied arcs */

            for( int i = 0; i < nunsatarcs; i++ )
            {
               const int tail = unsattails[i];

               if( gmark[tail] )
               {
                  shift++;
               }
               else
               {
                  const int a = unsatarcs[i];

                  assert(tailarr[a] > 0);
                  assert(rescap[a] > 0);

                  if( rescap[a] < min )
                     min = rescap[a];
                  if( shift )
                  {
                     unsattails[i - shift] = tail;
                     unsatarcs[i - shift] = a;
                  }
               }
            }

            assert(SCIPisLT(scip, min, FARAWAY));
            nunsatarcs -= shift;

            norgcutverts = ncutverts;

            /* 3. step: perform augmentation */

            /* create constraints/cuts ? */
            if( addcuts )
            {
               if( addconss )
               {
                  SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "da", 0, NULL, NULL,
                        1.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
               }
               else
               {
                  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "da", 1.0,
                        SCIPinfinity(scip), FALSE, FALSE, TRUE) );

                  SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
               }
            }

            shift = 0;

            /* update (dual) objective */
            dualobj += min;

            for( int i = 0; i < nunsatarcs; i++ )
            {
               const int a = unsatarcs[i];
               assert(a >= 0);

               if( addcuts )
               {
                  assert(vars != NULL);

                  if( addconss )
                     SCIP_CALL( SCIPaddCoefLinear(scip, cons, vars[edgearr[a]], 1.0) );
                  else
                     SCIP_CALL( SCIPaddVarToRow(scip, row, vars[edgearr[a]], 1.0) );
               }
               rescap[a] -= min;

               assert(SCIPisGE(scip, rescap[a], 0.0));

               if( rescap[a] <= DA_EPS )
               {
                  int tail = unsattails[i];

                  rescap[a] = 0.0;

                  assert(tail > 0);
                  assert(tailarr[a] > 0);

                  tailarr[a] *= -1;

                  if( active[tail - 1] >= 0 && is_active(active, tail - 1, v) )
                  {
                     assert(tail - 1 != v);
                     tailarr[a] = 0;
                     if( !isactive )
                     {
                        isactive = TRUE;
                        active[v] = tail;
                     }
                  }


                  if( !(gmark[tail])  )
                  {
                     assert(tail != 0);

                     gmark[tail] = TRUE;
                     tail--;
                     degsum += g->grad[tail];
                     cutverts[ncutverts++] = tail;
                  }

                  shift++;
               }
               else if( shift )
               {
                  unsattails[i - shift] = unsattails[i];
                  unsatarcs[i - shift] = a;
               }
            }

            if( addcuts )
            {
               if( addconss )
               {
                  SCIP_CALL( SCIPaddCons(scip, cons) );
                  SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               }
               else
               {
                  SCIP_Bool infeasible;

                  SCIP_CALL( SCIPflushRowExtensions(scip, row) );
                  SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
                  SCIP_CALL( SCIPreleaseRow(scip, &row) );

                  assert(!infeasible);
               }
            }

            if( isactive )
            {
               stacklength = 0;
               goto ENDOFLOOP;
            }
            nunsatarcs -= shift;

            continue;
         }
         else
         {
            SCIP_Bool insert = TRUE;

            if( is_pseudoroot )
            {
               int i = start[v];
               const int end = start[v + 1];

               assert(end - i == 2);

               for( ; i != end; i++ )
                  if( rescap[i] != 0.0 )
                     break;

               if( i == end )
               {
                  if( augmentingcomponent == -1 )
                     augmentingcomponent = v;

                  if( augmentingcomponent != v )
                     insert = FALSE;
               }
            }

            if( insert )
            {
               /* reinsert active vertex */
               gnodeact->dist = currscore;
               SCIP_CALL( SCIPpqueueInsert(pqueue, gnodeact) );
            }
         }

         ENDOFLOOP:

         for( int i = 0; i < ncutverts; i++ )
            gmark[cutverts[i] + 1] = FALSE;

         for( int i = 0; i < nnodes + 1; i++ )
         {
            assert(!gmark[i]);
         }

         break;
      } /* augmentation loop */
   } /* dual ascent loop */

   SCIPdebugMessage("DA: dualglobal: %f \n", dualobj);
   *objval = dualobj;

   for( int i = ncsredges; i < nedges; i++ )
   {
      edgearr[i] = i;
      rescap[i] = g->cost[i];
   }

   /* re-extend rescap array */
   for( int i = 0; i < ncsredges; i++ )
   {
      if( edgearr[i] != i  )
      {
         SCIP_Real bufferedval = rescap[i];
         int a = i;

         rescap[i] = g->cost[i];
         while( edgearr[a] != a )
         {
            const int shift = edgearr[a];
            const SCIP_Real min = rescap[shift];

            rescap[shift] = bufferedval;
            bufferedval = min;
            edgearr[a] = a;
            a = shift;
         }
      }
   }

#ifdef BITFIELDSARRAY
   SCIPfreeBufferArray(scip, &bitarr);
#endif

   SCIPfreeMemoryArray(scip, &stackarr);
   SCIPfreeMemoryArray(scip, &gmark);
   SCIPfreeMemoryArray(scip, &start);
   SCIPfreeMemoryArray(scip, &tailarr);
   SCIPfreeMemoryArray(scip, &edgearr);
   SCIPfreeMemoryArray(scip, &active);

   SCIPpqueueFree(&pqueue);
   SCIPfreeMemoryArrayNull(scip, &gnodearr);

   /* call Ascend-And-Prune? */
   if( ascendandprune )
   {
       SCIP_Bool success;

       SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, g, rescap, unsatarcs, root, &success, TRUE) );
   }

   SCIPfreeBufferArray(scip, &unsatarcs);
   SCIPfreeBufferArray(scip, &cutverts);

   if( redcost == NULL )
      SCIPfreeBufferArray(scip, &rescap);

   SCIPfreeBufferArray(scip, &unsattails);

   return SCIP_OKAY;
}

/** dual ascent heuristic for PCSPG and MWCSP */
SCIP_RETCODE dualascent_execPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            redcost,            /**< array to store reduced costs or NULL */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Bool             addcuts,            /**< should dual-ascent add Steiner cuts? */
   SCIP_Bool             ascendandprune,     /**< perform ascend-and-prune and add solution? */
   int                   nruns               /**< number of dual ascent runs */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_PQUEUE* pqueue;
   SCIP_VAR** vars;
   GRAPH* transgraph;
   SCIP_Real min;
   SCIP_Real prio1;
   SCIP_Real offset;
   SCIP_Real dualobj;
   SCIP_Real currscore;
   SCIP_Real maxdeviation;
   SCIP_Real* rescap;
   GNODE* gnodeact;
   GNODE** gnodearr;
   int s;
   int i;
   int k;
   int v;
   int a;
   int tail;
   int pnode;
   int shift;
   int root;
   int nnodes;
   int nterms;
   int nedges;
   int degsum;
   int ncutverts;
   int pseudoroot;
   int nunsatarcs;
   int stacklength;
   int norgcutverts;
   int* cutverts;
   int* stackarr;
   STP_Bool* origedge;
   int* unsatarcs;
   STP_Bool firstrun;
   STP_Bool* sat;
   STP_Bool* active;
   const SCIP_Bool addconss = (SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE);

   /* should currently not  be activated */
   assert(addconss || !addcuts);

   assert(g != NULL);
   assert(scip != NULL);
   assert(nruns >= 0);
   assert(objval != NULL);

   if( g->knots == 1 )
      return SCIP_OKAY;

   if( addcuts )
   {
      vars = SCIPprobdataGetVars(scip);
      assert(vars != NULL);
      if( !addconss )
      {
         conshdlr = SCIPfindConshdlr(scip, "stp");
         assert(conshdlr != NULL);
      }
   }
   else
   {
      vars = NULL;
   }

   root = g->source;
   degsum = 0;
   offset = 0.0;
   dualobj = 0.0;

   ncutverts = 0;
   norgcutverts = 0;
   maxdeviation = DEFAULT_DAMAXDEVIATION;

   SCIP_CALL( graph_transPcGetSap(scip, g, &transgraph, &offset) );

   nnodes = transgraph->knots;
   nedges = transgraph->edges;
   nterms = transgraph->terms;
   pseudoroot = nnodes - 1;

   if( redcost == NULL )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &rescap, nedges) );
   }
   else
   {
      rescap = redcost;
   }

   stacklength = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, &stackarr, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sat, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &active, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutverts, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &unsatarcs, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &origedge, nedges) );

   for( i = 0; i < nedges; i++ )
      if( !Is_term(transgraph->term[transgraph->tail[i]]) && transgraph->head[i] == pseudoroot )
         origedge[i] = FALSE;
      else if( transgraph->tail[i] == pseudoroot && !Is_term(transgraph->term[transgraph->head[i]])  )
         origedge[i] = FALSE;
      else
         origedge[i] = TRUE;

   for( i = 0; i < nterms - 1; i++ )
   {
      SCIP_CALL( SCIPallocBuffer(scip, &gnodearr[i]) ); /*lint !e866*/
   }

   SCIP_CALL( SCIPpqueueCreate( &pqueue, nnodes, 2.0, GNODECmpByDist) );

   k = 0;
   /* mark terminals as active, add all except root to pqueue */
   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(transgraph->term[i]) )
      {
         active[i] = TRUE;
         assert(transgraph->grad[i] > 0);
         if( i != root  )
         {
            gnodearr[k]->number = i;
            gnodearr[k]->dist = transgraph->grad[i];

            for( a = transgraph->inpbeg[i]; a != EAT_LAST; a = transgraph->ieat[a] )
               if( SCIPisEQ(scip, transgraph->cost[a], 0.0) )
                  break;

            if( a != EAT_LAST )
               gnodearr[k]->dist += transgraph->grad[transgraph->tail[a]] - 1;

            assert(gnodearr[k]->dist > 0);

            SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[k++]) );
         }
      }
      else
      {
         active[i] = FALSE;
      }
      transgraph->mark[i] = FALSE;
   }

   for( i = 0; i < nedges; i++ )
   {
      rescap[i] = transgraph->cost[i];
      if( SCIPisZero(scip, rescap[i]) )
         sat[i] = TRUE;
      else
         sat[i] = FALSE;
   }

   /* dual ascent loop */
   while( SCIPpqueueNElems(pqueue) > 0 && !SCIPisStopped(scip) )
   {
      /* get active vertex of minimum score */
      gnodeact = (GNODE*) SCIPpqueueRemove(pqueue);

      v = gnodeact->number;
      prio1 = gnodeact->dist;

      firstrun = TRUE;
      nunsatarcs = 0;

      /* perform augmentation as long as ... */
      for( ; ; )
      {
         assert(stacklength == 0);
         /* 1. step: BFS from v (or connected component) on saturated, incoming arcs */

         if( firstrun )
         {
            degsum = transgraph->grad[v];
            ncutverts = 0;
            firstrun = FALSE;
            nunsatarcs = 0;
            transgraph->mark[v] = TRUE;
            cutverts[ncutverts++] = v;
            stackarr[stacklength++] = v;
         }
         /* not in first processing of root component: */
         else
         {
            for( i = norgcutverts; i < ncutverts; i++ )
            {
               s = cutverts[i];
               assert(transgraph->mark[s]);
               if( active[s] )
               {
                  active[v] = FALSE;
                  stacklength = 0;
                  goto ENDOFLOOP;
               }

               stackarr[stacklength++] = s;
            }
         }

         while( stacklength )
         {
            pnode = stackarr[--stacklength];

            /* traverse incoming arcs */
            for( a = transgraph->inpbeg[pnode]; a != EAT_LAST; a = transgraph->ieat[a] )
            {
               tail = transgraph->tail[a];
               if( sat[a] )
               {
                  if( !transgraph->mark[tail] )
                  {
                     /* if an active vertex has been hit, break */
                     if( active[tail] )
                     {
                        active[v] = FALSE;
                        stacklength = 0;
                        goto ENDOFLOOP;
                     }

                     degsum += transgraph->grad[tail];
                     transgraph->mark[tail] = TRUE;
                     cutverts[ncutverts++] = tail;
                     stackarr[stacklength++] = tail;
                  }
               }
               else if( !transgraph->mark[tail] )
               {
                  unsatarcs[nunsatarcs++] = a;
               }
            }
         }

         currscore = degsum - (ncutverts - 1);

         assert(SCIPisGE(scip, currscore, prio1));

         /* augmentation criteria met? */
         if( SCIPisLE(scip, (currscore - prio1) / prio1, maxdeviation) || (SCIPpqueueNElems(pqueue) == 0) )
         {
            SCIP_Bool in = FALSE;
            SCIP_ROW* row;
            SCIP_CONS* cons = NULL;

            /* 2. pass: get minimum residual capacity among cut-arcs */

            /* adjust array of unsatisfied arcs */
            min = FARAWAY;
            shift = 0;

            for( i = 0; i < nunsatarcs; i++ )
            {
               a = unsatarcs[i];
               if( transgraph->mark[transgraph->tail[a]] )
               {
                  shift++;
               }
               else
               {

                  assert(!sat[a]);
                  if( SCIPisLT(scip, rescap[a], min) )
                     min = rescap[a];
                  if( shift != 0 )
                     unsatarcs[i - shift] = a;
               }
            }

            assert(SCIPisLT(scip, min, FARAWAY));
            nunsatarcs -= shift;

            if( nunsatarcs > 0)
               assert(!transgraph->mark[transgraph->tail[unsatarcs[nunsatarcs-1]]]);

            norgcutverts = ncutverts;


            /* 3. pass: perform augmentation */


            /* create constraint/row */

            if( addcuts )
            {
               if( addconss )
               {
                  SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "da", 0, NULL, NULL,
                        1.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );
               }
               else
               {
                  SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "da", 1.0, SCIPinfinity(scip), FALSE, FALSE, TRUE));
                  SCIP_CALL(SCIPcacheRowExtensions(scip, row));
               }
            }

            dualobj += min;
            for( i = 0; i < nunsatarcs; i++ )
            {
               a = unsatarcs[i];
               if( a == -1 )
                  continue;

               if( addcuts && origedge[a] )
               {
                  assert(vars != NULL);
                  assert(cons != NULL);

                  if( g->tail[a] == root && g->head[a] == v )
                     in = TRUE;

                  if( addconss )
                     SCIP_CALL( SCIPaddCoefLinear(scip, cons, vars[a], 1.0) );
                  else
                     SCIP_CALL( SCIPaddVarToRow(scip, row, vars[a], 1.0) );
               }
               rescap[a] -= min;

               assert(SCIPisGE(scip, rescap[a], 0.0));

               if( SCIPisEQ(scip, rescap[a], 0.0) )
               {
                  sat[a] = TRUE;
                  if( !(transgraph->mark[transgraph->tail[a]]) )
                  {
                     tail = transgraph->tail[a];
                     degsum += transgraph->grad[tail];
                     transgraph->mark[tail] = TRUE;
                     cutverts[ncutverts++] = tail;
                  }
               }
            }

            if( addcuts )
            {
               assert(vars != NULL);

               if( !in )
               {
                  for( i = g->outbeg[root]; i != EAT_LAST; i = g->oeat[i] )
                     if( g->head[i] == v )
                     {
                        if( addconss )
                           SCIP_CALL( SCIPaddCoefLinear(scip, cons, vars[i], 1.0) );
                        else
                           SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i], 1.0) );
                     }
               }

               if( addconss )
               {
                  SCIP_CALL( SCIPaddCons(scip, cons) );
                  SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               }
               else
               {
                  SCIP_Bool infeasible;
                  assert(row != NULL);

                  SCIP_CALL( SCIPflushRowExtensions(scip, row) );
                  SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
                  SCIP_CALL( SCIPreleaseRow(scip, &row) );

                  assert(!infeasible);
               }
            }

            continue;
         }
         else
         {
            /* reinsert active vertex */
            gnodeact->dist = currscore;
            SCIP_CALL( SCIPpqueueInsert(pqueue, gnodeact) );
         }

      ENDOFLOOP:

         for( i = 0; i < ncutverts; i++ )
            transgraph->mark[cutverts[i]] = FALSE;

         break;
      } /* augmentation loop */
   } /* dual ascent loop */


   *objval = dualobj + offset;
   SCIPdebugMessage("DA: dualglobal: %f \n", *objval + SCIPprobdataGetOffset(scip));

   /* call dual Ascend-And-Prune? */
   if( ascendandprune )
   {
      SCIP_Bool success;
      SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, g, rescap, unsatarcs, -1, &success, TRUE));
   }

   /* free memory */
   SCIPpqueueFree(&pqueue);

   for( i = nterms - 2; i >= 0; i-- )
      SCIPfreeBuffer(scip, &gnodearr[i]);

   SCIPfreeBufferArray(scip, &origedge);
   SCIPfreeBufferArray(scip, &unsatarcs);
   SCIPfreeBufferArray(scip, &cutverts);
   SCIPfreeBufferArray(scip, &gnodearr);
   SCIPfreeBufferArray(scip, &active);
   SCIPfreeBufferArray(scip, &sat);
   SCIPfreeBufferArray(scip, &stackarr);

   if( redcost == NULL )
      SCIPfreeBufferArray(scip, &rescap);

   graph_free(scip, &transgraph, TRUE);

   return SCIP_OKAY;
}


/** path based dual ascent heuristic */
SCIP_RETCODE dualascent_pathsPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          transgraph,         /**< transformed SAP graph */
   SCIP_Real* RESTRICT   redcost,            /**< array to store reduced costs */
   SCIP_Real*            objval,             /**< pointer to store (dual) objective value */
   const int*            result              /**< solution array or NULL */
)
{
   DAPATHS dapaths = { NULL, NULL, -1.0, -1 };

   assert(scip && transgraph && redcost && objval);

   SCIP_CALL( dapathsInit(scip, transgraph, &dapaths) );

   dapathsRunShortestPaths(transgraph, &dapaths);
   dapathsComputeRedCosts(transgraph, &dapaths, redcost, objval);

   dapathsFreeMembers(scip, &dapaths);

   assert(allTermsReachable(scip, transgraph, redcost));


   return SCIP_OKAY;
}



/**@} */
