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

/**@file   heur_tm.c
 * @brief  Shortest paths based primal heuristics for Steiner problems
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 * @author Michael Winkler
 *
 * This file implements several shortest paths based primal heuristics for Steiner problems, see
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" by
 * Gamrath, Koch, Maher, Rehfeldt and Shinano
 *
 * A list of all interface methods can be found in heur_tm.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>
#include "heur_tm.h"
#include "probdata_stp.h"
#include "portab.h"
#include "scip/misc.h"
#include <math.h>
#define HEUR_NAME             "TM"
#define HEUR_DESC             "shortest path based primal heuristics for Steiner trees"
#define HEUR_DISPCHAR         '+'
#define HEUR_PRIORITY         10000000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_EVALRUNS 25                  /**< number of runs */
#define DEFAULT_INITRUNS 100                 /**< number of initial runs */
#define DEFAULT_LEAFRUNS 15                  /**< number of runs at leafs */
#define DEFAULT_ROOTRUNS 50                  /**< number of runs at the root */
#define DEFAULT_DURINGLPFREQ 5               /**< frequency during LP solving */
#define DEFAULT_TYPE  0                      /**< heuristic to execute */
#define DEFAULT_RANDSEED 5                   /**< seed for pseudo-random functions */

#define AUTO        0
#define TM_SP       1
#define TM_VORONOI  2
#define TM_DIJKSTRA 3

#ifdef WITH_UG
extern
int getUgRank(void);
#endif

/*
 * Data structures
 */


/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          nlpiterations;      /**< number of total LP iterations*/
   SCIP_Longint          ncalls;             /**< number of total calls (of TM) */
   SCIP_Longint          nexecs;             /**< number of total executions (of TM) */
   SCIP_Real             hopfactor;          /**< edge multiplication factor for hop constrained problems */
   int                   stp_type;           /**< problem type */
   int                   evalruns;           /**< number of runs */
   int                   initruns;           /**< number of initial runs */
   int                   leafruns;           /**< number of runs at leafs */
   int                   rootruns;           /**< number of runs at the root */
   int                   duringlpfreq;       /**< frequency during LP solving */
   int                   type;               /**< Heuristic type: 0 automatic, 1 TM_SP, 2 TM_VORONOI, 3 TM_DIJKSTRA */
   int                   beststartnode;      /**< start node of the so far best found solution */
   unsigned int          randseed;           /**< seed value for random number generator */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   unsigned int          timing;             /**< timing for timing mask */
};

/*
 * Static methods
 */

/** information method for a parameter change of random seed */
static
SCIP_DECL_PARAMCHGD(paramChgdRandomseed)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int newrandseed;

   newrandseed = SCIPparamGetInt(param);

   heurdata = (SCIP_HEURDATA*)SCIPparamGetData(param);
   assert(heurdata != NULL);

   heurdata->randseed = (unsigned int)newrandseed;

   return SCIP_OKAY;
}


/** compute starting points among marked (w.r.t. g->mark) vertices for constructive heuristics */
void SCIPStpHeurTMCompStarts(
   GRAPH*                graph,              /**< graph data structure */
   int*                  starts,             /**< starting points array */
   int*                  runs               /**< pointer to number of runs */
   )
{
   int r;
   int k;
   int l;
   int root;
   int nruns;
   int nnodes;
   int nterms;
   int randval;

   assert(runs != NULL);
   assert(graph != NULL);
   assert(starts != NULL);

   nruns = *runs;
   root = graph->source;
   nnodes = graph->knots;
   nterms = graph->terms;

   r = 0;
   if( graph->mark[root] && nruns > 0 )
      starts[r++] = root;

   randval = nnodes - nterms;

   /* use non-isolated terminals as starting points for TM heuristic */
   for( k = 0; k < nnodes; k++ )
   {
      if( r >= nruns || r >= nterms )
         break;

      l = (k + randval) % nnodes;
      if( Is_term(graph->term[l]) && graph->mark[l] && l != root )
         starts[r++] = l;
   }

   /* fill empty slots randomly */
   for( k = 0; k < nnodes && r < nruns; k++ )
   {
      l = (k + randval) % nnodes;
      if( !Is_term(graph->term[l]) && graph->mark[l] )
         starts[r++] = l;
   }

   *runs = r;
}

/* prune the (rooted) prize collecting Steiner tree in such a way that all leaves are terminals */
SCIP_RETCODE SCIPStpHeurTMPrunePc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   int*                  result,             /**< ST edges (need to be set to UNKNOWN) */
   STP_Bool*             connected           /**< ST nodes */
   )
{
   PATH*  mst;
   int i;
   int j;
   int e1;
   int e2;
   int k1;
   int k2;
   int root;
   int count;
   int nnodes;
   SCIP_Bool rmw = (g->stp_type == STP_RMWCSP);
   SCIP_Bool rpc = (g->stp_type == STP_RPCSPG);
   nnodes = g->knots;
   root = g->source;

   if( rmw )
   {
      for( i = 0; i < nnodes; i++ )
      {
         if( connected[i] && g->mark[i] )
            g->mark[i] = TRUE;
         else
            g->mark[i] = FALSE;
      }
      if( !g->mark[root] )
      {
         int nedges = g->edges;
         for( i = 0; i < nedges; i++ )
            result[i] = CONNECT;
         return SCIP_OKAY;
      }
   }
   else
   {
      /* compute the MST, exclude all terminals */
      for( i = 0; i < nnodes; i++ )
      {
         if( connected[i] && !Is_term(g->term[i]) )
            g->mark[i] = TRUE;
         else
            g->mark[i] = FALSE;
      }
   }

   if( rpc )
   {
      g->mark[root] = TRUE;
   }
   else if( !rmw )
   {
      int a;
      for( a = g->outbeg[root]; a != EAT_LAST; a = g->oeat[a] )
      {
         i = g->head[a];
         if( !Is_term(g->term[i]) && connected[i] )
            break;
      }

      /* trivial solution? */
      if( a == EAT_LAST )
      {
         printf("trivial solution in pruning \n");
         for( a = g->outbeg[g->source]; a != EAT_LAST; a = g->oeat[a] )
         {
            i = g->head[a];
            if( Is_term(g->term[i]) )
            {
               assert(connected[i]);
               result[a] = CONNECT;
            }
         }
         return SCIP_OKAY;
      }

      assert(g->mark[i]);
      root = i;
   }
   assert(root >= 0);
   assert(root < nnodes);

   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nnodes) );
   graph_path_exec(scip, g, MST_MODE, root, cost, mst);

   for( i = 0; i < nnodes; i++ )
   {
      if( g->mark[i] && (mst[i].edge != -1) )
      {
         assert(g->path_state[i] == CONNECT);
         assert(g->head[mst[i].edge] == i);
         assert(result[mst[i].edge] == -1);
         result[mst[i].edge] = CONNECT;
      }
   }

   /* connect all terminals */
   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) && i != g->source )
      {
         if( rmw )
            if( g->mark[i] )
               continue;

         e1 = g->inpbeg[i];
         assert(e1 >= 0);
         e2 = g->ieat[e1];

         if( e2 == EAT_LAST )
         {
            result[e1] = CONNECT;
         }
         else
         {
            assert(e2 >= 0);

            assert(g->ieat[e2] == EAT_LAST);
            k1 = g->tail[e1];
            k2 = g->tail[e2];
            assert(k1 == g->source || k2 == g->source);
            if( k1 != g->source && g->path_state[k1] == CONNECT )
            {
               result[e1] = CONNECT;
            }
            else if( k2 != g->source && g->path_state[k2] == CONNECT )
            {
               result[e2] = CONNECT;
            }
            else if( k1 == g->source )
            {
               result[e1] = CONNECT;
            }
            else if( k2 == g->source )
            {
               result[e2] = CONNECT;
            }
         }
      }
      else if( i == root && !rpc && !rmw )
      {
         for( e1 = g->inpbeg[i]; e1 != EAT_LAST; e1 = g->ieat[e1] )
            if( g->tail[e1] == g->source )
               break;
         assert(e1 != EAT_LAST);
         result[e1] = CONNECT;
      }
   }

   /* prune */
   do
   {
      count = 0;

      for( i = nnodes - 1; i >= 0; --i )
      {
         if( !g->mark[i] || g->path_state[i] != CONNECT || Is_term(g->term[i]) )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
            if( result[j] == CONNECT )
               break;

         if( j == EAT_LAST )
         {
            /* there has to be exactly one incoming edge
             */
            for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            {
               if( result[j] == CONNECT )
               {
                  result[j]    = -1;
                  g->mark[i]   = FALSE;
                  connected[i] = FALSE;
                  count++;
                  break;
               }
            }
            assert(j != EAT_LAST);
         }
      }
   }
   while( count > 0 );

#ifndef NDEBUG
   /* make sure there is no unconnected vertex */
   for( i = 0; i < nnodes; i++ )
   {
      if( connected[i] && i != g->source )
      {
         for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            if( result[j] == CONNECT )
               break;

         assert(j != EAT_LAST);

      }
   }
#endif

   assert(graph_sol_valid(scip, g, result));
   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}



/** build (rooted) prize collecting Steiner tree in such a way that all leaves are terminals */
SCIP_RETCODE SCIPStpHeurTMBuildTreePcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   PATH*                 mst,                /**< path data structure array */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            objresult,          /**< pointer to store objective value of result */
   int*                  connected           /**< CONNECT/UNKNOWN */
   )
{
   SCIP_Real obj;
   int i;
   int j;
   int e1;
   int e2;
   int k1;
   int k2;
   int root;
   int count;
   int nnodes;
   int orgroot;

   assert(g != NULL);
   assert(mst != NULL);
   assert(scip != NULL);
   assert(cost != NULL);
   assert(connected != NULL);

   obj = 0.0;
   nnodes = g->knots;
   orgroot = g->source;

   /* compute the MST, exclude all terminals */
   for( i = nnodes - 1; i >= 0; --i )
   {
      if( connected[i] == CONNECT && !Is_term(g->term[i]) )
         g->mark[i] = TRUE;
      else
         g->mark[i] = FALSE;
   }

   if( g->stp_type == STP_RPCSPG )
   {
      root = orgroot;
      g->mark[root] = TRUE;
   }
   else
   {
      int a;
      for( a = g->outbeg[orgroot]; a != EAT_LAST; a = g->oeat[a] )
      {
         i = g->head[a];
         if( !Is_term(g->term[i]) && connected[i] == CONNECT )
            break;
      }

      /* trivial solution? */
      if( a == EAT_LAST )
      {
         for( i = 0; i < nnodes; i++ )
            mst[i].edge = UNKNOWN;

         printf("trivial solution in buildPcMwTree \n");
         for( a = g->outbeg[orgroot]; a != EAT_LAST; a = g->oeat[a] )
         {
            i = g->head[a];
            if( Is_term(g->term[i]) )
            {
               obj += cost[a];
               mst[i].edge = a;
            }
         }
         (*objresult) = obj;
         return SCIP_OKAY;
      }

      assert(g->mark[i]);
      root = i;
   }
   assert(root >= 0);
   assert(root < nnodes);

   graph_path_exec(scip, g, MST_MODE, root, cost, mst);

   /* connect all terminals */
   for( i = nnodes - 1; i >= 0; --i )
   {
      if( Is_term(g->term[i]) && i != orgroot )
      {
         e1 = g->inpbeg[i];
         assert(e1 >= 0);
         e2 = g->ieat[e1];

         if( e2 == EAT_LAST )
         {
            mst[i].edge = e1;
         }
         else
         {
            assert(e2 >= 0);

            assert(g->ieat[e2] == EAT_LAST);
            k1 = g->tail[e1];
            k2 = g->tail[e2];
            assert(k1 == orgroot || k2 == orgroot);

            if( k1 != orgroot && g->path_state[k1] == CONNECT )
            {
               mst[i].edge = e1;
            }
            else if( k2 != orgroot && g->path_state[k2] == CONNECT )
            {
               mst[i].edge = e2;
            }
            else if( k1 == orgroot )
            {
               mst[i].edge = e1;
            }
            else if( k2 == orgroot )
            {
               mst[i].edge = e2;
            }
         }
      }
      else if( i == root && g->stp_type != STP_RPCSPG )
      {
         for( e1 = g->inpbeg[i]; e1 != EAT_LAST; e1 = g->ieat[e1] )
            if( g->tail[e1] == orgroot )
               break;
         assert(e1 != EAT_LAST);
         mst[i].edge = e1;
      }
   }

   /* prune */
   do
   {
      count = 0;

      for( i = nnodes - 1; i >= 0; --i )
      {
         if( !g->mark[i] || g->path_state[i] != CONNECT || Is_term(g->term[i]) )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
         {
            e1 = mst[g->head[j]].edge;
            if( e1 == j )
               break;
         }

         if( j == EAT_LAST )
         {
            mst[i].edge = UNKNOWN;
            g->mark[i] = FALSE;
            connected[i] = UNKNOWN;
            count++;
            break;
         }
      }
   }
   while( count > 0 );

   for( i = nnodes - 1; i >= 0; --i )
      if( mst[i].edge >= 0 )
         obj += cost[mst[i].edge];

   (*objresult) = obj;

   return SCIP_OKAY;
}


/** prune a Steiner tree in such a way that all leaves are terminals */
SCIP_RETCODE SCIPStpHeurTMPrune(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   int                   layer,              /**< layer (usually 0) */
   int*                  result,             /**< ST edges, which need to be set to UNKNOWN */
   STP_Bool*             connected           /**< ST nodes */
   )
{
   PATH*  mst;
   int i;
   int j;
   int count;
   int nnodes;

   assert(g != NULL);
   assert(scip != NULL);
   assert(cost != NULL);
   assert(layer == 0);
   assert(result != NULL);
   assert(connected != NULL);

   j = 0;
   nnodes = g->knots;

   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nnodes) );

   /* compute the MST */
   for( i = nnodes - 1; i >= 0; --i )
   {
      if( connected[i] )
         j++;
      g->mark[i] = connected[i];
   }

   assert(g->source >= 0);
   assert(g->source < nnodes);
   assert(j >= g->terms);

   graph_path_exec(scip, g, MST_MODE, g->source, cost, mst);

   for( i = nnodes - 1; i >= 0; --i )
   {
      if( connected[i] && (mst[i].edge != -1) )
      {
         assert(g->head[mst[i].edge] == i);
         assert(result[mst[i].edge] == UNKNOWN);

         result[mst[i].edge] = layer;
      }
   }

   /* prune */
   do
   {
      SCIPdebug(fputc('C', stdout));
      SCIPdebug(fflush(stdout));

      count = 0;

      for( i = nnodes - 1; i >= 0; --i )
      {
         if( !g->mark[i] )
            continue;

         if( g->term[i] == layer )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
            if( result[j] == layer )
               break;

         if( j == EAT_LAST )
         {
            /* there has to be exactly one incoming edge
             */
            for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            {
               if( result[j] == layer )
               {
                  result[j]    = -1;
                  g->mark[i]   = FALSE;
                  connected[i] = FALSE;
                  count++;
                  break;
               }
            }
         }
      }
   }
   while( count > 0 );

   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}

/** build Steiner tree in such a way that all leaves are terminals */
SCIP_RETCODE SCIPStpHeurTMBuildTree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   PATH*                 mst,                /**< path data structure array */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            objresult,          /**< pointer to store objective value of result */
   int*                  connected           /**< CONNECT/UNKNOWN */
   )
{
   SCIP_Real obj;
   int i;
   int j;
   int e1;
   int root;
   int count;
   int nnodes;

   assert(g != NULL);
   assert(mst != NULL);
   assert(scip != NULL);
   assert(cost != NULL);
   assert(connected != NULL);

   obj = 0.0;
   root = g->source;
   nnodes = g->knots;

   /* compute the MST */
   for( i = nnodes - 1; i >= 0; --i )
      g->mark[i] = (connected[i] == CONNECT);

   graph_path_exec(scip, g, MST_MODE, root, cost, mst);

   /* prune */
   do
   {
      count = 0;

      for( i = nnodes - 1; i >= 0; --i )
      {
         if( !g->mark[i] || Is_term(g->term[i]) )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
         {
            e1 = mst[g->head[j]].edge;
            if( e1 == j )
               break;
         }

         if( j == EAT_LAST )
         {
            mst[i].edge = UNKNOWN;
            g->mark[i] = FALSE;
            connected[i] = UNKNOWN;
            count++;
            break;
         }
      }
   }
   while( count > 0 );

   for( i = nnodes - 1; i >= 0; --i )
   {
      if( mst[i].edge >= 0 )
         obj += cost[mst[i].edge];
      else if( Is_term(g->term[i]) && i != root )
      {
         obj = FARAWAY;
         break;
      }
   }

   *objresult = obj;

   return SCIP_OKAY;
}

/** prune a degree constrained Steiner tree in such a way that all leaves are terminals */
SCIP_RETCODE SCIPStpHeurTMBuildTreeDc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   int*                  result,             /**< ST edges */
   STP_Bool*             connected           /**< ST nodes (to be set) */
   )
{
   int* queue;
   int i;
   int j;
   int qsize;
   int count;
   int nnodes;
   assert(scip != NULL);
   assert(g != NULL);
   assert(result != NULL);
   assert(connected != NULL);

   nnodes = g->knots;

   /*
    * DFS until all terminals are reached
    */

   for( i = 0; i < nnodes; i++ )
      connected[i] = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &queue, nnodes) );

   qsize = 0;
   queue[qsize++] = g->source;
   connected[g->source] = TRUE;

   while( qsize )
   {
      const int node = queue[--qsize];
      for( int e = g->outbeg[node]; e != EAT_LAST; e = g->oeat[e] )
      {
         if( (result[e] == CONNECT || result[flipedge(e)] == CONNECT) && !(connected[g->head[e]]) )
         {
            i = g->head[e];
            result[e] = CONNECT;
            result[flipedge(e)] = UNKNOWN;
            connected[i] = TRUE;
            queue[qsize++] = i;
         }
      }
   }

   SCIPfreeBufferArray(scip, &queue);

   /* prune */
   do
   {
      SCIPdebug(fputc('C', stdout));
      SCIPdebug(fflush(stdout));

      count = 0;

      for( i = 0; i < nnodes; i++ )
      {
         if( !connected[i] )
            continue;

         if( Is_term(g->term[i]) )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
            if( result[j] == CONNECT )
               break;

         if( j == EAT_LAST )
         {
            /* there has to be exactly one incoming edge
             */
            for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            {
               if( result[j] == CONNECT )
               {
                  result[j]    = UNKNOWN;
                  connected[i] = FALSE;
                  count++;
                  break;
               }
            }
            assert(j != EAT_LAST);
         }
      }
   }
   while( count > 0 );

   return SCIP_OKAY;
}

/*
 *  local functions
 */

static
SCIP_RETCODE prune(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs for DHCSTP */
   int*                  result,             /**< ST edges */
   STP_Bool*             connected           /**< ST nodes */
   )
{
   const int nedges = g->edges;

   if( g->stp_type != STP_DHCSTP )
      for( int e = 0; e < nedges; e++ )
         result[e] = UNKNOWN;

   if( g->stp_type == STP_MWCSP || g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG || g->stp_type == STP_RMWCSP )
      SCIP_CALL( SCIPStpHeurTMPrunePc(scip, g, g->cost, result, connected) );
   else
      SCIP_CALL( SCIPStpHeurTMPrune(scip, g, (g->stp_type != STP_DHCSTP) ? g->cost : cost, 0, result, connected) );

   return SCIP_OKAY;
}

/** Dijkstra based shortest paths heuristic */
static
SCIP_RETCODE computeSteinerTreeDijk(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            dijkdist,           /**< distance array */
   int*                  result,             /**< solution array (on edges) */
   int*                  dijkedge,           /**< predecessor edge array */
   int                   start,              /**< start vertex*/
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   STP_Bool*             connected           /**< array marking all solution vertices*/
   )
{
   int k;
   int nnodes;

   assert(g != NULL);

   nnodes = g->knots;

   for( k = 0; k < nnodes; k++ )
      g->mark[k] = (g->grad[k] > 0);

   graph_path_st(scip, g, cost, dijkdist, dijkedge, start, randnumgen, connected);

   SCIP_CALL(prune(scip, g, cost, result, connected));

   return SCIP_OKAY;
}

/** Dijkstra based shortest paths heuristic for PCSTP and MWCSP */
static
SCIP_RETCODE computeSteinerTreeDijkPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            dijkdist,           /**< distance array */
   int*                  result,             /**< solution array (on edges) */
   int*                  dijkedge,           /**< predecessor edge array */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array marking all solution vertices*/
   )
{
   if( g->stp_type == STP_RMWCSP )
      graph_path_st_rmw(scip, g, cost, dijkdist, dijkedge, start, connected);
   else if( g->stp_type == STP_RPCSPG )
      graph_path_st_rpc(scip, g, cost, dijkdist, dijkedge, start, connected);
   else
      graph_path_st_pcmw(scip, g, cost, dijkdist, dijkedge, start, connected);

   SCIP_CALL(prune(scip, g, cost, result, connected));

   return SCIP_OKAY;
}

/** Dijkstra based shortest paths heuristic for PCSTP and MWCSP that computes tree spanning all positive
 * vertex weights and subsequently prunes this tree */
static
SCIP_RETCODE computeSteinerTreeDijkPcMwFull(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            dijkdist,           /**< distance array */
   int*                  result,             /**< solution array (on edges) */
   int*                  dijkedge,           /**< predecessor edge array */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array marking all solution vertices*/
   )
{
   graph_path_st_pcmw_full(scip, g, cost, dijkdist, dijkedge, start, connected);

   SCIP_CALL(prune(scip, g, cost, result, connected));

   return SCIP_OKAY;
}


/** shortest paths based heuristic */
static
SCIP_RETCODE computeSteinerTree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
   SCIP_Real**           pathdist,           /**< distance array */
   int                   start,              /**< start vertex */
   int*                  perm,               /**< permutation array (on nodes) */
   int*                  result,             /**< solution array (on edges) */
   int*                  cluster,            /**< array used to store current vertices of each subtree during construction */
   int**                 pathedge,           /**< predecessor edge array */
   STP_Bool*             connected,          /**< array marking all solution vertices */
   SCIP_RANDNUMGEN*      randnumgen          /**< random number generator */
   )
{
   SCIP_Real min;
   SCIP_Bool directed = (g->stp_type == STP_SAP || g->stp_type == STP_DHCSTP);
   int    k;
   int    e;
   int    i;
   int    j;
   int    l;
   int    z;
   int    old;
   int    root;
   int    csize;
   int    newval;
   int    nnodes;

   assert(scip      != NULL);
   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(cluster != NULL);
   assert(perm != NULL);

   root = g->source;
   csize = 0;
   nnodes = g->knots;

   assert(0 <= start && start < nnodes);

   SCIPdebugMessage("Heuristic: Start=%5d \n", start);

   cluster[csize++] = start;

   /* initialize arrays */
   for( i = 0; i < nnodes; i++ )
   {
      g->mark[i]   = (g->grad[i] > 0);
      connected[i] = FALSE;
      perm[i] = i;
   }

   connected[start] = TRUE;

   SCIPrandomPermuteIntArray(randnumgen, perm, 0, nnodes - 1);

   assert(graph_valid(g));

   /* CONSTCOND */
   for( ;; )
   {
      /* Find a terminal with minimal distance to current ST */
      min = FARAWAY;
      old = -1;
      newval = -1;

      /* time limit exceeded? */
      if( SCIPisStopped(scip) )
         return SCIP_OKAY;

      /* find shortest path from current tree to unconnected terminal */
      for( l = nnodes - 1; l >= 0; --l )
      {
         i = perm[l];
         if( !Is_term(g->term[i]) || connected[i] || !g->mark[i] || (directed && !connected[root] && i != root) )
            continue;

         z = SCIPrandomGetInt(randnumgen, 0, nnodes - 1);

         for( k = csize - 1; k >= 0; k-- )
         {
            j = cluster[(k + z) % csize];
            assert(i != j);
            assert(connected[j]);

            if( SCIPisLT(scip, pathdist[i][j], min) )
            {
               min = pathdist[i][j];
               newval = i;
               old = j;
            }
         }
      }

      /* no new path found? */
      if( newval == -1 )
         break;

      /* mark the new path */
      assert(old > -1);
      assert(newval > -1);
      assert(pathdist[newval] != NULL);
      assert(pathdist[newval][old] < FARAWAY);
      assert(g->term[newval] == 0);
      assert(!connected[newval]);
      assert(connected[old]);

      SCIPdebug(fputc('R', stdout));
      SCIPdebug(fflush(stdout));

      /* start from current tree */
      k = old;

      while( k != newval )
      {
         e = pathedge[newval][k];
         k = g->tail[e];
         if (!connected[k])
         {
            connected[k] = TRUE;
            cluster[csize++] = k;
         }
      }
   }

   /* prune the tree */
   SCIP_CALL( prune(scip, g, cost, result, connected) );

   return SCIP_OKAY;
}


/** heuristic for degree constrained STPs */
static
SCIP_RETCODE computeDegConsTree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
   SCIP_Real**           pathdist,           /**< distances from each terminal to all nodes */
   int                   start,              /**< start vertex */
   int*                  perm,               /**< permutation array */
   int*                  result,             /**< array to indicate whether an edge is in the solution */
   int*                  cluster,            /**< array for internal computations */
   int**                 pathedge,           /**< ancestor edges for shortest paths from each terminal to all nodes */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   STP_Bool*             connected,          /**< array to indicate whether a vertex is in the solution */
   STP_Bool*             solfound            /**< pointer to store whether a solution has been found */
   )
{
   SCIP_Real min;
   int    csize = 0;
   int    k;
   int    e;
   int    l;
   int    i;
   int    j;
   int    t;
   int    u;
   int    z;
   int    n;
   int    old;
   int    tldegcount;
   int    degcount;
   int    mindegsum;
   int    degmax;
   int    newval;
   int    nnodes;
   int    nterms;
   int    termcount;
   int*   degs;
   int*   maxdegs;
   assert(scip      != NULL);
   assert(g         != NULL);
   assert(g->maxdeg != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(cluster   != NULL);
   assert(perm   != NULL);

   z = 0;
   nterms = g->terms;
   nnodes = g->knots;
   mindegsum = 2;
   maxdegs = g->maxdeg;

   SCIPdebugMessage("Heuristic: Start=%5d ", start);

   SCIP_CALL( SCIPallocBufferArray(scip, &degs, nnodes) );

   cluster[csize++] = start;

   for( i = 0; i < nnodes; i++ )
   {
      degs[i] = 0;
      g->mark[i] = (g->grad[i] > 0);
      connected[i] = FALSE;
      perm[i] = i;
   }

   for( e = 0; e < g->edges; e++ )
   {
      assert(SCIPisGT(scip, cost[e], 0.0));
      assert(result[e] == UNKNOWN);
   }
   connected[start] = TRUE;
   tldegcount = MIN(g->grad[start], maxdegs[start]);
   SCIPrandomPermuteIntArray(randnumgen, perm, 0, nnodes - 1);

   if( Is_term(g->term[start]) )
      termcount = 1;
   else
      termcount = 0;

   for( n = 0; n < nnodes; n++ )
   {
      /* Find a terminal with minimal distance to the current ST */
      min = FARAWAY - 1;
      /* is the free degree sum at most one and are we not connecting the last terminal? */
      if( tldegcount <= 1 && termcount < nterms - 1)
         degmax = 1;
      else
         degmax = 0;
      old = -1;
      newval = -1;
      for( t = 0; t < nnodes; t++ )
      {
         i = perm[t];
         if( !Is_term(g->term[i]) || connected[i] || g->grad[i] == 0 )
            continue;

         z = SCIPrandomGetInt(randnumgen, 0, nnodes - 1);

         for( k = 0; k < csize; k++ )
         {
            j = cluster[(k + z) % csize];
            assert(i != j);
            assert(connected[j]);

            if( SCIPisLE(scip, pathdist[i][j], min) && degs[j] < maxdegs[j])
            {
               u = j;
               degcount = 1;
               while( u != i )
               {
                  u = g->tail[pathedge[i][u]];
                  if( !connected[u] )
                  {
                     if( (MIN(g->grad[u], maxdegs[u]) < 2 || Is_term(g->term[u])) && u != i )
                     {
                        degcount = -2;
                        break;
                     }
                     degcount += MIN(g->grad[u] - 2, maxdegs[u] - 2);
                  }
                  else
                  {
                     assert(u != i);
                     l = g->tail[pathedge[i][u]];
                     if( !connected[l] && degs[u] >= maxdegs[u] )
                     {
                        degcount = -2;
                        break;
                     }
                  }
               }
               if( degcount >= degmax || (degcount >= mindegsum && SCIPisLT(scip, pathdist[i][j], min)) )
               {
                  degmax = degcount;
                  min = pathdist[i][j];
                  newval = i;
                  old = j;
               }
            }
         }
      }

      if( newval == -1 )
      {
         j = UNKNOWN;
         for( k = 0; k < csize; k++ )
         {
            j = cluster[(k + z) % csize];
            if( degs[j] < maxdegs[j] )
               break;
         }
         if( j != UNKNOWN )
         {
            assert(k != csize);

            min = FARAWAY + 1;
            newval = UNKNOWN;
            for( e = g->outbeg[j]; e != EAT_LAST; e = g->oeat[e] )
            {
               u = g->head[e];
               if( !Is_term(g->term[u]) && !connected[u] && SCIPisGE(scip, min, cost[e]) )
               {
                  min = cost[e];
                  k = e;
                  newval = u;
               }
            }
            if( newval != UNKNOWN )
            {
               result[flipedge(k)] = CONNECT;
               degs[newval]++;
               degs[j]++;
               connected[newval] = TRUE;
               cluster[csize++] = newval;
               tldegcount += MIN(maxdegs[newval], g->grad[newval]) - 2;
               continue;
            }

         }
      }
      tldegcount += degmax - 1;
      /* break? */
      if( newval == -1 )
         break;
      /* Weg setzten
       */
      assert(old > -1);
      assert(newval > -1);
      assert(pathdist[newval] != NULL);
      assert(g->term[newval] == 0);
      assert(!connected[newval]);
      assert(connected[old]);

      SCIPdebug(fputc('R', stdout));
      SCIPdebug(fflush(stdout));

      /* mark new tree nodes/edges */
      k = old;
      if( Is_term(g->term[newval]) )
         termcount++;

      while(k != newval)
      {
         e = pathedge[newval][k];
         u = k;
         k = g->tail[e];

         if( !connected[k])
         {
            result[flipedge(e)] = CONNECT;
            degs[u]++;
            degs[k]++;
            connected[k] = TRUE;
            cluster[csize++] = k;
            if( k != newval )
               assert(!Is_term(g->term[k]));
         }
      }
      if( termcount == nterms )
         break;
      assert(degs[newval] == 1);
   }

   *solfound = TRUE;

   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) && !connected[i] )
      {
         *solfound = FALSE;
         break;
      }
   }

   if( *solfound )
   {
      /* prune the solution */
      SCIP_CALL( SCIPStpHeurTMBuildTreeDc(scip, g, result, connected) );

      for( t = 0; t < nnodes; t++ )
         if( degs[t] > maxdegs[t] )
            *solfound = FALSE;
   }

   SCIPfreeBufferArray(scip, &degs);
   return SCIP_OKAY;
}

/** Voronoi based shortest path heuristic */
static
SCIP_RETCODE computeSteinerTreeVnoi(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   SCIP_PQUEUE*          pqueue,             /**< priority queue */
   GNODE**               gnodearr,           /**< internal array */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
   SCIP_Real**           node_dist,          /**< internal array */
   int                   start,              /**< start vertex */
   int*                  result,             /**< array to indicate whether an edge is in the solution */
   int*                  vcount,             /**< internal array */
   int*                  nodenterms,         /**< internal array */
   int**                 node_base,          /**< internal array */
   int**                 node_edge,          /**< internal array */
   STP_Bool              firstrun,           /**< method called for the first time? (during one heuristic round) */
   STP_Bool*             connected           /**< array to indicate whether a vertex is in the solution */
   )
{
   int    k;
   int    i;
   int    j;
   int    best;
   int    term;
   int    count;
   int   nnodes;
   int   nterms;

   assert(scip      != NULL);
   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   nnodes = g->knots;
   nterms = g->terms;

   SCIPdebugMessage("TM_Polzin Heuristic: Start=%5d ", start);

   /* if the heuristic is called for the first time several data structures have to be set up */
   if( firstrun )
   {
      PATH* vnoi;
      SCIP_Real* vcost;
      int old;
      int oedge;
      int root = g->source;
      int   ntovisit;
      int   nneighbnodes;
      int   nneighbterms;
      int   nreachednodes;
      int*  state;
      int*  vbase;
      int*  terms;
      int*  tovisit;
      int*  reachednodes;
      STP_Bool* termsmark;
      STP_Bool* visited;
      int e;
      /* PHASE I: */
      for( i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);

      /* allocate memory needed in PHASE I */
      SCIP_CALL( SCIPallocBufferArray(scip, &terms, nterms) );
      SCIP_CALL( SCIPallocBufferArray(scip, &termsmark, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &visited, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &reachednodes, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tovisit, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vcost, nnodes) );

      j = 0;
      for( i = 0; i < nnodes; i++ )
      {
         visited[i] = FALSE;
         if( Is_term(g->term[i]) )
         {
            termsmark[i] = TRUE;
            terms[j++] = i;
         }
         else
         {
            termsmark[i] = FALSE;
         }
      }

      for( e = 0; e < g->edges; e++)
      {
         assert(SCIPisGE(scip, cost[e], 0.0));
         assert(SCIPisGE(scip, costrev[e], 0.0));
      }

      assert(j == nterms);
      graph_voronoi(scip, g, cost, costrev, termsmark, vbase, vnoi);
      state = g->path_state;

      for( i = 0; i < nnodes; i++ )
         if( Is_term(g->term[i]) )
            assert(vbase[i] == i);

      for( k = 0; k < nnodes; k++ )
      {
         connected[k] = FALSE;
         vcount[k] = 0;
         gnodearr[k]->number = k;
         if( !Is_term(g->term[k]) )
         {
            node_dist[k][0] = vnoi[k].dist;
            node_edge[k][0] = vnoi[k].edge;
            node_base[k][0] = vbase[k];
            nodenterms[k] = 1;
         }
         else
         {
            nodenterms[k] = 0;
            node_edge[k][0] = UNKNOWN;
            termsmark[k] = FALSE;
         }
         state[k] = UNKNOWN;
         vcost[k] = vnoi[k].dist;
         vnoi[k].dist = FARAWAY;
      }

      /* for each terminal: extend the voronoi regions until all neighbouring terminals have been visited */
      for( i = 0; i < nterms; i++ )
      {
         term = terms[i];
         nneighbterms = 0;
         nneighbnodes = 0;
         nreachednodes = 0;
         for( k = 0; k < nnodes; k++ )
            assert(termsmark[k] == FALSE);
         /* DFS (starting from terminal i) until the entire voronoi region has been visited */
         tovisit[0] = term;
         ntovisit = 1;
         visited[term] = TRUE;
         state[term] = CONNECT;
         while( ntovisit > 0 )
         {
            /* iterate all incident edges */
            old = tovisit[--ntovisit];

            for( oedge = g->outbeg[old]; oedge != EAT_LAST; oedge = g->oeat[oedge] )
            {
               k = g->head[oedge];

               /* is node k in the voronoi region of the i-th terminal ? */
               if( vbase[k] == term )
               {
                  if( !visited[k] )
                  {
                     state[k] = CONNECT;
                     assert(nnodes - (nneighbnodes + 1) > ntovisit);
                     tovisit[ntovisit++] = k;
                     visited[k] = TRUE;
                     reachednodes[nreachednodes++] = k;
                  }
               }
               else
               {
                  if( !visited[k] )
                  {
                     visited[k] = TRUE;
                     vnoi[k].dist = vcost[old] + ((vbase[k] == root)? cost[oedge] : costrev[oedge]);
                     vnoi[k].edge = oedge;

                     if( termsmark[vbase[k]] == FALSE )
                     {
                        termsmark[vbase[k]] = TRUE;
                        nneighbterms++;
                     }
                     assert(nnodes - (nneighbnodes + 1) > ntovisit - 1);
                     tovisit[nnodes - (++nneighbnodes)] = k;
                  }
                  else
                  {
                     /* if edge 'oedge' allows a shorter connection of node k, update */
                     if( SCIPisGT(scip, vnoi[k].dist, vcost[old] + ((vbase[k] == root)? cost[oedge] : costrev[oedge])) )
                     {
                        vnoi[k].dist = vcost[old] + ((vbase[k] == root)? cost[oedge] : costrev[oedge]);
                        vnoi[k].edge = oedge;
                     }
                  }
               }
            }
         }

         count = 0;
         for( j = 0; j < nneighbnodes; j++ )
         {
            assert(termsmark[vbase[tovisit[nnodes - j - 1]]]);
            heap_add(g->path_heap, state, &count, tovisit[nnodes - j - 1], vnoi);
         }
         SCIP_CALL( graph_voronoiExtend(scip, g, ((term == root)? cost : costrev), vnoi, node_dist, node_base, node_edge, termsmark, reachednodes, &nreachednodes, nodenterms,
               nneighbterms, term, nneighbnodes) );

         reachednodes[nreachednodes++] = term;

         for( j = 0; j < nreachednodes; j++ )
         {
            vnoi[reachednodes[j]].dist = FARAWAY;
            state[reachednodes[j]] = UNKNOWN;
            visited[reachednodes[j]] = FALSE;
         }

         for( j = 0; j < nneighbnodes; j++ )
         {
            vnoi[tovisit[nnodes - j - 1]].dist = FARAWAY;
            state[tovisit[nnodes - j - 1]] = UNKNOWN;
            visited[tovisit[nnodes - j - 1]] = FALSE;
         }
      }

      /* for each node v: sort the terminal arrays according to their distance to v */
      for( i = 0; i < nnodes && !SCIPisStopped(scip); i++ )
         SCIPsortRealIntInt(node_dist[i], node_base[i], node_edge[i], nodenterms[i]);

      /* free memory */
      SCIPfreeBufferArray(scip, &vcost);
      SCIPfreeBufferArray(scip, &tovisit);
      SCIPfreeBufferArray(scip, &vbase);
      SCIPfreeBufferArray(scip, &reachednodes);
      SCIPfreeBufferArray(scip, &visited);
      SCIPfreeBufferArray(scip, &vnoi);
      SCIPfreeBufferArray(scip, &termsmark);
      SCIPfreeBufferArray(scip, &terms);
   }

   /* PHASE II */
   else
   {
      for( k = 0; k < nnodes; k++ )
      {
         connected[k] = FALSE;
         vcount[k] = 0;
      }
   }

   connected[start] = TRUE;
   gnodearr[start]->dist = node_dist[start][0];
   SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[start]) );

   while( SCIPpqueueNElems(pqueue) > 0 && !SCIPisStopped(scip) )
   {
      best = ((GNODE*) SCIPpqueueRemove(pqueue))->number;

      term = node_base[best][vcount[best]];
      assert( Is_term(g->term[term]) );
      /* has the terminal already been connected? */
      if( !connected[term] )
      {
         /* connect the terminal */
         k = g->tail[node_edge[best][vcount[best]]];
         while( k != term )
         {
            j = 0;

            while( node_base[k][vcount[k] + j] != term )
               j++;

            assert(vcount[k] + j < nodenterms[k]);

            if( !connected[k] )
            {
               assert(vcount[k] == 0);

               connected[k] = TRUE;
               while( vcount[k] < nodenterms[k] && connected[node_base[k][vcount[k]]] )
               {
                  vcount[k]++;
                  j--;
               }

               if( vcount[k] < nodenterms[k] )
               {
                  gnodearr[k]->dist = node_dist[k][vcount[k]];
                  SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[k]) );
               }
            }

            assert( vcount[k] + j < nodenterms[k] );
            assert(node_base[k][vcount[k] + j] == term);
            k = g->tail[node_edge[k][vcount[k] + j]];
         }

         /* finally, connected the terminal */
         assert( k == term );
         assert( !connected[k] );
         connected[k] = TRUE;

         assert( vcount[k] == 0 );
         while( vcount[k] < nodenterms[k] && connected[node_base[k][vcount[k]]] )
         {
            vcount[k]++;
         }
         if( vcount[k] < nodenterms[k] )
         {
            gnodearr[k]->dist = node_dist[k][vcount[k]];
            SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[k]) );
         }
      }

      while( vcount[best] + 1 < nodenterms[best] )
      {
         if( !connected[node_base[best][++vcount[best]]] )
         {
            gnodearr[best]->dist = node_dist[best][vcount[best]];
            SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[best]) );
            break;
         }
      }
   }

   /* prune the ST, so that all leaves are terminals */
   SCIP_CALL( prune(scip, g, cost, result, connected) );

   return SCIP_OKAY;
}

/** execute shortest paths heuristic to obtain a Steiner tree */
SCIP_RETCODE SCIPStpHeurTMRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   int*                  starts,             /**< array containing start vertices (NULL to not provide any) */
   int*                  bestnewstart,       /**< pointer to the start vertex resulting in the best solution */
   int*                  best_result,        /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   int                   runs,               /**< number of runs */
   int                   bestincstart,       /**< best incumbent start vertex */
   SCIP_Real*            cost,               /**< arc costs */
   SCIP_Real*            costrev,            /**< reversed arc costs */
   SCIP_Real*            hopfactor,          /**< edge cost multiplicator for HC problems */
   SCIP_Real*            nodepriority,       /**< vertex priorities for vertices to be starting points (NULL for no priorities) */
   SCIP_Real             maxcost,            /**< maximal edge cost (only for HC) */
   SCIP_Bool*            success,            /**< pointer to store whether a solution could be found */
   SCIP_Bool             pcmwfull            /**< use full computation of tree (i.e. connect all terminals and prune), only for prize-collecting variants */
   )
{
   SCIP_PQUEUE* pqueue = NULL;
   SCIP_Longint nexecs;
   SCIP_Real obj;
   SCIP_Real min = FARAWAY;
   SCIP_Real* dijkdist = NULL;
   SCIP_Real** pathdist = NULL;
   SCIP_Real** node_dist = NULL;
   SCIP_Bool lsuccess;
   GNODE** gnodearr = NULL;
   int best;
   int k;
   int r;
   int e;
   int root;
   int mode;
   int nedges;
   int nnodes;
   int nterms;
   int* perm = NULL;
   int* start;
   int* vcount = NULL;
   int* RESTRICT result;
   int* cluster = NULL;
   int* dijkedge = NULL;
   int* nodenterms = NULL;
   int** pathedge = NULL;
   int** node_base = NULL;
   int** node_edge = NULL;
   SCIP_Bool startsgiven;
   STP_Bool* connected;

   assert(scip != NULL);
   assert(cost != NULL);
   assert(graph != NULL);
   assert(costrev != NULL);
   assert(best_result != NULL);

   if( heurdata == NULL )
   {
      assert(SCIPfindHeur(scip, "TM") != NULL);
      heurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));
   }

   best = bestincstart;
   root = graph->source;
   nnodes = graph->knots;
   nedges = graph->edges;
   nterms = graph->terms;
   nexecs = heurdata->nexecs;
   (*success) = FALSE;

   if( runs < 1 )
   {
      startsgiven = FALSE;
      runs = 1;
   }
   else
      startsgiven = (starts != NULL);

   for( e = 0; e < nedges; e++)
   {
      assert(SCIPisGE(scip, cost[e], 0.0));
      assert(SCIPisGE(scip, costrev[e], 0.0));
   }

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   assert(root >= 0);
   assert(nedges > 0);
   assert(nnodes > 0);
   assert(nterms > 0);
   assert(nexecs >= 0);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &connected, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &start, MIN(runs, nnodes)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );

   /* get user parameter */
   mode = heurdata->type;
   assert(mode == AUTO || mode == TM_SP || mode == TM_VORONOI || mode == TM_DIJKSTRA);

   if( graph->stp_type == STP_RPCSPG || graph->stp_type == STP_PCSPG || graph->stp_type == STP_MWCSP || graph->stp_type == STP_RMWCSP )
   {
      mode = TM_DIJKSTRA;
      graph_pc_2transcheck(graph);
   }
   else if( graph->stp_type == STP_DHCSTP )
   {
      mode = TM_SP;
   }
   else if( graph->stp_type == STP_DCSTP )
   {
      mode = TM_SP;
   }
   else if( graph->stp_type == STP_SAP )
   {
      mode = TM_SP;
   }
   else
   {
      if( mode == AUTO )
         mode = TM_DIJKSTRA;
   }

   /* allocate memory according to selected mode */
   if( mode == TM_DIJKSTRA || graph->stp_type == STP_DHCSTP )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &dijkdist, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &dijkedge, nnodes) );
   }
   if( mode == TM_VORONOI )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &nodenterms, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &node_base, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &node_dist, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &node_edge, nnodes) );

      for( k = 0; k < nnodes; k++ )
      {
         SCIP_CALL( SCIPallocBuffer(scip, &gnodearr[k]) ); /*lint !e866*/
         SCIP_CALL( SCIPallocBufferArray(scip, &node_base[k], nterms) ); /*lint !e866*/
         SCIP_CALL( SCIPallocBufferArray(scip, &node_dist[k], nterms) ); /*lint !e866*/
         SCIP_CALL( SCIPallocBufferArray(scip, &node_edge[k], nterms) ); /*lint !e866*/
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &vcount, nnodes) );
      SCIP_CALL( SCIPpqueueCreate( &pqueue, nnodes, 2.0, GNODECmpByDist) );
   }
   if( mode == TM_SP )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &pathdist, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes) );
      BMSclearMemoryArray(pathdist, nnodes);
      BMSclearMemoryArray(pathedge, nnodes);

      for( k = 0; k < nnodes; k++ )
      {
         graph->mark[k] = (graph->grad[k] > 0);

         if( Is_term(graph->term[k]) )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(pathdist[k]), nnodes) ); /*lint !e866*/
            SCIP_CALL( SCIPallocBufferArray(scip, &(pathedge[k]), nnodes) ); /*lint !e866*/
         }
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &cluster, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, nnodes) );
   }

   /* compute number of iterations and starting points for SPH */
   if( graph->stp_type == STP_PCSPG || graph->stp_type == STP_MWCSP || graph->stp_type == STP_RMWCSP || graph->stp_type == STP_RPCSPG )
   {
      if( runs > (nterms - 1) && graph->stp_type != STP_RMWCSP )
         runs = nterms - 1;
      else if( runs > (nterms) && graph->stp_type == STP_RMWCSP )
         runs = nterms;
   }
   else if( startsgiven )
   {
      for( k = 0; k < MIN(runs, nnodes); k++ )
         start[k] = starts[k];
   }
   else if( runs < nnodes || graph->stp_type == STP_DHCSTP )
   {
      if( best < 0 )
      {
         best = root;
      }
      else
      {
         int randint = SCIPrandomGetInt(heurdata->randnumgen, 0, 2);
         if( randint == 0 )
            best = -1;
         else if( randint == 1 )
            best = root;
      }
      r = 0;
      if( graph->stp_type == STP_DHCSTP )
         graph_path_execX(scip, graph, root, cost, dijkdist, dijkedge);

      /* allocate memory for permutation array */
      if( perm == NULL )
         SCIP_CALL( SCIPallocBufferArray(scip, &perm, nnodes) );

      /* are there no nodes are to be priorized a priori? */
      if( nodepriority == NULL )
      {
         for( k = 0; k < nnodes; k++ )
            perm[k] = k;
         SCIPrandomPermuteIntArray(heurdata->randnumgen, perm, 0, nnodes);

         /* use terminals (randomly permutated) as starting points for TM heuristic */
         for( k = 0; k < nnodes; k++ )
         {
            if( r >= runs || r >= nterms )
               break;

            if( Is_term(graph->term[perm[k]]) )
               start[r++] = perm[k];
         }

         /* fill empty slots */
         for( k = 0; k < nnodes && r < runs; k++ )
         {
            if( graph->stp_type == STP_DHCSTP )
            {
               assert(dijkdist != NULL);
               if( SCIPisGE(scip, dijkdist[perm[k]], BLOCKED) )
                  continue;
            }
            if( !Is_term(graph->term[perm[k]]) && graph->mark[k] )
               start[r++] = perm[k];
         }
      }
      else
      {
         SCIP_Real max = 0.0;
         int bbound = runs - runs / 3;

         for( k = 0; k < nnodes; k++ )
         {
            perm[k] = k;
            if( SCIPisLT(scip, max, nodepriority[k]) && Is_term(graph->term[k]) )
               max = nodepriority[k];
         }
         for( k = 0; k < nnodes; k++ )
         {
            if( Is_term(graph->term[k]) )
            {
               nodepriority[k] += SCIPrandomGetReal(heurdata->randnumgen, 0.0, max);
            }
            else if( SCIPisLE(scip, 1.0, nodepriority[k]) )
            {
               nodepriority[k] = nodepriority[k] * SCIPrandomGetReal(heurdata->randnumgen, 1.0, 2.0);
            }
         }

         SCIPsortRealInt(nodepriority, perm, nnodes);

         for( k = nnodes - 1; k >= 0; k-- )
         {
            if( r >= nterms || r >= bbound )
               break;

            if( Is_term(graph->term[perm[k]]) )
            {
               start[r++] = perm[k];
               perm[k] = -1;
            }
         }

         /* fill empty slots */
         for( k = nnodes - 1; k >= 0 && r < runs; k-- )
         {
            if( perm[k] == -1 )
               continue;
            if( graph->stp_type == STP_DHCSTP )
            {
               assert(dijkdist != NULL);
               if( SCIPisGE(scip, dijkdist[perm[k]], BLOCKED) )
                  continue;
            }
            if( graph->mark[k] )
               start[r++] = perm[k];
         }
      }
      /* not all slots filled? */
      if( r < runs )
         runs = r;

      if( best >= 0 )
      {
         /* check whether we have a already selected the best starting node */
         for( r = 0; r < runs; r++ )
            if( start[r] == best )
               break;

         /* no, we still have to */
         if( r == runs && runs > 0 )
            start[nexecs % runs] = best;
      }
   } /* STP_DHCSTP */
   else
   {
      runs = nnodes;
      for( k = 0; k < nnodes; k++ )
         start[k] = k;
   }

   /* perform SPH computations, differentiate between STP variants */
   if( graph->stp_type == STP_PCSPG || graph->stp_type == STP_RPCSPG || graph->stp_type == STP_MWCSP || graph->stp_type == STP_RMWCSP )
   {
      SCIP_Real* terminalprio;
      int* terminalperm;
      int t;
      const SCIP_Bool rmw = (graph->stp_type == STP_RMWCSP);
      const SCIP_Bool rpc = (graph->stp_type == STP_RPCSPG);

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &terminalperm, nterms) );
      SCIP_CALL( SCIPallocBufferArray(scip, &terminalprio, nterms) );

      /* initialize arrays */

      t = 0;
      if( nodepriority == NULL )
      {
         for( k = nnodes - 1; k >= 0; --k )
            if( Is_pterm(graph->term[k]) && graph->grad[k] > 0 )
            {
               assert(SCIPisGT(scip, graph->prize[k], 0.0));
               terminalperm[t] = k;
               terminalprio[t++] = SCIPrandomGetReal(heurdata->randnumgen, 0.0, graph->prize[k]);
            }
      }
      else
      {
         for( k = nnodes - 1; k >= 0; --k )
             if( Is_pterm(graph->term[k]) && graph->grad[k] > 0 )
             {
                assert(SCIPisGT(scip, graph->prize[k], 0.0));
                terminalperm[t] = k;
                terminalprio[t++] = SCIPrandomGetReal(heurdata->randnumgen, nodepriority[k] / 2.0, nodepriority[k]);
             }
      }

      if( rmw )
      {
         SCIP_Real max = 0.0;
         int head;
         for( k = t - 1; k >= 0; --k )
         {
            if( SCIPisGT(scip, terminalprio[k], max) )
               max = terminalprio[k];
         }

         for( k = nnodes - 1; k >= 0; --k )
            graph->mark[k] = (graph->grad[k] > 0);

         for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
         {
            if( SCIPisGT(scip, graph->cost[e], 0.0) && Is_term(graph->term[graph->head[e]]) )
            {
               head = graph->head[e];
               graph->mark[head] = FALSE;
               assert(graph->grad[head] == 2);
            }
         }

         for( k = nnodes - 1; k >= 0; --k )
         {
            if( Is_term(graph->term[k]) && graph->mark[k] )
            {
               assert(SCIPisGT(scip, graph->prize[k], 0.0));
               terminalperm[t] = k;
               terminalprio[t++] = SCIPrandomGetReal(heurdata->randnumgen, max / 2.0, 1.5 * max);
            }
         }
         assert(nterms == t);
         SCIPsortRealInt(terminalprio, terminalperm, nterms);
      }
      else
      {
         if( rpc )
         {
            terminalperm[t] = root;
            terminalprio[t++] = FARAWAY;
            runs++;

            SCIPsortRealInt(terminalprio, terminalperm, nterms);
         }
         else
         {
            SCIPsortRealInt(terminalprio, terminalperm, nterms - 1);
         }
      }

      if( best >= 0 && best < nnodes && Is_pterm(graph->term[best]) && SCIPrandomGetInt(heurdata->randnumgen, 0, 2) == 1 )
         terminalperm[nterms - 1] = best;

      /* local main loop */
      for( r = runs - 1; r >= 0; --r )
      {
         k = terminalperm[r];

         if( pcmwfull )
         {
            SCIP_CALL( computeSteinerTreeDijkPcMwFull(scip, graph, cost, dijkdist, result, dijkedge, k, connected) );
         }
         else
         {
            SCIP_CALL( computeSteinerTreeDijkPcMw(scip, graph, cost, dijkdist, result, dijkedge, k, connected) );
         }

         if( SCIPisStopped(scip) )
            break;

         /* compute objective value (wrt original costs) */
         obj = 0.0;
         for( e = nedges - 1; e >= 0; e-- ) /* todo: save result array as char put into computeSteinerTreeDijkPcMw */
            if( result[e] == CONNECT )
               obj += graph->cost[e];

         if( SCIPisLT(scip, obj, min) )
         {
            if( bestnewstart != NULL )
               *bestnewstart = k;
            min = obj;

            SCIPdebugMessage("\n Obj(run: %d, ncall: %d)=%.12e\n\n", r, (int) nexecs, obj);

            for( e = 0; e < nedges; e++ )
               best_result[e] = result[e];
            (*success) = TRUE;
         }
      }

      /* free memory */
      SCIPfreeBufferArray(scip, &terminalprio);
      SCIPfreeBufferArray(scip, &terminalperm);
   } /* STP_PCSPG or STP_(R)MWCSP */
   else
   {
      SCIP_Real* orgcost = NULL;
      int edgecount;
      STP_Bool solfound = FALSE;

      if( SCIPisLE(scip, maxcost, 0.0) )
         maxcost = 1.0;

      if( graph->stp_type == STP_DHCSTP )
      {
         SCIP_Real bestfactor = -1;

         assert(hopfactor != NULL);
         assert(SCIPisGT(scip, (*hopfactor), 0.0));

         SCIP_CALL( SCIPallocBufferArray(scip, &orgcost, nedges) );

         BMScopyMemoryArray(orgcost, cost, nedges);

         /* do a warm-up run */
         for( r = 0; r < 10; r++ )
         {
            for( e = 0; e < nedges; e++ )
            {
               if( (SCIPisLT(scip, cost[e], BLOCKED )) )
                  cost[e] = 1.0 + orgcost[e] / ((*hopfactor) * maxcost);
               result[e] = UNKNOWN;
            }

            SCIP_CALL( computeSteinerTreeDijk(scip, graph, cost, dijkdist, result, dijkedge, root, heurdata->randnumgen, connected) );

            obj = 0.0;
            edgecount = 0;

            for( e = 0; e < nedges; e++)
            {
               if( result[e] == CONNECT )
               {
                  obj += graph->cost[e];
                  edgecount++;
               }
            }
            lsuccess = FALSE;
            if( SCIPisLT(scip, obj, min) && edgecount <= graph->hoplimit )
            {
               min = obj;
               if( bestnewstart != NULL )
                  *bestnewstart = root;

               for( e = 0; e < nedges; e++ )
                  best_result[e] = result[e];
               (*success) = TRUE;
               lsuccess = TRUE;
               bestfactor = (*hopfactor);
            }

            if( !lsuccess || SCIPisGT(scip, fabs((double) edgecount - graph->hoplimit) / (double) graph->hoplimit, 0.05) )
            {
               if( !lsuccess )
               {
                  if( (*success) )
                  {
                     (*hopfactor) = (*hopfactor) * (1.0 + fabs((double) edgecount - graph->hoplimit) / (double) graph->hoplimit);
                  }
                  else
                  {
                     (*hopfactor) = (*hopfactor) * (1.0 + 3 * fabs((double) edgecount - graph->hoplimit) / (double) graph->hoplimit);
                     bestfactor = (*hopfactor);
                  }
               }
               else
               {
                  (*hopfactor) = (*hopfactor) / (1.0 + fabs((double) edgecount - graph->hoplimit) / (double) graph->hoplimit);
               }

               assert(SCIPisGT(scip, (*hopfactor), 0.0));
            }
            else
            {
               break;
            }
         }
         (*hopfactor) = bestfactor;

         for( e = 0; e < nedges; e++ )
            if( (SCIPisLT(scip, cost[e], BLOCKED )) )
               cost[e] = 1.0 + orgcost[e] / ((*hopfactor) * maxcost);
         for( e = 0; e < nedges; e++)
            costrev[e] = cost[flipedge(e)];
      }
      for( r = 0; r < runs; r++ )
      {
         assert(start[r] >= 0);
         assert(start[r] < nnodes);

         if( mode == TM_DIJKSTRA && graph->stp_type != STP_DCSTP )
         {
            SCIP_CALL( computeSteinerTreeDijk(scip, graph, cost, dijkdist, result, dijkedge, start[r],  heurdata->randnumgen, connected) );
         }
         else if( graph->stp_type == STP_DCSTP )
         {
            /* first run? */
            if( r == 0 )
            {
               assert(pathdist != NULL);
               assert(pathedge != NULL);

               for( k = 0; k < nnodes; k++ )
                  graph->mark[k] = (graph->grad[k] > 0);

               /* initialize shortest paths from all terminals */
               for( k = 0; k < nnodes; k++ )
               {
                  if( Is_term(graph->term[k]) )
                  {
                     if( root == k )
                        graph_path_execX(scip, graph, k, cost,  pathdist[k], pathedge[k]);
                     else
                        graph_path_execX(scip, graph, k, costrev, pathdist[k], pathedge[k]);
                  }
               }
            }
            for( e = 0; e < nedges; e++ )
               result[e] = UNKNOWN;

            SCIP_CALL( computeDegConsTree(scip, graph, cost, costrev, pathdist, start[r], perm, result, cluster, pathedge,  heurdata->randnumgen, connected, &solfound) );
         }
         else if( mode == TM_SP )
         {
            for( e = nedges - 1; e >= 0; --e )
               result[e] = UNKNOWN;

            if( r == 0 )
            {
               int i;
               assert(pathdist != NULL);
               assert(pathedge != NULL);

               for( i = 0; i < nnodes; i++ )
                  graph->mark[i] = (graph->grad[i] > 0);

               /* initialize shortest paths from all terminals */
               for( k = 0; k < nnodes; k++ )
               {
                  if( Is_term(graph->term[k]) )
                  {
                     if( root == k )
                        graph_path_execX(scip, graph, k, cost,  pathdist[k], pathedge[k]);
                     else
                        graph_path_execX(scip, graph, k, costrev, pathdist[k], pathedge[k]);
                  }
               }
            }

            SCIP_CALL( computeSteinerTree(scip, graph, cost, costrev, pathdist, start[r], perm, result, cluster, pathedge, connected, heurdata->randnumgen) );
         }
         else
         {
            SCIP_CALL( computeSteinerTreeVnoi(scip, graph, pqueue, gnodearr, cost, costrev, node_dist, start[r], result, vcount,
                  nodenterms, node_base, node_edge, (r == 0), connected) );
         }
         obj = 0.0;
         edgecount = 0;

         /* here another measure than in the do_(...) heuristics is being used */
         for( e = 0; e < nedges; e++)
         {
            if( result[e] >= 0 )
            {
               obj += graph->cost[e];
               edgecount++;
            }
         }

         SCIPdebugMessage(" Obj=%.12e\n", obj);

         if( SCIPisLT(scip, obj, min) && (graph->stp_type != STP_DCSTP || solfound) && !SCIPisStopped(scip) && r < runs )
         {
            if( graph->stp_type != STP_DHCSTP || edgecount <= graph->hoplimit )
            {
               min = obj;
               if( bestnewstart != NULL )
                  *bestnewstart = start[r];

               for( e = 0; e < nedges; e++ )
                  best_result[e] = result[e];
               (*success) = TRUE;
            }
         }

         /* time limit exceeded?*/
         if( SCIPisStopped(scip) )
            break;
      }

      if( graph->stp_type == STP_DHCSTP )
      {
         assert(orgcost != NULL);
         for( e = 0; e < nedges; e++ )
         {
            cost[e] = orgcost[e];
            costrev[e] = orgcost[flipedge(e)];
         }
         SCIPfreeBufferArray(scip, &orgcost);
      }
   }

   /* free allocated memory */
   SCIPfreeBufferArrayNull(scip, &perm);
   if( mode == TM_SP )
   {
      assert(pathedge != NULL);
      assert(pathdist != NULL);
      SCIPfreeBufferArray(scip, &cluster);
      for( k = nnodes - 1; k >= 0; k-- )
      {
         SCIPfreeBufferArrayNull(scip, &(pathedge[k]));
         SCIPfreeBufferArrayNull(scip, &(pathdist[k]));
      }

      SCIPfreeBufferArray(scip, &pathedge);
      SCIPfreeBufferArray(scip, &pathdist);
   }
   else if( mode == TM_VORONOI )
   {
      SCIPpqueueFree(&pqueue);

      SCIPfreeBufferArray(scip, &vcount);
      assert(node_edge != NULL);
      assert(node_dist != NULL);
      assert(node_base != NULL);
      assert(gnodearr != NULL);
      for( k = nnodes - 1; k >= 0; k-- )
      {
         SCIPfreeBufferArray(scip, &node_edge[k]);
         SCIPfreeBufferArray(scip, &node_dist[k]);
         SCIPfreeBufferArray(scip, &node_base[k]);
         SCIPfreeBuffer(scip, &gnodearr[k]);
      }
      SCIPfreeBufferArray(scip, &node_edge);
      SCIPfreeBufferArray(scip, &node_dist);
      SCIPfreeBufferArray(scip, &node_base);
      SCIPfreeBufferArray(scip, &gnodearr);
      SCIPfreeBufferArray(scip, &nodenterms);
   }
   if( mode == TM_DIJKSTRA || graph->stp_type == STP_DHCSTP )
   {
      SCIPfreeBufferArray(scip, &dijkedge);
      SCIPfreeBufferArray(scip, &dijkdist);
   }

   SCIPfreeBufferArray(scip, &result);
   SCIPfreeBufferArray(scip, &start);
   SCIPfreeBufferArray(scip, &connected);

   return SCIP_OKAY;
}

/** run shortest path heuristic, but bias edge costs towards best current LP solution */
SCIP_RETCODE SCIPStpHeurTMRunLP(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   SCIP_HEUR*            heur,               /**< heuristic or NULL */
   int*                  result,             /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   int                   runs,               /**< number of runs */
   SCIP_Real*            cost,               /**< arc costs (uninitialized) */
   SCIP_Real*            costrev,            /**< reversed arc costs (uninitialized) */
   SCIP_Bool*            success             /**< pointer to store whether a solution could be found */
   )
{
   SCIP_VAR** vars;
   SCIP_HEURDATA* heurdata;
   SCIP_Real* xval;
   SCIP_Real* nodepriority = NULL;
   SCIP_Real maxcost = 0.0;
   int beststart;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(result != NULL);
   assert(cost != NULL);
   assert(costrev != NULL);
   assert(success != NULL);

   assert(SCIPfindHeur(scip, "TM") != NULL);
   heurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));

   vars = SCIPprobdataGetVars(scip);
   assert(vars != NULL);

   /* LP was not solved */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      xval = NULL;
   }
   else
   {
      SCIP_SOL* sol = NULL;
      SCIP_CALL(SCIPcreateSol(scip, &sol, heur));

      /* copy the current LP solution to the working solution */
      SCIP_CALL(SCIPlinkLPSol(scip, sol));

      xval = SCIPprobdataGetXval(scip, sol);

      SCIP_CALL(SCIPfreeSol(scip, &sol));
   }

   /* set (edge) result array to default */
   for( int e = 0; e < nedges; e++ )
      result[e] = UNKNOWN;

   {
      const SCIP_Real randupper = SCIPrandomGetReal(heurdata->randnumgen, 1.1, 2.5);
      const SCIP_Real randlower = SCIPrandomGetReal(heurdata->randnumgen, 1.1, randupper);

      if( xval == NULL )
      {
         BMScopyMemoryArray(cost, graph->cost, nedges);

         /* hop constraint problem? */
         if( graph->stp_type == STP_DHCSTP )
         {
            for( int e = 0; e < nedges; e++ )
            {
               if( SCIPvarGetUbGlobal(vars[e]) < 0.5 )
                  cost[e] = BLOCKED;
               else if( SCIPisGT(scip, cost[e], maxcost) && SCIPisLT(scip, cost[e], FARAWAY) )
                  maxcost = cost[e];
            }
            for( int e = 0; e < nedges; e++ )
               costrev[e] = cost[flipedge(e)];
         }
         else
         {
            if( graph->stp_type != STP_PCSPG && graph->stp_type != STP_MWCSP )
            {
               SCIP_CALL(SCIPallocBufferArray(scip, &nodepriority, nnodes));
               for( int k = 0; k < nnodes; k++ )
               {
                  if( Is_term(graph->term[k]) )
                     nodepriority[k] = (SCIP_Real) graph->grad[k];
                  else
                     nodepriority[k] = SCIPrandomGetReal(heurdata->randnumgen, 0.0, 1.0);
               }
            }

            for( int e = 0; e < nedges; e += 2 )
            {
               if( SCIPvarGetUbGlobal(vars[e + 1]) < 0.5 )
               {
                  costrev[e] = BLOCKED;
                  cost[e + 1] = BLOCKED;
               }
               else
               {
                  costrev[e] = cost[e + 1];
               }

               if( SCIPvarGetUbGlobal(vars[e]) < 0.5 )
               {
                  costrev[e + 1] = BLOCKED;
                  cost[e] = BLOCKED;
               }
               else
               {
                  costrev[e + 1] = cost[e];
               }
            }
         }
      }
      else
      {
         SCIP_Bool partrand = FALSE;
         SCIP_Bool totalrand = FALSE;

         if( (heurdata->nlpiterations == SCIPgetNLPIterations(scip) && SCIPrandomGetInt(heurdata->randnumgen, 0, 5) != 1)
               || SCIPrandomGetInt(heurdata->randnumgen, 0, 15) == 5 )
            partrand = TRUE;

         if( !partrand && (heurdata->nlpiterations == SCIPgetNLPIterations(scip) ) )
            totalrand = TRUE;
         else if( graph->stp_type == STP_DCSTP && heurdata->ncalls != 1 && SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 1
               && (graph->maxdeg[graph->source] == 1 || SCIPrandomGetInt(heurdata->randnumgen, 0, 5) == 5) )
         {
            totalrand = TRUE;
            partrand = FALSE;
         }

         assert(nodepriority == NULL);

         SCIP_CALL(SCIPallocBufferArray(scip, &nodepriority, nnodes));

         if( graph->stp_type != STP_MWCSP && graph->stp_type != STP_RMWCSP )
         {
            for( int k = 0; k < nnodes; k++ )
               nodepriority[k] = 0.0;

            for( int e = 0; e < nedges; e++ )
            {
               nodepriority[graph->head[e]] += xval[e];
               nodepriority[graph->tail[e]] += xval[e];
            }
         }

         if( graph->stp_type == STP_DHCSTP )
         {
            for( int e = 0; e < nedges; e++ )
            {
               if( SCIPvarGetUbGlobal(vars[e]) < 0.5 )
               {
                  cost[e] = BLOCKED;
               }
               else
               {
                  if( totalrand )
                  {
                     const SCIP_Real randval = SCIPrandomGetReal(heurdata->randnumgen, randlower, randupper);
                     cost[e] = graph->cost[e] * randval;
                  }
                  else
                  {
                     cost[e] = ((1.0 - xval[e]) * graph->cost[e]);
                  }
               }
               if( partrand )
               {
                  const SCIP_Real randval = SCIPrandomGetReal(heurdata->randnumgen, randlower, randupper);
                  cost[e] = cost[e] * randval;
               }
               if( SCIPisLT(scip, cost[e], BLOCKED) && SCIPisGT(scip, cost[e], maxcost) )
                  maxcost = cost[e];
               assert(SCIPisGE(scip, cost[e], 0.0));
            }
            for( int e = 0; e < nedges; e++ )
               costrev[e] = cost[flipedge(e)];
         }
         else
         {
            /* swap costs; set a high cost if the variable is fixed to 0 */
            if( graph->stp_type == STP_MWCSP || graph->stp_type == STP_RMWCSP )
            {
               for( int e = 0; e < nnodes; e++ )
                  nodepriority[e] = 0.0;
               for( int e = 0; e < nedges; e++ )
                  nodepriority[graph->head[e]] += xval[e];

               for( int e = 0; e < nedges; e++ )
               {
                  if( graph->cost[e] >= FARAWAY )
                     cost[e] = graph->cost[e];

                  if( SCIPvarGetUbLocal(vars[e]) < 0.5 )
                     cost[e] = BLOCKED;
                  else
                     cost[e] = graph->cost[e] * (1.0 - MIN(1.0, nodepriority[graph->head[e]]));
               }

               for( int e = 0; e < nedges; e++ )
               {
                  nodepriority[graph->tail[e]] += xval[e];
                  costrev[flipedge(e)] = cost[e];
               }
            }
            else
            {
               for( int e = 0; e < nedges; e += 2 )
               {
                  const SCIP_Real randval = SCIPrandomGetReal(heurdata->randnumgen, randlower, randupper);

                  if( SCIPvarGetUbLocal(vars[e + 1]) < 0.5 )
                  {
                     costrev[e] = BLOCKED;
                     cost[e + 1] = BLOCKED;
                  }
                  else
                  {
                     if( totalrand )
                        costrev[e] = graph->cost[e + 1] * randval;
                     else
                        costrev[e] = ((1.0 - xval[e + 1]) * graph->cost[e + 1]);

                     if( partrand )
                        costrev[e] = costrev[e] * randval;

                     cost[e + 1] = costrev[e];
                  }

                  if( SCIPvarGetUbLocal(vars[e]) < 0.5 )
                  {
                     costrev[e + 1] = BLOCKED;
                     cost[e] = BLOCKED;
                  }
                  else
                  {
                     if( totalrand )
                        costrev[e + 1] = graph->cost[e] * randval;
                     else
                        costrev[e + 1] = ((1.0 - xval[e]) * graph->cost[e]);

                     if( partrand )
                        costrev[e + 1] = costrev[e + 1] * randval;
                     cost[e] = costrev[e + 1];
                  }
                  assert(SCIPisGE(scip, cost[e], 0.0));
                  assert(SCIPisGE(scip, costrev[e], 0.0));
               }
            }
         }
      }

      for( int e = 0; e < nedges; e++ )
      {
         if( SCIPisZero(scip, cost[e]) )
         {
            cost[e] = SCIPepsilon(scip) * 2.0;
            assert(!SCIPisZero(scip, cost[e]));
            assert(SCIPisZero(scip, costrev[flipedge(e)]));
            costrev[flipedge(e)] = cost[e];
         }
      }

      /* can we connect the network */
      SCIP_CALL( SCIPStpHeurTMRun(scip, heurdata, graph, NULL, &beststart, result, runs, heurdata->beststartnode,
            cost, costrev, &(heurdata->hopfactor), nodepriority, maxcost, success, FALSE) );
   }

   SCIPfreeBufferArrayNull(scip, &nodepriority);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */


/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyTM)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPStpIncludeHeurTM(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeTM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitTM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_PROBDATA* probdata;
   GRAPH* graph;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   if( graph == NULL )
   {
      heurdata->stp_type = STP_SPG;
      return SCIP_OKAY;
   }
   heurdata->stp_type = graph->stp_type;
   heurdata->beststartnode = -1;
   heurdata->ncalls = 0;
   heurdata->nlpiterations = -1;
   heurdata->nexecs = 0;

#ifdef WITH_UG
   heurdata->randseed += getUgRank();
#endif

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTM)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_PROBDATA* probdata;
   SCIP_HEURDATA* heurdata;
   GRAPH* graph;
   SCIP_Real* cost;
   SCIP_Real* costrev;
   int* soledges;
   int runs;
   int nedges;
   SCIP_Bool success = FALSE;

   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   runs = 0;

   /* set the runs, i.e. number of different starting points for the heuristic */
   if( heurtiming & SCIP_HEURTIMING_BEFORENODE )
   {
      if( SCIPgetDepth(scip) > 0 )
         return SCIP_OKAY;

      runs = heurdata->initruns;
   }
   else if( ((heurtiming & SCIP_HEURTIMING_DURINGLPLOOP) && (heurdata->ncalls % heurdata->duringlpfreq == 0)) || (heurtiming & SCIP_HEURTIMING_AFTERLPLOOP) )
   {
      if( graph->stp_type == STP_PCSPG || graph->stp_type == STP_MWCSP )
         runs = 2 * heurdata->evalruns;
      else
         runs = heurdata->evalruns;
   }
   else if( heurtiming & SCIP_HEURTIMING_AFTERNODE )
   {
      if( SCIPgetDepth(scip) == 0 )
         runs = heurdata->rootruns;
      else
         runs = heurdata->leafruns;
   }

   /* increase counter for number of (TM) calls */
   heurdata->ncalls++;

   if( runs == 0 )
      return SCIP_OKAY;

   heurdata->nexecs++;

   SCIPdebugMessage("Heuristic Start\n");

   /* get all variables (corresponding to the edges) */
   vars = SCIPprobdataGetVars(scip);
   if( vars == NULL )
      return SCIP_OKAY;

   assert(vars[0] != NULL);

   nedges = graph->edges;

   /* allocate memory */
   SCIP_CALL(SCIPallocBufferArray(scip, &cost, nedges));
   SCIP_CALL(SCIPallocBufferArray(scip, &costrev, nedges));
   SCIP_CALL(SCIPallocBufferArray(scip, &soledges, nedges));

   *result = SCIP_DIDNOTFIND;

   /* call the actual heuristic */
   SCIP_CALL( SCIPStpHeurTMRunLP(scip, graph, heur, soledges, runs, cost, costrev, &success) );

   if( success )
   {
      SCIP_Real* nval;
      const int nvars = SCIPprobdataGetNVars(scip);

      SCIP_CALL(SCIPallocBufferArray(scip, &nval, nvars));

      for( int v = 0; v < nvars; v++ )
         nval[v] = (soledges[v % nedges] == (v / nedges)) ? 1.0 : 0.0;

      SCIP_CALL( SCIPStpValidateSol(scip, graph, nval, &success) );
      if( success )
      {
         SCIP_SOL* sol = NULL;
         SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, sol, heur, &success) );

         if( success )
         {
            SCIPdebugMessage("TM solution added, value %f \n",
                  graph_sol_getObj(graph->cost, soledges, SCIPprobdataGetOffset(scip), nedges));

            *result = SCIP_FOUNDSOL;
         }
      }
      SCIPfreeBufferArray(scip, &nval);
   }

   heurdata->nlpiterations = SCIPgetNLPIterations(scip);
   SCIPfreeBufferArray(scip, &soledges);
   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the TM primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPStpIncludeHeurTM(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;
   char paramdesc[SCIP_MAXSTRLEN];

   /* create TM primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heur = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecTM, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyTM) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeTM) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitTM) );

   heurdata->ncalls = 0;
   heurdata->nlpiterations = -1;
   heurdata->nexecs = 0;
   heurdata->randseed = DEFAULT_RANDSEED;

   /* add TM primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/evalruns",
         "number of runs for eval",
         &heurdata->evalruns, FALSE, DEFAULT_EVALRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/randseed",
         "random seed for heuristic",
         NULL, FALSE, DEFAULT_RANDSEED, 1, INT_MAX, paramChgdRandomseed, (SCIP_PARAMDATA*)heurdata) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/initruns",
         "number of runs for init",
         &heurdata->initruns, FALSE, DEFAULT_INITRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/leafruns",
         "number of runs for leaf",
         &heurdata->leafruns, FALSE, DEFAULT_LEAFRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/rootruns",
         "number of runs for root",
         &heurdata->rootruns, FALSE, DEFAULT_ROOTRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/duringlpfreq",
         "frequency for calling heuristic during LP loop",
         &heurdata->duringlpfreq, FALSE, DEFAULT_DURINGLPFREQ, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/type",
         "Heuristic: 0 automatic, 1 TM_SP, 2 TM_VORONOI, 3 TM_DIJKSTRA",
         &heurdata->type, FALSE, DEFAULT_TYPE, 0, 3, NULL, NULL) );
   heurdata->hopfactor = DEFAULT_HOPFACTOR;

#ifdef WITH_UG
   heurdata->randseed += getUgRank();
#endif

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen, heurdata->randseed, TRUE) );

   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "timing when heuristic should be called (%u:BEFORENODE, %u:DURINGLPLOOP, %u:AFTERLPLOOP, %u:AFTERNODE)", SCIP_HEURTIMING_BEFORENODE, SCIP_HEURTIMING_DURINGLPLOOP, SCIP_HEURTIMING_AFTERLPLOOP, SCIP_HEURTIMING_AFTERNODE);
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/timing", paramdesc,
         (int*) &heurdata->timing, TRUE, (int) HEUR_TIMING, (int) SCIP_HEURTIMING_BEFORENODE, 2 * (int) SCIP_HEURTIMING_AFTERPSEUDONODE - 1, NULL, NULL) ); /*lint !e713*/

   return SCIP_OKAY;
}
