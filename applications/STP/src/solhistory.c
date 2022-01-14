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

/**@file   solhistory.c
 * @brief  includes methods working on the (reduction) history of solutions to Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file includes methods working on the (reduction) history of solutions to Steiner tree problems
 *
 * A list of all interface methods can be found in solhistory.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h) */
//#define SCIP_DEBUG

#include "solhistory.h"
#include "probdata_stp.h"
#include "portab.h"
#include "solstp.h"
#include "mst.h"
#include "shortestpath.h"


/*
 * Local methods
 */

/** updates */
static inline
void updateorgsol(
   const GRAPH*          graph,              /**< graph data structure */
   IDX*                  curr,               /**< head of solution edge list */
   STP_Bool* RESTRICT    orgnodes,           /**< array to mark whether a node is part of the original solution */
   STP_Bool* RESTRICT    orgedges,           /**< array to mark whether an edge is part of the original solution */
   int*                  nsolnodes,          /**< pointer to store the number of nodes in the original solution */
   int*                  nsoledges           /**< pointer to store the number of edges in the original solution */
)
{
   int e = 0;
   int k = 0;

   while( curr != NULL )
   {
      const int i = curr->index;
      assert(i >= 0);

      if( orgedges[i] == FALSE )
      {
         orgedges[i] = TRUE;
         e++;
      }

      if( !(orgnodes[graph->orgtail[i]]) )
      {
         k++;
         orgnodes[graph->orgtail[i]] = TRUE;;
      }

      if( !(orgnodes[graph->orghead[i]]) )
      {
         k++;
         orgnodes[graph->orghead[i]] = TRUE;
      }

      curr = curr->parent;
   }

   (*nsolnodes) += k;
   (*nsoledges) += e;
}


/** builds history */
static
void cleanHistory(
   const GRAPH*          g,                  /**< graph structure */
   SOLHISTORY*           solhistory          /**< the solution history */
   )
{
   STP_Bool* const orgnodes = solhistory->orgnodes_isInSol;
   STP_Bool* const orgedges = solhistory->orgedges_isInSol;
   const int norgnodes = solhistory->norgnodes;
   const int norgedges = solhistory->norgedges;

   assert(norgnodes == g->orgknots);
   assert(norgedges == g->orgedges);

   for( int k = 0; k < norgnodes; k++ )
      orgnodes[k] = FALSE;

   for( int e = 0; e < norgedges; e++ )
      orgedges[e] = FALSE;
}


/** builds history */
static
SCIP_RETCODE computeHistory(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             scipsol,            /**< solution */
   const GRAPH*          graph,              /**< graph structure */
   SOLHISTORY*           solhistory          /**< the solution history */
   )
{
   GRAPH* solgraph = NULL;
   SCIP_QUEUE* queue = NULL;
   SCIP_VAR** edgevars = SCIPprobdataGetVars(scip);
   STP_Bool* const orgnodes = solhistory->orgnodes_isInSol;
   STP_Bool* const orgedges = solhistory->orgedges_isInSol;
   int* nodechild = NULL;
   int* edgeancestor = NULL;
   int nsolnodes = 0;
   int nsoledges = 0;
   int norgnodes = solhistory->norgnodes;
   const int norgedges = solhistory->norgedges;

   assert(!graph_pc_isPcMw(graph));
   assert(graph->ancestors || graph->edges == 0);
   assert(edgevars || graph->edges == 0);

   /* iterate through the list of fixed edges */
   updateorgsol(graph, graph_getFixedges(graph), orgnodes, orgedges, &nsolnodes, &nsoledges);

   for( int e = 0; e < graph->edges; e++ )
   {
      if( !SCIPisZero(scip, SCIPgetSolVal(scip, scipsol, edgevars[e])) )
         /* iterate through the list of ancestors */
         updateorgsol(graph, graph_edge_getAncestors(graph, e), orgnodes, orgedges, &nsolnodes, &nsoledges);
   }

   if( nsolnodes == 0 )
   {
      assert(graph->terms == 1);

      solhistory->nsolnodes = 1;
      solhistory->nsoledges = 0;
      solhistory->norgnodes = norgnodes;
      solhistory->norgedges = norgedges;

      return SCIP_OKAY;
   }


   assert(nsolnodes > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &nodechild, norgnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgeancestor, 2 * nsoledges) );

   /* initialize new graph */
   SCIP_CALL( graph_init(scip, &solgraph, nsolnodes, 2 * nsoledges, 1) );

   /* add vertices to new graph */
   for( int k = 0; k < norgnodes; k++ )
   {
      if( orgnodes[k] )
      {
         nodechild[k] = solgraph->knots;
         graph_knot_add(solgraph, -1);
      }
      else
      {
         nodechild[k] = -1;
      }
   }

   assert(nsolnodes == solgraph->knots);
   assert(nodechild[graph->orgsource] >= 0);

   /* set root of new graph */
   solgraph->source = nodechild[graph->orgsource];

   /* add edges to new graph */
   for( int e = 0; e < norgedges; e++ )
   {
      int tail;
      int head;

      if( !orgedges[e] )
         continue;

      tail = nodechild[graph->orgtail[e]];
      head = nodechild[graph->orghead[e]];

      assert(tail >= 0);
      assert(head >= 0);

      edgeancestor[solgraph->edges] = e;
      edgeancestor[solgraph->edges + 1] = flipedge(e);
      graph_edge_add(scip, solgraph, tail, head, 1.0, 1.0);
   }

   for( int e = 0; e < norgedges; e++ )
      orgedges[e] = FALSE;

   for( int k = 0; k < norgnodes; k++ )
      orgnodes[k] = FALSE;

   SCIP_CALL( SCIPqueueCreate(&queue, nsolnodes, 1.1) );
   SCIP_CALL( SCIPqueueInsert(queue, &(solgraph->source)) );
   solgraph->mark[solgraph->source] = FALSE;

   nsolnodes = 1;

   /* BFS from root */
   while( !SCIPqueueIsEmpty(queue) )
   {
      const int k = *((int*) SCIPqueueRemove(queue));

      /* traverse outgoing arcs */
      for( int e = solgraph->outbeg[k]; e != EAT_LAST; e = solgraph->oeat[e] )
      {
         const int head = solgraph->head[e];

         /* vertex not labeled yet? */
         if( solgraph->mark[head] )
         {
            orgedges[edgeancestor[e]] = TRUE;
            solgraph->mark[head] = FALSE;
            nsolnodes++;
            SCIP_CALL(SCIPqueueInsert(queue, &(solgraph->head[e])));
         }
      }
   }

   nsoledges = nsolnodes - 1;

   SCIPqueueFree(&queue);

   graph_free(scip, &solgraph, TRUE);

   for( int e = 0; e < norgedges; e++ )
   {
      if( orgedges[e] )
      {
         orgnodes[graph->orgtail[e]] = TRUE;
         orgnodes[graph->orghead[e]] = TRUE;
      }
   }

   SCIPfreeBufferArray(scip, &edgeancestor);
   SCIPfreeBufferArray(scip, &nodechild);

   if( graph->stp_type == STP_GSTP )
   {
      assert(graph->norgmodelterms > 0);

      norgnodes -= graph->norgmodelterms;
      nsolnodes -= graph->norgmodelterms;
      nsoledges -= graph->norgmodelterms;
      assert(nsolnodes >= 0);
      assert(nsoledges >= 1);
   }

   solhistory->nsolnodes = nsolnodes;
   solhistory->nsoledges = nsoledges;
   solhistory->norgnodes = norgnodes;
   solhistory->norgedges = norgedges;

   return SCIP_OKAY;
}


/** builds history */
static
SCIP_RETCODE computeHistoryPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             scipsol,            /**< solution */
   const GRAPH*          graph,              /**< graph structure */
   SOLHISTORY*           solhistory          /**< the solution history */
   )
{
   SCIP_VAR** edgevars = SCIPprobdataGetVars(scip);
   int* solnodequeue;
   STP_Bool* const orgnodes = solhistory->orgnodes_isInSol;
   STP_Bool* const orgedges = solhistory->orgedges_isInSol;
   int nsolnodes = 0;
   int nsoledges = 0;
   int norgnodes = solhistory->norgnodes;
   const int norgedges = solhistory->norgedges;

   assert(graph->ancestors || graph->edges == 0);
   assert(edgevars || graph->edges == 0);
   assert(graph_pc_isPcMw(graph));
   assert(graph->source >= 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &solnodequeue, norgnodes) );

   /* cover RPCSPG/RMWCSP with single node solution */
   if( graph_pc_isRootedPcMw(graph) && graph->orgsource != -1 )
   {
      SCIPdebugMessage("graph->orgsource %d \n", graph->orgsource);
      solnodequeue[nsolnodes++] = graph->orgsource;
      orgnodes[graph->orgsource] = TRUE;
   }

   for( int e = 0; e <= graph->edges; e++ )
   {
      if( e == graph->edges || !SCIPisZero(scip, SCIPgetSolVal(scip, scipsol, edgevars[e])) )
      {
         IDX* curr;

         /* iterate through the list of ancestors/fixed edges */
         if( e < graph->edges )
            curr = graph_edge_getAncestors(graph, e);
         else
            curr = graph_getFixedges(graph);

         while (curr != NULL)
         {
            const int ancestoredge = curr->index;
            if( e < graph->edges && (graph->stp_type == STP_MWCSP || graph->stp_type == STP_RMWCSP || graph->stp_type == STP_BRMWCSP) )
            {
               if( !SCIPisZero(scip, SCIPgetSolVal(scip, scipsol, edgevars[flipedge(e)])) )
               {
                  curr = curr->parent;
                  continue;
               }
            }

            if( ancestoredge < graph->norgmodeledges )
            {
               if( !orgedges[ancestoredge] )
               {
                  orgedges[ancestoredge] = TRUE;
                  nsoledges++;
               }
               if( !orgnodes[graph->orgtail[ancestoredge]] )
               {
                  orgnodes[graph->orgtail[ancestoredge]] = TRUE;
                  solnodequeue[nsolnodes++] = graph->orgtail[ancestoredge];
               }
               if( !orgnodes[graph->orghead[ancestoredge]] )
               {
                  orgnodes[graph->orghead[ancestoredge]] = TRUE;
                  solnodequeue[nsolnodes++] = graph->orghead[ancestoredge];
               }
            }
            else if( graph->orghead[ancestoredge] < graph->norgmodelknots )
            {
               if( !orgnodes[graph->orghead[ancestoredge]])
               {
                  orgnodes[graph->orghead[ancestoredge]] = TRUE;
                  solnodequeue[nsolnodes++] = graph->orghead[ancestoredge];
               }
            }
            curr = curr->parent;
         }
      }
   }

   SCIP_CALL( solstp_markPcancestors(scip, graph->pcancestors, graph->orgtail, graph->orghead, norgnodes,
         orgnodes, orgedges, solnodequeue, &nsolnodes, &nsoledges ) );

   SCIPfreeBufferArray(scip, &solnodequeue);

   solhistory->nsolnodes = nsolnodes;
   solhistory->nsoledges = nsoledges;
   solhistory->norgnodes = norgnodes;
   solhistory->norgedges = norgedges;

   return SCIP_OKAY;
}



/*
 * Interface methods
 */



/** initializes */
SCIP_RETCODE solhistory_init(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   SOLHISTORY**          solhistory          /**< the solution history */
)
{
   SOLHISTORY* sol;

   assert(scip && graph);

   SCIP_CALL( SCIPallocMemory(scip, solhistory) );
   sol = *solhistory;

   sol->nsolnodes = 0;
   sol->nsoledges = 0;
   sol->norgedges = graph->orgedges;
   sol->norgnodes = graph->orgknots;
   assert(sol->norgedges >= 0);
   assert(sol->norgnodes >= 1);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(sol->orgnodes_isInSol), sol->norgnodes) );

   if( sol->norgedges > 0 )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(sol->orgedges_isInSol), sol->norgedges) );

   return SCIP_OKAY;
}


/** frees */
void solhistory_free(
   SCIP*                 scip,               /**< SCIP data structure */
   SOLHISTORY**          solhistory          /**< the solution history */
   )
{
   SOLHISTORY* sol;

   assert(scip && solhistory);
   sol = *solhistory;
   assert(sol);

   SCIPfreeMemoryArrayNull(scip, &(sol->orgedges_isInSol));
   SCIPfreeMemoryArray(scip, &(sol->orgnodes_isInSol));
   SCIPfreeMemory(scip, solhistory);
}


/** builds history */
SCIP_RETCODE solhistory_computeHistory(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             scipsol,            /**< solution */
   const GRAPH*          g,                  /**< graph structure */
   SOLHISTORY*           solhistory          /**< the solution history */
   )
{
   assert(scip && g && solhistory && scipsol);

   cleanHistory(g, solhistory);

   if( graph_pc_isPcMw(g) )
   {
      SCIP_CALL( computeHistoryPcMw(scip, scipsol, g, solhistory) );
   }
   else
   {
      SCIP_CALL( computeHistory(scip, scipsol, g, solhistory) );
   }

   return SCIP_OKAY;
}
