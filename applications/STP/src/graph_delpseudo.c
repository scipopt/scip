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

/**@file   graph_delpseudo.c
 * @brief  includes graph pseudo deletion methods for Steiner problems
 * @author Daniel Rehfeldt
 *
 * Graph node or edge pseudo-deletion methods (aka replacement methods) for Steiner problems
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h) */

#include "graph.h"
#include "portab.h"


#define STP_DELPSEUDO_NOEDGE    -1
#define STP_DELPSEUDO_SKIPEDGE  -2



/** Internal data for pseudo-deletion.
 *  Is used for both pseudo-elimination of vertex and edge.
 *  In the latter case, the adjacency information is w.r.t. the head of the edge. */
typedef struct pseudo_deletion
{
   SCIP_Real*            ecost;              /**< edge cost */
   SCIP_Real*            ecostrev;           /**< reverse edge cost */
   SCIP_Real*            ecostreal;          /**< reverse edge cost */
   SCIP_Real*            ecost_adapt;        /**< edge costs to adapt or NULL */
   SCIP_Real*            ecost_adaptrev;     /**< edge costs to adapt or NULL */
   int*                  incedge;            /**< incident edges */
   int*                  adjvert;            /**< adjacent vertices */
   int*                  neigbedge;          /**< neighboring edges array */
   SCIP_Real             vertexprize;        /**< prize for PC (either of head of edge, or of vertex itself) */
   int                   degree;             /**< degree of vertex to be deleted */
   int                   ancestorsnode;      /**< ancestor node for PC (either of head of edge, or of vertex itself) */
   int                   edge;               /**< edge for pseudo-deletion or UNKNOWN otherwise */
} DELPSEUDO;



/*
 * Local methods
 */

/** can edge in pseudo-elimination method be cut off? */
inline static
SCIP_Bool isCutoffEdge(
   SCIP*                 scip,               /**< SCIP data */
   const SCIP_Real*      cutoffs,            /**< cutoff values for each incident edge */
   const SCIP_Real*      cutoffsrev,         /**< revere cutoff values (or NULL if undirected) */
   const SCIP_Real*      ecost,              /**< edge cost*/
   const SCIP_Real*      ecostrev,           /**< reverse edge cost */
   SCIP_Real             prize,              /**< prize if PcMw, 0.0 otherwise */
   int                   edgeidx1,           /**< index of first edge to be checked (wrt provided arrays) */
   int                   edgeidx2,           /**< index of second edge to be checked (wrt provided arrays) */
   int                   cutoffidx           /**< index for cutoff array */
   )
{
   SCIP_Real newcost;

   assert(edgeidx1 != edgeidx2);

   if( cutoffs == NULL )
      return FALSE;

   newcost = ecostrev[edgeidx1] + ecost[edgeidx2] - prize;

   /* NOTE: don't replace by GT, to keep epsilon change valid! */
   if( !SCIPisGT(scip, newcost, cutoffs[cutoffidx]) )
      return FALSE;

   if( cutoffsrev != NULL )
   {
      const SCIP_Real newcostrev = ecost[edgeidx1] + ecostrev[edgeidx2];

      if( !SCIPisGT(scip, newcostrev, cutoffsrev[cutoffidx]) )
         return FALSE;
   }

   return TRUE;
}

#ifndef NDEBUG
/** in edge deletion mode? */
static
SCIP_Bool delPseudoIsEdgeDeletionMode(
   const DELPSEUDO*     delpseudo           /**< data */
)
{
   assert(delpseudo);
   assert(delpseudo->edge == UNKNOWN || delpseudo->edge >= 0);

   return (delpseudo->edge != UNKNOWN);
}
#endif


/** gets position of deletion edge in replacement arrays */
static inline
int delPseudoGetEdgePosition(
   const DELPSEUDO*     delpseudo           /**< data */
)
{
   const int edge_rev = flipedge(delpseudo->edge);
   const int degree = delpseudo->degree;
   const int* const incedges = delpseudo->incedge;
   int pos = -1;

   assert(edge_rev >= 0);
   assert(delPseudoIsEdgeDeletionMode(delpseudo));

   for( int i = 0; i < degree; i++ )
   {
      if( edge_rev == incedges[i] )
      {
         pos = i;
         break;
      }
   }

   assert(0 <= pos && pos < degree);

   return pos;
}


/** initializes */
static
SCIP_RETCODE delPseudoInit(
   SCIP*                 scip,               /**< SCIP data */
   const SCIP_Real*      edgecosts,          /**< edge costs for cutoff */
   const REDCOST*        redcostdata,        /**< reduced cost data for adaptation or NULL */
   int                   vertex,             /**< the vertex */
   GRAPH*                g,                  /**< graph */
   DELPSEUDO*            delpseudo           /**< data */
)
{
   SCIP_Real* RESTRICT ecost;
   SCIP_Real* RESTRICT ecostrev;
   SCIP_Real* RESTRICT ecostreal;
   SCIP_Real* RESTRICT ecostadapt = NULL;
   SCIP_Real* RESTRICT ecostadaptrev = NULL;
   int* RESTRICT incedge;
   int* RESTRICT adjvert;
   int* RESTRICT neigbedge;
   int edgecount = 0;
   const int redcost_nlevels = redcostdata ? redcosts_getNlevels(redcostdata) : -1;

   if( Is_term(g->term[vertex]) )
   {
      assert(graph_pc_isPcMw(g));
      assert(!graph_pc_knotHasMaxPrize(g, vertex) );
      assert(!delPseudoIsEdgeDeletionMode(delpseudo));

      delpseudo->vertexprize = g->prize[vertex];
      delpseudo->ancestorsnode = vertex;

      /* NOTE: for degree 3 the deletion is always possible.
       * Thus, we can already manipulate the graph here */
      graph_pc_termToNonTerm(scip, g, vertex);

      delpseudo->degree = g->grad[vertex];

      assert(delpseudo->vertexprize > 0.0);
      assert(delpseudo->degree == 3);
   }
   else
   {
      delpseudo->vertexprize = 0.0;
      delpseudo->ancestorsnode = -1;
      delpseudo->degree = g->grad[vertex];

      if( g->pcancestors && g->pcancestors[vertex] )
      {
         assert(graph_pc_isPcMw(g));
         assert(SCIPisZero(scip, g->prize[vertex]));

         delpseudo->ancestorsnode = vertex;
      }
   }

   if( redcostdata )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(ecostadapt), redcost_nlevels * STP_DELPSEUDO_MAXGRAD) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(ecostadaptrev), redcost_nlevels * STP_DELPSEUDO_MAXGRAD) );
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(ecost), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(ecostrev), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(ecostreal), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(incedge), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(adjvert), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(neigbedge), STP_DELPSEUDO_MAXNEDGES) );

   /* save the state of all incident edges */
   for( int e = g->outbeg[vertex]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int e_rev = flipedge(e);
      assert(e >= 0);

      if( redcostdata )
      {
         for( int i = 0; i < redcost_nlevels; i++ )
         {
            const SCIP_Real* const redcosts = redcosts_getEdgeCosts(redcostdata, i);
            const int position = i * STP_DELPSEUDO_MAXGRAD + edgecount;
            ecostadapt[position] = redcosts[e];
            ecostadaptrev[position] = redcosts[e_rev];
         }
      }

      incedge[edgecount] = e;
      ecostreal[edgecount] = g->cost[e];
      ecost[edgecount] = edgecosts[e];
      ecostrev[edgecount] = edgecosts[e_rev];
      adjvert[edgecount++] = g->head[e];

      assert(edgecount <= STP_DELPSEUDO_MAXGRAD);
   }

   assert(edgecount == delpseudo->degree);

   delpseudo->ecost = ecost;
   delpseudo->ecostrev = ecostrev;
   delpseudo->ecostreal = ecostreal;
   delpseudo->ecost_adapt = ecostadapt;
   delpseudo->ecost_adaptrev = ecostadaptrev;
   delpseudo->incedge = incedge;
   delpseudo->adjvert = adjvert;
   delpseudo->neigbedge = neigbedge;

   return SCIP_OKAY;
}


/** initializes for check */
static
SCIP_RETCODE delPseudoInitForCheck(
   SCIP*                 scip,               /**< SCIP data */
   const GRAPH*          g,                  /**< graph */
   const SCIP_Real*      edgecosts,          /**< edge costs for cutoff */
   int                   vertex,             /**< the vertex */
   DELPSEUDO*            delpseudo           /**< data */
)
{
   SCIP_Real* RESTRICT ecost;
   SCIP_Real* RESTRICT ecostrev;
   SCIP_Real* RESTRICT ecostreal;
   int* RESTRICT incedge;
   int* RESTRICT adjvert;
   int edgecount = 0;

   if( Is_term(g->term[vertex]) )
   {
      assert(graph_pc_isPcMw(g));
      assert(!graph_pc_knotHasMaxPrize(g, vertex) );
      assert(!delPseudoIsEdgeDeletionMode(delpseudo));

      delpseudo->vertexprize = g->prize[vertex];
      delpseudo->ancestorsnode = vertex;

      delpseudo->degree = 3;

      assert(delpseudo->vertexprize > 0.0);
      assert(delpseudo->degree == graph_pc_realDegree(g, vertex, FALSE));
   }
   else
   {
      delpseudo->vertexprize = 0.0;
      delpseudo->ancestorsnode = -1;
      delpseudo->degree = g->grad[vertex];

      if( g->pcancestors && g->pcancestors[vertex] )
      {
         assert(graph_pc_isPcMw(g));
         assert(SCIPisZero(scip, g->prize[vertex]));

         delpseudo->ancestorsnode = vertex;
      }
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(ecost), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(ecostrev), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(ecostreal), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(incedge), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(adjvert), STP_DELPSEUDO_MAXGRAD) );

   /* save the state of all incident edges */
   for( int e = g->outbeg[vertex]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int e_rev = flipedge(e);
      assert(e >= 0);

      incedge[edgecount] = e;
      ecostreal[edgecount] = g->cost[e];
      ecost[edgecount] = edgecosts[e];
      ecostrev[edgecount] = edgecosts[e_rev];
      adjvert[edgecount++] = g->head[e];

      assert(edgecount <= STP_DELPSEUDO_MAXGRAD);
   }

   assert(edgecount == delpseudo->degree);

   delpseudo->ecost = ecost;
   delpseudo->ecostrev = ecostrev;
   delpseudo->ecostreal = ecostreal;
   delpseudo->incedge = incedge;
   delpseudo->adjvert = adjvert;
   delpseudo->ecost_adapt = NULL;
   delpseudo->ecost_adaptrev = NULL;
   delpseudo->neigbedge = NULL;

   return SCIP_OKAY;
}


/** pseudo-eliminates vertex */
static
SCIP_RETCODE delPseudoDeleteVertex(
   SCIP*                 scip,               /**< SCIP data */
   int                   vertex,             /**< the vertex */
   GRAPH*                g,                  /**< graph */
   REDCOST*              redcostdata,        /**< reduced cost data for adaptation or NULL */
   DELPSEUDO*            delpseudo           /**< data */
)
{
   SINGLETONANS ancestors[STP_DELPSEUDO_MAXGRAD];
   const SCIP_Real* ecostreal = delpseudo->ecostreal;
   const SCIP_Real* ecost_adapt = delpseudo->ecost_adapt;
   const SCIP_Real* ecost_adaptrev = delpseudo->ecost_adaptrev;
   const int* incedge = delpseudo->incedge;
   const int* neigbedge = delpseudo->neigbedge;
   const int* adjvert = delpseudo->adjvert;
   const int degree = delpseudo->degree;
   const int nspareedges = degree; /* todo we might want to allow additional edges to be inserted */
   int edgecount = 0;
   int replacecount = 0;
   int pseudoancestor = -1;

   for( int i = 0; i < degree; i++ )
   {
      const int e = incedge[i];
      SCIP_CALL( graph_singletonAncestors_init(scip, g, e, &(ancestors[i])) );
   }

   assert(EQ(delpseudo->vertexprize, 0.0) || graph_pc_isPc(g));

   for( int i = 0; i < degree - 1; i++ )
   {
      for( int j = i + 1; j < degree; j++ )
      {
         const SCIP_Bool skipedge = (neigbedge[edgecount] == STP_DELPSEUDO_SKIPEDGE);

         /* do we need to insert edge at all? */
         if( !skipedge )
         {
            SCIP_Bool conflict;
            int newijedge;
            const SCIP_Real newijcost = ecostreal[i] + ecostreal[j] - delpseudo->vertexprize;
            const int oldincedge = incedge[(replacecount == nspareedges) ? replacecount - 1 : replacecount];
            const int oldijedge = neigbedge[edgecount];
#ifndef NDEBUG
            const int oldinctail = g->tail[oldincedge];
            const int oldinchead = g->head[oldincedge];
#endif
            assert(replacecount <= nspareedges);
            assert(replacecount < nspareedges || neigbedge[edgecount] != STP_DELPSEUDO_NOEDGE);

            SCIP_CALL( graph_edge_reinsert(scip, g, oldincedge, adjvert[i], adjvert[j], newijcost,
                  delpseudo->ancestorsnode, &(ancestors[i]), &(ancestors[j]), &newijedge, &conflict) );

            assert(!conflict);

            /* has a new edge been inserted (or the existing one been updated)? */
            if( newijedge >= 0 )
            {

               assert(g->tail[newijedge] == adjvert[i]);
               assert(g->head[newijedge] == adjvert[j]);

               if( redcostdata)
               {
                  const int redcost_nlevels = redcosts_getNlevels(redcostdata);
                  assert(ecost_adapt && ecost_adaptrev);

                  for( int level = 0; level < redcost_nlevels; level++ )
                  {
                     SCIP_Real* const redcosts = redcosts_getEdgeCosts(redcostdata, level);
                     const int offset = level * STP_DELPSEUDO_MAXGRAD;

                     redcosts[newijedge] = ecost_adaptrev[offset + i] + ecost_adapt[offset + j];
                     redcosts[flipedge(newijedge)] =  ecost_adaptrev[offset + j] + ecost_adapt[offset + i];
                  }
               //    graph_edge_printInfo(g, newijedge);
               //   printf("...with costs %f, %f \n", edgecosts_adapt[newijedge],  edgecosts_adapt[flipedge(newijedge)] );
               }

               if( pseudoancestor == -1 )
               {
                  graph_addPseudoAncestor(g, &pseudoancestor);
                  assert(pseudoancestor >= 0);
               }

               SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, newijedge, pseudoancestor, g) );
            }

            /* does no original edge exist? */
            if( oldijedge == STP_DELPSEUDO_NOEDGE )
            {
               replacecount++;
               assert(newijedge >= 0);
               assert(oldinctail != g->tail[oldincedge] || oldinchead != g->head[oldincedge]);
            }
            else
            {
               assert(newijedge == oldijedge || newijedge == - 1);
               assert(newijedge != oldincedge);
               assert(oldinctail == g->tail[oldincedge] && oldinchead == g->head[oldincedge]);
            }
         }

         edgecount++;
         assert(edgecount <= STP_DELPSEUDO_MAXNEDGES);
      }
   }

   /* delete remaining edges */
   graph_knot_del(scip, g, vertex, TRUE);

   for( int i = 0; i < degree; i++ )
      graph_singletonAncestors_freeMembers(scip, &(ancestors[i]));

   return SCIP_OKAY;
}


/** gets replacement edges; helper function for pseudo-elimination */
static
SCIP_RETCODE delPseudoGetReplaceEdges(
   SCIP*                 scip,               /**< SCIP data */
   const GRAPH*          g,                  /**< graph */
   const SCIP_Real*      cutoffs,            /**< cutoff values for each incident edge (or NULL) */
   const SCIP_Real*      cutoffsrev,         /**< reverse edge cutoff values (or NULL if undirected or non-existent) */
   DELPSEUDO*            delpseudo,          /**< data */
   SCIP_Bool*            success             /**< enough replace edges available?  */
)
{
   int* hasharr;
   int edgecount = 0;
   int replacecount = 0;
   const int degree = delpseudo->degree;
   const int nspareedges = degree;
   const int *incedge = delpseudo->incedge;
   const SCIP_Real *ecost = delpseudo->ecost;
   const SCIP_Real *ecostrev = delpseudo->ecostrev;
   const int *adjvert = delpseudo->adjvert;
   int* RESTRICT neigbedge = delpseudo->neigbedge;
   const SCIP_Real vertexprize = delpseudo->vertexprize;

   assert(scip && g && ecost && neigbedge);
   assert(degree >= 0 && degree <= STP_DELPSEUDO_MAXGRAD);

   *success = TRUE;

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &hasharr, graph_pseudoAncestorsGetHashArraySize(g->pseudoancestors)) );

   for( int i = 0; i < STP_DELPSEUDO_MAXNEDGES; i++ )
      neigbedge[i] = STP_DELPSEUDO_NOEDGE;

   for( int i = 0; i < degree - 1 && *success; i++ )
   {
      const int adjvertex = adjvert[i];
      const int iedge = incedge[i];

      graph_pseudoAncestors_hashEdge(g->pseudoancestors, iedge, hasharr);

      for( int j = i + 1; j < degree; j++ )
      {
         SCIP_Bool skipedge = isCutoffEdge(scip, cutoffs, cutoffsrev, ecost, ecostrev, vertexprize, i, j, edgecount);

         if( !skipedge )
         {
            const int jedge = incedge[j];
            skipedge = graph_pseudoAncestors_edgeIsHashed(g->pseudoancestors, jedge, hasharr);
         }

         edgecount++;
         assert(edgecount <= STP_DELPSEUDO_MAXNEDGES);

         /* can edge be discarded? */
         if( skipedge )
         {
            neigbedge[edgecount - 1] = STP_DELPSEUDO_SKIPEDGE;
         }
         else
         {
            int e;

            /* check whether edge already exists */
            for( e = g->outbeg[adjvertex]; e != EAT_LAST; e = g->oeat[e] )
            {
               if( g->head[e] == adjvert[j] )
               {
                  assert(e >= 0);
                  neigbedge[edgecount - 1] = e;
                  break;
               }
            }

            if( e != EAT_LAST )
               continue;

            /* not enough spare edges? */
            if( ++replacecount > nspareedges )
            {
               *success = FALSE;
               break;
            }
         }
      } /* inner neighbor loop */

      graph_pseudoAncestors_unhashEdge(g->pseudoancestors, iedge, hasharr);
   }  /* outer neighbor loop */

   SCIPfreeCleanBufferArray(scip, &hasharr);

   return SCIP_OKAY;
}


/** checks whether replacement is possible */
static
SCIP_RETCODE delPseudoCheckReplacement(
   SCIP*                 scip,               /**< SCIP data */
   const GRAPH*          g,                  /**< graph */
   const SCIP_Real*      cutoffs,            /**< cutoff values for each incident edge (or NULL) */
   const SCIP_Real*      cutoffsrev,         /**< reverse edge cutoff values (or NULL if undirected or non-existent) */
   DELPSEUDO*            delpseudo,          /**< data */
   SCIP_Bool*            success             /**< enough replace edges available?  */
)
{
   int* hasharr;
   int edgecount = 0;
   int replacecount = 0;
   const int degree = delpseudo->degree;
   const int nspareedges = degree;
   const int *incedge = delpseudo->incedge;
   const SCIP_Real *ecost = delpseudo->ecost;
   const SCIP_Real *ecostrev = delpseudo->ecostrev;
   const int *adjvert = delpseudo->adjvert;

   const SCIP_Real vertexprize = delpseudo->vertexprize;

   assert(scip && g);
   assert(incedge && ecost && ecostrev && adjvert);
   assert(degree >= 0 && degree <= STP_DELPSEUDO_MAXGRAD);

   *success = TRUE;

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &hasharr, graph_pseudoAncestorsGetHashArraySize(g->pseudoancestors)) );

   for( int i = 0; i < degree - 1 && *success; i++ )
   {
      const int adjvertex = adjvert[i];
      const int iedge = incedge[i];

      graph_pseudoAncestors_hashEdge(g->pseudoancestors, iedge, hasharr);

      for( int j = i + 1; j < degree; j++ )
      {
         SCIP_Bool skipedge = isCutoffEdge(scip, cutoffs, cutoffsrev, ecost, ecostrev, vertexprize, i, j, edgecount);

         if( !skipedge )
         {
            const int jedge = incedge[j];
            skipedge = graph_pseudoAncestors_edgeIsHashed(g->pseudoancestors, jedge, hasharr);
         }

         edgecount++;
         assert(edgecount <= STP_DELPSEUDO_MAXNEDGES);

         /* can edge not be discarded? */
         if( !skipedge )
         {
            int e;

            /* check whether edge already exists */
            for( e = g->outbeg[adjvertex]; e != EAT_LAST; e = g->oeat[e] )
            {
               if( g->head[e] == adjvert[j] )
               {
                  assert(e >= 0);
                  break;
               }
            }

            if( e != EAT_LAST )
               continue;

            /* not enough spare edges? */
            if( ++replacecount > nspareedges )
            {
               *success = FALSE;
               break;
            }
         }
      } /* inner neighbor loop */

      graph_pseudoAncestors_unhashEdge(g->pseudoancestors, iedge, hasharr);
   }  /* outer neighbor loop */

   SCIPfreeCleanBufferArray(scip, &hasharr);

   return SCIP_OKAY;
}


/** initializes for edge elimination */
static
SCIP_RETCODE delPseudoEdgeInit(
   SCIP*                 scip,               /**< SCIP data */
   const SCIP_Real*      edgecosts,          /**< edge costs for cutoff */
   const SCIP_Real*      edgecosts_adapt,    /**< edge costs that should be adapted */
   GRAPH*                g,                  /**< graph */
   DELPSEUDO*            delpseudo           /**< data */
)
{
   const int edge = delpseudo->edge;
   const int head = g->head[edge];
   assert(delPseudoIsEdgeDeletionMode(delpseudo));

   /* would need to change for REDCOST */
   assert(!edgecosts_adapt && "currently not supported!");

   SCIP_CALL( delPseudoInit(scip, edgecosts, NULL, head, g, delpseudo) );

   return SCIP_OKAY;
}


/** gets replacement edges for edge elimination; helper function for pseudo-elimination */
static
SCIP_RETCODE delPseudoEdgeGetReplaceEdges(
   SCIP*                 scip,               /**< SCIP data */
   const GRAPH*          g,                  /**< graph */
   const SCIP_Real*      cutoffs,            /**< cutoff values for each incident edge (or NULL) */
   const SCIP_Real*      cutoffsrev,         /**< reverse edge cutoff values (or NULL if undirected or non-existent) */
   DELPSEUDO*            delpseudo,          /**< data */
   SCIP_Bool*            success             /**< enough replace edges available?  */
)
{
   int* hasharr;
   int replacecount = 0;
   const int degree = delpseudo->degree;
   const int *incedge = delpseudo->incedge;
   const SCIP_Real *ecost = delpseudo->ecost;
   const SCIP_Real *ecostrev = delpseudo->ecostrev;
   const int *adjvert = delpseudo->adjvert;
   int* RESTRICT neigbedge = delpseudo->neigbedge;
   const SCIP_Real vertexprize = delpseudo->vertexprize;
   const int edge_pos = delPseudoGetEdgePosition(delpseudo);
   const int edge = delpseudo->edge;
   const int tail = g->tail[edge];

   assert(scip && g && ecost && neigbedge);
   assert(degree >= 0 && degree <= STP_DELPSEUDO_MAXGRAD);
   assert(delPseudoIsEdgeDeletionMode(delpseudo));
   assert(EQ(g->cost[edge], ecostrev[edge_pos]));
   assert(EQ(g->cost[flipedge(edge)], ecost[edge_pos]));

   *success = TRUE;

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &hasharr, graph_pseudoAncestorsGetHashArraySize(g->pseudoancestors)) );
   graph_pseudoAncestors_hashEdge(g->pseudoancestors, edge, hasharr);

   for( int i = 0; i < degree; i++ )
      neigbedge[i] = STP_DELPSEUDO_NOEDGE;

   for( int i = 0; i < degree; i++ )
   {
      if( i != edge_pos )
      {
         const int adjvertex = adjvert[i];
         const int iedge = incedge[i];
         SCIP_Bool skipedge = isCutoffEdge(scip, cutoffs, cutoffsrev, ecost, ecostrev, vertexprize, i, edge_pos, i);

         SCIPdebugMessage("in (degree=%d): %d->%d cutoff=%f \n", degree, tail, adjvertex, cutoffs[i]);

         assert((iedge / 2) != (edge / 2));

         if( !skipedge )
            skipedge = graph_pseudoAncestors_edgeIsHashed(g->pseudoancestors, iedge, hasharr);

         if( skipedge )
         {
            SCIPdebugMessage("skip edge %d->%d \n", tail, g->head[iedge]);
            neigbedge[i] = STP_DELPSEUDO_SKIPEDGE;
         }
         else
         {
            int e;

            /* check whether edge already exists */
            for( e = g->outbeg[adjvertex]; e != EAT_LAST; e = g->oeat[e] )
            {
               if( g->head[e] == tail )
               {
                  assert(e >= 0);
                  neigbedge[i] = e;
                  SCIPdebugMessage("edge exist: %d->%d  \n", g->tail[e], g->head[e]);
                  break;
               }
            }

            /* edge does not exist? */
            if( e == EAT_LAST )
            {
               /* more than one edge to be reinserted? */
               if( ++replacecount > 1 )
               {
                  *success = FALSE;
                  break;
               }
            }
         }
      }
      else
      {
         assert(EQ(cutoffs[i], -1.0));
      }
   }  /* neighbor loop */


   graph_pseudoAncestors_unhashEdge(g->pseudoancestors, edge, hasharr);
   SCIPfreeCleanBufferArray(scip, &hasharr);

   return SCIP_OKAY;
}


/** pseudo-eliminates edge */
static
SCIP_RETCODE delPseudoEdgeDeleteEdge(
   SCIP*                 scip,               /**< SCIP data */
   GRAPH*                g,                  /**< graph */
   SCIP_Real*            edgecosts_adapt,    /**< edge costs that should be adapted */
   DELPSEUDO*            delpseudo           /**< data */
)
{
   SINGLETONANS ancestors[STP_DELPSEUDO_MAXGRAD];
   const SCIP_Real* ecostreal = delpseudo->ecostreal;
   const SCIP_Real* ecost_adapt = delpseudo->ecost_adapt;
   const SCIP_Real* ecost_adaptrev = delpseudo->ecost_adaptrev;
   const int* incedge = delpseudo->incedge;
   const int* neigbedge = delpseudo->neigbedge;
   const int* adjvert = delpseudo->adjvert;
   const int degree = delpseudo->degree;
   const int edge_pos = delPseudoGetEdgePosition(delpseudo);
   const int edge = delpseudo->edge;
   const int tail = g->tail[edge];
   int replacecount = 0;
#ifndef NDEBUG
   const int head = g->head[edge];
#endif

   for( int i = 0; i < degree; i++ )
   {
      const int e = incedge[i];
      SCIP_CALL( graph_singletonAncestors_init(scip, g, e, &(ancestors[i])) );
   }

   assert(EQ(delpseudo->vertexprize, 0.0));

   for( int i = 0; i < degree; i++ )
   {
      if( i != edge_pos )
      {
         const SCIP_Bool skipedge = (neigbedge[i] == STP_DELPSEUDO_SKIPEDGE);

         if( !skipedge )
         {
            SCIP_Bool conflict;
            int newijedge;
            const SCIP_Real newijcost = ecostreal[i] + ecostreal[edge_pos] - delpseudo->vertexprize;
            const int oldincedge = (replacecount == 0)? flipedge(edge) : -1;
            const int oldAdjToTailEdge = neigbedge[i];

            assert(replacecount <= 1);
            assert(replacecount == 0 || oldAdjToTailEdge != STP_DELPSEUDO_NOEDGE);
            assert(adjvert[i] != head);

            SCIP_CALL( graph_edge_reinsert(scip, g, oldincedge, adjvert[i], tail, newijcost,
                  delpseudo->ancestorsnode, &(ancestors[i]), &(ancestors[edge_pos]), &newijedge, &conflict) );

            assert(!conflict);

            /* has a new edge been inserted (or the existing one been updated)? */
            if( newijedge >= 0 )
            {
               assert(g->tail[newijedge] == adjvert[i] && g->head[newijedge] == tail);

               if( edgecosts_adapt)
               {
                  assert(ecost_adapt && ecost_adaptrev);
                  edgecosts_adapt[newijedge] = ecost_adaptrev[i] + ecost_adapt[edge_pos];
                  edgecosts_adapt[flipedge(newijedge)] =  ecost_adaptrev[edge_pos] + ecost_adapt[i];
               }

#ifdef SCIP_DEBUG
               SCIPdebugMessage("pseudo-edge-elimination reinserted edge: ");
               graph_edge_printInfo(g, newijedge);
#endif

               // todo add conflicts
            }

            /* does no original edge exist? */
            if( oldAdjToTailEdge == STP_DELPSEUDO_NOEDGE )
            {
               replacecount++;
               assert(newijedge >= 0);
            }
            else
            {
               assert(newijedge == oldAdjToTailEdge || newijedge == -1);
               assert(newijedge != oldincedge || newijedge == -1);
            }
         }
      }
   }

   for( int i = 0; i < degree; i++ )
      graph_singletonAncestors_freeMembers(scip, &(ancestors[i]));

   if( replacecount == 0 )
   {
      assert(g->tail[edge] == tail && g->head[edge] == head);
#ifdef SCIP_DEBUG
      SCIPdebugMessage("pseudo-edge-elimination delete edge: ");
      graph_edge_printInfo(g, edge);
#endif
      graph_edge_del(scip, g, edge, TRUE);
   }
   else
   {
      assert(g->head[edge] != head);
   }

   return SCIP_OKAY;
}


/** frees data */
static
void delPseudoFreeData(
   SCIP*                 scip,               /**< SCIP data */
   DELPSEUDO*            delpseudo           /**< data */
)
{
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->neigbedge), STP_DELPSEUDO_MAXNEDGES);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->adjvert), STP_DELPSEUDO_MAXGRAD);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->incedge), STP_DELPSEUDO_MAXGRAD);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->ecostreal), STP_DELPSEUDO_MAXGRAD);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->ecostrev), STP_DELPSEUDO_MAXGRAD);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->ecost), STP_DELPSEUDO_MAXGRAD);

   if( delpseudo->ecost_adapt )
   {
      assert(delpseudo->ecost_adaptrev);
      SCIPfreeBufferArray(scip, &(delpseudo->ecost_adaptrev) );
      SCIPfreeBufferArray(scip, &(delpseudo->ecost_adapt) );
   }
}


/** frees data from check */
static
void delPseudoFreeDataForCheck(
   SCIP*                 scip,               /**< SCIP data */
   DELPSEUDO*            delpseudo           /**< data */
)
{
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->adjvert), STP_DELPSEUDO_MAXGRAD);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->incedge), STP_DELPSEUDO_MAXGRAD);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->ecostreal), STP_DELPSEUDO_MAXGRAD);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->ecostrev), STP_DELPSEUDO_MAXGRAD);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->ecost), STP_DELPSEUDO_MAXGRAD);

   assert(!delpseudo->ecost_adaptrev && !delpseudo->ecost_adapt && !delpseudo->neigbedge);
}


/** does the path replacement */
static inline
SCIP_RETCODE delPseudoPathCreatePseudoAncestorTuple(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   edge1,              /**< first edge */
   int                   edge2               /**< second edge */
   )
{
   int pseudoancestor;

   assert(edge1 / 2 != edge2 / 2);
   assert(graph_edge_isInRange(g, edge1));
   assert(graph_edge_isInRange(g, edge2));

   graph_addPseudoAncestor(g, &pseudoancestor);
   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, edge1, pseudoancestor, g) );
   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, edge2, pseudoancestor, g) );

   return SCIP_OKAY;
}

/** does the path replacement */
static
SCIP_RETCODE delPseudoPath(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   edge,               /**< the edge to be replaced, mind the orientation! */
   int                   edge_pathtail,      /**< tail edge of path */
   int                   edge_pathhead,      /**< head edge of path */
   SCIP_Real*            edgecosts_adapt     /**< costs to adapt or NULL */
   )
{
   const int old_tail = g->tail[edge];
   const int old_head = g->head[edge];
   const SCIP_Real edgecost_adapt
      = edgecosts_adapt ? (edgecosts_adapt[edge] + edgecosts_adapt[edge_pathtail] + edgecosts_adapt[edge_pathhead]) : -1.0;
   const SCIP_Real path_cost = g->cost[edge] + g->cost[edge_pathtail] + g->cost[edge_pathhead];
   const int path_tail = g->tail[edge_pathtail];
   const int path_head = g->head[edge_pathhead];
   const int newedge = graph_edge_redirect(scip, g, edge, path_tail, path_head, path_cost, FALSE, TRUE);

   /* is there a new edge? */
   if( newedge >= 0 )
   {
      SCIP_Bool conflict = FALSE;

#ifdef SCIP_DEBUG
      printf("have path replace edge: \n");
      graph_edge_printInfo(g, newedge);
#endif

      if( edgecost_adapt )
         edgecosts_adapt[newedge] = edgecost_adapt;

      if( !g->ancestors[newedge] || !g->ancestors[flipedge(newedge)] )
      {
         const int edge_even = Edge_even(newedge);
         assert(graph_typeIsUndirected(g));
         assert(edge_even == flipedge(newedge) || edge_even == edge);

         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[edge_even]), graph_edge_getAncestors(g, edge_pathtail), NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[edge_even]), graph_edge_getAncestors(g, edge_pathhead), NULL) );
      }
      else
      {
         assert(!graph_typeIsUndirected(g));
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[newedge]), g->ancestors[edge_pathtail], NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[newedge]), g->ancestors[edge_pathhead], NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[flipedge(newedge)]), g->ancestors[flipedge(edge_pathtail)], NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[flipedge(newedge)]), g->ancestors[flipedge(edge_pathhead)], NULL) );
      }

      SCIP_CALL( graph_pseudoAncestors_appendCopyEdge(scip, newedge, edge_pathtail, FALSE, g, &conflict) );
      assert(!conflict);

      SCIP_CALL( graph_pseudoAncestors_appendCopyEdge(scip, newedge, edge_pathhead, FALSE, g, &conflict) );
      assert(!conflict);

      /* create new ancestor relations */
      for( int e = g->outbeg[old_tail]; e != EAT_LAST; e = g->oeat[e] )
      {
         SCIP_CALL( delPseudoPathCreatePseudoAncestorTuple(scip, g, e, newedge) );
      }

      for( int e = g->outbeg[old_head]; e != EAT_LAST; e = g->oeat[e] )
      {
         SCIP_CALL( delPseudoPathCreatePseudoAncestorTuple(scip, g, e, newedge) );
      }
   }

   return SCIP_OKAY;
}



/*
 * Interface methods
 */



/** pseudo delete node, i.e. reconnect neighbors; maximum degree of STP_DELPSEUDO_MAXGRAD! */
SCIP_RETCODE graph_knot_delPseudo(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   const SCIP_Real*      cutoffcosts,        /**< edge costs for cutoff */
   const SCIP_Real*      cutoffs,            /**< cutoff values for each incident edge (or NULL) */
   const SCIP_Real*      cutoffsrev,         /**< reverse edge cutoff values (or NULL if undirected or non-existent) */
   int                   vertex,             /**< the vertex */
   REDCOST*              redcostdata,        /**< reduced cost data for adaptation or NULL */
   SCIP_Bool*            success             /**< has node been pseudo-eliminated? */
   )
{
   DELPSEUDO delpseudo = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -1.0, -1, -1, UNKNOWN};

   assert(scip && success && g);
   assert(vertex >= 0 && vertex < g->knots);
   assert(g->grad[vertex] <= STP_DELPSEUDO_MAXGRAD);
   assert(!delPseudoIsEdgeDeletionMode(&delpseudo));

#ifndef NDEBUG
   {
      int sum = 0;
      for( int i = 1; i < STP_DELPSEUDO_MAXGRAD; i++ )
         sum += i;
      assert(sum == STP_DELPSEUDO_MAXNEDGES);
   }
#endif

   *success = TRUE;

   if( g->grad[vertex] <= 1 )
      return SCIP_OKAY;

   SCIP_CALL( delPseudoInit(scip, cutoffcosts, redcostdata, vertex, g, &delpseudo) );
   SCIP_CALL( delPseudoGetReplaceEdges(scip, g, cutoffs, cutoffsrev, &delpseudo, success) );

   /* enough spare edges? */
   if( (*success) )
   {
      SCIP_CALL( delPseudoDeleteVertex(scip, vertex, g, redcostdata, &delpseudo) );
   }

   delPseudoFreeData(scip, &delpseudo);

   return SCIP_OKAY;
}


/** checks whether pseudo-deletion is possible; maximum degree of STP_DELPSEUDO_MAXGRAD! */
SCIP_RETCODE graph_knot_delPseudoCheckIfPossible(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   const SCIP_Real*      cutoffcosts,        /**< edge costs for cutoff */
   const SCIP_Real*      cutoffs,            /**< cutoff values for each incident edge (or NULL) */
   const SCIP_Real*      cutoffsrev,         /**< reverse edge cutoff values (or NULL if undirected or non-existent) */
   int                   vertex,             /**< the vertex */
   SCIP_Bool*            isPossible          /**< can vertex pseudo-eliminated? */
   )
{
   DELPSEUDO delpseudo = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -1.0, -1, -1, UNKNOWN};

   assert(scip && isPossible && g);
   assert(vertex >= 0 && vertex < g->knots);
   assert(g->grad[vertex] <= STP_DELPSEUDO_MAXGRAD);
   assert(!delPseudoIsEdgeDeletionMode(&delpseudo));

   *isPossible = TRUE;

   if( g->grad[vertex] <= 3 )
      return SCIP_OKAY;

   SCIP_CALL( delPseudoInitForCheck(scip, g, cutoffcosts, vertex, &delpseudo) );
   SCIP_CALL( delPseudoCheckReplacement(scip, g, cutoffs, cutoffsrev, &delpseudo, isPossible) );

   delPseudoFreeDataForCheck(scip, &delpseudo);

   return SCIP_OKAY;
}



/** Pseudo deletes edge, i.e. reconnects tail of edge with neighbors of head; maximum degree of STP_DELPSEUDO_MAXGRAD!
 *  The orientation of the edge is important! */
SCIP_RETCODE graph_edge_delPseudo(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   const SCIP_Real*      cutoffcosts,        /**< edge costs for cutoff */
   const SCIP_Real*      cutoffs,            /**< cutoff values for each incident edge (or NULL) */
   const SCIP_Real*      cutoffsrev,         /**< reverse edge cutoff values (or NULL if undirected or non-existent) */
   int                   edge,               /**< the edge, mind the orientation! */
   SCIP_Real*            edgecosts_adapt,    /**< costs to adapt or NULL */
   SCIP_Bool*            success             /**< has node been pseudo-eliminated? */
   )
{
   DELPSEUDO delpseudo = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -1.0, -1, -1, edge };

   assert(scip && success && g);
   assert(graph_edge_isInRange(g, edge));
   assert(4 <= g->grad[g->head[edge]]);
   assert(g->grad[g->head[edge]] <= STP_DELPSEUDO_MAXGRAD);
   assert(delPseudoIsEdgeDeletionMode(&delpseudo));

   *success = TRUE;

   SCIP_CALL( delPseudoEdgeInit(scip, cutoffcosts, edgecosts_adapt, g, &delpseudo) );
   SCIP_CALL( delPseudoEdgeGetReplaceEdges(scip, g, cutoffs, cutoffsrev, &delpseudo, success) );

   /* enough spare edges? */
   if( (*success) )
   {
      SCIP_CALL( delPseudoEdgeDeleteEdge(scip, g, edgecosts_adapt, &delpseudo) );
   }

   delPseudoFreeData(scip, &delpseudo);

   return SCIP_OKAY;
}


/** Path replacement of edge; path is given by three oriented edges:
 *  edge_pathtail -> edge -> edge_pathhead.
 *  Middle edges is replaced. */
SCIP_RETCODE graph_edge_delPseudoPath(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   edge,               /**< the edge to be replaced, mind the orientation! */
   int                   edge_pathtail,      /**< tail edge of path */
   int                   edge_pathhead,      /**< head edge of path */
   SCIP_Real*            edgecosts_adapt     /**< costs to adapt or NULL */
   )
{
   assert(scip && g);
   assert(graph_edge_isInRange(g, edge));
   assert(graph_edge_isInRange(g, edge_pathtail));
   assert(graph_edge_isInRange(g, edge_pathhead));
   assert(g->head[edge_pathtail] == g->tail[edge]);
   assert(g->head[edge] == g->tail[edge_pathhead]);

   SCIP_CALL( delPseudoPath(scip, g, edge, edge_pathtail, edge_pathhead, edgecosts_adapt) );

   return SCIP_OKAY;
}
