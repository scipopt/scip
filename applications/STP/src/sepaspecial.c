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

/**@file   sepaspecial.c
 * @brief  Separator for Steiner tree problem contraints beyond flow-balance-directed-cut constraints
 * @author Daniel Rehfeldt
 *
 * This file includes some special separator routines beyond the flow-balance directed cut formulation constraints.
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "sepaspecial.h"
#include "probdata_stp.h"
#include "portab.h"
#include "stpvector.h"
#include "prop_stp.h"
#include "scip/cons_linear.h"

#define PCIMPLICATIONS_ALLOC_FACTOR 4


/** pseudo ancestor cliques */
struct pseudoancestor_cliques
{
   SCIP_Real*            ancestors_lpweights;/**< of size nancestors */
   int*                  halfedges_starts;   /**< CSR like starts for accessing ancestors */
   int*                  halfedges_ancestors;/**< CSR like ancestors per edge */
   int                   halfnedges;         /**< |A| / 2 */
   int                   nancestors;         /**< total number of pseudo ancestors */
   int                   nfails;             /**< number of failures */
};


/** implications between potential terminals */
struct prize_collecting_implications
{
   int*                  pcimplstart;        /**< start for each proper potential terminal */
   int*                  pcimplverts;        /**< all vertices */
   int                   pcimplnppterms;     /**< number of proper potential terminals used */
};


/** cuts for implications between non-terminals and terminals */
struct vertex_terminal_implications
{
   STP_Vectype(int)      impverts;           /**< implications vertices */
   STP_Vectype(int)      imparcs;            /**< implications arcs (from non-terminals to terminals) */
};


/*
 * local methods
 */

/** returns incoming flow for given node */
static inline
SCIP_Real get_inflow(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      xval,               /**< edge values */
   int                   vert                /**< the vertex */
)
{
   double insum = 0.0;

   for( int e = g->inpbeg[vert]; e >= 0; e = g->ieat[e] )
      insum += xval[e];

   return insum;
}



/** initializes pseudo-ancestor CSR like map from edges to ancestors */
static
SCIP_RETCODE pacliquesBuildMap(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   PACLIQUES*            pacliques           /**< clique data */
)
{
   int ancestors_size = 0;
   int* RESTRICT halfedges_starts;
   int* RESTRICT halfedges_ancestors;
   const int nedges = graph_get_nEdges(g);

   assert(scip && pacliques && g);
   assert(pacliques->ancestors_lpweights && pacliques->halfedges_starts);
   assert(pacliques->nancestors > 0);

   /* count */
   for( int e = 0; e < nedges; e += 2 )
   {
      const int nancestors = graph_edge_nPseudoAncestors(g, e);

      assert(nancestors >= 0);
      assert(ancestors_size < INT_MAX - nancestors);
      assert(graph_edge_nPseudoAncestors(g, e) == graph_edge_nPseudoAncestors(g, e + 1));
      ancestors_size += nancestors;
   }

   /* NOTE: ancestor size might be 0 if all edges with ancestor have been deleted */
   assert(ancestors_size >= 0);
   printf("pseudo-ancestor cliques total number of ancestors: %d \n", ancestors_size);

   if( ancestors_size > 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pacliques->halfedges_ancestors), ancestors_size) );
   }
   else
   {
      pacliques->halfedges_ancestors = NULL;
   }

   halfedges_starts = pacliques->halfedges_starts;
   halfedges_ancestors = pacliques->halfedges_ancestors;
   halfedges_starts[0] = 0;

   /* insert */
   for( int e = 0; e < nedges; e += 2 )
   {
      const int e_half = e / 2;
      const int nancestors = graph_edge_nPseudoAncestors(g, e);

      halfedges_starts[e_half + 1] = halfedges_starts[e_half] + nancestors;

      if( nancestors > 0 )
      {
         const int* const ancestors = graph_edge_getPseudoAncestors(g, e);
         assert(halfedges_ancestors);

         BMScopyMemoryArray(halfedges_ancestors + halfedges_starts[e_half], ancestors, nancestors);
      }
   }
   assert(halfedges_starts[pacliques->halfnedges] == ancestors_size);

   return SCIP_OKAY;
}


/*
 * Interface methods
 */


/** initializes pseudo-ancestor cliques */
SCIP_RETCODE sepaspecial_pacliquesInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   PACLIQUES**           pacliques           /**< clique data */
)
{
   PACLIQUES* pac;
   const int nedges = graph_get_nEdges(g);

   assert(scip && g);

   SCIP_CALL( SCIPallocMemory(scip, pacliques) );
   pac = *pacliques;

   pac->halfnedges = nedges / 2;
   assert(pac->halfnedges > 0);
   pac->nancestors = graph_getNpseudoAncestors(g);
   assert(pac->nancestors >= 0);

   pac->nfails = 0;

   pac->halfedges_starts = NULL;

   if( pac->nancestors > 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pac->ancestors_lpweights), pac->nancestors) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pac->halfedges_starts), pac->halfnedges + 1) );

      SCIP_CALL( pacliquesBuildMap(scip, g, pac) );
   }
   else
   {
      pac->ancestors_lpweights = NULL;
      pac->halfedges_starts = NULL;
      pac->halfedges_ancestors = NULL;
   }

   return SCIP_OKAY;
}


/** frees */
void sepaspecial_pacliquesFree(
   SCIP*                 scip,               /**< SCIP data structure */
   PACLIQUES**           pacliques           /**< clique data */
)
{
   SCIPfreeMemoryArrayNull(scip, &((*pacliques)->halfedges_ancestors));
   SCIPfreeMemoryArrayNull(scip, &((*pacliques)->halfedges_starts));
   SCIPfreeMemoryArrayNull(scip, &((*pacliques)->ancestors_lpweights));

   SCIPfreeMemoryArray(scip, pacliques);
}



/** separates */
SCIP_RETCODE sepaspecial_pacliquesSeparate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   PACLIQUES*            pacliques,          /**< clique data */
   int                   maxcuts,            /**< maximal number of cuts */
   int*                  ncuts               /**< pointer to store number of cuts */
   )
{
   SCIP_ROW** rows;
   GRAPH* g = SCIPprobdataGetGraph2(scip);
   SCIP_VAR** vars;
   const int* halfedges_starts;
   const int* halfedges_ancestors;
   SCIP_Real* RESTRICT ancestorweights;
   const SCIP_Real* xval;
   const int nedges = graph_get_nEdges(g);
   int cutscount;
   const int nancestors = pacliques->nancestors;
   SCIP_Bool hasViolations = FALSE;

   assert(pacliques && conshdlr && ncuts);
   assert(maxcuts > 0);

   /* nothing to separate? */
   if( nancestors == 0 )
      return SCIP_OKAY;

   if( pacliques->nfails > 1 )
      return SCIP_OKAY;

   cutscount = 0;

   halfedges_starts = pacliques->halfedges_starts;
   halfedges_ancestors = pacliques->halfedges_ancestors;
   ancestorweights = pacliques->ancestors_lpweights;
   vars = SCIPprobdataGetVars(scip);
   xval = SCIPprobdataGetXval(scip, NULL);

   assert(halfedges_starts && ancestorweights);
   assert(vars && xval);

   for( int i = 0; i < nancestors; i++ )
   {
      ancestorweights[i] = 0.0;
   }

   /* add LP ancestor weights up */
   for( int e = 0; e < nedges; e += 2 )
   {
      const SCIP_Real edgeweight = xval[e] + xval[e + 1];

      if( GT(edgeweight, 0.0) )
      {
         const int e_half = e / 2;
         for( int i = halfedges_starts[e_half]; i != halfedges_starts[e_half + 1]; i++ )
         {
            const int ancestor = halfedges_ancestors[i];
            assert(0 <= ancestor && ancestor < nancestors);

            ancestorweights[ancestor] += edgeweight;

            if( !hasViolations && GT(ancestorweights[ancestor], 1.0) )
               hasViolations = TRUE;
         }
      }
   }

   if( !hasViolations )
   {
      printf("no violations found... \n");
      pacliques->nfails++;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &rows, nancestors) );

   for( int i = 0; i < nancestors; i++ )
   {
      if( SCIPisFeasGT(scip, ancestorweights[i], 1.0) )
      {
         printf("violation for ancestor %d (%f > 1) \n", i, ancestorweights[i]);

         SCIP_CALL(SCIPcreateEmptyRowConshdlr(scip, &(rows[i]), conshdlr, "pa-clique", -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE));
         SCIP_CALL(SCIPcacheRowExtensions(scip, rows[i]));
      }
      else
      {
         rows[i] = NULL;
      }
   }

   /* separation loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      const int e_half = e / 2;

      for( int i = halfedges_starts[e_half]; i != halfedges_starts[e_half + 1]; i++ )
      {
         const int ancestor = halfedges_ancestors[i];
         assert(0 <= ancestor && ancestor < nancestors);

         if( rows[ancestor] )
         {
            SCIP_CALL(SCIPaddVarToRow(scip, rows[ancestor], vars[e], 1.0));
            SCIP_CALL(SCIPaddVarToRow(scip, rows[ancestor], vars[e + 1], 1.0));
         }
      }
   }

   for( int i = 0; i < nancestors; i++ )
   {
      if( rows[i] )
      {
         SCIP_CALL(SCIPflushRowExtensions(scip, rows[i]));

         if( SCIPisCutEfficacious(scip, NULL, rows[i]) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL(SCIPaddRow(scip, rows[i], FALSE, &infeasible));
            assert(!infeasible);

            printf("added conflict-cut for ancestor %d \n", i);
#if ADDCUTSTOPOOL
         /* add cut to pool */
            if( !infeasible )
               SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif
            if( *ncuts + cutscount++ >= maxcuts )
            {
               SCIP_CALL(SCIPreleaseRow(scip, &(rows[i])));
               break;
            }
         }

         SCIP_CALL(SCIPreleaseRow(scip, &(rows[i])));
      }
   }

   SCIPfreeBufferArray(scip, &rows);

   if( cutscount == 0 )
      pacliques->nfails++;

   return SCIP_OKAY;
}


/** initializes (R)PC implications */
SCIP_RETCODE sepaspecial_pcimplicationsInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   PCIMPLICATION**       pcimplications      /**< implication data */
)
{
   PCIMPLICATION* pcimps;
   int* start;
   int* verts;
   int* termmark;
   int* visitlist;
   SCIP_Real* dist;
   STP_Bool* visited;
   int nspares;
   int termscount;
   int nimplications;
   const int nnodes = graph_get_nNodes(g);
   const int nppterms = graph_pc_nProperPotentialTerms(g);
   const int maxnimplications = PCIMPLICATIONS_ALLOC_FACTOR * g->edges;
   const int slotsize = ((nppterms == 0) ? 0 : maxnimplications / nppterms);

   assert(scip && g);
   assert(graph_pc_isPcMw(g) && !g->extended);

   SCIP_CALL( SCIPallocMemory(scip, pcimplications) );
   pcimps = *pcimplications;
   pcimps->pcimplnppterms = 0;

   if( nppterms == 0 )
   {
      pcimps->pcimplstart = NULL;
      pcimps->pcimplverts = NULL;
      return SCIP_OKAY;
   }

   assert(slotsize >= 1 && slotsize * nppterms <= maxnimplications);

   SCIP_CALL(SCIPallocBufferArray(scip, &dist, nnodes));
   SCIP_CALL(SCIPallocBufferArray(scip, &visited, nnodes));
   SCIP_CALL(SCIPallocBufferArray(scip, &visitlist, nnodes));
   SCIP_CALL(SCIPallocBufferArray(scip, &termmark, nnodes));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(pcimps->pcimplstart), nppterms + 1));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(pcimps->pcimplverts), maxnimplications));

   start = pcimps->pcimplstart;
   verts = pcimps->pcimplverts;
   pcimps->pcimplnppterms = nppterms;

   for( int i = 0; i < nnodes; i++ )
   {
      visited[i] = FALSE;
      g->path_state[i] = UNKNOWN;
      dist[i] = FARAWAY;
   }

   graph_pc_termMarkProper(g, termmark);

   start[0] = 0;
   nspares = 0;
   termscount = 0;
   nimplications = 0;

   /* main loop: initialize implication lists */
   for( int i = 0; i < nnodes; i++ )
   {
      int nvisits;
      int nadded;

      if( !Is_term(g->term[i]) || graph_pc_knotIsFixedTerm(g, i) || graph_pc_termIsNonLeafTerm(g, i) )
         continue;

      assert(i != g->source);
      assert(g->path_heap && g->path_state);

      (void) graph_sdWalksConnected(scip, g, termmark, g->cost, NULL, i, 1000, dist, visitlist, &nvisits, visited, TRUE);

      assert(nvisits >= 1 && visitlist[0] == i);
      assert(nspares >= 0);

      for( int j = 1; j < MIN(nvisits, slotsize + nspares + 1); j++ )
      {
         const int vert = visitlist[j];
         assert(nimplications < maxnimplications);

         verts[nimplications++] = vert;
      }

      nadded = nimplications - start[termscount];
      assert(nadded >= 0);

      if( nadded > slotsize )
         nspares -= nadded - slotsize;
      else
         nspares += slotsize - nadded;

      assert(termscount < nppterms);
      start[++termscount] = nimplications;
   }
   assert(termscount == nppterms);

#ifndef WITH_UG
   printf("number of implications %d \n", nimplications);
#endif

   SCIPfreeBufferArray(scip, &termmark);
   SCIPfreeBufferArray(scip, &visitlist);
   SCIPfreeBufferArray(scip, &visited);
   SCIPfreeBufferArray(scip, &dist);

   return SCIP_OKAY;
}


/** frees (R)PC implications */
void sepaspecial_pcimplicationsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   PCIMPLICATION**       pcimplications      /**< implication data */
)
{
   SCIPfreeMemoryArrayNull(scip, &((*pcimplications)->pcimplstart));
   SCIPfreeMemoryArrayNull(scip, &((*pcimplications)->pcimplverts));

   SCIPfreeMemoryArray(scip, pcimplications);
}


/** separates PCSPG/MWCS implications */
SCIP_RETCODE sepaspecial_pcimplicationsSeparate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   PCIMPLICATION*        pcimplications,     /**< implications */
   int                   maxcuts,            /**< maximal number of cuts */
   int*                  ncuts               /**< pointer to store number of cuts */
   )
{
   GRAPH* g = SCIPprobdataGetGraph2(scip);
   SCIP_Real* nodeinflow;
   SCIP_VAR** vars = SCIPprobdataGetVars(scip);
   SCIP_ROW* row = NULL;
   const SCIP_Real* xval = SCIPprobdataGetXval(scip, NULL);
   const int* verts;
   const int* start;
   const int nnodes = g->knots;
   int cutscount;
   int ptermcount;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(pcimplications != NULL);
   assert(g != NULL);
   assert(xval != NULL);

   /* nothing to separate? */
   if( pcimplications->pcimplnppterms == 0 )
   {
      assert(graph_pc_nNonFixedTerms(g) == 0);
      return SCIP_OKAY;
   }

   assert(graph_pc_nNonFixedTerms(g) > 0);

   assert(pcimplications );
   verts = pcimplications->pcimplverts;
   start = pcimplications->pcimplstart;
   assert(verts && start);

   SCIP_CALL( SCIPallocBufferArray(scip, &nodeinflow, nnodes) );

   /* initialize node sums */
   for( int i = 0; i < nnodes; i++ )
      nodeinflow[i] = get_inflow(g, xval, i);

   cutscount = 0;
   ptermcount = 0;

   assert(g->extended);

   /* main separation loop */
   for( int i = 0; i < nnodes; i++ )
   {
      int maxnode;
      SCIP_Real maxflow;
      const SCIP_Real inflow = nodeinflow[i];

      if( !Is_pseudoTerm(g->term[i]) )
         continue;

      ptermcount++;

      if( SCIPisFeasGE(scip, inflow, 1.0) )
         continue;

      maxnode = -1;
      maxflow = 0.0;
      for( int j = start[ptermcount - 1]; j < start[ptermcount]; j++ )
      {
         const int vert = verts[j];
         if( SCIPisFeasGT(scip, nodeinflow[vert], inflow) && nodeinflow[vert] > maxflow )
         {
            maxnode = vert;
            maxflow = nodeinflow[vert];
         }
      }

      /* separate? */
      if( maxnode >= 0 )
      {
#ifdef XXX_XXX
         SCIP_CALL(SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "pcimplicate", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE));
         SCIP_CALL(SCIPcacheRowExtensions(scip, row));

         for( int e = g->inpbeg[maxnode]; e != EAT_LAST; e = g->ieat[e] )
            SCIP_CALL(SCIPaddVarToRow(scip, row, vars[e], 1.0));

         for( int e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
            SCIP_CALL(SCIPaddVarToRow(scip, row, vars[e], -1.0));
#else
         {
            const int twinterm = graph_pc_getTwinTerm(g, i);
            const int rootedge = graph_pc_getRoot2PtermEdge(g, twinterm);
            assert(rootedge >= 0);

            SCIP_CALL(SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "pcimplicate", -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE));
            SCIP_CALL(SCIPcacheRowExtensions(scip, row));

            for( int e = g->inpbeg[maxnode]; e != EAT_LAST; e = g->ieat[e] )
               SCIP_CALL(SCIPaddVarToRow(scip, row, vars[e], 1.0));

            SCIP_CALL(SCIPaddVarToRow(scip, row, vars[rootedge], 1.0));
         }
#endif

         SCIP_CALL(SCIPflushRowExtensions(scip, row));

         if( SCIPisCutEfficacious(scip, NULL, row) )
         {
            SCIP_Bool infeasible;
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
#if ADDCUTSTOPOOL
         /* add cut to pool */
            if( !infeasible )
               SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

            if( *ncuts + cutscount++ >= maxcuts )
            {
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
               break;
            }
         }

         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }
   }
   assert((*ncuts + cutscount > maxcuts) || ptermcount == graph_pc_nProperPotentialTerms(g));

   *ncuts += cutscount;
   SCIPdebugMessage("PcImplication Separator: %d Inequalities added\n", cutscount);

   SCIPfreeBufferArray(scip, &nodeinflow);

   return SCIP_OKAY;
}

/** return number of proper potential terminals */
int sepaspecial_pcimplicationsGetNstarts(
   const PCIMPLICATION*  pcimp               /**< implications */
   )
{
   assert(pcimp);

   return pcimp->pcimplnppterms;
}


/** return CSR like starts */
const int* sepaspecial_pcimplicationsGetStarts(
   const PCIMPLICATION*  pcimp               /**< implications */
   )
{
   assert(pcimp);

   return pcimp->pcimplstart;
}


/** return CSR like vertices */
const int* sepaspecial_pcimplicationsGetVerts(
   const PCIMPLICATION*  pcimp               /**< implications */
   )
{
   assert(pcimp);

   return pcimp->pcimplverts;
}


/** initializes implications */
SCIP_RETCODE sepaspecial_vtimplicationsInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   VTIMPLICATION**       vtimplications      /**< implication data */
)
{
   VTIMPLICATION* vtimps;
   const int nnodes = graph_get_nNodes(g);

   assert(scip);

   SCIP_CALL( SCIPallocMemory(scip, vtimplications) );
   vtimps = *vtimplications;
   vtimps->impverts = NULL;
   vtimps->imparcs = NULL;

   // todo reactivate if it works well
   if( !graph_typeIsSpgLike(g) )
   {
      return SCIP_OKAY;
   }

   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) )
      {
         SCIP_Real min = FARAWAY;
         int minedge = -1;

         for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         {
            if( GE(g->cost[e], BLOCKED) )
            {
               minedge = -2;
               break;
            }

            if( LT(g->cost[e], min) )
            {
               min = g->cost[e];
               minedge = e;
            }
         }

         assert(minedge != -1);

         if( minedge == -2 )
            continue;

         assert(!Is_term(g->term[g->head[minedge]]));

         StpVecPushBack(scip, vtimps->impverts, g->head[minedge]);
         StpVecPushBack(scip, vtimps->imparcs, flipedge(minedge));
      }
   }

   return SCIP_OKAY;
}


/** frees implications */
void sepaspecial_vtimplicationsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   VTIMPLICATION**       vtimplications      /**< implication data */
)
{
   StpVecFree(scip, (*vtimplications)->imparcs);
   StpVecFree(scip, (*vtimplications)->impverts);

   SCIPfreeMemoryArray(scip, vtimplications);
}


/** separates implications */
SCIP_RETCODE sepaspecial_vtimplicationsSeparate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   VTIMPLICATION*        vtimplications,     /**< implication data */
   int                   maxcuts,            /**< maximal number of cuts */
   int*                  ncuts               /**< pointer to store number of cuts */
   )
{
   GRAPH* g = SCIPprobdataGetGraph2(scip);
   SCIP_VAR** vars = SCIPprobdataGetVars(scip);
   const SCIP_Real* xval = SCIPprobdataGetXval(scip, NULL);
   int cutscount = 0;
   int nimplications;

   assert(scip && conshdlr && vtimplications && ncuts);
   assert(maxcuts > 0);
   assert(StpVecGetSize(vtimplications->imparcs) == StpVecGetSize(vtimplications->impverts));

   nimplications = StpVecGetSize(vtimplications->imparcs);

   /* nothing to separate? */
   if( nimplications == 0)
   {
      return SCIP_OKAY;
   }

   for( int i = 0; i < nimplications; i++ )
   {
      const int impvert = vtimplications->impverts[i];
      const int imparc = vtimplications->imparcs[i];
      const int imparc_rev = flipedge(imparc);
      const SCIP_Real inflow = get_inflow(g, xval, impvert);

      assert(graph_edge_isInRange(g, imparc));
      assert(g->tail[imparc] == impvert);

      if( SCIPisFeasLT(scip, xval[imparc] + xval[imparc_rev], inflow) )
      {
         SCIP_ROW* row = NULL;

         SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "vtimplicate", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

         printf("found implication cut: %f < %f \n", xval[imparc] + xval[imparc_rev], inflow);

         for( int e = g->inpbeg[impvert]; e != EAT_LAST; e = g->ieat[e] )
         {
            if( e != imparc_rev )
            {
               SCIP_CALL(SCIPaddVarToRow(scip, row, vars[e], 1.0));
            }
            assert(e != imparc);
         }

         SCIP_CALL(SCIPaddVarToRow(scip, row, vars[imparc], -1.0));
         SCIP_CALL(SCIPflushRowExtensions(scip, row));

         if( SCIPisCutEfficacious(scip, NULL, row) )
         {
            SCIP_Bool infeasible;
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
            printf("...added implication cut! \n");

#if ADDCUTSTOPOOL
         /* add cut to pool */
            if( !infeasible )
               SCIP_CALL( SCIPaddPoolCut(scip, row) );
#endif

            if( *ncuts + cutscount++ >= maxcuts )
            {
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
               break;
            }
         }

         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }
   }

   *ncuts += cutscount;
   SCIPdebugMessage("Vertex-Terminal Implication Separator: %d Inequalities added\n", cutscount);

   return SCIP_OKAY;
}
