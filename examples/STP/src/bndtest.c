/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Function                                                      */
/*   File....: sdtest.c                                                      */
/*   Name....: Special Distance Test                                         */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "grph.h"
#include "portab.h"
#include "scip/scip.h"
#include "probdata_stp.h"

void get_close_terms(
      PATH**   path,
      double*  closetermsdist,
      int*     closetermshops,
      int*     closeterms,
      int*     vregion,
      int*     terms,
      int      termcount,
      int      node
      )
{
   int j;

   assert(closetermsdist   != NULL);
   assert(closeterms       != NULL);
   assert(vregion          != NULL);
   assert(terms            != NULL);

   closeterms[0] = vregion[node];
   closetermsdist[0] = path[vregion[node]][node].dist;
   closetermshops[0] = path[vregion[node]][node].hops;
   for( j = 1; j < 3; j++ )
   {
      closeterms[j] = -1;
      closetermsdist[j] = FARAWAY;
      closetermshops[j] = -1;
   }

   for( j = 0; j < termcount; j++ )
   {
      //printf("(%d, %f), ", terms[j], path[terms[j]][node].dist);
      if( j == vregion[node] )
         continue;

      if( LT(path[terms[j]][node].dist, closetermsdist[1]) ||
            (EQ(path[terms[j]][node].dist, closetermsdist[1]) && LT(path[terms[j]][node].hops, closetermshops[1])))
      {
         closeterms[2] = closeterms[1];
         closetermsdist[2] = closetermsdist[1];
         closetermshops[2] = closetermshops[1];

         closeterms[1] = terms[j];
         closetermsdist[1] = path[terms[j]][node].dist;
         closetermshops[1] = path[terms[j]][node].hops;
      }
      else if( LT(path[terms[j]][node].dist, closetermsdist[2]) ||
            (EQ(path[terms[j]][node].dist, closetermsdist[2]) && LT(path[terms[j]][node].hops, closetermshops[2])))
      {
         closeterms[2] = terms[j];
         closetermsdist[2] = path[terms[j]][node].dist;
         closetermshops[2] = path[terms[j]][node].hops;
      }
   }
   //printf("\n");
}

double compute_node_lb(
      double*  radius,
      double*  closetermsdist,
      int*     closetermshops,
      int*     closeterms,
      int*     radiushops,
      int      termcount,
      int      nodegrad,
      int      source,
      int      probtype,
      int*     hopsbound
      )
{
   double lowerbound;
   int j;

   assert(radius           != NULL);
   assert(closetermsdist   != NULL);
   assert(closetermshops   != NULL);
   assert(closeterms       != NULL);

   lowerbound = closetermsdist[0] + closetermsdist[1];
   (*hopsbound) = closetermshops[0] + closetermshops[1];
   for( j = 0; j < termcount - 2; j++ )
   {
      lowerbound += radius[j];
      (*hopsbound) += radiushops[j];
   }

   /* For the hop constrained problems we are not able to obtain a tighter bound given when the grad >= 3 condition */
   if( probtype == STP_HOP_CONS )
      return lowerbound;

#if 1
   if( nodegrad >= 3 )
   {
      lowerbound += closetermsdist[2];
      (*hopsbound) += closetermshops[2];
      lowerbound -= radius[termcount - 4];
      (*hopsbound) -= radiushops[termcount - 4];
   }
#endif

   return lowerbound;
}



SCIP_RETCODE bound_reduction(
      SCIP* scip,
      double bestbound,
      int* elimins
      )
{
   GRAPH* g;
   SCIP_PROBDATA* probdata;
   SCIP_VAR** vars;
   PATH**  path;
   PATH*   pathfromterm;
   double* distance;
   double* radius;
   double** termdist;
   int*    radiushops;
   int*    vregion;
   int*    heap;
   int*    state;
   int*    pred;
   int*    terms;
   int     termcount;
   int     i;
   int     j;
   int     k;
   int     e;
   int     closeterms[3] = {-1, -1, -1};
   double  closetermsdist[3] = {FARAWAY, FARAWAY, FARAWAY};
   int     closetermshops[3] = {-1, -1, -1};
   double  tempcost;
   double  lowerbound;
   int     hopsbound = 0;
   SCIP_Bool allfixed;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;

   SCIPdebugMessage("BoundBased-Reduction: ");
   fflush(stdout);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   g = SCIPprobdataGetGraph(probdata);
   assert(g != NULL);

   (*elimins) = 0;
   if( LT(g->greatestlb, SCIPgetDualbound(scip)))
      return SCIP_OKAY;

   /* getting the problem variables - the edge variables */
   vars = SCIPprobdataGetVars(scip);
   assert(vars != NULL);
   assert(vars[0] != NULL);

   /*
     graph_show(g);
   */
   path = malloc((size_t)g->knots * sizeof(PATH*));

   assert(path != NULL);

   pathfromterm = malloc((size_t)g->knots * sizeof(PATH));

   assert(pathfromterm != NULL);

   distance = malloc((size_t)g->knots * sizeof(double));
   radius = malloc((size_t)g->knots * sizeof(double));
   termdist = malloc((size_t)g->knots * sizeof(double*));

   assert(distance != NULL);
   assert(radius != NULL);
   assert(termdist != NULL);

   radiushops = malloc((size_t)g->knots * sizeof(int));
   vregion = malloc((size_t)g->knots * sizeof(int));

   assert(radiushops != NULL);
   assert(vregion != NULL);

   heap  = malloc((size_t)g->knots * sizeof(int));
   state = malloc((size_t)g->knots * sizeof(int));

   assert(heap != NULL);
   assert(state != NULL);

   pred = malloc((size_t)g->edges * sizeof(int));
   terms = malloc((size_t)g->terms * sizeof(int));

   assert(g->source[0] >= 0);

   termcount = 0;
   for(i = 0; i < g->knots; i++)
   {
      if( Is_term(g->term[i]) )
      {
         terms[termcount] = i;
         termcount++;
      }
      g->mark[i] = (g->grad[i] > 0);
      path[i] = NULL;
      termdist[i] = NULL;
      radiushops[i] = 0;
   }

   /* computing the voronoi regions inward to a node */
   if( g->stp_type == STP_HOP_CONS )
      voronoi_hop(g, g->cost, distance, radius, pathfromterm, vregion, heap, state, pred, radiushops);
   else
      voronoi_term(g, g->cost, distance, radius, termdist, pathfromterm, vregion, heap, state, pred, 1);

   /* sorting the radius values */
   SCIPsortRealInt(radius, radiushops, g->knots);

   /* computing the shortest paths from each terminal to every other node */
   if( g->stp_type == STP_HOP_CONS )
      calculate_distances(g, path, g->cost, BSP_MODE);
   else
      calculate_distances(g, path, g->cost, FSP_MODE);

   /* resetting the greatest lower bound */
   g->greatestlb = 0;

   for( i = 0; i < g->knots; i++ )
   {
      if( Is_term(g->term[i]) || (g->stp_type == STP_HOP_CONS && g->source[0] == vregion[i]) )
         continue;

      /* don't process the knot if it has already been eliminated */
      if( g->elimknots[i] < 0 )
         continue;

      /* knot is not reachable as a result of the variable fixings */
      if( vregion[i] < 0 )
         continue;

      //for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
         //printf("e: %d(%d), head: %d(%d,%d), tail: %d(%d,%d)\n", e, g->elimedges[e], g->head[e], g->term[g->head[e]],
               //g->elimknots[g->head[e]], g->tail[e], g->term[g->tail[e]], g->elimknots[g->tail[e]]);
      get_close_terms(path, closetermsdist, closetermshops, closeterms, vregion, terms, termcount, i);


      /* computing the lower bound for node i */
      lowerbound = compute_node_lb(radius, closetermsdist, closetermshops, closeterms, radiushops, termcount, g->grad[i],
            g->source[0], g->stp_type, &hopsbound);

      //printf("termdist: %f, %f, radius: %f\n", closetermsdist[0], closetermsdist[1],
            //lowerbound - (closetermsdist[0] + closetermsdist[1]));

      /* fixing variables in and out of the target node as a result of the bound violation */
      if( GT(lowerbound, bestbound) )
         printf("Knot (%d, %d) - lowerbound: %f, bestbound: %f, hopsbound: %d, hopslimit: %d\n", i, g->grad[i],
               lowerbound, bestbound, hopsbound, g->hoplimit);
      if( GT(lowerbound, bestbound) && (g->stp_type != STP_HOP_CONS || hopsbound < g->hoplimit) )
      {
         printf("=========================================================================================>>>\n");
         allfixed = TRUE;
         for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
         {
            /* strengthen upper bound */
            SCIP_CALL(SCIPtightenVarUbGlobal(scip, vars[j], 0.0, FALSE, &infeasible, &fixed));
            if( fixed || g->elimedges[j] == -1 )
               g->elimedges[j] = -1;
            else
               allfixed = FALSE;

            /* strengthen upper bound */
            SCIP_CALL(SCIPtightenVarUbGlobal(scip, vars[Edge_anti(j)], 0.0, FALSE, &infeasible, &fixed));
            if( fixed || g->elimedges[Edge_anti(j)] == -1 )
               g->elimedges[Edge_anti(j)] = -1;
            else
               allfixed = FALSE;

            printf("Fixed: %d, Infeasible: %d\n", fixed, infeasible);
         }

         if( allfixed )
         {
            g->elimknots[i] = -1;
            (*elimins)++;
         }

         /* computing the voronoi regions inward to a node */
         if( g->stp_type == STP_HOP_CONS )
            voronoi_hop(g, g->cost, distance, radius, pathfromterm, vregion, heap, state, pred, radiushops);
         else
            voronoi_term(g, g->cost, distance, radius, termdist, pathfromterm, vregion, heap, state, pred, 1);

         /* sorting the radius values */
         SCIPsortRealInt(radius, radiushops, g->knots);

         /* computing the shortest paths from each terminal to every other node */
         //if( g->stp_type == STP_HOP_CONS )
            //calculate_distances(g, path, g->cost, BSP_MODE);
         //else
            //calculate_distances(g, path, g->cost, FSP_MODE);
      }

      /* updating the greatest lower bound */
      if( LT(g->greatestlb, lowerbound) )
         g->greatestlb = lowerbound;

      /* don't process the knot if it has already been eliminated */
      if( g->elimknots[i] < 0 )
         continue;

      /* knot is not reachable as a result of the variable fixings */
      if( vregion[i] < 0 )
         continue;

      if( g->stp_type == STP_HOP_CONS )
         continue;
#if 1
      /* computing the lowerbound for the inclusion of a single edge */
      /* Resetting the lowerbound */
      lowerbound -= closetermsdist[1];
      if( g->grad[i] >= 3 )
      {
         lowerbound -= closetermsdist[2];
         lowerbound += radius[termcount - 4];
      }

      for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
      {
         if( g->elimedges[j] < 0 )
            continue;

         k = g->tail[j];
         tempcost = g->cost[j] + path[vregion[k]][k].dist;
         lowerbound += tempcost;

         //printf("Edge - lowerbound: %f, bestbound: %f\n", lowerbound, bestbound);
         if( GT(lowerbound, bestbound) )
         {
            /* strengthen upper bound */
            SCIP_CALL(SCIPtightenVarUbGlobal(scip, vars[j], 0.0, FALSE, &infeasible, &fixed));
            if( fixed )
               g->elimedges[j] = -1;

            /* strengthen upper bound */
            SCIP_CALL(SCIPtightenVarUbGlobal(scip, vars[Edge_anti(j)], 0.0, FALSE, &infeasible, &fixed));
            if( fixed )
               g->elimedges[Edge_anti(j)] = -1;

            (*elimins)++;

            /* computing the voronoi regions inward to a node */
            voronoi_term(g, g->cost, distance, radius, termdist, pathfromterm, vregion, heap, state, pred, 1);

            /* sorting the radius values */
            SCIPsortReal(radius, g->knots);

            /* computing the shortest paths from each terminal to every other node */
            //calculate_distances(g, path, g->cost, FSP_MODE);
            break;
         }

         /* updating the greatest lower bound */
         if( LT(g->greatestlb, lowerbound) )
            g->greatestlb = lowerbound;

         lowerbound -= tempcost;
      }
#endif


   }

   free(pred);
   free(vregion);
   free(radiushops);
   free(termdist);
   free(radius);
   free(distance);
   free(pathfromterm);
   free(path);

   assert(graph_valid(g));
   SCIPdebugMessage("bound_reduction: %d Knots deleted\n", (*elimins));
   printf("bound_reduction: %d Knots deleted\n", (*elimins));

   return SCIP_OKAY;
}
