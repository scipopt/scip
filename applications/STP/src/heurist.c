#ident "@(#) $Id: heurist.c,v 1.10 1999/12/01 12:21:27 bzfkocht Exp $"
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Function                                                      */
/*   File....: heurist.c                                                     */
/*   Name....: Steiner Tree Packing Heuristic                                */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "jack3.h"

/* Die Heuristic stoert sich nicht dran, wenn sie einzelne Wege nicht
 * routen kann, sondern erklaert das jeweilige Netz einfach fuer fertig.
 * Die Loesung muss also ueberprueft werden.
 */
static void do_heuristic(
   const GRAPH*  g,
   int           layer,
   int*          result,
   int           start,
   char*         connected,
   const double* cost,
   PATH**        path)
{
   PATH*  mst;
   int*   cluster;
   int    csize = 0;
   int    k;
   int    e;
   int    count;
   double min;
   int    i;
   int    j;
   int    old;
   int    new;

   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(path      != NULL);
   assert((layer > -1) && (layer < g->layers));

   printf("Heuristic: Start=%5d ", start);
   fflush(stdout);

   cluster = malloc((size_t)g->knots * sizeof(int));

   assert(cluster != NULL);

   for(i = 0; i < g->knots; i++)
   {
      g->mark[i]   = (g->grad[i] > 0);
      connected[i] = FALSE;
   }
   connected[start] = TRUE;
   cluster[csize++] = start;

   /* CONSTCOND */
   for(;;)
   {
      /* Suche das Terminal, das am dichtesten dran ist
       */
      min = FARAWAY;
      old = -1;
      new = -1;

      for(i = 0; i < g->knots; i++)
      {
         if (g->grad[i] == 0)
            continue;

         if (g->term[i] != layer)
            continue;

         if (connected[i])
            continue;

         /* Jetzt brauchen wir die Entfernungen.
          */
         if (path[i] == NULL)
         {
            path[i] = malloc((size_t)g->knots * sizeof(PATH));

            assert(path[i] != NULL);

            /* ! Ob das die guenstiges Richtung ist die Wege zu berechnen, wenn
             * ! die Kosten fuer Hin und Rueckweg unterschiedlich sind ist doch
             * ! sehr fraglich.
             * ! Koennte aber sein, weil wir die Kosten unten umgedreht haben.
             */
            graph_path_exec(g, FSP_MODE, i, cost, path[i]);
         }
         for(k = 0; k < csize; k++)
         {
            j = cluster[k];

            assert(i != j);
            assert(connected[j]);

            if (LT(path[i][j].dist, min))
            {
               min = path[i][j].dist;
               new = i;
               old = j;
            }
         }
      }
      /* Nichts mehr gefunden, also fertig
       */
      if (new == -1)
         break;

      /* Weg setzten
       */
      assert((old > -1) && (new > -1));
      assert(path[new] != NULL);
      assert(path[new][old].dist < FARAWAY);
      assert(g->term[new] == layer);
      assert(!connected[new]);
      assert(connected[old]);

      fputc('R', stdout);
      fflush(stdout);

      /*    printf("Connecting Knot %d-%d dist=%d\n", new, old, path[new][old].dist);
       */
      /* Gegen den Strom schwimmend alles markieren
       */
      k = old;

      while(k != new)
      {
         e = path[new][k].edge;
         k = g->tail[e];

         if (!connected[k])
         {
            connected[k] = TRUE;
            cluster[csize++] = k;
         }
      }
   }
   free(cluster);

   fputc('M', stdout);
   fflush(stdout);

   mst = malloc((size_t)g->knots * sizeof(PATH));

   assert(mst != NULL);

   /* MST berechnen
    */
   for(i = 0; i < g->knots; i++)
      g->mark[i] = connected[i];

   assert(g->source[layer] >= 0);
   assert(g->source[layer] <  g->knots);

   graph_path_exec(g, MST_MODE, g->source[layer], g->cost, mst);

   disp_path(g, mst);

   for(i = 0; i < g->knots; i++)
   {
      if (connected[i] && (mst[i].edge != -1))
      {
         assert(g->head[mst[i].edge] == i);
         assert(result[mst[i].edge] == -1);

         result[mst[i].edge] = layer;
      }
   }

   /* Baum beschneiden
    */
   do
   {
      fputc('C', stdout);
      fflush(stdout);

      count = 0;

      for(i = 0; i < g->knots; i++)
      {
         if (!g->mark[i])
            continue;

         if (g->term[i] == layer)
            continue;

         for(j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j])
            if (result[j] == layer)
               break;

         if (j == EAT_LAST)
         {
            /* Es muss genau eine eingehende Kante geben
             */
            for(j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j])
            {
               if (result[j] == layer)
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
   while(count > 0);

   free(mst);
}

static void do_layer(
   const GRAPH*  g,
   int           layer,
   int*          best_result,
   int           runs,
   const double* cost)
{
   static int best  = -1;

   PATH** path      = malloc((size_t)g->knots * sizeof(PATH *));
   char*  connected = malloc((size_t)g->knots * sizeof(char));
   int*   result    = malloc((size_t)g->edges * sizeof(int));
   int*   start     = malloc((size_t)g->knots * sizeof(int));
   int    i;
   int    k;
   int    t;
   double obj;
   double min       = FARAWAY;

   assert(g         != NULL);
   assert(result    != NULL);
   assert(path      != NULL);
   assert(connected != NULL);
   assert((layer > -1) && (layer < g->layers));

   /* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    * Patch um die heuristic nach einem restruct starten zu koennen
    * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    */
   if (best >= g->knots)
      best = -1;

   for(i = 0; i < g->knots; i++)
      path[i] = NULL;

   if (g->layers > 1)
      do_heuristic(g, layer, best_result, g->source[layer], connected, cost, path);
   else
   {
      runs = (runs > g->knots) ? g->knots : runs;
      best = (best <  0)       ? g->source[layer] : best;

      for(i = 0; i < g->knots; i++)
      {
         assert(g->grad[i] > 0);

         start[i] = i;
      }
      /* Wenn wir sowieso alle probieren, brauchen wir da nachfolgende
       * nicht machen.
       */
      if (runs < g->knots)
      {
         for(i = 0; i < runs; i++)
         {
            k        = rand() % g->knots;
            t        = start[i];
            start[i] = start[k];
            start[k] = t;
         }
         for(i = 0; i < runs; i++)
            if (start[i] == best)
               break;

         if (i == runs)
            start[0] = best;
      }

      for(i = 0; i < runs; i++)
      {
         /* Inkorrekt bei layers > 1 !
          */
         assert(g->layers == 1);

         for(k = 0; k < g->edges; k++)
            result[k] = -1;

         do_heuristic(g, layer, result, start[i], connected, cost, path);

         obj = 0.0;

         /* Hier wird mit anderem Mass gemessen als bei do_heuristic()
          */
         for(k = 0; k < g->edges; k++)
            obj += (result[k] > -1) ? g->cost[k] : 0.0;

         printf(" Obj=%.12e\n", obj);

         if (LT(obj, min))
         {
            min = obj;

            for(k = 0; k < g->edges; k++)
               best_result[k] = result[k];

            best = start[i];
         }
      }
   }
   printf("Freeing Memory\n");

   for(i = 0; i < g->knots; i++)
   {
      if (path[i] != NULL)
      {
         assert(g->term[i] == layer);
         free(path[i]);
      }
   }
   free(result);
   free(path);
   free(connected);
   free(start);
}

int heuristic(
   BACY*   b,
   double* xval,
   int     runs)
{
   GRAPH*  g;
   double* cost;
   int*    result;
   double* nval;
   double  pobj;
   int     retval  = FALSE;
   int     layer;
   int     i;

   assert(b       != NULL);
   assert(runs    >= 0);

   printf("Heuristic Start\n");

   g       = b->problem;
   cost    = malloc((size_t)g->edges * sizeof(double));
   result  = malloc((size_t)g->edges * sizeof(int));
   nval    = malloc((size_t)b->vars * sizeof(double));

   assert(g      != NULL);
   assert(cost   != NULL);
   assert(result != NULL);
   assert(nval   != NULL);

   for(i = 0; i < g->edges; i++)
      result[i] = -1;

   for(layer = 0; layer < g->layers; layer++)
   {
      if (xval == NULL)
      {
         for(i = 0; i < g->edges; i++)
            cost[i] = g->cost[i];
      }
      else
      {
         /* Kosten vertauschen
          */
         for(i = 0; i < g->edges; i += 2)
         {
            /* Permanent auf 0 fixierte sind keine gute Wahl
             */
#if 0
            if (b->state[layer * g->edges + i] | b->state[layer * g->edges + i + 1] | FOREVER)
            {
               /* Geht nicht ueber edge_hide(), weil die Fixierung ja
                * nur fuer diesen Layer ist.
                */
               cost[i]     = FARAWAY / 2;
               cost[i + 1] = FARAWAY / 2;
            }
            else
#endif
            {
               cost[i]     = ((1.0 - xval[layer * g->edges + i + 1]) * g->cost[i + 1]);
               cost[i + 1] = ((1.0 - xval[layer * g->edges + i    ]) * g->cost[i]);
            }
         }
      }
      /* Koennen wie das Netz verbinden ?
       */
      do_layer(g, layer, result, runs, cost);

      /* Weg besetzen
       */
      if (g->layers > 1)
         for(i = 0; i < g->edges; i += 2)
            if ((result[i] == layer) || (result[i + 1] == layer))
               graph_edge_hide(g, i);
   }
   if (g->layers > 1)
      graph_uncover(g);

   for(i = 0; i < b->vars; i++)
      nval[i] = (result[i % g->edges] == (i / g->edges)) ? 1.0 : 0.0;

   printf("Validation\n");

   if (validate(g, nval))
   {
      printf("Is Valid\n");

      pobj = 0.0;

      for(i = 0; i < b->vars; i++)
         pobj += g->cost[i % g->edges] * nval[i];

      printf("pobj=%.12e\n", pobj);

      if (pobj < b->pobj - EPSILON)
      {
         memcpy(b->xval, nval, (unsigned int)b->vars * sizeof(double));
         b->pobj = pobj;
         retval  = TRUE;

         printf("Heuristic: New Solution Value=%.12e\n", pobj);
      }
   }
   free(nval);
   free(result);
   free(cost);

   return(retval);
}
