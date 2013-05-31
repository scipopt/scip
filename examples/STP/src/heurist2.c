#ident "@(#) $Id: heurist2.c,v 1.6 1999/12/01 12:21:27 bzfkocht Exp $"
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Function                                                      */
/*   File....: heurist2.c                                                    */
/*   Name....: Steiner Tree Heuristic (Rayward-Smith)                        */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Der nachfolgend implementierte Algorithmus ist dem Artikel:
 *
 *       The computation of nearly minimal Steiner trees in graphs
 *
 * von V. J. Rayward-Smith entnommen.
 */
#define HEURIST2_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "jack3.h"

typedef struct component
{
   int    start;
   int    knot;
   double dist;
} COMP;

#define OUTPUT

#define COMP_MINI     0
#define COMP_LAST    -1
#define COMP_NONE    -2

static int comp_cmp1(
   const void* a,
   const void* b)
{
   const double da = ((const COMP*)a)->dist;
   const double db = ((const COMP*)b)->dist;

   if (EQ(da, db))
      return(0);

   return((da < db) ? -1 : 1);
}

static int comp_cmp2(
   const void* a,
   const void* b)
{
   const double da = ((const COMP*)a)->dist;
   const double db = ((const COMP*)b)->dist;

   if (EQ(da, db))
      return(0);

   return((db < da) ? -1 : 1);
}

inline static double evaluate_f(
   int         r,
   const COMP* tree)
{
   int    i;
   double ret = 0;

   assert(r    >= 2);
   assert(tree != NULL);

   for(i = 0; i < r; i++)
      ret += tree[i].dist;

   return(ret / (r - 1));
}

static double distance(
   const GRAPH* g,
   const PATH*  path,
   int          k,
   int          trees,
   COMP*        tree,
   const int*   next)
{
   int    i;
   int    j;
   double fval;

   assert(g       != NULL);
   assert(path    != NULL);
   assert(trees   >= 2);
   assert(trees   <= g->terms);
   assert(tree    != NULL);
   assert(next    != NULL);

   for(i = 0; i < trees; i++)
   {
      tree[i].dist = FARAWAY;

      for(j = tree[i].start; j >= COMP_MINI; j = next[j])
         if (GT(tree[i].dist, path[j].dist))
            tree[i].dist = path[j].dist;
   }
   qsort(tree, (size_t)trees, sizeof(COMP), comp_cmp1);

   i = 2;

   do
   {
      fval = evaluate_f(i, tree);
      i++;
   }
   while((i <= trees) && LT(tree[i].dist, fval));

   return(fval);
}

static int find_best(
   const GRAPH* g,
   PATH*        path,
   int          trees,
   COMP*        tree,
   const int*   next)
{
   double fval;
   double min_fval;
   int    min_knot;
   int    i;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(trees  <= g->terms);
   assert(tree   != NULL);
   assert(next   != NULL);

   min_fval = FARAWAY;
   min_knot = -1;

   for(i = 0; i < g->knots; i++)
   {
      graph_path_exec(g, FSP_MODE, i, g->cost, path);

      fval = distance(g, path, i, trees, tree, next);

      if (LT(fval, min_fval))
      {
         min_fval = fval;
         min_knot = i;
      }
   }
   assert(LT(min_fval, FARAWAY));
   assert(min_knot > -1);
   assert(min_knot < g->knots);

   return(min_knot);
}

static int connect2(
   const GRAPH* g,
   PATH*        path,
   int          knot,
   int          trees,
   COMP*        tree,
   int*         next)
{
   int i;
   int j;
   int e;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(knot   >= 0);
   assert(knot   <  g->knots);
   assert(trees  >= 2);
   assert(trees  <= g->terms);
   assert(tree   != NULL);
   assert(next   != NULL);

   graph_path_exec(g, FSP_MODE, knot, g->cost, path);

   for(i = 0; i < trees; i++)
   {
      tree[i].dist = FARAWAY;
      tree[i].knot = COMP_NONE;

      for(j = tree[i].start; j >= COMP_MINI; j = next[j])
      {
         if (GT(tree[i].dist, path[j].dist))
         {
            tree[i].dist = path[j].dist;
            tree[i].knot = j;
         }
      }
      assert(LT(tree[i].dist, FARAWAY));
      assert(tree[i].knot       >= COMP_MINI);
      assert(next[tree[i].knot] != COMP_NONE);
   }
   qsort(tree, (size_t)trees, sizeof(COMP), comp_cmp2);

   /* Beiden Baeume zusammen haengen
    */
   trees--;

   for(i = tree[trees - 1].start; next[i] >= COMP_MINI; i = next[i]);

   assert(next[i] == COMP_LAST);

   next[i] = tree[trees].start;

   /* Den Weg zwischen den letzten beiden Komponenten via 'knot' ermitteln
    */
   for(j = trees; j > trees - 2; j--)
   {
      i = tree[j].knot;

      assert(next[i] != COMP_NONE);

      while(i != knot)
      {
         e = path[i].edge;

         assert(e >= 0);
         assert(e < g->edges);

         i = g->tail[e];

         assert(i == g->tail[e]);
         assert(i >= 0);
         assert(i <  g->knots);

         if (next[i] == COMP_NONE)
         {
            assert(next[i] == COMP_NONE);

            /* Knoten zur Komponente hinzufuegen
             */
            next[i]               = tree[trees - 1].start;
            tree[trees - 1].start = i;
         }
      }
   }
   return(trees);
}

int heuristic2(
   BACY*   b,
   double* xval)
{
   /* xval wird momentan noch nicht gebraucht.
    */
   GRAPH*  g;
   PATH*   path;
   int     i;
   int     j;
   int     knot;
   int     count;
   int     trees  = 0;
   int     layer  = 0;
   int     retval = FALSE;
   double  pobj   = 0.0;
   COMP*   tree;
   int*    next;
   int*    result;
   double* nval;

   assert(b          != NULL);
   assert(b->problem != NULL);

   g = b->problem;

   assert(g->layers  == 1);

   puts("Heuristic (Rayward-Smith) Start");

   path   = malloc((size_t)g->knots * sizeof(PATH));
   tree   = malloc((size_t)g->terms * sizeof(COMP));
   next   = malloc((size_t)g->knots * sizeof(int));
   result = malloc((size_t)g->edges * sizeof(int));
   nval   = malloc((size_t)b->vars  * sizeof(double));

   assert(path   != NULL);
   assert(tree   != NULL);
   assert(next   != NULL);
   assert(result != NULL);
   assert(nval   != NULL);

   for(i = 0; i < g->knots; i++)
   {
      if (!Is_term(g->term[i]))
         next[i] = COMP_NONE;
      else
      {
         tree[trees].start = i;
         trees++;
         next[i] = COMP_LAST;
      }
   }

   for(i = 0; i < g->knots; i++)
      g->mark[i] = TRUE;

   /* Hauptschleife...
    */
   while(trees > 1)
   {
      fputc('X', stdout);
      fflush(stdout);

      knot = find_best(g, path, trees, tree, next);

      trees = connect2(g, path, knot, trees, tree, next);

#ifdef OUTPUT
      for(i = 0; i < trees; i++)
      {
         printf("Tree: %d = { ", i);
         for(j = tree[i].start; j >= COMP_MINI; j = next[j])
            printf("%d ", j);
         printf("}\n");
      }
      fputc('\n', stdout);
#endif
   }
   for(i = 0; i < g->knots; i++)
      g->mark[i] = (next[i] != COMP_NONE);

   graph_path_exec(g, MST_MODE, g->source[layer], g->cost, path);

   for(i = 0; i < g->edges; i++)
      result[i] = -1;

   for(i = 0; i < g->knots; i++)
   {
      if (g->mark[i] && (path[i].edge != -1))
      {
         assert(g->head[path[i].edge] == i);
         assert(result[path[i].edge] == -1);

         result[path[i].edge] = layer;
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
                  count++;
                  break;
               }
            }
            assert(j != EAT_LAST);
         }
      }
   }
   while(count > 0);

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
         memcpy(b->xval, nval, b->vars * sizeof(double));
         b->pobj = pobj;
         retval  = TRUE;

         printf("Heuristic: New Solution Value=%.12e\n", pobj);
      }
   }
   free(path);
   free(tree);
   free(next);
   free(result);
   free(nval);
/* free(cost);*/

   return(retval);
}
