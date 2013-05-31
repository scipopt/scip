/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Functions                                                     */
/*   File....: grphbase.c                                                    */
/*   Name....: Basic Graph Routines                                          */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*lint -esym(750,GRPHBASE_C) -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h)                                                   */

/* Ein Graph wird initialisiert, und dannach veraendert sich die Anzahl
 * seiner Knoten 'g->knots' und Kanten 'g->edges' nicht mehr kleiner.
 * Allerdings kann der Grad eines Knotens auf 0 zurueckgehen und eine
 * Kante als EAT_FREE gekennzeichnet werden. Wird dann 'graph_pack()'
 * aufgerufen, werden diese Knoten und Kanten nicht uebernommen.
 */
#define GRPHBASE_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "portab.h"
#include "grph.h"

GRAPH* graph_init(
   int ksize,
   int esize,
   int layers,
   int flags)
{
   GRAPH* p;
   int    i;

   assert(ksize > 0);
   assert(ksize < INT_MAX);
   assert(esize > 0);
   assert(esize < INT_MAX);
   assert(layers > 0);
   assert(layers < SHRT_MAX);

   p = malloc(sizeof(*p));

   assert(p != NULL);

   p->ksize  = ksize;
   p->knots  = 0;
   p->terms  = 0;
   p->flags  = flags;
   p->layers = layers;
   p->locals = malloc((size_t)layers * sizeof(int));
   p->source = malloc((size_t)layers * sizeof(int));

   p->term   = malloc((size_t)ksize * sizeof(int));
   p->mark   = malloc((size_t)ksize * sizeof(int));
   p->grad   = malloc((size_t)ksize * sizeof(int));
   p->inpbeg = malloc((size_t)ksize * sizeof(int));
   p->outbeg = malloc((size_t)ksize * sizeof(int));

   p->esize = esize;
   p->edges = 0;
   p->cost  = malloc((size_t)esize * sizeof(double));
   p->tail  = malloc((size_t)esize * sizeof(int));
   p->head  = malloc((size_t)esize * sizeof(int));

   p->ieat  = malloc((size_t)esize * sizeof(int));
   p->oeat  = malloc((size_t)esize * sizeof(int));

   p->xpos  = malloc((size_t)ksize * sizeof(int));
   p->ypos  = malloc((size_t)ksize * sizeof(int));

   assert(p->xpos   != NULL);
   assert(p->ypos   != NULL);
   assert(p->locals != NULL);
   assert(p->source != NULL);
   assert(p->term   != NULL);
   assert(p->mark   != NULL);
   assert(p->grad   != NULL);
   assert(p->inpbeg != NULL);
   assert(p->outbeg != NULL);
   assert(p->cost   != NULL);
   assert(p->tail   != NULL);
   assert(p->head   != NULL);
   assert(p->ieat   != NULL);
   assert(p->oeat   != NULL);

   for(i = 0; i < p->layers; i++)
   {
      p->source[i] = -1;
      p->locals[i] =  0;
   }
   return(p);
}

void graph_resize(
   GRAPH* p,
   int    ksize,
   int    esize,
   int    layers)
{
   int i;

   assert(p      != NULL);
   assert((ksize  < 0) || (ksize  >= p->knots));
   assert((esize  < 0) || (esize  >= p->edges));
   assert((layers < 0) || (layers >= p->layers));

   if ((layers > 0) && (layers != p->layers))
   {
      p->locals = realloc(p->locals, (size_t)layers * sizeof(int));
      p->source = realloc(p->source, (size_t)layers * sizeof(int));

      for(i = p->layers; i < layers; i++)
      {
         p->source[i] = -1;
         p->locals[i] =  0;
      }
      p->layers = layers;
   }
   if ((ksize > 0) && (ksize != p->ksize))
   {
      p->ksize  = ksize;
      p->term   = realloc(p->term,   (size_t)ksize * sizeof(int));
      p->mark   = realloc(p->mark,   (size_t)ksize * sizeof(int));
      p->grad   = realloc(p->grad,   (size_t)ksize * sizeof(int));
      p->inpbeg = realloc(p->inpbeg, (size_t)ksize * sizeof(int));
      p->outbeg = realloc(p->outbeg, (size_t)ksize * sizeof(int));
      p->xpos   = realloc(p->xpos,   (size_t)ksize * sizeof(int));
      p->ypos   = realloc(p->ypos,   (size_t)ksize * sizeof(int));
   }
   if ((esize > 0) && (esize != p->esize))
   {
      p->esize = esize;
      p->cost  = realloc(p->cost, (size_t)esize * sizeof(double));
      p->tail  = realloc(p->tail, (size_t)esize * sizeof(int));
      p->head  = realloc(p->head, (size_t)esize * sizeof(int));

      p->ieat  = realloc(p->ieat, (size_t)esize * sizeof(int));
      p->oeat  = realloc(p->oeat, (size_t)esize * sizeof(int));
   }
   assert(p->locals != NULL);
   assert(p->source != NULL);
   assert(p->term   != NULL);
   assert(p->mark   != NULL);
   assert(p->grad   != NULL);
   assert(p->inpbeg != NULL);
   assert(p->outbeg != NULL);
   assert(p->cost   != NULL);
   assert(p->tail   != NULL);
   assert(p->head   != NULL);
   assert(p->ieat   != NULL);
   assert(p->oeat   != NULL);
   assert(p->xpos   != NULL);
   assert(p->ypos   != NULL);
}

void graph_free(
   GRAPH* p)
{
   assert(p != NULL);

   free(p->locals);
   free(p->source);
   free(p->term);
   free(p->mark);
   free(p->grad);
   free(p->inpbeg);
   free(p->outbeg);
   free(p->cost);
   free(p->tail);
   free(p->head);
   free(p->ieat);
   free(p->oeat);
   free(p->xpos);
   free(p->ypos);

   free(p);
}

GRAPH* graph_copy(
   const GRAPH* p)
{
   GRAPH* g;

   assert(p != NULL);

   g = graph_init(p->ksize, p->esize, p->layers, p->flags);

   assert(g         != NULL);
   assert(g->xpos   != NULL);
   assert(g->ypos   != NULL);
   assert(g->locals != NULL);
   assert(g->source != NULL);
   assert(g->term   != NULL);
   assert(g->mark   != NULL);
   assert(g->grad   != NULL);
   assert(g->inpbeg != NULL);
   assert(g->outbeg != NULL);
   assert(g->cost   != NULL);
   assert(g->tail   != NULL);
   assert(g->head   != NULL);
   assert(g->ieat   != NULL);
   assert(g->oeat   != NULL);

   g->knots = p->knots;
   g->terms = p->terms;
   g->edges = p->edges;

   memcpy(g->locals, p->locals, p->layers * sizeof(*p->locals));
   memcpy(g->source, p->source, p->layers * sizeof(*p->source));
   memcpy(g->term,   p->term,   p->ksize  * sizeof(*p->term));
   memcpy(g->mark,   p->mark,   p->ksize  * sizeof(*p->mark));
   memcpy(g->grad,   p->grad,   p->ksize  * sizeof(*p->grad));
   memcpy(g->inpbeg, p->inpbeg, p->ksize  * sizeof(*p->inpbeg));
   memcpy(g->outbeg, p->outbeg, p->ksize  * sizeof(*p->outbeg));

   memcpy(g->cost,   p->cost,   p->esize  * sizeof(*p->cost));
   memcpy(g->tail,   p->tail,   p->esize  * sizeof(*p->tail));
   memcpy(g->head,   p->head,   p->esize  * sizeof(*p->head));
   memcpy(g->ieat,   p->ieat,   p->esize  * sizeof(*p->ieat));
   memcpy(g->oeat,   p->oeat,   p->esize  * sizeof(*p->oeat));

   memcpy(g->xpos,   p->xpos,   p->ksize  * sizeof(*p->xpos));
   memcpy(g->ypos,   p->ypos,   p->ksize  * sizeof(*p->ypos));

   assert(graph_valid(p));

   return g;
}

void graph_flags(
   GRAPH* p,
   int    flags)
{
   assert(p     != NULL);
   assert(flags >= 0);

   p->flags |= flags;
}

void graph_show(
   const GRAPH* p)
{
   int i;

   assert(p != NULL);

   for(i = 0; i < p->knots; i++)
      if (p->grad[i] > 0)
         (void)printf("Knot %d, term=%d, grad=%d, inpbeg=%d, outbeg=%d\n",
            i, p->term[i], p->grad[i], p->inpbeg[i], p->outbeg[i]);

   (void)fputc('\n', stdout);

   for(i = 0; i < p->edges; i++)
      if (p->ieat[i] != EAT_FREE)
         (void)printf("Edge %d, cost=%g, tail=%d, head=%d, ieat=%d, oeat=%d\n",
            i, p->cost[i], p->tail[i], p->head[i], p->ieat[i], p->oeat[i]);

   (void)fputc('\n', stdout);
}

void graph_ident(
   const GRAPH* p)
{
   int i;
   int ident = 0;

   assert(p != NULL);

   for(i = 0; i < p->knots; i++)
      ident += (i + 1) * (p->term[i] * 2 + p->grad[i] * 3
         + p->inpbeg[i] * 5 + p->outbeg[i] * 7);

   for(i = 0; i < p->edges; i++)
      ident += (i + 1) * ((int)p->cost[i] + p->tail[i]
         + p->head[i] + p->ieat[i] + p->oeat[i]);

   (void)printf("Graph Ident = %d\n", ident);
}

/* ARGSUSED */
void graph_knot_add(
   GRAPH* p,
   int    term,
   int    xpos,
   int    ypos)
{
   assert(p        != NULL);
   assert(p->ksize >  p->knots);
   assert(term     <  p->layers);

   p->term  [p->knots] = term;
   p->mark  [p->knots] = TRUE;
   p->grad  [p->knots] = 0;
   p->inpbeg[p->knots] = EAT_LAST;
   p->outbeg[p->knots] = EAT_LAST;
   p->xpos  [p->knots] = xpos;
   p->ypos  [p->knots] = ypos;

   if (Is_term(term))
   {
      p->terms++;
      p->locals[term]++;
   }
   p->knots++;
}

/* ARGSUSED */
void graph_knot_chg(
   GRAPH* p,
   int    knot,
   int    term,
   int    xpos,
   int    ypos)
{
   assert(p      != NULL);
   assert(knot   >= 0);
   assert(knot   < p->knots);
   assert(term   < p->layers);

   if (term != NO_CHANGE)
   {
      if (term != p->term[knot])
      {
         if (Is_term(p->term[knot]))
         {
            p->terms--;
            p->locals[p->term[knot]]--;
         }
         p->term[knot] = term;

         if (Is_term(p->term[knot]))
         {
            p->terms++;
            p->locals[p->term[knot]]++;
         }
      }
   }
   if (xpos >= 0)
      p->xpos[knot] = xpos;

   if (ypos >= 0)
      p->ypos[knot] = ypos;
}

void graph_knot_contract(
   GRAPH* p,
   int    t,
   int    s)
{
   typedef struct save_list
   {
      unsigned int mark;
      signed int   knot;
      double       cost;
   } SLIST;

   SLIST* slp;
   int    slc = 0;
   int    i;
   int    et;
   int    es;
   int    head;
   int    tail;

   assert(p          != NULL);
   assert(t          >= 0);
   assert(t          <  p->knots);
   assert(s          >= 0);
   assert(s          <  p->knots);
   assert(s          != t);
   assert(p->grad[s] >  0);
   assert(p->grad[t] >  0);
   assert(p->layers  == 1);
   assert(graph_valid(p));

   /* Feindliche Uebernahme des Terminals.
    */
   if (Is_term(p->term[s]))
      graph_knot_chg(p, t, p->term[s], -1, -1);

   /* Die Quelle darf nicht versiegen !
    */
   if (p->source[0] == s)
      p->source[0] = t;

   slp = malloc((size_t)p->grad[s] * sizeof(SLIST));

   assert(slp != NULL);

      /* Liste mit Kanten des aufzuloesenden Kontens merken
    */
   for(es = p->outbeg[s]; es != EAT_LAST; es = p->oeat[es])
   {
      assert(p->tail[es] == s);

      if (p->head[es] != t)
      {
         /*
         assert(EQ(p->cost[es], p->cost[Edge_anti(es)]));
         */
         slp[slc].mark = FALSE;
         slp[slc].knot = p->head[es];
         slp[slc].cost = p->cost[es];
         slc++;

         assert(slc < p->grad[s]);
      }
   }
   assert(slc == p->grad[s] - 1);

   /* Kantenliste durchgehen
    */
   for(i = 0; i < slc; i++)
   {
      /* Hat t schon eine Kante mit diesem Ziel ?
       */
      for(et = p->outbeg[t]; et != EAT_LAST; et = p->oeat[et])
         if (p->head[et] == slp[i].knot)
            break;

      /* Keine gefunden, Kante aus der Liste muss eingefuegt werden.
       */
      if (et == EAT_LAST)
          slp[i].mark = TRUE;
      else
      {
         /* Ja ist vorhanden !
          */
         assert(et != EAT_LAST);

         /* Muessen die Kosten korrigiert werden ?
          */
         if (GT(p->cost[et], slp[i].cost))
         {
            assert(EQ(p->cost[et], p->cost[Edge_anti(et)]));

            p->cost[et]            = slp[i].cost;
            p->cost[Edge_anti(et)] = slp[i].cost;
         }
      }
   }
   /* Einzufuegenden Kanten einfuegen
    */
   for(i = 0; i < slc; i++)
   {
      if (slp[i].mark)
      {
         es = p->outbeg[s];

         assert(es != EAT_LAST);

         graph_edge_del(p, es);

         head = slp[i].knot;
         tail = t;

         p->grad[head]++;
         p->grad[tail]++;

         p->cost[es]     = slp[i].cost;
         p->tail[es]     = tail;
         p->head[es]     = head;
         p->ieat[es]     = p->inpbeg[head];
         p->oeat[es]     = p->outbeg[tail];
         p->inpbeg[head] = es;
         p->outbeg[tail] = es;

         es = Edge_anti(es);

         p->cost[es]     = slp[i].cost;
         p->tail[es]     = head;
         p->head[es]     = tail;
         p->ieat[es]     = p->inpbeg[tail];
         p->oeat[es]     = p->outbeg[head];
         p->inpbeg[tail] = es;
         p->outbeg[head] = es;
      }
   }
   /* Alle uebrigen Kanten Loeschen.
    */
   while(p->outbeg[s] != EAT_LAST)
      graph_edge_del(p, p->outbeg[s]);

   free(slp);

   assert(p->grad[s]   == 0);
   assert(p->outbeg[s] == EAT_LAST);
   assert(p->inpbeg[s] == EAT_LAST);
   assert(graph_valid(p));
}

/* ARGSUSED */
void graph_edge_add(
   GRAPH* p,
   int    tail,
   int    head,
   double cost1,
   double cost2)
{
   int    e;

   assert(p      != NULL);
   assert(GE(cost1, 0));
   assert(GE(cost2, 0));
   assert(tail   >= 0);
   assert(tail   <  p->knots);
   assert(head   >= 0);
   assert(head   <  p->knots);

   assert(p->esize >= p->edges + 2);

   e = p->edges;

   p->grad[head]++;
   p->grad[tail]++;

   p->cost[e]      = cost1;
   p->tail[e]      = tail;
   p->head[e]      = head;
   p->ieat[e]      = p->inpbeg[head];
   p->oeat[e]      = p->outbeg[tail];
   p->inpbeg[head] = e;
   p->outbeg[tail] = e;

   e++;

   p->cost[e]      = cost2;
   p->tail[e]      = head;
   p->head[e]      = tail;
   p->ieat[e]      = p->inpbeg[tail];
   p->oeat[e]      = p->outbeg[head];
   p->inpbeg[tail] = e;
   p->outbeg[head] = e;

   p->edges += 2;
}

inline static void edge_remove(
   GRAPH* p,
   int    e)
{
   int    i;
   int    head;
   int    tail;

   assert(p          != NULL);
   assert(e          >= 0);
   assert(e          <  p->edges);

   head = p->head[e];
   tail = p->tail[e];

   if (p->inpbeg[head] == e)
      p->inpbeg[head] = p->ieat[e];
   else
   {
      for(i = p->inpbeg[head]; p->ieat[i] != e; i = p->ieat[i])
         assert(i >= 0);

      p->ieat[i] = p->ieat[e];
   }
   if (p->outbeg[tail] == e)
      p->outbeg[tail] = p->oeat[e];
   else
   {
      for(i = p->outbeg[tail]; p->oeat[i] != e; i = p->oeat[i])
         assert(i >= 0);

      p->oeat[i] = p->oeat[e];
   }
}

void graph_edge_del(
   GRAPH* p,
   int    e)
{
   assert(p          != NULL);
   assert(e          >= 0);
   assert(e          <  p->edges);

   /* Immer mit der ersten von beiden Kanten anfangen
    */
   e -= e % 2;

   assert(p->head[e] == p->tail[e + 1]);
   assert(p->tail[e] == p->head[e + 1]);

   p->grad[p->head[e]]--;
   p->grad[p->tail[e]]--;

   edge_remove(p, e);

   assert(p->ieat[e] != EAT_FREE);
   assert(p->ieat[e] != EAT_HIDE);
   assert(p->oeat[e] != EAT_FREE);
   assert(p->oeat[e] != EAT_HIDE);

   p->ieat[e] = EAT_FREE;
   p->oeat[e] = EAT_FREE;

   e++;

   edge_remove(p, e);

   assert(p->ieat[e] != EAT_FREE);
   assert(p->ieat[e] != EAT_HIDE);
   assert(p->oeat[e] != EAT_FREE);
   assert(p->oeat[e] != EAT_HIDE);

   p->ieat[e] = EAT_FREE;
   p->oeat[e] = EAT_FREE;
}

void graph_edge_hide(
   GRAPH* p,
   int    e)
{
   assert(p          != NULL);
   assert(e          >= 0);
   assert(e          <  p->edges);

   /* Immer mit der ersten von beiden Anfangen
    */
   e -= e % 2;

   assert(p->head[e] == p->tail[e + 1]);
   assert(p->tail[e] == p->head[e + 1]);

   p->grad[p->head[e]]--;
   p->grad[p->tail[e]]--;

   edge_remove(p, e);

   assert(p->ieat[e] != EAT_FREE);
   assert(p->ieat[e] != EAT_HIDE);
   assert(p->oeat[e] != EAT_FREE);
   assert(p->oeat[e] != EAT_HIDE);

   p->ieat[e] = EAT_HIDE;
   p->oeat[e] = EAT_HIDE;

   e++;

   edge_remove(p, e);

   assert(p->ieat[e] != EAT_FREE);
   assert(p->ieat[e] != EAT_HIDE);
   assert(p->oeat[e] != EAT_FREE);
   assert(p->oeat[e] != EAT_HIDE);

   p->ieat[e] = EAT_HIDE;
   p->oeat[e] = EAT_HIDE;
}

void graph_uncover(
   GRAPH* p)
{
   int head;
   int tail;
   int e;

   assert(p      != NULL);

   for(e = 0; e < p->edges; e++)
   {
      if (p->ieat[e] == EAT_HIDE)
      {
         assert(e % 2 == 0);
         assert(p->oeat[e] == EAT_HIDE);

	 head            = p->head[e];
	 tail            = p->tail[e];

	 p->grad[head]++;
	 p->grad[tail]++;

	 p->ieat[e]      = p->inpbeg[head];
	 p->oeat[e]      = p->outbeg[tail];
	 p->inpbeg[head] = e;
	 p->outbeg[tail] = e;

         e++;

         assert(p->ieat[e] == EAT_HIDE);
         assert(p->oeat[e] == EAT_HIDE);
         assert(p->head[e] == tail);
         assert(p->tail[e] == head);

	 head            = p->head[e];
	 tail            = p->tail[e];
	 p->ieat[e]      = p->inpbeg[head];
	 p->oeat[e]      = p->outbeg[tail];
	 p->inpbeg[head] = e;
	 p->outbeg[tail] = e;
      }
   }
}

GRAPH *graph_pack(
   GRAPH* p)
{
   const char* msg1   = " ok\nKnots: %d  Edges: %d  Terminals: %d\n";
   const char* msg2   = " vanished !\n";

   GRAPH* q;
   int*   new;
   int    knots = 0;
   int    edges = 0;
   int    i;
   int    l;

   assert(p      != NULL);
   assert(graph_valid(p));

   (void)printf("Packing Graph: ");

   new = malloc((size_t)p->knots * sizeof(new[0]));

   assert(new != NULL);

   /* Knoten zaehlen
    */
   for(i = 0; i < p->knots; i++)
   {
      new[i] = knots;

      /* Hat der Knoten noch Kanten ?
       */
      if (p->grad[i] > 0)
         knots++;
      else
         new[i] = -1;
   }

   /* Ist ueberhaupt noch ein Graph vorhanden ?
    */
   if (knots == 0)
   {
      free(new);
      graph_free(p);

      (void)printf(msg2);

      return(NULL);
   }

   /* Kanten zaehlen
    */
   for(i = 0; i < p->edges; i++)
   {
      if (p->oeat[i] != EAT_FREE)
      {
         assert(p->ieat[i] != EAT_FREE);
         edges++;
      }
   }

   q = graph_init(knots, edges, p->layers, p->flags);

   /* Knoten umladen
    */
   for(i = 0; i < p->knots; i++)
   {
      assert(p->term[i] < p->layers);

      if ((i % 100) == 0)
      {
         (void)fputc('k', stdout);
         (void)fflush(stdout);
      }

      if (p->grad[i] > 0)
         graph_knot_add(q, p->term[i], p->xpos[i], p->ypos[i]);
   }

   /* Kanten umladen
    */
   for(i = 0; i < p->edges; i += 2)
   {
      if ((i % 1000) == 0)
      {
         (void)fputc('e', stdout);
         (void)fflush(stdout);
      }

      if (p->ieat[i] == EAT_FREE)
      {
         assert(p->oeat[i]     == EAT_FREE);
         assert(p->ieat[i + 1] == EAT_FREE);
         assert(p->oeat[i + 1] == EAT_FREE);
         continue;
      }
      assert(p->ieat[i]      != EAT_FREE);
      assert(p->oeat[i]      != EAT_FREE);
      assert(p->ieat[i + 1]  != EAT_FREE);
      assert(p->oeat[i + 1]  != EAT_FREE);
      assert(new[p->tail[i]] >= 0);
      assert(new[p->head[i]] >= 0);

      graph_edge_add(q, new[p->tail[i]], new[p->head[i]],
         p->cost[i], p->cost[Edge_anti(i)]);
   }

   /* Wurzeln umladen
    */
   for(l = 0; l < q->layers; l++)
   {
      assert(q->term[new[p->source[l]]] == l);
      q->source[l] = new[p->source[l]];
   }

   free(new);

   graph_free(p);

#if 0
   for(l = 0; l < q->layers; l++)
      q->source[l] = -1;

   for(i = 0; i < q->knots; i++)
      if ((q->term[i] >= 0) && ((q->source[q->term[i]] < 0)
         || (q->grad[i] > q->grad[q->source[q->term[i]]])))
         q->source[q->term[i]] = i;
#endif
   assert(q->source[0] >= 0);

   (void)printf(msg1, q->knots, q->edges, q->terms);

   return(q);
}

void graph_trail(
   const GRAPH* p,
   int          i)
{
   int   k;

   assert(p      != NULL);
   assert(i      >= 0);
   assert(i      <  p->knots);

   if (!p->mark[i])
   {
      p->mark[i] = TRUE;

      for(k = p->outbeg[i]; k != EAT_LAST; k = p->oeat[k])
         if (!p->mark[p->head[k]])
            graph_trail(p, p->head[k]);
   }
}

int graph_valid(
   const GRAPH* p)
{
   const char* fehler1  = "*** Graph Validation Error: Head invalid, Knot %d, Edge %d, Tail=%d, Head=%d\n";
   const char* fehler2  = "*** Graph Validation Error: Tail invalid, Knot %d, Edge %d, Tail=%d, Head=%d\n";
   const char* fehler3  = "*** Graph Validation Error: Source invalid, Layer %d, Source %d, Terminal %d\n";
   const char* fehler4  = "*** Graph Validation Error: FREE invalid, Edge %d/%d\n";
   const char* fehler5  = "*** Graph Validation Error: Anti invalid, Edge %d/%d, Tail=%d/%d, Head=%d/%d\n";
   const char* fehler6  = "*** Graph Validation Error: Knot %d with Grad 0 has Edges\n";
   const char* fehler7  = "*** Graph Validation Error: Knot %d not connected\n";
   const char* fehler8  = "*** Graph Validation Error: Wrong locals count, Layer %d, count is %d, should be %d\n";
   const char* fehler9  = "*** Graph Validation Error: Wrong Terminal count, count is %d, should be %d\n";

   int    k;
   int    l;
   int    e;
   int    terms;
   int*   locals;

   assert(p      != NULL);

   terms  = p->terms;
   locals = malloc((size_t)p->layers * sizeof(int));

   assert(locals != NULL);

   for(l = 0; l < p->layers; l++)
      locals[l] = p->locals[l];

   for(k = 0; k < p->knots; k++)
   {
      if (Is_term(p->term[k]))
      {
         locals[p->term[k]]--;
         terms--;
      }
      for(e = p->inpbeg[k]; e != EAT_LAST; e = p->ieat[e])
         if (p->head[e] != k)
            break;

      if (e != EAT_LAST)
         return((void)fprintf(stderr, fehler1, k, e, p->tail[e], p->head[e]), FALSE);

      for(e = p->outbeg[k]; e != EAT_LAST; e = p->oeat[e])
         if (p->tail[e] != k)
            break;

      if (e != EAT_LAST)
         return((void)fprintf(stderr, fehler2, k, e, p->tail[e], p->head[e]), FALSE);
   }
   if (terms != 0)
      return((void)fprintf(stderr, fehler9, p->terms, p->terms - terms), FALSE);

   for(l = 0; l < p->layers; l++)
   {
      if (locals[l] != 0)
         return((void)fprintf(stderr, fehler8,
            l, p->locals[l], p->locals[l] - locals[l]), FALSE);

      if ((p->source[l] < 0)
       || (p->source[l] >= p->knots)
       || (p->term[p->source[l]] != l))
         return((void)fprintf(stderr, fehler3,
            l, p->source[l], p->term[p->source[l]]), FALSE);
   }
   free(locals);

   for(e = 0; e < p->edges; e += 2)
   {
      if ((p->ieat[e    ] == EAT_FREE) && (p->oeat[e    ] == EAT_FREE)
       && (p->ieat[e + 1] == EAT_FREE) && (p->oeat[e + 1] == EAT_FREE))
         continue;

      if ((p->ieat[e] == EAT_FREE) || (p->oeat[e] == EAT_FREE)
       || (p->ieat[e + 1] == EAT_FREE) || (p->oeat[e + 1] == EAT_FREE))
         return((void)fprintf(stderr, fehler4, e, e + 1), FALSE);

      if ((p->head[e] != p->tail[e + 1]) || (p->tail[e] != p->head[e + 1]))
         return((void)fprintf(stderr, fehler5,
            e, e + 1, p->head[e], p->tail[e + 1],
            p->tail[e], p->head[e + 1]), FALSE);

   }
   for(k = 0; k < p->knots; k++)
      p->mark[k] = FALSE;

   graph_trail(p, p->source[0]);

   for(k = 0; k < p->knots; k++)
   {
      if ((p->grad[k] == 0)
         && ((p->inpbeg[k] != EAT_LAST) || (p->outbeg[k] != EAT_LAST)))
         return((void)fprintf(stderr, fehler6, k), FALSE);

      if (!p->mark[k] && (p->grad[k] > 0))
         return((void)fprintf(stderr, fehler7, k), FALSE);
   }
   return(TRUE);
}
