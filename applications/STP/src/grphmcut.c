/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probdata_stp.c
 * @brief  Minimum cut routine for Steiner problems
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 * This file implements a graph minimum cut routine for Steiner problems. For more details see \ref MINCUT page.
 *
 * @page MINCUT Graph minimum cut routine
 *
 * The implemented algorithm is described in "A Faster Algorithm for Finding the Minimum Cut in a Graph" by Hao and Orlin.
 *
 * A list of all interface methods can be found in grph.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "portab.h"
#include "grph.h"

#define DEBUG        0        /* 0 = No, 1 = Validation, 2 = Show flow       */
#define STATIST      0

#ifdef NDEBUG
#undef STATIST
#undef DEBUG
#define STATIST      0
#define DEBUG        0
#endif

#define Q_LAST     -1         /* Last Element of Queue     */
#define Q_NMOQ     -2         /* Not a Member Of the Queue */

#if DEBUG
static int  is_valid(const GRAPH*, const int, const int, const int*, const int *);
static void show_flow(const GRAPH*, const int*, const int*);
#endif
#if STATIST
static void cut_statist(void);
static void cut_sum(const GRAPH*, const int*, const int*);
#endif

#if STATIST
static int    s_pushes = 0;
static int    n_pushes = 0;
static int    m_pushes = 0;
static int    x_pushes = 0;
static int    relabels = 0;
static int    s_sleeps = 0;
static int    m_sleeps = 0;
static int    searches = 0;
static int    cutsums  = 0;
#endif

/*---------------------------------------------------------------------------*/
/*--- Name     : GRAPH MINimumCUT INITialise                              ---*/
/*--- Function : Holt den Speicher fuer die Hilfsarrays die wir brauchen. ---*/
/*--- Parameter: Graph                                                    ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
SCIP_RETCODE graph_mincut_init(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                p                   /**< graph data structure */
     )
{
   assert(p    != NULL);
   assert(p->mincut_dist == NULL);
   assert(p->mincut_head == NULL);
   assert(p->mincut_numb == NULL);
   assert(p->mincut_prev == NULL);
   assert(p->mincut_next == NULL);
   assert(p->mincut_temp == NULL);
   assert(p->mincut_e    == NULL);
   assert(p->mincut_x    == NULL);
   assert(p->mincut_r    == NULL);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_dist), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_head), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_numb), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_prev), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_next), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_temp), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_e), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_x), p->edges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_r), p->edges) );

   return SCIP_OKAY;
}

/*---------------------------------------------------------------------------*/
/*--- Name     : GRAPH MINimumCUT EXIT                                    ---*/
/*--- Function : Gibt den Speicher fuer die Hilfsarrays wieder frei.      ---*/
/*--- Parameter: Keine                                                    ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
void graph_mincut_exit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                p                   /**< graph data structure */
     )
{
   assert(p->mincut_dist != NULL);
   assert(p->mincut_head != NULL);
   assert(p->mincut_numb != NULL);
   assert(p->mincut_prev != NULL);
   assert(p->mincut_next != NULL);
   assert(p->mincut_temp != NULL);
   assert(p->mincut_e    != NULL);
   assert(p->mincut_x    != NULL);
   assert(p->mincut_r    != NULL);

   SCIPfreeMemoryArray(scip, &(p->mincut_dist));
   SCIPfreeMemoryArray(scip, &(p->mincut_head));
   SCIPfreeMemoryArray(scip, &(p->mincut_numb));
   SCIPfreeMemoryArray(scip, &(p->mincut_prev));
   SCIPfreeMemoryArray(scip, &(p->mincut_next));
   SCIPfreeMemoryArray(scip, &(p->mincut_temp));
   SCIPfreeMemoryArray(scip, &(p->mincut_e));
   SCIPfreeMemoryArray(scip, &(p->mincut_x));
   SCIPfreeMemoryArray(scip, &(p->mincut_r));

#if STATIST
   cut_statist();
#endif
}

inline static void delete(
   const int knot,
   int*      q_dist,
   int*      q_head,
   int*      q_prev,
   int*      q_next)
{
   assert(knot         >= 0);
   assert(q_dist       != NULL);
   assert(q_head       != NULL);
   assert(q_prev       != NULL);
   assert(q_next       != NULL);
   assert(q_dist[knot] >  0);

   if (q_next[knot] != Q_NMOQ)
   {
      /* Etwa Erster ?
       */
      if (q_prev[knot] == Q_LAST)
      {
         assert(q_dist[knot]         >= 0);
         assert(q_head[q_dist[knot]] == knot);

         q_head[q_dist[knot]] = q_next[knot];
      }
      else
      {
         assert(q_prev[knot]         >= 0);
         assert(q_next[q_prev[knot]] == knot);

         q_next[q_prev[knot]] = q_next[knot];
      }

      /* Sind wir auch nicht letzter ?
       */
      if (q_next[knot] != Q_LAST)
      {
         assert(q_next[knot]         >= 0);
         assert(q_prev[q_next[knot]] == knot);

         q_prev[q_next[knot]] = q_prev[knot];
      }
      q_next[knot] = Q_NMOQ;
      q_prev[knot] = Q_NMOQ;
   }
   assert(q_next[knot] == Q_NMOQ);
   assert(q_prev[knot] == Q_NMOQ);
}

inline static int insert(
   const int knot,
   int*      q_dist,
   int*      q_head,
   int*      q_prev,
   int*      q_next,
   int       dmin)
{
   assert(knot         >= 0);
   assert(q_dist       != NULL);
   assert(q_head       != NULL);
   assert(q_prev       != NULL);
   assert(q_next       != NULL);
   assert(q_dist[knot] >  0);
   assert(dmin         >= 1);

   if (q_prev[knot] == Q_NMOQ)
   {
      q_prev[knot]         = Q_LAST;
      q_next[knot]         = q_head[q_dist[knot]];

      if (q_next[knot] != Q_LAST)
         q_prev[q_next[knot]] = knot;

      q_head[q_dist[knot]] = knot;

      if (q_dist[knot] < dmin)
         dmin = q_dist[knot];
   }
   assert(q_next[knot] != Q_NMOQ);
   assert(q_prev[knot] != Q_NMOQ);
   assert(q_dist[knot] >= dmin);

   return(dmin);
}

static int bfs(
   const GRAPH* p,
   const int    s,
   const int    t,
   int*         q_dist,
   int*         q_numb,
   int*         q_temp,
   int*         w)
{
   int          i;
   int          j;
   int          k;
   int          l;
   int          visited = 0;

   assert(q_temp != NULL);
   assert(q_dist != NULL);
   assert(q_numb != NULL);
   assert(w      != NULL);
   assert(s      >= 0);
   assert(s      < p->knots);
   assert(t      >= 0);
   assert(t      < p->knots);

   /* Beginnen wir bei der Senke
    * (Sie ist sich selbst unendlich nah).
    */
   q_temp[visited++] = t;
   q_dist[t]         = 0;

   /* Solange noch schon besuchte Knoten da sind, von denen aus nicht
    * versucht wurde weiter zu kommen:
    */
   for(j = 0; (j < visited) && (visited < p->knots); j++)
   {
      assert(visited < p->knots);
      assert(j       < visited);

      i = q_temp[j];

      assert(i         >= 0);
      assert(i         <  p->knots);
      assert(q_dist[i] >= 0);
      assert(q_dist[i] <  visited);
      assert(w[i]      == 0);

      assert((j == 0)           || (q_dist[i] >= q_dist[q_temp[j - 1]]));
      assert((j == visited - 1) || (q_dist[i] <= q_dist[q_temp[j + 1]]));

      /* Wo koennen wir den ueberall hin:
       */
      for(k = p->outbeg[i]; k != EAT_LAST; k = p->oeat[k])
      {
         /* Wo gehts hin ? Nach Knote l.
          */
         l = p->head[k];

         /* Waren wir da noch nicht ?
          */
         assert(!w[l] || (q_dist[l] >= 0));

         if (q_dist[l] < 0)
         {
            q_dist[l]         = q_dist[i] + 1;
            q_temp[visited++] = l;
            q_numb[q_dist[l]]++;

            assert(q_dist[l] < p->knots);
         }
      }
   }
   return(visited);
}


/*---------------------------------------------------------------------------*/
/*--- Name     : INITialise LABELS                                        ---*/
/*--- Function : Fuehrt eine BFS durch um die Distanzen der Knoten von    ---*/
/*---            der Senke zu ermitten. Kanten ohne Kapazitaet werden     ---*/
/*---            nicht beruecksichtigt, ebenso Knoten die nur ueber die   ---*/
/*---            Quelle zu erreichen sind.                                ---*/
/*--- Parameter: Graph, Quelle, Senke, Kantenkapazitaeten.                ---*/
/*--- Returns  : Anzahl der Aktiven (Erreichbaren) Knoten,                ---*/
/*---            Fuellt a[] mit den Nummern dieser Knoten und d[] mit den ---*/
/*---            Entfernungen der zugehoerigen Knoten zur Senke.          ---*/
/*---            a[] ist aufsteigend sortiert.                            ---*/
/*---------------------------------------------------------------------------*/
static void initialise(
   const GRAPH* p,
   const int    s,
   const int    t,
   const int*   capa,
   int*         q_dist,
   int*         q_numb,
   int*         q_head,
   int*         q_prev,
   int*         q_next,
   int*         q_temp,
   int*         excess,
   int*         transx,
   int*         residi,
   int*         w,
   int*         dmax)
{

   int i;
   int j;
   int k;

   assert(p      != NULL);
   assert(s      >= 0);
   assert(s      <  p->knots);
   assert(t      >= 0);
   assert(t      <  p->knots);
   assert(s      != t);
   assert(capa   != NULL);
   assert(w      != NULL);
   assert(q_head != NULL);
   assert(q_dist != NULL);
   assert(q_numb != NULL);
   assert(q_prev != NULL);
   assert(q_next != NULL);
   assert(q_temp != NULL);
   assert(excess != NULL);
   assert(transx != NULL);
   assert(residi != NULL);
   assert(p->mincut_r != NULL);
   assert(p->mincut_x != NULL);

   /* Knotenarrays Initialisieren
    */
   *dmax = 1;

   for(i = 0; i < p->knots; i++)
   {
      excess[i] =  0;
      w     [i] =  0;
      q_prev[i] =  Q_NMOQ;
      q_next[i] =  Q_NMOQ;
      q_head[i] =  Q_LAST;
      q_numb[i] =  0;
      q_dist[i] = -1;
   }
   /* Jetzt die Kantenarrays.
    */
   for(k = 0; k < p->edges; k++)
   {
      transx[k] = 0;
      residi[k] = capa[k];
   }
   /* Jetzt noch dist und numb.
    */
   (void)bfs(p, s, t, q_dist, q_numb, q_temp, w);

   /* Alles was wir nicht erreichen konnten schlafen legen.
    */
   for(i = 0; i < p->knots; i++)
      if (q_dist[i] < 0)
         w[i] = *dmax + 1;

   /* Quelle einschlaefern
    */
   w[s] = 1;    /* dmax */

   /* Falls wir die Quelle s nicht erreichen konnten sind wir fertig.
    */
   if (q_dist[s] < 0)
      return;

   assert(w[s]      >  0);
   assert(q_dist[s] >  0);
   assert(w[t]      == 0);
   assert(q_dist[t] == 0);

   /* Label der Quelle abziehen
    */
   q_numb[q_dist[s]]--;

   /* Von der Quelle alles wegschieben wofuer Kapazitaeten da sind.
    */
   for(k = p->outbeg[s]; k != EAT_LAST; k = p->oeat[k])
   {
      if (capa[k] == 0)
         continue;

      j         = p->head[k];
      transx[k] = capa[k];
      residi[k] = 0;                                         /* -= transx[k] */

      residi[Edge_anti(k)] += transx[k];     /* Ueberfluessig weil w[s] == 1 */
      excess[j]            += transx[k];

      if (j != t)
         (void)insert(j, q_dist, q_head, q_prev, q_next, 1);

      assert(w[j]                   == 0);
      assert(excess[j]              >  0);
      assert((j == t) || (q_next[j] != Q_NMOQ));
      assert((j == t) || (q_prev[j] != Q_NMOQ));

      assert((p->mincut_r)[k] + (p->mincut_r)[Edge_anti(k)] == capa[k] + capa[Edge_anti(k)]);
      assert((p->mincut_x)[k]                   >= 0);
      assert((p->mincut_x)[k]                   <= capa[k]);
      assert((p->mincut_r)[k]                   == capa[k] - (p->mincut_x)[k] + (p->mincut_x)[Edge_anti(k)]);
   }
#if DEBUG > 1
   show_flow(p, capa, w);
#endif
#if DEBUG > 0
   assert(is_valid(p, s, t, capa, w));
#endif
}

static void reinitialise(
   const GRAPH* p,
   const int    s,
   const int    t,
   const int*   capa,
   int*         q_dist,
   int*         q_numb,
   int*         q_head,
   int*         q_prev,
   int*         q_next,
   int*         q_temp,
   int*         excess,
   int*         transx,
   int*         residi,
   int*         w,
   int*         dmax)
{
   int wt;
   int i;
   int j;
   int k;
   int visited;

   assert(p      != NULL);
   assert(s      >= 0);
   assert(s      <  p->knots);
   assert(t      >= 0);
   assert(t      <  p->knots);
   assert(s      != t);
   assert(capa   != NULL);
   assert(w      != NULL);
   assert(q_head != NULL);
   assert(q_dist != NULL);
   assert(q_numb != NULL);
   assert(q_prev != NULL);
   assert(q_next != NULL);
   assert(q_temp != NULL);
   assert(excess != NULL);
   assert(transx != NULL);
   assert(residi != NULL);

   /* Knotenarrays Initialisieren
    */
   assert(w[s]);

   wt    = (w[t] == 0) ? p->knots + 1 : w[t];
   *dmax = 1;

   for(i = 0; i < p->knots; i++)
   {
      q_numb[i] = 0;
      q_head[i] = Q_LAST;
      q_next[i] = Q_NMOQ;
      q_prev[i] = Q_NMOQ;

      if ((w[i] == 0) || (w[i] >= wt))
      {
         w     [i] = 0;
         q_dist[i] = -1;
      }
      if (w[i] > *dmax)
         *dmax = w[i];
   }
   /* Jetzt noch dist und numb.
    */
   visited = bfs(p, s, t, q_dist, q_numb, q_temp, w);

   /* Alles was wir nicht erreichen konnten schlafen legen.
    */
   for(i = 0; i < p->knots; i++)
      if (q_dist[i] < 0)
         w[i] = *dmax + 1;

   /* Jetzt die Kantenarrays und ggf. e updaten.
    */
   for(k = 0; k < p->edges; k += 2)
   {
      i = p->head[k];
      j = p->tail[k];

      if (!w[i] && !w[j])
      {
         assert(w[s]);

         excess[i]    += transx[k + 1] - transx[k];
         excess[j]    += transx[k] - transx[k + 1];
         transx[k]     = 0;
         residi[k]     = capa[k];
         transx[k + 1] = 0;
         residi[k + 1] = capa[k + 1];
      }
   }
   assert(w[t]      == 0);
   assert(q_dist[t] == 0);
   assert(q_temp[0] == t);

   /* Jetzt noch die mit Excess einsortieren.
    */
   for(i = 1; i < visited; i++)
   {
      assert(w[q_temp[i]] == 0);
      assert(q_temp[i]    != s);
      assert(q_temp[i]    != t);

      if (excess[q_temp[i]] > 0)
         (void)insert(q_temp[i], q_dist, q_head, q_prev, q_next, 1);
   }
#if DEBUG > 1
   show_flow(p, capa, w);
#endif
#if DEBUG > 0
   assert(is_valid(p, s, t, capa, w));
#endif
}

/*---------------------------------------------------------------------------*/
/*--- Name     : GRAPH MINimumCUT EXECute                                 ---*/
/*--- Function : Fuehrt den Mincut Algorithmus durch und findet           ---*/
/*---            (hoffentlich) einen Minimalen (s,t) Schnitt im Graphen.  ---*/
/*--- Parameter: Graph, Quelle, Senke, Kantenkapazitaeten, Zustandsarray, ---*/
/*---            Flag um vorhandenen Fluss zu belassen.                   ---*/
/*--- Returns  : Nichts, fuellt aber w[] mit nicht Nulleintraegen fuer    ---*/
/*---            die Knoten, die auf der Quellenseite des Schnittes       ---*/
/*---            liegen und Null fuer die auf der Senkenseite.            ---*/
/*---------------------------------------------------------------------------*/
void graph_mincut_exec(
   GRAPH*       p,
   int          s,
   int          t,
   const int*   capa,
   int*         w,
   int          rerun)
{
   int    min_dist;
   int    min_capa;
   int    min_knot;
   int    min_arc;
   int    dmax;
   int    i;
   int    k;
   int    dmin = 1;
   int*   dist;
   int*   head;
   int*   numb;
   int*   prev;
   int*   next;
   int*   temp;
   int*   e;
   int*   x;
   int*   r;

   assert(p      != NULL);
   assert(s      >= 0);
   assert(s      <  p->knots);
   assert(t      >= 0);
   assert(t      <  p->knots);
   assert(s      != t);
   assert(capa   != NULL);
   assert(w      != NULL);
   assert(p->mincut_dist   != NULL);
   assert(p->mincut_numb   != NULL);
   assert(p->mincut_head   != NULL);
   assert(p->mincut_prev   != NULL);
   assert(p->mincut_next   != NULL);
   assert(p->mincut_temp   != NULL);
   assert(p->mincut_e      != NULL);
   assert(p->mincut_x      != NULL);
   assert(p->mincut_r      != NULL);

   dist = p->mincut_dist;
   head = p->mincut_head;
   numb = p->mincut_numb;
   prev = p->mincut_prev;
   next = p->mincut_next;
   temp = p->mincut_temp;
   e    = p->mincut_e;
   x    = p->mincut_x;
   r    = p->mincut_r;

#if DEBUG > 0
   (void)printf("graph_mincut_exec(p, s=%d, t=%d, capa, w, rerun=%d)\n",
      s, t, rerun);
#endif

   if (!rerun)
      initialise(p, s, t, capa, dist, numb, head, prev, next, temp, e, x, r, w, &dmax);
   else
      reinitialise(p, s, t, capa, dist, numb, head, prev, next, temp, e, x, r, w, &dmax);

   /* Solange wir nicht fertig sind ...
    */
   for(;;)
   {
#if DEBUG > 0
      assert(is_valid(p, s, t, capa, w));
#endif
#if STATIST
      searches++;
#endif

      /* Kein Knoten, keine Arbeit !
       */
      while((dmin < p->knots) && (head[dmin] == Q_LAST))
         dmin++;

      if (dmin == p->knots)
         break;

      /* Hole den naechsten Knoten mit Ueberschuss hervor.
       */
      i = head[dmin];

      assert(prev[i] == Q_LAST);
      assert(w[i]    == 0);
      assert(i       != t);
      assert(i       != s);

      /* Wenn es von i aus eine gangbare Kante gibt, versuche den
       * Ueberschuss wegzuschieben.
       * (Gangbar = hat Kapazitaet und ist in Bezug auf d[] eins abschuessig).
       */
      min_knot = -1;
      min_dist =  p->knots;
      min_capa =  0;
      min_arc  = -1;

      for(k = p->outbeg[i]; k != EAT_LAST; k = p->oeat[k])
      {
         /* Hat k keine Kapazitaet oder geht zu einem eingeschlafenen Knoten ?
          * --- Dann vergiss es !
          */
         if ((r[k] <= 0) || (w[p->head[k]]))
            continue;

         /* Fuhert k nach oben ?
          * --- Dann keine Chance was runterzuschieben, aber wir Merken uns
          *     den Koten, denn vielleicht steigen wird ja noch auf.
          */
         if (dist[i] != dist[p->head[k]] + 1)
         {
            assert(dist[i] <= dist[p->head[k]]);  /* Warum funktioniert das ? tut es das ? */

            if ((dist[p->head[k]] < min_dist)
             || ((dist[p->head[k]] == min_dist) && (r[k] > min_capa)))
            {
               min_knot = p->head[k];
               min_dist = dist[min_knot];
               min_capa = r[k];
               min_arc  = k;
            }
            continue;
         }

         /* Druecke delta := min{ e(i), r(i,j) } Einheiten Fluss von
          * Knoten i nach Knoten j.
          */
         assert(Min(e[i], r[k]) > 0);

         /* Mehr Kapazitaet als Ueberschuss ?
          */
         if (e[i] <= r[k])
         {
#if STATIST
            (e[i] == r[k]) ? x_pushes++ : n_pushes++;
#endif
            x[k]            += e[i];
            r[k]            -= e[i];
            r[Edge_anti(k)] += e[i];
            e[p->head[k]]   += e[i];
            e[i]             = 0;   /* -= e[i] */

            assert(e[p->head[k]] >  0);
            assert(w[p->head[k]] == 0);

            if (p->head[k] != t)
               dmin = insert(p->head[k], dist, head, prev, next, dmin);

            assert(r[k] + r[Edge_anti(k)] == capa[k] + capa[Edge_anti(k)]);
            assert(r[k]                   == capa[k] - x[k] + x[Edge_anti(k)]);

            /* Kein Ueberschuss uebrig, mit dem Knoten sind wir erstmal fertig.
             */
            break;
         }

         /* Mehr Ueberschuss als Kapazitaet.
          */
         r[Edge_anti(k)] += r[k];
         e[p->head[k]]   += r[k];
         e[i]            -= r[k];
         x[k]            += r[k];
         r[k]             = 0;      /* -= r[k] */

         if (p->head[k] != t)
            dmin = insert(p->head[k], dist, head, prev, next, dmin);

         assert(r[k] + r[Edge_anti(k)] == capa[k] + capa[Edge_anti(k)]);
         assert(r[k]                   == capa[k] - x[k] + x[Edge_anti(k)]);
         assert(e[i] > 0);

#if STATIST
         s_pushes++;
#endif
      }

      /* Kein Ueberschuss mehr, dann Knoten aus der aktiv Liste entfernen.
       */
      if (e[i] == 0)
      {
         delete(i, dist, head, prev, next);

         assert(prev[i] == Q_NMOQ);
         assert(next[i] == Q_NMOQ);

         continue;
      }
      /* Wir sind den Ueberschuss nicht losgeworden, also muss Knoten i
       * hochgelegt werden (relabel).
       */
      assert(numb[dist[i]] > 0);

      if (numb[dist[i]] == 1)
      {
         dmax++;

         assert(dmax <= p->knots);

         for(k = 0; k < p->knots; k++)
         {
            /* Kleinergleich um Knoten i auch zu erwischen.
             */
            if (!w[k] && (dist[i] <= dist[k]))
            {
               numb[dist[k]]--;
               w[k] = dmax;

               delete(k, dist, head, prev, next);

               assert(prev[k] == Q_NMOQ);
               assert(next[k] == Q_NMOQ);
            }
         }
#if STATIST
         m_sleeps++;
#endif
         continue;
      }
      /* Wir setzen d[i] um Eins hoeher als den Tiefsten drumherum der
       * aktiv ist und noch Kapazitaet auf der Leitung hat.
       */

      /* Keine Kapazitaeten mehr ?
       * --- Dann i einschlaefern.
       */
      if (min_knot == -1)
      {
         numb[dist[i]]--;
         w[i] = ++dmax;

         assert(dmax <= p->knots);

         delete(i, dist, head, prev, next);

         assert(prev[i] == Q_NMOQ);
         assert(next[i] == Q_NMOQ);
#if STATIST
         s_sleeps++;
#endif
      }
      else
      {
         assert(min_dist <  p->knots);
         assert(min_capa >  0);
         assert(min_knot >= 0);
         assert(min_arc  >= 0);

         delete(i, dist, head, prev, next);

         numb[dist[i]]--;

         dist[i] = min_dist + 1;

         (void)insert(i, dist, head, prev, next, dmin);

         numb[dist[i]]++;

         assert(dist[i] < p->knots);

         assert(min_capa         >  0);
         assert(min_capa         == r[min_arc]);
         assert(p->head[min_arc] == min_knot);
         assert(p->tail[min_arc] == i);
         assert(dist[i]          == dist[min_knot] + 1);
         assert(w[min_knot]      == 0);

         /* Wenn die Kante die wir uns gemerkt haben, genug Kapazitaet hat
          * um den Knoten zu erledigen, dann bringen wir es hinter uns.
          */
         if (e[i] <= min_capa)
         {
            x[min_arc]            += e[i];
            r[min_arc]            -= e[i];
            r[Edge_anti(min_arc)] += e[i];
            e[min_knot]           += e[i];
            e[i]                   = 0;   /* -= e[i] */

            if (p->head[min_arc] != t)
               dmin = insert(p->head[min_arc], dist, head, prev, next, dmin);

            delete(i, dist, head, prev, next);

            assert(r[min_arc] + r[Edge_anti(min_arc)] == capa[min_arc] + capa[Edge_anti(min_arc)]);
            assert(r[min_arc]                         >= 0);
            assert(r[min_arc]                         == capa[min_arc] - x[min_arc] + x[Edge_anti(min_arc)]);
#if STATIST
            m_pushes++;
#endif
         }
#if STATIST
         relabels++;
#endif
      }
   }
   assert(w[s]);
   assert(!w[t]);

#if STATIST
   cut_sum(p, capa, w);
#endif
#if DEBUG > 1
   show_flow(p, capa, w);
#endif
#if DEBUG > 0
   assert(is_valid(p, s, t, capa, w));
#endif
}

#if STATIST
static void cut_statist(void)
{
   (void)printf("Mincut Statistics:\n");
   (void)printf("Node-Searches=%d, Cut Sums=%d\n",
      searches, cutsums);
   (void)printf("S-Pushes=%d, N-Pushes=%d, X-Pushes=%d, M-Pushes=%d\n",
      s_pushes, n_pushes, x_pushes, m_pushes);
   (void)printf("Relabels=%d, S-Sleeps=%d, M-Sleeps=%d\n\n",
      relabels, s_sleeps, m_sleeps);

   s_pushes = 0;
   n_pushes = 0;
   x_pushes = 0;
   m_pushes = 0;
   relabels = 0;
   s_sleeps = 0;
   m_sleeps = 0;
   searches = 0;
   cutsums  = 0;
}

static void cut_sum(
   const GRAPH* p,
   const int*   capa,
   const int*   w)
{
   int          sum = 0;
   int          i;
   int          j;
   int          k;

   assert(p      != NULL);
   assert(capa   != NULL);
   assert(w      != NULL);

   for(k = 0; k < p->edges; k++)
   {
      i = p->head[k];
      j = p->tail[k];

      if ((w[i] && !w[j]) || (!w[i] && w[j]))
         sum += capa[k];
   }
#if DEBUG > 0
   (void)printf("Cut Sum=%d\n", sum);
#endif
   cutsums += sum;
}
#endif

#if DEBUG > 0
static int is_valid(
   const GRAPH* p,
   const int    s,
   const int    t,
   const int*   capa,
   const int*   w)
{
   int* e;
   int* r;
   int* x;
   int j;
   int k;

   assert(p      != NULL);
   assert(p->mincut_e != NULL);
   assert(p->mincut_r != NULL);
   assert(p->mincut_x != NULL);

   e = p->mincut_e;
   r = p->mincut_r;
   x = p->mincut_x;

   for(j = 0; j < p->knots; j++)
   {
#if 0
      if ((q[j] >= 0) && (a[q[j]] != j))
         return((void)fprintf(stderr, "Queue Error 1 Knot %d\n", j), FALSE);

      if (!w[j] && (q[j] < 0) && (e[j] > 0) && (j != t))
         return((void)fprintf(stderr, "Queue Error 2 Knot %d\n", j), FALSE);

      if (!w[j] && (q[j] >= 0) && (e[j] == 0))
         return((void)fprintf(stderr, "Queue Error 3 Knot %d\n", j), FALSE);

      if (w[j] && (q[j] >= 0))
         return((void)fprintf(stderr, "Queue Error 4 Knot %d\n", j), FALSE);
#endif
      if (e[j] < 0)
         return((void)fprintf(stderr, "Negativ Execess Knot %d\n", j), FALSE);

      if (p->mincut_dist[j] >= p->knots)
         return((void)fprintf(stderr, "Distance too big Knot %d\n", j), FALSE);

      /* Extended Dormacy Property
       */
      for(k = p->outbeg[j]; k != EAT_LAST; k = p->oeat[k])
      {
         if (r[k] > 0)
         {
            if ((w[j] && !w[p->head[k]]) || (w[j] && (w[j] < w[p->head[k]])))
            {
               (void)printf("k=%d r[k]=%d head=%d tail=%d w[h]=%d w[t]=%d\n",
                  k, r[k], p->head[k], p->tail[k], w[p->head[k]], w[p->tail[k]]);

               return((void)fprintf(stderr, "Extended Dormacy Violation Knot %d\n", j), FALSE);
            }
         }
      }
   }
   for(j = 0; j < p->edges; j++)
   {
      if (r[j] < 0)
         return((void)fprintf(stderr, "Negativ Residual Edge %d\n", j), FALSE);

      if (x[j] < 0)
         return((void)fprintf(stderr, "Negativ Flow Edge %d\n", j), FALSE);

      if (r[j] + r[Edge_anti(j)] != capa[j] + capa[Edge_anti(j)])
         return((void)fprintf(stderr, "Wrong Capacity Equation Edge %d\n", j), FALSE);

      if (r[j] != capa[j] - x[j] + x[Edge_anti(j)])
         return((void)fprintf(stderr, "Wrong Residual Equation Edge %d\n", j), FALSE);
   }
   return(TRUE);
}

static void show_flow(
   const GRAPH* p,
   const int*   capa,
   const int*   w)
{
   int          j;
   int*   head;
   int*   numb;
   int*   prev;
   int*   next;
   int*   e;
   int*   x;
   int*   r;

   assert(p != NULL);
   assert(w != NULL);
   assert(p->mincut_numb   != NULL);
   assert(p->mincut_head   != NULL);
   assert(p->mincut_prev   != NULL);
   assert(p->mincut_next   != NULL);
   assert(p->mincut_e      != NULL);
   assert(p->mincut_x      != NULL);
   assert(p->mincut_r      != NULL);

   head = p->mincut_head;
   numb = p->mincut_numb;
   prev = p->mincut_prev;
   next = p->mincut_next;
   e    = p->mincut_e;
   x    = p->mincut_x;
   r    = p->mincut_r;



   (void)printf("   ");
   for(j = 0; j < p->edges; j++)
      (void)printf("%6d ", j);
   (void)fputc('\n', stdout);

   (void)printf("ta:");
   for(j = 0; j < p->edges; j++)
      (void)printf("%6d ", p->tail[j]);
   (void)fputc('\n', stdout);

   (void)printf("he:");
   for(j = 0; j < p->edges; j++)
      (void)printf("%6d ", p->head[j]);
   (void)fputc('\n', stdout);

   (void)printf("x: ");
   for(j = 0; j < p->edges; j++)
      (void)printf("%6d ", x[j]);
   (void)fputc('\n', stdout);

   (void)printf("r: ");
   for(j = 0; j < p->edges; j++)
      (void)printf("%6d ", r[j]);
   (void)fputc('\n', stdout);

   (void)printf("ca:");
   for(j = 0; j < p->edges; j++)
      (void)printf("%6d ", capa[j]);
   (void)fputc('\n', stdout);

   (void)printf("w: ");
   for(j = 0; j < p->knots; j++)
      (void)printf("%2d ", w[j]);
   (void)fputc('\n', stdout);

   (void)printf("d: ");
   for(j = 0; j < p->knots; j++)
      (void)printf("%2d ", p->mincut_dist[j]);
   (void)fputc('\n', stdout);

   (void)printf("n: ");
   for(j = 0; j < p->knots; j++)
      (void)printf("%2d ", numb[j]);
   (void)fputc('\n', stdout);

   (void)printf("h: ");
   for(j = 0; j < p->knots; j++)
      (void)printf("%2d ", head[j]);
   (void)fputc('\n', stdout);

   (void)printf("p: ");
   for(j = 0; j < p->knots; j++)
      (void)printf("%2d ", prev[j]);
   (void)fputc('\n', stdout);

   (void)printf("n: ");
   for(j = 0; j < p->knots; j++)
      (void)printf("%2d ", next[j]);
   (void)fputc('\n', stdout);

   (void)printf("e: ");
   for(j = 0; j < p->knots; j++)
      (void)printf("%2d ", e[j]);
   (void)fputc('\n', stdout);
}

#endif /* DEBUG > 0 */
