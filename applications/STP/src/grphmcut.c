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

/**@file   grphmcut.c
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
#define CHECK        0

#ifndef RESTRICT
#define RESTRICT restrict
#endif

#ifdef NDEBUG
#undef STATIST
#undef DEBUG
#define STATIST      0
#define DEBUG        0
#endif

#define Q_NULL     -1         /* NULL element of queue/list */
#define GLOBALRELABEL_ADD  10 /* constant for global relabel heuristic */
#define GLOBALRELABEL_MULT 10 /* constant for global relabel heuristic */

/** remove first element from active (singly-linked) list */
#define listdelete(node,nodedist,next,headactive)\
{\
   assert(node >= 0);\
   assert(nodedist >= 0);\
   assert(headactive != NULL);\
   headactive[nodedist] = next[node];\
}

/** add element to active (singly-linked) list */
#define listinsert(node,nodedist,next,headactive,actmin,actmax)\
{\
   assert(node >= 0);\
   assert(next != NULL);\
   assert(headactive != NULL);\
   assert(nodedist >= 0);\
   next[node] = headactive[nodedist];\
   headactive[nodedist] = node;\
   if( nodedist < (actmin) )\
      actmin = nodedist;\
   if( nodedist > (actmax) )\
      actmax = nodedist;\
}
#if 0
/** add element to active (singly-linked) list */
#define listinsert(node,nodedist,next,headactive,actmin,actmax)\
{\
   assert(node >= 0);\
   assert(next != NULL);\
   assert(headactive != NULL);\
   assert(nodedist >= 0);\
   next[node] = headactive[nodedist];\
   headactive[nodedist] = node;\
   if( nodedist < (actmin) )\
      actmin = nodedist;\
   if( nodedist > (actmax) )\
      actmax = nodedist;\
}
#endif


/** remove first element from active (singly-linked) list */
#define activedelete(node,nodedist,next,headactive)\
{\
   assert(node >= 0);\
   assert(nodedist >= 0);\
   assert(headactive != NULL);\
   headactive[nodedist] = next[node];\
}

/** add element to active (singly-linked) list */
#define activeinsert(node,nodedist,next,headactive,actmin,actmax,glbmax)\
{\
   assert(node >= 0);\
   assert(next != NULL);\
   assert(headactive != NULL);\
   assert(nodedist >= 0);\
   next[node] = headactive[nodedist];\
   headactive[nodedist] = node;\
   if( nodedist < (actmin) )\
      actmin = nodedist;\
   if( nodedist > (actmax) )\
      actmax = nodedist;\
   if( (actmax) > (glbmax) )\
      glbmax = (actmax);\
}

/** remove element from inactive (doubly-linked) list */
#define inactivedelete(node,nodedist,next,prev,headinactive,nextnode,prevnode)\
{\
   assert(node >= 0);\
   assert(nodedist >= 0);\
   assert(next != NULL);\
   assert(prev != NULL);\
   assert(headinactive != NULL);\
   nextnode = next[node];\
   if( headinactive[nodedist] == node )\
   {\
      headinactive[nodedist] = nextnode;\
      if( nextnode >= 0 )\
         prev[nextnode] = Q_NULL;\
   }\
   else\
   {\
      prevnode = prev[node];\
      assert(prevnode >= 0);\
      assert(next[prevnode] == node);\
      next[prevnode] = nextnode;\
      if( nextnode >= 0 )\
         prev[nextnode] = prevnode;\
   }\
}

/** add element to inactive (doubly-linked) list */
#define inactiveinsert(node,nodedist,next,prev,headinactive,nextnode)\
{\
   assert(node >= 0);\
   assert(nodedist >= 0);\
   assert(next != NULL);\
   assert(prev != NULL);\
   assert(headinactive != NULL);\
   nextnode = headinactive[nodedist];\
   next[node] = nextnode;\
   prev[node] = Q_NULL;\
   if( nextnode >= 0 )\
      prev[nextnode] = node;\
   headinactive[nodedist] = node;\
}


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


#if 1

#if CHECK
/** checker */
static int is_valid_arr(
   const GRAPH* p,
   const int    s,
   const int    t,
   const int*   startarr,
   const int*   headarr,
   const int*   revedgearr,
   const int*   resarr,
   const int*   capa,
   const int*   w)
{
   int* e;
   int* r;
   int* dist;

   int j;
   int k;
   int end;
   int head;

   assert(p      != NULL);
   assert(p->mincut_e != NULL);
   assert(p->mincut_r != NULL);

   e = p->mincut_e;
   r = p->mincut_r;
   dist = p->mincut_dist;

   for(j = 0; j < p->knots; j++)
   {
      if( e[j] < 0 )
         return ((void) fprintf(stderr, "Negativ Execess Node %d\n", j), FALSE);

      if( dist[j] >= p->knots )
         return ((void) fprintf(stderr, "Distance too big Node %d dist: %d\n", j,  dist[j]), FALSE);

      /* extended dormacy property */
      for( k = startarr[j], end = startarr[j + 1]; k != end; k++ )
      {
         head = headarr[k];

         if( r[k] > 0 )
         {
            if( dist[j] > dist[head] + 1 && w[j] == 0 && w[head] == 0 )
            {
               printf("distance fail! %d->%d (%d)\n", j, head, k);
               return FALSE;
            }

            if( (w[j] && (w[j] < w[head])) )
            {
               printf("Extended Dormacy Violation %d %d \n", j, head);
            }

            if( w[j] && !w[head] )
            {
               printf("Simple Dormacy Violation %d %d \n", j, head);
            }

            if( (w[j] && !w[head]) || (w[j] && (w[j] < w[head])) )
            {
               (void) printf("e=%d r[e]=%d head=%d tail=%d w[h]=%d w[t]=%d\n",
                     k, r[k], head, p->tail[k], w[head],
                     w[p->tail[k]]);

               return ((void) fprintf(stderr,
                     "Dormacy Violation Knot %d\n", j), FALSE);
            }
         }

         if( r[k] < 0 )
            return ((void) fprintf(stderr, "Negativ Residual Edge %d\n", k), FALSE);
      }
   }

   return(TRUE);
}
#endif

/** initialize min cut arrays */
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

   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_dist), p->knots + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_head), p->knots + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_head_inact), p->knots + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_numb), p->knots + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_prev), p->knots + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_next), p->knots + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_temp), p->knots + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_e), p->knots + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_r), p->edges) );

   return SCIP_OKAY;
}

/** frees min cut arrays */
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
   assert(p->mincut_r    != NULL);

   SCIPfreeMemoryArray(scip, &(p->mincut_r));
   SCIPfreeMemoryArray(scip, &(p->mincut_e));
   SCIPfreeMemoryArray(scip, &(p->mincut_temp));
   SCIPfreeMemoryArray(scip, &(p->mincut_next));
   SCIPfreeMemoryArray(scip, &(p->mincut_prev));
   SCIPfreeMemoryArray(scip, &(p->mincut_numb));
   SCIPfreeMemoryArray(scip, &(p->mincut_head_inact));
   SCIPfreeMemoryArray(scip, &(p->mincut_head));
   SCIPfreeMemoryArray(scip, &(p->mincut_dist));

#if STATIST
   cut_statist();
#endif
}
#if 0
/** backwards bfs from t along arcs of cap > 0 and set distance labels according to distance from t */
static int bfs(
   const GRAPH* p,
   const int    s,
   const int    t,
   int*         dist,
   int*         temp,
   int*         w,
   const int*   capa,
   const int*   edgestart,
   const int*   edgearr,
   const int*   headarr
   )
{
   int          i;
   int          j;
   int          k;
   int          l;
   int          end;
   int          nnodes;
   int          dplus1;
   int          visited = 0;

   assert(temp != NULL);
   assert(dist != NULL);
   assert(w      != NULL);
   assert(s      >= 0);
   assert(s      < p->knots);
   assert(t      >= 0);
   assert(t      < p->knots);

   /* initialize */
   temp[visited++] = t;
   dist[t]         = 0;

   nnodes = p->knots;

   /* bfs loop */
   for( j = 0; j < visited; j++ )
   {
      assert(visited <= nnodes);

      i = temp[j];
      dplus1 = dist[i] + 1;

      assert(i         >= 0);
      assert(i         <  nnodes);
      assert(dist[i] >= 0);
      assert(dist[i] <  visited);
      assert(w[i]      == 0);

      for( l = edgestart[i], end = edgestart[i + 1]; l != end; l++ )
      {
         k = headarr[l];

         /* not visited yet? */
         if( dist[k] < 0 && w[k] == 0 )
         {
            assert(w[k] == 0);

            dist[k] = dplus1;
            temp[visited++] = k;
            assert(dist[k] < nnodes);
         }
      }
   } /* bfs loop */

   return (visited);
}
#endif

/** global relabel heuristic that sets distance of sink to zero and relabels all other nodes using backward bfs on residual
 * graph, starting from the sink.  */
static void globalrelabel(
   const GRAPH* p,
   const int    s,
   const int    t,
   int* RESTRICT dist,
   int* RESTRICT headactive,
   int* RESTRICT headinactive,
   int* RESTRICT edgecurr,
   int* RESTRICT next,
   int* RESTRICT prev,
   int* RESTRICT excess,
   int* RESTRICT residual,
   int* RESTRICT w,
   const int*   edgestart,
   const int*   edgearr,
   const int*   headarr,
   int*         actmin,
   int*         actmax,
   int*         glbmax,
   int*         dormmax
   )
{
   int i;
   int q;
   int k;
   int end;
   int nnodes;
   int currdist;
   int nextdist;
   int nextnode;
   int actmaxloc;
   int actminloc;
   int glbmaxloc;
   int dormmaxnext;
   SCIP_Bool hit;

   assert(p      != NULL);
   assert(s      >= 0);
   assert(s      <  p->knots);
   assert(t      >= 0);
   assert(t      <  p->knots);
   assert(s      != t);
   assert(w      != NULL);
   assert(headactive != NULL);
   assert(headinactive != NULL);
   assert(edgecurr != NULL);
   assert(prev != NULL);
   assert(next != NULL);
   assert(residual != NULL);
   assert(excess != NULL);
   assert(edgearr != NULL);
   assert(headarr != NULL);
   assert(edgecurr != NULL);
   assert(edgestart != NULL);
   assert(headactive != NULL);
   assert(headinactive != NULL);

   nnodes = p->knots;

   assert(*glbmax < nnodes);

   end = *glbmax;

   assert(w[t] == 0);

   for( i = 0; i <= end; i++ )
   {
      headactive[i] = Q_NULL;
      headinactive[i] = Q_NULL;
      if( w[i] == 0 )
         w[i] = -1;
   }

   for( i = 0; i < nnodes; i++ )
      assert(headactive[i] == Q_NULL && headinactive[i] == Q_NULL);

   for( i = end + 1; i < nnodes; i++ )
      if( w[i] == 0 )
         w[i] = -1;

   assert(w[t] == -1);

   w[t] = 0;
   dist[t] = 0;

   glbmaxloc = 0;
   actmaxloc = 0;
   actminloc = nnodes;

   /* add t to inactive list */
   inactiveinsert(t, 0, next, prev, headinactive, nextnode);

   /* bfs loop */
   for( currdist = 0; TRUE; currdist++ )
   {
      nextdist = currdist + 1;

      /* no more nodes with current distance? */
      if( (headactive[currdist] < 0) && (headinactive[currdist] < 0) )
         break;

      /* check inactive nodes */
      for( i = headinactive[currdist]; i >= 0; i = next[i] )
      {
         for( q = edgestart[i], end = edgestart[i + 1]; q != end; q++ )
         {
            if( residual[edgearr[q]] != 0 )
            {
               if( w[headarr[q]] < 0 )
               {
                  k = headarr[q];

                  w[k] = 0;

                  dist[k] = nextdist;
                  edgecurr[k] = edgestart[k];

                  if( excess[k] > 0 )
                  {
                     activeinsert(k, nextdist, next, headactive, actminloc, actmaxloc, glbmaxloc);
                  }
                  else
                  {
                     inactiveinsert(k, nextdist, next, prev, headinactive, nextnode);
                  }
               }
            }
         }
      }

      /* check active nodes */
      for( i = headactive[currdist]; i >= 0; i = next[i] )
      {
         for( q = edgestart[i], end = edgestart[i + 1]; q != end; q++ )
         {
            if( residual[edgearr[q]] != 0 )
            {
               k = headarr[q];

               if( w[k] < 0 )
               {
                  w[k] = 0;

                  dist[k] = nextdist;
                  edgecurr[k] = edgestart[k];

                  if( excess[k] > 0 )
                  {
                     activeinsert(k, nextdist, next, headactive, actminloc, actmaxloc, glbmaxloc);
                  }
                  else
                  {
                     inactiveinsert(k, nextdist, next, prev, headinactive, nextnode);
                  }
               }
            }
         }
      }
   }

   assert((headactive[currdist] == Q_NULL) && (headinactive[currdist] == Q_NULL));
   assert(currdist == 0 || (headactive[currdist - 1] != Q_NULL) || (headinactive[currdist - 1] != Q_NULL));

   /* set global maximum label */
   if( (currdist - 1) > (glbmaxloc) )
      glbmaxloc = currdist - 1;

   hit = FALSE;

   assert(w[t] == 0);

   dormmaxnext = (*dormmax) + 1;
   for( i = 0; i < nnodes; i++ )
   {
      if( w[i] < 0 )
      {
         w[i] = dormmaxnext;
         hit = TRUE;
      }
   }
   if( hit )
      *dormmax = dormmaxnext;

   *glbmax = glbmaxloc;
   *actmin = actminloc;
   *actmax = actmaxloc;

   return;
}


/** initialize data for the first max-flow run */
static void initialise(
   const GRAPH* p,
   const int    s,
   const int    t,
   const int    rootcutsize,
   const int*   rootcut,
   const int*   capa,
   int* RESTRICT dist,
   int* RESTRICT headactive,
   int* RESTRICT headinactive,
   int* RESTRICT edgecurr,
   int* RESTRICT next,
   int* RESTRICT prev,
   int* RESTRICT temp,
   int* RESTRICT excess,
   int* RESTRICT residual,
   int* RESTRICT w,
   const int*   edgestart,
   const int*   edgearr,
   const int*   headarr,
   int*         dormmax,
   int*         actmin,
   int*         actmax,
   int*         glbmax
   )
{
   int i;
#if 0
   int k;
   int e;
   int q;
   int end;
   int head;
   int nedges;
#endif
   int nnodes;
   int nextnode;
   int actmaxloc;
   int actminloc;
   int glbmaxloc;

   assert(p      != NULL);
   assert(s      >= 0);
   assert(s      <  p->knots);
   assert(t      >= 0);
   assert(t      <  p->knots);
   assert(s      != t);
   assert(capa   != NULL);
   assert(w      != NULL);
   assert(headactive != NULL);
   assert(headinactive != NULL);
   assert(prev != NULL);
   assert(next != NULL);
   assert(temp != NULL);
   assert(residual != NULL);
   assert(excess != NULL);
   assert(edgearr != NULL);
   assert(headarr != NULL);
   assert(edgecurr != NULL);
   assert(edgestart != NULL);
   assert(headactive != NULL);
   assert(headinactive != NULL);

   nnodes = p->knots;
#if 0
   nedges = p->edges;
#endif

   actmaxloc = 0;
   actminloc = nnodes;
   glbmaxloc = 0;

   *dormmax = 1;

   assert(w[s] > 0);
   assert(w[t] == 0);

   /* set distance labels, and add nodes to lists */
   for( i = nnodes - 1; i >= 0; i-- )
   {
      dist[i] = 1;

      /* in root cut? */
      if( w[i] > 0 )
         continue;

      /* sink? */
      if( i == t )
      {
         dist[i] = 0;
         inactiveinsert(i, 0, next, prev, headinactive, nextnode);
         continue;
      }

      assert(w[i] == 0);

      if( excess[i] > 0 )
      {
         /* add to active list */
         activeinsert(i, 1, next, headactive, actminloc, actmaxloc, glbmaxloc);
      }
      else
      {
         /* add to inactive list */
         inactiveinsert(i, 1, next, prev, headinactive, nextnode);
      }
   }

   glbmaxloc = 1;

#if 0 /* both */
   assert(nnodes <= p->edges);

   for( i = 0; i < nnodes; i++ )
   {
      dist[i] = -1;
      excess[i] = 0;
      edgecurr[i] = edgestart[i];
      headactive[i] = Q_NULL;
      headinactive[i] = Q_NULL;

      residual[i] = capa[i];
   }

   for( e = nnodes; e < nedges; e++ )
      residual[e] = capa[e];
#endif


#if 0 /* simple */

   /* push from source */
   for( q = edgestart[s], end = edgestart[s + 1]; q != end; q++ )
   {
      e = edgearr[q];

      if( residual[e] == 0 )
         continue;

      i = headarr[q];

      assert(w[i] == 0);
      assert(residual[e] == capa[e]);

      residual[flipedge(e)] += residual[e];
      excess[i] += residual[e];
      residual[e] = 0; /* -= residual[e] */
   }

   /* set distance labels, and add nodes to lists */
   for( i = 0; i < nnodes; i++ )
   {
      /* sink? */
      if( i == t )
      {
         dist[i] = 0;
         inactiveinsert(i, 0, next, prev, headinactive, nextnode);
         continue;
      }

      dist[i] = 1;

      /* source? */
      if( i == s )
      {
         assert(excess[s] == 0);

         /* put s (source) into lowest dormant set */
         w[i] = 1;    /* *dormmax */
         continue;
      }

      assert(w[i] == 0);

      if( excess[i] > 0 )
      {
         /* add to active list */
         activeinsert(i, 1, next, headactive, actminloc, actmaxloc, glbmaxloc);
      }
      else
      {
         /* add to inactive list */
         inactiveinsert(i, 1, next, prev, headinactive, nextnode);
      }
   }
   glbmaxloc = 1;
#endif

#if 0 /* with initial global relabel */
   /* backward bfs from t */
   (void)bfs(p, s, t, dist, temp, w, capa, edgestart, edgearr, headarr);

   /* put unreachable nodes to sleep */
   for( i = 0; i < nnodes; i++ )
      if( dist[i] < 0 )
      {
         if( (*dormmax) == 1 )
            *dormmax = 2;
         w[i] = 2;
      }

   /* put s (source) into lowest dormant set */
   w[s] = 1;    /* *dormmax */

   /* s not reached? */
   if( dist[s] < 0 )
      return;

   assert(w[s]      >  0);
   assert(dist[s] >  0);
   assert(w[t]      == 0);
   assert(dist[t] == 0);

   /* push from source */
   for( q = edgestart[s], end = edgestart[s + 1]; q != end; q++ )
   {
      e = edgearr[q];

      if( residual[e] == 0 )
         continue;

      i = headarr[q];

      assert(residual[e] == capa[e]);

      residual[flipedge(e)] += residual[e];
      excess[i] += residual[e];
      residual[e] = 0; /* -= residual[e] */
   }

   w[t] = 1;
   inactiveinsert(t, 0, next, prev, headinactive, nextnode);

   for( i = 0; i < nnodes; i++ )
   {
      /* including t and s */
      if( w[i] )
         continue;

      assert(w[i] == 0);
      assert(s != i);

      if( excess[i] > 0 )
      {
         /* add to active list */
         activeinsert(i, dist[i], next, headactive, actminloc, actmaxloc, glbmaxloc);
      }
      else
      {
         if( dist[i] > (glbmaxloc) )
          glbmaxloc = dist[i];

         /* add to inactive list */
         inactiveinsert(i, dist[i], next, prev, headinactive, nextnode);
      }
   }

   w[t] = 0;
#endif

   *glbmax = glbmaxloc;
   *actmax = actmaxloc;
   *actmin = actminloc;

#if DEBUG > 0
   assert(is_valid(p, s, t, capa, w));
#endif
}

/** initialize data for the repeated max-flow run */
static void reinitialise(
   const GRAPH* p,
   const int    s,
   const int    t,
   const int    rootcutsize,
   const int*   rootcut,
   const int*   capa,
   int* RESTRICT dist,
   int* RESTRICT headactive,
   int* RESTRICT headinactive,
   int* RESTRICT edgecurr,
   int* RESTRICT next,
   int* RESTRICT prev,
   int* RESTRICT temp,
   int* RESTRICT excess,
   int* RESTRICT residual,
   int* RESTRICT w,
   const int*   edgestart,
   const int*   edgearr,
   const int*   headarr,
   int*         dormmax,
   int*         actmin,
   int*         actmax,
   int*         glbmax
   )
{
   int i;
   int j;
   int k;
   int l;
   int end;
   int nnodes;
   int visited;
#if 0
   int a;
   int headnode;
#endif
   int nextnode;
   int actmaxloc;
   int actminloc;
   int glbmaxloc;
   int dormmaxloc;
   int dormmaxlocp1;
   SCIP_Bool hit;

   assert(p      != NULL);
   assert(s      >= 0);
   assert(s      <  p->knots);
   assert(t      >= 0);
   assert(t      <  p->knots);
   assert(s      != t);
   assert(capa   != NULL);
   assert(w      != NULL);
   assert(dist != NULL);
   assert(prev != NULL);
   assert(next != NULL);
   assert(temp != NULL);
   assert(residual != NULL);
   assert(excess != NULL);
   assert(edgearr != NULL);
   assert(headarr != NULL);
   assert(edgecurr != NULL);
   assert(edgestart != NULL);
   assert(headactive != NULL);
   assert(headinactive != NULL);

   /* initialize */
   nnodes = p->knots;

   assert(w[s] == 1);

   dormmaxloc = 1;

   actmaxloc = 0;
   actminloc = nnodes;
   glbmaxloc = 0;

   /* t already awake? */
   if( w[t] == 0 )
   {
      for( i = nnodes - 1; i >= 0; i-- )
      {
         headactive[i] = Q_NULL;
         headinactive[i] = Q_NULL;

         /* reset distance label for awake nodes and adapt incident edges */
         if( w[i] == 0 )
         {
            dist[i] = -1;
            edgecurr[i] = edgestart[i];
         }
         else if( w[i] > (dormmaxloc) )
         {
            dormmaxloc = w[i];
         }
      }
   }
   else
   {
      int wt = w[t];
      for( i = nnodes - 1; i >= 0; i-- )
      {
         headactive[i] = Q_NULL;
         headinactive[i] = Q_NULL;

         /* wake up nodes in higher or equal dormant layer */
         if( (w[i] == 0) || (w[i] >= wt) )
         {
            w[i] = 0;
            dist[i] = -1;
            edgecurr[i] = edgestart[i];
         }
         else if( w[i] > (dormmaxloc) )
         {
            dormmaxloc = w[i];
         }
      }
   }

   /* backward bfs from t */
   visited = 0;

   inactiveinsert(t, 0, next, prev, headinactive, nextnode);

   /* initialize */
   temp[visited++] = t;
   dist[t]         = 0;

   /* bfs loop */
   for( j = 0; j < visited; j++ )
   {
      assert(visited <= nnodes);

      i = temp[j];

      assert(i         >= 0);
      assert(i         <  nnodes);
      assert(dist[i] >= 0);
      assert(dist[i] <  visited);
      assert(w[i]      == 0);

      for( l = edgestart[i], end = edgestart[i + 1]; l != end; l++ )
      {
         k = headarr[l];

         /* not visited yet? */
         if( dist[k] < 0 && w[k] == 0 )
         {
            if( residual[edgearr[l]] > 0 )
            {
               dist[k] = dist[i] + 1;
               temp[visited++] = k;
               assert(dist[k] < nnodes);

               if( excess[k] > 0 )
               {
                  /* add to active list */
                  activeinsert(k, dist[k], next, headactive, actminloc, actmaxloc, glbmaxloc);
               }
               else
               {
                  if( dist[k] > (glbmaxloc) )
                     glbmaxloc = dist[k];

                  /* add to inactive list */
                  inactiveinsert(k, dist[k], next, prev, headinactive, nextnode);
               }
            }
         }
      }
   } /* bfs loop */

   hit = FALSE;
   dormmaxlocp1 = dormmaxloc + 1;
   for( i = nnodes - 1; i >= 0; i-- )
   {
      /* unreachable non-dormant node? */
      if( dist[i] < 0 && w[i] == 0 )
      {
         dist[i] = -1;
         w[i] = dormmaxlocp1;
         hit = TRUE;
      }
   }

   if( hit )
      dormmaxloc++;
#if 0
   for( i = nnodes - 1; i >= 0; i-- )
   {
      /* is the node awake? */
      if( w[i] == 0 )
      {
         /*
          * update excess and residual weights
          * */

         excess[i] = 0;

         /* iterate over outgoing arcs of i */
         for( j = edgestart[i], end = edgestart[i + 1]; j != end; j++ )
         {
            headnode = headarr[j];
            a = edgearr[j];

            if( w[headnode] )
            {
               assert(residual[flipedge(a)] == 0);
               excess[i] += capa[flipedge(a)];
            }
            else
            {
               residual[a] = capa[a];
            }
         }
      }
   }

   assert(w[t] == 0);
   assert(dist[t] == 0);
   assert(temp[0] == t);


   for( i = 1; i < visited && 0; i++ )
   {
      k = temp[i];
      assert(w[k] == 0);
      assert(k    != s);
      assert(k    != t);

      if( excess[k] > 0 )
      {
         /* add to active list */
         activeinsert(k, dist[k], next, headactive, actminloc, actmaxloc, glbmaxloc);
      }
      else
      {
         if( dist[k] > (glbmaxloc) )
            glbmaxloc = dist[k];

         /* add to inactive list */
         inactiveinsert(k, dist[k], next, prev, headinactive, nextnode);
      }
   }
#endif
   *dormmax = dormmaxloc;
   *actmin = actminloc;
   *actmax = actmaxloc;
   *glbmax = glbmaxloc;

#if DEBUG > 0
   assert(is_valid(p, s, t, capa, w));
#endif
}

/** finds a minimum s-t cut */
void graph_mincut_exec(
   const GRAPH*          p,
   const int             s,
   const int             t,
   const int             nnodesreal,
   const int             nedgesreal,
   const int             rootcutsize,
   const int*            rootcut,
   const int*            capa,
   int* RESTRICT         w,
   const int*            edgestart,
   const int*            edgearr,
   const int*            headarr,
   const SCIP_Bool       rerun
   )
{
   int    i;
   int    l;
   int    actmax;
   int    actmin;
   int    glbmax;
   int    dminus1;
   int    dormmax;
   int    mindist;
   int    minnode;
   int    minedgestart;
   int*   e;
   int*   r;
   int*   dist;
   int*   prev;
   int*   next;
   int*   temp;
   int*   edgecurr;
   int*   headactive;
   int*   headinactive;
   int j;
   int end;
   int nnodes;
   int headnode;
   int prevnode;
   int nextnode;
   int maxresval;
   int maxresedge;
   int maxresnode;
   int relabeltrigger;
   int relabelupdatebnd;

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
   assert(p->mincut_r      != NULL);

   dist = p->mincut_dist;
   headactive = p->mincut_head;
   headinactive = p->mincut_head_inact;
   edgecurr = p->mincut_numb;
   prev = p->mincut_prev;
   next = p->mincut_next;
   temp = p->mincut_temp;
   e    = p->mincut_e;
   r    = p->mincut_r;

   nnodes = p->knots;
   minnode = -1;
   maxresedge = -1;

   relabelupdatebnd = (GLOBALRELABEL_MULT * nnodesreal) + nedgesreal;
   relabeltrigger = 0;

   if( !rerun )
      initialise(p, s, t, rootcutsize, rootcut, capa, dist, headactive, headinactive, edgecurr, next, prev, temp, e, r, w, edgestart,
            edgearr, headarr, &dormmax, &actmin, &actmax, &glbmax);
   else
      reinitialise(p, s, t, rootcutsize, rootcut, capa, dist, headactive, headinactive, edgecurr, next, prev, temp, e, r, w, edgestart,
            edgearr, headarr, &dormmax, &actmin, &actmax, &glbmax);

   /* main loop: get highest label node */
   while( actmax >= actmin )
   {
      /* no active vertices with distance label == actmax? */
      if( headactive[actmax] < 0 )
      {
         actmax--;
         continue;
      }

      i = headactive[actmax];

      assert(actmax <= glbmax);
      assert(dist[i] == actmax);
      assert(w[i]    == 0);
      assert(i       != t);
      assert(i       != s);
      assert(e[i] > 0);

      /* delete i from active list */
      activedelete(i,actmax,next,headactive);

      end = edgestart[i + 1];

      assert(end > edgestart[i]);

      /* try to discharge i (repeatedly) */
      for( ;; )
      {
         /* iterate over outgoing arcs of i */
         for( j = edgecurr[i]; j != end; j++ )
         {
            /* non-saturated edge? */
            if( r[j] != 0 )
            {
               /* non-dormant head node? */
               if( w[headarr[j]] == 0 )
               {
                  headnode = headarr[j];

                  assert(w[headnode] == 0);
                  assert(headnode != s);

                  dminus1 = dist[i] - 1;

                  /* admissible arc? */
                  if( dist[headnode] == dminus1 )
                  {
                     assert(Min(e[i], r[j]) > 0);

                     /* is head node now active? */
                     if( e[headnode] == 0 )
                     {
                        if( headnode != t )
                        {
                           /* remove head node from inactive list */
                           inactivedelete(headnode, dminus1, next, prev, headinactive, nextnode, prevnode);

                           /* add head node to active list */
                           activeinsert(headnode, dminus1, next, headactive, actmin, actmax, glbmax);
                        }
                     }

                     /* not more residual capacity than excess? */
                     if( r[j] < e[i] )
                     {
                        e[i] -= r[j];
                        e[headnode] += r[j];
                        r[edgearr[j]] += r[j];
                        r[j] = 0; /* -= r[a] */

                        assert(e[i] > 0);
                     }
                     /* more residual capacity than excess */
                     else
                     {
                        r[edgearr[j]] += e[i];
                        r[j] -= e[i];
                        e[headnode] += e[i];
                        e[i] = 0;          /* -= e[i] */

                        /* excess vanished, so stop discharging */
                        break;
                     }
                  } /* admissible arc */
               } /*  non-dormant head node */
            } /* non-saturated edge */
         } /* all outgoing arcs */

         /* is there still excess on i? */
         if( j == end )
         {
            assert(e[i] > 0);

            /*
             * relabel i
             * */

            /* i only node of distance dist[i]? */
            if( (headactive[dist[i]] < 0) && (headinactive[dist[i]] < 0) )
            {
               /* put i into new dormant set */
               w[i] = ++dormmax;

               assert(dormmax <= nnodes);
               assert(dist[i] < nnodes);

               /* remove nodes of distance label >= dist[i] */
               for( j = dist[i] + 1; j <= glbmax; j++ )
               {
                  assert(headactive[j] == Q_NULL);

                  for( l = headinactive[j]; l >= 0; l = next[l] )
                  {
                     if( w[l] == 0 )
                        w[l] = dormmax;
                  }
                  headinactive[j] = Q_NULL;
               }

               actmax = dist[i] - 1;
               glbmax = actmax;

               break;
            }

            j = edgestart[i];

            assert(end > j);
            assert(end == edgestart[i + 1]);

            relabeltrigger += GLOBALRELABEL_ADD + (end - j);

            mindist = nnodes;
            maxresval = 0;

            /* todo only for debugging, could be deleted */
            maxresnode = -1;
            minedgestart = -1;

            /* find the first (!) minimum */
            for( ; j != end; j++ )
            {
               /* useable edge? */
               if( r[j] != 0 )
               {
                  /* non-dormant node? */
                  if( w[headarr[j]] == 0 )
                  {
#if 0
                     if( dist[headarr[j]] < mindist )
                     {
                        mindist = dist[headarr[j]];
                        minedgestart = j;
                     }
#endif
                     if( dist[headarr[j]] <= mindist )
                     {
                        if( dist[headarr[j]] < mindist )
                        {
                           mindist = dist[headarr[j]];
                           minedgestart = j;
                           maxresval = r[j];
                           maxresedge = j;
                           minnode =  headarr[j];
                           maxresnode = headarr[j];
                        }
                        else if( r[j] > maxresval )
                        {
                           maxresval = r[j];
                           maxresedge = j;
                           maxresnode = headarr[j];
                        }
                     }
                  }
               }
            }

            if( (++mindist) < nnodes )
            {
               assert(minedgestart >= 0);
               assert(r[minedgestart] > 0);

               dist[i] = mindist;
               //edgecurr[i] = minedgestart; // iff: #if 1

               if( glbmax < mindist )
                  glbmax = mindist;

               assert(maxresnode >= 0);
               assert(maxresnode == headarr[maxresedge]);
#if 1
               /* can we completely discharge i already? */
               if( r[minedgestart] >= e[i] )
               {
                  /* is node now active? */
                  if( e[minnode] == 0 )
                  {
                     if( minnode != t )
                     {
                        inactivedelete(minnode, dist[minnode], next, prev, headinactive, nextnode, prevnode);
                        activeinsert(minnode, dist[minnode], next, headactive, actmin, actmax, glbmax);
                     }
                  }

                  r[edgearr[minedgestart]] += e[i];
                  r[minedgestart] -= e[i];
                  e[minnode] += e[i];
                  e[i] = 0; /* -= e[i] */

                  inactiveinsert(i, mindist, next, prev, headinactive, nextnode);

                  edgecurr[i] = minedgestart;

                  /* excess vanished, so stop discharging */
                  break;
               }

               /* is node now active? */
               if( e[minnode] == 0 )
               {
                  if( minnode != t )
                  {
                     inactivedelete(minnode, dist[minnode], next, prev, headinactive, nextnode, prevnode);
                     activeinsert(minnode, dist[minnode], next, headactive, actmin, actmax, glbmax);
                  }
               }

               e[i] -= r[minedgestart];
               e[minnode] += r[minedgestart];
               r[edgearr[minedgestart]] += r[minedgestart];
               r[minedgestart] = 0; /* -= r[minedgestart] */

               assert(e[i] > 0);

               if( maxresval >= e[i] && minedgestart != maxresedge ) // todo
               {
                  /* is node now active? */
                  if( e[maxresnode] == 0 )
                  {
                     if( maxresnode != t )
                     {
                        /* remove head node from inactive list */
                        inactivedelete(maxresnode, dist[maxresnode], next, prev, headinactive, nextnode, prevnode);

                        /* add head node to active list */
                        activeinsert(maxresnode, dist[maxresnode], next, headactive, actmin, actmax, glbmax);
                     }
                  }

                  r[edgearr[maxresedge]] += e[i];
                  r[maxresedge] -= e[i];
                  e[maxresnode] += e[i];
                  e[i] = 0; /* -= e[i] */

                  inactiveinsert(i, mindist, next, prev, headinactive, nextnode);

                  edgecurr[i] = minedgestart;

                  /* excess vanished, so stop discharging */
                  break;
               }
               edgecurr[i] = minedgestart + 1;
#endif

               continue;
            }
            /* could not relabel i */
            else
            {
               /* put i into new dormant set */
               w[i] = ++dormmax;

               assert(dormmax <= nnodes);

               assert(!((headactive[dist[i]] == Q_NULL) && (headinactive[dist[i]] == Q_NULL)));

               if( relabeltrigger > relabelupdatebnd )
               {
                  if( actmax < actmin )
                     break;

                  /* execute global relabel heuristic */
                  globalrelabel(p, s, t, dist, headactive, headinactive, edgecurr, next, prev,
                        e, r, w, edgestart, edgearr, headarr, &actmin, &actmax, &glbmax, &dormmax);

                  relabeltrigger = 0;
               }

               break;
            }
         }
         /* no excess on i */
         else
         {
            assert(e[i] == 0);

            edgecurr[i] = j;

            inactiveinsert(i, dist[i], next, prev, headinactive, nextnode);

            break;
         }
      }

      if( relabeltrigger > relabelupdatebnd )
      {
         if( actmax < actmin )
            break;

         /* execute global relabel heuristic */
         globalrelabel(p, s, t, dist, headactive, headinactive, edgecurr, next, prev,
               e, r, w, edgestart, edgearr, headarr, &actmin, &actmax, &glbmax, &dormmax);

         relabeltrigger = 0;
      }
   }
   assert(w[s]);
   assert(!w[t]);

#if CHECK
   if( !is_valid_arr(p, s, t, edgestart, headarr, edgearr, r, capa, w) )
      printf("flow is not valid \n");
   assert(is_valid_arr(p, s, t, edgestart, headarr, edgearr, r, capa, w));

#endif


#if STATIST
   cut_sum(p, capa, w);
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
   int* dist;
#if 0
   int* x;
#endif
   int j;
   int k;

   assert(p      != NULL);
   assert(p->mincut_e != NULL);
   assert(p->mincut_r != NULL);
   assert(p->mincut_x != NULL);

   e = p->mincut_e;
   r = p->mincut_r;
   dist = p->mincut_dist;
#if 0
   x = p->mincut_x;
#endif

#if 0
   for (j = 0; j < p->knots; j++)
   {
      for (k = p->outbeg[j]; k != EAT_LAST; k = p->oeat[k])
      {
         if( r[k] > 0 )
         {
            if( dist[j] > dist[p->head[k]] + 1 )
            {
               printf("distance fail! %d->%d (%d)\n", j, p->head[k], k);
               return FALSE;
            }
         }
      }
   }
#endif

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
#if 0
      if (p->mincut_dist[j] >= p->knots)
         return((void)fprintf(stderr, "Distance too big Knot %d\n", j), FALSE);
#endif
      /* Extended Dormacy Property
       */
      for(k = p->outbeg[j]; k != EAT_LAST; k = p->oeat[k])
      {
         if (r[k] > 0)
         {
#if 1
            if( dist[j] > dist[p->head[k]] + 1 && w[j] == 0 &&  w[p->head[k]] == 0 )
            {
               printf("distance fail! %d->%d (%d)\n", j, p->head[k], k);
               return FALSE;
            }
#endif

            if((w[j] && (w[j] < w[p->head[k]])))
            {
               printf("fail %d %d\n", j, p->head[k]);
            }

            if( w[j] && !w[p->head[k]] )
            {
               printf("fail2 %d -> %d\n", j, p->head[k]);
               printf("w: %d -> w: %d \n", w[j], w[p->head[k]]);
               printf("edge %d \n", k);
            }

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

      if (r[j] + r[Edge_anti(j)] != capa[j] + capa[Edge_anti(j)] && w[j] == 0 && w[p->head[k]] == 0 )
         return((void)fprintf(stderr, "Wrong Capacity Equation Edge %d\n", j), FALSE);
#if 0
      if (x[j] < 0)
         return((void)fprintf(stderr, "Negativ Flow Edge %d\n", j), FALSE);

      if (r[j] != capa[j] - x[j] + x[Edge_anti(j)])
         return((void)fprintf(stderr, "Wrong Residual Equation Edge %d\n", j), FALSE);
#endif
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






#else
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

#if BLK
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->mincut_dist), p->knots) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->mincut_head), p->knots) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->mincut_numb), p->knots) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->mincut_prev), p->knots) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->mincut_next), p->knots) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->mincut_temp), p->knots) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->mincut_e), p->knots) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->mincut_x), p->edges) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->mincut_r), p->edges) );
#else
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_dist), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_head), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_numb), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_prev), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_next), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_temp), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_e), p->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_x), p->edges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mincut_r), p->edges) );
#endif

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

#if BLK
   SCIPfreeBlockMemoryArray(scip, &(p->mincut_r), p->edges);
   SCIPfreeBlockMemoryArray(scip, &(p->mincut_x), p->edges);
   SCIPfreeBlockMemoryArray(scip, &(p->mincut_e), p->knots);
   SCIPfreeBlockMemoryArray(scip, &(p->mincut_temp), p->knots);
   SCIPfreeBlockMemoryArray(scip, &(p->mincut_next), p->knots);
   SCIPfreeBlockMemoryArray(scip, &(p->mincut_prev), p->knots);
   SCIPfreeBlockMemoryArray(scip, &(p->mincut_numb), p->knots);
   SCIPfreeBlockMemoryArray(scip, &(p->mincut_head), p->knots);
   SCIPfreeBlockMemoryArray(scip, &(p->mincut_dist), p->knots);
#else
   SCIPfreeMemoryArray(scip, &(p->mincut_r));
   SCIPfreeMemoryArray(scip, &(p->mincut_x));
   SCIPfreeMemoryArray(scip, &(p->mincut_e));
   SCIPfreeMemoryArray(scip, &(p->mincut_temp));
   SCIPfreeMemoryArray(scip, &(p->mincut_next));
   SCIPfreeMemoryArray(scip, &(p->mincut_prev));
   SCIPfreeMemoryArray(scip, &(p->mincut_numb));
   SCIPfreeMemoryArray(scip, &(p->mincut_head));
   SCIPfreeMemoryArray(scip, &(p->mincut_dist));
#endif

#if STATIST
   cut_statist();
#endif
}

#if 0
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

   if (q_next[knot] != Q_NM)
   {
      /* Etwa Erster ?
       */
      if (q_prev[knot] == Q_NULL)
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
      if (q_next[knot] != Q_NULL)
      {
         assert(q_next[knot]         >= 0);
         assert(q_prev[q_next[knot]] == knot);

         q_prev[q_next[knot]] = q_prev[knot];
      }
      q_next[knot] = Q_NM;
      q_prev[knot] = Q_NM;
   }
   assert(q_next[knot] == Q_NM);
   assert(q_prev[knot] == Q_NM);
}

inline static int insert(
   const int knot,
   int*      q_dist,
   int*      q_head,
   int*      q_prev,
   int*      q_next,
   int       dmin,
   int*      dlmax)
{
   assert(knot         >= 0);
   assert(q_dist       != NULL);
   assert(q_head       != NULL);
   assert(q_prev       != NULL);
   assert(q_next       != NULL);
   assert(q_dist[knot] >  0);
   assert(dmin         >= 1);

   if( q_prev[knot] == Q_NM )
   {
      q_prev[knot] = Q_NULL;
      q_next[knot] = q_head[q_dist[knot]];

      if( q_next[knot] != Q_NULL )
         q_prev[q_next[knot]] = knot;

      q_head[q_dist[knot]] = knot;

      if( q_dist[knot] < dmin )
         dmin = q_dist[knot];
      if( q_dist[knot] > (*dlmax) )
         *dlmax = q_dist[knot];
   }
   assert(q_next[knot] != Q_NM);
   assert(q_prev[knot] != Q_NM);
   assert(q_dist[knot] >= dmin);

   return(dmin);
}
#endif

static int bfs(
   const GRAPH* p,
   const int    s,
   const int    t,
   int*         q_dist,
   int*         q_numb,
   int*         q_temp,
   int*         w,
   const int*   edgestart,
   const int*   edgearr,
   const int*   headarr)
{
   int          i;
   int          j;
   int          k;
   int          l;
   int          end;
   int          visited = 0;

   assert(q_temp != NULL);
   assert(q_dist != NULL);
   assert(q_numb != NULL);
   assert(w      != NULL);
   assert(s      >= 0);
   assert(s      < p->knots);
   assert(t      >= 0);
   assert(t      < p->knots);

   /* Beginnen wir bei der Senke */
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
      for( k = edgestart[i], end = edgestart[i + 1]; k != end; k++ )
      {
         l = headarr[k];

         /* Waren wir da noch nicht ?
          */
         assert(!w[l] || (q_dist[l] >= 0));

         if( q_dist[l] < 0 )
         {
            q_dist[l] = q_dist[i] + 1;
            q_temp[visited++] = l;
            q_numb[q_dist[l]]++;

            assert(q_dist[l] < p->knots);
         }
      }
   }
   return(visited);
}


/** global relabel heuristic that sets distance of sink to zero and relabels all other nodes using backward bfs on residual
 * graph, starting from the sink.  */
static void globalrelabel(
   const GRAPH* p,
   const int    s,
   const int    t,
   int*         q_dist,
   int*         q_numb,
   int*         q_head,
   int*         q_prev,
   int*         q_next,
   int*         q_temp,
   int*         excess,
   int*         transx,
   int*         residual,
   int*         w,
   int*         dormmax,
   int*         dlmax,
   int*         dmin)
{
   int i;
   int k;
   int l;
   int j;
   int hit;
   int dist1;
   int nnodes;
   int visited;

   assert(p      != NULL);
   assert(s      >= 0);
   assert(s      <  p->knots);
   assert(t      >= 0);
   assert(t      <  p->knots);
   assert(s      != t);
   assert(w      != NULL);
   assert(q_head != NULL);
   assert(q_dist != NULL);
   assert(q_numb != NULL);
   assert(q_prev != NULL);
   assert(q_next != NULL);
   assert(q_temp != NULL);
   assert(excess != NULL);
   assert(transx != NULL);
   assert(residual != NULL);

   nnodes = p->knots;

   assert(w[s]);

   /* backwards bfs starting from sink */

   for( i = 0; i < nnodes; i++ )
   {
      q_next[i] =  Q_NULL;
      q_head[i] =  Q_NULL;
      q_numb[i] =  0;

      if( w[i] == 0 )
         w[i] = -1;
   }

   *dlmax = 0;
   (*dmin) = nnodes;

   visited = 0;

   q_temp[visited++] = t;
   q_dist[t] = 0;
   w[t] = 0;

   for( j = 0; (j < visited) && (visited < nnodes); j++ )
   {
      assert(visited < nnodes);
      assert(j < visited);

      i = q_temp[j];

      assert(i >= 0);
      assert(i < nnodes);
      assert(q_dist[i] < nnodes);
      assert(w[i] == 0);

      dist1 = q_dist[i] + 1;

      assert(dist1 < nnodes);

      /* check all neighbors */
      for( k = p->outbeg[i]; k != EAT_LAST; k = p->oeat[k] )
      {
         l = p->head[k];

         /* is not node awake, allows for positive flow along edge, and has not been visited so far? */
         if( (w[l] < 0) && (residual[flipedge(k)] > 0) )
         {
            assert(l != s);
            w[l] = 0;
            q_temp[visited++] = l;

            q_dist[l] = dist1;
            q_numb[q_dist[l]]++;

            if( excess[l] > 0 )
            {
               /* renew the distance label */
               listinsert(l, dist1, q_next, q_head, (*dmin), (*dlmax));

               assert(q_dist[l] < nnodes);
            }
         }
      }
   }

   hit = FALSE;
   for( i = 0; i < nnodes; i++ )
      if( w[i] < 0 )
      {
         hit = TRUE;
         w[i] = *dormmax + 1;
      }

   if( hit )
      *dormmax = (*dormmax) + 1;

   return;
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
   int*         residual,
   int*         w,
   const int*   edgestart,
   const int*   edgearr,
   const int*   headarr,
   int*         dormmax,
   int*         dlmax)
{

   int i;
   int j;
   int k;
   int q;
   int end;
   int actminloc;

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
   assert(residual != NULL);
   assert(p->mincut_r != NULL);
   assert(p->mincut_x != NULL);

   /* Knotenarrays initialisieren
    */
   *dormmax = 1;
   *dlmax = -1;
   actminloc = p->knots;

   for(i = 0; i < p->knots; i++)
   {
      excess[i] =  0;
      w     [i] =  0;
      q_next[i] =  Q_NULL;
      q_head[i] =  Q_NULL;
      q_numb[i] =  0;
      q_dist[i] = -1;
   }
   /* Jetzt die Kantenarrays.
    */
   for(k = 0; k < p->edges; k++)
   {
      transx[k] = 0;
      residual[k] = capa[k];
   }
   /* Jetzt noch dist und numb.
    */
   (void)bfs(p, s, t, q_dist, q_numb, q_temp, w, edgestart, edgearr, headarr);

   /* Alles was wir nicht erreichen konnten schlafen legen.
    */
   for(i = 0; i < p->knots; i++)
      if (q_dist[i] < 0)
         w[i] = *dormmax + 1;

   /* put sink into lowest dormant set */
   w[s] = 1;    /* dormmax */

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
   for( q = edgestart[s], end = edgestart[s + 1]; q != end; q++ )
   {
      k = edgearr[q];

      if (capa[k] == 0)
         continue;

      j = headarr[q];

      transx[k] = capa[k];
      residual[k] = 0;                                         /* -= transx[k] */

      residual[Edge_anti(k)] += transx[k];     /* Ueberfluessig weil w[s] == 1 */
      excess[j]            += transx[k];

      if (j != t)
      {
         listinsert(j, q_dist[j], q_next, q_head, actminloc, (*dlmax));
      }

      assert(w[j]                   == 0);
      assert(excess[j]              >  0);

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
   int*         residual,
   int*         w,
   const int*   edgestart,
   const int*   edgearr,
   const int*   headarr,
   int*         dormmax,
   int*         dlmax)
{
   int wt;
   int i;
   int j;
   int k;
   int visited;
   int actminloc;

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
   assert(residual != NULL);

   /* initialize */
   assert(w[s]);

   wt    = (w[t] == 0) ? p->knots + 1 : w[t];
   actminloc = p->knots;
   *dormmax = 1;
   *dlmax = 1;

   for( i = 0; i < p->knots; i++ )
   {
      q_numb[i] = 0;
      q_head[i] = Q_NULL;
      q_next[i] = Q_NULL;

      /* wake up nodes in higher or equal dormant layer */
      if( (w[i] == 0) || (w[i] >= wt) )
      {
         w[i] = 0;
         q_dist[i] = -1;
      }
      else if( w[i] > *dormmax )
         *dormmax = w[i];
   }
   /* Jetzt noch dist und numb.
    */
   visited = bfs(p, s, t, q_dist, q_numb, q_temp, w, edgestart, edgearr, headarr);

   /* Alles was wir nicht erreichen konnten schlafen legen.
    */
   for( i = 0; i < p->knots; i++ )
   {
      if( q_dist[i] < 0 )
         w[i] = *dormmax + 1;
      else if( (w[i] == 0) && (q_dist[i] > *dlmax) )
         *dlmax = q_dist[i];
   }

   /* Jetzt die Kantenarrays und ggf. e updaten.
    */
   for( k = 0; k < p->edges; k += 2 )
   {
      i = p->head[k];
      j = p->tail[k];

      /* both tail and head awake? */
      if (!w[i] && !w[j])
      {
         assert(w[s]);

         excess[i]    += transx[k + 1] - transx[k];
         excess[j]    += transx[k] - transx[k + 1];
         transx[k]     = 0;
         residual[k]     = capa[k];
         transx[k + 1] = 0;
         residual[k + 1] = capa[k + 1];
      }
   }

   assert(w[t]      == 0);
   assert(q_dist[t] == 0);
   assert(q_temp[0] == t);

   /* Jetzt noch die mit Excess einsortieren.
    */
   for( i = 1; i < visited; i++ )
   {
      assert(w[q_temp[i]] == 0);
      assert(q_temp[i]    != s);
      assert(q_temp[i]    != t);

      if (excess[q_temp[i]] > 0)
      {
         listinsert(q_temp[i], q_dist[q_temp[i]], q_next, q_head, actminloc, *dlmax);
      }
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
   int          nedges,
   const int    rootcutsize,
   const int*   rootcut,
   const int*   capa,
   int*         w,
   const int*   edgestart,
   const int*   edgearr,
   const int*   headarr,
   int          rerun)
{
   int    min_dist;
   int    min_capa;
   int    min_knot;
   int    min_arc;
   int    dormmax;
   int    dlmax;
   int    i;
   int    k;
   int    di1;
   int    glbadd;
   int    nnodes;
   int    startind;
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
   int j;
   int end;
   int headnode;
   int relabelupdatebnd;
   int relabeltrigger;

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
   nnodes = p->knots;

   relabelupdatebnd =  (GLOBALRELABEL_MULT * nnodes) + nedges;
   relabeltrigger = 0;

   if( !rerun )
      initialise(p, s, t, capa, dist, numb, head, prev, next, temp, e, x, r, w, edgestart,
            edgearr, headarr, &dormmax, &dlmax);
   else
      reinitialise(p, s, t, capa, dist, numb, head, prev, next, temp, e, x, r, w, edgestart,
            edgearr, headarr, &dormmax, &dlmax);

   /* main loop */
   while( dlmax >= dmin )
   {
      if( (head[dlmax] == Q_NULL) )
      {
         dlmax--;
         continue;
      }

      /* get highest label node */
      i = head[dlmax];

      listdelete(i, dlmax, next, head);

      assert(dist[i] == dlmax);

      startind = edgestart[i];
      end = edgestart[i + 1];
      glbadd = (end - startind);

      /* discharge i */
      for( ;; )
      {
         assert(w[i]    == 0);
         assert(i       != t);
         assert(i       != s);
         assert(e[i] > 0);

         min_knot = -1;
         min_dist =  nnodes;
         min_capa =  0;
         min_arc  = -1;

         di1 = dist[i] - 1;

         /* traverse incident edges */
         for( j = startind; j != end; j++ )
         {
            /* res. cap. > 0? */
            if( r[edgearr[j]] > 0 )
            {
               /* head node awake? */
               if( w[headarr[j]] == 0 )
               {
                  k = edgearr[j];
                  headnode = headarr[j];

                  /* non-admissible edge? */
                  if( di1 != dist[headnode] )
                  {
                     assert(dist[i] <= dist[headnode]);

                     if( (dist[headnode] < min_dist) || ((dist[headnode] == min_dist) && (r[k] > min_capa)) )
                     {
                        min_knot = headnode;
                        min_dist = dist[min_knot];
                        min_capa = r[k];
                        min_arc = k;
                     }
                  }
                  else /* admissible edge */
                  {
                     assert(Min(e[i], r[k]) > 0);

                     assert(headnode == p->head[k]);

                     if( e[headnode] == 0 && headnode != t )
                     {
                        listinsert(headnode, dist[headnode], next, head, dmin, (dlmax));
                     }

                     if( e[i] <= r[k] )
                     {
   #if STATIST
                        (e[i] == r[k]) ? x_pushes++ : n_pushes++;
   #endif
                        x[k] += e[i];
                        r[k] -= e[i];
                        r[Edge_anti(k)] += e[i];
                        e[headnode] += e[i];
                        e[i] = 0; /* -= e[i] */

                        assert(e[headnode] > 0);
                        assert(w[headnode] == 0);

                        assert(r[k] + r[Edge_anti(k)] == capa[k] + capa[Edge_anti(k)]);
                        assert(r[k] == capa[k] - x[k] + x[Edge_anti(k)]);

                        break;
                     }

                     /* e[i] > r[k] */
                     r[Edge_anti(k)] += r[k];
                     e[headnode] += r[k];
                     e[i] -= r[k];
                     x[k] += r[k];
                     r[k] = 0; /* -= r[k] */

                     assert(r[k] + r[Edge_anti(k)] == capa[k] + capa[Edge_anti(k)]);
                     assert(r[k] == capa[k] - x[k] + x[Edge_anti(k)]);
                     assert(e[i] > 0);

                  } /* admissible edge */
               } /* head node awake? */
            } /* res. cap. > 0? */
   #if STATIST
            s_pushes++;
   #endif
         } /* traverse incident edges */

         /* e[i] == 0? */
         if( j != end )
         {
            assert(e[i] == 0);
            break;
         }

         assert(e[i] > 0);

         /*
          * relabel
          */
         assert(numb[dist[i]] > 0);

         if( numb[dist[i]] == 1 )
         {
            w[i] = ++dormmax;

            numb[dist[i]] = 0;


            assert(dormmax <= nnodes);

            for( k = 0; k < nnodes; k++ )
            {
               /* put all nodes with distance >= dist[i] (including i) into new dormant set */
               if (!w[k] && (dist[i] < dist[k]))
               {
                  numb[dist[k]]--;
                  w[k] = dormmax;
               }
            }
   #if STATIST
            m_sleeps++;
   #endif
            break;
         }

         /* no nodes reachable from i? */
         if( min_knot < 0 )
         {
            /* put i into new dormant set */
            numb[dist[i]]--;
            w[i] = ++dormmax;

            assert(dormmax <= nnodes);

   #if STATIST
            s_sleeps++;
   #endif
            break;
         }
         else
         {
            /* set distance label of i to minimum (distance + 1) among reachable, active adjacent nodes */
            assert(min_dist <  nnodes);
            assert(min_capa >  0);
            assert(min_knot >= 0);
            assert(min_arc  >= 0);

            relabeltrigger += GLOBALRELABEL_ADD + glbadd;

            numb[dist[i]]--;

            dist[i] = min_dist + 1;

            numb[dist[i]]++;

            assert(dist[i] < nnodes);

            assert(min_capa         >  0);
            assert(min_capa         == r[min_arc]);
            assert(p->head[min_arc] == min_knot);
            assert(p->tail[min_arc] == i);
            assert(dist[i]          == dist[min_knot] + 1);
            assert(w[min_knot]      == 0);

            if (e[i] <= min_capa)
            {
               if( (e[p->head[min_arc]] == 0) && (p->head[min_arc] != t) )
               {
                  listinsert(p->head[min_arc], dist[p->head[min_arc]], next, head, dmin, (dlmax));
               }

               x[min_arc]            += e[i];
               r[min_arc]            -= e[i];
               r[Edge_anti(min_arc)] += e[i];
               e[min_knot]           += e[i];
               e[i]                   = 0;   /* -= e[i] */

               assert(r[min_arc] + r[Edge_anti(min_arc)] == capa[min_arc] + capa[Edge_anti(min_arc)]);
               assert(r[min_arc]                         >= 0);
               assert(r[min_arc]                         == capa[min_arc] - x[min_arc] + x[Edge_anti(min_arc)]);
   #if STATIST
               m_pushes++;
   #endif

               if( relabeltrigger > relabelupdatebnd )
               {
                  /* start global relabel heuristic */
                  globalrelabel(p, s, t, dist, numb, head, prev, next, temp, e, x, r, w, &dormmax, &dlmax, &dmin);
                  relabeltrigger = 0;
               }

               break;
            }
   #if STATIST
            relabels++;
   #endif
            if( relabeltrigger > relabelupdatebnd )
            {
               /* start global relabel heuristic */
               globalrelabel(p, s, t, dist, numb, head, prev, next, temp, e, x, r, w, &dormmax, &dlmax, &dmin);
               relabeltrigger = 0;
               break;
            }
         }
      } /* discharge i */
   } /* main loop */

#if 0
   FILE *fptr;

   fptr = fopen("f2", "a");

   fprintf(fptr, "t: %d flow %d \n", t, e[t]);

   fclose(fptr);
#endif

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
#endif
