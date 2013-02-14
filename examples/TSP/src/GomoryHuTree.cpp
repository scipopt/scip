/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   GomoryHuTree.cpp
 * @brief  generator for global cuts in undirected graphs
 * @author Georg Skorobohatyj
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <stdio.h>
#include <assert.h>

#include "objscip/objscip.h"
#include "GomoryHuTree.h"

#define  EPS  1.0E-10

static GRAPHNODE **active;
static long *number;
static long max_dist, bound;
static SCIP_Bool co_check;

SCIP_Bool create_graph (int n, int m, GRAPH** gr)
{
   assert( gr != NULL );

   BMSallocMemory(gr);
   if( *gr == NULL )
      return FALSE;

   BMSallocMemoryArray( &(*gr)->nodes, n );
   if( (*gr)->nodes == NULL )
   {
      BMSfreeMemory(gr);
      return FALSE;
   }

   BMSallocMemoryArray( &(*gr)->edges, m );
   if( (*gr)->edges == NULL )
   {
      BMSfreeMemoryArray(&(*gr)->nodes);
      BMSfreeMemory(gr);
      return FALSE;
   }
   (*gr)->nuses = 1;
   (*gr)->nnodes = n;
   (*gr)->nedges = m/2;
   (*gr)->nedgesnonzero = m/2;
   return TRUE;
}

static 
void free_graph (GRAPH** gr)
{
   assert(gr != NULL);
   assert(*gr != NULL);
   assert((*gr)->nuses == 0);
   BMSfreeMemory(&(*gr)->nodes);
   BMSfreeMemory(&(*gr)->edges);
   BMSfreeMemory(gr);
}

void capture_graph (GRAPH* gr)
{
   assert(gr != NULL);
   gr->nuses++;
}

void release_graph (GRAPH** gr)
{
   assert(gr != NULL);
   assert(*gr != NULL);
   (*gr)->nuses--;
   if( (*gr)->nuses == 0 )
      free_graph(gr);
   *gr = NULL;
}

SCIP_Bool init_maxflow (long n)
{
   active = (GRAPHNODE **) malloc ((n+1L) * sizeof (GRAPHNODE *));
   /* holds stacks of active nodes arranged by distances */ 
   if ( active == (GRAPHNODE **) 0 )
   { 
      printf ("Unable to allocate memory\n");
      return FALSE;
   }
   number = (long *) malloc ((n+1L) * sizeof (long));
   /* counts occurences of node distances in set 
      of alive nodes, i.e. nodes not contained in
      the set of nodes disconnected from the sink */ 
   if ( number == (long *) 0 )
   { 
      printf ("Unable to allocate memory\n");
      return FALSE;
   }
   co_check = TRUE;
   return TRUE;

} /* init_maxflow */


void fini_maxflow ()
{
   free (active);
   free (number);

} /* fini_maxflow */


void global_relabel (GRAPH *gr, GRAPHNODE *tptr)
{ 
   /* breadth first search to get exact distance labels
      from sink with reordering of stack of active nodes */

   GRAPHNODE *front, *rear, *nptr, **ptr;
   GRAPHEDGE *eptr;
   long n, level, count, i;

   n = gr->nnodes;
   for ( nptr = &(gr->nodes[n-1L]); nptr >= gr->nodes; nptr-- )
   { 
      nptr->unmarked = TRUE;
      nptr->stack_link = NULL;
      nptr->scan_ptr = nptr->first_edge;
      while( nptr->scan_ptr != NULL && nptr->scan_ptr->cap <= EPS ) 
         nptr->scan_ptr = nptr->scan_ptr->next;
   }
   tptr->unmarked = FALSE;

   /* initialize stack of active nodes */
   for ( ptr = &(active[n]); ptr >= active; ptr-- )
      *ptr = NULL;
   for ( i = 0L; i <= n; i++ )
      number[i] = 0L;
   max_dist = 0L;
   count = 1L;     /* number of alive nodes */
   front = tptr;
   rear = front;

 bfs_next:
   level = rear->dist + 1L;
   eptr = rear->first_edge;
   while ( eptr != NULL )
   {
      nptr = eptr->adjac;
      if ( nptr->alive && nptr->unmarked && eptr->back->rcap > EPS ) 
      {
         assert(eptr->cap > EPS);
         assert(eptr->back->cap > EPS);
         nptr->unmarked = FALSE;
         nptr->dist = level;
         ++count;
         ++number[level];
         if ( nptr->excess > EPS )
         { 
            nptr->stack_link = active[level];
            active[level] = nptr;
            max_dist = level;
         }
         front->bfs_link = nptr;
         front = nptr;
      }
      eptr = eptr->next;
   }
   if ( front == rear )
      goto bfs_ready;
       
   rear = rear->bfs_link;
   goto bfs_next;

 bfs_ready: 

   if ( count < bound )
   {
      /* identify nodes that are marked alive but have
         not been reached by BFS and mark them as dead  */
      for ( nptr = &(gr->nodes[n-1L]); nptr >= gr->nodes; nptr-- )
      {
         if ( nptr->unmarked && nptr->alive )
         {
            nptr->dist = n;
            nptr->alive = FALSE;
         }
      }
      bound = count;
   }

} /* global_relabel */


double maxflow (GRAPH *gr, GRAPHNODE *s_ptr, GRAPHNODE *t_ptr)
{
   /* Determines maximum flow and minimum cut between nodes
      s (= *s_ptr) and t (= *t_ptr) in an undirected graph  

      References:
      ----------
      A. Goldberg/ E. Tarjan: "A New Approach to the
      Maximum Flow Problem", in Proc. 18th ACM Symp. 
      on Theory of Computing, 1986.
   */
   GRAPHNODE *aptr, *nptr, *q_front, *q_rear;
   GRAPHEDGE *eptr;
   long n, m, m0, level, i, n_discharge;
   double incre;
   long dmin;
   double cap;
   
   /* node ids range from 1 to n, node array indices  
      range from 0 to n-1                             */

   n = gr->nnodes;
   for ( nptr = &(gr->nodes[n-1L]); nptr >= gr->nodes; nptr-- )
   {  
      nptr->scan_ptr = nptr->first_edge;
      while( nptr->scan_ptr != NULL && nptr->scan_ptr->cap <= EPS ) 
         nptr->scan_ptr = nptr->scan_ptr->next;
      
      if ( nptr->scan_ptr == NULL )
      { 
         fprintf (stderr, "isolated node in input graph\n");
         return (FALSE);
      }

      nptr->excess = 0.0L;
      nptr->stack_link = NULL;
      nptr->alive = TRUE;
      nptr->unmarked = TRUE;
   }

   m = gr->nedgesnonzero;
   m0 = gr->nedges;
   for ( eptr = &(gr->edges[m-1L]); eptr >= gr->edges; eptr-- )
      eptr->rcap = eptr->cap;
   for ( eptr = &(gr->edges[m0+m-1L]); eptr >= &(gr->edges[m0]); eptr-- )
      eptr->rcap = eptr->cap;
      
   for ( i = n; i >= 0L; i-- )
   { 
      number[i] = 0L;
      active[i] = NULL;
   }
   t_ptr->dist = 0L;

   /* breadth first search to get exact distances 
      from sink and for test of graph connectivity */

   t_ptr->unmarked = FALSE; 
   q_front = t_ptr;
   q_rear = q_front;
 bfs_next:
   level = q_rear->dist + 1L;
   eptr = q_rear->first_edge;
   while ( eptr != NULL )
   {
      assert(eptr->back->cap == eptr->back->rcap);
      if ( eptr->adjac->unmarked && eptr->back->rcap > EPS )
      {
         nptr = eptr->adjac;
         nptr->unmarked = FALSE;
         nptr->dist = level;
         ++number[level];
         q_front->bfs_link = nptr;
         q_front = nptr;
      }
      eptr = eptr->next;
   }
   if ( q_rear == q_front )
      goto bfs_ready;

   q_rear = q_rear->bfs_link;
   goto bfs_next;

 bfs_ready:
   if ( co_check )
   { 
      co_check = FALSE;
      for ( nptr = &(gr->nodes[n-1]); nptr >= gr->nodes; --nptr )
      {
         if ( nptr->unmarked )
            return (-1.0L);
      }
   }

   s_ptr->dist = n; /* number[0] and number[n] not required */
   t_ptr->dist = 0L;
   t_ptr->excess = 1.0L;  /* to be subtracted again */


   /* initial preflow push from source node */

   max_dist = 0L;  /* = max_dist of active nodes */
   eptr = s_ptr->first_edge;
   while ( eptr != NULL )
   {
      if( eptr->cap > EPS ) 
      {
         nptr = eptr->adjac;
         cap = eptr->rcap;
         nptr->excess += cap;
         s_ptr->excess -= cap;
         eptr->back->rcap += cap;
         eptr->rcap = 0.0L;

         if ( nptr != t_ptr && nptr->excess <= cap + EPS ) 
         { 
            /* push node nptr onto stack for nptr->dist,
               but only once in case of double edges     */
            nptr->stack_link = active[nptr->dist];
            active[nptr->dist] = nptr;
            if ( nptr->dist > max_dist )
               max_dist = nptr->dist;
         }
      }
      eptr = eptr->next;
   }

   s_ptr->alive = FALSE;
   bound = n; 
   n_discharge = 0L;

   /* main loop */

   do
   { 
      /* get maximum distance active node */
      aptr = active[max_dist];
      while ( aptr != NULL )
      {
         active[max_dist] = aptr->stack_link;
         eptr = aptr->scan_ptr;
         assert(eptr != NULL);
         assert(eptr->back != NULL);
         assert(eptr->cap > EPS);

      edge_scan:  /* for current active node  */
         assert(eptr != NULL);
         assert(eptr->back != NULL);
         assert(eptr->cap > EPS);

         nptr = eptr->adjac;
         assert(nptr != NULL);
         if ( nptr->dist == aptr->dist - 1L && eptr->rcap > EPS ) 
         { 
            incre = aptr->excess;
            if ( incre <= eptr->rcap )
            { 
               /* perform a non saturating push */
               eptr->rcap -= incre;
               eptr->back->rcap += incre;
               aptr->excess = 0.0L;
               nptr->excess += incre;
               if ( nptr->excess <= incre + EPS )
               {
                  /* push nptr onto active stack */
                  nptr->stack_link = active[nptr->dist];
                  active[nptr->dist] = nptr;
               }
               assert(eptr->cap > EPS);
               aptr->scan_ptr = eptr;
               goto node_ready;
            }
            else
            {
               /* perform a saturating push */
               incre = eptr->rcap;
               eptr->back->rcap += incre; 
               aptr->excess -= incre;
               nptr->excess += incre;
               eptr->rcap = 0.0L;
               if ( nptr->excess <= incre + EPS )
               {
                  /* push nptr onto active stack */
                  nptr->stack_link = active[nptr->dist];
                  active[nptr->dist] = nptr;
               }
               if ( aptr->excess <= EPS )
               {
                  assert(eptr->cap > EPS);
                  aptr->scan_ptr = eptr;
                  goto node_ready;
               }
            }
         }

         /* go to the next non-zero edge */
         do 
         {
            eptr = eptr->next;
         }
         while( eptr != NULL && eptr->cap <= EPS );

         if ( eptr == NULL )
         {
            /* all incident arcs scanned, but node still
               has positive excess, check if for all nptr       
               nptr->dist != aptr->dist                  */

            if ( number[aptr->dist] == 1L )
            {
               /* put all nodes v with dist[v] >= dist[a] 
                  into the set of "dead" nodes since they
                  are disconnected from the sink          */
                     
               for ( nptr = &(gr->nodes[n-1L]); nptr >= gr->nodes; nptr-- )
               {
                  if ( nptr->alive && nptr->dist > aptr->dist )
                  { 
                     --number[nptr->dist];
                     active[nptr->dist] = NULL; 
                     nptr->alive = FALSE; 
                     nptr->dist = n;
                     --bound;
                  }
               }
               --number[aptr->dist];
               active[aptr->dist] = NULL;
               aptr->alive = FALSE;
               aptr->dist = n;
               --bound;
               goto node_ready;
            }
            else
            {
               /* determine new label value */
               dmin = n;
               aptr->scan_ptr = NULL;
               eptr = aptr->first_edge;
               while ( eptr != NULL )
               {
                  assert(eptr->rcap <= EPS || eptr->cap > EPS);
                  if ( eptr->adjac->dist < dmin && eptr->rcap > EPS )
                  { 
                     assert(eptr->cap > EPS);
                     dmin = eptr->adjac->dist;
                     if ( aptr->scan_ptr == NULL )
                        aptr->scan_ptr = eptr;
                  }
                  eptr = eptr->next;
               }
               
               if ( ++dmin < bound )
               { 
                  /* ordinary relabel operation */
                  --number[aptr->dist];
                  aptr->dist = dmin;
                  ++number[dmin];
                  max_dist = dmin;
                  eptr = aptr->scan_ptr;
                  assert(eptr != NULL);
                  assert(eptr->cap > EPS);
                  goto edge_scan;
               }
               else
               { 
                  aptr->alive = FALSE;
                  --number[aptr->dist];
                  aptr->dist = n;
                  --bound;
                  goto node_ready;
               }
            }
         }
         else
            goto edge_scan;
         

      node_ready: 
         ++n_discharge;
         if ( n_discharge == n )
         { 
            n_discharge = 0L;
            global_relabel (gr, t_ptr);
         }
         aptr = active[max_dist];
      } /* aptr != NULL */ 
      --max_dist;
   } 
   while ( max_dist > 0L );  

   return (t_ptr->excess - 1.0L);
}

SCIP_Bool nodeOnRootPath(GRAPH* gr, int i, int j)
{
   while( i != j && j != 0 )
   {
      j = gr->nodes[j].parent->id;
   }
   if( i == j ) 
      return TRUE;
   else
      return FALSE;
}

// constructs a list of cuts for a TSP relaxation polytope from a Gomory Hu Tree 
void constructCutList(GRAPH *gr, SCIP_Bool** cuts, int* ncuts, double minviol)
{
   int k = 0;
   for( int i = 1; i < gr->nnodes; i++ )
   {
      if( gr->nodes[i].mincap < 2.0 - minviol )
      {
         for( int j = 1 ; j < gr->nnodes; j++ )
            cuts[k][j] = nodeOnRootPath(gr, i, j);
         k++;
      }
   }
   *ncuts = k;
}

// if the non-zero-edges of a TSP relaxation induce a non-connected graph, 
// an according cut is generated, using information from BFS in method maxflow
void constructSingleCut(GRAPH *gr, SCIP_Bool** cuts)
{
   for( int i = 1; i < gr->nnodes; i++ )
      cuts[0][i]=gr->nodes[i].unmarked;
}

SCIP_Bool ghc_tree (GRAPH *gr, SCIP_Bool** cuts, int* ncuts, double minviol)
{
   /* Determines Gomory/Hu cut tree for input graph with
      capacitated edges, the tree structures is represented
      by parent pointers which are part of the node structure,
      the capacity of a tree edge is stored at the child node,
      the root of the cut tree is the first node in the list
      of graph nodes (&gr->nodes[0]). The implementation is
      described in [1].
     
      References:
      ----------
      1) D. Gusfield: "Very Simple Algorithms and Programs for
      All Pairs Network Flow Analysis", Computer Science Di-
      vision, University of California, Davis, 1987.
     
      2) R.E. Gomory and T.C. Hu: "Multi-Terminal Network Flows",
      SIAM J. Applied Math. 9 (1961), 551-570.

   */

   GRAPHNODE *nptr, *nptr1, *nptrn, *sptr, *tptr;
   long n;
   double maxfl;

   n = gr->nnodes;
   
   if ( ! init_maxflow (n) )
      return FALSE;

   nptr1 = gr->nodes;
   nptrn = &(gr->nodes[n-1L]);
   for ( nptr = nptrn; nptr >= nptr1; nptr-- )
      nptr->parent = nptr1;

   for ( sptr = &(gr->nodes[1L]); sptr <= nptrn; sptr++ )
   { 
      tptr = sptr->parent;
      maxfl = maxflow(gr, sptr, tptr);

      //maxfl < 0 <=> graph is not connected => generate cut
      if ( maxfl < 0L )
      {
         constructSingleCut(gr,cuts);
         *ncuts = 1;
	 fini_maxflow ();
         return TRUE;
      }
           
      sptr->mincap = maxfl;
      for ( nptr = &(gr->nodes[1L]); nptr <= nptrn; nptr++ )
         if ( nptr != sptr && ! nptr->alive && nptr->parent == tptr )
            nptr->parent = sptr;
      if ( ! tptr->parent->alive )
      { 
         sptr->parent = tptr->parent;
         tptr->parent = sptr;
         sptr->mincap = tptr->mincap;
         tptr->mincap = maxfl;
      }
   }
   fini_maxflow ();
   constructCutList(gr,cuts,ncuts,minviol);
   
   return TRUE;
} /* ghc_tree */
