/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   gminucut.cpp
 * @brief  generator for global minimum cuts in undirected graphs
 * @author Georg Skorobohatyj
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "gminucut.h"


Bool create_graph (int n, int m, GRAPH** gr)
{
   assert(gr != NULL);

   allocMemory(gr);
   if( gr == NULL )
      return FALSE;

   allocMemoryArray(&(*gr)->nodes, n);
   if( (*gr)->nodes == NULL )
   {
      freeMemory(gr);
      return FALSE;
   }

   allocMemoryArray(&(*gr)->edges, m);
   if( (*gr)->edges == NULL )
   {
      freeMemoryArray(&(*gr)->nodes);
      freeMemory(gr);
      return FALSE;
   }
   (*gr)->nuses = 1;
   (*gr)->nnodes = n;
   (*gr)->nedges = m;

   return TRUE;
}

static 
void free_graph (GRAPH** gr)
{
   assert(gr != NULL);
   assert(*gr != NULL);
   assert((*gr)->nuses == 0);
   freeMemory(&(*gr)->nodes);
   freeMemory(&(*gr)->edges);
   freeMemory(gr);
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
   if((*gr)->nuses == 0)
      free_graph(gr);
   *gr = NULL;
}

#define NO_MEM  fprintf(stderr,"Unable to allocate memory\n");
#define  EPS  1.0E-10
static long *number;

static
Bool initialize (GRAPH *gr)
{
   GRAPHNODE *nptr;
   long n, id;

   /* initial working set W contains nodes 2 to n,
      W represented as doubly linked list           */

   n = gr->nnodes;
   gr->nodes[n-1L].right_link = &(gr->nodes[1L]);
   gr->nodes[1L].left_link = &(gr->nodes[n-1L]);
   for (nptr = &(gr->nodes[n-1L]), id = n; id >= 1L; --nptr, --id)
   { if (nptr->first_edge == NULL)
      { if( n > 1 )
            fprintf (stderr,"Error in input graph - no incident edges for node %d\n",nptr->id);
         return (FALSE);
      }
      else
         nptr->scan_ptr = nptr->first_edge;
      nptr->excess = 0.0L;
      nptr->unmarked = TRUE;
      nptr->in_S = FALSE;
      nptr->in_W = TRUE;
      nptr->stack_link = NULL; 
      if (id > 2L)
         nptr->left_link = nptr - 1L;
      if (id < n)
         nptr->right_link = nptr + 1L;
   }
   return (TRUE);
}
static
long bfs0 (GRAPHNODE *t)
{
   long level, count;
   GRAPHNODE *q_rear, *q_front, *nptr;
   GRAPHEDGE *eptr; 

   t->dist = 0L;
   count = 1L;
   t->unmarked = FALSE;
   q_front = t;
   q_rear = q_front;

 bfs_next:
   level = q_rear->dist + 1L;
   eptr = q_rear->first_edge;
   do { if (eptr->adjac->unmarked && eptr->back->rcap > EPS)
      { nptr = eptr->adjac;
         nptr->unmarked = FALSE;
         nptr->dist = level;
         ++number[level];
         ++count;
         q_front->bfs_link = nptr;
         q_front = nptr;
      }
      eptr = eptr->next;
   }
   while (eptr != q_rear->first_edge);
   if (q_rear == q_front)
      goto bfs_ready;

   q_rear = q_rear->bfs_link;
   goto bfs_next;

 bfs_ready: ;
   return (count);
}

static
void bfs1 (GRAPHNODE *t)
{
   GRAPHNODE *q_front, *q_rear, *nptr;
   GRAPHEDGE *eptr;
   long level;

   t->unmarked = FALSE;
   t->dist = 0L;
   q_front = t;
   q_rear = q_front;

 bfs_next:
   level = q_rear->dist + 1L;
   eptr = q_rear->first_edge;
   do { if (eptr->adjac->in_W && eptr->adjac->unmarked && 
         eptr->back->rcap > EPS)
      { nptr = eptr->adjac;
         nptr->unmarked = FALSE;
         nptr->dist = level; 
         q_front->bfs_link = nptr;
         q_front = nptr;
      }
      eptr = eptr->next;
   }
   while (eptr != q_rear->first_edge);
   if (q_rear == q_front)
      goto bfs_ready; 

   q_rear = q_rear->bfs_link;
   goto bfs_next;

 bfs_ready:
   return;
}

Bool gmincut (GRAPH *gr, double *mincap, long *n_shore)
{
   /* Determines global minimum cut in an undirected graph,
      i.e. a cut of minimum capacity with respect to cuts
      between all pairs of nodes.

      References:
      ----------
      J. Hao/ J.B. Orlin: "A Faster Algorithm for Finding
      the Minimum Cut in a Graph", Proc. of the 3rd Annual
      ACM-SIAM Symposium on Discrete Algorithms, Orlando,
      Florida, 1992

   */
   GRAPHNODE *s_ptr, *t_ptr, *W_ptr, *w_end_ptr, 
      *aptr, *nptr, *nnptr;
   GRAPHNODE **dormant, **dor_ptr, **active, **ptr;
   GRAPHEDGE *eptr;
   long n, k, card_S, max_dor, max_dist;
   Bool found;
   double incre, cap;
   long adist, dmin, i;
  
   

   if (! initialize (gr))
      return (FALSE);

   *mincap = DBL_MAX; 
   n = gr->nnodes;
   dormant = (GRAPHNODE **) malloc (n * sizeof (GRAPHNODE *));
   if (dormant == (GRAPHNODE **) 0)
   { NO_MEM;
      return (FALSE);
   }

   active = (GRAPHNODE **) calloc (n+1L, sizeof (GRAPHNODE *));
   /* holds stacks of active nodes arranged by distances */ 
   if (active == (GRAPHNODE **) 0)
   { NO_MEM;
      return (FALSE);
   }

   number = (long *) calloc (n, sizeof (long));
   /* counts ocurrences of distances for nodes in W */
   if (number == (long *) 0)
   { NO_MEM;
      return (FALSE);
   }

   s_ptr = gr->nodes;
   s_ptr->in_S = TRUE;
   s_ptr->in_W = FALSE;
   card_S = 1;
   t_ptr = &(gr->nodes[n-1L]);

   /* breadth first search to get exact distances from first
      sink, exact distances used in test of graph connectivity */

   k = bfs0 (t_ptr);
   if (k < n)
   { /* input graph not connected */
      for (nptr = &(gr->nodes[n-1]); nptr >= gr->nodes; nptr--)
         nptr->partition = ! nptr->unmarked;
      *n_shore = k;
      *mincap = 0;
      return (TRUE);
   }

   number[0L] = 1L;
   W_ptr = &(gr->nodes[2L]);

   /* initialize set of dormant nodes */
   dormant[0L] = s_ptr;
   max_dor = 0L;
   s_ptr->left_link = s_ptr;
   s_ptr->right_link = s_ptr;

   /* initial preflow push from node s = gr->nodes[0] */

   max_dist = 0L;  /* = max_dist of active nodes */
   eptr = s_ptr->first_edge;
   do {
      nptr = eptr->adjac;
      cap = eptr->rcap;
      nptr->excess += cap;
      s_ptr->excess -= cap;
      eptr->back->rcap += cap;
      eptr->rcap = 0.0L;

      if (nptr != t_ptr && nptr->excess <= cap + EPS && cap > EPS)
      { /* push node nptr onto stack for nptr->dist,
           but only once in case of double dges      */
         assert(active[nptr->dist] == NULL || active[nptr->dist]->stack_link != nptr);
         nptr->stack_link = active[nptr->dist];
         active[nptr->dist] = nptr;
         if (nptr->dist > max_dist)
	    max_dist = nptr->dist;
      }
      
      eptr = eptr->next;
   }
   while (eptr != s_ptr->first_edge);

   /* main loop */

 next_cut:

   do { /* get maximum distance active node */
      
      aptr = active[max_dist];
     
      while (aptr != NULL)
      { /* remove node *aptr from stack */
         active[max_dist] = aptr->stack_link; 
         eptr = aptr->scan_ptr;

         /* node *aptr will not be put back onto stack again 
            in current mincut computation, either it is 
            processed until its excess becomes zero or
            else it will go into the set of dormant nodes */ 

      edge_scan:   /* for current active node */
         nptr = eptr->adjac;
     
         if (nptr->in_W && nptr->dist == aptr->dist-1L && eptr->rcap > EPS)
         {
            incre = aptr->excess;
            if (incre <= eptr->rcap)
            { /* perform a non saturating push */
               eptr->rcap -= incre;
               eptr->back->rcap += incre;
               aptr->excess = 0.0L;
               nptr->excess += incre;
               if (nptr != t_ptr && nptr->excess <= incre + EPS)
               {
                  assert(active[nptr->dist] == NULL || active[nptr->dist]->stack_link != nptr);
                  nptr->stack_link = active[nptr->dist];
                  active[nptr->dist] = nptr;
               }
               
               aptr->scan_ptr = eptr;
               goto node_ready;
            }
            else
            { /* perform a saturating push */
               incre = eptr->rcap;
               eptr->back->rcap += incre; 
               aptr->excess -= incre;
               nptr->excess += incre;
               eptr->rcap = 0.0L;
               if (nptr != t_ptr && nptr->excess <= incre + EPS)
               { 
                  assert(active[nptr->dist] == NULL || active[nptr->dist]->stack_link != nptr);
                  nptr->stack_link = active[nptr->dist];
                  active[nptr->dist] = nptr;
               }
               
               if (aptr->excess <= EPS)
               {
                  aptr->scan_ptr = eptr;
                  goto node_ready;
               }
            }
         }
         if (eptr->next == aptr->first_edge)
         { /* all admissable edges of current active node
              scanned, relabel or update set of dormant
              nodes now */
           
            adist = aptr->dist;
            if (number[adist] == 1L)
            {
               /* dist[j] != dist[i] for all j in W-{i},
                  extend dormant set by another layer */
               dor_ptr = &(dormant[++max_dor]);
               *dor_ptr = NULL;
               nptr = W_ptr;
               w_end_ptr = W_ptr->left_link;

            transfer:
               nnptr = nptr->right_link;
               if (nptr->dist >= adist)
               {
                  /* remove node nptr from set W */
                  nptr->left_link->right_link =
                     nptr->right_link;
                  nnptr->left_link = nptr->left_link; 
                  if (W_ptr == nptr)
                     W_ptr = nnptr;
                  /* W_ptr != NULL since t_ptr
                     is  contained in W       */
                  nptr->in_W = FALSE;
                  --number[nptr->dist];

                  /* clear stack for nptr->dist */
                  active[nptr->dist] = NULL;
                  nptr->scan_ptr = nptr->first_edge;

                  /* put node nptr into linked list
                     dormant[max_dor]              */
                  if (*dor_ptr == NULL)
                  { *dor_ptr = nptr;
                     nptr->right_link = nptr;
                     nptr->left_link = nptr;
                  }
                  else
                  { nptr->right_link = *dor_ptr;
                     nptr->left_link = (*dor_ptr)->left_link;
                     (*dor_ptr)->left_link = nptr;
                     nptr->left_link->right_link = nptr;
                  }
               }       
               if (nptr == w_end_ptr)
                  goto node_ready; 

               nptr = nnptr;
               goto transfer;

            } /* dist[j] != dist[i] for all j in W-{i} */
            else
            { 
               /* check if there is an edge (u, v), u=*aptr,
                  such that v in W and rcap(u,v) > 0    */
               eptr = aptr->first_edge;
               found = FALSE;
               do { if (eptr->adjac->in_W && eptr->rcap > EPS)
                  { found = TRUE;
                     dmin = eptr->adjac->dist;
                     break;
                  }
                  else
                     eptr = eptr->next;
               }
               while (eptr != aptr->first_edge);
               if (found)
               { aptr->scan_ptr = eptr;

                  /* get new distance label for *aptr */
                  while (eptr->next != aptr->first_edge)
                  { eptr = eptr->next;
                     if (eptr->adjac->in_W &&
                        eptr->adjac->dist < dmin && 
                        eptr->rcap > EPS) 
                        dmin = eptr->adjac->dist;
                  } 
                  --number[adist];
                  aptr->dist = dmin + 1L;
                  ++number[dmin+1];
                  max_dist = dmin + 1L;
                  eptr = aptr->scan_ptr;
                  goto edge_scan;
               }
               else
               { /* extend dormant set by another
                    layer containing node *aptr only,
                    remove aptr from W nodes first    */ 

                  aptr->in_W = FALSE;
                  --number[adist];
                  aptr->scan_ptr = aptr->first_edge;
                  aptr->left_link->right_link =
                     aptr->right_link;
                  aptr->right_link->left_link =
                     aptr->left_link;
                  if (W_ptr == aptr)
                     W_ptr = aptr->right_link;
                  dormant[++max_dor] = aptr;
                  aptr->right_link = aptr;
                  aptr->left_link = aptr;
                  goto node_ready;
               } 
            }
         } 
         else
         {
           
            eptr = eptr->next;
            goto edge_scan;
         }

      node_ready:
         aptr = active[max_dist];
      } /* aptr != NULL */
      --max_dist;
   }
   while (max_dist > 0L);

 check_min:

   if (*mincap > t_ptr->excess)
   { *mincap = t_ptr->excess;
      *n_shore = 0L;
      for (nptr = &(gr->nodes[n-1L]); nptr >= gr->nodes; --nptr)
         if (nptr->in_W)
            nptr->partition = FALSE;
         else
	 { nptr->partition = TRUE;
            ++(*n_shore);
         }
   }

   /* preparations for next cut step */

   /* delete t_ptr from W */
   t_ptr->in_W = FALSE;
   --number[t_ptr->dist];
   if (t_ptr->right_link == t_ptr)
      W_ptr = NULL;
   else
   { t_ptr->right_link->left_link = t_ptr->left_link;
      t_ptr->left_link->right_link = t_ptr->right_link;
      if (W_ptr == t_ptr)
         W_ptr = t_ptr->right_link;
   }

   /* put t_ptr into source set S and dormant[0] set */
   t_ptr->in_S = TRUE;
   ++card_S;
   t_ptr->right_link = dormant[0L]->right_link;
   t_ptr->left_link = dormant[0L];
   dormant[0L]->right_link->left_link = t_ptr;
   dormant[0L]->right_link = t_ptr;

   if (card_S == n)
      goto mincut_ready;
          
   /* saturate all arcs from *t_ptr to nodes not in S */
   eptr = t_ptr->first_edge;
   do { nptr = eptr->adjac;
      if (! nptr->in_S && eptr->rcap > EPS)
      { t_ptr->excess -= eptr->rcap;
         nptr->excess += eptr->rcap;
         eptr->back->rcap += eptr->rcap;
         eptr->rcap = 0.0L;
         nptr->scan_ptr = nptr->first_edge;
      }
      eptr = eptr->next;
   }
   while (eptr != t_ptr->first_edge);

   if (W_ptr == NULL)
   { /* set of W nodes empty, dormant[max_dor]
        taken as next set of W nodes
     */
      W_ptr = dormant[max_dor--];
      nptr = W_ptr;
      do { nptr->in_W = TRUE;
         nptr = nptr->right_link;
      }
      while (nptr != W_ptr);
   }

   /* get node from W with minimum distance as new sink */
   dmin = W_ptr->dist;
   t_ptr = W_ptr;
   nptr = t_ptr;
   while (nptr->right_link != W_ptr)
   { nptr = nptr->right_link;
      if (nptr->dist < dmin)
      { dmin = nptr->dist;
         t_ptr = nptr;
      }
   }

   /* breadth first search to get exact distances
      for nodes of W with respect to new sink,
      nodes of W with positive excess will be pushed
      onto stack of active nodes, not all nodes of W
      are reachable in the residual graph by breadth
      first search, however, all such nodes are put
      into another dormant set                       */

   if (W_ptr->right_link == W_ptr)
      /* only one node left in W = new sink */
      goto check_min; 

   for (i = n-1L, ptr = &(active[n-1L]); i >= 0L; i--, ptr--)
   { number[i] = 0L;  
      *ptr = NULL;
   }
   nptr = t_ptr;
   while (nptr->right_link != t_ptr)
   { nptr = nptr->right_link;
      nptr->unmarked = TRUE;
      nptr->scan_ptr = nptr->first_edge;
   }

   max_dist = 0L;
   number[0L] = 1L;

   bfs1 (t_ptr);

   /*   check next set W for nodes to be transferred 
	to another dormant set and for active nodes 
	to be pushed onto stack                      */

   dor_ptr = &(dormant[max_dor+1L]);
   *dor_ptr = NULL;
   nptr = W_ptr;
   w_end_ptr = W_ptr->left_link;

 check_W:
   nnptr = nptr->right_link;
   if (nptr->unmarked)
   {
      /* remove node *nptr from set W */
      nptr->in_W = FALSE;
      nptr->right_link->left_link = nptr->left_link;
      nptr->left_link->right_link = nptr->right_link;
      if (W_ptr == nptr)
         W_ptr = nnptr; 

      /* put node *nptr into new set dormant[dor_max+1] */
      if (*dor_ptr == NULL)
      { *dor_ptr = nptr;
	 nptr->right_link = nptr;
	 nptr->left_link = nptr;
      }
      else
      { nptr->right_link = (*dor_ptr)->right_link;
	 nptr->left_link = (*dor_ptr);
	 nptr->right_link->left_link = nptr;
	 (*dor_ptr)->right_link = nptr;
      }
   }
   else
      if (nptr != t_ptr)
      { 
         ++number[nptr->dist]; 
         if (nptr->excess > EPS)
         { 
            assert(active[nptr->dist] == NULL || active[nptr->dist]->stack_link != nptr);
            nptr->stack_link = active[nptr->dist];
            active[nptr->dist] = nptr;
            if (nptr->dist > max_dist)
               max_dist = nptr->dist;
         }
      }
   if (nptr == w_end_ptr)
      goto end_check;

   nptr = nnptr;
   goto check_W;

 end_check:
   if (*dor_ptr != NULL)
      ++max_dor;

   goto next_cut;

 mincut_ready:
   return (TRUE);
}
