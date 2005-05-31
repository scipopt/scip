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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                        This file is part of the program                   */
/*                    TCLIQUE --- Algorithm for Maximum Cliques              */
/*                                                                           */
/*    Copyright (C) 1996-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: tclique_branch.c,v 1.8 2005/05/31 17:20:23 bzfpfend Exp $"

/**@file   tclique_branch.c
 * @brief  branch and bound part of algorithm for maximum cliques
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "scip/def.h"
#include "scip/memory.h"
#include "scip/message.h"
#include "scip/tclique_graph.h"
#include "scip/tclique_branch.h"
#include "scip/tclique_coloring.h"



#define CHUNK_SIZE (64)
#define CLIQUETABLE_INITSIZE (1024)



/***************************
 * clique hash table methods
 ***************************/

typedef struct clique CLIQUE;           /**< single element of the clique table */

/** single element of the clique table */
struct clique
{
   int*             nodes;              /**< node number of the clique elements */
   int              nnodes;             /**< length of the clique */
};

typedef struct cliquetable CLIQUETABLE; /**< table for storing cliques */

/** table for storing cliques */
struct cliquetable
{
   CLIQUE**         cliques;            /**< elements of the table */
   int              cliquessize;        /**< size of the table */
   int              ncliques;           /**< number of cliques stored in the table */
};


/** creates a clique */
static
void createClique(
   CLIQUE**         clique,             /**< pointer to the clique */
   int*             nodes,              /**< nodes of the clique */
   int              nnodes              /**< number of nodes in the clique */
   )
{
   int i;

   assert(clique != NULL);

   ALLOC_ABORT( allocMemory(clique) );
   ALLOC_ABORT( allocMemoryArray(&(*clique)->nodes, nnodes) );

   /* sort the nodes into the clique's node array */
   for( i = 0; i < nnodes; ++i )
   {
      int node;
      int j;

      node = nodes[i];
      for( j = i; j > 0 && node < (*clique)->nodes[j-1]; --j )
         (*clique)->nodes[j] = (*clique)->nodes[j-1];
      (*clique)->nodes[j] = node;
   }
   (*clique)->nnodes = nnodes;
}

/** frees a clique */
static
void freeClique(
   CLIQUE**         clique              /**< pointer to the clique */
   )
{
   assert(clique != NULL);
   assert(*clique != NULL);

   freeMemoryArray(&(*clique)->nodes);
   freeMemory(clique);
}

/** checks, whether clique1 is a subset of clique2 and returns the following value:
 *   == 0 if clique1 == clique2, or clique1 is contained in clique2,
 *    < 0 if clique1 < clique2, and clique1 is not contained in clique2,
 *    > 0 if clique1 > clique2, and clique1 is not contained in clique2
 */
static
int compSubcliques(
   CLIQUE*          clique1,            /**< first clique to compare */
   CLIQUE*          clique2             /**< second clique to compare */
   )
{
   int pos1;
   int pos2;
   Bool clique2smaller;

   assert(clique1 != NULL);
   assert(clique2 != NULL);

   pos1 = 0;
   pos2 = 0;
   clique2smaller = FALSE;
   while( pos1 < clique1->nnodes && pos2 < clique2->nnodes )
   {
      if( clique1->nodes[pos1] < clique2->nodes[pos2] )
      {
         /* clique1 has an element not contained in clique2: clique1 is lex-smaller, if clique2 was not
          * detected earlier to be lex-smaller
          */
         return (clique2smaller ? +1 : -1);
      }
      else if( clique1->nodes[pos1] > clique2->nodes[pos2] )
      {
         /* clique2 has an element not contained in clique1: clique2 is lex-smaller, but probably clique1 is
          * contained in clique2
          */
         pos2++;
         clique2smaller = TRUE;
      }
      else
      {
         pos1++;
         pos2++;
      }
   }

   /* if clique1 has additional elements, it is not contained in clique2 */
   if( pos1 < clique1->nnodes )
      return (clique2smaller ? +1 : -1);

   /* clique1 is contained in clique2 */
   return 0;
}

#ifndef NDEBUG
/** performs an integrity check of the clique table */
static
void checkCliquetable(
   CLIQUETABLE*     cliquetable         /**< clique table */
   )
{
   int i;

   assert(cliquetable != NULL);

   for( i = 0; i < cliquetable->ncliques-1; ++i )
      assert(compSubcliques(cliquetable->cliques[i], cliquetable->cliques[i+1]) < 0);
}
#else
#define checkCliquetable(cliquetable) /**/
#endif

/** creates a table for storing cliques */
static
void createCliquetable(
   CLIQUETABLE**    cliquetable,        /**< pointer to store the clique table */
   int              tablesize           /**< initial size of the clique table */
   )
{
   assert(cliquetable != NULL);
   assert(tablesize > 0);

   ALLOC_ABORT( allocMemory(cliquetable) );
   ALLOC_ABORT( allocMemoryArray(&(*cliquetable)->cliques, tablesize) );
   (*cliquetable)->cliquessize = tablesize;
   (*cliquetable)->ncliques = 0;
}

/** clears the clique table and frees all inserted cliques */
static
void clearCliquetable(
   CLIQUETABLE*     cliquetable         /**< clique table */
   )
{
   int i;

   assert(cliquetable != NULL);

   /* free the cliques in the table */
   for( i = 0; i < cliquetable->ncliques; ++i )
      freeClique(&cliquetable->cliques[i]);

   cliquetable->ncliques = 0;
}

/** frees the table for storing cliques and all inserted cliques */
static
void freeCliquetable(
   CLIQUETABLE**    cliquetable         /**< pointer to the clique table */
   )
{
   assert(cliquetable != NULL);
   assert(*cliquetable != NULL);

   /* free the cliques in the table */
   clearCliquetable(*cliquetable);

   /* free the table data structure */
   freeMemoryArray(&(*cliquetable)->cliques);
   freeMemory(cliquetable);
}

/** ensures, that the clique table is able to store at least the given number of cliques */
static
void ensureCliquetableSize(
   CLIQUETABLE*     cliquetable,        /**< clique table */
   int              num                 /**< minimal number of cliques to store */
   )
{
   assert(cliquetable != NULL);

   if( num > cliquetable->cliquessize )
   {
      int newsize;

      newsize = 2*cliquetable->cliquessize;
      if( num > newsize )
         newsize = num;

      ALLOC_ABORT( reallocMemoryArray(&cliquetable->cliques, newsize) );
      cliquetable->cliquessize = newsize;
   }
   assert(cliquetable->cliquessize >= num);
}

#ifdef DEBUG
/** displayes clique table */
static
void printCliquetable(
   CLIQUETABLE*     cliquetable         /**< clique table */
   )
{
   int i;

   assert(cliquetable != NULL);

   printf("cliquetable (%d cliques):\n", cliquetable->ncliques);
   for( i = 0; i < cliquetable->ncliques; ++i )
   {
      int j;

      printf("%d:", i);
      for( j = 0; j < cliquetable->cliques[i]->nnodes; ++j )
         printf(" %d", cliquetable->cliques[i]->nodes[j]);
      printf("\n");
   }
}
#endif

/** searches the given clique in the clique table and returns whether it (or a stronger clique) exists */
static
Bool inCliquetable(
   CLIQUETABLE*     cliquetable,        /**< clique table */
   CLIQUE*          clique,             /**< clique to search in the table */
   int*             insertpos           /**< position where the clique should be inserted in the table */
   )
{
   int left;
   int right;
   int middle;
   int cmp;

   assert(cliquetable != NULL);
   assert(cliquetable->cliquessize > 0);
   assert(clique != NULL);
   assert(insertpos != NULL);

   /* perform a binary search on the clique table */
   left = 0;
   right = cliquetable->ncliques-1;
   while( left <= right )
   {
      middle = (left+right)/2;
      cmp = compSubcliques(clique, cliquetable->cliques[middle]);
      if( cmp > 0 )
         left = middle+1;
      else if( cmp < 0 )
         right = middle-1;
      else
      {
         *insertpos = middle;
         return TRUE;
      }
   }

   /* we found the correct insertion position for the clique, but it might be contained in a lex-smaller clique */
   *insertpos = left;
   for( middle = left-1; middle >= 0; --middle )
   {
      cmp = compSubcliques(clique, cliquetable->cliques[middle]);
      assert(cmp >= 0);
      if( cmp == 0 )
         return TRUE;
   }

   return FALSE;
}

/** inserts clique into clique table */
static
void insertClique(
   CLIQUETABLE*     cliquetable,        /**< clique table */
   CLIQUE*          clique,             /**< clique to search in the table */
   int              insertpos           /**< position to insert clique into table (returned by inCliquetable()) */
   )
{
   int i;

   assert(cliquetable != NULL);
   assert(clique != NULL);
   assert(0 <= insertpos && insertpos <= cliquetable->ncliques);

   /* allocate memory */
   ensureCliquetableSize(cliquetable, cliquetable->ncliques+1);
   
   /* insert clique into table */
   for( i = cliquetable->ncliques; i > insertpos; --i )
      cliquetable->cliques[i] = cliquetable->cliques[i-1];
   cliquetable->cliques[insertpos] = clique;
   cliquetable->ncliques++;

   /* check, whether the clique table is still sorted */
   checkCliquetable(cliquetable);

   debug(printCliquetable(cliquetable));
}




/****************************
 * clique calculation methods
 ****************************/

/** extends given clique by additional zero-weight nodes of the given node set */
static
void extendCliqueZeroWeight(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int*             Vzero,              /**< zero weighted nodes */
   int              nVzero,             /**< number of non-zero weighted nodes */
   int*             curcliquenodes,     /**< nodes of the clique */
   int*             ncurcliquenodes     /**< pointer to store number of nodes in the clique */
   )
{
   WEIGHT* weights;
   int* adjcount;
   int nnodes;
   int i;

   assert(Vzero != NULL);
   assert(curcliquenodes != NULL);
   assert(ncurcliquenodes != NULL);

   debugMessage("extending temporary clique (size %d) with zero-weighted nodes (nVzero=%d)\n", *ncurcliquenodes, nVzero);

   nnodes = tcliqueGetNNodes(tcliquedata);
   weights = tcliqueGetWeights(tcliquedata);
   assert(weights != NULL);

   /* for each node, count the number of adjacent clique nodes */
   ALLOC_ABORT( allocMemoryArray(&adjcount, nnodes) );
   clearMemoryArray(adjcount, nnodes);
   for( i = 0; i < *ncurcliquenodes; ++i )
   {
      int* currentadjedge;
      int* lastadjedge;
      int node;

      node = curcliquenodes[i];
      assert(weights[node] > 0);

      /* increase the degrees of the adjecent nodes */
      currentadjedge = tcliqueGetFirstAdjedge(tcliquedata, node);
      lastadjedge = tcliqueGetLastAdjedge(tcliquedata, node);
      for( ; currentadjedge <= lastadjedge; ++currentadjedge )
         adjcount[*currentadjedge]++;
   }

   /* put a zero-weighted node into the clique, if it is adjacent to all clique nodes */
   for( i = 0; i < nVzero; ++i )
   {
      int node;

      node = Vzero[i];
      if( adjcount[node] == *ncurcliquenodes && weights[node] == 0 )
      {
         int* currentadjedge;
         int* lastadjedge;

         /* put node into the clique */
         curcliquenodes[*ncurcliquenodes] = node;
         (*ncurcliquenodes)++;

         /* update the degrees of the adjecent nodes */
         currentadjedge = tcliqueGetFirstAdjedge(tcliquedata, node);
         lastadjedge = tcliqueGetLastAdjedge(tcliquedata, node);
         for( ; currentadjedge <= lastadjedge; ++currentadjedge )
            adjcount[*currentadjedge]++;

         debugMessage(" -> extended clique by zero-weighted node %d (new size: %d)\n", node, *ncurcliquenodes);
      }
   }

   freeMemoryArray(&adjcount);
}

/** calls user callback after a new solution was found, that is better than the current incumbent;
 *  the callback decides, whether this solution should be accepted as new incumbent, and whether the solution process
 *  should be stopped
 */
static
void newSolution(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   CLIQUETABLE*     cliquetable,        /**< clique table */
   int*             Vzero,              /**< zero weighted nodes */
   int              nVzero,             /**< number of non-zero weighted nodes */
   int*             curcliquenodes,     /**< nodes of the new clique */
   int              ncurcliquenodes,    /**< number of nodes in the new clique */
   WEIGHT           curcliqueweight,    /**< weight of the new clique */
   int*             maxcliquenodes,     /**< pointer to store nodes of the maximum weight clique */
   int*             nmaxcliquenodes,    /**< pointer to store number of nodes in the maximum weight clique */
   WEIGHT*          maxcliqueweight,    /**< pointer to store weight of the maximum weight clique */
   TCLIQUE_USRCALLBACK ((*usrcallback)),/**< user function to call on every new solution */
   void*            usrdata,            /**< user data to pass to user callback function */
   Bool*            stopsolving         /**< pointer to store whether the solving should be stopped */
   )
{
   CLIQUE* clique;
   int insertpos;
   Bool acceptsol;

   assert(cliquetable != NULL);
   assert(curcliquenodes != NULL);
   assert(maxcliquenodes != NULL);
   assert(nmaxcliquenodes != NULL);
   assert(maxcliqueweight != NULL);
   assert(curcliqueweight > *maxcliqueweight);
   assert(stopsolving != NULL);

   acceptsol = TRUE;
   *stopsolving = FALSE;
   clique = NULL;
   insertpos = 0;

   if( usrcallback != NULL )
   {
      /* check whether the clique is already stored in the table */
      if( cliquetable->ncliques > 0 )
      {
         createClique(&clique, curcliquenodes, ncurcliquenodes);
         acceptsol = !inCliquetable(cliquetable, clique, &insertpos);
      }
   }

   /* check, if this is a new clique */
   if( acceptsol )
   {
      /* extend the clique with the zero-weighted nodes */
      extendCliqueZeroWeight(tcliquedata, Vzero, nVzero, curcliquenodes, &ncurcliquenodes);

      if( usrcallback != NULL )
      {
         /* call user callback method */
         usrcallback(usrdata, curcliquenodes, ncurcliquenodes, curcliqueweight, maxcliqueweight, &acceptsol, stopsolving);

         /* if clique was accepted, clear the clique table; otherwise, insert it into the clique table, such that
          * the same or a weaker clique is not presented to the user again
          */
         if( acceptsol )
            clearCliquetable(cliquetable);
         else
         {
            /* if the clique was not yet created, do it now */
            if( clique == NULL )
            {
               assert(insertpos == 0);
               assert(cliquetable->ncliques == 0);
               createClique(&clique, curcliquenodes, ncurcliquenodes);
            }
            
            /* insert clique into clique table */
            insertClique(cliquetable, clique, insertpos);
            clique = NULL; /* the clique now belongs to the table */
         }
      }
   }

   /* free the clique, if it was created and not put into the clique table */
   if( clique != NULL )
      freeClique(&clique);

   if( acceptsol )
   {
      /* copy the solution to the incumbent */
      copyMemoryArray(maxcliquenodes, curcliquenodes, ncurcliquenodes);
      *nmaxcliquenodes = ncurcliquenodes;
      if( curcliqueweight > *maxcliqueweight )
         *maxcliqueweight = curcliqueweight;
   }

#ifdef DEBUG
   debugMessage(" -> clique %s (weight %d):", acceptsol ? "accepted" : "rejected", curcliqueweight);
   {
      int i;
      for( i = 0; i < ncurcliquenodes; ++i )
         printf(" %d", curcliquenodes[i]);
      printf("\n");
   }
#endif
}

/** tries to find a clique, if V has only one or two nodes */
static
void reduced( 
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int*             V,                  /**< non-zero weighted nodes for branching */
   int              nV,                 /**< number of non-zero weighted nodes for branching */
   WEIGHT*          apbound,            /**< apriori bound of nodes for branching */ 
   int*             tmpcliquenodes,     /**< buffer for storing the temporary clique */
   int*             ntmpcliquenodes,    /**< pointer to store number of nodes of the temporary clique */
   WEIGHT*          tmpcliqueweight     /**< pointer to store weight of the temporary clique */
   )
{
   WEIGHT* weights;
  
   assert(tcliquedata != NULL);
   assert(V != NULL);
   assert(0 <= nV && nV <= 2);
   assert(apbound != NULL);
   assert(tmpcliquenodes != NULL);
   assert(ntmpcliquenodes != NULL);
   assert(tmpcliqueweight != NULL);

   weights = tcliqueGetWeights(tcliquedata);
   assert(nV == 0 || weights[V[0]] > 0);
   assert(nV <= 1 || weights[V[1]] > 0);

   if( nV >= 1 )
      apbound[0] = weights[V[0]];
   if( nV >= 2 )
      apbound[1] = weights[V[1]];

   /* check if nodes are adjacent */ 
   if( nV >= 2 && tcliqueIsEdge(tcliquedata, V[0], V[1]) )
   {
      assert(tcliqueIsEdge(tcliquedata, V[1], V[0]));

      /* put nodes into clique */
      tmpcliquenodes[0] = V[0];
      tmpcliquenodes[1] = V[1];
      *ntmpcliquenodes = 2;
      *tmpcliqueweight = weights[V[0]] + weights[V[1]];
      apbound[0] += weights[V[1]];
   }
   else if( nV >= 2 && weights[V[1]] > weights[V[0]] )
   {
      /* put V[1] into clique */
      tmpcliquenodes[0] = V[1];
      *ntmpcliquenodes = 1;
      *tmpcliqueweight = weights[V[1]];
   }
   else if( nV >= 1 )
   {
      /* put V[0] into clique */
      tmpcliquenodes[0] = V[0];
      *ntmpcliquenodes = 1;
      *tmpcliqueweight = weights[V[0]];
   }
   else
   {
      *tmpcliqueweight = 0;
      *ntmpcliquenodes = 0;
   }
}

/** calculates upper bound on remaining subgraph, and heuristically generates a clique */
static
WEIGHT boundSubgraph(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   CHKMEM*          mem,                /**< block memory */
   int*             V,                  /**< non-zero weighted nodes for branching */
   int              nV,                 /**< number of non-zero weighted nodes for branching */
   NBC*             gsd,                /**< neighbour color information of all nodes */
   Bool*            iscolored,          /**< coloring status of all nodes */
   WEIGHT*          apbound,            /**< apriori bound of nodes for branching */ 
   int*             tmpcliquenodes,     /**< buffer for storing the temporary clique */
   int*             ntmpcliquenodes,    /**< pointer to store number of nodes of the temporary clique */
   WEIGHT*          tmpcliqueweight     /**< pointer to store weight of the temporary clique */
   )
{
   assert(tmpcliqueweight != NULL);

   /* check if we are in an easy case with at most 2 nodes left */
   if( nV <= 2 )
   {
      /* get 1- or 2-clique and bounds without coloring */
      reduced(tcliquedata, V, nV, apbound, tmpcliquenodes, ntmpcliquenodes, tmpcliqueweight);
      return *tmpcliqueweight;
   }
   else
   {
      /* color the graph induces by nodes of V to get an upper bound for the remaining subgraph */
      return tcliqueColoring(tcliquedata, mem, V, nV, gsd, iscolored, apbound, 
         tmpcliquenodes, ntmpcliquenodes, tmpcliqueweight);
   }
}

/** gets the index of the node of V with the maximum apriori bound; returns -1, if no node exists */
static
int getMaxApBoundIndex( 
   int              nV,  	        /**< number of nodes of V */
   WEIGHT*          apbound             /**< apriori bound of nodes of V */ 
   )
{
   WEIGHT maxapbound;
   int maxindex;
   int i;

   assert(apbound != NULL);

   maxapbound = 0;
   maxindex = -1;

   for( i = 0 ; i < nV; i++ )
   {
      assert(apbound[i] > 0);
      if( apbound[i] >= maxapbound )
      {
         maxapbound = apbound[i];
         maxindex = i;
      }      
   }

   return maxindex;
}

/** gets the index of the node of V with the maximum apriori bound, but ignores nodes with weights
 *  larger than the given maximal weight;
 *  returns -1 if no node with weight smaller or equal than maxweight is found
 */
static
int getMaxApBoundIndexNotMaxWeight( 
   int*             V,                  /**< non-zero weighted nodes for branching */
   int              nV,                 /**< number of non-zero weighted nodes for branching */
   WEIGHT*          apbound,            /**< apriori bound of nodes of V */ 
   WEIGHT*          weights,            /**< weights of nodes */
   WEIGHT           maxweight           /**< maximal weight of node to be candidate for selection */
   )
{
   WEIGHT maxapbound;
   int maxindex;
   int i;

   assert(apbound != NULL);

   maxapbound = 0;
   maxindex = -1;

   for( i = 0 ; i < nV; i++ )
   {
      assert(apbound[i] > 0);
      assert(weights[V[i]] > 0);

      if( apbound[i] >= maxapbound && weights[V[i]] <= maxweight )
      {
         maxapbound = apbound[i];
         maxindex = i;
      }      
   }
   
   return maxindex;
}

/** branches the searching tree, branching nodes are selected in decreasing order of their apriori bound, 
 *  returns whether to stop branching
 */
static
Bool branch( 
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   TCLIQUE_USRCALLBACK ((*usrcallback)),/**< user function to call on every new solution */
   void*            usrdata,            /**< user data to pass to user callback function */
   CHKMEM*          mem,                /**< block memory */
   CLIQUETABLE*     cliquetable,        /**< clique table */
   int              level,		/**< level of b&b tree */
   int*             V,                  /**< non-zero weighted nodes for branching */
   int              nV,                 /**< number of non-zero weighted nodes for branching */
   int*             Vzero,              /**< zero weighted nodes */
   int              nVzero,             /**< number of non-zero weighted nodes */
   NBC*             gsd,                /**< neighbour color information of all nodes */
   Bool*            iscolored,          /**< coloring status of all nodes */
   int*             K,                  /**< nodes from the b&b tree */ 
   WEIGHT           weightK,            /**< weight of the nodes from b&b tree */
   int*             maxcliquenodes,     /**< pointer to store nodes of the maximum weight clique */
   int*             nmaxcliquenodes,    /**< pointer to store number of nodes in the maximum weight clique */
   WEIGHT*          maxcliqueweight,    /**< pointer to store weight of the maximum weight clique */
   int*             curcliquenodes,     /**< pointer to store nodes of currenct clique */
   int*             ncurcliquenodes,    /**< pointer to store number of nodes in current clique */
   WEIGHT*          curcliqueweight,    /**< pointer to store weight of current clique */
   int*             tmpcliquenodes,     /**< buffer for storing the temporary clique */
   WEIGHT           maxfirstnodeweight, /**< maximum weight of branching nodes in level 0; 0 if not used 
                                         **  (for cliques with at least one fractional node) */
   int*             ntreenodes,         /**< pointer to store number of nodes of b&b tree */
   int              maxntreenodes       /**< maximum number of nodes of b&b tree */
   )
{
   Bool stopsolving;
   Bool isleaf;
   WEIGHT* apbound;
   WEIGHT* weights;
   WEIGHT subgraphweight;
   WEIGHT weightKold;
   WEIGHT tmpcliqueweight;
   int ntmpcliquenodes;
   int i;
   
   assert(tcliquedata != NULL);
   assert(mem != NULL);
   assert(V != NULL);
   assert(gsd != NULL);
   assert(iscolored != NULL);
   assert(K != NULL);
   assert(maxcliqueweight != NULL);
   assert(curcliquenodes != NULL);
   assert(ncurcliquenodes != NULL);
   assert(curcliqueweight != NULL);
   assert(ntreenodes != NULL);
   assert(maxfirstnodeweight >= 0);
   assert(*ntreenodes >= 0);
   assert(maxntreenodes >= 0);

   /* increase the number of nodes, and stop solving, if the node limit is exceeded */
   (*ntreenodes)++;
   debugMessage("(level %d, treenode %d) maxclique = %d, curclique = %d [mem=%lld (%lld), cliques=%d]\n",
      level, *ntreenodes, *maxcliqueweight, *curcliqueweight, 
      getChunkMemoryUsed(mem), getMemoryUsed(), cliquetable->ncliques);
   debugMessage(" -> current branching (weight %d):", weightK);
#ifdef DEBUG
   for( i = 0; i < level; ++i )
      printf(" %d", K[i]);
   printf("\n");
#endif
   debugMessage(" -> branching candidates:");
#ifdef DEBUG
   for( i = 0; i < nV; ++i )
      printf(" %d", V[i]);
   printf("\n");
#endif
   if( *ntreenodes > maxntreenodes )
      return TRUE;

   weights = tcliqueGetWeights(tcliquedata);
   stopsolving = FALSE;
   isleaf = TRUE;

   /* allocate temporary memory for a priori bounds */
   ALLOC_ABORT( allocMemoryArray(&apbound, nV) );
   clearMemoryArray(apbound, nV);

   /* use coloring relaxation to generate an upper bound for the current subtree and a heuristic solution */
   subgraphweight = boundSubgraph(tcliquedata, mem, V, nV, gsd, iscolored, apbound,
      tmpcliquenodes, &ntmpcliquenodes, &tmpcliqueweight);

#ifndef NDEBUG
   /* check correctness of V and apbound arrays */
   for( i = 0; i < nV; ++i )
   {
      assert(0 <= V[i] && V[i] < tcliqueGetNNodes(tcliquedata));
      assert(i == 0 || V[i-1] < V[i]);
      assert(apbound[i] >= 0);
      assert((apbound[i] == 0) == (weights[V[i]] == 0));
   }
#endif

   /* check, whether the heuristic solution is better than the current subtree's solution */
   if( weightK + tmpcliqueweight > *curcliqueweight )
   {
      /* install the newly generated clique as current clique */
      for( i = 0; i < level; ++i )
         curcliquenodes[i] = K[i];
      for( i = 0; i < ntmpcliquenodes; ++i )
         curcliquenodes[level+i] = tmpcliquenodes[i];
      *ncurcliquenodes = level + ntmpcliquenodes;
      *curcliqueweight = weightK + tmpcliqueweight;

#ifdef DEBUG
      debugMessage(" -> new current clique with weight %d at node %d in level %d:",
         *curcliqueweight, *ntreenodes, level);
      for( i = 0; i < *ncurcliquenodes; ++i )
         printf(" %d", curcliquenodes[i]);
      printf("\n");
#endif
   }

   /* discard subtree, if the upper bound is not better than the weight of the currently best clique;
    * if only 2 nodes are left, the maximal weighted clique was already calculated in boundSubgraph() and nothing
    * more has to be done
    */
   if( weightK + subgraphweight > *maxcliqueweight && nV > 2 )
   {
      int* Vcurrent;
      int nVcurrent;
      int nValive;
      int branchingnode;

      assert(nV > 0);

      /* process current subtree */
      level++;
   
      /* set up data structures */
      ALLOC_ABORT( allocMemoryArray(&Vcurrent, nV-1) );

      nValive = nV;
      weightKold = weightK;

      debugMessage("============================ branching level %d ===============================\n", level); 

      /* branch on the nodes of V by decreasing order of their apriori bound */
      while( !stopsolving && nValive > 0 )
      {
         int branchidx;

         /* get next branching node */
         if( level == 1 && maxfirstnodeweight > 0 ) 
            branchidx = getMaxApBoundIndexNotMaxWeight(V, nValive, apbound, weights, maxfirstnodeweight);
         else
            branchidx = getMaxApBoundIndex(nValive, apbound);
         if( branchidx < 0 )
            break;
         assert(0 <= branchidx && branchidx < nValive && nValive <= nV);
         assert(apbound[branchidx] > 0);
         assert(weights[V[branchidx]] > 0);

         /* test a priori bound */
         if( (weightKold + apbound[branchidx]) <= *maxcliqueweight )
            break;

         debugMessage("%d. branching in level %d: bidx=%d, node %d, weight %d, upperbound: %d+%d = %d, maxclique=%d\n", 
            nV-nValive+1, level, branchidx, V[branchidx], weights[V[branchidx]], weightKold, apbound[branchidx], 
            weightKold + apbound[branchidx], *maxcliqueweight); 

         /* because we branch on this node, the node is no leaf in the tree */
         isleaf = FALSE;

         /* update the set of nodes from the b&b tree 
          *   K = K & {branchingnode}
          */
         branchingnode = V[branchidx];
         K[level-1] = branchingnode;
         weightK = weightKold + weights[branchingnode];
 
         /* update the set of nodes for branching 
          *   V = V \ {branchingnode}
          */
         nValive--;
         for( i = branchidx; i < nValive; ++i )
         {
            V[i] = V[i+1];
            apbound[i] = apbound[i+1];
         }

         /* set the nodes for the next level of b&b tree 
          *   Vcurrent = nodes of V, that are adjacent to branchingnode
          */
         nVcurrent = tcliqueSelectAdjnodes(tcliquedata, branchingnode, V, nValive, Vcurrent);
      
         /* process the selected subtree */
         stopsolving = branch(tcliquedata, usrcallback, usrdata, mem, cliquetable, level, 
            Vcurrent, nVcurrent, Vzero, nVzero, gsd, iscolored, K, weightK,
            maxcliquenodes, nmaxcliquenodes, maxcliqueweight, 
            curcliquenodes, ncurcliquenodes, curcliqueweight, tmpcliquenodes,
            maxfirstnodeweight, ntreenodes, maxntreenodes);
      }
   
      debugMessage("========================== branching level %d end =============================\n\n", level); 

      /* free data structures */
      freeMemoryArray(&Vcurrent);
   }

   /* check, whether any branchings have been applied, or if this node is a leaf of the branching tree */
   if( isleaf )
   {
      /* the current clique is the best clique found on the path to this leaf
       * -> check, whether it is an improvement to the currently best clique
       */
      if( *curcliqueweight > *maxcliqueweight )
      {
         debugMessage("found clique of weight %d at node %d in level %d\n", *curcliqueweight, *ntreenodes, level);
         newSolution(tcliquedata, cliquetable, Vzero, nVzero,
            curcliquenodes, *ncurcliquenodes, *curcliqueweight, 
            maxcliquenodes, nmaxcliquenodes, maxcliqueweight, usrcallback, usrdata, &stopsolving);
      }

      /* discard the current clique */
      *ncurcliquenodes = 0;
      *curcliqueweight = 0;
   }

   /* free data structures */   
   freeMemoryArray(&apbound);

   return stopsolving;
}

/** finds maximum weight clique */
void tcliqueMaxClique(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   TCLIQUE_USRCALLBACK ((*usrcallback)),/**< user function to call on every new solution */
   void*            usrdata,            /**< user data to pass to user callback function */
   int*             maxcliquenodes,     /**< pointer to store nodes of the maximum weight clique */
   int*             nmaxcliquenodes,    /**< pointer to store number of nodes in the maximum weight clique */
   WEIGHT*          maxcliqueweight,    /**< pointer to store weight of the maximum weight clique */
   WEIGHT           maxfirstnodeweight, /**< maximum weight of branching nodes in level 0; 0 if not used 
                                         *   for cliques with at least one fractional node) */
   WEIGHT           minweight,          /**< lower bound for weight of generated cliques */
   int              maxntreenodes 	/**< maximum number of nodes of b&b tree */
   )
{
   CLIQUETABLE* cliquetable;
   WEIGHT* weights;
   int* K;	
   int* V;
   int* Vzero;
   int nnodes;
   int nV;
   int nVzero;
   int i;
   CHKMEM* mem;
   NBC* gsd;
   Bool* iscolored;
   int* curcliquenodes;
   int ncurcliquenodes;
   WEIGHT curcliqueweight;
   int* tmpcliquenodes;
   int ntreenodes;

   assert(tcliquedata != NULL);
   assert(maxcliquenodes != NULL);
   assert(nmaxcliquenodes != NULL);
   assert(maxcliqueweight != NULL);
   assert(maxntreenodes >= 0);

   debugMessage("calculating maximal weighted clique in graph (%d nodes, %d edges)\n", 
      tcliqueGetNNodes(tcliquedata), tcliqueGetNEdges(tcliquedata));

   nnodes = tcliqueGetNNodes(tcliquedata);

   /* set up data structures */
   createCliquetable(&cliquetable, CLIQUETABLE_INITSIZE);
   ALLOC_ABORT( allocMemoryArray(&K, nnodes) );
   ALLOC_ABORT( allocMemoryArray(&V, nnodes) );
   ALLOC_ABORT( allocMemoryArray(&Vzero, nnodes) );
   ALLOC_ABORT( allocMemoryArray(&gsd, nnodes) );
   ALLOC_ABORT( allocMemoryArray(&iscolored, nnodes) );
   ALLOC_ABORT( allocMemoryArray(&curcliquenodes, nnodes) );
   ALLOC_ABORT( allocMemoryArray(&tmpcliquenodes, nnodes) );

   /* set weight and number of nodes of maximum weighted clique */
   *nmaxcliquenodes = 0; 
   *maxcliqueweight = minweight-1;
   ncurcliquenodes = 0;
   curcliqueweight = 0;
   ntreenodes = 0;

   /* set up V and Vzero */
   weights = tcliqueGetWeights(tcliquedata);
   assert(weights != NULL);
   nV = 0;
   nVzero = 0;
   for( i = 0 ; i <  nnodes; i++ )
   {
      if( weights[i] == 0 )
      {
         Vzero[nVzero] = i;
         nVzero++;
      }
      else
      {
         V[nV] = i;
         nV++;
      }
   }

   /* initialize own memory allocator for coloring */ 
   mem = createChunkMemory(sizeof(LIST_ITV), CHUNK_SIZE, -1); 

   /* branch to find maximum weight clique */
   (void)branch(tcliquedata, usrcallback, usrdata, mem, cliquetable,
      0, V, nV, Vzero, nVzero, gsd, iscolored, K, 0, 
      maxcliquenodes, nmaxcliquenodes, maxcliqueweight, 
      curcliquenodes, &ncurcliquenodes, &curcliqueweight, tmpcliquenodes,
      maxfirstnodeweight, &ntreenodes, maxntreenodes);

   /* delete own memory allocator for coloring */ 
   destroyChunkMemory(&mem); 
   
   /* free data structures */
   freeMemoryArray(&tmpcliquenodes);
   freeMemoryArray(&curcliquenodes);
   freeMemoryArray(&iscolored); 
   freeMemoryArray(&gsd); 
   freeMemoryArray(&Vzero);
   freeMemoryArray(&V);
   freeMemoryArray(&K);
   freeCliquetable(&cliquetable);
}
