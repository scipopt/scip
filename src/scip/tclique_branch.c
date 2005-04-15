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
#pragma ident "@(#) $Id: tclique_branch.c,v 1.1 2005/04/15 11:46:54 bzfpfend Exp $"

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


/** tries to find a clique, if V has only one or two nodes */
static
void reduced( 
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              level,		/**< level of b&b tree */
   int*             V,                  /**< nodes for branching */ 
   int*             nV,		        /**< number of nodes for branching */
   int*             K,                  /**< nodes from the b&b tree (nK = level + 1) */ 
   WEIGHT           weightK,            /**< weight of the nodes from b&b tree */
   int*             mwc, 	        /**< pointer to store nodes of the maximum weight clique */
   int*             nmwc,	        /**< pointer to store number of nodes in the maximum weight clique */
   WEIGHT*          weightmwc           /**< pointer to store weight of the maximum weight clique */
   )
{
   WEIGHT* weights;
   WEIGHT weightclique;
   int i;
  
   assert(tcliquedata != NULL);
   assert(V != NULL);
   assert(*nV > 0 && *nV < 3);
   assert(K != NULL);
   assert(weightK >= 0);
   assert(mwc != NULL);
   assert(nmwc != NULL);
   assert(weightmwc != NULL);

   weights = tcliqueGetWeights(tcliquedata);
   weightclique = 0.;

   /* V has only two nodes */
   if( *nV == 2 )
   {
      /* checks if nodes are adjacent */ 
      if( tcliqueIsEdge(tcliquedata, V[0], V[1]) )
      {
         assert(tcliqueIsEdge(tcliquedata, V[1], V[0]));
         weightclique = weights[V[1]];
      }
      else
      {	
         assert(!tcliqueIsEdge(tcliquedata, V[1], V[0]));
         /* updates V 
          *   V = {node with maximum weight}
          */
         if( weights[V[1]] > weights[V[0]] )
            V[0] = V[1]; 
         *nV = 1;
      }
   }
   
   weightclique += weights[V[0]];
   
   /* improves global lower bound */
   if( (weightK + weightclique) > *weightmwc ) 
   {
      /* updates maximum weight clique 
       *   mwc = K & V
       */
      for( i = 0; i < level + 1; i++ )
         mwc[i] = K[i];
      for( i = 0; i < *nV; i++ )
         mwc[level+1+i] = V[i];
      *weightmwc = weightK + weightclique;
      *nmwc = level + 1 + *nV;
   }
}

/** gets the index of the node of V with the maximum apriori bound */  
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

   maxapbound = 0.;
   maxindex = -1;

   for( i = 0 ; i < nV; i++ )
   {
      if( apbound[i] >= maxapbound )
      {
         maxapbound = apbound[i];
         maxindex = i;
      }      
   }
   
   assert(maxindex > -1);

   return maxindex;
}

/** gets the index of the node of V with the maximum apriori bound, ignores nodes with weights larger than the
 *  given maximal weight; returns -1 if no node with weight smaller or equal than maxweight is found
 */
static
int getMaxApBoundIndexNotMaxWeight( 
   int              nV,  	        /**< number of nodes of V */
   WEIGHT*          apbound,            /**< apriori bound of nodes of V */ 
   WEIGHT*          weights,            /**< weights of nodes */
   WEIGHT           maxweight           /**< maximal weight of node to be candidate for selection */
   )
{
   WEIGHT maxapbound;
   int maxindex;
   int i;

   assert(apbound != NULL);

   maxapbound = 0.;
   maxindex = -1;

   for( i = 0 ; i < nV; i++ )
   {
      if( apbound[i] >= maxapbound && weights[i] <= maxweight )
      {
         maxapbound = apbound[i];
         maxindex = i;
      }      
   }
   
   return maxindex;
}

/** tries to improve the lower bound with the clique obtained in the coloring step. 
 *  also tries to improve the upper bound needed to decide, whether to discard the subproblem. 
 *  returns decision
 */
static
BOOL bound( 
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   CHKMEM*          mem,                /**< block memory */
   int              level,		/**< level of b&b tree */
   int*             V,                  /**< nodes for branching */ 
   int              nV,		        /**< number of nodes for branching */
   NBC*             gsd,                /**< neighbour color information of all nodes */
   BOOL*            iscolored,          /**< coloring status of all nodes */
   WEIGHT*          apbound,            /**< apriori bound of nodes for branching */ 
   int*             K,                  /**< nodes from the b&b tree (nK = level + 1) */ 
   WEIGHT           weightK,            /**< weight of the nodes from b&b tree */
   int*             mwc, 	        /**< pointer to store nodes of the maximum weight clique */
   int*             nmwc,	        /**< pointer to store number of nodes in the maximum weight clique */
   WEIGHT*          weightmwc,          /**< pointer to store weight of the maximum weight clique */
   int*             bestlevel		/**< pointer to store level of b&b tree where the current maximum weight clique was found */
   )
{
   WEIGHT weightclique;
   int* clique;	
   int nclique; 
   int ncolors;
   int i;
        
   assert(tcliquedata != NULL);
   assert(mem != NULL);
   assert(V != NULL);
   assert(gsd != NULL);
   assert(iscolored != NULL);
   assert(apbound != NULL);
   assert(K != NULL);
   assert(mwc != NULL);
   assert(nmwc != NULL);
   assert(weightmwc != NULL);
   assert(bestlevel != NULL);

   /* gets clique and bounds without coloring */
   if( nV < 3 )
   {
      reduced(tcliquedata, level, V, &nV, K, weightK, mwc, nmwc, weightmwc);
   }

   /* sets data structures */
   ALLOC_ABORT( allocMemoryArray(&clique, nV) );

   /* colors the graph induces by nodes of V */
   ncolors = tcliqueColoring(tcliquedata, mem, V, nV, gsd, iscolored, apbound, &clique, &nclique, &weightclique);

   /* improves global lower bound */
   if( (weightK + weightclique) > *weightmwc ) 
   {
      /* updates maximum weight clique 
       *   mwc = K & clique */
      for( i = 0; i < level + 1; i++ )
      {
         mwc[i] = K[i];
      }
      for( i = 0; i < nclique; i++ )
      {
         mwc[level+1+i] = clique[i];
      }
      *weightmwc = weightK + weightclique;
      *nmwc = level + 1 + nclique;
      *bestlevel = level+1;
   }
   
#ifdef DEBUG
   debugMessage("mwc=[ ");
   for( i = 0; i < *nmwc; i++ )
      printf("%d, ", mwc[i]);
   printf("] weightmwc=%d\n\n", *weightmwc);
#endif

   /* frees data structures */
   freeMemory(&clique);
    
   /* tests to discard subproblem */
   if( (weightK + ncolors) <= *weightmwc )
      return NO;
   
   return YES;
}

/** sorting function quicksort for sorting given set of nodes by increasing key */ 
static
void quicksortWeight(
   int**            nodes,              /**< given set of nodes */
   WEIGHT*          key,                /**< sorting key of nodes */
   int              pi,                 /**< index of first node in nodes */
   int              pf                  /**< index of last node in nodes */
   )
{
   int ii;
   int jj;
   int i;
   int j;
   int tmp;
   
   if( pi >= pf ) 
      return;
   
   ii = 1;
   jj = 0;
   i = pi;
   j = pf;

   while( i != j )
   {
      if( key[(*nodes)[i]] > key[(*nodes)[j]] )
      {
         tmp = (*nodes)[i];
         (*nodes)[i] = (*nodes)[j];
         (*nodes)[j] = tmp;
              
         ii = !ii;
         jj = !jj;
      }	
      i += ii;
      j -= jj;
           
   }

   if( i-pi > 1 )
      quicksortWeight(nodes, key, pi, i-1);

   if( pf-i > 1 )
      quicksortWeight(nodes, key, i+1, pf);
}

/** sorts a given set of nodes by increasing key */ 
static
void sortWeightInc(
   int              nnodes,             /**< number of given nodes */ 
   int**            nodes,              /**< given set of nodes */
   WEIGHT*          key                 /**< sorting key of nodes */
   )
{
   if( nnodes > 1 )
      quicksortWeight(nodes, key, 0, nnodes-1);
}

/** calls user callback after a new best solution was found */
static
void newSolution(
   int*             mwc, 	        /**< pointer to store nodes of the maximum weight clique */
   int              nmwc,	        /**< pointer to store number of nodes in the maximum weight clique */
   WEIGHT           weightmwc,          /**< weight of clique */
   int*             bestlevel,          /**< pointer to level where clique was found */
   TCLIQUE_USRCALLBACK ((*usrcallback)),/**< user function to call on every new solution */
   void*            usrdata,            /**< user data to pass to user callback function */
   BOOL*            stopsolving         /**< pointer to store whether the solving should be stopped */
   )
{
   assert(mwc != NULL);
   assert(bestlevel != NULL);
   assert(stopsolving != NULL);

   *stopsolving = NO;

   if( usrcallback != NULL )
      *stopsolving = usrcallback(mwc, nmwc, weightmwc, usrdata);

   *bestlevel = -1;
}

/** branches the searching tree, branching nodes are selected in decreasing order of their apriori bound, 
 *  returns whether to stop branching
 */
static
BOOL branch( 
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   TCLIQUE_USRCALLBACK ((*usrcallback)),/**< user function to call on every new solution */
   void*            usrdata,            /**< user data to pass to user callback function */
   CHKMEM*          mem,                /**< block memory */
   int*             level,		/**< level of b&b tree */
   int*             V,                  /**< nodes for branching */ 
   int              nV,		        /**< number of nodes for branching */
   NBC*             gsd,                /**< neighbour color information of all nodes */
   BOOL*            iscolored,          /**< coloring status of all nodes */
   WEIGHT*          apbound,            /**< apriori bound of nodes for branching */ 
   int*             K,                  /**< nodes from the b&b tree */ 
   WEIGHT*          weightK,            /**< weight of the nodes from b&b tree */
   int*             mwc, 	        /**< pointer to store nodes of the maximum weight clique */
   int*             nmwc,	        /**< pointer to store number of nodes in the maximum weight clique */
   WEIGHT*          weightmwc,          /**< pointer to store weight of the maximum weight clique */
   WEIGHT           maxweight,          /**< maximum weight of branching nodes in level 0; 0 if not used 
                                         ** (for cliques with at least one fractional node) */
   int*             bestlevel,    	/**< pointer to store level of b&b tree where the current maximum weight clique was found */
   int*             treenodes,    	/**< pointer to store number of nodes of b&b tree */
   int              maxtreenodes 	/**< maximum number of nodes of b&b tree */
   )
{
   BOOL stopsolving;
   BOOL stop;
   WEIGHT* apboundcurrent; 
   WEIGHT* weights;
   WEIGHT weightKold;
   int* Vcurrent;	
   int nVcurrent;
   int nValive;
   int branchingnode;
   int i;
   
   assert(tcliquedata != NULL);
   assert(mem != NULL);
   assert(V != NULL);
   assert(gsd != NULL);
   assert(iscolored != NULL);
   assert(apbound != NULL);
   assert(K != NULL);
   assert(weightK != NULL);
   assert(mwc != NULL);
   assert(nmwc != NULL);
   assert(weightmwc != NULL);
   assert(bestlevel != NULL);
   assert(treenodes != NULL);
   assert(maxweight >= 0);
   assert(*bestlevel >= -1);
   assert(*treenodes >= 0);
   assert(maxtreenodes >= 0);

   weights = tcliqueGetWeights(tcliquedata);

   /* sets data structures */
   ALLOC_ABORT( allocMemoryArray(&Vcurrent, nV-1) );
   ALLOC_ABORT( allocMemoryArray(&apboundcurrent, nV-1) );

   nVcurrent = nV;
   nValive = nV;

   weightKold = *weightK;
   (*level)++;
   stopsolving = NO;
   stop = NO;

   if( *treenodes > maxtreenodes )
      stop = YES;

   debugMessage("============================ branching level %d ===============================\n\n", *level); 

   /* branches on the nodes of V by decreasing order of their apriori bound */
   for( i = 0; !stop && i < nV; i++ )
   {
      int j;

      if( *level == 0 && maxweight > 0 ) 
      {
         j = getMaxApBoundIndexNotMaxWeight(nValive, apbound, weights, maxweight-1);
         if( j < 0 )
            break;
      }
      else
         j = getMaxApBoundIndex(nValive, apbound);
      assert(0 <= j && j < nV);
     
      /* tests a priori bound */
      if( (weightKold + apbound[j]) <= *weightmwc )
         break;

      branchingnode = V[j];
      (*treenodes)++;

#ifdef DEBUG
      {
         int k;
         debugMessage("apbound=[ ");
         for( k = 0; k < nValive; k++ )
            printf("%d:%d, ", k, apbound[k]);
         printf("]\n");
         debugMessage("branching node: vindex=%d vertex=%d)\n\n", j, V[j] ); 
      }
#endif      

      /* updates the set of nodes from the b&b tree 
       *   K = K & {branchingnode}
       */
      K[*level] = branchingnode;
      *weightK = weightKold + weights[branchingnode];
 
      /* update the set of nodes for branching 
       *   V = V \ {branchingnode}
       */
      nVcurrent--;
      nValive--;
      for( ; j < nValive; j++ )
      {
         V[j] = V[j+1];
         apbound[j] = apbound[j+1];
      }

      if( *treenodes > maxtreenodes )
      {
         stop = YES;
         break;
      }

      /* sets the nodes for the next level of b&b tree 
       *   V = nodes of V, that are adjacent to branchingnode
       */
      nVcurrent = tcliqueSelectAdjnodes(tcliquedata, branchingnode, V, nValive, Vcurrent);  
      
      /* branchingnode has no adjacent node in V */
      if( nVcurrent == 0) 
         continue;
      
      /* tests if branching on subproblem is usefull (uses coloring) 
       * and discrades subproblem if not */
      if( bound(tcliquedata, mem, *level, Vcurrent, nVcurrent, gsd, iscolored, apboundcurrent, K, *weightK, 
            mwc, nmwc, weightmwc, bestlevel) == NO )
      {
         if( *bestlevel == *level+1 )
         {
            newSolution(mwc, *nmwc, *weightmwc, bestlevel, usrcallback, usrdata, &stopsolving);
            debugMessage("treenodes=%d\n", *treenodes);

            /* stops branching on nodes of V (decision to stop made in this b&b node) */
            if( stopsolving )
            {
               stop = YES;
               break;
            }
         }
         continue ;
      }
      
      /* stops branching on nodes of V (decision to stop comes from a child b&b node) */
      if( branch(tcliquedata, usrcallback, usrdata, mem, level, Vcurrent, nVcurrent, gsd, iscolored, apboundcurrent, 
            K, weightK, mwc, nmwc, weightmwc, maxweight, bestlevel, treenodes, maxtreenodes)
         || *treenodes > maxtreenodes )
      {
         stop = YES;
         break;
      }
   }
   
   /* frees data structures */   
   freeMemory(&apboundcurrent);
   freeMemory(&Vcurrent);

   debugMessage("========================== branching level %d end =============================\n\n", *level); 

   /* check, if we found a new best solution */
   assert(*bestlevel <= *level);
   if( *bestlevel == *level )
   {
      newSolution(mwc, *nmwc, *weightmwc, bestlevel, usrcallback, usrdata, &stopsolving);
      debugMessage("treenodes=%d\n", *treenodes);
      if( stopsolving )
         stop = YES;
   }

   /* goes up in the searching tree */
   (*level)--; 
   *weightK = weightKold;

   debugMessage("treenodes=%d\n", *treenodes);

   return stop;
}

/** finds maximum weight clique */
void tcliqueMaxClique(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   TCLIQUE_USRCALLBACK ((*usrcallback)),/**< user function to call on every new solution */
   void*            usrdata,            /**< user data to pass to user callback function */
   int*             mwc, 	        /**< pointer to store nodes of the maximum weight clique */
   int*             nmwc,	        /**< pointer to store number of nodes in the maximum weight clique */
   WEIGHT*          weightmwc,          /**< pointer to store weight of the maximum weight clique */
   WEIGHT           maxfirstnodeweight, /**< maximum weight of branching nodes in level 0; 0 if not used 
                                         *   for cliques with at least one fractional node) */
   WEIGHT           minweight,          /**< lower bound for weight of generated cliques */
   int              maxtreenodes 	/**< maximum number of nodes of b&b tree */
   )
{
   WEIGHT* apbound;
   WEIGHT* weights;
   WEIGHT weightK;	
   WEIGHT weightsum;
   int* K;	
   int* V;
   int nV;
   int nnodes;
   int level;
   int i;
   CHKMEM* mem;
   NBC* gsd;
   BOOL* iscolored;
   int bestlevel;
   int treenodes;

   assert(tcliquedata != NULL);
   assert(mwc != NULL);
   assert(nmwc != NULL);
   assert(weightmwc != NULL);
   assert(maxtreenodes >= 0);

   nnodes = tcliqueGetNNodes(tcliquedata);
   level = -1;
   weightK = 0;

   /* sets data structures */
   ALLOC_ABORT( allocMemoryArray(&K, nnodes) );
   ALLOC_ABORT( allocMemoryArray(&V, nnodes) );
   ALLOC_ABORT( allocMemoryArray(&gsd, nnodes) );
   ALLOC_ABORT( allocMemoryArray(&iscolored, nnodes) );
   ALLOC_ABORT( allocMemoryArray(&apbound, nnodes) );

   /* set weight and number of nodes of maximum weighted clique */
   *nmwc = 0; 
   *weightmwc = minweight-1;
   bestlevel = -1;
   treenodes = 0;

   /* sets V */
   for( i = 0 ; i <  nnodes; i++ )
   {
      V[i] = i;
   }

   /* sets apriori bounds of all nodes 
    *    sort V, so that weights(v_i-1) <= weights(v_i) 
    *    apbound(v_i) = weights(v_0) + ... + weights(v_i)
    */

   /* sorts V by non decreasing weights */
   weights = tcliqueGetWeights(tcliquedata);
   sortWeightInc(nnodes, &V, weights);

   /* sets apbound of all nodes */
   weightsum = 0.;
   for( i = 0; i < nnodes; i++ )
   {
      weightsum += weights[V[i]];
      apbound[V[i]] = weightsum;
   }

   /* sorts V by increasing order of node index */ 
   for( i = 0; i < nnodes; i++ )
   {
      V[i] = i;
   }
   nV = nnodes;

   /* initialises own memory allocator for coloring */ 
   mem = createChunkMemory(sizeof(LIST_ITV), CHUNK_SIZE, -1); 

   /* branches to find maximum weight clique */
   branch(tcliquedata, usrcallback, usrdata, mem, &level, V, nV, gsd, iscolored, apbound, K, &weightK, 
      mwc, nmwc, weightmwc, maxfirstnodeweight, &bestlevel, &treenodes, maxtreenodes);

   /* deletes own memory allocator for coloring */ 
   destroyChunkMemory(&mem); 
   
   /* frees data structures */
   freeMemory(&apbound);
   freeMemory(&iscolored); 
   freeMemory(&gsd); 
   freeMemory(&V);
   freeMemory(&K);
}
