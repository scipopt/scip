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
#pragma ident "@(#) $Id: tclique_coloring.c,v 1.4 2005/04/29 13:50:50 bzfpfend Exp $"

/**@file   tclique_coloring.c
 * @brief  coloring part of algorithm for maximum cliques
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <memory.h>
#include <assert.h>
#include <stdlib.h>

#include "scip/def.h"
#include "scip/memory.h"
#include "scip/message.h"
#include "scip/tclique_graph.h"
#include "scip/tclique_coloring.h"



/** gets index of the uncolored node in a given array of nodes in V with maximum satdeg and positive weight;
 *  in case of a tie choose node with maximum weight;
 *  if no uncolored node with positive weight is found, -1 is returned
 */
static
int getMaxSatdegIndex(
   int*             V,                  /**< nodes in V */ 
   int              nV,		        /**< number of nodes in V */
   NBC*             gsd,                /**< neighbor color information of all nodes */
   BOOL*            iscolored,          /**< coloring status of all nodes */
   WEIGHT*          weights             /**< weight of nodes in grpah */
   )
{   
   WEIGHT maxweight;
   int maxsatdeg;
   int maxsatdegindex;
   int i;
	
   maxweight = -1;
   maxsatdeg = -1;
   maxsatdegindex = -1;
   
   assert(gsd != NULL);
   assert(iscolored != NULL);

   for( i = 0; i < nV; i++ )
   {
      WEIGHT weight;
      int satdeg;

      /* check only uncolored nodes */ 
      if( iscolored[i] ) 
         continue;

      /* check only nodes with positive weight */
      weight = weights[V[i]];
      assert(weight >= 0);
      if( weight <= 0 )
         continue;

      satdeg = gsd[i].satdeg;
      if( satdeg > maxsatdeg || (satdeg == maxsatdeg && weight > maxweight) )
      {
         maxsatdeg = satdeg;
         maxweight = weight;
         maxsatdegindex = i;
      }
   }

   return maxsatdegindex;	
}

/** gets index of the node in a given set of nodes with maximum weight */
static
int getMaxWeightIndex( 
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int*             nodes,              /**< given set of nodes (ordered by node index) */ 
   int              nnodes              /**< number of nodes in given set nodes */ 
   )
{
   WEIGHT* weights;
   WEIGHT maxweight;
   int maxweightindex;
   int i;
   
   assert(tcliquedata != NULL);
   assert(nnodes > 0);

   weights = tcliqueGetWeights(tcliquedata);

   maxweightindex = -1;
   maxweight = -1;

   /* try to improve maxweight */
   for( i = 0 ; i < nnodes; i++ )
   {
      assert(0 <= nodes[i] && nodes[i] < tcliqueGetNNodes(tcliquedata));
      assert(weights[nodes[i]] >= 0);
      if( weights[nodes[i]] > maxweight)
      {
         /* node has larger weight */
         maxweight = weights[nodes[i]];
         maxweightindex = i;
      }
   }
   assert(maxweightindex >= 0);

   return maxweightindex;
}

/** updates the neighbor colors information of a node: updates the list of neighbor color intervals 
 *  by making the union of the existing list and the given list of color intervals, and updates the saturation degree
 */
static
void updateNeighbor(
   CHKMEM*          mem,                /**< block memory */
   NBC*             pgsd,               /**< pointer to neighbor color information of node to update */
   LIST_ITV*        pnc                 /**< pointer to given list of color intervals */
   )
{
   LIST_ITV head;
   LIST_ITV* apciv;
   LIST_ITV* pciv;
   LIST_ITV* nciv;

   /* save the pointer to the first element of the list */
   head.next = pgsd->lcitv;
   apciv = &head;
   pciv = apciv->next;
   
   /* construct the union of the two intervals */
   while( (pnc != NULL) && (pciv != NULL) )
   {
      if( pnc->itv.inf < pciv->itv.inf ) 
      {	
         ALLOC_ABORT( allocChunkMemory(mem, &nciv) );
         nciv->itv = pnc->itv;
         nciv->next = pciv;
         apciv->next = nciv;
         apciv = nciv;
         
         pnc = pnc->next;	
      }
      else if( pnc->itv.inf <= pciv->itv.sup )
      {
         if( pnc->itv.sup > pciv->itv.sup )
            pciv->itv.sup = pnc->itv.sup;
         pnc = pnc->next;
      }
      else
      {
         apciv = pciv;
         pciv = pciv->next;
      }
   }

   while( pnc != NULL )
   {
      ALLOC_ABORT( allocChunkMemory(mem, &nciv) );
      nciv->itv = pnc->itv;
      nciv->next = NULL;
      
      apciv->next = nciv;
      apciv = nciv;
      
      pnc = pnc->next;
   }

   /* try to reduce the number of intervals */
   pgsd->satdeg = 0;
   apciv = head.next;
   while( (pciv = apciv->next) != NULL )
   {
      if( apciv->itv.sup < (pciv->itv.inf - 1) )
      {
         pgsd->satdeg += apciv->itv.sup - apciv->itv.inf + 1;
         apciv = apciv->next;
      }
      else
      {
         LIST_ITV* tmp;

         if( apciv->itv.sup < pciv->itv.sup )
            apciv->itv.sup = pciv->itv.sup;
         apciv->next = pciv->next;
       
         /* free data structure for created colorinterval */
         tmp = pciv->next; 
         freeChunkMemory(mem, &pciv); 
         pciv = tmp; 
      }
   }
   pgsd->satdeg += apciv->itv.sup - apciv->itv.inf + 1;
   
   /* updates the pointer to the first element of the list */		
   pgsd->lcitv = head.next;
}

/** colors the positive weighted nodes of a given set of nodes V with the lowest possible number of colors and 
 *  finds a clique in the graph induced by V, an upper bound and an apriori bound for further branching steps
 */
WEIGHT tcliqueColoring( 
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   CHKMEM*          mem,                /**< block memory */
   int*             V,                  /**< nodes for branching */ 
   int              nV,		        /**< number of nodes for branching */
   NBC*             gsd,                /**< neighbor color information of all nodes */
   BOOL*            iscolored,          /**< coloring status of all nodes */
   WEIGHT*          apbound,            /**< pointer to store apriori bound of nodes for branching */ 
   int*             clique,             /**< buffer for storing the clique */
   int*             nclique,            /**< pointer to store number of nodes in the clique */
   WEIGHT*          weightclique        /**< pointer to store the weight of the clique */
   )
{
   WEIGHT* weights;
   WEIGHT maxsatdegree; 
   WEIGHT range;
   BOOL growclique; 
   int node; 
   int nodeVindex;
   int i;     
   int j;
   LIST_ITV* colorinterval;
   LIST_ITV nwcitv;
   LIST_ITV* pnc;
   LIST_ITV* lcitv;
   LIST_ITV* item;
   LIST_ITV* tmpitem;
   int* workclique;
   int* currentclique;
   int ncurrentclique;
   int weightcurrentclique;
   int* currentadjedge;
   int* lastadjedge;

   assert(V != NULL);
   assert(nV > 0);
   assert(clique != NULL);
   assert(nclique != NULL);
   assert(weightclique != NULL);
   assert(gsd != NULL);
   assert(iscolored != NULL);

   weights = tcliqueGetWeights(tcliquedata);

   /* initialize maximum weight clique found so far */
   growclique = YES;
   *nclique = 0;
   *weightclique = 0;

   /* get node of V with maximum weight */
   nodeVindex = getMaxWeightIndex(tcliquedata, V, nV);
   node = V[nodeVindex];
   assert(0 <= node && node < tcliqueGetNNodes(tcliquedata));
   range = weights[node];
   assert(range >= 0);
   if( range == 0 ) /* there are no nodes with positive weight: nothing has to be done */
      return 0;

   /* set up data structures for coloring */
   clearMemoryArray(iscolored, nV); /* new-memory */
   clearMemoryArray(gsd, nV); /* new-memory */
   iscolored[nodeVindex] = YES;

   /* color the first node */
   debugMessage("---------------coloring-----------------\n");
   debugMessage("1. node choosen: vindex=%d, vertex=%d, satdeg=%d, range=%d)\n",
      nodeVindex, node, gsd[nodeVindex].satdeg, range);

   /* set apriori bound: apbound(v_i) = satdeg(v_i) + weight(v_i) */
   apbound[nodeVindex] = range;
   assert(apbound[nodeVindex] > 0);

   /* update maximum saturation degree: maxsatdeg = max { satdeg(v_i) + weight(v_i) | v_i in V } */
   maxsatdegree = range;

   debugMessage("-> updated neighbors:\n");

   /* set neighbor color of the adjacent nodes of node */
   currentadjedge = tcliqueGetFirstAdjedge(tcliquedata, node);
   lastadjedge = tcliqueGetLastAdjedge(tcliquedata, node);
   for( i = 0; i < nV; i++ )
   {
      assert(i == 0 || V[i] > V[i-1]);

      /* check if V[i] is contained in adjacency list of node started at position of V[i-1] 
       * (list is ordered by adjacent nodes)
       */
      for( ; currentadjedge <= lastadjedge; currentadjedge++ ) 
      {
         if( *currentadjedge >= V[i] )
         {
            if( *currentadjedge == V[i] )
            {
               assert(tcliqueIsEdge(tcliquedata, V[i], node));
               debugMessage("     nodeVindex=%d, node=%d, weight=%d, satdegold=%d  ->  ", 
                  i, V[i], weights[V[i]], gsd[i].satdeg); 
               
               /* sets satdeg for adjacent node */
               gsd[i].satdeg = range;
               
               /* creates new color interval [1,range] */
               ALLOC_ABORT( allocChunkMemory(mem, &colorinterval) );
               colorinterval->next = NULL;
               colorinterval->itv.inf = 1;
               colorinterval->itv.sup = range;
               
               /* colorinterval is the first added element of the list of neighborcolors of the adjacent node  */ 
               gsd[i].lcitv = colorinterval;

               debug(printf("satdegnew=%d, nbc=[%d,%d]\n", gsd[i].satdeg, gsd[i].lcitv->itv.inf, gsd[i].lcitv->itv.sup));
            }
            break;
         }
      }
   }

   /* set up data structures for the current clique */
   ALLOC_ABORT( allocMemoryArray(&currentclique, nV) );
   ncurrentclique = 0;
   weightcurrentclique = 0;
   workclique = clique;

   /* add node to the current clique */ 
   currentclique[0] = node; 
   ncurrentclique = 1; 
   weightcurrentclique = range; 
      
   /* color all other nodes of V */
   for( i = 0 ; i < nV-1; i++ )
   {
      assert((workclique == clique) != (currentclique == clique));

      /* selects the next uncolored node to color */
      nodeVindex = getMaxSatdegIndex(V, nV, gsd, iscolored, weights);
      if( nodeVindex == -1 ) /* no nodes left with positive weight */
         break;

      node = V[nodeVindex];
      assert(0 <= node && node < tcliqueGetNNodes(tcliquedata));
      range = weights[node];
      assert(range > 0);
      iscolored[nodeVindex] = YES;	

      debugMessage("%d. node choosen: vindex=%d, vertex=%d, satdeg=%d, range=%d, growclique=%d, weight=%d)\n",
         i+2, nodeVindex, node, gsd[nodeVindex].satdeg, range, growclique, weightcurrentclique);

      /* set apriori bound: apbound(v_i) = satdeg(v_i) + weight(v_i) */
      apbound[nodeVindex] = gsd[nodeVindex].satdeg + range;
      assert(apbound[nodeVindex] > 0);

      /* update maximum saturation degree: maxsatdeg = max { satdeg(v_i) + weight(v_i) | v_i in V } */
      if( maxsatdegree < apbound[nodeVindex] )
         maxsatdegree = apbound[nodeVindex];
      
      /* update clique */
      if( gsd[nodeVindex].satdeg == 0 )
      {
         /* current node is not adjacent to nodes of current clique, 
          * i.e. current clique can not be increased
          */
         debugMessage("current node not adjacend to current clique (weight:%d) -> starting new clique\n", 
            weightcurrentclique);

         /* check, if weight of current clique is larger than weight of maximum weight clique found so far */ 
         if( weightcurrentclique > *weightclique )
         {
            int* tmp;

            /* update maximum weight clique found so far */
            assert((workclique == clique) != (currentclique == clique));
            tmp = workclique;
            *weightclique = weightcurrentclique;
            *nclique = ncurrentclique;
            workclique = currentclique;
            currentclique = tmp;
            assert((workclique == clique) != (currentclique == clique));
         }
         weightcurrentclique = 0;
         ncurrentclique = 0;
         growclique = YES;
      }
      if( growclique )
      {
         /* check, if the current node is still adjacent to all nodes in the clique */
         if( gsd[nodeVindex].satdeg == weightcurrentclique )
         {
            assert(ncurrentclique < nV);
            currentclique[ncurrentclique] = node;
            ncurrentclique++; 
            weightcurrentclique += range;
#ifdef DEBUG
            {
               int k;
               printf("current clique (size:%d, weight:%d):", ncurrentclique, weightcurrentclique);
               for( k = 0; k < ncurrentclique; ++k )
                  printf(" %d", currentclique[k]);
               printf("\n");
            }
#endif
         }
         else
         {
            debugMessage("node satdeg: %d, clique weight: %d -> stop growing clique\n", 
               gsd[nodeVindex].satdeg, weightcurrentclique);
            growclique = NO;
         }
      }

      /* search for fitting color intervals for current node */
      pnc = &nwcitv;
      if( gsd[nodeVindex].lcitv == NULL )
      {
         /* current node has no colored neighbors yet: create new color interval [1,range] */
         ALLOC_ABORT( allocChunkMemory(mem, &colorinterval) );
         colorinterval->next = NULL;
         colorinterval->itv.inf = 1;
         colorinterval->itv.sup = range;
         
         /* add the new colorinterval [1, range] to the list of chosen colorintervals for node */
         pnc->next = colorinterval;
         pnc = colorinterval;
      }
      else
      {
         int tocolor;
         int dif;
         
         /* current node has colored neighbors */
         tocolor = range;
         lcitv = gsd[nodeVindex].lcitv;
         
         /* check, if first neighbor color interval [inf, sup] has inf > 1 */
         if( lcitv->itv.inf != 1 )
         {
            /* create new interval [1, min{range, inf}] */ 
            dif =  lcitv->itv.inf - 1 ;
            if( dif > tocolor )
               dif = tocolor;
            
            ALLOC_ABORT( allocChunkMemory(mem, &colorinterval) );
            colorinterval->next = NULL;
            colorinterval->itv.inf = 1;
            colorinterval->itv.sup = dif;

            tocolor -= dif;
            pnc->next = colorinterval;
            pnc = colorinterval;
         }

         /* as long as node is not colored with all colors, create new color interval by filling 
          * the gaps in the existing neighbor color intervals of the neighbors of node
          */
         while( tocolor > 0 )
         {	
            dif = tocolor;	
            
            ALLOC_ABORT( allocChunkMemory(mem, &colorinterval) );
            colorinterval->next = NULL;
            colorinterval->itv.inf = lcitv->itv.sup+1;			
            if( lcitv->next != NULL )
            {
               int min;

               min = lcitv->next->itv.inf - lcitv->itv.sup - 1;
          
               if( dif > min )  
                  dif = min;	
               lcitv = lcitv->next;
            }
            colorinterval->itv.sup = colorinterval->itv.inf + dif - 1;
            
            tocolor -= dif;
            pnc->next = colorinterval;
            pnc = colorinterval;
         }	
      }
      
      debugMessage("-> updated neighbors:\n"); 

      /* update saturation degree and neighbor colorintervals of all neighbors of node */
      currentadjedge = tcliqueGetFirstAdjedge(tcliquedata, node);
      lastadjedge = tcliqueGetLastAdjedge(tcliquedata, node);
      for( j = 0; j < nV; j++)
      {
         /* update only uncolored neighbors */
         if( iscolored[j] )
            continue;
      
         /* check if V[i] is contained in adjacency list of node started at position of V[i-1] 
          * (list is ordered by adjacent nodes)
          */
         for( ; currentadjedge <= lastadjedge; currentadjedge++ ) 
         {
            if( *currentadjedge >= V[j] )
            {
               if( *currentadjedge == V[j] )
               {
                  assert(tcliqueIsEdge(tcliquedata, V[j], node));
                  
                  debugMessage("     nodeVindex=%d, node=%d, weight=%d, satdegold=%d  ->  ", 
                     j, V[j], weights[V[j]], gsd[j].satdeg); 
                  updateNeighbor(mem, &gsd[j], nwcitv.next);
                  debug(printf("satdegnew=%d, nbc=[%d,%d]\n", 
                        gsd[j].satdeg, gsd[j].lcitv->itv.inf, gsd[j].lcitv->itv.sup));
               }
               break;
            }
         }
      }

      /* free data structure of created colorintervals */
      item = nwcitv.next;
      while( item != NULL )
      {
         tmpitem = item->next;                  
         freeChunkMemory(mem, &item);       
         item = tmpitem;                        
      }                                     

      /* free data structure of neighbor colorinterval of node just colored */
      item = gsd[nodeVindex].lcitv;
      while( item != NULL )
      {
         tmpitem = item->next;                  
         freeChunkMemory(mem, &item);       
         item = tmpitem;                        
      }                                     
   }
   assert((workclique == clique) != (currentclique == clique));
   
   /* update maximum weight clique found so far */
   if( weightcurrentclique > *weightclique )
   {
      int* tmp;
    
      tmp = workclique;
      *weightclique = weightcurrentclique;
      *nclique = ncurrentclique;
      workclique = currentclique;
      currentclique = tmp;
   }
   assert((workclique == clique) != (currentclique == clique));

   /* move the found clique to the provided clique pointer, if it is not the memory array */
   if( workclique != clique )
   {
      assert(clique == currentclique);
      assert(*nclique <= nV);
      copyMemoryArray(clique, workclique, *nclique);
      currentclique = workclique;
   }

   /* free data structures */
   freeMemoryArray(&currentclique);

   /* clear chunk memory */
   clearChunkMemory(mem);

   debugMessage("------------coloringend-----------------\n");

   return maxsatdegree;
}
