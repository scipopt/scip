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
#pragma ident "@(#) $Id: tclique_coloring.c,v 1.2 2005/04/25 14:34:08 bzfpfend Exp $"

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



/** gets index of the uncolored node in a given array of nodes in V with maximum satdeg. 
 *  in case of a tie choose node with maximum weight. V has to have uncolored nodes.
 */
static
int getMaxSatdegIndex(
   int*             V,                  /**< nodes in V */ 
   int              nV,		        /**< number of nodes in V */
   NBC*             gsd,                /**< neighbour color information of all nodes */
   BOOL*            iscolored,          /**< coloring status of all nodes */
   WEIGHT*          weights             /**< weight of nodes in grpah */
   )
{   
   WEIGHT maxweight;
   int maxsatdeg;
   int maxsatdegindex;
   int satdeg;
   int i;
	
   maxweight = -1.;
   maxsatdeg = -1.;
   maxsatdegindex = -1;
   
   assert(gsd != NULL);
   assert(iscolored != NULL);

   for( i = 0; i < nV; i++ )
   {
      /* checks only uncolored nodes */ 
      if( iscolored[i] ) 
         continue;
      
      satdeg = gsd[i].satdeg;
      if( satdeg  >= maxsatdeg )
      {
         /* tie: satdeg(v_i) = maxsatdeg */
         if( satdeg == maxsatdeg )
         {
            /* chooses node with maximum weight */
            if( weights[V[i]] > maxweight )
            { 	 
               maxweight = weights[V[i]];
               maxsatdegindex = i;
            }
         }
         /* satdeg(v_i) > maxsatdeg */
         else
         {
            maxsatdeg = satdeg;
            maxweight = weights[V[i]];
            maxsatdegindex = i;
         }
      }
   }
   
   assert(maxsatdegindex > -1);

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
   maxweight = -1.;

   /* tries to improve maxweight */
   for( i = 0 ; i < nnodes; i++ )
   {
      /* node has larger weight */
      if( weights[nodes[i]] > maxweight)
      {
         maxweight = weights[nodes[i]];
         maxweightindex = i;
      }
   }
   assert(maxweightindex >= 0);

   return maxweightindex;
}

/** updates the neighbour colors information of a node: updates the list of neighbour color intervals 
 *  by making the union of the existing list and the given list of color intervals, and updates the saturation degree
 */
static
void updateNeighbour(
   CHKMEM*          mem,                /**< block memory */
   NBC*             pgsd,               /**< pointer to neighbour color information of node to update */
   LIST_ITV*        pnc                 /**< pointer to given list of color intervals */
   )
{
   LIST_ITV head;
   LIST_ITV* apciv;
   LIST_ITV* pciv;
   LIST_ITV* nciv;

   /* saves the pointer to the first element of the list */
   head.next = pgsd->lcitv;
   apciv = &head;
   pciv = apciv->next;
   
   /* makes the union of the two intervals */
   while( (pnc != NULL) && (pciv != NULL) )
   {
      if( pnc->itv.inf < pciv->itv.inf ) 
      {	
         ALLOC_ABORT( allocChunkMemory(mem, &nciv) );
         nciv->next = NULL;
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
      nciv->next = NULL;
      nciv->itv = pnc->itv;
      nciv->next = NULL;
      
      apciv->next = nciv;
      apciv = nciv;
      
      pnc = pnc->next ;
   }

   /* tries to reduce the number of intervals */
   pgsd->satdeg = 0;
   apciv = head.next;
   while( (pciv = apciv->next) != NULL )
   {
      
      if( apciv->itv.sup < (pciv->itv.inf -1) )
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
       
         /* frees data structure for created colorinterval */
         tmp = pciv->next; 
         freeChunkMemory(mem, &pciv); 
         pciv = tmp; 
    
      }
   }
   pgsd->satdeg += apciv->itv.sup - apciv->itv.inf + 1;
   
   /* updates the pointer to the first element of the list */		
   pgsd->lcitv = head.next;
}

/** colors the nodes of a given set of nodes V with the lowest possible color and 
 *  finds a clique in the graph induced by V, an upper bound and an apriori bound for further branching steps
 */
int tcliqueColoring( 
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   CHKMEM*          mem,                /**< block memory */
   int*             V,                  /**< nodes for branching */ 
   int              nV,		        /**< number of nodes for branching */
   NBC*             gsd,                /**< neighbour color information of all nodes */
   BOOL*            iscolored,          /**< coloring status of all nodes */
   WEIGHT*          apbound,            /**< pointer to store apriori bound of nodes for branching */ 
   int**            clique,             /**< pointer to store the clique */           
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
   int* currentclique;
   int ncurrentclique;
   int weightcurrentclique;
   int* currentadjedge;
   int* lastadjedge;

   assert(clique != NULL);
   assert(*clique != NULL);
   assert(nclique != NULL);
   assert(weightclique != NULL);
   assert(gsd != NULL);
   assert(iscolored != NULL);

   weights = tcliqueGetWeights(tcliquedata);

   /* sets data structures for coloring */
   clearMemoryArray(iscolored, nV); /* new-memory */
   clearMemoryArray(gsd, nV); /* new-memory */
      
   /* initialise maximum weight clique found so far */
   growclique = YES;
   *nclique = 0;
   *weightclique = 0.;

   /* sets data structures for the current clique */
   ALLOC_ABORT( allocMemoryArray(&currentclique, nV) );

   ncurrentclique = 0;
   weightcurrentclique = 0.;
   
   /* -> colors the first node: */

   /* gets node of V with maximum weight */
   nodeVindex = getMaxWeightIndex(tcliquedata, V, nV);
   node = V[nodeVindex];
   iscolored[nodeVindex] = YES;
   range = weights[node];

   debugMessage("---------------coloring-----------------\n\n1. node choosen: vindex=%d vertex=%d)\n ",nodeVindex, node);

   /* sets apriori bound 
    *   apbound(v_i) = satdeg(v_i) + weight(v_i) */
   apbound[nodeVindex] = range;

   /* updates maximum saturation degree 
    *   maxsatdeg = max { satdeg(v_i) + weight(v_i) | v_i in V } */
   maxsatdegree = range;

   debugMessage("-> updated neighbours:\n");

   /* sets neighbourcolor of the adjacent nodes of node, 
    * if node has weight zero the colorinterval is empty and update of neighbours is not nes. */     
   if( range > 0 )
   {
      currentadjedge = tcliqueGetFirstAdjedge(tcliquedata, node);
      lastadjedge = tcliqueGetLastAdjedge(tcliquedata, node);
      for( i = 0; i < nV; i++ )
      {
         /* checks if V[i] is contained in adjacency list of node started at position of V[i-1] 
          * (list is ordered by adjacent nodes) */
         for( ; currentadjedge <= lastadjedge; currentadjedge++ ) 
         {
            if( *currentadjedge >= V[i] )
            {
               if( *currentadjedge == V[i] )
               {
                  assert(tcliqueIsEdge(tcliquedata, V[i], node));
                  debugMessage("       nodeVindex=%d, node=%d, satdegold=%d  ->  ", i, V[i], gsd[i].satdeg); 
               
                  /* sets satdeg for adjacent node */
                  gsd[i].satdeg = range;
               
                  /* creates new color interval [1,range] */
                  ALLOC_ABORT( allocChunkMemory(mem, &colorinterval) );
                  colorinterval->next = NULL;
                  colorinterval->itv.inf = 1;
                  colorinterval->itv.sup = range;
               
                  /* colorinterval is the first added element of the list of neighbourcolors of the adjacent node  */ 
                  gsd[i].lcitv = colorinterval;

#ifdef DEBUG
                  printf("satdegnew=%d, nbc=[%d,%d]\n", gsd[i].satdeg, gsd[i].lcitv->itv.inf, gsd[i].lcitv->itv.sup);
#endif
               }
               break;
            }
         }
      }
   }

   /* adds node to the current clique */ 
   currentclique[0] = node; 
   ncurrentclique = 1; 
   weightcurrentclique = range; 
      
   /* -> colors all other nodes of V: */
   for( i = 0 ; i < nV-1; i++ )
   {
      /* selects the next uncolored node to color */
      nodeVindex = getMaxSatdegIndex(V, nV, gsd, iscolored, weights);
      node = V[nodeVindex];
      iscolored[nodeVindex] = YES;	
      range = weights[node];

      debugMessage("\nn. node choosen: vindex=%d vertex=%d)\n",nodeVindex, node);

      /* sets apriori bound 
       *   apbound(v_i) = satdeg(v_i) + weight(v_i) */
      apbound[nodeVindex] = gsd[nodeVindex].satdeg + range;
      
      /* updates maximum saturation degree 
       *   maxsatdeg = max { satdeg(v_i) + weight(v_i) | v_i in V } */
      if( maxsatdegree < apbound[nodeVindex] )
         maxsatdegree = apbound[nodeVindex];
      
      /* only nodes with zero weight, that are not adjacent to current clique left 
       * (no coloring steps nessesary, clique with greater weight can't be found */ 
      if( range == 0 && gsd[nodeVindex].satdeg == 0 )
         break;

      /* updates clique: */
      /* current node is not adjacent to nodes of current clique, 
       * i.e. current clique can not be increased */
      if( gsd[nodeVindex].satdeg == 0 )
      {
         /* weight of current clique is larger than weight of maximum weight clique found so far */ 
         if( weightcurrentclique > *weightclique )
         {
            int* tmp;

            /* updates maximum weight clique found so far */
            *weightclique = weightcurrentclique;
            *nclique = ncurrentclique;
            tmp = currentclique;
            currentclique = *clique;
            *clique = tmp;
         }
         weightcurrentclique = 0;
         ncurrentclique = 0;
         growclique = YES;
      }
      if( growclique )
      {
         if( gsd[nodeVindex].satdeg == weightcurrentclique )
         {
            currentclique[ncurrentclique] = node;
            ncurrentclique++; 
            weightcurrentclique += range;
         }
         else
         {
            growclique = NO;
         }
      }

      /* searches for fitting color intervals for current node: */
      /* current node has no colored neighbours yet */
      pnc = &nwcitv;
      if( gsd[nodeVindex].lcitv == NULL )
      {
         /* creates new color interval [1,range] */
         ALLOC_ABORT( allocChunkMemory(mem, &colorinterval) );
         colorinterval->next = NULL;
         colorinterval->itv.inf = 1;
         colorinterval->itv.sup = range;
         
         /* adds the new colorinterval [1, range] to the list of chosen colorintervals for node */
         pnc->next = colorinterval;
         pnc = colorinterval;
      }
      /* current node has colored neighbours */
      else
      {
         int tocolor;
         int dif;
         
         tocolor = range;
         lcitv = gsd[nodeVindex].lcitv;
         
         /* first neighbour color interval [inf, sup] has inf > 1 */
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
          * the gaps in the existing neighbour color intervals of the neighbours of node */
         while( tocolor > 0 )
         {	
            dif = tocolor;	
            
            ALLOC_ABORT( allocChunkMemory(mem, &colorinterval) );
            colorinterval->next = NULL;
            colorinterval->itv.inf = lcitv->itv.sup+1;			
            if( lcitv->next != NULL )
            {
               int min;

               min = lcitv->next->itv.inf - lcitv->itv.sup-1;
          
               if( dif > min )  
                  dif = min;	
               lcitv = lcitv->next;
            }
            colorinterval->itv.sup = colorinterval->itv.inf + dif -1;
            
            tocolor -= dif;
            pnc->next = colorinterval;
            pnc = colorinterval;
         }	
      }
      
      debugMessage("-> updated neighbours:\n"); 

      /* updates saturation degree and neighbour colorintervals of all neighbours of node; 
       * if node has weight zero the colorinterval is empty and update of neighbours is not nes. */     
      if( range > 0.0 )
      {
         currentadjedge = tcliqueGetFirstAdjedge(tcliquedata, node);
         lastadjedge = tcliqueGetLastAdjedge(tcliquedata, node);
         for( j = 0; j < nV; j++)
         {
            /* updates only uncolored neighbours */
            if( iscolored[j] )
               continue;
      
            /* checks if V[i] is contained in adjacency list of node started at position of V[i-1] 
             * (list is ordered by adjacent nodes) */
            for( ; currentadjedge <= lastadjedge; currentadjedge++ ) 
            {
               if( *currentadjedge >= V[j] )
               {
                  if( *currentadjedge == V[j] )
                  {
                     assert(tcliqueIsEdge(tcliquedata, V[j], node));
                  
                     debugMessage("       nodeVindex=%d, node=%d, satdegold=%d  ->  ", j, V[j], gsd[j].satdeg); 
                  
                     updateNeighbour(mem, &gsd[j], nwcitv.next);
#ifdef DEBUG
                     printf("satdegnew=%d, nbc=[%d,%d]\n", gsd[j].satdeg, gsd[j].lcitv->itv.inf, gsd[j].lcitv->itv.sup); 
#endif   
                  }
                  break;
               }
            }
         }

         /* frees data structure of created colorintervals */
         item = nwcitv.next;
         while( item != NULL )
         {
            tmpitem = item->next;                  
            freeChunkMemory(mem, &item);       
            item = tmpitem;                        
         }                                     
      }

      /* frees data structure of neighbour colorinterval of node just colored; even if weight is zero */
      item = gsd[nodeVindex].lcitv;
      while( item != NULL )
      {
         tmpitem = item->next;                  
         freeChunkMemory(mem, &item);       
         item = tmpitem;                        
      }                                     
      
   }
   
   /* updates maximum weight clique found so far */
   if( weightcurrentclique > *weightclique )
   {
      int* tmp;
    
      *weightclique = weightcurrentclique;
      *nclique = ncurrentclique;
      tmp = currentclique;
      currentclique = *clique;
      *clique = tmp;
   }

   /* frees data structures */
   freeMemory(&currentclique);

   debugMessage("------------coloringend-----------------\n");

   return maxsatdegree;
}
