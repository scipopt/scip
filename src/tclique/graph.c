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
#pragma ident "@(#) $Id: graph.c,v 1.1 2005/03/10 17:11:17 bzfpfend Exp $"

/**@file   graph.c
 * @brief  tclique data part of algorithm for maximum cliques
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>

#include "scip/def.h"
#include "scip/memory.h"
#include "tclique/graph.h" 



#define EPSILON (0.000999999)       


/** creates tclique data structure */
void tcliqueCreate(
   TCLIQUEDATA**    tcliquedata         /**< pointer to store tclique data structure */
   )
{
   assert(tcliquedata != NULL);

   ALLOC_ABORT( allocMemory(tcliquedata) );

   (*tcliquedata)->nnodes = 0;
   (*tcliquedata)->nedges = 0;
   (*tcliquedata)->weights = NULL;
   (*tcliquedata)->degrees = NULL;
   (*tcliquedata)->adjnodes = NULL;
   (*tcliquedata)->adjedges = NULL;
   (*tcliquedata)->sizenodes = 0;
   (*tcliquedata)->sizeedges = 0;
}

/** frees tclique data structure */
void tcliqueFree(
   TCLIQUEDATA**    tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);

   if( *tcliquedata != NULL )
   {
      freeMemoryArray(&(*tcliquedata)->adjedges);
      freeMemoryArray(&(*tcliquedata)->adjnodes);
      freeMemoryArray(&(*tcliquedata)->degrees);
      freeMemoryArray(&(*tcliquedata)->weights);
      freeMemory(tcliquedata);
   }
}

/** ensures, that arrays concerning edges in tclique data structure can store at least num entries */
static
void tcliqueEnsureSizeEdges(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to tclique data structure */
   int              num                 /**< minimum number of entries concerning edges to store */
   )
{
   assert(tcliquedata != NULL);
   
   /* create tclique data structure, if not yet existing */
   if( *tcliquedata == NULL )
   {
      tcliqueCreate(tcliquedata);
   }
   assert(*tcliquedata != NULL);

   if( num > (*tcliquedata)->sizeedges )
   {
      int newsize;
      newsize = 2*(*tcliquedata)->sizeedges;
      if( newsize < num )
         newsize = num;

      ALLOC_ABORT( reallocMemoryArray(&(*tcliquedata)->adjnodes, newsize) );

      (*tcliquedata)->sizeedges = newsize;
   }

   assert(num <= (*tcliquedata)->sizeedges);
}

/** ensures, that arrays concerning nodes in tclique data structure can store at least num entries */
static
void tcliqueEnsureSizeNodes(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to tclique data structure */
   int              num                 /**< minimum number of entries concerning nodes to store */
   )
{
   assert(tcliquedata != NULL);
   
   /* create tclique data structure, if not yet existing */
   if( *tcliquedata == NULL )
   {
      tcliqueCreate(tcliquedata);
   }  
   assert(*tcliquedata != NULL);

   if( (*tcliquedata)->adjnodes == NULL ) 
   {
      tcliqueEnsureSizeEdges(tcliquedata, 1);
   }
   assert((*tcliquedata)->adjnodes != NULL);
   
   if( num > (*tcliquedata)->sizenodes )
   {
      int newsize;
      int i;

      newsize = 2*(*tcliquedata)->sizenodes;
      if( newsize < num )
         newsize = num;

      ALLOC_ABORT( reallocMemoryArray(&(*tcliquedata)->weights, newsize) );
      ALLOC_ABORT( reallocMemoryArray(&(*tcliquedata)->degrees, newsize) );
      ALLOC_ABORT( reallocMemoryArray(&(*tcliquedata)->adjedges, newsize) );

      for( i = (*tcliquedata)->sizenodes; i < newsize; i++ )
      {
         (*tcliquedata)->weights[i] = 0;
         (*tcliquedata)->degrees[i] = 0;
         (*tcliquedata)->adjedges[i].first = (*tcliquedata)->nedges;
         (*tcliquedata)->adjedges[i].last = (*tcliquedata)->adjedges[i].first;
         assert((*tcliquedata)->adjedges[i].first >= 0);
      }

      (*tcliquedata)->sizenodes = newsize;
   }
   assert(num <= (*tcliquedata)->sizenodes);
}


/** adds node to tclique data structure */
void tcliqueAddNode(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to tclique data structure */
   int              node,               /**< node to add */
   WEIGHT           weight              /**< weight of node to add (allready scaled) */
   )
{
   assert(weight >= 0);

   tcliqueEnsureSizeNodes(tcliquedata, node + 1);

   (*tcliquedata)->weights[node] = weight;
   
   assert((*tcliquedata)->degrees[node] == 0);
   assert((*tcliquedata)->adjedges[node].first <= (*tcliquedata)->nedges);
   assert((*tcliquedata)->adjedges[node].last == (*tcliquedata)->adjedges[node].first);
   ((*tcliquedata)->nnodes)++;
}

/** adds edge (node1, node2) to tclique data structure (node1 and node2 have to be contained in 
 *  tclique data structure) */
void tcliqueAddEdge(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to tclique data structure */
   int              node1,              /**< start node of edge to add */
   int              node2               /**< end node of edge to add */
   )
{
   int i;
   int firstadjedge;
   int lastadjedge;
   int currentadjedge;
   
   assert(tcliquedata != NULL);
   assert(*tcliquedata != NULL);
   assert(node1 >= 0);
   assert(node2 >= 0);
   assert(node1 != node2);
   assert((*tcliquedata)->weights[node1] >= 0);
   assert((*tcliquedata)->weights[node2] >= 0);
   assert((*tcliquedata)->degrees[node1] >= 0);
   assert((*tcliquedata)->degrees[node2] >= 0);

   tcliqueEnsureSizeEdges(tcliquedata, ((*tcliquedata)->nedges) + 1);

   assert((*tcliquedata)->adjedges[node1].first >= 0);
   assert((*tcliquedata)->adjedges[node2].first >= 0);
   assert((*tcliquedata)->adjedges[node1].first <= (*tcliquedata)->nedges);
   assert((*tcliquedata)->adjedges[node2].first <= (*tcliquedata)->nedges);

   /* updates adjnodes and adjedges for nodes after node1 */
   for( i = ((*tcliquedata)->sizenodes) - 1; i > node1; i-- )
   {
      if( (*tcliquedata)->degrees[i] > 0 )
      {
         firstadjedge = (*tcliquedata)->adjedges[i].first;
         lastadjedge = (*tcliquedata)->adjedges[i].last;
         currentadjedge = lastadjedge;
   
         assert(firstadjedge + (*tcliquedata)->degrees[i] == lastadjedge);

         for( ; currentadjedge > firstadjedge; currentadjedge-- )
         {
            if( currentadjedge != lastadjedge )
            {
               assert((*tcliquedata)->adjnodes[currentadjedge] > (*tcliquedata)->adjnodes[currentadjedge - 1]);
            }            
            (*tcliquedata)->adjnodes[currentadjedge] = (*tcliquedata)->adjnodes[currentadjedge - 1];
         }
         assert(currentadjedge == firstadjedge);
      }
      ((*tcliquedata)->adjedges[i].first)++;
      ((*tcliquedata)->adjedges[i].last)++;
   }
   assert(i == node1);
   
   /* updates adjnodes and adjedges for node1 */
   firstadjedge = (*tcliquedata)->adjedges[node1].first;
   lastadjedge = (*tcliquedata)->adjedges[node1].last;
   currentadjedge = lastadjedge;

   if( (*tcliquedata)->degrees[node1] > 0 ) 
   {
      assert(firstadjedge + (*tcliquedata)->degrees[node1] == lastadjedge);
      
      /* inserts node2 into adjnodes at the right position (adjnodes are sorted by startnodes and endnode) */
      for( ; currentadjedge > firstadjedge; currentadjedge-- )
      {
         if( (*tcliquedata)->adjnodes[currentadjedge - 1] > node2 ) 
         {
            (*tcliquedata)->adjnodes[currentadjedge] = (*tcliquedata)->adjnodes[currentadjedge - 1];
         }
         else
         {
            assert((*tcliquedata)->adjnodes[currentadjedge - 1] < node2);
            (*tcliquedata)->adjnodes[currentadjedge] = node2;
            break;
         }
      }
   }

   /* inserts node2 as the first adjacent node of node1 */ 
   if( (*tcliquedata)->degrees[node1] == 0 || currentadjedge == firstadjedge )
   {
      if( (*tcliquedata)->degrees[node1] > 0 )
      {
         assert((*tcliquedata)->adjnodes[firstadjedge + 1] > node2);
      }
      (*tcliquedata)->adjnodes[firstadjedge] = node2;
   }
   
   ((*tcliquedata)->adjedges[node1].last)++;
   ((*tcliquedata)->degrees[node1])++;
   ((*tcliquedata)->nedges)++;
}

/** changes weight of node in tclique data structure */
void tcliqueNodeChangeWeight(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to tclique data structure */
   int              node,               /**< node to set new weight */
   WEIGHT           weight              /**< new weight of node (allready scaled) */
   )
{
   assert(node < getNnodes(*tcliquedata));
      
   (*tcliquedata)->weights[node] = weight;
}

/** loads tclique data structure */
BOOL tcliqueLoadComplete(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to store tclique data structure */
   char*            filename,           /**< name of file with graph data */
   char*            probname            /**< name of problem */
   )
{
   FILE* file;
   float weight;
   int node1;
   int node2;
   int currentnode;
   int i;
   
   assert(tcliquedata != NULL);

   /* opens file */
   if( (file = fopen(filename, "rt")) == NULL )
   {
      if( (file = fopen("default.dat", "rt")) == NULL )
      {
         printf("\nCan't open file: %s", filename);
         return NO;
      }
   }
 
   ALLOC_ABORT( allocMemory(tcliquedata) );
 
   /* sets name of problem, number of nodes and number of edges in graph */
   fscanf(file, "%s", probname );
   fscanf(file, "%d", &(*tcliquedata)->nnodes );
   fscanf(file, "%d", &(*tcliquedata)->nedges );
   
   /* sets data structures for tclique */
   ALLOC_ABORT( allocMemoryArray(&(*tcliquedata)->weights, (*tcliquedata)->nnodes) );
   ALLOC_ABORT( allocMemoryArray(&(*tcliquedata)->degrees, (*tcliquedata)->nnodes) );
   ALLOC_ABORT( allocMemoryArray(&(*tcliquedata)->adjnodes, (*tcliquedata)->nedges) );
   ALLOC_ABORT( allocMemoryArray(&(*tcliquedata)->adjedges, (*tcliquedata)->nnodes) );

   /* sets weights of all nodes (scaled!) */
   for( i = 0; i < (*tcliquedata)->nnodes; i++ )
   {
      fscanf(file, "%f", &weight );
      (*tcliquedata)->weights[i] = (weight + EPSILON) * 1000;
   }

   /* sets adjacent edges and degree of all nodes */
   currentnode = -1;
   for( i = 0; i < (*tcliquedata)->nedges; i++ )
   {
      /* reads edge (node1, node2) */
      fscanf(file, "%d%d", &node1, &node2 );
      
      /* (node1, node2) is the first adjacent edge of node1 */
      if( node1 != currentnode )
      {
         currentnode = node1;
         (*tcliquedata)->degrees[currentnode] = 0;
         (*tcliquedata)->adjedges[currentnode].first = i;
         (*tcliquedata)->adjedges[currentnode].last = (*tcliquedata)->adjedges[currentnode].first;
      }
      (*tcliquedata)->degrees[currentnode]++;
      (*tcliquedata)->adjnodes[i] = node2;
      (*tcliquedata)->adjedges[currentnode].last++;
   }

   /* closes file */
   fclose(file);
   
   return YES;
}

/** loads tclique data structure */
BOOL tcliqueLoadStepwise(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to store tclique data structure */
   char*            filename,           /**< name of file with graph data */
   char*            probname            /**< name of problem */
   )
{
   FILE* file;
   float weight;
   int node1;
   int node2;
   int i;
   int nnodes;
   int nedges;

   assert(tcliquedata != NULL);
   assert(*tcliquedata == NULL);

   /* opens file */
   if( (file = fopen(filename, "rt")) == NULL )
   {
      if( (file = fopen("default.dat", "rt")) == NULL )
      {
         printf("\nCan't open file: %s", filename);
         return NO;
      }
   }
 
   /* sets name of problem, number of nodes and number of edges in graph */
   fscanf(file, "%s", probname );
   fscanf(file, "%d", &nnodes );
   fscanf(file, "%d", &nedges );

   /* sets weights of all nodes (scaled!) */
   for( i = 0; i < nnodes; i++ )
   {
      fscanf(file, "%f", &weight );
      tcliqueAddNode(tcliquedata, i, (weight + EPSILON) * 1000);
   }

   /* sets adjacent edges and degree of all nodes */
   for( i = 0; i < nedges; i++ )
   {
      /* reads edge (node1, node2) */
      fscanf(file, "%d%d", &node1, &node2 );
      tcliqueAddEdge(tcliquedata, node1, node2);
   }
   assert(nnodes == getNnodes(*tcliquedata));
   assert(nedges == getNedges(*tcliquedata));

   /* closes file */
   fclose(file);
   
   return YES;
}

/** gets number of nodes in the graph */
int getNnodes(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);
   
   return tcliquedata->nnodes;
}

/** gets number of edges in the graph */
int getNedges(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);
   
   return tcliquedata->nedges;
}

/** gets weight of nodes in the graph */
WEIGHT* getWeights(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);
   
   return tcliquedata->weights;
}

/** gets degree of nodes in graph */
int* getDegrees(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);
   
   return tcliquedata->degrees;
}

/** gets adjacent nodes of edges in graph */
int* getAdjnodes(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);
   
   return tcliquedata->adjnodes;
}

/** gets pointer to first and one after last adjacent edge of nodes in graph */
HEAD_ADJ* getAdjedges(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);
   
   return tcliquedata->adjedges;
}

/** gets pointer to first adjacent edge of given node in graph */
int* getFirstAdjedge(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node                /**< given node */
   )
{
   HEAD_ADJ* adjedges;
   int* adjnodes;

   assert(tcliquedata != NULL);
   assert(node >= 0);

   adjedges = getAdjedges(tcliquedata);
   assert(adjedges != NULL);
   assert(adjedges[node].first >= 0);
   assert(adjedges[node].first <= getNedges(tcliquedata));

   adjnodes = getAdjnodes(tcliquedata);
   assert(adjnodes != NULL);

   return &adjnodes[adjedges[node].first];
}

/** gets pointer to last adjacent edge of given node in graph */
int* getLastAdjedge(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node                /**< given node */
   )
{
   HEAD_ADJ* adjedges;
   int* adjnodes;
   int* degrees;

   assert(tcliquedata != NULL);
   assert(node >= 0);

   adjedges = getAdjedges(tcliquedata);
   degrees = getDegrees(tcliquedata);
   assert(adjedges != NULL);
   assert(degrees[node] == 0 || adjedges[node].last-1 >= 0);
   assert(adjedges[node].last-1 <= getNedges(tcliquedata));

   assert(adjedges[node].last - adjedges[node].first == degrees[node]);

   adjnodes = getAdjnodes(tcliquedata);
   assert(adjnodes != NULL);

   return &adjnodes[adjedges[node].last-1];
}

/** returns, whether the edge (node1, node2) is in the graph */
int isEdge(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node1,              /**< start node of edge */
   int              node2               /**< end node of edge */
   )
{
   int* currentadjedge;
   int* lastadjedge;
   int tmp;

   assert(tcliquedata != NULL);

   if( node1 < node2 )
   {
      tmp = node1;
      node1 = node2;
      node2 = tmp;
   }

   currentadjedge = getFirstAdjedge(tcliquedata, node1);
   lastadjedge = getLastAdjedge(tcliquedata, node1);
   
   if( *lastadjedge < node2 )
      return 0;

   /* checks if node2 is contained in adjacency list of node1 
    * (list is ordered by adjacent nodes) */
   while( currentadjedge <= lastadjedge ) 
   {
      if( *currentadjedge >= node2 )
      {
         if( *currentadjedge == node2 )
            return 1;
         else 
            break;
      }
      currentadjedge++;
   }

   return 0;
}

/* selects all nodes from a given set of nodes which are adjacent to a given node
 * and returns the number of selected nodes */
int selectAdjnodes( 
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node,               /**< given node */
   int*             nodes,              /**< given set of nodes (ordered by node index) */ 
   int              nnodes,             /**< number of nodes in given set nodes */ 
   int*             adjnodes            /**< pointer to store adjacent nodes of node, contained in nodes */ 
   )
{
   int nadjnodes;
   int* currentadjedge;
   int* lastadjedge;
   int i;

   assert(tcliquedata != NULL);
   assert(adjnodes != NULL);

   nadjnodes = 0;
   currentadjedge = getFirstAdjedge(tcliquedata, node);
   lastadjedge = getLastAdjedge(tcliquedata, node);

   /* checks for each node in given set nodes, if it is adjacent to given node 
    * (adjacent nodes are ordered by node index) */
   for( i = 0; i < nnodes; i++ )
   {
      for( ; currentadjedge <= lastadjedge; currentadjedge++ )
      {
         if( *currentadjedge >= nodes[i] )
         {
            /* current node is adjacent to given node */
            if( *currentadjedge == nodes[i] )
            {
               adjnodes[nadjnodes] = nodes[i]; 
               nadjnodes++;
            }
            break;
         } 
      }
   }
   
   return nadjnodes;
}

/* prints tclique data structure */
void printTcliquedata(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   int* weights;
   int* adjnodes;
   int* degrees;
   int i;

   assert(tcliquedata != NULL);

   degrees = getDegrees(tcliquedata);
   weights = getWeights(tcliquedata);
   adjnodes = getAdjnodes(tcliquedata);

   printf("tcliquedata:\n");
   printf("nnodes=%d, nedges=%d\n", getNnodes(tcliquedata), getNedges(tcliquedata));
   for( i = 0; i < getNnodes(tcliquedata); i++ )
   {
      int* currentadjedge;
      int* lastadjedge;

      printf("node %d: weight=%d, degree=%d, adjnodes=\n[ ", i, weights[i], degrees[i]);  

      currentadjedge = getFirstAdjedge(tcliquedata, i);
      lastadjedge = getLastAdjedge(tcliquedata, i);
      assert(lastadjedge + 1 - currentadjedge == degrees[i]);

      for( ; currentadjedge <= lastadjedge; currentadjedge++ )
      {
         printf("%d, ", *currentadjedge);
      }
      printf("]\n");
   }
}
