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
#pragma ident "@(#) $Id: tclique_graph.c,v 1.6 2005/05/31 17:20:24 bzfpfend Exp $"

/**@file   tclique_graph.c
 * @brief  tclique data part of algorithm for maximum cliques
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>

#include "scip/def.h"
#include "scip/memory.h"
#include "scip/tclique_graph.h" 



#define ALLOC_FALSE(x)  do                                                 \
   {                                                                    \
      if( NULL == (x) )                                                 \
      {                                                                 \
         printf("[%s:%d] ERROR: No memory in function call\n", __FILE__, __LINE__); \
         fflush(stdout);                                                \
         return FALSE;                                                     \
      }                                                                 \
   }                                                                    \
   while( FALSE )


struct _TCLIQUEDATA
{
   int              nnodes;		/**< number of nodes in graph */
   int              nedges;		/**< number of edges in graph */
   WEIGHT*          weights;	        /**< weight of nodes */
   int*             degrees;	        /**< degree of nodes */
   int*             adjnodes;	        /**< adjacent nodes of edges */
   HEAD_ADJ*        adjedges;           /**< pointer to first and one after last adjacent edge of nodes */
   int              sizenodes;		/**< size of arrays concerning nodes (weights, degrees and adjedges) */
   int              sizeedges;		/**< size of arrays concerning edges (adjnodes) */
   int*             cacheddegrees;      /**< number of adjacent cached edges for each node */
   int*             cachedorigs;        /**< origin nodes of cached edges */
   int*             cacheddests;        /**< destination nodes of cached edges */
   int              ncachededges;       /**< number of cached edges (not yet inserted in all data structures) */
   int              sizecachededges;    /**< size of arrays concerning cached edges */
}; 


/** creates tclique data structure */
Bool tcliqueCreate(
   TCLIQUEDATA**    tcliquedata         /**< pointer to store tclique data structure */
   )
{
   assert(tcliquedata != NULL);

   ALLOC_FALSE( allocMemory(tcliquedata) );

   (*tcliquedata)->nnodes = 0;
   (*tcliquedata)->nedges = 0;
   (*tcliquedata)->weights = NULL;
   (*tcliquedata)->degrees = NULL;
   (*tcliquedata)->adjnodes = NULL;
   (*tcliquedata)->adjedges = NULL;
   (*tcliquedata)->sizenodes = 0;
   (*tcliquedata)->sizeedges = 0;
   (*tcliquedata)->cacheddegrees = NULL;
   (*tcliquedata)->cachedorigs = NULL;
   (*tcliquedata)->cacheddests = NULL;
   (*tcliquedata)->ncachededges = 0;
   (*tcliquedata)->sizecachededges = 0;

   return TRUE;
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
      freeMemoryArrayNull(&(*tcliquedata)->cacheddegrees);
      freeMemoryArrayNull(&(*tcliquedata)->cachedorigs);
      freeMemoryArrayNull(&(*tcliquedata)->cacheddests);
      freeMemory(tcliquedata);
   }
}

/** ensures, that arrays concerning edges in tclique data structure can store at least num entries */
static
Bool tcliqueEnsureSizeEdges(
   TCLIQUEDATA*     tcliquedata,        /**< tclique data structure */
   int              num                 /**< minimum number of entries concerning edges to store */
   )
{
   assert(tcliquedata != NULL);
   
   if( num > tcliquedata->sizeedges )
   {
      int newsize;

      newsize = 2*tcliquedata->sizeedges;
      if( newsize < num )
         newsize = num;

      ALLOC_FALSE( reallocMemoryArray(&tcliquedata->adjnodes, newsize) );
      tcliquedata->sizeedges = newsize;
   }

   assert(num <= tcliquedata->sizeedges);

   return TRUE;
}

/** ensures, that arrays concerning cached edges in tclique data structure can store at least num entries */
static
Bool tcliqueEnsureSizeCachedEdges(
   TCLIQUEDATA*     tcliquedata,        /**< tclique data structure */
   int              num                 /**< minimum number of entries concerning cached edges to store */
   )
{
   assert(tcliquedata != NULL);
   
   if( num > tcliquedata->sizecachededges )
   {
      int newsize;

      newsize = 2*tcliquedata->sizecachededges;
      if( newsize < num )
         newsize = num;

      ALLOC_FALSE( reallocMemoryArray(&tcliquedata->cachedorigs, newsize) );
      ALLOC_FALSE( reallocMemoryArray(&tcliquedata->cacheddests, newsize) );
      tcliquedata->sizecachededges = newsize;
   }

   assert(num <= tcliquedata->sizecachededges);

   return TRUE;
}

/** ensures, that arrays concerning nodes in tclique data structure can store at least num entries */
static
Bool tcliqueEnsureSizeNodes(
   TCLIQUEDATA*     tcliquedata,        /**< tclique data structure */
   int              num                 /**< minimum number of entries concerning nodes to store */
   )
{
   assert(tcliquedata != NULL);
   
   if( !tcliqueEnsureSizeEdges(tcliquedata, 1) )
      return FALSE;
   assert(tcliquedata->adjnodes != NULL);
   
   if( num > tcliquedata->sizenodes )
   {
      int newsize;
      int i;

      newsize = 2*tcliquedata->sizenodes;
      if( newsize < num )
         newsize = num;

      ALLOC_FALSE( reallocMemoryArray(&tcliquedata->weights, newsize) );
      ALLOC_FALSE( reallocMemoryArray(&tcliquedata->degrees, newsize) );
      ALLOC_FALSE( reallocMemoryArray(&tcliquedata->adjedges, newsize) );

      for( i = tcliquedata->sizenodes; i < newsize; i++ )
      {
         tcliquedata->weights[i] = 0;
         tcliquedata->degrees[i] = 0;
         tcliquedata->adjedges[i].first = tcliquedata->nedges;
         tcliquedata->adjedges[i].last = tcliquedata->nedges;
      }

      if( tcliquedata->ncachededges > 0 )
      {
         assert(tcliquedata->cacheddegrees != NULL);
         ALLOC_FALSE( reallocMemoryArray(&tcliquedata->cacheddegrees, newsize) );
         for( i = tcliquedata->sizenodes; i < newsize; i++ )
            tcliquedata->cacheddegrees[i] = 0;
      }

      tcliquedata->sizenodes = newsize;
   }
   assert(num <= tcliquedata->sizenodes);

   return TRUE;
}


/** adds nodes up to the given node number to tclique data structure (intermediate nodes have weight 0) */
Bool tcliqueAddNode(
   TCLIQUEDATA*     tcliquedata,        /**< tclique data structure */
   int              node,               /**< node number to add */
   WEIGHT           weight              /**< weight of node to add */
   )
{
   assert(weight >= 0);

   if( !tcliqueEnsureSizeNodes(tcliquedata, node + 1) )
      return FALSE;

   tcliquedata->weights[node] = weight;
   
   assert(tcliquedata->degrees[node] == 0);
   assert(tcliquedata->adjedges[node].first <= tcliquedata->nedges);
   assert(tcliquedata->adjedges[node].last == tcliquedata->adjedges[node].first);
   tcliquedata->nnodes = MAX(tcliquedata->nnodes, node+1);

   return TRUE;
}

/** changes weight of node in tclique data structure */
void tcliqueChangeWeight(
   TCLIQUEDATA*     tcliquedata,        /**< tclique data structure */
   int              node,               /**< node to set new weight */
   WEIGHT           weight              /**< new weight of node (allready scaled) */
   )
{
   assert(0 <= node && node < tcliqueGetNNodes(tcliquedata));
   assert(weight >= 0);

   tcliquedata->weights[node] = weight;
}

/** adds edge (node1, node2) to tclique data structure (node1 and node2 have to be contained in 
 *  tclique data structure);
 *  new edges are cached, s.t. the graph data structures are not correct until a call to tcliqueFlush();
 *  you have to make sure, that no double edges are inserted
 */
Bool tcliqueAddEdge(
   TCLIQUEDATA*     tcliquedata,        /**< tclique data structure */
   int              node1,              /**< start node of edge to add */
   int              node2               /**< end node of edge to add */
   )
{
   assert(tcliquedata != NULL);
   assert(0 <= node1 && node1 < tcliquedata->nnodes);
   assert(0 <= node2 && node2 < tcliquedata->nnodes);
   assert(node1 != node2);

   if( !tcliqueEnsureSizeCachedEdges(tcliquedata, tcliquedata->ncachededges + 2) )
      return FALSE;

   /* make sure, the array for counting the cached node degrees exists */
   if( tcliquedata->ncachededges == 0 && tcliquedata->sizenodes > 0 )
   {
      assert(tcliquedata->cacheddegrees == NULL);
      ALLOC_FALSE( allocMemoryArray(&tcliquedata->cacheddegrees, tcliquedata->sizenodes) );
      clearMemoryArray(tcliquedata->cacheddegrees, tcliquedata->sizenodes);
   }
   assert(tcliquedata->cacheddegrees != NULL);

   /* just remember both new half edges in the cache; the full insertion is done later on demand */
   tcliquedata->cachedorigs[tcliquedata->ncachededges] = node1;
   tcliquedata->cacheddests[tcliquedata->ncachededges] = node2;
   tcliquedata->ncachededges++;
   tcliquedata->cachedorigs[tcliquedata->ncachededges] = node2;
   tcliquedata->cacheddests[tcliquedata->ncachededges] = node1;
   tcliquedata->ncachededges++;
   tcliquedata->cacheddegrees[node1]++;
   tcliquedata->cacheddegrees[node2]++;

   return TRUE;
}

/** inserts all cached edges into the data structures */
Bool tcliqueFlush(
   TCLIQUEDATA*     tcliquedata         /**< tclique data structure */
   )
{
   assert(tcliquedata != NULL);

   /* check, whether there are cached edges */
   if( tcliquedata->ncachededges > 0 )
   {
      int ninsertedholes;
      int pos;
      int n;
      int i;

      /* reallocate adjnodes array to be able to store all additional edges */
      if( !tcliqueEnsureSizeEdges(tcliquedata, tcliquedata->nedges + tcliquedata->ncachededges) )
         return FALSE;
      assert(tcliquedata->adjnodes != NULL);
      assert(tcliquedata->adjedges != NULL);

      /* move the old edges in the adjnodes array, s.t. there is enough free space for the additional edges */
      ninsertedholes = 0;
      pos = tcliquedata->nedges + tcliquedata->ncachededges - 1;
      for( n = tcliquedata->nnodes-1; ; --n ) /* no abort criterion, because at n == 0, the loop is break'ed */
      {
         int olddegree;

         assert(n >= 0);
         assert(tcliquedata->adjedges[n].last - tcliquedata->adjedges[n].first == tcliquedata->degrees[n]);

         /* increase the degree of the node */
         olddegree = tcliquedata->degrees[n];
         tcliquedata->degrees[n] += tcliquedata->cacheddegrees[n];

         /* skip space for new edges */
         pos -= tcliquedata->cacheddegrees[n];
         ninsertedholes += tcliquedata->cacheddegrees[n];
         assert(ninsertedholes <= tcliquedata->ncachededges);
         if( ninsertedholes == tcliquedata->ncachededges )
            break;
         assert(n > 0);

         /* move old edges */
         for( i = tcliquedata->adjedges[n].last - 1; i >= tcliquedata->adjedges[n].first; --i, --pos )
         {
            assert(0 <= i && i < pos && pos < tcliquedata->nedges + tcliquedata->ncachededges);
            tcliquedata->adjedges[pos] = tcliquedata->adjedges[i];
         }

         /* adjust the first and last edge pointers of the node */
         tcliquedata->adjedges[n].first = pos+1;
         tcliquedata->adjedges[n].last = pos+1 + olddegree;

         assert(n == tcliquedata->nnodes-1
            || tcliquedata->adjedges[n].first + tcliquedata->degrees[n] == tcliquedata->adjedges[n+1].first);
      }
      assert(ninsertedholes == tcliquedata->ncachededges - tcliquedata->nedges);
      assert(tcliquedata->adjedges[n].last == pos+1);
#ifndef NDEBUG
      for( --n; n >= 0; --n )
         assert(tcliquedata->cacheddegrees[n] == 0);
#endif

      /* insert the cached edges into the adjnodes array */
      for( i = 0; i < tcliquedata->ncachededges; ++i )
      {
         int dest;

         n = tcliquedata->cachedorigs[i];
         dest = tcliquedata->cacheddests[i];
         assert(0 <= n && n < tcliquedata->nnodes);
         assert(0 <= dest && dest < tcliquedata->nnodes);
         assert(tcliquedata->adjedges[n].last <= tcliquedata->nedges + tcliquedata->ncachededges);
         assert(n == tcliquedata->nnodes-1 || tcliquedata->adjedges[n].last <= tcliquedata->adjedges[n+1].first);
         assert(n == tcliquedata->nnodes-1
            || tcliquedata->adjedges[n].first + tcliquedata->degrees[n] == tcliquedata->adjedges[n+1].first);

         /* edges of each node must be sorted by increasing destination node number */
         for( pos = tcliquedata->adjedges[n].last;
              pos > tcliquedata->adjedges[n].first && dest < tcliquedata->adjnodes[pos-1]; --pos )
         {
            tcliquedata->adjnodes[pos] = tcliquedata->adjnodes[pos-1];
         }
         tcliquedata->adjnodes[pos] = dest;
         tcliquedata->adjedges[n].last++;

         assert(n == tcliquedata->nnodes-1 || tcliquedata->adjedges[n].last <= tcliquedata->adjedges[n+1].first);
      }

      /* update the number of edges */
      tcliquedata->nedges += tcliquedata->ncachededges;

      /* free the cache */
      freeMemoryArray(&tcliquedata->cacheddegrees);
      freeMemoryArray(&tcliquedata->cachedorigs);
      freeMemoryArray(&tcliquedata->cacheddests);
      tcliquedata->ncachededges = 0;
      tcliquedata->sizecachededges = 0;
   }

   /* the cache should now be freed */
   assert(tcliquedata->ncachededges == 0);
   assert(tcliquedata->sizecachededges == 0);
   assert(tcliquedata->cacheddegrees == NULL);
   assert(tcliquedata->cachedorigs == NULL);
   assert(tcliquedata->cacheddests == NULL);

#ifndef NDEBUG
   /* check integrity of the data structures */
   {
      int pos;
      int n;

      pos = 0;
      for( n = 0; n < tcliquedata->nnodes; ++n )
      {
         int i;

         assert(tcliquedata->adjedges[n].first == pos);
         assert(tcliquedata->adjedges[n].last == tcliquedata->adjedges[n].first + tcliquedata->degrees[n]);

         for( i = tcliquedata->adjedges[n].first; i < tcliquedata->adjedges[n].last-1; ++i )
         {
            assert(tcliquedata->adjnodes[i] < tcliquedata->adjnodes[i+1]);
         }
         pos = tcliquedata->adjedges[n].last;
      }
      assert(pos == tcliquedata->nedges);
   }   
#endif

   return TRUE;
}

/** loads tclique data structure from file */
Bool tcliqueLoadFile(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to store tclique data structure */
   const char*      filename,           /**< name of file with graph data */
   double           scaleval,           /**< value to scale weights (only integral part of scaled weights is considered) */
   char*            probname            /**< buffer to store the name of the problem */
   )
{
   FILE* file;
   double weight;
   int node1;
   int node2;
   int currentnode;
   int i;
   
   assert(tcliquedata != NULL);
   assert(scaleval > 0.0);

   /* open file */
   if( (file = fopen(filename, "r")) == NULL )
   {
      if( (file = fopen("default.dat", "r")) == NULL )
      {
         printf("\nCan't open file: %s", filename);
         return FALSE;
      }
   }
 
   if( !tcliqueCreate(tcliquedata) )
      return FALSE;
 
   /* set name of problem, number of nodes and number of edges in graph */
   fscanf(file, "%s", probname);
   fscanf(file, "%d", &(*tcliquedata)->nnodes);
   fscanf(file, "%d", &(*tcliquedata)->nedges);
   
   /* set data structures for tclique */
   ALLOC_FALSE( allocMemoryArray(&(*tcliquedata)->weights, (*tcliquedata)->nnodes) );
   ALLOC_FALSE( allocMemoryArray(&(*tcliquedata)->degrees, (*tcliquedata)->nnodes) );
   ALLOC_FALSE( allocMemoryArray(&(*tcliquedata)->adjnodes, (*tcliquedata)->nedges) );
   ALLOC_FALSE( allocMemoryArray(&(*tcliquedata)->adjedges, (*tcliquedata)->nnodes) );

   /* set weights of all nodes (scaled!) */
   for( i = 0; i < (*tcliquedata)->nnodes; i++ )
   {
      fscanf(file, "%lf", &weight);
      (*tcliquedata)->weights[i] = (WEIGHT)(weight * scaleval);
      assert((*tcliquedata)->weights[i] >= 0);
   }

   /* set adjacent edges and degree of all nodes */
   currentnode = -1;
   for( i = 0; i < (*tcliquedata)->nedges; i++ )
   {
      /* read edge (node1, node2) */
      fscanf(file, "%d%d", &node1, &node2);
      
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

   /* close file */
   fclose(file);
   
   return TRUE;
}

/** saves tclique data structure to file */
Bool tcliqueSaveFile(
   TCLIQUEDATA*     tcliquedata,        /**< tclique data structure */
   const char*      filename,           /**< name of file to create */
   double           scaleval,           /**< value to unscale weights with */
   const char*      probname            /**< name of the problem */
   )
{
   FILE* file;
   int i;
   int j;

   assert(tcliquedata != NULL);
   assert(scaleval > 0.0);

   /* create file */
   if( (file = fopen(filename, "w")) == NULL )
   {
      printf("\nCan't create file: %s", filename);
      return FALSE;
   }
 
   /* write name of problem, number of nodes and number of edges in graph */
   fprintf(file, "%s\n", probname);
   fprintf(file, "%d\n", tcliquedata->nnodes);
   fprintf(file, "%d\n", tcliquedata->nedges);
   
   /* write weights of all nodes (scaled!) */
   for( i = 0; i < tcliquedata->nnodes; i++ )
      fprintf(file, "%f\n", (double)tcliquedata->weights[i]/scaleval);

   /* write edges */
   for( i = 0; i < tcliquedata->nnodes; i++ )
   {
      for( j = tcliquedata->adjedges[i].first; j < tcliquedata->adjedges[i].last; j++ )
         fprintf(file, "%d %d\n", i, tcliquedata->adjnodes[j]);
   }

   /* close file */
   fclose(file);
   
   return TRUE;
}

/** gets number of nodes in the graph */
int tcliqueGetNNodes(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);

   return tcliquedata->nnodes;
}

/** gets number of edges in the graph */
int tcliqueGetNEdges(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);
   
   return tcliquedata->nedges + tcliquedata->ncachededges;
}

/** gets weight of nodes in the graph */
WEIGHT* tcliqueGetWeights(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);
   
   return tcliquedata->weights;
}

/** gets degree of nodes in graph */
int* tcliqueGetDegrees(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);
   assert(tcliquedata->ncachededges == 0);
   
   return tcliquedata->degrees;
}

/** gets adjacent nodes of edges in graph */
int* tcliqueGetAdjnodes(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);
   assert(tcliquedata->ncachededges == 0);
   
   return tcliquedata->adjnodes;
}

/** gets pointer to first and one after last adjacent edge of nodes in graph */
HEAD_ADJ* tcliqueGetAdjedges(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   assert(tcliquedata != NULL);
   assert(tcliquedata->ncachededges == 0);
   
   return tcliquedata->adjedges;
}

/** gets pointer to first adjacent edge of given node in graph */
int* tcliqueGetFirstAdjedge(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node                /**< given node */
   )
{
   HEAD_ADJ* adjedges;
   int* adjnodes;

   assert(tcliquedata != NULL);
   assert(tcliquedata->ncachededges == 0);
   assert(0 <= node && node < tcliquedata->nnodes);

   adjedges = tcliqueGetAdjedges(tcliquedata);
   assert(adjedges != NULL);
   assert(adjedges[node].first >= 0);
   assert(adjedges[node].first <= tcliqueGetNEdges(tcliquedata));

   adjnodes = tcliqueGetAdjnodes(tcliquedata);
   assert(adjnodes != NULL);

   return &adjnodes[adjedges[node].first];
}

/** gets pointer to last adjacent edge of given node in graph */
int* tcliqueGetLastAdjedge(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node                /**< given node */
   )
{
   HEAD_ADJ* adjedges;
   int* adjnodes;
   int* degrees;

   assert(tcliquedata != NULL);
   assert(tcliquedata->ncachededges == 0);
   assert(0 <= node && node < tcliquedata->nnodes);

   adjedges = tcliqueGetAdjedges(tcliquedata);
   degrees = tcliqueGetDegrees(tcliquedata);
   assert(adjedges != NULL);
   assert(degrees[node] == 0 || adjedges[node].last-1 >= 0);
   assert(adjedges[node].last-1 <= tcliqueGetNEdges(tcliquedata));

   assert(adjedges[node].last - adjedges[node].first == degrees[node]);

   adjnodes = tcliqueGetAdjnodes(tcliquedata);
   assert(adjnodes != NULL);

   return &adjnodes[adjedges[node].last-1];
}

/** returns, whether the edge (node1, node2) is in the graph */
int tcliqueIsEdge(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node1,              /**< start node of edge */
   int              node2               /**< end node of edge */
   )
{
   int* currentadjedge;
   int* lastadjedge;
   int tmp;

   assert(tcliquedata != NULL);
   assert(tcliquedata->ncachededges == 0);
   assert(0 <= node1 && node1 < tcliquedata->nnodes);
   assert(0 <= node2 && node2 < tcliquedata->nnodes);

   if( node1 < node2 )
   {
      tmp = node1;
      node1 = node2;
      node2 = tmp;
   }

   currentadjedge = tcliqueGetFirstAdjedge(tcliquedata, node1);
   lastadjedge = tcliqueGetLastAdjedge(tcliquedata, node1);
   
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
 * and returns the number of selected nodes
 */
int tcliqueSelectAdjnodes( 
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
   assert(tcliquedata->ncachededges == 0);
   assert(0 <= node && node < tcliquedata->nnodes);
   assert(nnodes == 0 || nodes != NULL);
   assert(adjnodes != NULL);

   nadjnodes = 0;
   currentadjedge = tcliqueGetFirstAdjedge(tcliquedata, node);
   lastadjedge = tcliqueGetLastAdjedge(tcliquedata, node);

   /* checks for each node in given set nodes, if it is adjacent to given node 
    * (adjacent nodes are ordered by node index)
    */
   for( i = 0; i < nnodes; i++ )
   {
      assert(0 <= nodes[i] && nodes[i] < tcliquedata->nnodes);
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
void tcliquePrintData(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   )
{
   int* weights;
   int* degrees;
   int i;

   assert(tcliquedata != NULL);
   assert(tcliquedata->ncachededges == 0);

   degrees = tcliqueGetDegrees(tcliquedata);
   weights = tcliqueGetWeights(tcliquedata);

   printf("nnodes=%d, nedges=%d\n", tcliqueGetNNodes(tcliquedata), tcliqueGetNEdges(tcliquedata));
   for( i = 0; i < tcliqueGetNNodes(tcliquedata); i++ )
   {
      int* currentadjedge;
      int* lastadjedge;

      printf("node %d: weight=%d, degree=%d, adjnodes=\n[ ", i, weights[i], degrees[i]);  

      currentadjedge = tcliqueGetFirstAdjedge(tcliquedata, i);
      lastadjedge = tcliqueGetLastAdjedge(tcliquedata, i);
      assert(lastadjedge + 1 - currentadjedge == degrees[i]);

      for( ; currentadjedge <= lastadjedge; currentadjedge++ )
      {
         printf("%d, ", *currentadjedge);
      }
      printf("]\n");
   }
}
