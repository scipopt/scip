/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   tree.c
 * @brief  branch-and-bound tree datastructures and operations
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "def.h"
#include "constraint.h"
#include "sort.h"
#include "tree.h"


struct Leaf                             /**< unsolved leaf of the tree */
{
   LPSTATE*         lpstate;            /**< pointer to LP state information */
   ROW**            addedRows;          /**< array with pointers to new rows added at this node into the LP */
   unsigned int     naddedRows:20;      /**< number of rows added at this node */
};

struct Fork                             /**< solved fork, where rows were only added to the father */
{
   ROW**            addedRows;          /**< array with pointers to new rows added at this node into the LP */
   unsigned int     naddedRows:20;      /**< number of rows added at this node */
   unsigned int     nchildren:16;       /**< number of children of this parent node */
};

struct Subroot                          /**< solved fork, where rows were added, deleted, or rearranged */
{
   ROW**            rows;               /**< array with pointers to the rows in the same order as in the LP */
   unsigned int     nrows:20;           /**< number of rows in the LP */
   unsigned int     nchildren:16;       /**< number of children of this parent node */
};

struct Node                             /**< node data structure */
{
   union
   {
      LEAF*         leaf;               /**< data for leaf nodes */
      FORK*         fork;               /**< data for fork nodes */
      SUBROOT*      subroot;            /**< data for subroot nodes */
   } data;
   NODE*            parent;             /**< parent node in the tree */
   CONSLIST*        conslist;           /**< full list of active constraints */
   VARDOMCHG*       vardomChg;          /**< list of domain changes at this node */
   double           lowerbound;         /**< lower (dual) LP bound of subtree */
   unsigned int     depth:16;           /**< depth in the tree */
   unsigned int     nodetype:2;         /**< type of node */
   unsigned int     active:1;           /**< is node in the path to the actual active node? */
};

struct Tree                             /**< branch and bound tree */
{
   NODE*            root;               /**< root node of the tree */
   PQUEUE*          leaves;             /**< leaves of the tree */
};


static
void captureParent(                     /**< increases number of children of the given parent node */
   NODE*            node                /**< parent node */
   )
{
   assert(node != NULL);

   switch( node->nodetype )
   {
   case SCIP_NODETYPE_LEAF:
      errorMessage("Node node is a leaf");
      abort();
   case SCIP_NODETYPE_FORK:
      assert(node->data.fork != NULL);
      node->data.fork->nchildren++;
      break;
   case SCIP_NODETYPE_SUBROOT:
      assert(node->data.subroot != NULL);
      node->data.subroot->nchildren++;
      break;
   default:
      errorMessage("Unknown node type");
      abort();
   }
}

static
void releaseParent(                     /**< decreases number of children, frees parent if no children left */
   MEM*             mem,                /**< block memory buffers */
   NODE*            node                /**< parent node */
   )
{
   assert(mem != NULL);
   assert(node != NULL);

   switch( node->nodetype )
   {
   case SCIP_NODETYPE_LEAF:
      errorMessage("Node node is a leaf");
      abort();
   case SCIP_NODETYPE_FORK:
      assert(node->data.fork != NULL);
      assert(node->data.fork->nchildren > 0);
      node->data.fork->nchildren--;
      if( node->data.fork->nchildren == 0 )
         SCIPfreeNode(mem, &node);
      break;
   case SCIP_NODETYPE_SUBROOT:
      assert(node->data.subroot != NULL);
      assert(node->data.subroot->nchildren > 0);
      node->data.subroot->nchildren--;
      if( node->data.subroot->nchildren == 0 )
         SCIPfreeNode(mem, &node);
      break;
   default:
      errorMessage("Unknown node type");
      abort();
   }
}

static
void assignParent(                      /**< assigns the node to be a child of the given parent node */
   NODE*            node,               /**< child node */
   NODE*            parent              /**< parent node (or NULL, if node is root) */
   )
{
   assert(node != NULL);
   assert(node->parent == NULL);

   node->parent = parent;
   if( parent != NULL )
      captureParent(parent);
}

static
void dismissParent(                     /**< releases the parent-child relation of the given child node */
   MEM*             mem,                /**< block memory buffers */
   NODE*            node                /**< child node */
   )
{
   assert(mem != NULL);
   assert(node != NULL);

   if( node->parent != NULL )
   {
      releaseParent(mem, node->parent);
      node->parent = NULL;
   }
}

static
void assignLPState(                     /**< assigns the given LP state to the node */
   LEAF*            leaf,               /**< leaf data of the node */
   LPSTATE*         lpstate             /**< LP state to be assigned */
   )
{
   assert(leaf != NULL);
   assert(leaf->lpstate == NULL);

   leaf->lpstate = lpstate;
   SCIPcaptureLPState(lpstate);
}

static
void dismissLPState(                    /**< releases the LP state from the given child node */
   MEM*             mem,                /**< block memory buffers */
   LEAF*            leaf                /**< leaf data of the node */
   )
{
   assert(mem != NULL);
   assert(leaf != NULL);
   assert(leaf->lpstate != NULL);

   SCIPreleaseLPState(mem, &(leaf->lpstate));
}


static
LEAF* createLeafData(                   /**< creates leaf data */
   MEM*             mem                 /**< block memory buffers */
   )
{
   LEAF* leaf;

   assert(mem != NULL);

   ALLOC_NULL( allocBlockMemory(mem->treemem, leaf) );
   leaf->lpstate = NULL;
   leaf->addedRows = NULL;
   leaf->naddedRows = 0;

   return leaf;
}

static
void freeLeafData(                      /**< frees leaf data */
   MEM*             mem,                /**< block memory buffers */
   LEAF**           leaf                /**< leaf data */
   )
{
   int r;

   assert(mem != NULL);
   assert(leaf != NULL);
   assert(*leaf != NULL);

   dismissLPState(mem, *leaf);
   for( r = 0; r < (*leaf)->naddedRows; ++r )
      SCIPreleaseRow(mem, &((*leaf)->addedRows[r]));
   
   freeBlockMemoryArrayNull(mem->treemem, (*leaf)->addedRows, (*leaf)->naddedRows);
   freeBlockMemory(mem->treemem, *leaf);
}

static
FORK* createForkData(                   /**< creates fork data */
   MEM*             mem                 /**< block memory buffers */
   )
{
   FORK* fork;

   assert(mem != NULL);

   ALLOC_NULL( allocBlockMemory(mem->treemem, fork) );
   fork->addedRows = NULL;
   fork->naddedRows = 0;
   fork->nchildren = 0;

   return fork;
}

static
void freeForkData(                      /**< frees fork data */
   MEM*             mem,                /**< block memory buffers */
   FORK**           fork                /**< fork data */
   )
{
   int r;

   assert(mem != NULL);
   assert(fork != NULL);
   assert(*fork != NULL);

   for( r = 0; r < (*fork)->naddedRows; ++r )
      SCIPreleaseRow(mem, &((*fork)->addedRows[r]));
   
   freeBlockMemoryArrayNull(mem->treemem, (*fork)->addedRows, (*fork)->naddedRows);
   freeBlockMemory(mem->treemem, *fork);
}

static
SUBROOT* createSubrootData(             /**< creates subroot data */
   MEM*             mem                 /**< block memory buffers */
   )
{
   SUBROOT* subroot;

   assert(mem != NULL);

   ALLOC_NULL( allocBlockMemory(mem->treemem, subroot) );
   subroot->rows = NULL;
   subroot->nrows = 0;
   subroot->nchildren = 0;

   return subroot;
}

static
void freeSubrootData(                   /**< frees subroot */
   MEM*             mem,                /**< block memory buffers */
   SUBROOT**        subroot             /**< subroot data */
   )
{
   int r;

   assert(mem != NULL);
   assert(subroot != NULL);
   assert(*subroot != NULL);

   for( r = 0; r < (*subroot)->nrows; ++r )
      SCIPreleaseRow(mem, &((*subroot)->rows[r]));

   freeBlockMemoryArrayNull(mem->treemem, (*subroot)->rows, (*subroot)->nrows);
   freeBlockMemory(mem->treemem, *subroot);
}

NODE* SCIPcreateLeaf(                   /**< creates a leaf node */
   MEM*             mem,                /**< block memory buffers */
   NODE*            parent,             /**< parent node in the tree */
   LPSTATE*         lpstate             /**< pointer to LP state information */
   )
{
   NODE* node;

   assert(mem != NULL);
   assert(lpstate != NULL);

   ALLOC_NULL( allocBlockMemory(mem->treemem, node) );
   ALLOC_NULL( node->data.leaf = createLeafData(mem) );
   node->parent = NULL;
   node->conslist = NULL;
   node->vardomChg = NULL;

   if( parent == NULL )
   {
      node->lowerbound = -SCIP_INFINITY;
      node->depth = 0;
   }
   else
   {
      node->lowerbound = parent->lowerbound;
      node->depth = parent->depth+1;
   }

   node->nodetype = SCIP_NODETYPE_LEAF;
   node->active = FALSE;

   assignParent(node, parent);
   assignLPState(node->data.leaf, lpstate);

   return node;
}

void SCIPfreeNode(                      /**< frees node */
   MEM*             mem,                /**< block memory buffers */
   NODE**           node                /**< node data */
   )
{
   assert(mem != NULL);
   assert(node != NULL);
   assert(*node != NULL);

   /* free nodetype specific data */
   switch((*node)->nodetype)
   {
   case SCIP_NODETYPE_LEAF:
      freeLeafData(mem, &((*node)->data.leaf));
      break;
   case SCIP_NODETYPE_FORK:
      freeForkData(mem, &((*node)->data.fork));
      break;
   case SCIP_NODETYPE_SUBROOT:
      freeSubrootData(mem, &((*node)->data.subroot));
      break;
   default:
      errorMessage("Unknown node type");
      break;
   }

   /* free common data */
   dismissParent(mem, *node);
   SCIPfreeConslist(mem, &((*node)->conslist));
   SCIPfreeVardomChg(mem, &((*node)->vardomChg));

   freeBlockMemory(mem->treemem, *node);
}

RETCODE SCIPleafToFork(                 /**< converts a leaf into a fork node */
   MEM*             mem,                /**< block memory buffers */
   NODE*            node                /**< node to convert */
   )
{
   FORK* fork;

   assert(mem != NULL);
   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_LEAF);
   assert(node->data.leaf != NULL);

   ALLOC_OKAY( fork = createForkData(mem) );
   fork->addedRows = node->data.leaf->addedRows;
   fork->naddedRows = node->data.leaf->naddedRows;
   
   freeLeafData(mem, &(node->data.leaf));

   node->nodetype = SCIP_NODETYPE_FORK;
   node->data.fork = fork;

   return SCIP_OKAY;
}

RETCODE SCIPforkToSubroot(              /**< converts a fork into a subroot node */
   MEM*             mem,                /**< block memory buffers */
   NODE*            node,               /**< node to convert */
   LP*              lp                  /**< actual LP data */
   )
{
   SUBROOT* subroot;

   assert(mem != NULL);
   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_FORK);
   assert(node->data.fork != NULL);
   assert(lp != NULL);

   ALLOC_OKAY( subroot = createSubrootData(mem) );
   ALLOC_OKAY( duplicateBlockMemoryArray(mem->treemem, subroot->rows, lp->rows, lp->nrows) );
   subroot->nrows = lp->nrows;
   subroot->nchildren = node->data.fork->nchildren;

   freeForkData(mem, &(node->data.fork));

   node->nodetype = SCIP_NODETYPE_SUBROOT;
   node->data.subroot = subroot;

   return SCIP_OKAY;
}
