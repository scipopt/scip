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


/** unsolved leaf of the tree */
struct Leaf
{
   LPSTATE*         lpstate;            /**< pointer to LP state information */
};

/** solved fork, where rows were only added to the father */
struct Fork
{
   COL**            addedCols;          /**< array with pointers to new columns added at this node into the LP */
   ROW**            addedRows;          /**< array with pointers to new rows added at this node into the LP */
   unsigned int     naddedCols:20;      /**< number of columns added at this node */
   unsigned int     naddedRows:20;      /**< number of rows added at this node */
   unsigned int     nchildren:16;       /**< number of children of this parent node */
};

/** solved fork, where rows were added, deleted, or rearranged */
struct Subroot
{
   COL**            cols;               /**< array with pointers to the columns in the same order as in the LP */
   ROW**            rows;               /**< array with pointers to the rows in the same order as in the LP */
   unsigned int     ncols:20;           /**< number of columns in the LP */
   unsigned int     nrows:20;           /**< number of rows in the LP */
   unsigned int     nchildren:16;       /**< number of children of this parent node */
};

/** node data structure */
struct Node
{
   union
   {
      LEAF*         leaf;               /**< data for leaf nodes */
      FORK*         fork;               /**< data for fork nodes */
      SUBROOT*      subroot;            /**< data for subroot nodes */
   } data;
   NODE*            parent;             /**< parent node in the tree */
   CONSLIST*        conslist;           /**< full list of active constraints */
   DOMCHG*          domchg;             /**< domain changes at this node or NULL */
   Real             lowerbound;         /**< lower (dual) LP bound of subtree */
   unsigned int     depth:16;           /**< depth in the tree */
   unsigned int     nodetype:2;         /**< type of node */
   unsigned int     active:1;           /**< is node in the path to the actual active node? */
};

/** branch and bound tree */
struct Tree
{
   NODE*            root;               /**< root node of the tree */
   PQUEUE*          leaves;             /**< leaves of the tree */
   NODE**           path;               /**< array of fork/subtree nodes storing the actual path from root to leaf */
   int*             pathnlpcols;        /**< array with number of LP columns for each problem in active path */
   int*             pathnlprows;        /**< array with number of LP rows for each problem in active path */
   int              pathlen;            /**< length of the actual path (== depth of the current node + 1) */
   int              pathsize;           /**< number of available slots in path arrays */
};




/*
 * Node methods
 */

DECL_SORTPTRCOMP(SCIPnodeCmpLowerbound) /**< node comparator for best lower bound */
{
   assert(elem1 != NULL);
   assert(elem2 != NULL);

   if( ((NODE*)elem1)->lowerbound < ((NODE*)elem2)->lowerbound )
      return -1;
   else if( ((NODE*)elem1)->lowerbound > ((NODE*)elem2)->lowerbound )
      return +1;
   else
      return 0;
}

static
void parentRelease(                     /**< decreases number of children, frees parent if no children left */
   NODE*            node,               /**< parent node */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(mem != NULL);
   assert(node != NULL);

   switch( node->nodetype )
   {
   case SCIP_NODETYPE_ACTNODE:
      errorMessage("Cannot release the parent-child relationship, if parent is the active node");
      abort();
   case SCIP_NODETYPE_LEAF:
      errorMessage("Node node is a leaf");
      abort();
   case SCIP_NODETYPE_FORK:
      assert(node->data.fork != NULL);
      assert(node->data.fork->nchildren > 0);
      node->data.fork->nchildren--;
      if( node->data.fork->nchildren == 0 )
         SCIPnodeFree(&node, mem, set);
      break;
   case SCIP_NODETYPE_SUBROOT:
      assert(node->data.subroot != NULL);
      assert(node->data.subroot->nchildren > 0);
      node->data.subroot->nchildren--;
      if( node->data.subroot->nchildren == 0 )
         SCIPnodeFree(&node, mem, set);
      break;
   default:
      errorMessage("Unknown node type");
      abort();
   }
}

static
void assignParent(                      /**< makes node a child of the given parent node, which must be the active node */
   NODE*            node,               /**< child node */
   NODE*            parent              /**< parent (= active) node (or NULL, if node is root) */
   )
{
   assert(node != NULL);
   assert(node->parent == NULL);
   assert(parent == NULL || parent->nodetype == SCIP_NODETYPE_ACTNODE);

   node->parent = parent;
   if( parent != NULL )
   {
      assert(parent->conslist != NULL); /* we need at least one constraint! */
      node->conslist = parent->conslist;
      node->lowerbound = parent->lowerbound;
      node->depth = parent->depth+1;
   }
}

static
void dismissParent(                     /**< releases the parent-child relation of the given child node */
   NODE*            node,               /**< child node */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(mem != NULL);
   assert(node != NULL);

   if( node->parent != NULL )
   {
      parentRelease(node->parent, mem, set);
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
   SCIPlpstateCapture(lpstate);
}

static
void dismissLPState(                    /**< releases the LP state from the given child node */
   LEAF*            leaf,               /**< leaf data of the node */
   MEM*             mem                 /**< block memory buffers */
   )
{
   assert(mem != NULL);
   assert(leaf != NULL);
   assert(leaf->lpstate != NULL);

   SCIPlpstateRelease(&(leaf->lpstate), mem);
}


static
LEAF* leafCreate(                       /**< creates leaf data */
   MEM*             mem                 /**< block memory buffers */
   )
{
   LEAF* leaf;

   assert(mem != NULL);

   ALLOC_NULL( allocBlockMemory(mem->treemem, leaf) );
   leaf->lpstate = NULL;

   return leaf;
}

static
void leafFree(                          /**< frees leaf data */
   LEAF**           leaf,               /**< leaf data */
   MEM*             mem                 /**< block memory buffers */
   )
{
   int r;

   assert(mem != NULL);
   assert(leaf != NULL);
   assert(*leaf != NULL);

   dismissLPState(*leaf, mem);
   freeBlockMemory(mem->treemem, *leaf);
}

static
FORK* forkCreate(                       /**< creates fork data */
   MEM*             mem                 /**< block memory buffers */
   )
{
   FORK* fork;

   assert(mem != NULL);

   ALLOC_NULL( allocBlockMemory(mem->treemem, fork) );
   fork->addedCols = NULL;
   fork->addedRows = NULL;
   fork->naddedCols = 0;
   fork->naddedRows = 0;
   fork->nchildren = 0;

   return fork;
}

static
void forkFree(                          /**< frees fork data */
   FORK**           fork,               /**< fork data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   Bool             releaseData         /**< should the columns and rows of the node's arrays be released? */
   )
{
   assert(mem != NULL);
   assert(fork != NULL);
   assert(*fork != NULL);

   if( releaseData )
   {
      int i;
      
      for( i = 0; i < (*fork)->naddedCols; ++i )
         SCIPcolRelease(&((*fork)->addedCols[i]), mem, set);
      for( i = 0; i < (*fork)->naddedRows; ++i )
         SCIProwRelease(&((*fork)->addedRows[i]), mem, set);
   }

   freeBlockMemoryArrayNull(mem->treemem, (*fork)->addedCols, (*fork)->naddedCols);
   freeBlockMemoryArrayNull(mem->treemem, (*fork)->addedRows, (*fork)->naddedRows);
   freeBlockMemory(mem->treemem, *fork);
}

static
SUBROOT* subrootCreate(                 /**< creates subroot data */
   MEM*             mem                 /**< block memory buffers */
   )
{
   SUBROOT* subroot;

   assert(mem != NULL);

   ALLOC_NULL( allocBlockMemory(mem->treemem, subroot) );
   subroot->cols = NULL;
   subroot->rows = NULL;
   subroot->ncols = 0;
   subroot->nrows = 0;
   subroot->nchildren = 0;

   return subroot;
}

static
void subrootFree(                       /**< frees subroot */
   SUBROOT**        subroot,            /**< subroot data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   Bool             releaseData         /**< should the columns and rows of the node's arrays be released? */
   )
{
   assert(mem != NULL);
   assert(subroot != NULL);
   assert(*subroot != NULL);

   if( releaseData )
   {
      int i;
      
      for( i = 0; i < (*subroot)->ncols; ++i )
         SCIPcolRelease(&((*subroot)->cols[i]), mem, set);
      for( i = 0; i < (*subroot)->nrows; ++i )
         SCIProwRelease(&((*subroot)->rows[i]), mem, set);
   }

   freeBlockMemoryArrayNull(mem->treemem, (*subroot)->cols, (*subroot)->ncols);
   freeBlockMemoryArrayNull(mem->treemem, (*subroot)->rows, (*subroot)->nrows);
   freeBlockMemory(mem->treemem, *subroot);
}

NODE* SCIPnodeCreate(                   /**< creates a leaf node */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   NODE*            parent,             /**< parent node in the tree */
   LPSTATE*         lpstate             /**< pointer to LP state information */
   )
{
   NODE* node;

   assert(mem != NULL);
   assert(lpstate != NULL);

   ALLOC_NULL( allocBlockMemory(mem->treemem, node) );
   ALLOC_NULL( node->data.leaf = leafCreate(mem) );
   node->parent = NULL;
   node->conslist = NULL;
   node->domchg = NULL;
   node->lowerbound = -SCIPinfinity(set);
   node->depth = 0;
   node->nodetype = SCIP_NODETYPE_LEAF;
   node->active = FALSE;

   assignParent(node, parent);
   assignLPState(node->data.leaf, lpstate);

   return node;
}

void SCIPnodeFree(                      /**< frees node */
   NODE**           node,               /**< node data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(mem != NULL);
   assert(node != NULL);
   assert(*node != NULL);

   /* free nodetype specific data */
   switch((*node)->nodetype)
   {
   case SCIP_NODETYPE_ACTNODE:
      break;
   case SCIP_NODETYPE_LEAF:
      leafFree(&((*node)->data.leaf), mem);
      break;
   case SCIP_NODETYPE_FORK:
      forkFree(&((*node)->data.fork), mem, set, TRUE);
      break;
   case SCIP_NODETYPE_SUBROOT:
      subrootFree(&((*node)->data.subroot), mem, set, TRUE);
      break;
   default:
      errorMessage("Unknown node type");
      break;
   }

   /* free common data */
   SCIPconslistFreePart(&((*node)->conslist), mem, (*node)->parent->conslist);
   SCIPdomchgFree(&((*node)->domchg), mem);
   dismissParent(*node, mem, set);

   freeBlockMemory(mem->treemem, *node);
}

static
RETCODE treeEnsurePathMem(              /**< resizes path array to be able to store at least num nodes */
   TREE*            tree,               /**< branch-and-bound tree */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in path */
   )
{
   assert(tree != NULL);
   assert(mem != NULL);
   assert(set != NULL);

   if( num > tree->pathsize )
   {
      int newsize;

      newsize = SCIPcalcPathGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(mem->treemem, tree->path, tree->pathsize, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(mem->treemem, tree->pathnlpcols, tree->pathsize, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(mem->treemem, tree->pathnlprows, tree->pathsize, newsize) );
      tree->pathsize = newsize;
   }
   assert(num <= tree->pathsize);

   return SCIP_OKAY;
}

static
RETCODE switchPath(                     /**< performs LP changes on the path from oldnode to node */
   NODE*            node,               /**< fork/subtree node to construct LP for */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   NODE*            oldnode             /**< old active node */
   )
{
   NODE* forknode;
   NODE* subroot;
   Bool samesubtree;
   Bool colsshrinked;
   Bool rowsshrinked;
   int forkdepth;
   int lpdepth;
   int correctcols;
   int correctrows;
   int d;

   assert(mem != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(tree != NULL);

   /* if oldnode == NULL, then node has to be the parent of the root node (NULL) and LP is empty */
   if( oldnode == NULL )
   {
      assert(node == NULL);
      assert(lp->ncols == 0);
      assert(lp->nrows == 0);
      assert(tree->pathlen == 0);
      return SCIP_OKAY;
   }

   assert(node != NULL && oldnode != NULL);
   assert(node->nodetype == SCIP_NODETYPE_FORK || node->nodetype == SCIP_NODETYPE_SUBROOT);
   assert(tree->pathlen == oldnode->depth+1);

   /* search the common fork node of the two nodes
    * detect the subroot of the new node, if existing
    * construct the new active path
    */
   CHECK_OKAY( treeEnsurePathMem(tree, mem, set, node->depth+1) );
   tree->pathlen = node->depth+1;
   forknode = node;
   subroot = NULL;
   while( forknode != NULL && !forknode->active )
   {
      assert(forknode->depth < tree->pathsize);
      tree->path[forknode->depth] = forknode;
      forknode->active = TRUE;
      if( subroot == NULL && forknode->nodetype == SCIP_NODETYPE_SUBROOT )
         subroot = forknode;
      assert(forknode->parent == NULL || forknode->parent->depth == forknode->depth-1);
      forknode = forknode->parent;
   }
   assert(forknode != NULL);
   assert(forknode->depth < tree->pathlen);
   assert(forknode->depth <= oldnode->depth);
   forkdepth = forknode->depth;

   /* undo the domain changes of the old subpath 
    * check, if we are in the same subtree
    * deactivate old subpath
    */
   samesubtree = (subroot == NULL);
   while( oldnode != NULL && oldnode != forknode )
   {
      assert(oldnode->active);
      CHECK_OKAY( SCIPlpUndoDomchg(lp, mem, set, oldnode->domchg) );
      oldnode->active = FALSE;
      assert(oldnode->nodetype == SCIP_NODETYPE_ACTNODE
         || oldnode->nodetype == SCIP_NODETYPE_FORK
         || oldnode->nodetype == SCIP_NODETYPE_SUBROOT);
      samesubtree &= (oldnode->nodetype != SCIP_NODETYPE_SUBROOT);
      oldnode = oldnode->parent;
   }
   assert(oldnode == forknode);

   /* construct the new LP */
   lpdepth = -1;
   correctcols = 0;
   correctrows = 0;
   colsshrinked = FALSE;
   rowsshrinked = FALSE;
   if( samesubtree )
   {
      /* same subtree */
      assert(forknode != NULL);
      assert(tree->path[forkdepth] == forknode);
      assert(0 <= tree->pathnlpcols[forkdepth] && tree->pathnlpcols[forkdepth] <= lp->ncols);
      assert(0 <= tree->pathnlprows[forkdepth] && tree->pathnlprows[forkdepth] <= lp->nrows);

      lpdepth = forkdepth;
      correctcols = tree->pathnlpcols[forkdepth];
      correctrows = tree->pathnlprows[forkdepth];
   }
   else if( subroot != NULL )
   {
      COL** cols;
      ROW** rows;
      int ncols;
      int nrows;

      /* different subtree with subroot */
      assert(subroot->nodetype == SCIP_NODETYPE_SUBROOT);
      assert(subroot->data.subroot != NULL);

      cols = subroot->data.subroot->cols;
      rows = subroot->data.subroot->rows;
      ncols = subroot->data.subroot->ncols;
      nrows = subroot->data.subroot->nrows;
      
      /* keep columns and rows, which wouldn't change position in new LP */
      for( correctcols = 0; correctcols < ncols && cols[correctcols]->lppos == correctcols; ++correctcols )
      {
         assert(correctcols < lp->ncols);
      }
      for( correctrows = 0; correctrows < nrows && rows[correctrows]->lppos == correctrows; ++correctrows )
      {
         assert(correctrows < lp->nrows);
      }

      /* shrink LP to the unchanged size, and extend LP with the rest of columns and rows in subroot data */
      if( correctcols < ncols )
      {
         CHECK_OKAY( SCIPlpShrinkCols(lp, correctcols) );
         colsshrinked = TRUE;
         for( ; correctcols < ncols; ++correctcols )
         {
            SCIPlpAddCol(lp, mem, set, cols[correctcols]);
         }
      }
      if( correctrows < nrows )
      {
         CHECK_OKAY( SCIPlpShrinkRows(lp, correctrows) );
         rowsshrinked = TRUE;
         for( ; correctrows < nrows; ++correctrows )
         {
            SCIPlpAddRow(lp, mem, set, rows[correctrows]);
         }
      }
      
      lpdepth = subroot->depth;
   }
   else
   {
      /* different subtree without subroot */
   }
   /* now, lpdepth is the depth until the LP is already constructed */
   assert(lpdepth <= node->depth);
   assert(lpdepth == -1 || (correctcols == tree->pathnlpcols[lpdepth] && correctrows == tree->pathnlprows[lpdepth]));
   assert(lpdepth >= 0 || (correctcols == 0 && correctrows == 0));

   /* add columns and rows of new path to LP
    * from here on, the path has to consist of fork nodes only
    */
   for( d = lpdepth+1; d <= node->depth; ++d )
   {
      COL** cols;
      ROW** rows;
      int ncols;
      int nrows;
      int col;
      int row;

      assert(d < tree->pathlen);
      assert(tree->path[d] != NULL);
      assert(tree->path[d]->nodetype == SCIP_NODETYPE_FORK);
      assert(tree->path[d]->data.fork != NULL);
      assert(tree->path[d]->active);
      assert(tree->path[d]->depth == d);

      cols = tree->path[d]->data.fork->addedCols;
      ncols = tree->path[d]->data.fork->naddedCols;
      col = 0;
      rows = tree->path[d]->data.fork->addedRows;
      nrows = tree->path[d]->data.fork->naddedRows;
      row = 0;

      if( !colsshrinked )
      {
         /* keep columns, which wouldn't change position in new LP */
         for( ; col < ncols && cols[col]->lppos == correctcols; ++col, ++correctcols )
         {
            assert(correctcols < lp->ncols);
         }

         /* if changed cols are left, shrink LP to the unchanged size */
         if( col < ncols )
         {
            CHECK_OKAY( SCIPlpShrinkCols(lp, correctcols) );
            colsshrinked = TRUE;
         }         
      }
      assert(colsshrinked || col == ncols);

      if( !rowsshrinked )
      {
         /* keep rows, which wouldn't change position in new LP */
         for( ; row < nrows && rows[row]->lppos == correctrows; ++row, ++correctrows )
         {
            assert(correctrows < lp->nrows);
         }

         /* if changed rows are left, shrink LP to the unchanged size */
         if( row < nrows )
         {
            CHECK_OKAY( SCIPlpShrinkRows(lp, correctrows) );
            rowsshrinked = TRUE;
         }         
      }
      assert(rowsshrinked || row == nrows);

      /* add additional columns and rows */
      for( ; col < ncols; ++col, ++correctcols )
      {
         assert(colsshrinked);
         CHECK_OKAY( SCIPlpAddCol(lp, mem, set, cols[col]) );
      }
      for( ; row < nrows; ++row, ++correctrows )
      {
         assert(rowsshrinked);
         CHECK_OKAY( SCIPlpAddRow(lp, mem, set, rows[row]) );
      }

      assert(tree->pathnlpcols[d] == correctcols);
      assert(tree->pathnlprows[d] == correctrows);
   }

   /* if there are additional columns or rows left, delete them */
   if( !colsshrinked )
   {
      CHECK_OKAY( SCIPlpShrinkCols(lp, correctcols) );
      colsshrinked = TRUE;
   }
   if( !rowsshrinked )
   {
      CHECK_OKAY( SCIPlpShrinkRows(lp, correctrows) );
      rowsshrinked = TRUE;
   }

   return SCIP_OKAY;
}

static
RETCODE leafToActnode(                  /**< converts a leaf node into an actnode node */
   NODE*            node,               /**< node to convert */
   MEM*             mem                 /**< block memory buffers */
   )
{
   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_LEAF);
   assert(mem != NULL);

   leafFree(&(node->data.leaf), mem);
   node->nodetype = SCIP_NODETYPE_ACTNODE;

   return SCIP_OKAY;
}

RETCODE SCIPnodeActivate(               /**< activates a leaf node */
   NODE*            node,               /**< leaf node to activate */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   NODE*            oldnode             /**< old active node */
   )
{
   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_LEAF);
   assert(node->data.leaf != NULL);
   assert(!node->active);
   assert(mem != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(tree != NULL);

   /* track the path from the old node to the new node's parent, and perform LP changes */
   CHECK_OKAY( switchPath(node->parent, mem, set, lp, tree, oldnode) );

   /* load LP state into solver */
   if( node->data.leaf->lpstate != NULL )
   {
      CHECK_OKAY( SCIPlpSetState(lp, mem, node->data.leaf->lpstate) );
   }

   /* convert node into the actual node */
   CHECK_OKAY( leafToActnode(node, mem) );
   node->active = TRUE;

   /* change domains of variables */
   if( node->domchg != NULL )
      SCIPlpApplyDomchg(lp, mem, set, node->domchg);

   /* insert node in path */
   CHECK_OKAY( SCIPtreeExtendPath(tree, mem, set, node) );

   /* flush changes to the LP */
   CHECK_OKAY( SCIPlpFlush(lp, mem, set) );

   return SCIP_OKAY;
}   

RETCODE SCIPactnodeToFork(              /**< converts the active node into a fork node */
   NODE*            node,               /**< node to convert */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   FORK* fork;

   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(node->active);
   assert(mem != NULL);
   assert(tree != NULL);
   assert(lp != NULL);

   ALLOC_OKAY( fork = forkCreate(mem) );
   
   fork->naddedCols = SCIPlpGetNumNewcols(lp);
   fork->naddedRows = SCIPlpGetNumNewrows(lp);

   if( fork->naddedCols > 0 )
   {
      /* copy the newly created columns to the fork's col array */
      ALLOC_OKAY( duplicateBlockMemoryArray(mem->treemem, fork->addedCols, SCIPlpGetNewcols(lp), fork->naddedCols) );
   }
   if( fork->naddedRows > 0 )
   {
      /* copy the newly created rows to the fork's row array */
      ALLOC_OKAY( duplicateBlockMemoryArray(mem->treemem, fork->addedRows, SCIPlpGetNewrows(lp), fork->naddedRows) );
   }
   
   node->nodetype = SCIP_NODETYPE_FORK;
   node->data.fork = fork;

   /* put node in the path of active forks and subroots */
   CHECK_OKAY( SCIPtreeExtendPath(tree, mem, set, node) );

   return SCIP_OKAY;
}

RETCODE SCIPactnodeToSubroot(           /**< converts the active node into a subroot node */
   NODE*            node,               /**< node to convert */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   SUBROOT* subroot;

   assert(mem != NULL);
   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(node->active);
   assert(set != NULL);
   assert(tree != NULL);
   assert(lp != NULL);

   ALLOC_OKAY( subroot = subrootCreate(mem) );

   subroot->ncols = lp->ncols;
   subroot->nrows = lp->nrows;
   subroot->nchildren = node->data.fork->nchildren;
   ALLOC_OKAY( duplicateBlockMemoryArray(mem->treemem, subroot->cols, lp->cols, subroot->ncols) );
   ALLOC_OKAY( duplicateBlockMemoryArray(mem->treemem, subroot->rows, lp->rows, subroot->nrows) );

   node->nodetype = SCIP_NODETYPE_SUBROOT;
   node->data.subroot = subroot;

   /* put node in the path of active forks and subroots */
   CHECK_OKAY( SCIPtreeExtendPath(tree, mem, set, node) );

   return SCIP_OKAY;
}


/*
 * Tree methods
 */

TREE* SCIPtreeCreate(                   /**< creates an initialized tree data structure */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   )
{
   TREE* tree;

   assert(mem != NULL);
   assert(set != NULL);
   assert(set->treeGrowInit >= 0);
   assert(set->treeGrowFac >= 1.0);

   ALLOC_NULL( allocBlockMemory(mem->treemem, tree) );
   tree->root = NULL;
   CHECK_NULL( SCIPpqueueInit(&(tree->leaves), set->treeGrowInit, set->treeGrowFac, set->nodecmp) );
   tree->pathlen = 0;
   tree->pathsize = SCIPcalcPathGrowSize(set, 0);
   ALLOC_NULL( allocBlockMemoryArray(mem->treemem, tree->path, tree->pathsize) );
   ALLOC_NULL( allocBlockMemoryArray(mem->treemem, tree->pathnlpcols, tree->pathsize) );
   ALLOC_NULL( allocBlockMemoryArray(mem->treemem, tree->pathnlprows, tree->pathsize) );
   
   return tree;
}

RETCODE SCIPtreeExtendPath(             /**< appends fork/subtree node to the path of active nodes and marks node active */
   TREE*            tree,               /**< branch-and-bound tree */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   NODE*            node                /**< fork or subtree node to append */
   )
{
   assert(tree != NULL);
   assert(mem != NULL);
   assert(set != NULL);
   assert(node != NULL);
   assert(node->depth == tree->pathlen);
   assert(node->nodetype == SCIP_NODETYPE_FORK || node->nodetype == SCIP_NODETYPE_SUBROOT);
   assert(!node->active);
   assert(node->depth == 0 || tree->path != NULL);
   assert(node->depth == 0 || tree->path[node->depth-1] == node->parent);

   CHECK_OKAY( treeEnsurePathMem(tree, mem, set, tree->pathlen+1) );
   tree->path[tree->pathlen] = node;
   switch( node->nodetype )
   {
   case SCIP_NODETYPE_FORK:
      assert(node->data.fork != NULL);
      if( tree->pathlen == 0 )
      {
         tree->pathnlpcols[0] = 0;
         tree->pathnlprows[0] = 0;
      }
      else
      {
         tree->pathnlpcols[tree->pathlen] = tree->pathnlpcols[tree->pathlen-1];
         tree->pathnlprows[tree->pathlen] = tree->pathnlprows[tree->pathlen-1];
      }
      tree->pathnlpcols[tree->pathlen] += node->data.fork->naddedCols;
      tree->pathnlprows[tree->pathlen] += node->data.fork->naddedRows;
      break;
   case SCIP_NODETYPE_SUBROOT:
      assert(node->data.subroot != NULL);
      tree->pathnlpcols[tree->pathlen] = node->data.subroot->ncols;
      tree->pathnlprows[tree->pathlen] = node->data.subroot->nrows;
      break;
   case SCIP_NODETYPE_LEAF:
      errorMessage("Cannot store leaf node in path of active forks");
      abort();
   case SCIP_NODETYPE_ACTNODE:
      errorMessage("Cannot store active node in path of active forks");
      abort();
   default:
      errorMessage("Unknown node type");
      abort();
   }
   tree->pathlen++;
   node->active = TRUE;

   return SCIP_OKAY;
}

void SCIPtreeShrinkPath(                /**< cuts off path of active nodes after given node and marks them inactive */
   TREE*            tree,               /**< branch-and-bound tree */
   NODE*            node                /**< last node in the path */
   )
{
   int i;

   assert(tree != NULL);
   assert(node != NULL);
   assert(node->depth < tree->pathlen);
   assert(node->active);
   assert(tree->path[node->depth] == node);

   for( i = node->depth+1; i < tree->pathlen; ++i )
   {
      assert(tree->path[i] != NULL);
      assert(tree->path[i]->depth == i);
      tree->path[i]->active = FALSE;
   }
   tree->pathlen = node->depth+1;

   assert(tree->pathlen <= tree->pathsize);
}

