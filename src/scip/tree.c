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


/** fork without LP solution, where only bounds and constraints have been changed */
struct Junction
{
   int              nchildren;          /**< number of children of this parent node */
};

/** fork with solved LP, where bounds and constraints have been changed, and rows and columns were added */
struct Fork
{
   COL**            addedCols;          /**< array with pointers to new columns added at this node into the LP */
   ROW**            addedRows;          /**< array with pointers to new rows added at this node into the LP */
   LPSTATE*         lpstate;            /**< LP state information */
   int              naddedCols;         /**< number of columns added at this node */
   int              naddedRows;         /**< number of rows added at this node */
   int              nchildren;          /**< number of children of this parent node */
   int              nlpstateref;        /**< number of times, the LP state is needed */   
};

/** fork with solved LP, where bounds and constraints have been changed, and rows and columns were removed and added */
struct Subroot
{
   COL**            cols;               /**< array with pointers to the columns in the same order as in the LP */
   ROW**            rows;               /**< array with pointers to the rows in the same order as in the LP */
   LPSTATE*         lpstate;            /**< LP state information */
   int              ncols;              /**< number of columns in the LP */
   int              nrows;              /**< number of rows in the LP */
   int              nchildren;          /**< number of children of this parent node */
   int              nlpstateref;        /**< number of times, the LP state is needed */   
};

/** node data structure */
struct Node
{
   union
   {
      JUNCTION*     junction;           /**< data for junction nodes */
      FORK*         fork;               /**< data for fork nodes */
      SUBROOT*      subroot;            /**< data for subroot nodes */
   } data;
   NODE*            parent;             /**< parent node in the tree */
   CONSLIST*        conslist;           /**< full list of active constraints */
   DOMCHG*          domchg;             /**< domain changes at this node or NULL */
   Real             lowerbound;         /**< lower (dual) LP bound of subtree */
   unsigned int     depth:16;           /**< depth in the tree */
   unsigned int     nodetype:3;         /**< type of node */
   unsigned int     active:1;           /**< is node in the path to the actual active node? */
};

/** branch and bound tree */
struct Tree
{
   NODE*            root;               /**< root node of the tree */
   PQUEUE*          leaves;             /**< leaves of the tree */
   NODE**           path;               /**< array of fork/subtree nodes storing the active path from root to leaf */
   NODE*            actLPFork;          /**< fork/subroot node defining the LP state of the active node */
   NODE*            actSubroot;         /**< root of the active subtree */
   NODE**           children;           /**< array with children of the active node */
   NODE**           siblings;           /**< array with siblings of the active node */
   DOMCHGDYN*       domchgdyn;          /**< domain changes of the active node */
   int*             pathnlpcols;        /**< array with number of LP columns for each problem in active path */
   int*             pathnlprows;        /**< array with number of LP rows for each problem in active path */
   int              pathlen;            /**< length of the actual path (== depth of the current node + 1) */
   int              pathsize;           /**< number of available slots in path arrays */
   int              correctLPDepth;     /**< depth to which current LP data corresponds to LP data of active path */
   int              childrensize;       /**< available slots in children vector */
   int              nchildren;          /**< actual number of children (number of used slots in children vector) */
   int              siblingssize;       /**< available slots in siblings vector */
   int              nsiblings;          /**< actual number of siblings (number of used slots in siblings vector) */
};




/*
 * dymanic memory arrays
 */

static
RETCODE treeEnsureChildrenMem(          /**< resizes children array to be able to store at least num nodes */
   TREE*            tree,               /**< branch-and-bound tree */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(tree != NULL);
   assert(mem != NULL);
   assert(set != NULL);

   if( num > tree->childrensize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(mem->treemem, tree->children, tree->childrensize, newsize) );
      tree->childrensize = newsize;
   }
   assert(num <= tree->childrensize);

   return SCIP_OKAY;
}

static
RETCODE treeEnsureSiblingsMem(          /**< resizes siblings array to be able to store at least num nodes */
   TREE*            tree,               /**< branch-and-bound tree */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(tree != NULL);
   assert(mem != NULL);
   assert(set != NULL);

   if( num > tree->siblingssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(mem->treemem, tree->siblings, tree->siblingssize, newsize) );
      tree->siblingssize = newsize;
   }
   assert(num <= tree->siblingssize);

   return SCIP_OKAY;
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
JUNCTION* junctionCreate(               /**< creates junction data */
   MEM*             mem,                /**< block memory buffers */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   JUNCTION* junction;

   assert(mem != NULL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);

   ALLOC_NULL( allocBlockMemory(mem->treemem, junction) );
   
   junction->nchildren = tree->nchildren;

   return junction;
}

static
RETCODE junctionFree(                   /**< frees junction data */
   JUNCTION**       junction,           /**< junction data */
   MEM*             mem                 /**< block memory buffers */
   )
{
   assert(junction != NULL);
   assert(*junction != NULL);
   assert((*junction)->nchildren == 0);

   freeBlockMemory(mem->treemem, *junction);

   return SCIP_OKAY;
}

static
FORK* forkCreate(                       /**< creates fork data */
   MEM*             mem,                /**< block memory buffers */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   FORK* fork;

   assert(mem != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(tree != NULL);
   assert(tree->nchildren > 0);

   ALLOC_NULL( allocBlockMemory(mem->treemem, fork) );

   CHECK_NULL( SCIPlpiGetState(lp->lpi, mem, &(fork->lpstate)) );
   fork->nlpstateref = 0;
   fork->addedCols = NULL;
   fork->addedRows = NULL;
   fork->naddedCols = SCIPlpGetNumNewcols(lp);
   fork->naddedRows = SCIPlpGetNumNewrows(lp);
   fork->nchildren = tree->nchildren;

   if( fork->naddedCols > 0 )
   {
      /* copy the newly created columns to the fork's col array */
      ALLOC_NULL( duplicateBlockMemoryArray(mem->treemem, fork->addedCols, SCIPlpGetNewcols(lp), fork->naddedCols) );
   }
   if( fork->naddedRows > 0 )
   {
      /* copy the newly created rows to the fork's row array */
      ALLOC_NULL( duplicateBlockMemoryArray(mem->treemem, fork->addedRows, SCIPlpGetNewrows(lp), fork->naddedRows) );
   }
   
   return fork;
}

static
RETCODE forkFree(                       /**< frees fork data */
   FORK**           fork,               /**< fork data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Bool             releaseData         /**< should the columns and rows of the node's arrays be released? */
   )
{
   assert(fork != NULL);
   assert(*fork != NULL);
   assert((*fork)->nchildren == 0);
   assert((*fork)->nlpstateref == 0);
   assert((*fork)->lpstate == NULL);
   assert(mem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( releaseData )
   {
      int i;
      
      for( i = 0; i < (*fork)->naddedCols; ++i )
         SCIPcolRelease(&((*fork)->addedCols[i]), mem, set, lp);
      for( i = 0; i < (*fork)->naddedRows; ++i )
         SCIProwRelease(&((*fork)->addedRows[i]), mem, set, lp);
   }

   freeBlockMemoryArrayNull(mem->treemem, (*fork)->addedCols, (*fork)->naddedCols);
   freeBlockMemoryArrayNull(mem->treemem, (*fork)->addedRows, (*fork)->naddedRows);
   freeBlockMemory(mem->treemem, *fork);

   return SCIP_OKAY;
}

static
RETCODE forkReleaseLPState(             /**< decreases the reference counter of the LP state in the fork */
   FORK*            fork,               /**< fork data */
   MEM*             mem,                /**< block memory buffers */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(fork != NULL);
   assert(fork->nlpstateref > 0);
   assert(fork->lpstate != NULL);
   assert(mem != NULL);
   assert(lp != NULL);

   fork->nlpstateref--;
   if( fork->nlpstateref == 0 )
      CHECK_OKAY( SCIPlpiFreeState(lp->lpi, mem, &(fork->lpstate)) );

   return SCIP_OKAY;
}

static
SUBROOT* subrootCreate(                 /**< creates subroot data */
   MEM*             mem,                /**< block memory buffers */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   SUBROOT* subroot;

   assert(mem != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(tree != NULL);
   assert(tree->nchildren > 0);

   ALLOC_NULL( allocBlockMemory(mem->treemem, subroot) );

   CHECK_NULL( SCIPlpiGetState(lp->lpi, mem, &(subroot->lpstate)) );
   subroot->nlpstateref = 0;
   subroot->ncols = lp->ncols;
   subroot->nrows = lp->nrows;
   subroot->nchildren = tree->nchildren;
   ALLOC_NULL( duplicateBlockMemoryArray(mem->treemem, subroot->cols, lp->cols, subroot->ncols) );
   ALLOC_NULL( duplicateBlockMemoryArray(mem->treemem, subroot->rows, lp->rows, subroot->nrows) );

   return subroot;
}

static
RETCODE subrootFree(                    /**< frees subroot */
   SUBROOT**        subroot,            /**< subroot data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Bool             releaseData         /**< should the columns and rows of the node's arrays be released? */
   )
{
   assert(subroot != NULL);
   assert(*subroot != NULL);
   assert((*subroot)->nchildren == 0);
   assert((*subroot)->nlpstateref == 0);
   assert(mem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   CHECK_OKAY( SCIPlpiFreeState(lp->lpi, mem, &((*subroot)->lpstate)) );

   if( releaseData )
   {
      int i;
      
      for( i = 0; i < (*subroot)->ncols; ++i )
         SCIPcolRelease(&((*subroot)->cols[i]), mem, set, lp);
      for( i = 0; i < (*subroot)->nrows; ++i )
         SCIProwRelease(&((*subroot)->rows[i]), mem, set, lp);
   }

   freeBlockMemoryArrayNull(mem->treemem, (*subroot)->cols, (*subroot)->ncols);
   freeBlockMemoryArrayNull(mem->treemem, (*subroot)->rows, (*subroot)->nrows);
   freeBlockMemory(mem->treemem, *subroot);

   return SCIP_OKAY;
}

static
RETCODE subrootReleaseLPState(          /**< decreases the reference counter of the LP state in the subroot */
   SUBROOT*         subroot,            /**< subroot data */
   MEM*             mem,                /**< block memory buffers */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(subroot != NULL);
   assert(subroot->nlpstateref > 0);
   assert(subroot->lpstate != NULL);
   assert(mem != NULL);
   assert(lp != NULL);

   subroot->nlpstateref--;
   if( subroot->nlpstateref == 0 )
      CHECK_OKAY( SCIPlpiFreeState(lp->lpi, mem, &(subroot->lpstate)) );
   
   return SCIP_OKAY;
}

static
RETCODE nodeAssignParent(               /**< makes node a child of the given parent node, which must be the active node */
   NODE*            node,               /**< child node */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   NODE*            parent              /**< parent (= active) node (or NULL, if node is root) */
   )
{
   assert(node != NULL);
   assert(node->parent == NULL);
   assert(mem != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(parent != NULL);
   assert(parent->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(parent->conslist != NULL); /* we need at least one constraint! */

   /* link node to parent */
   node->parent = parent;
   node->conslist = parent->conslist;
   node->lowerbound = parent->lowerbound;
   node->depth = parent->depth+1;

   /* register node in the childlist of the active (the parent) node */
   CHECK_OKAY( treeEnsureChildrenMem(tree, mem, set, tree->nchildren+1) );
   tree->children[tree->nchildren] = node;
   tree->nchildren++;

   return SCIP_OKAY;
}

static
void nodeReleaseParent(                 /**< decreases number of children of the parent, frees it if no children left */
   NODE*            node,               /**< child node */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   NODE* parent;
   Bool hasChildren = TRUE;

   assert(mem != NULL);
   assert(node != NULL);
   
   parent = node->parent;
   if( parent != NULL )
   {
      switch( parent->nodetype )
      {
      case SCIP_NODETYPE_LEAF:
         errorMessage("Leaf cannot be a parent node");
         abort();
      case SCIP_NODETYPE_ACTNODE:
         errorMessage("Cannot release the parent-child relationship, if parent is the active node");
         abort();
      case SCIP_NODETYPE_JUNCTION:
         assert(parent->data.junction != NULL);
         assert(parent->data.junction->nchildren > 0);
         parent->data.junction->nchildren--;
         hasChildren = (parent->data.junction->nchildren > 0);
         break;
      case SCIP_NODETYPE_FORK:
         assert(parent->data.fork != NULL);
         assert(parent->data.fork->nchildren > 0);
         parent->data.fork->nchildren--;
         hasChildren = (parent->data.fork->nchildren > 0);
         break;
      case SCIP_NODETYPE_SUBROOT:
         assert(parent->data.subroot != NULL);
         assert(parent->data.subroot->nchildren > 0);
         parent->data.subroot->nchildren--;
         hasChildren = (parent->data.subroot->nchildren > 0);
         break;
      default:
         errorMessage("Unknown node type");
         abort();
      }

      /* free parent, if it has no more children left and is not on the actual active path */
      if( !hasChildren && !parent->active )
         SCIPnodeFree(&node->parent, mem, set, lp);
   }
}

NODE* SCIPnodeCreate(                   /**< creates a child node of the active node */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   NODE* node;

   assert(mem != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->pathlen == 0 || tree->path != NULL);

   ALLOC_NULL( allocBlockMemory(mem->treemem, node) );
   node->parent = NULL;
   node->conslist = NULL;
   node->domchg = NULL;
   node->lowerbound = -SCIPinfinity(set);
   node->depth = 0;
   node->nodetype = SCIP_NODETYPE_LEAF;
   node->active = FALSE;
   
   if( tree->pathlen > 0 )
   {
      assert(tree->path[tree->pathlen-1]->nodetype == SCIP_NODETYPE_ACTNODE);
      CHECK_NULL( nodeAssignParent(node, mem, set, tree, tree->path[tree->pathlen-1]) );
   }

   return node;
}

RETCODE SCIPnodeFree(                   /**< frees node */
   NODE**           node,               /**< node data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(mem != NULL);
   assert(node != NULL);
   assert(*node != NULL);
   assert(!(*node)->active);

   /* free nodetype specific data */
   switch((*node)->nodetype)
   {
   case SCIP_NODETYPE_LEAF:
      break;
   case SCIP_NODETYPE_ACTNODE:
      break;
   case SCIP_NODETYPE_JUNCTION:
      CHECK_OKAY( junctionFree(&((*node)->data.junction), mem) );
      break;
   case SCIP_NODETYPE_FORK:
      CHECK_OKAY( forkFree(&((*node)->data.fork), mem, set, lp, TRUE) );
      break;
   case SCIP_NODETYPE_SUBROOT:
      CHECK_OKAY( subrootFree(&((*node)->data.subroot), mem, set, lp, TRUE) );
      break;
   default:
      errorMessage("Unknown node type");
      break;
   }

   /* free common data */
   SCIPconslistFreePart(&((*node)->conslist), mem, (*node)->parent->conslist);
   SCIPdomchgFree(&((*node)->domchg), mem);
   nodeReleaseParent(*node, mem, set, lp);

   freeBlockMemory(mem->treemem, *node);

   return SCIP_OKAY;
}

static
RETCODE nodeDeactivate(                 /**< informs node, that it is no longer on the active path */
   NODE**           node,               /**< node to deactivate */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   Bool hasChildren = TRUE;

   assert(node != NULL);
   assert(*node != NULL);
   assert((*node)->active);
   assert(mem != NULL);
   assert(set != NULL);
   assert(tree != NULL);

   (*node)->active = FALSE;

   switch( (*node)->nodetype )
   {
   case SCIP_NODETYPE_LEAF:
      errorMessage("Cannot deactivate leaf (which shouldn't be active)");
      abort();
   case SCIP_NODETYPE_ACTNODE:
      if( tree->nchildren > 0 )
      {
         errorMessage("Cannot deactivate active node with children");
         abort();
      }
      hasChildren = FALSE;
      break;
   case SCIP_NODETYPE_JUNCTION:
      assert((*node)->data.junction != NULL);
      assert((*node)->data.junction->nchildren > 0);
      hasChildren = ((*node)->data.junction->nchildren > 0);
      break;
   case SCIP_NODETYPE_FORK:
      assert((*node)->data.fork != NULL);
      assert((*node)->data.fork->nchildren > 0);
      hasChildren = ((*node)->data.fork->nchildren > 0);
      break;
   case SCIP_NODETYPE_SUBROOT:
      assert((*node)->data.subroot != NULL);
      assert((*node)->data.subroot->nchildren > 0);
      hasChildren = ((*node)->data.subroot->nchildren > 0);
      break;
   default:
      errorMessage("Unknown node type");
      abort();
   }

   /* free node, if it has no children */
   if( !hasChildren )
   {
      CHECK_OKAY( SCIPnodeFree(node, mem, set, lp) );
   }

   return SCIP_OKAY;
}



/*
 * Path Switching
 */

static
void treeUpdatePathLPSize(              /**< updates the LP sizes of the active path starting at the given depth */
   TREE*            tree,               /**< branch-and-bound tree */
   int              startdepth          /**< depth to start counting */
   )
{
   NODE* node;
   int ncols;
   int nrows;
   int i;

   assert(tree != NULL);
   assert(startdepth < tree->pathlen);

   if( startdepth == 0 )
   {
      ncols = 0;
      nrows = 0;
   }
   else
   {
      ncols = tree->pathnlpcols[startdepth-1];
      nrows = tree->pathnlprows[startdepth-1];
   }

   for( i = startdepth; i < tree->pathlen; ++i )
   {
      node = tree->path[i];
      assert(node != NULL);
      assert(node->active);
      assert(node->depth == i);
      
      switch( node->nodetype )
      {
      case SCIP_NODETYPE_LEAF:
         errorMessage("Leaf cannot be in the active path");
         abort();
      case SCIP_NODETYPE_ACTNODE:
         assert(i == tree->pathlen-1);
         break;
      case SCIP_NODETYPE_JUNCTION:
         assert(node->data.junction != NULL);
         break;
      case SCIP_NODETYPE_FORK:
         assert(node->data.fork != NULL);
         ncols += node->data.fork->naddedCols;
         nrows += node->data.fork->naddedRows;
         break;
      case SCIP_NODETYPE_SUBROOT:
         assert(node->data.subroot != NULL);
         ncols = node->data.subroot->ncols;
         nrows = node->data.subroot->nrows;
         break;
      default:
         errorMessage("Unknown node type");
         abort();
      }
      tree->pathnlpcols[i] = ncols;
      tree->pathnlprows[i] = nrows;
   }
}

static
RETCODE treeShrinkPath(                 /**< cuts off path of active nodes after given node, marks cutted nodes inactive,
                                           and undoes their domain changes */
   TREE*            tree,               /**< branch-and-bound tree */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   int              lastdepth           /**< depth of the last node in the shrinked path */
   )
{
   int i;

   assert(tree != NULL);
   assert(lastdepth >= -1);
   assert(lastdepth < tree->pathlen);

   for( i = tree->pathlen-1; i > lastdepth; --i )
   {
      assert(tree->path[i] != NULL);
      assert(tree->path[i]->depth == i);
      CHECK_OKAY( nodeDeactivate(&(tree->path[i]), mem, set, lp, tree) );
   }
   tree->pathlen = lastdepth+1;

   assert(tree->pathlen <= tree->pathsize);

   return SCIP_OKAY;
}

static
RETCODE treeSwitchPath(                 /**< switches the active path to end at the given node, applies domain changes */
   TREE*            tree,               /**< branch-and-bound tree */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   NODE*            node                /**< last node of the new active path (= new active node) */
   )
{
   NODE* commonfork;    /* common fork node */
   NODE* lpfork;        /* fork node defining the LP state of the new active node */
   NODE* subroot;       /* subroot of new active path */
   int commonforkdepth; /* depth of the common fork node */
   int i;

   assert(tree != NULL);
   assert(mem != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(node != NULL);
   assert(!node->active);
   assert(node->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(tree->actLPFork == NULL || tree->actLPFork->active);
   assert(tree->actSubroot == NULL || tree->actSubroot->active);

   /* find the common fork node, the new LP defining fork, and the new active subroot */
   commonfork = node;
   lpfork = NULL;
   subroot = NULL;
   while( commonfork != NULL && !commonfork->active )
   {
      if( lpfork == NULL && 
         (commonfork->nodetype == SCIP_NODETYPE_SUBROOT || commonfork->nodetype == SCIP_NODETYPE_SUBROOT) )
         lpfork = commonfork;
      if( subroot == NULL && commonfork->nodetype == SCIP_NODETYPE_SUBROOT )
         subroot = commonfork;
      commonfork = commonfork->parent;
   }
   if( commonfork == NULL )
      commonforkdepth = -1;
   else
      commonforkdepth = commonfork->depth;
   assert(lpfork == NULL || !lpfork->active);
   assert(subroot == NULL || !subroot->active);

   /* if not already found, continue searching the LP defining fork */
   if( lpfork == NULL )
   {
      if( tree->actLPFork != NULL && tree->actLPFork->depth >= commonforkdepth )
      {
         lpfork = commonfork;
         while( lpfork != NULL &&
            (lpfork->nodetype == SCIP_NODETYPE_SUBROOT || lpfork->nodetype == SCIP_NODETYPE_SUBROOT) )
            lpfork = lpfork->parent;
      }
      else
         lpfork = tree->actLPFork;
   }

   /* if not already found, continue searching the subroot */
   if( subroot == NULL )
   {
      if( tree->actSubroot != NULL && tree->actSubroot->depth >= commonforkdepth )
      {
         subroot = commonfork;
         while( subroot != NULL && subroot->nodetype != SCIP_NODETYPE_SUBROOT )
            subroot = subroot->parent;
      }
      else
         subroot = tree->actSubroot;
   }
   
   /* remember the depth of the common fork node for LP updates */
   if( subroot == tree->actSubroot )
   {
      /* we are in the same subtree: the LP is correct at most upto the fork depth */
      assert(subroot == NULL || subroot->active);
      tree->correctLPDepth = MIN(tree->correctLPDepth, commonforkdepth);
   }
   else
   {
      /* we are in a different subtree: the LP is completely incorrect */
      assert(subroot == NULL || !subroot->active);
      tree->correctLPDepth = -1;
   }
   
   /* undo the domain changes of the old active path */
   for( i = tree->pathlen-1; i > commonforkdepth; --i )
      CHECK_OKAY( SCIPlpUndoDomchg(lp, mem, set, tree->path[i]->domchg) );

   /* shrink active path to the common fork and deactivate the corresponding nodes */
   CHECK_OKAY( treeShrinkPath(tree, mem, set, lp, commonforkdepth) );
   assert(tree->pathlen == commonforkdepth+1);

   /* create the new active path */
   CHECK_OKAY( treeEnsurePathMem(tree, mem, set, node->depth+1) );
   tree->pathlen = node->depth+1;
   while( node != commonfork )
   {
      assert(node != NULL);
      tree->path[node->depth] = node;
      node->active = TRUE;
      node = node->parent;
   }

   /* count the new LP sizes of the path */
   treeUpdatePathLPSize(tree, commonforkdepth+1);

   /* apply domain changes of the new path */
   for( i = commonforkdepth+1; i < tree->pathlen; ++i )
      CHECK_OKAY( SCIPlpApplyDomchg(lp, mem, set, tree->path[i]->domchg) );

   /* remember LP defining fork and subroot */
   assert(subroot == NULL || lpfork != NULL);
   assert(subroot == NULL || subroot->depth <= lpfork->depth);
   tree->actLPFork = lpfork;
   tree->actSubroot = subroot;
   
   return SCIP_OKAY;
}

static
RETCODE subrootConstructLP(             /**< loads the subroot's LP data */
   NODE*            subroot,            /**< subroot node to construct LP for */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   COL** cols;
   ROW** rows;
   int ncols;
   int nrows;
   int c;
   int r;

   assert(subroot != NULL);
   assert(subroot->nodetype == SCIP_NODETYPE_SUBROOT);
   assert(subroot->data.subroot != NULL);
   assert(mem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   cols = subroot->data.subroot->cols;
   rows = subroot->data.subroot->rows;
   ncols = subroot->data.subroot->ncols;
   nrows = subroot->data.subroot->nrows;

   assert(ncols == 0 || cols != NULL);
   assert(nrows == 0 || rows != NULL);
   
   for( c = 0; c < ncols; ++c )
      CHECK_OKAY( SCIPlpAddCol(lp, mem, set, cols[c]) );
   for( r = 0; r < nrows; ++r )
      CHECK_OKAY( SCIPlpAddRow(lp, mem, set, rows[r]) );

   return SCIP_OKAY;
}
   
static
RETCODE forkAddLP(                      /**< loads the fork's additional LP data */
   NODE*            fork,               /**< fork node to construct additional LP for */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   COL** cols;
   ROW** rows;
   int ncols;
   int nrows;
   int c;
   int r;

   assert(fork != NULL);
   assert(fork->nodetype == SCIP_NODETYPE_FORK);
   assert(fork->data.fork != NULL);
   assert(mem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   cols = fork->data.fork->addedCols;
   rows = fork->data.fork->addedRows;
   ncols = fork->data.fork->naddedCols;
   nrows = fork->data.fork->naddedRows;

   assert(ncols == 0 || cols != NULL);
   assert(nrows == 0 || rows != NULL);
   
   for( c = 0; c < ncols; ++c )
      CHECK_OKAY( SCIPlpAddCol(lp, mem, set, cols[c]) );
   for( r = 0; r < nrows; ++r )
      CHECK_OKAY( SCIPlpAddRow(lp, mem, set, rows[r]) );

   return SCIP_OKAY;
}
   
RETCODE SCIPtreeLoadLP(                 /**< constructs the LP and loads LP state for fork/subroot of the active node */
   TREE*            tree,               /**< branch-and-bound tree */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   NODE* lpfork;
   NODE* pathnode;
   int lpforkdepth;
   int d;

   assert(mem != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->path != NULL);
   assert(tree->pathlen > 0);
   assert(tree->path[tree->pathlen-1]->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(lp != NULL);

   lpfork = tree->actLPFork;

   /* find out the lpfork's depth (or -1, if lpfork is NULL) */
   if( lpfork == NULL )
   {
      assert(tree->correctLPDepth == -1);
      assert(tree->actSubroot == NULL);
      lpforkdepth = -1;
   }
   else
   {
      assert(lpfork->nodetype == SCIP_NODETYPE_FORK || lpfork->nodetype == SCIP_NODETYPE_SUBROOT);
      assert(lpfork->active);
      assert(tree->path[lpfork->depth] == lpfork);
      lpforkdepth = lpfork->depth;
   }
   assert(lpforkdepth < tree->pathlen-1); /* lpfork must not be the last (the active) node of the active path */

   /* find out, if we are in the same subtree */
   if( tree->correctLPDepth >= 0 )
   {
      /* same subtree: shrink LP to the deepest node with correct LP */
      assert(lpfork != NULL);
      assert(tree->correctLPDepth <= lpforkdepth);
      CHECK_OKAY( SCIPlpShrinkCols(lp, tree->pathnlpcols[tree->correctLPDepth]) );
      CHECK_OKAY( SCIPlpShrinkRows(lp, tree->pathnlprows[tree->correctLPDepth]) );
   }
   else
   {
      /* other subtree: fill LP with the subroot LP data */
      CHECK_OKAY( SCIPlpClear(lp) );
      if( tree->actSubroot != NULL )
      {
         CHECK_OKAY( subrootConstructLP(tree->actSubroot, mem, set, lp) );
         tree->correctLPDepth = tree->actSubroot->depth; 
      }
   }

   assert(lpforkdepth < tree->pathlen);
   assert(tree->correctLPDepth <= lpforkdepth);

   /* add the missing columns and rows */
   for( d = tree->correctLPDepth+1; d <= lpforkdepth; ++d )
   {
      pathnode = tree->path[d];
      assert(pathnode != NULL);
      assert(pathnode->depth == d);
      assert(pathnode->nodetype == SCIP_NODETYPE_JUNCTION || pathnode->nodetype == SCIP_NODETYPE_FORK);
      if( pathnode->nodetype == SCIP_NODETYPE_FORK )
         CHECK_OKAY( forkAddLP(pathnode, mem, set, lp) );
   }
   tree->correctLPDepth = lpforkdepth;

   /* load LP state, if existing */
   if( lpfork != NULL )
   {
      LPSTATE* lpstate;

      if( lpfork->nodetype == SCIP_NODETYPE_FORK )
      {
         assert(lpfork->data.fork != NULL);
         lpstate = lpfork->data.fork->lpstate;
      }
      else
      {
         assert(lpfork->nodetype == SCIP_NODETYPE_SUBROOT);
         assert(lpfork->data.subroot != NULL);
         lpstate = lpfork->data.subroot->lpstate;
      }
      assert(lpstate != NULL);
      CHECK_OKAY( SCIPlpSetState(lp, mem, lpstate) );
   }
   
   return SCIP_OKAY;
}




/*
 * Node Conversion
 */

RETCODE SCIPnodeActivate(               /**< activates a leaf node */
   NODE*            node,               /**< leaf node to activate */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_LEAF);
   assert(!node->active);
   assert(mem != NULL);
   assert(tree != NULL);
   assert(lp != NULL);

   /* convert node into the actual node */
   node->nodetype = SCIP_NODETYPE_ACTNODE;

   /* track the path from the old active node to the new node, and perform domain changes */
   CHECK_OKAY( treeSwitchPath(tree, mem, set, lp, node) );
   assert(tree->pathlen > 0);
   assert(tree->path[tree->pathlen-1] == node);

   /* remember domain changes of the active node in the dynamic data */
   CHECK_OKAY( SCIPdomchgdynCopy(tree->domchgdyn, mem, set, node->domchg) );

   return SCIP_OKAY;
}   

RETCODE SCIPnodeToJunction(             /**< converts the active node into a junction node */
   NODE*            node,               /**< node to convert */
   MEM*             mem,                /**< block memory buffers */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   JUNCTION* junction;

   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(node->active);
   assert(tree != NULL);

   ALLOC_OKAY( junction = junctionCreate(mem, tree) );

   node->nodetype = SCIP_NODETYPE_JUNCTION;
   node->data.junction = junction;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, node->depth);

   return SCIP_OKAY;
}

RETCODE SCIPnodeToFork(                 /**< converts the active node into a fork node */
   NODE*            node,               /**< node to convert */
   MEM*             mem,                /**< block memory buffers */
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
   assert(lp->flushed);
   assert(lp->solved);

   ALLOC_OKAY( fork = forkCreate(mem, lp, tree) );
   
   node->nodetype = SCIP_NODETYPE_FORK;
   node->data.fork = fork;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, node->depth);

   return SCIP_OKAY;
}

RETCODE SCIPnodeToSubroot(              /**< converts the active node into a subroot node */
   NODE*            node,               /**< node to convert */
   MEM*             mem,                /**< block memory buffers */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   SUBROOT* subroot;

   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(node->active);
   assert(mem != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   ALLOC_OKAY( subroot = subrootCreate(mem, lp, tree) );

   node->nodetype = SCIP_NODETYPE_SUBROOT;
   node->data.subroot = subroot;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, node->depth);

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
   ALLOC_NULL( allocBlockMemoryArray(mem->treemem, tree->path, tree->pathsize) );
   ALLOC_NULL( allocBlockMemoryArray(mem->treemem, tree->pathnlpcols, tree->pathsize) );
   ALLOC_NULL( allocBlockMemoryArray(mem->treemem, tree->pathnlprows, tree->pathsize) );
   CHECK_NULL( SCIPpqueueInit(&(tree->leaves), set->treeGrowInit, set->treeGrowFac, set->nodecmp) );
   ALLOC_NULL( tree->domchgdyn = SCIPdomchgdynCreate(mem) );

   tree->root = NULL;
   tree->pathlen = 0;
   tree->pathsize = SCIPcalcPathGrowSize(set, 0);
   tree->actLPFork = NULL;
   tree->actSubroot = NULL;
   tree->correctLPDepth = -1;

   return tree;
}

































#if 0 // ???
   NODE* forknode;    /* common fork node */
   int forkdepth;     /* depth of the common fork node */
   int lpdepth;       /* depth to which LP is setup correctly */

   assert(node != NULL);
   assert(node->active);
   assert(node->nodetype == SCIP_NODETYPE_LPFORK || node->nodetype == SCIP_NODETYPE_SUBROOT);
   assert(mem != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(tree != NULL);

   forkdepth = -1;
   lpdepth = -1;
   correctcols = 0;
   correctrows = 0;
   colsshrinked = FALSE;
   rowsshrinked = FALSE;

   /* check, whether we are in the same subtree as the LPI defining node */
   if( tree->actSubroot == lp->lpiSubroot )
   {
      /* switch LPs of nodes from the same subtree */

      /* find the common fork node of LPI defining node and actual node */
      forknode = lp->actLPNode;
      assert(forknode == NULL || forknode->isActLPNode);
      while( forknode != NULL && !forknode->active )
         forknode = forknode->parent;
      if( forknode != NULL )
      {
         forkdepth = forknode->depth;

         assert(tree->path[forkdepth] == forknode);
         assert(0 <= tree->pathnlpcols[forkdepth] && tree->pathnlpcols[forkdepth] <= lp->ncols);
         assert(0 <= tree->pathnlprows[forkdepth] && tree->pathnlprows[forkdepth] <= lp->nrows);
         
         lpdepth = forkdepth;
         correctcols = tree->pathnlpcols[forkdepth];
         correctrows = tree->pathnlprows[forkdepth];
      }
   }
   else if( tree->actSubroot != NULL )
   {
      COL** cols;
      ROW** rows;
      int ncols;
      int nrows;

      /* switch LPs of nodes from different subtrees, where subroot of actual node is available */
      assert(tree->actSubroot->nodetype == SCIP_NODETYPE_SUBROOT);
      assert(tree->actSubroot->data.subroot != NULL);

      cols = tree->actSubroot->data.subroot->cols;
      rows = tree->actSubroot->data.subroot->rows;
      ncols = tree->actSubroot->data.subroot->ncols;
      nrows = tree->actSubroot->data.subroot->nrows;
      
      /* keep columns and rows, which wouldn't change position in new LP */
      for( correctcols = 0; correctcols < ncols && cols[correctcols]->lpipos == correctcols; ++correctcols )
      {
         assert(correctcols < lp->ncols);
      }
      for( correctrows = 0; correctrows < nrows && rows[correctrows]->lpipos == correctrows; ++correctrows )
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
      
      lpdepth = tree->actSubroot->depth;
   }
   else
   {
      /* switch LPs of nodes from different subtrees, where subroot of actual node is not available */
   }

   /* construct the new LP */
   if( samesubtree )
   {
      /* same subtree */
   }
   else if( subroot != NULL )
   {
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
   
}
#endif // ???


#if 0 // ???
static
RETCODE switchPath(                     /**< performs LP changes on the path from oldnode to node */
   NODE*            node,               /**< lpfork/fork/subtree node to construct LP for */
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
#endif // ???


#if 0 // ????????????
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
#endif // ???????????


