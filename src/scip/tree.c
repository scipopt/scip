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
   CONSLIST*        conslist;           /**< list of constraints created at this node */
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
   NODE*            actNode;            /**< active node */
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
 * dynamic memory arrays
 */

static
RETCODE treeEnsureChildrenMem(          /**< resizes children array to be able to store at least num nodes */
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(tree != NULL);
   assert(set != NULL);

   if( num > tree->childrensize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(tree->children, newsize) );
      tree->childrensize = newsize;
   }
   assert(num <= tree->childrensize);

   return SCIP_OKAY;
}

static
RETCODE treeEnsureSiblingsMem(          /**< resizes siblings array to be able to store at least num nodes */
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(tree != NULL);
   assert(set != NULL);

   if( num > tree->siblingssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(tree->siblings, newsize) );
      tree->siblingssize = newsize;
   }
   assert(num <= tree->siblingssize);

   return SCIP_OKAY;
}

static
RETCODE treeEnsurePathMem(              /**< resizes path array to be able to store at least num nodes */
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in path */
   )
{
   assert(tree != NULL);
   assert(set != NULL);

   if( num > tree->pathsize )
   {
      int newsize;

      newsize = SCIPcalcPathGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(tree->path, newsize) );
      ALLOC_OKAY( reallocMemoryArray(tree->pathnlpcols, newsize) );
      ALLOC_OKAY( reallocMemoryArray(tree->pathnlprows, newsize) );
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
RETCODE junctionCreate(                 /**< creates junction data */
   JUNCTION**       junction,           /**< pointer to junction data */
   MEMHDR*          memhdr,             /**< block memory */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(junction != NULL);
   assert(memhdr != NULL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);

   ALLOC_OKAY( allocBlockMemory(memhdr, *junction) );
   
   (*junction)->nchildren = tree->nchildren;

   return SCIP_OKAY;
}

static
RETCODE junctionFree(                   /**< frees junction data */
   JUNCTION**       junction,           /**< junction data */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(junction != NULL);
   assert(*junction != NULL);
   assert((*junction)->nchildren == 0);

   freeBlockMemory(memhdr, *junction);

   return SCIP_OKAY;
}

static
RETCODE forkCreate(                     /**< creates fork data */
   FORK**           fork,               /**< pointer to fork data */
   MEMHDR*          memhdr,             /**< block memory */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(fork != NULL);
   assert(memhdr != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(tree != NULL);
   assert(tree->nchildren > 0);

   ALLOC_OKAY( allocBlockMemory(memhdr, *fork) );

   CHECK_OKAY( SCIPlpiGetState(lp->lpi, memhdr, &((*fork)->lpstate)) );
   (*fork)->nlpstateref = 0;
   (*fork)->addedCols = NULL;
   (*fork)->addedRows = NULL;
   (*fork)->naddedCols = SCIPlpGetNumNewcols(lp);
   (*fork)->naddedRows = SCIPlpGetNumNewrows(lp);
   (*fork)->nchildren = tree->nchildren;

   if( (*fork)->naddedCols > 0 )
   {
      /* copy the newly created columns to the fork's col array */
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*fork)->addedCols, SCIPlpGetNewcols(lp), (*fork)->naddedCols) );
   }
   if( (*fork)->naddedRows > 0 )
   {
      /* copy the newly created rows to the fork's row array */
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*fork)->addedRows, SCIPlpGetNewrows(lp), (*fork)->naddedRows) );
   }
   
   return SCIP_OKAY;
}

static
RETCODE forkFree(                       /**< frees fork data */
   FORK**           fork,               /**< fork data */
   MEMHDR*          memhdr,             /**< block memory */
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
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( releaseData )
   {
      int i;
      
      for( i = 0; i < (*fork)->naddedCols; ++i )
         SCIPcolRelease(&((*fork)->addedCols[i]), memhdr, set, lp);
      for( i = 0; i < (*fork)->naddedRows; ++i )
         SCIProwRelease(&((*fork)->addedRows[i]), memhdr, set, lp);
   }

   freeBlockMemoryArrayNull(memhdr, (*fork)->addedCols, (*fork)->naddedCols);
   freeBlockMemoryArrayNull(memhdr, (*fork)->addedRows, (*fork)->naddedRows);
   freeBlockMemory(memhdr, *fork);

   return SCIP_OKAY;
}

static
RETCODE forkReleaseLPState(             /**< decreases the reference counter of the LP state in the fork */
   FORK*            fork,               /**< fork data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(fork != NULL);
   assert(fork->nlpstateref > 0);
   assert(fork->lpstate != NULL);
   assert(memhdr != NULL);
   assert(lp != NULL);

   fork->nlpstateref--;
   if( fork->nlpstateref == 0 )
      CHECK_OKAY( SCIPlpiFreeState(lp->lpi, memhdr, &(fork->lpstate)) );

   return SCIP_OKAY;
}

static
RETCODE subrootCreate(                  /**< creates subroot data */
   SUBROOT**        subroot,            /**< pointer to subroot data */
   MEMHDR*          memhdr,             /**< block memory */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(subroot != NULL);
   assert(memhdr != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(tree != NULL);
   assert(tree->nchildren > 0);

   ALLOC_OKAY( allocBlockMemory(memhdr, *subroot) );

   CHECK_OKAY( SCIPlpiGetState(lp->lpi, memhdr, &((*subroot)->lpstate)) );
   (*subroot)->nlpstateref = 0;
   (*subroot)->ncols = lp->ncols;
   (*subroot)->nrows = lp->nrows;
   (*subroot)->nchildren = tree->nchildren;
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*subroot)->cols, lp->cols, (*subroot)->ncols) );
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*subroot)->rows, lp->rows, (*subroot)->nrows) );

   return SCIP_OKAY;
}

static
RETCODE subrootFree(                    /**< frees subroot */
   SUBROOT**        subroot,            /**< subroot data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Bool             releaseData         /**< should the columns and rows of the node's arrays be released? */
   )
{
   assert(subroot != NULL);
   assert(*subroot != NULL);
   assert((*subroot)->nchildren == 0);
   assert((*subroot)->nlpstateref == 0);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   CHECK_OKAY( SCIPlpiFreeState(lp->lpi, memhdr, &((*subroot)->lpstate)) );

   if( releaseData )
   {
      int i;
      
      for( i = 0; i < (*subroot)->ncols; ++i )
         SCIPcolRelease(&((*subroot)->cols[i]), memhdr, set, lp);
      for( i = 0; i < (*subroot)->nrows; ++i )
         SCIProwRelease(&((*subroot)->rows[i]), memhdr, set, lp);
   }

   freeBlockMemoryArrayNull(memhdr, (*subroot)->cols, (*subroot)->ncols);
   freeBlockMemoryArrayNull(memhdr, (*subroot)->rows, (*subroot)->nrows);
   freeBlockMemory(memhdr, *subroot);

   return SCIP_OKAY;
}

static
RETCODE subrootReleaseLPState(          /**< decreases the reference counter of the LP state in the subroot */
   SUBROOT*         subroot,            /**< subroot data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(subroot != NULL);
   assert(subroot->nlpstateref > 0);
   assert(subroot->lpstate != NULL);
   assert(memhdr != NULL);
   assert(lp != NULL);

   subroot->nlpstateref--;
   if( subroot->nlpstateref == 0 )
      CHECK_OKAY( SCIPlpiFreeState(lp->lpi, memhdr, &(subroot->lpstate)) );
   
   return SCIP_OKAY;
}

static
RETCODE nodeAssignParent(               /**< makes node a child of the given parent node, which must be the active node */
   NODE*            node,               /**< child node */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   NODE*            parent              /**< parent (= active) node (or NULL, if node is root) */
   )
{
   assert(node != NULL);
   assert(node->parent == NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(parent != NULL);
   assert(parent->nodetype == SCIP_NODETYPE_ACTNODE);

   /* link node to parent */
   node->parent = parent;
   node->lowerbound = parent->lowerbound;
   node->depth = parent->depth+1;

   /* register node in the childlist of the active (the parent) node */
   CHECK_OKAY( treeEnsureChildrenMem(tree, set, tree->nchildren+1) );
   tree->children[tree->nchildren] = node;
   tree->nchildren++;

   return SCIP_OKAY;
}

static
void nodeReleaseParent(                 /**< decreases number of children of the parent, frees it if no children left */
   NODE*            node,               /**< child node */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   NODE* parent;
   Bool hasChildren = TRUE;

   assert(memhdr != NULL);
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
         SCIPnodeFree(&node->parent, memhdr, set, lp);
   }
}

RETCODE SCIPnodeCreate(                 /**< creates a child node of the active node */
   NODE**           node,               /**< pointer to node data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(node != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->pathlen == 0 || tree->path != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, *node) );
   (*node)->parent = NULL;
   (*node)->conslist = NULL;
   (*node)->domchg = NULL;
   (*node)->lowerbound = -SCIPinfinity(set);
   (*node)->depth = 0;
   (*node)->nodetype = SCIP_NODETYPE_LEAF;
   (*node)->active = FALSE;
   
   if( tree->pathlen > 0 )
   {
      assert(tree->path[tree->pathlen-1]->nodetype == SCIP_NODETYPE_ACTNODE);
      CHECK_OKAY( nodeAssignParent(*node, memhdr, set, tree, tree->path[tree->pathlen-1]) );
   }

   return SCIP_OKAY;
}

RETCODE SCIPnodeFree(                   /**< frees node */
   NODE**           node,               /**< node data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(memhdr != NULL);
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
      CHECK_OKAY( junctionFree(&((*node)->data.junction), memhdr) );
      break;
   case SCIP_NODETYPE_FORK:
      CHECK_OKAY( forkFree(&((*node)->data.fork), memhdr, set, lp, TRUE) );
      break;
   case SCIP_NODETYPE_SUBROOT:
      CHECK_OKAY( subrootFree(&((*node)->data.subroot), memhdr, set, lp, TRUE) );
      break;
   default:
      errorMessage("Unknown node type");
      break;
   }

   /* free common data */
   SCIPconslistFree(&((*node)->conslist), memhdr);
   SCIPdomchgFree(&((*node)->domchg), memhdr);
   nodeReleaseParent(*node, memhdr, set, lp);

   freeBlockMemory(memhdr, *node);

   return SCIP_OKAY;
}

static
RETCODE nodeDeactivate(                 /**< informs node, that it is no longer on the active path */
   NODE**           node,               /**< node to deactivate */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   Bool hasChildren = TRUE;

   assert(node != NULL);
   assert(*node != NULL);
   assert((*node)->active);
   assert(memhdr != NULL);
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
      CHECK_OKAY( SCIPnodeFree(node, memhdr, set, lp) );
   }

   return SCIP_OKAY;
}

RETCODE SCIPnodeAddConstraint(          /**< adds local constraint to the node */
   NODE*            node,               /**< node to add constraint to */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(node != NULL);

   CHECK_OKAY( SCIPconslistAdd(&(node->conslist), memhdr, cons) );

   return SCIP_OKAY;
}

RETCODE SCIPtreeAddLocalConstraint(     /**< adds local constraint to the active node */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(tree != NULL);
   assert(tree->actNode != NULL);

   CHECK_OKAY( SCIPnodeAddConstraint(tree->actNode, memhdr, cons) );

   return SCIP_OKAY;
}

RETCODE SCIPtreeAddGlobalConstraint(    /**< adds global constraint to the problem */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(tree != NULL);
   assert(tree->root != NULL);

   CHECK_OKAY( SCIPnodeAddConstraint(tree->root, memhdr, cons) );

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
   MEMHDR*          memhdr,             /**< block memory buffers */
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
      CHECK_OKAY( nodeDeactivate(&(tree->path[i]), memhdr, set, lp, tree) );
   }
   tree->pathlen = lastdepth+1;

   assert(tree->pathlen <= tree->pathsize);

   return SCIP_OKAY;
}

static
RETCODE treeSwitchPath(                 /**< switches the active path to end at the given node, applies domain changes */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory buffers */
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
   assert(memhdr != NULL);
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
      CHECK_OKAY( SCIPlpUndoDomchg(lp, set, tree->path[i]->domchg) );

   /* shrink active path to the common fork and deactivate the corresponding nodes */
   CHECK_OKAY( treeShrinkPath(tree, memhdr, set, lp, commonforkdepth) );
   assert(tree->pathlen == commonforkdepth+1);

   /* create the new active path */
   CHECK_OKAY( treeEnsurePathMem(tree, set, node->depth+1) );
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
      CHECK_OKAY( SCIPlpApplyDomchg(lp, set, tree->path[i]->domchg) );

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
   MEMHDR*          memhdr,             /**< block memory buffers */
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
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   cols = subroot->data.subroot->cols;
   rows = subroot->data.subroot->rows;
   ncols = subroot->data.subroot->ncols;
   nrows = subroot->data.subroot->nrows;

   assert(ncols == 0 || cols != NULL);
   assert(nrows == 0 || rows != NULL);
   
   for( c = 0; c < ncols; ++c )
      CHECK_OKAY( SCIPlpAddCol(lp, set, cols[c]) );
   for( r = 0; r < nrows; ++r )
      CHECK_OKAY( SCIPlpAddRow(lp, set, rows[r]) );

   return SCIP_OKAY;
}
   
static
RETCODE forkAddLP(                      /**< loads the fork's additional LP data */
   NODE*            fork,               /**< fork node to construct additional LP for */
   MEMHDR*          memhdr,             /**< block memory buffers */
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
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   cols = fork->data.fork->addedCols;
   rows = fork->data.fork->addedRows;
   ncols = fork->data.fork->naddedCols;
   nrows = fork->data.fork->naddedRows;

   assert(ncols == 0 || cols != NULL);
   assert(nrows == 0 || rows != NULL);
   
   for( c = 0; c < ncols; ++c )
      CHECK_OKAY( SCIPlpAddCol(lp, set, cols[c]) );
   for( r = 0; r < nrows; ++r )
      CHECK_OKAY( SCIPlpAddRow(lp, set, rows[r]) );

   return SCIP_OKAY;
}
   
RETCODE SCIPtreeLoadLP(                 /**< constructs the LP and loads LP state for fork/subroot of the active node */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   NODE* lpfork;
   NODE* pathnode;
   int lpforkdepth;
   int d;

   assert(memhdr != NULL);
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
         CHECK_OKAY( subrootConstructLP(tree->actSubroot, memhdr, set, lp) );
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
         CHECK_OKAY( forkAddLP(pathnode, memhdr, set, lp) );
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
      CHECK_OKAY( SCIPlpSetState(lp, memhdr, lpstate) );
   }
   
   return SCIP_OKAY;
}




/*
 * Node Conversion
 */

RETCODE SCIPnodeActivate(               /**< activates a leaf node */
   NODE*            node,               /**< leaf node to activate */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_LEAF);
   assert(!node->active);
   assert(memhdr != NULL);
   assert(tree != NULL);
   assert(lp != NULL);

   /* convert node into the actual node */
   node->nodetype = SCIP_NODETYPE_ACTNODE;

   /* track the path from the old active node to the new node, and perform domain changes */
   CHECK_OKAY( treeSwitchPath(tree, memhdr, set, lp, node) );
   assert(tree->pathlen > 0);
   assert(tree->path[tree->pathlen-1] == node);
   tree->actNode = node;

   /* remember domain changes of the active node in the dynamic data */
   CHECK_OKAY( SCIPdomchgdynCopy(tree->domchgdyn, set, node->domchg) );

   return SCIP_OKAY;
}   

RETCODE SCIPnodeToJunction(             /**< converts the active node into a junction node */
   NODE*            node,               /**< node to convert */
   MEMHDR*          memhdr,             /**< block memory buffers */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   JUNCTION* junction;

   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(node->active);
   assert(tree != NULL);

   CHECK_OKAY( junctionCreate(&junction, memhdr, tree) );

   node->nodetype = SCIP_NODETYPE_JUNCTION;
   node->data.junction = junction;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, node->depth);

   return SCIP_OKAY;
}

RETCODE SCIPnodeToFork(                 /**< converts the active node into a fork node */
   NODE*            node,               /**< node to convert */
   MEMHDR*          memhdr,             /**< block memory buffers */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   FORK* fork;

   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(node->active);
   assert(memhdr != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   CHECK_OKAY( forkCreate(&fork, memhdr, lp, tree) );
   
   node->nodetype = SCIP_NODETYPE_FORK;
   node->data.fork = fork;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, node->depth);

   return SCIP_OKAY;
}

RETCODE SCIPnodeToSubroot(              /**< converts the active node into a subroot node */
   NODE*            node,               /**< node to convert */
   MEMHDR*          memhdr,             /**< block memory buffers */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   SUBROOT* subroot;

   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(node->active);
   assert(memhdr != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   CHECK_OKAY( subrootCreate(&subroot, memhdr, lp, tree) );

   node->nodetype = SCIP_NODETYPE_SUBROOT;
   node->data.subroot = subroot;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, node->depth);

   return SCIP_OKAY;
}




/*
 * Tree methods
 */

RETCODE SCIPtreeCreate(                 /**< creates an initialized tree data structure */
   TREE**           tree,               /**< pointer to tree data structure */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(tree != NULL);
   assert(set != NULL);
   assert(set->treeGrowInit >= 0);
   assert(set->treeGrowFac >= 1.0);

   ALLOC_OKAY( allocMemory(*tree) );

   (*tree)->root = NULL;

   CHECK_OKAY( SCIPpqueueInit(&(*tree)->leaves, set->treeGrowInit, set->treeGrowFac, set->nodecmp) );

   (*tree)->path = NULL;
   (*tree)->actNode = NULL;
   (*tree)->actLPFork = NULL;
   (*tree)->actSubroot = NULL;
   (*tree)->children = NULL;
   (*tree)->siblings = NULL;

   CHECK_OKAY( SCIPdomchgdynCreate(&(*tree)->domchgdyn) );

   (*tree)->pathnlpcols = NULL;
   (*tree)->pathnlprows = NULL;
   (*tree)->pathlen = 0;
   (*tree)->pathsize = 0;
   (*tree)->correctLPDepth = -1;
   (*tree)->childrensize = 0;
   (*tree)->nchildren = 0;
   (*tree)->siblingssize = 0;
   (*tree)->nsiblings = 0;

   return SCIP_OKAY;
}

RETCODE SCIPtreeFree(                   /**< frees tree data structure */
   TREE**           tree                /**< pointer to tree data structure */
   )
{
   assert(tree != NULL);
   assert(*tree != NULL);

   SCIPpqueueFree(&(*tree)->leaves);
   SCIPdomchgdynFree(&(*tree)->domchgdyn);

   freeMemory(*tree);

   return SCIP_OKAY;
}
