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
#include "tree.h"


/*
 * dynamic memory arrays
 */

static
RETCODE treeEnsureChildrenMem(          /**< resizes children arrays to be able to store at least num nodes */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(tree != NULL);
   assert(set != NULL);

   if( num > tree->childrensize )
   {
      int newsize;
      int i;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(tree->children, newsize) );
      ALLOC_OKAY( reallocMemoryArray(tree->childrendomchg, newsize) );
      for( i = tree->childrensize; i < newsize; ++i )
      {
         CHECK_OKAY( SCIPdomchgdynCreate(&tree->childrendomchg[i], memhdr) );
      }
      tree->childrensize = newsize;
   }
   assert(num <= tree->childrensize);

   return SCIP_OKAY;
}

#if 0 /* not needed, because the siblings are always created by moving the children to the siblings array */
static
RETCODE treeEnsureSiblingsMem(          /**< resizes siblings array to be able to store at least num nodes */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(tree != NULL);
   assert(set != NULL);

   if( num > tree->siblingssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(tree->siblings, newsize) );
      ALLOC_OKAY( reallocMemoryArray(tree->siblingsdomchg, newsize) );
      for( i = tree->siblingssize; i < newsize; ++i )
      {
         CHECK_OKAY( SCIPdomchgdynCreate(&tree->siblingsdomchg[i], memhdr) );
      }
      tree->siblingssize = newsize;
   }
   assert(num <= tree->siblingssize);

   return SCIP_OKAY;
}
#endif

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

      newsize = SCIPsetCalcPathGrowSize(set, num);
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
void forkCaptureLPIState(               /**< increases the reference counter of the LP state in the fork */
   FORK*            fork,               /**< fork data */
   int              nuses               /**< number to add to the usage counter */
   )
{
   assert(fork != NULL);
   assert(fork->nlpistateref >= 0);
   assert(fork->lpistate != NULL);
   assert(nuses > 0);

   fork->nlpistateref += nuses;
   debugMessage("captured fork's LPI state %d times -> new nlpistateref=%d\n", nuses, fork->nlpistateref);
}

static
RETCODE forkReleaseLPIState(            /**< decreases the reference counter of the LP state in the fork */
   FORK*            fork,               /**< fork data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(fork != NULL);
   assert(fork->nlpistateref > 0);
   assert(fork->lpistate != NULL);
   assert(memhdr != NULL);
   assert(lp != NULL);

   fork->nlpistateref--;
   if( fork->nlpistateref == 0 )
      CHECK_OKAY( SCIPlpiFreeState(lp->lpi, memhdr, &(fork->lpistate)) );

   debugMessage("released fork's LPI state -> new nlpistateref=%d\n", fork->nlpistateref);
   return SCIP_OKAY;
}

static
void subrootCaptureLPIState(            /**< increases the reference counter of the LP state in the subroot */
   SUBROOT*         subroot,            /**< subroot data */
   int              nuses               /**< number to add to the usage counter */
   )
{
   assert(subroot != NULL);
   assert(subroot->nlpistateref >= 0);
   assert(subroot->lpistate != NULL);
   assert(nuses > 0);

   subroot->nlpistateref++;
   debugMessage("captured subroot's LPI state %d times -> new nlpistateref=%d\n", nuses, subroot->nlpistateref);
}

static
RETCODE subrootReleaseLPIState(         /**< decreases the reference counter of the LP state in the subroot */
   SUBROOT*         subroot,            /**< subroot data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(subroot != NULL);
   assert(subroot->nlpistateref > 0);
   assert(subroot->lpistate != NULL);
   assert(memhdr != NULL);
   assert(lp != NULL);

   subroot->nlpistateref--;
   if( subroot->nlpistateref == 0 )
      CHECK_OKAY( SCIPlpiFreeState(lp->lpi, memhdr, &(subroot->lpistate)) );
   
   debugMessage("released subroot's LPI state -> new nlpistateref=%d\n", subroot->nlpistateref);
   return SCIP_OKAY;
}

void SCIPnodeCaptureLPIState(           /**< increases the reference counter of the LP state in the fork or subroot node */
   NODE*            node,               /**< fork/subroot node */
   int              nuses               /**< number to add to the usage counter */
   )
{
   assert(node != NULL);

   debugMessage("capture %d times node's LPI state at depth %d\n", nuses, node->depth);
   switch( node->nodetype )
   {
   case SCIP_NODETYPE_FORK:
      forkCaptureLPIState(node->data.fork, nuses);
      break;
   case SCIP_NODETYPE_SUBROOT:
      subrootCaptureLPIState(node->data.subroot, nuses);
      break;
   default:
      errorMessage("node for capturing the LPI state is neither fork nor subroot");
      abort();
   }
}

RETCODE SCIPnodeReleaseLPIState(        /**< decreases the reference counter of the LP state in the fork or subroot node */
   NODE*            node,               /**< fork/subroot node */
   MEMHDR*          memhdr,             /**< block memory buffers */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(node != NULL);

   debugMessage("release node's LPI state at depth %d\n", node->depth);
   switch( node->nodetype )
   {
   case SCIP_NODETYPE_FORK:
      return forkReleaseLPIState(node->data.fork, memhdr, lp);
   case SCIP_NODETYPE_SUBROOT:
      return subrootReleaseLPIState(node->data.subroot, memhdr, lp);
   default:
      errorMessage("node for releasing the LPI state is neither fork nor subroot");
      return SCIP_INVALIDDATA;
   }
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

   /* increase the LPI state usage counter of the actual LP fork */
   if( tree->actlpfork != NULL )
      SCIPnodeCaptureLPIState(tree->actlpfork, tree->nchildren);

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
   const SET*       set,                /**< global SCIP settings */
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

   CHECK_OKAY( SCIPlpiGetState(lp->lpi, memhdr, set, &((*fork)->lpistate)) );
   (*fork)->nlpistateref = 0;
   (*fork)->addedcols = NULL;
   (*fork)->addedrows = NULL;
   (*fork)->naddedcols = SCIPlpGetNumNewcols(lp);
   (*fork)->naddedrows = SCIPlpGetNumNewrows(lp);
   (*fork)->nchildren = tree->nchildren;

   debugMessage("creating fork information with %d children (%d new cols, %d new rows)\n",
      (*fork)->nchildren, (*fork)->naddedcols, (*fork)->naddedrows);

   if( (*fork)->naddedcols > 0 )
   {
      /* copy the newly created columns to the fork's col array */
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*fork)->addedcols, SCIPlpGetNewcols(lp), (*fork)->naddedcols) );
   }
   if( (*fork)->naddedrows > 0 )
   {
      int i;
      
      /* copy the newly created rows to the fork's row array */
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*fork)->addedrows, SCIPlpGetNewrows(lp), (*fork)->naddedrows) );

      /* capture the added rows */
      for( i = 0; i < (*fork)->naddedrows; ++i )
         SCIProwCapture((*fork)->addedrows[i]);
   }

   /* capture the LPI state for the children */
   forkCaptureLPIState(*fork, tree->nchildren);
   
   return SCIP_OKAY;
}

static
RETCODE forkFree(                       /**< frees fork data */
   FORK**           fork,               /**< fork data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   int i;

   assert(fork != NULL);
   assert(*fork != NULL);
   assert((*fork)->nchildren == 0);
   assert((*fork)->nlpistateref == 0);
   assert((*fork)->lpistate == NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   /* release the added rows */
   for( i = 0; i < (*fork)->naddedrows; ++i )
   {
      CHECK_OKAY( SCIProwRelease(&(*fork)->addedrows[i], memhdr, set, lp) );
   }

   freeBlockMemoryArrayNull(memhdr, (*fork)->addedcols, (*fork)->naddedcols);
   freeBlockMemoryArrayNull(memhdr, (*fork)->addedrows, (*fork)->naddedrows);
   freeBlockMemory(memhdr, *fork);

   return SCIP_OKAY;
}

static
RETCODE subrootCreate(                  /**< creates subroot data */
   SUBROOT**        subroot,            /**< pointer to subroot data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   int i;
      
   assert(subroot != NULL);
   assert(memhdr != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(tree != NULL);
   assert(tree->nchildren > 0);

   ALLOC_OKAY( allocBlockMemory(memhdr, *subroot) );

   CHECK_OKAY( SCIPlpiGetState(lp->lpi, memhdr, set, &((*subroot)->lpistate)) );
   (*subroot)->nlpistateref = 0;
   (*subroot)->ncols = lp->ncols;
   (*subroot)->nrows = lp->nrows;
   (*subroot)->nchildren = tree->nchildren;
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*subroot)->cols, lp->cols, (*subroot)->ncols) );
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*subroot)->rows, lp->rows, (*subroot)->nrows) );

   /* capture the rows of the subroot */
   for( i = 0; i < (*subroot)->nrows; ++i )
      SCIProwCapture((*subroot)->rows[i]);

   /* capture the LPI state for the children */
   subrootCaptureLPIState(*subroot, tree->nchildren);
   
   return SCIP_OKAY;
}

static
RETCODE subrootFree(                    /**< frees subroot */
   SUBROOT**        subroot,            /**< subroot data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   int i;
      
   assert(subroot != NULL);
   assert(*subroot != NULL);
   assert((*subroot)->nchildren == 0);
   assert((*subroot)->nlpistateref == 0);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   CHECK_OKAY( SCIPlpiFreeState(lp->lpi, memhdr, &((*subroot)->lpistate)) );

   /* release the rows of the subroot */
   for( i = 0; i < (*subroot)->nrows; ++i )
   {
      CHECK_OKAY( SCIProwRelease(&(*subroot)->rows[i], memhdr, set, lp) );
   }

   freeBlockMemoryArrayNull(memhdr, (*subroot)->cols, (*subroot)->ncols);
   freeBlockMemoryArrayNull(memhdr, (*subroot)->rows, (*subroot)->nrows);
   freeBlockMemory(memhdr, *subroot);

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
   assert(node->nodetype == SCIP_NODETYPE_CHILD);
   assert(node->domchg == NULL);
   assert(node->data.child.arraypos == -1);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->actnode == parent);
   assert(parent == NULL || parent->nodetype == SCIP_NODETYPE_ACTNODE);

   /* link node to parent */
   node->parent = parent;
   if( parent != NULL )
   {
      node->lowerbound = parent->lowerbound;
      node->depth = parent->depth+1;
   }
   debugMessage("assigning parent %p to node %p in depth %d\n", parent, node, node->depth);

   /* register node in the childlist of the active (the parent) node */
   CHECK_OKAY( treeEnsureChildrenMem(tree, memhdr, set, tree->nchildren+1) );
   tree->children[tree->nchildren] = node;
   SCIPdomchgdynAttach(tree->childrendomchg[tree->nchildren], &node->domchg);
   node->data.child.arraypos = tree->nchildren;

   tree->nchildren++;

   return SCIP_OKAY;
}

static
void nodeReleaseParent(                 /**< decreases number of children of the parent, frees it if no children left */
   NODE*            node,               /**< child node */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   NODE* parent;
   Bool hasChildren = TRUE;

   assert(node != NULL);
   assert(memhdr != NULL);
   
   debugMessage("releasing parent-child relationship of node %p at depth %d of type %d with parent %p of type %d\n",
      node, node->depth, node->nodetype, node->parent, node->parent != NULL ? node->parent->nodetype : -1);
   parent = node->parent;
   if( parent != NULL )
   {
      switch( parent->nodetype )
      {
      case SCIP_NODETYPE_ACTNODE:
         errorMessage("Cannot release the parent-child relationship, if parent is the active node");
         abort();
      case SCIP_NODETYPE_SIBLING:
         errorMessage("Sibling cannot be a parent node");
         abort();
      case SCIP_NODETYPE_CHILD:
         errorMessage("Child cannot be a parent node");
         abort();
      case SCIP_NODETYPE_LEAF:
         errorMessage("Leaf cannot be a parent node");
         abort();
      case SCIP_NODETYPE_DEADEND:
         errorMessage("Deadend cannot be a parent node");
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
         SCIPnodeFree(&node->parent, memhdr, set, tree, lp);
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
   (*node)->lowerbound = -set->infinity;
   (*node)->depth = 0;
   (*node)->nodetype = SCIP_NODETYPE_CHILD;
   (*node)->data.child.arraypos = -1;
   (*node)->active = FALSE;
   
   if( tree->pathlen > 0 )
   {
      assert(tree->path[tree->pathlen-1]->nodetype == SCIP_NODETYPE_ACTNODE);
      CHECK_OKAY( nodeAssignParent(*node, memhdr, set, tree, tree->path[tree->pathlen-1]) );
   }
   else
   {
      /* we created the root node */
      assert(tree->actnode == NULL);
      CHECK_OKAY( nodeAssignParent(*node, memhdr, set, tree, NULL) );
   }

   debugMessage("created child node %p at depth %d\n", *node, (*node)->depth);
   return SCIP_OKAY;
}

static
void treeRemoveSibling(                 /**< removes given node from the siblings array */
   TREE*            tree,               /**< branch-and-bound tree */
   NODE*            sibling             /**< sibling node to remove */
   )
{
   DOMCHGDYN* domchgdyn;
   int delpos;

   assert(tree != NULL);
   assert(sibling != NULL);
   assert(sibling->nodetype == SCIP_NODETYPE_SIBLING);
   assert(sibling->data.sibling.arraypos >= 0 && sibling->data.sibling.arraypos < tree->nsiblings);
   assert(tree->siblings[sibling->data.sibling.arraypos] == sibling);
   assert(tree->siblings[tree->nsiblings-1]->nodetype == SCIP_NODETYPE_SIBLING);

   delpos = sibling->data.sibling.arraypos;

   /* switch domain change data of removed sibling and last sibling in array */
   domchgdyn = tree->siblingsdomchg[delpos];
   tree->siblingsdomchg[delpos] = tree->siblingsdomchg[tree->nsiblings-1];
   tree->siblingsdomchg[tree->nsiblings-1] = domchgdyn;

   /* move last sibling in array to position of removed sibling */
   tree->siblings[delpos] = tree->siblings[tree->nsiblings-1];
   tree->siblings[delpos]->data.sibling.arraypos = delpos;
   sibling->data.sibling.arraypos = -1;
   tree->nsiblings--;
}

static
void treeRemoveChild(                   /**< removes given node from the children array */
   TREE*            tree,               /**< branch-and-bound tree */
   NODE*            child               /**< child node to remove */
   )
{
   DOMCHGDYN* domchgdyn;
   int delpos;

   assert(tree != NULL);
   assert(child != NULL);
   assert(child->nodetype == SCIP_NODETYPE_CHILD);
   assert(child->data.child.arraypos >= 0 && child->data.child.arraypos < tree->nchildren);
   assert(tree->children[child->data.child.arraypos] == child);
   assert(tree->children[tree->nchildren-1]->nodetype == SCIP_NODETYPE_CHILD);

   delpos = child->data.child.arraypos;

   /* switch domain change data of removed child and last child in array */
   domchgdyn = tree->childrendomchg[delpos];
   tree->childrendomchg[delpos] = tree->childrendomchg[tree->nchildren-1];
   tree->childrendomchg[tree->nchildren-1] = domchgdyn;

   /* move last child in array to position of removed child */
   tree->children[delpos] = tree->children[tree->nchildren-1];
   tree->children[delpos]->data.child.arraypos = delpos;
   child->data.child.arraypos = -1;
   tree->nchildren--;
}

RETCODE SCIPnodeFree(                   /**< frees node */
   NODE**           node,               /**< node data */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(node != NULL);
   assert(*node != NULL);
   assert(!(*node)->active);
   assert(memhdr != NULL);
   assert(tree != NULL);

   debugMessage("free node %p at depth %d of type %d\n", *node, (*node)->depth, (*node)->nodetype);
   /* free nodetype specific data, and release no longer needed LPI states */
   switch((*node)->nodetype)
   {
   case SCIP_NODETYPE_ACTNODE:
      SCIPdomchgdynDiscard(tree->actnodedomchg, memhdr);
      /* the LPI state of the active node was already used and released when the node was activated */
      break;
   case SCIP_NODETYPE_SIBLING:
      assert((*node)->data.sibling.arraypos >= 0);
      assert((*node)->data.sibling.arraypos < tree->nsiblings);
      SCIPdomchgdynDiscard(tree->siblingsdomchg[(*node)->data.sibling.arraypos], memhdr);
      if( tree->actlpfork != NULL )
      {
         CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
      }
      treeRemoveSibling(tree, *node);
      break;
   case SCIP_NODETYPE_CHILD:
      assert((*node)->data.child.arraypos >= 0);
      assert((*node)->data.child.arraypos < tree->nchildren);
      SCIPdomchgdynDiscard(tree->childrendomchg[(*node)->data.child.arraypos], memhdr);
      /* The children capture the LPI state at the moment, where the active node is
       * converted into a junction, fork, or subroot, and a new node is activated.
       * At the same time, they become siblings or leaves, such that freeing a child
       * of the active node doesn't require to release the LPI state
       */
      treeRemoveChild(tree, *node);
      break;
   case SCIP_NODETYPE_LEAF:
      if( (*node)->data.leaf.lpfork != NULL )
      {
         CHECK_OKAY( SCIPnodeReleaseLPIState((*node)->data.leaf.lpfork, memhdr, lp) );
      }
      break;
   case SCIP_NODETYPE_DEADEND:
      break;
   case SCIP_NODETYPE_JUNCTION:
      CHECK_OKAY( junctionFree(&((*node)->data.junction), memhdr) );
      break;
   case SCIP_NODETYPE_FORK:
      CHECK_OKAY( forkFree(&((*node)->data.fork), memhdr, set, lp) );
      break;
   case SCIP_NODETYPE_SUBROOT:
      CHECK_OKAY( subrootFree(&((*node)->data.subroot), memhdr, set, lp) );
      break;
   default:
      errorMessage("Unknown node type");
      break;
   }

   /* free common data */
   SCIPconslistFree(&((*node)->conslist), memhdr, set);
   SCIPdomchgFree(&((*node)->domchg), memhdr);
   nodeReleaseParent(*node, memhdr, set, tree, lp);

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

   debugMessage("deactivate node %p at depth %d of type %d\n", *node, (*node)->depth, (*node)->nodetype);

   (*node)->active = FALSE;

   switch( (*node)->nodetype )
   {
   case SCIP_NODETYPE_ACTNODE:
      if( tree->nchildren > 0 )
      {
         errorMessage("Cannot deactivate active node with children");
         abort();
      }
      hasChildren = FALSE;
      break;
   case SCIP_NODETYPE_SIBLING:
      errorMessage("Cannot deactivate sibling (which shouldn't be active)");
      abort();
   case SCIP_NODETYPE_CHILD:
      errorMessage("Cannot deactivate child (which shouldn't be active)");
      abort();
   case SCIP_NODETYPE_LEAF:
      errorMessage("Cannot deactivate leaf (which shouldn't be active)");
      abort();
   case SCIP_NODETYPE_DEADEND:
      hasChildren = FALSE;
      break;
   case SCIP_NODETYPE_JUNCTION:
      assert((*node)->data.junction != NULL);
      hasChildren = ((*node)->data.junction->nchildren > 0);
      break;
   case SCIP_NODETYPE_FORK:
      assert((*node)->data.fork != NULL);
      hasChildren = ((*node)->data.fork->nchildren > 0);
      break;
   case SCIP_NODETYPE_SUBROOT:
      assert((*node)->data.subroot != NULL);
      hasChildren = ((*node)->data.subroot->nchildren > 0);
      break;
   default:
      errorMessage("Unknown node type");
      abort();
   }

   /* free node, if it has no children */
   if( !hasChildren )
   {
      CHECK_OKAY( SCIPnodeFree(node, memhdr, set, tree, lp) );
   }

   return SCIP_OKAY;
}

RETCODE SCIPnodeAddCons(                /**< adds local constraint to the node and captures it */
   NODE*            node,               /**< node to add constraint to */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(node != NULL);

   /* add the constraint to the node's constraint list and capture it */
   CHECK_OKAY( SCIPconslistAdd(&(node->conslist), memhdr, cons) );

   /* if the node is on the active path, add the constraint to the active constraints
    * of the constraint handler
    */
   if( node->active )
   {
      CHECK_OKAY( SCIPconsActivate(cons, set) );
   }

   return SCIP_OKAY;
}

RETCODE SCIPnodeAddBoundchg(            /**< adds bound change to active node, child or sibling of active node */
   NODE*            node,               /**< node to add bound change to */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   Real oldbound;

   assert(node != NULL);
   assert(var != NULL);
   
   debugMessage("adding boundchange at node in depth %d to variable <%s>: old bounds=[%g,%g], new %s bound: %g\n",
      node->depth, var->name, var->dom.lb, var->dom.ub, boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", newbound);
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
      oldbound = var->dom.lb;
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      oldbound = var->dom.ub;
   }

#ifndef NDEBUG
   if( SCIPsetIsEQ(set, newbound, oldbound) )
   {
      char s[255];
      sprintf(s, "variable's bound didn't change: var <%s>, oldbound=%g, newbound=%g", var->name, oldbound, newbound);
      errorMessage(s);
      return SCIP_INVALIDDATA;
   }
#endif

   switch( node->nodetype )
   {
   case SCIP_NODETYPE_ACTNODE:
      assert(tree->actnode == node);
      CHECK_OKAY( SCIPdomchgdynAddBoundchg(tree->actnodedomchg, memhdr, set, var, newbound, oldbound, boundtype) );
      CHECK_OKAY( SCIPvarChgBd(var, memhdr, set, lp, tree, newbound, boundtype) );
      return SCIP_OKAY;

   case SCIP_NODETYPE_SIBLING:
      assert(node->data.sibling.arraypos >= 0 && node->data.sibling.arraypos < tree->nsiblings);
      assert(tree->siblings[node->data.sibling.arraypos] == node);
      CHECK_OKAY( SCIPdomchgdynAddBoundchg(tree->siblingsdomchg[node->data.sibling.arraypos], memhdr, set,
                     var, newbound, oldbound, boundtype) );
      return SCIP_OKAY;

   case SCIP_NODETYPE_CHILD:
      assert(node->data.child.arraypos >= 0 && node->data.child.arraypos < tree->nchildren);
      assert(tree->children[node->data.child.arraypos] == node);
      CHECK_OKAY( SCIPdomchgdynAddBoundchg(tree->childrendomchg[node->data.child.arraypos], memhdr, set,
                     var, newbound, oldbound, boundtype) );
      return SCIP_OKAY;

   case SCIP_NODETYPE_LEAF:
   case SCIP_NODETYPE_DEADEND:
   case SCIP_NODETYPE_JUNCTION:
   case SCIP_NODETYPE_FORK:
   case SCIP_NODETYPE_SUBROOT:
      errorMessage("cannot add bound changes in nodes stored in the tree");
      return SCIP_INVALIDDATA;
   default:
      errorMessage("unknown node type");
      return SCIP_ERROR;
   }
}


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

NODETYPE SCIPnodeGetType(               /**< gets the type of the node */
   NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->nodetype;
}

int SCIPnodeGetDepth(                   /**< gets the depth of the node */
   NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->depth;
}

Real SCIPnodeGetLowerBound(             /**< gets the lower bound of the node */
   NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->lowerbound;
}

#endif



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
   assert(startdepth >= 0);
   assert(startdepth <= tree->pathlen);
   assert(tree->pathlen == 0 || startdepth < tree->pathlen);

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
      case SCIP_NODETYPE_ACTNODE:
         assert(i == tree->pathlen-1);
         break;
      case SCIP_NODETYPE_SIBLING:
         errorMessage("Sibling cannot be in the active path");
         abort();
      case SCIP_NODETYPE_CHILD:
         errorMessage("Child cannot be in the active path");
         abort();
      case SCIP_NODETYPE_LEAF:
         errorMessage("Leaf cannot be in the active path");
         abort();
      case SCIP_NODETYPE_DEADEND:
         errorMessage("Deadend cannot be in the active path");
         abort();
      case SCIP_NODETYPE_JUNCTION:
         assert(node->data.junction != NULL);
         break;
      case SCIP_NODETYPE_FORK:
         assert(node->data.fork != NULL);
         ncols += node->data.fork->naddedcols;
         nrows += node->data.fork->naddedrows;
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
RETCODE treeShrinkPath(                 /**< cuts off path of active nodes after given node, marks cutted nodes inactive */
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
   NODE*            node                /**< last node of the new active path (= new active node), or NULL */
   )
{
   NODE* commonfork;    /* common fork node */
   NODE* lpfork;        /* fork node defining the LP state of the new active node */
   NODE* subroot;       /* subroot of new active path */
   int nodedepth;       /* depth of the node, or -1 of node == NULL */
   int commonforkdepth; /* depth of the common fork node, or -1 if no common fork exists */
   int i;

   assert(tree != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(node == NULL || !node->active);
   assert(node == NULL || node->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(tree->actlpfork == NULL || tree->actlpfork->active);
   assert(tree->actsubroot == NULL || tree->actsubroot->active);

   /* get the node's depth */
   if( node == NULL )
      nodedepth = -1;
   else
      nodedepth = node->depth;
   debugMessage("switch path: nodedepth=%d\n", nodedepth);

   /* find the common fork node, the new LP defining fork, and the new active subroot */
   commonfork = node;
   lpfork = NULL;
   subroot = NULL;
   while( commonfork != NULL && !commonfork->active )
   {
      commonfork = commonfork->parent;
      if( commonfork != NULL )
      {
         if( lpfork == NULL
            && (commonfork->nodetype == SCIP_NODETYPE_FORK || commonfork->nodetype == SCIP_NODETYPE_SUBROOT) )
            lpfork = commonfork;
         if( subroot == NULL && commonfork->nodetype == SCIP_NODETYPE_SUBROOT )
            subroot = commonfork;
      }
   }
   if( commonfork == NULL )
      commonforkdepth = -1;
   else
      commonforkdepth = commonfork->depth;
   assert(lpfork == NULL || !lpfork->active || lpfork == commonfork);
   assert(subroot == NULL || !subroot->active || subroot == commonfork);
   debugMessage("switch path: commonforkdepth=%d\n", commonforkdepth);

   /* if not already found, continue searching the LP defining fork */
   if( lpfork == NULL )
   {
      if( tree->actlpfork != NULL && tree->actlpfork->depth > commonforkdepth )
      {
         /* actlpfork is not on the same active path as the new node: we have to search again */
         lpfork = commonfork;
         while( lpfork != NULL &&
            (lpfork->nodetype == SCIP_NODETYPE_FORK || lpfork->nodetype == SCIP_NODETYPE_SUBROOT) )
            lpfork = lpfork->parent;
      }
      else
      {
         /* actlpfork is on the same active path as the new node: old and new node have the same lpfork */
         lpfork = tree->actlpfork;
      }
   }
   debugMessage("switch path: lpforkdepth=%d\n", lpfork == NULL ? -1 : lpfork->depth);

   /* if not already found, continue searching the subroot */
   if( subroot == NULL )
   {
      if( tree->actsubroot != NULL && tree->actsubroot->depth >= commonforkdepth )
      {
         subroot = commonfork;
         while( subroot != NULL && subroot->nodetype != SCIP_NODETYPE_SUBROOT )
            subroot = subroot->parent;
      }
      else
         subroot = tree->actsubroot;
   }
   debugMessage("switch path: subrootdepth=%d\n", subroot == NULL ? -1 : subroot->depth);
   
   debugMessage("switch path: old correctlpdepth=%d\n", tree->correctlpdepth);
   /* remember the depth of the common fork node for LP updates */
   if( subroot == tree->actsubroot )
   {
      /* we are in the same subtree: the LP is correct at most upto the fork depth */
      assert(subroot == NULL || subroot->active);
      tree->correctlpdepth = MIN(tree->correctlpdepth, commonforkdepth);
   }
   else
   {
      /* we are in a different subtree: the LP is completely incorrect */
      assert(subroot == NULL || !subroot->active);
      tree->correctlpdepth = -1;
   }
   debugMessage("switch path: new correctlpdepth=%d\n", tree->correctlpdepth);

   debugMessage("switch path: pathlen=%d\n", tree->pathlen);   
   /* undo the domain changes of the old active path */
   for( i = tree->pathlen-1; i > commonforkdepth; --i )
   {
      debugMessage("switch path: undo domain changes in depth %d\n", i);
      CHECK_OKAY( SCIPdomchgUndo(tree->path[i]->domchg, memhdr, set, lp, tree) );
   }

   /* shrink active path to the common fork and deactivate the corresponding nodes */
   CHECK_OKAY( treeShrinkPath(tree, memhdr, set, lp, commonforkdepth) );
   assert(tree->pathlen == commonforkdepth+1);

   /* create the new active path */
   CHECK_OKAY( treeEnsurePathMem(tree, set, nodedepth+1) );
   tree->pathlen = nodedepth+1;
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
   {
      debugMessage("switch path: apply domain changes in depth %d\n", i);
      CHECK_OKAY( SCIPdomchgApply(tree->path[i]->domchg, memhdr, set, lp, tree) );
   }

   /* remember LP defining fork and subroot */
   assert(subroot == NULL || lpfork != NULL);
   assert(subroot == NULL || subroot->depth <= lpfork->depth);
   tree->actlpfork = lpfork;
   tree->actsubroot = subroot;
   
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

   cols = fork->data.fork->addedcols;
   rows = fork->data.fork->addedrows;
   ncols = fork->data.fork->naddedcols;
   nrows = fork->data.fork->naddedrows;

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

   assert(tree != NULL);
   assert(tree->path != NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnode->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(tree->pathlen > 0);
   assert(tree->path[tree->pathlen-1] == tree->actnode);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   lpfork = tree->actlpfork;

   /* find out the lpfork's depth (or -1, if lpfork is NULL) */
   if( lpfork == NULL )
   {
      assert(tree->correctlpdepth == -1);
      assert(tree->actsubroot == NULL);
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
   if( tree->correctlpdepth >= 0 )
   {
      /* same subtree: shrink LP to the deepest node with correct LP */
      assert(lpfork != NULL);
      assert(tree->correctlpdepth <= lpforkdepth);
      CHECK_OKAY( SCIPlpShrinkCols(lp, tree->pathnlpcols[tree->correctlpdepth]) );
      CHECK_OKAY( SCIPlpShrinkRows(lp, memhdr, set, tree->pathnlprows[tree->correctlpdepth]) );
   }
   else
   {
      /* other subtree: fill LP with the subroot LP data */
      CHECK_OKAY( SCIPlpClear(lp, memhdr, set) );
      if( tree->actsubroot != NULL )
      {
         CHECK_OKAY( subrootConstructLP(tree->actsubroot, memhdr, set, lp) );
         tree->correctlpdepth = tree->actsubroot->depth; 
      }
   }

   assert(lpforkdepth < tree->pathlen);
   assert(tree->correctlpdepth <= lpforkdepth);

   /* add the missing columns and rows */
   for( d = tree->correctlpdepth+1; d <= lpforkdepth; ++d )
   {
      pathnode = tree->path[d];
      assert(pathnode != NULL);
      assert(pathnode->depth == d);
      assert(pathnode->nodetype == SCIP_NODETYPE_JUNCTION || pathnode->nodetype == SCIP_NODETYPE_FORK);
      if( pathnode->nodetype == SCIP_NODETYPE_FORK )
      {
         CHECK_OKAY( forkAddLP(pathnode, memhdr, set, lp) );
      }
   }
   tree->correctlpdepth = lpforkdepth;

   /* load LP state, if existing */
   if( lpfork != NULL )
   {
      if( lpfork->nodetype == SCIP_NODETYPE_FORK )
      {
         assert(lpfork->data.fork != NULL);
         CHECK_OKAY( SCIPlpSetState(lp, memhdr, set, lpfork->data.fork->lpistate) );
         CHECK_OKAY( forkReleaseLPIState(lpfork->data.fork, memhdr, lp) );
      }
      else
      {
         assert(lpfork->nodetype == SCIP_NODETYPE_SUBROOT);
         assert(lpfork->data.subroot != NULL);
         CHECK_OKAY( SCIPlpSetState(lp, memhdr, set, lpfork->data.subroot->lpistate) );
         CHECK_OKAY( subrootReleaseLPIState(lpfork->data.subroot, memhdr, lp) );
      }
   }

   /* mark the LP's size, such that we know which rows and columns were added in the new node */
   SCIPlpMarkSize(lp);
   
   return SCIP_OKAY;
}




/*
 * Node Conversion
 */

static
RETCODE nodeToLeaf(                     /**< converts node into LEAF and puts it in the array on the node queue */
   NODE*            node,               /**< child or sibling node to convert */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   NODE* lpfork;

   assert(node->nodetype == SCIP_NODETYPE_SIBLING || node->nodetype == SCIP_NODETYPE_CHILD);
   
   /* find the lpfork for the node */
   lpfork = tree->actlpfork;
   if( node->parent != NULL )
   {
      if( node->parent->nodetype == SCIP_NODETYPE_FORK || node->parent->nodetype == SCIP_NODETYPE_SUBROOT )
         lpfork = node->parent;
   }
   
   /* convert node into leaf */
   debugMessage("convert node %p at depth %d to leaf with lpfork %p\n", node, node->depth, lpfork);
   node->nodetype = SCIP_NODETYPE_LEAF;
   node->data.leaf.lpfork = lpfork;
   
   /* insert leaf in node queue */
   CHECK_OKAY( SCIPnodepqInsert(tree->leaves, set, node) );

   return SCIP_OKAY;
}

static
RETCODE actnodeToDeadend(               /**< converts the active node into a deadend node */
   MEMHDR*          memhdr,             /**< block memory buffers */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(memhdr != NULL);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnode->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(tree->actnode->active);
   assert(tree->nchildren == 0);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   /* detach dynamic size data from domain change data of old active node */
   CHECK_OKAY( SCIPdomchgdynDetach(tree->actnodedomchg, memhdr) );

   tree->actnode->nodetype = SCIP_NODETYPE_DEADEND;

   return SCIP_OKAY;
}

static
RETCODE actnodeToJunction(              /**< converts the active node into a junction node */
   MEMHDR*          memhdr,             /**< block memory buffers */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   JUNCTION* junction;

   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnode->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(tree->actnode->active);

   /* detach dynamic size data from domain change data of old active node */
   CHECK_OKAY( SCIPdomchgdynDetach(tree->actnodedomchg, memhdr) );

   CHECK_OKAY( junctionCreate(&junction, memhdr, tree) );

   tree->actnode->nodetype = SCIP_NODETYPE_JUNCTION;
   tree->actnode->data.junction = junction;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, tree->actnode->depth);

   return SCIP_OKAY;
}

static
RETCODE actnodeToFork(                  /**< converts the active node into a fork node */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   FORK* fork;

   assert(memhdr != NULL);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnode->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(tree->actnode->active);
   assert(tree->nchildren > 0);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   /* detach dynamic size data from domain change data of old active node */
   CHECK_OKAY( SCIPdomchgdynDetach(tree->actnodedomchg, memhdr) );

   CHECK_OKAY( forkCreate(&fork, memhdr, set, lp, tree) );
   
   tree->actnode->nodetype = SCIP_NODETYPE_FORK;
   tree->actnode->data.fork = fork;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, tree->actnode->depth);

   return SCIP_OKAY;
}

static
RETCODE actnodeToSubroot(               /**< converts the active node into a subroot node */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   SUBROOT* subroot;

   assert(memhdr != NULL);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnode->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(tree->actnode->active);
   assert(tree->nchildren > 0);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   /* detach dynamic size data from domain change data of old active node */
   CHECK_OKAY( SCIPdomchgdynDetach(tree->actnodedomchg, memhdr) );

   CHECK_OKAY( subrootCreate(&subroot, memhdr, set, lp, tree) );

   tree->actnode->nodetype = SCIP_NODETYPE_SUBROOT;
   tree->actnode->data.subroot = subroot;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, tree->actnode->depth);

   return SCIP_OKAY;
}

static
RETCODE treeNodesToQueue(               /**< puts all nodes in the array on the node queue and makes them LEAFs */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   NODE**           nodes,              /**< array of nodes to put on the queue */
   DOMCHGDYN**      domchgdyn,          /**< array of dynamic domain changes */
   int*             nnodes              /**< pointer to number of nodes in the array */
   )
{
   NODE* node;
   int i;

   assert(tree != NULL);
   assert(set != NULL);
   assert(nnodes != NULL);
   assert(*nnodes == 0 || nodes != NULL);
   assert(*nnodes == 0 || domchgdyn != NULL);

   for( i = 0; i < *nnodes; ++i )
   {
      node = nodes[i];

      /* detach the dynamic size attachment of the domain change data to shrink the node's domain change data */
      CHECK_OKAY( SCIPdomchgdynDetach(domchgdyn[i], memhdr) );

      /* convert node to LEAF and put it into leaves queue */
      nodeToLeaf(node, set, tree);
   }
   *nnodes = 0;

   return SCIP_OKAY;
}

static
void treeChildrenToSiblings(            /**< converts children into siblings, clears children array */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   NODE** tmpnodes;
   DOMCHGDYN** tmpdomchg;
   int tmpnodessize;
   int i;

   assert(tree != NULL);
   assert(tree->nsiblings == 0);

   tmpnodes = tree->siblings;
   tmpdomchg = tree->siblingsdomchg;
   tmpnodessize = tree->siblingssize;

   tree->siblings = tree->children;
   tree->siblingsdomchg = tree->childrendomchg;
   tree->nsiblings = tree->nchildren;
   tree->siblingssize = tree->childrensize;

   tree->children = tmpnodes;
   tree->childrendomchg = tmpdomchg;
   tree->nchildren = 0;
   tree->childrensize = tmpnodessize;
   
   for( i = 0; i < tree->nsiblings; ++i )
   {
      assert(tree->siblings[i]->nodetype == SCIP_NODETYPE_CHILD);
      tree->siblings[i]->nodetype = SCIP_NODETYPE_SIBLING;

      /* because CHILD.arraypos and SIBLING.arraypos are on the same position, we do not have to copy it */
      assert(&(tree->siblings[i]->data.sibling.arraypos) == &(tree->siblings[i]->data.child.arraypos));
   }
}

RETCODE SCIPnodeActivate(               /**< activates a child, a sibling, or a leaf node */
   NODE*            node,               /**< leaf node to activate (or NULL to deactivate all nodes) */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(node == NULL
      || node->nodetype == SCIP_NODETYPE_SIBLING
      || node->nodetype == SCIP_NODETYPE_CHILD
      || node->nodetype == SCIP_NODETYPE_LEAF);
   assert(node == NULL || !node->active);
   assert(memhdr != NULL);
   assert(tree != NULL);
   assert(lp != NULL);

   /* convert the old active node into a fork node, if it has children */
   if( tree->actnode != NULL )
   {
      assert(tree->actnode->nodetype == SCIP_NODETYPE_ACTNODE);

      if( tree->nchildren > 0 )
      {
         /* convert old active node into a fork node */
         todoMessage("decide: old active node becomes junction, fork, or subroot");
         CHECK_OKAY( actnodeToFork(memhdr, set, tree, lp) );
      }
      else
      {
         CHECK_OKAY( actnodeToDeadend(memhdr, tree, lp) );
      }
   }

   /* set up the new lists of siblings and children */
   if( node == NULL )
   {
      /* move siblings to the queue, make them LEAFs */
      CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, tree->siblings, tree->siblingsdomchg, &tree->nsiblings) );

      /* move children to the queue, make them LEAFs */
      CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, tree->children, tree->childrendomchg, &tree->nchildren) );
   }
   else
   {
      DOMCHGDYN* domchgdyn;

      switch( node->nodetype )
      {
      case SCIP_NODETYPE_SIBLING:
         /* move children to the queue, make them LEAFs */
         CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, tree->children, tree->childrendomchg, &tree->nchildren) );

         /* switch domain change data of sibling and active node */
         domchgdyn = tree->actnodedomchg;
         tree->actnodedomchg = tree->siblingsdomchg[node->data.sibling.arraypos];
         tree->siblingsdomchg[node->data.sibling.arraypos] = domchgdyn;

         /* remove selected sibling from the siblings array */
         treeRemoveSibling(tree, node);
         
         break;
         
      case SCIP_NODETYPE_CHILD:
         /* move siblings to the queue, make them LEAFs */
         CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, tree->siblings, tree->siblingsdomchg, &tree->nsiblings) );

         /* switch domain change data of child and active node */
         domchgdyn = tree->actnodedomchg;
         tree->actnodedomchg = tree->childrendomchg[node->data.child.arraypos];
         tree->childrendomchg[node->data.child.arraypos] = domchgdyn;

         /* remove selected child from the children array */      
         treeRemoveChild(tree, node);
         
         /* move other children to the siblings array, make them SIBLINGs */
         treeChildrenToSiblings(tree);
         
         break;
         
      case SCIP_NODETYPE_LEAF:
         /* move siblings to the queue, make them LEAFs */
         CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, tree->siblings, tree->siblingsdomchg, &tree->nsiblings) );
         
         /* move children to the queue, make them LEAFs */
         CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, tree->children, tree->childrendomchg, &tree->nchildren) );

         /* attach dynamic size data to domain changes of the active node */
         SCIPdomchgdynAttach(tree->actnodedomchg, &node->domchg);

         /* remove node from the queue */
         if( node != SCIPnodepqRemove(tree->leaves, set) )
         {
            errorMessage("Selected node is a leaf, but not the first on the queue");
            return SCIP_INVALIDDATA;
         }
         break;
      default:
         errorMessage("Selected node is neither sibling, child, nor leaf");
         return SCIP_INVALIDDATA;
      }

      /* convert node into the active node */
      node->nodetype = SCIP_NODETYPE_ACTNODE;
   }

   /* track the path from the old active node to the new node, and perform domain changes */
   CHECK_OKAY( treeSwitchPath(tree, memhdr, set, lp, node) );
   assert(node == NULL || tree->pathlen > 0);
   assert(node != NULL || tree->pathlen == 0);
   assert(node == NULL || tree->path[tree->pathlen-1] == node);
   assert(tree->nchildren == 0);
   tree->actnode = node;

#ifndef NDEBUG
   {
      int i;
      for( i = 0; i < tree->nchildren; ++i )
         assert(SCIPdomchgdynGetDomchgPtr(tree->childrendomchg[i]) == &tree->children[i]->domchg);
      for( i = 0; i < tree->nsiblings; ++i )
         assert(SCIPdomchgdynGetDomchgPtr(tree->siblingsdomchg[i]) == &tree->siblings[i]->domchg);
      if( tree->actnode != NULL )
         assert(SCIPdomchgdynGetDomchgPtr(tree->actnodedomchg) == &tree->actnode->domchg);
      else
         assert(SCIPdomchgdynGetDomchgPtr(tree->actnodedomchg) == NULL);
   }
#endif
   
   return SCIP_OKAY;
}   




/*
 * Tree methods
 */

RETCODE SCIPtreeCreate(                 /**< creates an initialized tree data structure */
   TREE**           tree,               /**< pointer to tree data structure */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   LP*              lp,                 /**< actual LP data */
   PROB*            prob                /**< problem data */
   )
{
   VAR* var;
   int v;

   assert(tree != NULL);
   assert(set != NULL);
   assert(set->treeGrowInit >= 0);
   assert(set->treeGrowFac >= 1.0);
   assert(lp != NULL);
   assert(prob != NULL);

   ALLOC_OKAY( allocMemory(*tree) );

   (*tree)->root = NULL;

   CHECK_OKAY( SCIPnodepqCreate(&(*tree)->leaves, set) );

   (*tree)->path = NULL;
   (*tree)->actnode = NULL;
   (*tree)->actlpfork = NULL;
   (*tree)->actsubroot = NULL;
   (*tree)->children = NULL;
   (*tree)->siblings = NULL;

   CHECK_OKAY( SCIPdomchgdynCreate(&(*tree)->actnodedomchg, memhdr) );
   (*tree)->childrendomchg = NULL;
   (*tree)->siblingsdomchg = NULL;

   /* create root pseudo solution */
   CHECK_OKAY( SCIPsolCreate(&(*tree)->actpseudosol, memhdr, stat, NULL) );
   for( v = 0; v < prob->nvars; ++v )
   {
      var = prob->vars[v];
      assert(SCIPsolGetVal((*tree)->actpseudosol, var) == 0.0);
      CHECK_OKAY( SCIPsolSetVal((*tree)->actpseudosol, memhdr, set, var, SCIPvarGetBestBound(var)) );
   }
   
   (*tree)->pathnlpcols = NULL;
   (*tree)->pathnlprows = NULL;
   (*tree)->pathlen = 0;
   (*tree)->pathsize = 0;
   (*tree)->correctlpdepth = -1;
   (*tree)->childrensize = 0;
   (*tree)->nchildren = 0;
   (*tree)->siblingssize = 0;
   (*tree)->nsiblings = 0;

   /* create root node */
   CHECK_OKAY( SCIPnodeCreate(&(*tree)->root, memhdr, set, *tree) );

   return SCIP_OKAY;
}

RETCODE SCIPtreeFree(                   /**< frees tree data structure */
   TREE**           tree,               /**< pointer to tree data structure */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   int i;

   assert(tree != NULL);
   assert(*tree != NULL);
   assert((*tree)->nchildren == 0);
   assert((*tree)->nsiblings == 0);
   assert((*tree)->actnode == NULL);

   debugMessage("free tree\n");
#ifdef SCIP_BLOCKMEMORY
   /* we don't need to free the nodes and domain changes, because they are stored in block memory
    * we just need to free the queue data structure
    */
   /* SCIPnodepqDestroy(&(*tree)->leaves); */
   todoMessage("don't free the tree, but in some way release all variables and constraints");
#endif
   /* #else */
   /* free node queue */
   CHECK_OKAY( SCIPnodepqFree(&(*tree)->leaves, memhdr, set, *tree, lp) );
   
   /* free dynamic size domain change attachments */
   for( i = 0; i < (*tree)->childrensize; ++i )
      SCIPdomchgdynFree(&(*tree)->childrendomchg[i], memhdr);
   for( i = 0; i < (*tree)->siblingssize; ++i )
      SCIPdomchgdynFree(&(*tree)->siblingsdomchg[i], memhdr);
   SCIPdomchgdynFree(&(*tree)->actnodedomchg, memhdr);

   /* release actual pseudo solution */
   CHECK_OKAY( SCIPsolRelease(&(*tree)->actpseudosol, memhdr, set, lp) );
   /* #endif */

   /* free pointer arrays */
   freeMemoryArrayNull((*tree)->path);
   freeMemoryArrayNull((*tree)->children);
   freeMemoryArrayNull((*tree)->siblings);   
   freeMemoryArrayNull((*tree)->childrendomchg);
   freeMemoryArrayNull((*tree)->siblingsdomchg);   
   freeMemoryArrayNull((*tree)->pathnlpcols);
   freeMemoryArrayNull((*tree)->pathnlprows);

   freeMemory(*tree);

   return SCIP_OKAY;
}

RETCODE SCIPtreeCutoff(                 /**< cuts off nodes with lower bound not better than given upper bound */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             upperbound          /**< upper bound: all nodes with lowerbound >= upperbound are cut off */
   )
{
   NODE* node;
   int i;

   assert(tree != NULL);

   /* cut off leaf nodes in the queue */
   CHECK_OKAY( SCIPnodepqBound(tree->leaves, memhdr, set, tree, lp, upperbound) );

   /* cut off siblings: we have to loop backwards, because a removal leads to moving the last node in empty slot */
   for( i = tree->nsiblings-1; i >= 0; --i )
   {
      node = tree->siblings[i];
      if( SCIPsetIsGE(set, node->lowerbound, upperbound) )
      {
         debugMessage("cut off sibling %p in depth %d with lowerbound=%g at position %d\n", 
            node, node->depth, node->lowerbound, i);
         CHECK_OKAY( SCIPnodeFree(&node, memhdr, set, tree, lp) );
      }
   }

   /* cut off children: we have to loop backwards, because a removal leads to moving the last node in empty slot */
   for( i = tree->nchildren-1; i >= 0; --i )
   {
      node = tree->children[i];
      if( SCIPsetIsGE(set, node->lowerbound, upperbound) )
      {
         debugMessage("cut off child %p in depth %d with lowerbound=%g at position %d\n",
            node, node->depth, node->lowerbound, i);
         CHECK_OKAY( SCIPnodeFree(&node, memhdr, set, tree, lp) );
      }
   }

   return SCIP_OKAY;
}

RETCODE SCIPtreeAddLocalCons(           /**< adds local constraint to the active node and captures it */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(tree != NULL);
   assert(tree->actnode != NULL);

   /* add constraint to active node and capture it */
   CHECK_OKAY( SCIPnodeAddCons(tree->actnode, memhdr, set, cons) );

   return SCIP_OKAY;
}

RETCODE SCIPtreeAddGlobalCons(          /**< adds global constraint to the problem and captures it */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(tree != NULL);
   assert(tree->root != NULL);

   /* add constraint to root node and capture it */
   CHECK_OKAY( SCIPnodeAddCons(tree->root, memhdr, set, cons) );

   return SCIP_OKAY;
}

RETCODE SCIPtreeBoundChanged(           /**< notifies tree, that a bound of a variable changed */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   assert(var != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(tree != NULL);
      if( var->obj > 0.0 && boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         CHECK_OKAY( SCIPsolSetVal(tree->actpseudosol, memhdr, set, var, var->dom.lb) );
      }
      else if( var->obj < 0.0 && boundtype == SCIP_BOUNDTYPE_UPPER )
      {
         CHECK_OKAY( SCIPsolSetVal(tree->actpseudosol, memhdr, set, var, var->dom.ub) );
      }
      break;

   case SCIP_VARSTATUS_ORIGINAL:
   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_AGGREGATED:
   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("tree was informed of a bound change of a non-mutable variable");
      return SCIP_INVALIDDATA;
   default:
      errorMessage("unknown variable status");
      abort();
   }

   return SCIP_OKAY;
}

int SCIPtreeGetNLeaves(                 /**< gets number of leaves */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqLen(tree->leaves);
}
   
int SCIPtreeGetNNodes(                  /**< gets number of nodes (children + siblings + leaves + active) */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(tree != NULL);

   return tree->nchildren + tree->nsiblings + SCIPtreeGetNLeaves(tree) + (tree->actnode != NULL ? 1 : 0);
}

NODE* SCIPtreeGetBestChild(             /**< gets the best child of the active node */
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set                 /**< global SCIP settings */
   )
{
   SCIP* scip;
   NODESEL* nodesel;
   NODE* bestnode;
   int i;

   assert(tree != NULL);

   scip = set->scip;
   nodesel = set->nodesel;
   assert(scip != NULL);
   assert(nodesel != NULL);

   bestnode = NULL;
   for( i = 0; i < tree->nchildren; ++i )
   {
      if( bestnode == NULL || SCIPnodeselCompare(nodesel, scip, tree->children[i], bestnode) < 0 )
      {
         bestnode = tree->children[i];
      }
   }

   return bestnode;
}

NODE* SCIPtreeGetBestSibling(           /**< gets the best sibling of the active node */
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set                 /**< global SCIP settings */
   )
{
   SCIP* scip;
   NODESEL* nodesel;
   NODE* bestnode;
   int i;

   assert(tree != NULL);

   scip = set->scip;
   nodesel = set->nodesel;
   assert(scip != NULL);
   assert(nodesel != NULL);

   bestnode = NULL;
   for( i = 0; i < tree->nsiblings; ++i )
   {
      if( bestnode == NULL || SCIPnodeselCompare(nodesel, scip, tree->siblings[i], bestnode) < 0 )
      {
         bestnode = tree->siblings[i];
      }
   }
   
   return bestnode;
}

NODE* SCIPtreeGetBestLeaf(              /**< gets the best leaf from the node queue */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqFirst(tree->leaves);
}

NODE* SCIPtreeGetBestNode(              /**< gets the best node from the tree (child, sibling, or leaf) */
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set                 /**< global SCIP settings */
   )
{
   SCIP* scip;
   NODESEL* nodesel;
   NODE* bestchild;
   NODE* bestsibling;
   NODE* bestleaf;
   NODE* bestnode;

   assert(tree != NULL);
   assert(set != NULL);

   scip = set->scip;
   nodesel = set->nodesel;
   assert(scip != NULL);
   assert(nodesel != NULL);

   /* get the best child, sibling, and leaf */
   bestchild = SCIPtreeGetBestChild(tree, set);
   bestsibling = SCIPtreeGetBestSibling(tree, set);
   bestleaf = SCIPtreeGetBestLeaf(tree);

   /* return the best of the three */
   bestnode = bestchild;
   if( bestsibling != NULL && (bestnode == NULL || SCIPnodeselCompare(nodesel, scip, bestsibling, bestnode) < 0) )
      bestnode = bestsibling;
   if( bestleaf != NULL && (bestnode == NULL || SCIPnodeselCompare(nodesel, scip, bestleaf, bestnode) < 0) )
      bestnode = bestleaf;

   assert(SCIPtreeGetNLeaves(tree) == 0 || bestnode != NULL);

   return bestnode;
}

Real SCIPtreeGetLowerbound(             /**< gets the minimal lower bound of all nodes in the tree */
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set                 /**< global SCIP settings */
   )
{
   Real lowerbound;
   int i;

   assert(tree != NULL);
   assert(set != NULL);
   assert(set->nodesel != NULL);

   /* get the lower bound from the queue */
   lowerbound = SCIPnodepqGetLowerbound(tree->leaves, set);

   /* compare lower bound with active node */
   if( tree->actnode != NULL )
   {
      lowerbound = MIN(lowerbound, tree->actnode->lowerbound);
   }

   /* compare lower bound with siblings */
   for( i = 0; i < tree->nsiblings; ++i )
   {
      assert(tree->siblings[i] != NULL);
      lowerbound = MIN(lowerbound, tree->siblings[i]->lowerbound); 
   }

   /* compare lower bound with children */
   for( i = 0; i < tree->nchildren; ++i )
   {
      assert(tree->children[i] != NULL);
      lowerbound = MIN(lowerbound, tree->children[i]->lowerbound); 
   }

   return lowerbound;
}

Real SCIPtreeGetActLowerbound(          /**< gets the lower bound of the active node */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(tree != NULL);

   if( tree->actnode != NULL )
      return tree->actnode->lowerbound;
   else
      return SCIP_INVALID;
}

Real SCIPtreeGetAvgLowerbound(          /**< gets the average lower bound of all nodes in the tree */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   Real lowerboundsum;
   int nnodes;
   int i;

   assert(tree != NULL);

   /* get sum of lower bounds from nodes in the queue */
   lowerboundsum = SCIPnodepqGetLowerboundSum(tree->leaves);
   
   /* add lower bound of active node */
   if( tree->actnode != NULL )
   {
      lowerboundsum += tree->actnode->lowerbound;
   }

   /* add lower bounds of siblings */
   for( i = 0; i < tree->nsiblings; ++i )
   {
      assert(tree->siblings[i] != NULL);
      lowerboundsum += tree->siblings[i]->lowerbound;
   }

   /* add lower bounds of children */
   for( i = 0; i < tree->nchildren; ++i )
   {
      assert(tree->children[i] != NULL);
      lowerboundsum += tree->children[i]->lowerbound;
   }

   nnodes = SCIPtreeGetNNodes(tree);

   return nnodes == 0 ? 0.0 : lowerboundsum/nnodes;
}
