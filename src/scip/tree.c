/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
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

/** resizes children arrays to be able to store at least num nodes */
static
RETCODE treeEnsureChildrenMem(
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
      ALLOC_OKAY( reallocMemoryArray(&tree->children, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&tree->childrenconssetchg, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&tree->childrendomchg, newsize) );
      for( i = tree->childrensize; i < newsize; ++i )
      {
         CHECK_OKAY( SCIPconssetchgdynCreate(&tree->childrenconssetchg[i], memhdr) );
         CHECK_OKAY( SCIPdomchgdynCreate(&tree->childrendomchg[i], memhdr) );
      }
      tree->childrensize = newsize;
   }
   assert(num <= tree->childrensize);

   return SCIP_OKAY;
}

/** resizes path array to be able to store at least num nodes */
static
RETCODE treeEnsurePathMem(
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
      ALLOC_OKAY( reallocMemoryArray(&tree->path, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&tree->pathnlpcols, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&tree->pathnlprows, newsize) );
      tree->pathsize = newsize;
   }
   assert(num <= tree->pathsize);

   return SCIP_OKAY;
}




/*
 * Node methods
 */

/** node comparator for best lower bound */
DECL_SORTPTRCOMP(SCIPnodeCmpLowerbound)
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

/** increases the reference counter of the LP state in the fork */
static
void forkCaptureLPIState(
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

/** decreases the reference counter of the LP state in the fork */
static
RETCODE forkReleaseLPIState(
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

/** increases the reference counter of the LP state in the subroot */
static
void subrootCaptureLPIState(
   SUBROOT*         subroot,            /**< subroot data */
   int              nuses               /**< number to add to the usage counter */
   )
{
   assert(subroot != NULL);
   assert(subroot->nlpistateref >= 0);
   assert(subroot->lpistate != NULL);
   assert(nuses > 0);

   subroot->nlpistateref += nuses;
   debugMessage("captured subroot's LPI state %d times -> new nlpistateref=%d\n", nuses, subroot->nlpistateref);
}

/** decreases the reference counter of the LP state in the subroot */
static
RETCODE subrootReleaseLPIState(
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

/** increases the reference counter of the LP state in the fork or subroot node */
void SCIPnodeCaptureLPIState(
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

/** decreases the reference counter of the LP state in the fork or subroot node */
RETCODE SCIPnodeReleaseLPIState(
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

/** initializes junction data */
static
RETCODE junctionInit(
   JUNCTION*        junction,           /**< pointer to junction data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(junction != NULL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);

   junction->nchildren = tree->nchildren;

   /* increase the LPI state usage counter of the actual LP fork */
   if( tree->actlpfork != NULL )
      SCIPnodeCaptureLPIState(tree->actlpfork, tree->nchildren);

   return SCIP_OKAY;
}

/** creates fork data */
static
RETCODE forkCreate(
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

   ALLOC_OKAY( allocBlockMemory(memhdr, fork) );

   CHECK_OKAY( SCIPlpGetState(lp, memhdr, &((*fork)->lpistate)) );
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
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*fork)->addedcols, SCIPlpGetNewcols(lp), (*fork)->naddedcols) );
   }
   if( (*fork)->naddedrows > 0 )
   {
      int i;
      
      /* copy the newly created rows to the fork's row array */
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*fork)->addedrows, SCIPlpGetNewrows(lp), (*fork)->naddedrows) );

      /* capture the added rows */
      for( i = 0; i < (*fork)->naddedrows; ++i )
         SCIProwCapture((*fork)->addedrows[i]);
   }

   /* capture the LPI state for the children */
   forkCaptureLPIState(*fork, tree->nchildren);
   
   return SCIP_OKAY;
}

/** frees fork data */
static
RETCODE forkFree(
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

   freeBlockMemoryArrayNull(memhdr, &(*fork)->addedcols, (*fork)->naddedcols);
   freeBlockMemoryArrayNull(memhdr, &(*fork)->addedrows, (*fork)->naddedrows);
   freeBlockMemory(memhdr, fork);

   return SCIP_OKAY;
}

/** creates subroot data */
static
RETCODE subrootCreate(
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

   ALLOC_OKAY( allocBlockMemory(memhdr, subroot) );

   CHECK_OKAY( SCIPlpGetState(lp, memhdr, &((*subroot)->lpistate)) );
   (*subroot)->nlpistateref = 0;
   (*subroot)->ncols = lp->ncols;
   (*subroot)->nrows = lp->nrows;
   (*subroot)->nchildren = tree->nchildren;
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*subroot)->cols, lp->cols, (*subroot)->ncols) );
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*subroot)->rows, lp->rows, (*subroot)->nrows) );

   /* capture the rows of the subroot */
   for( i = 0; i < (*subroot)->nrows; ++i )
      SCIProwCapture((*subroot)->rows[i]);

   /* capture the LPI state for the children */
   subrootCaptureLPIState(*subroot, tree->nchildren);
   
   return SCIP_OKAY;
}

/** frees subroot */
static
RETCODE subrootFree(
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
   assert((*subroot)->lpistate == NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   /* release the rows of the subroot */
   for( i = 0; i < (*subroot)->nrows; ++i )
   {
      CHECK_OKAY( SCIProwRelease(&(*subroot)->rows[i], memhdr, set, lp) );
   }

   freeBlockMemoryArrayNull(memhdr, &(*subroot)->cols, (*subroot)->ncols);
   freeBlockMemoryArrayNull(memhdr, &(*subroot)->rows, (*subroot)->nrows);
   freeBlockMemory(memhdr, subroot);

   return SCIP_OKAY;
}

/** makes node a child of the given parent node, which must be the active node */
static
RETCODE nodeAssignParent(
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
   assert(node->conssetchg == NULL);
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
   SCIPconssetchgdynAttach(tree->childrenconssetchg[tree->nchildren], &node->conssetchg);
   SCIPdomchgdynAttach(tree->childrendomchg[tree->nchildren], &node->domchg);
   node->data.child.arraypos = tree->nchildren;

   tree->nchildren++;

   return SCIP_OKAY;
}

/** decreases number of children of the parent, frees it if no children left */
static
RETCODE nodeReleaseParent(
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
         assert(parent->active);
         /* nothing to do here, because tree->nchildren is updated in treeRemoveChild() */
         hasChildren = TRUE; /* don't kill the active node at this point */
         break;
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
         assert(parent->data.junction.nchildren > 0);
         parent->data.junction.nchildren--;
         hasChildren = (parent->data.junction.nchildren > 0);
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
      {
         CHECK_OKAY( SCIPnodeFree(&node->parent, memhdr, set, tree, lp) );
      }
   }

   return SCIP_OKAY;
}

/** creates a child node of the active node */
RETCODE SCIPnodeCreate(
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

   ALLOC_OKAY( allocBlockMemory(memhdr, node) );
   (*node)->parent = NULL;
   (*node)->conssetchg = NULL;
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

/** removes given node from the siblings array */
static
void treeRemoveSibling(
   TREE*            tree,               /**< branch-and-bound tree */
   NODE*            sibling             /**< sibling node to remove */
   )
{
   CONSSETCHGDYN* conssetchgdyn;
   DOMCHGDYN* domchgdyn;
   int delpos;

   assert(tree != NULL);
   assert(sibling != NULL);
   assert(sibling->nodetype == SCIP_NODETYPE_SIBLING);
   assert(sibling->data.sibling.arraypos >= 0 && sibling->data.sibling.arraypos < tree->nsiblings);
   assert(tree->siblings[sibling->data.sibling.arraypos] == sibling);
   assert(tree->siblings[tree->nsiblings-1]->nodetype == SCIP_NODETYPE_SIBLING);

   delpos = sibling->data.sibling.arraypos;

   /* switch constraint set change data of removed sibling and last sibling in array */
   conssetchgdyn = tree->siblingsconssetchg[delpos];
   tree->siblingsconssetchg[delpos] = tree->siblingsconssetchg[tree->nsiblings-1];
   tree->siblingsconssetchg[tree->nsiblings-1] = conssetchgdyn;

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

/** removes given node from the children array */
static
void treeRemoveChild(
   TREE*            tree,               /**< branch-and-bound tree */
   NODE*            child               /**< child node to remove */
   )
{
   CONSSETCHGDYN* conssetchgdyn;
   DOMCHGDYN* domchgdyn;
   int delpos;

   assert(tree != NULL);
   assert(child != NULL);
   assert(child->nodetype == SCIP_NODETYPE_CHILD);
   assert(child->data.child.arraypos >= 0 && child->data.child.arraypos < tree->nchildren);
   assert(tree->children[child->data.child.arraypos] == child);
   assert(tree->children[tree->nchildren-1]->nodetype == SCIP_NODETYPE_CHILD);

   delpos = child->data.child.arraypos;

   /* switch constraint set change data of removed child and last child in array */
   conssetchgdyn = tree->childrenconssetchg[delpos];
   tree->childrenconssetchg[delpos] = tree->childrenconssetchg[tree->nchildren-1];
   tree->childrenconssetchg[tree->nchildren-1] = conssetchgdyn;

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

/** frees node */
RETCODE SCIPnodeFree(
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
      CHECK_OKAY( SCIPconssetchgdynDiscard(tree->actnodeconssetchg, memhdr, set) );
      CHECK_OKAY( SCIPdomchgdynDiscard(tree->actnodedomchg, memhdr) );
      if( tree->actlpfork != NULL )
      {
         assert(tree->actlpfork->nodetype == SCIP_NODETYPE_FORK || tree->actlpfork->nodetype == SCIP_NODETYPE_SUBROOT);
         CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
      }
      break;
   case SCIP_NODETYPE_SIBLING:
      assert((*node)->data.sibling.arraypos >= 0);
      assert((*node)->data.sibling.arraypos < tree->nsiblings);
      CHECK_OKAY( SCIPconssetchgdynDiscard(tree->siblingsconssetchg[(*node)->data.sibling.arraypos], memhdr, set) );
      CHECK_OKAY( SCIPdomchgdynDiscard(tree->siblingsdomchg[(*node)->data.sibling.arraypos], memhdr) );
      if( tree->actlpfork != NULL )
      {
         assert(tree->actlpfork->nodetype == SCIP_NODETYPE_FORK || tree->actlpfork->nodetype == SCIP_NODETYPE_SUBROOT);
         CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
      }
      treeRemoveSibling(tree, *node);
      break;
   case SCIP_NODETYPE_CHILD:
      assert((*node)->data.child.arraypos >= 0);
      assert((*node)->data.child.arraypos < tree->nchildren);
      CHECK_OKAY( SCIPconssetchgdynDiscard(tree->childrenconssetchg[(*node)->data.child.arraypos], memhdr, set) );
      CHECK_OKAY( SCIPdomchgdynDiscard(tree->childrendomchg[(*node)->data.child.arraypos], memhdr) );
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
      /*CHECK_OKAY( junctionFree(&((*node)->data.junction), memhdr) );*/
      break;
   case SCIP_NODETYPE_FORK:
      CHECK_OKAY( forkFree(&((*node)->data.fork), memhdr, set, lp) );
      break;
   case SCIP_NODETYPE_SUBROOT:
      CHECK_OKAY( subrootFree(&((*node)->data.subroot), memhdr, set, lp) );
      break;
   default:
      errorMessage("Unknown node type");
      abort();
   }

   /* free common data */
   CHECK_OKAY( SCIPconssetchgFree(&(*node)->conssetchg, memhdr, set) );
   CHECK_OKAY( SCIPdomchgFree(&(*node)->domchg, memhdr) );
   CHECK_OKAY( nodeReleaseParent(*node, memhdr, set, tree, lp) );

   freeBlockMemory(memhdr, node);

   return SCIP_OKAY;
}

/** informs node, that it is no longer on the active path */
static
RETCODE nodeDeactivate(
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
      hasChildren = ((*node)->data.junction.nchildren > 0);
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

/** adds local constraint to the node and captures it */
RETCODE SCIPnodeAddCons(
   NODE*            node,               /**< node to add constraint to */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(node != NULL);
   assert(cons != NULL);
   assert(cons->node == NULL);
   assert(cons->arraypos == -1);

   switch( node->nodetype )
   {
   case SCIP_NODETYPE_ACTNODE:
      assert(tree->actnode == node);
      CHECK_OKAY( SCIPconssetchgdynAddAddedCons(tree->actnodeconssetchg, memhdr, set, node, cons) );
      CHECK_OKAY( SCIPconsActivate(cons, set) );
      break;

   case SCIP_NODETYPE_SIBLING:
      assert(node->data.sibling.arraypos >= 0 && node->data.sibling.arraypos < tree->nsiblings);
      assert(tree->siblings[node->data.sibling.arraypos] == node);
      CHECK_OKAY( SCIPconssetchgdynAddAddedCons(tree->siblingsconssetchg[node->data.sibling.arraypos],
                     memhdr, set, node, cons) );
      break;

   case SCIP_NODETYPE_CHILD:
      assert(node->data.child.arraypos >= 0 && node->data.child.arraypos < tree->nchildren);
      assert(tree->children[node->data.child.arraypos] == node);
      CHECK_OKAY( SCIPconssetchgdynAddAddedCons(tree->childrenconssetchg[node->data.child.arraypos], 
                     memhdr, set, node, cons) );
      break;

   case SCIP_NODETYPE_LEAF:
   case SCIP_NODETYPE_DEADEND:
   case SCIP_NODETYPE_JUNCTION:
   case SCIP_NODETYPE_FORK:
   case SCIP_NODETYPE_SUBROOT:
      errorMessage("cannot add constraints to nodes stored in the tree");
      return SCIP_INVALIDDATA;
   default:
      errorMessage("unknown node type");
      return SCIP_ERROR;
   }
   assert(node->conssetchg != NULL);
   assert(node->conssetchg->addedconss != NULL);
   assert(cons->node == node);
   assert(0 <= cons->arraypos && cons->arraypos < node->conssetchg->naddedconss);
   assert(node->conssetchg->addedconss[cons->arraypos] == cons);

   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities at the node, and captures constraint */
RETCODE SCIPnodeDisableCons(
   NODE*            node,               /**< node to add constraint to */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   CONS*            cons                /**< constraint to disable */
   )
{
   assert(node != NULL);
   assert(tree != NULL);
   assert(cons != NULL);

   switch( node->nodetype )
   {
   case SCIP_NODETYPE_ACTNODE:
      assert(node == tree->actnode);
      CHECK_OKAY( SCIPconssetchgdynAddDisabledCons(tree->actnodeconssetchg, memhdr, set, cons) );
      CHECK_OKAY( SCIPconsDisable(cons, set) );
      return SCIP_OKAY;

   case SCIP_NODETYPE_SIBLING:
      assert(node->data.sibling.arraypos >= 0 && node->data.sibling.arraypos < tree->nsiblings);
      assert(tree->siblings[node->data.sibling.arraypos] == node);
      CHECK_OKAY( SCIPconssetchgdynAddDisabledCons(tree->siblingsconssetchg[node->data.sibling.arraypos],
                     memhdr, set, cons) );
      return SCIP_OKAY;

   case SCIP_NODETYPE_CHILD:
      assert(node->data.child.arraypos >= 0 && node->data.child.arraypos < tree->nchildren);
      assert(tree->children[node->data.child.arraypos] == node);
      CHECK_OKAY( SCIPconssetchgdynAddDisabledCons(tree->childrenconssetchg[node->data.child.arraypos], 
                     memhdr, set, cons) );
      return SCIP_OKAY;

   case SCIP_NODETYPE_LEAF:
   case SCIP_NODETYPE_DEADEND:
   case SCIP_NODETYPE_JUNCTION:
   case SCIP_NODETYPE_FORK:
   case SCIP_NODETYPE_SUBROOT:
      errorMessage("cannot disable constraints in nodes stored in the tree");
      return SCIP_INVALIDDATA;
   default:
      errorMessage("unknown node type");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** adds bound change to active node, child or sibling of active node */
RETCODE SCIPnodeAddBoundchg(
   NODE*            node,               /**< node to add bound change to */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
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
   {
      oldbound = var->dom.lb;

      if( SCIPsetIsLE(set, newbound, oldbound) )
      {
         char s[255];
         sprintf(s, "variable's lower bound was not tightened: var <%s>, oldbound=%f, newbound=%f",
            var->name, oldbound, newbound);
         errorMessage(s);
         return SCIP_INVALIDDATA;
      }
   }
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      oldbound = var->dom.ub;

      if( SCIPsetIsGE(set, newbound, oldbound) )
      {
         char s[255];
         sprintf(s, "variable's upper bound was not tightened: var <%s>, oldbound=%f, newbound=%f",
            var->name, oldbound, newbound);
         errorMessage(s);
         return SCIP_INVALIDDATA;
      }
   }

   switch( node->nodetype )
   {
   case SCIP_NODETYPE_ACTNODE:
      assert(tree->actnode == node);
      CHECK_OKAY( SCIPdomchgdynAddBoundchg(tree->actnodedomchg, memhdr, set, var, newbound, oldbound, boundtype) );
      CHECK_OKAY( SCIPvarChgBd(var, memhdr, set, stat, lp, tree, branchcand, eventqueue, newbound, boundtype) );
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

      /* update the child's lower bound from changed pseudo solution */
      if( var->varstatus == SCIP_VARSTATUS_LOOSE || var->varstatus == SCIP_VARSTATUS_COLUMN )
      {
         Real pseudoobjval;

         pseudoobjval = tree->actpseudoobjval;
         if( var->obj > 0.0 && boundtype == SCIP_BOUNDTYPE_LOWER )
            pseudoobjval += (newbound - oldbound) * var->obj;
         else if( var->obj < 0.0 && boundtype == SCIP_BOUNDTYPE_UPPER )
            pseudoobjval += (newbound - oldbound) * var->obj;
         node->lowerbound = MAX(node->lowerbound, pseudoobjval);
      }

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

/** if given value is larger than the node's lower bound, sets the node's lower bound to the new value */
void SCIPnodeUpdateLowerBound(
   NODE*            node,               /**< node to update lower bound for */
   Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   )
{
   assert(node != NULL);

   node->lowerbound = MAX(node->lowerbound, newbound);
}


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets the type of the node */
NODETYPE SCIPnodeGetType(
   NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->nodetype;
}

/** gets the depth of the node */
int SCIPnodeGetDepth(
   NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->depth;
}

/** gets the lower bound of the node */
Real SCIPnodeGetLowerBound(
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

/** updates the LP sizes of the active path starting at the given depth */
static
void treeUpdatePathLPSize(
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

/** cuts off path of active nodes after given node, marks cutted nodes inactive */
static
RETCODE treeShrinkPath(
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

/** switches the active path to end at the given node, applies domain and constraint set changes */
static
RETCODE treeSwitchPath(
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   NODE*            node                /**< last node of the new active path (= new active node), or NULL */
   )
{
   NODE* commonfork;    /* common fork node */
   NODE* lpfork;        /* subroot/fork node defining the LP state of the new active node */
   NODE* subroot;       /* subroot of new active path */
   int nodedepth;       /* depth of the node, or -1 of node == NULL */
   int commonforkdepth; /* depth of the common subroot/fork/junction node, or -1 if no common fork exists */
   int lpforkdepth;     /* depth of the common subroot/fork node, or -1 if no common lp fork exists */
   int subrootdepth;    /* depth of the common subroot node, or -1 if no common subroot exists */
   int i;

   assert(tree != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(node == NULL || !node->active);
   assert(node == NULL || node->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(tree->actlpfork == NULL || tree->actlpfork->active);
   assert(tree->actlpfork == NULL
      || tree->actlpfork->nodetype == SCIP_NODETYPE_FORK || tree->actlpfork->nodetype == SCIP_NODETYPE_SUBROOT);
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
   commonforkdepth = (commonfork == NULL ? -1 : commonfork->depth);
   assert(lpfork == NULL || !lpfork->active || lpfork == commonfork);
   assert(subroot == NULL || !subroot->active || subroot == commonfork);
   debugMessage("switch path: commonforkdepth=%d\n", commonforkdepth);

   /* if not already found, continue searching the LP defining fork; it can not be deeper than the common fork */
   if( lpfork == NULL )
   {
      if( tree->actlpfork != NULL && tree->actlpfork->depth > commonforkdepth )
      {
         /* actlpfork is not on the same active path as the new node: we have to search again */
         lpfork = commonfork;
         while( lpfork != NULL && lpfork->nodetype != SCIP_NODETYPE_FORK && lpfork->nodetype != SCIP_NODETYPE_SUBROOT )
            lpfork = lpfork->parent;
      }
      else
      {
         /* actlpfork is on the same active path as the new node: old and new node have the same lpfork */
         lpfork = tree->actlpfork;
      }
      assert(lpfork == NULL || lpfork->depth <= commonforkdepth);
   }
   assert(lpfork == NULL || lpfork->nodetype == SCIP_NODETYPE_FORK || lpfork->nodetype == SCIP_NODETYPE_SUBROOT);
   lpforkdepth = (lpfork == NULL ? -1 : lpfork->depth);
   debugMessage("switch path: lpforkdepth=%d\n", lpforkdepth);

   /* if not already found, continue searching the subroot; it cannot be deeper than the LP fork and common fork */
   if( subroot == NULL )
   {
      if( tree->actsubroot != NULL && tree->actsubroot->depth > commonforkdepth )
      {
         if( lpforkdepth < commonforkdepth )
            subroot = lpfork;
         else
            subroot = commonfork;
         while( subroot != NULL && subroot->nodetype != SCIP_NODETYPE_SUBROOT )
            subroot = subroot->parent;
      }
      else
         subroot = tree->actsubroot;
   }
   assert(subroot == NULL || subroot->nodetype == SCIP_NODETYPE_SUBROOT);
   subrootdepth = (subroot == NULL ? -1 : subroot->depth);
   debugMessage("switch path: subrootdepth=%d\n", subrootdepth);
   assert(subrootdepth <= lpforkdepth);
   debugMessage("switch path: old correctlpdepth=%d\n", tree->correctlpdepth);

   /* remember the depth of the common fork node for LP updates */
   if( subroot == tree->actsubroot )
   {
      /* we are in the same subtree: the LP is correct at most upto the common fork depth */
      assert(subroot == NULL || subroot->active);
      tree->correctlpdepth = MIN(tree->correctlpdepth, commonforkdepth);
   }
   else
   {
      /* we are in a different subtree: the LP is completely incorrect */
      assert(tree->actsubroot != NULL);
      assert(subroot == NULL || !subroot->active || tree->actsubroot->depth > subrootdepth);
      tree->correctlpdepth = -1;
   }
   debugMessage("switch path: new correctlpdepth=%d\n", tree->correctlpdepth);
   debugMessage("switch path: pathlen=%d\n", tree->pathlen);   

   /* undo the domain and constraint set changes of the old active path */
   for( i = tree->pathlen-1; i > commonforkdepth; --i )
   {
      debugMessage("switch path: undo domain changes in depth %d\n", i);
      CHECK_OKAY( SCIPdomchgUndo(tree->path[i]->domchg, memhdr, set, stat, lp, tree, branchcand, eventqueue) );
      debugMessage("switch path: undo constraint set changed in depth %d\n", i);
      CHECK_OKAY( SCIPconssetchgUndo(tree->path[i]->conssetchg, set) );
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

   /* apply domain and constraint set changes of the new path */
   for( i = commonforkdepth+1; i < tree->pathlen; ++i )
   {
      debugMessage("switch path: apply constraint set changed in depth %d\n", i);
      CHECK_OKAY( SCIPconssetchgApply(tree->path[i]->conssetchg, memhdr, set) );
      debugMessage("switch path: apply domain changes in depth %d\n", i);
      CHECK_OKAY( SCIPdomchgApply(tree->path[i]->domchg, memhdr, set, stat, lp, tree, branchcand, eventqueue) );
   }

   /* remember LP defining fork and subroot */
   assert(subroot == NULL || lpfork != NULL);
   assert(subroot == NULL || subroot->depth <= lpfork->depth);
   tree->actlpfork = lpfork;
   tree->actsubroot = subroot;
   
   return SCIP_OKAY;
}

/** loads the subroot's LP data */
static
RETCODE subrootConstructLP(
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
   
/** loads the fork's additional LP data */
static
RETCODE forkAddLP(
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

#ifndef NDEBUG
/** checks validity of active path */
static
void treeCheckPath(
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   NODE* node;
   int ncols;
   int nrows;
   int d;
   char s[255];

   assert(tree != NULL);
   assert(tree->path != NULL);

   ncols = 0;
   nrows = 0;
   for( d = 0; d < tree->pathlen; ++d )
   {
      node = tree->path[d];
      assert(node != NULL);
      assert(node->depth == d);
      switch( node->nodetype )
      {
      case SCIP_NODETYPE_JUNCTION:
         break;
      case SCIP_NODETYPE_FORK:
         ncols += node->data.fork->naddedcols;
         nrows += node->data.fork->naddedrows;
         break;
      case SCIP_NODETYPE_SUBROOT:
         ncols = node->data.subroot->ncols;
         nrows = node->data.subroot->nrows;
         break;
      case SCIP_NODETYPE_ACTNODE:
         assert(d == tree->pathlen-1);
         break;
      default:
         sprintf(s, "node in depth %d on active path has to be of type FORK, SUBROOT, or ACTNODE, but is %d",
            d, node->nodetype);
         errorMessage(s);
         abort();
      }
      assert(tree->pathnlpcols[d] == ncols);
      assert(tree->pathnlprows[d] == nrows);
   }
}
#else
#define treeCheckPath(tree) /**/
#endif

/** constructs the LP and loads LP state for fork/subroot of the active node */
RETCODE SCIPtreeLoadLP(
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
   assert(tree->pathlen > 0);
   assert(tree->actnode != NULL);
   assert(tree->actnode->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(tree->actnode->depth == tree->pathlen-1);
   assert(tree->actnode == tree->path[tree->pathlen-1]);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   debugMessage("load LP for actual fork node %p at depth %d\n", 
      tree->actlpfork, tree->actlpfork == NULL ? -1 : tree->actlpfork->depth);
   debugMessage("-> old LP has %d cols and %d rows\n", lp->ncols, lp->nrows);
   debugMessage("-> correct LP has %d cols and %d rows\n", 
      tree->correctlpdepth >= 0 ? tree->pathnlpcols[tree->correctlpdepth] : 0,
      tree->correctlpdepth >= 0 ? tree->pathnlprows[tree->correctlpdepth] : 0);
   debugMessage("-> old correctlpdepth: %d\n", tree->correctlpdepth);

   treeCheckPath(tree);

   lpfork = tree->actlpfork;

   /* find out the lpfork's depth (or -1, if lpfork is NULL) */
   if( lpfork == NULL )
   {
      assert(tree->correctlpdepth == -1 || tree->pathnlpcols[tree->correctlpdepth] == 0);
      assert(tree->correctlpdepth == -1 || tree->pathnlprows[tree->correctlpdepth] == 0);
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
      assert(lpforkdepth == -1 || tree->pathnlpcols[tree->correctlpdepth] <= tree->pathnlpcols[lpforkdepth]);
      assert(lpforkdepth == -1 || tree->pathnlprows[tree->correctlpdepth] <= tree->pathnlprows[lpforkdepth]);
      assert(lpforkdepth >= 0 || tree->pathnlpcols[tree->correctlpdepth] == 0);
      assert(lpforkdepth >= 0 || tree->pathnlprows[tree->correctlpdepth] == 0);
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
   tree->correctlpdepth = MAX(tree->correctlpdepth, lpforkdepth);
   assert(lpforkdepth == -1 || tree->pathnlpcols[tree->correctlpdepth] == tree->pathnlpcols[lpforkdepth]);
   assert(lpforkdepth == -1 || tree->pathnlprows[tree->correctlpdepth] == tree->pathnlprows[lpforkdepth]);
   assert(lpforkdepth == -1 || lp->ncols == tree->pathnlpcols[lpforkdepth]);
   assert(lpforkdepth == -1 || lp->nrows == tree->pathnlprows[lpforkdepth]);
   assert(lpforkdepth >= 0 || lp->ncols == 0);
   assert(lpforkdepth >= 0 || lp->nrows == 0);

   /* load LP state, if existing */
   if( lpfork != NULL )
   {
      if( lpfork->nodetype == SCIP_NODETYPE_FORK )
      {
         assert(lpfork->data.fork != NULL);
         CHECK_OKAY( SCIPlpSetState(lp, memhdr, set, lpfork->data.fork->lpistate) );
      }
      else
      {
         assert(lpfork->nodetype == SCIP_NODETYPE_SUBROOT);
         assert(lpfork->data.subroot != NULL);
         CHECK_OKAY( SCIPlpSetState(lp, memhdr, set, lpfork->data.subroot->lpistate) );
      }
      assert(lp->primalfeasible);
      assert(lp->dualfeasible);

      /* check the path from LP fork to active node for domain changes (destroying primal feasibility of LP basis) */
      for( d = lpforkdepth; d < tree->actnode->depth && lp->primalfeasible; ++d )
      {
         assert(d < tree->pathlen);
         lp->primalfeasible &= (tree->path[d]->domchg == NULL || tree->path[d]->domchg->nboundchg == 0);
      }
   }

   /* mark the LP's size, such that we know which rows and columns were added in the new node */
   SCIPlpMarkSize(lp);

   debugMessage("-> new correctlpdepth: %d\n", tree->correctlpdepth);
   debugMessage("-> new LP has %d cols and %d rows, primalfeasible=%d, dualfeasible=%d\n", 
      lp->ncols, lp->nrows, lp->primalfeasible, lp->dualfeasible);

   return SCIP_OKAY;
}




/*
 * Node Conversion
 */

/** converts node into LEAF and puts it in the array on the node queue */
static
RETCODE nodeToLeaf(
   NODE*            node,               /**< child or sibling node to convert */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   NODE*            lpfork              /**< LP fork of the node */
   )
{
   assert(node->nodetype == SCIP_NODETYPE_SIBLING || node->nodetype == SCIP_NODETYPE_CHILD);
   assert(lpfork == NULL || lpfork->depth < node->depth);
   assert(lpfork == NULL || lpfork->active);
   assert(lpfork == NULL || lpfork->nodetype == SCIP_NODETYPE_FORK || lpfork->nodetype == SCIP_NODETYPE_SUBROOT);

#if 0
   /* find the lpfork for the node */
   lpfork = tree->actlpfork;
   if( node->parent != NULL )
   {
      if( node->parent->nodetype == SCIP_NODETYPE_FORK || node->parent->nodetype == SCIP_NODETYPE_SUBROOT )
         lpfork = node->parent;
   }
#endif
   
   /* convert node into leaf */
   debugMessage("convert node %p at depth %d to leaf with lpfork %p at depth %d\n",
      node, node->depth, lpfork, lpfork == NULL ? -1 : lpfork->depth);
   node->nodetype = SCIP_NODETYPE_LEAF;
   node->data.leaf.lpfork = lpfork;

   /* insert leaf in node queue */
   CHECK_OKAY( SCIPnodepqInsert(tree->leaves, set, node) );

   return SCIP_OKAY;
}

/** converts the active node into a deadend node */
static
RETCODE actnodeToDeadend(
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
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

   debugMessage("actnode to deadend at depth %d\n", tree->actnode->depth);

   /* detach dynamic size data from constraint set and domain change data of old active node */
   CHECK_OKAY( SCIPconssetchgdynDetach(tree->actnodeconssetchg, memhdr, set) );
   CHECK_OKAY( SCIPdomchgdynDetach(tree->actnodedomchg, memhdr) );

   tree->actnode->nodetype = SCIP_NODETYPE_DEADEND;

   /* release LPI state */
   if( tree->actlpfork != NULL )
   {
      CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
   }

   return SCIP_OKAY;
}

/** converts the active node into a junction node */
static
RETCODE actnodeToJunction(
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnode->nodetype == SCIP_NODETYPE_ACTNODE);
   assert(tree->actnode->active);

   debugMessage("actnode to junction at depth %d\n", tree->actnode->depth);

   /* detach dynamic size data from constraint set and domain change data of old active node */
   CHECK_OKAY( SCIPconssetchgdynDetach(tree->actnodeconssetchg, memhdr, set) );
   CHECK_OKAY( SCIPdomchgdynDetach(tree->actnodedomchg, memhdr) );

   tree->actnode->nodetype = SCIP_NODETYPE_JUNCTION;

   CHECK_OKAY( junctionInit(&tree->actnode->data.junction, tree) );

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, tree->actnode->depth);

   /* release LPI state */
   if( tree->actlpfork != NULL )
   {
      CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
   }

   return SCIP_OKAY;
}

/** converts the active node into a fork node */
static
RETCODE actnodeToFork(
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
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

   debugMessage("actnode %p to fork at depth %d\n", tree->actnode, tree->actnode->depth);

   /* clean up newly created part of LP to keep only necessary columns and rows */
   CHECK_OKAY( SCIPlpCleanupNew(lp, memhdr, set) );

   /* resolve LP after cleaning up */
   if( !lp->solved )
   {
      CHECK_OKAY( SCIPlpSolveDual(lp, memhdr, set, stat) );
   }
   assert(lp->solved);

   /* remember that this node is solved correctly */
   tree->correctlpdepth = tree->actnode->depth;

   /* detach dynamic size data from constraint set and domain change data of old active node */
   CHECK_OKAY( SCIPconssetchgdynDetach(tree->actnodeconssetchg, memhdr, set) );
   CHECK_OKAY( SCIPdomchgdynDetach(tree->actnodedomchg, memhdr) );

   CHECK_OKAY( forkCreate(&fork, memhdr, set, lp, tree) );
   
   tree->actnode->nodetype = SCIP_NODETYPE_FORK;
   tree->actnode->data.fork = fork;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, tree->actnode->depth);

   /* release LPI state */
   if( tree->actlpfork != NULL )
   {
      CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
   }

   /* set new actual LP fork */
   tree->actlpfork = tree->actnode;

   return SCIP_OKAY;
}

/** converts the active node into a subroot node */
static
RETCODE actnodeToSubroot(
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
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

   debugMessage("actnode %p to subroot at depth %d\n", tree->actnode, tree->actnode->depth);

   /* clean up whole LP to keep only necessary columns and rows */
   CHECK_OKAY( SCIPlpCleanupAll(lp, memhdr, set) );

   /* resolve LP after cleaning up */
   if( !lp->solved )
   {
      CHECK_OKAY( SCIPlpSolveDual(lp, memhdr, set, stat) );
   }
   assert(lp->solved);

   /* remember that this node is solved correctly */
   tree->correctlpdepth = tree->actnode->depth;

   /* detach dynamic size data from constraint set and domain change data of old active node */
   CHECK_OKAY( SCIPconssetchgdynDetach(tree->actnodeconssetchg, memhdr, set) );
   CHECK_OKAY( SCIPdomchgdynDetach(tree->actnodedomchg, memhdr) );

   CHECK_OKAY( subrootCreate(&subroot, memhdr, set, lp, tree) );

   tree->actnode->nodetype = SCIP_NODETYPE_SUBROOT;
   tree->actnode->data.subroot = subroot;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, tree->actnode->depth);

   /* release LPI state */
   if( tree->actlpfork != NULL )
   {
      CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
   }

   /* set new actual LP fork and actual subroot */
   tree->actlpfork = tree->actnode;
   tree->actsubroot = tree->actnode;

   return SCIP_OKAY;
}

/** puts all nodes in the array on the node queue and makes them LEAFs */
static
RETCODE treeNodesToQueue(
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   NODE**           nodes,              /**< array of nodes to put on the queue */
   CONSSETCHGDYN**  conssetchgdyn,      /**< array of dynamic constraint set changes */
   DOMCHGDYN**      domchgdyn,          /**< array of dynamic domain changes */
   int*             nnodes,             /**< pointer to number of nodes in the array */
   NODE*            lpfork              /**< LP fork of the nodes */
   )
{
   NODE* node;
   int i;

   assert(tree != NULL);
   assert(set != NULL);
   assert(nnodes != NULL);
   assert(*nnodes == 0 || nodes != NULL);
   assert(*nnodes == 0 || conssetchgdyn != NULL);
   assert(*nnodes == 0 || domchgdyn != NULL);

   for( i = 0; i < *nnodes; ++i )
   {
      node = nodes[i];

      /* detach dynamic size attachment to shrink the node's constraint set and domain change data */
      CHECK_OKAY( SCIPconssetchgdynDetach(conssetchgdyn[i], memhdr, set) );
      CHECK_OKAY( SCIPdomchgdynDetach(domchgdyn[i], memhdr) );

      /* convert node to LEAF and put it into leaves queue */
      nodeToLeaf(node, set, tree, lpfork);
   }
   *nnodes = 0;

   return SCIP_OKAY;
}

/** converts children into siblings, clears children array */
static
void treeChildrenToSiblings(
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   NODE** tmpnodes;
   CONSSETCHGDYN** tmpconssetchg;
   DOMCHGDYN** tmpdomchg;
   int tmpnodessize;
   int i;

   assert(tree != NULL);
   assert(tree->nsiblings == 0);

   tmpnodes = tree->siblings;
   tmpconssetchg = tree->siblingsconssetchg;
   tmpdomchg = tree->siblingsdomchg;
   tmpnodessize = tree->siblingssize;

   tree->siblings = tree->children;
   tree->siblingsconssetchg = tree->childrenconssetchg;
   tree->siblingsdomchg = tree->childrendomchg;
   tree->nsiblings = tree->nchildren;
   tree->siblingssize = tree->childrensize;

   tree->children = tmpnodes;
   tree->childrenconssetchg = tmpconssetchg;
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

/** activates a child, a sibling, or a leaf node */
RETCODE SCIPnodeActivate(
   NODE*            node,               /**< leaf node to activate (or NULL to deactivate all nodes) */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   NODE* oldlpfork;
   NODE* newlpfork;

   assert(node == NULL
      || node->nodetype == SCIP_NODETYPE_SIBLING
      || node->nodetype == SCIP_NODETYPE_CHILD
      || node->nodetype == SCIP_NODETYPE_LEAF);
   assert(node == NULL || !node->active);
   assert(memhdr != NULL);
   assert(tree != NULL);
   assert(lp != NULL);

   /* remember the actual LP fork, which is also the LP fork for the siblings of the old active node */
   oldlpfork = tree->actlpfork;

   /* deactivate old active node */
   if( tree->actnode != NULL )
   {
      assert(tree->actnode->nodetype == SCIP_NODETYPE_ACTNODE);

      /* convert the old active node into a fork node, if it has children, or a dead end otherwise */
      if( tree->nchildren > 0 )
      {
         if( tree->actnodehaslp )
         {
            todoMessage("decide: old active node becomes fork or subroot");
            if( tree->actnode->depth % 25 == 0 ) /* ????????????? */
            {
               /* convert old active node into a subroot node */
               CHECK_OKAY( actnodeToSubroot(memhdr, set, stat, tree, lp) );
            }
            else
            {
               /* convert old active node into a fork node */
               CHECK_OKAY( actnodeToFork(memhdr, set, stat, tree, lp) );
            }
         }
         else
         {
            /* convert old active node into junction */
            CHECK_OKAY( actnodeToJunction(memhdr, set, tree, lp) );
         }
      }
      else
      {
         CHECK_OKAY( actnodeToDeadend(memhdr, set, tree, lp) );
      }
   }

   /* now the old active node was converted to a subroot, fork, or junction; the first two cases made the old active
    * node the new LP fork, which is the correct LP fork for the child nodes of the old active node
    */
   newlpfork = tree->actlpfork;

   /* set up the new lists of siblings and children */
   if( node == NULL )
   {
      /* move siblings to the queue, make them LEAFs */
      CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, tree->siblings, tree->siblingsconssetchg, tree->siblingsdomchg,
                     &tree->nsiblings, oldlpfork) );

      /* move children to the queue, make them LEAFs */
      CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, tree->children, tree->childrenconssetchg, tree->childrendomchg,
                     &tree->nchildren, newlpfork) );
   }
   else
   {
      CONSSETCHGDYN* conssetchgdyn;
      DOMCHGDYN* domchgdyn;

      switch( node->nodetype )
      {
      case SCIP_NODETYPE_SIBLING:
         /* move children to the queue, make them LEAFs */
         CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, tree->children, tree->childrenconssetchg, tree->childrendomchg,
                        &tree->nchildren, newlpfork) );

         /* switch constraint set change data of sibling and active node */
         conssetchgdyn = tree->actnodeconssetchg;
         tree->actnodeconssetchg = tree->siblingsconssetchg[node->data.sibling.arraypos];
         tree->siblingsconssetchg[node->data.sibling.arraypos] = conssetchgdyn;

         /* switch domain change data of sibling and active node */
         domchgdyn = tree->actnodedomchg;
         tree->actnodedomchg = tree->siblingsdomchg[node->data.sibling.arraypos];
         tree->siblingsdomchg[node->data.sibling.arraypos] = domchgdyn;

         /* remove selected sibling from the siblings array */
         treeRemoveSibling(tree, node);

         /* reinstall the old LP fork, because this is the correct LP fork of the sibling */
         tree->actlpfork = oldlpfork;
         
         break;
         
      case SCIP_NODETYPE_CHILD:
         /* move siblings to the queue, make them LEAFs */
         CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, tree->siblings, tree->siblingsconssetchg, tree->siblingsdomchg,
                        &tree->nsiblings, oldlpfork) );

         /* switch constraint set change data of child and active node */
         conssetchgdyn = tree->actnodeconssetchg;
         tree->actnodeconssetchg = tree->childrenconssetchg[node->data.child.arraypos];
         tree->childrenconssetchg[node->data.child.arraypos] = conssetchgdyn;

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
         CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, tree->siblings, tree->siblingsconssetchg, tree->siblingsdomchg,
                        &tree->nsiblings, oldlpfork) );
         
         /* move children to the queue, make them LEAFs */
         CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, tree->children, tree->childrenconssetchg, tree->childrendomchg,
                        &tree->nchildren, newlpfork) );

         /* attach dynamic size data to constraint set and domain changes of the active node */
         SCIPconssetchgdynAttach(tree->actnodeconssetchg, &node->conssetchg);
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
   CHECK_OKAY( treeSwitchPath(tree, memhdr, set, stat, lp, branchcand, eventqueue, node) );
   assert(node == NULL || tree->pathlen > 0);
   assert(node != NULL || tree->pathlen == 0);
   assert(node == NULL || tree->path[tree->pathlen-1] == node);
   assert(tree->nchildren == 0);
   tree->actnode = node;

#ifndef NDEBUG
   {
      int i;
      for( i = 0; i < tree->nchildren; ++i )
      {
         assert(SCIPconssetchgdynGetConssetchgPtr(tree->childrenconssetchg[i]) == &tree->children[i]->conssetchg);
         assert(SCIPdomchgdynGetDomchgPtr(tree->childrendomchg[i]) == &tree->children[i]->domchg);
      }
      for( i = 0; i < tree->nsiblings; ++i )
      {
         assert(SCIPconssetchgdynGetConssetchgPtr(tree->siblingsconssetchg[i]) == &tree->siblings[i]->conssetchg);
         assert(SCIPdomchgdynGetDomchgPtr(tree->siblingsdomchg[i]) == &tree->siblings[i]->domchg);
      }
      if( tree->actnode != NULL )
      {
         assert(SCIPconssetchgdynGetConssetchgPtr(tree->actnodeconssetchg) == &tree->actnode->conssetchg);
         assert(SCIPdomchgdynGetDomchgPtr(tree->actnodedomchg) == &tree->actnode->domchg);
      }
      else
      {
         assert(SCIPconssetchgdynGetConssetchgPtr(tree->actnodeconssetchg) == NULL);
         assert(SCIPdomchgdynGetDomchgPtr(tree->actnodedomchg) == NULL);
      }
   }
#endif

   return SCIP_OKAY;
}   




/*
 * Tree methods
 */

/** creates an initialized tree data structure */
RETCODE SCIPtreeCreate(
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
   assert(lp != NULL);
   assert(prob != NULL);

   ALLOC_OKAY( allocMemory(tree) );

   (*tree)->root = NULL;

   CHECK_OKAY( SCIPnodepqCreate(&(*tree)->leaves, set) );

   (*tree)->path = NULL;
   (*tree)->actnode = NULL;
   (*tree)->actlpfork = NULL;
   (*tree)->actsubroot = NULL;
   (*tree)->children = NULL;
   (*tree)->siblings = NULL;

   CHECK_OKAY( SCIPconssetchgdynCreate(&(*tree)->actnodeconssetchg, memhdr) );
   (*tree)->childrenconssetchg = NULL;
   (*tree)->siblingsconssetchg = NULL;

   CHECK_OKAY( SCIPdomchgdynCreate(&(*tree)->actnodedomchg, memhdr) );
   (*tree)->childrendomchg = NULL;
   (*tree)->siblingsdomchg = NULL;

   /* calculate root pseudo solution value */
   (*tree)->actpseudoobjval = 0.0;
   for( v = 0; v < prob->nvars; ++v )
   {
      var = prob->vars[v];
      (*tree)->actpseudoobjval += SCIPvarGetPseudoSol(var) * var->obj;
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
   (*tree)->actnodehaslp = FALSE;

   /* create root node */
   CHECK_OKAY( SCIPnodeCreate(&(*tree)->root, memhdr, set, *tree) );
   assert((*tree)->nchildren == 1);
   /* move root to the queue, convert it to LEAF */
   CHECK_OKAY( treeNodesToQueue(*tree, memhdr, set, (*tree)->children, (*tree)->childrenconssetchg, (*tree)->childrendomchg,
                  &(*tree)->nchildren, NULL) );

   return SCIP_OKAY;
}

/** frees tree data structure */
RETCODE SCIPtreeFree(
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

   /* free node queue */
   CHECK_OKAY( SCIPnodepqFree(&(*tree)->leaves, memhdr, set, *tree, lp) );
   
   /* free dynamic size constraint set and domain change attachments */
   for( i = 0; i < (*tree)->childrensize; ++i )
   {
      SCIPconssetchgdynFree(&(*tree)->childrenconssetchg[i], memhdr);
      SCIPdomchgdynFree(&(*tree)->childrendomchg[i], memhdr);
   }
   for( i = 0; i < (*tree)->siblingssize; ++i )
   {
      SCIPconssetchgdynFree(&(*tree)->siblingsconssetchg[i], memhdr);
      SCIPdomchgdynFree(&(*tree)->siblingsdomchg[i], memhdr);
   }
   SCIPconssetchgdynFree(&(*tree)->actnodeconssetchg, memhdr);
   SCIPdomchgdynFree(&(*tree)->actnodedomchg, memhdr);

   /* free pointer arrays */
   freeMemoryArrayNull(&(*tree)->path);
   freeMemoryArrayNull(&(*tree)->children);
   freeMemoryArrayNull(&(*tree)->siblings);   
   freeMemoryArrayNull(&(*tree)->childrenconssetchg);
   freeMemoryArrayNull(&(*tree)->siblingsconssetchg);   
   freeMemoryArrayNull(&(*tree)->childrendomchg);
   freeMemoryArrayNull(&(*tree)->siblingsdomchg);   
   freeMemoryArrayNull(&(*tree)->pathnlpcols);
   freeMemoryArrayNull(&(*tree)->pathnlprows);

   freeMemory(tree);

   return SCIP_OKAY;
}

/** cuts off nodes with lower bound not better than given upper bound */
RETCODE SCIPtreeCutoff(
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

/** branches on a variable; if solution value x' is fractional, two child nodes are created
 *  (x <= floor(x'), x >= ceil(x')), if solution value is integral, three child nodes are created
 *  (x <= x'-1, x == x', x >= x'+1)
 */
RETCODE SCIPtreeBranchVar(
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             var                 /**< variable to branch on */
   )
{
   NODE* node;
   Real solval;

   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_LOOSE || var->varstatus == SCIP_VARSTATUS_COLUMN);
   assert(var->vartype == SCIP_VARTYPE_BINARY
      || var->vartype == SCIP_VARTYPE_INTEGER
      || var->vartype == SCIP_VARTYPE_IMPLINT);
   assert(SCIPsetIsIntegral(set, var->dom.lb));
   assert(SCIPsetIsIntegral(set, var->dom.ub));
   assert(!SCIPsetIsFixed(set, var->dom.lb, var->dom.ub));

   solval = SCIPvarGetSol(var, tree);
   assert(SCIPsetIsGE(set, solval, var->dom.lb));
   assert(SCIPsetIsLE(set, solval, var->dom.ub));
   
   if( SCIPsetIsIntegral(set, solval) )
   {
      Real fixval;

      /* create child nodes with x <= x'-1, x = x', and x >= x'+1 */
      fixval = SCIPsetCeil(set, solval);
      assert(SCIPsetIsEQ(set, SCIPsetCeil(set, solval), SCIPsetFloor(set, solval)));
      
      debugMessage("pseudo branch on variable <%s> with value %g\n", var->name, solval);
      
      /* create child node with x = x' */
      debugMessage(" -> creating child: <%s> == %g\n", var->name, fixval);
      CHECK_OKAY( SCIPnodeCreate(&node, memhdr, set, tree) );
      if( !SCIPsetIsEQ(set, var->dom.lb, fixval) )
      {
         CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, lp, tree, branchcand, eventqueue, var, fixval,
                        SCIP_BOUNDTYPE_LOWER) );
      }
      if( !SCIPsetIsEQ(set, var->dom.ub, fixval) )
      {
         CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, lp, tree, branchcand, eventqueue, var, fixval,
                        SCIP_BOUNDTYPE_UPPER) );
      }
      debugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
      
      /* create child node with x <= x'-1, if this would be feasible */
      if( SCIPsetIsGE(set, fixval-1, var->dom.lb) )
      {
         debugMessage(" -> creating child: <%s> <= %g\n", var->name, fixval-1);
         CHECK_OKAY( SCIPnodeCreate(&node, memhdr, set, tree) );
         CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, lp, tree, branchcand, eventqueue, var, fixval-1,
                        SCIP_BOUNDTYPE_UPPER) );
         debugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
      }
                  
      /* create child node with x >= x'+1, if this would be feasible */
      if( SCIPsetIsLE(set, fixval+1, var->dom.ub) )
      {
         debugMessage(" -> creating child: <%s> >= %g\n", var->name, fixval+1);
         CHECK_OKAY( SCIPnodeCreate(&node, memhdr, set, tree) );
         CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, lp, tree, branchcand, eventqueue, var, fixval+1,
                        SCIP_BOUNDTYPE_LOWER) );
         debugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
      }
   }
   else
   {   
      debugMessage("LP branch on variable <%s> with value %g\n", var->name, solval);
      
      /* create child node with x <= floor(x') */
      debugMessage(" -> creating child: <%s> <= %g\n", var->name, SCIPsetFloor(set, solval));
      CHECK_OKAY( SCIPnodeCreate(&node, memhdr, set, tree) );
      CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, lp, tree, branchcand, eventqueue, var,
                     SCIPsetFloor(set, solval), SCIP_BOUNDTYPE_UPPER) );
      debugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
      
      /* create child node with x >= ceil(x') */
      debugMessage(" -> creating child: <%s> >= %g\n", var->name, SCIPsetCeil(set, solval));
      CHECK_OKAY( SCIPnodeCreate(&node, memhdr, set, tree) );
      CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, lp, tree, branchcand, eventqueue, var,
                     SCIPsetCeil(set, solval), SCIP_BOUNDTYPE_LOWER) );
      debugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
   }

   return SCIP_OKAY;
}

/** notifies tree, that a bound of a variable changed */
RETCODE SCIPtreeBoundChanged(
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   Real             oldbound,           /**< old bound value */
   Real             newbound            /**< new bound value */
   )
{
   assert(var != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(tree != NULL);
      if( var->obj > 0.0 && boundtype == SCIP_BOUNDTYPE_LOWER )
         tree->actpseudoobjval += (newbound - oldbound) * var->obj;
      else if( var->obj < 0.0 && boundtype == SCIP_BOUNDTYPE_UPPER )
         tree->actpseudoobjval += (newbound - oldbound) * var->obj;
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

/** gets number of leaves */
int SCIPtreeGetNLeaves(
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqLen(tree->leaves);
}
   
/** gets number of nodes (children + siblings + leaves) */
int SCIPtreeGetNNodes(
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(tree != NULL);

   return tree->nchildren + tree->nsiblings + SCIPtreeGetNLeaves(tree);
}

/** gets the best child of the active node */
NODE* SCIPtreeGetBestChild(
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

/** gets the best sibling of the active node */
NODE* SCIPtreeGetBestSibling(
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

/** gets the best leaf from the node queue */
NODE* SCIPtreeGetBestLeaf(
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqFirst(tree->leaves);
}

/** gets the best node from the tree (child, sibling, or leaf) */
NODE* SCIPtreeGetBestNode(
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

/** gets the minimal lower bound of all nodes in the tree */
Real SCIPtreeGetLowerbound(
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

/** gets the lower bound of the active node */
Real SCIPtreeGetActLowerbound(
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(tree != NULL);

   if( tree->actnode != NULL )
      return tree->actnode->lowerbound;
   else
      return SCIP_INVALID;
}

/** gets the average lower bound of all nodes in the tree */
Real SCIPtreeGetAvgLowerbound(
   TREE*            tree,               /**< branch-and-bound tree */
   Real             upperbound          /**< global upper bound */
   )
{
   Real lowerboundsum;
   int nnodes;
   int i;

   assert(tree != NULL);

   /* get sum of lower bounds from nodes in the queue */
   lowerboundsum = SCIPnodepqGetLowerboundSum(tree->leaves);
   nnodes = SCIPtreeGetNLeaves(tree);

   /* add lower bound of active node */
   if( tree->actnode != NULL && tree->actnode->lowerbound < upperbound )
   {
      lowerboundsum += tree->actnode->lowerbound;
      nnodes++;
   }

   /* add lower bounds of siblings */
   for( i = 0; i < tree->nsiblings; ++i )
   {
      assert(tree->siblings[i] != NULL);
      lowerboundsum += tree->siblings[i]->lowerbound;
   }
   nnodes += tree->nsiblings;

   /* add lower bounds of children */
   for( i = 0; i < tree->nchildren; ++i )
   {
      assert(tree->children[i] != NULL);
      lowerboundsum += tree->children[i]->lowerbound;
   }
   nnodes += tree->nchildren;

   return nnodes == 0 ? 0.0 : lowerboundsum/nnodes;
}
