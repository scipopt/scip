/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: tree.c,v 1.110 2004/09/15 08:11:28 bzfpfend Exp $"

/**@file   tree.c
 * @brief  methods for branch and bound tree
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "def.h"
#include "message.h"
#include "set.h"
#include "stat.h"
#include "vbc.h"
#include "lp.h"
#include "var.h"
#include "tree.h"
#include "cons.h"
#include "nodesel.h"


#define MAXDEPTH   65535  /**< maximal depth level for nodes, must correspond to node data structure */


/*
 * dynamic memory arrays
 */

/** resizes children arrays to be able to store at least num nodes */
static
RETCODE treeEnsureChildrenMem(
   TREE*            tree,               /**< branch and bound tree */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(tree != NULL);
   assert(set != NULL);

   if( num > tree->childrensize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&tree->children, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&tree->childrenprio, newsize) );
      tree->childrensize = newsize;
   }
   assert(num <= tree->childrensize);

   return SCIP_OKAY;
}

/** resizes path array to be able to store at least num nodes */
static
RETCODE treeEnsurePathMem(
   TREE*            tree,               /**< branch and bound tree */
   SET*             set,                /**< global SCIP settings */
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
DECL_SORTPTRCOMP(SCIPnodeCompLowerbound)
{  /*lint --e{715}*/
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
   LP*              lp                  /**< current LP data */
   )
{
   assert(fork != NULL);
   assert(fork->nlpistateref > 0);
   assert(fork->lpistate != NULL);
   assert(memhdr != NULL);
   assert(lp != NULL);

   fork->nlpistateref--;
   if( fork->nlpistateref == 0 )
   {
      CHECK_OKAY( SCIPlpFreeState(lp, memhdr, &(fork->lpistate)) );
   }

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
   assert(subroot->lpistate != NULL || subroot->nrows == 0 || subroot->ncols == 0);
   assert(nuses > 0);

   subroot->nlpistateref += nuses;
   debugMessage("captured subroot's LPI state %d times -> new nlpistateref=%d\n", nuses, subroot->nlpistateref);
}

/** decreases the reference counter of the LP state in the subroot */
static
RETCODE subrootReleaseLPIState(
   SUBROOT*         subroot,            /**< subroot data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   LP*              lp                  /**< current LP data */
   )
{
   assert(subroot != NULL);
   assert(subroot->nlpistateref > 0);
   assert(subroot->lpistate != NULL || subroot->nrows == 0 || subroot->ncols == 0);
   assert(memhdr != NULL);
   assert(lp != NULL);

   subroot->nlpistateref--;
   if( subroot->nlpistateref == 0 )
   {
      CHECK_OKAY( SCIPlpFreeState(lp, memhdr, &(subroot->lpistate)) );
   }
   
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

   debugMessage("capture %d times LPI state of node %p at depth %d (current: %d)\n", nuses, node, node->depth, 
      SCIPnodeGetType(node) == SCIP_NODETYPE_FORK ? node->data.fork->nlpistateref : node->data.subroot->nlpistateref);
   switch( SCIPnodeGetType(node) )
   {  
   case SCIP_NODETYPE_FORK:
      forkCaptureLPIState(node->data.fork, nuses);
      break;
   case SCIP_NODETYPE_SUBROOT:
      subrootCaptureLPIState(node->data.subroot, nuses);
      break;
   default:
      errorMessage("node for capturing the LPI state is neither fork nor subroot\n");
      abort();
   }  /*lint !e788*/
}

/** decreases the reference counter of the LP state in the fork or subroot node */
RETCODE SCIPnodeReleaseLPIState(
   NODE*            node,               /**< fork/subroot node */
   MEMHDR*          memhdr,             /**< block memory buffers */
   LP*              lp                  /**< current LP data */
   )
{
   assert(node != NULL);

   debugMessage("release LPI state of node %p at depth %d (current: %d)\n", node, node->depth,
      SCIPnodeGetType(node) == SCIP_NODETYPE_FORK ? node->data.fork->nlpistateref : node->data.subroot->nlpistateref);
   switch( SCIPnodeGetType(node) )
   {  
   case SCIP_NODETYPE_FORK:
      return forkReleaseLPIState(node->data.fork, memhdr, lp);
   case SCIP_NODETYPE_SUBROOT:
      return subrootReleaseLPIState(node->data.subroot, memhdr, lp);
   default:
      errorMessage("node for releasing the LPI state is neither fork nor subroot\n");
      return SCIP_INVALIDDATA;
   }  /*lint !e788*/
}

/** initializes junction data */
static
RETCODE junctionInit(
   JUNCTION*        junction,           /**< pointer to junction data */
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(junction != NULL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);

   junction->nchildren = tree->nchildren;

   /* increase the LPI state usage counter of the current LP fork */
   if( tree->actlpfork != NULL )
      SCIPnodeCaptureLPIState(tree->actlpfork, tree->nchildren);

   return SCIP_OKAY;
}

/** creates fork data */
static
RETCODE forkCreate(
   FORK**           fork,               /**< pointer to fork data */
   MEMHDR*          memhdr,             /**< block memory */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   )
{
   assert(fork != NULL);
   assert(memhdr != NULL);
   assert(lp != NULL);
   assert(lp->solved);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);

   ALLOC_OKAY( allocBlockMemory(memhdr, fork) );

   CHECK_OKAY( SCIPlpGetState(lp, memhdr, &((*fork)->lpistate)) );
   (*fork)->nlpistateref = 0;
   (*fork)->addedcols = NULL;
   (*fork)->addedrows = NULL;
   (*fork)->naddedcols = SCIPlpGetNNewcols(lp);
   (*fork)->naddedrows = SCIPlpGetNNewrows(lp);
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
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
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
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   )
{
   int i;
      
   assert(subroot != NULL);
   assert(memhdr != NULL);
   assert(lp != NULL);
   assert(lp->solved);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);

   ALLOC_OKAY( allocBlockMemory(memhdr, subroot) );
   
   (*subroot)->nlpistateref = 0;
   (*subroot)->ncols = SCIPlpGetNCols(lp);
   (*subroot)->nrows = SCIPlpGetNRows(lp);
   (*subroot)->nchildren = tree->nchildren;
   CHECK_OKAY( SCIPlpGetState(lp, memhdr, &((*subroot)->lpistate)) );

   if( (*subroot)->ncols != 0 )
   {
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*subroot)->cols, SCIPlpGetCols(lp), (*subroot)->ncols) );
   }
   else
      (*subroot)->cols = NULL;
   if( (*subroot)->nrows != 0 )
   {
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*subroot)->rows, SCIPlpGetRows(lp), (*subroot)->nrows) );
   }
   else
      (*subroot)->rows = NULL;
  
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
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
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
   SET*             set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   NODE*            parent,             /**< parent (= active) node (or NULL, if node is root) */
   Real             nodeselprio         /**< node selection priority of child node */
   )
{
   assert(node != NULL);
   assert(node->parent == NULL);
   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);
   assert(node->conssetchg == NULL);
   assert(node->domchg == NULL);
   assert(SCIPsetIsInfinity(set, -node->lowerbound)); /* node was just created */
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->actnode == parent);
   assert(parent == NULL || SCIPnodeGetType(parent) == SCIP_NODETYPE_ACTNODE);

   /* link node to parent */
   node->parent = parent;
   if( parent != NULL )
   {
      node->lowerbound = parent->lowerbound;
      node->depth = parent->depth+1;
      if( parent->depth >= MAXDEPTH-1 )
      {
         errorMessage("maximal depth level exceeded\n");
         return SCIP_MAXDEPTHLEVEL;
      }
   }
   debugMessage("assigning parent %p to node %p in depth %d\n", parent, node, node->depth);

   /* register node in the childlist of the active (the parent) node */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD )
   {
      assert(node->data.child.arraypos == -1);
      CHECK_OKAY( treeEnsureChildrenMem(tree, set, tree->nchildren+1) );
      tree->children[tree->nchildren] = node;
      tree->childrenprio[tree->nchildren] = nodeselprio;
      node->data.child.arraypos = tree->nchildren;
      tree->nchildren++;
   }

   return SCIP_OKAY;
}

/** decreases number of children of the parent, frees it if no children left */
static
RETCODE nodeReleaseParent(
   NODE*            node,               /**< child node */
   MEMHDR*          memhdr,             /**< block memory buffer */
   SET*             set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   )
{
   NODE* parent;
   Bool hasChildren = TRUE;

   assert(node != NULL);
   assert(memhdr != NULL);
   
   debugMessage("releasing parent-child relationship of node %p at depth %d of type %d with parent %p of type %d\n",
      node, node->depth, SCIPnodeGetType(node), node->parent, 
      node->parent != NULL ? (int)(SCIPnodeGetType(node->parent)) : -1);
   parent = node->parent;
   if( parent != NULL )
   {
      switch( SCIPnodeGetType(parent) )
      {
      case SCIP_NODETYPE_ACTNODE:
         assert(parent->active);
         /* nothing to do here, because tree->nchildren is updated in treeRemoveChild() */
         hasChildren = TRUE; /* don't kill the active node at this point */
         break;
      case SCIP_NODETYPE_PROBINGNODE:
         errorMessage("probing node cannot be a parent node\n");
         abort();
      case SCIP_NODETYPE_SIBLING:
         errorMessage("sibling cannot be a parent node\n");
         abort();
      case SCIP_NODETYPE_CHILD:
         errorMessage("child cannot be a parent node\n");
         abort();
      case SCIP_NODETYPE_LEAF:
         errorMessage("leaf cannot be a parent node\n");
         abort();
      case SCIP_NODETYPE_DEADEND:
         errorMessage("deadend cannot be a parent node\n");
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
         errorMessage("unknown node type\n");
         abort();
      }

      /* free parent, if it has no more children left and is not on the current active path */
      if( !hasChildren && !parent->active )
      {
         CHECK_OKAY( SCIPnodeFree(&node->parent, memhdr, set, tree, lp) );
      }
   }

   return SCIP_OKAY;
}

/** creates a node data structure */
static
RETCODE nodeCreate(
   NODE**           node,               /**< pointer to node data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(node != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, node) );
   (*node)->parent = NULL;
   (*node)->conssetchg = NULL;
   (*node)->domchg = NULL;
   (*node)->lowerbound = -set->infinity;
   (*node)->depth = 0;
   (*node)->active = FALSE;
   (*node)->cutoff = FALSE;

   return SCIP_OKAY;
}

/** creates a child node of the active node */
RETCODE SCIPnodeCreateChild(
   NODE**           node,               /**< pointer to node data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   Real             nodeselprio         /**< node selection priority of new node */
   )
{
   assert(node != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(tree->pathlen == 0 || tree->path != NULL);
   assert((tree->pathlen == 0) == (tree->actnode == NULL));
   assert(tree->actnode == NULL || tree->actnode == tree->path[tree->pathlen-1]);
   assert(tree->actnode == NULL || SCIPnodeGetType(tree->actnode) == SCIP_NODETYPE_ACTNODE);

   stat->ncreatednodes++;

   /* create the node data structure */
   CHECK_OKAY( nodeCreate(node, memhdr, set) );

   /* mark node to be a child node */
   (*node)->nodetype = SCIP_NODETYPE_CHILD; /*lint !e641*/
   (*node)->data.child.arraypos = -1;

   /* make active node the parent of the new child */
   CHECK_OKAY( nodeAssignParent(*node, memhdr, set, tree, tree->actnode, nodeselprio) );

   /* output node creation to VBC file */
   CHECK_OKAY( SCIPvbcNewChild(stat->vbc, stat, *node) );

   debugMessage("created child node %p at depth %d\n", *node, (*node)->depth);

   return SCIP_OKAY;
}

/** creates a probing child node of the active node */
static
RETCODE nodeCreateProbingNode(
   NODE**           node,               /**< pointer to node data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(node != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->pathlen > 0);
   assert(tree->actnode != NULL);
   assert(tree->actnode == tree->path[tree->pathlen-1]);
   assert(SCIPnodeGetType(tree->actnode) == SCIP_NODETYPE_ACTNODE);

   /* create the node data structure */
   CHECK_OKAY( nodeCreate(node, memhdr, set) );

   /* mark node to be a probing child node */
   (*node)->nodetype = SCIP_NODETYPE_PROBINGNODE; /*lint !e641*/

   /* make active node the parent of the new probing child node */
   CHECK_OKAY( nodeAssignParent(*node, memhdr, set, tree, tree->actnode, 0.0) );

   debugMessage("created probing child node %p at depth %d\n", *node, (*node)->depth);

   return SCIP_OKAY;
}

/** removes given node from the siblings array */
static
void treeRemoveSibling(
   TREE*            tree,               /**< branch and bound tree */
   NODE*            sibling             /**< sibling node to remove */
   )
{
   int delpos;

   assert(tree != NULL);
   assert(sibling != NULL);
   assert(SCIPnodeGetType(sibling) == SCIP_NODETYPE_SIBLING);
   assert(sibling->data.sibling.arraypos >= 0 && sibling->data.sibling.arraypos < tree->nsiblings);
   assert(tree->siblings[sibling->data.sibling.arraypos] == sibling);
   assert(SCIPnodeGetType(tree->siblings[tree->nsiblings-1]) == SCIP_NODETYPE_SIBLING);

   delpos = sibling->data.sibling.arraypos;

   /* move last sibling in array to position of removed sibling */
   tree->siblings[delpos] = tree->siblings[tree->nsiblings-1];
   tree->siblingsprio[delpos] = tree->siblingsprio[tree->nsiblings-1];
   tree->siblings[delpos]->data.sibling.arraypos = delpos;
   sibling->data.sibling.arraypos = -1;
   tree->nsiblings--;
}

/** removes given node from the children array */
static
void treeRemoveChild(
   TREE*            tree,               /**< branch and bound tree */
   NODE*            child               /**< child node to remove */
   )
{
   int delpos;

   assert(tree != NULL);
   assert(child != NULL);
   assert(SCIPnodeGetType(child) == SCIP_NODETYPE_CHILD);
   assert(child->data.child.arraypos >= 0 && child->data.child.arraypos < tree->nchildren);
   assert(tree->children[child->data.child.arraypos] == child);
   assert(SCIPnodeGetType(tree->children[tree->nchildren-1]) == SCIP_NODETYPE_CHILD);

   delpos = child->data.child.arraypos;

   /* move last child in array to position of removed child */
   tree->children[delpos] = tree->children[tree->nchildren-1];
   tree->childrenprio[delpos] = tree->childrenprio[tree->nchildren-1];
   tree->children[delpos]->data.child.arraypos = delpos;
   child->data.child.arraypos = -1;
   tree->nchildren--;
}

/** frees node */
RETCODE SCIPnodeFree(
   NODE**           node,               /**< node data */
   MEMHDR*          memhdr,             /**< block memory buffer */
   SET*             set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   )
{
   Bool isroot;

   assert(node != NULL);
   assert(*node != NULL);
   assert(!(*node)->active);
   assert(memhdr != NULL);
   assert(tree != NULL);

   debugMessage("free node %p at depth %d of type %d\n", *node, (*node)->depth, SCIPnodeGetType(*node));

   /* free nodetype specific data, and release no longer needed LPI states */
   switch( SCIPnodeGetType(*node) )
   {
   case SCIP_NODETYPE_ACTNODE:
      assert(tree->actnode == *node);
      assert(tree->probingnode == NULL);
      if( tree->actlpfork != NULL )
      {
         assert(SCIPnodeGetType(tree->actlpfork) == SCIP_NODETYPE_FORK
            || SCIPnodeGetType(tree->actlpfork) == SCIP_NODETYPE_SUBROOT);
         CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
      }
      tree->actnode = NULL;
      break;
   case SCIP_NODETYPE_PROBINGNODE:
      assert(tree->probingnode == *node);
      break;
   case SCIP_NODETYPE_SIBLING:
      assert((*node)->data.sibling.arraypos >= 0);
      assert((*node)->data.sibling.arraypos < tree->nsiblings);
      assert(tree->siblings[(*node)->data.sibling.arraypos] == *node);
      if( tree->actlpfork != NULL )
      {
         assert(SCIPnodeGetType(tree->actlpfork) == SCIP_NODETYPE_FORK
            || SCIPnodeGetType(tree->actlpfork) == SCIP_NODETYPE_SUBROOT);
         CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
      }
      treeRemoveSibling(tree, *node);
      break;
   case SCIP_NODETYPE_CHILD:
      assert((*node)->data.child.arraypos >= 0);
      assert((*node)->data.child.arraypos < tree->nchildren);
      assert(tree->children[(*node)->data.child.arraypos] == *node);
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
   case SCIP_NODETYPE_JUNCTION:
      break;
   case SCIP_NODETYPE_FORK:
      CHECK_OKAY( forkFree(&((*node)->data.fork), memhdr, set, lp) );
      break;
   case SCIP_NODETYPE_SUBROOT:
      CHECK_OKAY( subrootFree(&((*node)->data.subroot), memhdr, set, lp) );
      break;
   default:
      errorMessage("unknown node type\n");
      abort();
   }

   /* check, if the node to be freed is the root node */
   isroot = ((*node)->depth == 0);

   /* free common data */
   CHECK_OKAY( SCIPconssetchgFree(&(*node)->conssetchg, memhdr, set) );
   CHECK_OKAY( SCIPdomchgFree(&(*node)->domchg, memhdr, set) );
   CHECK_OKAY( nodeReleaseParent(*node, memhdr, set, tree, lp) );

   freeBlockMemory(memhdr, node);

   /* delete the tree's root node pointer, if the freed node was the root */
   if( isroot )
      tree->root = NULL;

   return SCIP_OKAY;
}

/** cuts off node and whole sub tree from branch and bound tree */
void SCIPnodeCutoff(
   NODE*            node,               /**< node that should be cut off */
   SET*             set,                /**< global SCIP settings */
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(node != NULL);
   assert(set != NULL);
   assert(tree != NULL);

   debugMessage("cutting off node %p in depth %d\n", node, node->depth);

   node->cutoff = TRUE;
   node->lowerbound = set->infinity;
   tree->cutoffdepth = MIN(tree->cutoffdepth, node->depth);
}

/** informs node, that it is now on the active path and applies any domain and constraint set changes */
static
RETCODE nodeActivate(
   NODE*            node,               /**< node to activate */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   assert(node != NULL);
   assert(!node->active);

   debugMessage("activate node %p at depth %d of type %d\n", node, node->depth, SCIPnodeGetType(node));

   /* apply domain and constraint set changes */
   CHECK_OKAY( SCIPconssetchgApply(node->conssetchg, memhdr, set, stat, node->depth) );
   CHECK_OKAY( SCIPdomchgApply(node->domchg, memhdr, set, stat, lp, branchcand, eventqueue, node->depth) );

   /* mark node active */
   node->active = TRUE;

   return SCIP_OKAY;
}

/** informs node, that it is no longer on the active path and undoes any domain and constraint set changes;
 *  frees node, if it has no children
 */
static
RETCODE nodeDeactivate(
   NODE**           node,               /**< pointer to node to deactivate */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   Bool hasChildren;

   assert(node != NULL);
   assert(*node != NULL);
   assert((*node)->active);
   assert(tree != NULL);

   debugMessage("deactivate node %p at depth %d of type %d\n", *node, (*node)->depth, SCIPnodeGetType(*node));

   /* undo domain and constraint set changes */
   CHECK_OKAY( SCIPdomchgUndo((*node)->domchg, memhdr, set, stat, lp, branchcand, eventqueue) );
   CHECK_OKAY( SCIPconssetchgUndo((*node)->conssetchg, memhdr, set, stat) );

   /* mark node inactive */
   (*node)->active = FALSE;

   /* check, if node has any children */
   switch( SCIPnodeGetType(*node) )
   {
   case SCIP_NODETYPE_ACTNODE:
      if( tree->probingnode != NULL )
      {
         errorMessage("Cannot deactivate active node while in probing mode\n");
         abort();
      }
      hasChildren = (tree->nchildren > 0);
      break;
   case SCIP_NODETYPE_PROBINGNODE:
      hasChildren = FALSE;
      break;
   case SCIP_NODETYPE_SIBLING:
      errorMessage("Cannot deactivate sibling (which shouldn't be active)\n");
      abort();
   case SCIP_NODETYPE_CHILD:
      errorMessage("Cannot deactivate child (which shouldn't be active)\n");
      abort();
   case SCIP_NODETYPE_LEAF:
      errorMessage("Cannot deactivate leaf (which shouldn't be active)\n");
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
      errorMessage("unknown node type\n");
      abort();
   }

   /* free node, if it has no children */
   if( !hasChildren )
   {
      CHECK_OKAY( SCIPnodeFree(node, memhdr, set, tree, lp) );
   }

   return SCIP_OKAY;
}

/** adds constraint locally to the node and captures it; activates constraint, if node is active;
 *  if a local constraint is added to the root node, it is automatically upgraded into a global constraint
 */
RETCODE SCIPnodeAddCons(
   NODE*            node,               /**< node to add constraint to */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(node != NULL);
   assert(cons != NULL);

   /* if node is the root, mark constraint to be globally valid */
   if( node->depth == 0 )
   {
      assert(node == tree->root);
      cons->local = FALSE;
   }

   /* add constraint addition to the node's constraint set change data, and activate constraint if node is active */
   CHECK_OKAY( SCIPconssetchgAddAddedCons(&node->conssetchg, memhdr, set, stat, cons, node->depth, node->active) );
   assert(node->conssetchg != NULL);
   assert(node->conssetchg->addedconss != NULL);
   assert(!node->active || SCIPconsIsActive(cons));

   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities at the node, and captures constraint;
 *  disables constraint, if node is active
 */
RETCODE SCIPnodeDisableCons(
   NODE*            node,               /**< node to add constraint to */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   CONS*            cons                /**< constraint to disable */
   )
{
   assert(node != NULL);
   assert(tree != NULL);
   assert(cons != NULL);

   debugMessage("disabling constraint <%s> at node in depth %d\n", cons->name, node->depth);

   /* add constraint disabling to the node's constraint set change data */
   CHECK_OKAY( SCIPconssetchgAddDisabledCons(&node->conssetchg, memhdr, set, cons) );
   assert(node->conssetchg != NULL);
   assert(node->conssetchg->disabledconss != NULL);

   /* disable constraint, if node is active */
   if( node->active && cons->enabled && !cons->updatedisable )
   {
      CHECK_OKAY( SCIPconsDisable(cons, set, stat) );
   }

   return SCIP_OKAY;
}

/** adds bound change with inference information to active node, child of active node, or probing node;
 *  if possible, adjusts bound to integral value
 */
RETCODE SCIPnodeAddBoundinfer(
   NODE*            node,               /**< node to add bound change to */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   CONS*            infercons,          /**< constraint that deduced the bound change, or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   Bool             probingchange       /**< is the bound change a temporary setting due to probing? */
   )
{
   VAR* infervar;
   BOUNDTYPE inferboundtype;
   Real oldbound;

   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_ACTNODE
      || node->nodetype == SCIP_NODETYPE_PROBINGNODE
      || node->nodetype == SCIP_NODETYPE_CHILD);
   assert(tree != NULL);
   assert(var != NULL);
   assert(node->active || infercons == NULL);
   assert(node->nodetype == SCIP_NODETYPE_PROBINGNODE || !probingchange);

   debugMessage("adding boundchange at node in depth %d to variable <%s>: old bounds=[%g,%g], new %s bound: %g\n",
      node->depth, SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), 
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", newbound);

   /* remember variable as inference variable, and get corresponding active variable, bound and bound type */
   infervar = var;
   inferboundtype = boundtype;
   CHECK_OKAY( SCIPvarGetProbvarBound(&var, &newbound, &boundtype) );

   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      oldbound = SCIPvarGetLbLocal(var);

      /* adjust the new bound */
      SCIPvarAdjustLb(var, set, &newbound);

      if( SCIPsetIsLE(set, newbound, oldbound) )
      {
         errorMessage("variable's lower bound was not tightened: var <%s>, oldbound=%f, newbound=%f\n",
            SCIPvarGetName(var), oldbound, newbound);
         return SCIP_INVALIDDATA;
      }
   }
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      oldbound = SCIPvarGetUbLocal(var);

      /* adjust the new bound */
      SCIPvarAdjustUb(var, set, &newbound);

      if( SCIPsetIsGE(set, newbound, oldbound) )
      {
         errorMessage("variable's upper bound was not tightened: var <%s>, oldbound=%f, newbound=%f\n",
            SCIPvarGetName(var), oldbound, newbound);
         return SCIP_INVALIDDATA;
      }
   }
   
   debugMessage(" -> transformed to active variable <%s>: old bounds=[%g,%g], new %s bound: %g, obj: %g\n",
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", newbound, SCIPvarGetObj(var));

   stat->nboundchgs++;

   /* if the node is the root node: change local and global bound immediately */
   if( node->depth == 0 )
   {
      assert(node->active);
      assert(node->nodetype != SCIP_NODETYPE_PROBINGNODE);

      debugMessage(" -> bound change in root node: perform global bound change\n");
      CHECK_OKAY( SCIPvarChgBdGlobal(var, set, newbound, boundtype) );
      CHECK_OKAY( SCIPvarChgBdLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, newbound, boundtype) );

      stat->nrootboundchgs++;
      stat->nrootboundchgsrun++;

      return SCIP_OKAY;
   }

   /* if the node is a child, or the bound is a temporary probing bound
    *  - the bound change is a branching decision
    *  - the child's lower bound can be updated due to the changed pseudo solution
    * otherwise:
    *  - the bound change is an inference
    */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD || probingchange )
   {
      Real newpseudoobjval;
      Real lpsolval;

      assert(!node->active || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);

      /* get the solution value of variable in last solved LP on the active path:
       *  - if the LP was solved at the current node, the LP values of the columns are valid
       *  - if the last solved LP was the one in the current lpfork, the LP value in the columns are still valid
       *  - otherwise, the LP values are invalid
       */
      if( tree->actnodehaslp
         || (tree->actlpforklpcount == stat->lpcount && SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN) )
      {
         lpsolval = SCIPvarGetLPSol(var);
      }
      else
         lpsolval = SCIP_INVALID;

      /* remember the bound change as branching decision (infervar/infercons are not important: use NULL) */
      CHECK_OKAY( SCIPdomchgAddBoundchg(&node->domchg, memhdr, set, stat,
            var, newbound, boundtype, SCIP_BOUNDCHGTYPE_BRANCHING, 
            lpsolval, NULL, NULL, 0, inferboundtype) );
      
      /* update the child's lower bound */
      if( set->exactsolve )
         newpseudoobjval = SCIPlpGetModifiedProvedPseudoObjval(lp, set, var, oldbound, newbound, boundtype);
      else
         newpseudoobjval = SCIPlpGetModifiedPseudoObjval(lp, set, var, oldbound, newbound, boundtype);
      SCIPnodeUpdateLowerbound(node, stat, newpseudoobjval);

      /* update the branching history */
      CHECK_OKAY( SCIPvarIncNBranchings(var, stat, node->depth, 
                     boundtype == SCIP_BOUNDTYPE_LOWER ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS) );
   }
   else
   {
      /* remember the bound change as inference (lpsolval is not important: use 0.0) */
      CHECK_OKAY( SCIPdomchgAddBoundchg(&node->domchg, memhdr, set, stat,
            var, newbound, boundtype, SCIP_BOUNDCHGTYPE_INFERENCE, 
            0.0, infervar, infercons, inferinfo, inferboundtype) );

      /* update the inference history */
      if( stat->lastbranchvar != NULL )
      {
         CHECK_OKAY( SCIPvarIncNInferences(stat->lastbranchvar, stat, stat->lastbranchdir) );
      }
      /**@todo if last branching variable is unknown, retrieve it from the nodes' boundchg arrays */
   }

   assert(node->domchg != NULL);
   assert(node->domchg->domchgdyn.domchgtype == SCIP_DOMCHGTYPE_DYNAMIC); /*lint !e641*/
   assert(node->domchg->domchgdyn.boundchgs != NULL);
   assert(node->domchg->domchgdyn.nboundchgs > 0);
   assert(node->domchg->domchgdyn.boundchgs[node->domchg->domchgdyn.nboundchgs-1].var == var);
   assert(node->domchg->domchgdyn.boundchgs[node->domchg->domchgdyn.nboundchgs-1].newbound == newbound); /*lint !e777*/
   
   /* if node is active, apply the bound change immediately */
   if( node->active )
   {
      CHECK_OKAY( SCIPboundchgApply(&node->domchg->domchgdyn.boundchgs[node->domchg->domchgdyn.nboundchgs-1],
                     memhdr, set, stat, lp, branchcand, eventqueue, node->depth, node->domchg->domchgdyn.nboundchgs-1) );
   }

   return SCIP_OKAY;
}

/** adds bound change to active node, or child of active node; if possible, adjusts bound to integral value */
RETCODE SCIPnodeAddBoundchg(
   NODE*            node,               /**< node to add bound change to */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   Bool             probingchange       /**< is the bound change a temporary setting due to probing? */
   )
{
   CHECK_OKAY( SCIPnodeAddBoundinfer(node, memhdr, set, stat, tree, lp, branchcand, eventqueue, var, newbound, boundtype,
         NULL, 0, probingchange) );

   return SCIP_OKAY;
}

/** adds hole change to active node, or child of active node */
RETCODE SCIPnodeAddHolechg(
   NODE*            node,               /**< node to add bound change to */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   HOLELIST**       ptr,                /**< changed list pointer */
   HOLELIST*        newlist,            /**< new value of list pointer */
   HOLELIST*        oldlist             /**< old value of list pointer */
   )
{
   assert(node != NULL);
   assert(node->nodetype == SCIP_NODETYPE_ACTNODE
      || node->nodetype == SCIP_NODETYPE_PROBINGNODE
      || node->nodetype == SCIP_NODETYPE_CHILD);

   debugMessage("adding holechange at node in depth %d: changed pointer at %p from %p to %p\n",
      node->depth, ptr, newlist, oldlist);

   stat->nholechgs++;

   CHECK_OKAY( SCIPdomchgAddHolechg(&node->domchg, memhdr, set, stat, ptr, newlist, oldlist) );

   /**@todo apply hole change on active nodes and issue event */

   return SCIP_OKAY;
}

/** if given value is larger than the node's lower bound, sets the node's lower bound to the new value */
void SCIPnodeUpdateLowerbound(
   NODE*            node,               /**< node to update lower bound for */
   STAT*            stat,               /**< problem statistics */
   Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   )
{
   assert(node != NULL);
   assert(stat != NULL);

   if( newbound > node->lowerbound )
   {
      node->lowerbound = newbound;
      if( node->depth == 0 )
         stat->rootlowerbound = newbound;
   }
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

   return (NODETYPE)(node->nodetype);
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
Real SCIPnodeGetLowerbound(
   NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->lowerbound;
}

/** returns whether node is in the path to the current node */
Bool SCIPnodeIsActive(
   NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->active;
}

#endif




/*
 * Path Switching
 */

/** updates the LP sizes of the active path starting at the given depth */
static
void treeUpdatePathLPSize(
   TREE*            tree,               /**< branch and bound tree */
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
      assert((int)(node->depth) == i);
      
      switch( SCIPnodeGetType(node) )
      {
      case SCIP_NODETYPE_ACTNODE:
         assert(i == tree->pathlen-1
            || (i == tree->pathlen-2 && SCIPnodeGetType(tree->path[tree->pathlen-1]) == SCIP_NODETYPE_PROBINGNODE));
         break;
      case SCIP_NODETYPE_PROBINGNODE:
         assert(i == tree->pathlen-1);
         break;
      case SCIP_NODETYPE_SIBLING:
         errorMessage("Sibling cannot be in the active path\n");
         abort();
      case SCIP_NODETYPE_CHILD:
         errorMessage("Child cannot be in the active path\n");
         abort();
      case SCIP_NODETYPE_LEAF:
         errorMessage("Leaf cannot be in the active path\n");
         abort();
      case SCIP_NODETYPE_DEADEND:
         errorMessage("Deadend cannot be in the active path\n");
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
         errorMessage("unknown node type\n");
         abort();
      }
      tree->pathnlpcols[i] = ncols;
      tree->pathnlprows[i] = nrows;
   }
}

/** finds the common fork node, the new LP defining fork, and the new active subroot, if the path is switched to
 *  the given node
 */
static
void treeFindSwitchForks(
   TREE*            tree,               /**< branch and bound tree */
   NODE*            node,               /**< new active node, or NULL */
   NODE**           commonfork,         /**< pointer to store common fork node of old and new active node */
   NODE**           newlpfork,          /**< pointer to store the new LP defining fork node */
   NODE**           newsubroot,         /**< pointer to store the new subroot node */
   Bool*            cutoff              /**< pointer to store whether the given node can be cut off and no path switching
                                         *   should be performed */
   )
{
   NODE* fork;
   NODE* lpfork;
   NODE* subroot;

   assert(tree != NULL);
   assert(tree->root != NULL);
   assert((tree->actnode == NULL) == !tree->root->active);
   assert(tree->actlpfork == NULL || tree->actnode != NULL);
   assert(tree->actlpfork == NULL || tree->actlpfork->depth < tree->actnode->depth);
   assert(tree->actsubroot == NULL || tree->actlpfork != NULL);
   assert(tree->actsubroot == NULL || tree->actsubroot->depth <= tree->actlpfork->depth);
   assert(tree->cutoffdepth >= 0);
   assert(tree->cutoffdepth == INT_MAX || tree->cutoffdepth < tree->pathlen);
   assert(tree->cutoffdepth == INT_MAX || tree->path[tree->cutoffdepth]->cutoff);
   assert(commonfork != NULL);
   assert(newlpfork != NULL);
   assert(newsubroot != NULL);
   assert(cutoff != NULL);

   *commonfork = NULL;
   *newlpfork = NULL;
   *newsubroot = NULL;
   *cutoff = FALSE;

   /* if the new active node is NULL, there is no common fork node, and the new LP fork and subroot are NULL */
   if( node == NULL )
      return;

   *cutoff = node->cutoff;

   /* if the old active node is NULL, there is no common fork node, and we have to search the new LP fork and subroot */
   if( tree->actnode == NULL )
   {
      assert(!tree->root->active);
      assert(tree->pathlen == 0);
      assert(tree->cutoffdepth == INT_MAX);

      lpfork = node;
      while( SCIPnodeGetType(lpfork) != SCIP_NODETYPE_FORK && SCIPnodeGetType(lpfork) != SCIP_NODETYPE_SUBROOT )
      {
         lpfork = lpfork->parent;
         if( lpfork == NULL )
            return;
         if( lpfork->cutoff )
         {
            *cutoff = TRUE;
            return;
         }
      }
      *newlpfork = lpfork;

      subroot = lpfork;
      while( SCIPnodeGetType(subroot) != SCIP_NODETYPE_SUBROOT )
      {
         subroot = subroot->parent;
         if( subroot == NULL )
            return;
         if( subroot->cutoff )
         {
            *cutoff = TRUE;
            return;
         }
      }
      *newsubroot = subroot;

      fork = subroot;
      while( fork->parent != NULL )
      {
         fork = fork->parent;
         if( fork->cutoff )
         {
            *cutoff = TRUE;
            return;
         }
      }
      return;
   }

   /* find the common fork node, the new LP defining fork, and the new active subroot */
   fork = node;
   lpfork = NULL;
   subroot = NULL;
   while( !fork->active )
   {
      fork = fork->parent;
      assert(fork != NULL); /* because the root is active, there must be a common fork node */

      if( fork->cutoff )
      {
         *cutoff = TRUE;
         return;
      }
      if( lpfork == NULL
         && (SCIPnodeGetType(fork) == SCIP_NODETYPE_FORK || SCIPnodeGetType(fork) == SCIP_NODETYPE_SUBROOT) )
         lpfork = fork;
      if( subroot == NULL && SCIPnodeGetType(fork) == SCIP_NODETYPE_SUBROOT )
         subroot = fork;
   }
   assert(fork != NULL);
   assert(lpfork == NULL || !lpfork->active || lpfork == fork);
   assert(subroot == NULL || !subroot->active || subroot == fork);
   debugMessage("find switch forks: forkdepth=%d\n", fork->depth);

   /* if the common fork node is below the current cutoff depth, the cutoff node is a parent of the common fork
    * and thus a parent of the new active node, s.t. the new node can also be cut off
    */
   assert(fork->depth != tree->cutoffdepth);
   if( fork->depth > tree->cutoffdepth )
   {
#ifndef NDEBUG
      while( fork != NULL && !fork->cutoff )
         fork = fork->parent;
      assert(fork != NULL);
      assert(fork->depth >= tree->cutoffdepth);
#endif
      *cutoff = TRUE;
      return;
   }
   tree->cutoffdepth = INT_MAX;

   /* if not already found, continue searching the LP defining fork; it can not be deeper than the common fork */
   if( lpfork == NULL )
   {
      if( tree->actlpfork != NULL && (int)(tree->actlpfork->depth) > fork->depth )
      {
         /* actlpfork is not on the same active path as the new node: we have to continue searching */
         lpfork = fork;
         while( lpfork != NULL
            && SCIPnodeGetType(lpfork) != SCIP_NODETYPE_FORK
            && SCIPnodeGetType(lpfork) != SCIP_NODETYPE_SUBROOT )
         {
            assert(lpfork->active);
            lpfork = lpfork->parent;
         }
      }
      else
      {
         /* actlpfork is on the same active path as the new node: old and new node have the same lpfork */
         lpfork = tree->actlpfork;
      }
      assert(lpfork == NULL || (int)(lpfork->depth) <= fork->depth);
      assert(lpfork == NULL || lpfork->active);
   }
   assert(lpfork == NULL
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_FORK
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_SUBROOT);
   debugMessage("find switch forks: lpforkdepth=%d\n", lpfork == NULL ? -1 : (int)(lpfork->depth));

   /* if not already found, continue searching the subroot; it cannot be deeper than the LP fork and common fork */
   if( subroot == NULL )
   {
      if( tree->actsubroot != NULL && (int)(tree->actsubroot->depth) > fork->depth )
      {
         /* actsubroot is not on the same active path as the new node: we have to continue searching */
         if( lpfork != NULL && lpfork->depth < fork->depth )
            subroot = lpfork;
         else
            subroot = fork;
         while( subroot != NULL && SCIPnodeGetType(subroot) != SCIP_NODETYPE_SUBROOT )
         {
            assert(subroot->active);
            subroot = subroot->parent;
         }
      }
      else
         subroot = tree->actsubroot;
      assert(subroot == NULL || subroot->depth <= fork->depth);
      assert(subroot == NULL || subroot->active);
   }
   assert(subroot == NULL || SCIPnodeGetType(subroot) == SCIP_NODETYPE_SUBROOT);
   assert(subroot == NULL || lpfork != NULL);
   assert(subroot == NULL || subroot->depth <= lpfork->depth);
   debugMessage("find switch forks: subrootdepth=%d\n", subroot == NULL ? -1 : (int)(subroot->depth));

   *commonfork = fork;
   *newlpfork = lpfork;
   *newsubroot = subroot;

#ifndef NDEBUG
   while( fork != NULL )
   {
      assert(fork->active);
      assert(!fork->cutoff);
      fork = fork->parent;
   }
#endif
}

/** switches the active path to the new active node, applies domain and constraint set changes */
static
RETCODE treeSwitchPath(
   TREE*            tree,               /**< branch and bound tree */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   NODE*            fork                /**< common fork node of old and new active node */
   )
{
   NODE* actnode;       /* new active node of the tree */
   int actnodedepth;    /* depth of the new active node, or -1 if actnode == NULL */
   int forkdepth;       /* depth of the common subroot/fork/junction node, or -1 if no common fork exists */
   int i;

   assert(tree != NULL);
   assert(fork == NULL || (fork->active && !fork->cutoff));

   actnode = tree->actnode;
   assert(actnode == NULL || (!actnode->active && !actnode->cutoff));
   assert(actnode == NULL || SCIPnodeGetType(actnode) == SCIP_NODETYPE_ACTNODE);
   
   debugMessage("switch path: old pathlen=%d\n", tree->pathlen);   

   /* get the nodes' depth's */
   actnodedepth = (actnode != NULL ? actnode->depth : -1);
   forkdepth = (fork != NULL ? fork->depth : -1);
   assert(forkdepth <= actnodedepth);
   assert(forkdepth < tree->pathlen);

   /* undo the domain and constraint set changes of the old active path by deactivating the path's nodes */
   for( i = tree->pathlen-1; i > forkdepth; --i )
   {
      CHECK_OKAY( nodeDeactivate(&(tree->path[i]), memhdr, set, stat, tree, lp, branchcand, eventqueue) );
   }

   /* create the new active path */
   CHECK_OKAY( treeEnsurePathMem(tree, set, actnodedepth+1) );
   tree->pathlen = actnodedepth+1;
   while( actnode != fork )
   {
      assert(actnode != NULL);
      assert(!actnode->active);
      assert(!actnode->cutoff);
      tree->path[actnode->depth] = actnode;
      actnode = actnode->parent;
   }

   /* apply domain and constraint set changes of the new path by activating the path's nodes */
   for( i = forkdepth+1; i < tree->pathlen; ++i )
   {
      assert(!tree->path[i]->cutoff);
      CHECK_OKAY( nodeActivate(tree->path[i], memhdr, set, stat, lp, branchcand, eventqueue) );
   }

   /* count the new LP sizes of the path */
   treeUpdatePathLPSize(tree, forkdepth+1);

   debugMessage("switch path: new pathlen=%d\n", tree->pathlen);   

   return SCIP_OKAY;
}

/** loads the subroot's LP data */
static
RETCODE subrootConstructLP(
   NODE*            subroot,            /**< subroot node to construct LP for */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   COL** cols;
   ROW** rows;
   int ncols;
   int nrows;
   int c;
   int r;

   assert(subroot != NULL);
   assert(SCIPnodeGetType(subroot) == SCIP_NODETYPE_SUBROOT);
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
   {
      CHECK_OKAY( SCIPlpAddCol(lp, set, cols[c], subroot->depth) );
   }
   for( r = 0; r < nrows; ++r )
   {
      CHECK_OKAY( SCIPlpAddRow(lp, set, rows[r], subroot->depth) );
   }

   return SCIP_OKAY;
}
   
/** loads the fork's additional LP data */
static
RETCODE forkAddLP(
   NODE*            fork,               /**< fork node to construct additional LP for */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   COL** cols;
   ROW** rows;
   int ncols;
   int nrows;
   int c;
   int r;

   assert(fork != NULL);
   assert(SCIPnodeGetType(fork) == SCIP_NODETYPE_FORK);
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
   {
      CHECK_OKAY( SCIPlpAddCol(lp, set, cols[c], fork->depth) );
   }
   for( r = 0; r < nrows; ++r )
   {
      CHECK_OKAY( SCIPlpAddRow(lp, set, rows[r], fork->depth) );
   }

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** checks validity of active path */
static
void treeCheckPath(
   TREE*            tree                /**< branch and bound tree */
   )
{
   NODE* node;
   int ncols;
   int nrows;
   int d;

   assert(tree != NULL);
   assert(tree->path != NULL);

   ncols = 0;
   nrows = 0;
   for( d = 0; d < tree->pathlen; ++d )
   {
      node = tree->path[d];
      assert(node != NULL);
      assert((int)(node->depth) == d);
      switch( SCIPnodeGetType(node) )
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
         assert(d == tree->pathlen-1
            || (d == tree->pathlen-2 && SCIPnodeGetType(tree->path[tree->pathlen-1]) == SCIP_NODETYPE_PROBINGNODE));
         break;
      case SCIP_NODETYPE_PROBINGNODE:
         assert(d == tree->pathlen-1);
         break;
      default:
         errorMessage("node in depth %d on active path has to be of type FORK, SUBROOT, ACTNODE, or PROBINGNODE, but is %d\n",
            d, SCIPnodeGetType(node));
         abort();
      }  /*lint !e788*/
      assert(tree->pathnlpcols[d] == ncols);
      assert(tree->pathnlprows[d] == nrows);
   }
}
#else
#define treeCheckPath(tree) /**/
#endif

/** constructs the LP and loads LP state for fork/subroot of the active node */
RETCODE SCIPtreeLoadLP(
   TREE*            tree,               /**< branch and bound tree */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   LP*              lp                  /**< current LP data */
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
   assert(SCIPnodeGetType(tree->actnode) == SCIP_NODETYPE_ACTNODE);
   assert(SCIPnodeGetDepth(tree->actnode) == tree->pathlen-1);
   assert(tree->actnode == tree->path[tree->pathlen-1]);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   debugMessage("load LP for current fork node %p at depth %d\n", 
      tree->actlpfork, tree->actlpfork == NULL ? -1 : (int)(tree->actlpfork->depth));
   debugMessage("-> old LP has %d cols and %d rows\n", SCIPlpGetNCols(lp), SCIPlpGetNRows(lp));
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
      assert(SCIPnodeGetType(lpfork) == SCIP_NODETYPE_FORK || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_SUBROOT);
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
      assert((int)(pathnode->depth) == d);
      assert(SCIPnodeGetType(pathnode) == SCIP_NODETYPE_JUNCTION || SCIPnodeGetType(pathnode) == SCIP_NODETYPE_FORK);
      if( SCIPnodeGetType(pathnode) == SCIP_NODETYPE_FORK )
      {
         CHECK_OKAY( forkAddLP(pathnode, memhdr, set, lp) );
      }
   }
   tree->correctlpdepth = MAX(tree->correctlpdepth, lpforkdepth);
   assert(lpforkdepth == -1 || tree->pathnlpcols[tree->correctlpdepth] == tree->pathnlpcols[lpforkdepth]);
   assert(lpforkdepth == -1 || tree->pathnlprows[tree->correctlpdepth] == tree->pathnlprows[lpforkdepth]);
   assert(lpforkdepth == -1 || SCIPlpGetNCols(lp) == tree->pathnlpcols[lpforkdepth]);
   assert(lpforkdepth == -1 || SCIPlpGetNRows(lp) == tree->pathnlprows[lpforkdepth]);
   assert(lpforkdepth >= 0 || SCIPlpGetNCols(lp) == 0);
   assert(lpforkdepth >= 0 || SCIPlpGetNRows(lp) == 0);

   /* load LP state, if existing */
   if( lpfork != NULL )
   {
      if( tree->actlpforklpcount != stat->lpcount )
      {
         if( SCIPnodeGetType(lpfork) == SCIP_NODETYPE_FORK )
         {
            assert(lpfork->data.fork != NULL);
            CHECK_OKAY( SCIPlpSetState(lp, memhdr, set, lpfork->data.fork->lpistate) );
         }
         else
         {
            assert(SCIPnodeGetType(lpfork) == SCIP_NODETYPE_SUBROOT);
            assert(lpfork->data.subroot != NULL);
            CHECK_OKAY( SCIPlpSetState(lp, memhdr, set, lpfork->data.subroot->lpistate) );
         }
         assert(lp->primalfeasible);
         assert(lp->dualfeasible);
      }
      else
      {
         lp->primalfeasible = TRUE;
         lp->dualfeasible = TRUE;
      }

      /* check the path from LP fork to active node for domain changes (destroying primal feasibility of LP basis) */
      for( d = lpforkdepth; d < (int)(tree->actnode->depth) && lp->primalfeasible; ++d )
      {
         assert(d < tree->pathlen);
         lp->primalfeasible = lp->primalfeasible
            && (tree->path[d]->domchg == NULL || tree->path[d]->domchg->domchgbound.nboundchgs == 0);
      }
   }

   /* mark the LP's size, such that we know which rows and columns were added in the new node */
   SCIPlpMarkSize(lp);

   debugMessage("-> new correctlpdepth: %d\n", tree->correctlpdepth);
   debugMessage("-> new LP has %d cols and %d rows, primalfeasible=%d, dualfeasible=%d\n", 
      SCIPlpGetNCols(lp), SCIPlpGetNRows(lp), lp->primalfeasible, lp->dualfeasible);

   return SCIP_OKAY;
}




/*
 * Node Conversion
 */

/** converts node into LEAF and moves it into the array of the node queue
 *  if node's lower bound is greater or equal than the given upper bound, the node is deleted;
 *  otherwise, it is moved to the node queue; anyways, the given pointer is NULL after the call
 */
static
RETCODE nodeToLeaf(
   NODE**           node,               /**< pointer to child or sibling node to convert */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   NODE*            lpfork,             /**< LP fork of the node */
   Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   )
{
   assert(SCIPnodeGetType(*node) == SCIP_NODETYPE_SIBLING || SCIPnodeGetType(*node) == SCIP_NODETYPE_CHILD);
   assert(lpfork == NULL || lpfork->depth < (*node)->depth);
   assert(lpfork == NULL || lpfork->active);
   assert(lpfork == NULL
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_FORK
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_SUBROOT);

   /* convert node into leaf */
   debugMessage("convert node %p at depth %d to leaf with lpfork %p at depth %d\n",
      *node, (*node)->depth, lpfork, lpfork == NULL ? -1 : (int)(lpfork->depth));
   (*node)->nodetype = SCIP_NODETYPE_LEAF; /*lint !e641*/
   (*node)->data.leaf.lpfork = lpfork;
   assert(lpfork == NULL || lpfork->depth == (*node)->depth - 1); /*???????????????????????? WRONG !!*/

   /* if node is good enough to keep, put it on the node queue */
   if( SCIPsetIsLT(set, (*node)->lowerbound, cutoffbound) )
   {
      /* insert leaf in node queue */
      CHECK_OKAY( SCIPnodepqInsert(tree->leaves, set, *node) );
      
      /* make the domain change data static to save memory */
      CHECK_OKAY( SCIPdomchgMakeStatic(&(*node)->domchg, memhdr, set) );

      /* node is now member of the node queue: delete the pointer to forbid further access */
      *node = NULL;
   }
   else
   {
      /* delete node due to bound cut off */
      CHECK_OKAY( SCIPnodeFree(node, memhdr, set, tree, lp) );
   }
   assert(*node == NULL);

   return SCIP_OKAY;
}

/** converts the active node into a deadend node */
static
RETCODE actnodeToDeadend(
   MEMHDR*          memhdr,             /**< block memory buffers */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   )
{
   assert(memhdr != NULL);
   assert(tree != NULL);
   assert(tree->probingnode == NULL);
   assert(tree->actnode != NULL);
   assert(SCIPnodeGetType(tree->actnode) == SCIP_NODETYPE_ACTNODE);
   assert(tree->actnode->active);
   assert(tree->nchildren == 0);

   debugMessage("actnode to deadend %p at depth %d\n", tree->actnode, tree->actnode->depth);

   tree->actnode->nodetype = SCIP_NODETYPE_DEADEND; /*lint !e641*/

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
   SET*             set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);
   assert(tree->probingnode == NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnode->active);
   assert(SCIPnodeGetType(tree->actnode) == SCIP_NODETYPE_ACTNODE);

   debugMessage("actnode %p to junction at depth %d\n", tree->actnode, tree->actnode->depth);

   tree->actnode->nodetype = SCIP_NODETYPE_JUNCTION; /*lint !e641*/

   CHECK_OKAY( junctionInit(&tree->actnode->data.junction, tree) );

   /* release LPI state */
   if( tree->actlpfork != NULL )
   {
      CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
   }

   /* make the domain change data static to save memory */
   CHECK_OKAY( SCIPdomchgMakeStatic(&tree->actnode->domchg, memhdr, set) );

   return SCIP_OKAY;
}

/** converts the active node into a fork node */
static
RETCODE actnodeToFork(
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   )
{
   FORK* fork;
   Bool lperror;

   assert(memhdr != NULL);
   assert(tree != NULL);
   assert(tree->probingnode == NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnode->active);
   assert(SCIPnodeGetType(tree->actnode) == SCIP_NODETYPE_ACTNODE);
   assert(tree->nchildren > 0);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   debugMessage("actnode %p to fork at depth %d\n", tree->actnode, tree->actnode->depth);

   /* usually, the LP should be solved to optimality; otherwise, numerical troubles occured,
    * and we have to forget about the LP and transform the node into a junction (see below)
    */
   lperror = FALSE;
   if( lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      /* clean up newly created part of LP to keep only necessary columns and rows */
      CHECK_OKAY( SCIPlpCleanupNew(lp, memhdr, set, stat) );

      /* resolve LP after cleaning up */
      if( !lp->solved || !lp->flushed )
      {
         debugMessage("resolving LP after cleanup\n");
         CHECK_OKAY( SCIPlpSolve(lp, memhdr, set, stat, set->fastmip, FALSE, &lperror) );

         /* if the solution was valid before resolve, it is still valid (resolve didn't change the solution) */
         if( lp->validsollp == stat->lpcount-1 )
            lp->validsollp = stat->lpcount;
         if( lp->validfarkaslp == stat->lpcount-1 )
            lp->validfarkaslp = stat->lpcount;
      }
   }
   assert(lp->flushed);
   assert(lp->solved || lperror);

   /* There are two reasons, that the (reduced) LP is not solved to optimality:
    *  - The primal heuristics (called after the current node's LP was solved) found a new 
    *    solution, that is better than the current node's lower bound.
    *    (But in this case, all children should be cut off and the node should be converted
    *    into a deadend instead of a fork.)
    *  - Something numerically weird happened after cleaning up.
    * The only thing we can do, is to completely forget about the LP and treat the node as
    * if it was only a pseudo-solution node. Therefore we have to remove all additional
    * columns and rows from the LP and convert the node into a junction.
    * However, the node's lower bound is kept, thus automatically throwing away nodes that
    * were cut off due to a primal solution.
    */
   if( lperror || lp->lpsolstat != SCIP_LPSOLSTAT_OPTIMAL )
   {
      infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles: LP %d not optimal -- convert node into junction instead of fork\n", 
         stat->nnodes, stat->nlps);

      /* remove all additions to the LP at this node */
      CHECK_OKAY( SCIPlpShrinkCols(lp, SCIPlpGetNCols(lp) - SCIPlpGetNNewcols(lp)) );
      CHECK_OKAY( SCIPlpShrinkRows(lp, memhdr, set, SCIPlpGetNRows(lp) - SCIPlpGetNNewrows(lp)) );
      
      /* convert node into a junction */
      CHECK_OKAY( actnodeToJunction(memhdr, set, tree, lp) );
      
      return SCIP_OKAY;
   }
   assert(lp->solved);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL);

   /* create fork data */
   CHECK_OKAY( forkCreate(&fork, memhdr, tree, lp) );
   
   tree->actnode->nodetype = SCIP_NODETYPE_FORK; /*lint !e641*/
   tree->actnode->data.fork = fork;

   /* release LPI state */
   if( tree->actlpfork != NULL )
   {
      CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
   }

   /* make the domain change data static to save memory */
   CHECK_OKAY( SCIPdomchgMakeStatic(&tree->actnode->domchg, memhdr, set) );

   return SCIP_OKAY;
}

/** converts the active node into a subroot node */
static
RETCODE actnodeToSubroot(
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   )
{
   SUBROOT* subroot;
   Bool lperror;

   assert(memhdr != NULL);
   assert(tree != NULL);
   assert(tree->probingnode == NULL);
   assert(tree->actnode != NULL);
   assert(SCIPnodeGetType(tree->actnode) == SCIP_NODETYPE_ACTNODE);
   assert(tree->actnode->active);
   assert(tree->nchildren > 0);
   assert(lp != NULL);
   assert(lp->solved);

   debugMessage("actnode %p to subroot at depth %d\n", tree->actnode, tree->actnode->depth);

   /* usually, the LP should be solved to optimality; otherwise, numerical troubles occured,
    * and we have to forget about the LP and transform the node into a junction (see below)
    */
   lperror = FALSE;
   if( lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      /* clean up whole LP to keep only necessary columns and rows */
#if 0
      if( tree->actnode->depth == 0 )
      {
         CHECK_OKAY( SCIPlpCleanupAll(lp, memhdr, set) );
      }
      else
#endif
      {
         CHECK_OKAY( SCIPlpRemoveAllObsoletes(lp, memhdr, set, stat) );
      }

      /* resolve LP after cleaning up */
      if( !lp->solved || !lp->flushed )
      {
         CHECK_OKAY( SCIPlpSolve(lp, memhdr, set, stat, set->fastmip, FALSE, &lperror) );

         /* if the solution was valid before resolve, it is still valid (resolve didn't change the solution) */
         if( lp->validsollp == stat->lpcount-1 )
            lp->validsollp = stat->lpcount;
         if( lp->validfarkaslp == stat->lpcount-1 )
            lp->validfarkaslp = stat->lpcount;
      }
   }
   assert(lp->flushed);
   assert(lp->solved || lperror);

   /* There are two reasons, that the (reduced) LP is not solved to optimality:
    *  - The primal heuristics (called after the current node's LP was solved) found a new 
    *    solution, that is better than the current node's lower bound.
    *    (But in this case, all children should be cut off and the node should be converted
    *    into a deadend instead of a subroot.)
    *  - Something numerically weird happened after cleaning up.
    * The only thing we can do, is to completely forget about the LP and treat the node as
    * if it was only a pseudo-solution node. Therefore we have to remove all additional
    * columns and rows from the LP and convert the node into a junction.
    * However, the node's lower bound is kept, thus automatically throwing away nodes that
    * were cut off due to a primal solution.
    */
   if( lperror || lp->lpsolstat != SCIP_LPSOLSTAT_OPTIMAL )
   {
      infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles: LP %d not optimal -- convert node into junction instead of subroot\n", 
         stat->nnodes, stat->nlps);

      /* remove all additions to the LP at this node */
      CHECK_OKAY( SCIPlpShrinkCols(lp, SCIPlpGetNCols(lp) - SCIPlpGetNNewcols(lp)) );
      CHECK_OKAY( SCIPlpShrinkRows(lp, memhdr, set, SCIPlpGetNRows(lp) - SCIPlpGetNNewrows(lp)) );
      
      /* convert node into a junction */
      CHECK_OKAY( actnodeToJunction(memhdr, set, tree, lp) );
      
      return SCIP_OKAY;
   }
   assert(lp->solved);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL);

   /* create subroot data */
   CHECK_OKAY( subrootCreate(&subroot, memhdr, tree, lp) );

   tree->actnode->nodetype = SCIP_NODETYPE_SUBROOT; /*lint !e641*/
   tree->actnode->data.subroot = subroot;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, tree->actnode->depth);

   /* release LPI state */
   if( tree->actlpfork != NULL )
   {
      CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
   }

   /* make the domain change data static to save memory */
   CHECK_OKAY( SCIPdomchgMakeStatic(&tree->actnode->domchg, memhdr, set) );

   return SCIP_OKAY;
}

/** puts all nodes in the array on the node queue and makes them LEAFs */
static
RETCODE treeNodesToQueue(
   TREE*            tree,               /**< branch and bound tree */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   NODE**           nodes,              /**< array of nodes to put on the queue */
   int*             nnodes,             /**< pointer to number of nodes in the array */
   NODE*            lpfork,             /**< LP fork of the nodes */
   Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   )
{
   int i;

   assert(tree != NULL);
   assert(set != NULL);
   assert(nnodes != NULL);
   assert(*nnodes == 0 || nodes != NULL);

   for( i = 0; i < *nnodes; ++i )
   {
      /* convert node to LEAF and put it into leaves queue, or delete it if it's lower bound exceeds the cutoff bound */
      CHECK_OKAY( nodeToLeaf(&nodes[i], memhdr, set, tree, lp, lpfork, cutoffbound) );
      assert(nodes[i] == NULL);
   }
   *nnodes = 0;

   return SCIP_OKAY;
}

/** converts children into siblings, clears children array */
static
void treeChildrenToSiblings(
   TREE*            tree                /**< branch and bound tree */
   )
{
   NODE** tmpnodes;
   Real* tmpprios;
   int tmpnodessize;
   int i;

   assert(tree != NULL);
   assert(tree->nsiblings == 0);

   tmpnodes = tree->siblings;
   tmpprios = tree->siblingsprio;
   tmpnodessize = tree->siblingssize;

   tree->siblings = tree->children;
   tree->siblingsprio = tree->childrenprio;
   tree->nsiblings = tree->nchildren;
   tree->siblingssize = tree->childrensize;

   tree->children = tmpnodes;
   tree->childrenprio = tmpprios;
   tree->nchildren = 0;
   tree->childrensize = tmpnodessize;
   
   for( i = 0; i < tree->nsiblings; ++i )
   {
      assert(SCIPnodeGetType(tree->siblings[i]) == SCIP_NODETYPE_CHILD);
      tree->siblings[i]->nodetype = SCIP_NODETYPE_SIBLING; /*lint !e641*/

      /* because CHILD.arraypos and SIBLING.arraypos are on the same position, we do not have to copy it */
      assert(&(tree->siblings[i]->data.sibling.arraypos) == &(tree->siblings[i]->data.child.arraypos));
   }
}

/** activates a child, a sibling, or a leaf node */
RETCODE SCIPnodeActivate(
   NODE**           node,               /**< pointer to node to activate (or NULL to deactivate all nodes); the node
                                         *   is freed, if it was cut off due to a cut off subtree */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             cutoffbound,        /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   Bool*            cutoff              /**< pointer to store whether the given node can be cut off and no path switching
                                         *   should be performed */
   )
{
   NODE* fork;
   NODE* lpfork;
   NODE* subroot;
   NODE* childrenlpfork;

   assert(node != NULL);
   assert(*node == NULL
      || SCIPnodeGetType(*node) == SCIP_NODETYPE_SIBLING
      || SCIPnodeGetType(*node) == SCIP_NODETYPE_CHILD
      || SCIPnodeGetType(*node) == SCIP_NODETYPE_LEAF);
   assert(*node == NULL || !(*node)->active);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(tree->probingnode == NULL);
   assert(cutoff != NULL);

   /* find the common fork node, the new LP defining fork, and the new active subroot,
    * thereby checking, if the new node can be cut off
    */
   treeFindSwitchForks(tree, *node, &fork, &lpfork, &subroot, cutoff);
   debugMessage("activate node: actnodedepth=%d, forkdepth=%d, lpforkdepth=%d, subrootdepth=%d, cutoff=%d\n",
      *node != NULL ? (*node)->depth : -1, fork != NULL ? fork->depth : -1, lpfork != NULL ? lpfork->depth : -1,
      subroot != NULL ? subroot->depth : -1, *cutoff);

   /* free the new node, if it is located in a cut off subtree */
   if( *cutoff )
   {
      assert(*node != NULL);
      if( SCIPnodeGetType(*node) == SCIP_NODETYPE_LEAF )
      {
         CHECK_OKAY( SCIPnodepqRemove(tree->leaves, set, *node) );
      }
      CHECK_OKAY( SCIPnodeFree(node, memhdr, set, tree, lp) );
      stat->ndelayedcutoffs++;

      return SCIP_OKAY;
   }

   assert(fork == NULL || fork->active);
   assert(lpfork == NULL || fork != NULL);
   assert(subroot == NULL || lpfork != NULL);

   /* remember the depth of the common fork node for LP updates */
   debugMessage("activate node: old correctlpdepth=%d\n", tree->correctlpdepth);
   if( subroot == tree->actsubroot && fork != NULL )
   {
      /* we are in the same subtree: the LP is correct at most upto the common fork depth */
      assert(subroot == NULL || subroot->active);
      tree->correctlpdepth = MIN(tree->correctlpdepth, fork->depth);
   }
   else
   {
      /* we are in a different subtree: the LP is completely incorrect */
      assert(subroot == NULL || !subroot->active
         || (tree->actsubroot != NULL && (int)(tree->actsubroot->depth) > subroot->depth));
      tree->correctlpdepth = -1;
   }

   /* if the LP fork changed, the lpcount information for the new LP fork is unknown */
   if( lpfork != tree->actlpfork )
      tree->actlpforklpcount = -1;

   /* convert the old active node into a fork or subroot node, if it has children;
    * otherwise, convert it into a deadend, which will be freed later in treeSwitchPath()
    */
   childrenlpfork = NULL;
   if( tree->nchildren > 0 )
   {
      assert(tree->actnode != NULL);
      assert(SCIPnodeGetType(tree->actnode) == SCIP_NODETYPE_ACTNODE);

      if( tree->actnodehaslp )
      {
         /**@todo decide: old active node becomes fork or subroot */
         if( FALSE && tree->actnode->depth > 0 && tree->actnode->depth % 25 == 0 ) /* ????????????? */
         {
            /* convert old active node into a subroot node */
            CHECK_OKAY( actnodeToSubroot(memhdr, set, stat, tree, lp) );
            if( *node != NULL && SCIPnodeGetType(*node) == SCIP_NODETYPE_CHILD )
               subroot = tree->actnode;
         }
         else
         {
            /* convert old active node into a fork node */
            CHECK_OKAY( actnodeToFork(memhdr, set, stat, tree, lp) );
         }
         childrenlpfork = tree->actnode;

         /* update the path's LP size */
         tree->pathnlpcols[tree->actnode->depth] = lp->ncols;
         tree->pathnlprows[tree->actnode->depth] = lp->nrows;

         /* if a child of the old active node was selected as new active node, the old node becomes the new active
          * LP fork
          */
         if( *node != NULL && SCIPnodeGetType(*node) == SCIP_NODETYPE_CHILD )
         {
            lpfork = tree->actnode;
            tree->correctlpdepth = tree->actnode->depth;
            tree->actlpforklpcount = stat->lpcount;
         }
      }
      else
      {
         /* convert old active node into junction */
         CHECK_OKAY( actnodeToJunction(memhdr, set, tree, lp) );
         childrenlpfork = lpfork;
      }
   }
   else if( tree->actnode != NULL )
   {
      /* convert old active node into deadend */
      CHECK_OKAY( actnodeToDeadend(memhdr, tree, lp) );
   }
   debugMessage("activate node: new correctlpdepth=%d\n", tree->correctlpdepth);

   /* set up the new lists of siblings and children */
   if( *node == NULL )
   {
      /* move siblings to the queue, make them LEAFs */
      CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, lp, tree->siblings, &tree->nsiblings, tree->actlpfork,
            cutoffbound) );

      /* move children to the queue, make them LEAFs */
      CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, lp, tree->children, &tree->nchildren, childrenlpfork, cutoffbound) );
   }
   else
   {
      switch( SCIPnodeGetType(*node) )
      {  
      case SCIP_NODETYPE_SIBLING:
         /* move children to the queue, make them LEAFs */
         CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, lp, tree->children, &tree->nchildren, childrenlpfork,
               cutoffbound) );

         /* remove selected sibling from the siblings array */
         treeRemoveSibling(tree, *node);

         debugMessage("selected sibling node, lowerbound=%g, plungedepth=%d\n", (*node)->lowerbound, stat->plungedepth);
         break;
         
      case SCIP_NODETYPE_CHILD:
         /* move siblings to the queue, make them LEAFs */
         CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, lp, tree->siblings, &tree->nsiblings, tree->actlpfork,
               cutoffbound) );

         /* remove selected child from the children array */      
         treeRemoveChild(tree, *node);
         
         /* move other children to the siblings array, make them SIBLINGs */
         treeChildrenToSiblings(tree);
         
         stat->plungedepth++;
         debugMessage("selected child node, lowerbound=%g, plungedepth=%d\n", (*node)->lowerbound, stat->plungedepth);
         break;
         
      case SCIP_NODETYPE_LEAF:
         /* move siblings to the queue, make them LEAFs */
         CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, lp, tree->siblings, &tree->nsiblings, tree->actlpfork,
               cutoffbound) );
         
         /* move children to the queue, make them LEAFs */
         CHECK_OKAY( treeNodesToQueue(tree, memhdr, set, lp, tree->children, &tree->nchildren, childrenlpfork,
               cutoffbound) );

         /* remove node from the queue */
         CHECK_OKAY( SCIPnodepqRemove(tree->leaves, set, *node) );

         stat->plungedepth = 0;
         stat->nbacktracks++;
         debugMessage("selected leaf node, lowerbound=%g, plungedepth=%d\n", (*node)->lowerbound, stat->plungedepth);
         break;

      default:
         errorMessage("Selected node is neither sibling, child, nor leaf (nodetype=%d)\n", SCIPnodeGetType(*node));
         return SCIP_INVALIDDATA;
      }  /*lint !e788*/

      /* convert node into the active node */
      (*node)->nodetype = SCIP_NODETYPE_ACTNODE; /*lint !e641*/
   }
   assert(tree->nchildren == 0);
   
   /* set new active node, current LP fork, and subroot */
   assert(subroot == NULL || (lpfork != NULL && subroot->depth <= lpfork->depth));
   assert(lpfork == NULL || (*node != NULL && lpfork->depth < (*node)->depth));
   tree->actnode = *node;
   tree->actlpfork = lpfork;
   tree->actsubroot = subroot;

   /* track the path from the old active node to the new node, and perform domain changes */
   CHECK_OKAY( treeSwitchPath(tree, memhdr, set, stat, lp, branchcand, eventqueue, fork) );
   assert(tree->pathlen >= 0);
   assert((*node == NULL) == (tree->pathlen == 0));
   assert(*node == NULL || tree->path[tree->pathlen-1] == *node);

   return SCIP_OKAY;
}   




/*
 * Tree methods
 */

/** creates an initialized tree data structure */
RETCODE SCIPtreeCreate(
   TREE**           tree,               /**< pointer to tree data structure */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   NODESEL*         nodesel             /**< node selector to use for sorting leaves in the priority queue */
   )
{
   assert(tree != NULL);
   assert(set != NULL);

   ALLOC_OKAY( allocMemory(tree) );

   (*tree)->root = NULL;

   CHECK_OKAY( SCIPnodepqCreate(&(*tree)->leaves, set, nodesel) );

   (*tree)->path = NULL;
   (*tree)->actnode = NULL;
   (*tree)->actlpfork = NULL;
   (*tree)->actsubroot = NULL;
   (*tree)->children = NULL;
   (*tree)->siblings = NULL;
   (*tree)->probingnode = NULL;
   (*tree)->childrenprio = NULL;
   (*tree)->siblingsprio = NULL;
   (*tree)->pathnlpcols = NULL;
   (*tree)->pathnlprows = NULL;
   (*tree)->actlpforklpcount = -1;
   (*tree)->childrensize = 0;
   (*tree)->nchildren = 0;
   (*tree)->siblingssize = 0;
   (*tree)->nsiblings = 0;
   (*tree)->pathlen = 0;
   (*tree)->pathsize = 0;
   (*tree)->correctlpdepth = -1;
   (*tree)->cutoffdepth = INT_MAX;
   (*tree)->actnodehaslp = FALSE;
   (*tree)->cutoffdelayed = FALSE;

   /* create root node */
   CHECK_OKAY( SCIPnodeCreateChild(&(*tree)->root, memhdr, set, stat, *tree, 0.0) );
   assert((*tree)->nchildren == 1);

   /* move root to the queue, convert it to LEAF */
   CHECK_OKAY( treeNodesToQueue(*tree, memhdr, set, lp, (*tree)->children, &(*tree)->nchildren, NULL, set->infinity) );

   return SCIP_OKAY;
}

/** frees tree data structure */
RETCODE SCIPtreeFree(
   TREE**           tree,               /**< pointer to tree data structure */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);
   assert(*tree != NULL);
   assert((*tree)->nchildren == 0);
   assert((*tree)->nsiblings == 0);
   assert((*tree)->actnode == NULL);
   assert((*tree)->probingnode == NULL);

   debugMessage("free tree\n");

   /* free node queue */
   CHECK_OKAY( SCIPnodepqFree(&(*tree)->leaves, memhdr, set, *tree, lp) );
   
   /* free pointer arrays */
   freeMemoryArrayNull(&(*tree)->path);
   freeMemoryArrayNull(&(*tree)->children);
   freeMemoryArrayNull(&(*tree)->siblings);
   freeMemoryArrayNull(&(*tree)->childrenprio);
   freeMemoryArrayNull(&(*tree)->siblingsprio);
   freeMemoryArrayNull(&(*tree)->pathnlpcols);
   freeMemoryArrayNull(&(*tree)->pathnlprows);

   freeMemory(tree);

   return SCIP_OKAY;
}

/** returns the node selector associated with the given node priority queue */
NODESEL* SCIPtreeGetNodesel(
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqGetNodesel(tree->leaves);
}

/** sets the node selector used for sorting the nodes in the priority queue, and resorts the queue if necessary */
RETCODE SCIPtreeSetNodesel(
   TREE*            tree,               /**< branch and bound tree */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   NODESEL*         nodesel             /**< node selector to use for sorting the nodes in the queue */
   )
{
   assert(tree != NULL);
   assert(stat != NULL);

   if( SCIPnodepqGetNodesel(tree->leaves) != nodesel )
   {
      /* change the node selector used in the priority queue and resort the queue */
      CHECK_OKAY( SCIPnodepqSetNodesel(&tree->leaves, set, nodesel) );

      /* issue message */
      if( stat->nnodes > 0 )
      {
         infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
            "(node %lld) switching to node selector <%s>\n", stat->nnodes, SCIPnodeselGetName(nodesel));
      }
   }

   return SCIP_OKAY;
}

/** cuts off nodes with lower bound not better than given cutoff bound */
RETCODE SCIPtreeCutoff(
   TREE*            tree,               /**< branch and bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   )
{
   NODE* node;
   int i;

   assert(tree != NULL);
   assert(lp != NULL);

   /* if we are in diving mode, it is not allowed to cut off nodes, because this can lead to deleting LP rows which
    * would modify the currently unavailable (due to diving modifications) LP
    *  -> the cutoff must be delayed and executed after the diving ends
    */
   if( SCIPlpDiving(lp) )
   {
      tree->cutoffdelayed = TRUE;
      return SCIP_OKAY;
   }

   tree->cutoffdelayed = FALSE;

   /* cut off leaf nodes in the queue */
   CHECK_OKAY( SCIPnodepqBound(tree->leaves, memhdr, set, tree, lp, cutoffbound) );

   /* cut off siblings: we have to loop backwards, because a removal leads to moving the last node in empty slot */
   for( i = tree->nsiblings-1; i >= 0; --i )
   {
      node = tree->siblings[i];
      if( SCIPsetIsGE(set, node->lowerbound, cutoffbound) )
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
      if( SCIPsetIsGE(set, node->lowerbound, cutoffbound) )
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
   TREE*            tree,               /**< branch and bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             var                 /**< variable to branch on */
   )
{
   NODE* node;
   Real solval;
   Real rootsolval;

   assert(tree != NULL);
   assert(set != NULL);
   assert(var != NULL);

   /* get the corresponding active problem variable */
   var = SCIPvarGetProbvar(var);

   if( var == NULL )
   {
      errorMessage("cannot branch on a fixed variable\n");
      return SCIP_INVALIDDATA;
   }

   assert(SCIPvarGetProbindex(var) >= 0);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY
      || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER
      || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);
   assert(SCIPsetIsIntegral(set, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsIntegral(set, SCIPvarGetUbLocal(var)));
   assert(SCIPsetIsLT(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

   solval = SCIPvarGetSol(var, tree->actnodehaslp);
   rootsolval = SCIPvarGetRootSol(var);
   assert(SCIPsetIsGE(set, solval, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsLE(set, solval, SCIPvarGetUbLocal(var)));
   
   if( SCIPsetIsIntegral(set, solval) )
   {
      Real fixval;

      /* create child nodes with x <= x'-1, x = x', and x >= x'+1;
       * set the node selection priority in a way, s.t. a node is preferred whose branching goes in the same direction
       * as the deviation from the variable's root solution; evaluate x = x' first in any way
       */
      fixval = SCIPsetCeil(set, solval);
      assert(SCIPsetIsEQ(set, SCIPsetCeil(set, solval), SCIPsetFloor(set, solval)));
      
      debugMessage("pseudo branch on variable <%s> with value %g, priority %d\n", 
         SCIPvarGetName(var), solval, SCIPvarGetBranchPriority(var));
      
      /* create child node with x <= x'-1, if this would be feasible */
      if( SCIPsetIsGE(set, fixval-1, SCIPvarGetLbLocal(var)) )
      {
         debugMessage(" -> creating child: <%s> <= %g\n", SCIPvarGetName(var), fixval-1);
         CHECK_OKAY( SCIPnodeCreateChild(&node, memhdr, set, stat, tree, rootsolval - solval) );
         CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
               var, fixval-1, SCIP_BOUNDTYPE_UPPER, FALSE) );
         debugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
      }
                  
      /* create child node with x = x' */
      debugMessage(" -> creating child: <%s> == %g\n", SCIPvarGetName(var), fixval);
      CHECK_OKAY( SCIPnodeCreateChild(&node, memhdr, set, stat, tree, set->infinity) );
      if( !SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), fixval) )
      {
         CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
               var, fixval, SCIP_BOUNDTYPE_LOWER, FALSE) );
      }
      if( !SCIPsetIsEQ(set, SCIPvarGetUbLocal(var), fixval) )
      {
         CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
               var, fixval, SCIP_BOUNDTYPE_UPPER, FALSE) );
      }
      debugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
      
      /* create child node with x >= x'+1, if this would be feasible */
      if( SCIPsetIsLE(set, fixval+1, SCIPvarGetUbLocal(var)) )
      {
         debugMessage(" -> creating child: <%s> >= %g\n", SCIPvarGetName(var), fixval+1);
         CHECK_OKAY( SCIPnodeCreateChild(&node, memhdr, set, stat, tree, solval - rootsolval) );
         CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
               var, fixval+1, SCIP_BOUNDTYPE_LOWER, FALSE) );
         debugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
      }
   }
   else
   {   
      /* create child nodes with x <= floor(x'), and x >= ceil(x');
       * set the node selection priority in a way, s.t. a node is preferred whose branching goes in the same direction
       * as the deviation from the variable's root solution
       */
      debugMessage("LP branch on variable <%s> with value %g, priority %d\n", 
         SCIPvarGetName(var), solval, SCIPvarGetBranchPriority(var));
      
      /* create child node with x <= floor(x') */
      debugMessage(" -> creating child: <%s> <= %g\n", SCIPvarGetName(var), SCIPsetFloor(set, solval));
      CHECK_OKAY( SCIPnodeCreateChild(&node, memhdr, set, stat, tree, rootsolval - solval) );
      CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
            var, SCIPsetFloor(set, solval), SCIP_BOUNDTYPE_UPPER, FALSE) );
      debugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
      
      /* create child node with x >= ceil(x') */
      debugMessage(" -> creating child: <%s> >= %g\n", SCIPvarGetName(var), SCIPsetCeil(set, solval));
      CHECK_OKAY( SCIPnodeCreateChild(&node, memhdr, set, stat, tree, solval - rootsolval) );
      CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                     var, SCIPsetCeil(set, solval), SCIP_BOUNDTYPE_LOWER, FALSE) );
      debugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
   }

   return SCIP_OKAY;
}

/** switches to probing mode */
RETCODE SCIPtreeStartProbing(
   TREE*            tree,               /**< branch and bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(tree != NULL);
   assert(tree->probingnode == NULL);
   assert(tree->actnode != NULL);
   assert(tree->pathlen >= 1);
   assert(tree->path[tree->pathlen-1] == tree->actnode);
   assert(SCIPnodeGetType(tree->actnode) == SCIP_NODETYPE_ACTNODE);

   /* create temporary node for probing */
   CHECK_OKAY( nodeCreateProbingNode(&tree->probingnode, memhdr, set, tree) );
   assert(tree->probingnode != NULL);
   assert(SCIPnodeGetType(tree->probingnode) == SCIP_NODETYPE_PROBINGNODE);
   assert(SCIPnodeGetDepth(tree->probingnode) == tree->pathlen);

   /* create the new active path */
   CHECK_OKAY( treeEnsurePathMem(tree, set, tree->pathlen+1) );
   tree->probingnode->active = TRUE;
   tree->path[tree->pathlen] = tree->probingnode;
   tree->pathlen++;

   /* count the new LP sizes of the path */
   treeUpdatePathLPSize(tree, tree->pathlen-1);

   return SCIP_OKAY;
}

/** switches back from probing to normal operation mode, restores bounds of all variables and restores active constraints
 *  arrays of active node
 */
RETCODE SCIPtreeEndProbing(
   TREE*            tree,               /**< branch and bound tree */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   assert(tree != NULL);
   assert(tree->probingnode != NULL);
   assert(tree->actnode != NULL);
   assert(tree->probingnode->parent == tree->actnode);
   assert(SCIPnodeGetDepth(tree->probingnode) == SCIPnodeGetDepth(tree->actnode)+1);
   assert(tree->pathlen >= 2);
   assert(tree->path[tree->pathlen-1] == tree->probingnode);
   assert(tree->path[tree->pathlen-2] == tree->actnode);
   assert(SCIPnodeGetType(tree->probingnode) == SCIP_NODETYPE_PROBINGNODE);
   assert(SCIPnodeGetType(tree->actnode) == SCIP_NODETYPE_ACTNODE);

   /* undo the domain and constraint set changes of the temporary probing node */
   CHECK_OKAY( nodeDeactivate(&tree->probingnode, memhdr, set, stat, tree, lp, branchcand, eventqueue) );
   assert(tree->probingnode == NULL);

   /* remove the probing node from the active path */
   tree->pathlen--;
   assert(tree->pathlen == tree->actnode->depth+1);

   return SCIP_OKAY;
}


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets number of leaves */
int SCIPtreeGetNLeaves(
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqLen(tree->leaves);
}
   
/** gets number of nodes (children + siblings + leaves) */
int SCIPtreeGetNNodes(
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->nchildren + tree->nsiblings + SCIPtreeGetNLeaves(tree);
}

/** gets active node of the tree */
NODE* SCIPtreeGetActiveNode(
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->actnode;
}

/** gets depth of active node in the tree */
int SCIPtreeGetActiveDepth(
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->actnode != NULL);

   return tree->actnode != NULL ? tree->actnode->depth : -1;
}

/** returns, whether the LP was or is to be solved in the active node */
Bool SCIPtreeHasActiveNodeLP(
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->actnodehaslp;
}

/** sets mark to solve or to ignore the LP while processing the active node */
void SCIPtreeSetActiveNodeLP(
   TREE*            tree,               /**< branch and bound tree */
   Bool             solvelp             /**< should the LP be solved in active node? */
   )
{
   assert(tree != NULL);

   tree->actnodehaslp = solvelp;
}

/** gets current node of the tree */
NODE* SCIPtreeGetCurrentNode(
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->probingnode != NULL ? tree->probingnode : SCIPtreeGetActiveNode(tree);
}

/** gets depth of current node in the tree */
int SCIPtreeGetCurrentDepth(
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->probingnode != NULL ? tree->probingnode->depth : SCIPtreeGetActiveDepth(tree);
}

/** returns, whether the LP was or is to be solved in the current node */
Bool SCIPtreeHasCurrentNodeLP(
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->probingnode != NULL ? FALSE : SCIPtreeHasActiveNodeLP(tree);
}

/** returns whether the current node is a temporary probing node */
Bool SCIPtreeProbing(
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->probingnode == NULL || tree->probingnode->nodetype == SCIP_NODETYPE_PROBINGNODE);

   return (tree->probingnode != NULL);
}

/** returns the temporary probing node, or NULL if the current node is not the probing node */
NODE* SCIPtreeGetProbingNode(
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->probingnode == NULL || tree->probingnode->nodetype == SCIP_NODETYPE_PROBINGNODE);

   return tree->probingnode;
}

#endif


/** gets the best child of the active node w.r.t. the node selection priority assigned by the branching rule */
NODE* SCIPtreeGetPrioChild(
   TREE*            tree                /**< branch and bound tree */
   )
{
   NODE* bestnode;
   Real bestprio;
   int i;

   assert(tree != NULL);

   bestnode = NULL;
   bestprio = REAL_MIN;
   for( i = 0; i < tree->nchildren; ++i )
   {
      if( tree->childrenprio[i] > bestprio )
      {
         bestnode = tree->children[i];
         bestprio = tree->childrenprio[i];
      }
   }
   assert((tree->nchildren == 0) == (bestnode == NULL));

   return bestnode;
}

/** gets the best sibling of the active node w.r.t. the node selection priority assigned by the branching rule */
NODE* SCIPtreeGetPrioSibling(
   TREE*            tree                /**< branch and bound tree */
   )
{
   NODE* bestnode;
   Real bestprio;
   int i;

   assert(tree != NULL);

   bestnode = NULL;
   bestprio = REAL_MIN;
   for( i = 0; i < tree->nsiblings; ++i )
   {
      if( tree->siblingsprio[i] > bestprio )
      {
         bestnode = tree->siblings[i];
         bestprio = tree->siblingsprio[i];
      }
   }
   assert((tree->nsiblings == 0) == (bestnode == NULL));

   return bestnode;
}

/** gets the best child of the active node w.r.t. the node selection strategy */
NODE* SCIPtreeGetBestChild(
   TREE*            tree,               /**< branch and bound tree */
   SET*             set                 /**< global SCIP settings */
   )
{
   NODESEL* nodesel;
   NODE* bestnode;
   int i;

   assert(tree != NULL);

   nodesel = SCIPnodepqGetNodesel(tree->leaves);
   assert(nodesel != NULL);

   bestnode = NULL;
   for( i = 0; i < tree->nchildren; ++i )
   {
      if( bestnode == NULL || SCIPnodeselCompare(nodesel, set, tree->children[i], bestnode) < 0 )
      {
         bestnode = tree->children[i];
      }
   }

   return bestnode;
}

/** gets the best sibling of the active node w.r.t. the node selection strategy */
NODE* SCIPtreeGetBestSibling(
   TREE*            tree,               /**< branch and bound tree */
   SET*             set                 /**< global SCIP settings */
   )
{
   NODESEL* nodesel;
   NODE* bestnode;
   int i;

   assert(tree != NULL);

   nodesel = SCIPnodepqGetNodesel(tree->leaves);
   assert(nodesel != NULL);

   bestnode = NULL;
   for( i = 0; i < tree->nsiblings; ++i )
   {
      if( bestnode == NULL || SCIPnodeselCompare(nodesel, set, tree->siblings[i], bestnode) < 0 )
      {
         bestnode = tree->siblings[i];
      }
   }
   
   return bestnode;
}

/** gets the best leaf from the node queue w.r.t. the node selection strategy */
NODE* SCIPtreeGetBestLeaf(
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqFirst(tree->leaves);
}

/** gets the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy */
NODE* SCIPtreeGetBestNode(
   TREE*            tree,               /**< branch and bound tree */
   SET*             set                 /**< global SCIP settings */
   )
{
   NODESEL* nodesel;
   NODE* bestchild;
   NODE* bestsibling;
   NODE* bestleaf;
   NODE* bestnode;

   assert(tree != NULL);

   nodesel = SCIPnodepqGetNodesel(tree->leaves);
   assert(nodesel != NULL);

   /* get the best child, sibling, and leaf */
   bestchild = SCIPtreeGetBestChild(tree, set);
   bestsibling = SCIPtreeGetBestSibling(tree, set);
   bestleaf = SCIPtreeGetBestLeaf(tree);

   /* return the best of the three */
   bestnode = bestchild;
   if( bestsibling != NULL && (bestnode == NULL || SCIPnodeselCompare(nodesel, set, bestsibling, bestnode) < 0) )
      bestnode = bestsibling;
   if( bestleaf != NULL && (bestnode == NULL || SCIPnodeselCompare(nodesel, set, bestleaf, bestnode) < 0) )
      bestnode = bestleaf;

   assert(SCIPtreeGetNLeaves(tree) == 0 || bestnode != NULL);

   return bestnode;
}

/** gets the minimal lower bound of all nodes in the tree */
Real SCIPtreeGetLowerbound(
   TREE*            tree,               /**< branch and bound tree */
   SET*             set                 /**< global SCIP settings */
   )
{
   Real lowerbound;
   int i;

   assert(tree != NULL);
   assert(set != NULL);

   /* get the lower bound from the queue */
   lowerbound = SCIPnodepqGetLowerbound(tree->leaves, set);

   /* compare lower bound with children */
   for( i = 0; i < tree->nchildren; ++i )
   {
      assert(tree->children[i] != NULL);
      lowerbound = MIN(lowerbound, tree->children[i]->lowerbound); 
   }

   /* compare lower bound with siblings */
   for( i = 0; i < tree->nsiblings; ++i )
   {
      assert(tree->siblings[i] != NULL);
      lowerbound = MIN(lowerbound, tree->siblings[i]->lowerbound); 
   }

   /* compare lower bound with active node */
   if( tree->actnode != NULL )
   {
      lowerbound = MIN(lowerbound, tree->actnode->lowerbound);
   }

   return lowerbound;
}

/** gets the node with minimal lower bound of all nodes in the tree (child, sibling, or leaf) */
NODE* SCIPtreeGetLowerboundNode(
   TREE*            tree,               /**< branch and bound tree */
   SET*             set                 /**< global SCIP settings */
   )
{
   NODE* lowerboundnode;
   Real lowerbound;
   Real bestprio;
   int i;

   assert(tree != NULL);
   assert(set != NULL);

   /* get the lower bound from the queue */
   lowerboundnode = SCIPnodepqGetLowerboundNode(tree->leaves, set);
   lowerbound = lowerboundnode != NULL ? lowerboundnode->lowerbound : set->infinity;
   bestprio = -set->infinity;

   /* compare lower bound with children */
   for( i = 0; i < tree->nchildren; ++i )
   {
      assert(tree->children[i] != NULL);
      if( SCIPsetIsLE(set, tree->children[i]->lowerbound, lowerbound) )
      {
         if( SCIPsetIsLT(set, tree->children[i]->lowerbound, lowerbound) || tree->childrenprio[i] < bestprio )
         {
            lowerboundnode = tree->children[i]; 
            lowerbound = lowerboundnode->lowerbound; 
            bestprio = tree->childrenprio[i];
         }
      }
   }

   /* compare lower bound with siblings */
   for( i = 0; i < tree->nsiblings; ++i )
   {
      assert(tree->siblings[i] != NULL);
      if( SCIPsetIsLE(set, tree->siblings[i]->lowerbound, lowerbound) )
      {
         if( SCIPsetIsLT(set, tree->siblings[i]->lowerbound, lowerbound) || tree->siblingsprio[i] < bestprio )
         {
            lowerboundnode = tree->siblings[i]; 
            lowerbound = lowerboundnode->lowerbound; 
            bestprio = tree->siblingsprio[i];
         }
      }
   }

   return lowerboundnode;
}

/** gets the average lower bound of all nodes in the tree */
Real SCIPtreeGetAvgLowerbound(
   TREE*            tree,               /**< branch and bound tree */
   Real             cutoffbound         /**< global cutoff bound */
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
   if( tree->actnode != NULL && tree->actnode->lowerbound < cutoffbound )
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

