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
#pragma ident "@(#) $Id: tree.c,v 1.165 2006/01/03 11:36:49 bzfpfend Exp $"

/**@file   tree.c
 * @brief  methods for branch and bound tree
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/vbc.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/implics.h"
#include "scip/primal.h"
#include "scip/tree.h"
#include "scip/solve.h"
#include "scip/cons.h"
#include "scip/nodesel.h"
#include "scip/prop.h"
#include "scip/debug.h"


#define MAXDEPTH          65535  /**< maximal depth level for nodes; must correspond to node data structure */
#define MAXREPROPMARK       511  /**< maximal subtree repropagation marker; must correspond to node data structure */


/*
 * dynamic memory arrays
 */

/** resizes children arrays to be able to store at least num nodes */
static
SCIP_RETCODE treeEnsureChildrenMem(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of node slots in array */
   )
{
   assert(tree != NULL);
   assert(set != NULL);

   if( num > tree->childrensize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&tree->children, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&tree->childrenprio, newsize) );
      tree->childrensize = newsize;
   }
   assert(num <= tree->childrensize);

   return SCIP_OKAY;
}

/** resizes path array to be able to store at least num nodes */
static
SCIP_RETCODE treeEnsurePathMem(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of node slots in path */
   )
{
   assert(tree != NULL);
   assert(set != NULL);

   if( num > tree->pathsize )
   {
      int newsize;

      newsize = SCIPsetCalcPathGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&tree->path, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&tree->pathnlpcols, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&tree->pathnlprows, newsize) );
      tree->pathsize = newsize;
   }
   assert(num <= tree->pathsize);

   return SCIP_OKAY;
}




/*
 * Node methods
 */

/** node comparator for best lower bound */
SCIP_DECL_SORTPTRCOMP(SCIPnodeCompLowerbound)
{  /*lint --e{715}*/
   assert(elem1 != NULL);
   assert(elem2 != NULL);

   if( ((SCIP_NODE*)elem1)->lowerbound < ((SCIP_NODE*)elem2)->lowerbound )
      return -1;
   else if( ((SCIP_NODE*)elem1)->lowerbound > ((SCIP_NODE*)elem2)->lowerbound )
      return +1;
   else
      return 0;
}

/** increases the reference counter of the LP state in the fork */
static
void forkCaptureLPIState(
   SCIP_FORK*            fork,               /**< fork data */
   int                   nuses               /**< number to add to the usage counter */
   )
{
   assert(fork != NULL);
   assert(fork->nlpistateref >= 0);
   assert(nuses > 0);

   fork->nlpistateref += nuses;
   SCIPdebugMessage("captured LPI state of fork %p %d times -> new nlpistateref=%d\n", fork, nuses, fork->nlpistateref);
}

/** decreases the reference counter of the LP state in the fork */
static
SCIP_RETCODE forkReleaseLPIState(
   SCIP_FORK*            fork,               /**< fork data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(fork != NULL);
   assert(fork->nlpistateref > 0);
   assert(blkmem != NULL);
   assert(lp != NULL);

   fork->nlpistateref--;
   if( fork->nlpistateref == 0 )
   {
      SCIP_CALL( SCIPlpFreeState(lp, blkmem, &(fork->lpistate)) );
   }

   SCIPdebugMessage("released LPI state of fork %p -> new nlpistateref=%d\n", fork, fork->nlpistateref);

   return SCIP_OKAY;
}

/** increases the reference counter of the LP state in the subroot */
static
void subrootCaptureLPIState(
   SCIP_SUBROOT*         subroot,            /**< subroot data */
   int                   nuses               /**< number to add to the usage counter */
   )
{
   assert(subroot != NULL);
   assert(subroot->nlpistateref >= 0);
   assert(nuses > 0);

   subroot->nlpistateref += nuses;
   SCIPdebugMessage("captured LPI state of subroot %p %d times -> new nlpistateref=%d\n", 
      subroot, nuses, subroot->nlpistateref);
}

/** decreases the reference counter of the LP state in the subroot */
static
SCIP_RETCODE subrootReleaseLPIState(
   SCIP_SUBROOT*         subroot,            /**< subroot data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(subroot != NULL);
   assert(subroot->nlpistateref > 0);
   assert(blkmem != NULL);
   assert(lp != NULL);

   subroot->nlpistateref--;
   if( subroot->nlpistateref == 0 )
   {
      SCIP_CALL( SCIPlpFreeState(lp, blkmem, &(subroot->lpistate)) );
   }
   
   SCIPdebugMessage("released LPI state of subroot %p -> new nlpistateref=%d\n", subroot, subroot->nlpistateref);

   return SCIP_OKAY;
}

/** increases the reference counter of the LP state in the fork or subroot node */
void SCIPnodeCaptureLPIState(
   SCIP_NODE*            node,               /**< fork/subroot node */
   int                   nuses               /**< number to add to the usage counter */
   )
{
   assert(node != NULL);

   SCIPdebugMessage("capture %d times LPI state of node %p at depth %d (current: %d)\n", nuses, node, node->depth, 
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
      SCIPerrorMessage("node for capturing the LPI state is neither fork nor subroot\n");
      SCIPABORT();
   }  /*lint !e788*/
}

/** decreases the reference counter of the LP state in the fork or subroot node */
SCIP_RETCODE SCIPnodeReleaseLPIState(
   SCIP_NODE*            node,               /**< fork/subroot node */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(node != NULL);

   SCIPdebugMessage("release LPI state of node %p at depth %d (current: %d)\n", node, node->depth,
      SCIPnodeGetType(node) == SCIP_NODETYPE_FORK ? node->data.fork->nlpistateref : node->data.subroot->nlpistateref);
   switch( SCIPnodeGetType(node) )
   {  
   case SCIP_NODETYPE_FORK:
      return forkReleaseLPIState(node->data.fork, blkmem, lp);
   case SCIP_NODETYPE_SUBROOT:
      return subrootReleaseLPIState(node->data.subroot, blkmem, lp);
   default:
      SCIPerrorMessage("node for releasing the LPI state is neither fork nor subroot\n");
      return SCIP_INVALIDDATA;
   }  /*lint !e788*/
}

/** initializes junction data */
static
SCIP_RETCODE junctionInit(
   SCIP_JUNCTION*        junction,           /**< pointer to junction data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(junction != NULL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->focusnode != NULL);

   junction->nchildren = tree->nchildren;

   /* increase the LPI state usage counter of the current LP fork */
   if( tree->focuslpfork != NULL )
      SCIPnodeCaptureLPIState(tree->focuslpfork, tree->nchildren);

   return SCIP_OKAY;
}

/** creates fork data */
static
SCIP_RETCODE forkCreate(
   SCIP_FORK**           fork,               /**< pointer to fork data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(fork != NULL);
   assert(blkmem != NULL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->focusnode != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, fork) );

   SCIP_CALL( SCIPlpGetState(lp, blkmem, &((*fork)->lpistate)) );
   (*fork)->nlpistateref = 0;
   (*fork)->addedcols = NULL;
   (*fork)->addedrows = NULL;
   (*fork)->naddedcols = SCIPlpGetNNewcols(lp);
   (*fork)->naddedrows = SCIPlpGetNNewrows(lp);
   (*fork)->nchildren = tree->nchildren;

   SCIPdebugMessage("creating fork information with %d children (%d new cols, %d new rows)\n",
      (*fork)->nchildren, (*fork)->naddedcols, (*fork)->naddedrows);

   if( (*fork)->naddedcols > 0 )
   {
      /* copy the newly created columns to the fork's col array */
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*fork)->addedcols, SCIPlpGetNewcols(lp), (*fork)->naddedcols) );
   }
   if( (*fork)->naddedrows > 0 )
   {
      int i;
      
      /* copy the newly created rows to the fork's row array */
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*fork)->addedrows, SCIPlpGetNewrows(lp), (*fork)->naddedrows) );

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
SCIP_RETCODE forkFree(
   SCIP_FORK**           fork,               /**< fork data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(fork != NULL);
   assert(*fork != NULL);
   assert((*fork)->nchildren == 0);
   assert((*fork)->nlpistateref == 0);
   assert((*fork)->lpistate == NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   /* release the added rows */
   for( i = 0; i < (*fork)->naddedrows; ++i )
   {
      SCIP_CALL( SCIProwRelease(&(*fork)->addedrows[i], blkmem, set, lp) );
   }

   BMSfreeBlockMemoryArrayNull(blkmem, &(*fork)->addedcols, (*fork)->naddedcols);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*fork)->addedrows, (*fork)->naddedrows);
   BMSfreeBlockMemory(blkmem, fork);

   return SCIP_OKAY;
}

/** creates subroot data */
static
SCIP_RETCODE subrootCreate(
   SCIP_SUBROOT**        subroot,            /**< pointer to subroot data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;
      
   assert(subroot != NULL);
   assert(blkmem != NULL);
   assert(tree != NULL);
   assert(tree->nchildren > 0);
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->focusnode != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, subroot) );
   
   (*subroot)->nlpistateref = 0;
   (*subroot)->ncols = SCIPlpGetNCols(lp);
   (*subroot)->nrows = SCIPlpGetNRows(lp);
   (*subroot)->nchildren = tree->nchildren;
   SCIP_CALL( SCIPlpGetState(lp, blkmem, &((*subroot)->lpistate)) );

   if( (*subroot)->ncols != 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*subroot)->cols, SCIPlpGetCols(lp), (*subroot)->ncols) );
   }
   else
      (*subroot)->cols = NULL;
   if( (*subroot)->nrows != 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*subroot)->rows, SCIPlpGetRows(lp), (*subroot)->nrows) );
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
SCIP_RETCODE subrootFree(
   SCIP_SUBROOT**        subroot,            /**< subroot data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;
      
   assert(subroot != NULL);
   assert(*subroot != NULL);
   assert((*subroot)->nchildren == 0);
   assert((*subroot)->nlpistateref == 0);
   assert((*subroot)->lpistate == NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   /* release the rows of the subroot */
   for( i = 0; i < (*subroot)->nrows; ++i )
   {
      SCIP_CALL( SCIProwRelease(&(*subroot)->rows[i], blkmem, set, lp) );
   }

   BMSfreeBlockMemoryArrayNull(blkmem, &(*subroot)->cols, (*subroot)->ncols);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*subroot)->rows, (*subroot)->nrows);
   BMSfreeBlockMemory(blkmem, subroot);

   return SCIP_OKAY;
}

/** removes given sibling node from the siblings array */
static
void treeRemoveSibling(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            sibling             /**< sibling node to remove */
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

/** adds given child node to children array of focus node */
static
SCIP_RETCODE treeAddChild(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            child,              /**< child node to add */
   SCIP_Real             nodeselprio         /**< node selection priority of child node */
   )
{
   assert(tree != NULL);
   assert(child != NULL);
   assert(SCIPnodeGetType(child) == SCIP_NODETYPE_CHILD);
   assert(child->data.child.arraypos == -1);

   SCIP_CALL( treeEnsureChildrenMem(tree, set, tree->nchildren+1) );
   tree->children[tree->nchildren] = child;
   tree->childrenprio[tree->nchildren] = nodeselprio;
   child->data.child.arraypos = tree->nchildren;
   tree->nchildren++;

   return SCIP_OKAY;
}

/** removes given child node from the children array */
static
void treeRemoveChild(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            child               /**< child node to remove */
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

/** makes node a child of the given parent node, which must be the focus node; if the child is a probing node,
 *  the parent node can also be a refocused node or a probing node
 */
static
SCIP_RETCODE nodeAssignParent(
   SCIP_NODE*            node,               /**< child node */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            parent,             /**< parent (= focus) node (or NULL, if node is root) */
   SCIP_Real             nodeselprio         /**< node selection priority of child node */
   )
{
   assert(node != NULL);
   assert(node->parent == NULL);
   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);
   assert(node->conssetchg == NULL);
   assert(node->domchg == NULL);
   assert(SCIPsetIsInfinity(set, -node->lowerbound)); /* node was just created */
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1] == parent);
   assert(parent == tree->focusnode || SCIPnodeGetType(parent) == SCIP_NODETYPE_PROBINGNODE);
   assert(parent == NULL || SCIPnodeGetType(parent) == SCIP_NODETYPE_FOCUSNODE
      || (SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE
         && (SCIPnodeGetType(parent) == SCIP_NODETYPE_REFOCUSNODE
            || SCIPnodeGetType(parent) == SCIP_NODETYPE_PROBINGNODE)));

   /* link node to parent */
   node->parent = parent;
   if( parent != NULL )
   {
      node->lowerbound = parent->lowerbound;
      node->depth = parent->depth+1;
      if( parent->depth >= MAXDEPTH-1 )
      {
         SCIPerrorMessage("maximal depth level exceeded\n");
         return SCIP_MAXDEPTHLEVEL;
      }
   }
   SCIPdebugMessage("assigning parent %p to node %p at depth %d\n", parent, node, node->depth);

   /* register node in the childlist of the focus (the parent) node */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD )
   {
      assert(parent == NULL || SCIPnodeGetType(parent) == SCIP_NODETYPE_FOCUSNODE);
      SCIP_CALL( treeAddChild(tree, set, node, nodeselprio) );
   }

   node->priority = nodeselprio;

   return SCIP_OKAY;
}

/** decreases number of children of the parent, frees it if no children are left */
static
SCIP_RETCODE nodeReleaseParent(
   SCIP_NODE*            node,               /**< child node */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_NODE* parent;

   assert(node != NULL);
   assert(blkmem != NULL);
   
   SCIPdebugMessage("releasing parent-child relationship of node %p at depth %d of type %d with parent %p of type %d\n",
      node, node->depth, SCIPnodeGetType(node), node->parent,
      node->parent != NULL ? (int)SCIPnodeGetType(node->parent) : -1);
   parent = node->parent;
   if( parent != NULL )
   {
      SCIP_Bool freeParent;

      freeParent = FALSE;
      switch( SCIPnodeGetType(parent) )
      {
      case SCIP_NODETYPE_FOCUSNODE:
         assert(parent->active);
         assert(SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE
            || SCIPnodeGetType(node) == SCIP_NODETYPE_LEAF);
         if( SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD )
            treeRemoveChild(tree, node);
         freeParent = FALSE; /* don't kill the focus node at this point */
         break;
      case SCIP_NODETYPE_PROBINGNODE:
         assert(SCIPtreeProbing(tree));
         freeParent = FALSE; /* probing nodes have to be freed individually */
         break;
      case SCIP_NODETYPE_SIBLING:
         SCIPerrorMessage("sibling cannot be a parent node\n");
         return SCIP_INVALIDDATA;
      case SCIP_NODETYPE_CHILD:
         SCIPerrorMessage("child cannot be a parent node\n");
         return SCIP_INVALIDDATA;
      case SCIP_NODETYPE_LEAF:
         SCIPerrorMessage("leaf cannot be a parent node\n");
         return SCIP_INVALIDDATA;
      case SCIP_NODETYPE_DEADEND:
         SCIPerrorMessage("deadend cannot be a parent node\n");
         return SCIP_INVALIDDATA;
      case SCIP_NODETYPE_JUNCTION:
         assert(parent->data.junction.nchildren > 0);
         parent->data.junction.nchildren--;
         freeParent = (parent->data.junction.nchildren == 0); /* free parent if it has no more children */
         break;
      case SCIP_NODETYPE_FORK:
         assert(parent->data.fork != NULL);
         assert(parent->data.fork->nchildren > 0);
         parent->data.fork->nchildren--;
         freeParent = (parent->data.fork->nchildren == 0); /* free parent if it has no more children */
         break;
      case SCIP_NODETYPE_SUBROOT:
         assert(parent->data.subroot != NULL);
         assert(parent->data.subroot->nchildren > 0);
         parent->data.subroot->nchildren--;
         freeParent = (parent->data.subroot->nchildren == 0); /* free parent if it has no more children */
         break;
      case SCIP_NODETYPE_REFOCUSNODE:
         /* the only possible child a refocused node can have in its refocus state is the probing root node;
          * we don't want to free the refocused node, because we first have to convert it back to its original
          * type (where it possibly has children)
          */
         assert(SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);
         assert(!SCIPtreeProbing(tree));
         freeParent = FALSE;
         break;
      default:
         SCIPerrorMessage("unknown node type\n");
         return SCIP_INVALIDDATA;
      }

      /* free parent, if it is not on the current active path */
      if( freeParent && !parent->active )
      {
         SCIP_CALL( SCIPnodeFree(&node->parent, blkmem, set, tree, lp) );
      }
   }

   return SCIP_OKAY;
}

/** creates a node data structure */
static
SCIP_RETCODE nodeCreate(
   SCIP_NODE**           node,               /**< pointer to node data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(node != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, node) );
   (*node)->parent = NULL;
   (*node)->conssetchg = NULL;
   (*node)->domchg = NULL;
   (*node)->number = 0;
   (*node)->lowerbound = -SCIPsetInfinity(set);
   (*node)->priority = 0.0;
   (*node)->depth = 0;
   (*node)->active = FALSE;
   (*node)->cutoff = FALSE;
   (*node)->reprop = FALSE;
   (*node)->repropsubtreemark = 0;

   return SCIP_OKAY;
}

/** creates a child node of the focus node */
SCIP_RETCODE SCIPnodeCreateChild(
   SCIP_NODE**           node,               /**< pointer to node data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Real             nodeselprio         /**< node selection priority of new node */
   )
{
   assert(node != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->pathlen == 0 || tree->path != NULL);
   assert((tree->pathlen == 0) == (tree->focusnode == NULL));
   assert(tree->focusnode == NULL || tree->focusnode == tree->path[tree->pathlen-1]);
   assert(tree->focusnode == NULL || SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);

   stat->ncreatednodes++;
   stat->ncreatednodesrun++;

   /* create the node data structure */
   SCIP_CALL( nodeCreate(node, blkmem, set) );
   (*node)->number = stat->ncreatednodesrun;

   /* mark node to be a child node */
   (*node)->nodetype = SCIP_NODETYPE_CHILD; /*lint !e641*/
   (*node)->data.child.arraypos = -1;

   /* make focus node the parent of the new child */
   SCIP_CALL( nodeAssignParent(*node, blkmem, set, tree, tree->focusnode, nodeselprio) );
   
   /* output node creation to VBC file */
   SCIP_CALL( SCIPvbcNewChild(stat->vbc, stat, *node) );

   SCIPdebugMessage("created child node %p at depth %d\n", *node, (*node)->depth);

   return SCIP_OKAY;
}

/** frees node */
SCIP_RETCODE SCIPnodeFree(
   SCIP_NODE**           node,               /**< node data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Bool isroot;

   assert(node != NULL);
   assert(*node != NULL);
   assert(!(*node)->active);
   assert(blkmem != NULL);
   assert(tree != NULL);

   SCIPdebugMessage("free node %p at depth %d of type %d\n", *node, (*node)->depth, SCIPnodeGetType(*node));

   /* inform solution debugger, that the node has been freed */
   SCIP_CALL( SCIPdebugRemoveNode(*node) ); /*lint !e506 !e774*/

   /* free nodetype specific data, and release no longer needed LPI states */
   switch( SCIPnodeGetType(*node) )
   {
   case SCIP_NODETYPE_FOCUSNODE:
      assert(tree->focusnode == *node);
      assert(!SCIPtreeProbing(tree));
      SCIPerrorMessage("cannot free focus node - has to be converted into a dead end first\n");
      return SCIP_INVALIDDATA;
   case SCIP_NODETYPE_PROBINGNODE:
      assert(SCIPtreeProbing(tree));
      assert(SCIPnodeGetDepth(tree->probingroot) <= SCIPnodeGetDepth(*node));
      assert(SCIPnodeGetDepth(*node) > 0);
      break;
   case SCIP_NODETYPE_SIBLING:
      assert((*node)->data.sibling.arraypos >= 0);
      assert((*node)->data.sibling.arraypos < tree->nsiblings);
      assert(tree->siblings[(*node)->data.sibling.arraypos] == *node);
      if( tree->focuslpfork != NULL )
      {
         assert(SCIPnodeGetType(tree->focuslpfork) == SCIP_NODETYPE_FORK
            || SCIPnodeGetType(tree->focuslpfork) == SCIP_NODETYPE_SUBROOT);
         SCIP_CALL( SCIPnodeReleaseLPIState(tree->focuslpfork, blkmem, lp) );
      }
      treeRemoveSibling(tree, *node);
      break;
   case SCIP_NODETYPE_CHILD:
      assert((*node)->data.child.arraypos >= 0);
      assert((*node)->data.child.arraypos < tree->nchildren);
      assert(tree->children[(*node)->data.child.arraypos] == *node);
      /* The children capture the LPI state at the moment, where the focus node is
       * converted into a junction, fork, or subroot, and a new node is focused.
       * At the same time, they become siblings or leaves, such that freeing a child
       * of the focus node doesn't require to release the LPI state;
       * we don't need to call treeRemoveChild(), because this is done in nodeReleaseParent()
       */
      break;
   case SCIP_NODETYPE_LEAF:
      if( (*node)->data.leaf.lpfork != NULL )
      {
         SCIP_CALL( SCIPnodeReleaseLPIState((*node)->data.leaf.lpfork, blkmem, lp) );
      }
      break;
   case SCIP_NODETYPE_DEADEND:
   case SCIP_NODETYPE_JUNCTION:
      break;
   case SCIP_NODETYPE_FORK:
      SCIP_CALL( forkFree(&((*node)->data.fork), blkmem, set, lp) );
      break;
   case SCIP_NODETYPE_SUBROOT:
      SCIP_CALL( subrootFree(&((*node)->data.subroot), blkmem, set, lp) );
      break;
   case SCIP_NODETYPE_REFOCUSNODE:
      SCIPerrorMessage("cannot free node as long it is refocused\n");
      return SCIP_INVALIDDATA;
   default:
      SCIPerrorMessage("unknown node type\n");
      return SCIP_INVALIDDATA;
   }

   /* check, if the node to be freed is the root node */
   isroot = (SCIPnodeGetDepth(*node) == 0);

   /* free common data */
   SCIP_CALL( SCIPconssetchgFree(&(*node)->conssetchg, blkmem, set) );
   SCIP_CALL( SCIPdomchgFree(&(*node)->domchg, blkmem, set) );
   SCIP_CALL( nodeReleaseParent(*node, blkmem, set, tree, lp) );

   /* check, if the node is the current probing root */
   if( *node == tree->probingroot )
   {
      assert(SCIPnodeGetType(*node) == SCIP_NODETYPE_PROBINGNODE);
      tree->probingroot = NULL;
   }

   BMSfreeBlockMemory(blkmem, node);

   /* delete the tree's root node pointer, if the freed node was the root */
   if( isroot )
      tree->root = NULL;

   return SCIP_OKAY;
}

/** cuts off node and whole sub tree from branch and bound tree */
void SCIPnodeCutoff(
   SCIP_NODE*            node,               /**< node that should be cut off */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(node != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);

   node->cutoff = TRUE;
   node->lowerbound = SCIPsetInfinity(set);
   if( node->active )
      tree->cutoffdepth = MIN(tree->cutoffdepth, (int)node->depth);

   SCIPvbcCutoffNode(stat->vbc, stat, node);

   SCIPdebugMessage("cutting off %s node %p at depth %d (cutoffdepth: %d)\n", 
      node->active ? "active" : "inactive", node, node->depth, tree->cutoffdepth);
}

/** marks node, that propagation should be applied again the next time, a node of its subtree is focused */
void SCIPnodePropagateAgain(
   SCIP_NODE*            node,               /**< node that should be propagated again */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(node != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);

   node->reprop = TRUE;
   if( node->active )
      tree->repropdepth = MIN(tree->repropdepth, (int)node->depth);

   SCIPvbcMarkedRepropagateNode(stat->vbc, stat, node);

   SCIPdebugMessage("marked %s node %p at depth %d to be propagated again (repropdepth: %d)\n", 
      node->active ? "active" : "inactive", node, node->depth, tree->repropdepth);
}

/** marks node, that it is completely propagated in the current repropagation subtree level */
void SCIPnodeMarkPropagated(
   SCIP_NODE*            node,               /**< node that should be marked to be propagated */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(node != NULL);
   assert(tree != NULL);

   if( node->parent != NULL )
      node->repropsubtreemark = node->parent->repropsubtreemark;
   node->reprop = FALSE;

   /* if the node was the highest repropagation node in the path, update the repropdepth in the tree data */
   if( node->active && node->depth == tree->repropdepth )
   {
      do
      {
         assert(tree->repropdepth < tree->pathlen);
         assert(tree->path[tree->repropdepth]->active);
         assert(!tree->path[tree->repropdepth]->reprop);
         tree->repropdepth++;
      }
      while( tree->repropdepth < tree->pathlen && !tree->path[tree->repropdepth]->reprop );
      if( tree->repropdepth == tree->pathlen )
         tree->repropdepth = INT_MAX;
   }
}

/** moves the subtree repropagation counter to the next value */
static
void treeNextRepropsubtreecount(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   tree->repropsubtreecount++;
   tree->repropsubtreecount %= (MAXREPROPMARK+1);
}

/** applies propagation on the node, that was marked to be propagated again */
static
SCIP_RETCODE nodeRepropagate(
   SCIP_NODE*            node,               /**< node to apply propagation on */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_NODETYPE oldtype;
   SCIP_NODE* oldfocusnode;
   SCIP_NODE* oldfocuslpfork;
   SCIP_NODE* oldfocussubroot;
   int oldfocuslpforklpcount;
   int oldnchildren;
   int oldnsiblings;
   SCIP_Bool oldfocusnodehaslp;
   SCIP_Longint oldnboundchgs;
   SCIP_Bool initialreprop;
   SCIP_Bool clockisrunning;

   assert(node != NULL);
   assert((SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_FOCUSNODE
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_JUNCTION
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_FORK
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_SUBROOT);
   assert(node->active);
   assert(node->reprop || (node != NULL && node->repropsubtreemark != node->parent->repropsubtreemark));
   assert(stat != NULL);
   assert(tree != NULL);
   assert(SCIPeventqueueIsDelayed(eventqueue));
   assert(cutoff != NULL);

   SCIPdebugMessage("propagating again node %p at depth %d\n", node, node->depth);
   initialreprop = node->reprop;

   SCIPvbcRepropagatedNode(stat->vbc, stat, node);

   /* process the delayed events in order to flush the problem changes */
   SCIP_CALL( SCIPeventqueueProcess(eventqueue, blkmem, set, primal, lp, branchcand, eventfilter) );

   /* stop node activation timer */
   clockisrunning = SCIPclockIsRunning(stat->nodeactivationtime);
   if( clockisrunning )
      SCIPclockStop(stat->nodeactivationtime, set);

   /* mark the node refocused and temporarily install it as focus node */
   oldtype = (SCIP_NODETYPE)node->nodetype;
   oldfocusnode = tree->focusnode;
   oldfocuslpfork = tree->focuslpfork;
   oldfocussubroot = tree->focussubroot;
   oldfocuslpforklpcount = tree->focuslpforklpcount;
   oldnchildren = tree->nchildren;
   oldnsiblings = tree->nsiblings;
   oldfocusnodehaslp = tree->focusnodehaslp;
   node->nodetype = SCIP_NODETYPE_REFOCUSNODE; /*lint !e641*/
   tree->focusnode = node;
   tree->focuslpfork = NULL;
   tree->focussubroot = NULL;
   tree->focuslpforklpcount = -1;
   tree->nchildren = 0;
   tree->nsiblings = 0;
   tree->focusnodehaslp = FALSE;

   /* propagate the domains again */
   oldnboundchgs = stat->nboundchgs;
   SCIP_CALL( SCIPpropagateDomains(blkmem, set, stat, prob, tree, conflict, SCIPnodeGetDepth(node), 0, cutoff) );
   assert(!node->reprop);
   assert(node->parent == NULL || node->repropsubtreemark == node->parent->repropsubtreemark);
   assert((SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_REFOCUSNODE);
   assert(tree->focusnode == node);
   assert(tree->focuslpfork == NULL);
   assert(tree->focussubroot == NULL);
   assert(tree->focuslpforklpcount == -1);
   assert(tree->nchildren == 0);
   assert(tree->nsiblings == 0);
   assert(tree->focusnodehaslp == FALSE);
   assert(stat->nboundchgs >= oldnboundchgs);
   stat->nreprops++;
   stat->nrepropboundchgs += stat->nboundchgs - oldnboundchgs;

   SCIPdebugMessage("repropagation %"SCIP_LONGINT_FORMAT" at depth %d changed %"SCIP_LONGINT_FORMAT" bounds (total reprop bound changes: %"SCIP_LONGINT_FORMAT"), cutoff: %d\n",
      stat->nreprops, node->depth, stat->nboundchgs - oldnboundchgs, stat->nrepropboundchgs, *cutoff);

   /* if a propagation marked with the reprop flag was successful, we want to repropagate the whole subtree */
   /**@todo because repropsubtree is only a bit flag, we cannot mark a whole subtree a second time for
    *       repropagation; use a (small) part of the node's bits to be able to store larger numbers,
    *       and update tree->repropsubtreelevel with this number
    */
   if( initialreprop && !(*cutoff) && stat->nboundchgs > oldnboundchgs )
   {
      treeNextRepropsubtreecount(tree);
      node->repropsubtreemark = tree->repropsubtreecount; /*lint !e732*/
      SCIPdebugMessage("initial repropagation at depth %d changed %"SCIP_LONGINT_FORMAT" bounds -> repropagating subtree (new mark: %d)\n",
         node->depth, stat->nboundchgs - oldnboundchgs, tree->repropsubtreecount);
      assert((int)(node->repropsubtreemark) == tree->repropsubtreecount); /* bitfield must be large enough */
   }

   /* reset the node's type and reinstall the old focus node */
   node->nodetype = oldtype; /*lint !e641*/
   tree->focusnode = oldfocusnode;
   tree->focuslpfork = oldfocuslpfork;
   tree->focussubroot = oldfocussubroot;
   tree->focuslpforklpcount = oldfocuslpforklpcount;
   tree->nchildren = oldnchildren;
   tree->nsiblings = oldnsiblings;
   tree->focusnodehaslp = oldfocusnodehaslp;

   /* make the domain change data static again to save memory */
   if( (SCIP_NODETYPE)node->nodetype != SCIP_NODETYPE_FOCUSNODE )
   {
      SCIP_CALL( SCIPdomchgMakeStatic(&node->domchg, blkmem, set) );
   }

   /* start node activation timer again */
   if( clockisrunning )
      SCIPclockStart(stat->nodeactivationtime, set);

   /* delay events in path switching */
   SCIP_CALL( SCIPeventqueueDelay(eventqueue) );

   /* mark the node to be cut off if a cutoff was detected */
   if( *cutoff )
      SCIPnodeCutoff(node, set, stat, tree);

   return SCIP_OKAY;
}

/** informs node, that it is now on the active path and applies any domain and constraint set changes */
static
SCIP_RETCODE nodeActivate(
   SCIP_NODE*            node,               /**< node to activate */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   assert(node != NULL);
   assert(!node->active);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(cutoff != NULL);

   SCIPdebugMessage("activate node %p at depth %d of type %d (reprop subtree mark: %d)\n",
      node, node->depth, SCIPnodeGetType(node), node->repropsubtreemark);

   /* apply domain and constraint set changes */
   SCIP_CALL( SCIPconssetchgApply(node->conssetchg, blkmem, set, stat, node->depth) );
   SCIP_CALL( SCIPdomchgApply(node->domchg, blkmem, set, stat, lp, branchcand, eventqueue, node->depth, cutoff) );

   /* mark node active */
   node->active = TRUE;
   stat->nactivatednodes++;

   /* check if the domain change produced a cutoff */
   if( *cutoff )
   {
      /* try to repropagate the node to see, if the propagation also leads to a conflict and a conflict clause
       * could be generated; if propagation conflict analysis is turned off, repropagating the node makes no
       * sense, since it is already cut off
       */
      node->reprop = set->conf_useprop;

      /* mark the node to be cut off */
      SCIPnodeCutoff(node, set, stat, tree);
   }

   /* propagate node again, if the reprop flag is set; in the new focus node, no repropagation is necessary, because
    * the focus node is propagated anyways
    */
   if( SCIPnodeGetType(node) != SCIP_NODETYPE_FOCUSNODE
      && (node->reprop || (node->parent != NULL && node->repropsubtreemark != node->parent->repropsubtreemark)) )
   {
      SCIP_Bool propcutoff;

      SCIP_CALL( nodeRepropagate(node, blkmem, set, stat, prob, primal, tree, lp, branchcand, conflict, 
            eventfilter, eventqueue, &propcutoff) );
      *cutoff = *cutoff || propcutoff;
   }

   return SCIP_OKAY;
}

/** informs node, that it is no longer on the active path and undoes any domain and constraint set changes */
static
SCIP_RETCODE nodeDeactivate(
   SCIP_NODE*            node,               /**< node to deactivate */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   assert(node != NULL);
   assert(node->active);
   assert(tree != NULL);

   SCIPdebugMessage("deactivate node %p at depth %d of type %d  (reprop subtree mark: %d)\n",
      node, node->depth, SCIPnodeGetType(node), node->repropsubtreemark);

   /* undo domain and constraint set changes */
   SCIP_CALL( SCIPdomchgUndo(node->domchg, blkmem, set, stat, lp, branchcand, eventqueue) );
   SCIP_CALL( SCIPconssetchgUndo(node->conssetchg, blkmem, set, stat) );

   /* mark node inactive */
   node->active = FALSE;

   /* count number of deactivated nodes (ignoring probing switches) */
   if( !SCIPtreeProbing(tree) )
      stat->ndeactivatednodes++;

   return SCIP_OKAY;
}

/** adds constraint locally to the node and captures it; activates constraint, if node is active;
 *  if a local constraint is added to the root node, it is automatically upgraded into a global constraint
 */
SCIP_RETCODE SCIPnodeAddCons(
   SCIP_NODE*            node,               /**< node to add constraint to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   assert(node != NULL);
   assert(cons != NULL);
   assert(cons->validdepth <= node->depth);

   /* if node is the root, mark constraint to be globally valid */
   if( node->depth == 0 )
   {
      assert(node == tree->root);
      SCIPconsSetLocal(cons, FALSE);
   }

   /* add constraint addition to the node's constraint set change data, and activate constraint if node is active */
   SCIP_CALL( SCIPconssetchgAddAddedCons(&node->conssetchg, blkmem, set, stat, cons, node->depth, node->active) );
   assert(node->conssetchg != NULL);
   assert(node->conssetchg->addedconss != NULL);
   assert(!node->active || SCIPconsIsActive(cons));

   return SCIP_OKAY;
}

/** locally deletes constraint at the given node by disabling its separation, enforcing, and propagation capabilities
 *  at the node; captures constraint; disables constraint, if node is active
 */
SCIP_RETCODE SCIPnodeDelCons(
   SCIP_NODE*            node,               /**< node to add constraint to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_CONS*            cons                /**< constraint to locally delete */
   )
{
   assert(node != NULL);
   assert(tree != NULL);
   assert(cons != NULL);

   SCIPdebugMessage("disabling constraint <%s> at node at depth %d\n", cons->name, node->depth);

   /* add constraint disabling to the node's constraint set change data */
   SCIP_CALL( SCIPconssetchgAddDisabledCons(&node->conssetchg, blkmem, set, cons) );
   assert(node->conssetchg != NULL);
   assert(node->conssetchg->disabledconss != NULL);

   /* disable constraint, if node is active */
   if( node->active && cons->enabled && !cons->updatedisable )
   {
      SCIP_CALL( SCIPconsDisable(cons, set, stat) );
   }

   return SCIP_OKAY;
}

/** adds bound change with inference information to focus node, child of focus node, or probing node;
 *  if possible, adjusts bound to integral value;
 *  at most one of infercons and inferprop may be non-NULL
 */
SCIP_RETCODE SCIPnodeAddBoundinfer(
   SCIP_NODE*            node,               /**< node to add bound change to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change, or NULL */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             probingchange       /**< is the bound change a temporary setting due to probing? */
   )
{
   SCIP_VAR* infervar;
   SCIP_BOUNDTYPE inferboundtype;
   SCIP_Real oldlb;
   SCIP_Real oldub;
   SCIP_Real oldbound;

   assert(node != NULL);
   assert((SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_FOCUSNODE
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_PROBINGNODE
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_CHILD
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_REFOCUSNODE);
   assert(set != NULL);
   assert(tree != NULL);
   assert(var != NULL);
   assert(node->active || (infercons == NULL && inferprop == NULL));
   assert((SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_PROBINGNODE || !probingchange);

   SCIPdebugMessage("adding boundchange at node at depth %d to variable <%s>: old bounds=[%g,%g], new %s bound: %g (infer%s=<%s>, inferinfo=%d)\n",
      node->depth, SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), 
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", newbound,
      infercons != NULL ? "cons" : "prop", 
      infercons != NULL ? SCIPconsGetName(infercons) : (inferprop != NULL ? SCIPpropGetName(inferprop) : "-"), inferinfo);

   /* remember variable as inference variable, and get corresponding active variable, bound and bound type */
   infervar = var;
   inferboundtype = boundtype;
   SCIP_CALL( SCIPvarGetProbvarBound(&var, &newbound, &boundtype) );

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIPerrorMessage("cannot change bounds of multi-aggregated variale <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   oldlb = SCIPvarGetLbLocal(var);
   oldub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLT(set, oldlb, oldub));

   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      /* adjust the new lower bound */
      SCIPvarAdjustLb(var, set, &newbound);
      assert(SCIPsetIsGT(set, newbound, oldlb));
      assert(SCIPsetIsFeasLE(set, newbound, oldub));
      oldbound = oldlb;
      newbound = MIN(newbound, oldub);
   }
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);

      /* adjust the new upper bound */
      SCIPvarAdjustUb(var, set, &newbound);
      assert(SCIPsetIsLT(set, newbound, oldub));
      assert(SCIPsetIsFeasGE(set, newbound, oldlb));
      oldbound = oldub;
      newbound = MAX(newbound, oldlb);
   }
   
   SCIPdebugMessage(" -> transformed to active variable <%s>: old bounds=[%g,%g], new %s bound: %g, obj: %g\n",
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", newbound, SCIPvarGetObj(var));

   stat->nboundchgs++;

   /* if the node is the root node: change local and global bound immediately */
   if( node->depth == 0 )
   {
      assert(node->active);
      assert((SCIP_NODETYPE)node->nodetype != SCIP_NODETYPE_PROBINGNODE);
      assert(!probingchange);

      SCIPdebugMessage(" -> bound change in root node: perform global bound change\n");
      SCIP_CALL( SCIPvarChgBdGlobal(var, blkmem, set, stat, lp, branchcand, eventqueue, newbound, boundtype) );

      if( set->stage == SCIP_STAGE_SOLVING )
      {
         stat->nrootboundchgs++;
         stat->nrootboundchgsrun++;
      }

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
      SCIP_Real newpseudoobjval;
      SCIP_Real lpsolval;

      assert(!node->active || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);

      /* get the solution value of variable in last solved LP on the active path:
       *  - if the LP was solved at the current node, the LP values of the columns are valid
       *  - if the last solved LP was the one in the current lpfork, the LP value in the columns are still valid
       *  - otherwise, the LP values are invalid
       */
      if( SCIPtreeHasCurrentNodeLP(tree)
         || (tree->focuslpforklpcount == stat->lpcount && SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN) )
      {
         lpsolval = SCIPvarGetLPSol(var);
      }
      else
         lpsolval = SCIP_INVALID;

      /* remember the bound change as branching decision (infervar/infercons/inferprop are not important: use NULL) */
      SCIP_CALL( SCIPdomchgAddBoundchg(&node->domchg, blkmem, set, var, newbound, boundtype, SCIP_BOUNDCHGTYPE_BRANCHING, 
            lpsolval, NULL, NULL, NULL, 0, inferboundtype) );
      
      /* update the child's lower bound */
      if( set->misc_exactsolve )
         newpseudoobjval = SCIPlpGetModifiedProvedPseudoObjval(lp, set, var, oldbound, newbound, boundtype);
      else
         newpseudoobjval = SCIPlpGetModifiedPseudoObjval(lp, set, var, oldbound, newbound, boundtype);
      SCIPnodeUpdateLowerbound(node, stat, newpseudoobjval);
   }
   else
   {
      /* check the infered bound change on the debugging solution */
      SCIP_CALL( SCIPdebugCheckInference(blkmem, set, node, var, newbound, boundtype) ); /*lint !e506 !e774*/

      /* remember the bound change as inference (lpsolval is not important: use 0.0) */
      SCIP_CALL( SCIPdomchgAddBoundchg(&node->domchg, blkmem, set, var, newbound, boundtype,
            infercons != NULL ? SCIP_BOUNDCHGTYPE_CONSINFER : SCIP_BOUNDCHGTYPE_PROPINFER, 
            0.0, infervar, infercons, inferprop, inferinfo, inferboundtype) );
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
      SCIP_Bool cutoff;

      SCIP_CALL( SCIPboundchgApply(&node->domchg->domchgdyn.boundchgs[node->domchg->domchgdyn.nboundchgs-1],
            blkmem, set, stat, lp, branchcand, eventqueue, node->depth, node->domchg->domchgdyn.nboundchgs-1, &cutoff) );
      assert(node->domchg->domchgdyn.boundchgs[node->domchg->domchgdyn.nboundchgs-1].var == var);
      assert(!cutoff);
   }

   return SCIP_OKAY;
}

/** adds bound change to focus node, or child of focus node, or probing node;
 *  if possible, adjusts bound to integral value
 */
SCIP_RETCODE SCIPnodeAddBoundchg(
   SCIP_NODE*            node,               /**< node to add bound change to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   SCIP_Bool             probingchange       /**< is the bound change a temporary setting due to probing? */
   )
{
   SCIP_CALL( SCIPnodeAddBoundinfer(node, blkmem, set, stat, tree, lp, branchcand, eventqueue, var, newbound, boundtype,
         NULL, NULL, 0, probingchange) );

   return SCIP_OKAY;
}

/** adds hole change to focus node, or child of focus node */
SCIP_RETCODE SCIPnodeAddHolechg(
   SCIP_NODE*            node,               /**< node to add bound change to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_HOLELIST**       ptr,                /**< changed list pointer */
   SCIP_HOLELIST*        newlist,            /**< new value of list pointer */
   SCIP_HOLELIST*        oldlist             /**< old value of list pointer */
   )
{
   assert(node != NULL);
   assert((SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_FOCUSNODE
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_PROBINGNODE
      || (SCIP_NODETYPE)node->nodetype == SCIP_NODETYPE_CHILD);

   SCIPdebugMessage("adding holechange at node at depth %d: changed pointer at %p from %p to %p\n",
      node->depth, ptr, newlist, oldlist);

   stat->nholechgs++;

   SCIP_CALL( SCIPdomchgAddHolechg(&node->domchg, blkmem, set, ptr, newlist, oldlist) );

   /**@todo apply hole change on active nodes and issue event */

   return SCIP_OKAY;
}

/** if given value is larger than the node's lower bound, sets the node's lower bound to the new value */
void SCIPnodeUpdateLowerbound(
   SCIP_NODE*            node,               /**< node to update lower bound for */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
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

/** propagates implications of binary fixings at the given node triggered by the implication graph and the clique table */
SCIP_RETCODE SCIPnodePropagateImplics(
   SCIP_NODE*            node,               /**< node to propagate implications on */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   int nboundchgs;
   int i;

   assert(node != NULL);
   assert(SCIPnodeIsActive(node));
   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_FOCUSNODE
      || SCIPnodeGetType(node) == SCIP_NODETYPE_REFOCUSNODE
      || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);
   assert(cutoff != NULL);

   SCIPdebugMessage("implication graph propagation of node %p in depth %d\n", node, SCIPnodeGetDepth(node));

   *cutoff = FALSE;

   /* propagate all fixings of binary performed at this node */
   nboundchgs = SCIPdomchgGetNBoundchgs(node->domchg);
   for( i = 0; i < nboundchgs && !(*cutoff); ++i )
   {
      SCIP_BOUNDCHG* boundchg;
      SCIP_VAR* var;

      boundchg = SCIPdomchgGetBoundchg(node->domchg, i);
      var = SCIPboundchgGetVar(boundchg);
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      {
         SCIP_Bool varfixing;
         int nimpls;
         SCIP_VAR** implvars;
         SCIP_BOUNDTYPE* impltypes;
         SCIP_Real* implbounds;
         SCIP_CLIQUE** cliques;
         int ncliques;
         int j;

         varfixing = (SCIPboundchgGetBoundtype(boundchg) == SCIP_BOUNDTYPE_LOWER);
         nimpls = SCIPvarGetNImpls(var, varfixing);
         implvars = SCIPvarGetImplVars(var, varfixing);
         impltypes = SCIPvarGetImplTypes(var, varfixing);
         implbounds = SCIPvarGetImplBounds(var, varfixing);
   
         /* apply implications */
         for( j = 0; j < nimpls; ++j )
         {
            SCIP_Real lb;
            SCIP_Real ub;
                  
            /* check for infeasibility */
            lb = SCIPvarGetLbLocal(implvars[j]);
            ub = SCIPvarGetUbLocal(implvars[j]);
            if( impltypes[j] == SCIP_BOUNDTYPE_LOWER )
            {
               if( SCIPsetIsFeasGT(set, implbounds[j], ub) )
               {
                  *cutoff = TRUE;
                  return SCIP_OKAY;
               }
               if( SCIPsetIsFeasLE(set, implbounds[j], lb) )
                  continue;
            }
            else
            {
               if( SCIPsetIsFeasLT(set, implbounds[j], lb) )
               {
                  *cutoff = TRUE;
                  return SCIP_OKAY;
               }
               if( SCIPsetIsFeasGE(set, implbounds[j], ub) )
                  continue;
            }

            /* apply the implication */
            SCIP_CALL( SCIPnodeAddBoundinfer(node, blkmem, set, stat, tree, lp, branchcand, eventqueue,
                  implvars[j], implbounds[j], impltypes[j], NULL, NULL, 0, FALSE) );
         }

         /* apply cliques */
         ncliques = SCIPvarGetNCliques(var, varfixing);
         cliques = SCIPvarGetCliques(var, varfixing);
         for( j = 0; j < ncliques; ++j )
         {
            SCIP_VAR** vars;
            SCIP_Bool* values;
            int nvars;
            int k;

            nvars = SCIPcliqueGetNVars(cliques[j]);
            vars = SCIPcliqueGetVars(cliques[j]);
            values = SCIPcliqueGetValues(cliques[j]);
            for( k = 0; k < nvars; ++k )
            {
               SCIP_Real lb;
               SCIP_Real ub;
                  
               assert(SCIPvarGetType(vars[k]) == SCIP_VARTYPE_BINARY);

               if( vars[k] == var && values[k] == varfixing )
                  continue;

               /* check for infeasibility */
               lb = SCIPvarGetLbLocal(vars[k]);
               ub = SCIPvarGetUbLocal(vars[k]);
               if( values[k] == FALSE )
               {
                  if( ub < 0.5 )
                  {
                     *cutoff = TRUE;
                     return SCIP_OKAY;
                  }
                  if( lb > 0.5 )
                     continue;
               }
               else
               {
                  if( lb > 0.5 )
                  {
                     *cutoff = TRUE;
                     return SCIP_OKAY;
                  }
                  if( ub < 0.5 )
                     continue;
               }

               /* apply the clique implication */
               SCIP_CALL( SCIPnodeAddBoundinfer(node, blkmem, set, stat, tree, lp, branchcand, eventqueue,
                     vars[k], (SCIP_Real)(!values[k]), values[k] ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER,
                     NULL, NULL, 0, FALSE) );
            }
         }
      }
   }

   return SCIP_OKAY;
}




/*
 * Path Switching
 */

/** updates the LP sizes of the active path starting at the given depth */
static
void treeUpdatePathLPSize(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   startdepth          /**< depth to start counting */
   )
{
   SCIP_NODE* node;
   int ncols;
   int nrows;
   int i;

   assert(tree != NULL);
   assert(startdepth >= 0);
   assert(startdepth <= tree->pathlen);

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
      case SCIP_NODETYPE_FOCUSNODE:
         assert(i == tree->pathlen-1 || SCIPtreeProbing(tree));
         break;
      case SCIP_NODETYPE_PROBINGNODE:
         assert(SCIPtreeProbing(tree));
         break;
      case SCIP_NODETYPE_SIBLING:
         SCIPerrorMessage("sibling cannot be in the active path\n");
         SCIPABORT();
      case SCIP_NODETYPE_CHILD:
         SCIPerrorMessage("child cannot be in the active path\n");
         SCIPABORT();
      case SCIP_NODETYPE_LEAF:
         SCIPerrorMessage("leaf cannot be in the active path\n");
         SCIPABORT();
      case SCIP_NODETYPE_DEADEND:
         SCIPerrorMessage("deadend cannot be in the active path\n");
         SCIPABORT();
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
      case SCIP_NODETYPE_REFOCUSNODE:
         SCIPerrorMessage("node cannot be of type REFOCUSNODE at this point\n");
         SCIPABORT();
      default:
         SCIPerrorMessage("unknown node type\n");
         SCIPABORT();
      }
      tree->pathnlpcols[i] = ncols;
      tree->pathnlprows[i] = nrows;
   }
}

/** finds the common fork node, the new LP defining fork, and the new focus subroot, if the path is switched to
 *  the given node
 */
static
void treeFindSwitchForks(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            node,               /**< new focus node, or NULL */
   SCIP_NODE**           commonfork,         /**< pointer to store common fork node of old and new focus node */
   SCIP_NODE**           newlpfork,          /**< pointer to store the new LP defining fork node */
   SCIP_NODE**           newsubroot,         /**< pointer to store the new subroot node */
   SCIP_Bool*            cutoff              /**< pointer to store whether the given node can be cut off and no path switching
                                              *   should be performed */
   )
{
   SCIP_NODE* fork;
   SCIP_NODE* lpfork;
   SCIP_NODE* subroot;

   assert(tree != NULL);
   assert(tree->root != NULL);
   assert((tree->focusnode == NULL) == !tree->root->active);
   assert(tree->focuslpfork == NULL || tree->focusnode != NULL);
   assert(tree->focuslpfork == NULL || tree->focuslpfork->depth < tree->focusnode->depth);
   assert(tree->focussubroot == NULL || tree->focuslpfork != NULL);
   assert(tree->focussubroot == NULL || tree->focussubroot->depth <= tree->focuslpfork->depth);
   assert(tree->cutoffdepth >= 0);
   assert(tree->cutoffdepth == INT_MAX || tree->cutoffdepth < tree->pathlen);
   assert(tree->cutoffdepth == INT_MAX || tree->path[tree->cutoffdepth]->cutoff);
   assert(tree->repropdepth >= 0);
   assert(tree->repropdepth == INT_MAX || tree->repropdepth < tree->pathlen);
   assert(tree->repropdepth == INT_MAX || tree->path[tree->repropdepth]->reprop);
   assert(commonfork != NULL);
   assert(newlpfork != NULL);
   assert(newsubroot != NULL);
   assert(cutoff != NULL);

   *commonfork = NULL;
   *newlpfork = NULL;
   *newsubroot = NULL;
   *cutoff = FALSE;

   /* if the new focus node is NULL, there is no common fork node, and the new LP fork and subroot are NULL */
   if( node == NULL )
   {
      tree->cutoffdepth = INT_MAX;
      tree->repropdepth = INT_MAX;
      return;
   }

   /* check if the new node is marked to be cut off */
   if( node->cutoff )
   {
      *cutoff = TRUE;
      return;
   }

   /* if the old focus node is NULL, there is no common fork node, and we have to search the new LP fork and subroot */
   if( tree->focusnode == NULL )
   {
      assert(!tree->root->active);
      assert(tree->pathlen == 0);
      assert(tree->cutoffdepth == INT_MAX);
      assert(tree->repropdepth == INT_MAX);

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

   /* find the common fork node, the new LP defining fork, and the new focus subroot */
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
   SCIPdebugMessage("find switch forks: forkdepth=%d\n", fork->depth);

   /* if the common fork node is below the current cutoff depth, the cutoff node is an ancestor of the common fork
    * and thus an ancestor of the new focus node, s.t. the new node can also be cut off
    */
   assert((int)fork->depth != tree->cutoffdepth);
   if( (int)fork->depth > tree->cutoffdepth )
   {
#ifndef NDEBUG
      while( fork != NULL && !fork->cutoff )
         fork = fork->parent;
      assert(fork != NULL);
      assert((int)fork->depth >= tree->cutoffdepth);
#endif
      *cutoff = TRUE;
      return;
   }
   tree->cutoffdepth = INT_MAX;

   /* if not already found, continue searching the LP defining fork; it can not be deeper than the common fork */
   if( lpfork == NULL )
   {
      if( tree->focuslpfork != NULL && (int)(tree->focuslpfork->depth) > fork->depth )
      {
         /* focuslpfork is not on the same active path as the new node: we have to continue searching */
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
         /* focuslpfork is on the same active path as the new node: old and new node have the same lpfork */
         lpfork = tree->focuslpfork;
      }
      assert(lpfork == NULL || (int)(lpfork->depth) <= fork->depth);
      assert(lpfork == NULL || lpfork->active);
   }
   assert(lpfork == NULL
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_FORK
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_SUBROOT);
   SCIPdebugMessage("find switch forks: lpforkdepth=%d\n", lpfork == NULL ? -1 : (int)(lpfork->depth));

   /* if not already found, continue searching the subroot; it cannot be deeper than the LP fork and common fork */
   if( subroot == NULL )
   {
      if( tree->focussubroot != NULL && (int)(tree->focussubroot->depth) > fork->depth )
      {
         /* focussubroot is not on the same active path as the new node: we have to continue searching */
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
         subroot = tree->focussubroot;
      assert(subroot == NULL || subroot->depth <= fork->depth);
      assert(subroot == NULL || subroot->active);
   }
   assert(subroot == NULL || SCIPnodeGetType(subroot) == SCIP_NODETYPE_SUBROOT);
   assert(subroot == NULL || (lpfork != NULL && subroot->depth <= lpfork->depth));
   SCIPdebugMessage("find switch forks: subrootdepth=%d\n", subroot == NULL ? -1 : (int)(subroot->depth));

   /* if a node prior to the common fork should be repropagated, we select the node to be repropagated as common
    * fork in order to undo all bound changes up to this node, repropagate the node, and redo the bound changes
    * afterwards
    */
   if( (int)fork->depth > tree->repropdepth )
   {
      fork = tree->path[tree->repropdepth];
      assert(fork->active);
      assert(fork->reprop);
   }

   *commonfork = fork;
   *newlpfork = lpfork;
   *newsubroot = subroot;

#ifndef NDEBUG
   while( fork != NULL )
   {
      assert(fork->active);
      assert(!fork->cutoff);
      assert(fork->parent == NULL || !fork->parent->reprop);
      fork = fork->parent;
   }
#endif
   tree->repropdepth = INT_MAX;
}

/** switches the active path to the new focus node, applies domain and constraint set changes */
static
SCIP_RETCODE treeSwitchPath(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_NODE*            fork,               /**< common fork node of old and new focus node, or NULL */
   SCIP_NODE*            focusnode,          /**< new focus node, or NULL */
   SCIP_Bool*            cutoff              /**< pointer to store whether the new focus node can be cut off */
   )
{
   int focusnodedepth;  /* depth of the new focus node, or -1 if focusnode == NULL */
   int forkdepth;       /* depth of the common subroot/fork/junction node, or -1 if no common fork exists */
   int i;

   assert(tree != NULL);
   assert(fork == NULL || (fork->active && !fork->cutoff));
   assert(fork == NULL || focusnode != NULL);
   assert(focusnode == NULL || (!focusnode->active && !focusnode->cutoff));
   assert(focusnode == NULL || SCIPnodeGetType(focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(cutoff != NULL);

   *cutoff = FALSE;
   
   SCIPdebugMessage("switch path: old pathlen=%d\n", tree->pathlen);   

   /* get the nodes' depth's */
   focusnodedepth = (focusnode != NULL ? (int)focusnode->depth : -1);
   forkdepth = (fork != NULL ? (int)fork->depth : -1);
   assert(forkdepth <= focusnodedepth);
   assert(forkdepth < tree->pathlen);

   /* delay events in path switching */
   SCIP_CALL( SCIPeventqueueDelay(eventqueue) );
         
   /* undo the domain and constraint set changes of the old active path by deactivating the path's nodes */
   for( i = tree->pathlen-1; i > forkdepth; --i )
   {
      SCIP_CALL( nodeDeactivate(tree->path[i], blkmem, set, stat, tree, lp, branchcand, eventqueue) );
   }
   tree->pathlen = forkdepth+1;

   /* create the new active path */
   SCIP_CALL( treeEnsurePathMem(tree, set, focusnodedepth+1) );
   while( focusnode != fork )
   {
      assert(focusnode != NULL);
      assert(!focusnode->active);
      assert(!focusnode->cutoff);
      tree->path[focusnode->depth] = focusnode;
      focusnode = focusnode->parent;
   }

   /* propagate common fork again, if the reprop flag is set */
   if( fork != NULL && fork->reprop )
   {
      assert(tree->path[forkdepth] == fork);
      assert(fork->active);
      assert(!fork->cutoff);

      SCIP_CALL( nodeRepropagate(fork, blkmem, set, stat, prob, primal, tree, lp, branchcand, conflict, 
            eventfilter, eventqueue, cutoff) );
   }
   assert(fork != NULL || !(*cutoff));

   /* apply domain and constraint set changes of the new path by activating the path's nodes;
    * on the way, domain propagation might be applied again to the path's nodes, which can result in the cutoff of
    * the node (and its subtree)
    */
   for( i = forkdepth+1; i <= focusnodedepth && !(*cutoff); ++i )
   {
      assert(!tree->path[i]->cutoff);
      assert(tree->pathlen == i);

      /* activate the node, and apply domain propagation if the reprop flag is set */
      tree->pathlen++;
      SCIP_CALL( nodeActivate(tree->path[i], blkmem, set, stat, prob, primal, tree, lp, branchcand, conflict,
            eventfilter, eventqueue, cutoff) );
   }

   /* mark last node of path to be cut off, if a cutoff was found */
   if( *cutoff )
   {
      assert(tree->pathlen > 0);
      assert(tree->path[tree->pathlen-1]->active);
      SCIPnodeCutoff(tree->path[tree->pathlen-1], set, stat, tree);
   }

   /* count the new LP sizes of the path */
   treeUpdatePathLPSize(tree, forkdepth+1);

   /* process the delayed events */
   SCIP_CALL( SCIPeventqueueProcess(eventqueue, blkmem, set, primal, lp, branchcand, eventfilter) );

   SCIPdebugMessage("switch path: new pathlen=%d\n", tree->pathlen);   

   return SCIP_OKAY;
}

/** loads the subroot's LP data */
static
SCIP_RETCODE subrootConstructLP(
   SCIP_NODE*            subroot,            /**< subroot node to construct LP for */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_COL** cols;
   SCIP_ROW** rows;
   int ncols;
   int nrows;
   int c;
   int r;

   assert(subroot != NULL);
   assert(SCIPnodeGetType(subroot) == SCIP_NODETYPE_SUBROOT);
   assert(subroot->data.subroot != NULL);
   assert(blkmem != NULL);
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
      SCIP_CALL( SCIPlpAddCol(lp, set, cols[c], subroot->depth) );
   }
   for( r = 0; r < nrows; ++r )
   {
      SCIP_CALL( SCIPlpAddRow(lp, set, rows[r], subroot->depth) );
   }

   return SCIP_OKAY;
}
   
/** loads the fork's additional LP data */
static
SCIP_RETCODE forkAddLP(
   SCIP_NODE*            fork,               /**< fork node to construct additional LP for */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_COL** cols;
   SCIP_ROW** rows;
   int ncols;
   int nrows;
   int c;
   int r;

   assert(fork != NULL);
   assert(SCIPnodeGetType(fork) == SCIP_NODETYPE_FORK);
   assert(fork->data.fork != NULL);
   assert(blkmem != NULL);
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
      SCIP_CALL( SCIPlpAddCol(lp, set, cols[c], fork->depth) );
   }
   for( r = 0; r < nrows; ++r )
   {
      SCIP_CALL( SCIPlpAddRow(lp, set, rows[r], fork->depth) );
   }

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** checks validity of active path */
static
void treeCheckPath(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   SCIP_NODE* node;
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
      case SCIP_NODETYPE_FOCUSNODE:
      case SCIP_NODETYPE_REFOCUSNODE:
         assert(d == tree->pathlen-1 || SCIPtreeProbing(tree));
         break;
      case SCIP_NODETYPE_PROBINGNODE:
         assert(SCIPtreeProbing(tree));
         break;
      default:
         SCIPerrorMessage("node at depth %d on active path has to be of type FORK, SUBROOT, FOCUSNODE, REFOCUSNODE, or PROBINGNODE, but is %d\n",
            d, SCIPnodeGetType(node));
         SCIPABORT();
      }  /*lint !e788*/
      assert(tree->pathnlpcols[d] == ncols);
      assert(tree->pathnlprows[d] == nrows);
   }
}
#else
#define treeCheckPath(tree) /**/
#endif

/** constructs the LP and loads LP state for fork/subroot of the focus node */
SCIP_RETCODE SCIPtreeLoadLP(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool*            initroot            /**< pointer to store whether the root LP relaxation has to be initialized */
   )
{
   SCIP_NODE* lpfork;
   SCIP_NODE* pathnode;
   int lpforkdepth;
   int d;

   assert(tree != NULL);
   assert(!tree->focuslpconstructed);
   assert(tree->path != NULL);
   assert(tree->pathlen > 0);
   assert(tree->focusnode != NULL);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(SCIPnodeGetDepth(tree->focusnode) == tree->pathlen-1);
   assert(tree->focusnode == tree->path[tree->pathlen-1]);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(initroot != NULL);

   SCIPdebugMessage("load LP for current fork node %p at depth %d\n", 
      tree->focuslpfork, tree->focuslpfork == NULL ? -1 : (int)(tree->focuslpfork->depth));
   SCIPdebugMessage("-> old LP has %d cols and %d rows\n", SCIPlpGetNCols(lp), SCIPlpGetNRows(lp));
   SCIPdebugMessage("-> correct LP has %d cols and %d rows\n", 
      tree->correctlpdepth >= 0 ? tree->pathnlpcols[tree->correctlpdepth] : 0,
      tree->correctlpdepth >= 0 ? tree->pathnlprows[tree->correctlpdepth] : 0);
   SCIPdebugMessage("-> old correctlpdepth: %d\n", tree->correctlpdepth);

   treeCheckPath(tree);

   lpfork = tree->focuslpfork;

   /* find out the lpfork's depth (or -1, if lpfork is NULL) */
   if( lpfork == NULL )
   {
      assert(tree->correctlpdepth == -1 || tree->pathnlpcols[tree->correctlpdepth] == 0);
      assert(tree->correctlpdepth == -1 || tree->pathnlprows[tree->correctlpdepth] == 0);
      assert(tree->focussubroot == NULL);
      lpforkdepth = -1;
   }
   else
   {
      assert(SCIPnodeGetType(lpfork) == SCIP_NODETYPE_FORK || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_SUBROOT);
      assert(lpfork->active);
      assert(tree->path[lpfork->depth] == lpfork);
      lpforkdepth = lpfork->depth;
   }
   assert(lpforkdepth < tree->pathlen-1); /* lpfork must not be the last (the focus) node of the active path */

   /* find out, if we are in the same subtree */
   if( tree->correctlpdepth >= 0 )
   {
      /* same subtree: shrink LP to the deepest node with correct LP */
      assert(lpforkdepth == -1 || tree->pathnlpcols[tree->correctlpdepth] <= tree->pathnlpcols[lpforkdepth]);
      assert(lpforkdepth == -1 || tree->pathnlprows[tree->correctlpdepth] <= tree->pathnlprows[lpforkdepth]);
      assert(lpforkdepth >= 0 || tree->pathnlpcols[tree->correctlpdepth] == 0);
      assert(lpforkdepth >= 0 || tree->pathnlprows[tree->correctlpdepth] == 0);
      SCIP_CALL( SCIPlpShrinkCols(lp, tree->pathnlpcols[tree->correctlpdepth]) );
      SCIP_CALL( SCIPlpShrinkRows(lp, blkmem, set, tree->pathnlprows[tree->correctlpdepth]) );
   }
   else
   {
      /* other subtree: fill LP with the subroot LP data */
      SCIP_CALL( SCIPlpClear(lp, blkmem, set) );
      if( tree->focussubroot != NULL )
      {
         SCIP_CALL( subrootConstructLP(tree->focussubroot, blkmem, set, lp) );
         tree->correctlpdepth = tree->focussubroot->depth; 
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
         SCIP_CALL( forkAddLP(pathnode, blkmem, set, lp) );
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
      if( tree->focuslpforklpcount != stat->lpcount )
      {
         if( SCIPnodeGetType(lpfork) == SCIP_NODETYPE_FORK )
         {
            assert(lpfork->data.fork != NULL);
            SCIP_CALL( SCIPlpSetState(lp, blkmem, set, lpfork->data.fork->lpistate) );
         }
         else
         {
            assert(SCIPnodeGetType(lpfork) == SCIP_NODETYPE_SUBROOT);
            assert(lpfork->data.subroot != NULL);
            SCIP_CALL( SCIPlpSetState(lp, blkmem, set, lpfork->data.subroot->lpistate) );
         }
         assert(lp->primalfeasible);
         assert(lp->dualfeasible);
      }
      else
      {
         lp->primalfeasible = TRUE;
         lp->dualfeasible = TRUE;
      }

      /* check the path from LP fork to focus node for domain changes (destroying primal feasibility of LP basis) */
      for( d = lpforkdepth; d < (int)(tree->focusnode->depth) && lp->primalfeasible; ++d )
      {
         assert(d < tree->pathlen);
         lp->primalfeasible = lp->primalfeasible
            && (tree->path[d]->domchg == NULL || tree->path[d]->domchg->domchgbound.nboundchgs == 0);
      }
   }

   /* mark the LP's size, such that we know which rows and columns were added in the new node */
   SCIPlpMarkSize(lp);

   SCIPdebugMessage("-> new correctlpdepth: %d\n", tree->correctlpdepth);
   SCIPdebugMessage("-> new LP has %d cols and %d rows, primalfeasible=%d, dualfeasible=%d\n", 
      SCIPlpGetNCols(lp), SCIPlpGetNRows(lp), lp->primalfeasible, lp->dualfeasible);

   /* if the correct LP depth is still -1, the root LP relaxation has to be initialized */
   *initroot = (tree->correctlpdepth == -1);

   /* mark the LP of the focus node constructed */
   tree->focuslpconstructed = TRUE;

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
SCIP_RETCODE nodeToLeaf(
   SCIP_NODE**           node,               /**< pointer to child or sibling node to convert */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_NODE*            lpfork,             /**< LP fork of the node */
   SCIP_Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   )
{
   assert(SCIPnodeGetType(*node) == SCIP_NODETYPE_SIBLING || SCIPnodeGetType(*node) == SCIP_NODETYPE_CHILD);
   assert(stat != NULL);
   assert(lpfork == NULL || lpfork->depth < (*node)->depth);
   assert(lpfork == NULL || lpfork->active || SCIPsetIsGE(set, (*node)->lowerbound, cutoffbound));
   assert(lpfork == NULL
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_FORK
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_SUBROOT);

   /* convert node into leaf */
   SCIPdebugMessage("convert node %p at depth %d to leaf with lpfork %p at depth %d\n",
      *node, (*node)->depth, lpfork, lpfork == NULL ? -1 : (int)(lpfork->depth));
   (*node)->nodetype = SCIP_NODETYPE_LEAF; /*lint !e641*/
   (*node)->data.leaf.lpfork = lpfork;

#ifndef NDEBUG
   /* check, if the LP fork is the first node with LP information on the path back to the root */
   if( cutoffbound != SCIP_REAL_MIN ) /* if the node was cut off in SCIPnodeFocus(), the lpfork is invalid */
   {
      SCIP_NODE* pathnode;
      pathnode = (*node)->parent;
      while( pathnode != NULL && pathnode != lpfork )
      {
         assert(SCIPnodeGetType(pathnode) == SCIP_NODETYPE_JUNCTION);
         pathnode = pathnode->parent;
      }
      assert(pathnode == lpfork);
   }
#endif

   /* if node is good enough to keep, put it on the node queue */
   if( SCIPsetIsLT(set, (*node)->lowerbound, cutoffbound) )
   {
      /* insert leaf in node queue */
      SCIP_CALL( SCIPnodepqInsert(tree->leaves, set, *node) );
      
      /* make the domain change data static to save memory */
      SCIP_CALL( SCIPdomchgMakeStatic(&(*node)->domchg, blkmem, set) );

      /* node is now member of the node queue: delete the pointer to forbid further access */
      *node = NULL;
   }
   else
   {
      /* delete node due to bound cut off */
      SCIPvbcCutoffNode(stat->vbc, stat, *node);
      SCIP_CALL( SCIPnodeFree(node, blkmem, set, tree, lp) );
   }
   assert(*node == NULL);

   return SCIP_OKAY;
}

/** converts the focus node into a deadend node */
static
SCIP_RETCODE focusnodeToDeadend(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(blkmem != NULL);
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode != NULL);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(tree->nchildren == 0);

   SCIPdebugMessage("focusnode to deadend %p at depth %d\n", tree->focusnode, tree->focusnode->depth);

   tree->focusnode->nodetype = SCIP_NODETYPE_DEADEND; /*lint !e641*/

   /* release LPI state */
   if( tree->focuslpfork != NULL )
   {
      SCIP_CALL( SCIPnodeReleaseLPIState(tree->focuslpfork, blkmem, lp) );
   }

   return SCIP_OKAY;
}

/** converts the focus node into a junction node */
static
SCIP_RETCODE focusnodeToJunction(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode != NULL);
   assert(tree->focusnode->active); /* otherwise, no children could be created at the focus node */
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);

   SCIPdebugMessage("focusnode %p to junction at depth %d\n", tree->focusnode, tree->focusnode->depth);

   tree->focusnode->nodetype = SCIP_NODETYPE_JUNCTION; /*lint !e641*/

   SCIP_CALL( junctionInit(&tree->focusnode->data.junction, tree) );

   /* release LPI state */
   if( tree->focuslpfork != NULL )
   {
      SCIP_CALL( SCIPnodeReleaseLPIState(tree->focuslpfork, blkmem, lp) );
   }

   /* make the domain change data static to save memory */
   SCIP_CALL( SCIPdomchgMakeStatic(&tree->focusnode->domchg, blkmem, set) );

   return SCIP_OKAY;
}

/** converts the focus node into a fork node */
static
SCIP_RETCODE focusnodeToFork(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_FORK* fork;
   SCIP_Bool lperror;

   assert(blkmem != NULL);
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode != NULL);
   assert(tree->focusnode->active); /* otherwise, no children could be created at the focus node */
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(tree->nchildren > 0);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved || lp->resolvelperror);

   SCIPdebugMessage("focusnode %p to fork at depth %d\n", tree->focusnode, tree->focusnode->depth);

   /* usually, the LP should be solved to optimality; otherwise, numerical troubles occured,
    * and we have to forget about the LP and transform the node into a junction (see below)
    */
   lperror = FALSE;
   if( !lp->resolvelperror && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      /* clean up newly created part of LP to keep only necessary columns and rows */
      SCIP_CALL( SCIPlpCleanupNew(lp, blkmem, set, stat, (tree->focusnode->depth == 0)) );

      /* resolve LP after cleaning up */
      if( !lp->solved || !lp->flushed )
      {
         SCIPdebugMessage("resolving LP after cleanup\n");
         SCIP_CALL( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, FALSE, TRUE, &lperror) );
      }
   }
   assert(lp->flushed);
   assert(lp->solved || lperror || lp->resolvelperror);

   /* There are two reasons, that the (reduced) SCIP_LP is not solved to optimality:
    *  - The primal heuristics (called after the current node's LP was solved) found a new 
    *    solution, that is better than the current node's lower bound.
    *    (But in this case, all children should be cut off and the node should be converted
    *    into a deadend instead of a fork.)
    *  - Something numerically weird happened after cleaning up or after resolving a diving or probing LP.
    * The only thing we can do, is to completely forget about the LP and treat the node as
    * if it was only a pseudo-solution node. Therefore we have to remove all additional
    * columns and rows from the LP and convert the node into a junction.
    * However, the node's lower bound is kept, thus automatically throwing away nodes that
    * were cut off due to a primal solution.
    */
   if( lperror || lp->resolvelperror || SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %"SCIP_LONGINT_FORMAT") numerical troubles: LP %d not optimal -- convert node into junction instead of fork\n", 
         stat->nnodes, stat->nlps);

      /* remove all additions to the LP at this node */
      SCIP_CALL( SCIPlpShrinkCols(lp, SCIPlpGetNCols(lp) - SCIPlpGetNNewcols(lp)) );
      SCIP_CALL( SCIPlpShrinkRows(lp, blkmem, set, SCIPlpGetNRows(lp) - SCIPlpGetNNewrows(lp)) );
      
      /* convert node into a junction */
      SCIP_CALL( focusnodeToJunction(blkmem, set, tree, lp) );
      
      return SCIP_OKAY;
   }
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);

   /* create fork data */
   SCIP_CALL( forkCreate(&fork, blkmem, tree, lp) );
   
   tree->focusnode->nodetype = SCIP_NODETYPE_FORK; /*lint !e641*/
   tree->focusnode->data.fork = fork;

   /* release LPI state */
   if( tree->focuslpfork != NULL )
   {
      SCIP_CALL( SCIPnodeReleaseLPIState(tree->focuslpfork, blkmem, lp) );
   }

   /* make the domain change data static to save memory */
   SCIP_CALL( SCIPdomchgMakeStatic(&tree->focusnode->domchg, blkmem, set) );

   return SCIP_OKAY;
}

/** converts the focus node into a subroot node */
static
SCIP_RETCODE focusnodeToSubroot(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_SUBROOT* subroot;
   SCIP_Bool lperror;

   assert(blkmem != NULL);
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(tree->focusnode != NULL);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
   assert(tree->focusnode->active); /* otherwise, no children could be created at the focus node */
   assert(tree->nchildren > 0);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   SCIPdebugMessage("focusnode %p to subroot at depth %d\n", tree->focusnode, tree->focusnode->depth);

   /* usually, the LP should be solved to optimality; otherwise, numerical troubles occured,
    * and we have to forget about the LP and transform the node into a junction (see below)
    */
   lperror = FALSE;
   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      /* clean up whole LP to keep only necessary columns and rows */
#if 0
      if( tree->focusnode->depth == 0 )
      {
         SCIP_CALL( SCIPlpCleanupAll(lp, blkmem, set, stat, (tree->focusnode->depth == 0)) );
      }
      else
#endif
      {
         SCIP_CALL( SCIPlpRemoveAllObsoletes(lp, blkmem, set, stat) );
      }

      /* resolve LP after cleaning up */
      if( !lp->solved || !lp->flushed )
      {
         SCIPdebugMessage("resolving LP after cleanup\n");
         SCIP_CALL( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, FALSE, TRUE, &lperror) );
      }
   }
   assert(lp->flushed);
   assert(lp->solved || lperror);

   /* There are two reasons, that the (reduced) SCIP_LP is not solved to optimality:
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
   if( lperror || SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %"SCIP_LONGINT_FORMAT") numerical troubles: LP %d not optimal -- convert node into junction instead of subroot\n", 
         stat->nnodes, stat->nlps);

      /* remove all additions to the LP at this node */
      SCIP_CALL( SCIPlpShrinkCols(lp, SCIPlpGetNCols(lp) - SCIPlpGetNNewcols(lp)) );
      SCIP_CALL( SCIPlpShrinkRows(lp, blkmem, set, SCIPlpGetNRows(lp) - SCIPlpGetNNewrows(lp)) );
      
      /* convert node into a junction */
      SCIP_CALL( focusnodeToJunction(blkmem, set, tree, lp) );
      
      return SCIP_OKAY;
   }
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);

   /* create subroot data */
   SCIP_CALL( subrootCreate(&subroot, blkmem, tree, lp) );

   tree->focusnode->nodetype = SCIP_NODETYPE_SUBROOT; /*lint !e641*/
   tree->focusnode->data.subroot = subroot;

   /* update the LP column and row counter for the converted node */
   treeUpdatePathLPSize(tree, tree->focusnode->depth);

   /* release LPI state */
   if( tree->focuslpfork != NULL )
   {
      SCIP_CALL( SCIPnodeReleaseLPIState(tree->focuslpfork, blkmem, lp) );
   }

   /* make the domain change data static to save memory */
   SCIP_CALL( SCIPdomchgMakeStatic(&tree->focusnode->domchg, blkmem, set) );

   return SCIP_OKAY;
}

/** puts all nodes in the array on the node queue and makes them LEAFs */
static
SCIP_RETCODE treeNodesToQueue(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_NODE**           nodes,              /**< array of nodes to put on the queue */
   int*                  nnodes,             /**< pointer to number of nodes in the array */
   SCIP_NODE*            lpfork,             /**< LP fork of the nodes */
   SCIP_Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
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
      SCIP_CALL( nodeToLeaf(&nodes[i], blkmem, set, stat, tree, lp, lpfork, cutoffbound) );
      assert(nodes[i] == NULL);
   }
   *nnodes = 0;

   return SCIP_OKAY;
}

/** converts children into siblings, clears children array */
static
void treeChildrenToSiblings(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   SCIP_NODE** tmpnodes;
   SCIP_Real* tmpprios;
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

/** installs a child, a sibling, or a leaf node as the new focus node */
SCIP_RETCODE SCIPnodeFocus(
   SCIP_NODE**           node,               /**< pointer to node to focus (or NULL to remove focus); the node
                                              *   is freed, if it was cut off due to a cut off subtree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool*            cutoff              /**< pointer to store whether the given node can be cut off */
   )
{
   SCIP_NODE* oldfocusnode;
   SCIP_NODE* fork;
   SCIP_NODE* lpfork;
   SCIP_NODE* subroot;
   SCIP_NODE* childrenlpfork;
   int oldcutoffdepth;

   assert(node != NULL);
   assert(*node == NULL
      || SCIPnodeGetType(*node) == SCIP_NODETYPE_SIBLING
      || SCIPnodeGetType(*node) == SCIP_NODETYPE_CHILD
      || SCIPnodeGetType(*node) == SCIP_NODETYPE_LEAF);
   assert(*node == NULL || !(*node)->active);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(!SCIPtreeProbing(tree));
   assert(lp != NULL);
   assert(cutoff != NULL);

   SCIPdebugMessage("focussing node %p of type %d in depth %d\n",
      *node, *node != NULL ? (int)SCIPnodeGetType(*node) : 0, *node != NULL ? SCIPnodeGetDepth(*node) : 0);

   /* remember old cutoff depth in order to know, whether the children and siblings can be deleted */
   oldcutoffdepth = tree->cutoffdepth;

   /* find the common fork node, the new LP defining fork, and the new focus subroot,
    * thereby checking, if the new node can be cut off
    */
   treeFindSwitchForks(tree, *node, &fork, &lpfork, &subroot, cutoff);
   SCIPdebugMessage("focus node: focusnodedepth=%d, forkdepth=%d, lpforkdepth=%d, subrootdepth=%d, cutoff=%d\n",
      *node != NULL ? (*node)->depth : -1, fork != NULL ? fork->depth : -1, lpfork != NULL ? lpfork->depth : -1,
      subroot != NULL ? subroot->depth : -1, *cutoff);

   /* free the new node, if it is located in a cut off subtree */
   if( *cutoff )
   {
      assert(*node != NULL);
      assert(tree->cutoffdepth == oldcutoffdepth);
      if( SCIPnodeGetType(*node) == SCIP_NODETYPE_LEAF )
      {
         SCIP_CALL( SCIPnodepqRemove(tree->leaves, set, *node) );
      }
      SCIP_CALL( SCIPnodeFree(node, blkmem, set, tree, lp) );

      return SCIP_OKAY;
   }

   assert(tree->cutoffdepth == INT_MAX);
   assert(fork == NULL || fork->active);
   assert(lpfork == NULL || fork != NULL);
   assert(subroot == NULL || lpfork != NULL);

   /* remember the depth of the common fork node for LP updates */
   SCIPdebugMessage("focus node: old correctlpdepth=%d\n", tree->correctlpdepth);
   if( subroot == tree->focussubroot && fork != NULL && lpfork != NULL )
   {
      /* we are in the same subtree with valid LP fork: the LP is correct at most upto the common fork depth */
      assert(subroot == NULL || subroot->active);
      tree->correctlpdepth = MIN(tree->correctlpdepth, (int)fork->depth);
   }
   else
   {
      /* we are in a different subtree, or no valid LP fork exists: the LP is completely incorrect */
      assert(subroot == NULL || !subroot->active
         || (tree->focussubroot != NULL && (int)(tree->focussubroot->depth) > subroot->depth));
      tree->correctlpdepth = -1;
   }

   /* if the LP fork changed, the lpcount information for the new LP fork is unknown */
   if( lpfork != tree->focuslpfork )
      tree->focuslpforklpcount = -1;

   /* if the old focus node was cut off, we can delete its children;
    * if the old focus node's parent was cut off, we can also delete the focus node's siblings
    */
   if( tree->focusnode != NULL && oldcutoffdepth <= (int)tree->focusnode->depth )
   {
      SCIPdebugMessage("path to old focus node of depth %d was cut off at depth %d\n", 
         tree->focusnode->depth, oldcutoffdepth);

      /* delete the focus node's children by converting them to leaves with a cutoffbound of SCIP_REAL_MIN;
       * we cannot delete them directly, because in SCIPnodeFree(), the children array is changed, which is the
       * same array we would have to iterate over here;
       * the children don't have an LP fork, because the old focus node is not yet converted into a fork or subroot
       */
      SCIPdebugMessage(" -> deleting the %d children of the old focus node\n", tree->nchildren);
      SCIP_CALL( treeNodesToQueue(tree, blkmem, set, stat, lp, tree->children, &tree->nchildren, NULL, SCIP_REAL_MIN) );
      assert(tree->nchildren == 0);

      if( oldcutoffdepth < (int)tree->focusnode->depth )
      {
         /* delete the focus node's siblings by converting them to leaves with a cutoffbound of SCIP_REAL_MIN;
          * we cannot delete them directly, because in SCIPnodeFree(), the siblings array is changed, which is the
          * same array we would have to iterate over here;
          * the siblings have the same LP fork as the old focus node
          */
         SCIPdebugMessage(" -> deleting the %d siblings of the old focus node\n", tree->nsiblings);
         SCIP_CALL( treeNodesToQueue(tree, blkmem, set, stat, lp, tree->siblings, &tree->nsiblings, tree->focuslpfork,
               SCIP_REAL_MIN) );
         assert(tree->nsiblings == 0);
      }
   }

   /* convert the old focus node into a fork or subroot node, if it has children;
    * otherwise, convert it into a deadend, which will be freed later in treeSwitchPath()
    */
   childrenlpfork = tree->focuslpfork;
   if( tree->nchildren > 0 )
   {
      assert(tree->focusnode != NULL);
      assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE);
      assert(oldcutoffdepth == INT_MAX);

      if( tree->focusnodehaslp )
      {
         /**@todo decide: old focus node becomes fork or subroot */
         if( FALSE && tree->focusnode->depth > 0 && tree->focusnode->depth % 25 == 0 ) /* ????????????? */
         {
            /* convert old focus node into a subroot node */
            SCIP_CALL( focusnodeToSubroot(blkmem, set, stat, prob, tree, lp) );
            if( *node != NULL && SCIPnodeGetType(*node) == SCIP_NODETYPE_CHILD
               && SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_SUBROOT )
               subroot = tree->focusnode;
         }
         else
         {
            /* convert old focus node into a fork node */
            SCIP_CALL( focusnodeToFork(blkmem, set, stat, prob, tree, lp) );
         }

         /* update the path's LP size */
         tree->pathnlpcols[tree->focusnode->depth] = lp->ncols;
         tree->pathnlprows[tree->focusnode->depth] = lp->nrows;

         /* check, if the conversion into a subroot or fork was successful */
         if( SCIPnodeGetType(tree->focusnode) != SCIP_NODETYPE_JUNCTION )
         {
            assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_SUBROOT
               || SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FORK);
            childrenlpfork = tree->focusnode;

            /* if a child of the old focus node was selected as new focus node, the old node becomes the new focus
             * LP fork
             */
            if( *node != NULL && SCIPnodeGetType(*node) == SCIP_NODETYPE_CHILD )
            {
               lpfork = tree->focusnode;
               tree->correctlpdepth = tree->focusnode->depth;
               tree->focuslpforklpcount = stat->lpcount;
            }
         }
      }
      else
      {
         /* convert old focus node into junction */
         SCIP_CALL( focusnodeToJunction(blkmem, set, tree, lp) );
      }
   }
   else if( tree->focusnode != NULL )
   {
      /* convert old focus node into deadend */
      SCIP_CALL( focusnodeToDeadend(blkmem, tree, lp) );
   }
   assert(subroot == NULL || SCIPnodeGetType(subroot) == SCIP_NODETYPE_SUBROOT);
   assert(lpfork == NULL
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_SUBROOT
      || SCIPnodeGetType(lpfork) == SCIP_NODETYPE_FORK);
   assert(childrenlpfork == NULL
      || SCIPnodeGetType(childrenlpfork) == SCIP_NODETYPE_SUBROOT
      || SCIPnodeGetType(childrenlpfork) == SCIP_NODETYPE_FORK);
   SCIPdebugMessage("focus node: new correctlpdepth=%d\n", tree->correctlpdepth);
   
   /* set up the new lists of siblings and children */
   oldfocusnode = tree->focusnode;
   if( *node == NULL )
   {
      /* move siblings to the queue, make them LEAFs */
      SCIP_CALL( treeNodesToQueue(tree, blkmem, set, stat, lp, tree->siblings, &tree->nsiblings, tree->focuslpfork,
            primal->cutoffbound) );

      /* move children to the queue, make them LEAFs */
      SCIP_CALL( treeNodesToQueue(tree, blkmem, set, stat, lp, tree->children, &tree->nchildren, childrenlpfork, 
            primal->cutoffbound) );
   }
   else
   {
      SCIP_NODE* bestleaf;

      switch( SCIPnodeGetType(*node) )
      {  
      case SCIP_NODETYPE_SIBLING:
         /* reset plunging depth, if the selected node is better than all leaves */
         bestleaf = SCIPtreeGetBestLeaf(tree);
         if( bestleaf == NULL || SCIPnodepqCompare(tree->leaves, set, *node, bestleaf) <= 0 )
            stat->plungedepth = 0;

         /* move children to the queue, make them LEAFs */
         SCIP_CALL( treeNodesToQueue(tree, blkmem, set, stat, lp, tree->children, &tree->nchildren, childrenlpfork,
               primal->cutoffbound) );

         /* remove selected sibling from the siblings array */
         treeRemoveSibling(tree, *node);

         SCIPdebugMessage("selected sibling node, lowerbound=%g, plungedepth=%d\n", (*node)->lowerbound, stat->plungedepth);
         break;
         
      case SCIP_NODETYPE_CHILD:
         /* reset plunging depth, if the selected node is better than all leaves; otherwise, increase plunging depth */
         bestleaf = SCIPtreeGetBestLeaf(tree);
         if( bestleaf == NULL || SCIPnodepqCompare(tree->leaves, set, *node, bestleaf) <= 0 )
            stat->plungedepth = 0;
         else
            stat->plungedepth++;

         /* move siblings to the queue, make them LEAFs */
         SCIP_CALL( treeNodesToQueue(tree, blkmem, set, stat, lp, tree->siblings, &tree->nsiblings, tree->focuslpfork,
               primal->cutoffbound) );

         /* remove selected child from the children array */      
         treeRemoveChild(tree, *node);
         
         /* move remaining children to the siblings array, make them SIBLINGs */
         treeChildrenToSiblings(tree);
         
         SCIPdebugMessage("selected child node, lowerbound=%g, plungedepth=%d\n", (*node)->lowerbound, stat->plungedepth);
         break;
         
      case SCIP_NODETYPE_LEAF:
         /* move siblings to the queue, make them LEAFs */
         SCIP_CALL( treeNodesToQueue(tree, blkmem, set, stat, lp, tree->siblings, &tree->nsiblings, tree->focuslpfork,
               primal->cutoffbound) );
         
         /* move children to the queue, make them LEAFs */
         SCIP_CALL( treeNodesToQueue(tree, blkmem, set, stat, lp, tree->children, &tree->nchildren, childrenlpfork,
               primal->cutoffbound) );

         /* remove node from the queue */
         SCIP_CALL( SCIPnodepqRemove(tree->leaves, set, *node) );

         stat->plungedepth = 0;
         if( SCIPnodeGetDepth(*node) > 0 )
            stat->nbacktracks++;
         SCIPdebugMessage("selected leaf node, lowerbound=%g, plungedepth=%d\n", (*node)->lowerbound, stat->plungedepth);
         break;

      default:
         SCIPerrorMessage("selected node is neither sibling, child, nor leaf (nodetype=%d)\n", SCIPnodeGetType(*node));
         return SCIP_INVALIDDATA;
      }  /*lint !e788*/

      /* convert node into the focus node */
      (*node)->nodetype = SCIP_NODETYPE_FOCUSNODE; /*lint !e641*/
   }
   assert(tree->nchildren == 0);
   
   /* set new focus node, LP fork, and subroot */
   assert(subroot == NULL || (lpfork != NULL && subroot->depth <= lpfork->depth));
   assert(lpfork == NULL || (*node != NULL && lpfork->depth < (*node)->depth));
   tree->focusnode = *node;
   tree->focuslpfork = lpfork;
   tree->focussubroot = subroot;
   tree->focuslpconstructed = FALSE;
   lp->resolvelperror = FALSE;

   /* track the path from the old focus node to the new node, and perform domain and constraint set changes */
   SCIP_CALL( treeSwitchPath(tree, blkmem, set, stat, prob, primal, lp, branchcand, conflict, eventfilter, eventqueue, 
         fork, *node, cutoff) );
   assert(tree->pathlen >= 0);
   assert(*node != NULL || tree->pathlen == 0);
   assert(*node == NULL || tree->pathlen-1 <= (int)(*node)->depth);

   /* if the old focus node is a dead end (has no children), delete it */
   if( oldfocusnode != NULL && SCIPnodeGetType(oldfocusnode) == SCIP_NODETYPE_DEADEND )
   {
      SCIP_CALL( SCIPnodeFree(&oldfocusnode, blkmem, set, tree, lp) );
   }
   assert(*cutoff || SCIPtreeIsPathComplete(tree));

   return SCIP_OKAY;
}   




/*
 * Tree methods
 */

/** creates an initialized tree data structure */
SCIP_RETCODE SCIPtreeCreate(
   SCIP_TREE**           tree,               /**< pointer to tree data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODESEL*         nodesel             /**< node selector to use for sorting leaves in the priority queue */
   )
{
   assert(tree != NULL);

   SCIP_ALLOC( BMSallocMemory(tree) );

   (*tree)->root = NULL;

   SCIP_CALL( SCIPnodepqCreate(&(*tree)->leaves, set, nodesel) );

   (*tree)->path = NULL;
   (*tree)->focusnode = NULL;
   (*tree)->focuslpfork = NULL;
   (*tree)->focussubroot = NULL;
   (*tree)->children = NULL;
   (*tree)->siblings = NULL;
   (*tree)->probingroot = NULL;
   (*tree)->childrenprio = NULL;
   (*tree)->siblingsprio = NULL;
   (*tree)->pathnlpcols = NULL;
   (*tree)->pathnlprows = NULL;
   (*tree)->probinglpistate = NULL;
   (*tree)->focuslpforklpcount = -1;
   (*tree)->childrensize = 0;
   (*tree)->nchildren = 0;
   (*tree)->siblingssize = 0;
   (*tree)->nsiblings = 0;
   (*tree)->pathlen = 0;
   (*tree)->pathsize = 0;
   (*tree)->correctlpdepth = -1;
   (*tree)->cutoffdepth = INT_MAX;
   (*tree)->repropdepth = INT_MAX;
   (*tree)->repropsubtreecount = 0;
   (*tree)->focusnodehaslp = FALSE;
   (*tree)->probingnodehaslp = FALSE;
   (*tree)->focuslpconstructed = FALSE;
   (*tree)->cutoffdelayed = FALSE;
   (*tree)->probinglpwasflushed = FALSE;
   (*tree)->probinglpwassolved = FALSE;

   return SCIP_OKAY;
}

/** frees tree data structure */
SCIP_RETCODE SCIPtreeFree(
   SCIP_TREE**           tree,               /**< pointer to tree data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);
   assert(*tree != NULL);
   assert((*tree)->nchildren == 0);
   assert((*tree)->nsiblings == 0);
   assert((*tree)->focusnode == NULL);
   assert(!SCIPtreeProbing(*tree));

   SCIPdebugMessage("free tree\n");

   /* free node queue */
   SCIP_CALL( SCIPnodepqFree(&(*tree)->leaves, blkmem, set, *tree, lp) );
   
   /* free pointer arrays */
   BMSfreeMemoryArrayNull(&(*tree)->path);
   BMSfreeMemoryArrayNull(&(*tree)->children);
   BMSfreeMemoryArrayNull(&(*tree)->siblings);
   BMSfreeMemoryArrayNull(&(*tree)->childrenprio);
   BMSfreeMemoryArrayNull(&(*tree)->siblingsprio);
   BMSfreeMemoryArrayNull(&(*tree)->pathnlpcols);
   BMSfreeMemoryArrayNull(&(*tree)->pathnlprows);

   BMSfreeMemory(tree);

   return SCIP_OKAY;
}

/** clears and resets tree data structure and deletes all nodes */
SCIP_RETCODE SCIPtreeClear(
   SCIP_TREE*            tree,               /**< tree data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(tree->nsiblings == 0);
   assert(tree->focusnode == NULL);
   assert(!SCIPtreeProbing(tree));

   SCIPdebugMessage("clearing tree\n");

   /* clear node queue */
   SCIP_CALL( SCIPnodepqClear(tree->leaves, blkmem, set, tree, lp) );
   assert(tree->root == NULL);
   
   /* mark working arrays to be empty and reset data */
   tree->focuslpforklpcount = -1;
   tree->nchildren = 0;
   tree->nsiblings = 0;
   tree->pathlen = 0;
   tree->correctlpdepth = -1;
   tree->cutoffdepth = INT_MAX;
   tree->repropdepth = INT_MAX;
   tree->repropsubtreecount = 0;
   tree->focusnodehaslp = FALSE;
   tree->probingnodehaslp = FALSE;
   tree->cutoffdelayed = FALSE;
   tree->probinglpwasflushed = FALSE;
   tree->probinglpwassolved = FALSE;

   return SCIP_OKAY;
}

/** creates the root node of the tree and puts it into the leaves queue */
SCIP_RETCODE SCIPtreeCreateRoot(
   SCIP_TREE*            tree,               /**< tree data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(tree->nsiblings == 0);
   assert(tree->root == NULL);
   assert(tree->focusnode == NULL);
   assert(!SCIPtreeProbing(tree));

   /* create root node */
   SCIP_CALL( SCIPnodeCreateChild(&tree->root, blkmem, set, stat, tree, 0.0) );
   assert(tree->nchildren == 1);

#ifndef NDEBUG
   /* check, if the sizes in the data structures match the maximal numbers defined here */
   tree->root->depth = MAXDEPTH;
   tree->root->repropsubtreemark = MAXREPROPMARK;
   assert(tree->root->depth == MAXDEPTH);
   assert(tree->root->repropsubtreemark == MAXREPROPMARK);
   tree->root->depth++;             /* this should produce an overflow and reset the value to 0 */
   tree->root->repropsubtreemark++; /* this should produce an overflow and reset the value to 0 */
   assert(tree->root->depth == 0);
   assert((SCIP_NODETYPE)tree->root->nodetype == SCIP_NODETYPE_CHILD);
   assert(!tree->root->active);
   assert(!tree->root->cutoff);
   assert(!tree->root->reprop);
   assert(tree->root->repropsubtreemark == 0);
#endif

   /* move root to the queue, convert it to LEAF */
   SCIP_CALL( treeNodesToQueue(tree, blkmem, set, stat, lp, tree->children, &tree->nchildren, NULL, 
         SCIPsetInfinity(set)) );

   return SCIP_OKAY;
}

/** creates a temporary presolving root node of the tree and installs it as focus node */
SCIP_RETCODE SCIPtreeCreatePresolvingRoot(
   SCIP_TREE*            tree,               /**< tree data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_Bool cutoff;

   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(tree->nsiblings == 0);
   assert(tree->root == NULL);
   assert(tree->focusnode == NULL);
   assert(!SCIPtreeProbing(tree));

   /* create temporary presolving root node */
   SCIP_CALL( SCIPtreeCreateRoot(tree, blkmem, set, stat, lp) );
   assert(tree->root != NULL);

   /* install the temporary root node as focus node */
   SCIP_CALL( SCIPnodeFocus(&tree->root, blkmem, set, stat, prob, primal, tree, lp, branchcand, conflict,
         eventfilter, eventqueue, &cutoff) );
   assert(!cutoff);

   return SCIP_OKAY;
}

/** frees the temporary presolving root and resets tree data structure */
SCIP_RETCODE SCIPtreeFreePresolvingRoot(
   SCIP_TREE*            tree,               /**< tree data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_NODE* node;
   SCIP_Bool cutoff;

   assert(tree != NULL);
   assert(tree->root != NULL);
   assert(tree->focusnode == tree->root);
   assert(tree->pathlen == 1);

   /* unfocus the temporary root node */
   node = NULL;
   SCIP_CALL( SCIPnodeFocus(&node, blkmem, set, stat, prob, primal, tree, lp, branchcand, conflict, 
         eventfilter, eventqueue, &cutoff) );
   assert(!cutoff);
   assert(tree->root == NULL);
   assert(tree->focusnode == NULL);
   assert(tree->pathlen == 0);

   /** reset tree data structure */
   SCIP_CALL( SCIPtreeClear(tree, blkmem, set, lp) );

   return SCIP_OKAY;
}

/** returns the node selector associated with the given node priority queue */
SCIP_NODESEL* SCIPtreeGetNodesel(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqGetNodesel(tree->leaves);
}

/** sets the node selector used for sorting the nodes in the priority queue, and resorts the queue if necessary */
SCIP_RETCODE SCIPtreeSetNodesel(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODESEL*         nodesel             /**< node selector to use for sorting the nodes in the queue */
   )
{
   assert(tree != NULL);
   assert(stat != NULL);

   if( SCIPnodepqGetNodesel(tree->leaves) != nodesel )
   {
      /* change the node selector used in the priority queue and resort the queue */
      SCIP_CALL( SCIPnodepqSetNodesel(&tree->leaves, set, nodesel) );

      /* issue message */
      if( stat->nnodes > 0 )
      {
         SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %"SCIP_LONGINT_FORMAT") switching to node selector <%s>\n", stat->nnodes, SCIPnodeselGetName(nodesel));
      }
   }

   return SCIP_OKAY;
}

/** cuts off nodes with lower bound not better than given cutoff bound */
SCIP_RETCODE SCIPtreeCutoff(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   )
{
   SCIP_NODE* node;
   int i;

   assert(tree != NULL);
   assert(stat != NULL);
   assert(lp != NULL);

   /* if we are in diving mode, it is not allowed to cut off nodes, because this can lead to deleting LP rows which
    * would modify the currently unavailable (due to diving modifications) SCIP_LP
    *  -> the cutoff must be delayed and executed after the diving ends
    */
   if( SCIPlpDiving(lp) )
   {
      tree->cutoffdelayed = TRUE;
      return SCIP_OKAY;
   }

   tree->cutoffdelayed = FALSE;

   /* cut off leaf nodes in the queue */
   SCIP_CALL( SCIPnodepqBound(tree->leaves, blkmem, set, stat, tree, lp, cutoffbound) );

   /* cut off siblings: we have to loop backwards, because a removal leads to moving the last node in empty slot */
   for( i = tree->nsiblings-1; i >= 0; --i )
   {
      node = tree->siblings[i];
      if( SCIPsetIsGE(set, node->lowerbound, cutoffbound) )
      {
         SCIPdebugMessage("cut off sibling %p at depth %d with lowerbound=%g at position %d\n", 
            node, node->depth, node->lowerbound, i);
         SCIPvbcCutoffNode(stat->vbc, stat, node);
         SCIP_CALL( SCIPnodeFree(&node, blkmem, set, tree, lp) );
      }
   }

   /* cut off children: we have to loop backwards, because a removal leads to moving the last node in empty slot */
   for( i = tree->nchildren-1; i >= 0; --i )
   {
      node = tree->children[i];
      if( SCIPsetIsGE(set, node->lowerbound, cutoffbound) )
      {
         SCIPdebugMessage("cut off child %p at depth %d with lowerbound=%g at position %d\n",
            node, node->depth, node->lowerbound, i);
         SCIPvbcCutoffNode(stat->vbc, stat, node);
         SCIP_CALL( SCIPnodeFree(&node, blkmem, set, tree, lp) );
      }
   }

   return SCIP_OKAY;
}

/** branches on a variable; if solution value x' is fractional, two child nodes are created
 *  (x <= floor(x'), x >= ceil(x')), if solution value is integral, three child nodes are created
 *  (x <= x'-1, x == x', x >= x'+1)
 */
SCIP_RETCODE SCIPtreeBranchVar(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_BRANCHDIR        branchdir           /**< preferred branching direction; maybe overridden by user settings */
   )
{
   SCIP_NODE* node;
   SCIP_Real solval;
   SCIP_Real downprio;

   assert(tree != NULL);
   assert(set != NULL);
   assert(var != NULL);

   /* get the corresponding active problem variable */
   var = SCIPvarGetProbvar(var);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIPerrorMessage("cannot branch on fixed or multi-aggregated variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   assert(SCIPvarIsActive(var));
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY
      || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER
      || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);
   assert(SCIPsetIsFeasIntegral(set, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsFeasIntegral(set, SCIPvarGetUbLocal(var)));
   assert(SCIPsetIsLT(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

   solval = SCIPvarGetSol(var, tree->focusnodehaslp);
   assert(SCIPsetIsGE(set, solval, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsLE(set, solval, SCIPvarGetUbLocal(var)));

   /* choose preferred branching direction */
   switch( SCIPvarGetBranchDirection(var) )  /* first, look at the user settings for the specific variable */
   {
   case SCIP_BRANCHDIR_DOWNWARDS:
      downprio = 1.0;
      break;
   case SCIP_BRANCHDIR_UPWARDS:
      downprio = -1.0;
      break;
   case SCIP_BRANCHDIR_AUTO:
      switch( branchdir )                    /* second, look at what the caller wants */
      {
      case SCIP_BRANCHDIR_DOWNWARDS:
         downprio = 1.0;
         break;
      case SCIP_BRANCHDIR_UPWARDS:
         downprio = -1.0;
         break;
      case SCIP_BRANCHDIR_AUTO:
         if( SCIPtreeHasFocusNodeLP(tree) )
            downprio = SCIPvarGetRootSol(var) - solval;
         else
            downprio = (SCIP_Real)SCIPvarGetAvgInferences(var, stat, SCIP_BRANCHDIR_DOWNWARDS)
               - (SCIP_Real)SCIPvarGetAvgInferences(var, stat, SCIP_BRANCHDIR_UPWARDS);
         break;
      default:
         SCIPerrorMessage("invalid preferred automatic branching direction <%d> specified\n", SCIPvarGetBranchDirection(var));
         return SCIP_INVALIDDATA;
      }
      break;
   default:
      SCIPerrorMessage("invalid preferred branching direction <%d> of variable <%s>\n", 
         SCIPvarGetBranchDirection(var), SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   if( SCIPsetIsFeasIntegral(set, solval) )
   {
      SCIP_Real fixval;

      /* create child nodes with x <= x'-1, x = x', and x >= x'+1;
       * set the node selection priority in a way, s.t. a node is preferred whose branching goes in the same direction
       * as the deviation from the variable's root solution; evaluate x = x' first in any way
       */
      fixval = SCIPsetFeasCeil(set, solval);
      assert(SCIPsetIsEQ(set, SCIPsetFeasCeil(set, solval), SCIPsetFeasFloor(set, solval)));
      
      SCIPdebugMessage("pseudo branch on variable <%s> with value %g, priority %d\n", 
         SCIPvarGetName(var), solval, SCIPvarGetBranchPriority(var));
      
      /* create child node with x <= x'-1, if this would be feasible */
      if( SCIPsetIsGE(set, fixval-1, SCIPvarGetLbLocal(var)) )
      {
         SCIPdebugMessage(" -> creating child: <%s> <= %g\n", SCIPvarGetName(var), fixval-1);
         SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, downprio) );
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, tree, lp, branchcand, eventqueue, 
               var, fixval-1, SCIP_BOUNDTYPE_UPPER, FALSE) );
         SCIPdebugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
      }
                  
      /* create child node with x = x' */
      SCIPdebugMessage(" -> creating child: <%s> == %g\n", SCIPvarGetName(var), fixval);
      SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, 0.0) );
      if( !SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), fixval) )
      {
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, tree, lp, branchcand, eventqueue, 
               var, fixval, SCIP_BOUNDTYPE_LOWER, FALSE) );
      }
      if( !SCIPsetIsEQ(set, SCIPvarGetUbLocal(var), fixval) )
      {
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, tree, lp, branchcand, eventqueue, 
               var, fixval, SCIP_BOUNDTYPE_UPPER, FALSE) );
      }
      SCIPdebugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
      
      /* create child node with x >= x'+1, if this would be feasible */
      if( SCIPsetIsLE(set, fixval+1, SCIPvarGetUbLocal(var)) )
      {
         SCIPdebugMessage(" -> creating child: <%s> >= %g\n", SCIPvarGetName(var), fixval+1);
         SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, -downprio) );
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, tree, lp, branchcand, eventqueue, 
               var, fixval+1, SCIP_BOUNDTYPE_LOWER, FALSE) );
         SCIPdebugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
      }
   }
   else
   {   
      /* create child nodes with x <= floor(x'), and x >= ceil(x');
       * set the node selection priority in a way, s.t. a node is preferred whose branching goes in the same direction
       * as the deviation from the variable's root solution
       */
      SCIPdebugMessage("LP branch on variable <%s> with value %g, priority %d\n", 
         SCIPvarGetName(var), solval, SCIPvarGetBranchPriority(var));
      
      /* create child node with x <= floor(x') */
      SCIPdebugMessage(" -> creating child: <%s> <= %g\n", SCIPvarGetName(var), SCIPsetFeasFloor(set, solval));
      SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, downprio) );
      SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, tree, lp, branchcand, eventqueue, 
            var, SCIPsetFeasFloor(set, solval), SCIP_BOUNDTYPE_UPPER, FALSE) );
      SCIPdebugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
      
      /* create child node with x >= ceil(x') */
      SCIPdebugMessage(" -> creating child: <%s> >= %g\n", SCIPvarGetName(var), SCIPsetFeasCeil(set, solval));
      SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, -downprio) );
      SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, tree, lp, branchcand, eventqueue, 
            var, SCIPsetFeasCeil(set, solval), SCIP_BOUNDTYPE_LOWER, FALSE) );
      SCIPdebugMessage(" -> child's lowerbound: %g\n", node->lowerbound);
   }

   return SCIP_OKAY;
}

/** creates a probing child node of the current node, which must be the focus node, the current refocused node,
 *  or another probing node; if the current node is the focus or a refocused node, the created probing node is
 *  installed as probing root node
 */
static
SCIP_RETCODE treeCreateProbingNode(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_NODE* currentnode;
   SCIP_NODE* node;

   assert(tree != NULL);
   assert(SCIPtreeIsPathComplete(tree));
   assert(tree->pathlen > 0);
   assert(blkmem != NULL);
   assert(set != NULL);

   /* get the current node */
   currentnode = SCIPtreeGetCurrentNode(tree);
   assert(SCIPnodeGetType(currentnode) == SCIP_NODETYPE_FOCUSNODE
      || SCIPnodeGetType(currentnode) == SCIP_NODETYPE_REFOCUSNODE
      || SCIPnodeGetType(currentnode) == SCIP_NODETYPE_PROBINGNODE);
   assert((SCIPnodeGetType(currentnode) == SCIP_NODETYPE_PROBINGNODE) == SCIPtreeProbing(tree));

   /* create the node data structure */
   SCIP_CALL( nodeCreate(&node, blkmem, set) );
   assert(node != NULL);

   /* mark node to be a probing node */
   node->nodetype = SCIP_NODETYPE_PROBINGNODE; /*lint !e641*/

   /* make the current node the parent of the new probing node */
   SCIP_CALL( nodeAssignParent(node, blkmem, set, tree, currentnode, 0.0) );
   assert(SCIPnodeGetDepth(node) == tree->pathlen);

   /* check, if the node is the probing root node */
   if( tree->probingroot == NULL )
   {
      tree->probingroot = node;
      SCIPdebugMessage("created probing root node %p at depth %d\n", node, SCIPnodeGetDepth(node));
   }
   else
   {
      assert(SCIPnodeGetType(tree->probingroot) == SCIP_NODETYPE_PROBINGNODE);
      assert(SCIPnodeGetDepth(tree->probingroot) < SCIPnodeGetDepth(node));

      SCIPdebugMessage("created probing child node %p at depth %d, probing depth %d\n", 
         node, SCIPnodeGetDepth(node), SCIPnodeGetDepth(node) - SCIPnodeGetDepth(tree->probingroot));
   }

   /* create the new active path */
   SCIP_CALL( treeEnsurePathMem(tree, set, tree->pathlen+1) );
   node->active = TRUE;
   tree->path[tree->pathlen] = node;
   tree->pathlen++;

   /* count the new LP sizes of the path */
   treeUpdatePathLPSize(tree, tree->pathlen-1);

   /* the current probing node does not yet a solved LP */
   tree->probingnodehaslp = FALSE;

   return SCIP_OKAY;
}

/** switches to probing mode and creates a probing root */
SCIP_RETCODE SCIPtreeStartProbing(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);
   assert(tree->probinglpistate == NULL);
   assert(!SCIPtreeProbing(tree));
   assert(lp != NULL);

   SCIPdebugMessage("probing started in depth %d (LP flushed: %d, LP solved: %d, solstat: %d), probing root in depth %d\n",
      tree->pathlen-1, lp->flushed, lp->solved, SCIPlpGetSolstat(lp), tree->pathlen);

   /* inform LP about probing mode */
   SCIP_CALL( SCIPlpStartProbing(lp) );

   /* remember, whether the LP was flushed and solved */
   if( set->stage == SCIP_STAGE_SOLVING )
   {
      tree->probinglpwasflushed = lp->flushed;
      tree->probinglpwassolved = lp->solved;

      /* remember the LP state in order to restore the LP solution quickly after probing */
      if( lp->flushed && lp->solved )
      {
         SCIP_CALL( SCIPlpGetState(lp, blkmem, &tree->probinglpistate) );
      }
   }

   /* create temporary probing root node */
   SCIP_CALL( treeCreateProbingNode(tree, blkmem, set) );
   assert(SCIPtreeProbing(tree));

   return SCIP_OKAY;
}

/** creates a new probing child node in the probing path */
SCIP_RETCODE SCIPtreeCreateProbingNode(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(SCIPtreeProbing(tree));

   SCIPdebugMessage("new probing child in depth %d (probing depth: %d)\n",
      tree->pathlen, tree->pathlen-1 - SCIPnodeGetDepth(tree->probingroot));

   /* create temporary probing root node */
   SCIP_CALL( treeCreateProbingNode(tree, blkmem, set) );

   return SCIP_OKAY;
}

/** undoes all changes to the problem applied in probing up to the given probing depth */
static
SCIP_RETCODE treeBacktrackProbing(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   probingdepth        /**< probing depth of the node in the probing path that should be reactivated,
                                              *   -1 to even deactivate the probing root, thus exiting probing mode */
   )
{
   int newpathlen;

   assert(tree != NULL);
   assert(SCIPtreeProbing(tree));
   assert(tree->probingroot != NULL);
   assert(tree->focusnode != NULL);
   assert(SCIPnodeGetType(tree->probingroot) == SCIP_NODETYPE_PROBINGNODE);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE
      || SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_REFOCUSNODE);
   assert(tree->probingroot->parent == tree->focusnode);
   assert(SCIPnodeGetDepth(tree->probingroot) == SCIPnodeGetDepth(tree->focusnode)+1);
   assert(tree->pathlen >= 2);
   assert(SCIPnodeGetType(tree->path[tree->pathlen-1]) == SCIP_NODETYPE_PROBINGNODE);
   assert(-1 <= probingdepth && probingdepth <= SCIPtreeGetProbingDepth(tree));

   newpathlen = SCIPnodeGetDepth(tree->probingroot) + probingdepth + 1;
   assert(newpathlen >= 1); /* at least root node of the tree remains active */
   while( tree->pathlen > newpathlen )
   {
      assert(SCIPnodeGetType(tree->path[tree->pathlen-1]) == SCIP_NODETYPE_PROBINGNODE);
      assert(tree->pathlen-1 == SCIPnodeGetDepth(tree->path[tree->pathlen-1]));
      assert(tree->pathlen-1 >= SCIPnodeGetDepth(tree->probingroot));

      /* undo bound changes by deactivating the probing node */
      SCIP_CALL( nodeDeactivate(tree->path[tree->pathlen-1], blkmem, set, stat, tree, lp, branchcand, eventqueue) );

      /* free the probing node */
      SCIP_CALL( SCIPnodeFree(&tree->path[tree->pathlen-1], blkmem, set, tree, lp) );
      tree->pathlen--;
   }
   assert(tree->pathlen == newpathlen);

   /* if the highest cutoff or repropagation depth is inside the deleted part of the probing path,
    * reset them to infinity
    */
   if( tree->cutoffdepth >= tree->pathlen )
      tree->cutoffdepth = INT_MAX;
   if( tree->repropdepth >= tree->pathlen )
      tree->repropdepth = INT_MAX;

   return SCIP_OKAY;
}

/** undoes all changes to the problem applied in probing up to the given probing depth;
 *  the changes of the probing node of the given probing depth are the last ones that remain active;
 *  changes that were applied before calling SCIPtreeCreateProbingNode() cannot be undone
 */
SCIP_RETCODE SCIPtreeBacktrackProbing(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   probingdepth        /**< probing depth of the node in the probing path that should be reactivated */
   )
{
   assert(tree != NULL);
   assert(SCIPtreeProbing(tree));
   assert(0 <= probingdepth && probingdepth <= SCIPtreeGetProbingDepth(tree));

   /* undo the domain and constraint set changes and free the temporary probing nodes below the given probing depth */
   SCIP_CALL( treeBacktrackProbing(tree, blkmem, set, stat, lp, branchcand, eventqueue, probingdepth) );

   assert(SCIPtreeProbing(tree));
   assert(SCIPnodeGetType(SCIPtreeGetCurrentNode(tree)) == SCIP_NODETYPE_PROBINGNODE);

   return SCIP_OKAY;
}

/** switches back from probing to normal operation mode, frees all nodes on the probing path, restores bounds of all
 *  variables and restores active constraints arrays of focus node
 */
SCIP_RETCODE SCIPtreeEndProbing(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   assert(tree != NULL);
   assert(SCIPtreeProbing(tree));
   assert(tree->probingroot != NULL);
   assert(tree->focusnode != NULL);
   assert(SCIPnodeGetType(tree->probingroot) == SCIP_NODETYPE_PROBINGNODE);
   assert(SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_FOCUSNODE
      || SCIPnodeGetType(tree->focusnode) == SCIP_NODETYPE_REFOCUSNODE);
   assert(tree->probingroot->parent == tree->focusnode);
   assert(SCIPnodeGetDepth(tree->probingroot) == SCIPnodeGetDepth(tree->focusnode)+1);
   assert(tree->pathlen >= 2);
   assert(SCIPnodeGetType(tree->path[tree->pathlen-1]) == SCIP_NODETYPE_PROBINGNODE);
   assert(set != NULL);

   /* undo the domain and constraint set changes of the temporary probing nodes and free the probing nodes */
   SCIP_CALL( treeBacktrackProbing(tree, blkmem, set, stat, lp, branchcand, eventqueue, -1) );
   assert(SCIPtreeGetCurrentNode(tree) == tree->focusnode);
   assert(!SCIPtreeProbing(tree));

   /* if the LP was flushed before probing starts, flush it again */
   if( tree->probinglpwasflushed )
   {
      assert(set->stage == SCIP_STAGE_SOLVING);

      SCIP_CALL( SCIPlpFlush(lp, blkmem, set) );

      /* if the LP was solved before probing starts, solve it again to restore the LP solution */
      if( tree->probinglpwassolved )
      {
         SCIP_Bool lperror;
         
         /* reset the LP state before probing started */
         SCIP_CALL( SCIPlpSetState(lp, blkmem, set, tree->probinglpistate) );
         SCIP_CALL( SCIPlpFreeState(lp, blkmem, &tree->probinglpistate) );

         /* resolve LP to reset solution */
         SCIP_CALL( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, FALSE, FALSE, &lperror) );
         if( lperror )
         {
            SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "(node %"SCIP_LONGINT_FORMAT") unresolved numerical troubles while resolving LP %d after probing\n",
               stat->nnodes, stat->nlps);
            lp->resolvelperror = TRUE;
         }
      }
   }
   assert(tree->probinglpistate == NULL);
   tree->probinglpwasflushed = FALSE;
   tree->probinglpwassolved = FALSE;

   /* inform LP about end of probing mode */
   SCIP_CALL( SCIPlpEndProbing(lp) );

   SCIPdebugMessage("probing ended in depth %d (LP flushed: %d, solstat: %d)\n",
      tree->pathlen-1, lp->flushed, SCIPlpGetSolstat(lp));
   
   return SCIP_OKAY;
}

/** gets the best child of the focus node w.r.t. the node selection priority assigned by the branching rule */
SCIP_NODE* SCIPtreeGetPrioChild(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   SCIP_NODE* bestnode;
   SCIP_Real bestprio;
   int i;

   assert(tree != NULL);

   bestnode = NULL;
   bestprio = SCIP_REAL_MIN;
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

/** gets the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule */
SCIP_NODE* SCIPtreeGetPrioSibling(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   SCIP_NODE* bestnode;
   SCIP_Real bestprio;
   int i;

   assert(tree != NULL);

   bestnode = NULL;
   bestprio = SCIP_REAL_MIN;
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

/** gets the best child of the focus node w.r.t. the node selection strategy */
SCIP_NODE* SCIPtreeGetBestChild(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODE* bestnode;
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

/** gets the best sibling of the focus node w.r.t. the node selection strategy */
SCIP_NODE* SCIPtreeGetBestSibling(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODE* bestnode;
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
SCIP_NODE* SCIPtreeGetBestLeaf(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqFirst(tree->leaves);
}

/** gets the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy */
SCIP_NODE* SCIPtreeGetBestNode(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODE* bestchild;
   SCIP_NODE* bestsibling;
   SCIP_NODE* bestleaf;
   SCIP_NODE* bestnode;

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
SCIP_Real SCIPtreeGetLowerbound(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_Real lowerbound;
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

   /* compare lower bound with focus node */
   if( tree->focusnode != NULL )
   {
      lowerbound = MIN(lowerbound, tree->focusnode->lowerbound);
   }

   return lowerbound;
}

/** gets the node with minimal lower bound of all nodes in the tree (child, sibling, or leaf) */
SCIP_NODE* SCIPtreeGetLowerboundNode(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_NODE* lowerboundnode;
   SCIP_Real lowerbound;
   SCIP_Real bestprio;
   int i;

   assert(tree != NULL);
   assert(set != NULL);

   /* get the lower bound from the queue */
   lowerboundnode = SCIPnodepqGetLowerboundNode(tree->leaves, set);
   lowerbound = lowerboundnode != NULL ? lowerboundnode->lowerbound : SCIPsetInfinity(set);
   bestprio = -SCIPsetInfinity(set);

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
SCIP_Real SCIPtreeGetAvgLowerbound(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Real             cutoffbound         /**< global cutoff bound */
   )
{
   SCIP_Real lowerboundsum;
   int nnodes;
   int i;

   assert(tree != NULL);

   /* get sum of lower bounds from nodes in the queue */
   lowerboundsum = SCIPnodepqGetLowerboundSum(tree->leaves);
   nnodes = SCIPtreeGetNLeaves(tree);

   /* add lower bound of focus node */
   if( tree->focusnode != NULL && tree->focusnode->lowerbound < cutoffbound )
   {
      lowerboundsum += tree->focusnode->lowerbound;
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




/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPnodeGetType
#undef SCIPnodeGetNumber
#undef SCIPnodeGetDepth
#undef SCIPnodeGetLowerbound
#undef SCIPnodeGetPriority
#undef SCIPnodeIsActive
#undef SCIPnodeIsPropagatedAgain
#undef SCIPtreeGetNLeaves
#undef SCIPtreeGetNChildren
#undef SCIPtreeGetNSiblings
#undef SCIPtreeGetNNodes
#undef SCIPtreeIsPathComplete
#undef SCIPtreeProbing
#undef SCIPtreeGetProbingRoot
#undef SCIPtreeGetProbingDepth
#undef SCIPtreeGetFocusNode
#undef SCIPtreeGetFocusDepth
#undef SCIPtreeHasFocusNodeLP
#undef SCIPtreeSetFocusNodeLP
#undef SCIPtreeMarkProbingNodeHasLP
#undef SCIPtreeIsFocusNodeLPConstructed
#undef SCIPtreeGetCurrentNode
#undef SCIPtreeGetCurrentDepth
#undef SCIPtreeHasCurrentNodeLP

/** gets the type of the node */
SCIP_NODETYPE SCIPnodeGetType(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return (SCIP_NODETYPE)(node->nodetype);
}

/** gets successively assigned number of the node */
SCIP_Longint SCIPnodeGetNumber(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->number;
}

/** gets the depth of the node */
int SCIPnodeGetDepth(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->depth;
}

/** gets the lower bound of the node */
SCIP_Real SCIPnodeGetLowerbound(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->lowerbound;
}

/** gets the node selection priority of the node assigned by the branching rule */
SCIP_Real SCIPnodeGetPriority(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->priority;
}

/** returns whether node is in the path to the current node */
SCIP_Bool SCIPnodeIsActive(
   SCIP_NODE*            node                /**< node */
   )
{
   assert(node != NULL);

   return node->active;
}

/** returns whether the node is marked to be propagated again */
SCIP_Bool SCIPnodeIsPropagatedAgain(
   SCIP_NODE*            node                /**< node data */
   )
{
   assert(node != NULL);

   return node->reprop;
}

/** gets number of children */
int SCIPtreeGetNChildren(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->nchildren;
}

/** gets number of siblings */
int SCIPtreeGetNSiblings(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->nsiblings;
}

/** gets number of leaves */
int SCIPtreeGetNLeaves(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return SCIPnodepqLen(tree->leaves);
}
   
/** gets number of nodes (children + siblings + leaves) */
int SCIPtreeGetNNodes(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->nchildren + tree->nsiblings + SCIPtreeGetNLeaves(tree);
}

/** returns whether the active path goes completely down to the focus node */
SCIP_Bool SCIPtreeIsPathComplete(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->focusnode != NULL || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->focusnode != NULL);
   assert(tree->pathlen >= 2 || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1] != NULL);
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1]->depth == tree->pathlen-1);
   assert(tree->focusnode == NULL || (int)tree->focusnode->depth >= tree->pathlen
      || tree->path[tree->focusnode->depth] == tree->focusnode);

   return (tree->focusnode == NULL || (int)tree->focusnode->depth < tree->pathlen);
}

/** returns whether the current node is a temporary probing node */
SCIP_Bool SCIPtreeProbing(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->probingroot == NULL || (SCIP_NODETYPE)tree->probingroot->nodetype == SCIP_NODETYPE_PROBINGNODE);
   assert(tree->probingroot == NULL || tree->pathlen > SCIPnodeGetDepth(tree->probingroot));
   assert(tree->probingroot == NULL || tree->path[SCIPnodeGetDepth(tree->probingroot)] == tree->probingroot);

   return (tree->probingroot != NULL);
}

/** returns the temporary probing root node, or NULL if the we are not in probing mode */
SCIP_NODE* SCIPtreeGetProbingRoot(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->probingroot == NULL || (SCIP_NODETYPE)tree->probingroot->nodetype == SCIP_NODETYPE_PROBINGNODE);
   assert(tree->probingroot == NULL || tree->pathlen > SCIPnodeGetDepth(tree->probingroot));
   assert(tree->probingroot == NULL || tree->path[SCIPnodeGetDepth(tree->probingroot)] == tree->probingroot);

   return tree->probingroot;
}

/** gets focus node of the tree */
SCIP_NODE* SCIPtreeGetFocusNode(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->focusnode != NULL || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->focusnode != NULL);
   assert(tree->pathlen >= 2 || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1] != NULL);
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1]->depth == tree->pathlen-1);
   assert(tree->focusnode == NULL || (int)tree->focusnode->depth >= tree->pathlen
      || tree->path[tree->focusnode->depth] == tree->focusnode);

   return tree->focusnode;
}

/** gets depth of focus node in the tree */
int SCIPtreeGetFocusDepth(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->focusnode != NULL || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->focusnode != NULL);
   assert(tree->pathlen >= 2 || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1] != NULL);
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1]->depth == tree->pathlen-1);
   assert(tree->focusnode == NULL || (int)tree->focusnode->depth >= tree->pathlen
      || tree->path[tree->focusnode->depth] == tree->focusnode);

   return tree->focusnode != NULL ? (int)tree->focusnode->depth : -1;
}

/** returns, whether the LP was or is to be solved in the focus node */
SCIP_Bool SCIPtreeHasFocusNodeLP(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->focusnodehaslp;
}

/** sets mark to solve or to ignore the LP while processing the focus node */
void SCIPtreeSetFocusNodeLP(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             solvelp             /**< should the LP be solved in focus node? */
   )
{
   assert(tree != NULL);

   tree->focusnodehaslp = solvelp;
}

/** marks the probing node to have a solved LP relaxation */
void SCIPtreeMarkProbingNodeHasLP(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   tree->probingnodehaslp = TRUE;
}

/** returns whether the LP of the focus node is already constructed */
SCIP_Bool SCIPtreeIsFocusNodeLPConstructed(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);

   return tree->focuslpconstructed;
}

/** gets current node of the tree, i.e. the last node in the active path, or NULL if no current node exists */
SCIP_NODE* SCIPtreeGetCurrentNode(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->focusnode != NULL || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->focusnode != NULL);
   assert(tree->pathlen >= 2 || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1] != NULL);
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1]->depth == tree->pathlen-1);
   assert(tree->focusnode == NULL || (int)tree->focusnode->depth >= tree->pathlen
      || tree->path[tree->focusnode->depth] == tree->focusnode);

   return (tree->pathlen > 0 ? tree->path[tree->pathlen-1] : NULL);
}

/** gets depth of current node in the tree, i.e. the length of the active path minus 1, or -1 if no current node exists */
int SCIPtreeGetCurrentDepth(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(tree->focusnode != NULL || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->focusnode != NULL);
   assert(tree->pathlen >= 2 || !SCIPtreeProbing(tree));
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1] != NULL);
   assert(tree->pathlen == 0 || tree->path[tree->pathlen-1]->depth == tree->pathlen-1);
   assert(tree->focusnode == NULL || (int)tree->focusnode->depth >= tree->pathlen
      || tree->path[tree->focusnode->depth] == tree->focusnode);

   return tree->pathlen-1;
}

/** returns, whether the LP was or is to be solved in the current node */
SCIP_Bool SCIPtreeHasCurrentNodeLP(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(SCIPtreeIsPathComplete(tree));

   return SCIPtreeProbing(tree) ? tree->probingnodehaslp : SCIPtreeHasFocusNodeLP(tree);
}

/** returns the current probing depth, i.e. the number of probing sub nodes existing in the probing path */
int SCIPtreeGetProbingDepth(
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(tree != NULL);
   assert(SCIPtreeProbing(tree));

   return SCIPtreeGetCurrentDepth(tree) - SCIPnodeGetDepth(tree->probingroot);
}

