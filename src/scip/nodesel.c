/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nodesel.c,v 1.51 2006/09/17 01:58:42 bzfpfend Exp $"

/**@file   nodesel.c
 * @brief  methods for node selectors
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/vbc.h"
#include "scip/paramset.h"
#include "scip/tree.h"
#include "scip/scip.h"
#include "scip/nodesel.h"

#include "scip/struct_nodesel.h"



/* 
 * node priority queue methods
 */

#define PQ_PARENT(q) (((q)+1)/2-1)
#define PQ_LEFTCHILD(p) (2*(p)+1)
#define PQ_RIGHTCHILD(p) (2*(p)+2)

/** resizes node memory to hold at least the given number of nodes */
static
SCIP_RETCODE nodepqResize(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   minsize             /**< minimal number of storable nodes */
   )
{
   assert(nodepq != NULL);
   
   if( minsize <= nodepq->size )
      return SCIP_OKAY;

   nodepq->size = SCIPsetCalcTreeGrowSize(set, minsize);
   SCIP_ALLOC( BMSreallocMemoryArray(&nodepq->slots, nodepq->size) );

   return SCIP_OKAY;
}

/** updates the cached minimal lower bound of all nodes in the queue (used for node selection rules, that don't store
 *  the lowest bound node in the first slot of the queue)
 */
static
void nodepqUpdateLowerbound(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node to be inserted */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->nodesel != NULL);
   assert(!nodepq->nodesel->lowestboundfirst);
   assert(node != NULL);

   assert(nodepq->nlowerbounds > 0 || nodepq->lowerboundnode == NULL);
   SCIPdebugMessage("update queue's lower bound after adding node %p: nodebound=%g, queuebound=%g, nlowerbounds=%d, lowerboundnode=%p (bound=%g)\n",
      node, SCIPnodeGetLowerbound(node), nodepq->lowerbound, nodepq->nlowerbounds,
      nodepq->lowerboundnode, nodepq->lowerboundnode != NULL ? SCIPnodeGetLowerbound(nodepq->lowerboundnode) : 0.0);
   if( nodepq->validlowerbound )
   {
      assert(nodepq->lowerbound < SCIP_INVALID);
      if( SCIPsetIsLE(set, SCIPnodeGetLowerbound(node), nodepq->lowerbound) )
      {
         if( SCIPsetIsEQ(set, SCIPnodeGetLowerbound(node), nodepq->lowerbound) )
         {
            assert(nodepq->nlowerbounds >= 1);
            nodepq->nlowerbounds++;
         }
         else
         {
            nodepq->lowerboundnode = node;
            nodepq->lowerbound = SCIPnodeGetLowerbound(node);
            nodepq->nlowerbounds = 1;
         }
      }
   }
   SCIPdebugMessage(" -> new queuebound=%g, nlowerbounds=%d, lowerboundnode=%p\n",
      nodepq->lowerbound, nodepq->nlowerbounds, nodepq->lowerboundnode);

   assert(nodepq->nlowerbounds > 0 || nodepq->lowerboundnode == NULL);
}

/** calculates the minimal lower bound of all nodes in the queue (used for node selection rules, that don't store
 *  the lowest bound node in the first slot of the queue)
 */
static
void nodepqCalcLowerbound(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(nodepq != NULL);
   assert(nodepq->nodesel != NULL);
   assert(!nodepq->nodesel->lowestboundfirst);

   nodepq->validlowerbound = TRUE;
   nodepq->lowerboundnode = NULL;
   nodepq->lowerbound = SCIPsetInfinity(set);
   nodepq->nlowerbounds = 0;

   for( i = 0; i < nodepq->len; ++i )
      nodepqUpdateLowerbound(nodepq, set, nodepq->slots[i]);
}

/** creates node priority queue */
SCIP_RETCODE SCIPnodepqCreate(
   SCIP_NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODESEL*         nodesel             /**< node selector to use for sorting the nodes in the queue */
   )
{
   assert(nodepq != NULL);

   SCIP_ALLOC( BMSallocMemory(nodepq) );
   (*nodepq)->nodesel = nodesel;
   (*nodepq)->slots = NULL;
   (*nodepq)->len = 0;
   (*nodepq)->size = 0;
   (*nodepq)->lowerboundnode = NULL;
   (*nodepq)->lowerboundsum = 0.0;
   (*nodepq)->lowerbound = SCIPsetInfinity(set);
   (*nodepq)->nlowerbounds = 0;
   (*nodepq)->validlowerbound = TRUE;

   return SCIP_OKAY;
}

/** frees node priority queue, but not the data nodes themselves */
void SCIPnodepqDestroy(
   SCIP_NODEPQ**         nodepq              /**< pointer to a node priority queue */
   )
{
   assert(nodepq != NULL);
   assert(*nodepq != NULL);

   BMSfreeMemoryArrayNull(&(*nodepq)->slots);
   BMSfreeMemory(nodepq);
}

/** frees node priority queue and all nodes in the queue */
SCIP_RETCODE SCIPnodepqFree(
   SCIP_NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(nodepq != NULL);
   assert(*nodepq != NULL);

   /* free the nodes of the queue */
   SCIP_CALL( SCIPnodepqClear(*nodepq, blkmem, set, tree, lp) );
   
   /* free the queue data structure */
   SCIPnodepqDestroy(nodepq);

   return SCIP_OKAY;
}

/** deletes all nodes in the node priority queue */
SCIP_RETCODE SCIPnodepqClear(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(nodepq != NULL);

   /* free the nodes of the queue */
   for( i = 0; i < nodepq->len; ++i )
   {
      assert(nodepq->slots[i] != NULL);
      assert(SCIPnodeGetType(nodepq->slots[i]) == SCIP_NODETYPE_LEAF);
      SCIP_CALL( SCIPnodeFree(&nodepq->slots[i], blkmem, set, tree, lp) );
   }

   /* reset data */
   nodepq->len = 0;
   nodepq->lowerboundnode = NULL;
   nodepq->lowerboundsum = 0.0;
   nodepq->lowerbound = SCIPsetInfinity(set);
   nodepq->nlowerbounds = 0;
   nodepq->validlowerbound = TRUE;

   return SCIP_OKAY;
}

/** returns the node selector associated with the given node priority queue */
SCIP_NODESEL* SCIPnodepqGetNodesel(
   SCIP_NODEPQ*          nodepq              /**< node priority queue */
   )
{
   assert(nodepq != NULL);

   return nodepq->nodesel;
}

/** sets the node selector used for sorting the nodes in the queue, and resorts the queue if necessary */
SCIP_RETCODE SCIPnodepqSetNodesel(
   SCIP_NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODESEL*         nodesel             /**< node selector to use for sorting the nodes in the queue */
   )
{
   assert(nodepq != NULL);
   assert(*nodepq != NULL);
   assert((*nodepq)->len >= 0);
   assert(nodesel != NULL);
   assert(nodesel->nodeselcomp != NULL);

   if( (*nodepq)->nodesel != nodesel )
   {
      SCIP_NODEPQ* newnodepq;
      int i;

      /* create new node priority queue */
      SCIP_CALL( SCIPnodepqCreate(&newnodepq, set, nodesel) );
      
      /* resize the new node priority queue to be able to store all nodes */
      SCIP_CALL( nodepqResize(newnodepq, set, (*nodepq)->len) );
      
      /* insert all nodes in the new node priority queue */
      for( i = 0; i < (*nodepq)->len; ++i )
      {
         SCIP_CALL( SCIPnodepqInsert(newnodepq, set, (*nodepq)->slots[i]) );
      }
      
      /* destroy the old node priority queue without freeing the nodes */
      SCIPnodepqDestroy(nodepq);
      
      /* use the new node priority queue */
      *nodepq = newnodepq;
   }

   return SCIP_OKAY;
}

/** compares two nodes; returns -1/0/+1 if node1 better/equal/worse than node2 */
int SCIPnodepqCompare(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node1,              /**< first node to compare */
   SCIP_NODE*            node2               /**< second node to compare */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->nodesel != NULL);
   assert(nodepq->nodesel->nodeselcomp != NULL);
   assert(set != NULL);

   return SCIPnodeselCompare(nodepq->nodesel, set, node1, node2);
}

/** inserts node into node priority queue */
SCIP_RETCODE SCIPnodepqInsert(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node to be inserted */
   )
{
   SCIP_NODESEL* nodesel;
   int pos;

   assert(nodepq != NULL);
   assert(nodepq->len >= 0);
   assert(set != NULL);
   assert(node != NULL);

   nodesel = nodepq->nodesel;
   assert(nodesel != NULL);
   assert(nodesel->nodeselcomp != NULL);

   SCIP_CALL( nodepqResize(nodepq, set, nodepq->len+1) );

   /* insert node as leaf in the tree, move it towards the root as long it is better than its parent */
   pos = nodepq->len;
   nodepq->len++;
   nodepq->lowerboundsum += SCIPnodeGetLowerbound(node);
   while( pos > 0 && nodesel->nodeselcomp(set->scip, nodesel, node, nodepq->slots[PQ_PARENT(pos)]) < 0 )
   {
      nodepq->slots[pos] = nodepq->slots[PQ_PARENT(pos)];
      pos = PQ_PARENT(pos);
   }
   nodepq->slots[pos] = node;

   if( !nodesel->lowestboundfirst )
   {
      /* update the minimal lower bound */
      nodepqUpdateLowerbound(nodepq, set, node);
   }

   return SCIP_OKAY;
}

/** deletes node at given position from the node priority queue; returns TRUE, if the parent fell down to the
 *  free position
 */
static
SCIP_Bool nodepqDelPos(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   rempos              /**< queue position of node to remove */
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODE* lastnode;
   int freepos;
   int childpos;
   int parentpos;
   int brotherpos;
   SCIP_Bool parentfelldown;

   assert(nodepq != NULL);
   assert(nodepq->len > 0);
   assert(set != NULL);
   assert(0 <= rempos && rempos < nodepq->len);

   nodesel = nodepq->nodesel;
   assert(nodesel != NULL);
   assert(nodesel->nodeselcomp != NULL);

   if( !nodesel->lowestboundfirst )
   {
      assert(nodepq->nlowerbounds > 0 || nodepq->lowerboundnode == NULL);

      /* update the minimal lower bound */
      if( nodepq->nlowerbounds > 0 )
      {
         SCIP_NODE* node;
         
         node = nodepq->slots[rempos];
         assert(SCIPsetIsGE(set, SCIPnodeGetLowerbound(node), nodepq->lowerbound));
         
         SCIPdebugMessage("update queue's lower bound after removal of node %p: nodebound=%g, queuebound=%g, nlowerbounds=%d, lowerboundnode=%p (bound=%g)\n",
            node, SCIPnodeGetLowerbound(node), nodepq->lowerbound, nodepq->nlowerbounds, 
            nodepq->lowerboundnode, nodepq->lowerboundnode != NULL ? SCIPnodeGetLowerbound(nodepq->lowerboundnode) : 0.0);
         if( SCIPsetIsEQ(set, SCIPnodeGetLowerbound(node), nodepq->lowerbound) )
         {
            nodepq->nlowerbounds--;
            if( nodepq->nlowerbounds == 0 )
            {
               nodepq->validlowerbound = FALSE;
               nodepq->lowerbound = SCIP_INVALID;
            }
         }
         if( node == nodepq->lowerboundnode )
            nodepq->lowerboundnode = NULL;
         SCIPdebugMessage(" -> new queuebound=%g, nlowerbounds=%d, lowerboundnode=%p\n",
            nodepq->lowerbound, nodepq->nlowerbounds, nodepq->lowerboundnode);
      }
      assert(nodepq->nlowerbounds > 0 || nodepq->lowerboundnode == NULL);
   }

   /* remove node of the tree and get a free slot,
    * if the removed node was the last node of the queue
    *  - do nothing
    * if the last node of the queue is better than the parent of the removed node:
    *  - move the parent to the free slot, until the last node can be placed in the free slot
    * if the last node of the queue is not better than the parent of the free slot:
    *  - move the better child to the free slot until the last node can be placed in the free slot
    */
   nodepq->lowerboundsum -= SCIPnodeGetLowerbound(nodepq->slots[rempos]);
   freepos = rempos;
   lastnode = nodepq->slots[nodepq->len-1];
   nodepq->len--;

   if( freepos == nodepq->len )
      return FALSE;
   assert(freepos < nodepq->len);

   /* try to move parents downwards to insert last node */
   parentfelldown = FALSE;
   parentpos = PQ_PARENT(freepos);
   while( freepos > 0 && nodesel->nodeselcomp(set->scip, nodesel, lastnode, nodepq->slots[parentpos]) < 0 )
   {
      nodepq->slots[freepos] = nodepq->slots[parentpos];
      freepos = parentpos;
      parentpos = PQ_PARENT(freepos);
      parentfelldown = TRUE;
   }
   if( !parentfelldown )
   {
      /* downward moving of parents was not successful -> move children upwards */
      while( freepos <= PQ_PARENT(nodepq->len-1) ) /* as long as free slot has children... */
      {
         /* select the better child of free slot */
         childpos = PQ_LEFTCHILD(freepos);
         assert(childpos < nodepq->len);
         brotherpos = PQ_RIGHTCHILD(freepos);
         if( brotherpos < nodepq->len
            && nodesel->nodeselcomp(set->scip, nodesel, nodepq->slots[brotherpos], nodepq->slots[childpos]) < 0 )
            childpos = brotherpos;
         /* exit search loop if better child is not better than last node */
         if( nodesel->nodeselcomp(set->scip, nodesel, lastnode, nodepq->slots[childpos]) <= 0 )
            break;
         /* move better child upwards, free slot is now the better child's slot */
         nodepq->slots[freepos] = nodepq->slots[childpos];
         freepos = childpos;
      }
   }
   assert(0 <= freepos && freepos < nodepq->len);
   assert(!parentfelldown || PQ_LEFTCHILD(freepos) < nodepq->len);
   nodepq->slots[freepos] = lastnode;

   return parentfelldown;
}

/** returns the position of given node in the priority queue, or -1 if not existing */
static
int nodepqFindNode(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node to find */
   )
{
   int pos;

   assert(nodepq != NULL);
   assert(nodepq->len >= 0);
   assert(set != NULL);
   assert(node != NULL);

   /* search the node in the queue */
   for( pos = 0; pos < nodepq->len && node != nodepq->slots[pos]; ++pos )
   {}

   if( pos == nodepq->len )
      pos = -1;
   
   return pos;
}

/** removes node from the node priority queue */
SCIP_RETCODE SCIPnodepqRemove(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node to remove */
   )
{
   int pos;

   pos = nodepqFindNode(nodepq, set, node);
   if( pos == -1 )
   {
      SCIPerrorMessage("node doesn't exist in node priority queue\n");
      return SCIP_INVALIDDATA;
   }

   (void)nodepqDelPos(nodepq, set, pos);

   return SCIP_OKAY;
}

/** returns the best node of the queue without removing it */
SCIP_NODE* SCIPnodepqFirst(
   const SCIP_NODEPQ*    nodepq              /**< node priority queue */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->len >= 0);

   if( nodepq->len == 0 )
      return NULL;

   assert(nodepq->slots[0] != NULL);

   return nodepq->slots[0];
}

/** returns the nodes array of the queue */
SCIP_NODE** SCIPnodepqNodes(
   const SCIP_NODEPQ*    nodepq              /**< node priority queue */
   )
{
   assert(nodepq != NULL);

   return nodepq->slots;
}

/** returns the number of nodes stored in the node priority queue */
int SCIPnodepqLen(
   const SCIP_NODEPQ*    nodepq              /**< node priority queue */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->len >= 0);

   return nodepq->len;
}

/** gets the minimal lower bound of all nodes in the queue */
SCIP_Real SCIPnodepqGetLowerbound(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->nodesel != NULL);
   assert(set != NULL);

   if( nodepq->nodesel->lowestboundfirst )
   {
      /* the node selector's compare method sorts the minimal lower bound to the front */
      if( nodepq->len > 0 )
      {
         assert(nodepq->slots[0] != NULL);
         return SCIPnodeGetLowerbound(nodepq->slots[0]);
      }
      else
         return SCIPsetInfinity(set);
   }
   else
   {
      /* we use bookkeeping to remember the lowest bound */

      /* if the cached lower bound is invalid, calculate it */
      if( !nodepq->validlowerbound )
         nodepqCalcLowerbound(nodepq, set);

      assert(nodepq->validlowerbound);
      assert(nodepq->lowerbound < SCIP_INVALID);

      return nodepq->lowerbound;
   }
}

/** gets the node with minimal lower bound of all nodes in the queue */
SCIP_NODE* SCIPnodepqGetLowerboundNode(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->nodesel != NULL);
   assert(set != NULL);

   if( nodepq->nodesel->lowestboundfirst )
   {
      /* the node selector's compare method sorts the minimal lower bound to the front */
      if( nodepq->len > 0 )
      {
         assert(nodepq->slots[0] != NULL);
         return nodepq->slots[0];
      }
      else
         return NULL;
   }
   else
   {
      /* we use bookkeeping to remember the lowest bound */

      /* if the cached lower bound node is invalid, calculate it */
      if( !nodepq->validlowerbound || nodepq->lowerboundnode == NULL )
         nodepqCalcLowerbound(nodepq, set);
      
      assert(nodepq->validlowerbound);
      assert(nodepq->lowerbound < SCIP_INVALID);
      assert(nodepq->lowerbound == SCIPsetInfinity(set) || nodepq->lowerboundnode != NULL); /*lint !e777*/

      return nodepq->lowerboundnode;
   }
}

/** gets the sum of lower bounds of all nodes in the queue */
SCIP_Real SCIPnodepqGetLowerboundSum(
   SCIP_NODEPQ*          nodepq              /**< node priority queue */
   )
{
   assert(nodepq != NULL);

   return nodepq->lowerboundsum;
}

/** free all nodes from the queue that are cut off by the given upper bound */
SCIP_RETCODE SCIPnodepqBound(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   )
{
   SCIP_NODE* node;
   int pos;
   SCIP_Bool parentfelldown;

   assert(nodepq != NULL);

   SCIPdebugMessage("bounding node queue of length %d with cutoffbound=%g\n", nodepq->len, cutoffbound);
   pos = nodepq->len-1;
   while( pos >= 0 )
   {
      assert(pos < nodepq->len);
      node = nodepq->slots[pos];
      assert(node != NULL);
      assert(SCIPnodeGetType(node) == SCIP_NODETYPE_LEAF);
      if( SCIPsetIsGE(set, SCIPnodeGetLowerbound(node), cutoffbound) )
      {
         SCIPdebugMessage("free node in slot %d (len=%d) at depth %d with lowerbound=%g\n",
            pos, nodepq->len, SCIPnodeGetDepth(node), SCIPnodeGetLowerbound(node));

         /* cut off node; because we looped from back to front, the existing children of the node must have a smaller
          * lower bound than the cut off value
          */
         assert(PQ_LEFTCHILD(pos) >= nodepq->len
            || SCIPsetIsLT(set, SCIPnodeGetLowerbound(nodepq->slots[PQ_LEFTCHILD(pos)]), cutoffbound));
         assert(PQ_RIGHTCHILD(pos) >= nodepq->len
            || SCIPsetIsLT(set, SCIPnodeGetLowerbound(nodepq->slots[PQ_RIGHTCHILD(pos)]), cutoffbound));

         /* free the slot in the node PQ */
         parentfelldown = nodepqDelPos(nodepq, set, pos);

         /* - if the slot was occupied by the parent, we have to check this slot (the parent) again; unfortunately,
          *   we will check the node which occupied the parent's slot again, even though it cannot be cut off;
          * - otherwise, the slot was the last slot or it was occupied by a node with a position greater than
          *   the current position; this node was already checked and we can decrease the position
          */
         if( !parentfelldown )
            pos--;

         SCIPvbcCutoffNode(stat->vbc, stat, node);

         /* free memory of the node */
         SCIP_CALL( SCIPnodeFree(&node, blkmem, set, tree, lp) );
      }
      else
         pos--;
   }
   SCIPdebugMessage(" -> bounded node queue has length %d\n", nodepq->len);

   return SCIP_OKAY;
}




/*
 * node selector methods 
 */

/** method to call, when the standard mode priority of a node selector was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdNodeselStdPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetNodeselStdPriority() to mark the nodesels unsorted */
   SCIP_CALL( SCIPsetNodeselStdPriority(scip, (SCIP_NODESEL*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** method to call, when the memory saving mode priority of a node selector was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdNodeselMemsavePriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetNodeselMemsavePriority() to mark the nodesels unsorted */
   SCIP_CALL( SCIPsetNodeselMemsavePriority(scip, (SCIP_NODESEL*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a node selector */
SCIP_RETCODE SCIPnodeselCreate(
   SCIP_NODESEL**        nodesel,            /**< pointer to store node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of node selector */
   const char*           desc,               /**< description of node selector */
   int                   stdpriority,        /**< priority of the node selector in standard mode */
   int                   memsavepriority,    /**< priority of the node selector in memory saving mode */
   SCIP_Bool             lowestboundfirst,   /**< does node comparison sorts w.r.t. lower bound as primal criterion? */
   SCIP_DECL_NODESELFREE ((*nodeselfree)),   /**< destructor of node selector */
   SCIP_DECL_NODESELINIT ((*nodeselinit)),   /**< initialize node selector */
   SCIP_DECL_NODESELEXIT ((*nodeselexit)),   /**< deinitialize node selector */
   SCIP_DECL_NODESELINITSOL((*nodeselinitsol)),/**< solving process initialization method of node selector */
   SCIP_DECL_NODESELEXITSOL((*nodeselexitsol)),/**< solving process deinitialization method of node selector */
   SCIP_DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   SCIP_DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   SCIP_NODESELDATA*     nodeseldata         /**< node selector data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(nodesel != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(nodeselselect != NULL);
   assert(nodeselcomp != NULL);

   SCIP_ALLOC( BMSallocMemory(nodesel) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*nodesel)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*nodesel)->desc, desc, strlen(desc)+1) );
   (*nodesel)->stdpriority = stdpriority;
   (*nodesel)->memsavepriority = memsavepriority;
   (*nodesel)->nodeselfree = nodeselfree;
   (*nodesel)->nodeselinit = nodeselinit;
   (*nodesel)->nodeselexit = nodeselexit;
   (*nodesel)->nodeselinitsol = nodeselinitsol;
   (*nodesel)->nodeselexitsol = nodeselexitsol;
   (*nodesel)->nodeselselect = nodeselselect;
   (*nodesel)->nodeselcomp = nodeselcomp;
   (*nodesel)->nodeseldata = nodeseldata;
   (*nodesel)->lowestboundfirst = lowestboundfirst;
   (*nodesel)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "nodeselection/%s/stdpriority", name);
   sprintf(paramdesc, "priority of branching rule <%s> in standard mode", name);
   SCIP_CALL( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
                  &(*nodesel)->stdpriority, stdpriority, INT_MIN/4, INT_MAX/4, 
                  paramChgdNodeselStdPriority, (SCIP_PARAMDATA*)(*nodesel)) ); /*lint !e740*/

   sprintf(paramname, "nodeselection/%s/memsavepriority", name);
   sprintf(paramdesc, "priority of branching rule <%s> in memory saving mode", name);
   SCIP_CALL( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
                  &(*nodesel)->memsavepriority, memsavepriority, INT_MIN/4, INT_MAX/4, 
                  paramChgdNodeselMemsavePriority, (SCIP_PARAMDATA*)(*nodesel)) ); /*lint !e740*/

   return SCIP_OKAY;
}
   
/** frees memory of node selector */
SCIP_RETCODE SCIPnodeselFree(
   SCIP_NODESEL**        nodesel,            /**< pointer to node selector data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodesel != NULL);
   assert(*nodesel != NULL);
   assert(!(*nodesel)->initialized);
   assert(set != NULL);

   /* call destructor of node selector */
   if( (*nodesel)->nodeselfree != NULL )
   {
      SCIP_CALL( (*nodesel)->nodeselfree(set->scip, *nodesel) );
   }

   BMSfreeMemoryArray(&(*nodesel)->name);
   BMSfreeMemoryArray(&(*nodesel)->desc);
   BMSfreeMemory(nodesel);

   return SCIP_OKAY;
}

/** initializes node selector */
SCIP_RETCODE SCIPnodeselInit(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);

   if( nodesel->initialized )
   {
      SCIPerrorMessage("node selector <%s> already initialized", nodesel->name);
      return SCIP_INVALIDCALL;
   }

   if( nodesel->nodeselinit != NULL )
   {
      SCIP_CALL( nodesel->nodeselinit(set->scip, nodesel) );
   }
   nodesel->initialized = TRUE;

   return SCIP_OKAY;
}

/** deinitializes node selector */
SCIP_RETCODE SCIPnodeselExit(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);

   if( !nodesel->initialized )
   {
      SCIPerrorMessage("node selector <%s> not initialized", nodesel->name);
      return SCIP_INVALIDCALL;
   }

   if( nodesel->nodeselexit != NULL )
   {
      SCIP_CALL( nodesel->nodeselexit(set->scip, nodesel) );
   }
   nodesel->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs node selector that the branch and bound process is being started */
SCIP_RETCODE SCIPnodeselInitsol(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);

   /* call solving process initialization method of node selector */
   if( nodesel->nodeselinitsol != NULL )
   {
      SCIP_CALL( nodesel->nodeselinitsol(set->scip, nodesel) );
   }

   return SCIP_OKAY;
}

/** informs node selector that the branch and bound process data is being freed */
SCIP_RETCODE SCIPnodeselExitsol(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of node selector */
   if( nodesel->nodeselexitsol != NULL )
   {
      SCIP_CALL( nodesel->nodeselexitsol(set->scip, nodesel) );
   }

   return SCIP_OKAY;
}

/** select next node to be processed */
SCIP_RETCODE SCIPnodeselSelect(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE**           selnode             /**< pointer to store node to be processed next */
   )
{
   assert(nodesel != NULL);
   assert(nodesel->nodeselselect != NULL);
   assert(set != NULL);
   assert(selnode != NULL);

   SCIP_CALL( nodesel->nodeselselect(set->scip, nodesel, selnode) );

   return SCIP_OKAY;
}

/** compares two nodes; returns -1/0/+1 if node1 better/equal/worse than node2 */
int SCIPnodeselCompare(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node1,              /**< first node to compare */
   SCIP_NODE*            node2               /**< second node to compare */
   )
{
   assert(nodesel != NULL);
   assert(nodesel->nodeselcomp != NULL);
   assert(set != NULL);
   assert(node1 != NULL);
   assert(node2 != NULL);

   return nodesel->nodeselcomp(set->scip, nodesel, node1, node2);
}

/** gets name of node selector */
const char* SCIPnodeselGetName(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->name;
}

/** gets description of node selector */
const char* SCIPnodeselGetDesc(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->desc;
}

/** gets priority of node selector in standard mode */
int SCIPnodeselGetStdPriority(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->stdpriority;
}

/** gets priority of node selector in standard mode */
void SCIPnodeselSetStdPriority(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the node selector */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);

   nodesel->stdpriority = priority;
   set->nodesel = NULL;
}

/** gets priority of node selector in memory saving mode */
int SCIPnodeselGetMemsavePriority(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->memsavepriority;
}

/** sets priority of node selector in memory saving mode */
void SCIPnodeselSetMemsavePriority(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the node selector */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);
   
   nodesel->memsavepriority = priority;
   set->nodesel = NULL;
}

/** gets user data of node selector */
SCIP_NODESELDATA* SCIPnodeselGetData(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->nodeseldata;
}

/** sets user data of node selector; user has to free old data in advance! */
void SCIPnodeselSetData(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_NODESELDATA*     nodeseldata         /**< new node selector user data */
   )
{
   assert(nodesel != NULL);

   nodesel->nodeseldata = nodeseldata;
}

/** is node selector initialized? */
SCIP_Bool SCIPnodeselIsInitialized(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->initialized;
}

