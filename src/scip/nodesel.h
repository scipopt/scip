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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nodesel.h,v 1.35 2005/01/21 09:16:57 bzfpfend Exp $"

/**@file   nodesel.h
 * @brief  internal methods for node selectors and node priority queues
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __NODESEL_H__
#define __NODESEL_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_lp.h"
#include "type_tree.h"
#include "type_scip.h"
#include "pub_nodesel.h"



/* 
 * node priority queue methods
 */

/** creates node priority queue */
extern
RETCODE SCIPnodepqCreate(
   NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   SET*             set,                /**< global SCIP settings */
   NODESEL*         nodesel             /**< node selector to use for sorting the nodes in the queue */
   );

/** frees node priority queue, but not the data nodes themselves */
extern
void SCIPnodepqDestroy(
   NODEPQ**         nodepq              /**< pointer to a node priority queue */
   );

/** frees node priority queue and all nodes in the queue */
extern
RETCODE SCIPnodepqFree(
   NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   );

/** deletes all nodes in the node priority queue */
extern
RETCODE SCIPnodepqClear(
   NODEPQ*          nodepq,             /**< node priority queue */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   );

/** returns the node selector associated with the given node priority queue */
extern
NODESEL* SCIPnodepqGetNodesel(
   NODEPQ*          nodepq              /**< node priority queue */
   );

/** sets the node selector used for sorting the nodes in the queue, and resorts the queue if necessary */
extern
RETCODE SCIPnodepqSetNodesel(
   NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   SET*             set,                /**< global SCIP settings */
   NODESEL*         nodesel             /**< node selector to use for sorting the nodes in the queue */
   );

/** compares two nodes; returns -1/0/+1 if node1 better/equal/worse than node2 */
extern
int SCIPnodepqCompare(
   NODEPQ*          nodepq,             /**< node priority queue */
   SET*             set,                /**< global SCIP settings */
   NODE*            node1,              /**< first node to compare */
   NODE*            node2               /**< second node to compare */
   );

/** inserts node into node priority queue */
extern
RETCODE SCIPnodepqInsert(
   NODEPQ*          nodepq,             /**< node priority queue */
   SET*             set,                /**< global SCIP settings */
   NODE*            node                /**< node to be inserted */
   );

/** removes node from the node priority queue */
extern
RETCODE SCIPnodepqRemove(
   NODEPQ*          nodepq,             /**< node priority queue */
   SET*             set,                /**< global SCIP settings */
   NODE*            node                /**< node to remove */
   );

/** returns the best node of the queue without removing it */
extern
NODE* SCIPnodepqFirst(
   const NODEPQ*    nodepq              /**< node priority queue */
   );

/** returns the nodes array of the queue */
extern
NODE** SCIPnodepqNodes(
   const NODEPQ*    nodepq              /**< node priority queue */
   );

/** returns the number of nodes stored in the node priority queue */
extern
int SCIPnodepqLen(
   const NODEPQ*    nodepq              /**< node priority queue */
   );

/** gets the minimal lower bound of all nodes in the queue */
extern
Real SCIPnodepqGetLowerbound(
   NODEPQ*          nodepq,             /**< node priority queue */
   SET*             set                 /**< global SCIP settings */
   );

/** gets the node with minimal lower bound of all nodes in the queue */
extern
NODE* SCIPnodepqGetLowerboundNode(
   NODEPQ*          nodepq,             /**< node priority queue */
   SET*             set                 /**< global SCIP settings */
   );

/** gets the sum of lower bounds of all nodes in the queue */
extern
Real SCIPnodepqGetLowerboundSum(
   NODEPQ*          nodepq              /**< node priority queue */
   );

/** free all nodes from the queue that are cut off by the given upper bound */
extern
RETCODE SCIPnodepqBound(
   NODEPQ*          nodepq,             /**< node priority queue */
   MEMHDR*          memhdr,             /**< block memory buffer */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   );




/*
 * node selector methods 
 */

/** creates a node selector */
extern
RETCODE SCIPnodeselCreate(
   NODESEL**        nodesel,            /**< pointer to store node selector */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of node selector */
   const char*      desc,               /**< description of node selector */
   int              stdpriority,        /**< priority of the node selector in standard mode */
   int              memsavepriority,    /**< priority of the node selector in memory saving mode */
   Bool             lowestboundfirst,   /**< does node comparison sorts w.r.t. lower bound as primal criterion? */
   DECL_NODESELFREE ((*nodeselfree)),   /**< destructor of node selector */
   DECL_NODESELINIT ((*nodeselinit)),   /**< initialize node selector */
   DECL_NODESELEXIT ((*nodeselexit)),   /**< deinitialize node selector */
   DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   NODESELDATA*     nodeseldata         /**< node selector data */
   );

/** frees memory of node selector */
extern
RETCODE SCIPnodeselFree(
   NODESEL**        nodesel,            /**< pointer to node selector data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** initializes node selector */
extern
RETCODE SCIPnodeselInit(
   NODESEL*         nodesel,            /**< node selector */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** deinitializes node selector */
extern
RETCODE SCIPnodeselExit(
   NODESEL*         nodesel,            /**< node selector */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** select next node to be processed */
extern
RETCODE SCIPnodeselSelect(
   NODESEL*         nodesel,            /**< node selector */
   SET*             set,                /**< global SCIP settings */
   NODE**           selnode             /**< pointer to store node to be processed next */
   );

/** compares two nodes; returns -1/0/+1 if node1 better/equal/worse than node2 */
extern
int SCIPnodeselCompare(
   NODESEL*         nodesel,            /**< node selector */
   SET*             set,                /**< global SCIP settings */
   NODE*            node1,              /**< first node to compare */
   NODE*            node2               /**< second node to compare */
   );

/** sets priority of node selector in standard mode */
extern
void SCIPnodeselSetStdPriority(
   NODESEL*         nodesel,            /**< node selector */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the node selector */
   );

/** sets priority of node selector in memory saving mode */
extern
void SCIPnodeselSetMemsavePriority(
   NODESEL*         nodesel,            /**< node selector */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the node selector */
   );


#endif
