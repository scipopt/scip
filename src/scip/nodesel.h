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

/**@file   nodesel.h
 * @brief  methods and datastructures for node selectors
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __NODESEL_H__
#define __NODESEL_H__


typedef struct NodePQ NODEPQ;           /**< node priority queue */
typedef struct Nodesel NODESEL;         /**< node selector data structure */
typedef struct NodeselData NODESELDATA; /**< node selector specific data */


/** destructor of node selector to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define DECL_NODESELFREE(x) RETCODE x (SCIP* scip, NODESEL* nodesel)

/** initialization method of node selector (called when problem solving starts)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define DECL_NODESELINIT(x) RETCODE x (SCIP* scip, NODESEL* nodesel)

/** deinitialization method of node selector (called when problem solving exits)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define DECL_NODESELEXIT(x) RETCODE x (SCIP* scip, NODESEL* nodesel)

/** node selection method of node selector
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 *  - selnode         : pointer to store the selected node
 *
 *  possible return values for *selnode:
 *  - NULL    : problem is solved, because tree is empty
 *  - non-NULL: node to be solved next
 */
#define DECL_NODESELSELECT(x) RETCODE x (SCIP* scip, NODESEL* nodesel, NODE** selnode)

/** node comparison method of node selector
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 *
 *  possible return values:
 *  - value < 0: node1 comes before (is better than) node2
 *  - value = 0: both nodes are equally good
 *  - value > 0: node2 comes after (is worse than) node2
 */
#define DECL_NODESELCOMP(x) int x (SCIP* scip, NODESEL* nodesel, NODE* node1, NODE* node2)




#include "scip.h"
#include "retcode.h"
#include "set.h"
#include "tree.h"
#include "lp.h"



/* 
 * node priority queue methods
 */

/** creates node priority queue */
extern
RETCODE SCIPnodepqCreate(
   NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   const SET*       set                 /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   );

/** inserts node into node priority queue */
extern
RETCODE SCIPnodepqInsert(
   NODEPQ*          nodepq,             /**< node priority queue */
   const SET*       set,                /**< global SCIP settings */
   NODE*            node                /**< node to be inserted */
   );

/** removes node from the node priority queue */
extern
RETCODE SCIPnodepqRemove(
   NODEPQ*          nodepq,             /**< node priority queue */
   const SET*       set,                /**< global SCIP settings */
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
   const SET*       set                 /**< global SCIP settings */
   );

/** gets the node with minimal lower bound of all nodes in the queue */
extern
NODE* SCIPnodepqGetLowerboundNode(
   NODEPQ*          nodepq,             /**< node priority queue */
   const SET*       set                 /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   Real             upperbound          /**< upper bound: all nodes with lowerbound >= upperbound are cut off */
   );

/** resorts the priority queue (necessary for changes in node selector) */
extern
RETCODE SCIPnodepqResort(
   NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );




/*
 * node selector methods 
 */

/** creates a node selector */
extern
RETCODE SCIPnodeselCreate(
   NODESEL**        nodesel,            /**< pointer to store node selector */
   const char*      name,               /**< name of node selector */
   const char*      desc,               /**< description of node selector */
   DECL_NODESELFREE ((*nodeselfree)),   /**< destructor of node selector */
   DECL_NODESELINIT ((*nodeselinit)),   /**< initialize node selector */
   DECL_NODESELEXIT ((*nodeselexit)),   /**< deinitialize node selector */
   DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   NODESELDATA*     nodeseldata,        /**< node selector data */
   Bool             lowestboundfirst    /**< does node comparison sorts w.r.t. lower bound as primal criterion? */
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
   const SET*       set,                /**< global SCIP settings */
   NODE**           selnode             /**< pointer to store node to be processed next */
   );

/** compares two nodes; returns -1/0/+1 if node1 better/equal/worse than node2 */
extern
int SCIPnodeselCompare(
   NODESEL*         nodesel,            /**< node selector */
   const SET*       set,                /**< global SCIP settings */
   NODE*            node1,              /**< first node to compare */
   NODE*            node2               /**< second node to compare */
   );

/** gets name of node selector */
extern
const char* SCIPnodeselGetName(
   NODESEL*         nodesel             /**< node selector */
   );

/** gets user data of node selector */
extern
NODESELDATA* SCIPnodeselGetData(
   NODESEL*         nodesel             /**< node selector */
   );

/** sets user data of node selector; user has to free old data in advance! */
extern
void SCIPnodeselSetData(
   NODESEL*         nodesel,            /**< node selector */
   NODESELDATA*     nodeseldata         /**< new node selector user data */
   );

/** is node selector initialized? */
extern
Bool SCIPnodeselIsInitialized(
   NODESEL*         nodesel             /**< node selector */
   );


#endif
