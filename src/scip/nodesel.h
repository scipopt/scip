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

/**@file   nodesel.h
 * @brief  datastructures and methods for node selectors
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __NODESEL_H__
#define __NODESEL_H__


typedef struct NodePQ NODEPQ;           /**< node priority queue */
typedef struct Nodesel NODESEL;         /**< node selector data structure */
typedef struct NodeselData NODESELDATA; /**< node selector specific data */


/** initialization method of node selector
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_NODESELINIT(x) RETCODE x (NODESEL* nodesel, SCIP* scip)

/** deinitialization method of node selector
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_NODESELEXIT(x) RETCODE x (NODESEL* nodesel, SCIP* scip)

/** node selection method of node selector
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 *  possible return values for *selnode:
 *    NULL    : problem is solved, because tree is empty
 *    non-NULL: node to be solved next
 */
#define DECL_NODESELSLCT(x) RETCODE x (NODESEL* nodesel, SCIP* scip, NODE** selnode)

/** node comparison method of node selector
 *  possible return values:
 *    < 0: node1 comes before (is better than) node2
 *    = 0: both nodes have the same value
 *    > 0: node2 comes after (is worse than) node2
 */
#define DECL_NODESELCOMP(x) int x (NODESEL* nodesel, SCIP* scip, NODE* node1, NODE* node2)




#include "scip.h"
#include "retcode.h"
#include "set.h"
#include "tree.h"
#include "lp.h"


extern
RETCODE SCIPnodepqCreate(               /**< creates node priority queue */
   NODEPQ**         nodepq              /**< pointer to a node priority queue */
   );

extern
void SCIPnodepqDestroy(                 /**< frees node priority queue, but not the data nodes themselves */
   NODEPQ**         nodepq              /**< pointer to a node priority queue */
   );

extern
RETCODE SCIPnodepqFree(                 /**< frees node priority queue and all nodes in the queue */
   NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPnodepqInsert(               /**< inserts node into node priority queue */
   NODEPQ*          nodepq,             /**< pointer to a node priority queue */
   const SET*       set,                /**< global SCIP settings */
   NODE*            node                /**< node to be inserted */
   );

extern
NODE* SCIPnodepqRemove(                 /**< removes and returns best node from the node priority queue */
   NODEPQ*          nodepq,             /**< pointer to a node priority queue */
   const SET*       set                 /**< global SCIP settings */
   );

extern
NODE* SCIPnodepqFirst(                  /**< returns the best node of the queue without removing it */
   const NODEPQ*    nodepq              /**< pointer to a node priority queue */
   );

extern
int SCIPnodepqLen(                      /**< returns the number of nodes stored in the node priority queue */
   const NODEPQ*    nodepq              /**< pointer to a node priority queue */
   );

extern
Real SCIPnodepqGetLowerbound(           /**< gets the minimal lower bound of all nodes in the queue */
   NODEPQ*          nodepq,             /**< pointer to a node priority queue */
   const SET*       set                 /**< global SCIP settings */
   );

extern
Real SCIPnodepqGetLowerboundSum(        /**< gets the sum of lower bounds of all nodes in the queue */
   NODEPQ*          nodepq              /**< pointer to a node priority queue */
   );

extern
RETCODE SCIPnodepqBound(                /**< free all nodes from the queue that are cut off by the given upper bound */
   NODEPQ*          nodepq,             /**< pointer to a node priority queue */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   Real             upperbound          /**< upper bound: all nodes with lowerbound >= upperbound are cut off */
   );

extern
RETCODE SCIPnodeselCreate(              /**< creates a node selector */
   NODESEL**        nodesel,            /**< pointer to store node selector */
   const char*      name,               /**< name of node selector */
   const char*      desc,               /**< description of node selector */
   DECL_NODESELINIT((*nodeselinit)),    /**< initialise node selector */
   DECL_NODESELEXIT((*nodeselexit)),    /**< deinitialise node selector */
   DECL_NODESELSLCT((*nodeselslct)),    /**< node selection method */
   DECL_NODESELCOMP((*nodeselcomp)),    /**< node comparison method */
   NODESELDATA*     nodeseldata,        /**< node selector data */
   Bool             lowestboundfirst    /**< does node comparison sorts w.r.t. lower bound as primal criterion? */
   );

extern
RETCODE SCIPnodeselFree(                /**< frees memory of node selector */
   NODESEL**        nodesel             /**< pointer to node selector data structure */
   );

extern
RETCODE SCIPnodeselInit(                /**< initializes node selector */
   NODESEL*         nodesel,            /**< node selector */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPnodeselExit(                /**< deinitializes node selector */
   NODESEL*         nodesel,            /**< node selector */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPnodeselSelect(              /**< select next node to be processed */
   NODESEL*         nodesel,            /**< node selector */
   SCIP*            scip,               /**< SCIP data structure */   
   NODE**           selnode             /**< pointer to store node to be processed next */
   );

extern
int SCIPnodeselCompare(                 /**< compares two nodes; returns -1/0/+1 if node1 better/equal/worse than node2 */
   NODESEL*         nodesel,            /**< node selector */
   SCIP*            scip,               /**< SCIP data structure */   
   NODE*            node1,              /**< first node to compare */
   NODE*            node2               /**< second node to compare */
   );

extern
const char* SCIPnodeselGetName(         /**< gets name of node selector */
   NODESEL*         nodesel             /**< node selector */
   );

extern
Bool SCIPnodeselIsInitialized(          /**< is node selector initialized? */
   NODESEL*         nodesel             /**< node selector */
   );


#endif
