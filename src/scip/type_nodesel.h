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
#pragma ident "@(#) $Id: type_nodesel.h,v 1.9 2005/05/31 17:20:25 bzfpfend Exp $"

/**@file   type_nodesel.h
 * @brief  type definitions for node selectors
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_NODESEL_H__
#define __TYPE_NODESEL_H__


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

/** initialization method of node selector (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define DECL_NODESELINIT(x) RETCODE x (SCIP* scip, NODESEL* nodesel)

/** deinitialization method of node selector (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define DECL_NODESELEXIT(x) RETCODE x (SCIP* scip, NODESEL* nodesel)

/** solving process initialization method of node selector (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The node selector may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define DECL_NODESELINITSOL(x) RETCODE x (SCIP* scip, NODESEL* nodesel)

/** solving process deinitialization method of node selector (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The node selector should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define DECL_NODESELEXITSOL(x) RETCODE x (SCIP* scip, NODESEL* nodesel)

/** node selection method of node selector
 *
 *  This method is called to select the next leaf of the branch and bound tree to be processed.
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
 *  This method is called to compare two nodes regarding their order in the node priority queue.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 *  - node1           : first node to compare
 *  - node2           : second node to compare
 *
 *  possible return values:
 *  - value < 0: node1 comes before (is better than) node2
 *  - value = 0: both nodes are equally good
 *  - value > 0: node2 comes after (is worse than) node2
 */
#define DECL_NODESELCOMP(x) int x (SCIP* scip, NODESEL* nodesel, NODE* node1, NODE* node2)



#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_tree.h"
#include "scip/type_scip.h"


#endif
