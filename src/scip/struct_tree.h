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
#pragma ident "@(#) $Id: struct_tree.h,v 1.9 2004/04/27 15:50:05 bzfpfend Exp $"

/**@file   struct_tree.h
 * @brief  datastructures for branch and bound tree
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_TREE_H__
#define __STRUCT_TREE_H__


#include "def.h"
#include "type_lpi.h"
#include "type_lp.h"
#include "type_var.h"
#include "type_tree.h"
#include "type_cons.h"
#include "type_nodesel.h"



/** child information (should not exceed the size of a pointer) */
struct Child
{
   int              arraypos;           /**< position of node in the children array */
};

/** sibling information (should not exceed the size of a pointer) */
struct Sibling
{
   int              arraypos;           /**< position of node in the siblings array */
};

/** leaf information (should not exceed the size of a pointer) */
struct Leaf
{
   NODE*            lpfork;             /**< fork/subroot node defining the LP state of the leaf */
};

/** fork without LP solution, where only bounds and constraints have been changed */
struct Junction
{
   int              nchildren;          /**< number of children of this parent node */
};

/** fork with solved LP, where bounds and constraints have been changed, and rows and columns were added */
struct Fork
{
   COL**            addedcols;          /**< array with pointers to new columns added at this node into the LP */
   ROW**            addedrows;          /**< array with pointers to new rows added at this node into the LP */
   LPISTATE*        lpistate;           /**< LP state information */
   int              naddedcols;         /**< number of columns added at this node */
   int              naddedrows;         /**< number of rows added at this node */
   int              nchildren;          /**< number of children of this parent node */
   int              nlpistateref;       /**< number of times, the LP state is needed */   
};

/** fork with solved LP, where bounds and constraints have been changed, and rows and columns were removed and added */
struct Subroot
{
   COL**            cols;               /**< array with pointers to the columns in the same order as in the LP */
   ROW**            rows;               /**< array with pointers to the rows in the same order as in the LP */
   LPISTATE*        lpistate;           /**< LP state information */
   int              ncols;              /**< number of columns in the LP */
   int              nrows;              /**< number of rows in the LP */
   int              nchildren;          /**< number of children of this parent node */
   int              nlpistateref;       /**< number of times, the LP state is needed */   
};

/** node data structure */
struct Node
{
   Real             lowerbound;         /**< lower (dual) LP bound of subtree */
   union
   {
      SIBLING       sibling;            /**< data for sibling nodes */
      CHILD         child;              /**< data for child nodes */
      LEAF          leaf;               /**< data for leaf nodes */
      JUNCTION      junction;           /**< data for junction nodes */
      FORK*         fork;               /**< data for fork nodes */
      SUBROOT*      subroot;            /**< data for subroot nodes */
   } data;
   NODE*            parent;             /**< parent node in the tree */
   CONSSETCHG*      conssetchg;         /**< constraint set changes at this node or NULL */
   DOMCHG*          domchg;             /**< domain changes at this node or NULL */
   unsigned int     depth:16;           /**< depth in the tree */
   unsigned int     nodetype:3;         /**< type of node */
   unsigned int     active:1;           /**< is node in the path to the current active node? */
};

/** branch and bound tree */
struct Tree
{
   NODE*            root;               /**< root node of the tree */
   NODEPQ*          leaves;             /**< leaves of the tree */
   NODE**           path;               /**< array of fork/subtree nodes storing the active path from root to leaf */
   NODE*            actnode;            /**< currently active node */
   NODE*            actlpfork;          /**< fork/subroot node defining the LP state of the active node */
   NODE*            actsubroot;         /**< root of the active subtree */
   NODE**           children;           /**< array with children of the active node */
   NODE**           siblings;           /**< array with siblings of the active node */
   int*             pathnlpcols;        /**< array with number of LP columns for each problem in active path */
   int*             pathnlprows;        /**< array with number of LP rows for each problem in active path */
   int              actlpforklpcount;   /**< LP number of last solved LP in current LP fork, or -1 if unknown */
   int              childrensize;       /**< available slots in children vector */
   int              nchildren;          /**< current number of children (number of used slots in children vector) */
   int              siblingssize;       /**< available slots in siblings vector */
   int              nsiblings;          /**< current number of siblings (number of used slots in siblings vector) */
   int              pathlen;            /**< length of the current path (== depth of the current node + 1) */
   int              pathsize;           /**< number of available slots in path arrays */
   int              correctlpdepth;     /**< depth to which current LP data corresponds to LP data of active path */
   Bool             actnodehaslp;       /**< is LP being processed in the currently active node? */
   Bool             cutoffdelayed;      /**< the treeCutoff() call was delayed because of diving and has to be executed */
};


#endif
