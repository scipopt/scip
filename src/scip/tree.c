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

#include "def.h"
#include "constraint.h"
#include "sort.h"
#include "tree.h"


struct Leaf
{
   BASIS*           basis;              /**< pointer to LP basis information */
};

struct Fork
{
   unsigned int     nchildren:16;       /**< number of children of this parent node */
};

struct Node
{
   NODE*            parent;             /**< parent node in the tree */
   ROWLIST*         rowlist;            /**< active LP rows */
   CONSLIST*        conslist;           /**< active constraints */
   ROWLIST*         insertedRows;       /**< rows inserted at this node into the LP */
   ROWLIST*         removedRows;        /**< rows removed at this node from the LP */   
   DOMAINCHG*       domainChg;          /**< list of domain changes */
   union
   {
      LEAF*         leaf;               /**< data for leaf nodes */
      FORK*         fork;               /**< data for fork nodes */
   } data;
   double           lowerbound;         /**< lower (dual) LP bound of subtree */
   unsigned int     depth:16;           /**< depth in the tree */
   unsigned int     nodetype:1;         /**< type of node */
};

struct Tree
{
   NODE*            root;               /**< root node of the tree */
   PQUEUE*          leaves;             /**< leaves of the tree */
};


NODE* SCIPcreateLeaf(                   /**< creates a leaf node */
   NODE*            parent,             /**< parent node in the tree */
   BASIS*           basis               /**< pointer to LP basis information */
   )
{
   NODE* node;

   CHECK_NULL( allocMemory(node) );
   node->parent = parent;
   node->rowlist = NULL;
   node->conslist = NULL;
   node->insertedRows = NULL;
   node->removedRows = NULL;
   node->domainChg = NULL;
   CHECK_NULL( allocMemory(node->data.leaf) );
   node->data.leaf->basis = basis;
   SCIPuseBasis(basis);
   node->nodetype = SCIP_NODETYPE_LEAF;
   if( parent == NULL )
   {
      node->depth = 0;
      node->lowerbound = -SCIP_INFINITY;
   }
   else
   {
      node->depth = parent->depth+1;
      node->lowerbound = parent->lowerbound;
   }

   return node;
}
