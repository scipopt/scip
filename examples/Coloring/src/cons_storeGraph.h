/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_storeGraph.h
 * @brief  constraint handler for storing the graph at each node of the tree
 * @author Gerald Gamrath
 */

#ifndef CONSSTOREGRAPH_H
#define CONSSTOREGRAPH_H

#include "scip/scip.h"
#include "tclique/tclique.h"

#ifdef __cplusplus
extern "C" {
#endif

/* type of storeGraph constraint: differ, same or root */
enum COLOR_ConsType
{
   COLOR_CONSTYPE_DIFFER = 0,  /* constraint representing the branching decision differ(i,j) */
   COLOR_CONSTYPE_SAME   = 1,  /* constraint representing the branching decision same(i,j) */
   COLOR_CONSTYPE_ROOT   = 2   /* constraint created for the root, is created automatically */
};
typedef enum COLOR_ConsType COLOR_CONSTYPE;

/** returns the store graph constraint of the current node, needs the pointer to the constraint handler */
extern
SCIP_CONS* COLORconsGetActiveStoreGraphConsFromHandler(
   SCIP_CONSHDLR*        conshdlr            /**< constaint handler for store-graph constraints */
   );


/** returns the store graph constraint of the current node, needs only the pointer to scip */
extern
SCIP_CONS* COLORconsGetActiveStoreGraphCons(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** returns array of representatives of all nodes */
extern
int* COLORconsGetRepresentatives(
   SCIP*                 scip                 /**< SCIP data structure */
   );


/** returns the current graph */
extern
TCLIQUE_GRAPH* COLORconsGetCurrentGraph(
   SCIP*                 scip            /**< SCIP data structure */
   );


/** returns the complementary graph */
TCLIQUE_GRAPH* COLORconsGetComplementaryGraph(
   SCIP*                 scip            /**< SCIP data structure */
   );


/** returns representative of the union which contains a given node */
extern
int COLORconsGetRepresentative(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node                /**< the node, for wich the representative is searched */
   );


/** returns array of all unions, a union is saved in the array at the position of its representative */
extern
void COLORconsGetUnions(
   SCIP*                 scip,               /**< SCIP data structure */
   int***                unions,             /**< output: array containing array which contains nodes in the union */
   int**                 lengths             /**< output: lengths of the unions */
   );


/** returns the union which has a given node as representative */
extern
void COLORconsGetUnion(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 unionSet,           /**< output: array containig nodes in the union */
   int*                  length,             /**< output: length of the union */
   int                   node                /**< the node, whose union we want to get */
   );


/** creates the handler for graph storing constraints and includes it in SCIP */
extern
SCIP_RETCODE COLORincludeConshdlrStoreGraph(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a storeGraph constraint, uses knowledge of the B&B-father*/
extern
SCIP_RETCODE COLORcreateConsStoreGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONS*            fatherconstraint,   /**< constraint in B&B-father */
   int                   type,               /**< type of the constraint: ROOT for root-constraint, else SAME or DIFFER */
   int                   node1,              /**< the first node of the constraint or -1 if root-constraint */
   int                   node2,              /**< the second node of the constraint or -1 if root-constraint */
   SCIP_NODE*            stickingnode        /**< the B&B-tree node at which the constraint will be sticking */     
   );


/** returns the stack and the number of elements on it */
extern
void COLORconsGetStack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   );

#ifdef __cplusplus
}
#endif

#endif
