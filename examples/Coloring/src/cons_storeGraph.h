/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
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
 *
 * This file implements the constraints that are used for the branching in the coloring algorithm.
 *
 * For each node in the branch-and-bound tree, a constraint of this type is created, which stores
 * all restrictions related to that branch-and-bound node.
 *
 * First of all, it stores the type of the constraint ("same" or "differ", the root has type root)
 * and the two nodes in the graph on which this restriction is applied.  When the branch-and-bound
 * node corresponding to the constraint is examined for the first time, the constraint creates a
 * graph that takes into account all the restrictions, which are active at this node.
 * At the root, this is the original (preprocessed) graph.  At any other branch-and-bound node, it
 * takes the graph of the constraint related to the branch-and-bound father of the current node and
 * modifies it so that all restrictions up to this node are respected.  Since the graph in the
 * branch-and-bound father respects all restrictions on the path to that node, only the last
 * requirement, the one saved at the current branch-and-bound node, must be added.
 * This is done as follows: Adding a DIFFER(v,w) constraint is easy, since it suffices to add
 * an edge between v and w. For a SAME(v,w) constraint, the original idea is to collapse the nodes v
 * and w into one single vertex. Since this is not possible in the tclique-graph data structure, we
 * introduce new edges in the graph, so that v and w have the same neighborhood.  Hence, in the
 * pricing routine, each new stable set will either contain both nodes or none of them, since we
 * create (inclusion-) maximal sets.
 *
 * This does of course not hold for sets created in a higher level of the branch-and-bound tree or
 * in another subtree. In order to forbid all of these sets, which do not fulfill the current
 * restrictions, a propagation is started when the node is entered the first time and repeated
 * later, if the node is reentered after the creation of new variables in another subtree. The
 * propagation simply fixes to 0 all variables representing a stable set that does not
 * fulfill the restriction at the current node.
 *
 * The information about all fusions of nodes (caused by the SAME() operation) is stored, so that the nodes
 * constituting a union can be accessed easily. Each union has a representative and a set of nodes, whereas
 * each node knows the representative of the union it belongs to. At the beginning, each node forms its own
 * union and therefore each node also represents this union, consisting of only this node.  Later on, some
 * nodes represent unions of several nodes, while other nodes are part of a union which they do not represent,
 * so they have another node as representative. The representatives of the nodes are returned by the methods
 * COLORconsGetRepresentative() / COLORconsGetRepresentatives(), the union represented by a node is returned
 * by COLORconsGetUnion(), the array of unions, indexed by the representing node, is returned by
 * COLORconsGetUnions().
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
