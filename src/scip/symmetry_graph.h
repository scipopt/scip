/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   symmetry_graph.h
 * @ingroup PUBLICCOREAPI
 * @brief  methods for dealing with symmetry detection graphs
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SYMMETRY_GRAPH_H__
#define __SCIP_SYMMETRY_GRAPH_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"
#include <symmetry/type_symmetry.h>

#ifdef __cplusplus
extern "C" {
#endif


/**@addtogroup PublicSymmetryGraphMethods
 *
 * @{
 */

/** creates and initializes a symmetry detection graph with memory for the given number of nodes and edges
 *
 *  @note at some point, the graph needs to be freed!
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateSymgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SYMTYPE           symtype,            /**< type of symmetries encoded in graph */
   SYM_GRAPH**           graph,              /**< pointer to hold symmetry detection graph */
   SCIP_VAR**            symvars,            /**< variables used in symmetry detection */
   int                   nsymvars,           /**< number of variables used in symmetry detection */
   int                   nopnodes,           /**< number of operator nodes */
   int                   nvalnodes,          /**< number of value nodes */
   int                   nconsnodes,         /**< number of constraint nodes */
   int                   nedges              /**< number of edges */
   );

/** frees a symmetry detection graph */
SCIP_EXPORT
SCIP_RETCODE SCIPfreeSymgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH**           graph               /**< pointer to symmetry detection graph */
   );

/** copies an existing graph and changes variable nodes according to a permutation */
SCIP_EXPORT
SCIP_RETCODE SCIPcopySymgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH**           graph,              /**< pointer to hold copy of symmetry detection graph */
   SYM_GRAPH*            origgraph,          /**< graph to be copied */
   int*                  perm,               /**< permutation of variables */
   SYM_SPEC              fixedtype           /**< variable types that must be fixed by symmetries */
   );

/** adds a symmetry detection graph for a linear constraint to existing graph
 *
 *  For permutation symmetries, a constraint node is added that is connected to all
 *  variable nodes in the constraint. Edges are colored according to the variable coefficients.
 *  For signed permutation symmetries, also edges connecting the constraint node and
 *  the negated variable nodes are added, these edges are colored by the negative coefficients.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPextendPermsymDetectionGraphLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_VAR**            vars,               /**< variable array of linear constraint */
   SCIP_Real*            vals,               /**< coefficients of linear constraint */
   int                   nvars,              /**< number of variables in linear constraint */
   SCIP_CONS*            cons,               /**< constraint for which we encode symmetries */
   SCIP_Real             lhs,                /**< left-hand side of constraint */
   SCIP_Real             rhs,                /**< right-hand side of constraint */
   SCIP_Bool*            success             /**< pointer to store whether graph could be built */
   );

/** adds nodes and edges corresponding to the aggregation of a variable to a symmetry detection graph
 *
 *  For permutation symmetries, the root node is connected with all variable nodes in the aggregation.
 *  Edges are colored according to the variable coefficients.
 *  For signed permutation symmetries, also edges connecting the root node and the negated variable
 *  nodes are added, these edges are colored by the negative coefficients.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddSymgraphVarAggregation(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   rootidx,            /**< index of root node of the aggregation */
   SCIP_VAR**            vars,               /**< array of variables in aggregation */
   SCIP_Real*            vals,               /**< coefficients of variables */
   int                   nvars,              /**< number of variables in aggregation */
   SCIP_Real             constant            /**< constant of aggregation */
   );

/*
 * methods for adding and accessing nodes
 */

/** adds an operator node to a symmetry detection graph and returns its node index */
SCIP_EXPORT
SCIP_RETCODE SCIPaddSymgraphOpnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   op,                 /**< int associated with operator of node */
   int*                  nodeidx             /**< pointer to hold index of created node */
   );

/** adds a value node to a symmetry detection graph and returns its node index */
SCIP_EXPORT
SCIP_RETCODE SCIPaddSymgraphValnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_Real             val,                /**< value of node */
   int*                  nodeidx             /**< pointer to hold index of created node */
   );

/** adds a constraint node to a symmetry detection graph and returns its node index */
SCIP_EXPORT
SCIP_RETCODE SCIPaddSymgraphConsnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_CONS*            cons,               /**< constraint of node */
   SCIP_Real             lhs,                /**< left-hand side of node */
   SCIP_Real             rhs,                /**< right-hand side of node */
   int*                  nodeidx             /**< pointer to hold index of created node */
   );

/** returns the (hypothetical) node index of a variable */
SCIP_EXPORT
int SCIPgetSymgraphVarnodeidx(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_VAR*             var                 /**< variable */
   );

/** returns the (hypothetical) node index of a negated variable */
SCIP_EXPORT
int SCIPgetSymgraphNegatedVarnodeidx(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_VAR*             var                 /**< variable */
   );

/** updates the lhs of a constraint node */
SCIP_EXPORT
SCIP_RETCODE SCIPupdateSymgraphLhs(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   nodeidx,            /**< index of constraint node in graph */
   SCIP_Real             newlhs              /**< new left-hand side of node */
   );

/** updates the rhs of a constraint node */
SCIP_EXPORT
SCIP_RETCODE SCIPupdateSymgraphRhs(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   nodeidx,            /**< index of constraint node in graph */
   SCIP_Real             newrhs              /**< new reft-hand side of node */
   );

/** registers a variable node (corresponding to active variable) to be fixed by symmetry */
SCIP_EXPORT
SCIP_RETCODE SCIPfixSymgraphVarnode(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_VAR*             var                 /**< active variable that needs to be fixed */
   );

/*
 * methods for adding edges
 */

/** adds an edge to a symmetry detection graph */
SCIP_EXPORT
SCIP_RETCODE SCIPaddSymgraphEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   first,              /**< first node index of edge */
   int                   second,             /**< second node index of edge */
   SCIP_Bool             hasval,             /**< whether the edge has a value */
   SCIP_Real             val                 /**< value of the edge (is it has a value) */
   );

/*
 * methods to compute colors
 */

/** computes colors of nodes and edges */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeSymgraphColors(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SYM_SPEC              fixedtype           /**< variable types that must be fixed by symmetries */
   );

/*
 * general methods
 */

/** returns the type of symmetries encoded in graph */
SCIP_EXPORT
SYM_SYMTYPE SCIPgetSymgraphSymtype(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   );

/** returns the variables in a symmetry detection graph */
SCIP_EXPORT
SCIP_VAR** SCIPgetSymgraphVars(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   );

/** returns the number of variables in a symmetry detection graph */
SCIP_EXPORT
int SCIPgetSymgraphNVars(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   );

/** returns the number of constraint nodes in a symmetry detection graph */
SCIP_EXPORT
int SCIPgetSymgraphNConsnodes(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   );

/** returns the number of non-variable nodes in a graph */
SCIP_EXPORT
int SCIPgetSymgraphNNodes(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   );

/** returns the number of edges in a graph */
SCIP_EXPORT
int SCIPgetSymgraphNEdges(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   );

/** return the first node index of an edge */
SCIP_EXPORT
int SCIPgetSymgraphEdgeFirst(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   edgeidx             /**< index of edge */
   );

/** return the second node index of an edge */
SCIP_EXPORT
int SCIPgetSymgraphEdgeSecond(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   edgeidx             /**< index of edge */
   );

/** returns the color of a variable node */
SCIP_EXPORT
int SCIPgetSymgraphVarnodeColor(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   nodeidx             /**< index of variable node */
   );

/** returns the type of a node */
SCIP_EXPORT
SYM_NODETYPE SCIPgetSymgraphNodeType(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   nodeidx             /**< index of node */
   );

/** returns the color of a non-variable node */
SCIP_EXPORT
int SCIPgetSymgraphNodeColor(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   nodeidx             /**< index of node */
   );

/** returns whether an edge is colored */
SCIP_EXPORT
SCIP_Bool SCIPisSymgraphEdgeColored(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   edgeidx             /**< index of edge */
   );

/** returns color of an edge */
SCIP_EXPORT
int SCIPgetSymgraphEdgeColor(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   edgeidx             /**< index of edge */
   );

/** returns the number of unique variable colors in the graph */
SCIP_EXPORT
int SCIPgetSymgraphNVarcolors(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   );

/** returns whether the graph has a unique edge type */
SCIP_EXPORT
SCIP_Bool SCIPhasGraphUniqueEdgetype(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   );

/** creates consnodeperm array for symmetry detection graph
 *
 *  @note @p colors of symmetry detection graph must have been computed
 */
SCIP_EXPORT
SCIP_RETCODE SCIPallocateSymgraphConsnodeperm(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   );

/** creates consnodeperm array for symmetry detection graph
 *
 *  @note @p colors of symmetry detection graph must have been computed
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateSymgraphConsnodeperm(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   );

/** returns consnodeperm array for symmetry detection graph
 *
 *  @note @p colors of symmetry detection graph must have been computed
 */
SCIP_EXPORT
int* SCIPgetSymgraphConsnodeperm(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   );

/** frees consnodeperm array for symmetry detection graph */
SCIP_EXPORT
SCIP_RETCODE SCIPfreeSymgraphConsnodeperm(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   );

/** Transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant.
 *
 *  For permutation symmetries, active variables as encoded in SCIP are used. For signed permutation symmetries,
 *  active variables are shifted such that their domain is centered at 0 (if both their upper and lower bounds
 *  are finite).
 *
 *  @note @p constant needs to be initialized!
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetSymActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SYMTYPE           symtype,            /**< type of symmetries for which variables are required */
   SCIP_VAR***           vars,               /**< pointer to vars array to get active variables for */
   SCIP_Real**           scalars,            /**< pointer to scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c */
   SCIP_Bool             transformed         /**< transformed constraint? */
   );

/** frees symmetry information of an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPfreeSymDataExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_EXPRDATA**        symdata             /**< symmetry information of an expression */
   );

/** returns number of coefficients from exprdata */
SCIP_EXPORT
int SCIPgetSymExprdataNConstants(
   SYM_EXPRDATA*         symdata             /**< symmetry information of an expression */
   );

/** returns number of coefficients form exprdata */
SCIP_EXPORT
SCIP_Real* SCIPgetSymExprdataConstants(
   SYM_EXPRDATA*         symdata             /**< symmetry information of an expression */
   );

/** gets coefficient of expression from parent expression */
SCIP_EXPORT
SCIP_RETCODE SCIPgetCoefSymData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression for which coefficient needs to be found */
   SCIP_EXPR*            parentexpr,         /**< parent of expr */
   SCIP_Real*            coef,               /**< buffer to store coefficient */
   SCIP_Bool*            success             /**< whether a coefficient is found */
   );


/** @} */

#ifdef __cplusplus
}
#endif

#endif
