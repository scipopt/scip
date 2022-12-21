/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
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

/**@file   struct_symmetry.h
 * @brief  structs for symmetry computations
 * @author Marc Pfetsch
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SYMMETRY_H_
#define __SCIP_STRUCT_SYMMETRY_H_

#include "scip/scip.h"
#include "symmetry/type_symmetry.h"
#include "scip/type_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data of variables that are considered to be equivalent */
struct SYM_Vartype
{
   SCIP_Real             obj;                /**< objective of variable */
   SCIP_Real             lb;                 /**< lower bound of variable */
   SCIP_Real             ub;                 /**< upper bound of variable */
   SCIP_VARTYPE          type;               /**< type of variable */
   int                   nconss;             /**< number of conss a variable is contained in */
   int                   color;              /**< store color */
};

/** data of operators that are considered to be equivalent */
struct SYM_Optype
{
   SCIP_EXPR*            expr;               /**< the underlying expression */
   int                   level;              /**< level of operator in its expression tree */
   int                   color;              /**< store color */
};

/** data of constants that are considered to be equivalent */
struct SYM_Consttype
{
   SCIP_Real             value;              /**< value of constant */
   int                   color;              /**< store color */
};

/** data of coefficients that are considered to be equivalent */
struct SYM_Rhstype
{
   SCIP_Real             lhs;                /**< value of left-hand-side */
   SCIP_Real             rhs;                /**< value of right-hand-side */
   int                   color;              /**< store color */
};

/** data to compare a graph */
struct SYM_Compgraph
{
   int**                 adjacentnodecolors;
   int**                 incidentedgecolors;
   int*                  nneighbors;
   int*                  nodecolors;
};

/** data for symmetry group computation on linear constraints */
struct SYM_Matrixdata
{
   SCIP_Real*            matcoef;            /**< nonzero coefficients appearing in the matrix */
   SCIP_Real*            rhscoef;            /**< rhs coefficients */
   SYM_RHSSENSE*         rhssense;           /**< sense of rhs */
   int*                  matrhsidx;          /**< indices of rhs corresponding to matrix entries */
   int*                  matvaridx;          /**< indices of variables for matrix entries */
   int*                  matidx;             /**< indices in mat(rhs/var)idx array corresponding to matrix coefficients */
   int*                  rhsidx;             /**< indices in rhstype array corresponding to rhs coefficients */
   int*                  permvarcolors;      /**< array for storing the colors of the individual variables */
   int*                  matcoefcolors;      /**< array for storing the colors of all matrix coefficients */
   int*                  rhscoefcolors;      /**< array for storing the colors of all rhs coefficients */
   SCIP_VAR**            permvars;           /**< variables on which permutations act */
   int                   npermvars;          /**< number of variables for permutations */
   int                   nmatcoef;           /**< number of coefficients in matrix */
   int                   nrhscoef;           /**< number of coefficients in rhs */
   int                   nmaxmatcoef;        /**< maximal number of matrix coefficients (will be increase on demand) */
   int                   nuniquevars;        /**< number of unique variable types */
   int                   nuniquerhs;         /**< number of unique rhs types */
   int                   nuniquemat;         /**< number of unique matrix coefficients */
};

/** data for symmetry group computation on nonlinear constraints */
struct SYM_Exprdata
{
   int                   nuniqueconstants;   /**< number of unique constants */
   int                   nuniqueoperators;   /**< number of unique operators */
   int                   nuniquecoefs;       /**< number of unique coefficients */
};

/** data for reflection symmetry group computation */
struct SYM_Reflsymdata
{
   SYM_NODETYPE*         trees;              /**< array to encode all constraints as consecutive expression trees
                                              *   (trees are encoded by DFS-traversal, coefficients appear before
                                              *    the corresponding variables) */
   SCIP_VAR**            treevars;           /**< array of unique variables that appear in expression trees */
   int*                  treebegins;         /**< array containing begin positions of new tree in trees */
   SYM_CONSTYPE*         treeconstype;       /**< array of constraint types stores in trees */
   SCIP_Real*            treerhs;            /**< right-hand side coefficients of trees */
   int*                  treeparentidx;      /**< array assigning each position in trees the position of its parent
                                              *   (or -1 in case the position corresponds to the root of a tree) */
   int*                  treevaridx;         /**< indices of variables in expression trees (order according t trees) */
   SCIP_Real*            treecoefs;          /**< var coefficients in expression trees (order according to trees) */
   SCIP_Real*            treevals;           /**< numerical values in expression trees (order according to trees) */
   SCIP_EXPRHDLR**       treeops;            /**< operators used in expression trees (order according to trees) */
   int*                  treemap;            /**< maps position in trees array to the corresponding position in
                                              *   treecoefs/treevals/treeops/treevaridx (depending on node type) */
   int*                  rhsidx;             /**< maps index of treerhs to index of corresponding tree */
   int*                  varidx;             /**< maps index of treevaridx to position in trees */
   int*                  coefidx;            /**< maps index of treecoefs to position in trees */
   int*                  validx;             /**< maps index of treevals to position in trees */
   int*                  opsidx;             /**< maps index of treeops to position in trees */
   int                   ntrees;             /**< number of elements in trees */
   int                   ntreevars;          /**< number of elements in treevars */
   int                   ntreevaridx;        /**< number of elements in treevaridx */
   int                   ntreerhs;           /**< number of elements in treerhs */
   int                   ntreecoefs;         /**< number of elements in treecoefs */
   int                   ntreevals;          /**< number of elements in treevals */
   int                   ntreeops;           /**< number of elements in treeops */
   int                   maxntrees;          /**< maximum number of elements that fit into trees */
   int                   maxntreerhs;        /**< maximum number of elements that fit into treerhs */
   int                   maxntreevaridx;     /**< maximum number of elements that fit into treevaridx */
   int                   maxntreecoefs;      /**< maximum number of elements that fit into treecoefs */
   int                   maxntreevals;       /**< maximum number of elements that fit into treevals */
   int                   maxntreeops;        /**< maximum number of elements that fit into treeops */
   int                   maxntreebegins;     /**< maximum number of elements that fit into treebegins */
   int*                  varcolors;          /**< array to store colors of individual variables */
   int*                  invvarcolors;       /**< array to store colors of individual negated variables */
   int*                  coefcolors;         /**< array to store colors of individual coefficients */
   int*                  invcoefcolors;      /**< array to store inverse colors of individual coefficients */
   int*                  opscolors;          /**< array to store colors of individual operators */
   int*                  valcolors;          /**< array to store colors of individual values */
   int*                  rhscolors;          /**< array to store colors of individual right-hand sides */
   int                   nuniquevars;        /**< number of unique variable types */
   int                   nuniquecoefs;       /**< number of unique coefficients */
   int                   nuniqueops;         /**< number of unique operators */
   int                   nuniquevals;        /**< number of unique values */
   int                   nuniquerhs;         /**< number of unique right-hand sides */
};

/** information about a constraint used in symmetry computation */
struct SYM_Consinfo
{
   SCIP_CONSHDLR*        conshdlr;           /**< pointer to the constraint handler of a constraint */
   SCIP_CONS*            cons;               /**< pointer to the constraint */
   SCIP_Real             lhs;                /**< left-hand side of constraint */
   SCIP_Real             rhs;                /**< right-hand side of constraint */
};

/** data to encode a node of a symmetry detection graph */
struct SYM_Node
{
   int                   id;                 /**< ID of the node (must be unique in the graph) */
   SYM_NODETYPE          nodetype;           /**< type of the node */
   SCIP_EXPRHDLR*        op;                 /**< operator encoded by the node
                                              *   (if nodetype is SYM_NODETYPE_OPERATOR) */
   SCIP_VAR*             var;                /**< variable encoded by the node
                                              *   (if nodetype is SYM_NODETYPE_VAR) */
   int                   varidx;             /**< index of variable encoded by the node
                                              *   (if nodetype is SYM_NODETYPE_VAR) */
   SCIP_Real             value;              /**< index of value encoded by the node
                                              *   (if nodetype is SYM_NODETYPE_VAL) */
   SCIP_Bool             hasinfo;            /**< whether the node encodes information about a constraint */
   SYM_CONSINFO*         consinfo;           /**< pointer to information about constraint (or NULL) */
   int                   computedcolor;      /**< color computed for symmetry detection */
};

/** data to encode an edge of a symmetry detection graph */
struct SYM_Edge
{
   SYM_NODE*             first;              /**< pointer to the first node of an edge */
   SYM_NODE*             second;             /**< pointer to the second node of an edge */
   SCIP_Bool             iscolored;          /**< whether the edge is colored */
   SCIP_Real             color;              /**< color of the edge (if it is colored) */
   int                   computedcolor;      /**< color computed for symmetry detection */
};

/** data to encode a symmetry detection graph */
struct SYM_Graph
{
   SYM_NODE**            nodes;              /**< array of nodes in the graph */
   int                   nnodes;             /**< number of nodes encoded in nodes */
   int                   maxnnodes;          /**< maximum number of nodes that can be hold by nodes */
   SYM_EDGE**            edges;              /**< array of edges in the graph */
   int                   nedges;             /**< number of edges encoded in edged */
   int                   maxnedges;          /**< maximum number of edges that can be hold by edges */
};

#ifdef __cplusplus
}
#endif

#endif
