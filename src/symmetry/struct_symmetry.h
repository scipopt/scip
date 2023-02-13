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

/** data to encode a symmetry detection graph */
struct SYM_Graph
{
   /* information about nodes and node arrays */
   int                   nnodes;             /**< number of nodes in graph */
   int                   maxnnodes;          /**< maximum number of entries in node-based arrays */
   int                   nopnodes;           /**< number of operator nodes in graph */
   int                   maxnopnodes;        /**< maximum number of entries in operator-based arrays */
   int                   nvalnodes;          /**< number of value nodes in graph */
   int                   maxnvalnodes;       /**< maximum number of entries in value-based arrays */
   int                   nconsnodes;         /**< number of constraint nodes */
   int                   maxnconsnodes;      /**< maximum number of constraint-based arrays */
   int                   nvarcolors;         /**< number of variable colors */
   SCIP_Bool             islocked;           /**< whether graph is locked, i.e., cannot be modified anymore
                                              *   (computing colors will lock the graph to avoid incosistencies) */

   /* node-based arrays */
   SYM_NODETYPE*         nodetypes;          /**< array storing each node's type */
   int*                  nodeinfopos;        /**< array assigning each node the position in the corresponding
                                              *   containing its information (operator, variable, or value) */

   /* information-based arrays */
   int*                  ops;                /**< operators corresponding to nodes in graph */
   SCIP_Real*            vals;               /**< values corresponding to nodes in graph */
   SCIP_CONS**           conss;              /**< constraints corresponding to cons nodes */
   SCIP_Real*            lhs;                /**< array of left-hand sides for cons nodes */
   SCIP_Real*            rhs;                /**< array of right-hand sides for cons nodes */

   /* information about edges and edge arrays */
   int                   nedges;             /**< number of edges in graph */
   int                   maxnedges;          /**< maximum number of entries in edge-based arrays */

   /* edge-based arrays; negative entries in edgefirst and edgesecond correspond to variables */
   int*                  edgefirst;          /**< array of first nodes of edges */
   int*                  edgesecond;         /**< array of second nodes of edges */
   SCIP_Real*            edgevals;           /**< array assigning each edge a value (SCIPinfinity if unassigned) */

   /* information about variables */
   SCIP_VAR**            symvars;            /**< variables on which symmetries act */
   int                   nsymvars;           /**< number of variables in symvars */

   /* arrays of colors used for symmetry detection */
   int*                  varcolors;          /**< variable colors for symmetry detection */
   int*                  opcolors;           /**< operator colors for symmetry detection */
   int*                  valcolors;          /**< value colors for symmetry detection */
   int*                  conscolors;         /**< constraint colors for symmetry detection */
   int*                  edgecolors;         /**< edge colors used for symmetry detection (-1 uncolored) */
   SCIP_Bool             uniqueedgetype;     /**< whether all edges are equivalent */
};

/** (additional) data used to encode an expression, which is not encoded as another expression */
struct SYM_ExprData
{
   SCIP_Real*            constants;          /**< constants used in an expression */
   int                   nconstants;         /**< number of constants */
   SCIP_Real*            coefficients;       /**< coefficients of children */
   int                   ncoefficients;      /**< number of coefficients */
   SCIP_EXPR**           children;           /**< children of expression with a coefficient */
};

#ifdef __cplusplus
}
#endif

#endif
