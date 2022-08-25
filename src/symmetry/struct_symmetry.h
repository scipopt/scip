/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_symmetry.h
 * @brief  structs for symmetry computations
 * @author Marc Pfetsch
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
   SCIP_Real*            treerhs;            /**< right-hand side coefficients of trees */
   int*                  treeparentidx;      /**< array assigning each position in trees the position of its parent
                                              *   (or -1 in case the position corresponds to the root of a tree) */
   int*                  treevaridx;         /**< indices of variables in expression trees (order according to trees) */
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
};

#ifdef __cplusplus
}
#endif

#endif
