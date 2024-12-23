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

/** data to encode a symmetry detection graph */
struct SYM_Graph
{
   /* general information about graph */
   SYM_SYMTYPE           symtype;            /**< type of symmetries encoded in graph */
   SCIP_Bool             islocked;           /**< whether graph is locked, i.e., cannot be modified anymore
                                              *   (computing colors will lock the graph to avoid inconsistencies) */
   SCIP_Real             infinity;           /**< values as least as large as this are regarded as infinite */

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

   /* node-based arrays */
   SYM_NODETYPE*         nodetypes;          /**< array storing each node's type */
   int*                  nodeinfopos;        /**< array assigning each node the position in the corresponding array
                                              *   containing its information (operator, variable, or value) */
   int*                  consnodeperm;       /**< array to hold permutation to sort constraint nodes
                                              *   (graph needs to be locked to avoid inconsistencies) */

   /* information-based arrays */
   int*                  ops;                /**< operators corresponding to nodes in graph */
   SCIP_Real*            vals;               /**< values corresponding to nodes in graph */
   SCIP_CONS**           conss;              /**< constraints corresponding to cons nodes */
   SCIP_Real*            lhs;                /**< array of left-hand sides for cons nodes */
   SCIP_Real*            rhs;                /**< array of right-hand sides for cons nodes */

   /* information about edges and edge arrays */
   int                   nedges;             /**< number of edges in graph */
   int                   maxnedges;          /**< maximum number of entries in edge-based arrays */

   /* edge-based arrays; negative entries in edge-first and edge-second correspond to variables */
   int*                  edgefirst;          /**< array of first nodes of edges */
   int*                  edgesecond;         /**< array of second nodes of edges */
   SCIP_Real*            edgevals;           /**< array assigning each edge a value (SCIPinfinity if unassigned) */

   /* information about variables */
   SCIP_VAR**            symvars;            /**< variables on which symmetries act */
   int                   nsymvars;           /**< number of variables in symvars */
   SCIP_Bool*            isfixedvar;         /**< whether a variable needs to be fixed */

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
