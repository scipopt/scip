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

/**@file   symmetry.h
 * @ingroup PUBLICCOREAPI
 * @brief  methods for handling symmetries
 * @author Christopher Hojny
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SYMMETRY_H__
#define __SCIP_SYMMETRY_H__

#include "scip/def.h"
#include "scip/pub_misc.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"
#include <symmetry/type_symmetry.h>

#ifdef __cplusplus
extern "C" {
#endif


/**@addtogroup PublicSymmetryMethods
 *
 * @{
 */


/** compute non-trivial orbits of symmetry group
 *
 *  The non-tivial orbits of the group action are stored in the array orbits of length npermvars. This array contains
 *  the indices of variables from the permvars array such that variables that are contained in the same orbit appear
 *  consecutively in the orbits array. The variables of the i-th orbit have indices
 *  orbits[orbitbegins[i]], ... , orbits[orbitbegins[i + 1] - 1].
 *  Note that the description of the orbits ends at orbitbegins[norbits] - 1.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeOrbitsSym(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR**            permvars,           /**< variables considered in a permutation array */
   int                   npermvars,          /**< length of a permutation array */
   int**                 perms,              /**< matrix containing in each row a permutation of the symmetry group */
   int                   nperms,             /**< number of permutations encoded in perms */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int*                  norbits             /**< pointer to number of orbits currently stored in orbits */
   );


/** compute non-trivial orbits of symmetry group using filtered generators
 *
 *  The non-trivial orbits of the group action are stored in the array orbits of length npermvars. This array contains
 *  the indices of variables from the permvars array such that variables that are contained in the same orbit appear
 *  consecutively in the orbits array. The variables of the i-th orbit have indices
 *  orbits[orbitbegins[i]], ... , orbits[orbitbegins[i + 1] - 1].
 *  Note that the description of the orbits ends at orbitbegins[norbits] - 1.
 *
 *  Only permutations that are not inactive (as marked by @p inactiveperms) are used. Thus, one can use this array to
 *  filter out permutations.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeOrbitsFilterSym(
   SCIP*                 scip,               /**< SCIP instance */
   int                   npermvars,          /**< length of a permutation array */
   int**                 permstrans,         /**< transposed matrix containing in each column a
                                              *   permutation of the symmetry group */
   int                   nperms,             /**< number of permutations encoded in perms */
   SCIP_Shortbool*       inactiveperms,      /**< array to store whether permutations are inactive */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int*                  norbits,            /**< pointer to number of orbits currently stored in orbits */
   int*                  components,         /**< array containing the indices of permutations sorted by components */
   int*                  componentbegins,    /**< array containing in i-th position the first position of
                                              *   component i in components array */
   int*                  vartocomponent,     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   unsigned*             componentblocked,   /**< array to store which symmetry methods have been used on a component
                                              *   using the same bitset information as for misc/usesymmetry */
   int                   ncomponents,        /**< number of components of symmetry group */
   int                   nmovedpermvars      /**< number of variables moved by any permutation in a symmetry component
                                              *   that is handled by orbital fixing */
   );

/** compute non-trivial orbits of symmetry group
 *
 *  The non-tivial orbits of the group action are stored in the array orbits of length npermvars. This array contains
 *  the indices of variables from the permvars array such that variables that are contained in the same orbit appear
 *  consecutively in the orbits array. The variables of the i-th orbit have indices
 *  orbits[orbitbegins[i]], ... , orbits[orbitbegins[i + 1] - 1].
 *  Note that the description of the orbits ends at orbitbegins[norbits] - 1.
 *
 *  This function is adapted from SCIPcomputeOrbitsFilterSym().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeOrbitsComponentsSym(
   SCIP*                 scip,               /**< SCIP instance */
   int                   npermvars,          /**< length of a permutation array */
   int**                 permstrans,         /**< transposed matrix containing in each column a permutation of the symmetry group */
   int                   nperms,             /**< number of permutations encoded in perms */
   int*                  components,         /**< array containing the indices of permutations sorted by components */
   int*                  componentbegins,    /**< array containing in i-th position the first position of component i in components array */
   int*                  vartocomponent,     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   int                   ncomponents,        /**< number of components of symmetry group */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int*                  norbits,            /**< pointer to number of orbits currently stored in orbits */
   int*                  varorbitmap         /**< array for storing the orbits for each variable */
   );

/** Compute orbit of a given variable and store it in @p orbit. The first entry of the orbit will
 *  be the given variable index and the rest is filled with the remaining variables excluding
 *  the ones specified in @p ignoredvars.
 *
 *  @pre orbit is an initialized array of size propdata->npermvars
 *  @pre at least one of @p perms and @p permstrans should not be NULL
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeOrbitVar(
   SCIP*                 scip,               /**< SCIP instance */
   int                   npermvars,          /**< number of variables in permvars */
   int**                 perms,              /**< the generators of the permutation group (or NULL) */
   int**                 permstrans,         /**< the transposed matrix of generators (or NULL) */
   int*                  components,         /**< the components of the permutation group */
   int*                  componentbegins,    /**< array containing the starting index of each component */
   SCIP_Shortbool*       ignoredvars,        /**< array indicating which variables should be ignored */
   SCIP_Shortbool*       varfound,           /**< bitmap to mark which variables have been added (or NULL) */
   int                   varidx,             /**< index of variable for which the orbit is requested */
   int                   component,          /**< component that var is in */
   int *                 orbit,              /**< array in which the orbit should be stored */
   int*                  orbitsize           /**< buffer to store the size of the orbit */
   );

/** Checks whether a permutation is a composition of 2-cycles and in this case determine the number of overall
 *  2-cycles and binary 2-cycles. It is a composition of 2-cycles iff @p ntwocyclesperm > 0 upon termination.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPisInvolutionPerm(
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< array of variables perm is acting on */
   int                   nvars,              /**< number of variables */
   int*                  ntwocyclesperm,     /**< pointer to store number of 2-cycles */
   int*                  nbincyclesperm,     /**< pointer to store number of binary cycles */
   SCIP_Bool             earlytermination    /**< whether we terminate early if not all affected variables are binary */
   );

/** determine number of variables affected by symmetry group */
SCIP_EXPORT
SCIP_RETCODE SCIPdetermineNVarsAffectedSym(
   SCIP*                 scip,               /**< SCIP instance */
   int**                 perms,              /**< permutations */
   int                   nperms,             /**< number of permutations in perms */
   SCIP_VAR**            permvars,           /**< variables corresponding to permutations */
   int                   npermvars,          /**< number of permvars in perms */
   int*                  nvarsaffected       /**< pointer to store number of all affected variables */
   );

/** compute components of symmetry group */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeComponentsSym(
   SCIP*                 scip,               /**< SCIP instance */
   int**                 perms,              /**< permutation generators as
                                              *   (either nperms x npermvars or npermvars x nperms) matrix */
   int                   nperms,             /**< number of permutations */
   SCIP_VAR**            permvars,           /**< variables on which permutations act */
   int                   npermvars,          /**< number of variables for permutations */
   SCIP_Bool             transposed,         /**< transposed permutation generators as (npermvars x nperms) matrix */
   int**                 components,         /**< array containing the indices of permutations sorted by components */
   int**                 componentbegins,    /**< array containing in i-th position the first position of
                                              *   component i in components array */
   int**                 vartocomponent,     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   unsigned**            componentblocked,   /**< array to store which symmetry methods have been used on a component
                                              *   using the same bitset information as for misc/usesymmetry */
   int*                  ncomponents         /**< pointer to store number of components of symmetry group */
   );

/** Given a matrix with nrows and \#perms + 1 columns whose first nfilledcols columns contain entries of variables, this routine
 *  checks whether the 2-cycles of perm intersect each row of column coltoextend in exactly one position. In this case,
 *  we add one column to the suborbitope of the first nfilledcols columns.
 *
 *  @pre Every non-trivial cycle of perm is a 2-cycle.
 *  @pre perm has nrows many 2-cycles
 */
SCIP_EXPORT
SCIP_RETCODE SCIPextendSubOrbitope(
   int**                 suborbitope,        /**< matrix containing suborbitope entries */
   int                   nrows,              /**< number of rows of suborbitope */
   int                   nfilledcols,        /**< number of columns of suborbitope which are filled with entries */
   int                   coltoextend,        /**< index of column that should be extended by perm */
   int*                  perm,               /**< permutation */
   SCIP_Bool             leftextension,      /**< whether we extend the suborbitope to the left */
   int**                 nusedelems,         /**< pointer to array storing how often an element was used in the orbitope */
   SCIP_VAR**            permvars,           /**< permutation vars array */
   SCIP_Shortbool*       rowisbinary,        /**< array encoding whether variables in an orbitope row are binary */
   SCIP_Bool*            success,            /**< pointer to store whether extension was successful */
   SCIP_Bool*            infeasible          /**< pointer to store if the number of intersecting cycles is too small */
   );

/** generate variable matrix for orbitope constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPgenerateOrbitopeVarsMatrix(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR****          vars,               /**< pointer to matrix of orbitope variables */
   int                   nrows,              /**< number of rows of orbitope */
   int                   ncols,              /**< number of columns of orbitope */
   SCIP_VAR**            permvars,           /**< superset of variables that are contained in orbitope */
   int                   npermvars,          /**< number of variables in permvars array */
   int**                 orbitopevaridx,     /**< permuted index table of variables in permvars that are contained in orbitope */
   int*                  columnorder,        /**< permutation to reorder column of orbitopevaridx */
   int*                  nusedelems,         /**< array storing how often an element was used in the orbitope */
   SCIP_Shortbool*       rowisbinary,        /**< array encoding whether a row contains only binary variables */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the potential orbitope is not an orbitope */
   SCIP_Bool             storelexorder,      /**< whether the lexicographic order induced by the orbitope shall be stored */
   int**                 lexorder,           /**< pointer to array storing the lexorder (or NULL) */
   int*                  nvarsorder,         /**< pointer to store number of variables in lexorder (or NULL) */
   int*                  maxnvarsorder       /**< pointer to store maximum number of variables in lexorder (or NULL) */
   );

/** checks whether an orbitope is a packing or partitioning orbitope */
SCIP_RETCODE SCIPisPackingPartitioningOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< variable matrix of orbitope constraint */
   int                   nrows,              /**< pointer to number of rows of variable matrix */
   int                   ncols,              /**< number of columns of variable matrix */
   SCIP_Bool**           pprows,             /**< pointer to store which rows are are contained in
                                              *   packing/partitioning constraints or NULL if not needed */
   int*                  npprows,            /**< pointer to store how many rows are contained
                                              *   in packing/partitioning constraints or NULL if not needed */
   SCIP_ORBITOPETYPE*    type                /**< pointer to store type of orbitope constraint after strengthening */
   );

/** Transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant.
 *
 *  @note @p constant needs to be initialized!
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to vars array to get active variables for */
   SCIP_Real**           scalars,            /**< pointer to scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c */
   SCIP_Bool             transformed         /**< transformed constraint? */
   );

/** creates permutation symmetry detection graph for linear constraints
 *
 *  @note at some point, the graph needs to be freed!
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreatePermsymDetectionGraphLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH**           graph,              /**< pointer to hold symmetry detection graph */
   SCIP_VAR**            vars,               /**< variable array of linear constraint */
   SCIP_Real*            vals,               /**< coefficients of linear constraint */
   int                   nvars,              /**< number of variables in linear constraint */
   SCIP_Real             lhs,                /**< left-hand side of linear constraint */
   SCIP_Real             rhs,                /**< right-hand side of linear constraint */
   SCIP_CONS*            cons,               /**< constraint for which we encode symmetries */
   SCIP_Bool*            success             /**< pointer to store whether graph could be built */
   );

/** adds nodes and edges corresponding to the aggregation of a variable to a symmetry detection graph */
SCIP_EXPORT
SCIP_RETCODE SCIPaddSymgraphVarAggegration(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   rootidx,            /**< index of root of the aggegration */
   int*                  idx,                /**< pointer to store index of lastly added node
                                              *   (initialized with index of first node to be added) */
   SCIP_VAR**            vars,               /**< array of variables in aggregation */
   SCIP_Real*            vals,               /**< coefficients of variables */
   int                   nvars,              /**< number of variables in aggregation */
   SCIP_Real             constant            /**< constant of aggregation */
   );

/** frees symmetry information of an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPfreeSymdataExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_EXPRDATA2**       symdata             /**< symmetry information of an expression */
   );

/** gets coefficient of expression from parent expression */
SCIP_EXPORT
SCIP_RETCODE SCIPgetCoefSymdata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression for which coefficient needs to be found */
   SCIP_EXPR*            parentexpr,         /**< parent of expr */
   SCIP_Real*            coef,               /**< buffer to store coefficient */
   SCIP_Bool*            success             /**< whether a coefficient is found */
   );

/** creates and initializes a symmetry detection graph with memory for the given number of nodes and edges */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateSymgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH**           graph,              /**< pointer to hold symmetry detection graph */
   int                   nopnodes,           /**< number of operator nodes */
   int                   nvarnodes,          /**< number of variable nodes */
   int                   nvalnodes,          /**< number of value nodes */
   int                   nedges              /**< number of edges */
   );

/** frees a symmetry detection graph */
SCIP_EXPORT
SCIP_RETCODE SCIPfreeSymgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH**           graph               /**< pointer to hold symmetry detection graph */
   );

/** adds an operator node to a symmetry detection graph */
SCIP_EXPORT
SCIP_RETCODE SCIPaddSymgraphOpnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   SCIP_EXPRHDLR*        op                  /**< operator of node to be added */
   );

/** adds a variable node to a symmetry detection graph */
SCIP_EXPORT
SCIP_RETCODE SCIPaddSymgraphVarnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   SCIP_VAR*             var                 /**< variable of node to be added */
   );

/** adds a value node to a symmetry detection graph */
SCIP_EXPORT
SCIP_RETCODE SCIPaddSymgraphValnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   SCIP_Real             val                 /**< value of node to be added */
   );

/** adds a rhs node to a symmetry detection graph */
SCIP_RETCODE SCIPaddSymgraphRhsnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   SCIP_CONS*            cons,               /**< constraint associated with the node */
   SCIP_Real             lhs,                /**< lhs associated with the node */
   SCIP_Real             rhs                 /**< rhs associated with the node */
   );

/** adds edge to a symmetry detection graph */
SCIP_EXPORT
SCIP_RETCODE SCIPaddSymgraphEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   first,              /**< node index of first node in edge */
   int                   second,             /**< node index of second node in edge */
   SCIP_Bool             iscolored,          /**< whether the edge is colored */
   SCIP_Real             color               /**< color of edge (if it is colored) */
   );

/** returns type of a node */
SCIP_EXPORT
SYM_NODETYPE SCIPgetSymgraphNodeType(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   nodeidx             /**< index of node whose type needs to be returned */
   );

/** returns operator of an operator node */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetSymgraphNodeOperator(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   nodeidx             /**< index of node whose operator needs to be returned */
   );

/** returns variable of a variable node */
SCIP_EXPORT
SCIP_VAR* SCIPgetSymgraphNodeVar(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   nodeidx             /**< index of node whose variable needs to be returned */
   );

/** returns value of a value node */
SCIP_EXPORT
SCIP_Real SCIPgetSymgraphNodeVal(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   nodeidx             /**< index of node whose value needs to be returned */
   );

/** returns node color used for symmetry detection */
SCIP_EXPORT
int SCIPgetSymgraphNodeSymcolor(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   nodeidx             /**< index of node whose color needs to be returned */
   );

/** returns node color of rhs node (or -1 if node does not exist) */
SCIP_EXPORT
int SCIPgetSymgraphRhsnodeSymcolor(
   SYM_GRAPH*            graph               /**< pointer to symmetry detection graph */
   );

/** returns node color used for symmetry detection */
SCIP_EXPORT
SCIP_Bool SCIPhasSymgraphNodeSymcolor(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   nodeidx             /**< index of node whose color needs to be returned */
   );

/** sets node color used for symmetry detection */
SCIP_EXPORT
SCIP_RETCODE SCIPsetSymgraphNodeSymcolor(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   nodeidx,            /**< index of node whose color needs to be returned */
   int                   color               /**< color used for symmetry detection */
   );

/** sets rhs node color used for symmetry detection */
SCIP_EXPORT
SCIP_RETCODE SCIPsetSymgraphRhscolor(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   color               /**< color used for symmetry detection */
   );

/** gets first and second node of an edge */
SCIP_EXPORT
SCIP_RETCODE SCIPgetSymgraphEdge(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   edgeidx,            /**< index of edge that needs to be returned */
   int*                  first,              /**< buffer to store first index of node in edge */
   int*                  second              /**< buffer to store second index of node in edge */
   );

/** returns whether an edge is colored */
SCIP_EXPORT
SCIP_RETCODE SCIPisSymgraphEdgeColored(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   edgeidx             /**< index of edge */
   );

/** returns the color of an edge (or infinity if not colored) */
SCIP_EXPORT
SCIP_Real SCIPgetSymgraphEdgeColor(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   edgeidx             /**< index of edge */
   );

/** returns edge color used for symmetry detection */
SCIP_EXPORT
int SCIPgetSymgraphEdgeSymcolor(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   edgeidx             /**< index of edge whose color needs to be returned */
   );

/** sets edge color used for symmetry detection */
SCIP_EXPORT
SCIP_RETCODE SCIPsetSymgraphEdgeSymcolor(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   edgeidx,            /**< index of edge whose color needs to be returned */
   int                   color               /**< color used for symmetry detection */
   );

/** returns whether edge colors for symmetry detection have been computed */
SCIP_EXPORT
SCIP_Bool SCIPhasSymgraphEdgeSymcolor(
   SYM_GRAPH*            graph               /**< pointer to symmetry detection graph */
   );

/** gets constraint information from a symmetry detection graph */
SCIP_EXPORT
SCIP_RETCODE SCIPgetSymgraphConsinfo(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   SCIP_Bool*            hasconsinfo,        /**< buffer to store whether graph has constraint information */
   SCIP_CONS**           cons,               /**< buffer to store pointer to constraint (if graph has information) */
   SCIP_Real*            lhs,                /**< buffer to store lhs (if graph has information) */
   SCIP_Real*            rhs                 /**< buffer to store rhs (if graph has information) */
   );

/** return the number of nodes in a symmetry detection graph */
SCIP_EXPORT
int SCIPgetSymgraphNnodes(
   SYM_GRAPH*            graph               /**< pointer to symmetry detection graph */
   );

/** return the number of operator nodes in a symmetry detection graph */
SCIP_EXPORT
int SCIPgetSymgraphNopnodes(
   SYM_GRAPH*            graph               /**< pointer to symmetry detection graph */
   );

/** return the number of variable nodes in a symmetry detection graph */
SCIP_EXPORT
int SCIPgetSymgraphNvarnodes(
   SYM_GRAPH*            graph               /**< pointer to symmetry detection graph */
   );

/** return the number of value nodes in a symmetry detection graph */
SCIP_EXPORT
int SCIPgetSymgraphNvalnodes(
   SYM_GRAPH*            graph               /**< pointer to symmetry detection graph */
   );

/** return whether the graph has a rhs node */
SCIP_EXPORT
SCIP_Bool SCIPhasSymgraphRhsnode(
   SYM_GRAPH*            graph               /**< pointer to symmetry detection graph */
   );

/** return the number of edges in a symmetry detection graph */
SCIP_EXPORT
int SCIPgetSymgraphNedges(
   SYM_GRAPH*            graph               /**< pointer to symmetry detection graph */
   );

/** returns whether a node of a symmetry detection graph is an operator node */
SCIP_EXPORT
SCIP_Bool SCIPisSymgraphNodeOpnode(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   nodeidx             /**< index of node in graph */
   );

/** returns whether a node of a symmetry detection graph is a variable node */
SCIP_EXPORT
SCIP_Bool SCIPisSymgraphNodeVarnode(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   nodeidx             /**< index of node in graph */
   );

/** returns whether a node of a symmetry detection graph is a value node */
SCIP_EXPORT
SCIP_Bool SCIPisSymgraphNodeValnode(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   nodeidx             /**< index of node in graph */
   );

/** returns whether a node of a symmetry detection graph is a rhs node */
SCIP_EXPORT
SCIP_Bool SCIPisSymgraphNodeRhsnode(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   nodeidx             /**< index of node in graph */
   );

/** returns the variable of a variable node */
SCIP_EXPORT
SCIP_VAR* SCIPgetSymgraphVarnodeVar(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   nodeidx             /**< index of node in graph */
   );

/** returns the first node index of an edge */
SCIP_EXPORT
int SCIPgetSymgraphEdgeFirst(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   edgeidx             /**< index of edge in graph */
   );

/** returns the second node index of an edge */
SCIP_EXPORT
int SCIPgetSymgraphEdgeSecond(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   int                   edgeidx             /**< index of edge in graph */
   );

/** updates lhs constraint information of a symmetry detection graph */
SCIP_EXPORT
SCIP_RETCODE SCIPupdateSymgraphLhs(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   SCIP_Real             lhs                 /**< lhs */
   );

/** updates rhs constraint information of a symmetry detection graph */
SCIP_EXPORT
SCIP_RETCODE SCIPupdateSymgraphRhs(
   SYM_GRAPH*            graph,              /**< pointer to symmetry detection graph */
   SCIP_Real             rhs                 /**< rhs */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
