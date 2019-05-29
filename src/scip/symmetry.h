/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
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
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicCutMethods
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
SCIP_RETCODE SCIPcomputeGroupOrbitsSymbreak(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR**            permvars,           /**< variables considered by symbreak presolver */
   int                   npermvars,          /**< length of a permutation array */
   int**                 perms,              /**< matrix containing in each row a permutation of the symmetry group */
   int                   nperms,             /**< number of permutations encoded in perms */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int*                  norbits             /**< pointer to number of orbits currently stored in orbits */
   );


/** compute non-trivial orbits of symmetry group using filtered generators
 *
 *  The non-tivial orbits of the group action are stored in the array orbits of length npermvars. This array contains
 *  the indices of variables from the permvars array such that variables that are contained in the same orbit appear
 *  consecutively in the orbits array. The variables of the i-th orbit have indices
 *  orbits[orbitbegins[i]], ... , orbits[orbitbegins[i + 1] - 1].
 *  Note that the description of the orbits ends at orbitbegins[norbits] - 1.
 */
SCIP_EXPORT
SCIP_RETCODE computeGroupOrbitsFilter(
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
   SCIP_Shortbool*       componentblocked,   /**< array to store whether a component is blocked to be considered by
                                              *   further symmetry handling techniques */
   int                   ncomponents,        /**< number of components of symmetry group */
   int                   nmovedpermvars      /**< number of variables moved by any permutation in a symmetry component
                                              *   that is handled by OF */
   );


/** check whether a permutation is a composition of 2-cycles of binary variables and in this case determines the number of 2-cycles */
SCIP_EXPORT
SCIP_RETCODE getPermProperties(
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< array of variables perm is acting on */
   int                   nvars,              /**< number of variables */
   SCIP_Bool*            iscompoftwocycles,  /**< pointer to store whether permutation is a composition of 2-cycles */
   int*                  ntwocyclesperm,     /**< pointer to store number of 2-cycles */
   SCIP_Bool*            allvarsbinary       /**< pointer to store whether perm is acting on binary variables only */
   );


/** determine whether binary variable is affected (and potentially compute number of affected variables) */
SCIP_RETCODE SCIPdetermineBinvarAffected(
   SCIP*                 scip,               /**< SCIP instance */
   int**                 perms,              /**< permutations */
   int                   nperms,             /**< number of permutations in perms */
   SCIP_VAR**            permvars,           /**< variables corresponding to permutations */
   int                   npermvars,          /**< number of permvars in perms */
   SCIP_Bool             completestatistic,  /**< whether number of affected vars should be computed */
   int*                  naffected           /**< pointer to store number of affected vars */
   );


/** Given a matrix with nrows and \#perms + 1 columns whose first nfilledcols columns contain entries of variables, this routine
 *  checks whether the 2-cycles of perm intersect each row of column coltoextend in exactly one position. In this case,
 *  we add one column to the suborbitope of the first nfilledcols columns.
 *
 *  @pre Every non-trivial cycle of perm is a 2-cycle.
 *  @pre perm has nrows many 2-cycles
 */
SCIP_EXPORT
SCIP_RETCODE extendSubOrbitope(
   int**                 suborbitope,        /**< matrix containing suborbitope entries */
   int                   nrows,              /**< number of rows of suborbitope */
   int                   nfilledcols,        /**< number of columns of suborbitope which are filled with entries */
   int                   coltoextend,        /**< index of column that should be extended by perm */
   int*                  perm,               /**< permutation */
   SCIP_Bool             leftextension,      /**< whether we extend the suborbitope to the left */
   int**                 nusedelems,         /**< pointer to array storing how often an element was used in the orbitope */
   SCIP_Bool*            success,            /**< pointer to store whether extension was successful */
   SCIP_Bool*            infeasible          /**< pointer to store if the number of intersecting cycles is too small */
   );


/** generate variable matrix for orbitope constraint handler */
SCIP_EXPORT
SCIP_RETCODE generateOrbitopeVarsMatrix(
   SCIP_VAR****          vars,               /**< pointer to matrix of orbitope variables */
   int                   nrows,              /**< number of rows of orbitope */
   int                   ncols,              /**< number of columns of orbitope */
   SCIP_VAR**            permvars,           /**< superset of variables that are contained in orbitope */
   int                   npermvars,          /**< number of variables in permvars array */
   int**                 orbitopevaridx,     /**< permuted index table of variables in permvars that are contained in orbitope */
   int*                  columnorder,        /**< permutation to reorder column of orbitopevaridx */
   int*                  nusedelems,         /**< array storing how often an element was used in the orbitope */
   SCIP_Bool*            infeasible          /**< pointer to store whether the potential orbitope is not an orbitope */
   );


/* @} */

#ifdef __cplusplus
}
#endif

#endif
