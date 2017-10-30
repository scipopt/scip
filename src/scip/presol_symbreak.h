/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_symbreak.h
 * @ingroup PRESOLVERS
 * @brief  presolver for adding symmetry breaking constraints
 * @author Marc Pfetsch
 * @author Thomas Rehn
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_SYMBREAK_H__
#define __SCIP_PRESOL_SYMBREAK_H__

#include "scip/scip.h"

/** creates the symmetry breaking presolver and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePresolSymbreak(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** compute non-trivial orbits of symmetry group
 *
 *  The orbits of the group action are stored in the array orbits of length npermvars. This array contains
 *  the indices of variables from the permvars array such that variables that are contained in the same
 *  orbit appear consecutively in the orbits array. To detect when a new orbit starts, we use the array
 *  orbitbegins:
 *    The variables of the i-th orbit have indices orbits[orbitbegins[i]], ... , orbits[orbitbegins[i + 1] - 1]
 */
EXTERN
SCIP_RETCODE computeGroupOrbits(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR**            permvars,           /**< variables considered by symbreak presolver */
   int                   npermvars,          /**< length of a permutation array */
   int**                 perms,              /**< matrix containing in each row a permutation of the symmetry group */
   int                   nperms,             /**< number of permutations encoded in perms */
   SCIP_Bool*            activeperms,        /**< array for marking active permutations (or NULL) */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int*                  norbits,            /**< pointer to number of orbits currently stored in orbits */
   int*                  nvarsinorbits       /**< pointer to number of variables contained in non-trivial orbits */
   );

#ifdef __cplusplus
}
#endif

#endif
