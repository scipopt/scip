/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_symmetry.h
 * @ingroup PRESOLVERS
 * @brief  presolver for storing symmetry information about current problem
 * @author Marc Pfetsch
 * @author Thomas Rehn
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_SYMMETRY_H_
#define __SCIP_PRESOL_SYMMETRY_H_

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <symmetry/type_symmetry.h>

/** include symmetry presolver */
EXTERN
SCIP_RETCODE SCIPincludePresolSymmetry(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** return symmetry group generators */
EXTERN
SCIP_RETCODE SCIPgetGeneratorsSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SPEC              symspecrequire,     /**< symmetry specification for which we need to compute symmetries */
   SYM_SPEC              symspecrequirefixed,/**< symmetry specification of variables which must be fixed by symmetries */
   SCIP_Bool             recompute,          /**< Have symmetries already been computed? */
   SCIP_Bool             transposedperms,    /**< whether permutations should be computed/stored in transposed form (i.e., npermvars x nperms) */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   int*                  nperms,             /**< pointer to store number of permutations */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars or transposed) matrix */
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of group size (or NULL) */
   SCIP_Bool*            binvaraffected,     /**< pointer to store whether binary variables are affected */
   int**                 components,         /**< pointer to store components of symmetry group (or NULL) */
   int**                 componentbegins,    /**< pointer to store begin positions of components in components array (or NULL) */
   int**                 vartocomponent,     /**< pointer to store assignment from variable to its component (or NULL) */
   int*                  ncomponents         /**< pointer to store number of components (or NULL) */
   );

/** return objective coefficients of permuted variables at time of symmetry computation */
EXTERN
SCIP_RETCODE SCIPgetPermvarsObjSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           permvarsobj         /**< pointer to store objective coefficients of permuted variables (NULL if not available) */
   );

/* block component of symmetry group to be considered by symmetry handling routines */
EXTERN
SCIP_RETCODE SCIPblockSymmetryComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   i                   /**< index of component to block */
   );

#ifdef __cplusplus
}
#endif

#endif
