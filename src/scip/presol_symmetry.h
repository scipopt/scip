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

/**@file   presol_symmetry.h
 * @ingroup PRESOLVERS
 * @brief  presovler for storing symmetry information about current problem
 * @author Marc Pfetsch
 * @author Thomas Rehn
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_SYMMETRY_H_
#define __SCIP_PRESOL_SYMMETRY_H_

#include <scip/scip.h>

#include <symmetry/type_symmetry.h>

/** include symmetry presolver */
EXTERN
SCIP_RETCODE SCIPincludePresolSymmetry(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** return symmetry group generators */
EXTERN
SCIP_RETCODE SCIPgetSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< symmetry presolver */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   int*                  nperms,             /**< pointer to store number of permutations */
   int***                perms               /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   );

/** specify symmetry type for which we need symmetries */
EXTERN
void SYMsetSpecRequirement(
   SCIP_PRESOL*          presol,             /**< symmetry presolver */
   SYM_SPEC              type                /**< variable types the callee is interested in */
   );

/** specify symmetry type which symmetry group must fix */
EXTERN
void SYMsetSpecRequirementFixed(
   SCIP_PRESOL*          presol,             /**< symmetry presolver */
   SYM_SPEC              fixedtype           /**< variable types that callee wants to have fixed */
   );

/** whether symmetry should be computed for after presolving */
EXTERN
SCIP_Bool SYMcomputeSymmetryPresolved(
   SCIP_PRESOL*          presol              /**< symmetry presolver */
   );

#endif
