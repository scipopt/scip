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

/**@file   cons_symmetry.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for storing symmetry information about current problem
 * @author Marc Pfetsch
 * @author Thomas Rehn
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_SYMMETRY_H_
#define __SCIP_CONS_SYMMETRY_H_

#include <scip/scip.h>

#include <symmetry/type_symmetry.h>

/** include symmetry constraint handler */
extern
SCIP_RETCODE SCIPincludeConshdlrSymmetry(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** create symmetries constraint */
extern
SCIP_RETCODE SCIPcreateConsSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name                /**< name of constraint */
   );

/** return symmetry group generators */
extern
SCIP_RETCODE SCIPgetSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< symmetry constraint handler */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   SCIP_HASHMAP**        permvarmap,         /**< map of variables to indices in permvars array */
   int*                  nperms,             /**< pointer to store number of permutations */
   int***                perms               /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   );

/** specify symmetry type for which we need symmetries */
extern
void SYMsetSpecRequirement(
   SCIP_CONSHDLR*        conshdlr,           /**< symmetry constraint handler */
   SYM_SPEC              type                /**< variable types the callee is interested in */
   );

/** specify symmetry type which symmetry group must fix */
extern
void SYMsetSpecRequirementFixed(
   SCIP_CONSHDLR*        conshdlr,           /**< symmetry constraint handler */
   SYM_SPEC              fixedtype           /**< variable types that callee wants to have fixed */
   );

/** whether symmetry should be computed for presolved system */
extern
SCIP_Bool SYMdetectSymmetryPresolved(
   SCIP_CONSHDLR*        conshdlr            /**< symmetry constraint handler */
   );

#endif
