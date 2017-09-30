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

/**@file   compute_symmetry.h
 * @brief  interface for symmetry computations
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_COMPUTE_SYMMETRY_H_
#define __SCIP_COMPUTE_SYMMETRY_H_

#include <scip/scip.h>

#include <symmetry/type_symmetry.h>

#ifdef __cplusplus
extern "C" {
#endif

/** return whether symmetry can be computed */
EXTERN
SCIP_Bool SYMcanComputeSymmetry(void);

/** return name of external program used to compute generators */
EXTERN
const char* SYMsymmetryGetName(void);

/** return description of external program used to compute generators */
EXTERN
const char* SYMsymmetryGetDesc(void);

/** compute generators of symmetry group */
EXTERN
SCIP_RETCODE SYMcomputeSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   SCIP_Bool             local,              /**< Use local variable bounds? */
   int                   npermvars,          /**< number of variables for permutations */
   SCIP_VAR**            permvars,           /**< variables on which permutations act */
   SYM_MATRIXDATA*       matrixdata,         /**< data for MIP matrix */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms               /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   );

#ifdef __cplusplus
}
#endif

#endif
