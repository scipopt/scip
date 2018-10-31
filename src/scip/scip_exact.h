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

/**@file   scip_sol.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for solutions
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_EXACT_H__
#define __SCIP_SCIP_EXACT_H__


#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_heur.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"
#include "scip/type_var.h"

/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons).
 * Additionally, the internal "set.h" is included, such that the defines in set.h are
 * available in optimized mode.
 */
#ifdef NDEBUG
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"
#include "scip/var.h"
#include "scip/cons.h"
#include "scip/solve.h"
#include "scip/debug.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** returns whether the solution process is arithmetically exact, i.e., not subject to roundoff errors
 *
 *  @note This feature is not supported yet!
 *
 *  @return Returns TRUE if \SCIP is exact solving mode, otherwise FALSE
 */
EXTERN
SCIP_Bool SCIPisExactSolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* CERT: we probably want to use an FP approximation since we do not use bounding method 'n' */
/** returns whether the floating point problem should be a relaxation of the original problem (instead of an approximation);
 *  only relevant for solving the problem provably correct
 */
SCIP_Bool SCIPuseFPRelaxation(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns which method is used for computing truely valid dual bounds at the nodes ('n'eumaier and shcherbina,
 *  'v'erify LP basis, 'r'epair LP basis, 'p'roject and scale, 'e'xact LP,'i'nterval neumaier and shcherbina,
 *  e'x'act neumaier and shcherbina, 'a'utomatic); only relevant for solving the problem provably correct
 */
char SCIPdualBoundMethod(
   SCIP*                 scip                /**< SCIP data structure */
   );


#ifdef __cplusplus
}
#endif

#endif