
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_lpexact.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for the LP relaxation, rows and columns
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_LPEXACT_H__
#define __SCIP_SCIP_LPEXACT_H__


#include "lpi/type_lpi.h"
#include "scip/def.h"
#include "scip/rational.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_lpexact.h"
#include "scip/type_misc.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sepa.h"
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

/** creates and captures an LP row without any coefficients from a constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateEmptyRowConsExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT**       rowexact,           /**< pointer to row */
   SCIP_ROW*             fprow,              /**< pointer to fp-row that corresponds to this row */
   SCIP_Rational*        lhs,                /**< left hand side of row */
   SCIP_Rational*        rhs                 /**< right hand side of row */
   );

/** decreases usage counter of LP row, and frees memory if necessary
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_EXPORT
SCIP_RETCODE SCIPreleaseRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT**       row                 /**< pointer to LP row */
   );


/** resolves variables to columns and adds them with the coefficients to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @attention If a coefficients absolute value is below the SCIP epsilon tolerance, the variable with its value is not added.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddVarsToRowEx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   int                   nvars,              /**< number of variables to add to the row */
   SCIP_VAR**            vars,               /**< problem variables to add */
   SCIP_Rational**       vals                /**< values of coefficients */
   );

/** returns the activity of a row in the last LP or pseudo solution
 *
 *  @return the activity of a row in the last LP or pseudo solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
void SCIPgetRowActivityExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_Rational*        result              /**< result pointer */
   );

/** returns the feasibility of a row in the last LP or pseudo solution
 *
 *  @return the feasibility of a row in the last LP or pseudo solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
void SCIPgetRowFeasibilityExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_Rational*        result              /**< result pointer */
   );

/** returns the activity of a row for the given primal solution
 *
 *  @return the activitiy of a row for the given primal solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
void SCIPgetRowSolActivityExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             useexact,           /**< true if sol should be considered instead of sol */
   SCIP_Rational*        result              /**< result pointer */
   );

/** returns the feasibility of a row for the given primal solution
 *
 *  @return the feasibility of a row for the given primal solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
void SCIPgetRowSolFeasibilityExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Rational*        result              /**< result pointer */
   );

/** output exact row to file stream via the message handler system
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_EXPORT
SCIP_RETCODE SCIPprintRowex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** returns whether the exact lp was solved */
SCIP_EXPORT
SCIP_Bool SCIPlpExactIsSolved(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif