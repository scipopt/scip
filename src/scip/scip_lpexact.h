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
   SCIP_ROW*             fprowrhs,           /**< rhs-part of fp-relaxation of this row if necessary, NULL otherwise */
   SCIP_Rational*        lhs,                /**< left hand side of row */
   SCIP_Rational*        rhs,                /**< right hand side of row */
   SCIP_Bool             isfprelaxable       /**< is it possible to create an fp relaxation of this row? */
   );

/** creates and captures an exact LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT**       row,                /**< pointer to row */
   SCIP_ROW*             fprow,              /**< corresponding fp approximation/relaxation */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COLEXACT**       cols,               /**< array with columns of row entries */
   SCIP_Rational**       vals,               /**< array with coefficients of row entries */
   SCIP_Rational*        lhs,                /**< left hand side of row */
   SCIP_Rational*        rhs,                /**< right hand side of row */
   SCIP_Bool             isfprelaxable       /**< is it possible to make fp-relaxation of this row */
   );

/** creates and captures an exact LP row from an existing fp row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateRowExactFromRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             fprow
                /**< corresponding fp approximation/relaxation */
   );

/** generates two fprows that are a relaxation of the exact row wrt the lhs/rhs, respectively
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgenerateFpRowsFromRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< SCIP exact row */
   SCIP_ROW*             rowlhs,             /**< fp row-relaxation wrt lhs */
   SCIP_ROW*             rowrhs,             /**< fp row-relaxation wrt rhs */
   SCIP_Bool*            onerowrelax,        /**< is one row enough to represent the exact row */
   SCIP_Bool*            hasfprelax          /**< is it possible to generate relaxations at all for this row? */
   );

/** increases usage counter of exact LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcaptureRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row                 /**< row to capture */
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

/** changes left hand side of exact LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgRowExactLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_Rational*        lhs                 /**< new left hand side */
   );

/** changes right hand side of exact LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgRowExactRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_Rational*        rhs                 /**< new right hand side */
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
SCIP_RETCODE SCIPaddVarsToRowExact(
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

/** returns the activity of a row for the given primal solution with running error analysis
 *
 *  @return the activitiy of a row for the given primal solution and the error bound of the activity; returns true on success
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_Bool SCIPgetRowSolActivityWithErrorboundExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            activity,           /**< the approximate activity */
   SCIP_Real*            errorbound          /**< the error bound */
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
SCIP_RETCODE SCIPgetRowSolActivityExact(
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
SCIP_RETCODE SCIPgetRowSolFeasibilityExact(
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
SCIP_RETCODE SCIPprintRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** returns whether the exact lp was solved */
SCIP_EXPORT
SCIP_Bool SCIPlpExactIsSolved(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets solution status of current exact LP
  *
  *  @return the solution status of current exact LP.
  *
  *  @pre This method can be called if @p scip is in one of the following stages:
  *       - \ref SCIP_STAGE_SOLVING
  *
  *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
  */
SCIP_EXPORT
SCIP_LPSOLSTAT SCIPgetLPExactSolstat(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets objective value of current exact LP (which is the sum of column and loose objective value)
 *
 *  @return the objective value of current LP (which is the sum of column and loose objective value).
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This method returns the objective value of the current LP solution, which might be primal or dual infeasible
 *        if a limit was hit during solving. It must not be used as a dual bound if the exact LP solution status
 *        returned by SCIPgetLPExactSolstat() is SCIP_LPSOLSTAT_ITERLIMIT or SCIP_LPSOLSTAT_TIMELIMIT.
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
void SCIPgetLPExactObjval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Rational*        result              /**< result pointer */
   );

/** changes variable's lower bound in current exact dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgVarLbExactDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Rational*        newbound            /**< new value for bound */
   );

/** changes variable's upper bound in current exact dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgVarUbExactDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Rational*        newbound            /**< new value for bound */
   );

/** solves the exact LP of the current dive; no separation or pricing is applied
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsolveExactDiveLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool*            cutoff              /**< pointer to store whether the diving LP was infeasible or the objective
                                              *   limit was reached (or NULL, if not needed) */
   );

/** initiates exact LP diving, making methods SCIPchgVarObjExactDive(), SCIPchgVarLbExactDive(), and SCIPchgVarUbExactDive() available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note In parallel to exact LP diving, this method also starts the regular LP diving mode by calling SCIPstartDive().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPstartExactDive(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** checks if exact diving mode is possible at this point in time
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note In parallel to exact LP diving, this method also starts the regular LP diving mode by calling SCIPstartDive().
 */
SCIP_EXPORT
SCIP_Bool SCIPisExactDivePossible(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether we are in exact diving mode
 *
 *  @return whether we are in exact diving mode.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_Bool SCIPinExactDive(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** quits exact LP diving and resets bounds and objective values of columns to the current node's values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPendExactDive(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** writes current exact LP to a file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPwriteLPexact(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   );

#ifdef __cplusplus
}
#endif

#endif
