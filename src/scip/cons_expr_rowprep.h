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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_rowprep.h
 * @brief  linear inequalities in preparation
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_ROWPREP_H__
#define __SCIP_CONS_EXPR_ROWPREP_H__

#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sepa.h"
#include "scip/type_sol.h"
#include "scip/type_var.h"
#include "scip/type_cons_expr.h"


/** creates a SCIP_ROWPREP datastructure
 *
 * Initial row represents 0 <= 0.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        rowprep,            /**< buffer to store pointer to rowprep */
   SCIP_SIDETYPE         sidetype,           /**< whether cut will be or lower-equal or larger-equal type */
   SCIP_Bool             local               /**< whether cut will be valid only locally */
);

/** frees a SCIP_ROWPREP datastructure */
SCIP_EXPORT
void SCIPfreeRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        rowprep             /**< pointer that stores pointer to rowprep */
);

/** creates a copy of a SCIP_ROWPREP datastructure */
SCIP_EXPORT
SCIP_RETCODE SCIPcopyRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        target,             /**< buffer to store pointer of rowprep copy */
   SCIP_ROWPREP*         source              /**< rowprep to copy */
);

/** ensures that rowprep has space for at least given number of additional terms
 *
 * Useful when knowing in advance how many terms will be added.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPensureRowprepSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   int                   size                /**< number of additional terms for which to alloc space in rowprep */
);

/** prints a rowprep */
SCIP_EXPORT
void SCIPprintRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be printed */
   FILE*                 file                /**< file to print to, or NULL for stdout */
);

/** prints a rowprep and values in solution */
SCIP_EXPORT
void SCIPprintRowprepSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be printed */
   SCIP_SOL*             sol,                /**< solution for activity */
   FILE*                 file                /**< file to print to, or NULL for stdout */
);

/** adds a term coef*var to a rowprep */
SCIP_EXPORT
SCIP_RETCODE SCIPaddRowprepTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_Real             coef                /**< coefficient to add */
);

/** adds several terms coef*var to a rowprep */
SCIP_EXPORT
SCIP_RETCODE SCIPaddRowprepTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   int                   nvars,              /**< number of terms to add */
   SCIP_VAR**            vars,               /**< variables to add */
   SCIP_Real*            coefs               /**< coefficients to add */
);

/** adds constant value to side of rowprep */
SCIP_EXPORT
void SCIPaddRowprepSide(
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_Real             side                /**< constant value to be added to side */
);

/** adds constant term to rowprep
 *
 * Substracts constant from side.
 */
SCIP_EXPORT
void SCIPaddRowprepConstant(
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_Real             constant            /**< constant value to be added */
);

#ifdef NDEBUG
#define SCIPaddRowprepSide(rowprep, sideadd)  ((rowprep)->side += (sideadd))
#define SCIPaddRowprepConstant(rowprep, constant)  SCIPaddRowprepSide(rowprep, -(constant))
#endif

/** computes violation of cut in a given solution
 *
 * Can return whether the violation value is reliable from a float-point accuracy point of view.
 * The value will not be deemed reliable when its calculation involved the subtraction of large numbers.
 * To be precise, the violation of an inequality \f$ \sum_i a_ix_i \leq b \f$ in a solution \f$x^*\f$ is deemed
 * reliable if \f$ |\sum_i a_ix^*_i - b| >= 2^-50 \max (|b|, \max_i |a_ix^*_i|) \f$.
 */
SCIP_EXPORT
SCIP_Real SCIPgetRowprepViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_SOL*             sol,                /**< solution or NULL for LP solution */
   SCIP_Bool*            reliable            /**< buffer to store whether computed violation is reliable (numerically), or NULL if not of interest */
);

/** computes violation of cut in a given solution and reports whether that value seem numerically reliable
 *
 * @see SCIPgetRowprepViolation
 */
SCIP_EXPORT
SCIP_Bool SCIPisRowprepViolationReliable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_SOL*             sol                 /**< solution or NULL for LP solution */
);

/** Merge terms that use same variable and eliminate zero coefficients.
 *
 * Removes a variable if its bounds have a relative difference of below epsilon.
 * Local bounds are checked for local rows, otherwise global bounds are used.
 * If the bounds are not absolute equal, the bound that relaxes the row is used.
 *
 * Terms are sorted by variable (@see SCIPvarComp) after return.
 */
SCIP_EXPORT
void SCIPmergeRowprepTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep             /**< rowprep to be cleaned up */
);

/* Cleans up and attempts to improve rowprep
 *
 * Drops small or large coefficients if coefrange is too large, if this can be done by relaxing the cut.
 * Scales coefficients up to reach minimal violation, if possible.
 * Scaling is omitted if violation is very small (ROWPREP_SCALEUP_VIOLNONZERO) or
 * maximal coefficient would become huge (ROWPREP_SCALEUP_MAXMAXCOEF).
 * Scales coefficients and side down if they are large and if the minimal violation is still reached.
 * Rounds coefficients close to integral values to integrals, if this can be done by relaxing the cut.
 * Rounds side within epsilon of 0 to 0.0 or +/-1.1*epsilon, whichever relaxes the cut least.
 *
 * After return, the terms in the rowprep will be sorted by absolute value of coefficient, in decreasing order.
 * Thus, the coef.range can be obtained via REALABS(rowprep->coefs[0]) / REALABS(rowprep->coefs[rowprep->nvars-1]) (if nvars>0).
 *
 * success is set to TRUE if and only if the rowprep satisfies the following:
 * - the coef.range is below maxcoefrange
 * - the violation is at least minviol
 * - the violation is reliable or minviol == 0
 * - the absolute value of coefficients are below SCIPs value of infinity
 * - the absolute value of the side is below SCIPs value of infinity
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcleanupRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be cleaned */
   SCIP_SOL*             sol,                /**< solution that we try to cut off, or NULL for LP solution */
   SCIP_Real             maxcoefrange,       /**< maximal allowed coefficients range */
   SCIP_Real             minviol,            /**< minimal absolute violation the row should achieve (w.r.t. sol) */
   SCIP_Real*            viol,               /**< buffer to store absolute violation of cleaned up cut in sol, or NULL if not of interest */
   SCIP_Bool*            success             /**< buffer to store whether cut cleanup was successful, or NULL if not of interest */
);


/* Cleans up and attempts to improve rowprep
 *
 * Drops small or large coefficients if coefrange is too large, if this can be done by relaxing the cut.
 * Scales coefficients and side to have maximal coefficient in [1/maxcoefbound,maxcoefbound].
 * Rounds coefficients close to integral values to integrals, if this can be done by relaxing the cut.
 * Rounds side within epsilon of 0 to 0.0 or +/-1.1*epsilon, whichever relaxes the cut least.
 *
 * After return, the terms in the rowprep will be sorted by absolute value of coefficient, in decreasing order.
 * Thus, the coef.range can be obtained via REALABS(rowprep->coefs[0]) / REALABS(rowprep->coefs[rowprep->nvars-1]) (if nvars>0).
 *
 * success is set to TRUE if and only if the rowprep satisfies the following:
 * - the coef.range is below maxcoefrange
 * - the absolute value of coefficients are below SCIPs value of infinity
 * - the absolute value of the side is below SCIPs value of infinity
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcleanupRowprep2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be cleaned */
   SCIP_SOL*             sol,                /**< solution that we try to cut off, or NULL for LP solution */
   SCIP_Real             maxcoefrange,       /**< maximal allowed coefficients range */
   SCIP_Real             maxcoefbound,       /**< bound on absolute value of largest coefficient */
   SCIP_Bool*            success             /**< buffer to store whether cut cleanup was successful, or NULL if not of interest */
   );

/** Scales up a rowprep to increase coefficients/sides that are within epsilon to an integer value, if possible.
 *
 * Computes the minimal fractionality of all fractional coefficients and the side of the rowprep.
 * If this fractionality is below epsilon, the rowprep is scaled up such that the fractionality exceeds epsilon,
 * if this will not put any coefficient or side above SCIPhugeValue.
 *
 * This does not relax the rowprep.
 * *success is set to TRUE if the resulting rowprep can be turned into a SCIP_ROW, that is,
 * all coefs and the side is below SCIPinfinity and fractionalities are above epsilon.
 * If *success is set to FALSE, then the rowprep will not have been modified.
 *
 * @return The applied scaling factor, if *success is set to TRUE.
 */
SCIP_EXPORT
SCIP_Real SCIPscaleupRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be cleaned */
   SCIP_Real             minscaleup,         /**< minimal factor by which to scale up row, or <= 1.0 if to be ignored */
   SCIP_Bool*            success             /**< buffer to store whether rowprep could be turned into SCIP_ROW without loss, or NULL if not of interest */
);

/** scales a rowprep
 *
 * @return Exponent of actually applied scaling factor, if written as 2^x.
 */
SCIP_EXPORT
int SCIPscaleRowprep(
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be scaled */
   SCIP_Real             factor              /**< suggested scale factor */
);

/** generates a SCIP_ROW from a rowprep */
SCIP_EXPORT
SCIP_RETCODE SCIPgetRowprepRowConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store pointer to new row */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be turned into a row */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
);

/** generates a SCIP_ROW from a rowprep */
SCIP_EXPORT
SCIP_RETCODE SCIPgetRowprepRowCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store pointer to new row */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be turned into a row */
   SCIP_CONS*            cons                /**< constraint */
);

/** generates a SCIP_ROW from a rowprep */
SCIP_EXPORT
SCIP_RETCODE SCIPgetRowprepRowSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store pointer to new row */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be turned into a row */
   SCIP_SEPA*            sepa                /**< separator */
);

#endif /* __SCIP_CONS_EXPR_ROWPREP_H__ */
