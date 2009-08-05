/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_exactlp.h,v 1.1.2.2 2009/08/05 10:10:27 bzfwolte Exp $"

/**@file   cons_exactlp.h
 * @brief  constraint handler for exactlp constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXACTLP_H__
#define __SCIP_CONS_EXACTLP_H__


#include <mpfr.h> /* mpfr.h has to be included before gmp.h */ /* todo: only necessary because of gmp<->fp functions (which maybe move) ?????? */
#include <gmp.h>
#include "scip/scip.h"
#include "scip/lpiex.h"


/** returns value treated as negative infinite in exactlp constraint handler */
extern
const mpq_t* negInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< exactlp constraint handler data */
   );

/** checks if value is treated as positive infinite in exactlp constraint handler */
extern
const mpq_t* posInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< exactlp constraint handler data */
   );

/** checks if value is treated as negative infinite in exactlp constraint handler */
extern
SCIP_Bool isNegInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   mpq_t                 val                 /**< value to be compared against infinity */
   );

/** checks if value is treated as positive infinite in exactlp constraint handler */
extern
SCIP_Bool isPosInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   mpq_t                 val                 /**< value to be compared against infinity */
   );

/** returns whether given rational number can be stored as FP number without roundinf errors */
extern
SCIP_Bool mpqIsReal(
   SCIP*                 scip,               /**< SCIP data structure */
   mpq_t                 val                 /**< given rational number */
   );

/** converts given rational number into an FP number; uses given rounding mode during conversion 
 * (should be used to construct an FP relaxation of a constraint) 
 */
extern
SCIP_Real mpqGetRealRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   mpq_t                 val,                /**< given rational number */
   mp_rnd_t              roundmode           /**< rounding mode to be used for the conversion */
   );

/** converts given rational number into an FP number; uses default rounding mode during conversion 
 * (should be used to construct an FP approximation of a constraint) 
 */
extern
SCIP_Real mpqGetRealApprox(
   SCIP*                 scip,               /**< SCIP data structure */
   mpq_t                 val                 /**< given rational number */
   );

/** creates the handler for exactlp constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrExactlp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a exactlp constraint */
extern
SCIP_RETCODE SCIPcreateConsExactlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   nvars,              /**< number of variables */
   mpq_t*                obj,                /**< objective function values of variables */
   mpq_t*                lb,                 /**< lower bounds of variables */
   mpq_t*                ub,                 /**< upper bounds of variables */
   int                   nconss,              /**< number of constraints */
   mpq_t*                lhs,                /**< left hand sides of constraints */
   mpq_t*                rhs,                /**< right hand sides of constraints */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each variable in ind- and val-array */
   int*                  ind,                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val,                /**< values of nonzero constraint matrix entries */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** returns value treated as positive infinite in exactlp constraint handler */
void SCIPpositivInfinityExactlp(
   SCIP_CONSHDLR*        conshdlr,           /**< exactlp constraint handler */
   mpq_t                 posinfinity         /**< pointer to store positive infinity */
   );

/** returns value treated as negative infinite in exactlp constraint handler */
void SCIPnegativInfinityExactlp(
   SCIP_CONSHDLR*        conshdlr,           /**< exactlp constraint handler */
   mpq_t                 neginfinity         /**< pointer to store negative infinity */
   );

#endif
