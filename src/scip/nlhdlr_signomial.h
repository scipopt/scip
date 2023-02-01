/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_signomial.h
 * @ingroup NLHDLRS
 * @brief  signomial nonlinear handler
 * @author Liding Xu
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLHDLR_SIGNOMIAL_H__
#define __SCIP_NLHDLR_SIGNOMIAL_H__

#include "scip/scip.h"
#include "scip/pub_nlhdlr.h"

#ifdef __cplusplus
extern "C" {
#endif


/** includes signomial nonlinear handler to nonlinear constraint handler
 *
 * @ingroup NlhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrSignomial(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**@addtogroup NLHDLRS
 * @{
 *
 * @name Signomial nonlinear handler
 *
 * This nonlinear handler detects and collects signomial terms and provides specialized propagation and estimation functionality.
 *
 * @{
 */

SCIP_EXPORT
SCIP_EXPR** SCIPgetExprsSignomial(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   );

/** returns an array of nonlinear handler expressions data of expressions that have been detected by the signomial nonlinear handler */
SCIP_EXPORT
SCIP_NLHDLREXPRDATA** SCIPgetExprsdataSignomial(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   );

/** returns the total number of expressions that have been detected by the signomial nonlinear handler */
SCIP_EXPORT
int SCIPgetNExprsSignomial(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   );

/** adds a globally valid inequality of the form \f$\text{xcoef}\cdot x \leq \text{ycoef} \cdot y + \text{constant}\f$ to a product expression of the form \f$x\cdot y\f$ */
SCIP_EXPORT
SCIP_RETCODE SCIPaddIneqSignomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_EXPR*            expr,               /**< product expression */
   SCIP_Real             xcoef,              /**< x coefficient */
   SCIP_Real             ycoef,              /**< y coefficient */
   SCIP_Real             constant,           /**< constant part */
   SCIP_Bool*            success             /**< buffer to store whether inequality has been accepted */
   );



/** @}
  * @}
  */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLHDLR_SIGNOMIAL_H__ */
