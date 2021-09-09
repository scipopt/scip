/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_quotient.h
 * @ingroup NLHDLRS
 * @brief  quotient nonlinear handler
 * @author Benjamin Mueller
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLHDLR_QUOTIENT_H__
#define __SCIP_NLHDLR_QUOTIENT_H__

#include "scip/scip.h"
#include "scip/pub_nlhdlr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes quotient nonlinear handler in nonlinear constraint handler
 *
 * @ingroup NlhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrQuotient(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**@addtogroup NLHDLRS
 * @{
 *
 * @name Quotient nonlinear handler
 *
 * This nonlinear handler detects quotients and provides specialized propagation and estimation functionality.
 *
 * @{
 */

/** @}
  * @}
  */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLHDLR_QUOTIENT_H__ */
