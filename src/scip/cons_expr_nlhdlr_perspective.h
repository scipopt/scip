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

/**@file   cons_expr_nlhdlr_erspective.h
 * @brief  PERSPECTIVE nonlinear handler
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_NLHDLR_PERSPECTIVE_H__
#define __SCIP_CONS_EXPR_NLHDLR_PERSPECTIVE_H__

#include "scip/scip.h"
#include "scip/type_cons_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data structure to store information of a semicontinuous variable
 */
struct SCIP_SCVarData
{
   SCIP_VAR*             bvar;           /**< the binary variable on which the variable domain depends */
   SCIP_Real             val0;            /**< var value when bvar = 0 */
};
typedef struct SCIP_SCVarData SCIP_SCVARDATA;

/** includes perspective nonlinear handler to consexpr */
EXTERN
SCIP_RETCODE SCIPincludeConsExprNlhdlrPerspective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_NLHDLR_DEFAULT_H__ */
