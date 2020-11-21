/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_abs.h
 * @brief  absolute expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPR_ABS_H__
#define __SCIP_EXPR_ABS_H__


#include "scip/scip.h"
#include "scip/type_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates a absolute expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprAbs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_DECL_EXPR_OWNERDATACREATE((*ownerdatacreate)), /**< function to call to create ownerdata */
   SCIP_EXPR_OWNERDATACREATEDATA* ownerdatacreatedata  /**< data to pass to ownerdatacreate */
   );

/** creates the handler for absolute expression and includes it into the expression constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprHdlrAbs(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPR_ABS_H__ */
