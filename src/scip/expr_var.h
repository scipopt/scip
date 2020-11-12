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

/**@file   expr_var.h
 * @brief  variable expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPR_VAR_H__
#define __SCIP_EXPR_VAR_H__


#include "scip/scip.h"
#include "scip/type_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for variable expression and includes it into SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprHdlrVar(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a variable expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_VAR*             var                 /**< variable to be stored */
   );

#if !1  // FIXME move into cons_nonlinear
/** registers event handler to catch variable events on variable
 *
 * Additionally, the given constraint is stored in the data of the variable-expression.
 * When an event occurs, all stored constraints are notified.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcatchConsExprExprVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   expr,               /**< variable expression */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< expr constraint */
   );

/** unregisters event handler to catch variable events on variable
 *
 * The given constraint is removed from the constraints array in the data of the variable-expression.
 * If this was the last constraint, then the event handler is unregistered for this variable.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdropConsExprExprVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   expr,               /**< variable expression */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< expr constraint */
   );

/** returns whether the variable events on variable are catched */
SCIP_EXPORT
SCIP_Bool SCIPisConsExprExprVarEventCatched(
   SCIP_EXPR*   expr                /**< variable expression */
   );

/** gives number of constraints for which the expression catches bound change events on the variable */
SCIP_EXPORT
int SCIPgetConsExprExprVarNConss(
   SCIP_EXPR*   expr                /**< variable expression */
   );

/** gives constraints for which the expression catches bound change events on the variable */
SCIP_EXPORT
SCIP_CONS** SCIPgetConsExprExprVarConss(
   SCIP_EXPR*   expr                /**< variable expression */
   );
#endif

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPR_VAR_H__ */
