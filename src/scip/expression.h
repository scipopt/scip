/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: expression.h,v 1.2 2010/05/06 16:12:56 bzfviger Exp $"

/**@file   expression.h
 * @brief  more methods for expressions and expression trees
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 *
 * This file contains methods for handling and manipulating expressions and expression trees
 * that are SCIP specific and thus not included in nlpi/expression.*
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPRESSION_H__
#define __SCIP_EXPRESSION_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_var.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"
#include "nlpi/expression.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns variables of expression tree */
SCIP_VAR** SCIPexprtreeGetVars(
   SCIP_EXPRTREE*        tree                /**< expression tree */
);

/** stores array of variables in expression tree */
SCIP_RETCODE SCIPexprtreeSetVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars                /**< variables */
);

/** evaluates an expression tree for a primal solution or LP solution */
SCIP_RETCODE SCIPexprtreeEvalSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_SOL*             sol,                /**< a solution, or NULL for current LP solution */
   SCIP_Real*            val                 /**< buffer to store value */
);

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPRESSION_H_ */
