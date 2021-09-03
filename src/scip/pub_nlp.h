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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_nlp.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for NLP management
 * @author Thorsten Gellermann
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_NLP_H__
#define __SCIP_PUB_NLP_H__

#include <stdio.h>

#include "scip/def.h"
#include "scip/type_message.h"
#include "blockmemshell/memory.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_nlp.h"
#include "scip/type_var.h"
#include "scip/type_sol.h"
#include "scip/type_expr.h"
#include "scip/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicNLRowMethods
 *
 * @{
 */

/** gets constant */
SCIP_EXPORT
SCIP_Real SCIPnlrowGetConstant(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets number of variables of linear part */
SCIP_EXPORT
int SCIPnlrowGetNLinearVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets array with variables of linear part */
SCIP_EXPORT
SCIP_VAR** SCIPnlrowGetLinearVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets array with coefficients in linear part */
SCIP_EXPORT
SCIP_Real* SCIPnlrowGetLinearCoefs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets expression */
SCIP_EXPORT
SCIP_EXPR* SCIPnlrowGetExpr(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** returns the left hand side of a nonlinear row */
SCIP_EXPORT
SCIP_Real SCIPnlrowGetLhs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** returns the right hand side of a nonlinear row */
SCIP_EXPORT
SCIP_Real SCIPnlrowGetRhs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** returns the curvature of a nonlinear row */
SCIP_EXPORT
SCIP_EXPRCURV SCIPnlrowGetCurvature(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** sets the curvature of a nonlinear row */
SCIP_EXPORT
void SCIPnlrowSetCurvature(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_EXPRCURV         curvature           /**< curvature of NLP row */
   );

/** returns the name of a nonlinear row */
SCIP_EXPORT
const char* SCIPnlrowGetName(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets position of a nonlinear row in current NLP, or -1 if not in NLP */
SCIP_EXPORT
int SCIPnlrowGetNLPPos(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** returns TRUE iff row is member of current NLP */
SCIP_EXPORT
SCIP_Bool SCIPnlrowIsInNLP(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets the dual NLP solution of a nlrow
 *
 * for a ranged constraint, the dual value is positive if the right hand side is active and negative if the left hand side is active
 */
SCIP_EXPORT
SCIP_Real SCIPnlrowGetDualsol(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_PUB_NLP_H__ */
