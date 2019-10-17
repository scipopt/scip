/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   var.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for  data of problem variables
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_VAREX_H__
#define __SCIP_PUB_VAREX_H__


#include "blockmemshell/memory.h"
#include "scip/def.h"
#include "scip/rational.h"
#include "scip/type_branch.h"
#include "scip/type_cons.h"
#include "scip/type_event.h"
#include "scip/type_history.h"
#include "scip/type_implics.h"
#include "scip/type_lp.h"
#include "scip/type_message.h"
#include "scip/type_misc.h"
#include "scip/type_primal.h"
#include "scip/type_prob.h"
#include "scip/type_prop.h"
#include "scip/type_relax.h"
#include "scip/type_reopt.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_set.h"
#include "scip/type_sol.h"
#include "scip/type_stat.h"
#include "scip/type_tree.h"
#include "scip/type_var.h"

#ifndef NDEBUG
#include "scip/struct_var.h"
#else
#include "scip/event.h"
#include "scip/pub_history.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
/*
 * getter methods for  varbounds
 */

/** gets original lower bound of original problem variable (i.e. the bound set in problem creation) */
SCIP_Rational* SCIPvarGetLbOriginalExact(
   SCIP_VAR*             var                 /**< original problem variable */
   );

/** gets original upper bound of original problem variable (i.e. the bound set in problem creation) */
SCIP_Rational* SCIPvarGetUbOriginalExact(
   SCIP_VAR*             var                 /**< original problem variable */
   );

/** gets global lower bound of variable */
SCIP_Rational* SCIPvarGetLbGlobalExact(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets global upper bound of variable */
SCIP_Rational* SCIPvarGetUbGlobalExact(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets best global bound of variable with respect to the objective function */
SCIP_Rational* SCIPvarGetBestBoundGlobalExact(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets worst global bound of variable with respect to the objective function */
SCIP_Rational* SCIPvarGetWorstBoundGlobalExact(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets current lower bound of variable */
SCIP_Rational* SCIPvarGetLbLocalExact(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets current upper bound of variable */
SCIP_Rational* SCIPvarGetUbLocalExact(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets best local bound of variable with respect to the objective function */
SCIP_Rational* SCIPvarGetBestBoundLocalExact(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets worst local bound of variable with respect to the objective function */
SCIP_Rational* SCIPvarGetWorstBoundLocalExact(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets type (lower or upper) of best bound of variable with respect to the objective function */
SCIP_BOUNDTYPE SCIPvarGetBestBoundTypeExact(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets type (lower or upper) of worst bound of variable with respect to the objective function */
SCIP_BOUNDTYPE SCIPvarGetWorstBoundTypeExact(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets objective function value of variable */
EXTERN
SCIP_Rational* SCIPvarGetObjExact(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** retransforms given variable, scalar and constant to the corresponding original variable, scalar
 *  and constant, if possible; if the retransformation is impossible, NULL is returned as variable
 */
EXTERN
SCIP_RETCODE SCIPvarGetOrigvarSumExact(
   SCIP_VAR**            var,                /**< pointer to problem variable x in sum a*x + c */
   SCIP_Rational*        scalar,             /**< pointer to scalar a in sum a*x + c */
   SCIP_Rational*        constant            /**< pointer to constant c in sum a*x + c */
   );

#ifdef __cplusplus
}
#endif

#endif