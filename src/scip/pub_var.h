/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_var.h,v 1.14 2004/04/30 11:16:25 bzfpfend Exp $"

/**@file   pub_var.h
 * @brief  public methods for problem variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_VAR_H__
#define __PUB_VAR_H__


#include "def.h"
#include "type_retcode.h"
#include "type_history.h"
#include "type_var.h"
#include "type_cons.h"

#ifdef NDEBUG
#include "struct_var.h"
#include "history.h"
#endif



/*
 * methods for variables 
 */

/** increases lock number for rounding down by one; tells variable, that rounding its value down will make the
 *  solution infeasible
 */
extern
void SCIPvarLockDown(
   VAR*             var                 /**< problem variable */
   );

/** increases lock number for rounding up by one; tells variable, that rounding its value up will make the
 *  solution infeasible
 */
extern
void SCIPvarLockUp(
   VAR*             var                 /**< problem variable */
   );

/** increases lock number for rounding down and up by one; tells variable, that rounding value in either direction will
 *  make the solution infeasible
 */
extern
void SCIPvarLockBoth(
   VAR*             var                 /**< problem variable */
   );

/** declares that rounding down the given variable would destroy the feasibility of the given constraint;
 *  locks the roundings of the variable corresponding to the lock status of the constraint and its negation
 */
extern
void SCIPvarLockDownCons(
   VAR*             var,                /**< problem variable */
   CONS*            cons                /**< constraint */
   );

/** declares that rounding up the given variable would destroy the feasibility of the given constraint;
 *  locks the roundings of the variable corresponding to the lock status of the constraint and its negation
 */
extern
void SCIPvarLockUpCons(
   VAR*             var,                /**< problem variable */
   CONS*            cons                /**< constraint */
   );

/** declares that rounding the given variable in any direction would destroy the feasibility of the given constraint;
 *  locks the roundings of the variable corresponding to the lock status of the constraint and its negation
 */
extern
void SCIPvarLockBothCons(
   VAR*             var,                /**< problem variable */
   CONS*            cons                /**< constraint */
   );

/** increases lock number for roundings of variable; tells variable, that rounding value in a direction set to
 *  a positive value will make the solution infeasible
 */
extern
void SCIPvarLock(
   VAR*             var,                /**< problem variable */
   int              nlocksdown,         /**< increase in number of rounding down locks */
   int              nlocksup            /**< increase in number of rounding up locks */
   );

/** decreases lock number for rounding down by one; cancels a prior SCIPvarLockDown() */
extern
void SCIPvarUnlockDown(
   VAR*             var                 /**< problem variable */
   );

/** decreases lock number for rounding up by one; cancels a prior SCIPvarLockUp() */
extern
void SCIPvarUnlockUp(
   VAR*             var                 /**< problem variable */
   );

/** decreases lock number for rounding down and up by one; cancels a prior SCIPvarLockBoth() */
extern
void SCIPvarUnlockBoth(
   VAR*             var                 /**< problem variable */
   );

/** declares that rounding down the given variable would no longer destroy the feasibility of the given constraint;
 *  unlocks the roundings of the variable corresponding to the lock status of the constraint and its negation
 */
extern
void SCIPvarUnlockDownCons(
   VAR*             var,                /**< problem variable */
   CONS*            cons                /**< constraint */
   );

/** declares that rounding up the given variable would no longer destroy the feasibility of the given constraint;
 *  unlocks the roundings of the variable corresponding to the lock status of the constraint and its negation
 */
extern
void SCIPvarUnlockUpCons(
   VAR*             var,                /**< problem variable */
   CONS*            cons                /**< constraint */
   );

/** declares that rounding the given variable in any direction would no longer destroy the feasibility of the given
 *  constraint; unlocks the roundings of the variable corresponding to the lock status of the constraint and its negation
 */
extern
void SCIPvarUnlockBothCons(
   VAR*             var,                /**< problem variable */
   CONS*            cons                /**< constraint */
   );

/** decreases lock number for roundings of variable; cancels a prior call to SCIPvarLock() */
extern
void SCIPvarUnlock(
   VAR*             var,                /**< problem variable */
   int              nunlocksdown,       /**< decrease in number of rounding down locks */
   int              nunlocksup          /**< decrease in number of rounding up locks */
   );

/** gets number of locks for rounding down */
extern
int SCIPvarGetNLocksDown(
   VAR*             var                 /**< problem variable */
   );

/** gets number of locks for rounding up */
extern
int SCIPvarGetNLocksUp(
   VAR*             var                 /**< problem variable */
   );

/** is it possible, to round variable down and stay feasible? */
extern
Bool SCIPvarMayRoundDown(
   VAR*             var                 /**< problem variable */
   );

/** is it possible, to round variable up and stay feasible? */
extern
Bool SCIPvarMayRoundUp(
   VAR*             var                 /**< problem variable */
   );

/** compares the index of two variables, returns -1 if first is smaller than, and +1 if first is greater than second
 *  variable index; returns 0 if both indices are equal, which means both variables are equal
 */
extern
int SCIPvarCmp(
   VAR*             var1,               /**< first problem variable */
   VAR*             var2                /**< second problem variable */
   );

/** gets corresponding active problem variable of a variable; returns NULL for fixed variables */
extern
VAR* SCIPvarGetProbvar(
   VAR*             var                 /**< problem variable */
   );

/** transforms given variable, boundtype and bound to the corresponding active variable values */
extern
RETCODE SCIPvarGetProbvarBound(
   VAR**            var,                /**< pointer to problem variable */
   Real*            bound,              /**< pointer to bound value to transform */
   BOUNDTYPE*       boundtype           /**< pointer to type of bound: lower or upper bound */
   );

/** transforms given variable, scalar and constant to the corresponding active variable, scalar and constant;
 *  if the variable resolves to a fixed variable, the returned variable will be NULL, "scalar" will be 0.0 and
 *  the value of the sum will be stored in "constant"
 */
extern
RETCODE SCIPvarGetProbvarSum(
   VAR**            var,                /**< pointer to problem variable x in sum a*x + c */
   Real*            scalar,             /**< pointer to scalar a in sum a*x + c */
   Real*            constant            /**< pointer to constant c in sum a*x + c */
   );

/** returns the number of times, a bound of the variable was changed in given direction due to branching */
extern
Longint SCIPvarGetNBranchings(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** returns the number of inferences branching on this variable in given direction triggered */
extern
Longint SCIPvarGetNInferences(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** returns the average depth of bound changes in given direction due to branching on the variable */
extern
Real SCIPvarGetAvgBranchdepth(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction */
   );


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** get name of variable */
extern
const char* SCIPvarGetName(
   VAR*             var                 /**< problem variable */
   );

/** gets status of variable */
extern
VARSTATUS SCIPvarGetStatus(
   VAR*             var                 /**< problem variable */
   );

/** returns whether the variable belongs to the original problem */
extern
Bool SCIPvarIsOriginal(
   VAR*             var                 /**< problem variable */
   );

/** returns whether the variable belongs to the transformed problem */
extern
Bool SCIPvarIsTransformed(
   VAR*             var                 /**< problem variable */
   );

/** returns whether the variable was created by negation of a different variable */
extern
Bool SCIPvarIsNegated(
   VAR*             var                 /**< problem variable */
   );

/** gets type of variable */
extern
VARTYPE SCIPvarGetType(
   VAR*             var                 /**< problem variable */
   );

/** returns whether variable's column should be present in the initial root LP */
extern
Bool SCIPvarIsInitial(
   VAR*             var                 /**< problem variable */
   );

/** returns whether variable's column is removeable from the LP (due to aging or cleanup) */
extern
Bool SCIPvarIsRemoveable(
   VAR*             var                 /**< problem variable */
   );

/** returns whether variable is an active (neither fixed nor aggregated) variable */
extern
Bool SCIPvarIsActive(
   VAR*             var                 /**< problem variable */
   );

/** gets unique index of variable */
extern
int SCIPvarGetIndex(
   VAR*             var                 /**< problem variable */
   );

/** gets position of variable in problem, or -1 if variable is not active */
extern
int SCIPvarGetProbindex(
   VAR*             var                 /**< problem variable */
   );

/** gets transformed variable of ORIGINAL variable */
extern
VAR* SCIPvarGetTransVar(
   VAR*             var                 /**< problem variable */
   );

/** gets column of COLUMN variable */
extern
COL* SCIPvarGetCol(
   VAR*             var                 /**< problem variable */
   );

/** gets aggregation variable y of an aggregated variable x = a*y + c */
extern
VAR* SCIPvarGetAggrVar(
   VAR*             var                 /**< problem variable */
   );

/** gets aggregation scalar a of an aggregated variable x = a*y + c */
extern
Real SCIPvarGetAggrScalar(
   VAR*             var                 /**< problem variable */
   );

/** gets aggregation constant c of an aggregated variable x = a*y + c */
extern
Real SCIPvarGetAggrConstant(
   VAR*             var                 /**< problem variable */
   );

/** gets number n of aggregation variables of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
extern
int SCIPvarGetMultaggrNVars(
   VAR*             var                 /**< problem variable */
   );

/** gets vector of aggregation variables y of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
extern
VAR** SCIPvarGetMultaggrVars(
   VAR*             var                 /**< problem variable */
   );

/** gets vector of aggregation scalars a of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
extern
Real* SCIPvarGetMultaggrScalars(
   VAR*             var                 /**< problem variable */
   );

/** gets aggregation constant c of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
extern
Real SCIPvarGetMultaggrConstant(
   VAR*             var                 /**< problem variable */
   );

/** gets the negation of the given variable; may return NULL, if no negation is existing yet */
extern
VAR* SCIPvarGetNegatedVar(
   VAR*             var                 /**< negated problem variable */
   );

/** gets the negation variable x of a negated variable x' = offset - x */
extern
VAR* SCIPvarGetNegationVar(
   VAR*             var                 /**< negated problem variable */
   );

/** gets the negation offset of a negated variable x' = offset - x */
extern
Real SCIPvarGetNegationConstant(
   VAR*             var                 /**< negated problem variable */
   );

/** gets objective function value of variable */
extern
Real SCIPvarGetObj(
   VAR*             var                 /**< problem variable */
   );

/** gets global lower bound of variable */
extern
Real SCIPvarGetLbGlobal(
   VAR*             var                 /**< problem variable */
   );

/** gets global upper bound of variable */
extern
Real SCIPvarGetUbGlobal(
   VAR*             var                 /**< problem variable */
   );

/** gets current lower bound of variable */
extern
Real SCIPvarGetLbLocal(
   VAR*             var                 /**< problem variable */
   );

/** gets current upper bound of variable */
extern
Real SCIPvarGetUbLocal(
   VAR*             var                 /**< problem variable */
   );

/** gets inference variable of variable (variable that was assigned: parent of var, or var itself), or NULL */
extern
VAR* SCIPvarGetInferVar(
   VAR*             var                 /**< problem variable */
   );

/** gets inference constraint of variable (constraint that deduced the current assignment), or NULL */
extern
CONS* SCIPvarGetInferCons(
   VAR*             var                 /**< problem variable */
   );

/** gets user information for inference to help resolving the conflict */
extern
int SCIPvarGetInferInfo(
   VAR*             var                 /**< problem variable */
   );

/** gets depth level, where the binary variable was fixed, or -1 if unfixed */
extern
int SCIPvarGetFixDepth(
   VAR*             var                 /**< problem variable */
   );

/** gets fixing index of the variable in the depth level, where the binary variable was fixed, or -1 if unfixed or
 *  fixed during preprocessing
 */
extern
int SCIPvarGetFixIndex(
   VAR*             var                 /**< problem variable */
   );

/** returns TRUE iff first variable was fixed earlier than second variable */
extern
Bool SCIPvarWasFixedEarlier(
   VAR*             var1,               /**< first problem variable */
   VAR*             var2                /**< second problem variable */
   );

/** gets type of bound change of fixed binary variable (fixed due to branching or due to inference) */
extern
BOUNDCHGTYPE SCIPvarGetBoundchgType(
   VAR*             var                 /**< problem variable */
   );

/** gets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
extern
Real SCIPvarGetBranchFactor(
   VAR*             var                 /**< problem variable */
   );

/** gets the branch priority of the variable; variables with higher priority should always be prefered to variables
 *  with lower priority
 */
extern
int SCIPvarGetBranchPriority(
   VAR*             var                 /**< problem variable */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPvarGetName(var)             (var)->name
#define SCIPvarGetStatus(var)           (VARSTATUS)((var)->varstatus)
#define SCIPvarIsOriginal(var)          ((var)->varstatus == SCIP_VARSTATUS_ORIGINAL \
      || ((var)->varstatus == SCIP_VARSTATUS_NEGATED && (var)->negatedvar->varstatus == SCIP_VARSTATUS_ORIGINAL))
#define SCIPvarIsTransformed(var)       ((var)->varstatus != SCIP_VARSTATUS_ORIGINAL \
      && ((var)->varstatus != SCIP_VARSTATUS_NEGATED || (var)->negatedvar->varstatus != SCIP_VARSTATUS_ORIGINAL))
#define SCIPvarIsNegated(var)           ((var)->varstatus == SCIP_VARSTATUS_NEGATED)
#define SCIPvarGetType(var)             ((VARTYPE)((var)->vartype))
#define SCIPvarIsInitial(var)           (var)->initial
#define SCIPvarIsRemoveable(var)        (var)->removeable
#define SCIPvarIsActive(var)            ((var)->probindex >= 0)
#define SCIPvarGetIndex(var)            (var)->index
#define SCIPvarGetProbindex(var)        (var)->probindex
#define SCIPvarGetTransVar(var)         (var)->data.transvar
#define SCIPvarGetCol(var)              (var)->data.col
#define SCIPvarGetAggrVar(var)          (var)->data.aggregate.var
#define SCIPvarGetAggrScalar(var)       (var)->data.aggregate.scalar
#define SCIPvarGetAggrConstant(var)     (var)->data.aggregate.constant
#define SCIPvarGetMultaggrNVars(var)    (var)->data.multaggr.nvars
#define SCIPvarGetMultaggrVars(var)     (var)->data.multaggr.vars
#define SCIPvarGetMultaggrScalars(var)  (var)->data.multaggr.scalars
#define SCIPvarGetMultaggrConstant(var) (var)->data.multaggr.constant
#define SCIPvarGetNegatedVar(var)       (var)->negatedvar
#define SCIPvarGetNegationVar(var)      (var)->negatedvar
#define SCIPvarGetNegationConstant(var) (var)->data.negate.constant
#define SCIPvarGetObj(var)              (var)->obj
#define SCIPvarGetLbGlobal(var)         (var)->glbdom.lb
#define SCIPvarGetUbGlobal(var)         (var)->glbdom.ub
#define SCIPvarGetLbLocal(var)          (var)->locdom.lb
#define SCIPvarGetUbLocal(var)          (var)->locdom.ub
#define SCIPvarGetInferVar(var)         (var)->infervar
#define SCIPvarGetInferCons(var)        (var)->infercons
#define SCIPvarGetInferInfo(var)        (var)->inferinfo
#define SCIPvarGetFixDepth(var)         (var)->fixdepth
#define SCIPvarGetFixIndex(var)         (var)->fixindex
#define SCIPvarWasFixedEarlier(var1,var2) ((var1)->fixdepth >= 0                          \
                                            && ((var2)->fixdepth == -1                    \
                                              || (var1)->fixdepth < (var2)->fixdepth      \
                                              || ((var1)->fixdepth == (var2)->fixdepth    \
                                                && (var1)->fixindex < (var2)->fixindex)))
#define SCIPvarGetBoundchgType(var)     (var)->boundchgtype
#define SCIPvarGetBranchFactor(var)     (var)->branchfactor
#define SCIPvarGetBranchPriority(var)   (var)->branchpriority

#endif

/** gets best local bound of variable with respect to the objective function */
extern
Real SCIPvarGetBestBound(
   VAR*             var                 /**< problem variable */
   );

/** gets worst local bound of variable with respect to the objective function */
extern
Real SCIPvarGetWorstBound(
   VAR*             var                 /**< problem variable */
   );

/** gets type (lower or upper) of best bound of variable with respect to the objective function */
extern
BOUNDTYPE SCIPvarGetBestBoundType(
   VAR*             var                 /**< problem variable */
   );

/** gets type (lower or upper) of worst bound of variable with respect to the objective function */
extern
BOUNDTYPE SCIPvarGetWorstBoundType(
   VAR*             var                 /**< problem variable */
   );

/** gets primal LP solution value of variable */
extern
Real SCIPvarGetLPSol(
   VAR*             var                 /**< problem variable */
   );

/** gets pseudo solution value of variable at current node */
extern
Real SCIPvarGetPseudoSol(
   VAR*             var                 /**< problem variable */
   );

/** gets current LP or pseudo solution value of variable */
extern
Real SCIPvarGetSol(
   VAR*             var,                /**< problem variable */
   Bool             getlpval            /**< should the LP solution value be returned? */
   );


#endif
