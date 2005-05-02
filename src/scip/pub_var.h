/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_var.h,v 1.43 2005/05/02 11:42:55 bzfpfend Exp $"

/**@file   pub_var.h
 * @brief  public methods for problem variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_VAR_H__
#define __PUB_VAR_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_history.h"
#include "scip/type_var.h"
#include "scip/type_cons.h"

#ifdef NDEBUG
#include "scip/struct_var.h"
#include "scip/history.h"
#endif



/*
 * methods for variables 
 */

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
int SCIPvarCompare(
   VAR*             var1,               /**< first problem variable */
   VAR*             var2                /**< second problem variable */
   );

/** gets corresponding active, fixed, or multi-aggregated problem variable of a variable */
extern
VAR* SCIPvarGetProbvar(
   VAR*             var                 /**< problem variable */
   );

/** gets corresponding active, fixed, or multi-aggregated problem variable of a binary variable and updates the given
 *  negation status
 */
extern
RETCODE SCIPvarGetProbvarBinary(
   VAR**            var,                /**< pointer to binary problem variable */
   Bool*            negated             /**< pointer to update the negation status */
   );

/** transforms given variable, boundtype and bound to the corresponding active, fixed, or multi-aggregated variable
 *  values
 */
extern
RETCODE SCIPvarGetProbvarBound(
   VAR**            var,                /**< pointer to problem variable */
   Real*            bound,              /**< pointer to bound value to transform */
   BOUNDTYPE*       boundtype           /**< pointer to type of bound: lower or upper bound */
   );

/** transforms given variable, scalar and constant to the corresponding active, fixed, or multi-aggregated variable,
 *  scalar and constant;
 *  if the variable resolves to a fixed variable, "scalar" will be 0.0 and the value of the sum will be stored
 *  in "constant"
 */
extern
RETCODE SCIPvarGetProbvarSum(
   VAR**            var,                /**< pointer to problem variable x in sum a*x + c */
   Real*            scalar,             /**< pointer to scalar a in sum a*x + c */
   Real*            constant            /**< pointer to constant c in sum a*x + c */
   );

/** retransforms given variable, scalar and constant to the corresponding original variable, scalar and constant,
 *  if possible;
 *  if the retransformation is impossible, NULL is returned as variable
 */
extern
RETCODE SCIPvarGetOrigvarSum(
   VAR**            var,                /**< pointer to problem variable x in sum a*x + c */
   Real*            scalar,             /**< pointer to scalar a in sum a*x + c */
   Real*            constant            /**< pointer to constant c in sum a*x + c */
   );

/** returns the number of times, a bound of the variable was changed in given direction due to branching */
extern
Longint SCIPvarGetNBranchings(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the number of times, a bound of the variable was changed in given direction due to branching
 *  in the current run
 */
extern
Longint SCIPvarGetNBranchingsCurrentRun(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the number of inferences branching on this variable in given direction triggered */
extern
Longint SCIPvarGetNInferences(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the number of inferences branching on this variable in given direction triggered
 *  in the current run
 */
extern
Longint SCIPvarGetNInferencesCurrentRun(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the number of cutoffs branching on this variable in given direction produced */
extern
Longint SCIPvarGetNCutoffs(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the number of cutoffs branching on this variable in given direction produced in the current run */
extern
Longint SCIPvarGetNCutoffsCurrentRun(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average depth of bound changes in given direction due to branching on the variable */
extern
Real SCIPvarGetAvgBranchdepth(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average depth of bound changes in given direction due to branching on the variable
 *  in the current run
 */
extern
Real SCIPvarGetAvgBranchdepthCurrentRun(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
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

/** returns the user data of the variable */
extern
VARDATA* SCIPvarGetData(
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

/** returns whether variable is of integral type (binary, integer, or implicit integer) */
extern
Bool SCIPvarIsIntegral(
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

/** returns whether the variable is a COLUMN variable that is member of the current LP */
extern
Bool SCIPvarIsInLP(
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

/** gets original lower bound of original problem variable (i.e. the bound set in problem creation) */
extern
Real SCIPvarGetLbOriginal(
   VAR*             var                 /**< original problem variable */
   );

/** gets original upper bound of original problem variable (i.e. the bound set in problem creation) */
extern
Real SCIPvarGetUbOriginal(
   VAR*             var                 /**< original problem variable */
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

/** gets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
extern
Real SCIPvarGetBranchFactor(
   VAR*             var                 /**< problem variable */
   );

/** gets the branch priority of the variable; variables with higher priority should always be preferred to variables
 *  with lower priority
 */
extern
int SCIPvarGetBranchPriority(
   VAR*             var                 /**< problem variable */
   );

/** gets the preferred branch direction of the variable (downwards, upwards, or auto) */
extern
BRANCHDIR SCIPvarGetBranchDirection(
   VAR*             var                 /**< problem variable */
   );

/** gets number of variable lower bounds x >= b_i*z_i + d_i of given variable x */
extern
int SCIPvarGetNVlbs(
   VAR*             var                 /**< problem variable */
   );

/** gets array with bounding variables z_i in variable lower bounds x >= b_i*z_i + d_i of given variable x */
extern
VAR** SCIPvarGetVlbVars(
   VAR*             var                 /**< problem variable */
   );

/** gets array with bounding coefficients b_i in variable lower bounds x >= b_i*z_i + d_i of given variable x */
extern
Real* SCIPvarGetVlbCoefs(
   VAR*             var                 /**< problem variable */
   );

/** gets array with bounding constants d_i in variable lower bounds x >= b_i*z_i + d_i of given variable x */
extern
Real* SCIPvarGetVlbConstants(
   VAR*             var                 /**< problem variable */
   );

/** gets number of variable upper bounds x <= b_i*z_i + d_i of given variable x */
extern
int SCIPvarGetNVubs(
   VAR*             var                 /**< problem variable */
   );

/** gets array with bounding variables z_i in variable upper bounds x <= b_i*z_i + d_i of given variable x */
extern
VAR** SCIPvarGetVubVars(
   VAR*             var                 /**< problem variable */
   );

/** gets array with bounding coefficients b_i in variable upper bounds x <= b_i*z_i + d_i of given variable x */
extern
Real* SCIPvarGetVubCoefs(
   VAR*             var                 /**< problem variable */
   );

/** gets array with bounding constants d_i in variable upper bounds x <= b_i*z_i + d_i of given variable x */
extern
Real* SCIPvarGetVubConstants(
   VAR*             var                 /**< problem variable */
   );

/** gets number of implications  y <= b or y >= b for x == 0 or x == 1 of given variable x, 
 *  there are no implications for nonbinary variable x
 */
extern
int SCIPvarGetNImpls(
   VAR*             var,                /**< problem variable */
   Bool             varfixing           /**< FALSE for implications for x == 0, TRUE for x == 1 */
   );

/** gets number of implications  y <= 0 or y >= 1 for x == 0 or x == 1 of given variable x with binary y, 
 *  there are no implications for nonbinary variable x
 */
extern
int SCIPvarGetNBinImpls(
   VAR*             var,                /**< problem variable */
   Bool             varfixing           /**< FALSE for implications for x == 0, TRUE for x == 1 */
   );

/** gets array with implication variables y of implications  y <= b or y >= b for x == 0 or x == 1 of given variable x,  
 *  there are no implications for nonbinary variable x
 */
extern
VAR** SCIPvarGetImplVars(
   VAR*             var,                /**< problem variable */
   Bool             varfixing           /**< FALSE for implications for x == 0, TRUE for x == 1 */
   );

/** gets array with implication types of implications  y <= b or y >= b for x == 0 or x == 1 of given variable x
 *  (SCIP_BOUNDTYPE_UPPER if y <= b, SCIP_BOUNDTYPE_LOWER if y >= b), 
 *  there are no implications for nonbinary variable x
 */
extern
BOUNDTYPE* SCIPvarGetImplTypes(
   VAR*             var,                /**< problem variable */
   Bool             varfixing           /**< FALSE for implications for x == 0, TRUE for x == 1 */
   );

/** gets array with implication bounds b of implications  y <= b or y >= b for x == 0 or x == 1 of given variable x,  
 *  there are no implications for nonbinary variable x
 */
extern
Real* SCIPvarGetImplBounds(
   VAR*             var,                /**< problem variable */
   Bool             varfixing           /**< FALSE for implications for x == 0, TRUE for x == 1 */
   );

/** gets array with unique ids of implications  y <= b or y >= b for x == 0 or x == 1 of given variable x,  
 *  there are no implications for nonbinary variable x
 */
extern
int* SCIPvarGetImplIds(
   VAR*             var,                /**< problem variable */
   Bool             varfixing           /**< FALSE for implications for x == 0, TRUE for x == 1 */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPvarGetName(var)             (var)->name
#define SCIPvarGetData(var)             (var)->vardata
#define SCIPvarGetStatus(var)           (VARSTATUS)((var)->varstatus)
#define SCIPvarIsOriginal(var)          ((var)->varstatus == SCIP_VARSTATUS_ORIGINAL \
      || ((var)->varstatus == SCIP_VARSTATUS_NEGATED && (var)->negatedvar->varstatus == SCIP_VARSTATUS_ORIGINAL))
#define SCIPvarIsTransformed(var)       ((var)->varstatus != SCIP_VARSTATUS_ORIGINAL \
      && ((var)->varstatus != SCIP_VARSTATUS_NEGATED || (var)->negatedvar->varstatus != SCIP_VARSTATUS_ORIGINAL))
#define SCIPvarIsNegated(var)           ((var)->varstatus == SCIP_VARSTATUS_NEGATED)
#define SCIPvarGetType(var)             ((VARTYPE)((var)->vartype))
#define SCIPvarIsIntegral(var)          ((var)->vartype != SCIP_VARTYPE_CONTINUOUS)
#define SCIPvarIsInitial(var)           (var)->initial
#define SCIPvarIsRemoveable(var)        (var)->removeable
#define SCIPvarIsActive(var)            ((var)->probindex >= 0)
#define SCIPvarGetIndex(var)            (var)->index
#define SCIPvarGetProbindex(var)        (var)->probindex
#define SCIPvarGetTransVar(var)         (var)->data.original.transvar
#define SCIPvarGetCol(var)              (var)->data.col
#define SCIPvarIsInLP(var)              ((var)->varstatus == SCIP_VARSTATUS_COLUMN && SCIPcolIsInLP((var)->data.col))
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
#define SCIPvarGetLbOriginal(var)       ((var)->varstatus == SCIP_VARSTATUS_ORIGINAL \
      ? (var)->data.original.origdom.lb                                 \
      : (var)->data.negate.constant - (var)->negatedvar->data.original.origdom.ub)
#define SCIPvarGetUbOriginal(var)       ((var)->varstatus == SCIP_VARSTATUS_ORIGINAL \
      ? (var)->data.original.origdom.ub                                 \
      : (var)->data.negate.constant - (var)->negatedvar->data.original.origdom.lb)
#define SCIPvarGetLbGlobal(var)         (var)->glbdom.lb
#define SCIPvarGetUbGlobal(var)         (var)->glbdom.ub
#define SCIPvarGetLbLocal(var)          (var)->locdom.lb
#define SCIPvarGetUbLocal(var)          (var)->locdom.ub
#define SCIPvarGetBranchFactor(var)     (var)->branchfactor
#define SCIPvarGetBranchPriority(var)   (var)->branchpriority
#define SCIPvarGetBranchDirection(var)  (var)->branchdirection
#define SCIPvarGetNVlbs(var)            ((var)->vlbs != NULL ? (var)->vlbs->len : 0)
#define SCIPvarGetVlbVars(var)          ((var)->vlbs != NULL ? (var)->vlbs->vars : NULL)
#define SCIPvarGetVlbCoefs(var)         ((var)->vlbs != NULL ? (var)->vlbs->coefs : NULL)
#define SCIPvarGetVlbConstants(var)     ((var)->vlbs != NULL ? (var)->vlbs->constants : NULL)
#define SCIPvarGetNVubs(var)            ((var)->vubs != NULL ? (var)->vubs->len : 0)
#define SCIPvarGetVubVars(var)          ((var)->vubs != NULL ? (var)->vubs->vars : NULL)
#define SCIPvarGetVubCoefs(var)         ((var)->vubs != NULL ? (var)->vubs->coefs : NULL)
#define SCIPvarGetVubConstants(var)     ((var)->vubs != NULL ? (var)->vubs->constants : NULL)
#define SCIPvarGetNImpls(var, fix)      ((var)->implics != NULL ? (var)->implics->nimpls[fix] : 0)
#define SCIPvarGetNBinImpls(var, fix)   ((var)->implics != NULL ? (var)->implics->nbinimpls[fix] : 0)
#define SCIPvarGetImplVars(var, fix)    ((var)->implics != NULL ? (var)->implics->implvars[fix] : NULL)
#define SCIPvarGetImplTypes(var, fix)   ((var)->implics != NULL ? (var)->implics->impltypes[fix] : NULL)
#define SCIPvarGetImplBounds(var, fix)  ((var)->implics != NULL ? (var)->implics->implbounds[fix] : NULL)
#define SCIPvarGetImplIds(var, fix)     ((var)->implics != NULL ? (var)->implics->implids[fix] : NULL)
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

/** returns the solution of the variable in the root node's relaxation, if the root relaxation is not yet completely
 *  solved, zero is returned
 */
extern
Real SCIPvarGetRootSol(
   VAR*             var                 /**< problem variable */
   );

/** returns a weighted average solution value of the variable in all feasible primal solutions found so far */
extern
Real SCIPvarGetAvgSol(
   VAR*             var                 /**< problem variable */
   );

/** returns the bound change information for the last lower bound change on given active problem variable before or
 *  after the bound change with the given index was applied;
 *  returns NULL, if no change to the lower bound was applied up to this point of time
 */
extern
BDCHGINFO* SCIPvarGetLbchgInfo(
   VAR*             var,                /**< active problem variable */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   );

/** returns the bound change information for the last upper bound change on given active problem variable before or
 *  after the bound change with the given index was applied;
 *  returns NULL, if no change to the upper bound was applied up to this point of time
 */
extern
BDCHGINFO* SCIPvarGetUbchgInfo(
   VAR*             var,                /**< active problem variable */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   );

/** returns the bound change information for the last lower or upper bound change on given active problem variable
 *  before or after the bound change with the given index was applied;
 *  returns NULL, if no change to the lower/upper bound was applied up to this point of time
 */
extern
BDCHGINFO* SCIPvarGetBdchgInfo(
   VAR*             var,                /**< active problem variable */
   BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   );

/** returns lower bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
extern
Real SCIPvarGetLbAtIndex(
   VAR*             var,                /**< problem variable */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   );

/** returns upper bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
extern
Real SCIPvarGetUbAtIndex(
   VAR*             var,                /**< problem variable */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   );

/** returns lower or upper bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
extern
Real SCIPvarGetBdAtIndex(
   VAR*             var,                /**< problem variable */
   BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   );

/** returns whether the binary variable was fixed at the time given by the bound change index */
extern
Bool SCIPvarWasFixedAtIndex(
   VAR*             var,                /**< problem variable */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   );

/** returns the last bound change index, at which the bounds of the given variable were tightened */
extern
BDCHGIDX* SCIPvarGetLastBdchgIndex(
   VAR*             var                 /**< problem variable */
   );

/** returns the last depth level, at which the bounds of the given variable were tightened;
 *  returns -2, if the variable's bounds are still the global bounds
 *  returns -1, if the variable was fixed in presolving
 */
extern
int SCIPvarGetLastBdchgDepth(
   VAR*             var                 /**< problem variable */
   );

/** returns whether the first binary variable was fixed earlier than the second one;
 *  returns FALSE, if the first variable is not fixed, and returns TRUE, if the first variable is fixed, but the
 *  second one is not fixed
 */
extern
Bool SCIPvarWasFixedEarlier(
   VAR*             var1,               /**< first binary variable */
   VAR*             var2                /**< second binary variable */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns whether first bound change index belongs to an earlier applied bound change than second one;
 *  if a bound change index is NULL, the bound change index represents the current time, i.e. the time after the
 *  last bound change was applied to the current node
 */
extern
Bool SCIPbdchgidxIsEarlier(
   BDCHGIDX*        bdchgidx1,          /**< first bound change index, or NULL */
   BDCHGIDX*        bdchgidx2           /**< second bound change index, or NULL */
   );

/** returns whether first bound change index belongs to an earlier applied bound change than second one */
extern
Bool SCIPbdchgidxIsEarlierNonNull(
   BDCHGIDX*        bdchgidx1,          /**< first bound change index */
   BDCHGIDX*        bdchgidx2           /**< second bound change index */
   );

/** returns old bound that was overwritten for given bound change information */
extern
Real SCIPbdchginfoGetOldbound(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returns new bound installed for given bound change information */
extern
Real SCIPbdchginfoGetNewbound(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returns variable that belongs to the given bound change information */
extern
VAR* SCIPbdchginfoGetVar(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returns whether the bound change information belongs to a branching decision or a deduction */
extern
BOUNDCHGTYPE SCIPbdchginfoGetChgtype(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returns whether the bound change information belongs to a lower or upper bound change */
extern
BOUNDTYPE SCIPbdchginfoGetBoundtype(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returs depth level of given bound change information */
extern
int SCIPbdchginfoGetDepth(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returs bound change position in its depth level of given bound change information */
extern
int SCIPbdchginfoGetPos(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returs bound change index of given bound change information */
extern
BDCHGIDX* SCIPbdchginfoGetIdx(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returs inference variable of given bound change information */
extern
VAR* SCIPbdchginfoGetInferVar(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returs inference constraint of given bound change information */
extern
CONS* SCIPbdchginfoGetInferCons(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returs inference propagator of given bound change information, or NULL if no propagator was responsible */
extern
PROP* SCIPbdchginfoGetInferProp(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returs inference user information of given bound change information */
extern
int SCIPbdchginfoGetInferInfo(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returs inference bound of inference variable of given bound change information */
extern
BOUNDTYPE SCIPbdchginfoGetInferBoundtype(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

/** returns whether the bound change has an inference reason (constraint or propagator), that can be resolved */
extern
Bool SCIPbdchginfoHasInferenceReason(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPbdchgidxIsEarlierNonNull(idx1,idx2)                         \
   ((idx1)->depth < (idx2)->depth || ((idx1)->depth == (idx2)->depth && (idx1)->pos < (idx2)->pos))
#define SCIPbdchgidxIsEarlier(idx1,idx2)                                \
   ((idx1) != NULL && ((idx2) == NULL || SCIPbdchgidxIsEarlierNonNull(idx1, idx2)))
#define SCIPbdchginfoGetOldbound(bdchginfo)       (bdchginfo)->oldbound
#define SCIPbdchginfoGetNewbound(bdchginfo)       (bdchginfo)->newbound
#define SCIPbdchginfoGetVar(bdchginfo)            (bdchginfo)->var
#define SCIPbdchginfoGetChgtype(bdchginfo)        (BOUNDCHGTYPE)((bdchginfo)->boundchgtype)
#define SCIPbdchginfoGetBoundtype(bdchginfo)      (BOUNDTYPE)((bdchginfo)->boundtype)
#define SCIPbdchginfoGetDepth(bdchginfo)          (bdchginfo)->bdchgidx.depth
#define SCIPbdchginfoGetPos(bdchginfo)            (bdchginfo)->bdchgidx.pos
#define SCIPbdchginfoGetIdx(bdchginfo)            (&(bdchginfo)->bdchgidx)
#define SCIPbdchginfoGetInferVar(bdchginfo)       (bdchginfo)->inferencedata.var
#define SCIPbdchginfoGetInferCons(bdchginfo)      (bdchginfo)->inferencedata.reason.cons
#define SCIPbdchginfoGetInferProp(bdchginfo)      (bdchginfo)->inferencedata.reason.prop
#define SCIPbdchginfoGetInferInfo(bdchginfo)      (bdchginfo)->inferencedata.info
#define SCIPbdchginfoGetInferBoundtype(bdchginfo) (BOUNDTYPE)((bdchginfo)->inferboundtype)
#define SCIPbdchginfoHasInferenceReason(bdchginfo)                      \
   (((bdchginfo)->boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER)          \
      || ((bdchginfo)->boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER && (bdchginfo)->inferencedata.reason.prop != NULL))

#endif

#endif
