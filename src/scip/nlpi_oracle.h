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

/**@file   nlpi_oracle.h
 * @brief  methods to store an NLP and request function, gradient, and hessian values
 * @author Stefan Vigerske
 * 
 * Not a full NLPI, but implements a part of an NLPI that takes care of the problem storage.
 * Is used to complete NLPIs for solvers that need only function, gradient, and hessian evaluation methods.
 */

#ifndef __SCIP_NLPI_ORACLE_H__
#define __SCIP_NLPI_ORACLE_H__

#include "scip/scip.h"
#include "scip/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_NlpiOracle SCIP_NLPIORACLE; /**< NLPI oracle data structure */

/** creates an NLPIORACLE data structure */
extern
SCIP_RETCODE SCIPnlpiOracleCreate(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE**     oracle              /**< pointer to store NLPIORACLE data structure */
   );

/** initializes NLPI oracle */
extern
SCIP_RETCODE SCIPnlpiOracleInit(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle              /**< NLPIORACLE data structure */
   );

/** frees an NLPIORACLE data structure */
extern 
SCIP_RETCODE SCIPnlpiOracleFree(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE**     oracle              /**< pointer to NLPIORACLE data structure */
   );

/** adds variables */
extern
SCIP_RETCODE SCIPnlpiOracleAddVars(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nvars,              /**< number of variables to add */
   const SCIP_Real*      lbs,                /**< array with lower bounds of new variables, or NULL if all -infinity */
   const SCIP_Real*      ubs,                /**< array with upper bounds of new variables, or NULL if all +infinity */
   const char**          varnames            /**< array with names of new variables, or NULL if no names should be stored */
   );

/** adds constraints 
 * 
 *  linear coefficients: row(=constraint) oriented matrix;
 *  quadratic coefficiens: row oriented matrix for each constraint
 */
extern
SCIP_RETCODE SCIPnlpiOracleAddConstraints(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nconss,             /**< number of constraints to add */
   const SCIP_Real*      lhss,               /**< array with left-hand sides of constraints, or NULL if all -infinity */
   const SCIP_Real*      rhss,               /**< array with right-hand sides of constraints, or NULL if all +infinity */
   const int*            nlininds,           /**< number of linear coefficients for each constraint, may be NULL in case of no linear part */
   int* const*           lininds,            /**< indices of variables for linear coefficients for each constraint, may be NULL in case of no linear part */
   SCIP_Real* const*     linvals,            /**< values of linear coefficient for each constraint, may be NULL in case of no linear part */
   const int*            nquadrows,          /**< NULL if no quadratic parts, otherwise nquadrows[.] gives the number of columns in the matrix 
                                              *   of the quadratic part, or 0 if no quadratic part in this constraint */
   int* const*           quadrowidxs,        /**< NULL if no quadratic parts, otherwise quadrowidx[.] gives the indices of variables 
                                              *   for which a quadratic part is specified, or NULL if no quadratic part in this constraint */
   int* const*           quadoffsets,        /**< NULL if no quadratic parts, otherwise quadoffsets[.] gives start index of each rows quadratic coefficients 
                                              *   in quadidxs[.] and quadvals[.] (quadoffsets[.][nvars] gives length of quadidxs[.] and quadvals[.]), 
                                              *   or NULL if no quadratic part in this constraint */
   int* const*           quadidxs,           /**< NULL if no quadratic parts, otherwise quadidxs[.] gives column indices for quadratic part, 
                                              *   or NULL if no quadratic part in this constraint */
   SCIP_Real* const*     quadvals,           /**< NULL if no quadratic parts, otherwise quadvals[.] gives coefficients in quadratic part, 
                                              *   or NULL if no quadratic part in this constraint */
   int* const*           exprvaridxs,        /**< NULL if no nonquadratic parts, otherwise epxrvaridxs[.] maps variable indices in expression tree to indices in nlp */
   SCIP_EXPRTREE* const* exprtrees,          /**< NULL if no nonquadratic parts, otherwise exprtrees[.] gives nonquadratic part, 
                                              *   or NULL if no nonquadratic part in this constraint */
   const char**          consnames           /**< names of new constraints, or NULL if no names should be stored */
   );

/** sets or overwrites objective, a minization problem is expected
 * 
 *  May change sparsity pattern.
 */
extern
SCIP_RETCODE SCIPnlpiOracleSetObjective(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real       constant,           /**< constant part of objective */
   int                   nlin,               /**< number of linear variable coefficients */ 
   const int*            lininds,            /**< indices of linear variables, or NULL if no linear part */
   const SCIP_Real*      linvals,            /**< coefficients of linear variables, or NULL if no linear part */
   int                   nquadrows,          /**< number of columns in the matrix of the quadratic part */
   const int*            quadrowidxs,        /**< indices of variables appearing in quadratic part, or NULL if no quadratic part */
   const int*            quadoffsets,        /**< start index of each rows quadratic coefficients in quadidxs and quadvals 
                                              *   (quadoffsets[.][nvars] gives length of quadidxs and quadvals), or NULL if no quadratic part */
   const int*            quadidxs,           /**< column indices in quadratic part, or NULL if no quadratic part */ 
   const SCIP_Real*      quadvals,           /**< coefficients in quadratic part, or NULL if no quadratic part */
   const int*            exprvaridxs,        /**< maps variable indices in expression tree to indices in nlp, or NULL if no nonquadratic part */
   const SCIP_EXPRTREE*  exprtree            /**< expression tree of nonquadratic part, or NULL if no nonquadratic part */
   );

/** change variable bounds */
extern
SCIP_RETCODE SCIPnlpiOracleChgVarBounds(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< array with indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< array with new lower bounds, or NULL if all should be -infty */
   const SCIP_Real*      ubs                 /**< array with new upper bounds, or NULL if all should be +infty */
   );

/** change constraint bounds */
extern
SCIP_RETCODE SCIPnlpiOracleChgConsBounds(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             nconss,             /**< number of constraints to change bounds */
   const int*            indices,            /**< array with indices of constraints to change bounds */
   const SCIP_Real*      lhss,               /**< array with new left-hand sides, or NULL if all should be -infty */
   const SCIP_Real*      rhss                /**< array with new right-hand sides, or NULL if all should be +infty */
   );

/** deletes a set of variables */
extern
SCIP_RETCODE SCIPnlpiOracleDelVarSet(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< array with deletion status of vars in input (1 if var should be deleted, 0 if not); 
                                              *   new position of var in output (-1 if var was deleted) */
   );

/** deletes a set of constraints */
extern
SCIP_RETCODE SCIPnlpiOracleDelConsSet(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< array with deletion status of rows in input (1 if row should be deleted, 0 if not); 
                                              *   new position of row in output (-1 if row was deleted) */
   );

/** changes (or adds) linear coefficients in one constraint or objective
 */
extern
SCIP_RETCODE SCIPnlpiOracleChgLinearCoefs(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where linear coefficients should be changed, or -1 for objective */
   int                   nentries,           /**< number of coefficients to change */
   const int*            varidxs,            /**< array with indices of variables which coefficients should be changed */
   const SCIP_Real*      newcoefs            /**< array with new coefficients of variables */
   );

/** changes (or adds) coefficients in the quadratic part of one constraint or objective
 */
extern
SCIP_RETCODE SCIPnlpiOracleChgQuadCoefs(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where quadratic coefficients should be changed, or -1 for objective */
   const int             nentries,           /**< number of coefficients to change */
   const int*            rowidxs,            /**< array with row indices of quadratic matrix entries for which new values are provided */
   const int*            colidxs,            /**< array with column indices of quadratic matrix entries for which new values are provided */
   SCIP_Real*            newcoefs            /**< array with new quadratic coefficients */ 
   );

/** gives the current number of variables */
extern
int SCIPnlpiOracleGetNVars(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the current number of constraints */
extern
int SCIPnlpiOracleGetNConstraints(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the variables lower bounds */
extern
const SCIP_Real* SCIPnlpiOracleGetVarLbs(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the variables upper bounds */
extern
const SCIP_Real* SCIPnlpiOracleGetVarUbs(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** Gives maximum degree of a variable w.r.t. objective and all constraints.
 *  The degree of a variable is the degree of the summand where it appears in, and is infinity for nonpolynomial terms.
 */ 
extern
int SCIPnlpiOracleGetVarDegree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   varidx              /**< the variable for which the degree is returned */
   );

/** Gives maximum degree of all variables w.r.t. objective and all constraints.
 *  The degree of a variable is the degree of the summand where it appears in, and is infinity for nonpolynomial terms.
 */ 
extern
int* SCIPnlpiOracleGetVarDegrees(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the constraints left-hand sides */
extern
const SCIP_Real* SCIPnlpiOracleGetConstraintLhss(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the constraints right-hand sides */
extern
const SCIP_Real* SCIPnlpiOracleGetConstraintRhss(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** Gives maximum degree of a constraints.
 *  The degree of a constraint is the maximal degree of all summands which appear in it, and is infinity for nonpolynomial terms.
 */ 
extern
int SCIPnlpiOracleGetConstraintDegree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< index of constraint for which the degree is requested */
   );

/** evaluates the objective function in a given point */
extern
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveValue(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            objval              /**< pointer to store objective value */  
   );

/** evaluates one constraint function in a given point */
extern
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValue(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             considx,            /**< index of constraint to evaluate */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            conval              /**< pointer to store constraint value */  
   );

/** evaluates all constraint functions in a given point */
extern
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValues(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            convals             /**< pointer to store constraint values */  
   );

/** computes the objective gradient in a given point */
extern
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveGradient(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             newx,               /**< indicates whether the function has not been evaluated for this point before */
   SCIP_Real*            objval,             /**< pointer to buffer objective value */
   SCIP_Real*            objgrad             /**< pointer to buffer (dense) objective gradient */  
   );

/** computes a constraints gradient in a given point */
extern
SCIP_RETCODE SCIPnlpiOracleEvalConstraintGradient(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             considx,            /**< index of constraint to compute gradient for */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             newx,               /**< indicates whether the function has not been evaluated for this point before */ 
   SCIP_Real*            conval,             /**< pointer to store constraint value */
   SCIP_Real*            congrad             /**< pointer to store (dense) constraint gradient */  
   );

/** Gets sparsity pattern (rowwise) of Jacobian matrix.
 * 
 *  Note that internal data is returned in *offset and *col, thus the user does not need to allocate memory there.
 *  Adding or deleting constraints destroys the sparsity structure and make another call to this function necessary. 
 */
extern
SCIP_RETCODE SCIPnlpiOracleGetJacobianSparsity(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row, 
                                              *   offsets[nconss] gives length of col, can be NULL */
   );

/** Evaluates the Jacobi matrix in a given point.
 * 
 *  The values in the Jacobi matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetJacobianSparsity.
 *  The user need to call SCIPnlpiOracleGetJacobianSparsity at least ones before using this function. 
 */
extern
SCIP_RETCODE SCIPnlpiOracleEvalJacobian(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             newx,               /**< indicates whether some function has not been evaluated for this point before */
   SCIP_Real*            convals,            /**< pointer to store constraint values, can be NULL */ 
   SCIP_Real*            jacobi              /**< pointer to store sparse jacobian values */  
   );

/** Gets sparsity pattern of the Hessian matrix of the Lagrangian.
 * 
 *  Note that internal data is returned in *offset and *col, thus the user does not need to allocate memory there.
 *  Adding or deleting variables, objective, or constraints may destroy the sparsity structure and make another call to this function necessary.
 *  Only elements of the lower left triangle and the diagonal are counted.
 */
extern
SCIP_RETCODE SCIPnlpiOracleGetHessianLagSparsity(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row, 
                                              *   offsets[nconss] gives length of col, can be NULL */
   );

/** Evaluates the Hessian matrix of the Lagrangian in a given point.
 * 
 *  The values in the Hessian matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetHessianLagSparsity.
 *  The user must call SCIPnlpiOracleGetHessianLagSparsity at least ones before using this function. 
 *  Only elements of the lower left triangle and the diagonal are computed.
 */
extern
SCIP_RETCODE SCIPnlpiOracleEvalHessianLag(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             newx,               /**< indicates whether some function has not been evaluated for this point before */
   SCIP_Real             objfactor,          /**< weight for objective function */
   const SCIP_Real*      lambdas,            /**< array with weights (Lagrangian multipliers) for the constraints */ 
   SCIP_Real*            hessian             /**< pointer to store sparse hessian values */  
   );

/** prints the problem to a file. */
extern
SCIP_RETCODE SCIPnlpiOraclePrintProblem(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   );

/** prints the problem to a file in GAMS format */
extern
SCIP_RETCODE SCIPnlpiOraclePrintProblemGams(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real*            initval,            /**< starting point values for variables or NULL */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   );

#ifdef __cplusplus
}
#endif

#endif
