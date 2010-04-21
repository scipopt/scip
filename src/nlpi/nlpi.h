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
#pragma ident "@(#) $Id: nlpi.h,v 1.2 2010/04/21 14:21:14 bzfviger Exp $"

/**@file   nlpi.h
 * @brief  internal methods for NLPI solver interfaces
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_H__
#define __SCIP_NLPI_H__

#include "nlpi/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates an NLP solver interface */
extern
SCIP_RETCODE SCIPnlpiCreate(
   SCIP_NLPI**                     nlpi,                        /**< pointer to NLP interface data structure */
   const char*                     name,                        /**< name of NLP interface */
   const char*                     description,                 /**< description of NLP interface */
   int                             priority,                    /**< priority of NLP interface */
   SCIP_DECL_NLPIFREE              ((*nlpifree)),               /**< free NLPI user data */
   SCIP_DECL_NLPIGETSOLVERPOINTER  ((*nlpigetsolverpointer)),   /**< get solver pointer */
   SCIP_DECL_NLPICREATEPROBLEM     ((*nlpicreateproblem)),      /**< create a new problem instance */
   SCIP_DECL_NLPIFREEPROBLEM       ((*nlpifreeproblem)),        /**< free a problem instance */
   SCIP_DECL_NLPIGETPROBLEMPOINTER ((*nlpigetproblempointer)),  /**< get problem pointer */
   SCIP_DECL_NLPIADDVARS           ((*nlpiaddvars)),            /**< add variables */
   SCIP_DECL_NLPIADDCONSTRAINTS    ((*nlpiaddconstraints)),     /**< add constraints */
   SCIP_DECL_NLPISETOBJECTIVE      ((*nlpisetobjective)),       /**< set objective */
   SCIP_DECL_NLPICHGVARBOUNDS      ((*nlpichgvarbounds)),       /**< change variable bounds */
   SCIP_DECL_NLPICHGCONSBOUNDS     ((*nlpichgconsbounds)),      /**< change constraint bounds */
   SCIP_DECL_NLPIDELVARSET         ((*nlpidelvarset)),          /**< delete a set of constraints */
   SCIP_DECL_NLPIDELCONSSET        ((*nlpidelconsset)),         /**< delete a set of constraints */
   SCIP_DECL_NLPICHGLINEARCOEFS    ((*nlpichglinearcoef)),      /**< change one coefficient  in linear part */
   SCIP_DECL_NLPICHGQUADCOEFS      ((*nlpichgquadcoef)),        /**< change one coefficient  in quadratic part */
   SCIP_DECL_NLPICHGNONLINCOEF     ((*nlpichgnonlincoef)),      /**< change one parameter in nonlinear expressions */
   SCIP_DECL_NLPISETINITIALGUESS   ((*nlpisetinitialguess)),    /**< set initial guess for primal variables */
   SCIP_DECL_NLPISOLVE             ((*nlpisolve)),              /**< solve NLP */
   SCIP_DECL_NLPIGETSOLSTAT        ((*nlpigetsolstat)),         /**< get solution status */
   SCIP_DECL_NLPIGETTERMSTAT       ((*nlpigettermstat)),        /**< get termination status */
   SCIP_DECL_NLPIGETSOLUTION       ((*nlpigetsolution)),        /**< get solution */
   SCIP_DECL_NLPIGETSTATISTICS     ((*nlpigetstatistics)),      /**< get solve statistics */
   SCIP_DECL_NLPIGETWARMSTARTSIZE  ((*nlpigetwarmstartsize)),   /**< get size for warmstart object buffer */
   SCIP_DECL_NLPIGETWARMSTARTMEMO  ((*nlpigetwarmstartmemo)),   /**< get warmstart object */
   SCIP_DECL_NLPISETWARMSTARTMEMO  ((*nlpisetwarmstartmemo)),   /**< set warmstart object */
   SCIP_DECL_NLPIGETINTPAR         ((*nlpigetintpar)),          /**< get value of integer parameter */
   SCIP_DECL_NLPISETINTPAR         ((*nlpisetintpar)),          /**< set value of integer parameter */
   SCIP_DECL_NLPIGETREALPAR        ((*nlpigetrealpar)),         /**< get value of floating point parameter */
   SCIP_DECL_NLPISETREALPAR        ((*nlpisetrealpar)),         /**< set value of floating point parameter */
   SCIP_DECL_NLPIGETSTRINGPAR      ((*nlpigetstringpar)),       /**< get value of string parameter */
   SCIP_DECL_NLPISETSTRINGPAR      ((*nlpisetstringpar)),       /**< set value of string parameter */
   SCIP_NLPIDATA*                  nlpidata                     /**< NLP interface local data */
);

/** frees NLPI user data */
extern 
SCIP_RETCODE SCIPnlpiFree(
   SCIP_NLPI**           nlpi                /**< pointer to NLPI data structure */
);

/** gets pointer for NLP solver
 * @return void pointer to solver
 */
extern
void* SCIPnlpiGetSolverPointer(
   SCIP_NLPI*            nlpi                /**< pointer to NLPI datastructure */
);

/** creates a problem instance */
extern
SCIP_RETCODE SCIPnlpiCreateProblem(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM**    problem,            /**< pointer to store problem data */
   const char*           name                /**< name of problem, can be NULL */
);

/** frees a problem instance */
extern
SCIP_RETCODE SCIPnlpiFreeProblem(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM**    problem             /**< pointer where problem data is stored */
);

/** gets pointer to solver-internal problem instance
 * @return void pointer to problem instance
 */
extern
void* SCIPnlpiGetProblemPointer(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer where problem data is stored */
);

/** add variables to nlpi */
extern
SCIP_RETCODE SCIPnlpiAddVars(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI data structure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int                   nvars,              /**< number of variables */               
   const SCIP_Real*      lbs,                /**< lower bounds of variables, can be NULL if -infinity */
   const SCIP_Real*      ubs,                /**< upper bounds of variables, can be NULL if +infinity */
   const char**          varnames            /**< varnames names of variables, can be NULL */
);

/** add constraints to nlpi */
extern
SCIP_RETCODE SCIPnlpiAddConstraints(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI data structure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int                   nconss,             /**< number of added constraints */
   const SCIP_Real*      lhss,               /**< left hand sides of constraints */
   const SCIP_Real*      rhss,               /**< right hand sides of constraints */
   const int*            nlininds,           /**< number of linear coefficients for each constraint, may be NULL in case of no linear part */
   int* const*           lininds,            /**< indices of variables for linear coefficients for each constraint, may be NULL in case of no linear part */
   SCIP_Real* const*     linvals,            /**< values of linear coefficient for each constraint, may be NULL in case of no linear part */
   const int*            nquadrows,          /**< number of columns in matrix of quadratic part for each constraint, may be
                                              * NULL in case of no quadratic part in any constraint */
   int* const*           quadrowidxs,        /**< indices of variables for which a quadratic part is specified, may be NULL
                                              * in case of no quadratic part in any constraint */
   int* const*           quadoffsets,        /**< start index of each rows quadratic coefficients in quadinds[.] and quadvals[.],
                                              * indices are given w.r.t. quadrowidxs., i.e., quadoffsets[.][i] gives the start
                                              * index of row quadrowidxs[.][i] in quadvals[.], quadoffsets[.][nquadrows[.]] gives
                                              * length of quadinds[.] and quadvals[.], entry of array may be NULL in case of no
                                              * quadratic part, may be NULL in case of no quadratic part in any constraint */
   int* const*           quadinds,           /**< column indices w.r.t. quadrowidxs, i.e., quadrowidxs[quadinds[.][i]] gives the
                                              * index of the variable corresponding to entry i, entry of array may be NULL in
                                              * case of no quadratic part, may be NULL in case of no quadratic part in any constraint */
   SCIP_Real* const*     quadvals,           /**< coefficient values, entry of array may be NULL in case of no quadratic part,
                                              * may be NULL in case of no quadratic part in any constraint */
   int* const*           exprvaridxs,        /**< indices of variables in expression tree, maps variable indices in expression
                                              * tree to indices in nlp, entry of array may be NULL in case of no expression
                                              * tree, may be NULL in case of no expression tree in any constraint */
   SCIP_EXPRTREE* const* exprtrees,          /**< exprtrees expression tree for nonquadratic part of constraints, entry of
                                              * array may be NULL in case of no nonquadratic part, may be NULL in case of no
                                              * nonquadratic part in any constraint */
   const char**          names               /**< names of constraints, may be NULL or entries may be NULL */
);


/** sets or overwrites objective, a minimization problem is expected */
extern
SCIP_RETCODE SCIPnlpiSetObjective(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int                   nlins,              /**< number of linear variables */
   const int*            lininds,            /**< variable indices, may be NULL in case of no linear part */
   const SCIP_Real*      linvals,            /**< coefficient values, may be NULL in case of no linear part */
   int                   nquadcols,          /**< number of columns in matrix of quadratic part */
   const int*            quadcols,           /**< indices of variables for which a quadratic part is specified, may be NULL in
                                              * case of no quadratic part */
   const int*            quadoffsets,        /**< start index of each rows quadratic coefficients in quadinds and quadvals,
                                              * quadoffsets[.][nquadcols] gives length of quadinds and quadvals, may be NULL in
                                              * case of no quadratic part */
   const int*            quadinds,           /**< column indices, may be NULL in case of no quadratic part */
   const SCIP_Real*      quadvals,           /**< coefficient values, may be NULL in case of no quadratic part */
   const int*            exprvaridxs,        /**< indices of variables in expression tree, maps variable indices in expression
                                              * tree to indices in nlp, may be NULL in case of no expression tree */
   const SCIP_EXPRTREE*  exprtree,           /**< expression tree for nonquadratic part of objective function, may be NULL in
                                              * case of no nonquadratic part */
   const SCIP_Real       constant            /**< objective value offset*/
);

/** change variable bounds */
extern
SCIP_RETCODE SCIPnlpiChgVarBounds(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   const int             nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds */
   const SCIP_Real*      ubs                 /**< new upper bounds */
);

/** change constraint bounds */
extern
SCIP_RETCODE SCIPnlpiChgConsBounds(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   const int             nconss,             /**< number of constraints to change bounds */
   const int*            indices,            /**< indices of constraints to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds */
   const SCIP_Real*      ubs                 /**< new upper bounds */
);

/** delete a set of variables */
extern
SCIP_RETCODE SCIPnlpiDelVarSet(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int*                  dstats              /**< deletion status of vars; 1 if var should be deleted, 0 if not; afterwards -1
                                              * if var was deleted */
);

/** delete a set of constraints */
extern
SCIP_RETCODE SCIPnlpiDelConsSet(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int*                  dstats              /**< deletion status of rows; 1 if row should be deleted, 0 if not; afterwards -1
                                              * if row was deleted */
);

/** changes or adds one linear coefficient in a constraint or objective */
extern
SCIP_RETCODE SCIPnlpiChgLinearCoefs(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   const int             idx,                /**< index of constraint or -1 for objective */
   int                   nvals,              /**< number of values in linear constraint */
   const int*            varidxs,            /**< indices of variable */
   const SCIP_Real*      vals                /**< new values for coefficient */
);

/** changes or adds one coefficient in the quadratic part of a constraint or objective */
extern
SCIP_RETCODE SCIPnlpiChgQuadCoefs(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   const int             idx,                /**< index of constraint or -1 for objective */
   const int             nentries,           /**< nentries number of values in quadratic constraint */
   const int*            rows,               /**< row offset containing modified indices */
   const int*            cols,               /**< cols containing modified indices to the corresponding row offset */
   SCIP_Real*            values              /**< coefficients corresponding to same indices as used when constraint/objective
                                              * was constructed */
);

/** change the value of one parameter in the nonlinear part */
extern
SCIP_RETCODE SCIPnlpiChgNonlinCoef(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   const int             idxcons,            /**< index of constraint or -1 for objective */
   const int             idxparam,           /**< index of parameter */
   SCIP_Real             value               /**< new value for nonlinear parameter */
);

/** sets initial guess for primal variables */
extern
SCIP_RETCODE SCIPnlpiSetInitialGuess(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_Real*            values              /**< initial starting solution, or NULL to clear previous starting solution */
);

/** tries to solve NLP */
extern
SCIP_RETCODE SCIPnlpiSolve(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
);

/** gives solution status
 * @return solution status */
extern
SCIP_NLPSOLSTAT SCIPnlpiGetSolstat(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
);

/** gives termination reason
 * @return termination status */
extern
SCIP_NLPTERMSTAT SCIPnlpiGetTermstat(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
);

/** gives primal solution */
extern
SCIP_RETCODE SCIPnlpiGetSolution(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_Real**           primalvalues        /**< pointer to store primal values */
);

/** gives solve statistics */
extern
SCIP_RETCODE SCIPnlpiGetStatistics(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
);

/** gives required size of a buffer to store a warmstart object */
extern
SCIP_RETCODE SCIPnlpiGetWarmstartSize(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   size_t*               size                /**< pointer to store required size for warmstart buffer */
);

/** stores warmstart information in buffer */
extern
SCIP_RETCODE SCIPnlpiGetWarmstartMemo(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   void*                 buffer              /**< memory to store warmstart information */
);

/** sets warmstart information in solver */
extern
SCIP_RETCODE SCIPnlpiSetWarmstartMemo(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   void*                 buffer              /**< warmstart information */
);

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP */
extern
SCIP_RETCODE SCIPnlpiGetIntPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
);

/** sets integer parameter of NLP */
extern
SCIP_RETCODE SCIPnlpiSetIntPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
);

/** gets floating point parameter of NLP
 * if problem is NULL and type == SCIP_NLPPAR_INFINITY, then gets solver-wide value for infinity */
extern
SCIP_RETCODE SCIPnlpiGetRealPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure, can be NULL only if type == SCIP_NLPPAR_INFINITY */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
);

/** sets floating point parameter of NLP
 * if problem is NULL and type == SCIP_NLPPAR_INFINITY, then sets solver-wide value for infinity */
extern
SCIP_RETCODE SCIPnlpiSetRealPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure, can be NULL only if type == SCIP_NLPPAR_INFINITY */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
);

/** gets string parameter of NLP */
extern
SCIP_RETCODE SCIPnlpiGetStringPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char**          sval                /**< pointer to store the parameter value, the user must not modify the string */
);

/** sets string parameter of NLP */
extern
SCIP_RETCODE SCIPnlpiSetStringPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char*           sval                /**< parameter value */
);

/** gets data of an NLPI */
SCIP_NLPIDATA* SCIPnlpiGetData(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
);

/** gets NLP solver name */
const char* SCIPnlpiGetName(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
);

/** creates an NLP statistics structure */
SCIP_RETCODE SCIPnlpStatisticsCreate(
   SCIP_NLPSTATISTICS**  statistics          /**< pointer where to store NLP statistics structure */
);

/** frees an NLP statistics structure */
void SCIPnlpStatisticsFree(
   SCIP_NLPSTATISTICS**  statistics          /**< pointer where to store NLP statistics structure */
);

/** gets the number of iterations from an NLP statistics structure */
int SCIPnlpStatisticsGetNIterations(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
);

/** gets the total time from an NLP statistics structure */
SCIP_Real SCIPnlpStatisticsGetTotalTime(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
);

/** sets the number of iterations in an NLP statistics structure */
void SCIPnlpStatisticsSetNIterations(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   int                   niterations         /**< number of iterations to store */
);

/** sets the total time in an NLP statistics structure */
void SCIPnlpStatisticsSetTotalTime(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   SCIP_Real             totaltime           /**< solution time to store */
);

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_H__ */
