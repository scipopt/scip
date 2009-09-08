/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nlpi.h,v 1.4 2009/09/08 17:00:08 bzfviger Exp $"

/**@file   nlpi.h
 * @brief  internal methods for NLPI solver interfaces
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_H__
#define __SCIP_NLPI_H__

#include "scip/scip.h"
#include "scip/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates an NLP solver interface */
extern
SCIP_RETCODE SCIPnlpiCreate(
   SCIP*                           scip,                        /**< pointer to SCIP */
   SCIP_NLPI**                     nlpi,                        /**< pointer to NLP interface data structure */
   const char*                     name,                        /**< name of NLP interface */
   const char*                     description,                 /**< description of NLP interface */
   int                             priority,                    /**< priority of NLP interface */
   SCIP_DECL_NLPIINIT              ((*nlpiinit)),               /**< initialize NLPI user data */
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
   SCIP_DECL_NLPIGETSOLVERPOINTER  ((*nlpigetsolverpointer)),   /**< get solver pointer */
   SCIP_DECL_NLPIGETINTPAR         ((*nlpigetintpar)),          /**< get value of integer parameter */
   SCIP_DECL_NLPISETINTPAR         ((*nlpisetintpar)),          /**< set value of integer parameter */
   SCIP_DECL_NLPIGETREALPAR        ((*nlpigetrealpar)),         /**< get value of floating point parameter */
   SCIP_DECL_NLPISETREALPAR        ((*nlpisetrealpar)),         /**< set value of floating point parameter */
   SCIP_DECL_NLPIFREE              ((*nlpifree)),               /**< free NLPI user data */
   SCIP_NLPIDATA*                  nlpidata                     /**< NLP interface local data */
);

/** initializes an NLP interface structure
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - name problem name
 */
SCIP_DECL_NLPIINIT( SCIPnlpiInit );

/** frees NLPI user data */
extern 
SCIP_RETCODE SCIPnlpiFree(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPI**           nlpi                /**< pointer to NLPI data structure */
);

/** add variables
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - nvars number of variables 
 *  - lb lower bounds of variables
 *  - ub upper bounds of variables
 *  - type types of variables, saying NULL means all are continuous
 *  - varnames names of variables, can be NULL
 */
SCIP_DECL_NLPIADDVARS( SCIPnlpiAddVars );

/** add constraints 
 * linear coefficients: row(=constraint) oriented matrix
 * quadratic coefficiens: row oriented matrix for each constraint
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - ncons number of added constraints
 *  - linoffset start index of each constraints linear coefficients in linind and linval
 *    length: ncons + 1, linoffset[ncons] gives length of linind and linval
 *    may be NULL in case of no linear part
 *  - linind variable indices
 *    may be NULL in case of no linear part
 *  - linval coefficient values
 *    may be NULL in case of no linear part
 *  - nquadrows number of columns in matrix of quadratic part for each constraint
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadrowidx indices of variables for which a quadratic part is specified
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadoffset start index of each rows quadratic coefficients in quadind[.] and quadval[.]
 *    indices are given w.r.t. quadrowidx., i.e., quadoffset[.][i] gives the start index of row quadrowidx[.][i] in quadval[.]
 *    quadoffset[.][nquadrows[.]] gives length of quadind[.] and quadval[.]
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadind column indices w.r.t. quadrowidx, i.e., quadrowidx[quadind[.][i]] gives the index of the variable corresponding to entry i
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadval coefficient values
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - exprvaridx indices of variables in expression tree, maps variable indices in expression tree to indices in nlp
 *    entry of array may be NULL in case of no expression tree
 *    may be NULL in case of no expression tree in any constraint
 *  - exprtree expression tree for nonquadratic part of constraints
 *    entry of array may be NULL in case of no nonquadratic part
 *    may be NULL in case of no nonquadratic part in any constraint
 *  - names of constraints, may be NULL or entries may be NULL
 */
SCIP_DECL_NLPIADDCONSTRAINTS( SCIPnlpiAddConstraints );

/** sets or overwrites objective, a minimization problem is expected
 *  May change sparsity pattern.
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - linind variable indices
 *    may be NULL in case of no linear part
 *  - linval coefficient values
 *    may be NULL in case of no linear part
 *  - nquadcols number of columns in matrix of quadratic part
 *  - quadcols indices of variables for which a quadratic part is specified
 *    may be NULL in case of no quadratic part
 *  - quadoffset start index of each rows quadratic coefficients in quadind and quadval
 *    quadoffset[.][nquadcols] gives length of quadind and quadval
 *    may be NULL in case of no quadratic part
 *  - quadind column indices
 *    may be NULL in case of no quadratic part
 *  - quadval coefficient values
 *    may be NULL in case of no quadratic part
 *  - exprvaridx indices of variables in expression tree, maps variable indices in expression tree to indices in nlp
 *    may be NULL in case of no expression tree
 *  - exprtree expression tree for nonquadratic part of constraints
 *    may be NULL in case of no nonquadratic part
 *  - objective values offset
 */
SCIP_DECL_NLPISETOBJECTIVE( SCIPnlpiSetObjective );

/** change variable bounds
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - nvars number of variables to change bounds
 *  - indices indices of variables to change bounds
 *  - lb new lower bounds
 *  - ub new upper bounds
 */
SCIP_DECL_NLPICHGVARBOUNDS( SCIPnlpiChgVarBounds );

/** change constraint bounds
 *
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - ncons number of constraints to change bounds
 *  - indices indices of constraints to change bounds
 *  - lb new lower bounds
 *  - ub new upper bounds
 */
SCIP_DECL_NLPICHGCONSBOUNDS( SCIPnlpiChgConsBounds );

/** delete a set of variables
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - dstat deletion status of vars; 1 if var should be deleted, 0 if not
 * 
 * output:
 *  - dstat new position of var, -1 if var was deleted
 */
SCIP_DECL_NLPIDELVARSET( SCIPnlpiDelVarSet );

/** delete a set of constraints
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - dstat deletion status of rows; 1 if row should be deleted, 0 if not
 * 
 * output:
 *  - dstat new position of row, -1 if row was deleted
 */
SCIP_DECL_NLPIDELCONSSET( SCIPnlpiDelConsSet );

/** change one linear coefficient in a constraint or objective
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - cons   index of constraint or -1 for objective
 *  - nvals number of values in linear constraint
 *  - varidx index of variable
 *  - value  new value for coefficient
 * 
 * return: Error if coefficient did not exist before
 */
SCIP_DECL_NLPICHGLINEARCOEFS( SCIPnlpiChgLinearCoefs );
  
/** change one coefficient in the quadratic part of a constraint or objective
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - cons index of constraint or -1 for objective
 *  - row row offset containing modified indices
 *  - col cols containing modified indices to the corresponding row offset
 *  - values coefficients corresponding to same indices as used when constraint/objective was constructed
 * 
 * return: Error if coefficient did not exist before
 */
SCIP_DECL_NLPICHGQUADCOEFS( SCIPnlpiChgQuadCoefs );

/** change a parameter constant in the nonlinear part
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - cons index of constraint or -1 for objective
 *  - idx index of parameter
 *  - value new value for nonlinear parameter
 * 
 * return: Error if parameter does not exist
 */
SCIP_DECL_NLPICHGNONLINCOEF( SCIPnlpiChgNonlinCoef );

/** sets initial guess for primal variables
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - values initial starting solution
 */
SCIP_DECL_NLPISETINITIALGUESS( SCIPnlpiSetInitialGuess );

/** tries to solve NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 */
SCIP_DECL_NLPISOLVE( SCIPnlpiSolve );

/** gives solution status
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 * 
 * return: Solution Status
 */
SCIP_DECL_NLPIGETSOLSTAT( SCIPnlpiGetSolstat );

/** gives termination reason
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 * 
 * return: Termination Status
 */
SCIP_DECL_NLPIGETTERMSTAT( SCIPnlpiGetTermstat );

/** gives primal solution
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - primalvalues buffer to store pointer to primal values
 * 
 * output:
 *  - primalvalues primal values of solution
 */
SCIP_DECL_NLPIGETSOLUTION( SCIPnlpiGetSolution );

/** gives solve statistics
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - statistics buffer to store statistics
 * 
 * output:
 *  - statistics solve statistics
 */
SCIP_DECL_NLPIGETSTATISTICS( SCIPnlpiGetStatistics );

/** gives required size of a buffer to store a warmstart object
 * 
 *  input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - size pointer to store required size for warmstart buffer
 * 
 * output:
 *  - size required size for warmstart buffer
 */
SCIP_DECL_NLPIGETWARMSTARTSIZE( SCIPnlpiGetWarmstartSize );

/** stores warmstart information in buffer
 * 
 * required size of buffer should have been obtained by SCIPnlpiGetWarmstartSize before
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - buffer memory to store warmstart information
 * 
 * output:
 *  - buffer warmstart information in solver specific data structure
 */
SCIP_DECL_NLPIGETWARMSTARTMEMO( SCIPnlpiGetWarmstartMemo );

/** sets warmstart information in solver
 * 
 * write warmstart to buffer
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - buffer warmstart information
 */
SCIP_DECL_NLPISETWARMSTARTMEMO( SCIPnlpiSetWarmstartMemo );

/** gets pointer for NLP solver
 * 
 *  to do dirty stuff
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  
 * return: void pointer to solver
 */
SCIP_DECL_NLPIGETSOLVERPOINTER( SCIPnlpiGetSolverPointer );

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi NLP interface structure
 *  - type parameter number
 *  - ival buffer to store the parameter value
 * 
 * output:
 *  - ival parameter value
 */
SCIP_DECL_NLPIGETINTPAR( SCIPnlpiGetIntPar );

/** sets integer parameter of NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi NLP interface structure
 *  - type parameter number
 *  - ival parameter value
 */
SCIP_DECL_NLPISETINTPAR( SCIPnlpiSetIntPar );

/** gets floating point parameter of NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi NLP interface structure
 *  - type parameter number
 *  - dval buffer to store the parameter value
 * 
 * output:
 *  - dval parameter value
 */
SCIP_DECL_NLPIGETREALPAR( SCIPnlpiGetRealPar );

/** sets floating point parameter of NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi NLP interface structure
 *  - type parameter number
 *  - dval parameter value
 */
SCIP_DECL_NLPISETREALPAR( SCIPnlpiSetRealPar );

/** get nlpi data
 */
SCIP_NLPIDATA* SCIPnlpiGetNlpiData(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
);

/** gets NLP solver name
 */
const char* SCIPnlpiGetName(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
);

/** Creates an NLP statistics structure.
 */
SCIP_RETCODE SCIPnlpStatisticsCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPSTATISTICS**  statistics          /**< pointer where to store NLP statistics structure */
);

/** Frees an NLP statistics structure.
 */
void SCIPnlpStatisticsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPSTATISTICS**  statistics          /**< pointer where to store NLP statistics structure */
);

/** Gets the number of iterations from an NLP statistics structure.
 */
int SCIPnlpStatisticsGetNIterations(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
);

/** Gets the total time from an NLP statistics structure.
 */
SCIP_Real SCIPnlpStatisticsGetTotalTime(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
);

/** Sets the number of iterations in an NLP statistics structure.
 */
void SCIPnlpStatisticsSetNIterations(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   int                   niterations         /**< number of iterations to store */
);

/** Sets the total time in an NLP statistics structure.
 */
void SCIPnlpStatisticsSetTotalTime(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   SCIP_Real             totaltime           /**< solution time to store */
);

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_H__ */
