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
#pragma ident "@(#) $Id: nlpi.c,v 1.7 2009/09/08 17:05:23 bzfviger Exp $"

/**@file   nlpi.c
 * @brief  methods for handling nlp interface
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>

#include "scip/scip.h"
#include "scip/nlpi.h"
#include "scip/struct_nlpi.h"

/** creates an NLP solver interface */
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
   SCIP_DECL_NLPICHGLINEARCOEFS    ((*nlpichglinearcoefs)),     /**< change one coefficient  in linear part */
   SCIP_DECL_NLPICHGQUADCOEFS      ((*nlpichgquadcoefs)),       /**< change one coefficient  in quadratic part */
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
)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nlpi != NULL);

   assert(nlpiinit             != NULL);
   assert(nlpifree             != NULL);
   assert(nlpiaddvars          != NULL);
   assert(nlpiaddconstraints   != NULL);
   assert(nlpisetobjective     != NULL);
   assert(nlpichgvarbounds     != NULL);
   assert(nlpichgconsbounds    != NULL);
   assert(nlpidelconsset       != NULL);
   assert(nlpichglinearcoefs   != NULL);
   assert(nlpichgquadcoefs     != NULL);
   assert(nlpichgnonlincoef    != NULL);
   assert(nlpisetinitialguess  != NULL);
   assert(nlpisolve            != NULL);
   assert(nlpigetsolstat       != NULL);
   assert(nlpigettermstat      != NULL);
   assert(nlpigetsolution      != NULL);
   assert(nlpigetstatistics    != NULL);
   assert(nlpigetwarmstartsize != NULL);
   assert(nlpigetwarmstartmemo != NULL);
   assert(nlpisetwarmstartmemo != NULL);
   assert(nlpigetsolverpointer != NULL);
   assert(nlpigetintpar        != NULL);
   assert(nlpisetintpar        != NULL);
   assert(nlpigetrealpar       != NULL);
   assert(nlpisetrealpar       != NULL);

   SCIP_CALL( SCIPallocMemory( scip, nlpi ) );
   (*nlpi) -> nlpiinit             = nlpiinit;
   (*nlpi) -> nlpifree             = nlpifree;
   (*nlpi) -> nlpiaddvars          = nlpiaddvars;
   (*nlpi) -> nlpiaddconstraints   = nlpiaddconstraints;
   (*nlpi) -> nlpisetobjective     = nlpisetobjective;
   (*nlpi) -> nlpichgvarbounds     = nlpichgvarbounds;
   (*nlpi) -> nlpichgconsbounds    = nlpichgconsbounds;
   (*nlpi) -> nlpidelvarset        = nlpidelvarset;
   (*nlpi) -> nlpidelconsset       = nlpidelconsset;
   (*nlpi) -> nlpichglinearcoefs   = nlpichglinearcoefs;
   (*nlpi) -> nlpichgquadcoefs     = nlpichgquadcoefs;
   (*nlpi) -> nlpichgnonlincoef    = nlpichgnonlincoef;
   (*nlpi) -> nlpisetinitialguess  = nlpisetinitialguess;
   (*nlpi) -> nlpisolve            = nlpisolve;
   (*nlpi) -> nlpigetsolstat       = nlpigetsolstat;
   (*nlpi) -> nlpigettermstat      = nlpigettermstat;
   (*nlpi) -> nlpigetsolution      = nlpigetsolution;
   (*nlpi) -> nlpigetstatistics    = nlpigetstatistics;
   (*nlpi) -> nlpigetwarmstartsize = nlpigetwarmstartsize;
   (*nlpi) -> nlpigetwarmstartmemo = nlpigetwarmstartmemo;
   (*nlpi) -> nlpisetwarmstartmemo = nlpisetwarmstartmemo;
   (*nlpi) -> nlpigetsolverpointer = nlpigetsolverpointer;
   (*nlpi) -> nlpigetintpar        = nlpigetintpar;
   (*nlpi) -> nlpisetintpar        = nlpisetintpar;
   (*nlpi) -> nlpigetrealpar       = nlpigetrealpar;
   (*nlpi) -> nlpisetrealpar       = nlpisetrealpar;
   (*nlpi) -> nlpidata             = nlpidata;

   return SCIP_OKAY;
}

/** frees NLPI user data */
SCIP_RETCODE SCIPnlpiFree(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPI**           nlpi                /**< pointer to NLPI data structure */
   )
{
   assert(nlpi  != NULL);
   assert(*nlpi != NULL);
   assert(scip  != NULL);

   SCIP_CALL( (*(*nlpi)->nlpifree)(scip, (*nlpi)->nlpidata) );
   SCIPfreeMemory( scip, nlpi );

   return SCIP_OKAY;
}

/** initializes an NLP interface structure
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - name problem name
 */
SCIP_DECL_NLPIINIT( SCIPnlpiInit )
{
   assert(nlpi != NULL);
   
   if( nlpi -> nlpiinit != NULL )
   {
      SCIP_CALL( (*nlpi->nlpiinit)(scip, nlpi, name) );
   }

   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIADDVARS( SCIPnlpiAddVars )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpiaddvars)(scip, nlpi, nvars, lb, ub, type, varnames) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIADDCONSTRAINTS( SCIPnlpiAddConstraints )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpiaddconstraints)(scip, nlpi, ncons, lhs, rhs, linoffset, linind, linval,
      nquadrows, quadrowidx, quadoffset, quadind, quadval, exprvaridx, exprtree, names) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPISETOBJECTIVE( SCIPnlpiSetObjective )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetobjective)(scip, nlpi, nlin, linind, linval, nquadcols, quadcols, 
      quadoffset, quadind, quadval, exprvaridx, exprtree, constant) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPICHGVARBOUNDS( SCIPnlpiChgVarBounds )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpichgvarbounds)(scip, nlpi, nvars, indices, lb, ub) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPICHGCONSBOUNDS( SCIPnlpiChgConsBounds )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpichgconsbounds)(scip, nlpi, ncons, indices, lb, ub) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIDELVARSET( SCIPnlpiDelVarSet )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpidelvarset)(scip, nlpi, dstat) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIDELCONSSET( SCIPnlpiDelConsSet )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpidelconsset)(scip, nlpi, dstat) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPICHGLINEARCOEFS( SCIPnlpiChgLinearCoefs )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpichglinearcoefs)(scip, nlpi, cons, nvals, varidx, value) );
   
   return SCIP_OKAY;
}
  
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
SCIP_DECL_NLPICHGQUADCOEFS( SCIPnlpiChgQuadCoefs )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpichgquadcoefs)(scip, nlpi, cons, nentries, row, col, value) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPICHGNONLINCOEF( SCIPnlpiChgNonlinCoef )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpichgnonlincoef)(scip, nlpi, cons, idx, value) );
   
   return SCIP_OKAY;
}

/** sets initial guess for primal variables
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - values initial starting solution
 */
SCIP_DECL_NLPISETINITIALGUESS( SCIPnlpiSetInitialGuess )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetinitialguess)(scip, nlpi, values) );
   
   return SCIP_OKAY;
}

/** tries to solve NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 */
SCIP_DECL_NLPISOLVE( SCIPnlpiSolve )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpisolve)(scip, nlpi) );
   
   return SCIP_OKAY;
}

/** gives solution status
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 * 
 * return: Solution Status
 */
SCIP_DECL_NLPIGETSOLSTAT( SCIPnlpiGetSolstat )
{
  assert(nlpi != NULL);
  
  return (*nlpi->nlpigetsolstat)(scip, nlpi);
}

/** gives termination reason
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 * 
 * return: Termination Status
 */
SCIP_DECL_NLPIGETTERMSTAT( SCIPnlpiGetTermstat )
{
   assert(nlpi != NULL);
   
   return (*nlpi->nlpigettermstat)(scip, nlpi);
}

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
SCIP_DECL_NLPIGETSOLUTION( SCIPnlpiGetSolution )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetsolution)(scip, nlpi, primalvalues) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETSTATISTICS( SCIPnlpiGetStatistics )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetstatistics)(scip, nlpi, statistics) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETWARMSTARTSIZE( SCIPnlpiGetWarmstartSize )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetwarmstartsize)(scip, nlpi, size) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETWARMSTARTMEMO( SCIPnlpiGetWarmstartMemo )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetwarmstartmemo)(scip, nlpi, buffer) );
   
   return SCIP_OKAY;
}

/** sets warmstart information in solver
 * 
 * write warmstart to buffer
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - buffer warmstart information
 */
SCIP_DECL_NLPISETWARMSTARTMEMO( SCIPnlpiSetWarmstartMemo )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetwarmstartmemo)(scip, nlpi, buffer) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETSOLVERPOINTER( SCIPnlpiGetSolverPointer )
{
   assert(nlpi != NULL);
   
   return (*nlpi->nlpigetsolverpointer)(scip, nlpi);
}

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
SCIP_DECL_NLPIGETINTPAR( SCIPnlpiGetIntPar )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetintpar)(scip, nlpi, type, ival) );
   
   return SCIP_OKAY;
}

/** sets integer parameter of NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi NLP interface structure
 *  - type parameter number
 *  - ival parameter value
 */
SCIP_DECL_NLPISETINTPAR( SCIPnlpiSetIntPar )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetintpar)(scip, nlpi, type, ival) );
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETREALPAR( SCIPnlpiGetRealPar )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetrealpar)(scip, nlpi, type, dval) );
   
   return SCIP_OKAY;
}

/** sets floating point parameter of NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi NLP interface structure
 *  - type parameter number
 *  - dval parameter value
 */
SCIP_DECL_NLPISETREALPAR( SCIPnlpiSetRealPar )
{
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetrealpar)(scip, nlpi, type, dval) );
   
   return SCIP_OKAY;
}

/** get nlpi data
 */
SCIP_NLPIDATA* SCIPnlpiGetNlpiData(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
)
{
   assert(nlpi != NULL);
   
   return nlpi->nlpidata;
}

/** gets NLP solver name
 */
const char* SCIPnlpiGetName(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
)
{
   assert(nlpi != NULL);
   
   return nlpi->name;
}

/** Creates an NLP statistics structure.
 */
SCIP_RETCODE SCIPnlpStatisticsCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPSTATISTICS**  statistics          /**< pointer where to store NLP statistics structure */
)
{
   assert(scip != NULL);
   assert(statistics != NULL);
   
   SCIP_CALL( SCIPallocMemory(scip, statistics) );
   assert(*statistics != NULL);
   
   (*statistics)->niterations = -1;
   (*statistics)->totaltime   = -1.0;
   
   return SCIP_OKAY;
}

/** Frees an NLP statistics structure.
 */
void SCIPnlpStatisticsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPSTATISTICS**  statistics          /**< pointer where to store NLP statistics structure */
)
{
   assert(scip != NULL);
   assert(statistics != NULL);
   
   SCIPfreeMemory(scip, statistics);
   assert(*statistics == NULL);
}

/** Gets the number of iterations from an NLP statistics structure.
 */
int SCIPnlpStatisticsGetNIterations(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
)
{
   assert(statistics != NULL);
   
   return statistics->niterations;
}

/** Gets the total time from an NLP statistics structure.
 */
SCIP_Real SCIPnlpStatisticsGetTotalTime(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
)
{
   assert(statistics != NULL);
   
   return statistics->totaltime;
}

/** Sets the number of iterations in an NLP statistics structure.
 */
void SCIPnlpStatisticsSetNIterations(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   int                   niterations         /**< number of iterations to store */
)
{
   assert(statistics != NULL);
   statistics->niterations = niterations;
}

/** Sets the total time in an NLP statistics structure.
 */
void SCIPnlpStatisticsSetTotalTime(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   SCIP_Real             totaltime           /**< solution time to store */
)
{
   assert(statistics != NULL);
   statistics->totaltime = totaltime;
}
