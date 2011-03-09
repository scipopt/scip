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

   assert(nlpiinit != NULL);
   assert(nlpifree != NULL);
   assert(nlpiaddvars != NULL);
   assert(nlpiaddconstraints != NULL);
   assert(nlpisetobjective != NULL);
   assert(nlpichgvarbounds != NULL);
   assert(nlpichgconsbounds != NULL);
   assert(nlpidelconsset != NULL);
   assert(nlpichglinearcoefs != NULL);
   assert(nlpichgquadcoefs != NULL);
   assert(nlpichgnonlincoef != NULL);
   assert(nlpisetinitialguess != NULL);
   assert(nlpisolve != NULL);
   assert(nlpigetsolstat != NULL);
   assert(nlpigettermstat != NULL);
   assert(nlpigetsolution != NULL);
   assert(nlpigetstatistics != NULL);
   assert(nlpigetwarmstartsize != NULL);
   assert(nlpigetwarmstartmemo != NULL);
   assert(nlpisetwarmstartmemo != NULL);
   assert(nlpigetsolverpointer != NULL);
   assert(nlpigetintpar != NULL);
   assert(nlpisetintpar != NULL);
   assert(nlpigetrealpar != NULL);
   assert(nlpisetrealpar != NULL);

   SCIP_CALL( SCIPallocMemory(scip, nlpi) );
   (*nlpi)->nlpiinit = nlpiinit;
   (*nlpi)->nlpifree = nlpifree;
   (*nlpi)->nlpiaddvars = nlpiaddvars;
   (*nlpi)->nlpiaddconstraints = nlpiaddconstraints;
   (*nlpi)->nlpisetobjective = nlpisetobjective;
   (*nlpi)->nlpichgvarbounds = nlpichgvarbounds;
   (*nlpi)->nlpichgconsbounds = nlpichgconsbounds;
   (*nlpi)->nlpidelvarset = nlpidelvarset;
   (*nlpi)->nlpidelconsset = nlpidelconsset;
   (*nlpi)->nlpichglinearcoefs = nlpichglinearcoefs;
   (*nlpi)->nlpichgquadcoefs = nlpichgquadcoefs;
   (*nlpi)->nlpichgnonlincoef = nlpichgnonlincoef;
   (*nlpi)->nlpisetinitialguess = nlpisetinitialguess;
   (*nlpi)->nlpisolve = nlpisolve;
   (*nlpi)->nlpigetsolstat = nlpigetsolstat;
   (*nlpi)->nlpigettermstat = nlpigettermstat;
   (*nlpi)->nlpigetsolution = nlpigetsolution;
   (*nlpi)->nlpigetstatistics = nlpigetstatistics;
   (*nlpi)->nlpigetwarmstartsize = nlpigetwarmstartsize;
   (*nlpi)->nlpigetwarmstartmemo = nlpigetwarmstartmemo;
   (*nlpi)->nlpisetwarmstartmemo = nlpisetwarmstartmemo;
   (*nlpi)->nlpigetsolverpointer = nlpigetsolverpointer;
   (*nlpi)->nlpigetintpar = nlpigetintpar;
   (*nlpi)->nlpisetintpar = nlpisetintpar;
   (*nlpi)->nlpigetrealpar = nlpigetrealpar;
   (*nlpi)->nlpisetrealpar = nlpisetrealpar;
   (*nlpi)->nlpidata = nlpidata;

   return SCIP_OKAY;
}

/** frees NLPI user data */
SCIP_RETCODE SCIPnlpiFree(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPI**           nlpi                /**< pointer to NLPI data structure */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(*nlpi != NULL);

   SCIP_CALL( (*(*nlpi)->nlpifree)(scip, (*nlpi)->nlpidata) );
   SCIPfreeMemory(scip, nlpi);

   return SCIP_OKAY;
}

/** initializes an NLP interface structure */
SCIP_RETCODE SCIPnlpiInit(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI data structure */
   const char*           name                /**< problem name */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   if( nlpi->nlpiinit != NULL )
   {
      SCIP_CALL( (*nlpi->nlpiinit)(scip, nlpi, name) );
   }

   return SCIP_OKAY;
}

/** add variables to nlpi */
SCIP_RETCODE SCIPnlpiAddVars(
   SCIP*                 scip,               /**< pointer to SCIP */                 
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI data structure */
   int                   nvars,              /**< number of variables */               
   const SCIP_Real*      lbs,                /**< lower bounds of variables */
   const SCIP_Real*      ubs,                /**< ubs upper bounds of variables */
   SCIP_VARTYPE*         types,              /**< types of variables, saying NULL means all are continuous */
   const char**          varnames            /**< varnames names of variables, can be NULL */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(lbs != NULL);
   assert(ubs != NULL);
   
   SCIP_CALL( (*nlpi->nlpiaddvars)(scip, nlpi, nvars, lbs, ubs, types, varnames) );
   
   return SCIP_OKAY;
}

/** add constraints to nlpi */
SCIP_RETCODE SCIPnlpiAddConstraints(
   SCIP*                 scip,               /**< pointer to SCIP */                 
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI data structure */
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
   int* const*           quadoffsets,        /**< start index of each rows quadratic coefficients in quadind[.] and quadval[.],
                                              * indices are given w.r.t. quadrowidx., i.e., quadoffset[.][i] gives the start
                                              * index of row quadrowidx[.][i] in quadval[.], quadoffset[.][nquadrows[.]] gives
                                              * length of quadind[.] and quadval[.], entry of array may be NULL in case of no
                                              * quadratic part, may be NULL in case of no quadratic part in any constraint */
   int* const*           quadinds,           /**< column indices w.r.t. quadrowidx, i.e., quadrowidx[quadind[.][i]] gives the
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
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpiaddconstraints)(scip, nlpi, nconss, lhss, rhss, nlininds, lininds, linvals,
      nquadrows, quadrowidxs, quadoffsets, quadinds, quadvals, exprvaridxs, exprtrees, names) );
   
   return SCIP_OKAY;
}

/** sets or overwrites objective, a minimization problem is expected */
SCIP_RETCODE SCIPnlpiSetObjective(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   int                   nlins,              /**< number of linear variables */
   const int*            lininds,            /**< variable indices, may be NULL in case of no linear part */
   const SCIP_Real*      linvals,            /**< coefficient values, may be NULL in case of no linear part */
   int                   nquadcols,          /**< number of columns in matrix of quadratic part */
   const int*            quadcols,           /**< indices of variables for which a quadratic part is specified, may be NULL in
                                              * case of no quadratic part */
   const int*            quadoffsets,        /**< start index of each rows quadratic coefficients in quadind and quadval,
                                              * quadoffset[.][nquadcols] gives length of quadind and quadval, may be NULL in
                                              * case of no quadratic part */
   const int*            quadinds,           /**< column indices, may be NULL in case of no quadratic part */
   const SCIP_Real*      quadvals,           /**< coefficient values, may be NULL in case of no quadratic part */
   const int*            exprvaridxs,        /**< indices of variables in expression tree, maps variable indices in expression
                                              * tree to indices in nlp, may be NULL in case of no expression tree */
   const SCIP_EXPRTREE*  exprtree,           /**< expression tree for nonquadratic part of objective function, may be NULL in
                                              * case of no nonquadratic part */
   const SCIP_Real       constant            /**< objective value offset*/
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetobjective)(scip, nlpi, nlins, lininds, linvals, nquadcols, quadcols, 
      quadoffsets, quadinds, quadvals, exprvaridxs, exprtree, constant) );
   
   return SCIP_OKAY;
}

/** change variable bounds */
SCIP_RETCODE SCIPnlpiChgVarBounds(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   const int             nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds */
   const SCIP_Real*      ubs                 /**< new upper bounds */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpichgvarbounds)(scip, nlpi, nvars, indices, lbs, ubs) );
   
   return SCIP_OKAY;
}

/** change constraint bounds */
SCIP_RETCODE SCIPnlpiChgConsBounds(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   const int             nconss,              /**< number of constraints to change bounds */
   const int*            indices,            /**< indices of constraints to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds */
   const SCIP_Real*      ubs                 /**< new upper bounds */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpichgconsbounds)(scip, nlpi, nconss, indices, lbs, ubs) );
   
   return SCIP_OKAY;
}

/** delete a set of variables */
SCIP_RETCODE SCIPnlpiDelVarSet(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   int*                  dstats              /**< deletion status of vars; 1 if var should be deleted, 0 if not; afterwards -1
                                              * if var was deleted */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpidelvarset)(scip, nlpi, dstats) );
   
   return SCIP_OKAY;
}

/** delete a set of constraints */
SCIP_RETCODE SCIPnlpiDelConsSet(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   int*                  dstats              /**< deletion status of rows; 1 if row should be deleted, 0 if not; afterwards -1
                                              * if row was deleted */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpidelconsset)(scip, nlpi, dstats) );
   
   return SCIP_OKAY;
}

/** change one linear coefficient in a constraint or objective; returns: Error if coefficient did not exist before */
SCIP_RETCODE SCIPnlpiChgLinearCoefs(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   const int             idx,                /**< index of constraint or -1 for objective */
   int                   nvals,              /**< number of values in linear constraint */
   const int*            varidxs,            /**< indices of variable */
   const SCIP_Real*      vals                /**< new values for coefficient */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpichglinearcoefs)(scip, nlpi, idx, nvals, varidxs, vals) );
   
   return SCIP_OKAY;
}
  
/** change one coefficient in the quadratic part of a constraint or objective; return: Error if coefficient did not exist before */
SCIP_RETCODE SCIPnlpiChgQuadCoefs(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   const int             idx,                /**< index of constraint or -1 for objective */
   const int             nentries,           /**< nentries number of values in quadratic constraint */
   const int*            rows,               /**< row offset containing modified indices */
   const int*            cols,               /**< cols containing modified indices to the corresponding row offset */
   SCIP_Real*            values              /**< coefficients corresponding to same indices as used when constraint/objective
                                              * was constructed */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpichgquadcoefs)(scip, nlpi, idx, nentries, rows, cols, values) );
   
   return SCIP_OKAY;
}

/** change one coefficient in the nonlinear part; return: Error if parameter does not exist */
SCIP_RETCODE SCIPnlpiChgNonlinCoef(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   const int             considx,            /**< index of constraint or -1 for objective */
   const int             paramidx,           /**< index of parameter */
   SCIP_Real             value               /**< new value for nonlinear parameter */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpichgnonlincoef)(scip, nlpi, considx, paramidx, value) );
   
   return SCIP_OKAY;
}

/** sets initial guess for primal variables */
SCIP_RETCODE SCIPnlpiSetInitialGuess(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_Real*            values              /**< initial starting solution */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetinitialguess)(scip, nlpi, values) );
   
   return SCIP_OKAY;
}

/** tries to solve NLP */
SCIP_RETCODE SCIPnlpiSolve(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi                /**< pointer to NLPI datastructure */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpisolve)(scip, nlpi) );
   
   return SCIP_OKAY;
}

/** gives solution status, return: Solution Status */
SCIP_NLPSOLSTAT SCIPnlpiGetSolstat(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi                /**< pointer to NLPI datastructure */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   return (*nlpi->nlpigetsolstat)(scip, nlpi);
}

/** gives termination reason; return: Termination Status */
SCIP_NLPTERMSTAT SCIPnlpiGetTermstat(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi                /**< pointer to NLPI datastructure */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   return (*nlpi->nlpigettermstat)(scip, nlpi);
}

/** gives primal solution */
SCIP_RETCODE SCIPnlpiGetSolution(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_Real**           primalvalues        /**< pointer to store primal values */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetsolution)(scip, nlpi, primalvalues) );
   
   return SCIP_OKAY;
}

/** gives solve statistics */
SCIP_RETCODE SCIPnlpiGetStatistics(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetstatistics)(scip, nlpi, statistics) );
   
   return SCIP_OKAY;
}

/** gives required size of a buffer to store a warmstart object */
SCIP_RETCODE SCIPnlpiGetWarmstartSize(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   size_t*               size                /**< pointer to store required size for warmstart buffer */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetwarmstartsize)(scip, nlpi, size) );
   
   return SCIP_OKAY;
}

/** stores warmstart information in buffer */
SCIP_RETCODE SCIPnlpiGetWarmstartMemo(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   void*                 buffer              /**< memory to store warmstart information */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetwarmstartmemo)(scip, nlpi, buffer) );
   
   return SCIP_OKAY;
}

/** sets warmstart information in solver */
SCIP_RETCODE SCIPnlpiSetWarmstartMemo(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   void*                 buffer              /**< warmstart information */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetwarmstartmemo)(scip, nlpi, buffer) );
   
   return SCIP_OKAY;
}

/** gets pointer for NLP solver, return: void pointer to solver */
void* SCIPnlpiGetSolverPointer(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi                /**< pointer to NLPI datastructure */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   return (*nlpi->nlpigetsolverpointer)(scip, nlpi);
}

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP */
SCIP_RETCODE SCIPnlpiGetIntPar(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(ival != NULL);
      
   SCIP_CALL( (*nlpi->nlpigetintpar)(scip, nlpi, type, ival) );
   
   return SCIP_OKAY;
}

/** sets integer parameter of NLP */
SCIP_RETCODE SCIPnlpiSetIntPar(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetintpar)(scip, nlpi, type, ival) );
   
   return SCIP_OKAY;
}

/** gets floating point parameter of NLP */
SCIP_RETCODE SCIPnlpiGetRealPar(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(dval != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetrealpar)(scip, nlpi, type, dval) );
   
   return SCIP_OKAY;
}

/** sets floating point parameter of NLP */
SCIP_RETCODE SCIPnlpiSetRealPar(
   SCIP*                 scip,               /**< pointer to SCIP */              
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetrealpar)(scip, nlpi, type, dval) );
   
   return SCIP_OKAY;
}

/** get nlpi data */
SCIP_NLPIDATA* SCIPnlpiGetNlpiData(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
)
{
   assert(nlpi != NULL);
   
   return nlpi->nlpidata;
}

/** gets NLP solver name */
const char* SCIPnlpiGetName(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
)
{
   assert(nlpi != NULL);
   
   return nlpi->name;
}

/** Creates an NLP statistics structure. */
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
   (*statistics)->totaltime = -1.0;
   
   return SCIP_OKAY;
}

/** Frees an NLP statistics structure. */
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

/** Gets the number of iterations from an NLP statistics structure. */
int SCIPnlpStatisticsGetNIterations(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
)
{
   assert(statistics != NULL);
   
   return statistics->niterations;
}

/** Gets the total time from an NLP statistics structure. */
SCIP_Real SCIPnlpStatisticsGetTotalTime(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
)
{
   assert(statistics != NULL);
   
   return statistics->totaltime;
}

/** Sets the number of iterations in an NLP statistics structure. */
void SCIPnlpStatisticsSetNIterations(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   int                   niterations         /**< number of iterations to store */
)
{
   assert(statistics != NULL);
   statistics->niterations = niterations;
}

/** Sets the total time in an NLP statistics structure. */
void SCIPnlpStatisticsSetTotalTime(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   SCIP_Real             totaltime           /**< solution time to store */
)
{
   assert(statistics != NULL);
   statistics->totaltime = totaltime;
}
