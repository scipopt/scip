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
#pragma ident "@(#) $Id: nlpi.c,v 1.3 2010/05/03 15:23:57 bzfviger Exp $"

/**@file   nlpi.c
 * @brief  methods for handling nlp interface
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "nlpi/nlpi.h"
#include "nlpi/struct_nlpi.h"
#include "blockmemshell/memory.h"

/** compares two NLPIs w.r.t. their priority */
SCIP_DECL_SORTPTRCOMP(SCIPnlpiComp)
{  /*lint --e{715}*/
   return ((SCIP_NLPI*)elem2)->priority - ((SCIP_NLPI*)elem1)->priority;
}

/** creates an NLP solver interface */
SCIP_RETCODE SCIPnlpiCreate(
   SCIP_NLPI**                     nlpi,                        /**< pointer to NLP interface data structure */
   const char*                     name,                        /**< name of NLP interface */
   const char*                     description,                 /**< description of NLP interface */
   int                             priority,                    /**< priority of NLP interface */
   SCIP_DECL_NLPICOPY              ((*nlpicopy)),               /**< copying an NLPI */
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
   SCIP_DECL_NLPIGETINTPAR         ((*nlpigetintpar)),          /**< get value of integer parameter */
   SCIP_DECL_NLPISETINTPAR         ((*nlpisetintpar)),          /**< set value of integer parameter */
   SCIP_DECL_NLPIGETREALPAR        ((*nlpigetrealpar)),         /**< get value of floating point parameter */
   SCIP_DECL_NLPISETREALPAR        ((*nlpisetrealpar)),         /**< set value of floating point parameter */
   SCIP_DECL_NLPIGETSTRINGPAR      ((*nlpigetstringpar)),       /**< get value of string parameter */
   SCIP_DECL_NLPISETSTRINGPAR      ((*nlpisetstringpar)),       /**< set value of string parameter */
   SCIP_NLPIDATA*                  nlpidata                     /**< NLP interface local data */
)
{  /*lint --e{715}*/
   assert(nlpi != NULL);

   assert(name != NULL);
   assert(description != NULL);
   assert(nlpicopy != NULL);
   assert(nlpifree != NULL);
   assert(nlpigetsolverpointer != NULL);
   assert(nlpicreateproblem != NULL);
   assert(nlpifreeproblem != NULL);
   assert(nlpigetproblempointer != NULL);
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
   assert(nlpigetintpar != NULL);
   assert(nlpisetintpar != NULL);
   assert(nlpigetrealpar != NULL);
   assert(nlpisetrealpar != NULL);
   assert(nlpigetstringpar != NULL);
   assert(nlpisetstringpar != NULL);

   if( BMSallocMemory(nlpi) == NULL )
      return SCIP_NOMEMORY;
   
   (*nlpi)->name = strdup(name);
   (*nlpi)->description = strdup(description);
   (*nlpi)->priority = priority;
   (*nlpi)->nlpicopy = nlpicopy;
   (*nlpi)->nlpifree = nlpifree;
   (*nlpi)->nlpigetsolverpointer = nlpigetsolverpointer;
   (*nlpi)->nlpicreateproblem = nlpicreateproblem;
   (*nlpi)->nlpifreeproblem = nlpifreeproblem;
   (*nlpi)->nlpigetproblempointer = nlpigetproblempointer;
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
   (*nlpi)->nlpigetintpar = nlpigetintpar;
   (*nlpi)->nlpisetintpar = nlpisetintpar;
   (*nlpi)->nlpigetrealpar = nlpigetrealpar;
   (*nlpi)->nlpisetrealpar = nlpisetrealpar;
   (*nlpi)->nlpigetstringpar = nlpigetstringpar;
   (*nlpi)->nlpisetstringpar = nlpisetstringpar;
   (*nlpi)->nlpidata = nlpidata;

   return SCIP_OKAY;
}

/** copies an NLPI */
SCIP_RETCODE SCIPnlpiCopy(
   SCIP_NLPI*            sourcenlpi,         /**< pointer to NLPI data structure to copy */
   SCIP_NLPI**           targetnlpi          /**< buffer to store pointer to copied NLPI data structure */
)
{
   assert(sourcenlpi != NULL);
   assert(targetnlpi != NULL);

   SCIP_CALL( (*sourcenlpi->nlpicopy)(sourcenlpi, targetnlpi) );

   return SCIP_OKAY;
}

/** frees NLPI user data */
SCIP_RETCODE SCIPnlpiFree(
   SCIP_NLPI**           nlpi                /**< pointer to NLPI data structure */
   )
{
   assert(nlpi  != NULL);
   assert(*nlpi != NULL);

   SCIP_CALL( (*(*nlpi)->nlpifree)((*nlpi)) );
   free((*nlpi)->name);
   free((*nlpi)->description);
   BMSfreeMemory(nlpi);
   
   assert(*nlpi == NULL);

   return SCIP_OKAY;
}

/** gets pointer for NLP solver
 * @return void pointer to solver
 */
void* SCIPnlpiGetSolverPointer(
   SCIP_NLPI*            nlpi                /**< pointer to NLPI datastructure */
   )
{
   assert(nlpi != NULL);
   
   return (*nlpi->nlpigetsolverpointer)(nlpi);
}

/** creates a problem instance */
SCIP_RETCODE SCIPnlpiCreateProblem(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM**    problem,            /**< pointer to store problem data */
   const char*           name                /**< name of problem, can be NULL */
)
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   return (*nlpi->nlpicreateproblem)(nlpi, problem, name);
}

/** frees a problem instance */
SCIP_RETCODE SCIPnlpiFreeProblem(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM**    problem             /**< pointer where problem data is stored */
)
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   return (*nlpi->nlpifreeproblem)(nlpi, problem);
}

/** gets pointer to solver-internal problem instance
 * @return void pointer to problem instance
 */
void* SCIPnlpiGetProblemPointer(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer where problem data is stored */
)
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   return (*nlpi->nlpigetproblempointer)(nlpi, problem);
}

/** add variables to nlpi */
SCIP_RETCODE SCIPnlpiAddVars(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI data structure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int                   nvars,              /**< number of variables */               
   const SCIP_Real*      lbs,                /**< lower bounds of variables, can be NULL if -infinity */
   const SCIP_Real*      ubs,                /**< ubs upper bounds of variables, can be NULL if +infinity */
   const char**          varnames            /**< varnames names of variables, can be NULL */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpiaddvars)(nlpi, problem, nvars, lbs, ubs, varnames) );
   
   return SCIP_OKAY;
}

/** add constraints to nlpi */
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
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpiaddconstraints)(nlpi, problem, nconss, lhss, rhss, nlininds, lininds, linvals,
      nquadrows, quadrowidxs, quadoffsets, quadinds, quadvals, exprvaridxs, exprtrees, names) );
   
   return SCIP_OKAY;
}

/** sets or overwrites objective, a minimization problem is expected */
SCIP_RETCODE SCIPnlpiSetObjective(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
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
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetobjective)(nlpi, problem, nlins, lininds, linvals, nquadcols, quadcols, 
      quadoffsets, quadinds, quadvals, exprvaridxs, exprtree, constant) );
   
   return SCIP_OKAY;
}

/** change variable bounds */
SCIP_RETCODE SCIPnlpiChgVarBounds(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   const int             nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds */
   const SCIP_Real*      ubs                 /**< new upper bounds */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpichgvarbounds)(nlpi, problem, nvars, indices, lbs, ubs) );
   
   return SCIP_OKAY;
}

/** change constraint bounds */
SCIP_RETCODE SCIPnlpiChgConsBounds(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   const int             nconss,              /**< number of constraints to change bounds */
   const int*            indices,            /**< indices of constraints to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds */
   const SCIP_Real*      ubs                 /**< new upper bounds */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpichgconsbounds)(nlpi, problem, nconss, indices, lbs, ubs) );
   
   return SCIP_OKAY;
}

/** delete a set of variables */
SCIP_RETCODE SCIPnlpiDelVarSet(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int*                  dstats              /**< deletion status of vars; 1 if var should be deleted, 0 if not; afterwards -1
                                              * if var was deleted */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpidelvarset)(nlpi, problem, dstats) );
   
   return SCIP_OKAY;
}

/** delete a set of constraints */
SCIP_RETCODE SCIPnlpiDelConsSet(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int*                  dstats              /**< deletion status of rows; 1 if row should be deleted, 0 if not; afterwards -1
                                              * if row was deleted */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpidelconsset)(nlpi, problem, dstats) );
   
   return SCIP_OKAY;
}

/** changes or adds one linear coefficient in a constraint or objective */
SCIP_RETCODE SCIPnlpiChgLinearCoefs(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   const int             idx,                /**< index of constraint or -1 for objective */
   int                   nvals,              /**< number of values in linear constraint */
   const int*            varidxs,            /**< indices of variable */
   const SCIP_Real*      vals                /**< new values for coefficient */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpichglinearcoefs)(nlpi, problem, idx, nvals, varidxs, vals) );
   
   return SCIP_OKAY;
}
  
/** changes or adds one coefficient in the quadratic part of a constraint or objective */
SCIP_RETCODE SCIPnlpiChgQuadCoefs(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   const int             idx,                /**< index of constraint or -1 for objective */
   const int             nentries,           /**< nentries number of values in quadratic constraint */
   const int*            rows,               /**< row offset containing modified indices */
   const int*            cols,               /**< cols containing modified indices to the corresponding row offset */
   SCIP_Real*            values              /**< coefficients corresponding to same indices as used when constraint/objective
                                              * was constructed */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpichgquadcoefs)(nlpi, problem, idx, nentries, rows, cols, values) );
   
   return SCIP_OKAY;
}

/** change the value of one parameter in the nonlinear part */
SCIP_RETCODE SCIPnlpiChgNonlinCoef(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   const int             considx,            /**< index of constraint or -1 for objective */
   const int             paramidx,           /**< index of parameter */
   SCIP_Real             value               /**< new value for nonlinear parameter */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpichgnonlincoef)(nlpi, problem, considx, paramidx, value) );
   
   return SCIP_OKAY;
}

/** sets initial guess for primal variables */
SCIP_RETCODE SCIPnlpiSetInitialGuess(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_Real*            values              /**< initial starting solution */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetinitialguess)(nlpi, problem, values) );
   
   return SCIP_OKAY;
}

/** tries to solve NLP */
SCIP_RETCODE SCIPnlpiSolve(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpisolve)(nlpi, problem) );
   
   return SCIP_OKAY;
}

/** gives solution status, return: Solution Status */
SCIP_NLPSOLSTAT SCIPnlpiGetSolstat(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   return (*nlpi->nlpigetsolstat)(nlpi, problem);
}

/** gives termination reason; return: Termination Status */
SCIP_NLPTERMSTAT SCIPnlpiGetTermstat(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   return (*nlpi->nlpigettermstat)(nlpi, problem);
}

/** gives primal solution */
SCIP_RETCODE SCIPnlpiGetSolution(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_Real**           primalvalues        /**< pointer to store primal values */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetsolution)(nlpi, problem, primalvalues) );
   
   return SCIP_OKAY;
}

/** gives solve statistics */
SCIP_RETCODE SCIPnlpiGetStatistics(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetstatistics)(nlpi, problem, statistics) );
   
   return SCIP_OKAY;
}

/** gives required size of a buffer to store a warmstart object */
SCIP_RETCODE SCIPnlpiGetWarmstartSize(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   size_t*               size                /**< pointer to store required size for warmstart buffer */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetwarmstartsize)(nlpi, problem, size) );
   
   return SCIP_OKAY;
}

/** stores warmstart information in buffer */
SCIP_RETCODE SCIPnlpiGetWarmstartMemo(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   void*                 buffer              /**< memory to store warmstart information */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetwarmstartmemo)(nlpi, problem, buffer) );
   
   return SCIP_OKAY;
}

/** sets warmstart information in solver */
SCIP_RETCODE SCIPnlpiSetWarmstartMemo(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   void*                 buffer              /**< warmstart information */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetwarmstartmemo)(nlpi, problem, buffer) );
   
   return SCIP_OKAY;
}

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP */
SCIP_RETCODE SCIPnlpiGetIntPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   assert(ival    != NULL);

   SCIP_CALL( (*nlpi->nlpigetintpar)(nlpi, problem, type, ival) );
   
   return SCIP_OKAY;
}

/** sets integer parameter of NLP */
SCIP_RETCODE SCIPnlpiSetIntPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetintpar)(nlpi, problem, type, ival) );
   
   return SCIP_OKAY;
}

/** gets floating point parameter of NLP
 * if problem is NULL and type == SCIP_NLPPAR_INFINITY, then gets solver-wide value for infinity */
SCIP_RETCODE SCIPnlpiGetRealPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure, can be NULL only if type == SCIP_NLPPAR_INFINITY */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   assert(dval    != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetrealpar)(nlpi, problem, type, dval) );
   
   return SCIP_OKAY;
}

/** sets floating point parameter of NLP
 * if problem is NULL and type == SCIP_NLPPAR_INFINITY, then sets solver-wide value for infinity */
SCIP_RETCODE SCIPnlpiSetRealPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure, can be NULL only if type == SCIP_NLPPAR_INFINITY */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   assert(nlpi    != NULL);
   assert(problem != NULL || type == SCIP_NLPPAR_INFINITY);
   
   SCIP_CALL( (*nlpi->nlpisetrealpar)(nlpi, problem, type, dval) );
   
   return SCIP_OKAY;
}

/** gets string parameter of NLP */
SCIP_RETCODE SCIPnlpiGetStringPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char**          sval                /**< pointer to store the parameter value, the user must not modify the string */
)
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   assert(sval    != NULL);
   
   SCIP_CALL( (*nlpi->nlpigetstringpar)(nlpi, problem, type, sval) );
   
   return SCIP_OKAY;
}

/** sets string parameter of NLP */
SCIP_RETCODE SCIPnlpiSetStringPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char*           sval                /**< parameter value */
)
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   SCIP_CALL( (*nlpi->nlpisetstringpar)(nlpi, problem, type, sval) );
   
   return SCIP_OKAY;
}

/** gets data of an NLPI */
SCIP_NLPIDATA* SCIPnlpiGetData(
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

/** gets NLP solver descriptions */
const char* SCIPnlpiGetDesc(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
)
{
   assert(nlpi != NULL);

   return nlpi->description;
}

/** gets NLP solver priority */
int SCIPnlpiGetPriority(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
)
{
   assert(nlpi != NULL);

   return nlpi->priority;
}

/** sets NLP solver priority */
void SCIPnlpiSetPriority(
   SCIP_NLPI*            nlpi,               /**< NLP interface structure */
   int                   priority            /**< new priority of NLPI */
)
{
   assert(nlpi != NULL);

   nlpi->priority = priority;
}

/** creates an NLP statistics structure */
SCIP_RETCODE SCIPnlpStatisticsCreate(
   SCIP_NLPSTATISTICS**  statistics          /**< pointer where to store NLP statistics structure */
)
{
   assert(statistics != NULL);
   
   if( BMSallocMemory(statistics) == NULL )
      return SCIP_NOMEMORY;
   assert(*statistics != NULL);
   
   (*statistics)->niterations = -1;
   (*statistics)->totaltime = -1.0;
   
   return SCIP_OKAY;
}

/** frees an NLP statistics structure */
void SCIPnlpStatisticsFree(
   SCIP_NLPSTATISTICS**  statistics          /**< pointer where to store NLP statistics structure */
)
{
   assert(statistics != NULL);
   
   BMSfreeMemory(statistics);
   
   assert(*statistics == NULL);
}

/** gets the number of iterations from an NLP statistics structure */
int SCIPnlpStatisticsGetNIterations(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
)
{
   assert(statistics != NULL);
   
   return statistics->niterations;
}

/** gets the total time from an NLP statistics structure */
SCIP_Real SCIPnlpStatisticsGetTotalTime(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
)
{
   assert(statistics != NULL);
   
   return statistics->totaltime;
}

/** sets the number of iterations in an NLP statistics structure */
void SCIPnlpStatisticsSetNIterations(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   int                   niterations         /**< number of iterations to store */
)
{
   assert(statistics != NULL);
   statistics->niterations = niterations;
}

/** sets the total time in an NLP statistics structure */
void SCIPnlpStatisticsSetTotalTime(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   SCIP_Real             totaltime           /**< solution time to store */
)
{
   assert(statistics != NULL);
   statistics->totaltime = totaltime;
}
