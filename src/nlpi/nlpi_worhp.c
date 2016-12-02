/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    nlpi_worhp.c
 * @ingroup NLPIS
 * @brief   WORHP NLP interface
 * @author  Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "nlpi/nlpi_worhp.h"
#include "nlpi/nlpi.h"
#include "nlpi/nlpioracle.h"
#include "nlpi/exprinterpret.h"
#include "scip/interrupt.h"
#include "scip/pub_misc.h"

#include <stdio.h>
#include <stdlib.h>

#include "worhp/worhp.h"

#define NLPI_NAME              "worhp"                      /* short concise name of solver */
#define NLPI_DESC              "WORHP interface"            /* description of solver */
#define NLPI_PRIORITY          -1                           /* priority of NLP solver */

/*
 * Data structures
 */

/* TODO: fill in the necessary NLP solver interface data */

struct SCIP_NlpiData
{
   BMS_BLKMEM*                 blkmem;       /**< block memory */
   SCIP_MESSAGEHDLR*           messagehdlr;  /**< message handler */
   SCIP_Real                   infinity;     /**< initial value for infinity */
};

/* TODO: fill in the necessary NLP problem instance data */

struct SCIP_NlpiProblem
{
   SCIP_NLPIORACLE*            oracle;       /**< Oracle-helper to store and evaluate NLP */
   BMS_BLKMEM*                 blkmem;       /**< block memory */

   SCIP_Real*                  lastsol;      /**< solution from last run, if available */
   SCIP_NLPTERMSTAT            lasttermstat; /**< termination status from last run */
   SCIP_NLPSOLSTAT             lastsolstat;  /**< solution status from last run */
   SCIP_Real                   lastsolinfeas;/**< infeasibility (constraint violation) of solution stored in lastsol */
   SCIP_Real                   lasttime;     /**< time spend in last run */
   int                         lastsolsize;  /**< size of lastsol array */
   int                         lastniter;    /**< number of iterations in last run */

   /* Worhp data structures (not used yet) */
   OptVar*                     opt;          /**< Worhp variables  */
   Workspace*                  wsp;          /**< Worhp working space */
   Params*                     par;          /**< Worhp parameters */
   Control*                    cnt;          /**< Worhp control */

};


/*
 * Local methods
 */

/** clears the last solution information */
static
void invalidateSolution(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(problem != NULL);
   assert(problem->blkmem != NULL);

   BMSfreeBlockMemoryArrayNull(problem->blkmem, &(problem->lastsol), problem->lastsolsize);

   problem->lastsol = NULL;
   problem->lastsolsize = 0;
   problem->lastsolinfeas = SCIP_INVALID;
   problem->lastsolstat = SCIP_NLPSOLSTAT_UNKNOWN;
   problem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
}

/** evaluate last Worhp run */
static
SCIP_RETCODE evaluateWorhpRun(
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler (might be NULL) */
   )
{
   assert(problem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->par != NULL);
   assert(problem->cnt != NULL);

  /* initialization error */
  if (problem->cnt->status == InitError)
  {
     SCIPmessagePrintWarning(messagehdlr, "Worhp failed because of initialization error!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_MEMERR;
  }
  /* data error */
  else if (problem->cnt->status == DataError)
  {
     SCIPmessagePrintWarning(messagehdlr, "Worhp failed because of data error!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
  }
  /* license error */
  else if (problem->cnt->status == LicenseError)
  {
     SCIPmessagePrintWarning(messagehdlr, "Worhp failed because of license error!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_LICERR;
  }

  /* evaluation errors */
  else if (problem->cnt->status == FunctionErrorF)
  {
     SCIPmessagePrintWarning(messagehdlr, "Worhp failed because of an invalid number in objective function evaluation!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_EVALERR;
  }
  else if (problem->cnt->status == FunctionErrorG)
  {
     SCIPmessagePrintWarning(messagehdlr, "Worhp failed because of an invalid number in constraint evaluation!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_EVALERR;
  }
  else if (problem->cnt->status == FunctionErrorDF)
  {
     SCIPmessagePrintWarning(messagehdlr, "Worhp failed because of an invalid number in gradient of objective function evaluation!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_EVALERR;
  }
  else if (problem->cnt->status == FunctionErrorDG)
  {
     SCIPmessagePrintWarning(messagehdlr, "Worhp failed because of an invalid number in jacobian of constraint evaluation!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_EVALERR;
  }
  else if (problem->cnt->status == FunctionErrorHM)
  {
     SCIPmessagePrintWarning(messagehdlr, "Worhp failed because of an invalid number in hessian evaluation!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_EVALERR;
  }

  /*
   * catch unsuccessful states
   */

  /* numerical errors during solution of NLP */
  else if( (problem->cnt->status == QPerror) || (problem->cnt->status == MinimumStepsize) || (problem->cnt->status == TooBig) )
  {
     SCIPmessagePrintWarning(messagehdlr, "Worhp failed because of an numerical error during optimization!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_NUMERR;
  }
  /* maximal number of calls or iteration */
  else if( (problem->cnt->status == MaxCalls) || (problem->cnt->status == MaxIter) )
  {
     SCIPmessagePrintWarning(messagehdlr, "Worhp failed because maximal number of calls or iterations is reached!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_ITLIM;
  }
  /* time limit reached */
  else if( problem->cnt->status == Timeout )
  {
     SCIPmessagePrintWarning(messagehdlr, "Worhp failed because time limit is reached!\n");
     invalidateSolution(problem);
     problem->lastsolstat = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_TILIM;
  }
  /* infeasible stationary point found */
  else if( problem->cnt->status == LocalInfeas )
  {
     SCIPmessagePrintWarning(messagehdlr, "Worhp failed because of convergence against infeasible stationary point!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }

  /*
   * catch successful states
   */

  /* everything went fine */
  else if( problem->cnt->status == OptimalSolution )
  {
     SCIPdebugMessage("Worhp terminated successfully at a local optimum!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_LOCOPT;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* feasible and no further progress */
  else if( problem->cnt->status == LowPassFilterOptimal )
  {
     SCIPdebugMessage("Worhp terminated at feasible solution without further progress!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* feasible and in feasibility mode, i.e., optimality not required */
  else if( problem->cnt->status == FeasibleSolution )
  {
     SCIPdebugMessage("Worhp terminated at feasible solution, optimality was not required!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* acceptable solution found, but stopped due to limit or error */
  else if( problem->cnt->status == AcceptableSolution )
  {
     SCIPdebugMessage("Worhp terminated at acceptable solution due to limit or error!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* previously acceptable solution was found, but stopped due to limit or error */
  else if( problem->cnt->status == AcceptablePrevious )
  {
     SCIPdebugMessage("Worhp previously found acceptable solution but terminated due to limit or error!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* acceptable solution found, and no further progress */
  else if( problem->cnt->status == LowPassFilterAcceptable )
  {
     SCIPdebugMessage("Worhp found acceptable solution but terminated due no further progress!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* acceptable solution found, but search direction is small or zero */
  else if( (problem->cnt->status == SearchDirectionZero) || (problem->cnt->status == SearchDirectionSmall) )
  {
     SCIPdebugMessage("Worhp found acceptable solution but search direction is small or zero!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* acceptable solution found, but not optimal */
  else if( (problem->cnt->status == FritzJohn) || (problem->cnt->status == NotDiffable) || (problem->cnt->status == Unbounded) )
  {
     SCIPdebugMessage("Worhp found acceptable solution but terminated perhaps due to nondifferentiability, unboundedness or at Fritz John point!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }

  /*
   * catch other states
   */
  else
  {
     SCIPerrorMessage("Worhp returned with unknown solution status %d\n", problem->cnt->status);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
     return SCIP_ERROR;
  }

  /* store solution */
  if( problem->lastsol == NULL )
  {
     SCIP_ALLOC( BMSduplicateBlockMemoryArray(problem->blkmem, &problem->lastsol, problem->opt->X, problem->opt->n) );
     problem->lastsolsize = problem->opt->n;
  }
  else
  {
     BMScopyMemoryArray(problem->lastsol, problem->opt->X, problem->opt->n);
  }

  return SCIP_OKAY;
}


/** evaluates objective function and store the result in the corresponding WORHP data fields */
static
SCIP_RETCODE userF(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   SCIP_Real objval;

   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->blkmem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n = SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m = SCIPnlpiOracleGetNConstraints(problem->oracle));

   SCIP_CALL( SCIPnlpiOracleEvalObjectiveValue(problem->oracle, problem->opt->X, &objval) );
   problem->opt->F *= problem->wsp->ScaleObj * objval;

   return SCIP_OKAY;
}

/** evaluates objective function and store the result in the corresponding WORHP data fields */
static
SCIP_RETCODE userG(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->blkmem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n = SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m = SCIPnlpiOracleGetNConstraints(problem->oracle));

   SCIP_CALL( SCIPnlpiOracleEvalConstraintValues(problem->oracle, problem->opt->X, problem->opt->G) );

   return SCIP_OKAY;
}

/** computes objective gradient and store the result in the corresponding WORHP data fields */
static
SCIP_RETCODE userDF(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   SCIP_Real objval;

   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->blkmem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n = SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m = SCIPnlpiOracleGetNConstraints(problem->oracle));

   /* TODO this needs to be changed if we store the gradient of the objective function in a sparse format */
   SCIP_CALL( SCIPnlpiOracleEvalObjectiveGradient(problem->oracle, problem->opt->X, TRUE, &objval,
         problem->wsp->DF.val) );

   /* scale gradient if necessary */
   if( problem->wsp->ScaleObj != 1.0 )
   {
      int i;
      for( i = 0; i < problem->opt->n; ++i )
         problem->wsp->DF.val[i] *= problem->wsp->ScaleObj;
   }

   return SCIP_OKAY;
}

/** computes jacobian matrix and store the result in the corresponding WORHP data fields */
static
SCIP_RETCODE userDG(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->blkmem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n = SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m = SCIPnlpiOracleGetNConstraints(problem->oracle));

   SCIP_CALL( SCIPnlpiOracleEvalJacobian(problem->oracle, problem->opt->X, TRUE, NULL, problem->wsp->DG.val) );

   return SCIP_OKAY;
}

/** computes hessian matrix and store the result in the corresponding WORHP data fields */
static
SCIP_RETCODE userHM(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->blkmem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n = SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m = SCIPnlpiOracleGetNConstraints(problem->oracle));

   return SCIP_OKAY;
}


/*
 * Callback methods of NLP solver interface
 */

/* TODO: Implement all necessary NLP interface methods. The methods with an #if 0 ... #else #define ... are optional
 * (currently, all methods are required) */

/** copy method of NLP interface (called when SCIP copies plugins)
 *
 * input:
 *  - blkmem block memory in target SCIP
 *  - sourcenlpi the NLP interface to copy
 *  - targetnlpi buffer to store pointer to copy of NLP interface
 */
static
SCIP_DECL_NLPICOPY( nlpiCopyWorhp )
{
   SCIP_NLPIDATA* sourcedata;
   SCIP_NLPIDATA* targetdata;

   assert(sourcenlpi != NULL);
   assert(targetnlpi != NULL);

   sourcedata = SCIPnlpiGetData(sourcenlpi);
   assert(sourcedata != NULL);

   SCIP_CALL( SCIPcreateNlpSolverWorhp(blkmem, targetnlpi) );
   assert(*targetnlpi != NULL);

   SCIP_CALL( SCIPnlpiSetRealPar(*targetnlpi, NULL, SCIP_NLPPAR_INFINITY, sourcedata->infinity) );
   SCIP_CALL( SCIPnlpiSetMessageHdlr(*targetnlpi, sourcedata->messagehdlr) );

   targetdata = SCIPnlpiGetData(*targetnlpi);
   assert(targetdata != NULL);

   targetdata->blkmem = sourcedata->blkmem;
   targetdata->messagehdlr = sourcedata->messagehdlr;
   targetdata->infinity = sourcedata->infinity;

   return SCIP_OKAY;
}  /*lint !e715*/

/** destructor of NLP interface to free nlpi data
 *
 * input:
 *  - nlpi datastructure for solver interface
 */
static
SCIP_DECL_NLPIFREE( nlpiFreeWorhp )
{
   SCIP_NLPIDATA* data;

   assert(nlpi != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);
   assert(data->blkmem != NULL);

   BMSfreeMemory(&data);
   data = NULL;

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets pointer for NLP solver
 *
 *  to do dirty stuff
 *
 * input:
 *  - nlpi datastructure for solver interface
 *
 * return: void pointer to solver
 */
static
SCIP_DECL_NLPIGETSOLVERPOINTER(nlpiGetSolverPointerWorhp)
{
   assert(nlpi != NULL);

   return NULL;
}  /*lint !e715*/

/** creates a problem instance
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem pointer to store the problem data
 *  - name name of problem, can be NULL
 */
static
SCIP_DECL_NLPICREATEPROBLEM(nlpiCreateProblemWorhp)
{
   SCIP_NLPIDATA* data;

   assert(nlpi    != NULL);
   assert(problem != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(data->blkmem, problem) );
   if( *problem == NULL )
      return SCIP_NOMEMORY;

   /* initialize data for oracle */
   SCIP_CALL( SCIPnlpiOracleCreate(data->blkmem, &(*problem)->oracle) );
   SCIP_CALL( SCIPnlpiOracleSetInfinity((*problem)->oracle, data->infinity) );
   SCIP_CALL( SCIPnlpiOracleSetProblemName((*problem)->oracle, name) );

   /* allocate memory for WORHP data */
   SCIP_ALLOC( BMSallocBlockMemory(data->blkmem, &(*problem)->opt) );
   SCIP_ALLOC( BMSallocBlockMemory(data->blkmem, &(*problem)->wsp) );
   SCIP_ALLOC( BMSallocBlockMemory(data->blkmem, &(*problem)->par) );
   SCIP_ALLOC( BMSallocBlockMemory(data->blkmem, &(*problem)->cnt) );
   /* WorhpPreInit((*problem)->opt, (*problem)->wsp, (*problem)->par, (*problem)->cnt); */

   return SCIP_OKAY;
}  /*lint !e715*/

/** free a problem instance
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem pointer where problem data is stored
 */
static
SCIP_DECL_NLPIFREEPROBLEM(nlpiFreeProblemWorhp)
{
   SCIP_NLPIDATA* data;

   assert(nlpi     != NULL);
   assert(problem  != NULL);
   assert(*problem != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   if( (*problem)->opt != NULL )
   {
      assert((*problem)->wsp != NULL);
      assert((*problem)->par != NULL);
      assert((*problem)->cnt != NULL);

      /* free WORHP data */
      /* WorhpFree((*problem)->opt, (*problem)->wsp, (*problem)->par, (*problem)->cnt); */
      BMSfreeBlockMemory(data->blkmem, &(*problem)->cnt);
      BMSfreeBlockMemory(data->blkmem, &(*problem)->par);
      BMSfreeBlockMemory(data->blkmem, &(*problem)->wsp);
      BMSfreeBlockMemory(data->blkmem, &(*problem)->opt);
   }

   if( (*problem)->oracle != NULL )
   {
      SCIP_CALL( SCIPnlpiOracleFree(&(*problem)->oracle) );
   }

   BMSfreeBlockMemory(data->blkmem, problem);
   *problem = NULL;

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets pointer to solver-internal problem instance
 *
 *  to do dirty stuff
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 * return: void pointer to problem instance
 */
static
SCIP_DECL_NLPIGETPROBLEMPOINTER(nlpiGetProblemPointerWorhp)
{
   assert(nlpi    != NULL);
   assert(problem != NULL);

   return NULL;
}  /*lint !e715*/

/** add variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nvars number of variables
 *  - lbs lower bounds of variables, can be NULL if -infinity
 *  - ubs upper bounds of variables, can be NULL if +infinity
 *  - varnames names of variables, can be NULL
 */
static
SCIP_DECL_NLPIADDVARS( nlpiAddVarsWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleAddVars(problem->oracle, nvars, lbs, ubs, varnames) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/


/** add constraints
 * quadratic coefficiens: row oriented matrix for each constraint
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - ncons number of added constraints
 *  - lhss left hand sides of constraints
 *  - rhss right hand sides of constraints
 *  - nlininds number of linear coefficients for each constraint
 *    may be NULL in case of no linear part
 *  - lininds indices of variables for linear coefficients for each constraint
 *    may be NULL in case of no linear part
 *  - linvals values of linear coefficient for each constraint
 *    may be NULL in case of no linear part
 *  - nquadrows number of columns in matrix of quadratic part for each constraint
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadrowidxs indices of variables for which a quadratic part is specified
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadoffsets start index of each rows quadratic coefficients in quadinds[.] and quadvals[.]
 *    indices are given w.r.t. quadrowidxs., i.e., quadoffsets[.][i] gives the start index of row quadrowidxs[.][i] in quadvals[.]
 *    quadoffsets[.][nquadrows[.]] gives length of quadinds[.] and quadvals[.]
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadinds column indices w.r.t. quadrowidxs, i.e., quadrowidxs[quadinds[.][i]] gives the index of the variable corresponding
 *    to entry i, entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadvals coefficient values
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - exprvaridxs indices of variables in expression tree, maps variable indices in expression tree to indices in nlp
 *    entry of array may be NULL in case of no expression tree
 *    may be NULL in case of no expression tree in any constraint
 *  - exprtrees expression tree for nonquadratic part of constraints
 *    entry of array may be NULL in case of no nonquadratic part
 *    may be NULL in case of no nonquadratic part in any constraint
 *  - names of constraints, may be NULL or entries may be NULL
 */
static
SCIP_DECL_NLPIADDCONSTRAINTS( nlpiAddConstraintsWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleAddConstraints(problem->oracle,
         ncons, lhss, rhss,
         nlininds, lininds, linvals,
         nquadelems, quadelems,
         exprvaridxs, exprtrees, names) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets or overwrites objective, a minimization problem is expected
 *  May change sparsity pattern.
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nlins number of linear variables
 *  - lininds variable indices
 *    may be NULL in case of no linear part
 *  - linvals coefficient values
 *    may be NULL in case of no linear part
 *  - nquadcols number of columns in matrix of quadratic part
 *  - quadcols indices of variables for which a quadratic part is specified
 *    may be NULL in case of no quadratic part
 *  - quadoffsets start index of each rows quadratic coefficients in quadinds and quadvals
 *    quadoffsets[.][nquadcols] gives length of quadinds and quadvals
 *    may be NULL in case of no quadratic part
 *  - quadinds column indices
 *    may be NULL in case of no quadratic part
 *  - quadvals coefficient values
 *    may be NULL in case of no quadratic part
 *  - exprvaridxs indices of variables in expression tree, maps variable indices in expression tree to indices in nlp
 *    may be NULL in case of no expression tree
 *  - exprtree expression tree for nonquadratic part of objective function
 *    may be NULL in case of no nonquadratic part
 *  - constant objective value offset
 */
static
SCIP_DECL_NLPISETOBJECTIVE( nlpiSetObjectiveWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleSetObjective(problem->oracle,
         constant, nlins, lininds, linvals,
         nquadelems, quadelems,
         exprvaridxs, exprtree) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change variable bounds
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nvars number of variables to change bounds
 *  - indices indices of variables to change bounds
 *  - lbs new lower bounds
 *  - ubs new upper bounds
 */
static
SCIP_DECL_NLPICHGVARBOUNDS( nlpiChgVarBoundsWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgVarBounds(problem->oracle, nvars, indices, lbs, ubs) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change constraint bounds
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nconss number of constraints to change sides
 *  - indices indices of constraints to change sides
 *  - lhss new left hand sides
 *  - rhss new right hand sides
 */
static
SCIP_DECL_NLPICHGCONSSIDES( nlpiChgConsSidesWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgConsSides(problem->oracle, nconss, indices, lhss, rhss) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** delete a set of variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - dstats deletion status of vars; 1 if var should be deleted, 0 if not
 *
 * output:
 *  - dstats new position of var, -1 if var was deleted
 */
static
SCIP_DECL_NLPIDELVARSET( nlpiDelVarSetWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelVarSet(problem->oracle, dstats) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** delete a set of constraints
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - dstats deletion status of rows; 1 if row should be deleted, 0 if not
 *
 * output:
 *  - dstats new position of row, -1 if row was deleted
 */
static
SCIP_DECL_NLPIDELCONSSET( nlpiDelConstraintSetWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelConsSet(problem->oracle, dstats) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** changes (or adds) linear coefficients in a constraint or objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idx index of constraint or -1 for objective
 *  - nvals number of values in linear constraint to change
 *  - varidxs indices of variables which coefficient to change
 *  - vals new values for coefficients
 */
static
SCIP_DECL_NLPICHGLINEARCOEFS( nlpiChgLinearCoefsWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgLinearCoefs(problem->oracle, idx, nvals, varidxs, vals) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** changes (or adds) coefficients in the quadratic part of a constraint or objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idx index of constraint or -1 for objective
 *  - nentries number of entries in quadratic matrix to change
 *  - rows row indices of entries in quadratic matrix where values should be changed
 *  - cols column indices of entries in quadratic matrix where values should be changed
 *  - values new values for entries in quadratic matrix
 */
static
SCIP_DECL_NLPICHGQUADCOEFS( nlpiChgQuadraticCoefsWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgQuadCoefs(problem->oracle, idx, nquadelems, quadelems) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** replaces the expression tree of a constraint or objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idxcons index of constraint or -1 for objective
 *  - exprvaridxs indices of variables in expression tree, maps variable indices in expression tree to indices in nlp, or NULL
 *  - exprtree new expression tree for constraint or objective, or NULL to only remove previous tree
 */
static
SCIP_DECL_NLPICHGEXPRTREE( nlpiChgExprtreeWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgExprtree(problem->oracle, idxcons, exprvaridxs, exprtree) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change one coefficient in the nonlinear part
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idxcons index of constraint or -1 for objective
 *  - idxparam index of parameter
 *  - value new value for nonlinear parameter
 *
 * return: Error if parameter does not exist
 */
static
SCIP_DECL_NLPICHGNONLINCOEF( nlpiChgNonlinCoefWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgExprParam(problem->oracle, idxcons, idxparam, value) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change the constant offset in the objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - objconstant new value for objective constant
 */
static
SCIP_DECL_NLPICHGOBJCONSTANT( nlpiChgObjConstantWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgObjConstant(problem->oracle, objconstant) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets initial guess for primal variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues initial primal values for variables, or NULL to clear previous values
 *  - consdualvalues initial dual values for constraints, or NULL to clear previous values
 *  - varlbdualvalues  initial dual values for variable lower bounds, or NULL to clear previous values
 *  - varubdualvalues  initial dual values for variable upper bounds, or NULL to clear previous values
 */
static
SCIP_DECL_NLPISETINITIALGUESS( nlpiSetInitialGuessWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** tries to solve NLP
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 */
static
SCIP_DECL_NLPISOLVE( nlpiSolveWorhp )
{
   SCIP_NLPIDATA* nlpidata;
   const SCIP_Real* lbs;
   const SCIP_Real* ubs;
   const int* offset;
   const int* cols;
   int i;

   OptVar* opt = problem->opt;
   Workspace* wsp = problem->wsp;
   Params* par = problem->par;
   Control* cnt = problem->cnt;

   assert(nlpi != NULL);
   assert(problem != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   /* properly zeros everything */
   WorhpPreInit(opt, wsp, par, cnt);
   par->Infty = SCIP_DEFAULT_INFINITY;

   /* set problem dimensions */
   opt->n = SCIPnlpiOracleGetNVars(problem->oracle);
   opt->m = SCIPnlpiOracleGetNConstraints(problem->oracle);
   SCIPdebugMessage("nvars %d nconss %d\n", opt->n, opt->m);

   /* assume that objective function is dense TODO use sparse representation */
   wsp->DF.nnz = opt->n;

   /* get number of non-zero entries in Jacobian */
   SCIP_CALL( SCIPnlpiOracleGetJacobianSparsity(problem->oracle, &offset, NULL) );
   wsp->DG.nnz = offset[opt->m];
   SCIPdebugMessage("nnonz jacobian %d\n", wsp->DG.nnz);

   /* get number of non-zero entries in Hessian */
   SCIP_CALL( SCIPnlpiOracleGetHessianLagSparsity(problem->oracle, &offset, NULL) );
   wsp->HM.nnz = offset[opt->n];
   SCIPdebugMessage("nnonz hessian %d\n", wsp->HM.nnz);

   /* initialize data in WORHP */
   WorhpInit(opt, wsp, par, cnt);
   if (cnt->status != FirstCall)
   {
      SCIPerrorMessage("Initialisation failed.\n");
      return SCIP_ERROR;
   }

   /* set variable bounds */
   lbs = SCIPnlpiOracleGetVarLbs(problem->oracle);
   ubs = SCIPnlpiOracleGetVarUbs(problem->oracle);
   for( i = 0; i < opt->n; ++i )
   {
      opt->XL[i] = lbs[i];
      opt->XU[i] = ubs[i];
      SCIPdebugMessage("bounds %d [%g,%g]\n", i, opt->XL[i], opt->XU[i]);
   }

   /* set constraint sides */
   for( i = 0; i < opt->m; ++i )
   {
      opt->GL[i] = SCIPnlpiOracleGetConstraintLhs(problem->oracle, i);
      opt->GU[i] = SCIPnlpiOracleGetConstraintRhs(problem->oracle, i);
      SCIPdebugMessage("sides %d [%g,%g]\n", i, opt->GL[i], opt->GU[i]);
   }

   /* set column indices of objective function; note that indices go from 1 to n */
   if( wsp->DF.NeedStructure )
   {
      SCIPdebugMessage("column indices of objective function:");
      for( i = 0; i < opt->n; ++i )
      {
         wsp->DF.row[i] = i + 1;
         SCIPdebugMessage(" %d", wsp->DF.row[i]);
      }
      SCIPdebugMessage("\n");
   }

   /* set column and row indices of non-zero entries in Jacobian matrix */
   if( wsp->DG.NeedStructure )
   {
      int nnonz;
      int j;

      SCIP_CALL( SCIPnlpiOracleGetJacobianSparsity(problem->oracle, &offset, &cols) );
      assert(offset[opt->m] == wsp->DG.nnz);
      SCIPdebugMessage("column and row indices of jacobian:\n");

      nnonz = 0;
      j = offset[0];
      for( i = 0; i < opt->m; ++i )
      {
         for( ; j < offset[i+1]; ++j )
         {
            /* note that column and row indices are shifted by one */
            wsp->DG.row[nnonz] = i + 1;
            wsp->DG.col[nnonz] = cols[j] + 1;
            SCIPdebugMessage("  entry %d: (row,col) = (%d,%d)\n", nnonz, wsp->DG.row[nnonz], wsp->DG.col[nnonz]);
            ++nnonz;
         }
      }
      assert(nnonz == wsp->DG.nnz);

      /* sort arrays such that the entries are in column-major order */
      SCIPsortIntInt(wsp->DG.col, wsp->DG.row, nnonz);
   }

   /* set column and row indices of non-zero entries in hessian matrix */
   if( wsp->HM.NeedStructure )
   {
      int nnonz;
      int j;

      SCIP_CALL( SCIPnlpiOracleGetHessianLagSparsity(problem->oracle, &offset, &cols) );
      assert(offset[opt->n] == wsp->HM.nnz);
      SCIPdebugMessage("column and row indices of hessian:\n");

      nnonz = 0;
      j = offset[0];
      for( i = 0; i < opt->n; ++i )
      {
         for( ; j < offset[i+1]; ++j )
         {
            /* note that column and row indices are shifted by one */
            wsp->HM.row[nnonz] = i + 1;
            wsp->HM.col[nnonz] = cols[j] + 1;
            SCIPdebugMessage("  entry %d: (row,col) = (%d,%d)\n", nnonz, wsp->HM.row[nnonz], wsp->HM.col[nnonz]);
            ++nnonz;
         }
      }
      assert(nnonz == wsp->HM.nnz);

      /* sort arrays such that the entries are in column-major order */
      SCIPsortIntInt(wsp->HM.col, wsp->HM.row, nnonz);
   }

   /* TODO set a start point */

   /*
    * WORHP Reverse Communication loop.
    * In every iteration poll GetUserAction for the requested action, i.e. one
    * of {callWorhp, iterOutput, evalF, evalG, evalDF, evalDG, evalHM, fidif}.
    *
    * Make sure to reset the requested user action afterwards by calling
    * DoneUserAction, except for 'callWorhp' and 'fidif'.
    */
   while( cnt->status < TerminateSuccess && cnt->status > TerminateError )
   {
      /*
       * WORHP's main routine.
       * Do not manually reset callWorhp, this is only done by the FD routines.
       */
      if( GetUserAction(cnt, callWorhp) )
      {
         Worhp(opt, wsp, par, cnt);
         /* No DoneUserAction! */
      }

      /*
       * Show iteration output.
       * The call to IterationOutput() may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, iterOutput) )
      {
         IterationOutput(opt, wsp, par, cnt);
         DoneUserAction(cnt, iterOutput);
      }

      /*
       * Evaluate the objective function.
       * The call to UserF may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalF) )
      {
         SCIP_CALL( userF(problem) );
         DoneUserAction(cnt, evalF);
      }

      /*
       * Evaluate the constraints.
       * The call to UserG may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalG) )
      {
         SCIP_CALL( userG(problem) );
         DoneUserAction(cnt, evalG);
      }

      /*
       * Evaluate the gradient of the objective function.
       * The call to UserDF may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalDF) )
      {
         SCIP_CALL( userDF(problem) );
         DoneUserAction(cnt, evalDF);
      }

      /*
       * Evaluate the Jacobian of the constraints.
       * The call to UserDG may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalDG) )
      {
         SCIP_CALL( userDG(problem) );
         DoneUserAction(cnt, evalDG);
      }

      /*
       * Evaluate the Hessian matrix of the Lagrange function (L = f + mu*g)
       * The call to UserHM may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalHM) )
      {
         SCIP_CALL( userHM(problem) );
         DoneUserAction(cnt, evalHM);
      }

      /*
       * Use finite differences with RC to determine derivatives
       * Do not reset fidif, this is done by the FD routine.
       */
      if( GetUserAction(cnt, fidif) )
      {
         WorhpFidif(opt, wsp, par, cnt);
         /* No DoneUserAction! */
      }
   }

   /* interpret Worhp result */
   SCIP_CALL( evaluateWorhpRun(problem, nlpidata->messagehdlr) );

   /* free memory allocated in WORHP */
   WorhpFree(opt, wsp, par, cnt);

   return SCIP_OKAY;
}  /*lint !e715*/

/** gives solution status
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 * return: Solution Status
 */
static
SCIP_DECL_NLPIGETSOLSTAT( nlpiGetSolstatWorhp )
{
   /* TODO */

   return SCIP_NLPSOLSTAT_UNKNOWN;  /*lint !e527*/
}  /*lint !e715*/

/** gives termination reason
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 * return: Termination Status
 */
static
SCIP_DECL_NLPIGETTERMSTAT( nlpiGetTermstatWorhp )
{
   /* TODO */

   return SCIP_NLPTERMSTAT_OTHER;  /*lint !e527*/
}  /*lint !e715*/

/** gives primal and dual solution values
 *
 * solver can return NULL in dual values if not available
 * but if solver provides dual values for one side of variable bounds, then it must also provide those for the other side
 *
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues buffer to store pointer to array to primal values, or NULL if not needed
 *  - consdualvalues buffer to store pointer to array to dual values of constraints, or NULL if not needed
 *  - varlbdualvalues buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed
 *  - varubdualvalues buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed
 */
static
SCIP_DECL_NLPIGETSOLUTION( nlpiGetSolutionWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives solve statistics
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - statistics pointer to store statistics
 *
 * output:
 *  - statistics solve statistics
 */
static
SCIP_DECL_NLPIGETSTATISTICS( nlpiGetStatisticsWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives required size of a buffer to store a warmstart object
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - size pointer to store required size for warmstart buffer
 *
 * output:
 *  - size required size for warmstart buffer
 */
static
SCIP_DECL_NLPIGETWARMSTARTSIZE( nlpiGetWarmstartSizeWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** stores warmstart information in buffer
 *
 * required size of buffer should have been obtained by SCIPnlpiGetWarmstartSize before
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - buffer memory to store warmstart information
 *
 * output:
 *  - buffer warmstart information in solver specific data structure
 */
static
SCIP_DECL_NLPIGETWARMSTARTMEMO( nlpiGetWarmstartMemoWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets warmstart information in solver
 *
 * write warmstart to buffer
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - buffer warmstart information
 */
static
SCIP_DECL_NLPISETWARMSTARTMEMO( nlpiSetWarmstartMemoWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gets integer parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival pointer to store the parameter value
 *
 * output:
 *  - ival parameter value
 */
static
SCIP_DECL_NLPIGETINTPAR( nlpiGetIntParWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets integer parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival parameter value
 */
static
SCIP_DECL_NLPISETINTPAR( nlpiSetIntParWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gets floating point parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval pointer to store the parameter value
 *
 * output:
 *  - dval parameter value
 */
static
SCIP_DECL_NLPIGETREALPAR( nlpiGetRealParWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets floating point parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval parameter value
 */
static
SCIP_DECL_NLPISETREALPAR( nlpiSetRealParWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gets string parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - sval pointer to store the string value, the user must not modify the string
 *
 * output:
 *  - sval parameter value
 */
static
SCIP_DECL_NLPIGETSTRINGPAR( nlpiGetStringParWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets string parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - sval parameter value
 */
static
SCIP_DECL_NLPISETSTRINGPAR( nlpiSetStringParWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets message handler for message output
 *
 * input:
 *  - nlpi NLP interface structure
 *  - messagehdlr SCIP message handler, or NULL to suppress all output
 */
static
SCIP_DECL_NLPISETMESSAGEHDLR( nlpiSetMessageHdlrWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/*
 * NLP solver interface specific interface methods
 */

/** create solver interface for Worhp solver */
SCIP_RETCODE SCIPcreateNlpSolverWorhp(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi                /**< pointer to buffer for nlpi address */
   )
{
   SCIP_NLPIDATA* nlpidata;

   assert(blkmem != NULL);
   assert(nlpi   != NULL);

   /* create WORHP solver interface data */
   BMSallocMemory(&nlpidata);
   BMSclearMemory(nlpidata);

   nlpidata->blkmem = blkmem;
   nlpidata->infinity = SCIP_DEFAULT_INFINITY;

   /* checks the version of the library and header files */
   CHECK_WORHP_VERSION

   /* create solver interface */
   SCIP_CALL( SCIPnlpiCreate(nlpi,
         NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
         nlpiCopyWorhp, nlpiFreeWorhp, nlpiGetSolverPointerWorhp,
         nlpiCreateProblemWorhp, nlpiFreeProblemWorhp, nlpiGetProblemPointerWorhp,
         nlpiAddVarsWorhp, nlpiAddConstraintsWorhp, nlpiSetObjectiveWorhp,
         nlpiChgVarBoundsWorhp, nlpiChgConsSidesWorhp, nlpiDelVarSetWorhp, nlpiDelConstraintSetWorhp,
         nlpiChgLinearCoefsWorhp, nlpiChgQuadraticCoefsWorhp, nlpiChgExprtreeWorhp, nlpiChgNonlinCoefWorhp,
         nlpiChgObjConstantWorhp, nlpiSetInitialGuessWorhp, nlpiSolveWorhp, nlpiGetSolstatWorhp, nlpiGetTermstatWorhp,
         nlpiGetSolutionWorhp, nlpiGetStatisticsWorhp,
         nlpiGetWarmstartSizeWorhp, nlpiGetWarmstartMemoWorhp, nlpiSetWarmstartMemoWorhp,
         nlpiGetIntParWorhp, nlpiSetIntParWorhp, nlpiGetRealParWorhp, nlpiSetRealParWorhp, nlpiGetStringParWorhp, nlpiSetStringParWorhp,
         nlpiSetMessageHdlrWorhp,
         nlpidata) );

   return SCIP_OKAY;
}

/** gets string that identifies Worhp (version number) */
const char* SCIPgetSolverNameWorhp(void)
{
   /* TODO maybe you can add this macro to worhp_version.h */
   /* return "WORHP " WORHP_VERSION; */
   return "WORHP 1.9.1";
}

/** gets string that describes Worhp (version number) */
extern
const char* SCIPgetSolverDescWorhp(void)
{
   /* TODO */
   return "Sequential Quadratic Programming developed at Research Institute Steinbeis (www.worhp.de)";
}

/** returns whether Worhp is available, i.e., whether it has been linked in */
extern
SCIP_Bool SCIPisWorhpAvailableWorhp(void)
{
   return TRUE;
}
