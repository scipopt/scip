/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nlpi_xxx.c,v 1.8 2010/09/15 14:57:34 bzfviger Exp $"

/**@file    nlpi_xxx.c
 * @ingroup NLPIS
 * @brief   XXX NLP interface
 * @author  you
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "nlpi/nlpi_xxx.h"
#include "nlpi/nlpi.h"

#define NLPI_NAME              "xxx"                       /* short concise name of solver */
#define NLPI_DESC              "solver interface template" /* description of solver */
#define NLPI_PRIORITY          0                           /* priority of NLP solver */

/*
 * Data structures
 */

/* TODO: fill in the necessary NLP solver interface data */

struct SCIP_NlpiData
{
};

/* TODO: fill in the necessary NLP problem instance data */

struct SCIP_NlpiProblem
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of NLP solver interface
 */

/* TODO: Implement all necessary NLP interface methods. The methods with an #if 0 ... #else #define ... are optional
 * (currently, all methods are required) */

/** copy method of NLP interface (called when SCIP copies plugins)
 *
 * input:
 *  - sourcenlpi the NLP interface to copy
 *  - targetnlpi buffer to store pointer to copy of NLP interface
 */
static
SCIP_DECL_NLPICOPY( nlpiCopyXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;
}

/** destructor of NLP interface to free nlpi data
 * 
 * input:
 *  - nlpi datastructure for solver interface
 */
static
SCIP_DECL_NLPIFREE( nlpiFreeXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETSOLVERPOINTER(nlpiGetSolverPointerXxx)
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

/** creates a problem instance
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem pointer to store the problem data
 *  - name name of problem, can be NULL
 */
static
SCIP_DECL_NLPICREATEPROBLEM(nlpiCreateProblemXxx)
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

/** free a problem instance
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem pointer where problem data is stored 
 */
static
SCIP_DECL_NLPIFREEPROBLEM(nlpiFreeProblemXxx)
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETPROBLEMPOINTER(nlpiGetProblemPointerXxx)
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return NULL;
}

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
SCIP_DECL_NLPIADDVARS( nlpiAddVarsXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}


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
SCIP_DECL_NLPIADDCONSTRAINTS( nlpiAddConstraintsXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPISETOBJECTIVE( nlpiSetObjectiveXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPICHGVARBOUNDS( nlpiChgVarBoundsXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPICHGCONSSIDES( nlpiChgConsSidesXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIDELVARSET( nlpiDelVarSetXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIDELCONSSET( nlpiDelConstraintSetXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPICHGLINEARCOEFS( nlpiChgLinearCoefsXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPICHGQUADCOEFS( nlpiChgQuadraticCoefsXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPICHGEXPRTREE( nlpiChgExprtreeXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;
}

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
SCIP_DECL_NLPICHGNONLINCOEF( nlpiChgNonlinCoefXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

/** change the constant offset in the objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - objconstant new value for objective constant
 */
static
SCIP_DECL_NLPICHGOBJCONSTANT( nlpiChgObjConstantXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

/** sets initial guess for primal variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - values initial starting solution, or NULL to clear previous starting solution
 */
static
SCIP_DECL_NLPISETINITIALGUESS( nlpiSetInitialGuessXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

/** tries to solve NLP
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 */
static
SCIP_DECL_NLPISOLVE( nlpiSolveXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

/** gives solution status
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 * 
 * return: Solution Status
 */
static
SCIP_DECL_NLPIGETSOLSTAT( nlpiGetSolstatXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_NLPSOLSTAT_UNKNOWN;
}

/** gives termination reason
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 * 
 * return: Termination Status
 */
static
SCIP_DECL_NLPIGETTERMSTAT( nlpiGetTermstatXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_NLPITERMSTAT_OTHER;
}

/** gives primal solution
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues pointer to store primal values
 * 
 * output:
 *  - primalvalues primal values of solution
 */
static
SCIP_DECL_NLPIGETSOLUTION( nlpiGetSolutionXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETSTATISTICS( nlpiGetStatisticsXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETWARMSTARTSIZE( nlpiGetWarmstartSizeXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETWARMSTARTMEMO( nlpiGetWarmstartMemoXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPISETWARMSTARTMEMO( nlpiSetWarmstartMemoXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETINTPAR( nlpiGetIntParXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

/** sets integer parameter of NLP
 * 
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival parameter value
 */
static
SCIP_DECL_NLPISETINTPAR( nlpiSetIntParXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETREALPAR( nlpiGetRealParXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}
 
/** sets floating point parameter of NLP
 * 
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval parameter value
 */
static
SCIP_DECL_NLPISETREALPAR( nlpiSetRealParXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

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
SCIP_DECL_NLPIGETSTRINGPAR( nlpiGetStringParXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

/** sets string parameter of NLP
 * 
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - sval parameter value
 */
static
SCIP_DECL_NLPISETSTRINGPAR( nlpiSetStringParXxx )
{
   SCIPerrorMessage("method of xxx nonlinear solver is not implemented\n");
   SCIPABORT();
   
   return SCIP_OKAY;
}

/*
 * NLP solver interface specific interface methods
 */

/** create solver interface for Xxx solver */
SCIP_RETCODE SCIPcreateNlpSolverXxx(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi                /**< pointer to buffer for nlpi address */
)
{
   SCIP_NLPIDATA* nlpidata;
   
   assert(blkmem != NULL);
   assert(nlpi   != NULL);
  
   /* create xxx solver interface data */
   nlpidata = NULL;
   /* TODO: (optional) create solver interface specific data here */

   /* create solver interface */
   SCIP_CALL( SCIPnlpiCreate(nlpi,
      NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
      nlpiCopyXxx, nlpiFreeXxx, nlpiGetSolverPointerXxx,
      nlpiCreateProblemXxx, nlpiFreeProblemXxx, nlpiGetProblemPointerXxx,
      nlpiAddVarsXxx, nlpiAddConstraintsXxx, nlpiSetObjectiveXxx, 
      nlpiChgVarBoundsXxx, nlpiChgConsSidesXxx, nlpiDelVarSetXxx, nlpiDelConstraintSetXxx,
      nlpiChgLinearCoefsXxx, nlpiChgQuadraticCoefsXxx, nlpiChgExprtreeXxx, nlpiChgNonlinCoefXxx,
      nlpiChgObjConstantXxx, nlpiSetInitialGuessXxx, nlpiSolveXxx, nlpiGetSolstatXxx, nlpiGetTermstatXxx,
      nlpiGetSolutionXxx, nlpiGetStatisticsXxx,
      nlpiGetWarmstartSizeXxx, nlpiGetWarmstartMemoXxx, nlpiSetWarmstartMemoXxx,
      nlpiGetIntParXxx, nlpiSetIntParXxx, nlpiGetRealParXxx, nlpiSetRealParXxx, nlpiGetStringParXxx, nlpiSetStringParXxx,
      nlpidata) );

   return SCIP_OKAY;
}
