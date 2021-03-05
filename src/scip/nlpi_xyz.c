/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    nlpi_xyz.c
 * @ingroup NLPIS
 * @brief   XYZ NLP interface
 * @author  you
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/nlpi_xyz.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_nlp.h"
#include "scip/scip_nlpi.h"

#define NLPI_NAME              "xyz"                       /* short concise name of solver */
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

/** copy method of NLP interface (called when SCIP copies plugins) */
static
SCIP_DECL_NLPICOPY(nlpiCopyXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** destructor of NLP interface to free nlpi data */
static
SCIP_DECL_NLPIFREE(nlpiFreeXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gets pointer for NLP solver */
static
SCIP_DECL_NLPIGETSOLVERPOINTER(nlpiGetSolverPointerXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return NULL;  /*lint !e527*/
}  /*lint !e715*/

/** creates a problem instance */
static
SCIP_DECL_NLPICREATEPROBLEM(nlpiCreateProblemXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** free a problem instance */
static
SCIP_DECL_NLPIFREEPROBLEM(nlpiFreeProblemXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gets pointer to solver-internal problem instance */
static
SCIP_DECL_NLPIGETPROBLEMPOINTER(nlpiGetProblemPointerXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return NULL;  /*lint !e527*/
}  /*lint !e715*/

/** add variables */
static 
SCIP_DECL_NLPIADDVARS(nlpiAddVarsXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/


/** add constraints */
static
SCIP_DECL_NLPIADDCONSTRAINTS(nlpiAddConstraintsXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets or overwrites objective, a minimization problem is expected */
static
SCIP_DECL_NLPISETOBJECTIVE(nlpiSetObjectiveXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change variable bounds */
static
SCIP_DECL_NLPICHGVARBOUNDS(nlpiChgVarBoundsXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change constraint bounds */
static
SCIP_DECL_NLPICHGCONSSIDES(nlpiChgConsSidesXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** delete a set of variables */
static
SCIP_DECL_NLPIDELVARSET(nlpiDelVarSetXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** delete a set of constraints */
static
SCIP_DECL_NLPIDELCONSSET(nlpiDelConstraintSetXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** changes (or adds) linear coefficients in a constraint or objective */
static
SCIP_DECL_NLPICHGLINEARCOEFS(nlpiChgLinearCoefsXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** replaces the expression tree of a constraint or objective */
static
SCIP_DECL_NLPICHGEXPR(nlpiChgExprXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change the constant offset in the objective */
static
SCIP_DECL_NLPICHGOBJCONSTANT(nlpiChgObjConstantXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets initial guess for primal variables */
static
SCIP_DECL_NLPISETINITIALGUESS(nlpiSetInitialGuessXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** tries to solve NLP */
static
SCIP_DECL_NLPISOLVE(nlpiSolveXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives solution status */
static
SCIP_DECL_NLPIGETSOLSTAT(nlpiGetSolstatXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_NLPSOLSTAT_UNKNOWN;  /*lint !e527*/
}  /*lint !e715*/

/** gives termination reason */
static
SCIP_DECL_NLPIGETTERMSTAT(nlpiGetTermstatXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_NLPTERMSTAT_OTHER;  /*lint !e527*/
}  /*lint !e715*/

/** gives primal and dual solution values */
static
SCIP_DECL_NLPIGETSOLUTION(nlpiGetSolutionXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives solve statistics */
static
SCIP_DECL_NLPIGETSTATISTICS(nlpiGetStatisticsXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives required size of a buffer to store a warmstart object */
static
SCIP_DECL_NLPIGETWARMSTARTSIZE(nlpiGetWarmstartSizeXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** stores warmstart information in buffer */
static
SCIP_DECL_NLPIGETWARMSTARTMEMO(nlpiGetWarmstartMemoXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets warmstart information in solver */
static
SCIP_DECL_NLPISETWARMSTARTMEMO(nlpiSetWarmstartMemoXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gets integer parameter of NLP */
static
SCIP_DECL_NLPIGETINTPAR(nlpiGetIntParXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets integer parameter of NLP */
static
SCIP_DECL_NLPISETINTPAR(nlpiSetIntParXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gets floating point parameter of NLP */
static
SCIP_DECL_NLPIGETREALPAR(nlpiGetRealParXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets floating point parameter of NLP */
static
SCIP_DECL_NLPISETREALPAR(nlpiSetRealParXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gets string parameter of NLP */
static
SCIP_DECL_NLPIGETSTRINGPAR(nlpiGetStringParXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets string parameter of NLP */
static
SCIP_DECL_NLPISETSTRINGPAR(nlpiSetStringParXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/*
 * NLP solver interface specific interface methods
 */

/** create solver interface for Xyz solver and includes it into SCIP */
SCIP_RETCODE SCIPincludeNlpSolverXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLPI* nlpi;
   SCIP_NLPIDATA* nlpidata;

   /* create xyz solver interface data */
   nlpidata = NULL;
   /* TODO: (optional) create solver interface specific data here */

   /* create solver interface */
   SCIP_CALL( SCIPnlpiCreate(&nlpi,
         NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
         nlpiCopyXyz, nlpiFreeXyz, nlpiGetSolverPointerXyz,
         nlpiCreateProblemXyz, nlpiFreeProblemXyz, nlpiGetProblemPointerXyz,
         nlpiAddVarsXyz, nlpiAddConstraintsXyz, nlpiSetObjectiveXyz,
         nlpiChgVarBoundsXyz, nlpiChgConsSidesXyz, nlpiDelVarSetXyz, nlpiDelConstraintSetXyz,
         nlpiChgLinearCoefsXyz, nlpiChgExprXyz,
         nlpiChgObjConstantXyz, nlpiSetInitialGuessXyz, nlpiSolveXyz, nlpiGetSolstatXyz, nlpiGetTermstatXyz,
         nlpiGetSolutionXyz, nlpiGetStatisticsXyz,
         nlpiGetWarmstartSizeXyz, nlpiGetWarmstartMemoXyz, nlpiSetWarmstartMemoXyz,
         nlpiGetIntParXyz, nlpiSetIntParXyz, nlpiGetRealParXyz, nlpiSetRealParXyz, nlpiGetStringParXyz, nlpiSetStringParXyz,
         nlpidata) );
   assert(nlpi != NULL);

   /* include NLPI into SCIP */
   SCIP_CALL( SCIPincludeNlpi(scip, nlpi) );

   return SCIP_OKAY;
}
