/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* NLP interface for the CONOPT solver */

/**@file    nlpi_conopt.c
 * @ingroup DEFPLUGINS_NLPI
 * @brief   CONOPT NLP interface
 * @author  you
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/nlpi_conopt.h"
#include "scip/nlpioracle.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_nlp.h"
#include "scip/scip_nlpi.h"
#include "scip/pub_message.h"

#include "coiheader.h"
#include "scip_message.h"

#define NLPI_NAME              "conopt"                    /**< short concise name of solver */
#define NLPI_DESC              "solver interface template" /**< description of solver */
#define NLPI_PRIORITY          0                           /**< priority of NLP solver */

#define SCIP_DEBUG

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
   SCIP*                 scip;               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle;             /**< Oracle-helper to store and evaluate NLP */

   SCIP_Bool             firstrun;
   SCIP_NLPSOLSTAT       solstat;            /**< solution status from last NLP solve */
   SCIP_NLPTERMSTAT      termstat;           /**< termination status from last NLP solve */
   SCIP_Real             solvetime;          /**< time spend for last NLP solve */
   int                   niterations;        /**< number of iterations for last NLP solve */
   SCIP_Real             objval;             /**< objective value from last run */

   coiHandle_t           CntVect;            /**< pointer to CONOPT Control Vector */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** Implementations of CONOPT callbacks */

/* TODO implement solution routine */

/* TODO complete this */
int COI_CALLCONV Tut_ReadMatrix(double LOWER[], double CURR[], double UPPER[], int VSTA[], int TYPE[], double RHS[],
      int ESTA[], int COLSTA[], int ROWNO[], double VALUE[], int NLFLAG[], int NUMVAR, int NUMCON, int NUMNZ, void* USRMEM)
{
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;
   SCIP_NLPIORACLE* oracle;

   assert(problem != NULL);
   assert(problem->scip != NULL);

   oracle = problem->oracle;
   assert(oracle != NULL);

   lbs = SCIPnlpiOracleGetVarLbs(oracle);
   ubs = SCIPnlpiOracleGetVarUbs(oracle);

   for( int i = 0; i < NUMVAR; i++ )
   {
      LOWER[i] = lbs[i];
      UPPER[i] = ubs[i];
      CURR[i] = 0.0; /* TODO get proper initial values */
   }

   /* TODO check if VSTA is needed - seems to not be the case */

   for ( int i = 0; i < NUMCON; i++ )
   {
      SCIP_Real lhs = SCIPnlpiOracleGetConstraintLhs(oracle, i);
      SCIP_Real rhs = SCIPnlpiOracleGetConstraintRhs(oracle, i);

      assert(!SCIPisInfinity(problem->scip, lhs) || !SCIPisInfinity(problem->scip, rhs));

      if( !SCIPisInfinity(problem->scip, lhs) && !SCIPisInfinity(problem->scip, rhs) )
      {
         TYPE[i] = 0; /* an equality or ranged row modelled as equality */

         if( !SCIPisEQ(problem->scip, lhs, rhs) )
         {
            /* range constraint lhs <= g(x) <= rhs:
             * reformulate as g(x) - s = 0 and lhs <= s <= rhs */
            RHS[i] = 0.0;
            /* TODO introduce slack */
         }
         else
            RHS[i] = lhs;
      }
      else if( !SCIPisInfinity(problem->scip, lhs) )
      {
         TYPE[i] = 1;
         RHS[i] = lhs;
      }
      else
      {
         TYPE[i] = 2;
         RHS[i] = rhs;
      }
   }

   /* TODO Jacobian information */

   return 0;
}

/* callback for CONOPT's standard output */
static int COI_CALLCONV Std_Message(int SMSG, int DMSG, int NMSG, char* MSGV[], void* USRMEM)
{
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;

   assert(problem != NULL);
   assert(problem->scip != NULL);

   for( int i = 0; i < SMSG; i++ )
      SCIPinfoMessage(problem->scip, NULL, "%s\n", MSGV[i]);

   return 0;
}

/* callback for CONOPT's standard error output */
int COI_CALLCONV Std_ErrMsg(int ROWNO, int COLNO, int POSNO, const char* MSG, void* USRMEM)
{
   if( ROWNO == -1 )
      SCIPerrorMessage("Variable %d : ", COLNO);
   else if( COLNO == -1 )
      SCIPerrorMessage("Constraint %d : ", ROWNO);
   else
      SCIPerrorMessage("Variable %d appearing in constraint %d : ", COLNO, ROWNO);
   SCIPerrorMessage("%s\n", MSG);

   return 0;
}

/* callback for CONOPT to report the solving statuses */
int COI_CALLCONV Std_Status(int MODSTA, int SOLSTA, int ITER, double OBJVAL, void* USRMEM)
{
   SCIP* scip;
   SCIP_NLPIPROBLEM* problem = (SCIP_NLPIPROBLEM*)USRMEM;

   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   SCIPdebugMsg(scip, "CONOPT has finished optimizing\n");
   SCIPdebugMsg(scip, "Iteration count = %8d\n", ITER);
   SCIPdebugMsg(scip, "Objective value = %10f\n", OBJVAL);

   problem->niterations = ITER;
   problem->objval = OBJVAL;

   switch( MODSTA )
   {
      case 1:
         SCIPdebugMsg(scip, "NLP problem solved to global optimality\n");
         problem->solstat = SCIP_NLPSOLSTAT_GLOBOPT;
         break;
      case 2:
         SCIPdebugMsg(scip, "NLP problem solved to local optimality\n");
         problem->solstat = SCIP_NLPSOLSTAT_LOCOPT;
         break;
      case 3:
         SCIPdebugMsg(scip, "NLP problem unbounded\n");
         problem->solstat = SCIP_NLPSOLSTAT_UNBOUNDED;
         break;
      case 4:
         SCIPdebugMsg(scip, "NLP problem infeasible\n");
         problem->solstat = SCIP_NLPSOLSTAT_GLOBINFEASIBLE;
         break;
      case 5:
         SCIPdebugMsg(scip, "NLP problem locally infeasible\n");
         problem->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
         break;
      case 6: /* intermediate infeasible */
      case 7: /* intermediate non-optimal */
      case 12: /* unknown error */
      case 13:
         SCIPdebugMsg(scip, "NLP problem status unknown (CONOPT status %d)\n", MODSTA);
         problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         break;
      default:
         SCIPerrorMessage("CONOPT returned an unexpected solution status %d\n", MODSTA);
         problem->solstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   }

   switch( SOLSTA )
   {
      case 1:
         problem->termstat = SCIP_NLPTERMSTAT_OKAY;
         SCIPdebugMsg(scip, "CONOPT terminated with normal status.\n");
         break;
      case 2:
         problem->termstat = SCIP_NLPTERMSTAT_ITERLIMIT;
         SCIPdebugMsg(scip, "CONOPT terminated due to an iteration limit.\n");
         break;
      case 3:
         problem->termstat = SCIP_NLPTERMSTAT_TIMELIMIT;
         SCIPdebugMsg(scip, "CONOPT terminated due to a time limit.\n");
         break;
      case 5:
         problem->termstat = SCIP_NLPTERMSTAT_EVALERROR;
         SCIPdebugMsg(scip, "CONOPT terminated due to evaluation errors.\n");
         break;
      case 8:
         problem->termstat = SCIP_NLPTERMSTAT_INTERRUPT;
         SCIPdebugMsg(scip, "CONOPT interrupted by user.\n");
         break;
      case 4: /* terminated by solver */
      case 6: /* unknown */
      case 9: /* error: setup failure */
      case 10: /* error: solver failure */
      case 11: /* error: internal solver error */
         SCIPdebugMsg(scip, "CONOPT terminated with status %d\n", SOLSTA);
         problem->termstat = SCIP_NLPTERMSTAT_OTHER;
         break;
      default:
         SCIPerrorMessage("CONOPT returned an unexpected termination status %d\n", SOLSTA);
         problem->termstat = SCIP_NLPTERMSTAT_OTHER;
   }

   return 0;
}

/* NLPI local methods */

/** sets the solstat and termstat to unknown and other, resp. */
static
void invalidateSolution(
   SCIP_NLPIPROBLEM*     problem             /**< data structure of problem */
   )
{
   assert(problem != NULL);

   problem->solstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   problem->termstat = SCIP_NLPTERMSTAT_OTHER;
}

static SCIP_RETCODE initConopt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to NLPI problem structure */
)
{
   int COI_Error = 0; /* CONOPT error counter */
   int nrangeconss = 0;
   int nconss;
   const int* jacoffset;
   const int* hessoffset;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   nconss = SCIPnlpiOracleGetNConstraints(problem->oracle);

   /* count range constraints: because CONOPT doesn't support them directly, will need to add a slack variable for each ranged constraint */
   for( int i = 0; i < SCIPnlpiOracleGetNConstraints(problem->oracle); i++ )
   {
      SCIP_Real lhs = SCIPnlpiOracleGetConstraintLhs(problem->oracle, i);
      SCIP_Real rhs = SCIPnlpiOracleGetConstraintRhs(problem->oracle, i);

      if( !SCIPisEQ(scip, lhs, rhs) )
         nrangeconss++;
   }

   /* inform CONOPT about problem sizes */

   COI_Error += COIDEF_NumVar(problem->CntVect, SCIPnlpiOracleGetNVars(problem->oracle) + nrangeconss);
   COI_Error += COIDEF_NumCon(problem->CntVect, nconss);

   /* jacobian information */
   SCIP_CALL( SCIPnlpiOracleGetJacobianSparsity(scip, problem->oracle, &jacoffset, NULL) );
   COI_Error += COIDEF_NumNz(problem->CntVect, jacoffset[nconss] + nrangeconss); /* each slack var adds a Jacobian nnz */
   // COI_Error += COIDEF_NumNlNz(problem->CntVect, ); TODO find out how to go about nonlinear jacobian entries

   /* hessian information */
   SCIP_CALL( SCIPnlpiOracleGetHessianLagSparsity(scip, problem->oracle, &hessoffset, NULL) );
   COI_Error += COIDEF_NumHess(problem->CntVect, hessoffset[nconss]);

   /* tell conopt to minimise the objective (the oracle always gives us a minimisation problem) */
   COI_Error += COIDEF_OptDir(problem->CntVect, -1);

   /* oracle gives objective as a constraint, hence use ObjCon (not ObjVar) here */
   COI_Error += COIDEF_ObjCon(problem->CntVect, nconss); /* TODO: decide which index to use for the objective constraint; nconss? */

   /* register callback routines */

   COI_Error += COIDEF_Message(problem->CntVect, &Std_Message);
   COI_Error += COIDEF_ErrMsg(problem->CntVect, &Std_ErrMsg);
   COI_Error += COIDEF_Status(problem->CntVect, &Std_Status);
   // COI_Error +=  COIDEF_Solution  ( CntVect, &Std_Solution );  /* Register the callback Solution    */
   // COI_Error +=  COIDEF_ReadMatrix( CntVect, &Tut_ReadMatrix); /* Register the callback ReadMatrix  */
   // COI_Error +=  COIDEF_FDEval    ( CntVect, &Tut_FDEval);     /* Register the callback FDEval      */

   /* pass the problem pointer to CONOPT, so that it may be used in CONOPT callbacks */
   COI_Error += COIDEF_UsrMem(problem->CntVect, (void*)problem);

   if( COI_Error )
      SCIPinfoMessage(scip, NULL, "Error %d encountered during adding variables", COI_Error);

   return SCIP_OKAY;
}

static void freeConopt(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
)
{

}

static void updateConopt(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
)
{

}


/*
 * Callback methods of NLP solver interface
 */

/* TODO: Implement all necessary NLP interface methods. The methods with an #ifdef SCIP_DISABLED_CODE ... #else #define ... are optional */

#ifdef SCIP_DISABLED_CODE
/** copy method of NLP interface (called when SCIP copies plugins) */
static
SCIP_DECL_NLPICOPY(nlpiCopyXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/
#else
#define nlpiCopyConopt NULL
#endif

/** destructor of NLP interface to free nlpi data */
static
SCIP_DECL_NLPIFREE(nlpiFreeConopt)
{
   assert(nlpi != NULL);
   assert(nlpidata != NULL);
   assert(*nlpidata != NULL);

   SCIPfreeBlockMemory(scip, nlpidata);
   assert(*nlpidata == NULL);

   return SCIP_OKAY;
}  /*lint !e715*/

#ifdef SCIP_DISABLED_CODE
/** gets pointer for NLP solver */
static
SCIP_DECL_NLPIGETSOLVERPOINTER(nlpiGetSolverPointerXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return NULL;  /*lint !e527*/
}  /*lint !e715*/
#else
#define nlpiGetSolverPointerConopt NULL
#endif

/** create a problem instance */
static
SCIP_DECL_NLPICREATEPROBLEM(nlpiCreateProblemConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, problem) );
   assert(*problem != NULL);

   (*problem)->firstrun = TRUE;
   (*problem)->scip = scip;

   /* initialize oracle */
   SCIP_CALL( SCIPnlpiOracleCreate(scip, &(*problem)->oracle) );
   SCIP_CALL( SCIPnlpiOracleSetProblemName(scip, (*problem)->oracle, name) );

   coiCreate(&((*problem)->CntVect)); /* TODO handle errors here? */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** free a problem instance */
static
SCIP_DECL_NLPIFREEPROBLEM(nlpiFreeProblemConopt)
{
   assert(nlpi     != NULL);
   assert(problem  != NULL);
   assert(*problem != NULL);

   if( (*problem)->oracle != NULL )
   {
      SCIP_CALL( SCIPnlpiOracleFree(scip, &(*problem)->oracle) );
   }

   coiFree(&((*problem)->CntVect));

   SCIPfreeBlockMemory(scip, problem);
   *problem = NULL;

   return SCIP_OKAY;
}  /*lint !e715*/

#ifdef SCIP_DISABLED_CODE
/** gets pointer to solver-internal problem instance */
static
SCIP_DECL_NLPIGETPROBLEMPOINTER(nlpiGetProblemPointerXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return NULL;  /*lint !e527*/
}  /*lint !e715*/
#else
#define nlpiGetProblemPointerConopt NULL
#endif

/** add variables */
static
SCIP_DECL_NLPIADDVARS(nlpiAddVarsConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleAddVars(scip, problem->oracle, nvars, lbs, ubs, varnames) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/


/** add constraints */
static
SCIP_DECL_NLPIADDCONSTRAINTS(nlpiAddConstraintsConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleAddConstraints(scip, problem->oracle, nconss, lhss, rhss,
         nlininds, lininds, linvals, exprs, names) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets or overwrites objective, a minimization problem is expected */
static
SCIP_DECL_NLPISETOBJECTIVE(nlpiSetObjectiveConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   /* TODO do we need to reset CONOPT? */
   if( expr != NULL || SCIPnlpiOracleIsConstraintNonlinear(problem->oracle, -1) )
      problem->firstrun = TRUE;

   SCIP_CALL( SCIPnlpiOracleSetObjective(scip, problem->oracle, constant, nlins, lininds, linvals, expr) );

   invalidateSolution(problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change variable bounds */
static
SCIP_DECL_NLPICHGVARBOUNDS(nlpiChgVarBoundsConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   /* TODO do this on the conopt side */

   SCIP_CALL( SCIPnlpiOracleChgVarBounds(scip, problem->oracle, nvars, indices, lbs, ubs) );

   invalidateSolution(problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change constraint bounds */
static
SCIP_DECL_NLPICHGCONSSIDES(nlpiChgConsSidesConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   /* TODO do this on the conopt side */

   SCIP_CALL( SCIPnlpiOracleChgConsSides(scip, problem->oracle, nconss, indices, lhss, rhss) );

   invalidateSolution(problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** delete a set of variables */
static
SCIP_DECL_NLPIDELVARSET(nlpiDelVarSetConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelVarSet(scip, problem->oracle, dstats) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** delete a set of constraints */
static
SCIP_DECL_NLPIDELCONSSET(nlpiDelConstraintSetConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelConsSet(scip, problem->oracle, dstats) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** changes (or adds) linear coefficients in a constraint or objective */
static
SCIP_DECL_NLPICHGLINEARCOEFS(nlpiChgLinearCoefsConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgLinearCoefs(scip, problem->oracle, idx, nvals, varidxs, vals) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** replaces the expression tree of a constraint or objective */
static
SCIP_DECL_NLPICHGEXPR(nlpiChgExprConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgExpr(scip, problem->oracle, idxcons, expr) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change the constant offset in the objective */
static
SCIP_DECL_NLPICHGOBJCONSTANT(nlpiChgObjConstantConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgObjConstant(scip, problem->oracle, objconstant) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/* TODO */
#ifdef SCIP_DISABLED_CODE
/** sets initial guess */
static
SCIP_DECL_NLPISETINITIALGUESS(nlpiSetInitialGuessXyz)
{
   SCIPerrorMessage("method of xyz nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/
#else
#define nlpiSetInitialGuessConopt NULL
#endif

/** try to solve NLP
 *
 * Note that SCIP will already have reset a timelimit of SCIP_REAL_MAX to the time remaining for the SCIP solve in SCIPnlpiSolve().
 */
static
SCIP_DECL_NLPISOLVE(nlpiSolveConopt)
{
   int COI_Error = 0; /* CONOPT error counter */

   assert(nlpi != NULL);
   assert(problem != NULL);

   SCIPdebugMsg(scip, "solve with parameters " SCIP_NLPPARAM_PRINT(param));

   SCIP_CALL( SCIPnlpiOracleResetEvalTime(scip, problem->oracle) );

   if( param.timelimit == 0.0 )
   {
      /* there is nothing we can do if we are not given any time */
      problem->niterations = 0;
      problem->solvetime = 0.0;
      problem->termstat = SCIP_NLPTERMSTAT_TIMELIMIT;
      problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;

      return SCIP_OKAY;
   }

   problem->niterations = -1;
   problem->solvetime  = -1.0;

   /* TODO handle param.verblevel */

   // /* if warmstart parameter is disabled, then we will not warmstart */
   // if( !param.warmstart )
      // problem->warmstart = FALSE;

   /* initialize Conopt data if necessary */
   if( problem->firstrun )
   { /* TODO implement the functions */
      freeConopt(problem);
      initConopt(scip, nlpi, problem);
      problem->firstrun = FALSE;
   }
   else
   {
      updateConopt(problem);
   }

   /* TODO set parameters */

   /* TODO set initial guess (if available) */



   COI_Error = COI_Solve(problem->CntVect); /* optimize */
   if( COI_Error )
   {
      SCIPinfoMessage(scip, NULL, "Errors encountered during solution, %d", COI_Error);
      exit(1);
   }

   /* TODO interpret conopt result (update solstat, termstat, print some debug message) */

   /* TODO handle some statistics */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives solution status */
static
SCIP_DECL_NLPIGETSOLSTAT(nlpiGetSolstatConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   return problem->solstat;
}  /*lint !e715*/

/** gives termination reason */
static
SCIP_DECL_NLPIGETTERMSTAT(nlpiGetTermstatConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   return problem->termstat;
}  /*lint !e715*/

/** gives primal and dual solution values */
static
SCIP_DECL_NLPIGETSOLUTION(nlpiGetSolutionConopt)
{
   SCIPerrorMessage("method of conopt nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives solve statistics */
static
SCIP_DECL_NLPIGETSTATISTICS(nlpiGetStatisticsConopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(statistics != NULL);

   statistics->niterations = problem->niterations;
   statistics->totaltime = problem->solvetime;
   statistics->evaltime = SCIPnlpiOracleGetEvalTime(scip, problem->oracle);
   // statistics->consviol = problem->wsp->FeasOrigMax;
   statistics->boundviol = 0.0;

   return SCIP_OKAY;
}  /*lint !e715*/

/*
 * NLP solver interface specific interface methods
 */

/** create solver interface for Conopt solver and includes it into SCIP */
SCIP_RETCODE SCIPincludeNlpSolverConopt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLPIDATA* nlpidata;

   /* create Conopt solver interface data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nlpidata) );

   /* create and include solver interface */
   SCIP_CALL( SCIPincludeNlpi(scip,
         NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
         nlpiCopyConopt, nlpiFreeConopt, nlpiGetSolverPointerConopt,
         nlpiCreateProblemConopt, nlpiFreeProblemConopt, nlpiGetProblemPointerConopt,
         nlpiAddVarsConopt, nlpiAddConstraintsConopt, nlpiSetObjectiveConopt,
         nlpiChgVarBoundsConopt, nlpiChgConsSidesConopt, nlpiDelVarSetConopt, nlpiDelConstraintSetConopt,
         nlpiChgLinearCoefsConopt, nlpiChgExprConopt, nlpiChgObjConstantConopt,
         nlpiSetInitialGuessConopt, nlpiSolveConopt,
         nlpiGetSolstatConopt, nlpiGetTermstatConopt, nlpiGetSolutionConopt, nlpiGetStatisticsConopt,
         nlpidata) );

   return SCIP_OKAY;
}
