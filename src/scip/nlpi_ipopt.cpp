/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nlpi_ipopt.cpp,v 1.4 2009/08/09 15:49:58 bzfviger Exp $"

/**@file    nlpi_ipopt.cpp
 * @brief   Ipopt NLP interface
 * @ingroup NLPINTERFACES
 * @author  Stefan Vigerske
 */

/* @TODO warm starts
 * @TODO ScipJournal to redirect Ipopt output 
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/nlpi_ipopt.h"
#include "scip/nlpi.h"
#include "scip/nlpi_oracle.h"

#include <new>      /* for std::bad_alloc */
#include <cstring>  /* for memcpy */

#include "IpIpoptApplication.hpp"
namespace Ipopt
{
   class IpoptNLP;
   class IpoptData;
}
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpSolveStatistics.hpp"

using namespace Ipopt;

#define NLPI_NAME              "Ipopt"           /**< short concise name of solver */
#define NLPI_DESC              "Ipopt interface" /**< description of solver */
#define NLPI_TYPE              "IP"              /**< solver type */
#define NLPI_PRIORITY          0                 /**< priority */

#ifdef SCIP_DEBUG
#define DEFAULT_PRINTLEVEL     J_WARNING         /**< default print level of Ipopt */
#else
#define DEFAULT_PRINTLEVEL     J_STRONGWARNING   /**< default print level of Ipopt */
#endif
#define DEFAULT_MAXITER        3000              /**< default iteration limit for Ipopt */

class ScipNLP;

struct SCIP_NlpiData
{
public:
   SCIP_NLPIORACLE*            oracle;
   
   SmartPtr<IpoptApplication>  ipopt;
   SmartPtr<ScipNLP>           nlp;
   
   SCIP_Bool                   firstrun;        /**< whether the next NLP solve will be the first one (with the current problem structure) */
   SCIP_Real*                  initguess;       /**< initial values for primal variables, or NULL if not known */
   
   SCIP_NLPSOLSTAT             lastsolstat;     /**< solution status from last run */
   SCIP_NLPITERMSTAT           lasttermstat;    /**< termination status from last run */
   SCIP_Real*                  lastsol;         /**< solution from last run, if available */
   int                         lastniter;       /**< number of iterations in last run */
   SCIP_Real                   lasttime;        /**< time spend in last run */

   SCIP_NlpiData()
   : oracle(NULL),
     firstrun(TRUE), initguess(NULL),
     lastsolstat(SCIP_NLPSOLSTAT_UNKNOWN), lasttermstat(SCIP_NLPITERMSTAT_OTHER), lastsol(NULL),
     lastniter(-1), lasttime(-1.0)
   { }
};

class ScipNLP : public TNLP
{
private:
   SCIP*             scip;
   SCIP_NLPIDATA*    nlpidata;

public:

   ScipNLP(SCIP* scip_ = NULL, SCIP_NLPIDATA* nlpidata_ = NULL)
   : scip(scip_), nlpidata(nlpidata_)
   { }

   ~ScipNLP() { }
   
   void setSCIP(SCIP* scip_)
   {
      assert(scip_ != NULL);
      scip = scip_;
   }

   void setNLPIDATA(SCIP_NLPIDATA* nlpidata_)
   {
      assert(nlpidata_ != NULL);
      nlpidata = nlpidata_;
   }


   /** Method to return some info about the nlp */
   bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style);

   /** Method to return the bounds for my problem */
   bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u);

   /** Method to return the starting point for the algorithm */
   bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda);

   /** Method to return the variables linearity. */
   bool get_variables_linearity(Index n, LinearityType* var_types);

   /** Method to return the constraint linearity. */
   bool get_constraints_linearity(Index m, LinearityType* const_types);

   /** Method to return the objective value */
   bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

   /** Method to return the gradient of the objective */
   bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

   /** Method to return the constraint residuals */
   bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

   /** Method to return:
    *   1) The structure of the jacobian (if "values" is NULL)
    *   2) The values of the jacobian (if "values" is not NULL)
    */
   bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values);

   /** Method to return:
    *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
    *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
    */
   bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values);

   /** Method called by the solver at each iteration.
    * Checks whether ^C was hit.
    */
   bool intermediate_callback (AlgorithmMode mode, Index iter, Number obj_value, Number inf_pr, Number inf_du, Number mu, Number d_norm, Number regularization_size, Number alpha_du, Number alpha_pr, Index ls_trials, const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq);

   /** This method is called when the algorithm is complete so the TNLP can store/write the solution.
    */
   void finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g, const Number* lambda, Number obj_value, const IpoptData* data, IpoptCalculatedQuantities* cq);
};

static
void SCIPnlpiIpoptInvalidateSolution(SCIP* scip, SCIP_NLPI* nlpi)
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   
   SCIPfreeMemoryArrayNull(scip, &data->lastsol);
   data->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   data->lasttermstat = SCIP_NLPITERMSTAT_OTHER;
}

/** initializes an NLP interface structure
 * Input:
 *  - nlpi datastructure for solver interface
 *  - name problem name
 */
static
SCIP_DECL_NLPIINIT( nlpiInitIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);

   assert(data->oracle == NULL);
   SCIP_CALL( SCIPnlpiOracleCreate(scip, &data->oracle) );
   SCIP_CALL( SCIPnlpiOracleInit(scip, data->oracle) );
   
   try
   {
      data->ipopt = new IpoptApplication();
      if (IsNull(data->ipopt))
         throw std::bad_alloc();
   }
   catch (std::bad_alloc)
   {
      SCIPerrorMessage("Not enough memory to allocate IpoptApplication.\n");
      return SCIP_NOMEMORY;
   }

   try
   {
      data->nlp = new ScipNLP(scip, data);
      if (IsNull(data->nlp))
         throw std::bad_alloc();
   }
   catch (std::bad_alloc)
   {
      SCIPerrorMessage("Not enough memory to allocate ScipNLP.\n");
      return SCIP_NOMEMORY;
   }
   
   data->ipopt->Options()->SetIntegerValue("print_level", DEFAULT_PRINTLEVEL);
   /* data->ipopt->Options()->SetStringValue("print_timing_statistics", "yes"); */
   data->ipopt->Options()->SetStringValue("mu_strategy", "adaptive");
   data->ipopt->Options()->SetStringValue("expect_infeasible_problem", "yes");
   /* it seem to be better to let Ipopt relax bounds a bit to ensure that a relative interior exists;
    * however, if we relax the bounds too much, then the solutions tend to be slightly infeasible */
   data->ipopt->Options()->SetNumericValue("tol", SCIPfeastol(scip)/2);
   data->ipopt->Options()->SetNumericValue("bound_relax_factor", SCIPfeastol(scip)/2);
   data->ipopt->Options()->SetNumericValue("constr_viol_tol", 0.75*SCIPfeastol(scip));
   data->ipopt->Options()->SetIntegerValue("max_iter", DEFAULT_MAXITER);
   data->ipopt->Options()->SetNumericValue("nlp_lower_bound_inf", -SCIPinfinity(scip), false);
   data->ipopt->Options()->SetNumericValue("nlp_upper_bound_inf",  SCIPinfinity(scip), false);
   data->ipopt->Options()->SetNumericValue("diverging_iterates_tol", SCIPinfinity(scip), false);
   /* data->ipopt->Options()->SetStringValue("dependency_detector", "ma28"); */
   /* data->ipopt->Options()->SetStringValue("hessian_approximation", "limited-memory"); */
#ifdef SCIP_DEBUG
   data->ipopt->Options()->SetStringValue("derivative_test", "second-order");
#endif

   data->ipopt->Initialize(""); // do not read default option file

   return SCIP_OKAY;
}

/** frees nlpi solver data */
static
SCIP_DECL_NLPIFREE( nlpiFreeIpopt )
{
   assert(scip != NULL);
   assert(data != NULL);
  
   if (data->oracle)
      SCIP_CALL( SCIPnlpiOracleFree(scip, &data->oracle) );
   
   SCIPfreeMemoryArrayNull(scip, &data->initguess);
   SCIPfreeMemoryArrayNull(scip, &data->lastsol);
   
   delete data;

   return SCIP_OKAY;
}

/** add constraints */
static
SCIP_DECL_NLPIADDVARS( nlpiAddVarsIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleAddVars(scip, data->oracle, nvars, lb, ub, varnames) );
   
   data->firstrun = TRUE;
   SCIPfreeMemoryArrayNull(scip, &data->initguess);
   SCIPnlpiIpoptInvalidateSolution(scip, nlpi);
   
   return SCIP_OKAY;
}

/** add restrictions to the solver */
static
SCIP_DECL_NLPIADDCONSTRAINTS( nlpiAddConstraintsIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);
 
   SCIP_CALL( SCIPnlpiOracleAddConstraints(scip, data->oracle,
      ncons, lhs, rhs,
      linoffset, linind, linval,
      nquadrows, quadrowidx, quadoffset, quadind, quadval,
      exprvaridx, exprtree, names) );

   data->firstrun = TRUE;
   SCIPnlpiIpoptInvalidateSolution(scip, nlpi);

   return SCIP_OKAY;
}

/** Overwrite objective, may change sparsity pattern */
static
SCIP_DECL_NLPISETOBJECTIVE( nlpiSetObjectiveIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleSetObjective(scip, data->oracle,
      constant, nlin, linind, linval,
      nquadcols, quadcols, quadoffset, quadind, quadval,
      exprvaridx, exprtree) );

   data->firstrun = TRUE;
   SCIPnlpiIpoptInvalidateSolution(scip, nlpi);

   return SCIP_OKAY;
}

/** change variable bounds */
static
SCIP_DECL_NLPICHGVARBOUNDS( nlpiChgVarBoundsIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);
    
   SCIP_CALL( SCIPnlpiOracleChgVarBounds(scip, data->oracle, nvars, indices, lb, ub) );

   SCIPnlpiIpoptInvalidateSolution(scip, nlpi);

   return SCIP_OKAY;
}

/** change constraint bounds */
static
SCIP_DECL_NLPICHGCONSBOUNDS( nlpiChgConsBoundsIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);
   
   int* realindices;
   SCIP_CALL( SCIPallocBufferArray(scip, &realindices, ncons) );
   for (int i = 0; i < ncons; ++i)
      realindices[i] = indices[i];
   
   SCIP_CALL( SCIPnlpiOracleChgConsBounds(scip, data->oracle, ncons, realindices, lb, ub) );
   
   SCIPfreeBufferArray(scip, &realindices);

   SCIPnlpiIpoptInvalidateSolution(scip, nlpi);

   return SCIP_OKAY;
}

/** delete a set of constraints */
static
SCIP_DECL_NLPIDELVARSET( nlpiDelVarSetIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleDelVarSet(scip, data->oracle, dstat) );

   data->firstrun = TRUE;
   SCIPfreeMemoryArrayNull(scip, &data->initguess); // @TODO keep initguess for remaining variables 

   SCIPnlpiIpoptInvalidateSolution(scip, nlpi);

   return SCIP_OKAY;
}

/** Removes a row set from problem */
static
SCIP_DECL_NLPIDELCONSSET( nlpiDelConstraintSetIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleDelConsSet(scip, data->oracle, dstat) );

   data->firstrun = TRUE;

   SCIPnlpiIpoptInvalidateSolution(scip, nlpi);

   return SCIP_OKAY;
}

/** Changes coefficient in linear part of restriction */
static
SCIP_DECL_NLPICHGLINEARCOEFS( nlpiChgLinearCoefsIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleChgLinearCoefs(scip, data->oracle, cons, nvals, varidx, value) );
   SCIPnlpiIpoptInvalidateSolution(scip, nlpi);

   return SCIP_OKAY;
}

/** Changes coefficient in quadratic part of restriction */
static
SCIP_DECL_NLPICHGQUADCOEFS( nlpiChgQuadraticCoefsIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);
   SCIP_CALL( SCIPnlpiOracleChgQuadCoefs(scip, data->oracle, cons, nentries, row, col, value) );
   SCIPnlpiIpoptInvalidateSolution(scip, nlpi);

   return SCIP_OKAY;
}

/** Change modifiable constant of an expression */
static
SCIP_DECL_NLPICHGNONLINCOEF( nlpiChgNonlinCoefIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
/*    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);
*/ 
   SCIPerrorMessage("ChgNonlinCoef method of Ipopt nonlinear solver is not implemented\n");

   SCIPnlpiIpoptInvalidateSolution(scip, nlpi);

   return SCIP_ERROR;
}

/** Set initial estimate as a ,,good solution'' */
static
SCIP_DECL_NLPISETINITIALGUESS( nlpiSetInitialGuessIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);

   if (!data->initguess)
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &data->initguess, values, SCIPnlpiOracleGetNVars(data->oracle)) );
   else
      memcpy(data->initguess, values, SCIPnlpiOracleGetNVars(data->oracle) * sizeof(SCIP_Real));

   return SCIP_OKAY;
}

/** Solve specified problem */
static
SCIP_DECL_NLPISOLVE( nlpiSolveIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPnlpiOraclePrintProblem(scip, data->oracle, NULL) );
#endif
/*
   FILE* gmsfile = fopen("nlpiprob.gms", "w");
   SCIP_CALL( SCIPnlpiOraclePrintProblemGams(scip, data->oracle, data->initguess, gmsfile) );
   fclose(gmsfile);
*/ 
   assert(IsValid(data->ipopt));
   assert(IsValid(data->nlp));
   
   data->nlp->setSCIP(scip);
   data->nlp->setNLPIDATA(data);
   
   data->lastniter = -1;
   data->lasttime  = -1.0;
   
   ApplicationReturnStatus status;
   try
   {
      if (data->firstrun)
         status = data->ipopt->OptimizeTNLP(GetRawPtr(data->nlp));
      else
         status = data->ipopt->ReOptimizeTNLP(GetRawPtr(data->nlp));
      
      // catch the very bad status codes
      switch (status) {
         case Invalid_Problem_Definition:
         case Invalid_Option:
         case Unrecoverable_Exception:
         case NonIpopt_Exception_Thrown:
         case Internal_Error:
            SCIPerrorMessage("Ipopt returned with application return status %d\n", status);
            return SCIP_ERROR;
         case Insufficient_Memory:
            SCIPerrorMessage("Ipopt returned with status \"Insufficient Memory\"\n");
            return SCIP_NOMEMORY;
         case Invalid_Number_Detected:
            SCIPwarningMessage("Ipopt failed because of an invalid number in function or derivative value\n");
            SCIPnlpiIpoptInvalidateSolution(scip, nlpi);
            data->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
            data->lasttermstat = SCIP_NLPITERMSTAT_EVALERR;
         default: ;
      }

      SmartPtr<SolveStatistics> stats = data->ipopt->Statistics();
      if (IsValid(stats))
      {
         data->lastniter = stats->IterationCount();
         data->lasttime  = stats->TotalCPUTime();
      }
   }
   catch (IpoptException except)
   {
      SCIPerrorMessage("Ipopt returned with exception: %s\n", except.Message().c_str());
      return SCIP_ERROR;
   }
   
   data->firstrun = FALSE;

   return SCIP_OKAY;
}

/** Determine termination status of solution */
static
SCIP_DECL_NLPIGETSOLSTAT( nlpiGetSolstatIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);

   return data->lastsolstat;
}

/** Determine termination status of solution */
static
SCIP_DECL_NLPIGETTERMSTAT( nlpiGetSoltermIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);

   return data->lasttermstat;
}

/** Get last found solution, or NULL if no point is available. */
static
SCIP_DECL_NLPIGETSOLUTION( nlpiGetSolutionIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(primalvalues != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   
   *primalvalues = data->lastsol;
   
   return SCIP_OKAY;
}

/** Get last found solution, or NULL if no point is available. */
static
SCIP_DECL_NLPIGETSTATISTICS( nlpiGetStatisticsIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);

   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);

   SCIPnlpStatisticsSetNIterations(statistics, data->lastniter);
   SCIPnlpStatisticsSetTotalTime  (statistics, data->lasttime);
   
   return SCIP_OKAY;
}

/** Determine size of warmstart information */
static
SCIP_DECL_NLPIGETWARMSTARTSIZE( nlpiGetWarmstatSizeIpopt )
{
   SCIPerrorMessage("method of Ipopt nonlinear solver is not implemented\n");
   SCIPABORT();
   return SCIP_OKAY;
}

/* Write warmstart information to buffer */
static
SCIP_DECL_NLPIGETWARMSTARTMEMO( nlpiGetWarmstatMemoIpopt )
{
   SCIPerrorMessage("method of Ipopt nonlinear solver is not implemented\n");
   SCIPABORT();
   return SCIP_OKAY;
}

/** Write back warmstartmemo to solver */
static
SCIP_DECL_NLPISETWARMSTARTMEMO( nlpiSetWarmstatMemoIpopt )
{
   SCIPerrorMessage("method of Ipopt nonlinear solver is not implemented\n");
   SCIPABORT();
   return SCIP_OKAY;
}

/** get pointer to solver interface */
static
SCIP_DECL_NLPIGETSOLVERPOINTER( nlpiGetSolverPointerIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);

   return GetRawPtr(data->ipopt);
}

/** gets integer parameter of NLP solver */
static
SCIP_DECL_NLPIGETINTPAR( nlpiGetIntparIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(ival != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(IsValid(data->ipopt));

   //@TODO try-catch block for Ipopt exceptions
   switch (type)
   {
      case SCIP_NLPPAR_FROMSCRATCH:
      {
         *ival = 1;
         break;
      }
         
      case SCIP_NLPPAR_VERBLEVEL:
      {
         int printlevel;
         data->ipopt->Options()->GetIntegerValue("print_level", printlevel, "");
         if (printlevel <= J_STRONGWARNING)
            *ival = 0;
         else if (printlevel >= J_DETAILED)
            *ival = 2;
         else /* J_SUMMARY or J_WARNING or J_ITERSUMMARY */
            *ival = 1;
         break;
      }
      
      case SCIP_NLPPAR_FEASTOL:
      {
         SCIPerrorMessage("feasibility tolerance parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_RELOBJTOL:
      {
         SCIPerrorMessage("relative objective tolerance parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
      
      case SCIP_NLPPAR_LOBJLIM:
      {
         SCIPerrorMessage("objective limit parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
      
      case SCIP_NLPPAR_INFINITY:
      {
         SCIPerrorMessage("infinity parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
      
      case SCIP_NLPPAR_ITLIM:
      {
         data->ipopt->Options()->GetIntegerValue("max_iter", *ival, "");
         break;
      }

      case SCIP_NLPPAR_TILIM:
      {
         SCIPerrorMessage("time limit parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
      
      default:
      {
         SCIPerrorMessage("Parameter %d not known to Ipopt interface.\n", type);
         return SCIP_PARAMETERUNKNOWN;
      }
   }

   return SCIP_OKAY;
}

/** sets integer parameter of NLP solver */
static
SCIP_DECL_NLPISETINTPAR( nlpiSetIntparIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(IsValid(data->ipopt));

   switch (type)
   {
      case SCIP_NLPPAR_FROMSCRATCH:
      {
         if (ival == 0 || ival == 1)
         {
            SCIPwarningMessage("from scratch parameter not supported by Ipopt interface yet. Ignored.\n");
         }
         else
         {
            SCIPerrorMessage("Value %d for parameter from scratch out of range {0, 1}\n", ival);
            return SCIP_PARAMETERWRONGVAL;
         }
         break;
      }
         
      case SCIP_NLPPAR_VERBLEVEL:
      {
         switch (ival)
         {
            case 0:
               data->ipopt->Options()->SetIntegerValue("print_level", J_STRONGWARNING);
               break;
            case 1:
               data->ipopt->Options()->SetIntegerValue("print_level", J_ITERSUMMARY);
               break;
            case 2:
               data->ipopt->Options()->SetIntegerValue("print_level", J_DETAILED);
               break;
            default:
               SCIPerrorMessage("Value %d for parameter from verbosity level out of range {0, 1, 2}\n", ival);
               return SCIP_PARAMETERWRONGVAL;
         }
         break;
      }
      
      case SCIP_NLPPAR_FEASTOL:
      {
         SCIPerrorMessage("feasibility tolerance parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_RELOBJTOL:
      {
         SCIPerrorMessage("relative objective tolerance parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
      
      case SCIP_NLPPAR_LOBJLIM:
      {
         SCIPerrorMessage("objective limit parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
      
      case SCIP_NLPPAR_INFINITY:
      {
         SCIPerrorMessage("infinity parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
      
      case SCIP_NLPPAR_ITLIM:
      {
         if (ival >= 0)
         {
            data->ipopt->Options()->SetIntegerValue("max_iter", ival);
         }
         else
         {
            SCIPerrorMessage("Value %d for parameter iteration limit is negative\n", ival);
            return SCIP_PARAMETERWRONGVAL;
         }
         break;
      }

      case SCIP_NLPPAR_TILIM:
      {
         SCIPerrorMessage("time limit parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
      
      default:
      {
         SCIPerrorMessage("Parameter %d not known to Ipopt interface.\n", type);
         return SCIP_PARAMETERUNKNOWN;
      }
   }

   return SCIP_OKAY;
}

/** gets floating point parameter of NLP */
static
SCIP_DECL_NLPIGETREALPAR( nlpiGetRealParIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
   assert(dval != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(IsValid(data->ipopt));

   switch (type)
   {
      case SCIP_NLPPAR_FROMSCRATCH:
      {
         SCIPerrorMessage("from scratch parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
         
      case SCIP_NLPPAR_VERBLEVEL:
      {
         SCIPerrorMessage("verbosity level parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
      
      case SCIP_NLPPAR_FEASTOL:
      {
         data->ipopt->Options()->GetNumericValue("constr_viol_tol", *dval, "");
         break;
      }

      case SCIP_NLPPAR_RELOBJTOL:
      {
         data->ipopt->Options()->GetNumericValue("dual_inf_tol", *dval, "");
         break;
      }
      
      case SCIP_NLPPAR_LOBJLIM:
      {
         *dval = -SCIPinfinity(scip);
         break;
      }
      
      case SCIP_NLPPAR_INFINITY:
      {
         *dval = SCIPinfinity(scip);
         break;
      }
      
      case SCIP_NLPPAR_ITLIM:
      {
         SCIPerrorMessage("iteration limit parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_TILIM:
      {
         data->ipopt->Options()->GetNumericValue("max_cpu_time", *dval, "");
         break;
      }
      
      default:
      {
         SCIPerrorMessage("Parameter %d not known to Ipopt interface.\n", type);
         return SCIP_PARAMETERUNKNOWN;
      }
   }

   return SCIP_OKAY;
}
 
/** sets floating point parameter of NLP solver*/
static
SCIP_DECL_NLPISETREALPAR( nlpiSetRealparIpopt )
{
   assert(scip != NULL);
   assert(nlpi != NULL);
    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(IsValid(data->ipopt));

   switch (type)
   {
      case SCIP_NLPPAR_FROMSCRATCH:
      {
         SCIPerrorMessage("from scratch parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
         
      case SCIP_NLPPAR_VERBLEVEL:
      {
         SCIPerrorMessage("verbosity level parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
      
      case SCIP_NLPPAR_FEASTOL:
      {
         if (dval >= 0)
         {
            data->ipopt->Options()->SetNumericValue("constr_viol_tol", dval);
            /* Let's think that when the user wants to set the feas. tolerance below the ipopt default of the bound_relax_factor,
             * then (s)he has problem to have SCIP accept a solution found by Ipopt.
             * Thus, we turn off the bound_relax_factor completely.
             */ 
            data->ipopt->Options()->SetNumericValue("bound_relax_factor", dval < 1e-8 ? 0. : dval/2.);
         }
         else
         {
            SCIPerrorMessage("Value %g for parameter feasibility tolerance is negative\n", dval);
            return SCIP_PARAMETERWRONGVAL;
         }
         break;
      }

      case SCIP_NLPPAR_RELOBJTOL:
      {
         if (dval >= 0)
         {
            data->ipopt->Options()->SetNumericValue("dual_inf_tol", dval);
         }
         else
         {
            SCIPerrorMessage("Value %g for parameter relative objective tolerance is negative\n", dval);
            return SCIP_PARAMETERWRONGVAL;
         }
         break;
      }
      
      case SCIP_NLPPAR_LOBJLIM:
      {
         SCIPwarningMessage("Parameter lower objective limit not supported by Ipopt interface yet. Ignored.\n");
         break;
      }
      
      case SCIP_NLPPAR_INFINITY:
      {
         data->ipopt->Options()->SetNumericValue("diverging_iterates_tol", dval);
         break;
      }
      
      case SCIP_NLPPAR_ITLIM:
      {
         SCIPerrorMessage("iteration limit parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_TILIM:
      {
         if (dval >= 0)
         {
            data->ipopt->Options()->SetNumericValue("max_cpu_time", dval);
         }
         else
         {
            SCIPerrorMessage("Value %g for parameter time limit is negative\n", dval);
            return SCIP_PARAMETERWRONGVAL;
         }
         break;
      }
      
      default:
      {
         SCIPerrorMessage("Parameter %d not known to Ipopt interface.\n", type);
         return SCIP_PARAMETERUNKNOWN;
      }
   }

   return SCIP_OKAY;
}


/** create solver interface for Ipopt solver */
SCIP_RETCODE
SCIPcreateNlpSolverIpopt(
   SCIP*       scip,  /**< central scip datastructure */
   SCIP_NLPI** nlpi   /**< pointer to buffer for nlpi address */
)
{
   assert( scip && nlpi );
   
   SCIP_NlpiData* nlpidata;
   try
   {
      nlpidata = new SCIP_NlpiData();
      if (nlpidata == NULL)
         throw std::bad_alloc();
   }
   catch (std::bad_alloc)
   {
      SCIPerrorMessage("Not enough memory to allocate SCIP_NlpiData structure in Ipopt interface.\n");
      return SCIP_NOMEMORY;
   }
  
   SCIP_CALL( SCIPnlpiCreate( scip, nlpi,
      NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
      nlpiInitIpopt, nlpiAddVarsIpopt, nlpiAddConstraintsIpopt, nlpiSetObjectiveIpopt, 
      nlpiChgVarBoundsIpopt, nlpiChgConsBoundsIpopt, nlpiDelVarSetIpopt, nlpiDelConstraintSetIpopt,
      nlpiChgLinearCoefsIpopt, nlpiChgQuadraticCoefsIpopt,
      nlpiChgNonlinCoefIpopt, nlpiSetInitialGuessIpopt,
      nlpiSolveIpopt, nlpiGetSolstatIpopt, nlpiGetSoltermIpopt, nlpiGetSolutionIpopt, nlpiGetStatisticsIpopt,
      nlpiGetWarmstatSizeIpopt, nlpiGetWarmstatMemoIpopt,
      nlpiSetWarmstatMemoIpopt, nlpiGetSolverPointerIpopt, nlpiGetIntparIpopt,
      nlpiSetIntparIpopt, nlpiGetRealParIpopt,
      nlpiSetRealparIpopt, nlpiFreeIpopt, nlpidata ) );

   return SCIP_OKAY;
}


/** Method to return some info about the nlp */
bool ScipNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)
{
   assert(scip != NULL);
   assert(nlpidata != NULL);
   assert(nlpidata->oracle != NULL);
   
   n = SCIPnlpiOracleGetNVars(nlpidata->oracle);
   m = SCIPnlpiOracleGetNConstraints(nlpidata->oracle);
   
   const int* offset;
   SCIP_RETCODE retcode = SCIPnlpiOracleGetJacobianSparsity(scip, nlpidata->oracle, &offset, NULL);
   if (retcode != SCIP_OKAY)
      return false;
   assert(offset != NULL);
   nnz_jac_g = offset[m];

   retcode = SCIPnlpiOracleGetHessianLagSparsity(scip, nlpidata->oracle, &offset, NULL);
   if (retcode != SCIP_OKAY)
      return false;
   assert(offset != NULL);
   nnz_h_lag = offset[n];
   
   index_style = TNLP::C_STYLE;
   
   return true;
}

/** Method to return the bounds for my problem */
bool ScipNLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
   assert(scip != NULL);
   assert(nlpidata != NULL);
   assert(nlpidata->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpidata->oracle));
   assert(m == SCIPnlpiOracleGetNConstraints(nlpidata->oracle));
   
   assert(SCIPnlpiOracleGetVarLb(nlpidata->oracle) != NULL);
   assert(SCIPnlpiOracleGetVarUb(nlpidata->oracle) != NULL);
   assert(SCIPnlpiOracleGetConstraintsLhs(nlpidata->oracle) != NULL);
   assert(SCIPnlpiOracleGetConstraintsRhs(nlpidata->oracle) != NULL);
   
   memcpy(x_l, SCIPnlpiOracleGetVarLb(nlpidata->oracle), n * sizeof(SCIP_Real));
   memcpy(x_u, SCIPnlpiOracleGetVarUb(nlpidata->oracle), n * sizeof(SCIP_Real));
   memcpy(g_l, SCIPnlpiOracleGetConstraintsLhs(nlpidata->oracle), m * sizeof(SCIP_Real));
   memcpy(g_u, SCIPnlpiOracleGetConstraintsRhs(nlpidata->oracle), m * sizeof(SCIP_Real));

   return true;
}

/** Method to return the starting point for the algorithm */
bool ScipNLP::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda)
{
   assert(scip != NULL);
   assert(nlpidata != NULL);
   assert(nlpidata->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpidata->oracle));
   assert(m == SCIPnlpiOracleGetNConstraints(nlpidata->oracle));

   if (init_x)
   {
      if (nlpidata->initguess)
         memcpy(x, nlpidata->initguess, n * sizeof(SCIP_Real));
      else
      {
         SCIPwarningMessage("Ipopt started without intial primal values.\n");
         return false; // do not have initial guess, this will make Ipopt fail; @TODO: should we make up some point?
      }
   }
   if (init_z || init_lambda)
      return false;

   return true;
}

/** Method to return the variables linearity. */
bool ScipNLP::get_variables_linearity(Index n, LinearityType* var_types)
{
   assert(scip != NULL);
   assert(nlpidata != NULL);
   assert(nlpidata->oracle != NULL);

   assert(n == SCIPnlpiOracleGetNVars(nlpidata->oracle));
   
   for (int i = 0; i < n; ++i)
   {
      var_types[i] = SCIPnlpiOracleGetVarDegree(nlpidata->oracle, i) <= 1 ? LINEAR : NON_LINEAR;
//      printf("var %d is %s\n", i, var_types[i] == LINEAR ? "linear" : "nonlinear");
   }
   
   return true;
}

/** Method to return the constraint linearity. */
bool ScipNLP::get_constraints_linearity(Index m, LinearityType* const_types)
{
   assert(scip != NULL);
   assert(nlpidata != NULL);
   assert(nlpidata->oracle != NULL);

   assert(m == SCIPnlpiOracleGetNConstraints(nlpidata->oracle));
   
   for (int i = 0; i < m; ++i)
   {
      const_types[i] = SCIPnlpiOracleGetConstraintDegree(nlpidata->oracle, i) <= 1 ? LINEAR : NON_LINEAR;
//      printf("con %d is %s\n", i, const_types[i] == LINEAR ? "linear" : "nonlinear");
   }
   
   return true;
}

/** Method to return the objective value */
bool ScipNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
   assert(scip != NULL);
   assert(nlpidata != NULL);
   assert(nlpidata->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpidata->oracle));
   
   SCIP_RETCODE retcode = SCIPnlpiOracleEvalObjectiveValue(scip, nlpidata->oracle, x, &obj_value);

   return retcode == SCIP_OKAY ? true : false;
}

/** Method to return the gradient of the objective */
bool ScipNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
   assert(scip != NULL);
   assert(nlpidata != NULL);
   assert(nlpidata->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpidata->oracle));
   
   SCIP_Real dummy;
   SCIP_RETCODE retcode = SCIPnlpiOracleEvalObjectiveGradient(scip, nlpidata->oracle, x, new_x, &dummy, grad_f);

   return retcode == SCIP_OKAY ? true : false;
}

/** Method to return the constraint residuals */
bool ScipNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
   assert(scip != NULL);
   assert(nlpidata != NULL);
   assert(nlpidata->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpidata->oracle));
   
   SCIP_RETCODE retcode = SCIPnlpiOracleEvalConstraintValues(scip, nlpidata->oracle, x, g);

   return retcode == SCIP_OKAY ? true : false;
}

/** Method to return:
 *   1) The structure of the jacobian (if "values" is NULL)
 *   2) The values of the jacobian (if "values" is not NULL)
 */
bool ScipNLP::eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values)
{
   assert(scip != NULL);
   assert(nlpidata != NULL);
   assert(nlpidata->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpidata->oracle));
   assert(m == SCIPnlpiOracleGetNConstraints(nlpidata->oracle));
   
   if (values == NULL)
   { /* Ipopt wants to know sparsity structure */
      assert(iRow != NULL);
      assert(jCol != NULL);
      
      const int* jacoffset;
      const int* jaccol;
      SCIP_RETCODE retcode = SCIPnlpiOracleGetJacobianSparsity(scip, nlpidata->oracle, &jacoffset, &jaccol);
      if (retcode != SCIP_OKAY)
         return false;
      
      assert(jacoffset[0] == 0);
      assert(jacoffset[m] == nele_jac);
      int j = jacoffset[0];
      for (int i = 0; i < m; ++i)
         for (; j < jacoffset[i+1]; ++j)
            iRow[j] = i;
      
      memcpy(jCol, jaccol, nele_jac * sizeof(int));
   }
   else
   {
      SCIP_RETCODE retcode = SCIPnlpiOracleEvalJacobian(scip, nlpidata->oracle, x, new_x, NULL, values);
      if (retcode != SCIP_OKAY)
         return false;
   }

   return true;
}

/** Method to return:
 *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
 *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
 */
bool ScipNLP::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values)
{
   assert(scip != NULL);
   assert(nlpidata != NULL);
   assert(nlpidata->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpidata->oracle));
   assert(m == SCIPnlpiOracleGetNConstraints(nlpidata->oracle));
   
   if (values == NULL)
   { /* Ipopt wants to know sparsity structure */
      assert(iRow != NULL);
      assert(jCol != NULL);
     
      const int* heslagoffset;
      const int* heslagcol;
      SCIP_RETCODE retcode = SCIPnlpiOracleGetHessianLagSparsity(scip, nlpidata->oracle, &heslagoffset, &heslagcol);
      if (retcode != SCIP_OKAY)
         return false;
     
      assert(heslagoffset[0] == 0);
      assert(heslagoffset[n] == nele_hess);
      int j = heslagoffset[0];
      for (int i = 0; i < n; ++i)
         for (; j < heslagoffset[i+1]; ++j)
            iRow[j] = i;
     
      memcpy(jCol, heslagcol, nele_hess * sizeof(int));
   }
   else
   {
      SCIP_RETCODE retcode = SCIPnlpiOracleEvalHessianLag(scip, nlpidata->oracle, x, new_x, obj_factor, lambda, values);
      if (retcode != SCIP_OKAY)
         return false;
   }
   
   return true;
}

/** Method called by the solver at each iteration.
 * Checks whether ^C was hit.
 */
bool ScipNLP::intermediate_callback(AlgorithmMode mode, Index iter, Number obj_value, Number inf_pr, Number inf_du, Number mu, Number d_norm, Number regularization_size, Number alpha_du, Number alpha_pr, Index ls_trials, const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq)
{
   assert(scip != NULL);
   
   return !SCIPpressedCtrlC(scip);
}

/** This method is called when the algorithm is complete so the TNLP can store/write the solution.
 */
void ScipNLP::finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g, const Number* lambda, Number obj_value, const IpoptData* data, IpoptCalculatedQuantities* cq)
{
   assert(scip != NULL);
   assert(nlpidata != NULL);
   assert(nlpidata->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpidata->oracle));
   assert(m == SCIPnlpiOracleGetNConstraints(nlpidata->oracle));
   
   bool check_feasibility = false; // whether we should check x for feasibility, if not NULL
   switch (status)
   {
      case SUCCESS:
      case STOP_AT_ACCEPTABLE_POINT:
         nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_LOCOPT;
         nlpidata->lasttermstat = SCIP_NLPITERMSTAT_OKAY;
         assert(x != NULL);
         break;
         
      case FEASIBLE_POINT_FOUND:
         nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
         nlpidata->lasttermstat = SCIP_NLPITERMSTAT_OKAY;
         assert(x != NULL);
         break;
         
      case MAXITER_EXCEEDED:
         check_feasibility = true;
         nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
         nlpidata->lasttermstat = SCIP_NLPITERMSTAT_ITLIM;
         break;
         
      case CPUTIME_EXCEEDED:
         check_feasibility = true;
         nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
         nlpidata->lasttermstat = SCIP_NLPITERMSTAT_TILIM;
         break;
         
      case STOP_AT_TINY_STEP:
      case RESTORATION_FAILURE:
      case ERROR_IN_STEP_COMPUTATION:
         check_feasibility = true;
         nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
         nlpidata->lasttermstat = SCIP_NLPITERMSTAT_NUMERR;
         break;

      case LOCAL_INFEASIBILITY:
         nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
         nlpidata->lasttermstat = SCIP_NLPITERMSTAT_OKAY;
         break;

      case DIVERGING_ITERATES:
         nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_UNBOUNDED;
         nlpidata->lasttermstat = SCIP_NLPITERMSTAT_UOBJLIM;
         break;

      case INVALID_NUMBER_DETECTED:
         nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
         nlpidata->lasttermstat = SCIP_NLPITERMSTAT_EVALERR;
         break;

      case USER_REQUESTED_STOP:
         nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
         nlpidata->lasttermstat = SCIP_NLPITERMSTAT_TILIM;
         break;

      case TOO_FEW_DEGREES_OF_FREEDOM:
      case INTERNAL_ERROR:
         nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
         nlpidata->lasttermstat = SCIP_NLPITERMSTAT_OTHER;
         break;

      default:
         SCIPerrorMessage("Ipopt returned with unknown solution status %d\n", status);
         nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
         nlpidata->lasttermstat = SCIP_NLPITERMSTAT_OTHER;
         break;
   }

   if (x)
   {
      if (!nlpidata->lastsol)
      {
         SCIP_RETCODE retcode = SCIPduplicateMemoryArray(scip, &nlpidata->lastsol, x, n);
         if (retcode != SCIP_OKAY)
         {
            nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
            nlpidata->lasttermstat = retcode == SCIP_NOMEMORY ? SCIP_NLPITERMSTAT_MEMERR : SCIP_NLPITERMSTAT_OTHER;
            return;
         }
      }
      else
         memcpy(nlpidata->lastsol, x, n * sizeof(SCIP_Real));
      
      if (check_feasibility && cq)
      {
         Number constrviol = cq->curr_constraint_violation();
         Number constrvioltol;
         nlpidata->ipopt->Options()->GetNumericValue("constr_viol_tol", constrvioltol, "");
         if (constrviol <= constrvioltol)
         {
            nlpidata->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
         }
         else
         {
            nlpidata->ipopt->Options()->GetNumericValue("acceptable_constr_viol_tol", constrvioltol, "");
            if (constrviol <= constrvioltol)
               nlpidata->lastsolstat = SCIP_NLPSOLSTAT_FEASIBLE;
         }
      }
   }
}
