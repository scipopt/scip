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
#pragma ident "@(#) $Id: nlpi_ipopt.cpp,v 1.2 2010/04/21 14:21:14 bzfviger Exp $"

/**@file    nlpi_ipopt.cpp
 * @ingroup NLPIS
 * @brief   Ipopt NLP interface
 * @author  Stefan Vigerske
 *
 * @todo warm starts
 * @todo use new_x: Ipopt sets new_x = false if any function has been evaluated for the current x already, while oracle allows new_x to be false only if the current function has been evaluated for the current x before
 * @todo using SCIPinterrupted() in ScipNLP::intermediate_callback creates a dependency to the scip library, which we probably want to avoid if using the nlpi standalone
 *       also do we really want to stop if *some* SCIP instance was interrupted? sometimes we really want to, since NLPs solver can "run crazy" 
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "nlpi/nlpi_ipopt.h"
#include "nlpi/nlpi.h"
#include "nlpi/nlpioracle.h"
#include "scip/interrupt.h"

#include <new>      /* for std::bad_alloc */

#include "IpIpoptApplication.hpp"
namespace Ipopt
{
   class IpoptNLP;
   class IpoptData;
}
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpSolveStatistics.hpp"
#include "IpJournalist.hpp"

using namespace Ipopt;

#define NLPI_NAME          "Ipopt"           /**< short concise name of solver */
#define NLPI_DESC          "Ipopt interface" /**< description of solver */
#define NLPI_TYPE          "IP"              /**< solver type */
#define NLPI_PRIORITY      0                 /**< priority */

#ifdef SCIP_DEBUG
#define DEFAULT_PRINTLEVEL J_WARNING         /**< default print level of Ipopt */
#else
#define DEFAULT_PRINTLEVEL J_STRONGWARNING   /**< default print level of Ipopt */
#endif
#define DEFAULT_MAXITER    3000              /**< default iteration limit for Ipopt */

class ScipNLP;

struct SCIP_NlpiData
{
   BMS_BLKMEM*                 blkmem;       /**< block memory */
   SCIP_Real                   infinity;     /**< initial value for infinity */
};

struct SCIP_NlpiProblem
{
public:
   SCIP_NLPIORACLE*            oracle;       /**< Oracle-helper to store and evaluate NLP */
   
   SmartPtr<IpoptApplication>  ipopt;        /**< Ipopt application */
   SmartPtr<ScipNLP>           nlp;          /**< NLP in Ipopt form */
   std::string                 optfile;      /**< name of options file */
   
   SCIP_Bool                   firstrun;     /**< whether the next NLP solve will be the first one (with the current problem structure) */
   SCIP_Real*                  initguess;    /**< initial values for primal variables, or NULL if not known */
   
   SCIP_NLPSOLSTAT             lastsolstat;  /**< solution status from last run */
   SCIP_NLPTERMSTAT            lasttermstat; /**< termination status from last run */
   SCIP_Real*                  lastsol;      /**< solution from last run, if available */
   int                         lastniter;    /**< number of iterations in last run */
   SCIP_Real                   lasttime;     /**< time spend in last run */

   SCIP_NlpiProblem()
   : oracle(NULL),
     firstrun(TRUE), initguess(NULL),
     lastsolstat(SCIP_NLPSOLSTAT_UNKNOWN), lasttermstat(SCIP_NLPTERMSTAT_OTHER), lastsol(NULL),
     lastniter(-1), lasttime(-1.0)
   { }
};

/** TNLP implementation for SCIPs NLP */
class ScipNLP : public TNLP
{
private:
   SCIP_NLPIPROBLEM*     nlpiproblem;        /**< NLPI problem data */

public:
   /** constructor */
   ScipNLP(
      SCIP_NLPIPROBLEM*  nlpiproblem_ = NULL /**< NLPI problem data */
   )
   : nlpiproblem(nlpiproblem_)
   { }

   /** destructor */
   ~ScipNLP() { }
   
   /** sets NLPI data structure */
   void setNLPIPROBLEM(SCIP_NLPIPROBLEM* nlpiproblem_)
   {
      assert(nlpiproblem_ != NULL);
      nlpiproblem = nlpiproblem_;
   }

   /** Method to return some info about the nlp */
   bool get_nlp_info(
      Index&             n,                  /**< place to store number of variables */ 
      Index&             m,                  /**< place to store number of constraints */ 
      Index&             nnz_jac_g,          /**< place to store number of nonzeros in jacobian */
      Index&             nnz_h_lag,          /**< place to store number of nonzeros in hessian */
      IndexStyleEnum&    index_style         /**< place to store used index style (0-based or 1-based) */
   );

   /** Method to return the bounds for my problem */
   bool get_bounds_info(
      Index              n,                  /**< number of variables */ 
      Number*            x_l,                /**< buffer to store lower bounds on variables */
      Number*            x_u,                /**< buffer to store upper bounds on variables */
      Index              m,                  /**< number of constraints */
      Number*            g_l,                /**< buffer to store lower bounds on constraints */
      Number*            g_u                 /**< buffer to store lower bounds on constraints */
   );

   /** Method to return the starting point for the algorithm */
   bool get_starting_point(
      Index              n,                  /**< number of variables */ 
      bool               init_x,             /**< whether initial values for primal values are requested */ 
      Number*            x,                  /**< buffer to store initial primal values */
      bool               init_z,             /**< whether initial values for dual values of variable bounds are requested */  
      Number*            z_L,                /**< buffer to store dual values for variable lower bounds */
      Number*            z_U,                /**< buffer to store dual values for variable upper bounds */
      Index              m,                  /**< number of constraints */
      bool               init_lambda,        /**< whether initial values for dual values of constraints are required */
      Number*            lambda              /**< buffer to store dual values of constraints */
   );

   /** Method to return the variables linearity. */
   bool get_variables_linearity(
      Index              n,                  /**< number of variables */ 
      LinearityType*     var_types           /**< buffer to store linearity types of variables */
   );

   /** Method to return the constraint linearity. */
   bool get_constraints_linearity(
      Index              m,                  /**< number of constraints */
      LinearityType*     const_types         /**< buffer to store linearity types of constraints */
   );

   /** Method to return the objective value */
   bool eval_f(
      Index              n,                  /**< number of variables */ 
      const Number*      x,                  /**< point to evaluate */ 
      bool               new_x,              /**< whether some function evaluation method has been called for this point before */
      Number&            obj_value           /**< place to store objective function value */
   );

   /** Method to return the gradient of the objective */
   bool eval_grad_f(
      Index              n,                  /**< number of variables */ 
      const Number*      x,                  /**< point to evaluate */ 
      bool               new_x,              /**< whether some function evaluation method has been called for this point before */
      Number*            grad_f              /**< buffer to store objective gradient */
   );

   /** Method to return the constraint residuals */
   bool eval_g(
      Index              n,                  /**< number of variables */ 
      const Number*      x,                  /**< point to evaluate */ 
      bool               new_x,              /**< whether some function evaluation method has been called for this point before */
      Index              m,                  /**< number of constraints */
      Number*            g                   /**< buffer to store constraint function values */
   );

   /** Method to return:
    *   1) The structure of the jacobian (if "values" is NULL)
    *   2) The values of the jacobian (if "values" is not NULL)
    */
   bool eval_jac_g(
      Index              n,                  /**< number of variables */ 
      const Number*      x,                  /**< point to evaluate */ 
      bool               new_x,              /**< whether some function evaluation method has been called for this point before */
      Index              m,                  /**< number of constraints */
      Index              nele_jac,           /**< number of nonzero entries in jacobian */ 
      Index*             iRow,               /**< buffer to store row indices of nonzero jacobian entries, or NULL if values 
                                              * are requested */
      Index*             jCol,               /**< buffer to store column indices of nonzero jacobian entries, or NULL if values
                                              * are requested */                  
      Number*            values              /**< buffer to store values of nonzero jacobian entries, or NULL if structure is
                                              * requested */
   );

   /** Method to return:
    *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
    *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
    */
   bool eval_h(
      Index              n,                  /**< number of variables */ 
      const Number*      x,                  /**< point to evaluate */ 
      bool               new_x,              /**< whether some function evaluation method has been called for this point before */
      Number             obj_factor,         /**< weight for objective function */ 
      Index              m,                  /**< number of constraints */
      const Number*      lambda,             /**< weights for constraint functions */ 
      bool               new_lambda,         /**< whether the hessian has been evaluated for these values of lambda before */
      Index              nele_hess,          /**< number of nonzero entries in hessian */
      Index*             iRow,               /**< buffer to store row indices of nonzero hessian entries, or NULL if values
                                              * are requested */
      Index*             jCol,               /**< buffer to store column indices of nonzero hessian entries, or NULL if values
                                              * are requested */                  
      Number*            values              /**< buffer to store values of nonzero hessian entries, or NULL if structure is requested */
   );
   
   /** Method called by the solver at each iteration.
    * 
    * Checks whether Ctrl-C was hit.
    */
   bool intermediate_callback(
      AlgorithmMode      mode,               /**< current mode of algorithm */
      Index              iter,               /**< current iteration number */
      Number             obj_value,          /**< current objective value */
      Number             inf_pr,             /**< current primal infeasibility */
      Number             inf_du,             /**< current dual infeasibility */
      Number             mu,                 /**< current barrier parameter */
      Number             d_norm,             /**< current gradient norm */
      Number             regularization_size,/**< current size of regularization */
      Number             alpha_du,           /**< current dual alpha */
      Number             alpha_pr,           /**< current primal alpha */
      Index              ls_trials,          /**< current number of linesearch trials */
      const IpoptData*   ip_data,            /**< pointer to Ipopt Data */
      IpoptCalculatedQuantities* ip_cq       /**< pointer to current calculated quantities */
   );

   /** This method is called when the algorithm is complete so the TNLP can store/write the solution. */
   void finalize_solution(
      SolverReturn       status,             /**< solve and solution status */ 
      Index              n,                  /**< number of variables */ 
      const Number*      x,                  /**< primal solution values */ 
      const Number*      z_L,                /**< dual values of variable lower bounds */
      const Number*      z_U,                /**< dual values of variable upper bounds */
      Index              m,                  /**< number of constraints */ 
      const Number*      g,                  /**< values of constraints */ 
      const Number*      lambda,             /**< dual values of constraints */ 
      Number             obj_value,          /**< objective function value */ 
      const IpoptData*   data,               /**< pointer to Ipopt Data */ 
      IpoptCalculatedQuantities* cq          /**< pointer to calculated quantities */
   );
};

/** A particular Ipopt::Journal implementation that uses the SCIP message routines for output.
 */
class ScipJournal : public Ipopt::Journal {
public:
  ScipJournal(const char* name, Ipopt::EJournalLevel default_level)
  : Ipopt::Journal(name, default_level)
  { }

  ~ScipJournal() { }

protected:
  void PrintImpl(Ipopt::EJournalCategory category, Ipopt::EJournalLevel level, const char* str)
  {
     SCIPmessagePrintInfo(str);
  }

  void PrintfImpl(Ipopt::EJournalCategory category, Ipopt::EJournalLevel level, const char* pformat, va_list ap)
  {
     SCIPmessageVPrintInfo(pformat, ap);
  }

  void FlushBufferImpl() { }
};

/** clears the last solution arrays and sets the solstat and termstat to unknown and other, resp. */
static
void SCIPnlpiIpoptInvalidateSolution(
   SCIP_NLPIPROBLEM*     problem             /**< data structure of problem */
)
{
   assert(problem != NULL);
   
   BMSfreeMemoryArrayNull(&problem->lastsol);
   problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   problem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
}

/** destructor of NLP interface to free nlpi data
 * 
 * input:
 *  - nlpi datastructure for solver interface
 */
static
SCIP_DECL_NLPIFREE(nlpiFreeIpopt)
{
   SCIP_NLPIDATA* data;

   assert(nlpi != NULL);
   
   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);
   
   BMSfreeMemory(&data);
   
   assert(data == NULL);
   
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
SCIP_DECL_NLPIGETSOLVERPOINTER(nlpiGetSolverPointerIpopt)
{
   assert(nlpi != NULL);

   return NULL;
}

/** creates a problem instance
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem pointer to store the problem data
 *  - name name of problem, can be NULL
 */
static
SCIP_DECL_NLPICREATEPROBLEM(nlpiCreateProblemIpopt)
{
   SCIP_NLPIDATA* data;

   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);
   
   *problem = new SCIP_NLPIPROBLEM;
   if( *problem == NULL )
      return SCIP_NOMEMORY;

   SCIP_CALL( SCIPnlpiOracleCreate(data->blkmem, &(*problem)->oracle) );
   SCIP_CALL( SCIPnlpiOracleSetInfinity((*problem)->oracle, data->infinity) );
   
   try
   {
      /* initialize IPOPT without default journal */
      (*problem)->ipopt = new IpoptApplication(false);
      if( IsNull((*problem)->ipopt) )
         throw std::bad_alloc();
      
      /* plugin our journal to get output through SCIP message handler */
      SmartPtr<Journal> jrnl = new ScipJournal("console", J_ITERSUMMARY);
      if( IsNull(jrnl) )
         throw std::bad_alloc();
      jrnl->SetPrintLevel(J_DBG, J_NONE);
      if( !(*problem)->ipopt->Jnlst()->AddJournal(jrnl) )
      {
         SCIPerrorMessage("Failed to register ScipJournal for IPOPT output.");
      }

      /* initialize Ipopt/SCIP NLP interface */
      (*problem)->nlp = new ScipNLP(*problem);
      if( IsNull((*problem)->nlp) )
         throw std::bad_alloc();
   }
   catch( std::bad_alloc )
   {
      SCIPerrorMessage("Not enough memory to initialize Ipopt.\n");
      return SCIP_NOMEMORY;
   }
   
   (*problem)->ipopt->Options()->SetIntegerValue("print_level", DEFAULT_PRINTLEVEL);
   /* (*problem)->ipopt->Options()->SetStringValue("print_timing_statistics", "yes"); */
   (*problem)->ipopt->Options()->SetStringValue("mu_strategy", "adaptive");
   (*problem)->ipopt->Options()->SetStringValue("expect_infeasible_problem", "yes");
   /* it seem to be better to let Ipopt relax bounds a bit to ensure that a relative interior exists;
    * however, if we relax the bounds too much, then the solutions tend to be slightly infeasible */
   (*problem)->ipopt->Options()->SetNumericValue("tol", SCIP_DEFAULT_FEASTOL/2);
   (*problem)->ipopt->Options()->SetNumericValue("bound_relax_factor", SCIP_DEFAULT_FEASTOL/2);
   (*problem)->ipopt->Options()->SetNumericValue("constr_viol_tol", 0.75*SCIP_DEFAULT_FEASTOL);
   (*problem)->ipopt->Options()->SetIntegerValue("max_iter", DEFAULT_MAXITER);
   (*problem)->ipopt->Options()->SetNumericValue("nlp_lower_bound_inf", -data->infinity, false);
   (*problem)->ipopt->Options()->SetNumericValue("nlp_upper_bound_inf",  data->infinity, false);
   (*problem)->ipopt->Options()->SetNumericValue("diverging_iterates_tol", data->infinity, false);
   /* (*problem)->ipopt->Options()->SetStringValue("dependency_detector", "ma28"); */
   /* (*problem)->ipopt->Options()->SetStringValue("hessian_approximation", "limited-memory"); */
#ifdef SCIP_DEBUG
   (*problem)->ipopt->Options()->SetStringValue("derivative_test", "second-order");
#endif

   if( (*problem)->ipopt->Initialize((*problem)->optfile) != Solve_Succeeded )
   {
      SCIPerrorMessage("Error during initialization of Ipopt using optionfile \"%s\"\n", (*problem)->optfile.c_str());
      return SCIP_ERROR;
   }

   /* TODO store the problem name somewhere and use in print */
   
   return SCIP_OKAY;
}

/** free a problem instance
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem pointer where problem data is stored 
 */
static
SCIP_DECL_NLPIFREEPROBLEM(nlpiFreeProblemIpopt)
{
   assert(nlpi     != NULL);
   assert(problem  != NULL);
   assert(*problem != NULL);
  
   if( (*problem)->oracle != NULL )
   {
      SCIP_CALL( SCIPnlpiOracleFree(&(*problem)->oracle) );
   }
   
   BMSfreeMemoryArrayNull(&(*problem)->initguess);
   BMSfreeMemoryArrayNull(&(*problem)->lastsol);
   
   delete *problem;
   *problem = NULL;

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
SCIP_DECL_NLPIGETPROBLEMPOINTER(nlpiGetProblemPointerIpopt)
{
   assert(nlpi    != NULL);
   assert(problem != NULL);
   
   return GetRawPtr(problem->nlp);
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
SCIP_DECL_NLPIADDVARS(nlpiAddVarsIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleAddVars(problem->oracle, nvars, lbs, ubs, varnames) );
   
   problem->firstrun = TRUE;
   BMSfreeMemoryArrayNull(&problem->initguess);
   SCIPnlpiIpoptInvalidateSolution(problem);
   
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
 *  - linoffsets start index of each constraints linear coefficients in lininds and linvals
 *    length: ncons + 1, linoffsets[ncons] gives length of lininds and linvals
 *    may be NULL in case of no linear part
 *  - lininds variable indices
 *    may be NULL in case of no linear part
 *  - linvals coefficient values
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
SCIP_DECL_NLPIADDCONSTRAINTS(nlpiAddConstraintsIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);
 
   SCIP_CALL( SCIPnlpiOracleAddConstraints(problem->oracle,
      ncons, lhss, rhss,
      nlininds, lininds, linvals,
      nquadrows, quadrowidxs, quadoffsets, quadinds, quadvals,
      exprvaridxs, exprtrees, names) );

   problem->firstrun = TRUE;
   SCIPnlpiIpoptInvalidateSolution(problem);

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
SCIP_DECL_NLPISETOBJECTIVE(nlpiSetObjectiveIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleSetObjective(problem->oracle,
      constant, nlins, lininds, linvals,
      nquadcols, quadcols, quadoffsets, quadinds, quadvals,
      exprvaridxs, exprtree) );

   problem->firstrun = TRUE;
   SCIPnlpiIpoptInvalidateSolution(problem);

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
SCIP_DECL_NLPICHGVARBOUNDS(nlpiChgVarBoundsIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);
    
   SCIP_CALL( SCIPnlpiOracleChgVarBounds(problem->oracle, nvars, indices, lbs, ubs) );

   SCIPnlpiIpoptInvalidateSolution(problem);

   return SCIP_OKAY;
}

/** change constraint bounds
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - ncons number of constraints to change bounds
 *  - indices indices of constraints to change bounds
 *  - lbs new lower bounds
 *  - ubs new upper bounds
 */
static
SCIP_DECL_NLPICHGCONSBOUNDS(nlpiChgConsBoundsIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);
      
   SCIP_CALL( SCIPnlpiOracleChgConsBounds(problem->oracle, ncons, indices, lbs, ubs) );
   
   SCIPnlpiIpoptInvalidateSolution(problem);

   return SCIP_OKAY;
}

/** delete a set of variables
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nlpi datastructure for solver interface
 *  - dstats deletion status of vars; 1 if var should be deleted, 0 if not
 * 
 * output:
 *  - dstats new position of var, -1 if var was deleted
 */
static
SCIP_DECL_NLPIDELVARSET(nlpiDelVarSetIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleDelVarSet(problem->oracle, dstats) );

   problem->firstrun = TRUE;
   BMSfreeMemoryArrayNull(&problem->initguess); // @TODO keep initguess for remaining variables 

   SCIPnlpiIpoptInvalidateSolution(problem);

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
SCIP_DECL_NLPIDELCONSSET(nlpiDelConstraintSetIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleDelConsSet(problem->oracle, dstats) );

   problem->firstrun = TRUE;

   SCIPnlpiIpoptInvalidateSolution(problem);

   return SCIP_OKAY;
}

/** change one linear coefficient in a constraint or objective
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idx index of constraint or -1 for objective
 *  - nvals number of values in linear constraint
 *  - varidxs indices of variable
 *  - vals new values for coefficient
 *
 * return: Error if coefficient did not exist before
 */
static
SCIP_DECL_NLPICHGLINEARCOEFS(nlpiChgLinearCoefsIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleChgLinearCoefs(problem->oracle, idx, nvals, varidxs, vals) );
   SCIPnlpiIpoptInvalidateSolution(problem);

   return SCIP_OKAY;
}

/** change one coefficient in the quadratic part of a constraint or objective
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idx index of constraint or -1 for objective
 *  - nentries number of values in quadratic constraint
 *  - rows row offset containing modified indices
 *  - cols cols containing modified indices to the corresponding row offset
 *  - values coefficients corresponding to same indices as used when constraint/objective was constructed
 *
 * return: Error if coefficient did not exist before
 */
static
SCIP_DECL_NLPICHGQUADCOEFS(nlpiChgQuadraticCoefsIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleChgQuadCoefs(problem->oracle, idx, nentries, rows, cols, values) );
   SCIPnlpiIpoptInvalidateSolution(problem);

   return SCIP_OKAY;
}

/** change the value of one parameter in the nonlinear part
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - considx index of constraint or -1 for objective
 *  - paramidx index of parameter
 *  - value new value for nonlinear parameter
 * 
 * return: Error if parameter does not exist
 */
static
SCIP_DECL_NLPICHGNONLINCOEF(nlpiChgNonlinCoefIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);
/*    
   SCIP_NLPIDATA* data = SCIPnlpiGetNlpiData(nlpi);
   assert(data != NULL);
   assert(data->oracle != NULL);
*/ 
   SCIPerrorMessage("ChgNonlinCoef method of Ipopt nonlinear solver is not implemented\n");

   SCIPnlpiIpoptInvalidateSolution(problem);

   return SCIP_ERROR;
}

/** sets initial guess for primal variables
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - values initial starting solution, or NULL to clear previous starting solution
 */
static
SCIP_DECL_NLPISETINITIALGUESS(nlpiSetInitialGuessIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   if( values != NULL )
   {
      if( !problem->initguess )
      {
         if( BMSduplicateMemoryArray(&problem->initguess, values, SCIPnlpiOracleGetNVars(problem->oracle)) == NULL )
            return SCIP_NOMEMORY;
      }
      else
      {
         BMScopyMemoryArray(problem->initguess, values, SCIPnlpiOracleGetNVars(problem->oracle));
      }
   }
   else
   {
      BMSfreeMemoryArrayNull(&problem->initguess);
   }

   return SCIP_OKAY;
}

/** tries to solve NLP
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 */
static
SCIP_DECL_NLPISOLVE(nlpiSolveIpopt)
{
   ApplicationReturnStatus status;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);
    
#ifdef SCIP_DEBUG
//   SCIP_CALL( SCIPnlpiOraclePrintProblem(problem->oracle, NULL) );
#endif
   
   assert(IsValid(problem->ipopt));
   assert(IsValid(problem->nlp));
   
   problem->nlp->setNLPIPROBLEM(problem);
   
   problem->lastniter = -1;
   problem->lasttime  = -1.0;
   
   try
   {
      SmartPtr<SolveStatistics> stats;

      if( problem->firstrun )
         status = problem->ipopt->OptimizeTNLP(GetRawPtr(problem->nlp));
      else
         status = problem->ipopt->ReOptimizeTNLP(GetRawPtr(problem->nlp));
      
      // catch the very bad status codes
      switch( status ) {
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
            SCIPnlpiIpoptInvalidateSolution(problem);
            problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
            problem->lasttermstat = SCIP_NLPTERMSTAT_EVALERR;
         default: ;
      }

      stats = problem->ipopt->Statistics();
      if( IsValid(stats) )
      {
         problem->lastniter = stats->IterationCount();
         problem->lasttime  = stats->TotalCPUTime();
      }
   }
   catch( IpoptException except )
   {
      SCIPerrorMessage("Ipopt returned with exception: %s\n", except.Message().c_str());
      return SCIP_ERROR;
   }
   
   problem->firstrun = FALSE;

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
SCIP_DECL_NLPIGETSOLSTAT(nlpiGetSolstatIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   return problem->lastsolstat;
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
SCIP_DECL_NLPIGETTERMSTAT(nlpiGetSoltermIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   return problem->lasttermstat;
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
SCIP_DECL_NLPIGETSOLUTION(nlpiGetSolutionIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(primalvalues != NULL);
    
   *primalvalues = problem->lastsol;
   
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
SCIP_DECL_NLPIGETSTATISTICS(nlpiGetStatisticsIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   SCIPnlpStatisticsSetNIterations(statistics, problem->lastniter);
   SCIPnlpStatisticsSetTotalTime  (statistics, problem->lasttime);
   
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
SCIP_DECL_NLPIGETWARMSTARTSIZE(nlpiGetWarmstatSizeIpopt)
{
   SCIPerrorMessage("method of Ipopt nonlinear solver is not implemented\n");
   return SCIP_ERROR;
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
SCIP_DECL_NLPIGETWARMSTARTMEMO(nlpiGetWarmstatMemoIpopt)
{
   SCIPerrorMessage("method of Ipopt nonlinear solver is not implemented\n");
   return SCIP_ERROR;
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
SCIP_DECL_NLPISETWARMSTARTMEMO(nlpiSetWarmstatMemoIpopt)
{
   SCIPerrorMessage("method of Ipopt nonlinear solver is not implemented\n");
   SCIPABORT();
   return SCIP_OKAY;
}

/** gets integer parameter of NLP
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival pointer to store the parameter value
 * 
 * output:
 *  - ival parameter value
 */
static
SCIP_DECL_NLPIGETINTPAR(nlpiGetIntParIpopt)
{
   assert(nlpi != NULL);
   assert(ival != NULL);
   assert(problem != NULL);
   assert(IsValid(problem->ipopt));

   //@TODO try-catch block for Ipopt exceptions
   switch( type )
   {
      case SCIP_NLPPAR_FROMSCRATCH:
      {
         *ival = 1;
         break;
      }
         
      case SCIP_NLPPAR_VERBLEVEL:
      {
         int printlevel;
         problem->ipopt->Options()->GetIntegerValue("print_level", printlevel, "");
         if( printlevel <= J_STRONGWARNING )
            *ival = 0;
         else if( printlevel >= J_DETAILED )
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
         problem->ipopt->Options()->GetIntegerValue("max_iter", *ival, "");
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

/** sets integer parameter of NLP
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival parameter value
 */
static
SCIP_DECL_NLPISETINTPAR(nlpiSetIntParIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(IsValid(problem->ipopt));

   switch( type )
   {
      case SCIP_NLPPAR_FROMSCRATCH:
      {
         if( ival == 0 || ival == 1 )
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
         switch( ival )
         {
            case 0:
               problem->ipopt->Options()->SetIntegerValue("print_level", J_STRONGWARNING);
               break;
            case 1:
               problem->ipopt->Options()->SetIntegerValue("print_level", J_ITERSUMMARY);
               break;
            case 2:
               problem->ipopt->Options()->SetIntegerValue("print_level", J_DETAILED);
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
         if( ival >= 0 )
         {
            problem->ipopt->Options()->SetIntegerValue("max_iter", ival);
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

/** gets floating point parameter of NLP
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval pointer to store the parameter value
 * 
 * output:
 *  - dval parameter value
 */
static
SCIP_DECL_NLPIGETREALPAR(nlpiGetRealParIpopt)
{
   assert(nlpi != NULL);
   assert(dval != NULL);
   assert(type == SCIP_NLPPAR_INFINITY || problem != NULL);
   assert(type == SCIP_NLPPAR_INFINITY || IsValid(problem->ipopt));

   switch( type )
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
         problem->ipopt->Options()->GetNumericValue("constr_viol_tol", *dval, "");
         break;
      }

      case SCIP_NLPPAR_RELOBJTOL:
      {
         problem->ipopt->Options()->GetNumericValue("dual_inf_tol", *dval, "");
         break;
      }
      
      case SCIP_NLPPAR_LOBJLIM:
      {
         *dval = -SCIPnlpiOracleGetInfinity(problem->oracle);
         break;
      }
      
      case SCIP_NLPPAR_INFINITY:
      {
         if( problem )
         {
            *dval = SCIPnlpiOracleGetInfinity(problem->oracle);
         }
         else
         {
            SCIP_NLPIDATA* data = SCIPnlpiGetData(nlpi);
            assert(data != NULL);
            *dval = data->infinity;
         }
         break;
      }
      
      case SCIP_NLPPAR_ITLIM:
      {
         SCIPerrorMessage("iteration limit parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_TILIM:
      {
         problem->ipopt->Options()->GetNumericValue("max_cpu_time", *dval, "");
         break;
      }

      case SCIP_NLPPAR_OPTFILE:
      {
         SCIPerrorMessage("option file parameter is of type string.\n");
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
 
/** sets floating point parameter of NLP
 * 
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval parameter value
 */
static
SCIP_DECL_NLPISETREALPAR(nlpiSetRealParIpopt)
{
   assert(nlpi != NULL);
   assert(type == SCIP_NLPPAR_INFINITY || problem != NULL);
   assert(type == SCIP_NLPPAR_INFINITY || IsValid(problem->ipopt));

   switch( type )
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
         if( dval >= 0 )
         {
            problem->ipopt->Options()->SetNumericValue("constr_viol_tol", dval);
            /* Let's think that when the user wants to set the feas. tolerance below the ipopt default of the bound_relax_factor,
             * then (s)he has problem to have SCIP accept a solution found by Ipopt.
             * Thus, we turn off the bound_relax_factor completely.
             */ 
            problem->ipopt->Options()->SetNumericValue("bound_relax_factor", dval < 1e-8 ? 0. : dval/2.);
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
         if( dval >= 0 )
         {
            problem->ipopt->Options()->SetNumericValue("dual_inf_tol", dval);
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
         if( dval < 0.0 )
            return SCIP_PARAMETERWRONGVAL;
         if( problem )
         {
            problem->ipopt->Options()->SetNumericValue("diverging_iterates_tol", dval);
            problem->ipopt->Options()->SetNumericValue("nlp_lower_bound_inf", -dval);
            problem->ipopt->Options()->SetNumericValue("nlp_upper_bound_inf",  dval);
            SCIPnlpiOracleSetInfinity(problem->oracle, dval);
         }
         else
         {
            SCIP_NLPIDATA* data = SCIPnlpiGetData(nlpi);
            assert(data != NULL);
            data->infinity = dval;
         }
         break;
      }
      
      case SCIP_NLPPAR_ITLIM:
      {
         SCIPerrorMessage("iteration limit parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_TILIM:
      {
         if( dval >= 0 )
         {
            problem->ipopt->Options()->SetNumericValue("max_cpu_time", dval);
         }
         else
         {
            SCIPerrorMessage("Value %g for parameter time limit is negative\n", dval);
            return SCIP_PARAMETERWRONGVAL;
         }
         break;
      }

      case SCIP_NLPPAR_OPTFILE:
      {
         SCIPerrorMessage("option file parameter is of type string.\n");
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
SCIP_DECL_NLPIGETSTRINGPAR( nlpiGetStringParIpopt )
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   switch( type )
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
         SCIPerrorMessage("feasibility tolerance parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_RELOBJTOL:
      {
         SCIPerrorMessage("objective tolerance parameter is of type real.\n");
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
         SCIPerrorMessage("iteration limit parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_TILIM:
      {
         SCIPerrorMessage("time limit parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
      
      case SCIP_NLPPAR_OPTFILE:
      {
         if( !problem->optfile.empty() )
            *sval = problem->optfile.c_str();
         else
            *sval = NULL;
         return SCIP_OKAY;
      }
      
      default:
      {
         SCIPerrorMessage("Parameter %d not known to Ipopt interface.\n", type);
         return SCIP_PARAMETERUNKNOWN;
      }
   }
   
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
SCIP_DECL_NLPISETSTRINGPAR( nlpiSetStringParIpopt )
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   switch( type )
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
         SCIPerrorMessage("feasibility tolerance parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_RELOBJTOL:
      {
         SCIPerrorMessage("objective tolerance parameter is of type real.\n");
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
         SCIPerrorMessage("iteration limit parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_TILIM:
      {
         SCIPerrorMessage("time limit parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }
      
      case SCIP_NLPPAR_OPTFILE:
      {
         if( sval != NULL )
            problem->optfile = sval;
         else
            problem->optfile.clear();
         
         if( problem->ipopt->Initialize(problem->optfile) != Solve_Succeeded )
         {
            SCIPerrorMessage("Error initializing Ipopt using optionfile \"%s\"\n", problem->optfile.c_str());
            return SCIP_ERROR;
         }
         problem->firstrun = TRUE;
         
         return SCIP_OKAY;
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
SCIP_RETCODE SCIPcreateNlpSolverIpopt(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi                /**< pointer to buffer for nlpi address */
   )
{
   SCIP_NLPIDATA* nlpidata;

   assert(blkmem != NULL);
   assert(nlpi   != NULL);
   
   if( BMSallocMemory(&nlpidata) == NULL )
      return SCIP_NOMEMORY;
   nlpidata->blkmem = blkmem;
   nlpidata->infinity = SCIP_DEFAULT_INFINITY;
   
   SCIP_CALL( SCIPnlpiCreate(nlpi,
      NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
      nlpiFreeIpopt, nlpiGetSolverPointerIpopt, nlpiCreateProblemIpopt, nlpiFreeProblemIpopt, nlpiGetProblemPointerIpopt,
      nlpiAddVarsIpopt, nlpiAddConstraintsIpopt, nlpiSetObjectiveIpopt, 
      nlpiChgVarBoundsIpopt, nlpiChgConsBoundsIpopt, nlpiDelVarSetIpopt, nlpiDelConstraintSetIpopt,
      nlpiChgLinearCoefsIpopt, nlpiChgQuadraticCoefsIpopt,
      nlpiChgNonlinCoefIpopt, nlpiSetInitialGuessIpopt,
      nlpiSolveIpopt, nlpiGetSolstatIpopt, nlpiGetSoltermIpopt, nlpiGetSolutionIpopt, nlpiGetStatisticsIpopt,
      nlpiGetWarmstatSizeIpopt, nlpiGetWarmstatMemoIpopt,
      nlpiSetWarmstatMemoIpopt, nlpiGetIntParIpopt,
      nlpiSetIntParIpopt, nlpiGetRealParIpopt, nlpiSetRealParIpopt, nlpiGetStringParIpopt, nlpiSetStringParIpopt,
      nlpidata) );

   return SCIP_OKAY;
}

/** Method to return some info about the nlp */
bool ScipNLP::get_nlp_info(
   Index&             n,                  /**< place to store number of variables */ 
   Index&             m,                  /**< place to store number of constraints */ 
   Index&             nnz_jac_g,          /**< place to store number of nonzeros in jacobian */
   Index&             nnz_h_lag,          /**< place to store number of nonzeros in hessian */
   IndexStyleEnum&    index_style         /**< place to store used index style (0-based or 1-based) */
   )
{
   const int* offset;
   SCIP_RETCODE retcode;

   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);
   
   n = SCIPnlpiOracleGetNVars(nlpiproblem->oracle);
   m = SCIPnlpiOracleGetNConstraints(nlpiproblem->oracle);
   
   retcode = SCIPnlpiOracleGetJacobianSparsity(nlpiproblem->oracle, &offset, NULL);
   if( retcode != SCIP_OKAY )
      return false;
   assert(offset != NULL);
   nnz_jac_g = offset[m];

   retcode = SCIPnlpiOracleGetHessianLagSparsity(nlpiproblem->oracle, &offset, NULL);
   if( retcode != SCIP_OKAY )
      return false;
   assert(offset != NULL);
   nnz_h_lag = offset[n];
   
   index_style = TNLP::C_STYLE;
   
   return true;
}

/** Method to return the bounds for my problem */
bool ScipNLP::get_bounds_info(
   Index              n,                  /**< number of variables */ 
   Number*            x_l,                /**< buffer to store lower bounds on variables */
   Number*            x_u,                /**< buffer to store upper bounds on variables */
   Index              m,                  /**< number of constraints */
   Number*            g_l,                /**< buffer to store lower bounds on constraints */
   Number*            g_u                 /**< buffer to store lower bounds on constraints */
   )
{
   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));
   assert(m == SCIPnlpiOracleGetNConstraints(nlpiproblem->oracle));
   
   assert(SCIPnlpiOracleGetVarLbs(nlpiproblem->oracle) != NULL);
   assert(SCIPnlpiOracleGetVarUbs(nlpiproblem->oracle) != NULL);
   assert(SCIPnlpiOracleGetConstraintLhss(nlpiproblem->oracle) != NULL);
   assert(SCIPnlpiOracleGetConstraintRhss(nlpiproblem->oracle) != NULL);
   
   BMScopyMemoryArray(x_l, SCIPnlpiOracleGetVarLbs(nlpiproblem->oracle), n);
   BMScopyMemoryArray(x_u, SCIPnlpiOracleGetVarUbs(nlpiproblem->oracle), n);
   BMScopyMemoryArray(g_l, SCIPnlpiOracleGetConstraintLhss(nlpiproblem->oracle), m);
   BMScopyMemoryArray(g_u, SCIPnlpiOracleGetConstraintRhss(nlpiproblem->oracle), m);

   return true;
}

/** Method to return the starting point for the algorithm */
bool ScipNLP::get_starting_point(
   Index              n,                  /**< number of variables */ 
   bool               init_x,             /**< whether initial values for primal values are requested */ 
   Number*            x,                  /**< buffer to store initial primal values */
   bool               init_z,             /**< whether initial values for dual values of variable bounds are requested */  
   Number*            z_L,                /**< buffer to store dual values for variable lower bounds */
   Number*            z_U,                /**< buffer to store dual values for variable upper bounds */
   Index              m,                  /**< number of constraints */
   bool               init_lambda,        /**< whether initial values for dual values of constraints are required */
   Number*            lambda              /**< buffer to store dual values of constraints */
   )
{
   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));
   assert(m == SCIPnlpiOracleGetNConstraints(nlpiproblem->oracle));

   if( init_x )
   {
      if( nlpiproblem->initguess )
      {
         BMScopyMemoryArray(x, nlpiproblem->initguess, n);
      }
      else
      {
         SCIP_Real lb, ub;
         SCIPdebugMessage("Ipopt started without intial primal values; make up starting guess by projecting 0 onto variable bounds\n");
         for( int i = 0; i < n; ++i )
         {
            lb = SCIPnlpiOracleGetVarLbs(nlpiproblem->oracle)[i];
            ub = SCIPnlpiOracleGetVarUbs(nlpiproblem->oracle)[i];
            if( lb > 0.0 )
               x[i] = lb;
            else if( ub < 0.0 )
               x[i] = ub;
            else
               x[i] = 0.0;
         }
      }
   }
   if( init_z || init_lambda )
      return false;

   return true;
}

/** Method to return the variables linearity. */
bool ScipNLP::get_variables_linearity(
   Index              n,                  /**< number of variables */ 
   LinearityType*     var_types           /**< buffer to store linearity types of variables */
   )
{
   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);

   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));
   
   for( int i = 0; i < n; ++i )
      var_types[i] = (SCIPnlpiOracleGetVarDegree(nlpiproblem->oracle, i) <= 1 ? LINEAR : NON_LINEAR);
   
   return true;
}

/** Method to return the constraint linearity. */
bool ScipNLP::get_constraints_linearity(
   Index              m,                  /**< number of constraints */
   LinearityType*     const_types         /**< buffer to store linearity types of constraints */
   )
{
   int i;

   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);

   assert(m == SCIPnlpiOracleGetNConstraints(nlpiproblem->oracle));
   
   for( i = 0; i < m; ++i )
      const_types[i] = (SCIPnlpiOracleGetConstraintDegree(nlpiproblem->oracle, i) <= 1 ? LINEAR : NON_LINEAR);
   
   return true;
}

/** Method to return the objective value */
bool ScipNLP::eval_f(
   Index              n,                  /**< number of variables */ 
   const Number*      x,                  /**< point to evaluate */ 
   bool               new_x,              /**< whether some function evaluation method has been called for this point before */
   Number&            obj_value           /**< place to store objective function value */
   )
{
   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));

   return (SCIPnlpiOracleEvalObjectiveValue(nlpiproblem->oracle, x, &obj_value) == SCIP_OKAY ? true : false);
}

/** Method to return the gradient of the objective */
bool ScipNLP::eval_grad_f(
   Index              n,                  /**< number of variables */ 
   const Number*      x,                  /**< point to evaluate */ 
   bool               new_x,              /**< whether some function evaluation method has been called for this point before */
   Number*            grad_f              /**< buffer to store objective gradient */
   )
{
   SCIP_Real dummy;

   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));

   return (SCIPnlpiOracleEvalObjectiveGradient(nlpiproblem->oracle, x, TRUE, &dummy, grad_f) == SCIP_OKAY ? true : false);
}

/** Method to return the constraint residuals */
bool ScipNLP::eval_g(
   Index              n,                  /**< number of variables */ 
   const Number*      x,                  /**< point to evaluate */ 
   bool               new_x,              /**< whether some function evaluation method has been called for this point before */
   Index              m,                  /**< number of constraints */
   Number*            g                   /**< buffer to store constraint function values */
   )
{
   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));

   return (SCIPnlpiOracleEvalConstraintValues(nlpiproblem->oracle, x, g) == SCIP_OKAY ? true : false);
}

/** Method to return:
 *   1) The structure of the jacobian (if "values" is NULL)
 *   2) The values of the jacobian (if "values" is not NULL)
 */
bool ScipNLP::eval_jac_g(
   Index              n,                  /**< number of variables */ 
   const Number*      x,                  /**< point to evaluate */ 
   bool               new_x,              /**< whether some function evaluation method has been called for this point before */
   Index              m,                  /**< number of constraints */
   Index              nele_jac,           /**< number of nonzero entries in jacobian */ 
   Index*             iRow,               /**< buffer to store row indices of nonzero jacobian entries, or NULL if values are requested */
   Index*             jCol,               /**< buffer to store column indices of nonzero jacobian entries, or NULL if values are requested */                  
   Number*            values              /**< buffer to store values of nonzero jacobian entries, or NULL if structure is requested */
   )
{
   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));
   assert(m == SCIPnlpiOracleGetNConstraints(nlpiproblem->oracle));
   
   if( values == NULL )
   { /* Ipopt wants to know sparsity structure */
      const int* jacoffset;
      const int* jaccol;
      int j;
      int i;

      assert(iRow != NULL);
      assert(jCol != NULL);
      
      if( SCIPnlpiOracleGetJacobianSparsity(nlpiproblem->oracle, &jacoffset, &jaccol) != SCIP_OKAY )
         return false;
      
      assert(jacoffset[0] == 0);
      assert(jacoffset[m] == nele_jac);
      j = jacoffset[0];
      for( i = 0; i < m; ++i )
         for( ; j < jacoffset[i+1]; ++j )
            iRow[j] = i;
      
      BMScopyMemoryArray(jCol, jaccol, nele_jac);
   }
   else
   {
      if( SCIPnlpiOracleEvalJacobian(nlpiproblem->oracle, x, TRUE, NULL, values) != SCIP_OKAY )
         return false;
   }

   return true;
}

/** Method to return:
 *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
 *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
 */
bool ScipNLP::eval_h(
   Index              n,                  /**< number of variables */ 
   const Number*      x,                  /**< point to evaluate */ 
   bool               new_x,              /**< whether some function evaluation method has been called for this point before */
   Number             obj_factor,         /**< weight for objective function */ 
   Index              m,                  /**< number of constraints */
   const Number*      lambda,             /**< weights for constraint functions */ 
   bool               new_lambda,         /**< whether the hessian has been evaluated for these values of lambda before */
   Index              nele_hess,          /**< number of nonzero entries in hessian */
   Index*             iRow,               /**< buffer to store row indices of nonzero hessian entries, or NULL if values are requested */
   Index*             jCol,               /**< buffer to store column indices of nonzero hessian entries, or NULL if values are requested */                  
   Number*            values              /**< buffer to store values of nonzero hessian entries, or NULL if structure is requested */
)
{
   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));
   assert(m == SCIPnlpiOracleGetNConstraints(nlpiproblem->oracle));
   
   if( values == NULL )
   { /* Ipopt wants to know sparsity structure */
      const int* heslagoffset;
      const int* heslagcol;
      int j;
      int i;

      assert(iRow != NULL);
      assert(jCol != NULL);
     
      if( SCIPnlpiOracleGetHessianLagSparsity(nlpiproblem->oracle, &heslagoffset, &heslagcol) != SCIP_OKAY )
         return false;
     
      assert(heslagoffset[0] == 0);
      assert(heslagoffset[n] == nele_hess);
      j = heslagoffset[0];
      for( i = 0; i < n; ++i )
         for( ; j < heslagoffset[i+1]; ++j )
            iRow[j] = i;
     
      BMScopyMemoryArray(jCol, heslagcol, nele_hess);
   }
   else
   {
      if( SCIPnlpiOracleEvalHessianLag(nlpiproblem->oracle, x, TRUE, obj_factor, lambda, values) != SCIP_OKAY )
         return false;
   }
   
   return true;
}

/** Method called by the solver at each iteration.
 * 
 * Checks whether Ctrl-C was hit.
 */
bool ScipNLP::intermediate_callback(
   AlgorithmMode      mode,               /**< current mode of algorithm */
   Index              iter,               /**< current iteration number */
   Number             obj_value,          /**< current objective value */
   Number             inf_pr,             /**< current primal infeasibility */
   Number             inf_du,             /**< current dual infeasibility */
   Number             mu,                 /**< current barrier parameter */
   Number             d_norm,             /**< current gradient norm */
   Number             regularization_size,/**< current size of regularization */
   Number             alpha_du,           /**< current dual alpha */
   Number             alpha_pr,           /**< current primal alpha */
   Index              ls_trials,          /**< current number of linesearch trials */
   const IpoptData*   ip_data,            /**< pointer to Ipopt Data */
   IpoptCalculatedQuantities* ip_cq       /**< pointer to current calculated quantities */
)
{
   return (SCIPinterrupted() == FALSE);
}

/** This method is called when the algorithm is complete so the TNLP can store/write the solution.
 */
void ScipNLP::finalize_solution(
   SolverReturn       status,             /**< solve and solution status */ 
   Index              n,                  /**< number of variables */ 
   const Number*      x,                  /**< primal solution values */ 
   const Number*      z_L,                /**< dual values of variable lower bounds */
   const Number*      z_U,                /**< dual values of variable upper bounds */
   Index              m,                  /**< number of constraints */ 
   const Number*      g,                  /**< values of constraints */ 
   const Number*      lambda,             /**< dual values of constraints */ 
   Number             obj_value,          /**< objective function value */ 
   const IpoptData*   data,               /**< pointer to Ipopt Data */ 
   IpoptCalculatedQuantities* cq          /**< pointer to calculated quantities */
)
{
   assert(nlpiproblem         != NULL);
   assert(nlpiproblem->oracle != NULL);
   
   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));
   assert(m == SCIPnlpiOracleGetNConstraints(nlpiproblem->oracle));
   
   bool check_feasibility = false; // whether we should check x for feasibility, if not NULL
   switch( status )
   {
      case SUCCESS:
      case STOP_AT_ACCEPTABLE_POINT:
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_LOCOPT;
         nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
         assert(x != NULL);
         break;
         
      case FEASIBLE_POINT_FOUND:
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
         nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
         assert(x != NULL);
         break;
         
      case MAXITER_EXCEEDED:
         check_feasibility = true;
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
         nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_ITLIM;
         break;
         
      case CPUTIME_EXCEEDED:
         check_feasibility = true;
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
         nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_TILIM;
         break;
         
      case STOP_AT_TINY_STEP:
      case RESTORATION_FAILURE:
      case ERROR_IN_STEP_COMPUTATION:
         check_feasibility = true;
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
         nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_NUMERR;
         break;

      case LOCAL_INFEASIBILITY:
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
         nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
         break;

      case DIVERGING_ITERATES:
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNBOUNDED;
         nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_UOBJLIM;
         break;

      case INVALID_NUMBER_DETECTED:
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
         nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_EVALERR;
         break;

      case USER_REQUESTED_STOP:
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
         nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_TILIM;
         break;

      case TOO_FEW_DEGREES_OF_FREEDOM:
      case INTERNAL_ERROR:
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
         nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
         break;

      default:
         SCIPerrorMessage("Ipopt returned with unknown solution status %d\n", status);
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
         nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
         break;
   }

   if( x != NULL )
   {
      if( nlpiproblem->lastsol == NULL )
      {
         if( BMSduplicateMemoryArray(&nlpiproblem->lastsol, x, n) == NULL )
         {
            nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
            nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_MEMERR;
            return;
         }
      }
      else
      {
         BMScopyMemoryArray(nlpiproblem->lastsol, x, n);
      }
      
      if( check_feasibility && cq != NULL )
      {
         Number constrviol;
         Number constrvioltol;

         constrviol = cq->curr_constraint_violation();

         nlpiproblem->ipopt->Options()->GetNumericValue("constr_viol_tol", constrvioltol, "");
         if( constrviol <= constrvioltol )
         {
            nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
         }
         else
         {
            nlpiproblem->ipopt->Options()->GetNumericValue("acceptable_constr_viol_tol", constrvioltol, "");
            if( constrviol <= constrvioltol )
               nlpiproblem->lastsolstat = SCIP_NLPSOLSTAT_FEASIBLE;
         }
      }
   }
}
