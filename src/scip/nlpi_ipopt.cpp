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

/**@file    nlpi_ipopt.cpp
 * @ingroup NLPIS
 * @brief   Ipopt NLP interface
 * @author  Stefan Vigerske
 * @author  Benjamin Müller
 *
 * @todo warm starts
 * @todo use new_x: Ipopt sets new_x = false if any function has been evaluated for the current x already, while oracle allows new_x to be false only if the current function has been evaluated for the current x before
 * @todo influence output by SCIP verblevel, too, e.g., print strong warnings if SCIP verblevel is full; but currently we have no access to SCIP verblevel
 * @todo if too few degrees of freedom, solve a slack-minimization problem instead?
 *
 * This file can only be compiled if Ipopt is available.
 * Otherwise, to resolve public functions, use nlpi_ipopt_dummy.c.
 * Since the dummy code is C instead of C++, it has been moved into a separate file.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/nlpi_ipopt.h"

#include "scip/nlpioracle.h"
#include "scip/exprinterpret.h"
#include "scip/interrupt.h"
#include "scip/scip_nlpi.h"
#include "scip/scip_nlp.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_message.h"
#include "scip/scip_general.h"
#include "scip/scip_numerics.h"
#include "scip/pub_misc.h"

#include <new>      /* for std::bad_alloc */
#include <sstream>

/* fallback to non-thread-safe version if C++ is too old to have mutex */
#if __cplusplus < 201103L && defined(SCIP_THREADSAFE)
#undef SCIP_THREADSAFE
#endif

#ifdef SCIP_THREADSAFE
#include <mutex>
static std::mutex solve_mutex;  /*lint !e1756*/
#endif

/* turn off some lint warnings for file */
/*lint --e{1540,750,3701}*/

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wshadow"
#endif
#define IPOPT_DEPRECATED  // to avoid warnings about using functions that became deprecated in Ipopt 3.14
#include "IpoptConfig.h"
#include "IpIpoptApplication.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpSolveStatistics.hpp"
#include "IpJournalist.hpp"
#include "IpIpoptData.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpLapack.hpp"
#ifdef __GNUC__
#pragma GCC diagnostic warning "-Wshadow"
#endif

#if (IPOPT_VERSION_MAJOR < 3 || (IPOPT_VERSION_MAJOR == 3 && IPOPT_VERSION_MINOR < 12))
#error "The Ipopt interface requires at least 3.12."
#endif

using namespace Ipopt;

#define NLPI_NAME          "ipopt"           /**< short concise name of solver */
#define NLPI_DESC          "Ipopt interface" /**< description of solver */
#define NLPI_PRIORITY      1000              /**< priority */

#ifdef SCIP_DEBUG
#define DEFAULT_PRINTLEVEL J_ITERSUMMARY     /**< default print level of Ipopt */
#else
#define DEFAULT_PRINTLEVEL J_ERROR           /**< default print level of Ipopt */
#endif
#define DEFAULT_MAXITER    3000              /**< default iteration limit for Ipopt */

#define MAXPERTURB         0.01              /**< maximal perturbation of bounds in starting point heuristic */
#define FEASTOLFACTOR      0.05              /**< factor for user-given feasibility tolerance to get feasibility tolerance that is actually passed to Ipopt */

#define DEFAULT_RANDSEED   71                /**< initial random seed */


/* Convergence check (see ScipNLP::intermediate_callback)
 *
 * If the fastfail option is enabled, then we stop Ipopt if the reduction in
 * primal infeasibility is not sufficient for a consecutive number of iterations.
 * With the parameters as given below, we require Ipopt to
 * - not increase the primal infeasibility after 5 iterations
 * - reduce the primal infeasibility by at least 50% within 10 iterations
 * - reduce the primal infeasibility by at least 90% within 30 iterations
 * The targets are updated once they are reached and the limit on allowed iterations to reach the new target is reset.
 *
 * In certain situations, it is allowed to exceed an iteration limit:
 * - If we are in the first 10 (convcheck_startiter) iterations.
 * - If we are within 10 (convcheck_startiter) iterations after the restoration phase ended.
 *   The reason for this is that during feasibility restoration phase Ipopt aims completely on
 *   reducing constraint violation, completely forgetting the objective function.
 *   When returning from feasibility restoration and considering the original objective again,
 *   it is unlikely that Ipopt will continue to decrease primal infeasibility, since it may now target on
 *   more on optimality again. Thus, we do not check convergence for a number of iterations.
 * - If the target on dual infeasibility reduction has been achieved, we are below twice the iteration limit, and
 *   we are not in restoration mode.
 *   The reason for this is that if Ipopt makes good progress towards optimality,
 *   we want to allow some more iterations where primal infeasibility is not reduced.
 *   However, in restoration mode, dual infeasibility does not correspond to the original problem and
 *   the complete aim is to restore primal infeasibility.
 */
static const int convcheck_nchecks                         = 3;                 /**< number of convergence checks */
static const int convcheck_startiter                       = 10;                /**< iteration where to start convergence checking */
static const int convcheck_maxiter[convcheck_nchecks]      = { 5,   15,  30 };  /**< maximal number of iterations to achieve each convergence check */
static const SCIP_Real convcheck_minred[convcheck_nchecks] = { 1.0, 0.5, 0.1 }; /**< minimal required infeasibility reduction in each convergence check */

class ScipNLP;

struct SCIP_NlpiData
{
public:
   std::string                 defoptions;   /**< modified default options for Ipopt */

   /** constructor */
   explicit SCIP_NlpiData() { }
};

struct SCIP_NlpiProblem
{
public:
   SCIP_NLPIORACLE*            oracle;       /**< Oracle-helper to store and evaluate NLP */

   SmartPtr<IpoptApplication>  ipopt;        /**< Ipopt application */
   SmartPtr<ScipNLP>           nlp;          /**< NLP in Ipopt form */
   std::string                 optfile;      /**< name of options file */
   bool                        storeintermediate;/**< whether to store intermediate solutions */
   bool                        fastfail;     /**< whether to stop Ipopt if convergence seems slow */

   SCIP_Bool                   firstrun;     /**< whether the next NLP solve will be the first one (with the current problem structure) */
   SCIP_Real*                  initguess;    /**< initial values for primal variables, or NULL if not known */

   SCIP_NLPSOLSTAT             lastsolstat;  /**< solution status from last run */
   SCIP_NLPTERMSTAT            lasttermstat; /**< termination status from last run */
   SCIP_Real*                  lastsolprimals; /**< primal solution values from last run, if available */
   SCIP_Real*                  lastsoldualcons; /**< dual solution values of constraints from last run, if available */
   SCIP_Real*                  lastsoldualvarlb; /**< dual solution values of variable lower bounds from last run, if available */
   SCIP_Real*                  lastsoldualvarub; /**< dual solution values of variable upper bounds from last run, if available */
   SCIP_Real                   lastsolinfeas;/**< infeasibility (constraint violation) of solution stored in lastsolprimals */
   int                         lastniter;    /**< number of iterations in last run */
   SCIP_Real                   lasttime;     /**< time spend in last run */

   /** constructor */
   SCIP_NlpiProblem()
      : oracle(NULL),
        storeintermediate(false), fastfail(false),
        firstrun(TRUE), initguess(NULL),
        lastsolstat(SCIP_NLPSOLSTAT_UNKNOWN), lasttermstat(SCIP_NLPTERMSTAT_OTHER),
        lastsolprimals(NULL), lastsoldualcons(NULL), lastsoldualvarlb(NULL), lastsoldualvarub(NULL),
        lastsolinfeas(0.0), lastniter(-1), lasttime(-1.0)
   { }
};

/** TNLP implementation for SCIPs NLP */
class ScipNLP : public TNLP
{
private:
   SCIP_NLPIPROBLEM*     nlpiproblem;        /**< NLPI problem data */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP*                 scip;               /**< SCIP data structure */

   SCIP_Real             conv_prtarget[convcheck_nchecks]; /**< target primal infeasibility for each convergence check */
   SCIP_Real             conv_dutarget[convcheck_nchecks]; /**< target dual infeasibility for each convergence check */
   int                   conv_iterlim[convcheck_nchecks];  /**< iteration number where target primal infeasibility should to be achieved */
   int                   conv_lastrestoiter;               /**< last iteration number in restoration mode, or -1 if none */

public:
   bool                  approxhessian;      /**< do we tell Ipopt to approximate the hessian? (may also be false if user set to approx. hessian via option file) */

   // cppcheck-suppress uninitMemberVar
   /** constructor */
   ScipNLP(
      SCIP_NLPIPROBLEM*  nlpiproblem_ = NULL,/**< NLPI problem data */
      SCIP*              scip_ = NULL        /**< SCIP data structure */
      )
      : nlpiproblem(nlpiproblem_), randnumgen(NULL), scip(scip_), conv_lastrestoiter(-1), approxhessian(false)
   {
      assert(scip != NULL);
      SCIP_CALL_ABORT_QUIET( SCIPcreateRandom(scip, &randnumgen, DEFAULT_RANDSEED, TRUE) );
   }

   /** destructor */
   ~ScipNLP()
   { /*lint --e{1540}*/
      assert(randnumgen != NULL);
      SCIPfreeRandom(scip, &randnumgen);
   }

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

   /** Method to return the number of nonlinear variables. */
   Index get_number_of_nonlinear_variables();

   /** Method to return the indices of the nonlinear variables */
   bool get_list_of_nonlinear_variables(
      Index              num_nonlin_vars,    /**< number of nonlinear variables */
      Index*             pos_nonlin_vars     /**< array to fill with indices of nonlinear variables */
      );

   /** Method to return metadata about variables and constraints */
   bool get_var_con_metadata(
      Index              n,                  /**< number of variables */
      StringMetaDataMapType& var_string_md,  /**< variable meta data of string type */
      IntegerMetaDataMapType& var_integer_md,/**< variable meta data of integer type */
      NumericMetaDataMapType& var_numeric_md,/**< variable meta data of numeric type */
      Index              m,                  /**< number of constraints */
      StringMetaDataMapType& con_string_md,  /**< constraint meta data of string type */
      IntegerMetaDataMapType& con_integer_md,/**< constraint meta data of integer type */
      NumericMetaDataMapType& con_numeric_md /**< constraint meta data of numeric type */
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

/** A particular Ipopt::Journal implementation that uses the SCIP message routines for output. */
class ScipJournal : public Ipopt::Journal {
private:
   SCIP*                 scip;               /**< SCIP data structure */

public:
   ScipJournal(
      const char*          name,             /**< name of journal */
      Ipopt::EJournalLevel default_level,    /**< default verbosity level */
      SCIP*                scip_             /**< SCIP data structure */
      )
      : Ipopt::Journal(name, default_level),
        scip(scip_)
   { }

   ~ScipJournal() { }

protected:
   /*lint -e{715}*/
   void PrintImpl(
      Ipopt::EJournalCategory category,      /**< category of message */
      Ipopt::EJournalLevel    level,         /**< verbosity level of message */
      const char*             str            /**< message to print */
      )
   {  /*lint --e{715} */
      if( level == J_ERROR )
      {
         SCIPmessagePrintError("%s", str);
      }
      else
      {
         SCIPinfoMessage(scip, NULL, "%s", str);
      }
   }

   /*lint -e{715}*/
   void PrintfImpl(
      Ipopt::EJournalCategory category,      /**< category of message */
      Ipopt::EJournalLevel    level,         /**< verbosity level of message */
      const char*             pformat,       /**< message printing format */
      va_list                 ap             /**< arguments of message */
      )
   {  /*lint --e{715} */
      if( level == J_ERROR )
      {
         SCIPmessageVPrintError(pformat, ap);
      }
      else
      {
         SCIPmessageVPrintInfo(SCIPgetMessagehdlr(scip), pformat, ap);
      }
   }

   void FlushBufferImpl() { }
};

/** clears the last solution arrays and sets the solstat and termstat to unknown and other, resp. */
static
void invalidateSolution(
   SCIP_NLPIPROBLEM*     problem             /**< data structure of problem */
   )
{
   assert(problem != NULL);

   BMSfreeMemoryArrayNull(&problem->lastsolprimals);
   BMSfreeMemoryArrayNull(&problem->lastsoldualcons);
   BMSfreeMemoryArrayNull(&problem->lastsoldualvarlb);
   BMSfreeMemoryArrayNull(&problem->lastsoldualvarub);
   problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   problem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
   problem->lastsolinfeas = SCIP_INVALID;
}

/** sets feasibility tolerance parameters in Ipopt
 *
 * Sets tol and constr_viol_tol to FEASTOLFACTOR*feastol and acceptable_tol and acceptable_viol_tol to feastol.
 * Since the users and Ipopts conception of feasibility may differ, we let Ipopt try to compute solutions
 * that are more accurate (w.r.t. constraint violation) than requested by the user.
 * Only if Ipopt has problems to achieve this accuracy, we also accept solutions that are accurate w.r.t. feastol only.
 * The additional effort for computing a more accurate solution should be small if one can assume fast convergence when close to a local minimizer.
 */
static
void setFeastol(
   SCIP_NLPIPROBLEM* nlpiproblem,
   SCIP_Real         feastol
   )
{
   assert(nlpiproblem != NULL);

   (void) nlpiproblem->ipopt->Options()->SetNumericValue("tol", FEASTOLFACTOR * feastol);
   (void) nlpiproblem->ipopt->Options()->SetNumericValue("constr_viol_tol", FEASTOLFACTOR * feastol);

   (void) nlpiproblem->ipopt->Options()->SetNumericValue("acceptable_tol", feastol);
   (void) nlpiproblem->ipopt->Options()->SetNumericValue("acceptable_constr_viol_tol", feastol);

   /* It seem to be better to let Ipopt relax bounds a bit to ensure that a relative interior exists.
    * However, if we relax the bounds too much, then the solutions tend to be slightly infeasible.
    * If the user wants to set a tight feasibility tolerance, then (s)he has probably difficulties to compute accurate enough solutions.
    * Thus, we turn off the bound_relax_factor completely if it would be below its default value of 1e-8.
    */
   (void) nlpiproblem->ipopt->Options()->SetNumericValue("bound_relax_factor", feastol < 1e-8/FEASTOLFACTOR ? 0.0 : FEASTOLFACTOR * feastol);
}

/** sets optimality tolerance parameters in Ipopt
 *
 * Sets dual_inf_tol and compl_inf_tol to opttol.
 * We leave acceptable_dual_inf_tol and acceptable_compl_inf_tol untouched, which means that if Ipopt has convergence problems, then
 * it can stop with a solution that is still feasible (see setFeastol), but essentially without a proof of local optimality.
 * Note, that we report the solution as feasible only if Ipopt stopped on an "acceptable point" (see ScipNLP::finalize_solution).
 *
 * Note, that parameters tol and acceptable_tol (set in setFeastol) also apply.
 */
static
void setOpttol(
   SCIP_NLPIPROBLEM* nlpiproblem,
   SCIP_Real         opttol
   )
{
   assert(nlpiproblem != NULL);

   (void) nlpiproblem->ipopt->Options()->SetNumericValue("dual_inf_tol", opttol);
   (void) nlpiproblem->ipopt->Options()->SetNumericValue("compl_inf_tol", opttol);
}

/** copy method of NLP interface (called when SCIP copies plugins)
 *
 *  input:
 *  - blkmem block memory of target SCIP
 *  - sourcenlpi the NLP interface to copy
 *  - targetnlpi buffer to store pointer to copy of NLP interface
 */
static
SCIP_DECL_NLPICOPY(nlpiCopyIpopt)
{
   SCIP_NLPI* targetnlpi;
   SCIP_NLPIDATA* sourcedata;
   SCIP_NLPIDATA* targetdata;

   assert(sourcenlpi != NULL);

   SCIP_CALL( SCIPincludeNlpSolverIpopt(scip) );

   targetnlpi = SCIPfindNlpi(scip, NLPI_NAME);
   assert(targetnlpi != NULL);

   sourcedata = SCIPnlpiGetData(sourcenlpi);
   assert(sourcedata != NULL);

   targetdata = SCIPnlpiGetData(targetnlpi);
   assert(targetdata != NULL);

   targetdata->defoptions = sourcedata->defoptions;

   return SCIP_OKAY;
}

/** destructor of NLP interface to free nlpi data
 * 
 * input:
 *  - nlpi datastructure for solver interface
 */
static
SCIP_DECL_NLPIFREE(nlpiFreeIpopt)
{
   assert(nlpidata != NULL);

   delete *nlpidata;
   *nlpidata = NULL;

   return SCIP_OKAY;
}

/** gets pointer for NLP solver to do dirty stuff
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *
 *  return: void pointer to solver
 */
static
SCIP_DECL_NLPIGETSOLVERPOINTER(nlpiGetSolverPointerIpopt)
{
   assert(nlpi != NULL);

   return NULL;
}

/** creates a problem instance
 *
 *  input:
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

   SCIP_ALLOC( *problem = new SCIP_NLPIPROBLEM ); /*lint !e774*/

   SCIP_CALL( SCIPnlpiOracleCreate(scip, &(*problem)->oracle) );
   SCIP_CALL( SCIPnlpiOracleSetProblemName(scip, (*problem)->oracle, name) );

   try
   {
      /* initialize IPOPT without default journal */
      (*problem)->ipopt = new IpoptApplication(false);
      if( IsNull((*problem)->ipopt) )
         throw std::bad_alloc();

      /* plugin our journal to get output through SCIP message handler */
      SmartPtr<Journal> jrnl = new ScipJournal("console", J_ITERSUMMARY, scip);
      if( IsNull(jrnl) )
         throw std::bad_alloc();
      jrnl->SetPrintLevel(J_DBG, J_NONE);
      if( !(*problem)->ipopt->Jnlst()->AddJournal(jrnl) )
      {
         SCIPerrorMessage("Failed to register ScipJournal for IPOPT output.");
      }

      /* initialize Ipopt/SCIP NLP interface */
      (*problem)->nlp = new ScipNLP(*problem, scip);
      if( IsNull((*problem)->nlp) )
         throw std::bad_alloc();
   }
   catch( const std::bad_alloc& )
   {
      SCIPerrorMessage("Not enough memory to initialize Ipopt.\n");
      return SCIP_NOMEMORY;
   }

   /* modify Ipopt's default settings to what we believe is appropriate */
   (*problem)->ipopt->RegOptions()->AddStringOption2("store_intermediate", "whether to store the most feasible intermediate solutions", "no", "yes", "", "no", "", "useful when Ipopt looses a once found feasible solution and then terminates with an infeasible point");
   (void) (*problem)->ipopt->Options()->SetIntegerValue("print_level", DEFAULT_PRINTLEVEL);
   /* (*problem)->ipopt->Options()->SetStringValue("print_timing_statistics", "yes"); */
#ifdef SCIP_DEBUG
   (void) (*problem)->ipopt->Options()->SetStringValue("print_user_options", "yes");
#endif
   (void) (*problem)->ipopt->Options()->SetStringValue("mu_strategy", "adaptive");
   (void) (*problem)->ipopt->Options()->SetIntegerValue("max_iter", DEFAULT_MAXITER);
   (void) (*problem)->ipopt->Options()->SetNumericValue("nlp_lower_bound_inf", -SCIPinfinity(scip), false);
   (void) (*problem)->ipopt->Options()->SetNumericValue("nlp_upper_bound_inf",  SCIPinfinity(scip), false);
   (void) (*problem)->ipopt->Options()->SetNumericValue("diverging_iterates_tol", SCIPinfinity(scip), false);
   /* (void) (*problem)->ipopt->Options()->SetStringValue("dependency_detector", "ma28"); */
   /* if the expression interpreter does not give hessians, tell Ipopt to approximate hessian */
   setFeastol(*problem, SCIP_DEFAULT_FEASTOL);
   setOpttol(*problem, SCIP_DEFAULT_DUALFEASTOL);

   /* apply user's given modifications to Ipopt's default settings */
   if( data->defoptions.length() > 0 )
   {
      std::istringstream is(data->defoptions);

#if (IPOPT_VERSION_MAJOR > 3) || (IPOPT_VERSION_MAJOR == 3 && IPOPT_VERSION_MINOR > 12) || (IPOPT_VERSION_MAJOR == 3 && IPOPT_VERSION_MINOR == 12 && IPOPT_VERSION_RELEASE >= 5)
      if( !(*problem)->ipopt->Options()->ReadFromStream(*(*problem)->ipopt->Jnlst(), is, true) )
#else
      if( !(*problem)->ipopt->Options()->ReadFromStream(*(*problem)->ipopt->Jnlst(), is) )
#endif
      {
         SCIPerrorMessage("Error when modifying Ipopt options using options string\n%s\n", data->defoptions.c_str());
         return SCIP_ERROR;
      }
   }

   /* apply user's given options file (this one is NLPI problem specific) */
   if( (*problem)->ipopt->Initialize((*problem)->optfile) != Solve_Succeeded )
   {
      SCIPerrorMessage("Error during initialization of Ipopt using optionfile \"%s\"\n", (*problem)->optfile.c_str());
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** free a problem instance
 *
 *  input:
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
      SCIP_CALL( SCIPnlpiOracleFree(scip, &(*problem)->oracle) );
   }

   BMSfreeMemoryArrayNull(&(*problem)->initguess);
   BMSfreeMemoryArrayNull(&(*problem)->lastsolprimals);
   BMSfreeMemoryArrayNull(&(*problem)->lastsoldualcons);
   BMSfreeMemoryArrayNull(&(*problem)->lastsoldualvarlb);
   BMSfreeMemoryArrayNull(&(*problem)->lastsoldualvarub);

   delete *problem;
   *problem = NULL;

   return SCIP_OKAY;
}

/** gets pointer to solver-internal problem instance to do dirty stuff
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 *  return: void pointer to problem instance
 */
static
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

   SCIP_CALL( SCIPnlpiOracleAddVars(scip, problem->oracle, nvars, lbs, ubs, varnames) );

   problem->firstrun = TRUE;
   BMSfreeMemoryArrayNull(&problem->initguess);
   invalidateSolution(problem);

   return SCIP_OKAY;
}

/** add constraints
 *
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
 *  - nquadelems number of quadratic elements for each constraint
 *    may be NULL in case of no quadratic part
 *  - quadelems quadratic elements for each constraint
 *    may be NULL in case of no quadratic part
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

   SCIP_CALL( SCIPnlpiOracleAddConstraints(scip, problem->oracle, nconss, lhss, rhss, nlininds, lininds, linvals, exprs, names) );

   problem->firstrun = TRUE;
   invalidateSolution(problem);

   return SCIP_OKAY;
}

/** sets or overwrites objective, a minimization problem is expected
 *
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
 *  - nquadelems number of elements in matrix of quadratic part
 *  - quadelems elements of quadratic part
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

   /* We pass the objective gradient in dense form to Ipopt, so if the sparsity of that gradient changes, we do not need to reset Ipopt (firstrun=TRUE).
    * However, if the sparsity of the Hessian matrix of the objective changes, then the sparsity pattern of the Hessian of the Lagrangian may change.
    * Thus, reset Ipopt if the objective was and/or becomes nonlinear, but leave firstrun untouched if it was and stays linear.
    */
   if( expr != NULL || SCIPnlpiOracleGetConstraintDegree(problem->oracle, -1) > 1 )
      problem->firstrun = TRUE;

   SCIP_CALL( SCIPnlpiOracleSetObjective(scip, problem->oracle, constant, nlins, lininds, linvals, expr) );

   invalidateSolution(problem);

   return SCIP_OKAY;
}

/** change variable bounds
 *
 *  input:
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
   int i;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   /* If some variable is fixed or unfixed, then better don't reoptimize the NLP in the next solve.
    * Calling Optimize instead of ReOptimize should remove fixed variables from the problem that is solved by Ipopt.
    * This way, the variable fixing is satisfied exactly in a solution, see also #1254.
    */
   for( i = 0; i < nvars && !problem->firstrun; ++i )
      if( (SCIPnlpiOracleGetVarLbs(problem->oracle)[indices[i]] == SCIPnlpiOracleGetVarUbs(problem->oracle)[indices[i]]) != (lbs[i] == ubs[i]) )  /*lint !e777*/
         problem->firstrun = TRUE;

   SCIP_CALL( SCIPnlpiOracleChgVarBounds(scip, problem->oracle, nvars, indices, lbs, ubs) );

   invalidateSolution(problem);

   return SCIP_OKAY;
}

/** change constraint bounds
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nconss number of constraints to change sides
 *  - indices indices of constraints to change sides
 *  - lhss new left hand sides
 *  - rhss new right hand sides
 */
static
SCIP_DECL_NLPICHGCONSSIDES(nlpiChgConsSidesIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgConsSides(scip, problem->oracle, nconss, indices, lhss, rhss) );

   invalidateSolution(problem);

   return SCIP_OKAY;
}

/** delete a set of variables
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nlpi datastructure for solver interface
 *  - dstats deletion status of vars; 1 if var should be deleted, 0 if not
 *
 *  output:
 *  - dstats new position of var, -1 if var was deleted
 */
static
SCIP_DECL_NLPIDELVARSET(nlpiDelVarSetIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelVarSet(scip, problem->oracle, dstats) );

   problem->firstrun = TRUE;
   BMSfreeMemoryArrayNull(&problem->initguess); // @TODO keep initguess for remaining variables 

   invalidateSolution(problem);

   return SCIP_OKAY;
}

/** delete a set of constraints
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - dstats deletion status of rows; 1 if row should be deleted, 0 if not
 *
 *  output:
 *  - dstats new position of row, -1 if row was deleted
 */
static
SCIP_DECL_NLPIDELCONSSET(nlpiDelConstraintSetIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelConsSet(scip, problem->oracle, dstats) );

   problem->firstrun = TRUE;

   invalidateSolution(problem);

   return SCIP_OKAY;
}

/** change one linear coefficient in a constraint or objective
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idx index of constraint or -1 for objective
 *  - nvals number of values in linear constraint
 *  - varidxs indices of variable
 *  - vals new values for coefficient
 *
 *  return: Error if coefficient did not exist before
 */
static
SCIP_DECL_NLPICHGLINEARCOEFS(nlpiChgLinearCoefsIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgLinearCoefs(scip, problem->oracle, idx, nvals, varidxs, vals) );
   invalidateSolution(problem);

   return SCIP_OKAY;
}

/** replaces the expression tree of a constraint or objective
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idxcons index of constraint or -1 for objective
 *  - exprtree new expression tree for constraint or objective, or NULL to only remove previous tree
 */
static
SCIP_DECL_NLPICHGEXPR(nlpiChgExprIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgExpr(scip, problem->oracle, idxcons, expr) );
   invalidateSolution(problem);

   return SCIP_OKAY;
}

/** change the constant offset in the objective
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - objconstant new value for objective constant
 */
static
SCIP_DECL_NLPICHGOBJCONSTANT(nlpiChgObjConstantIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgObjConstant(scip, problem->oracle, objconstant) );

   return SCIP_OKAY;
}

/** sets initial guess for primal variables
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues initial primal values for variables, or NULL to clear previous values
 *  - consdualvalues initial dual values for constraints, or NULL to clear previous values
 *  - varlbdualvalues  initial dual values for variable lower bounds, or NULL to clear previous values
 *  - varubdualvalues  initial dual values for variable upper bounds, or NULL to clear previous values
 */
static
SCIP_DECL_NLPISETINITIALGUESS(nlpiSetInitialGuessIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   if( primalvalues != NULL )
   {
      if( !problem->initguess )
      {
         if( BMSduplicateMemoryArray(&problem->initguess, primalvalues, SCIPnlpiOracleGetNVars(problem->oracle)) == NULL )
            return SCIP_NOMEMORY;
      }
      else
      {
         BMScopyMemoryArray(problem->initguess, primalvalues, SCIPnlpiOracleGetNVars(problem->oracle));
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
 *  input:
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

   assert(IsValid(problem->ipopt));
   assert(IsValid(problem->nlp));

   problem->nlp->setNLPIPROBLEM(problem);

   problem->lastniter = -1;
   problem->lasttime  = -1.0;
   problem->lastsolinfeas = SCIP_INVALID;

   try
   {
      SmartPtr<SolveStatistics> stats;

#ifdef SCIP_THREADSAFE
      /* lock solve_mutex if Ipopt is going to use Mumps as linear solver
       * unlocking will happen in the destructor of guard, which is called when this block is left
       */
      std::unique_lock<std::mutex> guard(solve_mutex, std::defer_lock);  /*lint !e{728}*/
      std::string linsolver;
      (void) problem->ipopt->Options()->GetStringValue("linear_solver", linsolver, "");
      if( linsolver == "mumps" )
         guard.lock();
#endif

      if( problem->firstrun )
      {
         SCIP_EXPRINTCAPABILITY cap;

         cap = SCIPexprintGetCapability() & SCIPnlpiOracleGetEvalCapability(scip, problem->oracle);

         /* if the expression interpreter or some user expression do not support function values and gradients and Hessians,
          * change NLP parameters or give an error
          */
         if( (cap & (SCIP_EXPRINTCAPABILITY_FUNCVALUE | SCIP_EXPRINTCAPABILITY_GRADIENT | SCIP_EXPRINTCAPABILITY_HESSIAN)) != (SCIP_EXPRINTCAPABILITY_FUNCVALUE | SCIP_EXPRINTCAPABILITY_GRADIENT | SCIP_EXPRINTCAPABILITY_HESSIAN) )
         {
            /* @todo could enable Jacobian approximation in Ipopt */
            if( !(SCIPexprintGetCapability() & SCIP_EXPRINTCAPABILITY_FUNCVALUE) ||
                !(SCIPexprintGetCapability() & SCIP_EXPRINTCAPABILITY_GRADIENT) )
            {
               SCIPerrorMessage("Do not have expression interpreter that can compute function values and gradients. Cannot solve NLP with Ipopt.\n");
               problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
               problem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
               return SCIP_OKAY;
            }

            /* enable Hessian approximation if we are nonquadratic and the expression interpreter or user expression do not support Hessians */
            if( !(cap & SCIP_EXPRINTCAPABILITY_HESSIAN) )
            {
               (void) problem->ipopt->Options()->SetStringValue("hessian_approximation", "limited-memory");
               problem->nlp->approxhessian = true;
            }
            else
               problem->nlp->approxhessian = false;
         }

#ifdef SCIP_DEBUG
         // enabling the derivative tester leads to calling TNLP::get_starting_point() twice, which has undesired consequences if the starting point is generated randomly
         // so we don't enable derivative tester if SCIP_DEBUG is defined for now
         // problem->ipopt->Options()->SetStringValue("derivative_test", problem->nlp->approxhessian ? "first-order" : "second-order");
#endif

         status = problem->ipopt->OptimizeTNLP(GetRawPtr(problem->nlp));
      }
      else
      {
         status = problem->ipopt->ReOptimizeTNLP(GetRawPtr(problem->nlp));
      }

      // catch the very bad status codes
      switch( status ) {
         case Invalid_Problem_Definition:
         case Invalid_Option:
         case Unrecoverable_Exception:
         case NonIpopt_Exception_Thrown:
            SCIPerrorMessage("Ipopt returned with application return status %d\n", status);
            return SCIP_ERROR;
         case Internal_Error:
            // could be a fail in the linear solver
            SCIPerrorMessage("Ipopt returned with status \"Internal Error\"\n");
            invalidateSolution(problem);
            problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
            problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
            break;
         case Insufficient_Memory:
            SCIPerrorMessage("Ipopt returned with status \"Insufficient Memory\"\n");
            return SCIP_NOMEMORY;
         case Invalid_Number_Detected:
            SCIPdebugMsg(scip, "Ipopt failed because of an invalid number in function or derivative value\n");
            invalidateSolution(problem);
            problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
            problem->lasttermstat = SCIP_NLPTERMSTAT_EVALERR;
            break;
         default: ;
      }

      stats = problem->ipopt->Statistics();
      if( IsValid(stats) )
      {
         problem->lastniter = stats->IterationCount();
         problem->lasttime  = stats->TotalCPUTime();
      }
      else
      {
         /* Ipopt does not provide access to the statistics when all variables have been fixed */
         problem->lastniter = 0;
         problem->lasttime  = 0.0;
      }
   }
   catch( IpoptException& except )
   {
      SCIPerrorMessage("Ipopt returned with exception: %s\n", except.Message().c_str());
      return SCIP_ERROR;
   }

   problem->firstrun = FALSE;

   return SCIP_OKAY;
}

/** gives solution status
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 *  return: Solution Status
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
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 *  return: Termination Status
 */
static
SCIP_DECL_NLPIGETTERMSTAT(nlpiGetTermstatIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   return problem->lasttermstat;
}

/** gives primal and dual solution values
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues buffer to store pointer to array to primal values, or NULL if not needed
 *  - consdualvalues buffer to store pointer to array to dual values of constraints, or NULL if not needed
 *  - varlbdualvalues buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed
 *  - varubdualvalues buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed
 *  - objval buffer store the objective value, or NULL if not needed
 */
static
SCIP_DECL_NLPIGETSOLUTION(nlpiGetSolutionIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   if( primalvalues != NULL )
      *primalvalues = problem->lastsolprimals;

   if( consdualvalues != NULL )
      *consdualvalues = problem->lastsoldualcons;

   if( varlbdualvalues != NULL )
      *varlbdualvalues = problem->lastsoldualvarlb;

   if( varubdualvalues != NULL )
      *varubdualvalues = problem->lastsoldualvarub;

   if( objval != NULL )
   {
      if( problem->lastsolprimals != NULL )
      {
         /* TODO store last solution value instead of reevaluating the objective function */
         SCIP_CALL( SCIPnlpiOracleEvalObjectiveValue(scip, problem->oracle, problem->lastsolprimals, objval) );
      }
      else
         *objval = SCIP_INVALID;
   }

   return SCIP_OKAY;
}

/** gives solve statistics
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - statistics pointer to store statistics
 *
 *  output:
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
 *  output:
 *  - size required size for warmstart buffer
 */
static
SCIP_DECL_NLPIGETWARMSTARTSIZE(nlpiGetWarmstartSizeIpopt)
{
   SCIPerrorMessage("method of Ipopt nonlinear solver is not implemented\n");
   return SCIP_ERROR;
}

/** stores warmstart information in buffer
 *
 *  Required size of buffer should have been obtained by SCIPnlpiGetWarmstartSize before.
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - buffer memory to store warmstart information
 *
 *  output:
 *  - buffer warmstart information in solver specific data structure
 */
static
SCIP_DECL_NLPIGETWARMSTARTMEMO(nlpiGetWarmstartMemoIpopt)
{
   SCIPerrorMessage("method of Ipopt nonlinear solver is not implemented\n");
   return SCIP_ERROR;
}

/** sets warmstart information in solver
 *
 *  Write warmstart to buffer.
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - buffer warmstart information
 */
static
SCIP_DECL_NLPISETWARMSTARTMEMO(nlpiSetWarmstartMemoIpopt)
{
   SCIPerrorMessage("method of Ipopt nonlinear solver is not implemented\n");
   SCIPABORT();
   return SCIP_OKAY;
}

/** gets integer parameter of NLP
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival pointer to store the parameter value
 *
 *  output:
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
      (void) problem->ipopt->Options()->GetIntegerValue("print_level", printlevel, "");
      if( printlevel <= J_ERROR )
         *ival = 0;
      else if( printlevel >= J_DETAILED )
         *ival = printlevel - J_ITERSUMMARY + 1;
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

   case SCIP_NLPPAR_ITLIM:
   {
      (void) problem->ipopt->Options()->GetIntegerValue("max_iter", *ival, "");
      break;
   }

   case SCIP_NLPPAR_TILIM:
   {
      SCIPerrorMessage("time limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_OPTFILE:
   {
      SCIPerrorMessage("optfile parameter is of type string.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FASTFAIL:
   {
      *ival = problem->fastfail ? 1 : 0;
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

/** sets integer parameter of NLP
 *
 *  input:
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
         SCIPwarningMessage(scip, "from scratch parameter not supported by Ipopt interface yet. Ignored.\n");
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
         (void) problem->ipopt->Options()->SetIntegerValue("print_level", J_ERROR);
         break;
      case 1:
         (void) problem->ipopt->Options()->SetIntegerValue("print_level", J_ITERSUMMARY);
         break;
      case 2:
         (void) problem->ipopt->Options()->SetIntegerValue("print_level", J_DETAILED);
         break;
      default:
         if( ival > 2 )
         {
            (void) problem->ipopt->Options()->SetIntegerValue("print_level", MIN(J_ITERSUMMARY + (ival-1), J_ALL));
            break;
         }
         else
         {
            SCIPerrorMessage("Value %d for parameter from verbosity level out of range {0, 1, 2}\n", ival);
            return SCIP_PARAMETERWRONGVAL;
         }
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

   case SCIP_NLPPAR_ITLIM:
   {
      if( ival >= 0 )
      {
         (void) problem->ipopt->Options()->SetIntegerValue("max_iter", ival);
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

   case SCIP_NLPPAR_OPTFILE:
   {
      SCIPerrorMessage("optfile parameter is of type string.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FASTFAIL:
   {
      if( ival == 0 || ival == 1 )
      {
         problem->fastfail = (bool)ival;
         problem->storeintermediate = (bool)ival;
      }
      else
      {
         SCIPerrorMessage("Value %d for parameter fastfail out of range {0, 1}\n", ival);
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

/** gets floating point parameter of NLP
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - dval pointer to store the parameter value
 *
 *  output:
 *  - dval parameter value
 */
static
SCIP_DECL_NLPIGETREALPAR(nlpiGetRealParIpopt)
{
   assert(nlpi != NULL);
   assert(dval != NULL);
   assert(problem != NULL);
   assert(IsValid(problem->ipopt));

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
      (void) problem->ipopt->Options()->GetNumericValue("acceptable_constr_viol_tol", *dval, "");
      break;
   }

   case SCIP_NLPPAR_RELOBJTOL:
   {
      (void) problem->ipopt->Options()->GetNumericValue("dual_inf_tol", *dval, "");
      break;
   }

   case SCIP_NLPPAR_LOBJLIM:
   {
      *dval = -SCIPinfinity(scip);
      break;
   }

   case SCIP_NLPPAR_ITLIM:
   {
      SCIPerrorMessage("iteration limit parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_TILIM:
   {
      (void) problem->ipopt->Options()->GetNumericValue("max_cpu_time", *dval, "");
      break;
   }

   case SCIP_NLPPAR_OPTFILE:
   {
      SCIPerrorMessage("option file parameter is of type string.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FASTFAIL:
   {
      SCIPerrorMessage("fastfail parameter is of type int.\n");
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
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - dval parameter value
 */
static
SCIP_DECL_NLPISETREALPAR(nlpiSetRealParIpopt)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(IsValid(problem->ipopt));

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
         setFeastol(problem, dval);
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
         setOpttol(problem, dval);
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
      SCIPwarningMessage(scip, "Parameter lower objective limit not supported by Ipopt interface yet. Ignored.\n");
      break;
   }

   case SCIP_NLPPAR_ITLIM:
   {
      SCIPerrorMessage("iteration limit parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_TILIM:
   {
      if( dval >= 0.0 )
      {
         /* Ipopt doesn't like a setting of exactly 0 for the max_cpu_time, so increase as little as possible in that case */
         (void) problem->ipopt->Options()->SetNumericValue("max_cpu_time", MAX(dval, DBL_MIN));
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

   case SCIP_NLPPAR_FASTFAIL:
   {
      SCIPerrorMessage("fastfail parameter is of type int.\n");
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
 *  input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - sval pointer to store the string value, the user must not modify the string
 *
 *  output:
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

   case SCIP_NLPPAR_FASTFAIL:
   {
      SCIPerrorMessage("fastfail parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   default:
   {
      SCIPerrorMessage("Parameter %d not known to Ipopt interface.\n", type);
      return SCIP_PARAMETERUNKNOWN;
   }
   }

   return SCIP_OKAY; /*lint !e527*/
}

/** sets string parameter of NLP
 *
 *  input:
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
      (void) problem->ipopt->Options()->GetBoolValue("store_intermediate", problem->storeintermediate, "");
      problem->firstrun = TRUE;

      return SCIP_OKAY;
   }

   case SCIP_NLPPAR_FASTFAIL:
   {
      SCIPerrorMessage("fastfail parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   default:
   {
      SCIPerrorMessage("Parameter %d not known to Ipopt interface.\n", type);
      return SCIP_PARAMETERUNKNOWN;
   }
   }

   return SCIP_OKAY; /*lint !e527*/
}

/** create solver interface for Ipopt solver and includes it into SCIP, if Ipopt is available */
SCIP_RETCODE SCIPincludeNlpSolverIpopt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLPIDATA* nlpidata;

   assert(scip != NULL);

   SCIP_ALLOC( nlpidata = new SCIP_NLPIDATA() ); /*lint !e774*/

   SCIP_CALL( SCIPincludeNlpi(scip,
         NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
         nlpiCopyIpopt, nlpiFreeIpopt, nlpiGetSolverPointerIpopt,
         nlpiCreateProblemIpopt, nlpiFreeProblemIpopt, nlpiGetProblemPointerIpopt,
         nlpiAddVarsIpopt, nlpiAddConstraintsIpopt, nlpiSetObjectiveIpopt,
         nlpiChgVarBoundsIpopt, nlpiChgConsSidesIpopt, nlpiDelVarSetIpopt, nlpiDelConstraintSetIpopt,
         nlpiChgLinearCoefsIpopt, nlpiChgExprIpopt,
         nlpiChgObjConstantIpopt, nlpiSetInitialGuessIpopt, nlpiSolveIpopt, nlpiGetSolstatIpopt, nlpiGetTermstatIpopt,
         nlpiGetSolutionIpopt, nlpiGetStatisticsIpopt,
         nlpiGetWarmstartSizeIpopt, nlpiGetWarmstartMemoIpopt, nlpiSetWarmstartMemoIpopt,
         nlpiGetIntParIpopt, nlpiSetIntParIpopt, nlpiGetRealParIpopt, nlpiSetRealParIpopt,
         nlpiGetStringParIpopt, nlpiSetStringParIpopt,
         nlpidata) );

   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameIpopt(), SCIPgetSolverDescIpopt()) );

   return SCIP_OKAY;
}  /*lint !e429 */

/** gets string that identifies Ipopt (version number) */
const char* SCIPgetSolverNameIpopt(void)
{
   return "Ipopt " IPOPT_VERSION;
}

/** gets string that describes Ipopt */
const char* SCIPgetSolverDescIpopt(void)
{
   return "Interior Point Optimizer developed by A. Waechter et.al. (www.coin-or.org/Ipopt)";
}

/** returns whether Ipopt is available, i.e., whether it has been linked in */
SCIP_Bool SCIPisIpoptAvailableIpopt(void)
{
   return TRUE;
}

/** gives a pointer to the IpoptApplication object stored in Ipopt-NLPI's NLPI problem data structure */
void* SCIPgetIpoptApplicationPointerIpopt(
   SCIP_NLPIPROBLEM*     nlpiproblem         /**< NLP problem of Ipopt-NLPI */
   )
{
   assert(nlpiproblem != NULL);

   return (void*)GetRawPtr(nlpiproblem->ipopt);
}

/** gives a pointer to the NLPIORACLE object stored in Ipopt-NLPI's NLPI problem data structure */
void* SCIPgetNlpiOracleIpopt(
   SCIP_NLPIPROBLEM*     nlpiproblem         /**< NLP problem of Ipopt-NLPI */
   )
{
   assert(nlpiproblem != NULL);

   return nlpiproblem->oracle;
}

/** sets modified default settings that are used when setting up an Ipopt problem
 *
 *  Do not forget to add a newline after the last option in optionsstring.
 */
void SCIPsetModifiedDefaultSettingsIpopt(
   SCIP_NLPI*            nlpi,               /**< Ipopt NLP interface */
   const char*           optionsstring,      /**< string with options as in Ipopt options file */
   SCIP_Bool             append              /**< whether to append to modified default settings or to overwrite */
   )
{
   SCIP_NLPIDATA* data;

   assert(nlpi != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   if( append )
      data->defoptions += optionsstring;
   else
      data->defoptions = optionsstring;
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

   retcode = SCIPnlpiOracleGetJacobianSparsity(scip, nlpiproblem->oracle, &offset, NULL);
   if( retcode != SCIP_OKAY )
      return false;
   assert(offset != NULL);
   nnz_jac_g = offset[m];

   if( !approxhessian )
   {
      retcode = SCIPnlpiOracleGetHessianLagSparsity(scip, nlpiproblem->oracle, &offset, NULL);
      if( retcode != SCIP_OKAY )
         return false;
      assert(offset != NULL);
      nnz_h_lag = offset[n];
   }
   else
   {
      nnz_h_lag = 0;
   }

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

   assert(SCIPnlpiOracleGetVarLbs(nlpiproblem->oracle) != NULL || n == 0);
   assert(SCIPnlpiOracleGetVarUbs(nlpiproblem->oracle) != NULL || n == 0);

   BMScopyMemoryArray(x_l, SCIPnlpiOracleGetVarLbs(nlpiproblem->oracle), n);
   BMScopyMemoryArray(x_u, SCIPnlpiOracleGetVarUbs(nlpiproblem->oracle), n);
#ifndef NDEBUG
   for( int i = 0; i < n; ++i )
      assert(x_l[i] <= x_u[i]);
#endif

   /* Ipopt performs better when unused variables do not appear, which we can achieve by fixing them,
    * since Ipopts TNLPAdapter will hide them from Ipopts NLP. In the dual solution, bound multipliers (z_L, z_U)
    * for fixed variables have value 0.0.
    */
   for( int i = 0; i < n; ++i )
   {
      int vardegree;
      if( SCIPnlpiOracleGetVarDegree(scip, nlpiproblem->oracle, i, &vardegree) != SCIP_OKAY )
         return false;
      if( vardegree == 0 )
      {
         SCIPdebugMsg(scip, "fix unused variable x%d [%g,%g] to 0.0 or bound\n", i, x_l[i], x_u[i]);
         assert(x_l[i] <= x_u[i]);
         x_l[i] = x_u[i] = MAX(MIN(x_u[i], 0.0), x_l[i]);
      }
   }

   for( int i = 0; i < m; ++i )
   {
      g_l[i] = SCIPnlpiOracleGetConstraintLhs(nlpiproblem->oracle, i);
      g_u[i] = SCIPnlpiOracleGetConstraintRhs(nlpiproblem->oracle, i);
      assert(g_l[i] <= g_u[i]);
   }

   return true;
}

/** Method to return the starting point for the algorithm */  /*lint -e{715}*/
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
{  /*lint --e{715} */
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

         SCIPdebugMsg(scip, "Ipopt started without initial primal values; make up starting guess by projecting 0 onto variable bounds\n");

         for( int i = 0; i < n; ++i )
         {
            lb = SCIPnlpiOracleGetVarLbs(nlpiproblem->oracle)[i];
            ub = SCIPnlpiOracleGetVarUbs(nlpiproblem->oracle)[i];
            if( lb > 0.0 )
               x[i] = SCIPrandomGetReal(randnumgen, lb, lb + MAXPERTURB*MIN(1.0, ub-lb));
            else if( ub < 0.0 )
               x[i] = SCIPrandomGetReal(randnumgen, ub - MAXPERTURB*MIN(1.0, ub-lb), ub);
            else
               x[i] = SCIPrandomGetReal(randnumgen, MAX(lb, -MAXPERTURB*MIN(1.0, ub-lb)), MIN(ub, MAXPERTURB*MIN(1.0, ub-lb)));
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
   {
      int vardegree;
      if( SCIPnlpiOracleGetVarDegree(scip, nlpiproblem->oracle, i, &vardegree) != SCIP_OKAY )
         return false;
      var_types[i] = vardegree <= 1 ? LINEAR : NON_LINEAR;
   }

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

/** Method to return the number of nonlinear variables. */
Index ScipNLP::get_number_of_nonlinear_variables()
{
   int count;
   int n;

   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);

   n = SCIPnlpiOracleGetNVars(nlpiproblem->oracle);

   count = 0;
   for( int i = 0; i < n; ++i )
   {
      int vardegree;
      if( SCIPnlpiOracleGetVarDegree(scip, nlpiproblem->oracle, i, &vardegree) != FALSE )
         return -1;  // TODO what is the right way to report an error back to Ipopt here?
      if( vardegree > 1 )
         ++count;
   }

   return count;
}

/** Method to return the indices of the nonlinear variables */
bool ScipNLP::get_list_of_nonlinear_variables(
   Index              num_nonlin_vars,    /**< number of nonlinear variables */
   Index*             pos_nonlin_vars     /**< array to fill with indices of nonlinear variables */
   )
{
   int count;
   int n;

   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);

   n = SCIPnlpiOracleGetNVars(nlpiproblem->oracle);

   count = 0;
   for( int i = 0; i < n; ++i )
   {
      int vardegree;
      if( SCIPnlpiOracleGetVarDegree(scip, nlpiproblem->oracle, i, &vardegree) != SCIP_OKAY )
         return false;
      if( vardegree > 1 )
      {
         assert(count < num_nonlin_vars);
         pos_nonlin_vars[count++] = i;
      }
   }

   assert(count == num_nonlin_vars);

   return true;
}

/** Method to return metadata about variables and constraints */  /*lint -e{715}*/
bool ScipNLP::get_var_con_metadata(
   Index              n,                  /**< number of variables */
   StringMetaDataMapType& var_string_md,  /**< variable meta data of string type */
   IntegerMetaDataMapType& var_integer_md,/**< variable meta data of integer type */
   NumericMetaDataMapType& var_numeric_md,/**< variable meta data of numeric type */
   Index              m,                  /**< number of constraints */
   StringMetaDataMapType& con_string_md,  /**< constraint meta data of string type */
   IntegerMetaDataMapType& con_integer_md,/**< constraint meta data of integer type */
   NumericMetaDataMapType& con_numeric_md /**< constraint meta data of numeric type */
   )
{ /*lint --e{715}*/
   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);
   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));
   assert(m == SCIPnlpiOracleGetNConstraints(nlpiproblem->oracle));

   char** varnames = SCIPnlpiOracleGetVarNames(nlpiproblem->oracle);
   if( varnames != NULL )
   {
      std::vector<std::string>& varnamesvec(var_string_md["idx_names"]);
      varnamesvec.reserve((size_t)n);
      for( int i = 0; i < n; ++i )
      {
         if( varnames[i] != NULL )
         {
            varnamesvec.push_back(varnames[i]);  /*lint !e3701*/
         }
         else
         {
            char buffer[20];
            (void) sprintf(buffer, "nlpivar%8d", i);
            varnamesvec.push_back(buffer);
         }
      }
   }

   std::vector<std::string>& consnamesvec(con_string_md["idx_names"]);
   consnamesvec.reserve((size_t)m);
   for( int i = 0; i < m; ++i )
   {
      if( SCIPnlpiOracleGetConstraintName(nlpiproblem->oracle, i) != NULL )
      {
         consnamesvec.push_back(SCIPnlpiOracleGetConstraintName(nlpiproblem->oracle, i));
      }
      else
      {
         char buffer[20];
         (void) sprintf(buffer, "nlpicons%8d", i);
         consnamesvec.push_back(buffer);  /*lint !e3701*/
      }
   }

   return true;
}

/** Method to return the objective value */  /*lint -e{715}*/
bool ScipNLP::eval_f(
   Index              n,                  /**< number of variables */ 
   const Number*      x,                  /**< point to evaluate */ 
   bool               new_x,              /**< whether some function evaluation method has been called for this point before */
   Number&            obj_value           /**< place to store objective function value */
   )
{ /*lint --e{715}*/
   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);

   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));

   return SCIPnlpiOracleEvalObjectiveValue(scip, nlpiproblem->oracle, x, &obj_value) == SCIP_OKAY;
}

/** Method to return the gradient of the objective */  /*lint -e{715}*/
bool ScipNLP::eval_grad_f(
   Index              n,                  /**< number of variables */ 
   const Number*      x,                  /**< point to evaluate */ 
   bool               new_x,              /**< whether some function evaluation method has been called for this point before */
   Number*            grad_f              /**< buffer to store objective gradient */
   )
{ /*lint --e{715}*/
   SCIP_Real dummy;

   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);

   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));

   return SCIPnlpiOracleEvalObjectiveGradient(scip, nlpiproblem->oracle, x, TRUE, &dummy, grad_f) == SCIP_OKAY;
}

/** Method to return the constraint residuals */  /*lint -e{715}*/
bool ScipNLP::eval_g(
   Index              n,                  /**< number of variables */ 
   const Number*      x,                  /**< point to evaluate */ 
   bool               new_x,              /**< whether some function evaluation method has been called for this point before */
   Index              m,                  /**< number of constraints */
   Number*            g                   /**< buffer to store constraint function values */
   )
{ /*lint --e{715}*/
   assert(nlpiproblem != NULL);
   assert(nlpiproblem->oracle != NULL);

   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));

   return SCIPnlpiOracleEvalConstraintValues(scip, nlpiproblem->oracle, x, g) == SCIP_OKAY;
}

/** Method to return:
 *   1) The structure of the jacobian (if "values" is NULL)
 *   2) The values of the jacobian (if "values" is not NULL)
 */  /*lint -e{715}*/
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
{ /*lint --e{715}*/
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

      if( SCIPnlpiOracleGetJacobianSparsity(scip, nlpiproblem->oracle, &jacoffset, &jaccol) != SCIP_OKAY )
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
      if( SCIPnlpiOracleEvalJacobian(scip, nlpiproblem->oracle, x, TRUE, NULL, values) != SCIP_OKAY )
         return false;
   }

   return true;
}

/** Method to return:
 *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
 *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
 */   /*lint -e{715}*/
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
{  /*lint --e{715}*/
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

      if( SCIPnlpiOracleGetHessianLagSparsity(scip, nlpiproblem->oracle, &heslagoffset, &heslagcol) != SCIP_OKAY )
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
      if( SCIPnlpiOracleEvalHessianLag(scip, nlpiproblem->oracle, x, TRUE, obj_factor, lambda, values) != SCIP_OKAY )
         return false;
   }

   return true;
}

/** Method called by the solver at each iteration.
 * 
 * Checks whether Ctrl-C was hit.
 */   /*lint -e{715}*/
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
{  /*lint --e{715}*/
   if( nlpiproblem->storeintermediate && mode == RegularMode && inf_pr < nlpiproblem->lastsolinfeas )
   {
      Ipopt::TNLPAdapter* tnlp_adapter;

      tnlp_adapter = NULL;
      if( ip_cq != NULL )
      {
         Ipopt::OrigIpoptNLP* orignlp;

         orignlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
         if( orignlp != NULL )
            tnlp_adapter = dynamic_cast<TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
      }

      if( tnlp_adapter != NULL && ip_data != NULL && IsValid(ip_data->curr()) )
      {
         SCIPdebugMsg(scip, "update lastsol: inf_pr old = %g -> new = %g\n", nlpiproblem->lastsolinfeas, inf_pr);

         if( nlpiproblem->lastsolprimals == NULL )
         {
            assert(nlpiproblem->lastsoldualcons == NULL);
            assert(nlpiproblem->lastsoldualvarlb == NULL);
            assert(nlpiproblem->lastsoldualvarub == NULL);
            if( BMSallocMemoryArray(&nlpiproblem->lastsolprimals, SCIPnlpiOracleGetNVars(nlpiproblem->oracle)) == NULL ||
               BMSallocMemoryArray(&nlpiproblem->lastsoldualcons, SCIPnlpiOracleGetNConstraints(nlpiproblem->oracle)) == NULL ||
               BMSallocMemoryArray(&nlpiproblem->lastsoldualvarlb, SCIPnlpiOracleGetNVars(nlpiproblem->oracle)) == NULL ||
               BMSallocMemoryArray(&nlpiproblem->lastsoldualvarub, SCIPnlpiOracleGetNVars(nlpiproblem->oracle)) == NULL )
            {
               SCIPerrorMessage("out-of-memory in ScipNLP::intermediate_callback()\n");
               return true;
            }
         }

         assert(IsValid(ip_data->curr()->x()));
         tnlp_adapter->ResortX(*ip_data->curr()->x(), nlpiproblem->lastsolprimals);
         nlpiproblem->lastsolinfeas = inf_pr;

         assert(IsValid(ip_data->curr()->y_c()));
         assert(IsValid(ip_data->curr()->y_d()));
         tnlp_adapter->ResortG(*ip_data->curr()->y_c(), *ip_data->curr()->y_d(), nlpiproblem->lastsoldualcons);

         // need to clear arrays first because ResortBnds only sets values for non-fixed variables
         BMSclearMemoryArray(nlpiproblem->lastsoldualvarlb, SCIPnlpiOracleGetNVars(nlpiproblem->oracle));
         BMSclearMemoryArray(nlpiproblem->lastsoldualvarub, SCIPnlpiOracleGetNVars(nlpiproblem->oracle));
         assert(IsValid(ip_data->curr()->z_L()));
         assert(IsValid(ip_data->curr()->z_U()));
         tnlp_adapter->ResortBnds(*ip_data->curr()->z_L(), nlpiproblem->lastsoldualvarlb, *ip_data->curr()->z_U(), nlpiproblem->lastsoldualvarub);

      }
   }

   /* do convergence test if fastfail is enabled */
   if( nlpiproblem->fastfail )
   {
      int i;

      if( iter == 0 )
      {
         conv_lastrestoiter = -1;
      }
      else if( mode == RestorationPhaseMode )
      {
         conv_lastrestoiter = iter;
      }
      else if( conv_lastrestoiter == iter-1 )
      {
         /* just switched back from restoration mode, reset dual reduction targets */
         for( i = 0; i < convcheck_nchecks; ++i )
            conv_dutarget[i] = convcheck_minred[i] * inf_du;
      }

      if( iter == convcheck_startiter )
      {
         /* define initial targets and iteration limits */
         for( i = 0; i < convcheck_nchecks; ++i )
         {
            conv_prtarget[i] = convcheck_minred[i] * inf_pr;
            conv_dutarget[i] = convcheck_minred[i] * inf_du;
            conv_iterlim[i] = iter + convcheck_maxiter[i];
         }
      }
      else if( iter > convcheck_startiter )
      {
         /* check if we should stop */
         for( i = 0; i < convcheck_nchecks; ++i )
         {
            if( inf_pr <= conv_prtarget[i] )
            {
               /* sufficient reduction w.r.t. primal infeasibility target
                * reset target w.r.t. current infeasibilities
                */
               conv_prtarget[i] = convcheck_minred[i] * inf_pr;
               conv_dutarget[i] = convcheck_minred[i] * inf_du;
               conv_iterlim[i] = iter + convcheck_maxiter[i];
            }
            else if( iter >= conv_iterlim[i] )
            {
               /* we hit a limit, should we really stop? */
               SCIPdebugMsg(scip, "convcheck %d: inf_pr = %e > target %e; inf_du = %e target %e: ",
                  i, inf_pr, conv_prtarget[i], inf_du, conv_dutarget[i]);
               if( mode == RegularMode && iter <= conv_lastrestoiter + convcheck_startiter )
               {
                  /* if we returned from feasibility restoration recently, we allow some more iterations,
                   * because Ipopt may go for optimality for some iterations, at the costs of infeasibility
                   */
                  SCIPdebugPrintf("continue, because restoration phase only %d iters ago\n", iter - conv_lastrestoiter);
               }
               else if( mode == RegularMode && inf_du <= conv_dutarget[i] && iter < conv_iterlim[i] + convcheck_maxiter[i] )
               {
                  /* if dual reduction is sufficient, we allow for twice the number of iterations to reach primal infeas reduction */
                  SCIPdebugPrintf("continue, because dual infeas. red. sufficient and only %d iters above limit\n", iter - conv_iterlim[i]);
               }
               else
               {
                  SCIPdebugPrintf("abort\n");
                  return false;
               }
            }
         }
      }
   }

   return (SCIPinterrupted() == FALSE);
}

/** This method is called when the algorithm is complete so the TNLP can store/write the solution. */  /*lint -e{715}*/
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
{ /*lint --e{715}*/
   assert(nlpiproblem         != NULL);
   assert(nlpiproblem->oracle != NULL);

   assert(n == SCIPnlpiOracleGetNVars(nlpiproblem->oracle));
   assert(m == SCIPnlpiOracleGetNConstraints(nlpiproblem->oracle));

   bool check_feasibility = false; // whether we should check x for feasibility, if not NULL
   switch( status )
   {
   case SUCCESS:
      nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_LOCOPT;
      nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      assert(x != NULL);
      break;

   case STOP_AT_ACCEPTABLE_POINT:
      /* if stop at acceptable point, then dual infeasibility can be arbitrary large, so claim only feasibility */
   case FEASIBLE_POINT_FOUND:
      nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
      nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      assert(x != NULL);
      break;

   case MAXITER_EXCEEDED:
      check_feasibility = true;
      nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_ITLIM;
      break;

   case CPUTIME_EXCEEDED:
      check_feasibility = true;
      nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_TILIM;
      break;

   case STOP_AT_TINY_STEP:
   case RESTORATION_FAILURE:
   case ERROR_IN_STEP_COMPUTATION:
      check_feasibility = true;
      nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_NUMERR;
      break;

   case LOCAL_INFEASIBILITY:
      /* still check feasibility, since we let Ipopt solve with higher tolerance than actually required */
      check_feasibility = true;
      nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
      nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;

   case DIVERGING_ITERATES:
      nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNBOUNDED;
      nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;

   case INVALID_NUMBER_DETECTED:
      nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_EVALERR;
      break;

   case USER_REQUESTED_STOP:
      check_feasibility = true;
      nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;

   case TOO_FEW_DEGREES_OF_FREEDOM:
   case INTERNAL_ERROR:
   case INVALID_OPTION:
      nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
      break;

   case OUT_OF_MEMORY:
      nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_MEMERR;
      break;

   default:
      SCIPerrorMessage("Ipopt returned with unknown solution status %d\n", status);
      nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
      break;
   }

   /* if Ipopt reports its solution as locally infeasible or we don't know feasibility, then report the intermediate point with lowest constraint violation, if available */
   if( (x == NULL || nlpiproblem->lastsolstat == SCIP_NLPSOLSTAT_LOCINFEASIBLE || nlpiproblem->lastsolstat == SCIP_NLPSOLSTAT_UNKNOWN) && nlpiproblem->lastsolinfeas != SCIP_INVALID )  /*lint !e777*/
   {
      /* if infeasibility of lastsol is not invalid, then lastsol values should exist */
      assert(nlpiproblem->lastsolprimals != NULL);
      assert(nlpiproblem->lastsoldualcons != NULL);
      assert(nlpiproblem->lastsoldualvarlb != NULL);
      assert(nlpiproblem->lastsoldualvarub != NULL);

      /* check if lastsol is feasible */
      Number constrvioltol;
      (void) nlpiproblem->ipopt->Options()->GetNumericValue("acceptable_constr_viol_tol", constrvioltol, "");
      if( nlpiproblem->lastsolinfeas <= constrvioltol )
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
      else
         nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;

      SCIPdebugMsg(scip, "drop Ipopt's final point and report intermediate locally %sfeasible solution with infeas %g instead (acceptable: %g)\n",
         nlpiproblem->lastsolstat == SCIP_NLPSOLSTAT_LOCINFEASIBLE ? "in" : "", nlpiproblem->lastsolinfeas, constrvioltol);
   }
   else
   {
      assert(x != NULL);
      assert(lambda != NULL);
      assert(z_L != NULL);
      assert(z_U != NULL);

      if( nlpiproblem->lastsolprimals == NULL )
      {
         assert(nlpiproblem->lastsoldualcons == NULL);
         assert(nlpiproblem->lastsoldualvarlb == NULL);
         assert(nlpiproblem->lastsoldualvarub == NULL);
         BMSallocMemoryArray(&nlpiproblem->lastsolprimals,   n);
         BMSallocMemoryArray(&nlpiproblem->lastsoldualcons,  m);
         BMSallocMemoryArray(&nlpiproblem->lastsoldualvarlb, n);
         BMSallocMemoryArray(&nlpiproblem->lastsoldualvarub, n);

         if( nlpiproblem->lastsolprimals == NULL || nlpiproblem->lastsoldualcons == NULL ||
            nlpiproblem->lastsoldualvarlb == NULL || nlpiproblem->lastsoldualvarub == NULL )
         {
            nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
            nlpiproblem->lasttermstat = SCIP_NLPTERMSTAT_MEMERR;
            return;
         }
      }

      BMScopyMemoryArray(nlpiproblem->lastsolprimals, x, n);
      BMScopyMemoryArray(nlpiproblem->lastsoldualcons, lambda, m);
      BMScopyMemoryArray(nlpiproblem->lastsoldualvarlb, z_L, n);
      BMScopyMemoryArray(nlpiproblem->lastsoldualvarub, z_U, n);

      if( check_feasibility && cq != NULL )
      {
         Number constrviol;
         Number constrvioltol;

         constrviol = cq->curr_constraint_violation();

         (void) nlpiproblem->ipopt->Options()->GetNumericValue("acceptable_constr_viol_tol", constrvioltol, "");
         if( constrviol <= constrvioltol )
            nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
         else if( nlpiproblem->lastsolstat != SCIP_NLPSOLSTAT_LOCINFEASIBLE )
            nlpiproblem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      }
   }

   if( nlpiproblem->lastsolstat == SCIP_NLPSOLSTAT_LOCINFEASIBLE )
   {
      assert(lambda != NULL);
      SCIP_Real tol;
      (void) nlpiproblem->ipopt->Options()->GetNumericValue("tol", tol, "");

      // Jakobs paper ZR_20-20 says we should have lambda*g(x) + mu*h(x) > 0
      //   if the NLP is min f(x) s.t. g(x) <= 0, h(x) = 0
      // we check this here and change solution status to unknown if the test fails
      bool infreasonable = true;
      SCIP_Real infproof = 0.0;
      for( int i = 0; i < m && infreasonable; ++i )
      {
         if( fabs(lambda[i]) < tol )
            continue;
         SCIP_Real side;
         if( lambda[i] < 0.0 )
         {
            // lhs <= g(x) should be active
            // in the NLP above, this should be lhs - g(x) <= 0 with negated dual
            // so this contributes -lambda*(lhs-g(x)) = lambda*(g(x)-side)
            side = SCIPnlpiOracleGetConstraintLhs(nlpiproblem->oracle, i);
            if( SCIPisInfinity(scip, -side) )
            {
               SCIPdebugMessage("inconsistent dual, lambda = %g, but lhs = %g\n", lambda[i], side);
               infreasonable = false;
            }
         }
         else
         {
            // g(x) <= rhs should be active
            // in the NLP above, this should be g(x) - rhs <= 0
            // so this contributes lambda*(g(x)-rhs)
            side = SCIPnlpiOracleGetConstraintRhs(nlpiproblem->oracle, i);
            if( SCIPisInfinity(scip, side) )
            {
               SCIPdebugMessage("inconsistent dual, lambda = %g, but rhs = %g\n", lambda[i], side);
               infreasonable = false;
            }
         }

         // g(x) <= 0
         infproof += lambda[i] * (g[i] - side);
         // SCIPdebugMessage("cons %d lambda %g, slack %g\n", i, lambda[i], g[i] - side);
      }
      if( infreasonable )
      {
         SCIPdebugMessage("infproof = %g should be positive to be valid\n", infproof);
         if( infproof <= 0.0 )
            infreasonable = false;
      }

      if( !infreasonable )
      {
         // change status to say we don't know
         nlpiproblem->lastsolstat = SCIP_NLPSOLSTAT_UNKNOWN;
      }
   }
}

/** Calls Lapacks Dsyev routine to compute eigenvalues and eigenvectors of a dense matrix.
 *
 *  It's here, because we use Ipopt's interface to Lapack.
 */
SCIP_RETCODE LapackDsyev(
   SCIP_Bool             computeeigenvectors,/**< should also eigenvectors should be computed ? */
   int                   N,                  /**< dimension */
   SCIP_Real*            a,                  /**< matrix data on input (size N*N); eigenvectors on output if computeeigenvectors == TRUE */
   SCIP_Real*            w                   /**< buffer to store eigenvalues (size N) */
   )
{
   int info;

   IpLapackDsyev((bool)computeeigenvectors, N, a, N, w, info);

   if( info != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEV. INFO = %d\n", info);
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** solves a linear problem of the form Ax = b for a regular matrix 3*3 A */
static
SCIP_RETCODE SCIPsolveLinearProb3(
   SCIP_Real*            A,                  /**< matrix data on input (size 3*3); filled column-wise */
   SCIP_Real*            b,                  /**< right hand side vector (size 3) */
   SCIP_Real*            x,                  /**< buffer to store solution (size 3) */
   SCIP_Bool*            success             /**< pointer to store if the solving routine was successful */
   )
{
   SCIP_Real Acopy[9];
   SCIP_Real bcopy[3];
   int pivotcopy[3];
   const int N = 3;
   int info;

   assert(A != NULL);
   assert(b != NULL);
   assert(x != NULL);
   assert(success != NULL);

   BMScopyMemoryArray(Acopy, A, N*N);
   BMScopyMemoryArray(bcopy, b, N);

   /* compute the LU factorization */
   IpLapackDgetrf(N, Acopy, pivotcopy, N, info);

   if( info != 0 )
   {
      SCIPdebugMessage("There was an error when calling Dgetrf. INFO = %d\n", info);
      *success = FALSE;
   }
   else
   {
      *success = TRUE;

      /* solve linear problem */
      IpLapackDgetrs(N, 1, Acopy, N, pivotcopy, bcopy, N);

      /* copy the solution */
      BMScopyMemoryArray(x, bcopy, N);
   }

   return SCIP_OKAY;
}

/** solves a linear problem of the form Ax = b for a regular matrix A
 *
 *  Calls Lapacks IpLapackDgetrf routine to calculate a LU factorization and uses this factorization to solve
 *  the linear problem Ax = b.
 *  It's here, because Ipopt is linked against Lapack.
 */
SCIP_RETCODE SCIPsolveLinearProb(
   int                   N,                  /**< dimension */
   SCIP_Real*            A,                  /**< matrix data on input (size N*N); filled column-wise */
   SCIP_Real*            b,                  /**< right hand side vector (size N) */
   SCIP_Real*            x,                  /**< buffer to store solution (size N) */
   SCIP_Bool*            success             /**< pointer to store if the solving routine was successful */
   )
{
   SCIP_Real* Acopy;
   SCIP_Real* bcopy;
   int* pivotcopy;
   int info;

   assert(N > 0);
   assert(A != NULL);
   assert(b != NULL);
   assert(x != NULL);
   assert(success != NULL);

   /* call SCIPsolveLinearProb3() for performance reasons */
   if( N == 3 )
   {
      SCIP_CALL( SCIPsolveLinearProb3(A, b, x, success) );
      return SCIP_OKAY;
   }

   Acopy = NULL;
   bcopy = NULL;
   pivotcopy = NULL;

   SCIP_ALLOC( BMSduplicateMemoryArray(&Acopy, A, N*N) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&bcopy, b, N) );
   SCIP_ALLOC( BMSallocMemoryArray(&pivotcopy, N) );

   /* compute the LU factorization */
   IpLapackDgetrf(N, Acopy, pivotcopy, N, info);

   if( info != 0 )
   {
      SCIPdebugMessage("There was an error when calling Dgetrf. INFO = %d\n", info);
      *success = FALSE;
   }
   else
   {
      *success = TRUE;

      /* solve linear problem */
      IpLapackDgetrs(N, 1, Acopy, N, pivotcopy, bcopy, N);

      /* copy the solution */
      BMScopyMemoryArray(x, bcopy, N);
   }

   BMSfreeMemoryArray(&pivotcopy);
   BMSfreeMemoryArray(&bcopy);
   BMSfreeMemoryArray(&Acopy);

   return SCIP_OKAY;
}
