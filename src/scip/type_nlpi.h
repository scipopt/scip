/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_nlpi.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for NLP solver interfaces
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_NLPI_H__
#define __SCIP_TYPE_NLPI_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_expr.h"
#include "scip/type_nlp.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Nlpi          SCIP_NLPI;          /**< NLP solver interface */
typedef struct SCIP_NlpiData      SCIP_NLPIDATA;      /**< locally defined NLP solver interface data */
typedef struct SCIP_NlpiProblem   SCIP_NLPIPROBLEM;   /**< locally defined NLP solver interface data for a specific problem instance */

/** NLP solver parameter */
enum SCIP_NlpParam
{
   SCIP_NLPPAR_FROMSCRATCH    =  0,      /**< solver should start from scratch at next call?: 0 no, 1 yes (int) */
   SCIP_NLPPAR_VERBLEVEL      =  1,      /**< verbosity level of output of NLP solver to the screen: 0 off, 1 normal, 2 debug, > 2 more debug (int) */
   SCIP_NLPPAR_FEASTOL        =  2,      /**< feasibility tolerance for primal variables and slacks (real) */
   SCIP_NLPPAR_RELOBJTOL      =  3,      /**< relative objective tolerance (real) */
   SCIP_NLPPAR_LOBJLIM        =  4,      /**< lower objective limit (cutoff) (real) */
   SCIP_NLPPAR_ITLIM          =  5,      /**< NLP iteration limit (int) */
   SCIP_NLPPAR_TILIM          =  6,      /**< NLP time limit (real) */
   SCIP_NLPPAR_OPTFILE        =  7,      /**< name of a solver specific option file (string) */
   SCIP_NLPPAR_FASTFAIL       =  8       /**< should the NLP solver stop early if convergence is slow?: 0 no, 1 yes (int) */
};
typedef enum SCIP_NlpParam SCIP_NLPPARAM;  /**< NLP solver parameter */

/** NLP solution status */
enum SCIP_NlpSolStat
{
   SCIP_NLPSOLSTAT_GLOBOPT        = 0,    /**< solved to global optimality */
   SCIP_NLPSOLSTAT_LOCOPT         = 1,    /**< solved to local optimality */
   SCIP_NLPSOLSTAT_FEASIBLE       = 2,    /**< feasible solution found */
   SCIP_NLPSOLSTAT_LOCINFEASIBLE  = 3,    /**< solution found is local infeasible */
   SCIP_NLPSOLSTAT_GLOBINFEASIBLE = 4,    /**< problem is proven infeasible */
   SCIP_NLPSOLSTAT_UNBOUNDED      = 5,    /**< problem is unbounded */
   SCIP_NLPSOLSTAT_UNKNOWN        = 6     /**< unknown solution status (e.g., problem not solved yet) */
};
typedef enum SCIP_NlpSolStat SCIP_NLPSOLSTAT;      /** NLP solution status */

/** NLP solver termination status */
enum SCIP_NlpTermStat
{
   SCIP_NLPTERMSTAT_OKAY          = 0,    /**< terminated successfully */
   SCIP_NLPTERMSTAT_TILIM         = 1,    /**< time limit exceeded */
   SCIP_NLPTERMSTAT_ITLIM         = 2,    /**< iteration limit exceeded */
   SCIP_NLPTERMSTAT_LOBJLIM       = 3,    /**< lower objective limit reached */
   SCIP_NLPTERMSTAT_NUMERR        = 5,    /**< stopped on numerical error */
   SCIP_NLPTERMSTAT_EVALERR       = 6,    /**< stopped on function evaluation error */
   SCIP_NLPTERMSTAT_MEMERR        = 7,    /**< memory exceeded */
   SCIP_NLPTERMSTAT_LICERR        = 8,    /**< licence error */
   SCIP_NLPTERMSTAT_OTHER         = 9     /**< other error (= this should never happen) */
};
typedef enum SCIP_NlpTermStat SCIP_NLPTERMSTAT;  /** NLP solver termination status */

/** copy method of NLP interface (called when SCIP copies plugins)
 *
 * input:
 *  - scip target SCIP where to include copy of NLPI
 *  - sourcenlpi the NLP interface to copy
 */
#define SCIP_DECL_NLPICOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* sourcenlpi)

/** destructor of NLP interface to free nlpi data
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - nlpidata NLPI data to free
 */
#define SCIP_DECL_NLPIFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIDATA** nlpidata)

/** gets pointer to solver-internal NLP solver
 * 
 *  to do dirty stuff
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  
 * return: void pointer to solver
 */
#define SCIP_DECL_NLPIGETSOLVERPOINTER(x) void* x (SCIP* scip, SCIP_NLPI* nlpi)

/** creates a problem instance
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem pointer to store the problem data
 *  - name name of problem, can be NULL
 */
#define SCIP_DECL_NLPICREATEPROBLEM(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM** problem, const char* name)

/** free a problem instance
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem pointer where problem data is stored 
 */
#define SCIP_DECL_NLPIFREEPROBLEM(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM** problem)

/** gets pointer to solver-internal problem instance
 * 
 *  to do dirty stuff
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  
 * return: void pointer to problem instance
 */
#define SCIP_DECL_NLPIGETPROBLEMPOINTER(x) void* x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem)

/** add variables
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nvars number of variables 
 *  - lbs lower bounds of variables, can be NULL if -infinity
 *  - ubs upper bounds of variables, can be NULL if +infinity
 *  - varnames names of variables, can be NULL
 */
#define SCIP_DECL_NLPIADDVARS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, int nvars, const SCIP_Real* lbs, \
      const SCIP_Real* ubs, const char** varnames)

/** add constraints
 * quadratic coefficiens: row oriented matrix for each constraint
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - ncons number of added constraints
 *  - lhss left hand sides of constraints, can be NULL if -infinity
 *  - rhss right hand sides of constraints, can be NULL if +infinity
 *  - nlininds number of linear coefficients for each constraint
 *    may be NULL in case of no linear part
 *  - lininds indices of variables for linear coefficients for each constraint
 *    may be NULL in case of no linear part
 *  - linvals values of linear coefficient for each constraint
 *    may be NULL in case of no linear part
 *  - exprs expressions for nonlinear part of constraints
 *    entry of array may be NULL in case of no nonlinear part
 *    may be NULL in case of no nonlinear part in any constraint
 *  - names of constraints, may be NULL or entries may be NULL
 */
#define SCIP_DECL_NLPIADDCONSTRAINTS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, int nconss, const SCIP_Real* lhss, \
      const SCIP_Real* rhss, const int* nlininds, int* const* lininds, SCIP_Real* const* linvals, \
      SCIP_EXPR** exprs, const char** names)

/** sets or overwrites objective, a minimization problem is expected
 *  May change sparsity pattern.
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nlins number of linear variables
 *  - lininds variable indices
 *    may be NULL in case of no linear part
 *  - linvals coefficient values
 *    may be NULL in case of no linear part
 *  - expr expression for nonlinear part of objective function
 *    may be NULL in case of no nonlinear part
 *  - constant objective value offset
 */
#define SCIP_DECL_NLPISETOBJECTIVE(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, int nlins, const int* lininds, \
      const SCIP_Real* linvals, SCIP_EXPR* expr, \
      const SCIP_Real constant)

/** change variable bounds
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nvars number of variables to change bounds
 *  - indices indices of variables to change bounds
 *  - lbs new lower bounds
 *  - ubs new upper bounds
 */
#define SCIP_DECL_NLPICHGVARBOUNDS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, const int nvars, const int* indices, \
      const SCIP_Real* lbs, const SCIP_Real* ubs)

/** change constraint sides
 *
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nconss number of constraints to change sides
 *  - indices indices of constraints to change sides
 *  - lhss new left hand sides
 *  - rhss new right hand sides
 */
#define SCIP_DECL_NLPICHGCONSSIDES(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, int nconss, const int* indices, \
      const SCIP_Real* lhss, const SCIP_Real* rhss)

/** delete a set of variables
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - dstats deletion status of vars; 1 if var should be deleted, 0 if not
 *  - size of the dstats array
 * 
 * output:
 *  - dstats new position of var, -1 if var was deleted
 */
#define SCIP_DECL_NLPIDELVARSET(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, int* dstats, int dstatssize)

/** delete a set of constraints
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - dstats deletion status of rows; 1 if row should be deleted, 0 if not
 *  - size of the dstats array
 * 
 * output:
 *  - dstats new position of row, -1 if row was deleted
 */
#define SCIP_DECL_NLPIDELCONSSET(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, int* dstats, int dstatssize)

/** changes (or adds) linear coefficients in a constraint or objective
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idx index of constraint or -1 for objective
 *  - nvals number of values in linear constraint to change
 *  - varidxs indices of variables which coefficient to change
 *  - vals new values for coefficients
 */
#define SCIP_DECL_NLPICHGLINEARCOEFS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, int idx, int nvals, \
      const int* varidxs, const SCIP_Real* vals)

/** replaces the expression of a constraint or objective
 *
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idxcons index of constraint or -1 for objective
 *  - expr new expression for constraint or objective, or NULL to only remove previous tree
 */
#define SCIP_DECL_NLPICHGEXPR(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, int idxcons, SCIP_EXPR* expr)

/** change the constant offset in the objective
 *
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - objconstant new value for objective constant
 */
#define SCIP_DECL_NLPICHGOBJCONSTANT(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, SCIP_Real objconstant)

/** sets initial guess for primal variables
 *
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues initial primal values for variables, or NULL to clear previous values
 *  - consdualvalues initial dual values for constraints, or NULL to clear previous values
 *  - varlbdualvalues  initial dual values for variable lower bounds, or NULL to clear previous values
 *  - varubdualvalues  initial dual values for variable upper bounds, or NULL to clear previous values
 */
#define SCIP_DECL_NLPISETINITIALGUESS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, SCIP_Real* primalvalues, \
      SCIP_Real* consdualvalues, SCIP_Real* varlbdualvalues, SCIP_Real* varubdualvalues)

/** tries to solve NLP
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 */
#define SCIP_DECL_NLPISOLVE(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem)

/** gives solution status
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 * 
 * return: Solution Status
 */
#define SCIP_DECL_NLPIGETSOLSTAT(x) SCIP_NLPSOLSTAT x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem)

/** gives termination reason
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 * 
 * return: Termination Status
 */
#define SCIP_DECL_NLPIGETTERMSTAT(x) SCIP_NLPTERMSTAT x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem)

/** gives primal and dual solution values
 * 
 * solver can return NULL in dual values if not available
 * but if solver provides dual values for one side of variable bounds, then it must also provide those for the other side
 *
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 *
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues buffer to store pointer to array to primal values, or NULL if not needed
 *  - consdualvalues buffer to store pointer to array to dual values of constraints, or NULL if not needed
 *  - varlbdualvalues buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed
 *  - varubdualvalues buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed
 *  - objval pointer to store the objective value, or NULL if not needed
 */
#define SCIP_DECL_NLPIGETSOLUTION(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, SCIP_Real** primalvalues, \
      SCIP_Real** consdualvalues, SCIP_Real** varlbdualvalues, SCIP_Real** varubdualvalues, SCIP_Real* objval)

/** gives solve statistics
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - statistics pointer to store statistics
 * 
 * output:
 *  - statistics solve statistics
 */
#define SCIP_DECL_NLPIGETSTATISTICS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, SCIP_NLPSTATISTICS* statistics)

/** gives required size of a buffer to store a warmstart object
 * 
 *  input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - size pointer to store required size for warmstart buffer
 * 
 * output:
 *  - size required size for warmstart buffer
 */
#define SCIP_DECL_NLPIGETWARMSTARTSIZE(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, size_t* size)

/** stores warmstart information in buffer
 * 
 * required size of buffer should have been obtained by SCIPnlpiGetWarmstartSize before
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - buffer memory to store warmstart information
 * 
 * output:
 *  - buffer warmstart information in solver specific data structure
 */
#define SCIP_DECL_NLPIGETWARMSTARTMEMO(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, void* buffer)

/** sets warmstart information in solver
 * 
 * write warmstart to buffer
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - buffer warmstart information
 */
#define SCIP_DECL_NLPISETWARMSTARTMEMO(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, void* buffer)

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival pointer to store the parameter value
 * 
 * output:
 *  - ival parameter value
 */
#define SCIP_DECL_NLPIGETINTPAR(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, SCIP_NLPPARAM type, int* ival)

/** sets integer parameter of NLP
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival parameter value
 */
#define SCIP_DECL_NLPISETINTPAR(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, SCIP_NLPPARAM type, int ival)

/** gets floating point parameter of NLP
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval pointer to store the parameter value
 * 
 * output:
 *  - dval parameter value
 */
#define SCIP_DECL_NLPIGETREALPAR(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, SCIP_NLPPARAM type, SCIP_Real* dval)

/** sets floating point parameter of NLP
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval parameter value
 */
#define SCIP_DECL_NLPISETREALPAR(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, SCIP_NLPPARAM type, SCIP_Real dval)

/** gets string parameter of NLP
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - sval pointer to store the string value, the user must not modify the string
 * 
 * output:
 *  - sval parameter value
 */
#define SCIP_DECL_NLPIGETSTRINGPAR(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, SCIP_NLPPARAM type, const char** sval)

/** sets string parameter of NLP
 * 
 * input:
 *  - scip SCIP data structure
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - sval parameter value
 */
#define SCIP_DECL_NLPISETSTRINGPAR(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, SCIP_NLPPARAM type, const char* sval)

/**@} */
#ifdef __cplusplus
}
#endif

#endif /*__SCIP_TYPE_NLPI_H__ */
