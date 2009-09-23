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
#pragma ident "@(#) $Id: type_nlpi.h,v 1.11 2009/09/23 08:54:57 bzfheinz Exp $"

/**@file   type_nlpi.h
 * @ingroup TYPEDEFINITONS
 * @brief  type definitions for specific NLP solvers interface
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 *
 * @note This is a non-final version of the NLPI specification. It is likely to change soon. Use with care!
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_NLPI_H__
#define __SCIP_TYPE_NLPI_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_ExprTree      SCIP_EXPRTREE;      /**< expression tree for general nonlinear functions (maybe in SCIP 1.3) */
typedef struct SCIP_Nlpi          SCIP_NLPI;          /**< NLP solver interface */
typedef struct SCIP_NlpiData      SCIP_NLPIDATA;      /**< locally defined NLP solver interface data */
typedef struct SCIP_NlpStatistics SCIP_NLPSTATISTICS; /**< NLP solve statistics */

/** NLP solver parameter */
enum SCIP_NlpParam
{
   SCIP_NLPPAR_FROMSCRATCH    =  0,      /**< solver should start from scratch at next call?: 0 no, 1 yes (int) */
   SCIP_NLPPAR_VERBLEVEL      =  1,      /**< verbosity level of output of NLP solver to the screen: 0 off, 1 normal, 2 debug (int) */
   SCIP_NLPPAR_FEASTOL        =  2,      /**< feasibility tolerance for primal variables and slacks (real) */
   SCIP_NLPPAR_RELOBJTOL      =  3,      /**< relative objective tolerance (real) */
   SCIP_NLPPAR_LOBJLIM        =  4,      /**< lower objective limit (cutoff) (real) */
   SCIP_NLPPAR_INFINITY       =  5,      /**< value for infinity used to decide unbounded sides, unbounded variable and constraint bounds, and upper objective limit (real) */
   SCIP_NLPPAR_ITLIM          =  6,      /**< NLP iteration limit (int) */
   SCIP_NLPPAR_TILIM          =  7       /**< NLP time limit (real) */
};
typedef enum SCIP_NlpParam SCIP_NLPIPARAM;  /**< NLP solver parameter */

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
enum SCIP_NlpiTermStat
{
   SCIP_NLPITERMSTAT_OKAY          = 0,    /**< terminated successfully */
   SCIP_NLPITERMSTAT_TILIM         = 1,    /**< time limit exceeded */
   SCIP_NLPITERMSTAT_ITLIM         = 2,    /**< iteration limit exceeded */
   SCIP_NLPITERMSTAT_LOBJLIM       = 3,    /**< lower objective limit reached */
   SCIP_NLPITERMSTAT_UOBJLIM       = 4,    /**< upper objective limit (= infinity) reached */
   SCIP_NLPITERMSTAT_NUMERR        = 5,    /**< stopped on numerical error */
   SCIP_NLPITERMSTAT_EVALERR       = 6,    /**< stopped on function evaluation error */
   SCIP_NLPITERMSTAT_MEMERR        = 7,    /**< memory exceeded */
   SCIP_NLPITERMSTAT_LICERR        = 8,    /**< licence error */
   SCIP_NLPITERMSTAT_OTHER         = 9     /**< other error (= this should never happen) */
};
typedef enum SCIP_NlpiTermStat SCIP_NLPITERMSTAT;  /** NLP solver termination status */

/** initializes an NLP interface structure
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - name problem name
 */
#define SCIP_DECL_NLPIINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, const char* name)

/** destructor of NLP interface to free nlpi data
 * 
 * input:
 *  - scip SCIP datastructure
 *  - data nlpi data of solver interface
 */
#define SCIP_DECL_NLPIFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPIDATA* data)

/** add variables
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - nvars number of variables 
 *  - lbs lower bounds of variables
 *  - ubs upper bounds of variables
 *  - types types of variables, saying NULL means all are continuous
 *  - varnames names of variables, can be NULL
 */
#define SCIP_DECL_NLPIADDVARS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, int nvars, const SCIP_Real* lbs, \
   const SCIP_Real* ubs, SCIP_VARTYPE* types, const char** varnames)

/** add constraints 
 * linear coefficients: row(=constraint) oriented matrix
 * quadratic coefficiens: row oriented matrix for each constraint
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
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
#define SCIP_DECL_NLPIADDCONSTRAINTS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, int ncons, const SCIP_Real* lhss, \
   const SCIP_Real* rhss, const int* linoffsets, const int* lininds, const SCIP_Real* linvals, const int* nquadrows, \
   int* const* quadrowidxs, int* const* quadoffsets, int* const* quadinds, SCIP_Real* const* quadvals, \
   int* const* exprvaridxs, SCIP_EXPRTREE* const* exprtrees, const char** names)

/** sets or overwrites objective, a minimization problem is expected
 *  May change sparsity pattern.
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
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
#define SCIP_DECL_NLPISETOBJECTIVE(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, int nlins, const int* lininds, \
   const SCIP_Real* linvals, int nquadcols, const int* quadcols, const int* quadoffsets, const int* quadinds, \
   const SCIP_Real* quadvals, const int* exprvaridxs, const SCIP_EXPRTREE* exprtree, const SCIP_Real constant)

/** change variable bounds
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - nvars number of variables to change bounds
 *  - indices indices of variables to change bounds
 *  - lbs new lower bounds
 *  - ubs new upper bounds
 */
#define SCIP_DECL_NLPICHGVARBOUNDS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, const int nvars, const int* indices, \
   const SCIP_Real* lbs, const SCIP_Real* ubs)

/** change constraint bounds
 *
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - ncons number of constraints to change bounds
 *  - indices indices of constraints to change bounds
 *  - lbs new lower bounds
 *  - ubs new upper bounds
 */
#define SCIP_DECL_NLPICHGCONSBOUNDS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, const int ncons, const int* indices, \
   const SCIP_Real* lbs, const SCIP_Real* ubs)

/** delete a set of variables
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - dstats deletion status of vars; 1 if var should be deleted, 0 if not
 * 
 * output:
 *  - dstats new position of var, -1 if var was deleted
 */
#define SCIP_DECL_NLPIDELVARSET(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, int* dstats)

/** delete a set of constraints
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - dstats deletion status of rows; 1 if row should be deleted, 0 if not
 * 
 * output:
 *  - dstats new position of row, -1 if row was deleted
 */
#define SCIP_DECL_NLPIDELCONSSET(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, int* dstats)

/** change one linear coefficient in a constraint or objective
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - idx index of constraint or -1 for objective
 *  - nvals number of values in linear constraint
 *  - varidxs indices of variable
 *  - vals new values for coefficient
 * 
 * return: Error if coefficient did not exist before
 */
#define SCIP_DECL_NLPICHGLINEARCOEFS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, const int idx, int nvals, \
   const int* varidxs, const SCIP_Real* vals)

/** change one coefficient in the quadratic part of a constraint or objective
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - idx index of constraint or -1 for objective
 *  - nentries number of values in quadratic constraint
 *  - rows row offset containing modified indices
 *  - cols cols containing modified indices to the corresponding row offset
 *  - values coefficients corresponding to same indices as used when constraint/objective was constructed
 * 
 * return: Error if coefficient did not exist before
 */
#define SCIP_DECL_NLPICHGQUADCOEFS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, const int idx, const int nentries, \
   const int* rows, const int* cols, SCIP_Real* values)

/** change one coefficient in the nonlinear part
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - idxcons index of constraint or -1 for objective
 *  - idxparam index of parameter
 *  - value new value for nonlinear parameter
 * 
 * return: Error if parameter does not exist
 */
#define SCIP_DECL_NLPICHGNONLINCOEF(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, const int idxcons, const int idxparam, \
   SCIP_Real value)

/** sets initial guess for primal variables
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - values initial starting solution
 */
#define SCIP_DECL_NLPISETINITIALGUESS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_Real* values)

/** tries to solve NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 */
#define SCIP_DECL_NLPISOLVE(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi)

/** gives solution status
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 * 
 * return: Solution Status
 */
#define SCIP_DECL_NLPIGETSOLSTAT(x) SCIP_NLPSOLSTAT x (SCIP* scip, SCIP_NLPI* nlpi)

/** gives termination reason
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 * 
 * return: Termination Status
 */
#define SCIP_DECL_NLPIGETTERMSTAT(x) SCIP_NLPITERMSTAT x (SCIP* scip, SCIP_NLPI* nlpi)

/** gives primal solution
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - primalvalues pointer to store primal values
 * 
 * output:
 *  - primalvalues primal values of solution
 */
#define SCIP_DECL_NLPIGETSOLUTION(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_Real** primalvalues)

/** gives solve statistics
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - statistics pointer to store statistics
 * 
 * output:
 *  - statistics solve statistics
 */
#define SCIP_DECL_NLPIGETSTATISTICS(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPSTATISTICS* statistics)

/** gives required size of a buffer to store a warmstart object
 * 
 *  input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - size pointer to store required size for warmstart buffer
 * 
 * output:
 *  - size required size for warmstart buffer
 */
#define SCIP_DECL_NLPIGETWARMSTARTSIZE(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, size_t* size)

/** stores warmstart information in buffer
 * 
 * required size of buffer should have been obtained by SCIPnlpiGetWarmstartSize before
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - buffer memory to store warmstart information
 * 
 * output:
 *  - buffer warmstart information in solver specific data structure
 */
#define SCIP_DECL_NLPIGETWARMSTARTMEMO(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, void* buffer)

/** sets warmstart information in solver
 * 
 * write warmstart to buffer
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  - buffer warmstart information
 */
#define SCIP_DECL_NLPISETWARMSTARTMEMO(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, void* buffer)

/** gets pointer for NLP solver
 * 
 *  to do dirty stuff
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 *  
 * return: void pointer to solver
 */
#define SCIP_DECL_NLPIGETSOLVERPOINTER(x) void* x (SCIP* scip, SCIP_NLPI* nlpi)

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi NLP interface structure
 *  - type parameter number
 *  - ival pointer to store the parameter value
 * 
 * output:
 *  - ival parameter value
 */
#define SCIP_DECL_NLPIGETINTPAR(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPARAM type, int* ival)

/** sets integer parameter of NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi NLP interface structure
 *  - type parameter number
 *  - ival parameter value
 */
#define SCIP_DECL_NLPISETINTPAR(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPARAM type, int ival)

/** gets floating point parameter of NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi NLP interface structure
 *  - type parameter number
 *  - dval pointer to store the parameter value
 * 
 * output:
 *  - dval parameter value
 */
#define SCIP_DECL_NLPIGETREALPAR(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPARAM type, SCIP_Real* dval)

/** sets floating point parameter of NLP
 * 
 * input:
 *  - scip SCIP datastructure
 *  - nlpi NLP interface structure
 *  - type parameter number
 *  - dval parameter value
 */
#define SCIP_DECL_NLPISETREALPAR(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPARAM type, SCIP_Real dval)

/**@} */
#ifdef __cplusplus
}
#endif

#endif /*__SCIP_TYPE_NLPI_H__ */
