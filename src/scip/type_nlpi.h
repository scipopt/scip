/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_nlpi.h,v 1.2 2009/07/28 13:10:15 bzfviger Exp $"

/**@file   type_nlpi.h
 * @brief  type definitions for specific NLP solvers interface
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 *
 * @note This is a non-final version of the NLPI specification. It is likely to change soon. Use with care!
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_NLPI_H__
#define __SCIP_TYPE_NLPI_H__

#ifdef WITH_NL
#include "type_expression.h"
#else
typedef struct SCIP_ExprTree      SCIP_EXPRTREE;      /**< expression tree for general nonlinear functions (maybe in SCIP 1.3) */
#endif

typedef struct SCIP_Nlpi          SCIP_NLPI;          /**< NLP solver interface */
typedef struct SCIP_NlpiData      SCIP_NLPIDATA;      /**< locally defined NLP solver interface data */
typedef struct SCIP_NlpStatistics SCIP_NLPSTATISTICS; /**< NLP solve statistics */

/** NLP solver parameter */
enum
SCIP_NlpParam
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
typedef enum SCIP_NlpParam SCIP_NLPIPARAM;

/** NLP solution status */
enum
SCIP_NlpSolStat
{
   SCIP_NLPSOLSTAT_GLOBOPT        = 0,    /**< solved to global optimality */
   SCIP_NLPSOLSTAT_LOCOPT         = 1,    /**< solved to local  optimality */
   SCIP_NLPSOLSTAT_FEASIBLE       = 2,    /**< feasible solution found */
   SCIP_NLPSOLSTAT_LOCINFEASIBLE  = 3,    /**< solution found is local infeasible */
   SCIP_NLPSOLSTAT_GLOBINFEASIBLE = 4,    /**< problem is proven infeasible */
   SCIP_NLPSOLSTAT_UNBOUNDED      = 5,    /**< problem is unbounded */
   SCIP_NLPSOLSTAT_UNKNOWN        = 6     /**< unknown solution status (e.g., problem not solved yet) */
};
typedef enum SCIP_NlpSolStat SCIP_NLPSOLSTAT;

/** NLP termination status */
enum
SCIP_NlpiTermStat
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
typedef enum SCIP_NlpiTermStat SCIP_NLPITERMSTAT;

/** initializes an NLP interface structure
 * Input:
 *  - nlpi datastructure for solver interface
 *  - name problem name
 */
#define SCIP_DECL_NLPIINIT( x ) \
SCIP_RETCODE \
x ( \
   SCIP*        scip, \
   SCIP_NLPI*   nlpi, \
   const char*  name \
   )

/** destructor of NLP interface to free user data
 * Input:
 *  - scip SCIP datastructure
 *  - nlpi datastructure for solver interface
 */
#define SCIP_DECL_NLPIFREE( x ) \
SCIP_RETCODE \
x ( \
   SCIP*           scip, \
   SCIP_NLPIDATA*  data \
   )

/** add variables
 * Input:
 *  - nlpi datastructure for solver interface
 *  - nvars number of variables 
 *  - lb lower bounds of variables
 *  - ub upper bounds of variables
 *  - type types of variables, saying NULL means all are continuous
 *  - varnames names of variables, can be NULL
 */
#define SCIP_DECL_NLPIADDVARS( x )  \
SCIP_RETCODE \
x ( \
   SCIP*            scip, \
   SCIP_NLPI*       nlpi, \
   int              nvars, \
   const SCIP_Real* lb, \
   const SCIP_Real* ub, \
   SCIP_VARTYPE*    type, \
   const char**     varnames \
   )

/* vorsicht namenskonflikt */
/** add constraints 
 * linear coefficients: row(=constraint) oriented matrix
 * quadratic coefficiens: row oriented matrix for each constraint
 * Input:
 *  - nlpi datastructure for solver interface
 *  - ncons number of added constraints
 *  - linoffset start index of each constraints linear coefficients in linind and linval
 *    length: ncons + 1, linoffset[ncons] gives length of linind and linval
 *    may be NULL in case of no linear part
 *  - linind variable indices
 *    may be NULL in case of no linear part
 *  - linval coefficient values
 *    may be NULL in case of no linear part
 *  - nquadrows number of columns in matrix of quadratic part for each constraint
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadrowidx indices of variables for which a quadratic part is specified
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadoffset start index of each rows quadratic coefficients in quadind[.] and quadval[.]
 *    indices are given w.r.t. quadrowidx., i.e., quadoffset[.][i] gives the start index of row quadrowidx[.][i] in quadval[.]
 *    quadoffset[.][nquadrows[.]] gives length of quadind[.] and quadval[.]
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadind column indices w.r.t. quadrowidx, i.e., quadrowidx[quadind[.][i]] gives the index of the variable corresponding to entry i
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadval coefficient values
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - exprtree expression tree for nonquadratic part of constraints
 *    entry of array may be NULL in case of no nonquadratic part
 *    may be NULL in case of no nonquadratic part in any constraint
 *  - names of constraints, may be NULL or entries may be NULL
 *  - cls solvers classification of constraint
 */
#define SCIP_DECL_NLPIADDCONSTRAINTS( x )  \
SCIP_RETCODE \
x ( \
   SCIP*                       scip, \
   SCIP_NLPI*                  nlpi, \
   int                         ncons, \
   const SCIP_Real*            lhs, \
   const SCIP_Real*            rhs, \
   const int*                  linoffset, \
   const int*                  linind, \
   const SCIP_Real*            linval, \
   const int*                  nquadrows, \
   int* const*                 quadrowidx, \
   int* const*                 quadoffset, \
   int* const*                 quadind, \
   SCIP_Real* const*           quadval, \
   SCIP_EXPRTREE* const*       exprtree, \
   const char**                names \
   )

/** sets or overwrites objective, a minization problem is expected
 *  May change sparsity pattern.
 * Input:
 *  - nlpi datastructure for solver interface
 *  - linind variable indices
 *    may be NULL in case of no linear part
 *  - linval coefficient values
 *    may be NULL in case of no linear part
 *  - nquadcols number of columns in matrix of quadratic part
 *  - quadcols indices of variables for which a quadratic part is specified
 *    may be NULL in case of no quadratic part
 *  - quadoffset start index of each rows quadratic coefficients in quadind and quadval
 *    quadoffset[.][nquadcols] gives length of quadind and quadval
 *    may be NULL in case of no quadratic part
 *  - quadind column indices
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadval coefficient values
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - exprtree expression tree for nonquadratic part of constraints
 *    entry of array may be NULL in case of no nonquadratic part
 *    may be NULL in case of no nonquadratic part in any constraint
 *  - objective values offset
 */
#define SCIP_DECL_NLPISETOBJECTIVE( x ) \
SCIP_RETCODE \
x ( \
   SCIP*                      scip, \
   SCIP_NLPI*                 nlpi, \
   int                        nlin, \
   const int*                 linind, \
   const SCIP_Real*           linval, \
   int                        nquadcols, \
   const int*                 quadcols, \
   const int*                 quadoffset, \
   const int*                 quadind, \
   const SCIP_Real*           quadval, \
   const SCIP_EXPRTREE*       exprtree, \
   const SCIP_Real            constant \
   )

/** change variable bounds
 * Input:
 *  - nlpi datastructure for solver interface
 *  - nvars number of variables to change bounds
 *  - indices of variables to change bounds
 *  - new lower bounds
 *  - new upper bounds
 */
#define SCIP_DECL_NLPICHGVARBOUNDS( x ) \
SCIP_RETCODE \
x ( \
   SCIP*                 scip, \
   SCIP_NLPI*            nlpi, \
   const int             nvars, \
   const int*            indices, \
   const SCIP_Real*      lb, \
   const SCIP_Real*      ub \
   )


/** change constraint bounds
 * Input:
 *  - nlpi datastructure for solver interface
 *  - ncons number of constraints to change bounds
 *  - indices of constraints to change bounds
 *  - new lower bounds
 *  - new upper bounds
 */
#define SCIP_DECL_NLPICHGCONSBOUNDS( x ) \
SCIP_RETCODE \
x ( \
   SCIP*                 scip, \
   SCIP_NLPI*            nlpi, \
   const int             ncons, \
   const int*            indices, \
   const SCIP_Real*      lb, \
   const SCIP_Real*      ub \
   )

/** delete a set of variables
 * Input:
 *  - nlpi datastructure for solver interface
 *  - dstat input: deletion status of vars; 1 if var should be deleted, 0 if not
 *          output: new position of var, -1 if var was deleted
 */
#define SCIP_DECL_NLPIDELVARSET( x ) \
SCIP_RETCODE \
x( \
   SCIP*          scip, \
   SCIP_NLPI*     nlpi, \
   int*           dstat \
   )

/** delete a set of constraints
 * Input:
 *  - nlpi datastructure for solver interface
 *  - dstat input: deletion status of rows; 1 if row should be deleted, 0 if not
 *          output: new position of row, -1 if row was deleted
 */
#define SCIP_DECL_NLPIDELCONSSET( x ) \
SCIP_RETCODE \
x( \
   SCIP*          scip, \
   SCIP_NLPI*     nlpi, \
   int*           dstat \
   )

/** change one linear coefficient in a constraint or objective
 * Input:
 *  - nlpi datastructure for solver interface
 *  - cons   index of constraint or -1 for objective
 *  - nvals number of values in linear constraint
 *  - varidx index of variable
 *  - value  new value for coefficient
 * Return: Error if coefficient did not exist before
 */
#define SCIP_DECL_NLPICHGLINEARCOEFS( x ) \
SCIP_RETCODE \
x( \
   SCIP*                scip, \
   SCIP_NLPI*           nlpi, \
   const int            cons, \
   int                  nvals, \
   const int            *varidx, \
   const SCIP_Real      *value \
   )

# if 0
/** change all coefficients in linear part of a constraint or objective
 * Input:
 *  - nlpi datastructure for solver interface
 *  - cons   index of constraint or -1 for objective
 *  - ncoeff number of coefficients to set (should be same as initially given)
 *  - values coefficients corresponding to same indices as used when constraint/objective was constructed
 */
#define SCIP_DECL_NLPICHGALLLINEARCOEF( x ) \
SCIP_RETCODE \
x ( \
   SCIP*            scip, \
   SCIP_NLPI*       nlpi, \
   const int        cons, \
   const int        ncoeff, \
   const SCIP_Real* values \
   );

#endif

/** change one coefficient in the quadratic part of a constraint or objective
 * Input:
 *  - nlpi datastructure for solver interface
 *  - cons index of constraint or -1 for objective
 *  - row row offset containing modified indices
 *  - col cols containing modified indices to the corresponding row offset
 *  - values coefficients corresponding to same indices as used when constraint/objective was constructed
 * Return: Error if coefficient did not exist before
 */
#define SCIP_DECL_NLPICHGQUADCOEFS( x ) \
SCIP_RETCODE \
x ( \
   SCIP*              scip, \
   SCIP_NLPI*         nlpi, \
   const int          cons, \
   const int          nentries, \
   const int*         row, \
   const int*         col, \
   SCIP_Real*         value \
   )

#if 0
/** change all coefficients in quadratic part of a constraint or objective
 * Input:
 *  - nlpi datastructure for solver interface
 *  - cons   index of constraint or -1 for objective
 *  - ncoeff number of coefficients to set (should be same as initially given)
 *  - values coefficients corresponding to same indices as used when constraint/objective was constructed
 */
#define SCIP_DECL_NLPICHGALLQUADCOEF( x ) \
SCIP_RETCODE \
x ( \
   SCIP*            scip, \
   SCIP_NLPI*       nlpi, \
   const int        cons, \
   const int        ncoeff, \
   const SCIP_Real* values \
   );
#endif

/** change a parameter constant in the nonlinear part
 * to discuss
 * Input:
 *  - nlpi datastructure for solver interface
 *  - cons index of constraint or -1 for objective
 *  - idx index of parameter
 *  - value new value for nonlinear parameter
 * Return: Error if parameter does not exist
 */
#define SCIP_DECL_NLPICHGNONLINCOEF( x ) \
SCIP_RETCODE \
x ( \
   SCIP*               scip, \
   SCIP_NLPI*          nlpi, \
   const int           cons, \
   const int           idx, \
   SCIP_Real           value \
   )

/** sets initial guess for primal variables
 * Input:
 *  - nlpi datastructure for solver interface
 *  - values initial starting solution
 */
#define SCIP_DECL_NLPISETINITIALGUESS( x ) \
SCIP_RETCODE \
x ( \
   SCIP*       scip, \
   SCIP_NLPI*  nlpi, \
   SCIP_Real*  values \
   )

/** tries to solve NLP
 *  - nlpi datastructure for solver interface
 */
#define SCIP_DECL_NLPISOLVE( x ) \
SCIP_RETCODE \
x ( \
   SCIP*       scip, \
   SCIP_NLPI*  nlpi \
   )

/** gives solution status
 *  - nlpi datastructure for solver interface
 *  Retrun: Solution Status
 */
#define SCIP_DECL_NLPIGETSOLSTAT( x ) \
SCIP_NLPSOLSTAT \
x ( \
   SCIP*         scip, \
   SCIP_NLPI*    nlpi \
   )

/** gives termination reason
 *  Input:
 *   - nlpi datastructure for solver interface
 *  Retrun: Termination Status
 */
#define SCIP_DECL_NLPIGETTERMSTAT( x ) \
SCIP_NLPITERMSTAT \
x ( \
   SCIP*        scip, \
   SCIP_NLPI*   nlpi \
   )

/** gives primal solution
 *  Input:
 *   - nlpi datastructure for solver interface
 *   - buffer to store pointer to primal values
 */
#define SCIP_DECL_NLPIGETSOLUTION( x ) \
SCIP_RETCODE \
x ( \
   SCIP*         scip, \
   SCIP_NLPI*    nlpi, \
   SCIP_Real**   primalvalues \
   )

/** gives solve statistics
 *  Input:
 *   - nlpi datastructure for solver interface
 *   - buffer to store statistics
 */
#define SCIP_DECL_NLPIGETSTATISTICS( x ) \
SCIP_RETCODE \
x ( \
   SCIP*               scip, \
   SCIP_NLPI*          nlpi, \
   SCIP_NLPSTATISTICS* statistics \
   )

/** gives required size of a buffer to store a warmstart object
 *  Input:
 *   - nlpi datastructure for solver interface
 *   - size pointer to store required size for warmstart buffer
 */
#define SCIP_DECL_NLPIGETWARMSTARTSIZE( x ) \
SCIP_RETCODE \
x ( \
   SCIP*       scip, \
   SCIP_NLPI*  nlpi, \
   size_t*     size \
   )

/** stores warmstart information in buffer
 * required size of buffer should have been obtained by SCIPnlpiGetWarmstartSize before
 *  Input:
 *   - nlpi datastructure for solver interface
 *   - buffer to store warmstart information
 */
#define SCIP_DECL_NLPIGETWARMSTARTMEMO( x ) \
SCIP_RETCODE \
x ( \
   SCIP*       scip, \
   SCIP_NLPI*  nlpi, \
   void*       buffer \
   )

/** sets warmstart information in solver
 * write warmstart to buffer
 *  Input:
 *   - nlpi datastructure for solver interface
 *   - buffer to store warmstart information

 */
#define SCIP_DECL_NLPISETWARMSTARTMEMO( x ) \
SCIP_RETCODE \
x ( \
   SCIP*       scip, \
   SCIP_NLPI*  nlpi, \
   void*       buffer \
   )

/** gets pointer for NLP solver
 *  to do dirty stuff
 *  Input:
 *   - nlpi datastructure for solver interface
 *  Return: void pointer to solve interface
 */
#define SCIP_DECL_NLPIGETSOLVERPOINTER( x ) \
void* \
x ( \
   SCIP*       scip, \
   SCIP_NLPI*  nlpi \
   )

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP
 * input:
 * - nlpi  NLP interface structure
 * - type  parameter number
 * - ival  buffer to store the parameter value
*/
#define SCIP_DECL_NLPIGETINTPAR( x ) \
SCIP_RETCODE \
x ( \
   SCIP*                 scip, \
   SCIP_NLPI*            nlpi, \
   SCIP_NLPIPARAM        type, \
   int*                  ival \
   )

/** sets integer parameter of NLP
 * input:
 * - nlpi  NLP interface structure
 * - type  parameter number
 * - ival  parameter value
 */
#define SCIP_DECL_NLPISETINTPAR( x ) \
SCIP_RETCODE \
x ( \
   SCIP*            scip, \
   SCIP_NLPI*       nlpi, \
   SCIP_NLPIPARAM   type, \
   int              ival \
   )


/** gets floating point parameter of NLP
 * input:
 * - nlpi  NLP interface structure
 * - type  parameter number
 * - dval  buffer to store the parameter value
 */
#define SCIP_DECL_NLPIGETREALPAR( x ) \
SCIP_RETCODE \
x ( \
   SCIP*                 scip, \
   SCIP_NLPI*            nlpi, \
   SCIP_NLPIPARAM        type, \
   SCIP_Real*            dval \
   )

/** sets floating point parameter of NLP
 * input:
 * - nlpi  NLP interface structure
 * - type  parameter number
 * - dval  parameter value
 */
#define SCIP_DECL_NLPISETREALPAR( x ) \
SCIP_RETCODE \
x ( \
   SCIP*                  scip, \
   SCIP_NLPI*             nlpi, \
   SCIP_NLPIPARAM         type, \
   SCIP_Real              dval \
   )

/**@} */

#endif /*__SCIP_TYPE_NLPI_H__ */
