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

/**@file   nlpi.c
 * @ingroup OTHER_CFILES
 * @brief  methods for handling nlp interface
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/nlpi.h"
#include "scip/pub_message.h"
#include "scip/pub_nlpi.h"
#include "scip/struct_nlpi.h"
#include "scip/struct_set.h"

/** compares two NLPIs w.r.t. their priority */
SCIP_DECL_SORTPTRCOMP(SCIPnlpiComp)
{  /*lint --e{715}*/
   return ((SCIP_NLPI*)elem2)->priority - ((SCIP_NLPI*)elem1)->priority;
}

/** creates an NLP solver interface */
SCIP_RETCODE SCIPnlpiCreate(
   SCIP_NLPI**                     nlpi,                        /**< pointer to NLP interface data structure */
   const char*                     name,                        /**< name of NLP interface */
   const char*                     description,                 /**< description of NLP interface */
   int                             priority,                    /**< priority of NLP interface */
   SCIP_DECL_NLPICOPY              ((*nlpicopy)),               /**< copying of NLPI, can be NULL */
   SCIP_DECL_NLPIFREE              ((*nlpifree)),               /**< free NLPI user data */
   SCIP_DECL_NLPIGETSOLVERPOINTER  ((*nlpigetsolverpointer)),   /**< get solver pointer, can be NULL */
   SCIP_DECL_NLPICREATEPROBLEM     ((*nlpicreateproblem)),      /**< create a new problem instance */
   SCIP_DECL_NLPIFREEPROBLEM       ((*nlpifreeproblem)),        /**< free a problem instance */
   SCIP_DECL_NLPIGETPROBLEMPOINTER ((*nlpigetproblempointer)),  /**< get problem pointer, can be NULL */
   SCIP_DECL_NLPIADDVARS           ((*nlpiaddvars)),            /**< add variables */
   SCIP_DECL_NLPIADDCONSTRAINTS    ((*nlpiaddconstraints)),     /**< add constraints */
   SCIP_DECL_NLPISETOBJECTIVE      ((*nlpisetobjective)),       /**< set objective */
   SCIP_DECL_NLPICHGVARBOUNDS      ((*nlpichgvarbounds)),       /**< change variable bounds */
   SCIP_DECL_NLPICHGCONSSIDES      ((*nlpichgconssides)),       /**< change constraint sides */
   SCIP_DECL_NLPIDELVARSET         ((*nlpidelvarset)),          /**< delete a set of constraints */
   SCIP_DECL_NLPIDELCONSSET        ((*nlpidelconsset)),         /**< delete a set of constraints */
   SCIP_DECL_NLPICHGLINEARCOEFS    ((*nlpichglinearcoefs)),     /**< change coefficients in linear part of a constraint or objective */
   SCIP_DECL_NLPICHGEXPR           ((*nlpichgexpr)),            /**< change nonlinear expression a constraint or objective */
   SCIP_DECL_NLPICHGOBJCONSTANT    ((*nlpichgobjconstant)),     /**< change the constant offset in the objective */
   SCIP_DECL_NLPISETINITIALGUESS   ((*nlpisetinitialguess)),    /**< set initial guess, can be NULL */
   SCIP_DECL_NLPISOLVE             ((*nlpisolve)),              /**< solve NLP */
   SCIP_DECL_NLPIGETSOLSTAT        ((*nlpigetsolstat)),         /**< get solution status */
   SCIP_DECL_NLPIGETTERMSTAT       ((*nlpigettermstat)),        /**< get termination status */
   SCIP_DECL_NLPIGETSOLUTION       ((*nlpigetsolution)),        /**< get solution */
   SCIP_DECL_NLPIGETSTATISTICS     ((*nlpigetstatistics)),      /**< get solve statistics */
   SCIP_DECL_NLPIGETINTPAR         ((*nlpigetintpar)),          /**< get value of integer parameter */
   SCIP_DECL_NLPISETINTPAR         ((*nlpisetintpar)),          /**< set value of integer parameter */
   SCIP_DECL_NLPIGETREALPAR        ((*nlpigetrealpar)),         /**< get value of floating point parameter */
   SCIP_DECL_NLPISETREALPAR        ((*nlpisetrealpar)),         /**< set value of floating point parameter */
   SCIP_DECL_NLPIGETSTRINGPAR      ((*nlpigetstringpar)),       /**< get value of string parameter */
   SCIP_DECL_NLPISETSTRINGPAR      ((*nlpisetstringpar)),       /**< set value of string parameter */
   SCIP_NLPIDATA*                  nlpidata                     /**< NLP interface local data */
   )
{  /*lint --e{715}*/
   assert(nlpi != NULL);

   assert(name != NULL);
   assert(description != NULL);
   assert(nlpifree != NULL);
   assert(nlpicreateproblem != NULL);
   assert(nlpifreeproblem != NULL);
   assert(nlpiaddvars != NULL);
   assert(nlpiaddconstraints != NULL);
   assert(nlpisetobjective != NULL);
   assert(nlpichgvarbounds != NULL);
   assert(nlpichgconssides != NULL);
   assert(nlpidelconsset != NULL);
   assert(nlpichglinearcoefs != NULL);
   assert(nlpichgobjconstant != NULL);
   assert(nlpisolve != NULL);
   assert(nlpigetsolstat != NULL);
   assert(nlpigettermstat != NULL);
   assert(nlpigetsolution != NULL);
   assert(nlpigetstatistics != NULL);
   assert(nlpigetintpar != NULL);
   assert(nlpisetintpar != NULL);
   assert(nlpigetrealpar != NULL);
   assert(nlpisetrealpar != NULL);
   assert(nlpigetstringpar != NULL);
   assert(nlpisetstringpar != NULL);

   SCIP_ALLOC( BMSallocMemory(nlpi) );

   SCIP_ALLOC( BMSduplicateMemoryArray(&(*nlpi)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*nlpi)->description, description, strlen(description)+1) );
   (*nlpi)->priority = priority;
   (*nlpi)->nlpicopy = nlpicopy;
   (*nlpi)->nlpifree = nlpifree;
   (*nlpi)->nlpigetsolverpointer = nlpigetsolverpointer;
   (*nlpi)->nlpicreateproblem = nlpicreateproblem;
   (*nlpi)->nlpifreeproblem = nlpifreeproblem;
   (*nlpi)->nlpigetproblempointer = nlpigetproblempointer;
   (*nlpi)->nlpiaddvars = nlpiaddvars;
   (*nlpi)->nlpiaddconstraints = nlpiaddconstraints;
   (*nlpi)->nlpisetobjective = nlpisetobjective;
   (*nlpi)->nlpichgvarbounds = nlpichgvarbounds;
   (*nlpi)->nlpichgconssides = nlpichgconssides;
   (*nlpi)->nlpidelvarset = nlpidelvarset;
   (*nlpi)->nlpidelconsset = nlpidelconsset;
   (*nlpi)->nlpichglinearcoefs = nlpichglinearcoefs;
   (*nlpi)->nlpichgobjconstant = nlpichgobjconstant;
   (*nlpi)->nlpisetinitialguess = nlpisetinitialguess;
   (*nlpi)->nlpisolve = nlpisolve;
   (*nlpi)->nlpigetsolstat = nlpigetsolstat;
   (*nlpi)->nlpigettermstat = nlpigettermstat;
   (*nlpi)->nlpigetsolution = nlpigetsolution;
   (*nlpi)->nlpigetstatistics = nlpigetstatistics;
   (*nlpi)->nlpigetintpar = nlpigetintpar;
   (*nlpi)->nlpisetintpar = nlpisetintpar;
   (*nlpi)->nlpigetrealpar = nlpigetrealpar;
   (*nlpi)->nlpisetrealpar = nlpisetrealpar;
   (*nlpi)->nlpigetstringpar = nlpigetstringpar;
   (*nlpi)->nlpisetstringpar = nlpisetstringpar;
   (*nlpi)->nlpidata = nlpidata;

   return SCIP_OKAY;
}

/** copies an NLPI and includes it into another SCIP instance */
SCIP_RETCODE SCIPnlpiCopyInclude(
   SCIP_NLPI*            sourcenlpi,         /**< the NLP interface to copy */
   SCIP_SET*             targetset           /**< global SCIP settings where to include copy */
   )
{
   assert(sourcenlpi != NULL);
   assert(targetset != NULL);

   if( sourcenlpi->nlpicopy != NULL )
   {
      SCIP_CALL( sourcenlpi->nlpicopy(targetset->scip, sourcenlpi) );
   }

   return SCIP_OKAY;
}

/** frees NLPI */
SCIP_RETCODE SCIPnlpiFree(
   SCIP_NLPI**           nlpi,               /**< pointer to NLPI data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nlpi  != NULL);
   assert(*nlpi != NULL);
   assert(set   != NULL);

   if( (*nlpi)->nlpifree != NULL )
   {
      SCIP_CALL( (*nlpi)->nlpifree(set->scip, *nlpi, &(*nlpi)->nlpidata) );
      assert((*nlpi)->nlpidata == NULL);
   }
   BMSfreeMemoryArray(&(*nlpi)->name);
   BMSfreeMemoryArray(&(*nlpi)->description);
   BMSfreeMemory(nlpi);

   assert(*nlpi == NULL);

   return SCIP_OKAY;
}

/** gets pointer for NLP solver */
void* SCIPnlpiGetSolverPointer(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi                /**< solver interface */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);

   if( nlpi->nlpigetsolverpointer != NULL )
      return nlpi->nlpigetsolverpointer(set->scip, nlpi);
   else
      return NULL;
}

/** creates a problem instance */
SCIP_RETCODE SCIPnlpiCreateProblem(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM**    problem,            /**< problem pointer to store the problem data */
   const char*           name                /**< name of problem, can be NULL */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpicreateproblem != NULL);
   assert(problem != NULL);

   return nlpi->nlpicreateproblem(set->scip, nlpi, problem, name);
}

/** frees a problem instance */
SCIP_RETCODE SCIPnlpiFreeProblem(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM**    problem             /**< pointer where problem instance is stored */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpifreeproblem != NULL);
   assert(problem != NULL);

   return nlpi->nlpifreeproblem(set->scip, nlpi, problem);
}

/** gets pointer to solver-internal problem instance */
void* SCIPnlpiGetProblemPointer(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem             /**< problem instance */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(problem != NULL);

   if( nlpi->nlpigetproblempointer != NULL )
      return nlpi->nlpigetproblempointer(set->scip, nlpi, problem);
   else
      return NULL;
}

/** add variables to nlpi */
SCIP_RETCODE SCIPnlpiAddVars(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      lbs,                /**< lower bounds of variables, can be NULL if -infinity */
   const SCIP_Real*      ubs,                /**< upper bounds of variables, can be NULL if +infinity */
   const char**          varnames            /**< names of variables, can be NULL */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpiaddvars != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpiaddvars(set->scip, nlpi, problem, nvars, lbs, ubs, varnames) );

   return SCIP_OKAY;
}

/** add constraints to nlpi */
SCIP_RETCODE SCIPnlpiAddConstraints(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   nconss,             /**< number of constraints */
   const SCIP_Real*      lhss,               /**< left hand sides of constraints, can be NULL if -infinity */
   const SCIP_Real*      rhss,               /**< right hand sides of constraints, can be NULL if +infinity */
   const int*            nlininds,           /**< number of linear coefficients for each constraint, may be NULL in case of no linear part */
   int* const*           lininds,            /**< indices of variables for linear coefficients for each constraint, may be NULL in case of no linear part */
   SCIP_Real* const*     linvals,            /**< values of linear coefficient for each constraint, may be NULL in case of no linear part */
   SCIP_EXPR**           exprs,              /**< expressions for nonlinear part of constraints, entry of array may be NULL in case of no nonlinear part, may be NULL in case of no nonlinear part in any constraint */
   const char**          names               /**< names of constraints, may be NULL or entries may be NULL */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpiaddconstraints != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpiaddconstraints(set->scip, nlpi, problem, nconss, lhss, rhss, nlininds, lininds, linvals, exprs, names) );

   return SCIP_OKAY;
}

/** sets or overwrites objective, a minimization problem is expected */
SCIP_RETCODE SCIPnlpiSetObjective(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   nlins,              /**< number of linear variables */
   const int*            lininds,            /**< variable indices, may be NULL in case of no linear part */
   const SCIP_Real*      linvals,            /**< coefficient values, may be NULL in case of no linear part */
   SCIP_EXPR*            expr,               /**< expression for nonlinear part of objective function, may be NULL in case of no nonlinear part */
   const SCIP_Real       constant            /**< objective value offset */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpisetobjective != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpisetobjective(set->scip, nlpi, problem, nlins, lininds, linvals, expr, constant) );

   return SCIP_OKAY;
}

/** change variable bounds */
SCIP_RETCODE SCIPnlpiChgVarBounds(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   const int             nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds */
   const SCIP_Real*      ubs                 /**< new upper bounds */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpichgvarbounds != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpichgvarbounds(set->scip, nlpi, problem, nvars, indices, lbs, ubs) );

   return SCIP_OKAY;
}

/** change constraint sides */
SCIP_RETCODE SCIPnlpiChgConsSides(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   nconss,             /**< number of constraints to change sides */
   const int*            indices,            /**< indices of constraints to change sides */
   const SCIP_Real*      lhss,               /**< new left hand sides */
   const SCIP_Real*      rhss                /**< new right hand sides */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpichgconssides != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpichgconssides(set->scip, nlpi, problem, nconss, indices, lhss, rhss) );

   return SCIP_OKAY;
}

/** delete a set of variables */
SCIP_RETCODE SCIPnlpiDelVarSet(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int*                  dstats,             /**< deletion status of vars; 1 if var should be deleted, 0 if not */
   int                   dstatssize          /**< size of the dstats array */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpidelvarset != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpidelvarset(set->scip, nlpi, problem, dstats, dstatssize) );

   return SCIP_OKAY;
}

/** delete a set of constraints */
SCIP_RETCODE SCIPnlpiDelConsSet(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int*                  dstats,             /**< deletion status of constraints; 1 if constraint should be deleted, 0 if not */
   int                   dstatssize          /**< size of the dstats array */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpidelconsset != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpidelconsset(set->scip, nlpi, problem, dstats, dstatssize) );

   return SCIP_OKAY;
}

/** changes or adds linear coefficients in a constraint or objective */
SCIP_RETCODE SCIPnlpiChgLinearCoefs(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   idx,                /**< index of constraint or -1 for objective */
   int                   nvals,              /**< number of values in linear constraint to change */
   const int*            varidxs,            /**< indices of variables which coefficient to change */
   const SCIP_Real*      vals                /**< new values for coefficients */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpichglinearcoefs != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpichglinearcoefs(set->scip, nlpi, problem, idx, nvals, varidxs, vals) );

   return SCIP_OKAY;
}

/** change the expression in the nonlinear part */
SCIP_RETCODE SCIPnlpiChgExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   idxcons,            /**< index of constraint or -1 for objective */
   SCIP_EXPR*            expr                /**< new expression for constraint or objective, or NULL to only remove previous tree */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpichgexpr != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpichgexpr(set->scip, nlpi, problem, idxcons, expr) );

   return SCIP_OKAY;
}

/** change the constant offset in the objective */
SCIP_RETCODE SCIPnlpiChgObjConstant(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_Real             objconstant         /**< new value for objective constant */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpichgobjconstant != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpichgobjconstant(set->scip, nlpi, problem, objconstant) );

   return SCIP_OKAY;
}

/** sets initial guess for primal variables */
SCIP_RETCODE SCIPnlpiSetInitialGuess(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_Real*            primalvalues,       /**< initial primal values for variables, or NULL to clear previous values */
   SCIP_Real*            consdualvalues,     /**< initial dual values for constraints, or NULL to clear previous values */
   SCIP_Real*            varlbdualvalues,    /**< initial dual values for variable lower bounds, or NULL to clear previous values */
   SCIP_Real*            varubdualvalues     /**< initial dual values for variable upper bounds, or NULL to clear previous values */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(problem != NULL);

   if( nlpi->nlpisetinitialguess != NULL )
   {
      SCIP_CALL( nlpi->nlpisetinitialguess(set->scip, nlpi, problem, primalvalues, consdualvalues, varlbdualvalues, varubdualvalues) );
   }

   return SCIP_OKAY;
}

/** tries to solve NLP */
SCIP_RETCODE SCIPnlpiSolve(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem             /**< problem instance */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpisolve != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpisolve(set->scip, nlpi, problem) );

   return SCIP_OKAY;
}

/** gives solution status */
SCIP_NLPSOLSTAT SCIPnlpiGetSolstat(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem             /**< problem instance */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpigetsolstat != NULL);
   assert(problem != NULL);

   return nlpi->nlpigetsolstat(set->scip, nlpi, problem);
}

/** gives termination reason */
SCIP_NLPTERMSTAT SCIPnlpiGetTermstat(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem             /**< problem instance */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpigettermstat != NULL);
   assert(problem != NULL);

   return nlpi->nlpigettermstat(set->scip, nlpi, problem);
}

/** gives primal and dual solution
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 */
SCIP_RETCODE SCIPnlpiGetSolution(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_Real**           primalvalues,       /**< buffer to store pointer to array to primal values, or NULL if not needed */
   SCIP_Real**           consdualvalues,     /**< buffer to store pointer to array to dual values of constraints, or NULL if not needed */
   SCIP_Real**           varlbdualvalues,    /**< buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed */
   SCIP_Real**           varubdualvalues,    /**< buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed */
   SCIP_Real*            objval              /**< pointer to store the objective value, or NULL if not needed */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpigetsolution != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpigetsolution(set->scip, nlpi, problem, primalvalues, consdualvalues, varlbdualvalues, varubdualvalues, objval) );

   return SCIP_OKAY;
}

/** gives solve statistics */
SCIP_RETCODE SCIPnlpiGetStatistics(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpigetstatistics != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpigetstatistics(set->scip, nlpi, problem, statistics) );

   return SCIP_OKAY;
}

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP */
SCIP_RETCODE SCIPnlpiGetIntPar(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpigetintpar != NULL);
   assert(problem != NULL);
   assert(ival != NULL);

   SCIP_CALL( nlpi->nlpigetintpar(set->scip, nlpi, problem, type, ival) );

   return SCIP_OKAY;
}

/** sets integer parameter of NLP */
SCIP_RETCODE SCIPnlpiSetIntPar(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpisetintpar != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpisetintpar(set->scip, nlpi, problem, type, ival) );

   return SCIP_OKAY;
}

/** gets floating point parameter of NLP */
SCIP_RETCODE SCIPnlpiGetRealPar(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpigetrealpar != NULL);
   assert(problem != NULL);
   assert(dval != NULL);

   SCIP_CALL( nlpi->nlpigetrealpar(set->scip, nlpi, problem, type, dval) );

   return SCIP_OKAY;
}

/** sets floating point parameter of NLP */
SCIP_RETCODE SCIPnlpiSetRealPar(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpisetrealpar != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpisetrealpar(set->scip, nlpi, problem, type, dval) );

   return SCIP_OKAY;
}

/** gets string parameter of NLP */
SCIP_RETCODE SCIPnlpiGetStringPar(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char**          sval                /**< pointer to store the string value, the user must not modify the string */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpigetstringpar != NULL);
   assert(problem != NULL);
   assert(sval != NULL);

   SCIP_CALL( nlpi->nlpigetstringpar(set->scip, nlpi, problem, type, sval) );

   return SCIP_OKAY;
}

/** sets string parameter of NLP */
SCIP_RETCODE SCIPnlpiSetStringPar(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char*           sval                /**< parameter value */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpisetstringpar != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpisetstringpar(set->scip, nlpi, problem, type, sval) );

   return SCIP_OKAY;
}
/**@} */

/** gets data of an NLPI */
SCIP_NLPIDATA* SCIPnlpiGetData(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);

   return nlpi->nlpidata;
}

/** gets NLP solver name */
const char* SCIPnlpiGetName(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);

   return nlpi->name;
}

/** gets NLP solver descriptions */
const char* SCIPnlpiGetDesc(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);

   return nlpi->description;
}

/** gets NLP solver priority */
int SCIPnlpiGetPriority(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);

   return nlpi->priority;
}

/** sets NLP solver priority */
void SCIPnlpiSetPriority(
   SCIP_NLPI*            nlpi,               /**< NLP interface structure */
   int                   priority            /**< new priority of NLPI */
   )
{
   assert(nlpi != NULL);

   nlpi->priority = priority;
}
