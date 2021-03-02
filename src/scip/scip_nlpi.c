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

/**@file   scip_nlpi.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for NLP interfaces
 * @author Stefan Vigerske
 *
 * @todo check SCIP_STAGE_* switches
 * @todo allow for optional callbacks
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/scip_param.h"
#include "scip/scip_nlpi.h"
#include "scip/debug.h"
#include "scip/nlpi.h"
#include "scip/paramset.h"
#include "scip/set.h"
#include "scip/struct_scip.h"


/** method to call, when the priority of an NLPI was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdNlpiPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetSetPriorityNlpi() to mark the nlpis unsorted */
   SCIP_CALL( SCIPsetNlpiPriority(scip, (SCIP_NLPI*)paramdata, SCIPparamGetInt(param)) );

   return SCIP_OKAY;
}

/** create a NLPI and includes it into SCIP */
SCIP_RETCODE SCIPincludeNlpi(
   SCIP*                           scip,                        /**< SCIP data structure */
   const char*                     name,                        /**< name of NLP interface */
   const char*                     description,                 /**< description of NLP interface */
   int                             priority,                    /**< priority of NLP interface */
   SCIP_DECL_NLPICOPY              ((*nlpicopy)),               /**< copying an NLPI */
   SCIP_DECL_NLPIFREE              ((*nlpifree)),               /**< free NLPI user data */
   SCIP_DECL_NLPIGETSOLVERPOINTER  ((*nlpigetsolverpointer)),   /**< get solver pointer */
   SCIP_DECL_NLPICREATEPROBLEM     ((*nlpicreateproblem)),      /**< create a new problem instance */
   SCIP_DECL_NLPIFREEPROBLEM       ((*nlpifreeproblem)),        /**< free a problem instance */
   SCIP_DECL_NLPIGETPROBLEMPOINTER ((*nlpigetproblempointer)),  /**< get problem pointer */
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
   SCIP_DECL_NLPISETINITIALGUESS   ((*nlpisetinitialguess)),    /**< set initial guess for primal variables */
   SCIP_DECL_NLPISOLVE             ((*nlpisolve)),              /**< solve NLP */
   SCIP_DECL_NLPIGETSOLSTAT        ((*nlpigetsolstat)),         /**< get solution status */
   SCIP_DECL_NLPIGETTERMSTAT       ((*nlpigettermstat)),        /**< get termination status */
   SCIP_DECL_NLPIGETSOLUTION       ((*nlpigetsolution)),        /**< get solution */
   SCIP_DECL_NLPIGETSTATISTICS     ((*nlpigetstatistics)),      /**< get solve statistics */
   SCIP_DECL_NLPIGETWARMSTARTSIZE  ((*nlpigetwarmstartsize)),   /**< get size for warmstart object buffer */
   SCIP_DECL_NLPIGETWARMSTARTMEMO  ((*nlpigetwarmstartmemo)),   /**< get warmstart object */
   SCIP_DECL_NLPISETWARMSTARTMEMO  ((*nlpisetwarmstartmemo)),   /**< set warmstart object */
   SCIP_DECL_NLPIGETINTPAR         ((*nlpigetintpar)),          /**< get value of integer parameter */
   SCIP_DECL_NLPISETINTPAR         ((*nlpisetintpar)),          /**< set value of integer parameter */
   SCIP_DECL_NLPIGETREALPAR        ((*nlpigetrealpar)),         /**< get value of floating point parameter */
   SCIP_DECL_NLPISETREALPAR        ((*nlpisetrealpar)),         /**< set value of floating point parameter */
   SCIP_DECL_NLPIGETSTRINGPAR      ((*nlpigetstringpar)),       /**< get value of string parameter */
   SCIP_DECL_NLPISETSTRINGPAR      ((*nlpisetstringpar)),       /**< set value of string parameter */
   SCIP_NLPIDATA*                  nlpidata                     /**< NLP interface local data */
   )
{
   SCIP_NLPI* nlpi = NULL;
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeNlpi", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether NLPI of given name is already present */
   if( SCIPfindNlpi(scip, name) != NULL )
   {
      SCIPerrorMessage("NLPI <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPnlpiCreate(&nlpi, name, description, priority,
      nlpicopy, nlpifree, nlpigetsolverpointer,
      nlpicreateproblem, nlpifreeproblem, nlpigetproblempointer,
      nlpiaddvars, nlpiaddconstraints, nlpisetobjective, nlpichgvarbounds, nlpichgconssides, nlpidelvarset, nlpidelconsset, nlpichglinearcoefs, nlpichgexpr, nlpichgobjconstant,
      nlpisetinitialguess, nlpisolve, nlpigetsolstat, nlpigettermstat, nlpigetsolution, nlpigetstatistics,
      nlpigetwarmstartsize, nlpigetwarmstartmemo, nlpisetwarmstartmemo,
      nlpigetintpar, nlpisetintpar, nlpigetrealpar, nlpisetrealpar, nlpigetstringpar, nlpisetstringpar,
      nlpidata) );
   assert(nlpi != NULL);

   SCIP_CALL( SCIPsetIncludeNlpi(scip->set, nlpi) );

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "nlpi/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of NLPI <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, paramname, paramdesc,
         NULL, FALSE, SCIPnlpiGetPriority(nlpi), INT_MIN/4, INT_MAX/4,
         paramChgdNlpiPriority, (SCIP_PARAMDATA*)nlpi) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** returns the NLPI of the given name, or NULL if not existing */
SCIP_NLPI* SCIPfindNlpi(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of NLPI */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindNlpi(scip->set, name);
}

/** returns the array of currently available NLPIs (sorted by priority) */
SCIP_NLPI** SCIPgetNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortNlpis(scip->set);

   return scip->set->nlpis;
}

/** returns the number of currently available NLPIs */
int SCIPgetNNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nnlpis;
}

/** sets the priority of an NLPI */
SCIP_RETCODE SCIPsetNlpiPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< NLPI */
   int                   priority            /**< new priority of the NLPI */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSetPriorityNlpi(scip->set, nlpi, priority);

   return SCIP_OKAY;
}

/** gets pointer for NLP solver
 * @return void pointer to solver
 */
SCIP_DECL_NLPIGETSOLVERPOINTER(SCIPgetNlpiSolverPointer)
{
   assert(scip != NULL);

   return SCIPnlpiGetSolverPointer(scip->set, nlpi);
}

/** creates a problem instance */
SCIP_DECL_NLPICREATEPROBLEM(SCIPcreateNlpiProblem)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiCreateProblem(scip->set, nlpi, problem, name) );

   return SCIP_OKAY;
}

/** frees a problem instance */
SCIP_DECL_NLPIFREEPROBLEM(SCIPfreeNlpiProblem)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiFreeProblem(scip->set, nlpi, problem) );

   return SCIP_OKAY;
}

/** gets pointer to solver-internal problem instance
 * @return void pointer to problem instance
 */
SCIP_DECL_NLPIGETPROBLEMPOINTER(SCIPgetNlpiProblemPointer)
{
   assert(scip != NULL);

   return SCIPnlpiGetProblemPointer(scip->set, nlpi, problem);
}

/** add variables to nlpi */
SCIP_DECL_NLPIADDVARS(SCIPaddNlpiVars)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiAddVars(scip->set, nlpi, problem, nvars, lbs, ubs, varnames) );

   return SCIP_OKAY;
}

/** add constraints to nlpi */
SCIP_DECL_NLPIADDCONSTRAINTS(SCIPaddNlpiConstraints)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiAddConstraints(scip->set, nlpi, problem, nconss, lhss, rhss, nlininds, lininds, linvals, exprs, names) );

   return SCIP_OKAY;
}

/** sets or overwrites objective, a minimization problem is expected */
SCIP_DECL_NLPISETOBJECTIVE(SCIPsetNlpiObjective)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiSetObjective(scip->set, nlpi, problem, nlins, lininds, linvals, expr, constant) );

   return SCIP_OKAY;
}

/** change variable bounds */
SCIP_DECL_NLPICHGVARBOUNDS(SCIPchgNlpiVarBounds)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiChgVarBounds(scip->set, nlpi, problem, nvars, indices, lbs, ubs) );

   return SCIP_OKAY;
}

/** change constraint sides */
SCIP_DECL_NLPICHGCONSSIDES(SCIPchgNlpiConsSides)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiChgConsSides(scip->set, nlpi, problem, nconss, indices, lhss, rhss) );

   return SCIP_OKAY;
}

/** delete a set of variables */
SCIP_DECL_NLPIDELVARSET(SCIPdelNlpiVarSet)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiDelVarSet(scip->set, nlpi, problem, dstats, dstatssize) );

   return SCIP_OKAY;
}

/** delete a set of constraints */
SCIP_DECL_NLPIDELCONSSET(SCIPdelNlpiConsSet)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiDelConsSet(scip->set, nlpi, problem, dstats, dstatssize) );

   return SCIP_OKAY;
}

/** changes or adds linear coefficients in a constraint or objective */
SCIP_DECL_NLPICHGLINEARCOEFS(SCIPchgNlpiLinearCoefs)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiChgLinearCoefs(scip->set, nlpi, problem, idx, nvals, varidxs, vals) );

   return SCIP_OKAY;
}

/** change the expression in the nonlinear part */
SCIP_DECL_NLPICHGEXPR(SCIPchgNlpiExpr)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiChgExpr(scip->set, nlpi, problem, idxcons, expr) );

   return SCIP_OKAY;
}

/** change the constant offset in the objective */
SCIP_DECL_NLPICHGOBJCONSTANT(SCIPchgNlpiObjConstant)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiChgObjConstant(scip->set, nlpi, problem, objconstant) );

   return SCIP_OKAY;
}

/** sets initial guess for primal variables */
SCIP_DECL_NLPISETINITIALGUESS(SCIPsetNlpiInitialGuess)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiSetInitialGuess(scip->set, nlpi, problem, primalvalues, consdualvalues, varlbdualvalues, varubdualvalues) );

   return SCIP_OKAY;
}

/** tries to solve NLP */
SCIP_DECL_NLPISOLVE(SCIPsolveNlpi)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiSolve(scip->set, nlpi, problem) );

   return SCIP_OKAY;
}

/** gives solution status */
SCIP_DECL_NLPIGETSOLSTAT(SCIPgetNlpiSolstat)
{
   assert(scip != NULL);

   return SCIPnlpiGetSolstat(scip->set, nlpi, problem);
}

/** gives termination reason */
SCIP_DECL_NLPIGETTERMSTAT(SCIPgetNlpiTermstat)
{
   assert(scip != NULL);

   return SCIPnlpiGetTermstat(scip->set, nlpi, problem);
}

/** gives primal and dual solution
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 */
SCIP_DECL_NLPIGETSOLUTION(SCIPgetNlpiSolution)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiGetSolution(scip->set, nlpi, problem, primalvalues, consdualvalues, varlbdualvalues, varubdualvalues, objval) );

   return SCIP_OKAY;
}

/** gives solve statistics */
SCIP_DECL_NLPIGETSTATISTICS(SCIPgetNlpiStatistics)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiGetStatistics(scip->set, nlpi, problem, statistics) );

   return SCIP_OKAY;
}

/** gives required size of a buffer to store a warmstart object */
SCIP_DECL_NLPIGETWARMSTARTSIZE(SCIPgetNlpiWarmstartSize)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiGetWarmstartSize(scip->set, nlpi, problem, size) );

   return SCIP_OKAY;
}

/** stores warmstart information in buffer */
SCIP_DECL_NLPIGETWARMSTARTMEMO(SCIPgetNlpiWarmstartMemo)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiGetWarmstartMemo(scip->set, nlpi, problem, buffer) );

   return SCIP_OKAY;
}

/** sets warmstart information in solver */
SCIP_DECL_NLPISETWARMSTARTMEMO(SCIPsetNlpiWarmstartMemo)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiSetWarmstartMemo(scip->set, nlpi, problem, buffer) );

   return SCIP_OKAY;
}

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP */
SCIP_DECL_NLPIGETINTPAR(SCIPgetNlpiIntPar)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiGetIntPar(scip->set, nlpi, problem, type, ival) );

   return SCIP_OKAY;
}

/** sets integer parameter of NLP */
SCIP_DECL_NLPISETINTPAR(SCIPsetNlpiIntPar)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiSetIntPar(scip->set, nlpi, problem, type, ival) );

   return SCIP_OKAY;
}

/** gets floating point parameter of NLP
 * if problem is NULL and type == SCIP_NLPPAR_INFINITY, then gets solver-wide value for infinity */
SCIP_DECL_NLPIGETREALPAR(SCIPgetNlpiRealPar)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiGetRealPar(scip->set, nlpi, problem, type, dval) );

   return SCIP_OKAY;
}

/** sets floating point parameter of NLP
 * if problem is NULL and type == SCIP_NLPPAR_INFINITY, then sets solver-wide value for infinity */
SCIP_DECL_NLPISETREALPAR(SCIPsetNlpiRealPar)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiSetRealPar(scip->set, nlpi, problem, type, dval) );

   return SCIP_OKAY;
}

/** gets string parameter of NLP */
SCIP_DECL_NLPIGETSTRINGPAR(SCIPgetNlpiStringPar)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiGetStringPar(scip->set, nlpi, problem, type, sval) );

   return SCIP_OKAY;
}

/** sets string parameter of NLP */
SCIP_DECL_NLPISETSTRINGPAR(SCIPsetNlpiStringPar)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiSetStringPar(scip->set, nlpi, problem, type, sval) );

   return SCIP_OKAY;
}
