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

/**@file   scip_nlpi.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for NLPI solver interfaces
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_NLPI_H__
#define __SCIP_SCIP_NLPI_H__

#include "scip/type_nlpi.h"
#include "scip/type_misc.h"
#include "scip/type_lp.h"
#include "blockmemshell/memory.h"
#include "scip/pub_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicNLPIInterfaceMethods
 *
 * @{
 */

/** includes an NLPI in SCIP */
SCIP_EXPORT
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
   );

/** returns the NLPI of the given name, or NULL if not existing */
SCIP_EXPORT
SCIP_NLPI* SCIPfindNlpi(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of NLPI */
   );

/** returns the array of currently available NLPIs (sorted by priority) */
SCIP_EXPORT
SCIP_NLPI** SCIPgetNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available NLPIs */
SCIP_EXPORT
int SCIPgetNNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of an NLPI */
SCIP_EXPORT
SCIP_RETCODE SCIPsetNlpiPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< NLPI */
   int                   priority            /**< new priority of the NLPI */
   );

/** gets pointer for NLP solver
 * @return void pointer to solver
 */
SCIP_EXPORT
SCIP_DECL_NLPIGETSOLVERPOINTER(SCIPgetNlpiSolverPointer);

/** creates an empty problem instance */
SCIP_EXPORT
SCIP_DECL_NLPICREATEPROBLEM(SCIPcreateNlpiProblem);

/** frees a problem instance */
SCIP_EXPORT
SCIP_DECL_NLPIFREEPROBLEM(SCIPfreeNlpiProblem);

/** gets pointer to solver-internal problem instance
 * @return void pointer to problem instance
 */
SCIP_EXPORT
SCIP_DECL_NLPIGETPROBLEMPOINTER(SCIPgetNlpiProblemPointer);

/** add variables to nlpi */
SCIP_EXPORT
SCIP_DECL_NLPIADDVARS(SCIPaddNlpiVars);

/** add constraints to nlpi */
SCIP_EXPORT
SCIP_DECL_NLPIADDCONSTRAINTS(SCIPaddNlpiConstraints);

/** sets or overwrites objective, a minimization problem is expected */
SCIP_EXPORT
SCIP_DECL_NLPISETOBJECTIVE(SCIPsetNlpiObjective);

/** change variable bounds */
SCIP_EXPORT
SCIP_DECL_NLPICHGVARBOUNDS(SCIPchgNlpiVarBounds);

/** change constraint sides */
SCIP_EXPORT
SCIP_DECL_NLPICHGCONSSIDES(SCIPchgNlpiConsSides);

/** delete a set of variables */
SCIP_EXPORT
SCIP_DECL_NLPIDELVARSET(SCIPdelNlpiVarSet);

/** delete a set of constraints */
SCIP_EXPORT
SCIP_DECL_NLPIDELCONSSET(SCIPdelNlpiConsSet);

/** changes or adds linear coefficients in a constraint or objective */
SCIP_EXPORT
SCIP_DECL_NLPICHGLINEARCOEFS(SCIPchgNlpiLinearCoefs);

/** change the expression in the nonlinear part */
SCIP_EXPORT
SCIP_DECL_NLPICHGEXPR(SCIPchgNlpiExpr);

/** change the constant offset in the objective */
SCIP_EXPORT
SCIP_DECL_NLPICHGOBJCONSTANT(SCIPchgNlpiObjConstant);

/** sets initial guess for primal variables */
SCIP_EXPORT
SCIP_DECL_NLPISETINITIALGUESS(SCIPsetNlpiInitialGuess);

/** tries to solve NLP */
SCIP_EXPORT
SCIP_DECL_NLPISOLVE(SCIPsolveNlpi);

/** gives solution status */
SCIP_EXPORT
SCIP_DECL_NLPIGETSOLSTAT(SCIPgetNlpiSolstat);

/** gives termination reason */
SCIP_EXPORT
SCIP_DECL_NLPIGETTERMSTAT(SCIPgetNlpiTermstat);

/** gives primal and dual solution
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 */
SCIP_EXPORT
SCIP_DECL_NLPIGETSOLUTION(SCIPgetNlpiSolution);

/** gives solve statistics */
SCIP_EXPORT
SCIP_DECL_NLPIGETSTATISTICS(SCIPgetNlpiStatistics);

/** gives required size of a buffer to store a warmstart object */
SCIP_EXPORT
SCIP_DECL_NLPIGETWARMSTARTSIZE(SCIPgetNlpiWarmstartSize);

/** stores warmstart information in buffer */
SCIP_EXPORT
SCIP_DECL_NLPIGETWARMSTARTMEMO(SCIPgetNlpiWarmstartMemo);

/** sets warmstart information in solver */
SCIP_EXPORT
SCIP_DECL_NLPISETWARMSTARTMEMO(SCIPsetNlpiWarmstartMemo);

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP */
SCIP_EXPORT
SCIP_DECL_NLPIGETINTPAR(SCIPgetNlpiIntPar);

/** sets integer parameter of NLP */
SCIP_EXPORT
SCIP_DECL_NLPISETINTPAR(SCIPsetNlpiIntPar);

/** gets floating point parameter of NLP
 * if problem is NULL and type == SCIP_NLPPAR_INFINITY, then gets solver-wide value for infinity */
SCIP_EXPORT
SCIP_DECL_NLPIGETREALPAR(SCIPgetNlpiRealPar);

/** sets floating point parameter of NLP
 * if problem is NULL and type == SCIP_NLPPAR_INFINITY, then sets solver-wide value for infinity */
SCIP_EXPORT
SCIP_DECL_NLPISETREALPAR(SCIPsetNlpiRealPar);

/** gets string parameter of NLP */
SCIP_EXPORT
SCIP_DECL_NLPIGETSTRINGPAR(SCIPgetNlpiStringPar);

/** sets string parameter of NLP */
SCIP_EXPORT
SCIP_DECL_NLPISETSTRINGPAR(SCIPsetNlpiStringPar);

/** @} */

/**@name Convenience methods to setup and update an NLPI problem using NLROWS
 *
 * These methods can be used, for example, to create a NLPI problem that contains only the convex rows of the SCIP NLP relaxation.
 * @{
 */

/** creates a NLPI problem from given nonlinear rows
 *
 * The function computes for each variable the number of non-linear occurrences and stores it in the nlscore array.
 *
 * @note the first row corresponds always to the cutoff row (even if cutoffbound is SCIPinfinity(scip))
 **/
SCIP_EXPORT
SCIP_RETCODE SCIPcreateNlpiProblemFromNlRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM**    nlpiprob,           /**< buffer to store pointer to created nlpi problem */
   const char*           name,               /**< name to give to problem */
   SCIP_NLROW**          nlrows,             /**< nonlinear rows */
   int                   nnlrows,            /**< number of nonlinear rows */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpiprob */
   SCIP_HASHMAP*         nlrow2idx,          /**< empty hash map to store mapping between variables and indices in nlpiprob, can be NULL */
   SCIP_Real*            nlscore,            /**< array to store the score of each nonlinear variable (NULL if not needed) */
   SCIP_Real             cutoffbound,        /**< cutoff bound */
   SCIP_Bool             setobj,             /**< whether the objective function should be set to one of the SCIP problem */
   SCIP_Bool             onlyconvex          /**< filter only for convex constraints */
   );

/** updates variable bounds and the cutoff row in a NLPI problem
 *
 * The NLPI problem must have been setup by SCIPcreateNlpiProblemFromNlRows().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPupdateNlpiProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem representing the convex NLP relaxation */
   SCIP_HASHMAP*         var2nlpiidx,        /**< mapping between variables and nlpi indices */
   SCIP_VAR**            nlpivars,           /**< array containing all variables of the nlpi */
   int                   nlpinvars,          /**< total number of nlpi variables */
   SCIP_Real             cutoffbound         /**< new cutoff bound */
   );

/** adds SCIP_ROWs to a NLPI problem */
SCIP_EXPORT
SCIP_RETCODE SCIPaddNlpiProblemRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpiprob */
   SCIP_ROW**            rows,               /**< rows to add */
   int                   nrows               /**< number of rows to add */
   );

/** adds SCIP_NLROWs to a NLPI problem */
SCIP_EXPORT
SCIP_RETCODE SCIPaddNlpiProblemNlRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpiprob */
   SCIP_NLROW**          nlrows,             /**< rows to add */
   int                   nnlrows             /**< number of rows to add */
   );

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_SCIP_NLPI_H__ */
