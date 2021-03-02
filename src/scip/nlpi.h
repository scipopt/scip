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

/**@file   nlpi.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for NLP solver interfaces
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_H__
#define __SCIP_NLPI_H__

#include "scip/type_nlpi.h"
#include "scip/type_misc.h"
#include "blockmemshell/memory.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup NLPIS
 *
 * @{
 */

/** creates an NLP solver interface */
SCIP_RETCODE SCIPnlpiCreate(
   SCIP_NLPI**                     nlpi,                        /**< pointer to NLP interface data structure */
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

/** copies an NLPI */
SCIP_DECL_NLPICOPY(SCIPnlpiCopy);

/** frees NLPI */
SCIP_RETCODE SCIPnlpiFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI**           nlpi                /**< pointer to NLPI data structure */
);

/** gets pointer for NLP solver
 * @return void pointer to solver
 */
SCIP_DECL_NLPIGETSOLVERPOINTER(SCIPnlpiGetSolverPointer);

/** creates a problem instance */
SCIP_DECL_NLPICREATEPROBLEM(SCIPnlpiCreateProblem);

/** frees a problem instance */
SCIP_DECL_NLPIFREEPROBLEM(SCIPnlpiFreeProblem);

/** gets pointer to solver-internal problem instance
 * @return void pointer to problem instance
 */
SCIP_DECL_NLPIGETPROBLEMPOINTER(SCIPnlpiGetProblemPointer);

/** add variables to nlpi */
SCIP_DECL_NLPIADDVARS(SCIPnlpiAddVars);

/** add constraints to nlpi */
SCIP_DECL_NLPIADDCONSTRAINTS(SCIPnlpiAddConstraints);

/** sets or overwrites objective, a minimization problem is expected */
SCIP_DECL_NLPISETOBJECTIVE(SCIPnlpiSetObjective);

/** change variable bounds */
SCIP_DECL_NLPICHGVARBOUNDS(SCIPnlpiChgVarBounds);

/** change constraint sides */
SCIP_DECL_NLPICHGCONSSIDES(SCIPnlpiChgConsSides);

/** delete a set of variables */
SCIP_DECL_NLPIDELVARSET(SCIPnlpiDelVarSet);

/** delete a set of constraints */
SCIP_DECL_NLPIDELCONSSET(SCIPnlpiDelConsSet);

/** changes or adds linear coefficients in a constraint or objective */
SCIP_DECL_NLPICHGLINEARCOEFS(SCIPnlpiChgLinearCoefs);

/** change the expression in the nonlinear part */
SCIP_DECL_NLPICHGEXPR(SCIPnlpiChgExpr);

/** change the constant offset in the objective */
SCIP_DECL_NLPICHGOBJCONSTANT(SCIPnlpiChgObjConstant);

/** sets initial guess for primal variables */
SCIP_DECL_NLPISETINITIALGUESS(SCIPnlpiSetInitialGuess);

/** tries to solve NLP */
SCIP_DECL_NLPISOLVE(SCIPnlpiSolve);
   
/** gives solution status */
SCIP_DECL_NLPIGETSOLSTAT(SCIPnlpiGetSolstat);

/** gives termination reason */
SCIP_DECL_NLPIGETTERMSTAT(SCIPnlpiGetTermstat);

/** gives primal and dual solution
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 */
SCIP_DECL_NLPIGETSOLUTION(SCIPnlpiGetSolution);

/** gives solve statistics */
SCIP_DECL_NLPIGETSTATISTICS(SCIPnlpiGetStatistics);

/** gives required size of a buffer to store a warmstart object */
SCIP_DECL_NLPIGETWARMSTARTSIZE(SCIPnlpiGetWarmstartSize);

/** stores warmstart information in buffer */
SCIP_DECL_NLPIGETWARMSTARTMEMO(SCIPnlpiGetWarmstartMemo);

/** sets warmstart information in solver */
SCIP_DECL_NLPISETWARMSTARTMEMO(SCIPnlpiSetWarmstartMemo);

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP */
SCIP_DECL_NLPIGETINTPAR(SCIPnlpiGetIntPar);
   
/** sets integer parameter of NLP */
SCIP_DECL_NLPISETINTPAR(SCIPnlpiSetIntPar);

/** gets floating point parameter of NLP */
SCIP_DECL_NLPIGETREALPAR(SCIPnlpiGetRealPar);

/** sets floating point parameter of NLP */
SCIP_DECL_NLPISETREALPAR(SCIPnlpiSetRealPar);

/** gets string parameter of NLP */
SCIP_DECL_NLPIGETSTRINGPAR(SCIPnlpiGetStringPar);

/** sets string parameter of NLP */
SCIP_DECL_NLPISETSTRINGPAR(SCIPnlpiSetStringPar);

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_H__ */
