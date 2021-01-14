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

#include "scip/pub_message.h"
#include "nlpi/nlpi.h"
#include "nlpi/struct_nlpi.h"
#include "blockmemshell/memory.h"

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
{  /*lint --e{715}*/
   assert(nlpi != NULL);

   assert(name != NULL);
   assert(description != NULL);
   assert(nlpicopy != NULL);
   assert(nlpifree != NULL);
   assert(nlpigetsolverpointer != NULL);
   assert(nlpicreateproblem != NULL);
   assert(nlpifreeproblem != NULL);
   assert(nlpigetproblempointer != NULL);
   assert(nlpiaddvars != NULL);
   assert(nlpiaddconstraints != NULL);
   assert(nlpisetobjective != NULL);
   assert(nlpichgvarbounds != NULL);
   assert(nlpichgconssides != NULL);
   assert(nlpidelconsset != NULL);
   assert(nlpichglinearcoefs != NULL);
   assert(nlpichgobjconstant != NULL);
   assert(nlpisetinitialguess != NULL);
   assert(nlpisolve != NULL);
   assert(nlpigetsolstat != NULL);
   assert(nlpigettermstat != NULL);
   assert(nlpigetsolution != NULL);
   assert(nlpigetstatistics != NULL);
   assert(nlpigetwarmstartsize != NULL);
   assert(nlpigetwarmstartmemo != NULL);
   assert(nlpisetwarmstartmemo != NULL);
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
   (*nlpi)->nlpigetwarmstartsize = nlpigetwarmstartsize;
   (*nlpi)->nlpigetwarmstartmemo = nlpigetwarmstartmemo;
   (*nlpi)->nlpisetwarmstartmemo = nlpisetwarmstartmemo;
   (*nlpi)->nlpigetintpar = nlpigetintpar;
   (*nlpi)->nlpisetintpar = nlpisetintpar;
   (*nlpi)->nlpigetrealpar = nlpigetrealpar;
   (*nlpi)->nlpisetrealpar = nlpisetrealpar;
   (*nlpi)->nlpigetstringpar = nlpigetstringpar;
   (*nlpi)->nlpisetstringpar = nlpisetstringpar;
   (*nlpi)->nlpidata = nlpidata;

   return SCIP_OKAY;
}

/** copies an NLPI */
SCIP_DECL_NLPICOPY(SCIPnlpiCopy)
{
   assert(sourcenlpi != NULL);
   assert(sourcenlpi->nlpicopy != NULL);

   SCIP_CALL( sourcenlpi->nlpicopy(scip, sourcenlpi) );

   return SCIP_OKAY;
}

/** frees NLPI */
SCIP_RETCODE SCIPnlpiFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI**           nlpi                /**< pointer to NLPI data structure */
)
{
   assert(nlpi  != NULL);
   assert(*nlpi != NULL);

   if( (*nlpi)->nlpifree != NULL )
   {
      SCIP_CALL( (*nlpi)->nlpifree(scip, *nlpi, &(*nlpi)->nlpidata) );
      assert((*nlpi)->nlpidata == NULL);
   }
   BMSfreeMemoryArray(&(*nlpi)->name);
   BMSfreeMemoryArray(&(*nlpi)->description);
   BMSfreeMemory(nlpi);

   assert(*nlpi == NULL);

   return SCIP_OKAY;
}

/** gets pointer for NLP solver
 * @return void pointer to solver
 */
SCIP_DECL_NLPIGETSOLVERPOINTER(SCIPnlpiGetSolverPointer)
{
   assert(nlpi != NULL);
   assert(nlpi->nlpigetsolverpointer != NULL);

   return nlpi->nlpigetsolverpointer(scip, nlpi);
}

/** creates a problem instance */
SCIP_DECL_NLPICREATEPROBLEM(SCIPnlpiCreateProblem)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpicreateproblem != NULL);
   assert(problem != NULL);

   return nlpi->nlpicreateproblem(scip, nlpi, problem, name);
}

/** frees a problem instance */
SCIP_DECL_NLPIFREEPROBLEM(SCIPnlpiFreeProblem)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpifreeproblem != NULL);
   assert(problem != NULL);

   return nlpi->nlpifreeproblem(scip, nlpi, problem);
}

/** gets pointer to solver-internal problem instance
 * @return void pointer to problem instance
 */
SCIP_DECL_NLPIGETPROBLEMPOINTER(SCIPnlpiGetProblemPointer)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpigetproblempointer != NULL);
   assert(problem != NULL);

   return nlpi->nlpigetproblempointer(scip, nlpi, problem);
}

/** add variables to nlpi */
SCIP_DECL_NLPIADDVARS(SCIPnlpiAddVars)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpiaddvars != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpiaddvars(scip, nlpi, problem, nvars, lbs, ubs, varnames) );

   return SCIP_OKAY;
}

/** add constraints to nlpi */
SCIP_DECL_NLPIADDCONSTRAINTS(SCIPnlpiAddConstraints)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpiaddconstraints != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpiaddconstraints(scip, nlpi, problem, nconss, lhss, rhss, nlininds, lininds, linvals, exprs, names) );

   return SCIP_OKAY;
}

/** sets or overwrites objective, a minimization problem is expected */
SCIP_DECL_NLPISETOBJECTIVE(SCIPnlpiSetObjective)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpisetobjective != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpisetobjective(scip, nlpi, problem, nlins, lininds, linvals, expr, constant) );

   return SCIP_OKAY;
}

/** change variable bounds */
SCIP_DECL_NLPICHGVARBOUNDS(SCIPnlpiChgVarBounds)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpichgvarbounds != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpichgvarbounds(scip, nlpi, problem, nvars, indices, lbs, ubs) );

   return SCIP_OKAY;
}

/** change constraint bounds */
SCIP_DECL_NLPICHGCONSSIDES(SCIPnlpiChgConsSides)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpichgconssides != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpichgconssides(scip, nlpi, problem, nconss, indices, lhss, rhss) );

   return SCIP_OKAY;
}

/** delete a set of variables */
SCIP_DECL_NLPIDELVARSET(SCIPnlpiDelVarSet)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpidelvarset != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpidelvarset(scip, nlpi, problem, dstats, dstatssize) );

   return SCIP_OKAY;
}

/** delete a set of constraints */
SCIP_DECL_NLPIDELCONSSET(SCIPnlpiDelConsSet)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpidelconsset != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpidelconsset(scip, nlpi, problem, dstats, dstatssize) );

   return SCIP_OKAY;
}

/** changes or adds linear coefficients in a constraint or objective */
SCIP_DECL_NLPICHGLINEARCOEFS(SCIPnlpiChgLinearCoefs)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpichglinearcoefs != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpichglinearcoefs(scip, nlpi, problem, idx, nvals, varidxs, vals) );

   return SCIP_OKAY;
}

/** change the expression in the nonlinear part */
SCIP_DECL_NLPICHGEXPR(SCIPnlpiChgExpr)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpichgexpr != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpichgexpr(scip, nlpi, problem, idxcons, expr) );

   return SCIP_OKAY;
}

/** change the constant offset in the objective */
SCIP_DECL_NLPICHGOBJCONSTANT(SCIPnlpiChgObjConstant)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpichgobjconstant != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpichgobjconstant(scip, nlpi, problem, objconstant) );

   return SCIP_OKAY;
}

/** sets initial guess for primal variables */
SCIP_DECL_NLPISETINITIALGUESS(SCIPnlpiSetInitialGuess)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpisetinitialguess != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpisetinitialguess(scip, nlpi, problem, primalvalues, consdualvalues, varlbdualvalues, varubdualvalues) );

   return SCIP_OKAY;
}

/** tries to solve NLP */
SCIP_DECL_NLPISOLVE(SCIPnlpiSolve)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpisolve != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpisolve(scip, nlpi, problem) );

   return SCIP_OKAY;
}

/** gives solution status, return: Solution Status */
SCIP_DECL_NLPIGETSOLSTAT(SCIPnlpiGetSolstat)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpigetsolstat != NULL);
   assert(problem != NULL);

   return nlpi->nlpigetsolstat(scip, nlpi, problem);
}

/** gives termination reason; return: Termination Status */
SCIP_DECL_NLPIGETTERMSTAT(SCIPnlpiGetTermstat)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpigettermstat != NULL);
   assert(problem != NULL);

   return nlpi->nlpigettermstat(scip, nlpi, problem);
}

/** gives primal and dual solution
  * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
  */
SCIP_DECL_NLPIGETSOLUTION(SCIPnlpiGetSolution)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpigetsolution != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpigetsolution(scip, nlpi, problem, primalvalues, consdualvalues, varlbdualvalues, varubdualvalues, objval) );

   return SCIP_OKAY;
}

/** gives solve statistics */
SCIP_DECL_NLPIGETSTATISTICS(SCIPnlpiGetStatistics)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpigetstatistics != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpigetstatistics(scip, nlpi, problem, statistics) );

   return SCIP_OKAY;
}

/** gives required size of a buffer to store a warmstart object */
SCIP_DECL_NLPIGETWARMSTARTSIZE(SCIPnlpiGetWarmstartSize)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpigetwarmstartsize != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpigetwarmstartsize(scip, nlpi, problem, size) );

   return SCIP_OKAY;
}

/** stores warmstart information in buffer */
SCIP_DECL_NLPIGETWARMSTARTMEMO(SCIPnlpiGetWarmstartMemo)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpigetwarmstartmemo != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpigetwarmstartmemo(scip, nlpi, problem, buffer) );

   return SCIP_OKAY;
}

/** sets warmstart information in solver */
SCIP_DECL_NLPISETWARMSTARTMEMO(SCIPnlpiSetWarmstartMemo)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpisetwarmstartmemo != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpisetwarmstartmemo(scip, nlpi, problem, buffer) );

   return SCIP_OKAY;
}

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP */
SCIP_DECL_NLPIGETINTPAR(SCIPnlpiGetIntPar)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpigetintpar != NULL);
   assert(problem != NULL);
   assert(ival    != NULL);

   SCIP_CALL( nlpi->nlpigetintpar(scip, nlpi, problem, type, ival) );

   return SCIP_OKAY;
}

/** sets integer parameter of NLP */
SCIP_DECL_NLPISETINTPAR(SCIPnlpiSetIntPar)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpisetintpar != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpisetintpar(scip, nlpi, problem, type, ival) );

   return SCIP_OKAY;
}

/** gets floating point parameter of NLP */
SCIP_DECL_NLPIGETREALPAR(SCIPnlpiGetRealPar)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpigetrealpar != NULL);
   assert(problem != NULL);
   assert(dval    != NULL);

   SCIP_CALL( nlpi->nlpigetrealpar(scip, nlpi, problem, type, dval) );

   return SCIP_OKAY;
}

/** sets floating point parameter of NLP */
SCIP_DECL_NLPISETREALPAR(SCIPnlpiSetRealPar)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpisetrealpar != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpisetrealpar(scip, nlpi, problem, type, dval) );

   return SCIP_OKAY;
}

/** gets string parameter of NLP */
SCIP_DECL_NLPIGETSTRINGPAR(SCIPnlpiGetStringPar)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpigetstringpar != NULL);
   assert(problem != NULL);
   assert(sval    != NULL);

   SCIP_CALL( nlpi->nlpigetstringpar(scip, nlpi, problem, type, sval) );

   return SCIP_OKAY;
}

/** sets string parameter of NLP */
SCIP_DECL_NLPISETSTRINGPAR(SCIPnlpiSetStringPar)
{
   assert(nlpi    != NULL);
   assert(nlpi->nlpisetstringpar != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpisetstringpar(scip, nlpi, problem, type, sval) );

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

/** creates an NLP statistics structure */
SCIP_RETCODE SCIPnlpStatisticsCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NLPSTATISTICS**  statistics          /**< pointer where to store NLP statistics structure */
   )
{
   assert(blkmem != NULL);
   assert(statistics != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, statistics) );

   (*statistics)->niterations = -1;
   (*statistics)->totaltime = -1.0;

   return SCIP_OKAY;
}

/** frees an NLP statistics structure */
void SCIPnlpStatisticsFree(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NLPSTATISTICS**  statistics          /**< pointer where to store NLP statistics structure */
   )
{
   assert(blkmem != NULL);
   assert(statistics != NULL);
   assert(*statistics != NULL);

   BMSfreeBlockMemory(blkmem, statistics);

   assert(*statistics == NULL);
}

/** gets the number of iterations from an NLP statistics structure */
int SCIPnlpStatisticsGetNIterations(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
   )
{
   assert(statistics != NULL);

   return statistics->niterations;
}

/** gets the total time from an NLP statistics structure */
SCIP_Real SCIPnlpStatisticsGetTotalTime(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
   )
{
   assert(statistics != NULL);

   return statistics->totaltime;
}

/** sets the number of iterations in an NLP statistics structure */
void SCIPnlpStatisticsSetNIterations(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   int                   niterations         /**< number of iterations to store */
   )
{
   assert(statistics != NULL);
   statistics->niterations = niterations;
}

/** sets the total time in an NLP statistics structure */
void SCIPnlpStatisticsSetTotalTime(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   SCIP_Real             totaltime           /**< solution time to store */
   )
{
   assert(statistics != NULL);
   statistics->totaltime = totaltime;
}
