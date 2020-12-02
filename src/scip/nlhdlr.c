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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr.c
 * @brief  functions for nonlinearity handlers of nonlinear constraint handler
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#include <assert.h>

#include "scip/intervalarith.h"
#include "scip/pub_nlhdlr.h"
#include "scip/nlhdlr.h"
#include "scip/struct_nlhdlr.h"
#include "scip/scip_timing.h"

/* nlhdlr public API functions from pub_nlhdlr.h */

/** set the copy handler callback of a nonlinear handler */
void SCIPnlhdlrSetCopyHdlr(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRCOPYHDLR((*copy))         /**< copy callback (can be NULL) */
)
{
   assert(nlhdlr != NULL);

   nlhdlr->copyhdlr = copy;
}

/** set the nonlinear handler callback to free the nonlinear handler data */
void SCIPnlhdlrSetFreeHdlrData(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRFREEHDLRDATA((*freehdlrdata)) /**< handler free callback (can be NULL) */
)
{
   assert(nlhdlr != NULL);

   nlhdlr->freehdlrdata = freehdlrdata;
}

/** set the nonlinear handler callback to free expression specific data of nonlinear handler */
void SCIPnlhdlrSetFreeExprData(
   SCIP_NLHDLR*          nlhdlr,            /**< nonlinear handler */
   SCIP_DECL_NLHDLRFREEEXPRDATA((*freeexprdata)) /**< nonlinear handler expression data free callback (can be NULL if data does not need to be freed) */
)
{
   assert(nlhdlr != NULL);

   nlhdlr->freeexprdata = freeexprdata;
}

/** set the initialization and deinitialization callback of a nonlinear handler */
void SCIPnlhdlrSetInitExit(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRINIT((*init)),            /**< initialization callback (can be NULL) */
   SCIP_DECL_NLHDLREXIT((*exit_))            /**< deinitialization callback (can be NULL) */
)
{
   assert(nlhdlr != NULL);

   nlhdlr->init = init;
   nlhdlr->exit = exit_;
}

/** set the propagation callbacks of a nonlinear handler */
void SCIPnlhdlrSetProp(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRINTEVAL((*inteval)),      /**< interval evaluation callback (can be NULL) */
   SCIP_DECL_NLHDLRREVERSEPROP((*reverseprop)) /**< reverse propagation callback (can be NULL) */
)
{
   assert(nlhdlr != NULL);

   nlhdlr->inteval = inteval;
   nlhdlr->reverseprop = reverseprop;
}

/** set the enforcement callbacks of a nonlinear handler */
void SCIPnlhdlrSetSepa(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRINITSEPA((*initsepa)),    /**< separation initialization callback (can be NULL) */
   SCIP_DECL_NLHDLRENFO((*enfo)),            /**< enforcement callback (can be NULL if estimate is not NULL) */
   SCIP_DECL_NLHDLRESTIMATE((*estimate)),    /**< estimation callback (can be NULL if sepa is not NULL) */
   SCIP_DECL_NLHDLREXITSEPA((*exitsepa))     /**< separation deinitialization callback (can be NULL) */
)
{
   assert(nlhdlr != NULL);
   assert(enfo != NULL || estimate != NULL);

   nlhdlr->initsepa = initsepa;
   nlhdlr->enfo = enfo;
   nlhdlr->estimate = estimate;
   nlhdlr->exitsepa = exitsepa;
}

/** gives name of nonlinear handler */
const char* SCIPnlhdlrGetName(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
)
{
   assert(nlhdlr != NULL);

   return nlhdlr->name;
}

/** gives description of nonlinear handler, can be NULL */
const char* SCIPnlhdlrGetDesc(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
)
{
   assert(nlhdlr != NULL);

   return nlhdlr->desc;
}

/** gives detection priority of nonlinear handler */
int SCIPnlhdlrGetDetectPriority(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
)
{
   assert(nlhdlr != NULL);

   return nlhdlr->detectpriority;
}

/** gives enforcement priority of nonlinear handler */
int SCIPnlhdlrGetEnfoPriority(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
)
{
   assert(nlhdlr != NULL);

   return nlhdlr->enfopriority;
}

/** gives handler data of nonlinear handler */
SCIP_NLHDLRDATA* SCIPnlhdlrGetData(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
)
{
   assert(nlhdlr != NULL);

   return nlhdlr->data;
}

/** returns whether nonlinear handler implements the interval evaluation callback */
SCIP_Bool SCIPnlhdlrHasIntEval(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
)
{
   return nlhdlr->inteval != NULL;
}

/** returns whether nonlinear handler implements the reverse propagation callback */
SCIP_Bool SCIPnlhdlrHasReverseProp(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
)
{
   return nlhdlr->reverseprop != NULL;
}

/** returns whether nonlinear handler implements the separation initialization callback */
SCIP_Bool SCIPnlhdlrHasInitSepa(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
)
{
   return nlhdlr->initsepa != NULL;
}

/** returns whether nonlinear handler implements the separation deinitialization callback */
SCIP_Bool SCIPnlhdlrHasExitSepa(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
)
{
   return nlhdlr->exitsepa != NULL;
}

/** returns whether nonlinear handler implements the enforcement callback */
SCIP_Bool SCIPnlhdlrHasEnfo(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
)
{
   return nlhdlr->enfo != NULL;
}

/** returns whether nonlinear handler implements the estimator callback */
SCIP_Bool SCIPnlhdlrHasEstimate(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
)
{
   return nlhdlr->estimate != NULL;
}

/* nlhdlr private API functions from pub_nlhdlr.h */

/** call the detect callback of a nonlinear handler */
SCIP_DECL_NLHDLRDETECT(SCIPnlhdlrDetect)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->detect != NULL);
   assert(nlhdlr->detecttime != NULL);
   assert(participating != NULL);

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->detecttime) );
   SCIP_CALL( nlhdlr->detect(scip, conshdlr, nlhdlr, expr, cons, enforcing, participating, nlhdlrexprdata) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->detecttime) );

   if( *participating != SCIP_NLHDLR_METHOD_NONE )
   {
      ++nlhdlr->ndetections;
      ++nlhdlr->ndetectionslast;
   }

   return SCIP_OKAY;
}

/** call the auxiliary evaluation callback of a nonlinear handler */
SCIP_DECL_NLHDLREVALAUX(SCIPnlhdlrEvalaux)
{
   assert(nlhdlr != NULL);
   assert(nlhdlr->evalaux != NULL);

   SCIP_CALL( nlhdlr->evalaux(scip, nlhdlr, expr, nlhdlrexprdata, auxvalue, sol) );

   return SCIP_OKAY;
}

/** calls the interval evaluation callback of a nonlinear handler */
SCIP_DECL_NLHDLRINTEVAL(SCIPnlhdlrInteval)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->intevaltime != NULL);

   if( nlhdlr->inteval != NULL )
   {
      SCIP_CALL( SCIPstartClock(scip, nlhdlr->intevaltime) );
      SCIP_CALL( nlhdlr->inteval(scip, nlhdlr, expr, nlhdlrexprdata, interval, intevalvar, intevalvardata) );
      SCIP_CALL( SCIPstopClock(scip, nlhdlr->intevaltime) );

      ++nlhdlr->nintevalcalls;
   }

   return SCIP_OKAY;
}

/** calls the reverse propagation callback of a nonlinear handler */
SCIP_DECL_NLHDLRREVERSEPROP(SCIPnlhdlrReverseprop)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->proptime != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

   if( nlhdlr->reverseprop == NULL )
   {
      *infeasible = FALSE;
      *nreductions = 0;

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->proptime) );
   SCIP_CALL( nlhdlr->reverseprop(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, bounds, infeasible, nreductions) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->proptime) );

   /* update statistics */
   nlhdlr->ndomreds += *nreductions;
   if( *infeasible )
      ++nlhdlr->ncutoffs;
   ++nlhdlr->npropcalls;

   return SCIP_OKAY;
}

/** calls the separation initialization callback of a nonlinear handler */
SCIP_DECL_NLHDLRINITSEPA(SCIPnlhdlrInitsepa)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);
   assert(infeasible != NULL);

   if( nlhdlr->initsepa == NULL )
   {
      *infeasible = FALSE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( nlhdlr->initsepa(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, overestimate, underestimate, infeasible) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );

   ++nlhdlr->nenfocalls;
   if( *infeasible )
      ++nlhdlr->ncutoffs;

   return SCIP_OKAY;
}

/** calls the separation deinitialization callback of a nonlinear handler */
SCIP_DECL_NLHDLREXITSEPA(SCIPnlhdlrExitsepa)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);

   if( nlhdlr->exitsepa != NULL )
   {
      SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
      SCIP_CALL( nlhdlr->exitsepa(scip, nlhdlr, expr, nlhdlrexprdata) );
      SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );
   }

   return SCIP_OKAY;
}

/** calls the enforcement callback of a nonlinear handler */
SCIP_DECL_NLHDLRENFO(SCIPnlhdlrEnfo)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);
   assert(result != NULL);

   if( nlhdlr->enfo == NULL )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

#ifndef NDEBUG
   /* check that auxvalue is correct by reevaluating */
   {
      SCIP_Real auxvaluetest;
      SCIP_CALL( SCIPnlhdlrEvalaux(scip, nlhdlr, expr, nlhdlrexprdata, &auxvaluetest, sol) );
      assert(auxvalue == auxvaluetest);  /* we should get EXACTLY the same value from calling evalaux with the same solution as before */  /*lint !e777*/
   }
#endif

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( nlhdlr->enfo(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate, allowweakcuts, separated, addbranchscores, result) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );

   /* update statistics */
   ++nlhdlr->nenfocalls;
   switch( *result )
   {
      case SCIP_SEPARATED :
         ++nlhdlr->nseparated;
         break;
      case SCIP_BRANCHED:
         ++nlhdlr->nbranchscores;
         break;
      case SCIP_CUTOFF:
         ++nlhdlr->ncutoffs;
         break;
      case SCIP_REDUCEDDOM:
         ++nlhdlr->ndomreds;
         break;
      default: ;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** calls the estimator callback of a nonlinear handler */
SCIP_DECL_NLHDLRESTIMATE(SCIPnlhdlrEstimate)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);
   assert(success != NULL);
   assert(addedbranchscores != NULL);

   if( nlhdlr->estimate == NULL )
   {
      *success = FALSE;
      *addedbranchscores = FALSE;
      return SCIP_OKAY;
   }

#ifndef NDEBUG
   /* check that auxvalue is correct by reevaluating */
   {
      SCIP_Real auxvaluetest;
      SCIP_CALL( SCIPnlhdlrEvalaux(scip, nlhdlr, expr, nlhdlrexprdata, &auxvaluetest, sol) );
      assert(auxvalue == auxvaluetest);  /* we should get EXACTLY the same value from calling evalaux with the same solution as before */  /*lint !e777*/
   }
#endif

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( nlhdlr->estimate(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate, targetvalue,
         rowpreps, success, addbranchscores, addedbranchscores) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );

   /* update statistics */
   ++nlhdlr->nenfocalls;

   return SCIP_OKAY;
}
