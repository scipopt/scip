/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_lj.h
 * @brief  Lennard-Jones potential nonlinear handler
 * @author Stefan Vigerske
 */

#include "nlhdlr_lj.h"
#include "scip/cons_nonlinear.h"
#include "scip/pub_misc_rowprep.h"


/* fundamental nonlinear handler properties */
#define NLHDLR_NAME               "lj"
#define NLHDLR_DESC               "handler of Lennard-Jones potential function"
#define NLHDLR_DETECTPRIORITY     0
#define NLHDLR_ENFOPRIORITY       0

/** minimum of r^(-6) - r^(-3); 1.26 */
#define RMINIMUM pow(2, 1.0/3.0)
/** inflection point of r^(-6) - r^(-3); 1.52 */
#define RINFLECT pow(3.5, 1.0/3.0)

/** maximal number of iterations for Newton algorithm */
#define NEWTONITER 20

/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_NlhdlrExprData
{
   SCIP_EXPR*          r;                    /**< r argument */
   SCIP_EXPR*          p;                    /**< p argument */
   SCIP_Real           rcoef;                /**< coefficient for r-term: 1 or -1 */
   SCIP_Real           pcoef;                /**< coefficient of p */
};

/*
 * Local methods
 */

/** Function evaluation callback for Newton algorithm.
 *
 * To find point such that f(rub) - f(point) = f'(point) * (rub - point), we use Newton's algorithm to solve
 * g(point) := f(rub) - f(point) - f'(point) * (rub - point).
 * This function evaluates g(point).
 */
static
SCIP_DECL_NEWTONEVAL(rootfunc)
{
   SCIP_Real rub;
   SCIP_Real val;

   assert(nparams == 1);
   assert(point > 1e-9);

   rub = params[0];

   val = pow(rub, -6.0) - pow(rub, -3.0) - pow(point, -6.0) + pow(point, -3.0);
   val -= (-6.0 * pow(point, -7.0) + 3.0 * pow(point, -4.0)) * (rub - point);

   return val;
}

/** Gradient evaluation callback for Newton algorithm.
 *
 * Evaluates g'(point).
 */
static
SCIP_DECL_NEWTONEVAL(rootfuncderiv)
{
   SCIP_Real rub;

   assert(nparams == 1);
   assert(point > 1e-9);

   rub = params[0];

   /* derivative of f(rub) - f(point) - f'(point) * (rub - point) w.r.t. point
    * -f'(point) - f''(point) (rub - point) + f'(point) = - f''(point) (rub - point)
    */
   return - (42.0 * pow(point, -8.0) - 12.0 * pow(point, -5.0)) * (rub - point);
}

/** Compute underestimator for r^(-6) - r^(-3). */
static
SCIP_Bool underestimator(
   SCIP*                 scip,
   SCIP_Real*            slope,
   SCIP_Real*            constant,
   SCIP_Real             rval,
   SCIP_Real             rlb,
   SCIP_Real             rub
)
{
   SCIP_Real r;

   /* too close or negative of zero */
   if( SCIPisNegative(scip, rval) )
      return FALSE;

   if( SCIPisLE(scip, rval, RMINIMUM) || SCIPisLE(scip, rub, RINFLECT) )
   {
      /* rub <= RINFLECT: function is convex on [rlb,rub]
       * rval <= RMINIMUM: function is convex at rval, and slope is negative
       * underestimator is f(r) + f'(rval) (r - rval)
       */
      *slope = -6.0 * pow(rval, -7.0) + 3.0 * pow(rval, -4.0);
      *constant = pow(rval, -6.0) - pow(rval, -3.0) - *slope * rval;

      assert(rval > RMINIMUM || SCIPisNegative(scip, *slope));

      return TRUE;
   }

   if( SCIPisHugeValue(scip, rub) )
      return FALSE;

   if( SCIPisGE(scip, rlb, RINFLECT) )
   {
      if( SCIPisRelEQ(scip, rlb, rub) )
         return FALSE;

      /* function is concave on [rlb, rub]
       * secant f(rlb) + (f(rub) - f(rlb)) / (rub - rlb) * (r - rlb) */
      *slope = (pow(rub, -6.0) - pow(rub, -3.0) - pow(rlb, -6.0) + pow(rlb, -3.0)) / (rub - rlb);
      *constant = pow(rlb, -6.0) - pow(rlb, -3.0) - *slope * rlb;

      return TRUE;
   }

   /* find r such that the secant from r to rub has slope f'(r) */
   r = SCIPcalcRootNewton(rootfunc, rootfuncderiv, &rub, 1, (RMINIMUM + RINFLECT) / 2.0, SCIPepsilon(scip), NEWTONITER);

   SCIPdebugMsg(scip, "newton returned r = %.15g for rub = %.15g\n", r, rub);

   if( r == SCIP_INVALID )
      return FALSE;

   if( rlb <= r )
   {
      /* if rlb <= r: underestimator is tangent at r (== secant from r to rub) */
      *slope = -6.0 * pow(r, -7.0) + 3.0 * pow(r, -4.0);
      *constant = pow(r, -6.0) - pow(r, -3.0) - *slope * r;
   }
   else
   {
      /* for rlb > r, secant between rlb and rub */
      *slope = (pow(rub, -6.0) - pow(rub, -3.0) - pow(rlb, -6.0) + pow(rlb, -3.0)) / (rub - rlb);
      *constant = pow(rlb, -6.0) - pow(rlb, -3.0) - *slope * rlb;
   }

   return TRUE;
}

/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler copy callback */
static
SCIP_DECL_NLHDLRCOPYHDLR(nlhdlrCopyhdlrLJ)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(sourcenlhdlr != NULL);

   SCIP_STRINGEQ( SCIPnlhdlrGetName(sourcenlhdlr), NLHDLR_NAME, SCIP_INVALIDCALL );

   SCIP_CALL( SCIPincludeNlhdlrLJ(targetscip) );

   return SCIP_OKAY;
}

/** callback to free expression specific data */
static
SCIP_DECL_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataLJ)
{  /*lint --e{715}*/
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** callback to detect structure in expression tree */
static
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectLJ)
{ /*lint --e{715}*/
   SCIP_EXPR* r = NULL;
   SCIP_EXPR* p = NULL;
   SCIP_Real pcoef = 0.0;
   SCIP_Bool six = FALSE;
   SCIP_Bool three = FALSE;
   SCIP_Bool dobelow = FALSE;
   SCIP_Bool doabove = FALSE;
   int i;

   /* if no under- or overestimators are needed, then do nothing */
   if( (*enforcing & SCIP_NLHDLR_METHOD_SEPABOTH) == SCIP_NLHDLR_METHOD_SEPABOTH )
      return SCIP_OKAY;

   /* look for an expression r^(-6) - r^(-3) - p */

   if( !SCIPisExprSum(scip, expr) )
      return SCIP_OKAY;

   if( SCIPexprGetNChildren(expr) != 3 )
      return SCIP_OKAY;

   for( i = 0; i < 3; ++i )
   {
      SCIP_EXPR* child;
      SCIP_Real coef;

      child = SCIPexprGetChildren(expr)[i];
      assert(child != NULL);

      coef = SCIPgetCoefsExprSum(expr)[i];

      if( SCIPisExprPower(scip, child) && SCIPgetExponentExprPow(child) == -6.0 && !six )
      {
         /* our first r^(-6) */
         if( coef == 1.0 && !doabove )
            dobelow = TRUE;
         else if( coef == -1.0 && !dobelow )
            doabove = TRUE;
         else /* coef not 1 or -1 or inconsistent sign */
            return SCIP_OKAY;

         if( r == NULL )
            r = SCIPexprGetChildren(child)[0];
         else if( r != SCIPexprGetChildren(child)[0] )
            return SCIP_OKAY;

         six = TRUE;
         continue;
      }

      if( SCIPisExprPower(scip, child) && SCIPgetExponentExprPow(child) == -3.0 && !three )
      {
         /* our first -r^(-3) */
         if( coef == -1.0 && !doabove )
            dobelow = TRUE;
         else if( coef == 1.0 && !dobelow )
            doabove = TRUE;
         else /* coef not 1 or -1 or inconsistent sign */
            return SCIP_OKAY;

         if( r == NULL )
            r = SCIPexprGetChildren(child)[0];
         else if( r != SCIPexprGetChildren(child)[0] )
            return SCIP_OKAY;

         three = TRUE;
         continue;
      }

      if( p == NULL )
      {
         p = child;
         pcoef = coef;
         continue;
      }

      return SCIP_OKAY;
   }

   if( !three || !six || p == NULL )
      return SCIP_OKAY;
   assert(r != NULL);

   SCIPdebugMsg(scip, "detected LJ Potential in expr ");
#ifdef SCIP_DEBUG
   SCIPprintExpr(scip, expr, NULL);
   SCIPinfoMessage(scip, NULL, " -> sepa %s\n", dobelow ? "below" : "above");
#endif

   assert(dobelow != doabove);
   if( dobelow )
   {
      /* if someone already provides good underestimators, then skip */
      if( *enforcing & SCIP_NLHDLR_METHOD_SEPABELOW )
         return SCIP_OKAY;
      *participating = SCIP_NLHDLR_METHOD_SEPABELOW;
   }
   else
   {
      /* if someone already provides good overestimators, then skip */
      if( *enforcing & SCIP_NLHDLR_METHOD_SEPAABOVE )
         return SCIP_OKAY;
      *participating = SCIP_NLHDLR_METHOD_SEPAABOVE;
   }

   /* our method is enforcing */
   *enforcing |= *participating;

   /* we need an auxiliary variable for r and will use its activity for the under- or overestimator */
   SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, r, TRUE, FALSE, dobelow, doabove) );
   /* we need an (auxiliary) variable for p */
   SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, p, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPallocBlockMemory(scip, nlhdlrexprdata) );
   (*nlhdlrexprdata)->r = r;
   (*nlhdlrexprdata)->p = p;
   (*nlhdlrexprdata)->rcoef = dobelow ? 1.0 : -1.0;
   (*nlhdlrexprdata)->pcoef = pcoef;

   return SCIP_OKAY;
}

/** auxiliary evaluation callback of nonlinear handler */
static
SCIP_DECL_NLHDLREVALAUX(nlhdlrEvalauxLJ)
{ /*lint --e{715}*/
   SCIP_VAR* raux;
   SCIP_VAR* paux;
   SCIP_Real rval;
   SCIP_Real pval;

   assert(nlhdlrexprdata != NULL);

   raux = SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->r);
   paux = SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->p);
   assert(raux != NULL);
   assert(paux != NULL);

   /* get solution values of the auxiliary variables */
   rval = SCIPgetSolVal(scip, sol, raux);
   pval = SCIPgetSolVal(scip, sol, paux);

   /* evaluate expression w.r.t. the values of the auxiliary variables
    * return SCIP_INVALID if r is zero
    */
   if( rval == 0.0 )
      *auxvalue = SCIP_INVALID;
   else
      *auxvalue = nlhdlrexprdata->rcoef * (pow(rval, -6.0) - pow(rval, -3.0)) + nlhdlrexprdata->pcoef * pval;

   return SCIP_OKAY;
}

/** separation initialization method of a nonlinear handler (called during CONSINITLP) */
#if 0
static
SCIP_DECL_NLHDLRINITSEPA(nlhdlrInitSepaLJ)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of lj nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrInitSepaLJ NULL
#endif

/** nonlinear handler under/overestimation callback */
static
SCIP_DECL_NLHDLRESTIMATE(nlhdlrEstimateLJ)
{ /*lint --e{715}*/
   SCIP_ROWPREP* rowprep;
   SCIP_VAR* rvar;
   SCIP_Real rval;
   SCIP_Real rlb;
   SCIP_Real rub;
   SCIP_Real slope;
   SCIP_Real constant;

   *success = FALSE;

   if( overestimate != (nlhdlrexprdata->rcoef == -1.0) )
      return SCIP_OKAY;

   rvar = SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->r);
   assert(rvar != NULL);

   rval = SCIPgetSolVal(scip, sol, rvar);

   rlb = SCIPvarGetLbLocal(rvar);
   rub = SCIPvarGetUbLocal(rvar);

   SCIPdebugMsg(scip, "estimate at r=%g in [%g,%g]\n", rval, rlb, rub);

   *success = underestimator(scip, &slope, &constant, rval, rlb, rub);

   if( !*success )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "-> %g * r %+g\n", slope, constant);

   /* cuts are globally valid if generated left at RMININUM, or if the function is convex w.r.t. global bounds */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, (rval > RMINIMUM) && (SCIPvarGetUbGlobal(rvar) > RINFLECT)) );
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, rvar, slope * nlhdlrexprdata->rcoef) );
   SCIProwprepAddConstant(rowprep, constant * nlhdlrexprdata->rcoef);
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->p), nlhdlrexprdata->pcoef) );

   SCIP_CALL( SCIPsetPtrarrayVal(scip, rowpreps, 0, rowprep) );

   /* add branching scores if requested */
   if( addbranchscores && rval > RMINIMUM )
   {
      SCIP_Real violation;

      /* compute violation w.r.t. the auxiliary variable(s) */
#ifndef BRSCORE_ABSVIOL
      SCIP_CALL( SCIPgetExprRelAuxViolationNonlinear(scip, expr, auxvalue, sol, &violation, NULL, NULL) );
#else
      SCIP_CALL( SCIPgetExprAbsAuxViolationNonlinear(scip, expr, auxvalue, sol, &violation, NULL, NULL) );
#endif
      assert(violation > 0.0);  /* there should be a violation if we were called to enforce */

      SCIP_CALL( SCIPaddExprsViolScoreNonlinear(scip, &nlhdlrexprdata->r, 1, violation, sol, addedbranchscores) );
   }

   return SCIP_OKAY;
}

/** nonlinear handler solution linearization callback */
#if !1
static
SCIP_DECL_NLHDLRSOLLINEARIZE(nlhdlrSollinearizeLJ)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of lj nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrSollinearizeLJ NULL
#endif


/** nonlinear handler interval evaluation callback */
#if 0
static
SCIP_DECL_NLHDLRINTEVAL(nlhdlrIntevalLJ)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of lj nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrIntevalLJ NULL
#endif


/** nonlinear handler callback for reverse propagation */
#if 0
static
SCIP_DECL_NLHDLRREVERSEPROP(nlhdlrReversepropLJ)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of lj nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrReversepropLJ NULL
#endif


/*
 * nonlinear handler specific interface methods
 */

/** includes LJ nonlinear handler in nonlinear constraint handler */
SCIP_RETCODE SCIPincludeNlhdlrLJ(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLHDLR* nlhdlr;

   assert(scip != NULL);

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY,
      NLHDLR_ENFOPRIORITY, nlhdlrDetectLJ, nlhdlrEvalauxLJ, NULL) );
   assert(nlhdlr != NULL);

   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrLJ);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrFreeExprDataLJ);
   SCIPnlhdlrSetSepa(nlhdlr, nlhdlrInitSepaLJ, NULL, nlhdlrEstimateLJ, NULL);
   SCIPnlhdlrSetProp(nlhdlr, nlhdlrIntevalLJ, nlhdlrReversepropLJ);

   return SCIP_OKAY;
}
