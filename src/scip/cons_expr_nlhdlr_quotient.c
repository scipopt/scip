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

/**@file   cons_expr_nlhdlr_quotient.h
 * @brief  quotient nonlinear handler
 * @author Benjamin Mueller
 * @author Fabian Wegscheider
 *
 * @todo implement INITSEPA
 */

#include <string.h>

#include "scip/cons_expr_nlhdlr_quotient.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME         "quotient"
#define NLHDLR_DESC         "quotient handler for quotient expressions"
#define NLHDLR_PRIORITY     0

/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_VAR*             nomvar;             /**< variable of the nominator */
   SCIP_Real             nomcoef;            /**< coefficient of the nominator */
   SCIP_Real             nomconst;           /**< constant of the nominator */
   SCIP_VAR*             denomvar;           /**< variable of the denominator */
   SCIP_Real             denomcoef;          /**< coefficient of the denominator */
   SCIP_Real             denomconst;         /**< constant of the denominator */
   SCIP_Real             constant;           /**< constant */
};

/*
 * Local methods
 */

/** helper method to create nonlinear handler expression data */
static
SCIP_RETCODE exprdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata, /**< nonlinear handler expression data */
   SCIP_VAR*             nomvar,             /**< variable of the nominator */
   SCIP_Real             nomcoef,            /**< coefficient of the nominator */
   SCIP_Real             nomconst,           /**< constant of the nominator */
   SCIP_VAR*             denomvar,           /**< variable of the denominator */
   SCIP_Real             denomcoef,          /**< coefficient of the denominator */
   SCIP_Real             denomconst,         /**< constant of the denominator */
   SCIP_Real             constant            /**< constant */
   )
{
   assert(nlhdlrexprdata != NULL);
   assert(nomvar != NULL);
   assert(denomvar != NULL);
   assert(!SCIPisZero(scip, nomcoef));
   assert(!SCIPisZero(scip, denomcoef));

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, nlhdlrexprdata) );

   /* store values */
   (*nlhdlrexprdata)->nomvar = nomvar;
   (*nlhdlrexprdata)->nomcoef = nomcoef;
   (*nlhdlrexprdata)->nomconst = nomconst;
   (*nlhdlrexprdata)->denomvar = denomvar;
   (*nlhdlrexprdata)->denomcoef = denomcoef;
   (*nlhdlrexprdata)->denomconst = denomconst;
   (*nlhdlrexprdata)->constant = constant;

   /* capture variables */
   SCIP_CALL( SCIPcaptureVar(scip, nomvar) );
   SCIP_CALL( SCIPcaptureVar(scip, denomvar) );

   return SCIP_OKAY;
}

/** helper method to free nonlinear handler expression data */
static
SCIP_RETCODE exprdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata  /**< nonlinear handler expression data */
   )
{
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);
   assert((*nlhdlrexprdata)->nomvar != NULL);
   assert((*nlhdlrexprdata)->denomvar != NULL);

   /* release variables */
   SCIP_CALL( SCIPreleaseVar(scip, &(*nlhdlrexprdata)->denomvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &(*nlhdlrexprdata)->nomvar) );

   /* free expression data of nonlinear handler */
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** helper method to detect whether an expression is of the form a*x + b */
static
SCIP_Bool isExprUnivariateLinear(
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef,               /**< pointer to store the coefficient */
   SCIP_Real*            constant            /**< pointer to store the constant */
   )
{
   assert(expr != NULL);
   assert(conshdlr != NULL);
   assert(coef != NULL);
   assert(constant != NULL);

   *var = NULL;
   *coef = 0.0;
   *constant = 0.0;

   /* expression is a variable, i.e., a = 1, b = 0 */
   if( SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrVar(conshdlr) )
   {
      *var = SCIPgetConsExprExprVarVar(expr);
      *coef = 1.0;
      *constant = 0.0;
      return TRUE;
   }
   /* expression is a sum; check whether it consists only of one variable expression */
   else if( SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr) && SCIPgetConsExprExprNChildren(expr) == 1 )
   {
      SCIP_CONSEXPR_EXPR* child = SCIPgetConsExprExprChildren(expr)[0];
      assert(child != NULL);

      /* child must be a variable */
      if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrVar(conshdlr) )
      {
         *var = SCIPgetConsExprExprVarVar(child);
         *coef = SCIPgetConsExprExprSumCoefs(expr)[0];
         *constant = SCIPgetConsExprExprSumConstant(expr);
         return TRUE;
      }
   }

   return FALSE;
}

/** helper method to detect an expression of the form (a*x + b) / (c*y + d) + e; due to the expansion of products,
  * there are two types of expressions that can be detected:
  *
  * 1. prod(f(x), pow(g(y),-1))
  * 2. sum(prod(f(x),pow(g(y),-1)), pow(g(y),-1))
  *
  * @todo At the moment quotients like xy / z are not detected, because they are turned into a product expression
  * with three children, i,e., x * y * (1 / z).
  */
static
SCIP_RETCODE detectExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata, /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            success             /**< pointer to store whether nonlinear handler should be called for this expression */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* prodhdlr;
   SCIP_CONSEXPR_EXPRHDLR* sumhdlr;
   SCIP_CONSEXPR_EXPRHDLR* powhdlr;
   SCIP_CONSEXPR_EXPR** children;
   SCIP_CONSEXPR_EXPR* denomexpr = NULL;
   SCIP_CONSEXPR_EXPR* nomexpr = NULL;
   SCIP_VAR* x = NULL;
   SCIP_VAR* y = NULL;
   SCIP_Real a, b, c, d, e;
   SCIP_Real nomfac = 1.0;
   SCIP_Real nomconst = 0.0;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);

   *success = FALSE;
   a = 0.0;
   b = 0.0;
   c = 0.0;
   d = 0.0;
   e = 0.0;

   /* possible structures only have two children */
   if( SCIPgetConsExprExprNChildren(expr) != 2 )
      return SCIP_OKAY;

   /* collect expression handlers */
   prodhdlr = SCIPgetConsExprExprHdlrProduct(conshdlr);
   sumhdlr = SCIPgetConsExprExprHdlrSum(conshdlr);
   powhdlr = SCIPgetConsExprExprHdlrPower(conshdlr);

   /* expression must be either a product or a sum */
   if( SCIPgetConsExprExprHdlr(expr) != prodhdlr && SCIPgetConsExprExprHdlr(expr) != sumhdlr )
      return SCIP_OKAY;

   children = SCIPgetConsExprExprChildren(expr);
   assert(children != NULL);

   /* case: prod(f(x), pow(g(y),-1)) */
   if( SCIPgetConsExprExprHdlr(expr) == prodhdlr )
   {
      if( SCIPgetConsExprExprHdlr(children[0]) == powhdlr && SCIPgetConsExprExprPowExponent(children[0]) == -1 )
      {
         denomexpr = SCIPgetConsExprExprChildren(children[0])[0];
         nomexpr = children[1];
      }
      else if( SCIPgetConsExprExprHdlr(children[1]) == powhdlr && SCIPgetConsExprExprPowExponent(children[1]) == -1 )
      {
         denomexpr = SCIPgetConsExprExprChildren(children[1])[0];
         nomexpr = children[0];
      }

      /* remember to scale the nominator by the coefficient stored in the product expression */
      nomfac = SCIPgetConsExprExprProductCoef(expr);
   }
   /* case: sum(prod(f(x),pow(g(y),-1)), pow(g(y),-1)) */
   else
   {
      SCIP_Real* sumcoefs;

      assert(SCIPgetConsExprExprHdlr(expr) == sumhdlr);
      sumcoefs = SCIPgetConsExprExprSumCoefs(expr);

      /* children[0] is 1/g(y) and children[1] is a product of f(x) and 1/g(y) */
      if( SCIPgetConsExprExprHdlr(children[0]) == powhdlr && SCIPgetConsExprExprPowExponent(children[0]) == -1
         && SCIPgetConsExprExprHdlr(children[1]) == prodhdlr && SCIPgetConsExprExprNChildren(children[1]) == 2 )
      {
         SCIP_Real prodcoef = SCIPgetConsExprExprProductCoef(children[1]);

         if( children[0] == SCIPgetConsExprExprChildren(children[1])[0] )
         {
            denomexpr = SCIPgetConsExprExprChildren(children[0])[0];
            nomexpr = SCIPgetConsExprExprChildren(children[1])[1];
         }
         else if( children[0] == SCIPgetConsExprExprChildren(children[1])[1] )
         {
            denomexpr = SCIPgetConsExprExprChildren(children[0])[0];
            nomexpr = SCIPgetConsExprExprChildren(children[1])[0];
         }

         /* remember scalar and constant for nominator */
         nomfac = sumcoefs[1] * prodcoef;
         nomconst = sumcoefs[0];
      }
      /* children[1] is 1/g(y) and children[0] is a product of f(x) and 1/g(y) */
      else if( SCIPgetConsExprExprHdlr(children[1]) == powhdlr && SCIPgetConsExprExprPowExponent(children[1]) == -1
         && SCIPgetConsExprExprHdlr(children[0]) == prodhdlr && SCIPgetConsExprExprNChildren(children[0]) == 2 )
      {
         SCIP_Real prodcoef = SCIPgetConsExprExprProductCoef(children[0]);

         if( children[1] == SCIPgetConsExprExprChildren(children[0])[0] )
         {
            denomexpr = SCIPgetConsExprExprChildren(children[1])[0];
            nomexpr = SCIPgetConsExprExprChildren(children[0])[1];
         }
         else if( children[1] == SCIPgetConsExprExprChildren(children[0])[1] )
         {
            denomexpr = SCIPgetConsExprExprChildren(children[1])[0];
            nomexpr = SCIPgetConsExprExprChildren(children[0])[0];
         }

         /* remember scalar and constant for nominator */
         nomfac = sumcoefs[0] * prodcoef;
         nomconst = sumcoefs[1];
      }

      /* remember the constant of the sum expression */
      e = SCIPgetConsExprExprSumConstant(expr);
   }

   if( denomexpr != NULL && nomexpr != NULL )
   {
      /* nominator and denominator are univariate linear functions -> no auxiliary variables are needed */
      if( isExprUnivariateLinear(nomexpr, conshdlr, &x, &a, &b)
         && isExprUnivariateLinear(denomexpr, conshdlr, &y, &c, &d) )
      {
         SCIPdebugMsg(scip, "detected nominator (%g * %s + %g) and denominator (%g * %s + %g) to be univariate and linear\n",
            a, SCIPvarGetName(x), b, c, SCIPvarGetName(y), d);

         /* during presolving, it only makes sense to detect the quotient if both variables are the same */
         *success = (SCIPgetStage(scip) == SCIP_STAGE_SOLVING) || (x == y);

         /* if the variables are different and it is not of the form x / y, add auxiliary variables */
         if( *success && x != y && (a != 0.0 || b != 0.0 || c != 0.0 || d != 0.0) )
         {
            SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, nomexpr, &x) );
            a = 1.0;
            b = 0.0;

            SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, denomexpr, &y) );
            c = 1.0;
            d = 0.0;
         }
      }
      /* create auxiliary variables if we are in the solving stage */
      else if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      {
         assert(x == NULL);
         assert(y == NULL);

         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, nomexpr, &x) );
         a = 1.0;
         b = 0.0;

#ifdef SCIP_DEBUG
         SCIPinfoMessage(scip, NULL, "Expression for nominator: ");
         SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, nomexpr, NULL) );
         SCIPinfoMessage(scip, NULL, " is not univariate and linear -> add auxiliary variable %s\n", SCIPvarGetName(x));
#endif

         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, denomexpr, &y) );
         c = 1.0;
         d = 0.0;

#ifdef SCIP_DEBUG
         SCIPinfoMessage(scip, NULL, "Expression for denominator: ");
         SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, denomexpr, NULL) );
         SCIPinfoMessage(scip, NULL, " is not univariate and linear -> add auxiliary variable %s\n", SCIPvarGetName(y));
#endif
         *success = TRUE;
      }
   }

   /* create nonlinear handler expression data */
   if( *success )
   {
      assert(x != NULL);
      assert(y != NULL);
      assert(a != 0.0);
      assert(c != 0.0);

      a = nomfac * a;
      b = nomfac * b + nomconst;

      SCIPdebug( SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, expr, NULL) ); )
      SCIPdebug( SCIPinfoMessage(scip, NULL, "\n") );
      SCIPdebugMsg(scip, "detected quotient expression (%g * %s + %g) / (%g * %s + %g) + %g\n", a, SCIPvarGetName(x),
         b, c, SCIPvarGetName(y), d, e);
      SCIP_CALL( exprdataCreate(scip, nlhdlrexprdata, x, a, b, y, c, d, e) );
   }

   return SCIP_OKAY;
}

/** helper method to compute interval for (a x + b) / (c x + d) + e */
static
SCIP_INTERVAL intEvalQuotient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTERVAL         bnds,               /**< bounds on x */
   SCIP_Real             a,                  /**< coefficient in nominator */
   SCIP_Real             b,                  /**< constant in nominator */
   SCIP_Real             c,                  /**< coefficient in denominator */
   SCIP_Real             d,                  /**< constant in denominator */
   SCIP_Real             e                   /**< constant */
   )
{
   SCIP_INTERVAL result;
   SCIP_INTERVAL denominterval;
   SCIP_Real infeval;
   SCIP_Real supeval;

   assert(scip != NULL);

   /* return empty interval if the domain of x is empty */
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, bnds) )
   {
      SCIPintervalSetEmpty(&result);
      return result;
   }

   /* compute bounds for denominator */
   SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &denominterval, bnds, c);
   SCIPintervalAddScalar(SCIP_INTERVAL_INFINITY, &denominterval, denominterval, d);

   /* there is no useful interval if 0 is in the interior of the interval of the denominator */
   if( SCIPintervalGetInf(denominterval) < 0.0 && SCIPintervalGetSup(denominterval) > 0.0 )
   {
      SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &result);
      return result;
   }

   assert(!SCIPisZero(scip, c));

   if( SCIPisInfinity(scip, -bnds.inf) )
      infeval = a / c;
   else
      infeval = (a * bnds.inf + b) / (c * bnds.inf + d) + e;

   if( SCIPisInfinity(scip, bnds.sup) )
      supeval = a / c;
   else
      supeval = (a * bnds.sup + b) / (c * bnds.sup + d) + e;

   /* f(x) = (a x + b) / (c x + d) + e implies f'(x) = (a d - b c) / (d + c x)^2 */
   if( a*d - b*c > 0.0 ) /* monotone increasing */
   {
      SCIPintervalSetBounds(&result, infeval, supeval);
   }
   else if( a*d - b*c < 0.0 ) /* monotone decreasing */
   {
      SCIPintervalSetBounds(&result, supeval, infeval);
   }
   else /* a d = b c implies that f(x) = b / d + e, i.e., f is constant */
   {
      assert(a*d - b*c == 0.0);

      SCIPintervalSet(&result, b / d + e);
   }

   return result;
}

/** helper method to compute reverse propagation for (a x + b) / (c x + d) + e */
static
SCIP_INTERVAL reversepropQuotient(
   SCIP_INTERVAL         bnds,               /**< bounds on (a x + b) / (c x + d) + e */
   SCIP_Real             a,                  /**< coefficient in nominator */
   SCIP_Real             b,                  /**< constant in nominator */
   SCIP_Real             c,                  /**< coefficient in denominator */
   SCIP_Real             d,                  /**< constant in denominator */
   SCIP_Real             e                   /**< constant */
   )
{
   SCIP_INTERVAL result;
   SCIP_Real infpropval;
   SCIP_Real suppropval;

   /* return empty interval if the domain of the expression is empty */
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, bnds) )
   {
      SCIPintervalSetEmpty(&result);
      return result;
   }

   /* substract constrant from bounds of the expression */
   SCIPintervalSubScalar(SCIP_INTERVAL_INFINITY, &bnds, bnds, e);

   /* if the expression is constant or the limit lies inside the domain, nothing can be propagated */
   if( a*d - b*c == 0.0 || (bnds.inf < a / c && bnds.sup > a / c) )
   {
      SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &result);
      return result;
   }

   infpropval = (d * bnds.inf - b) / (a - c * bnds.inf);
   suppropval = (d * bnds.sup - b) / (a - c * bnds.sup);

   /* f(x) = (a x + b) / (c x + d) + e implies f'(x) = (a d - b c) / (d + c x)^2 */
   if( a*d - b*c > 0.0 ) /* monotone increasing */
   {
      assert(infpropval <= suppropval);
      SCIPintervalSetBounds(&result, infpropval, suppropval);
   }
   else /* monotone decreasing */
   {
      assert(a*d - b*c < 0.0);
      assert(suppropval <= infpropval);
      SCIPintervalSetBounds(&result, suppropval, infpropval);
   }

   return result;
}

/** prepares a rowprep from given data; the generated estimator is always locally valid
 *
 *  @note the constant is moved to the left- or right-hand side
 */
static
SCIP_RETCODE createRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP *        rowprep,            /**< a rowprep where to store the estimator */
   const char*           name,               /**< name of type of cut */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_Real*            coefs,              /**< coefficients */
   SCIP_Real             constant,           /**< constant */
   int                   nlinvars            /**< total number of variables */
   )
{
   int i;

   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(coefs != NULL);
   assert(vars != NULL);

   /* create rowprep */
   (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, name);
   SCIPaddRowprepSide(rowprep, -constant);
   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, nlinvars + 1) );

   /* add coefficients */
   for( i = 0; i < nlinvars; ++i )
   {
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, vars[i], coefs[i]) );
   }

   return SCIP_OKAY;
}

/** separates a given point in the univariate case (ax + b) / (cx + d) + e
 *
 *  Depending on the reference point, the estimator is a tangent or a secant on the graph.
 *  It depends on whether we are under- or overestimating, whether we are on the left or
 *  on the right side of the singularity at -d/c, and whether it is the monotone increasing
 *  (ab - cd > 0) or decreasing part (ab - cd < 0). Together, there are 8 cases:
 *
 *  mon. incr. + overestimate + left hand side  -->  secant
 *  mon. incr. + overestimate + right hand side -->  tangent
 *  mon. incr. + understimate + left hand side  -->  tangent
 *  mon. incr. + understimate + right hand side -->  secant
 *  mon. decr. + overestimate + left hand side  -->  tangent
 *  mon. decr. + overestimate + right hand side -->  secant
 *  mon. decr. + understimate + left hand side  -->  secant
 *  mon. decr. + understimate + right hand side -->  tangent
 *
 * @TODO: check whether the formula is correct
 */
static
SCIP_RETCODE sepaUnivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution point (or NULL for the LP solution) */
   SCIP_VAR*             x,                  /**< argument variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_Real             a,                  /**< coefficient in nominator */
   SCIP_Real             b,                  /**< constant in nominator */
   SCIP_Real             c,                  /**< coefficient in denominator */
   SCIP_Real             d,                  /**< constant in denominator */
   SCIP_Real             e,                  /**< constant */
   SCIP_Bool             overestimate,       /**< whether the expression should be overestimated */
   SCIP_ROWPREP*         rowprep,            /**< a rowprep where to store the estimator */
   SCIP_Bool*            success             /**< buffer to store whether separation was successful */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_INTERVAL bnds;
   SCIP_Real singularity;
   SCIP_Real linconst;
   SCIP_Real lincoef;
   SCIP_Bool monincreasing;
   SCIP_Bool isinleftpart;

   assert(scip != NULL);
   assert(success != NULL);
   assert(rowprep != NULL);
   assert(x != NULL);
   assert(c != 0.0);

   *success = FALSE;

   bnds.inf = SCIPvarGetLbLocal(x);
   bnds.sup = SCIPvarGetUbLocal(x);
   singularity = -d / c;

   /* if 0 is in the denom interval, estimation is not possible */
   if( SCIPisLE(scip, bnds.inf, singularity) && SCIPisGE(scip, bnds.sup, singularity) )
      return SCIP_OKAY;

   isinleftpart = (bnds.sup < singularity);
   monincreasing = (a * b - c * d > 0.0);

   /* this encodes the 8 cases explained above */
   if( monincreasing == (overestimate == isinleftpart) )
   {
      SCIP_Real lbeval;
      SCIP_Real ubeval;

      /* if one of the bounds is infinite, secant cannot be computed */
      if( SCIPisInfinity(scip, -bnds.inf) || SCIPisInfinity(scip, bnds.sup) )
         return SCIP_OKAY;

      lbeval = (a * bnds.inf + b) / (c * bnds.inf + d) + e;
      ubeval = (a * bnds.sup + b) / (c * bnds.sup + d) + e;

      /* compute coefficient and constant of linear estimator */
      lincoef = (ubeval - lbeval) / (bnds.sup - bnds.inf);
      linconst = ubeval - lincoef * bnds.sup;
   }
   else
   {
      SCIP_Real solvarval;
      SCIP_Real soleval;

      solvarval = SCIPgetSolVal(scip, sol, x);
      soleval = (a * solvarval + b) / (c * solvarval + d) + e;

      /* compute coefficient and constant of linear estimator */
      lincoef = (a * d - b * c) / SQR(d + c * solvarval);
      linconst = soleval - lincoef * solvarval;
   }

   /* avoid huge values in the cut */
   if( SCIPisHugeValue(scip, ABS(lincoef)) || SCIPisHugeValue(scip, ABS(linconst)) )
      return SCIP_OKAY;

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "quot_%s_%lld", SCIPvarGetName(x), SCIPgetNLPs(scip));

   /* create rowprep */
   SCIP_CALL( createRowprep(scip, rowprep, name, &x, &lincoef, linconst, 1) );

   /* clean up rowprep */
   SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, sol, SCIP_CONSEXPR_CUTMAXRANGE, SCIPinfinity(scip), success) );

   return SCIP_OKAY;
}

/** helper method to compute a gradient cut for h^c(x,y) at a given reference point */
static
void hcGradCut(
   SCIP_Real             lbx,                /**< lower bound of x */
   SCIP_Real             ubx,                /**< upper bound of x */
   SCIP_Real             solx,               /**< solution value of x */
   SCIP_Real             soly,               /**< solution value of y */
   SCIP_Real*            coefx,              /**< pointer to store the coefficient of x */
   SCIP_Real*            coefy,              /**< pointer to store the coefficient of y */
   SCIP_Real*            constant            /**< pointer to store the constant */
   )
{
   SCIP_Real tmp1;
   SCIP_Real tmp2;

   assert(lbx >= 0.0);
   assert(lbx <= ubx);
   assert(soly > 0.0);
   assert(coefx != NULL);
   assert(coefy != NULL);
   assert(constant != NULL);

   tmp1 = SQRT(lbx * ubx) + solx;
   tmp2 = SQR(SQRT(lbx) + SQRT(ubx)) * soly;
   assert(tmp2 > 0.0);

   *coefx = 2.0 * tmp1 / tmp2;
   *coefy = -SQR(tmp1) / (tmp2 * soly);
   *constant = 2.0 * SQRT(lbx * ubx) * tmp1 / tmp2;
}

/** separates a given point in the bivariate case x/y <=/>= z (auxvar)
 *
 *  There are the following cases for y > 0:
 *
 *    1. lbx > 0 or ubx < 0
 *       a) underestimation: use z >= x/y >= (1/y) * ( (x + sqrt(lbx * ubx)) / (sqrt(lbx) + sqrt(ubx)) )
 *                           and build gradient cut
 *       b) overestimation:  use z <= 1/(lby*uby) * min{uby*x - lbx*y + lbx*lby, lby*x - ubx*y + ubx*uby}
 *
 *    2. lbx <= 0 and ubx >= 0
 *          use McCormick for x <=/>= y * z and transform resulting linear inequality
 *          x <=/>= ay + bz + c to z >=/<= (x - ay - c) / b. The direction of the inequality
 *          is preserved, since we assume y > 0.
 *
 *    If y < 0, swap and negate its bounds and compute the respective opposite estimator (and negate it).
 *    If 0 is in the interval of y, nothing is possible.
 */
static
SCIP_RETCODE sepaBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution point (or NULL for the LP solution) */
   SCIP_VAR*             x,                  /**< nominator variable */
   SCIP_VAR*             y,                  /**< denominator variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_Bool             overestimate,       /**< whether the expression should be overestimated */
   SCIP_ROWPREP*         rowprep,            /**< a rowprep where to store the estimator */
   SCIP_Bool*            success             /**< buffer to store whether separation was successful */
   )
{
   SCIP_VAR* linvars[2];
   SCIP_Real lincoefs[2];
   char name[SCIP_MAXSTRLEN];
   SCIP_Real lbx;
   SCIP_Real ubx;
   SCIP_Real lby;
   SCIP_Real uby;
   SCIP_Real solx;
   SCIP_Real soly;
   SCIP_Real linconst;
   SCIP_Bool xisnonnegative;
   SCIP_Bool yispositive;

   assert(scip != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   *success = FALSE;

   lbx = SCIPvarGetLbLocal(x);
   ubx = SCIPvarGetUbLocal(x);
   lby = SCIPvarGetLbLocal(y);
   uby = SCIPvarGetUbLocal(y);

   /* if 0 is in the interior of [lby,uby], no estimator is possible */
   if( SCIPisLT(scip, lby, 0.0) && SCIPisGT(scip, uby, 0.0) )
      return SCIP_OKAY;

   solx = SCIPgetSolVal(scip, sol, x);
   soly = SCIPgetSolVal(scip, sol, y);

   yispositive = SCIPisGT(scip, lby, 0.0);

   /* if y is not positive, swap and negate its bounds */
   if( !yispositive )
   {
      SCIP_Real tmp;

      tmp = uby;
      uby = -lby;
      lby = -tmp;
      soly = -soly;
   }

   /* case 1: 0 is not in the interior of [lbx,ubx] */
   if( SCIPisGE(scip, lbx, 0.0) || SCIPisLE(scip, ubx, 0.0) )
   {
      xisnonnegative = SCIPisGE(scip, lbx, 0.0);

      /* if x is not non-negative, swap and negate its bounds */
      if( !xisnonnegative )
      {
         SCIP_Real tmp;

         tmp = ubx;
         ubx = -lbx;
         lbx = -tmp;
         solx = -solx;
      }

      assert(SCIPisGE(scip, lbx, 0.0));
      assert(SCIPisGT(scip, lby, 0.0));
      assert(SCIPisLE(scip, lbx, solx) && SCIPisLE(scip, solx, ubx));
      assert(SCIPisLE(scip, lby, soly) && SCIPisLE(scip, soly, uby));

      /* case 1a: underestimating the original or overestimating the negated expression */
      if( overestimate != (xisnonnegative == yispositive) )
      {
         (void) hcGradCut(lbx, ubx, solx, soly, &lincoefs[0], &lincoefs[1], &linconst);
      }
      /* case 1b: overestimating the original or underestimating the negated expression */
      else
      {
         if( uby * solx - lbx * soly + lbx * lby <= lby * solx - ubx * soly + ubx * uby )
         {
            lincoefs[0] = 1.0 / lby;
            lincoefs[1] = -lbx / (lby * uby);
            linconst = lbx / uby;
         }
         else
         {
            lincoefs[0] = 1.0 / uby;
            lincoefs[1] = -ubx / (lby * uby);
            linconst = ubx / lby;
         }
      }

      /* avoid huge values in the estimator */
      if( SCIPisHugeValue(scip, ABS(lincoefs[0])) || SCIPisHugeValue(scip, ABS(lincoefs[1]))
         || SCIPisHugeValue(scip, ABS(linconst)) )
      {
         return SCIP_OKAY;
      }

      /* negate coefficients when considering the negated expression */
      if( xisnonnegative != yispositive )
      {
         lincoefs[0] = -lincoefs[0];
         lincoefs[1] = -lincoefs[1];
         linconst = -linconst;
      }
   }
   /* case 2: 0 is in the interior of [lbx,ubx] */
   else
   {
      SCIP_Real mccoefy = 0.0;
      SCIP_Real mccoefaux = 0.0;

      linconst = 0.0;

      SCIPaddBilinMcCormick(scip, 1.0, SCIPvarGetLbLocal(auxvar), SCIPvarGetUbLocal(auxvar),
         SCIPgetSolVal(scip, sol, auxvar), lby, uby, soly, overestimate,
         &mccoefaux, &mccoefy, &linconst, success);

      /* the McCormick coefficients of auxvar is always lby or uby, so it has to be > 0 */
      assert(SCIPisGT(scip, mccoefaux, 0.0));

      if( !(*success) )
         return SCIP_OKAY;

      lincoefs[0] = 1.0 / mccoefaux;
      lincoefs[1] = -mccoefy / mccoefaux;
      linconst = -linconst / mccoefaux;
   }

   linvars[0] = x;
   linvars[1] = y;

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "quot_%s_%s_%lld", SCIPvarGetName(x), SCIPvarGetName(y),
      SCIPgetNLPs(scip));

   /* create rowprep */
   SCIP_CALL( createRowprep(scip, rowprep, name, linvars, lincoefs, linconst, 2) );

   /* clean up rowprep */
   SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, sol, SCIP_CONSEXPR_CUTMAXRANGE, SCIPinfinity(scip), success) );

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler copy callback */
static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrCopyhdlrQuotient)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPgetConsExprNlhdlrName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConsExprNlhdlrQuotient(targetscip, targetconsexprhdlr) );

   return SCIP_OKAY;
}


/** callback to free expression specific data */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataQuotient)
{  /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   /* free expression data of nonlinear handler */
   SCIP_CALL( exprdataFree(scip, nlhdlrexprdata) );

   return SCIP_OKAY;
}


/** callback to detect structure in expression tree */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectQuotient)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);

   /* call detection routine */
   SCIP_CALL( detectExpr(scip, conshdlr, expr, nlhdlrexprdata, success) );

   return SCIP_OKAY;
}


/** auxiliary evaluation callback of nonlinear handler */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalauxQuotient)
{ /*lint --e{715}*/
   SCIP_Real solvalx;
   SCIP_Real solvaly;
   SCIP_Real nomval;
   SCIP_Real denomval;

   assert(expr != NULL);
   assert(auxvalue != NULL);

   solvalx = SCIPgetSolVal(scip, sol, nlhdlrexprdata->nomvar);
   solvaly = SCIPgetSolVal(scip, sol, nlhdlrexprdata->denomvar);
   nomval = nlhdlrexprdata->nomcoef *  solvalx + nlhdlrexprdata->nomconst;
   denomval = nlhdlrexprdata->denomcoef *  solvaly + nlhdlrexprdata->denomconst;

   /* return SCIP_INVALID if the denominator evaluates to zero */
   *auxvalue = (denomval != 0.0) ? nlhdlrexprdata->constant + nomval / denomval : SCIP_INVALID;

   return SCIP_OKAY;
}


/** nonlinear handler under/overestimation callback
 *
 * @todo which of the paramters did I not use, but have to be taken into consideration?
*/
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimateQuotient)
{ /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);

   if( nlhdlrexprdata->nomvar == nlhdlrexprdata->denomvar )
   {
      SCIP_CALL( sepaUnivariate(scip, sol, nlhdlrexprdata->nomvar, SCIPgetConsExprExprAuxVar(expr),
            nlhdlrexprdata->nomcoef, nlhdlrexprdata->nomconst, nlhdlrexprdata->denomcoef,
            nlhdlrexprdata->denomconst, nlhdlrexprdata->constant, overestimate, rowprep, success) );
   }
   else
   {
      assert(nlhdlrexprdata->nomcoef == 1.0);
      assert(nlhdlrexprdata->nomconst == 0.0);
      assert(nlhdlrexprdata->denomcoef == 1.0);
      assert(nlhdlrexprdata->denomconst == 0.0);

      SCIP_CALL( sepaBivariate(scip, sol, nlhdlrexprdata->nomvar, nlhdlrexprdata->denomvar,
         SCIPgetConsExprExprAuxVar(expr), overestimate, rowprep, success) );
   }

   assert(!(*success) || rowprep != NULL);

   return SCIP_OKAY;
}


/** nonlinear handler interval evaluation callback */
static
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(nlhdlrIntevalQuotient)
{ /*lint --e{715}*/
   SCIP_INTERVAL varbnds;
   SCIP_INTERVAL tmp;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nomvar != NULL);
   assert(nlhdlrexprdata->denomvar != NULL);

   /* it is not possible to compute tighter intervals if both variables are different */
   if( nlhdlrexprdata->nomvar != nlhdlrexprdata->denomvar )
      return SCIP_OKAY;

   SCIPintervalSetBounds(&varbnds, SCIPvarGetLbLocal(nlhdlrexprdata->nomvar),
      SCIPvarGetUbLocal(nlhdlrexprdata->nomvar));

   tmp = intEvalQuotient(scip, varbnds, nlhdlrexprdata->nomcoef, nlhdlrexprdata->nomconst,
      nlhdlrexprdata->denomcoef, nlhdlrexprdata->denomconst, nlhdlrexprdata->constant);

   /* intersect intervals if we have learned a tighter interval */
   if( SCIPisGT(scip, tmp.inf, (*interval).inf) || SCIPisLT(scip, tmp.sup, (*interval).sup) )
      SCIPintervalIntersect(interval, *interval, tmp);

   return SCIP_OKAY;
}


/** nonlinear handler callback for reverse propagation */
static
SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(nlhdlrReversepropQuotient)
{ /*lint --e{715}*/
   SCIP_INTERVAL exprbounds;
   SCIP_INTERVAL result;
   SCIP_Real varlb;
   SCIP_Real varub;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nomvar != NULL);
   assert(nlhdlrexprdata->denomvar != NULL);

   /* it is not possible to compute tighter intervals if both variables are different */
   if( nlhdlrexprdata->nomvar != nlhdlrexprdata->denomvar )
      return SCIP_OKAY;

   exprbounds = SCIPgetConsExprExprActivity(scip, expr);
   varlb = SCIPvarGetLbLocal(nlhdlrexprdata->nomvar);
   varub = SCIPvarGetUbLocal(nlhdlrexprdata->nomvar);
   SCIPdebugMsg(scip, "call reverse propagation for expression (%g x + %g) / (%g y + %g) + %g bounds [%g,%g]\n",
      nlhdlrexprdata->nomcoef, nlhdlrexprdata->nomconst, nlhdlrexprdata->denomcoef, nlhdlrexprdata->denomconst,
      nlhdlrexprdata->constant, exprbounds.inf, exprbounds.sup);

   result = reversepropQuotient(exprbounds, nlhdlrexprdata->nomcoef, nlhdlrexprdata->nomconst,
      nlhdlrexprdata->denomcoef, nlhdlrexprdata->denomconst, nlhdlrexprdata->constant);

   if( SCIPisLT(scip, varlb, result.inf) || SCIPisGT(scip, varub, result.sup) )
   {
      SCIP_INTERVAL varbnds;
      SCIP_Bool tightened;

      /* if force=TRUE, take the bound strengthening tolerance into account */
      if( !force && !SCIPisLbBetter(scip, result.inf, varlb, varub) && !SCIPisUbBetter(scip, result.sup, varlb, varub) )
         return SCIP_OKAY;

      SCIPintervalSetBounds(&varbnds, varlb, varub);
      SCIPintervalIntersect(&result, result, varbnds);

      /* tighten bounds of x */
      SCIPdebugMsg(scip, "try to tighten bounds of %s: [%g,%g] -> [%g,%g]\n",
         SCIPvarGetName(nlhdlrexprdata->nomvar), varlb, varub, result.inf, result.sup);

      SCIP_CALL( SCIPtightenVarLb(scip, nlhdlrexprdata->nomvar, result.inf, force,
         infeasible, &tightened) );

      if( tightened )
         ++(*nreductions);

      if( !(*infeasible) )
      {
         SCIP_CALL( SCIPtightenVarUb(scip, nlhdlrexprdata->nomvar, result.sup, force,
            infeasible, &tightened) );

         if( tightened )
            ++(*nreductions);
      }
   }

   return SCIP_OKAY;
}


/*
 * nonlinear handler specific interface methods
 */

/** includes Quotient nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeConsExprNlhdlrQuotient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_CONSEXPR_NLHDLR* nlhdlr;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   /* create nonlinear handler data */
   nlhdlrdata = NULL;

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME,
      NLHDLR_DESC, NLHDLR_PRIORITY, nlhdlrDetectQuotient, nlhdlrEvalauxQuotient, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrQuotient);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrFreeExprDataQuotient);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, NULL, nlhdlrEstimateQuotient, NULL);
   SCIPsetConsExprNlhdlrProp(scip, nlhdlr, nlhdlrIntevalQuotient, nlhdlrReversepropQuotient);

   return SCIP_OKAY;
}
