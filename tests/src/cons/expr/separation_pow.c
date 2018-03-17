/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   separation_prod.c
 * @brief  tests separation of power expressions
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_pow.c"
#include "separation.h"

/* test estimateTangent */
Test(estimation, tangent, .description = "test computation of tangent")
{
   SCIP_Real exponent;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&scip) );

   for( exponent = -3.0; exponent <= 3.0; exponent += 0.5 )
   {
      if( exponent == 0.0 )
         continue;

      for( xref = -2.0; xref <= 2.0; xref += 1.0 )
      {
         /* skip negative reference points when exponent is fractional */
         if( xref < 0.0 && !EPSISINT(exponent, 0.0) )
            continue;

         /* skip zero reference point when exponent is negative */
         if( xref == 0.0 && exponent < 0.0 )
            continue;

         success = FALSE;
         constant = DBL_MAX;
         slope = DBL_MAX;

         estimateTangent(scip, exponent, xref, &constant, &slope, &success);

         /* x^p -> x0^p + p*x0^{p-1} (x-x0) */

         /* estimateTangent must fail iff xref is 0 and exponent < 1 (infinite gradient in reference point) */
         cr_assert(success != (xref == 0.0 && exponent < 1.0));

         if( success )
         {
            cr_assert(SCIPisEQ(scip, slope, exponent * pow(xref, exponent-1.0)));
            cr_assert(SCIPisEQ(scip, constant, pow(xref, exponent) - slope * xref));
         }
      }
   }

   SCIP_CALL( SCIPfree(&scip) );
}

/* test estimateSecant */
Test(estimation, secant, .description = "test computation of secant")
{
   SCIP_Real exponent;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&scip) );

   for( exponent = -3.0; exponent <= 3.0; exponent += 0.5 )
   {
      if( exponent == 0.0 )
         continue;

      for( xlb = -2.0; xlb <= 2.0; xlb += 1.0 )
      {
         for( xub = xlb + 1.0; xub <= 3.0; xub += 1.0 )
         {
            /* skip negative lower bound when exponent is fractional */
            if( xlb < 0.0 && !EPSISINT(exponent, 0.0) )
               continue;

            success = FALSE;
            constant = DBL_MAX;
            slope = DBL_MAX;

            estimateSecant(scip, exponent, xlb, xub, &constant, &slope, &success);

            /* x^p -> xlb^p + (xub^p - xlb^p) / (xub - xlb) * (x - xlb) */

            /* estimateSecant must fail iff xlb or xub is 0 and exponent < 0 (pole at boundary) */
            cr_assert(success != ((xlb == 0.0 || xub == 0.0) && exponent < 0.0));

            if( success )
            {
               cr_assert(SCIPisEQ(scip, slope, (pow(xub, exponent) - pow(xlb, exponent)) / (xub - xlb)));
               cr_assert(SCIPisEQ(scip, constant, pow(xlb, exponent) - slope * xlb));
            }
         }
      }
   }

   /* do one more test where cancellation is likely
    * cancellation when computing slope occurs, e.g., when xub^exponent - xlb^exponent is too small
    * in double precision, with xlb = 1 and xub = 1 + 2*SCIPepsilon (estimateSecant forbids SCIPisEQ(xlb,xub)), this means
    *     (1+2*SCIPepsilon)^exponent - 1 < DBL_EPSILON
    * <-> 1+2*SCIPepsilon < (1+DBL_EPSILON)^(1/exponent)
    * <-> log(1+2*SCIPepsilon) < 1/exponent * log(1+DBL_EPSILON)
    * <-> exponent < log(1+DBL_EPSILON) / log(1+2*SCIPepsilon)
    */
   xlb = 1.0;
   xub = 1.0 + 2 * SCIPepsilon(scip);
   exponent = log(1+DBL_EPSILON) / log(xub) / 2.0;
   cr_assert(exponent > 0.0);  /* exponent is about 1e-7, so we look at a very very flat power function */
   cr_assert(xlb < xub);

   /* in double precision, xlb^exponent looks the same as xub^exponent */
   cr_assert_eq(pow(xlb, exponent), pow(xub, exponent));

   estimateSecant(scip, exponent, xlb, xub, &constant, &slope, &success);

   /* estimateSecant should either fail or produce a positive slope */
   cr_assert(!success || (slope > 0.0));


   SCIP_CALL( SCIPfree(&scip) );
}

/* test estimateParabola */
Test(estimation, parabola, .description = "test computation of parabola estimators")
{
   SCIP_Real exponent;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Bool islocal;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&scip) );

   for( exponent = 1.5; exponent <= 4.0; exponent += 0.5 )
   {
      /* if exponent not even, then start at 0 (otherwise not parabola) */
      for( xref = EPSISINT(exponent/2.0, 0.0) ? -2.0 : 0.0; xref <= 2.0; xref += 1.5 )
      {
         success = FALSE;
         islocal = TRUE;
         constant = DBL_MAX;
         slope = DBL_MAX;

         /* check underestimator (-> tangent) */
         estimateParabola(scip, exponent, FALSE, xref, xref+1.0, xref, &constant, &slope, &islocal, &success);

         cr_assert(success);
         cr_assert(!islocal);
         cr_assert(SCIPisEQ(scip, constant + slope * xref, pow(xref, exponent)));  /* should touch in reference point */
         cr_assert(SCIPisLE(scip, constant + slope * (xref+1.0), pow(xref+1.0, exponent)));  /* should be underestimating in xref+1 */

         /* check overestimator (-> secant) */
         xlb = xref;
         for( xub = xlb + 1.0; xub <= xlb + 2.0; xub += 1.0 )
         {
            success = FALSE;
            islocal = FALSE;
            constant = DBL_MAX;
            slope = DBL_MAX;

            estimateParabola(scip, exponent, TRUE, xlb, xub, (xlb + xub)/2.0, &constant, &slope, &islocal, &success);

            cr_assert(success);
            cr_assert(islocal);
            cr_assert(SCIPisEQ(scip, constant + slope * xlb, pow(xlb, exponent)));  /* should touch at bounds */
            cr_assert(SCIPisEQ(scip, constant + slope * xub, pow(xub, exponent)));  /* should touch at bounds */
            cr_assert(SCIPisGE(scip, constant + slope * (xlb + xub)/2.0, pow((xlb + xub)/2.0, exponent)));  /* should be overestimating in middle point */
         }
      }
   }

   SCIP_CALL( SCIPfree(&scip) );
}

/* test estimateHyperbolaPositive */
Test(estimation, hyperbolaPositive, .description = "test computation of estimators for positive hyperbola")
{
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Bool islocal;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&scip) );

   /* x^(-2) on [-infty,+infty] */
   success = FALSE;
   constant = slope = 5.0;
   estimateHyperbolaPositive(scip, -2.0, FALSE, -SCIPinfinity(scip), SCIPinfinity(scip), -0.5, &constant, &slope, &islocal, &success);
   cr_assert(success); /* underestimator == 0 */
   cr_assert(!islocal);
   cr_assert_eq(constant, 0.0);
   cr_assert_eq(slope, 0.0);

   /* x^(-2) on [-infty,1.0] */
   success = TRUE;
   estimateHyperbolaPositive(scip, -2.0, FALSE, -SCIPinfinity(scip), 1.0, -0.5, &constant, &slope, &islocal, &success);
   cr_assert(!success); /* underestimator not implemented yet */

   /* x^(-2) on [-1.0,infty] */
   success = TRUE;
   estimateHyperbolaPositive(scip, -2.0, FALSE, -1.0, SCIPinfinity(scip), -0.5, &constant, &slope, &islocal, &success);
   cr_assert(!success); /* underestimator not implemented yet */

   /* x^(-2) on [-1,1] */
   success = TRUE;
   estimateHyperbolaPositive(scip, -2.0, FALSE, -1.0, 1.0, -0.5, &constant, &slope, &islocal, &success);
   cr_assert(!success); /* underestimator not implemented yet */

   success = TRUE;
   estimateHyperbolaPositive(scip, -2.0, TRUE, -1.0, 1.0, -0.5, &constant, &slope, &islocal, &success);
   cr_assert(!success); /* overestimator does not exist (or equals infty) */

   /* x^(-2) on [-2,-1] -> underestimator = tangent, overestimator = secant */
   success = FALSE;
   xref = -1.5;
   estimateHyperbolaPositive(scip, -2.0, FALSE, -2.0, -1.0, xref, &constant, &slope, &islocal, &success);
   cr_assert(success);
   cr_assert(!islocal);
   cr_assert(SCIPisEQ(scip, slope, -2.0 * pow(xref, -3.0)));  /* exponent * xref^(exponent-1) */
   cr_assert(SCIPisEQ(scip, constant, pow(xref, -2.0) - slope * xref));

   success = FALSE;
   estimateHyperbolaPositive(scip, -2.0, TRUE, -2.0, -1.0, xref, &constant, &slope, &islocal, &success);
   cr_assert(success);
   cr_assert(islocal);
   cr_assert(SCIPisEQ(scip, slope, (pow(-1.0, -2.0) - pow(-2.0, -2.0))));
   cr_assert(SCIPisEQ(scip, constant, pow(-2.0, -2.0) - slope * (-2.0)));

   /* x^(-2) on [1, 2] -> underestimator = tangent, overestimator = secant */
   success = FALSE;
   xref = 1.5;
   estimateHyperbolaPositive(scip, -2.0, FALSE, 1.0, 2.0, xref, &constant, &slope, &islocal, &success);
   cr_assert(success);
   cr_assert(!islocal);
   cr_assert(SCIPisEQ(scip, slope, -2.0 * pow(xref, -3.0)));  /* exponent * xref^(exponent-1) */
   cr_assert(SCIPisEQ(scip, constant, pow(xref, -2.0) - slope * xref));

   success = FALSE;
   estimateHyperbolaPositive(scip, -2.0, TRUE, 1.0, 2.0, xref, &constant, &slope, &islocal, &success);
   cr_assert(success);
   cr_assert(islocal);
   cr_assert(SCIPisEQ(scip, slope, -0.75));  /* (2^(-2) - 1^(-2)) / (2-1) */
   cr_assert(SCIPisEQ(scip, constant, 1.75)); /* 1^(-2) - slope * 1 */

   SCIP_CALL( SCIPfree(&scip) );
}

/* test estimateHyperbolaMixed */
Test(estimation, hyperbolaMixed, .description = "test computation of estimators for mixed-sign hyperbola")
{
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Bool islocal;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&scip) );

   /* x^(-3) on [-1.0,1.0] */
   success = TRUE;
   estimateHyperbolaMixed(scip, -3.0, FALSE, -1.0, 1.0, -0.5, &constant, &slope, &islocal, &success);
   cr_assert(!success); /* underestimator does not exist (pole in domain) */

   success = TRUE;
   estimateHyperbolaMixed(scip, -3.0, TRUE, -1.0, 1.0, -0.5, &constant, &slope, &islocal, &success);
   cr_assert(!success); /* overestimator does not exist (pole in domain) */

   /* x^(-3) on [-1.0,0.0] -> underestimator does not exist (upper bound is pole); overestimator is tangent */
   success = TRUE;
   xref = -0.5;
   estimateHyperbolaMixed(scip, -3.0, FALSE, -1.0, 0.0, xref, &constant, &slope, &islocal, &success);
   cr_assert(!success); /*  */

   success = FALSE;
   constant = slope = 5.0;
   estimateHyperbolaMixed(scip, -3.0, TRUE, -1.0, 0.0, xref, &constant, &slope, &islocal, &success);
   cr_assert(success);
   cr_assert(!islocal);
   cr_assert(SCIPisEQ(scip, slope, -3.0 * pow(xref, -4.0)));
   cr_assert(SCIPisEQ(scip, constant, pow(xref, -3.0) - slope * xref));

   /* x^(-3) on [-2.0,-1.0] -> underestimator is secant */
   success = FALSE;
   constant = slope = 5.0;
   estimateHyperbolaMixed(scip, -3.0, FALSE, -2.0, -1.0, xref, &constant, &slope, &islocal, &success);
   cr_assert(success); /* underestimator does not exist (upper bound is pole) */
   cr_assert(islocal);
   cr_assert(SCIPisEQ(scip, slope, -1.0 - pow(-2.0, -3.0)));
   cr_assert(SCIPisEQ(scip, constant, pow(-2.0, -3.0) - slope * (-2.0)));

   /* x^(-3) on [-infty,-1.0] -> underestimator does not exist */
   estimateHyperbolaMixed(scip, -3.0, FALSE, -SCIPinfinity(scip), -1.0, xref, &constant, &slope, &islocal, &success);
   cr_assert(!success);

   SCIP_CALL( SCIPfree(&scip) );
}

Test(separation, convexsquare, .init = setup, .fini = teardown,
   .description = "test separation for a convex square expression"
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROW* cut;
   int i;

   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &expr, xexpr, 2.0) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   /*
    * compute cut for w = x^2 with x* = 1.0 and w* = -5.0
    * this should result in an gradient cut at x*
    */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -5.0) );

   cut = NULL;
   SCIP_CALL( separatePointPow(scip,  conshdlr, expr, sol, 0.0, FALSE, &cut) );
   cr_assert(cut != NULL);

   cr_assert_eq(SCIProwGetNNonz(cut), 2);
   cr_assert_eq(SCIProwGetLhs(cut), -SCIPinfinity(scip));
   cr_assert_eq(SCIProwGetRhs(cut), 1.0);

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_assert_eq(coef, 2.0);
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_assert_eq(coef, -1.0);
      else
         cr_assert(FALSE, "found an unknown variable");
   }
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   /*
    * compute cut for w = x^2 with x* = 1.0 and w* = 5.0
    * this should result in a secant cut
    */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 5.0) );

   cut = NULL;
   SCIP_CALL( separatePointPow(scip,  conshdlr, expr, sol, 0.0, TRUE, &cut) );
   cr_assert(cut != NULL);
   cr_assert_eq(SCIProwGetNNonz(cut), 2);
   cr_assert_eq(SCIProwGetLhs(cut), -5.0);
   cr_assert_eq(SCIProwGetRhs(cut), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_assert_eq(coef, 4.0);
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_assert_eq(coef, -1.0);
      else
         cr_assert(FALSE, "found an unknown variable");
   }
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
