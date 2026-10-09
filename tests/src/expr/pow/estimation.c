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

/**@file   estimation.c
 * @brief  tests estimation of power and signed power expressions
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/expr_pow.c"
#include "../estimation.h"

/* test computeTangent */
/** @brief test computation of tangent */
void test_estimation_tangent(void)
{
   SCIP* localscip;
   SCIP_Real exponent;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Bool success;
   unsigned int signpower;

   SCIP_CALL( SCIPcreate(&localscip) );

   for( exponent = -3.0; exponent <= 3.0; exponent += 0.5 )
   {
      if( exponent == 0.0 )
         continue;

      for( signpower = 0; signpower <= 1; ++signpower )
      {
         for( xref = -2.0; xref <= 2.0; xref += 1.0 )
         {
            /* skip negative reference points when exponent is fractional and not signpower */
            if( xref < 0.0 && !EPSISINT(exponent, 0.0) && !signpower )
               continue;

            /* skip zero reference point when exponent is negative */
            if( xref == 0.0 && exponent < 0.0 )
               continue;

            success = FALSE;
            constant = DBL_MAX;
            slope = DBL_MAX;

            computeTangent(localscip, signpower, exponent, xref, &constant, &slope, &success);

            /* normal: x^p -> x0^p + p*x0^{p-1} (x-x0)
             * signpower with x0 < 0: x^p -> -(-x0)^p + p*(-x0)^{p-1} (x-x0)
             */

            /* computeTangent must fail iff xref is 0 and exponent < 1 (infinite gradient in reference point) */
            TEST_ASSERT(success != (xref == 0.0 && exponent < 1.0));

            if( success )
            {
               if( !signpower )
               {
                  TEST_ASSERT(SCIPisEQ(localscip, slope, exponent * pow(xref, exponent-1.0)));
                  TEST_ASSERT(SCIPisEQ(localscip, constant, pow(xref, exponent) - slope * xref));
               }
               else
               {
                  TEST_ASSERT(SCIPisEQ(localscip, slope, exponent * pow(REALABS(xref), exponent-1.0)));
                  TEST_ASSERT(SCIPisEQ(localscip, constant, SIGN(xref) * pow(REALABS(xref), exponent) - slope * xref));
               }
            }
         }
      }
   }

   SCIP_CALL( SCIPfree(&localscip) );
}

/* test computeSecant */
/** @brief test computation of secant */
void test_estimation_secant(void)
{
   SCIP* localscip;
   SCIP_Real exponent;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Bool success;
   unsigned int signpower;

   SCIP_CALL( SCIPcreate(&localscip) );

   for( exponent = -3.0; exponent <= 3.0; exponent += 0.5 )
   {
      if( exponent == 0.0 || exponent == 1.0 )
         continue;

      for( signpower = 0; signpower <= 1; ++signpower )
      {
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

               computeSecant(localscip, signpower, exponent, xlb, xub, &constant, &slope, &success);

               /* f(x) -> f(xlb) + (f(xub) - f(xlb)) / (xub - xlb) * (x - xlb) */

               /* computeSecant must fail iff xlb or xub is 0 and exponent < 0 (pole at boundary) */
               TEST_ASSERT(success != ((xlb == 0.0 || xub == 0.0) && exponent < 0.0));

               if( success )
               {
                  if( !signpower )
                  {
                     TEST_ASSERT(SCIPisEQ(localscip, slope, (pow(xub, exponent) - pow(xlb, exponent)) / (xub - xlb)));
                     TEST_ASSERT(SCIPisEQ(localscip, constant, pow(xlb, exponent) - slope * xlb));
                  }
                  else
                  {
                     TEST_ASSERT(SCIPisEQ(localscip, slope, (SIGN(xub) * pow(REALABS(xub), exponent) - SIGN(xlb) * pow(REALABS(xlb), exponent)) / (xub - xlb)));
                     TEST_ASSERT(SCIPisEQ(localscip, constant, SIGN(xlb) * pow(REALABS(xlb), exponent) - slope * xlb));
                  }
               }
            }
         }
      }
   }

   /* do one more test where cancellation is likely
    * cancellation when computing slope occurs, e.g., when xub^exponent - xlb^exponent is too small
    * in double precision, with xlb = 1 and xub = 1 + 2*SCIPepsilon, this means
    *     (1+2*SCIPepsilon)^exponent - 1 < DBL_EPSILON
    * <-> 1+2*SCIPepsilon < (1+DBL_EPSILON)^(1/exponent)
    * <-> log(1+2*SCIPepsilon) < 1/exponent * log(1+DBL_EPSILON)
    * <-> exponent < log(1+DBL_EPSILON) / log(1+2*SCIPepsilon)
    */
   xlb = 1.0;
   xub = 1.0 + 2 * SCIPepsilon(localscip);
   exponent = log(1+DBL_EPSILON) / log(xub) / 2.0;
   TEST_ASSERT(exponent > 0.0);  /* exponent is about 1e-7, so we look at a very very flat power function */
   TEST_ASSERT(xlb < xub);

   /* in double precision, xlb^exponent looks the same as xub^exponent */
   /* TEST_ASSERT_EQUAL(pow(xlb, exponent), pow(xub, exponent)); */ /* assert fails only on some architectures */

   computeSecant(localscip, FALSE, exponent, xlb, xub, &constant, &slope, &success);

   /* computeSecant should either fail or produce a positive slope */
   TEST_ASSERT(!success || (slope > 0.0));


   /* do one more test where cancellation is even more likely, but is circumvented in computeSecant
    * similar to above, but with xlb = 1 and xub = 1 + 0.5*SCIPepsilon  (computeSecant() checks SCIPisEQ(xlb,xub))
    */
   xlb = 1.0;
   xub = 1.0 + 0.5 * SCIPepsilon(localscip);
   exponent = log(1+DBL_EPSILON) / log(xub) / 2.0;
   TEST_ASSERT(exponent > 0.0);
   TEST_ASSERT(xlb < xub);

   computeSecant(localscip, FALSE, exponent, xlb, xub, &constant, &slope, &success);

   /* in double precision, xlb^exponent looks the same as xub^exponent */
   TEST_ASSERT_EQUAL(pow(xlb, exponent), pow(xub, exponent)); /* assert fails only on some architectures? */

   /* computeSecant should not fail but produce a positive slope */
   TEST_ASSERT(slope > 0.0);
   TEST_ASSERT(SCIPisEQ(localscip, constant, pow(xlb, exponent) - slope * xlb));

   SCIP_CALL( SCIPfree(&localscip) );
}

/* test estimateParabola */
/** @brief test computation of parabola estimators */
void test_estimation_parabola(void)
{
   SCIP* localscip;
   SCIP_Real exponent;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Bool islocal;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&localscip) );

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
         estimateParabola(localscip, exponent, FALSE, xref, xref+1.0, xref, &constant, &slope, &islocal, &success);

         TEST_ASSERT(success);
         TEST_ASSERT(!islocal);
         TEST_ASSERT(SCIPisEQ(localscip, constant + slope * xref, pow(xref, exponent)));  /* should touch in reference point */
         TEST_ASSERT(SCIPisLE(localscip, constant + slope * (xref+1.0), pow(xref+1.0, exponent)));  /* should be underestimating in xref+1 */

         /* check overestimator (-> secant) */
         xlb = xref;
         for( xub = xlb + 1.0; xub <= xlb + 2.0; xub += 1.0 )
         {
            success = FALSE;
            islocal = FALSE;
            constant = DBL_MAX;
            slope = DBL_MAX;

            estimateParabola(localscip, exponent, TRUE, xlb, xub, (xlb + xub)/2.0, &constant, &slope, &islocal, &success);

            TEST_ASSERT(success);
            TEST_ASSERT(islocal);
            TEST_ASSERT(SCIPisEQ(localscip, constant + slope * xlb, pow(xlb, exponent)));  /* should touch at bounds */
            TEST_ASSERT(SCIPisEQ(localscip, constant + slope * xub, pow(xub, exponent)));  /* should touch at bounds */
            TEST_ASSERT(SCIPisGE(localscip, constant + slope * (xlb + xub)/2.0, pow((xlb + xub)/2.0, exponent)));  /* should be overestimating in middle point */
         }
      }
   }

   SCIP_CALL( SCIPfree(&localscip) );
}

/* test computeSignpowerRoot */
/** @brief test calculation of roots for signpower estimators */
void test_estimation_signpower_root(void)
{
   SCIP* localscip;
   SCIP_Real exponent;
   SCIP_Real root;

   SCIP_CALL( SCIPcreate(&localscip) );

   /* try integer exponents, includes lookup table */
   for( exponent = 2.0; exponent < 20.0; exponent += 1.0 )
   {
      SCIP_CALL( computeSignpowerRoot(localscip, &root, exponent) );
      TEST_ASSERT(root > 0.0);
      TEST_ASSERT(root < 1.0);
      /* check that root is a root of (n-1) y^n + n y^(n-1) - 1 */
      TEST_ASSERT(SCIPisEQ(localscip, (exponent-1) * pow(root, exponent) + exponent * pow(root, exponent - 1.0), 1.0));
   }

   /* try some rational exponents and also bigger ones (for exponent 95, Newton fails, but that is crazy anyway) */
   for( exponent = 1.1; exponent < 70.0; exponent *= 1.5 )
   {
      SCIP_CALL( computeSignpowerRoot(localscip, &root, exponent) );
      TEST_ASSERT(root > 0.0);
      TEST_ASSERT(root < 1.0);
      /* check that root is a root of (n-1) y^n + n y^(n-1) - 1 */
      TEST_ASSERT(SCIPisEQ(localscip, (exponent-1) * pow(root, exponent) + exponent * pow(root, exponent - 1.0), 1.0));
   }

   /* try a special rational exponent (has a lookup) */
   exponent = 1.852;
   SCIP_CALL( computeSignpowerRoot(localscip, &root, exponent) );
   TEST_ASSERT(root > 0.0);
   TEST_ASSERT(root < 1.0);
   /* check that root is a root of (n-1) y^n + n y^(n-1) - 1 */
   TEST_ASSERT(SCIPisEQ(localscip, (exponent-1) * pow(root, exponent) + exponent * pow(root, exponent - 1.0), 1.0));

   SCIP_CALL( SCIPfree(&localscip) );
}

/* test estimateSignedpower */
/** @brief test computation of signpower estimators */
void test_estimation_signpower(void)
{
   SCIP* localscip;
   SCIP_Real exponent;
   SCIP_Real root;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Bool islocal;
   SCIP_Bool branchcand;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&localscip) );

   for( exponent = 3.0; exponent <= 5.0; exponent += 2.0 )
   {
      /* later I want this loop to also cover even or rational exponents */

      SCIP_CALL( computeSignpowerRoot(localscip, &root, exponent) );

      /* on [-10,-5] and [-10,0], we should get secants (underestimator) and tangents (overestimator) */
      xlb = -10.0;
      for( xub = -5; xub <= 0.0; xub += 5.0 )
      {
         xref = (xlb + xub) / 2.0;

         success = FALSE;
         islocal = FALSE;
         branchcand = TRUE;
         slope = constant = -5;
         estimateSignedpower(localscip, exponent, root, FALSE, xlb, xub, xref, xlb, xub, &constant, &slope, &islocal, &branchcand, &success);
         TEST_ASSERT(success);
         TEST_ASSERT(islocal);
         TEST_ASSERT(branchcand);
         TEST_ASSERT(SCIPisEQ(localscip, -pow(-xlb, exponent), constant + slope * xlb));
         TEST_ASSERT(SCIPisEQ(localscip, -pow(-xub, exponent), constant + slope * xub));

         success = FALSE;
         islocal = TRUE;
         branchcand = TRUE;
         slope = constant = -5;
         estimateSignedpower(localscip, exponent, root, TRUE, xlb, xub, xref, xlb, xub, &constant, &slope, &islocal, &branchcand, &success);
         TEST_ASSERT(success);
         TEST_ASSERT(!islocal);
         TEST_ASSERT(!branchcand);
         TEST_ASSERT(SCIPisEQ(localscip, slope, exponent * pow(-xref, exponent - 1.0)));
         TEST_ASSERT(SCIPisEQ(localscip, -pow(-xref, exponent), constant + slope * xref));

         /* if global upper bound is small enough (< -xref/root), then overestimator should still be global */
         success = FALSE;
         islocal = TRUE;
         branchcand = TRUE;
         estimateSignedpower(localscip, exponent, root, TRUE, xlb, xub, xref, xlb, - xref/root / 2.0 , &constant, &slope, &islocal, &branchcand, &success);
         TEST_ASSERT(success);
         TEST_ASSERT(!islocal);
         TEST_ASSERT(!branchcand);

         /* if global upper bound is too large (> -xref/root), then overestimator is only locally valid */
         success = FALSE;
         islocal = FALSE;
         branchcand = TRUE;
         estimateSignedpower(localscip, exponent, root, TRUE, xlb, xub, xref, xlb, - xref/root * 2.0 , &constant, &slope, &islocal, &branchcand, &success);
         TEST_ASSERT(success);
         TEST_ASSERT(islocal);
         TEST_ASSERT(!branchcand);
      }

      /* on [-10,10] it gets more interesting */
      xub = 10.0;
      for( xref = xlb; xref <= xub; xref += 2.0 )
      {
         /* underestimator is secant for xref < -xlb * root, otherwise tangent */
         success = FALSE;
         islocal = !(xref < -xlb * root);
         branchcand = TRUE;
         slope = constant = -5;
         estimateSignedpower(localscip, exponent, root, FALSE, xlb, xub, xref, xlb, xub, &constant, &slope, &islocal, &branchcand, &success);
         TEST_ASSERT(success);
         if( xref < -xlb * root )
         {
            /* expect secant between xlb and -xlb*root */
            TEST_ASSERT(islocal);
            TEST_ASSERT(branchcand);
            TEST_ASSERT(SCIPisEQ(localscip, -pow(-xlb, exponent), constant + slope * xlb));
            TEST_ASSERT(SCIPisEQ(localscip, pow(-xlb*root, exponent), constant + slope * (-xlb*root)));
         }
         else
         {
            /* expect tangent */
            TEST_ASSERT(!islocal);
            TEST_ASSERT(!branchcand);
            TEST_ASSERT(SCIPisEQ(localscip, slope, exponent * pow(xref, exponent - 1.0)));
            TEST_ASSERT(SCIPisEQ(localscip, constant, pow(xref, exponent) - slope * xref));
         }

         /* overestimator is secant for xref > -xub * root, otherwise tangent */
         success = FALSE;
         islocal = !(xref > -xub * root);
         branchcand = TRUE;
         slope = constant = -5;
         estimateSignedpower(localscip, exponent, root, TRUE, xlb, xub, xref, xlb, xub, &constant, &slope, &islocal, &branchcand, &success);
         TEST_ASSERT(success);
         if( xref > -xub * root )
         {
            /* expect secant between -xub*root and xub */
            TEST_ASSERT(islocal);
            TEST_ASSERT(branchcand);
            TEST_ASSERT(SCIPisEQ(localscip, -pow(xub*root, exponent), constant + slope * (-xub*root)));
            TEST_ASSERT(SCIPisEQ(localscip, pow(xub, exponent), constant + slope * xub));
         }
         else
         {
            /* expect tangent */
            TEST_ASSERT(!islocal);
            TEST_ASSERT(!branchcand);
            TEST_ASSERT(SCIPisEQ(localscip, slope, exponent * pow(xref, exponent - 1.0)));
            TEST_ASSERT(SCIPisEQ(localscip, constant, pow(xref, exponent) - slope * xref));
         }
      }
   }

   SCIP_CALL( SCIPfree(&localscip) );
}

/* test computeHyperbolaRoot */
/** @brief test calculation of roots for positive hyperbola estimators */
void test_estimation_hyperbola_root(void)
{
   SCIP* localscip;
   SCIP_Real exponent;
   SCIP_Real root;

   SCIP_CALL( SCIPcreate(&localscip) );

   /* try odd negative integer exponents (for exponent -42, Newton fails, but that is crazy anyway) */
   for( exponent = -2.0; exponent > -40.0; exponent -= 2.0 )
   {
      SCIP_CALL( computeHyperbolaRoot(localscip, &root, exponent) );
      TEST_ASSERT(root < 0.0);
      /* check that root is a root of (n-1) y^n - n y^(n-1) + 1 */
      TEST_ASSERT(SCIPisZero(localscip, (exponent-1) * pow(root, exponent) - exponent * pow(root, exponent - 1.0) + 1.0));
   }

   SCIP_CALL( SCIPfree(&localscip) );
}

/* test estimateHyperbolaPositive */
/** @brief test computation of estimators for positive hyperbola */
void test_estimation_hyperbolaPositive(void)
{
   SCIP* localscip;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Real root;
   SCIP_Bool islocal;
   SCIP_Bool branchcand;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&localscip) );

   /* compute root for exponent -2 */
   SCIP_CALL( computeHyperbolaRoot(localscip, &root, -2.0) );

   /* x^(-2) on [-infty,+infty] */
   success = FALSE;
   islocal = TRUE;
   branchcand = TRUE;
   constant = slope = 5.0;
   estimateHyperbolaPositive(localscip, -2.0, root, FALSE, -SCIPinfinity(localscip), SCIPinfinity(localscip), -0.5, -SCIPinfinity(localscip), SCIPinfinity(localscip), &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success); /* underestimator == 0 */
   TEST_ASSERT(!islocal);
   TEST_ASSERT(branchcand);
   TEST_ASSERT_EQUAL(constant, 0.0);
   TEST_ASSERT_EQUAL(slope, 0.0);

   /* x^(-2) on [-1,1]; underestimator is secant between -1 and 1 */
   success = FALSE;
   islocal = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(localscip, -2.0, root, FALSE, -1.0, 1.0, -0.5, -1.0, 1.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(islocal);
   TEST_ASSERT(branchcand);
   TEST_ASSERT(constant == 1.0);
   TEST_ASSERT(slope == 0.0);

   /* x^(-2) on [-1.0,infty]; underestimator is secant between -1 and 2 for xref = -0.5 (< 2) */
   success = FALSE;
   islocal = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(localscip, -2.0, root, FALSE, -1.0, SCIPinfinity(localscip), -0.5, -1.0, SCIPinfinity(localscip), &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(islocal);
   TEST_ASSERT(branchcand);
   TEST_ASSERT(SCIPisEQ(localscip, 1, constant + slope * (-1))); /* touch at -1, (-1)^(-2) = 1 */
   TEST_ASSERT(SCIPisEQ(localscip, 0.25, constant + slope * 2)); /* touch at 2, 2^(-2) = 0.25 */

   /* x^(-2) on [-1.0,infty]; underestimator is tangent for xref > 2 */
   success = FALSE;
   islocal = FALSE;
   branchcand = TRUE;
   xref = 4.0;
   estimateHyperbolaPositive(localscip, -2.0, root, FALSE, -1.0, SCIPinfinity(localscip), xref, -1.0, SCIPinfinity(localscip), &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(!islocal); /* the tangent is also globally valid, since global bounds equal local bounds here */
   TEST_ASSERT(!branchcand);
   TEST_ASSERT(SCIPisEQ(localscip, slope, -2.0 * pow(xref, -3.0)));  /* slope should be gradient at xref */
   TEST_ASSERT(SCIPisEQ(localscip, pow(xref, -2.0), constant + slope * xref)); /* touch at xref */

   /* x^(-2) on [-infty,1.0]; underestimator is secant between -2 and 1 */
   success = TRUE;
   islocal = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(localscip, -2.0, root, FALSE, -SCIPinfinity(localscip), 1.0, -0.5, -SCIPinfinity(localscip), 1.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(islocal);
   TEST_ASSERT(branchcand);
   TEST_ASSERT(SCIPisEQ(localscip, 0.25, constant + slope * (-2))); /* touch at -2, (-2)^(-2) = 0.25 */
   TEST_ASSERT(SCIPisEQ(localscip, 1, constant + slope * 1)); /* touch at 1, 1^(-2) = 1 */

   success = TRUE;
   branchcand = TRUE;
   estimateHyperbolaPositive(localscip, -2.0, SCIP_INVALID, TRUE, -1.0, 1.0, -0.5, -1.0, 1.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(!success); /* overestimator does not exist (or equals infty) */
   TEST_ASSERT(branchcand);

   /* x^(-2) on [-2,-1] -> underestimator = tangent, overestimator = secant */
   success = FALSE;
   islocal = TRUE;
   branchcand = TRUE;
   xref = -1.5;
   estimateHyperbolaPositive(localscip, -2.0, SCIP_INVALID, FALSE, -2.0, -1.0, xref, -2.0, -1.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(!islocal);
   TEST_ASSERT(!branchcand);
   TEST_ASSERT(SCIPisEQ(localscip, slope, -2.0 * pow(xref, -3.0)));  /* exponent * xref^(exponent-1) */
   TEST_ASSERT(SCIPisEQ(localscip, constant, pow(xref, -2.0) - slope * xref));

   success = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(localscip, -2.0, root, FALSE, -2.0, -1.0, xref, -2.0, 2.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(islocal);  /* if global domain is [-2,2], then the tangent is not globally valid if xref > -xubglobal = -2 */
   TEST_ASSERT(!branchcand);  /* but branching will not change the tangent */

   success = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(localscip, -2.0, root, FALSE, -2.0, -1.0, xref, -2.0, 0.5, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(!islocal);  /* if global domain is [-2,0.5], then the tangent is globally valid, since xref = -1.5 < xubglobal*root = 0.5*(-2) = -1 */
   TEST_ASSERT(!branchcand);

   success = FALSE;
   islocal = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(localscip, -2.0, SCIP_INVALID, TRUE, -2.0, -1.0, xref, -2.0, -1.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(islocal);
   TEST_ASSERT(branchcand);
   TEST_ASSERT(SCIPisEQ(localscip, slope, (pow(-1.0, -2.0) - pow(-2.0, -2.0))));
   TEST_ASSERT(SCIPisEQ(localscip, constant, pow(-2.0, -2.0) - slope * (-2.0)));

   /* x^(-2) on [1, 2] -> underestimator = tangent, overestimator = secant */
   success = FALSE;
   islocal = TRUE;
   branchcand = TRUE;
   xref = 1.5;
   estimateHyperbolaPositive(localscip, -2.0, SCIP_INVALID, FALSE, 1.0, 2.0, xref, 1.0, 2.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(!islocal);
   TEST_ASSERT(!branchcand);
   TEST_ASSERT(SCIPisEQ(localscip, slope, -2.0 * pow(xref, -3.0)));  /* exponent * xref^(exponent-1) */
   TEST_ASSERT(SCIPisEQ(localscip, constant, pow(xref, -2.0) - slope * xref));

   success = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(localscip, -2.0, root, FALSE, 1.0, 2.0, xref, -1.0, 2.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(islocal);  /* if global domain is [-1,2], then tangent is not globally valid */
   TEST_ASSERT(!branchcand);

   success = FALSE;
   islocal = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(localscip, -2.0, SCIP_INVALID, TRUE, 1.0, 2.0, xref, 1.0, 2.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(islocal);
   TEST_ASSERT(branchcand);
   TEST_ASSERT(SCIPisEQ(localscip, slope, -0.75));  /* (2^(-2) - 1^(-2)) / (2-1) */
   TEST_ASSERT(SCIPisEQ(localscip, constant, 1.75)); /* 1^(-2) - slope * 1 */

   SCIP_CALL( SCIPfree(&localscip) );
}

/* test estimateHyperbolaMixed */
/** @brief test computation of estimators for mixed-sign hyperbola */
void test_estimation_hyperbolaMixed(void)
{
   SCIP* localscip;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Bool islocal;
   SCIP_Bool branchcand;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&localscip) );

   /* x^(-3) on [-1.0,1.0] */
   success = TRUE;
   branchcand = TRUE;
   estimateHyperbolaMixed(localscip, -3.0, FALSE, -1.0, 1.0, -0.5, -1.0, 1.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(!success); /* underestimator does not exist (pole in domain) */
   TEST_ASSERT(branchcand);

   success = TRUE;
   branchcand = TRUE;
   estimateHyperbolaMixed(localscip, -3.0, TRUE, -1.0, 1.0, -0.5, -1.0, 1.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(!success); /* overestimator does not exist (pole in domain) */
   TEST_ASSERT(branchcand);

   /* x^(-3) on [-1.0,0.0] -> underestimator does not exist (upper bound is pole); overestimator is tangent */
   success = TRUE;
   branchcand = TRUE;
   xref = -0.5;
   estimateHyperbolaMixed(localscip, -3.0, FALSE, -1.0, 0.0, xref, -1.0, 0.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(!success);
   TEST_ASSERT(branchcand);

   success = FALSE;
   islocal = TRUE;
   branchcand = TRUE;
   constant = slope = 5.0;
   estimateHyperbolaMixed(localscip, -3.0, TRUE, -1.0, 0.0, xref, -1.0, 0.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(!islocal);
   TEST_ASSERT(!branchcand);
   TEST_ASSERT(SCIPisEQ(localscip, slope, -3.0 * pow(xref, -4.0)));
   TEST_ASSERT(SCIPisEQ(localscip, constant, pow(xref, -3.0) - slope * xref));

   success = FALSE;
   branchcand = TRUE;
   estimateHyperbolaMixed(localscip, -3.0, TRUE, -1.0, 0.0, xref, -1.0, 1.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(islocal);  /* if global domain is [-1,1], then tangent is not globally valid */
   TEST_ASSERT(!branchcand);

   /* x^(-3) on [-2.0,-1.0] -> underestimator is secant */
   success = FALSE;
   islocal = FALSE;
   branchcand = TRUE;
   constant = slope = 5.0;
   estimateHyperbolaMixed(localscip, -3.0, FALSE, -2.0, -1.0, xref, -2.0, 0.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(success); /* underestimator does not exist (upper bound is pole) */
   TEST_ASSERT(islocal);
   TEST_ASSERT(branchcand);
   TEST_ASSERT(SCIPisEQ(localscip, slope, -1.0 - pow(-2.0, -3.0)));
   TEST_ASSERT(SCIPisEQ(localscip, constant, pow(-2.0, -3.0) - slope * (-2.0)));

   /* x^(-3) on [-infty,-1.0] -> underestimator does not exist */
   success = TRUE;
   branchcand = TRUE;
   estimateHyperbolaMixed(localscip, -3.0, FALSE, -SCIPinfinity(localscip), -1.0, xref, -SCIPinfinity(localscip), -1.0, &constant, &slope, &islocal, &branchcand, &success);
   TEST_ASSERT(!success);
   TEST_ASSERT(branchcand);

   SCIP_CALL( SCIPfree(&localscip) );
}

/* test SCIPestimateRoot */
/** @brief test computation of estimators for roots (<1) */
void test_estimation_root(void)
{
   SCIP* localscip;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Bool islocal;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&localscip) );

   /* x^0.25 on [0.0,infty] -> underestimator does not exist */
   success = TRUE;
   SCIPestimateRoot(localscip, 0.25, FALSE, 0.0, SCIPinfinity(localscip), 0.5, &constant, &slope, &islocal, &success);
   TEST_ASSERT(!success);

   /* x^0.25 on [0.0,16.0] -> underestimator is secant; overestimator is tangent */
   xref = 4.0;
   success = FALSE;
   islocal = FALSE;
   constant = slope = -5.0;
   SCIPestimateRoot(localscip, 0.25, FALSE, 0.0, 16.0, xref, &constant, &slope, &islocal, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(islocal);
   TEST_ASSERT(SCIPisEQ(localscip, slope, 2.0/16.0));
   TEST_ASSERT(SCIPisEQ(localscip, constant, 0.0));

   success = FALSE;
   islocal = TRUE;
   SCIPestimateRoot(localscip, 0.25, TRUE, 0.0, 16.0, xref, &constant, &slope, &islocal, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(!islocal);
   TEST_ASSERT(SCIPisEQ(localscip, slope, 0.25 * pow(xref, -0.75)));
   TEST_ASSERT(SCIPisEQ(localscip, constant, pow(xref, 0.25) - slope * xref));

   /* if reference point at 0.0, then tangent will still be computed, but it will not touch at 0.0 */
   success = FALSE;
   islocal = TRUE;
   SCIPestimateRoot(localscip, 0.25, TRUE, 0.0, 16.0, 0.0, &constant, &slope, &islocal, &success);
   TEST_ASSERT(success);
   TEST_ASSERT(!islocal);
   TEST_ASSERT(constant != 0.0);

   /* if reference point at 0.0 and bounds on x are very small, then no estimator is computed */
   SCIPestimateRoot(localscip, 0.25, TRUE, 0.0, SCIPepsilon(localscip), 0.0, &constant, &slope, &islocal, &success);
   TEST_ASSERT(!success);

   SCIP_CALL( SCIPfree(&localscip) );
}

/** @brief test separation for a convex square expression */
void test_estimation_convexsquare(void)
{
   SCIP_EXPR* expr;
   SCIP_Real xval;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Bool islocal;
   SCIP_Bool branchcand;
   SCIP_Bool success;
   SCIP_INTERVAL bnd;

   SCIP_CALL( SCIPcreateExprPow(scip, &expr, xexpr, 2.0, NULL, NULL) );
   SCIP_CALL( SCIPevalExprActivity(scip, expr) );
   bnd = SCIPexprGetActivity(xexpr);

   /*
    * compute underestimator for x^2 with x* = 1.0
    * this should result in an gradient estimator
    */
   xval = 1.0;

   branchcand = TRUE;
   SCIP_CALL( estimatePow(scip, expr, &bnd, &bnd, &xval, FALSE, SCIPinfinity(scip), &slope, &constant, &islocal, &success, &branchcand) );

   TEST_ASSERT(success);
   SOFT_ASSERT_DOUBLE_WITHIN(constant, -1.0, SCIPepsilon(scip));
   SOFT_ASSERT_DOUBLE_WITHIN(slope, 2.0, SCIPepsilon(scip));
   SOFT_ASSERT(!islocal);
   SOFT_ASSERT(!branchcand);

   /*
    * compute overestimator for x^2 with x* = 1.0
    * this should result in a secant estimator
    */
   xval = 1.0;

   branchcand = TRUE;
   SCIP_CALL( estimatePow(scip, expr, &bnd, &bnd, &xval, TRUE, -SCIPinfinity(scip), &slope, &constant, &islocal, &success, &branchcand) );
   TEST_ASSERT(success);
   SOFT_ASSERT_DOUBLE_WITHIN(constant, 5.0, SCIPepsilon(scip));
   SOFT_ASSERT_DOUBLE_WITHIN(slope, 4.0, SCIPepsilon(scip));
   SOFT_ASSERT(islocal);
   SOFT_ASSERT(branchcand);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_estimation_tangent);
   RUN_TEST(test_estimation_secant);
   RUN_TEST(test_estimation_parabola);
   RUN_TEST(test_estimation_signpower_root);
   RUN_TEST(test_estimation_signpower);
   RUN_TEST(test_estimation_hyperbola_root);
   RUN_TEST(test_estimation_hyperbolaPositive);
   RUN_TEST(test_estimation_hyperbolaMixed);
   RUN_TEST(test_estimation_root);
   RUN_TEST(test_estimation_convexsquare);
   return UNITY_END();
}
