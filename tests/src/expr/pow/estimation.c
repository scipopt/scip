/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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
Test(estimation, tangent, .description = "test computation of tangent")
{
   SCIP_Real exponent;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Bool success;
   unsigned int signpower;

   SCIP_CALL( SCIPcreate(&scip) );

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

            computeTangent(scip, signpower, exponent, xref, &constant, &slope, &success);

            /* normal: x^p -> x0^p + p*x0^{p-1} (x-x0)
             * signpower with x0 < 0: x^p -> -(-x0)^p + p*(-x0)^{p-1} (x-x0)
             */

            /* computeTangent must fail iff xref is 0 and exponent < 1 (infinite gradient in reference point) */
            cr_assert(success != (xref == 0.0 && exponent < 1.0));

            if( success )
            {
               if( !signpower )
               {
                  cr_assert(SCIPisEQ(scip, slope, exponent * pow(xref, exponent-1.0)));
                  cr_assert(SCIPisEQ(scip, constant, pow(xref, exponent) - slope * xref));
               }
               else
               {
                  cr_assert(SCIPisEQ(scip, slope, exponent * pow(REALABS(xref), exponent-1.0)));
                  cr_assert(SCIPisEQ(scip, constant, SIGN(xref) * pow(REALABS(xref), exponent) - slope * xref));
               }
            }
         }
      }
   }

   SCIP_CALL( SCIPfree(&scip) );
}

/* test computeSecant */
Test(estimation, secant, .description = "test computation of secant")
{
   SCIP_Real exponent;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Bool success;
   unsigned int signpower;

   SCIP_CALL( SCIPcreate(&scip) );

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

               computeSecant(scip, signpower, exponent, xlb, xub, &constant, &slope, &success);

               /* f(x) -> f(xlb) + (f(xub) - f(xlb)) / (xub - xlb) * (x - xlb) */

               /* computeSecant must fail iff xlb or xub is 0 and exponent < 0 (pole at boundary) */
               cr_assert(success != ((xlb == 0.0 || xub == 0.0) && exponent < 0.0));

               if( success )
               {
                  if( !signpower )
                  {
                     cr_assert(SCIPisEQ(scip, slope, (pow(xub, exponent) - pow(xlb, exponent)) / (xub - xlb)));
                     cr_assert(SCIPisEQ(scip, constant, pow(xlb, exponent) - slope * xlb));
                  }
                  else
                  {
                     cr_assert(SCIPisEQ(scip, slope, (SIGN(xub) * pow(REALABS(xub), exponent) - SIGN(xlb) * pow(REALABS(xlb), exponent)) / (xub - xlb)));
                     cr_assert(SCIPisEQ(scip, constant, SIGN(xlb) * pow(REALABS(xlb), exponent) - slope * xlb));
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
   xub = 1.0 + 2 * SCIPepsilon(scip);
   exponent = log(1+DBL_EPSILON) / log(xub) / 2.0;
   cr_assert(exponent > 0.0);  /* exponent is about 1e-7, so we look at a very very flat power function */
   cr_assert(xlb < xub);

   /* in double precision, xlb^exponent looks the same as xub^exponent */
   /* cr_assert_eq(pow(xlb, exponent), pow(xub, exponent)); */ /* assert fails only on some architectures */

   computeSecant(scip, FALSE, exponent, xlb, xub, &constant, &slope, &success);

   /* computeSecant should either fail or produce a positive slope */
   cr_assert(!success || (slope > 0.0));


   /* do one more test where cancellation is even more likely, but is circumvented in computeSecant
    * similar to above, but with xlb = 1 and xub = 1 + 0.5*SCIPepsilon  (computeSecant() checks SCIPisEQ(xlb,xub))
    */
   xlb = 1.0;
   xub = 1.0 + 0.5 * SCIPepsilon(scip);
   exponent = log(1+DBL_EPSILON) / log(xub) / 2.0;
   cr_assert(exponent > 0.0);
   cr_assert(xlb < xub);

   computeSecant(scip, FALSE, exponent, xlb, xub, &constant, &slope, &success);

   /* in double precision, xlb^exponent looks the same as xub^exponent */
   cr_assert_eq(pow(xlb, exponent), pow(xub, exponent)); /* assert fails only on some architectures? */

   /* computeSecant should not fail but produce a positive slope */
   cr_assert(slope > 0.0);
   cr_assert(SCIPisEQ(scip, constant, pow(xlb, exponent) - slope * xlb));

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

/* test computeSignpowerRoot */
Test(estimation, signpower_root, .description = "test calculation of roots for signpower estimators")
{
   SCIP_Real exponent;
   SCIP_Real root;

   SCIP_CALL( SCIPcreate(&scip) );

   /* try integer exponents, includes lookup table */
   for( exponent = 2.0; exponent < 20.0; exponent += 1.0 )
   {
      SCIP_CALL( computeSignpowerRoot(scip, &root, exponent) );
      cr_assert(root > 0.0);
      cr_assert(root < 1.0);
      /* check that root is a root of (n-1) y^n + n y^(n-1) - 1 */
      cr_assert(SCIPisEQ(scip, (exponent-1) * pow(root, exponent) + exponent * pow(root, exponent - 1.0), 1.0));
   }

   /* try some rational exponents and also bigger ones (for exponent 95, Newton fails, but that is crazy anyway) */
   for( exponent = 1.1; exponent < 70.0; exponent *= 1.5 )
   {
      SCIP_CALL( computeSignpowerRoot(scip, &root, exponent) );
      cr_assert(root > 0.0);
      cr_assert(root < 1.0);
      /* check that root is a root of (n-1) y^n + n y^(n-1) - 1 */
      cr_assert(SCIPisEQ(scip, (exponent-1) * pow(root, exponent) + exponent * pow(root, exponent - 1.0), 1.0));
   }

   /* try a special rational exponent (has a lookup) */
   exponent = 1.852;
   SCIP_CALL( computeSignpowerRoot(scip, &root, exponent) );
   cr_assert(root > 0.0);
   cr_assert(root < 1.0);
   /* check that root is a root of (n-1) y^n + n y^(n-1) - 1 */
   cr_assert(SCIPisEQ(scip, (exponent-1) * pow(root, exponent) + exponent * pow(root, exponent - 1.0), 1.0));

   SCIP_CALL( SCIPfree(&scip) );
}

/* test estimateSignedpower */
Test(estimation, signpower, .description = "test computation of signpower estimators")
{
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

   SCIP_CALL( SCIPcreate(&scip) );

   for( exponent = 3.0; exponent <= 5.0; exponent += 2.0 )
   {
      /* later I want this loop to also cover even or rational exponents */

      SCIP_CALL( computeSignpowerRoot(scip, &root, exponent) );

      /* on [-10,-5] and [-10,0], we should get secants (underestimator) and tangents (overestimator) */
      xlb = -10.0;
      for( xub = -5; xub <= 0.0; xub += 5.0 )
      {
         xref = (xlb + xub) / 2.0;

         success = FALSE;
         islocal = FALSE;
         branchcand = TRUE;
         slope = constant = -5;
         estimateSignedpower(scip, exponent, root, FALSE, xlb, xub, xref, xlb, xub, &constant, &slope, &islocal, &branchcand, &success);
         cr_assert(success);
         cr_assert(islocal);
         cr_assert(branchcand);
         cr_assert(SCIPisEQ(scip, -pow(-xlb, exponent), constant + slope * xlb));
         cr_assert(SCIPisEQ(scip, -pow(-xub, exponent), constant + slope * xub));

         success = FALSE;
         islocal = TRUE;
         branchcand = TRUE;
         slope = constant = -5;
         estimateSignedpower(scip, exponent, root, TRUE, xlb, xub, xref, xlb, xub, &constant, &slope, &islocal, &branchcand, &success);
         cr_assert(success);
         cr_assert(!islocal);
         cr_assert(!branchcand);
         cr_assert(SCIPisEQ(scip, slope, exponent * pow(-xref, exponent - 1.0)));
         cr_assert(SCIPisEQ(scip, -pow(-xref, exponent), constant + slope * xref));

         /* if global upper bound is small enough (< -xref/root), then overestimator should still be global */
         success = FALSE;
         islocal = TRUE;
         branchcand = TRUE;
         estimateSignedpower(scip, exponent, root, TRUE, xlb, xub, xref, xlb, - xref/root / 2.0 , &constant, &slope, &islocal, &branchcand, &success);
         cr_assert(success);
         cr_assert(!islocal);
         cr_assert(!branchcand);

         /* if global upper bound is too large (> -xref/root), then overestimator is only locally valid */
         success = FALSE;
         islocal = FALSE;
         branchcand = TRUE;
         estimateSignedpower(scip, exponent, root, TRUE, xlb, xub, xref, xlb, - xref/root * 2.0 , &constant, &slope, &islocal, &branchcand, &success);
         cr_assert(success);
         cr_assert(islocal);
         cr_assert(!branchcand);
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
         estimateSignedpower(scip, exponent, root, FALSE, xlb, xub, xref, xlb, xub, &constant, &slope, &islocal, &branchcand, &success);
         cr_assert(success);
         if( xref < -xlb * root )
         {
            /* expect secant between xlb and -xlb*root */
            cr_assert(islocal);
            cr_assert(branchcand);
            cr_assert(SCIPisEQ(scip, -pow(-xlb, exponent), constant + slope * xlb));
            cr_assert(SCIPisEQ(scip, pow(-xlb*root, exponent), constant + slope * (-xlb*root)));
         }
         else
         {
            /* expect tangent */
            cr_assert(!islocal);
            cr_assert(!branchcand);
            cr_assert(SCIPisEQ(scip, slope, exponent * pow(xref, exponent - 1.0)));
            cr_assert(SCIPisEQ(scip, constant, pow(xref, exponent) - slope * xref));
         }

         /* overestimator is secant for xref > -xub * root, otherwise tangent */
         success = FALSE;
         islocal = !(xref > -xub * root);
         branchcand = TRUE;
         slope = constant = -5;
         estimateSignedpower(scip, exponent, root, TRUE, xlb, xub, xref, xlb, xub, &constant, &slope, &islocal, &branchcand, &success);
         cr_assert(success);
         if( xref > -xub * root )
         {
            /* expect secant between -xub*root and xub */
            cr_assert(islocal);
            cr_assert(branchcand);
            cr_assert(SCIPisEQ(scip, -pow(xub*root, exponent), constant + slope * (-xub*root)));
            cr_assert(SCIPisEQ(scip, pow(xub, exponent), constant + slope * xub));
         }
         else
         {
            /* expect tangent */
            cr_assert(!islocal);
            cr_assert(!branchcand);
            cr_assert(SCIPisEQ(scip, slope, exponent * pow(xref, exponent - 1.0)));
            cr_assert(SCIPisEQ(scip, constant, pow(xref, exponent) - slope * xref));
         }
      }
   }

   SCIP_CALL( SCIPfree(&scip) );
}

/* test computeHyperbolaRoot */
Test(estimation, hyperbola_root, .description = "test calculation of roots for positive hyperbola estimators")
{
   SCIP_Real exponent;
   SCIP_Real root;

   SCIP_CALL( SCIPcreate(&scip) );

   /* try odd negative integer exponents (for exponent -42, Newton fails, but that is crazy anyway) */
   for( exponent = -2.0; exponent > -40.0; exponent -= 2.0 )
   {
      SCIP_CALL( computeHyperbolaRoot(scip, &root, exponent) );
      cr_assert(root < 0.0);
      /* check that root is a root of (n-1) y^n - n y^(n-1) + 1 */
      cr_assert(SCIPisZero(scip, (exponent-1) * pow(root, exponent) - exponent * pow(root, exponent - 1.0) + 1.0));
   }

   SCIP_CALL( SCIPfree(&scip) );
}

/* test estimateHyperbolaPositive */
Test(estimation, hyperbolaPositive, .description = "test computation of estimators for positive hyperbola")
{
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Real root;
   SCIP_Bool islocal;
   SCIP_Bool branchcand;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&scip) );

   /* compute root for exponent -2 */
   SCIP_CALL( computeHyperbolaRoot(scip, &root, -2.0) );

   /* x^(-2) on [-infty,+infty] */
   success = FALSE;
   islocal = TRUE;
   branchcand = TRUE;
   constant = slope = 5.0;
   estimateHyperbolaPositive(scip, -2.0, root, FALSE, -SCIPinfinity(scip), SCIPinfinity(scip), -0.5, -SCIPinfinity(scip), SCIPinfinity(scip), &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success); /* underestimator == 0 */
   cr_assert(!islocal);
   cr_assert(branchcand);
   cr_assert_eq(constant, 0.0);
   cr_assert_eq(slope, 0.0);

   /* x^(-2) on [-1,1]; underestimator is secant between -1 and 1 */
   success = FALSE;
   islocal = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(scip, -2.0, root, FALSE, -1.0, 1.0, -0.5, -1.0, 1.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(islocal);
   cr_assert(branchcand);
   cr_assert(constant == 1.0);
   cr_assert(slope == 0.0);

   /* x^(-2) on [-1.0,infty]; underestimator is secant between -1 and 2 for xref = -0.5 (< 2) */
   success = FALSE;
   islocal = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(scip, -2.0, root, FALSE, -1.0, SCIPinfinity(scip), -0.5, -1.0, SCIPinfinity(scip), &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(islocal);
   cr_assert(branchcand);
   cr_assert(SCIPisEQ(scip, 1, constant + slope * (-1))); /* touch at -1, (-1)^(-2) = 1 */
   cr_assert(SCIPisEQ(scip, 0.25, constant + slope * 2)); /* touch at 2, 2^(-2) = 0.25 */

   /* x^(-2) on [-1.0,infty]; underestimator is tangent for xref > 2 */
   success = FALSE;
   islocal = FALSE;
   branchcand = TRUE;
   xref = 4.0;
   estimateHyperbolaPositive(scip, -2.0, root, FALSE, -1.0, SCIPinfinity(scip), xref, -1.0, SCIPinfinity(scip), &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(!islocal); /* the tangent is also globally valid, since global bounds equal local bounds here */
   cr_assert(!branchcand);
   cr_assert(SCIPisEQ(scip, slope, -2.0 * pow(xref, -3.0)));  /* slope should be gradient at xref */
   cr_assert(SCIPisEQ(scip, pow(xref, -2.0), constant + slope * xref)); /* touch at xref */

   /* x^(-2) on [-infty,1.0]; underestimator is secant between -2 and 1 */
   success = TRUE;
   islocal = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(scip, -2.0, root, FALSE, -SCIPinfinity(scip), 1.0, -0.5, -SCIPinfinity(scip), 1.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(islocal);
   cr_assert(branchcand);
   cr_assert(SCIPisEQ(scip, 0.25, constant + slope * (-2))); /* touch at -2, (-2)^(-2) = 0.25 */
   cr_assert(SCIPisEQ(scip, 1, constant + slope * 1)); /* touch at 1, 1^(-2) = 1 */

   success = TRUE;
   branchcand = TRUE;
   estimateHyperbolaPositive(scip, -2.0, SCIP_INVALID, TRUE, -1.0, 1.0, -0.5, -1.0, 1.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(!success); /* overestimator does not exist (or equals infty) */
   cr_assert(branchcand);

   /* x^(-2) on [-2,-1] -> underestimator = tangent, overestimator = secant */
   success = FALSE;
   islocal = TRUE;
   branchcand = TRUE;
   xref = -1.5;
   estimateHyperbolaPositive(scip, -2.0, SCIP_INVALID, FALSE, -2.0, -1.0, xref, -2.0, -1.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(!islocal);
   cr_assert(!branchcand);
   cr_assert(SCIPisEQ(scip, slope, -2.0 * pow(xref, -3.0)));  /* exponent * xref^(exponent-1) */
   cr_assert(SCIPisEQ(scip, constant, pow(xref, -2.0) - slope * xref));

   success = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(scip, -2.0, root, FALSE, -2.0, -1.0, xref, -2.0, 2.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(islocal);  /* if global domain is [-2,2], then the tangent is not globally valid if xref > -xubglobal = -2 */
   cr_assert(!branchcand);  /* but branching will not change the tangent */

   success = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(scip, -2.0, root, FALSE, -2.0, -1.0, xref, -2.0, 0.5, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(!islocal);  /* if global domain is [-2,0.5], then the tangent is globally valid, since xref = -1.5 < xubglobal*root = 0.5*(-2) = -1 */
   cr_assert(!branchcand);

   success = FALSE;
   islocal = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(scip, -2.0, SCIP_INVALID, TRUE, -2.0, -1.0, xref, -2.0, -1.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(islocal);
   cr_assert(branchcand);
   cr_assert(SCIPisEQ(scip, slope, (pow(-1.0, -2.0) - pow(-2.0, -2.0))));
   cr_assert(SCIPisEQ(scip, constant, pow(-2.0, -2.0) - slope * (-2.0)));

   /* x^(-2) on [1, 2] -> underestimator = tangent, overestimator = secant */
   success = FALSE;
   islocal = TRUE;
   branchcand = TRUE;
   xref = 1.5;
   estimateHyperbolaPositive(scip, -2.0, SCIP_INVALID, FALSE, 1.0, 2.0, xref, 1.0, 2.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(!islocal);
   cr_assert(!branchcand);
   cr_assert(SCIPisEQ(scip, slope, -2.0 * pow(xref, -3.0)));  /* exponent * xref^(exponent-1) */
   cr_assert(SCIPisEQ(scip, constant, pow(xref, -2.0) - slope * xref));

   success = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(scip, -2.0, root, FALSE, 1.0, 2.0, xref, -1.0, 2.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(islocal);  /* if global domain is [-1,2], then tangent is not globally valid */
   cr_assert(!branchcand);

   success = FALSE;
   islocal = FALSE;
   branchcand = TRUE;
   estimateHyperbolaPositive(scip, -2.0, SCIP_INVALID, TRUE, 1.0, 2.0, xref, 1.0, 2.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(islocal);
   cr_assert(branchcand);
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
   SCIP_Bool branchcand;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&scip) );

   /* x^(-3) on [-1.0,1.0] */
   success = TRUE;
   branchcand = TRUE;
   estimateHyperbolaMixed(scip, -3.0, FALSE, -1.0, 1.0, -0.5, -1.0, 1.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(!success); /* underestimator does not exist (pole in domain) */
   cr_assert(branchcand);

   success = TRUE;
   branchcand = TRUE;
   estimateHyperbolaMixed(scip, -3.0, TRUE, -1.0, 1.0, -0.5, -1.0, 1.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(!success); /* overestimator does not exist (pole in domain) */
   cr_assert(branchcand);

   /* x^(-3) on [-1.0,0.0] -> underestimator does not exist (upper bound is pole); overestimator is tangent */
   success = TRUE;
   branchcand = TRUE;
   xref = -0.5;
   estimateHyperbolaMixed(scip, -3.0, FALSE, -1.0, 0.0, xref, -1.0, 0.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(!success);
   cr_assert(branchcand);

   success = FALSE;
   islocal = TRUE;
   branchcand = TRUE;
   constant = slope = 5.0;
   estimateHyperbolaMixed(scip, -3.0, TRUE, -1.0, 0.0, xref, -1.0, 0.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(!islocal);
   cr_assert(!branchcand);
   cr_assert(SCIPisEQ(scip, slope, -3.0 * pow(xref, -4.0)));
   cr_assert(SCIPisEQ(scip, constant, pow(xref, -3.0) - slope * xref));

   success = FALSE;
   branchcand = TRUE;
   estimateHyperbolaMixed(scip, -3.0, TRUE, -1.0, 0.0, xref, -1.0, 1.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success);
   cr_assert(islocal);  /* if global domain is [-1,1], then tangent is not globally valid */
   cr_assert(!branchcand);

   /* x^(-3) on [-2.0,-1.0] -> underestimator is secant */
   success = FALSE;
   islocal = FALSE;
   branchcand = TRUE;
   constant = slope = 5.0;
   estimateHyperbolaMixed(scip, -3.0, FALSE, -2.0, -1.0, xref, -2.0, 0.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(success); /* underestimator does not exist (upper bound is pole) */
   cr_assert(islocal);
   cr_assert(branchcand);
   cr_assert(SCIPisEQ(scip, slope, -1.0 - pow(-2.0, -3.0)));
   cr_assert(SCIPisEQ(scip, constant, pow(-2.0, -3.0) - slope * (-2.0)));

   /* x^(-3) on [-infty,-1.0] -> underestimator does not exist */
   success = TRUE;
   branchcand = TRUE;
   estimateHyperbolaMixed(scip, -3.0, FALSE, -SCIPinfinity(scip), -1.0, xref, -SCIPinfinity(scip), -1.0, &constant, &slope, &islocal, &branchcand, &success);
   cr_assert(!success);
   cr_assert(branchcand);

   SCIP_CALL( SCIPfree(&scip) );
}

/* test estimateRoot */
Test(estimation, root, .description = "test computation of estimators for roots (<1)")
{
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Real xref;
   SCIP_Bool islocal;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreate(&scip) );

   /* x^0.25 on [0.0,infty] -> underestimator does not exist */
   success = TRUE;
   estimateRoot(scip, 0.25, FALSE, 0.0, SCIPinfinity(scip), 0.5, &constant, &slope, &islocal, &success);
   cr_assert(!success);

   /* x^0.25 on [0.0,16.0] -> underestimator is secant; overestimator is tangent */
   xref = 4.0;
   success = FALSE;
   islocal = FALSE;
   constant = slope = -5.0;
   estimateRoot(scip, 0.25, FALSE, 0.0, 16.0, xref, &constant, &slope, &islocal, &success);
   cr_assert(success);
   cr_assert(islocal);
   cr_assert(SCIPisEQ(scip, slope, 2.0/16.0));
   cr_assert(SCIPisEQ(scip, constant, 0.0));

   success = FALSE;
   islocal = TRUE;
   estimateRoot(scip, 0.25, TRUE, 0.0, 16.0, xref, &constant, &slope, &islocal, &success);
   cr_assert(success);
   cr_assert(!islocal);
   cr_assert(SCIPisEQ(scip, slope, 0.25 * pow(xref, -0.75)));
   cr_assert(SCIPisEQ(scip, constant, pow(xref, 0.25) - slope * xref));

   /* if reference point at 0.0, then tangent will still be computed, but it will not touch at 0.0 */
   success = FALSE;
   islocal = TRUE;
   estimateRoot(scip, 0.25, TRUE, 0.0, 16.0, 0.0, &constant, &slope, &islocal, &success);
   cr_assert(success);
   cr_assert(!islocal);
   cr_assert(constant != 0.0);

   /* if reference point at 0.0 and bounds on x are very small, then no estimator is computed */
   estimateRoot(scip, 0.25, TRUE, 0.0, SCIPepsilon(scip), 0.0, &constant, &slope, &islocal, &success);
   cr_assert(!success);

   SCIP_CALL( SCIPfree(&scip) );
}

Test(estimation, convexsquare, .init = setup, .fini = teardown,
   .description = "test separation for a convex square expression"
   )
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

   cr_assert(success);
   cr_assert_float_eq(constant, -1.0, SCIPepsilon(scip));
   cr_assert_float_eq(slope, 2.0, SCIPepsilon(scip));
   cr_assert(!islocal);
   cr_assert(!branchcand);

   /*
    * compute overestimator for x^2 with x* = 1.0
    * this should result in a secant estimator
    */
   xval = 1.0;

   branchcand = TRUE;
   SCIP_CALL( estimatePow(scip, expr, &bnd, &bnd, &xval, TRUE, -SCIPinfinity(scip), &slope, &constant, &islocal, &success, &branchcand) );
   cr_assert(success);
   cr_assert_float_eq(constant, 5.0, SCIPepsilon(scip));
   cr_assert_float_eq(slope, 4.0, SCIPepsilon(scip));
   cr_assert(islocal);
   cr_assert(branchcand);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
