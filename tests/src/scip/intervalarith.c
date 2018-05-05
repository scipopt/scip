/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   intervalarith.c
 * @brief  unit tests to check that bug1861 doesn't appear
 * @author Stefan Vigerske
 */

#include "scip/intervalarith.h"
#include "include/scip_test.h"
#include <stdio.h>

Test(intervalarith, issue1861)
{
   SCIP_Real             infinity;
   SCIP_INTERVAL         resultant;
   SCIP_Real             ax;
   SCIP_Real             ay;
   SCIP_Real             axy;
   SCIP_Real             bx;
   SCIP_Real             by;
   SCIP_INTERVAL         rhs;
   SCIP_INTERVAL         xbnds;
   SCIP_INTERVAL         ybnds;

   infinity = 1e43;
   ax = 12;
   ay = 27;
   axy = -36;
   bx = -32;
   by = 48;
   SCIPintervalSetBounds(&rhs, -21.333333336466665741681936197, -21.3333333341333322152877371991);
   SCIPintervalSetBounds(&xbnds, -10.0382009139778674011722614523, 0.0);
   SCIPintervalSetBounds(&ybnds, -7.93110393120117596055251851794, -0.537524704319278456843278490851);

   // the quadratic equation
   // ax*x^2 + ay*y^2 + axy*x*y + bx*x + by*y \in rhs
   // does not have a solution for x \in xbnds and y \in ybnds
   // however, relaxing rhs very slightly gives a solution
   // SCIPintervalSolveBivariateQuadExpressionAllScalar() is expected to either
   // - figure out that there is no solution, or
   // - return compute some bounds because it decided to relax the bounds a bit
   // However, it is not supposed to throw an assert (#1861).

   SCIPintervalSolveBivariateQuadExpressionAllScalar(
      infinity,           /**< value for infinity in interval arithmetics */
      &resultant,          /**< buffer where to store result of operation */
      ax,                 /**< square coefficient of x */
      ay,                 /**< square coefficient of y */
      axy,                /**< bilinear coefficients */
      bx,                 /**< linear coefficient of x */
      by,                 /**< linear coefficient of y */
      rhs,                /**< right-hand-side of equation */
      xbnds,              /**< bounds on x */
      ybnds               /**< bounds on y */
      );
 }


Test(intervalarith, issue2250)
{
   SCIP_Real             infinity;
   SCIP_INTERVAL         resultant;
   SCIP_Real             ax;
   SCIP_Real             ay;
   SCIP_Real             axy;
   SCIP_Real             bx;
   SCIP_Real             by;
   SCIP_INTERVAL         rhs;
   SCIP_INTERVAL         xbnds;
   SCIP_INTERVAL         ybnds;

   infinity = 1e43;
   ax = 1;
   ay = 0;
   axy = 1;
   bx = -2;
   by = -6;
   SCIPintervalSetBounds(&rhs, -infinity, 0.0);
   SCIPintervalSetBounds(&xbnds, 0.0, infinity);
   SCIPintervalSetBounds(&ybnds, 0.0, infinity);

   /* x=y=4 is feasible for this equation, so x=4 should be part of the solution of this equation */
   cr_assert(ax*4*4 + bx*4 + ay*4*4 + by*4 + axy*4*4 <= 1e-12);

   SCIPintervalSolveBivariateQuadExpressionAllScalar(
      infinity,           /**< value for infinity in interval arithmetics */
      &resultant,         /**< buffer where to store result of operation */
      ax,                 /**< square coefficient of x */
      ay,                 /**< square coefficient of y */
      axy,                /**< bilinear coefficients */
      bx,                 /**< linear coefficient of x */
      by,                 /**< linear coefficient of y */
      rhs,                /**< right-hand-side of equation */
      xbnds,              /**< bounds on x */
      ybnds               /**< bounds on y */
      );

   cr_assert(resultant.inf <= 4.0);
   cr_assert(resultant.sup >= 4.0);
}
