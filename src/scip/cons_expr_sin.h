/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_sin.h
 * @brief  handler for sin expressions
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_SIN_H__
#define __SCIP_CONS_EXPR_SIN_H__


#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_sin.h"

#define NEWTON_NITERATIONS    100
#define NEWTON_PRECISION      1e-12

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for sin expressions and includes it into the expression constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeConsExprExprHdlrSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** creates a sin expression */
EXTERN
SCIP_RETCODE SCIPcreateConsExprExprSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPR*   child               /**< single child */
   );

/** evaluates the function a*x + b - sin(x) for some coefficient a and constant b at a given point p
 *  the constants a and b are expected to be stored in that order in params
 */
static
SCIP_DECL_NEWTONEVAL(function1)
{
   assert(params != NULL);
   assert(nparams == 2);

   return params[0]*point + params[1] - SIN(point);
}

/** evaluates the derivative of a*x + b - sin(x) for some coefficient a and constant b at a given point p
 *  the constants a and b are expected to be stored in that order in params
 */
static
SCIP_DECL_NEWTONEVAL(derivative1)
{
   assert(params != NULL);
   assert(nparams == 2);

   return params[0] - COS(point);
}

/** evaluates the function sin(x) + (alpha - x)*cos(x) - sin(alpha) for some constant alpha at a given point p
 *  the constant alpha is expected to be stored in params
 * */
static
SCIP_DECL_NEWTONEVAL(function2)
{
   assert(params != NULL);
   assert(nparams == 1);

   return SIN(point) + (params[0] - point) * COS(point) - SIN(params[0]);
}

/** evaluates the derivative of sin(x) + (alpha - x)*cos(x) - sin(alpha) for some constant alpha at a given point p
 *  the constant alpha is expected to be stored in params
 */
static
SCIP_DECL_NEWTONEVAL(derivative2)
{
   assert(params != NULL);
   assert(nparams == 1);

   return (point - params[0]) * SIN(point);
}

/** helper function to compute the secant if it is a feasible underestimating cut
 *  returns true if the cut was computed successfully
 */
static
SCIP_Bool computeSecantSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Real             lb,                 /**< lower bound of argument variable */
   SCIP_Real             ub                  /**< upper bound of argument variable */
)
{
   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);
   assert(lb < ub);

   /* if range is too big, secant is infeasible */
   if( ub - lb >= M_PI )
      return FALSE;

   /* if bounds are not within positive bay, cut is infeasible or not underestimating */
   if( SIN(lb) < 0.0 || SIN(ub) < 0.0  || (SIN(lb) == 0.0 && COS(lb) < 0.0) )
      return FALSE;

   *lincoef = (SIN(ub) - SIN(lb)) / (ub - lb);
   *linconst = SIN(ub) - (*lincoef) * ub;

   return TRUE;
}

/** helper function to compute the tangent at lower bound if it is a feasible underestimating cut
 *  returns true if the cut was computed successfully
 */
static
SCIP_Bool computeLeftTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Real             lb                  /**< lower bound of argument variable */
)
{
   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);

   if( SCIPisInfinity(scip, lb) )
      return FALSE;

   /* left tangent is only feasible and underestimating in [pi, 1.5*pi) *2kpi */
   if( SIN(lb) > 0.0 || COS(lb) >= 0.0 )
      return FALSE;

   *lincoef = COS(lb);
   *linconst = SIN(lb) - (*lincoef) * lb;

   return TRUE;
}

/** helper function to compute the tangent at upper bound if it is a feasible underestimating cut
 *  returns true if the cut was computed successfully
 */
static
SCIP_Bool computeRightTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Real             ub                  /**< upper bound of argument variable */
)
{
   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);

   if( SCIPisInfinity(scip, ub) )
      return FALSE;

   /* left tangent is only feasible and underestimating in (1.5*pi, 2*pi] *2kpi */
   if( SIN(ub) > 0.0 || COS(ub) <= 0.0 )
      return FALSE;

   *lincoef = COS(ub);
   *linconst = SIN(ub) - (*lincoef) * ub;

   return TRUE;
}

/** helper function to compute the tangent at solution point if it is a feasible underestimating cut
 *  returns true if the cut was computed successfully
 */
static
SCIP_Bool computeSolTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Real             lb,                 /**< lower bound of argument variable */
   SCIP_Real             ub,                 /**< upper bound of argument variable */
   SCIP_Real             solpoint            /**< solution point to be separated */
)
{
   SCIP_Real params[2];
   SCIP_Real startingpoints[3];
   SCIP_Real solpointmodpi;
   SCIP_Real intersection;
   int i;

   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);

   /* tangent is only underestimating in negative bay */
   if( SIN(solpoint) > 0.0 )
      return FALSE;

   /* compute solution point mod pi */
   solpointmodpi = fmod(solpoint, M_PI);
   if( solpoint < 0.0 )
      solpointmodpi += M_PI;


   /* if the point is too far away from the bounds or is at a multiple of pi, the tangent is infeasible */
   if( SCIPisGE(scip, solpoint - lb, 2*M_PI) || SCIPisGE(scip, ub - solpoint, 2*M_PI)
      || SCIPisZero(scip, solpointmodpi) )
      return FALSE;

   params[0] = COS(solpoint);
   params[1] = SIN(solpoint) - params[0] * solpoint;

   /* choose starting points for newton procedure */
   if( SCIPisGT(scip, solpointmodpi, 0.5*M_PI) )
   {
      startingpoints[0] = solpoint + (M_PI - solpointmodpi) + 0.5*M_PI;
      startingpoints[1] = startingpoints[0] + 0.5*M_PI;
      startingpoints[2] = startingpoints[1] + 0.5*M_PI;
   }
   else
   {
      startingpoints[0] = solpoint - solpointmodpi - 0.5*M_PI;
      startingpoints[1] = startingpoints[0] - 0.5*M_PI;
      startingpoints[2] = startingpoints[1] - 0.5*M_PI;
   }

   /* use newton procedure to test if cut is valid */
   for( i = 0; i < 3; ++i)
   {
      intersection = SCIPcomputeRootNewton(function1, derivative1, params, 2, startingpoints[i], NEWTON_PRECISION,
         NEWTON_NITERATIONS);

      if( intersection != SCIP_INVALID && !SCIPisEQ(scip, intersection, solpoint) ) /*lint !e777*/
         break;
   }

   /* if newton failed or intersection point lies within bounds, cut is not feasible */
   if( intersection == SCIP_INVALID || (intersection >= lb && intersection <= ub) ) /*lint !e777*/
      return FALSE;

   *lincoef = params[0];
   *linconst = params[1];

   return TRUE;
}

/** helper function to compute the tangent at some other point that goes through (lb,sin(lb)) and is underestimating
 *  returns true if the cut was computed successfully
 */
static
SCIP_Bool computeLeftMidTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Bool*            issecant,           /**< buffer to store whether cut is actually a secant */
   SCIP_Real             lb,                 /**< lower bound of argument variable */
   SCIP_Real             ub                  /**< upper bound of argument variable */
)
{

   SCIP_Real lbmodpi;
   SCIP_Real tangentpoint;
   SCIP_Real startingpoint;

   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);
   assert(lb < ub);

   *issecant = FALSE;

   if( SCIPisInfinity(scip, lb) )
      return FALSE;

   /* compute shifted bounds for case evaluation */
   lbmodpi = fmod(lb, M_PI);
   if( lb < 0.0 )
      lbmodpi += M_PI;

   /* choose starting point for newton procedure */
   if( COS(lb) < 0.0 )
   {
      /* in [pi/2,pi] underestimating doesn't work; otherwise, take the midpoint of possible area */
      if( SIN(lb) <= 0.0 )
         return FALSE;
      else
         startingpoint = lb + 1.25*M_PI - lbmodpi;
   }
   else
   {
      /* in ascending area, take the midpoint of the possible area in descending part */
      if( SIN(lb) <= 0.0 )
         startingpoint = lb + 2.25*M_PI - lbmodpi;
      else
         startingpoint = lb + 1.25*M_PI - lbmodpi;
   }

   /* use newton procedure to find the point where the tangent intersects sine at lower bound */
   tangentpoint = SCIPcomputeRootNewton(function2, derivative2, &lb, 1, startingpoint, NEWTON_PRECISION,
      NEWTON_NITERATIONS);

   /* if newton procedure failed, no cut is added */
   if( tangentpoint == SCIP_INVALID ) /*lint !e777*/
      return FALSE;

   /* if the computed point lies outside the bounds, it is shifted to upper bound */
   if( SCIPisGE(scip, tangentpoint, ub) )
   {
      tangentpoint = ub;

      /* check whether cut is still underestimating */
      if( SIN(0.5 * (ub + lb)) < SIN(lb) + 0.5*(SIN(ub) - SIN(lb)) / (ub - lb) )
         return FALSE;

      *issecant = TRUE;
   }

   /* compute secant between lower bound and connection point */
   *lincoef = (SIN(tangentpoint) - SIN(lb)) / (tangentpoint - lb);
   *linconst = SIN(lb) - (*lincoef) * lb;

   return TRUE;
}

/** helper function to compute the tangent at some other point that goes through (ub,sin(ub)) and is underestimating
 *  returns true if the cut was computed successfully
 */
static
SCIP_Bool computeRightMidTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Bool*            issecant,           /**< buffer to stroe whether cut is actually a secant */
   SCIP_Real             lb,                 /**< lower bound of argument variable */
   SCIP_Real             ub                  /**< upper bound of argument variable */
)
{

   SCIP_Real ubmodpi;
   SCIP_Real tangentpoint;
   SCIP_Real startingpoint;

   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);
   assert(lb < ub);

   *issecant = FALSE;

   if( SCIPisInfinity(scip, ub) )
      return FALSE;

   /* compute shifted bounds for case evaluation */
   ubmodpi = fmod(ub, M_PI);
   if( ub < 0.0 )
      ubmodpi += M_PI;

   /* choose starting point for newton procedure */
   if( COS(ub) > 0.0 )
   {
      /* in [pi/2,pi] underestimating doesn't work; othherwise, take the midpoint of possible area */
      if( SIN(ub) <= 0.0 )
         return FALSE;
      else
         startingpoint = ub - 0.25*M_PI - ubmodpi;
   }
   else
   {
      /* in ascending area, take the midpoint of the possible area in descending part */
      if( SIN(ub) <= 0.0 )
         startingpoint = ub - 1.25*M_PI - ubmodpi;
      else
         startingpoint = ub - 0.25*M_PI - ubmodpi;
   }

   /* use newton procedure to find the point where the tangent intersects sine at lower bound */
   tangentpoint = SCIPcomputeRootNewton(function2, derivative2, &ub, 1, startingpoint, NEWTON_PRECISION,
      NEWTON_NITERATIONS);

   /* if newton procedure failed, no cut is added */
   if( tangentpoint == SCIP_INVALID ) /*lint !e777*/
      return FALSE;

   /* if the computed point lies outside the bounds, it is shifted to upper bound */
   if( SCIPisLE(scip, tangentpoint, lb) )
   {
      tangentpoint = lb;

      /* check whether cut is still underestimating */
      if( SIN(0.5 * (ub + lb)) < SIN(lb) + 0.5*(SIN(ub) - SIN(lb)) / (ub - lb) )
         return FALSE;

      *issecant = TRUE;
   }

   /* compute secant between lower bound and connection point */
   *lincoef = (SIN(tangentpoint) - SIN(ub)) / (tangentpoint - ub);
   *linconst = SIN(ub) - (*lincoef) * ub;

   return TRUE;
}

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_SIN_H__ */
