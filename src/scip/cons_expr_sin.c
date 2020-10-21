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

/**@file   cons_expr_sin.c
 * @brief  handler for sine expressions
 * @author Fabian Wegscheider
 *
 * The estimator/separator code always computes underestimators for sin(x).
 * For overestimator or cos(x), we first reduce to underestimators of sin(x).
 *
 * Overestimator for sin(x):
 *   Assume that a*y+b <= sin(y) for y in [-ub,-lb].
 *   Then we have a*(-y)-b >= -sin(y) = sin(-y) for y in [-ub,-lb].
 *   Thus, a*x-b >= sin(x) for x in [lb,ub].
 *
 * Underestimator for cos(x):
 *   Assume that a*y+b <= sin(y) for y in [lb+pi/2,ub+pi/2].
 *   Then we have a*(x+pi/2) + b <= sin(x+pi/2) = cos(x) for x in [lb,ub].
 *   Thus, a*x + (b+a*pi/2) <= cos(x) for x in [lb,ub].
 *
 * Overestimator for cos(x):
 *   Assume that a*z+b <= sin(z) for z in [-(ub+pi/2),-(lb+pi/2)].
 *   Then, a*y-b >= sin(y) for y in [lb+pi/2,ub+pi/2].
 *   Then, a*x-b+a*pi/2 >= cos(x) for x in [lb,ub].
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define _USE_MATH_DEFINES   /* to get M_PI on Windows */  /*lint !750 */

#include <string.h>
#include <math.h>
#include "scip/cons_expr_sin.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_rowprep.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2         1.57079632679489661923
#endif
#ifndef M_PI_4
#define M_PI_4         0.785398163397448309616
#endif

/* fundamental expression handler properties */
#define EXPRHDLR_NAME         "sin"
#define EXPRHDLR_DESC         "sine expression"
#define EXPRHDLR_PRECEDENCE   91000
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(82457.0)

#define MAXCHILDABSVAL        1e+6                       /**< maximum absolute value that is accepted for propagation */

/*
 * Data structures
 */

/*
 * Local methods
 */

/** evaluates the function a*x + b - sin(x) for some coefficient a and constant b at a given point p
 *  the constants a and b are expected to be stored in that order in params
 */
static
SCIP_DECL_NEWTONEVAL(function1)
{  /*lint --e{715}*/
   assert(params != NULL);
   assert(nparams == 2);

   return params[0]*point + params[1] - SIN(point);
}

/** evaluates the derivative of a*x + b - sin(x) for some coefficient a and constant b at a given point p
 *  the constants a and b are expected to be stored in that order in params
 */
static
SCIP_DECL_NEWTONEVAL(derivative1)
{  /*lint --e{715}*/
   assert(params != NULL);
   assert(nparams == 2);

   return params[0] - COS(point);
}

/** evaluates the function sin(x) + (alpha - x)*cos(x) - sin(alpha) for some constant alpha at a given point p
 *  the constant alpha is expected to be stored in params
 */
static
SCIP_DECL_NEWTONEVAL(function2)
{  /*lint --e{715}*/
   assert(params != NULL);
   assert(nparams == 1);

   return SIN(point) + (params[0] - point) * COS(point) - SIN(params[0]);
}

/** evaluates the derivative of sin(x) + (alpha - x)*cos(x) - sin(alpha) for some constant alpha at a given point p
 *  the constant alpha is expected to be stored in params
 */
static
SCIP_DECL_NEWTONEVAL(derivative2)
{  /*lint --e{715}*/
   assert(params != NULL);
   assert(nparams == 1);

   return (point - params[0]) * SIN(point);
}

/** helper function to compute the secant if it is a valid underestimator
 *  returns true if the estimator was computed successfully
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

   /* if range is too big, secant is not underestimating */
   if( ub - lb >= M_PI )
      return FALSE;

   /* if bounds are not within positive bay, secant is not underestimating */
   if( SIN(lb) < 0.0 || SIN(ub) < 0.0  || (SIN(lb) == 0.0 && COS(lb) < 0.0) )
      return FALSE;

   *lincoef = (SIN(ub) - SIN(lb)) / (ub - lb);
   *linconst = SIN(ub) - (*lincoef) * ub;

   return TRUE;
}

/** helper function to compute the tangent at lower bound if it is underestimating
 *  returns true if the underestimator was computed successfully
 */
static
SCIP_Bool computeLeftTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of tangent */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of tangent */
   SCIP_Real             lb                  /**< lower bound of argument variable */
   )
{
   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);

   if( SCIPisInfinity(scip, -lb) )
      return FALSE;

   /* left tangent is only underestimating in [pi, 1.5*pi) *2kpi */
   if( SIN(lb) > 0.0 || COS(lb) >= 0.0 )
      return FALSE;

   *lincoef = COS(lb);
   *linconst = SIN(lb) - (*lincoef) * lb;

   return TRUE;
}

/** helper function to compute the tangent at upper bound if it is an underestimator
 *  returns true if the underestimator was computed successfully
 */
static
SCIP_Bool computeRightTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of tangent */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of tangent */
   SCIP_Real             ub                  /**< upper bound of argument variable */
   )
{
   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);

   if( SCIPisInfinity(scip, ub) )
      return FALSE;

   /* left tangent is only underestimating in (1.5*pi, 2*pi] *2kpi */
   if( SIN(ub) > 0.0 || COS(ub) <= 0.0 )
      return FALSE;

   *lincoef = COS(ub);
   *linconst = SIN(ub) - (*lincoef) * ub;

   return TRUE;
}

/** helper function to compute the tangent at solution point if it is an underestimator
 *  returns true if the underestimator was computed successfully
 */
static
SCIP_Bool computeSolTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of tangent */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of tangent */
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

   /* if the point is too far away from the bounds or is at a multiple of pi, then tangent is not underestimating */
   if( SCIPisGE(scip, solpoint - lb, 2*M_PI) || SCIPisGE(scip, ub - solpoint, 2*M_PI)
      || SCIPisZero(scip, solpointmodpi) )
      return FALSE;

   params[0] = COS(solpoint);
   params[1] = SIN(solpoint) - params[0] * solpoint;

   /* choose starting points for Newton procedure */
   if( SCIPisGT(scip, solpointmodpi, M_PI_2) )
   {
      startingpoints[0] = solpoint + (M_PI - solpointmodpi) + M_PI_2;
      startingpoints[1] = startingpoints[0] + M_PI_2;
      startingpoints[2] = startingpoints[1] + M_PI_2;
   }
   else
   {
      startingpoints[0] = solpoint - solpointmodpi - M_PI_2;
      startingpoints[1] = startingpoints[0] - M_PI_2;
      startingpoints[2] = startingpoints[1] - M_PI_2;
   }

   /* use Newton procedure to test if cut is valid */
   for( i = 0; i < 3; ++i )
   {
      intersection = SCIPcomputeRootNewton(function1, derivative1, params, 2, startingpoints[i], NEWTON_PRECISION,
         NEWTON_NITERATIONS);

      if( intersection != SCIP_INVALID && !SCIPisEQ(scip, intersection, solpoint) ) /*lint !e777*/
         break;
   }

   /* if Newton failed or intersection point lies within bounds, underestimator is not valid */
   if( intersection == SCIP_INVALID || (intersection >= lb && intersection <= ub) ) /*lint !e777*/
      return FALSE;

   *lincoef = params[0];
   *linconst = params[1];

   return TRUE;
}

/** helper function to compute the tangent at some other point that goes through (lb,sin(lb)) and is underestimating
 *  returns true if the underestimator was computed successfully
 */
static
SCIP_Bool computeLeftMidTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of tangent */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of tangent */
   SCIP_Bool*            issecant,           /**< buffer to store whether underestimator is actually a secant */
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

   if( SCIPisInfinity(scip, -lb) )
      return FALSE;

   /* compute shifted bounds for case evaluation */
   lbmodpi = fmod(lb, M_PI);
   if( lb < 0.0 )
      lbmodpi += M_PI;

   /* choose starting point for Newton procedure */
   if( COS(lb) < 0.0 )
   {
      /* in [pi/2,pi] underestimating doesn't work; otherwise, take the midpoint of possible area */
      if( SCIPisLE(scip, SIN(lb), 0.0) )
         return FALSE;
      else
         startingpoint = lb + 1.25*M_PI - lbmodpi;
   }
   else
   {
      /* in ascending area, take the midpoint of the possible area in descending part */
      if( SCIPisLT(scip, SIN(lb), 0.0) )
         startingpoint = lb + 2.25*M_PI - lbmodpi;
      else
         startingpoint = lb + 1.25*M_PI - lbmodpi;
   }

   /* use Newton procedure to find the point where the tangent intersects sine at lower bound */
   tangentpoint = SCIPcomputeRootNewton(function2, derivative2, &lb, 1, startingpoint, NEWTON_PRECISION,
      NEWTON_NITERATIONS);

   /* if Newton procedure failed, no cut is added */
   if( tangentpoint == SCIP_INVALID ) /*lint !e777*/
      return FALSE;

   /* if the computed point lies outside the bounds, it is shifted to upper bound */
   if( SCIPisGE(scip, tangentpoint, ub) )
   {
      tangentpoint = ub;

      /* check whether affine function is still underestimating */
      if( SCIPisLE(scip, SIN(0.5 * (ub + lb)), SIN(lb) + 0.5*(SIN(ub) - SIN(lb))) )
         return FALSE;

      *issecant = TRUE;
   }

   if( SCIPisEQ(scip, tangentpoint, lb) )  /*lint !e777 */
      return FALSE;

   /* compute secant between lower bound and connection point */
   *lincoef = (SIN(tangentpoint) - SIN(lb)) / (tangentpoint - lb);
   *linconst = SIN(lb) - (*lincoef) * lb;

   /* if the bounds are to close to each other, it's possible that the underestimator is not valid */
   if( *lincoef >= COS(lb) )
      return FALSE;

   SCIPdebugMsg(scip, "leftmidtangent: %g + %g*x <= sin(x) on [%g,%g]\n", *linconst, *lincoef, lb, ub);

   return TRUE;
}

/** helper function to compute the tangent at some other point that goes through (ub,sin(ub)) and is underestimating
 *  returns true if the underestimator was computed successfully
 */
static
SCIP_Bool computeRightMidTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of tangent */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of tangent */
   SCIP_Bool*            issecant,           /**< buffer to store whether underestimator is actually a secant */
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

   /* choose starting point for Newton procedure */
   if( COS(ub) > 0.0 )
   {
      /* in [3*pi/2,2*pi] underestimating doesn't work; otherwise, take the midpoint of possible area */
      if( SCIPisLE(scip, SIN(ub), 0.0) )
         return FALSE;
      else
         startingpoint = ub - M_PI_4 - ubmodpi;
   }
   else
   {
      /* in descending area, take the midpoint of the possible area in ascending part */
      if( SCIPisLE(scip, SIN(ub), 0.0) )
         startingpoint = ub - 1.25*M_PI - ubmodpi;
      else
         startingpoint = ub - M_PI_4 - ubmodpi;
   }

   /* use Newton procedure to find the point where the tangent intersects sine at lower bound */
   tangentpoint = SCIPcomputeRootNewton(function2, derivative2, &ub, 1, startingpoint, NEWTON_PRECISION,
      NEWTON_NITERATIONS);

   /* if Newton procedure failed, no underestimator is found */
   if( tangentpoint == SCIP_INVALID ) /*lint !e777*/
      return FALSE;

   /* if the computed point lies outside the bounds, it is shifted to upper bound */
   if( SCIPisLE(scip, tangentpoint, lb) )
   {
      tangentpoint = lb;

      /* check whether affine function is still underestimating */
      if( SCIPisLE(scip, SIN(0.5 * (ub + lb)), SIN(lb) + 0.5*(SIN(ub) - SIN(lb))) )
         return FALSE;

      *issecant = TRUE;
   }

   if( SCIPisEQ(scip, tangentpoint, ub) )  /*lint !e777 */
      return FALSE;

   /* compute secant between lower bound and connection point */
   *lincoef = (SIN(tangentpoint) - SIN(ub)) / (tangentpoint - ub);
   *linconst = SIN(ub) - (*lincoef) * ub;

   /* if the bounds are to close to each other, it's possible that the underestimator is not valid */
   if( *lincoef <= COS(lb) )
      return FALSE;

   return TRUE;
}

/** sets up a rowprep from given data */
static
SCIP_RETCODE assembleRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        rowprep,            /**< buffer to store rowprep */
   const char*           name,               /**< name of type of cut */
   SCIP_Bool             iscos,              /**< whether we are doing a cut for cosine instead of sine */
   SCIP_Bool             underestimate,      /**< whether underestimating */
   SCIP_Real             linconst,           /**< constant term */
   SCIP_Real             lincoef,            /**< coefficient of childvar */
   SCIP_VAR*             childvar,           /**< child var */
   SCIP_VAR*             auxvar              /**< auxiliary variable */
   )
{
   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(childvar != NULL);
   assert(auxvar != NULL);

   /* for overestimators, mirror back */
   if( !underestimate )
      linconst *= -1.0;

   /* further, for cos expressions, the estimator needs to be shifted back to match original bounds */
   if( iscos )
      linconst += lincoef * M_PI_2;

   SCIP_CALL( SCIPcreateRowprep(scip, rowprep, underestimate ? SCIP_SIDETYPE_RIGHT : SCIP_SIDETYPE_LEFT, TRUE) );
   (void) SCIPsnprintf((*rowprep)->name, SCIP_MAXSTRLEN, "%s_%s_%s_%lld", iscos ? "cos" : "sin", name,
      SCIPvarGetName(childvar), SCIPgetNLPs(scip));

   SCIPaddRowprepConstant(*rowprep, linconst);

   SCIP_CALL( SCIPensureRowprepSize(scip, *rowprep, 2) );
   SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, auxvar, -1.0) );
   SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, childvar, lincoef) );

   return SCIP_OKAY;
}


/** helper function to compute the new interval for child in reverse propagation */
SCIP_RETCODE SCIPcomputeRevPropIntervalSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTERVAL         parentbounds,       /**< bounds for sine expression */
   SCIP_INTERVAL         childbounds,        /**< bounds for child expression */
   SCIP_INTERVAL*        newbounds           /**< buffer to store new child bounds */
   )
{
   SCIP_Real newinf = childbounds.inf;
   SCIP_Real newsup = childbounds.sup;

   /* if the absolute values of the bounds are too large, skip reverse propagation
    * TODO: if bounds are close but too large, shift them to [0,2pi] and do the computation there
    */
   if( ABS(newinf) > MAXCHILDABSVAL || ABS(newsup) > MAXCHILDABSVAL )
   {
      SCIPintervalSetBounds(newbounds, newinf, newsup);
      return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip, -newinf) )
   {
      /* l(x) and u(x) are lower/upper bound of child, l(s) and u(s) are lower/upper bound of sin expr
       *
       * if sin(l(x)) < l(s), we are looking for k minimal s.t. a + 2k*pi > l(x) where a = asin(l(s))
       * then the new lower bound is a + 2k*pi
       */
      if( SCIPisLT(scip, SIN(newinf), parentbounds.inf) )
      {
         SCIP_Real a = ASIN(parentbounds.inf);
         int k = (int) ceil((newinf - a) / (2.0*M_PI));
         newinf = a + 2.0*M_PI * k;
      }

      /* if sin(l(x)) > u(s), we are looking for k minimal s.t. pi - a + 2k*pi > l(x) where a = asin(u(s))
       * then the new lower bound is pi - a + 2k*pi
       */
      else if( SCIPisGT(scip, SIN(newinf), parentbounds.sup) )
      {
         SCIP_Real a = ASIN(parentbounds.sup);
         int k = (int) ceil((newinf + a) / (2.0*M_PI) - 0.5);
         newinf = M_PI * (2.0*k + 1.0) - a;
      }

      assert(newinf >= childbounds.inf);
      assert(SCIPisFeasGE(scip, SIN(newinf), parentbounds.inf));
      assert(SCIPisFeasLE(scip, SIN(newinf), parentbounds.sup));
   }

   if( !SCIPisInfinity(scip, newsup) )
   {
      /* if sin(u(x)) > u(s), we are looking for k minimal s.t. a + 2k*pi > u(x) - 2*pi where a = asin(u(s))
       * then the new upper bound is a + 2k*pi
       */
      if ( SCIPisGT(scip, SIN(newsup), parentbounds.sup) )
      {
         SCIP_Real a = ASIN(parentbounds.sup);
         int k = (int) ceil((newsup - a ) / (2.0*M_PI)) - 1;
         newsup = a + 2.0*M_PI * k;
      }

      /* if sin(u(x)) < l(s), we are looking for k minimal s.t. pi - a + 2k*pi > l(x) - 2*pi where a = asin(l(s))
       * then the new upper bound is pi - a + 2k*pi
       */
      if( SCIPisLT(scip, SIN(newsup), parentbounds.inf) )
      {
         SCIP_Real a = ASIN(parentbounds.inf);
         int k = (int) ceil((newsup + a) / (2.0*M_PI) - 0.5) - 1;
         newsup = M_PI * (2.0*k + 1.0) - a;
      }

      assert(newsup <= childbounds.sup);
      assert(SCIPisFeasGE(scip, SIN(newsup), parentbounds.inf));
      assert(SCIPisFeasLE(scip, SIN(newsup), parentbounds.sup));
   }

   /* if the new interval is invalid, the old one was already invalid */
   if( newinf <= newsup )
      SCIPintervalSetBounds(newbounds, newinf, newsup);
   else
      SCIPintervalSetEmpty(newbounds);

   return SCIP_OKAY;
}

/** helper function to compute coefficients and constant term of a linear estimator at a given point
 *
 *  The function will try to compute the following estimators in that order:
 *  - soltangent: tangent at specified refpoint
 *  - secant: secant between the points (lb,sin(lb)) and (ub,sin(ub))
 *  - lmidtangent: tangent at some other point that goes through (lb,sin(lb))
 *  - rmidtangent: tangent at some other point that goes through (ub,sin(ub))
 *
 *  They are ordered such that a successful computation for one of them cannot be improved by following ones in terms
 *  of value at the reference point
 */
SCIP_Bool SCIPcomputeEstimatorsTrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sin or cos expression */
   SCIP_Real*            lincoef,            /**< buffer to store the linear coefficient */
   SCIP_Real*            linconst,           /**< buffer to store the constant term */
   SCIP_Real             refpoint,           /**< point at which to underestimate (can be SCIP_INVALID) */
   SCIP_Real             childlb,            /**< lower bound of child variable */
   SCIP_Real             childub,            /**< upper bound of child variable */
   SCIP_Bool             underestimate       /**< whether the estimator should be underestimating */
   )
{
   SCIP_Bool success;
   SCIP_Bool iscos;
   SCIP_Bool issecant;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "sin") == 0 || strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "cos") == 0);
   assert(SCIPisLE(scip, childlb, childub));

   /* if child is essentially constant, then there should be no point in estimation */
   if( SCIPisEQ(scip, childlb, childub) ) /* @todo maybe return a constant estimator? */
      return FALSE;

   iscos = strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "cos") == 0;

   /* for cos expressions, the bounds have to be shifted before and after computation */
   if( iscos )
   {
      childlb += M_PI_2;
      childub += M_PI_2;
      refpoint += M_PI_2;
   }

   if( !underestimate )
   {
      SCIP_Real tmp = childlb;
      childlb = -childub;
      childub = -tmp;
      refpoint *= -1;
   }

   /* try out tangent at solution point */
   success = computeSolTangentSin(scip, lincoef, linconst, childlb, childub, refpoint);

   /* otherwise, try out secant */
   if( !success )
      success = computeSecantSin(scip, lincoef, linconst, childlb, childub);

   /* otherwise, try left middle tangent, that is tangent at some other point which goes through (lb,sin(lb)) */
   if( !success )
      success = computeLeftMidTangentSin(scip, lincoef, linconst, &issecant, childlb, childub);

   /* otherwise, try right middle tangent, that is tangent at some other point which goes through (ub,sin(ub)) */
   if( !success )
      success = computeRightMidTangentSin(scip, lincoef, linconst, &issecant, childlb, childub);

   if( !success )
      return FALSE;

   /* for overestimators, mirror back */
   if( !underestimate )
      (*linconst) *= -1.0;

   /* for cos expressions, shift back */
   if( iscos )
      (*linconst) += (*lincoef) * M_PI_2;

   return TRUE;
}

/** helper function to create initial cuts for sine and cosine separation
 *
 *  The following 5 cuts can be generated:
 *  - secant: secant between the points (lb,sin(lb)) and (ub,sin(ub))
 *  - ltangent/rtangent: tangents at the points (lb,sin(lb)) or (ub,sin(ub))
 *  - lmidtangent/rmidtangent: tangent at some other point that goes through (lb,sin(lb)) or (ub,sin(ub))
 *
 *  If one of the computations fails or turns out to be irrelevant, the respective argument pointer is set to NULL.
 */
SCIP_RETCODE SCIPcomputeInitialCutsTrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sin or cos expression */
   SCIP_ROWPREP**        secant,             /**< pointer to store the secant */
   SCIP_ROWPREP**        ltangent,           /**< pointer to store the left tangent */
   SCIP_ROWPREP**        rtangent,           /**< pointer to store the right tangent */
   SCIP_ROWPREP**        lmidtangent,        /**< pointer to store the left middle tangent */
   SCIP_ROWPREP**        rmidtangent,        /**< pointer to store the right middle tangent */
   SCIP_Real             childlb,            /**< lower bound of child variable */
   SCIP_Real             childub,            /**< upper bound of child variable */
   SCIP_Bool             underestimate       /**< whether the cuts should be underestimating */
   )
{
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* auxvar;
   SCIP_VAR* childvar;
   SCIP_Real lincoef;
   SCIP_Real linconst;
   SCIP_Bool success;
   SCIP_Bool iscos;
   SCIP_Bool issecant;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "sin") == 0 || strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "cos") == 0);
   assert(SCIPisLE(scip, childlb, childub));

   assert(secant != NULL);
   assert(ltangent != NULL);
   assert(rtangent != NULL);
   assert(lmidtangent != NULL);
   assert(rmidtangent != NULL);

   /* caller must ensure that variable is not already fixed */
   assert(!SCIPisEQ(scip, childlb, childub));

   /* get expression data */
   auxvar = SCIPgetConsExprExprAuxVar(expr);
   assert(auxvar != NULL);
   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   childvar = SCIPgetConsExprExprAuxVar(child);
   assert(childvar != NULL);

   iscos = strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "cos") == 0;

   /* for cos expressions, the bounds have to be shifted before and after computation */
   if( iscos )
   {
      childlb += M_PI_2;
      childub += M_PI_2;
   }

   /*
    * Compute all initial cuts
    * For each linear equation z = a*x + b with bounds [lb,ub] the parameters can be computed by:
    *
    * a = cos(x^)    and     b = sin(x^) - a * x^        where x^ is any known point in [lb,ub]
    *
    * and the resulting cut is       a*x - z <=/>= -b           depending on over-/underestimation
    */

   /* compute secant between lower and upper bound */
   *secant = NULL;

   if( underestimate )
      success = computeSecantSin(scip, &lincoef, &linconst, childlb, childub);
   else
      success = computeSecantSin(scip, &lincoef, &linconst, -childub, -childlb);

   if( success )
   {
      SCIP_CALL( assembleRowprep(scip, secant, "secant", iscos, underestimate, linconst, lincoef, childvar, auxvar) );
   }

   /* compute tangent at lower bound */
   *ltangent = NULL;

   if( underestimate )
      success = computeLeftTangentSin(scip, &lincoef, &linconst, childlb);
   else
      success = computeRightTangentSin(scip, &lincoef, &linconst, -childlb);

   if( success )
   {
      SCIP_CALL( assembleRowprep(scip, ltangent, "ltangent", iscos, underestimate, linconst, lincoef, childvar, auxvar) );
   }

   /* compute tangent at upper bound */
   *rtangent = NULL;

   if( underestimate )
      success = computeRightTangentSin(scip, &lincoef, &linconst, childub);
   else
      success = computeLeftTangentSin(scip, &lincoef, &linconst, -childub);

   if( success )
   {
      SCIP_CALL( assembleRowprep(scip, rtangent, "rtangent", iscos, underestimate, linconst, lincoef, childvar, auxvar) );
   }

   /* compute left middle tangent, that is tangent at some other point which goes through (lb,sin(lb))
    * if secant is feasible, this cut can never beat it so don't compute it
    */
   *lmidtangent = NULL;

   if( *secant == NULL )
   {
      if( underestimate )
         success = computeLeftMidTangentSin(scip, &lincoef, &linconst, &issecant, childlb, childub);
      else
         success = computeRightMidTangentSin(scip, &lincoef, &linconst, &issecant, -childub, -childlb);

      if( success )
      {
         /* if the cut connects bounds, it is stored in secant */
         SCIP_CALL( assembleRowprep(scip, issecant ? secant : lmidtangent, "lmidtangent", iscos, underestimate, linconst, lincoef, childvar, auxvar) );
      }
   }

   /* compute right middle tangent, that is tangent at some other point which goes through (ub,sin(ub))
    * if secant or soltangent are feasible, this cut can never beat them
    */
   *rmidtangent = NULL;

   if( *secant == NULL )
   {
      if( underestimate )
         success = computeRightMidTangentSin(scip, &lincoef, &linconst, &issecant, childlb, childub);
      else
         success = computeLeftMidTangentSin(scip, &lincoef, &linconst, &issecant, -childub, -childlb);

      if( success )
      {
         /* if the cut connects bounds, it is stored in secant */
         SCIP_CALL( assembleRowprep(scip, issecant ? secant : rmidtangent, "rmidtangent", iscos, underestimate, linconst, lincoef, childvar, auxvar) );
      }
   }

   return SCIP_OKAY;
}

/* helper function that computes the curvature of a sine expression for given bounds and curvature of child */
SCIP_EXPRCURV SCIPcomputeCurvatureSin(
   SCIP_EXPRCURV         childcurvature,     /**< curvature of child */
   SCIP_Real             lb,                 /**< lower bound of child */
   SCIP_Real             ub                  /**< upper bound of child */
   )
{
   SCIP_Real lbsin = SIN(lb);
   SCIP_Real ubsin = SIN(ub);
   SCIP_Real lbcos = COS(lb);
   SCIP_Real ubcos = COS(ub);

   /* curvature can only be determined if bounds lie within one bay*/
   if( (ub - lb <= M_PI) && (lbsin * ubsin >= 0.0) )
   {
      /* special case that both sin(ub) and sin(lb) are 0 (i.e. ub - lb = pi) */
      if( lbsin == 0.0 && ubsin == 0.0 )
      {
         if( childcurvature == SCIP_EXPRCURV_LINEAR )
            return (fmod(lb, 2.0*M_PI) == 0.0) ? SCIP_EXPRCURV_CONCAVE : SCIP_EXPRCURV_CONVEX;
      }

      /* if sine is monotone on the interval, the curvature depends on the child curvature and on the segment */
      else if( lbcos * ubcos >= 0.0 )
      {
         /* on [0, pi/2], sine is concave iff child is concave */
         if( lbsin >= 0.0 && lbcos >= 0.0 && ((int)(childcurvature & SCIP_EXPRCURV_CONCAVE) != 0))
            return SCIP_EXPRCURV_CONCAVE;

         /* on [pi/2, pi], sine is concave iff child is convex */
         if( lbsin >= 0.0 && lbcos <= 0.0 && ((int)(childcurvature & SCIP_EXPRCURV_CONVEX) != 0))
            return SCIP_EXPRCURV_CONCAVE;

         /* on [pi, 3pi/2], sine is convex iff child is concave */
         if( lbsin <= 0.0 && lbcos <= 0.0 && ((int)(childcurvature & SCIP_EXPRCURV_CONCAVE) != 0))
            return SCIP_EXPRCURV_CONVEX;

         /* on [3pi/2, 2pi], sine is convex iff child is convex */
         if( lbsin <= 0.0 && lbcos >= 0.0 && ((int)(childcurvature & SCIP_EXPRCURV_CONVEX) != 0))
            return SCIP_EXPRCURV_CONVEX;
      }

      /* otherwise, we can only say something if the child is linear */
      else if( childcurvature == SCIP_EXPRCURV_LINEAR )
         return (lbsin >= 0.0 && ubsin >= 0.0) ? SCIP_EXPRCURV_CONCAVE : SCIP_EXPRCURV_CONVEX;
   }

   return SCIP_EXPRCURV_UNKNOWN;
}

/*
 * Callback methods of expression handler
 */

/** expression handler copy callback */
static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrSin)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrSin(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

/** simplifies a sin expression
 * Evaluates the sine value function when its child is a value expression
 * TODO: add further simplifications
 */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifySin)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);

   /* check for value expression */
   if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrValue(conshdlr) )
   {
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr,
            SIN(SCIPgetConsExprExprValueValue(child))) );
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureConsExprExpr(*simplifiedexpr);
   }

   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_CONSEXPR_EXPRPARSE(parseSin)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseConsExprExpr(scip, consexprhdlr, string, endstring, &childexpr) );
   assert(childexpr != NULL);

   /* create sine expression */
   SCIP_CALL( SCIPcreateConsExprExprSin(scip, consexprhdlr, expr, childexpr) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the sine expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPREVAL(evalSin)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   *val = SIN(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffSin)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;

   assert(expr != NULL);
   assert(childidx == 0);
   assert(SCIPgetConsExprExprValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child)), "val") != 0);

   *val = COS(SCIPgetConsExprExprValue(child));

   return SCIP_OKAY;
}

/** derivative evaluation callback:
 * computes <gradient, children.dot>
 * if expr is SIN(child), then computes
 * COS(child) dot(child)
 */
static
SCIP_DECL_CONSEXPR_EXPRFWDIFF(fwdiffSin)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child)), "val") != 0);
   assert(SCIPgetConsExprExprDot(child) != SCIP_INVALID); /*lint !e777*/

   *dot = cos(SCIPgetConsExprExprValue(child)) * SCIPgetConsExprExprDot(child);

   return SCIP_OKAY;
}

/** expression backward forward derivative evaluation callback
 * computes partial/partial child ( <gradient, children.dot> )
 * if expr is SIN(child), then computes
 * -SIN(child) dot(child)
 * */
static
SCIP_DECL_CONSEXPR_EXPRBWFWDIFF(bwfwdiffSin)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprValue(expr) != SCIP_INVALID); /*lint !e777*/
   assert(childidx == 0);

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child)), "val") != 0);
   assert(SCIPgetConsExprExprDot(child) != SCIP_INVALID); /*lint !e777*/

   *bardot = -sin(SCIPgetConsExprExprValue(child)) * SCIPgetConsExprExprDot(child);

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalSin)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   childinterval = SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(expr)[0]);
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      SCIPintervalSetEmpty(interval);
   else
      SCIPintervalSin(SCIP_INTERVAL_INFINITY, interval, childinterval);

   return SCIP_OKAY;
}

/** separation initialization callback */
static
SCIP_DECL_CONSEXPR_EXPRINITSEPA(initSepaSin)
{  /*lint --e{715}*/
   SCIP_VAR* childvar;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Bool success;

   SCIP_ROWPREP* cuts[5];   /* 0: secant, 1: left tangent, 2: right tangent, 3: left mid tangent, 4: right mid tangent */
   int i;

   *infeasible = FALSE;

   childvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[0]);
   assert(childvar != NULL);

   childlb = SCIPvarGetLbLocal(childvar);
   childub = SCIPvarGetUbLocal(childvar);

   /* no need for cut if child is fixed */
   if( SCIPisRelEQ(scip, childlb, childub) )
      return SCIP_OKAY;

   /* compute underestimating cuts */
   if( underestimate )
   {
      SCIP_CALL(SCIPcomputeInitialCutsTrig(scip, conshdlr, expr, &cuts[0], &cuts[1], &cuts[2], &cuts[3], &cuts[4],
            childlb, childub, TRUE) );

      for( i = 0; i < 5; ++i)
      {
         /* only the cuts which could be created are added */
         if( !*infeasible && cuts[i] != NULL )
         {
            SCIP_CALL( SCIPcleanupRowprep(scip, cuts[i], NULL, SCIP_CONSEXPR_CUTMAXRANGE, 0.0, NULL, &success) );

            if( success && cuts[i]->nvars == 2 )
            {
               /* make a SCIP_ROW and add to LP */
               SCIP_ROW* row;

               SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, cuts[i], cons) );
               SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }

            SCIPfreeRowprep(scip, &cuts[i]);
         }
      }
   }

   /* compute overestimating cuts */
   if( overestimate && !*infeasible )
   {
      SCIP_CALL(SCIPcomputeInitialCutsTrig(scip, conshdlr, expr, &cuts[0], &cuts[1], &cuts[2], &cuts[3], &cuts[4],
            childlb, childub, FALSE) );

      for( i = 0; i < 5; ++i )
      {
         /* only the cuts which could be created are added */
         if( !*infeasible && cuts[i] != NULL )
         {
            SCIP_CALL( SCIPcleanupRowprep(scip, cuts[i], NULL, SCIP_CONSEXPR_CUTMAXRANGE, 0.0, NULL, &success) );

            if( success && cuts[i]->nvars == 2 )
            {
               /* make a SCIP_ROW and add to LP */
               SCIP_ROW* row;

               SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, cuts[i], cons) );
               SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }

            SCIPfreeRowprep(scip, &cuts[i]);
         }
      }
   }

   return SCIP_OKAY;
}

/** expression estimator callback */
static
SCIP_DECL_CONSEXPR_EXPRESTIMATE(estimateSin)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* childvar;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(coefs != NULL);
   assert(constant != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);
   assert(success != NULL);

   /* get expression data */
   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   childvar = SCIPgetConsExprExprAuxVar(child);
   assert(childvar != NULL);

   *success = SCIPcomputeEstimatorsTrig(scip, conshdlr, expr, coefs, constant, SCIPgetSolVal(scip, sol, childvar),
      SCIPvarGetLbLocal(childvar), SCIPvarGetUbLocal(childvar), !overestimate);
   *islocal = TRUE;  /* TODO there are cases where cuts would be globally valid */

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_EXPRREVERSEPROP(reversepropSin)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_INTERVAL newbounds;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(nreductions != NULL);
   assert(SCIPintervalGetInf(bounds) >= -1.0);
   assert(SCIPintervalGetSup(bounds) <= 1.0);

   *nreductions = 0;

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);

   newbounds = SCIPgetConsExprExprBounds(scip, conshdlr, child);
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, newbounds) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* compute the new child interval */
   SCIP_CALL( SCIPcomputeRevPropIntervalSin(scip, bounds, newbounds, &newbounds) );

   /* try to tighten the bounds of the child node */
   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, conshdlr, child, newbounds, infeasible, nreductions) );

   return SCIP_OKAY;
}

/** sin hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashSin)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= childrenhashes[0];

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_CONSEXPR_EXPRCURVATURE(curvatureSin)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_INTERVAL childinterval;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(childcurv != NULL);
   assert(success != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, child, &childinterval, FALSE) );

   /* TODO rewrite SCIPcomputeCurvatureSin so it provides the reverse operation */
   *success = TRUE;
   if( SCIPcomputeCurvatureSin(SCIP_EXPRCURV_CONVEX, childinterval.inf, childinterval.sup) == exprcurvature )
      *childcurv = SCIP_EXPRCURV_CONVEX;
   else if( SCIPcomputeCurvatureSin(SCIP_EXPRCURV_CONCAVE, childinterval.inf, childinterval.sup) == exprcurvature )
      *childcurv = SCIP_EXPRCURV_CONCAVE;
   if( SCIPcomputeCurvatureSin(SCIP_EXPRCURV_LINEAR, childinterval.inf, childinterval.sup) == exprcurvature )
      *childcurv = SCIP_EXPRCURV_LINEAR;
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_CONSEXPR_EXPRMONOTONICITY(monotonicitySin)
{  /*lint --e{715}*/
   SCIP_INTERVAL interval;
   SCIP_Real inf;
   SCIP_Real sup;
   int k;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(childidx == 0);

   assert(SCIPgetConsExprExprChildren(expr)[0] != NULL);
   SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, SCIPgetConsExprExprChildren(expr)[0], &interval, FALSE) );

   *result = SCIP_MONOTONE_UNKNOWN;
   inf = SCIPintervalGetInf(interval);
   sup = SCIPintervalGetSup(interval);

   /* expression is not monotone because the interval is too large */
   if( sup - inf > M_PI )
      return SCIP_OKAY;

   /* compute k s.t. PI * (2k+1) / 2 <= interval.inf <= PI * (2k+3) / 2 */
   k = (int)floor(inf/M_PI - 0.5);
   assert(M_PI * (2.0*k + 1.0) / 2.0 <= inf);
   assert(M_PI * (2.0*k + 3.0) / 2.0 >= inf);

   /* check whether [inf,sup] are in containing in an interval for which the sine function is monotone */
   if( M_PI * (2.0*k + 3.0) / 2.0 <= sup )
      *result = ((k % 2 + 2) % 2) == 1 ? SCIP_MONOTONE_INC : SCIP_MONOTONE_DEC;

   return SCIP_OKAY;
}

/** creates the handler for sin expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   /* include expression handler */
   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC,
         EXPRHDLR_PRECEDENCE, evalSin, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrSin, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifySin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrParse(scip, consexprhdlr, exprhdlr, parseSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, initSepaSin, NULL, estimateSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrDiff(scip, consexprhdlr, exprhdlr, bwdiffSin, fwdiffSin, bwfwdiffSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCurvature(scip, consexprhdlr, exprhdlr, curvatureSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrMonotonicity(scip, consexprhdlr, exprhdlr, monotonicitySin) );

   return SCIP_OKAY;
}

/** creates a sin expression */
SCIP_RETCODE SCIPcreateConsExprExprSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPR*   child               /**< single child */
   )
{
   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindConsExprExprHdlr(consexprhdlr, EXPRHDLR_NAME) != NULL);

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPfindConsExprExprHdlr(consexprhdlr, EXPRHDLR_NAME), NULL, 1,
         &child) );

   return SCIP_OKAY;
}
