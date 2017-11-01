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

/**@file   cons_expr_sin.c
 * @brief  handler for sin expressions
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define _USE_MATH_DEFINES   /* to get M_PI on Windows */  /*lint !750 */

#include <string.h>
#include <math.h>
#include "scip/cons_expr_sin.h"
#include "cons_expr_value.h"

/* fundamental expression handler properties */
#define EXPRHDLR_NAME         "sin"
#define EXPRHDLR_DESC         "sine expression"
#define EXPRHDLR_PRECEDENCE   91000
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(82457.0)

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
 */
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
SCIP_Bool SCIPcomputeSecantSin(
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
SCIP_Bool SCIPcomputeLeftTangentSin(
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
SCIP_Bool SCIPcomputeRightTangentSin(
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
SCIP_Bool SCIPcomputeSolTangentSin(
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
SCIP_Bool SCIPcomputeLeftMidTangentSin(
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
      if( SIN(lb) < 0.0 )
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
SCIP_Bool SCIPcomputeRightMidTangentSin(
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
      if( SIN(ub) < 0.0 )
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

/** helper function to create cuts for point- or initial separation
 *
 *  A total of 6 different cuts can be generated. All except soltangent are independent of a specific solution and
 *  use only the bounds of the child variable. If their pointers are passed with NULL, the respective coputation
 *  is not performed at all. If one of the computations failes or turns out to be irrelevant, the respective argument
 *  pointer is set to NULL.
 */
static
SCIP_RETCODE computeCutsSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_ROW**            secant,             /**< pointer to store the secant */
   SCIP_ROW**            ltangent,           /**< pointer to store the left tangent */
   SCIP_ROW**            rtangent,           /**< pointer to store the right tangent */
   SCIP_ROW**            lmidtangent,        /**< pointer to store the left middle tangent */
   SCIP_ROW**            rmidtangent,        /**< pointer to store the right middle tangent */
   SCIP_ROW**            soltangent,         /**< pointer to store the solution tangent */
   SCIP_Bool             underestimate       /**< whether the cuts should be underestimating */
   )
{
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* auxvar;
   SCIP_VAR* childvar;
   SCIP_Real lincoef;
   SCIP_Real linconst;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool success;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);

   /* get expression data */
   auxvar = SCIPgetConsExprExprLinearizationVar(expr);
   assert(auxvar != NULL);
   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   childvar = SCIPgetConsExprExprLinearizationVar(child);
   assert(childvar != NULL);
   childlb = SCIPvarGetLbLocal(childvar);
   childub = SCIPvarGetUbLocal(childvar);
   assert(SCIPisLE(scip, childlb, childub));

   /* if variable is fixed, it does not make sense to add cuts */
   if( SCIPisEQ(scip, childlb, childub) )
      return SCIP_OKAY;

   /*
    * Compute all cuts that where specified upon call.
    * For each linear equation z = a*x + b with bounds [lb,ub] the parameters can be computed by:
    *
    * a = cos(x^)    and     b = sin(x^) - a * x^        where x^ is any known point in [lb,ub]
    *
    * and the resulting cut is       a*x - z <=/>= -b           depending on over-/underestimation
    */

   /* compute secant between lower and upper bound */
   if( secant != NULL )
   {
      *secant = NULL;

      if( underestimate )
         success = SCIPcomputeSecantSin(scip, &lincoef, &linconst, childlb, childub);
      else
         success = SCIPcomputeSecantSin(scip, &lincoef, &linconst, -childub, -childlb);

      if( success )
      {
         lhs = underestimate ? -SCIPinfinity(scip) : linconst;
         rhs = underestimate ? -linconst : SCIPinfinity(scip);

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sin_secant_%s", SCIPvarGetName(childvar));

         SCIP_CALL( SCIPcreateEmptyRowCons(scip, secant, conshdlr, name, lhs, rhs,
            TRUE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddVarToRow(scip, *secant, auxvar, -1.0) );
         SCIP_CALL( SCIPaddVarToRow(scip, *secant, childvar, lincoef) );
      }
   }

   /* compute tangent at lower bound */
   if( ltangent != NULL )
   {
      *ltangent = NULL;

      if( underestimate )
         success = SCIPcomputeLeftTangentSin(scip, &lincoef, &linconst, childlb);
      else
         success = SCIPcomputeRightTangentSin(scip, &lincoef, &linconst, -childlb);

      if( success )
      {
         lhs = underestimate ? -SCIPinfinity(scip) : linconst;
         rhs = underestimate ? -linconst : SCIPinfinity(scip);

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sin_ltangent_%s", SCIPvarGetName(childvar));

         SCIP_CALL( SCIPcreateEmptyRowCons(scip, ltangent, conshdlr, name, lhs, rhs,
            TRUE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddVarToRow(scip, *ltangent, auxvar, -1.0) );
         SCIP_CALL( SCIPaddVarToRow(scip, *ltangent, childvar, lincoef) );
      }
   }

   /* compute tangent at upper bound */
   if( rtangent != NULL )
   {
      *rtangent = NULL;

      if( underestimate )
         success = SCIPcomputeRightTangentSin(scip, &lincoef, &linconst, childub);
      else
         success = SCIPcomputeLeftTangentSin(scip, &lincoef, &linconst, -childub);

      if( success )
      {
         lhs = underestimate ? -SCIPinfinity(scip) : linconst;
         rhs = underestimate ? -linconst : SCIPinfinity(scip);

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sin_rtangent_%s", SCIPvarGetName(childvar));

         SCIP_CALL( SCIPcreateEmptyRowCons(scip, rtangent, conshdlr, name, lhs, rhs,
            TRUE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddVarToRow(scip, *rtangent, auxvar, -1.0) );
         SCIP_CALL( SCIPaddVarToRow(scip, *rtangent, childvar, lincoef) );
      }
   }

   /* compute tangent at solution point */
   if( soltangent != NULL )
   {
      SCIP_Real refpoint;

      *soltangent = NULL;
      refpoint = SCIPgetSolVal(scip, sol, childvar);

      if( underestimate )
         success = SCIPcomputeSolTangentSin(scip, &lincoef, &linconst, childlb, childub, refpoint);
      else
         success = SCIPcomputeSolTangentSin(scip, &lincoef, &linconst, -childub, -childlb, -refpoint);

      if( success )
      {
         lhs = underestimate ? -SCIPinfinity(scip) : linconst;
         rhs = underestimate ? -linconst : SCIPinfinity(scip);

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sin_soltangent_%s", SCIPvarGetName(childvar));

         SCIP_CALL( SCIPcreateEmptyRowCons(scip, soltangent, conshdlr, name, lhs, rhs,
            TRUE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddVarToRow(scip, *soltangent, auxvar, -1.0) );
         SCIP_CALL( SCIPaddVarToRow(scip, *soltangent, childvar, lincoef) );
      }
   }


   /* compute left middle tangent, that is tangent at some other point which goes through (lb,sin(lb))
    * if secant or soltangent are feasible, this cut can never beat them
    */
   if( lmidtangent != NULL && (secant == NULL || *secant == NULL) && (soltangent == NULL || *soltangent == NULL) )
   {
      SCIP_ROW** cutbuffer;
      SCIP_Bool issecant;

      *lmidtangent = NULL;

      if( underestimate )
         success = SCIPcomputeLeftMidTangentSin(scip, &lincoef, &linconst, &issecant, childlb, childub);
      else
         success = SCIPcomputeRightMidTangentSin(scip, &lincoef, &linconst, &issecant, -childub, -childlb);

      /* if the cut connects bounds it is stored in secant */
      cutbuffer = (issecant && secant != NULL) ? secant : lmidtangent;

      if( success )
      {
         lhs = underestimate ? -SCIPinfinity(scip) : linconst;
         rhs = underestimate ? -linconst : SCIPinfinity(scip);

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sin_lmidtangent_%s", SCIPvarGetName(childvar));

         SCIP_CALL( SCIPcreateEmptyRowCons(scip, cutbuffer, conshdlr, name, lhs, rhs,
            TRUE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddVarToRow(scip, *cutbuffer, auxvar, -1.0) );
         SCIP_CALL( SCIPaddVarToRow(scip, *cutbuffer, childvar, lincoef) );
      }
   }
   else if( lmidtangent != NULL )
      *lmidtangent = NULL;

   /* compute right middle tangent, that is tangent at some other point which goes through (ub,sin(ub))
    * if secant or soltangent are feasible, this cut can never beat them
    */
   if( rmidtangent != NULL && (secant == NULL || *secant == NULL) && (soltangent == NULL || *soltangent == NULL) )
   {
      SCIP_ROW** cutbuffer;
      SCIP_Bool issecant;

      *rmidtangent = NULL;

      if( underestimate )
         success = SCIPcomputeRightMidTangentSin(scip, &lincoef, &linconst, &issecant, childlb, childub);
      else
         success = SCIPcomputeLeftMidTangentSin(scip, &lincoef, &linconst, &issecant, -childub, -childlb);

      /* if the cut connects bounds it is stored in secant */
      cutbuffer = (issecant && secant != NULL) ? secant : rmidtangent;

      if( success )
      {
         lhs = underestimate ? -SCIPinfinity(scip) : linconst;
         rhs = underestimate ? -linconst : SCIPinfinity(scip);

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sin_rmidtangent_%s", SCIPvarGetName(childvar));

         SCIP_CALL( SCIPcreateEmptyRowCons(scip, cutbuffer, conshdlr, name, lhs, rhs,
            TRUE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddVarToRow(scip, *cutbuffer, auxvar, -1.0) );
         SCIP_CALL( SCIPaddVarToRow(scip, *cutbuffer, childvar, lincoef) );
      }
   }
   else if( rmidtangent != NULL )
      *rmidtangent = NULL;

   return SCIP_OKAY;
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
   SCIP_CONSHDLR* conshdlr;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

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

/** expression data copy callback */
static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataSin)
{  /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetexprdata != NULL);
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);
   assert(SCIPgetConsExprExprData(sourceexpr) == NULL);

   *targetexprdata = NULL;

   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataSin)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

/** expression print callback */
static
SCIP_DECL_CONSEXPR_EXPRPRINT(printSin)
{  /*lint --e{715}*/
   assert(expr != NULL);

   switch( stage )
   {
   case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
   {
      /* print function with opening parenthesis */
      SCIPinfoMessage(scip, file, "%s(", EXPRHDLR_NAME);
      break;
   }

   case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
   {
      assert(SCIPgetConsExprExprWalkCurrentChild(expr) == 0);
      break;
   }

   case SCIP_CONSEXPREXPRWALK_LEAVEEXPR :
   {
      /* print closing parenthesis */
      SCIPinfoMessage(scip, file, ")");
      break;
   }

   case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
   default: ;
   }

   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_CONSEXPR_EXPRPARSE(parseSin)
{
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
   assert(idx >= 0 && idx < SCIPgetConsExprExprNChildren(expr));
   assert(SCIPgetConsExprExprValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPgetConsExprExprChildren(expr)[idx];
   assert(child != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child)), "val") != 0);

   *val = COS(SCIPgetConsExprExprValue(child));

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalSin)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   childinterval = SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[0]);
   assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), childinterval));

   SCIPintervalSin(SCIPinfinity(scip), interval, childinterval);

   return SCIP_OKAY;
}

/** separation initialization callback */
static
SCIP_DECL_CONSEXPR_EXPRINITSEPA(initSepaSin)
{
   SCIP_ROW* cuts[5];   /* 0: secant, 1: left tangent, 2: right tangent, 3: left mid tangent, 4: right mid tangent */
   int i;

   *infeasible = FALSE;

   /* compute underestimating cuts */
   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, NULL, &cuts[0], &cuts[1], &cuts[2], &cuts[3], &cuts[4], NULL, TRUE) );

   for( i = 0; i < 5; ++i)
   {
      /* only the cuts which could be created are added */
      if( !*infeasible && cuts[i] != NULL )
      {
         SCIP_CALL( SCIPmassageConsExprExprCut(scip, &cuts[i], NULL, -SCIPinfinity(scip)) );

         if( cuts[i] != NULL ) {
            SCIP_CALL( SCIPaddCut(scip, NULL, cuts[i], FALSE, infeasible) );
         }
      }

      /* release the row */
      if( cuts[i] != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &cuts[i]) );
      }
   }

   /* compute overestimating cuts */
   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, NULL, &cuts[0], &cuts[1], &cuts[2], &cuts[3], &cuts[4], NULL, FALSE) );

   for( i = 0; i < 5; ++i ) {
      /* only the cuts which could be created are added */
      if (!*infeasible && cuts[i] != NULL) {
         SCIP_CALL(SCIPmassageConsExprExprCut(scip, &cuts[i], NULL, -SCIPinfinity(scip)));

         if (cuts[i] != NULL) {
            SCIP_CALL(SCIPaddCut(scip, NULL, cuts[i], FALSE, infeasible));
         }
      }

      /* release the row */
      if (cuts[i] != NULL) {
         SCIP_CALL(SCIPreleaseRow(scip, &cuts[i]));
      }
   }

   return SCIP_OKAY;
}

/** expression separation callback */
static
SCIP_DECL_CONSEXPR_EXPRSEPA(sepaSin)
{
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* auxvar;
   SCIP_VAR* childvar;
   SCIP_ROW* cuts[4] = {NULL, NULL, NULL, NULL};
   SCIP_Real violation;
   SCIP_Bool underestimate;
   SCIP_Bool infeasible;
   int i;

   /* get expression data */
   auxvar = SCIPgetConsExprExprLinearizationVar(expr);
   assert(auxvar != NULL);
   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   childvar = SCIPgetConsExprExprLinearizationVar(child);
   assert(childvar != NULL);

   infeasible = FALSE;
   *ncuts = 0;
   *result = SCIP_DIDNOTFIND;

   /* compute the violation; this determines whether we need to over- or underestimate */
   violation = exp(SCIPgetSolVal(scip, sol, childvar)) - SCIPgetSolVal(scip, sol, auxvar);

   /* check if there is a violation */
   if( SCIPisEQ(scip, violation, 0.0) )
      return SCIP_OKAY;

   /* determine if we need to under- or overestimate */
   underestimate = SCIPisGT(scip, violation, 0.0);

   /* compute all possible inequalities; the resulting cuts are stored in the cuts array
    *
    *  - cuts[0] = secant
    *  - cuts[1] = secant connecting (lb,sin(lbx)) with left tangent point
    *  - cuts[1] = secant connecting (ub,sin(ubx)) with right tangent point
    *  - cuts[3] = solution tangent (for convex / concave segments that globally under- / overestimate)
    */
   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &cuts[0], NULL, NULL, &cuts[1], &cuts[2], &cuts[3],
      underestimate) );

   for( i = 0; i < 4; ++i )
   {
      if( cuts[i] == NULL )
         continue;

      SCIP_CALL( SCIPmassageConsExprExprCut(scip, &cuts[i], sol, minviolation) );

      if( cuts[i] == NULL )
         continue;

      violation = -SCIPgetRowSolFeasibility(scip, cuts[i], sol);
      if( SCIPisGE(scip, violation, minviolation) )
      {
         SCIP_CALL( SCIPaddCut(scip, sol, cuts[i], FALSE, &infeasible) );
         *ncuts += 1;

         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            break;
         }
         else
            *result = SCIP_SEPARATED;
      }

      /* release the secant */
      if( cuts[i] != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &cuts[i]) );
      }
   }

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_REVERSEPROP(reversepropSin)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_INTERVAL interval;
   SCIP_INTERVAL childbound;
   SCIP_Real newinf;
   SCIP_Real newsup;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(nreductions != NULL);
   assert(SCIPintervalGetInf(SCIPgetConsExprExprInterval(expr)) >= -1.0);
   assert(SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr)) <= 1.0);

   *nreductions = 0;

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);

   childbound = SCIPgetConsExprExprInterval(child);
   newinf = childbound.inf;
   newsup = childbound.sup;

   interval = SCIPgetConsExprExprInterval(expr);

   if( !SCIPisInfinity(scip, -newinf) )
   {
      /* l(x) and u(x) are lower/upper bound of child, l(s) and u(s) are lower/upper bound of sin expr
       *
       * if sin(l(x)) < l(s), we are looking for k minimal s.t. a + 2k*pi > l(x) where a = asin(l(s))
       * then the new lower bound is a + 2k*pi
       */
      if( SCIPisLT(scip, SIN(newinf), interval.inf) )
      {
         SCIP_Real a = ASIN(interval.inf);
         int k = (int) ceil((newinf - a) / (2.0*M_PI));
         newinf = a + 2.0*M_PI * k;
      }

      /* if sin(l(x)) > u(s), we are looking for k minimal s.t. pi - a + 2k*pi > l(x) where a = asin(u(s))
       * then the new lower bound is pi - a + 2k*pi
       */
      else if( SCIPisGT(scip, SIN(newinf), interval.sup) )
      {
         SCIP_Real a = ASIN(interval.sup);
         int k = (int) ceil((newinf + a) / (2.0*M_PI) - 0.5);
         newinf = M_PI * (2.0*k + 1.0) - a;
      }
   }

   if( !SCIPisInfinity(scip, newsup) )
   {
      /* if sin(u(x)) > u(s), we are looking for k minimal s.t. a + 2k*pi > u(x) - 2*pi where a = asin(u(s))
       * then the new upper bound is a + 2k*pi
       */
      if ( SCIPisGT(scip, SIN(newsup), interval.sup) )
      {
         SCIP_Real a = ASIN(interval.sup);
         int k = (int) ceil((newsup - a ) / (2.0*M_PI)) - 1;
         newsup = a + 2.0*M_PI * k;
      }

      /* if sin(u(x)) < l(s), we are looking for k minimal s.t. pi - a + 2k*pi > l(x) - 2*pi where a = asin(l(s))
       * then the new upper bound is pi - a + 2k*pi
       */
      if( SCIPisLT(scip, SIN(newsup), interval.inf) )
      {
         SCIP_Real a = ASIN(interval.inf);
         int k = (int) ceil((newsup + a) / (2.0*M_PI) - 0.5) - 1;
         newsup = M_PI * (2.0*k + 1.0) - a;
      }
   }

   assert(newinf >= childbound.inf);
   assert(newsup <= childbound.sup);
   assert(newinf <= newsup);
   assert(SCIPisGE(scip, SIN(newinf), interval.inf));
   assert(SCIPisLE(scip, SIN(newinf), interval.sup));
   assert(SCIPisGE(scip, SIN(newsup), interval.inf));
   assert(SCIPisLE(scip, SIN(newsup), interval.sup));

   SCIPintervalSetBounds(&childbound, newinf, newsup);

   /* try to tighten the bounds of the child node */
   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, child, childbound, force, infeasible, nreductions) );

   return SCIP_OKAY;
}

/** sin hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashSin)
{  /*lint --e{715}*/
   unsigned int childhash;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(expr2key != NULL);
   assert(hashkey != NULL);

   *hashkey = EXPRHDLR_HASHKEY;

   assert(SCIPhashmapExists(expr2key, (void*) SCIPgetConsExprExprChildren(expr)[0]));
   childhash = (unsigned int)(size_t) SCIPhashmapGetImage(expr2key, SCIPgetConsExprExprChildren(expr)[0]);

   *hashkey ^= childhash;

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
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataSin, freedataSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifySin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrParse(scip, consexprhdlr, exprhdlr, parseSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrInitSepa(scip, consexprhdlr, exprhdlr, initSepaSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, sepaSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashSin) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffSin) );

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
