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
#define NEWTON_NITERATIONS    100
#define NEWTON_PRECISION      1e-10

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
   SCIP_Bool             overestimate,       /**< whether overestimating cuts are wanted */
   SCIP_Bool             underestimate       /**< whether underestimating cuts are wanted */
   )
{
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* auxvar;                     /* auxiliary variable for sine expression */
   SCIP_VAR* childvar;                   /* auxiliary variable for child expression */
   SCIP_Real params[2];                  /* array to store linear coefficient and constant for cuts */
   SCIP_Real childlb;                    /* lower bound of x */
   SCIP_Real childub;                    /* upper bound of x */
   SCIP_Real shiftedlbhalf;              /* lb mod pi */
   SCIP_Real shiftedubhalf;              /* ub mod pi */
   SCIP_Bool isoverestimating;           /* whether a computed cut is overestimating */

   SCIP_Real refpoint;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(overestimate || underestimate);

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
   {
      return SCIP_OKAY;
   }

   /* compute shifted bounds for case evaluation */
   shiftedlbhalf = fmod(childlb, M_PI);
   shiftedubhalf = fmod(childub, M_PI);
   if( childlb < 0.0 )
   {
      shiftedlbhalf += M_PI;
   }
   if( childub < 0.0 )
   {
      shiftedubhalf += M_PI;
   }

   /*
    * Compute all cuts that where specified upon call.
    * For each linear equation z = a*x + b with bounds [lb,ub] the parameters can be computed by:
    *
    * a = cos(x^)    and     b = sin(x^) - a * x^        where x^ is any known point in [lb,ub]
    *
    * and the resulting cut is       a*x - z <=/>= -b           depending on bay
    */

   /* compute secant between lower and upper bound */
   if( secant != NULL )
   {
      *secant = NULL;

      /* if bounds are not within one bay, secant does not make sense */
      if( childub - childlb <= M_PI - shiftedlbhalf )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sin_secant_%s", SCIPvarGetName(childvar));

         params[0] = (SIN(childub) - SIN(childlb)) / (childub - childlb);
         params[1] = SIN(childub) - params[0] * childub;

         /* inside the negative bay the cut is overestimating */
         isoverestimating = (SIN(childlb) < 0.0 || (SCIPisZero(scip, SIN(childlb)) && COS(childlb) < 0.0));

         /* inside the negative bay the cut is overestimating */
         if( isoverestimating && overestimate )
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, secant, conshdlr, name, -params[1], SCIPinfinity(scip),
                  TRUE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddVarToRow(scip, *secant, auxvar, -1.0) );
            SCIP_CALL( SCIPaddVarToRow(scip, *secant, childvar, params[0]) );
         }
         else if( !isoverestimating && underestimate )
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, secant, conshdlr, name, -SCIPinfinity(scip), -params[1],
                  TRUE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddVarToRow(scip, *secant, auxvar, -1.0) );
            SCIP_CALL( SCIPaddVarToRow(scip, *secant, childvar, params[0]) );
         }

      }
   }

   /* compute tangent at lower bound */
   if( ltangent != NULL )
   {
      *ltangent = NULL;

      if( shiftedlbhalf < 0.5*M_PI )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sin_ltangent_%s", SCIPvarGetName(childvar));

         params[0] = COS(childlb);
         params[1] = SIN(childlb) - params[0] * childlb;

         /* inside the positive bay the cut ist overestimating */
         isoverestimating = (SIN(childlb) > 0.0 || (SCIPisZero(scip, SIN(childlb)) && COS(childlb) > 0.0));

         if( isoverestimating && overestimate )
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, ltangent, conshdlr, name, -params[1], SCIPinfinity(scip),
                  TRUE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddVarToRow(scip, *ltangent, auxvar, -1.0) );
            SCIP_CALL( SCIPaddVarToRow(scip, *ltangent, childvar, params[0]) );
         }
         else if( !isoverestimating && underestimate )
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, ltangent, conshdlr, name, -SCIPinfinity(scip), -params[1],
                  TRUE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddVarToRow(scip, *ltangent, auxvar, -1.0) );
            SCIP_CALL( SCIPaddVarToRow(scip, *ltangent, childvar, params[0]) );
         }
      }
   }

   /* compute tangent at upper bound */
   if( rtangent != NULL )
   {
      *rtangent = NULL;

      if( shiftedubhalf > 0.5*M_PI || shiftedubhalf == 0.0 )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sin_rtangent_%s", SCIPvarGetName(childvar));

         params[0] = COS(childub);
         params[1] = SIN(childub) - params[0] * childub;

         /* inside the positive bay the cut ist overestimating */
         isoverestimating = (SIN(childub) > 0.0 || (SCIPisZero(scip, SIN(childlb)) && COS(childlb) < 0.0));

         if( isoverestimating && overestimate )
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, rtangent, conshdlr, name, -params[1], SCIPinfinity(scip),
                  TRUE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddVarToRow(scip, *rtangent, auxvar, -1.0) );
            SCIP_CALL( SCIPaddVarToRow(scip, *rtangent, childvar, params[0]) );

         }
         else if( !isoverestimating && underestimate )
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, rtangent, conshdlr, name, -SCIPinfinity(scip), -params[1],
                  TRUE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddVarToRow(scip, *rtangent, auxvar, -1.0) );
            SCIP_CALL( SCIPaddVarToRow(scip, *rtangent, childvar, params[0]) );

         }
      }
   }

   /* compute tangent at solution point */
   if( soltangent != NULL )
   {
      SCIP_Real startingpoints[3];
      SCIP_Real shiftedpointhalf;
      SCIP_Real intersection;
      int i;

      *soltangent = NULL;

      refpoint = SCIPgetSolVal(scip, sol, childvar);

      /* compute refpoint mod pi */
      shiftedpointhalf = fmod(refpoint, M_PI);
      if( refpoint < 0.0 )
      {
         shiftedpointhalf += M_PI;
      }

      /* if the point is too far away from the bounds or is at a multiple of pi, the cut will not help */
      if( SCIPisLE(scip, refpoint - childlb, 2*M_PI) && SCIPisLE(scip, childub - refpoint, 2*M_PI)
         && !SCIPisZero(scip, shiftedpointhalf) )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sin_soltangent_%s", SCIPvarGetName(childvar));

         params[0] = COS(refpoint);
         params[1] = SIN(refpoint) - params[0] * refpoint;

         /* choose starting points for newton procedure */
         if( SCIPisGT(scip, shiftedpointhalf, 0.5*M_PI) )
         {
            startingpoints[0] = refpoint + (M_PI - shiftedpointhalf) + 0.5*M_PI;
            startingpoints[1] = startingpoints[0] + 0.5*M_PI;
            startingpoints[2] = startingpoints[1] + 0.5*M_PI;
         }
         else
         {
            startingpoints[0] = refpoint - shiftedpointhalf - 0.5*M_PI;
            startingpoints[1] = startingpoints[0] - 0.5*M_PI;
            startingpoints[2] = startingpoints[1] - 0.5*M_PI;
         }

         /* use newton procedure to test if cut is valid */
         for( i = 0; i < 3; ++i)
         {
            intersection = SCIPcomputeRootNewton(function1, derivative1, params, 2, startingpoints[i], NEWTON_PRECISION,
               NEWTON_NITERATIONS);

            if( intersection != SCIP_INVALID && intersection != refpoint)
               break;
         }

         if( intersection != SCIP_INVALID && (intersection <= childlb || intersection >= childub) ) /*lint !e777*/
         {
            isoverestimating = (SIN(refpoint) > 0.0);

            if( isoverestimating && overestimate )
            {
               SCIP_CALL( SCIPcreateEmptyRowCons(scip, soltangent, conshdlr, name, -params[1], SCIPinfinity(scip),
                  TRUE, FALSE, FALSE) );

               SCIP_CALL( SCIPaddVarToRow(scip, *soltangent, auxvar, -1.0) );
               SCIP_CALL( SCIPaddVarToRow(scip, *soltangent, childvar, params[0]) );
            }
            else if( !isoverestimating && underestimate )
            {
               SCIP_CALL( SCIPcreateEmptyRowCons(scip, soltangent, conshdlr, name, -SCIPinfinity(scip), -params[1],
                  TRUE, FALSE, FALSE) );

               SCIP_CALL( SCIPaddVarToRow(scip, *soltangent, auxvar, -1.0) );
               SCIP_CALL( SCIPaddVarToRow(scip, *soltangent, childvar, params[0]) );
            }
         }
      }
   }

   /* compute left middle tangent, that is tangent at some other point which goes through (lb,sin(lb))
    * if secant or soltangent are feasible, this cut can never beat them
    */
   if( lmidtangent != NULL && (secant == NULL || *secant == NULL) && (soltangent == NULL || *soltangent == NULL) )
   {
      SCIP_Real tangentpoint;
      SCIP_Real connectionpoint;
      SCIP_Real testpoint;

      *lmidtangent = NULL;

      /* use newton procedure to find the point where the tangent intersects sine at lower bound */
      tangentpoint = SCIPcomputeRootNewton(function2, derivative2, &childlb, 1, childlb + 2*M_PI, NEWTON_PRECISION,
         NEWTON_NITERATIONS);

      /* if newton procedure failed, no cut is added */
      if( tangentpoint != SCIP_INVALID ) /*lint !e777*/
      {
         /* if the computed point lies outside the bounds, it is shifted to upper bound */
         connectionpoint = SCIPisGE(scip, tangentpoint, childub) ? childub : tangentpoint;

         /* compute secant between lower bound and connection point */
         params[0] = (SIN(connectionpoint) - SIN(childlb)) / (connectionpoint - childlb );
         params[1] = SIN(childlb) - params[0] * childlb;

         /* determine whether the cut is under- or overestimating */
         testpoint = 0.5 * (connectionpoint - childlb);
         isoverestimating = SCIPisGT(scip, params[0] * testpoint + params[1] - SIN(testpoint),  0.0);

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sin_lmidtangent_%s", SCIPvarGetName(childvar));

         if( isoverestimating && overestimate )
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, lmidtangent, conshdlr, name, -params[1], SCIPinfinity(scip),
                  TRUE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddVarToRow(scip, *lmidtangent, auxvar, -1.0) );
            SCIP_CALL( SCIPaddVarToRow(scip, *lmidtangent, childvar, params[0]) );
         }
         else if( !isoverestimating && underestimate )
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, lmidtangent, conshdlr, name, -SCIPinfinity(scip), -params[1],
                  TRUE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddVarToRow(scip, *lmidtangent, auxvar, -1.0) );
            SCIP_CALL( SCIPaddVarToRow(scip, *lmidtangent, childvar, params[0]) );
         }
      }
   }
   else if( lmidtangent != NULL )
      *lmidtangent = NULL;

   /* compute right middle tangent, that is tangent at some other point which goes through (ub,sin(ub))
    * if secant or soltangent are feasible, this cut can never beat them
    */
   if( rmidtangent != NULL && (secant == NULL || *secant != NULL) && (soltangent == NULL || *soltangent == NULL) )
   {
      SCIP_Real tangentpoint;
      SCIP_Real connectionpoint;
      SCIP_Real testpoint;

      *rmidtangent = NULL;

      /* use newton procedure to find the point where the tangent intersects sine at upper bound */
      tangentpoint = SCIPcomputeRootNewton(function2, derivative2, &childub, 1, childub - 2*M_PI, NEWTON_PRECISION,
         NEWTON_NITERATIONS);

      /* if newton procedure failed, no cut is added */
      if( tangentpoint != SCIP_INVALID ) /*lint !e777*/
      {
         /* if the computed point lies outside the bounds, it is shifted to lower bound */
         connectionpoint = SCIPisLE(scip, tangentpoint, childlb) ? childlb : tangentpoint;

         /* compute secant between upper bound and connection point */
         params[0] = (SIN(connectionpoint) - SIN(childub)) / (connectionpoint - childub );
         params[1] = SIN(childub) - params[0] * childub;

         /* determine whether the cut is under- or overestimating */
         testpoint = 0.5 * (childub - connectionpoint);
         isoverestimating = SCIPisGT(scip, params[0] * testpoint + params[1] - SIN(testpoint),  0.0);

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sin_rmidtangent_%s", SCIPvarGetName(childvar));

         if( isoverestimating && overestimate )
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, rmidtangent, conshdlr, name, -params[1], SCIPinfinity(scip),
               TRUE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddVarToRow(scip, *rmidtangent, auxvar, -1.0) );
            SCIP_CALL( SCIPaddVarToRow(scip, *rmidtangent, childvar, params[0]) );
         }
         else if( !isoverestimating && underestimate )
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, rmidtangent, conshdlr, name, -SCIPinfinity(scip), -params[1],
               TRUE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddVarToRow(scip, *rmidtangent, auxvar, -1.0) );
            SCIP_CALL( SCIPaddVarToRow(scip, *rmidtangent, childvar, params[0]) );
         }
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

   /* release child expression since it has been captured by the absolute expression */
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

   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, NULL, &cuts[0], &cuts[1], &cuts[2], &cuts[3], &cuts[4], NULL,
         TRUE, TRUE) );

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
   SCIP_Bool overestimate;
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
   overestimate = SCIPisLT(scip, violation, 0.0);

   /* compute all possible inequalities; the resulting cuts are stored in the cuts array
    *
    *  - cuts[0] = secant
    *  - cuts[1] = secant connecting (lb,sin(lbx)) with left tangent point
    *  - cuts[1] = secant connecting (ub,sin(ubx)) with right tangent point
    *  - cuts[3] = solution tangent (for convex / concave segments that globally under- / overestimate)
    */
   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &cuts[0], NULL, NULL, &cuts[1], &cuts[2], &cuts[3],
         overestimate, !overestimate) );

   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL( SCIPmassageConsExprExprCut(scip, &cuts[i], sol, minviolation) );

      if( cuts[i] == NULL || SCIProwIsInLP(cuts[i]) )
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
         assert(newinf >= childbound.inf);
      }

      /* if sin(l(x)) > u(s), we are looking for k minimal s.t. pi - a + 2k*pi > l(x) where a = asin(u(s))
       * then the new lower bound is pi - a + 2k*pi
       */
      else if( SCIPisGT(scip, SIN(newinf), interval.sup) )
      {
         SCIP_Real a = ASIN(interval.sup);
         int k = (int) ceil((newinf + a) / (2.0*M_PI) - 0.5);
         newinf = M_PI * (2.0*k + 1.0) - a;
         assert(newinf >= childbound.inf);
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
         assert(newsup <= childbound.sup);
      }

      /* if sin(u(x)) < l(s), we are looking for k minimal s.t. pi - a + 2k*pi > l(x) - 2*pi where a = asin(l(s))
       * then the new upper bound is pi - a + 2k*pi
       */
      if( SCIPisLT(scip, SIN(newsup), interval.inf) )
      {
         SCIP_Real a = ASIN(interval.inf);
         int k = (int) ceil((newsup + a) / (2*M_PI) - 0.5) - 1;
         newsup = M_PI * (2.0*k + 1.0) - a;
         assert(newsup <= childbound.sup);
      }
   }

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
