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

/**@file   separation_sin.c
 * @brief  tests separation of sin()
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_sin.c"
#include "separation.h"

static
void checkCut(SCIP_ROW* cut, SCIP_VAR* childvar, SCIP_Real lhs, SCIP_Real rhs, SCIP_Real varcoef)
{
   SCIP_VAR* var;
   SCIP_Real coef;
   int i;

   cr_assert(cut != NULL);
   cr_expect_eq(SCIProwGetNNonz(cut), 2);
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(cut), lhs));
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(cut), rhs));

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];

      if( var == SCIPvarGetTransVar(childvar) )
         cr_expect(SCIPisEQ(scip, coef, varcoef));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }
}

/* tests for interval [-1,5] */
Test(separation, sinus_x, .init = setup, .fini = teardown,
   .description = "test separation for a sinus expression in large range"
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROW* secant;
   SCIP_ROW* ltangent;
   SCIP_ROW* rtangent;
   SCIP_ROW* lmidtangent;
   SCIP_ROW* rmidtangent;
   SCIP_ROW* soltangent;
   SCIP_Real newtonpoint;
   SCIP_Real refpoint;
   SCIP_Real childlb;
   SCIP_Real childub;

   secant = NULL;
   ltangent = NULL;
   rtangent = NULL;
   lmidtangent = NULL;
   rmidtangent = NULL;
   soltangent = NULL;

   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr, xexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   refpoint = SCIP_INVALID;
   childlb = SCIPvarGetLbLocal(x);
   childub = SCIPvarGetUbLocal(x);

   /*
    * test initial overestimation
    */

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
      refpoint, childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(ltangent == NULL);
   cr_expect(rtangent == NULL);

   /* check lmidtangent */
   newtonpoint = 0.4936608602;
   checkCut(lmidtangent, x, -COS(newtonpoint) - SIN(-1), SCIPinfinity(scip), COS(newtonpoint));

   /* check rmidtangent */
   newtonpoint = 2.2544608804;
   checkCut(rmidtangent, x, 5 * COS(newtonpoint) - SIN(5), SCIPinfinity(scip), COS(newtonpoint));

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &lmidtangent) );
   SCIP_CALL( SCIPreleaseRow(scip, &rmidtangent) );

   /*
    * test initial underestimation
    */

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
         refpoint, childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(ltangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check rtangent */
   checkCut(rtangent, x, -SCIPinfinity(scip), 5 * COS(5) - SIN(5), COS(5));

   /* check lmidtangent */
   newtonpoint = 4.6845658560;
   checkCut(lmidtangent, x, -SCIPinfinity(scip), -COS(newtonpoint) - SIN(-1), COS(newtonpoint));

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &rtangent) );
   SCIP_CALL( SCIPreleaseRow(scip, &lmidtangent) );

   /*
    * test solution overestimation
    */

   /* create solution that induces overestimating cut */
   refpoint = 1.5;
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 2.0) );

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent */
   checkCut(soltangent, x, 1.5 * COS(1.5) - SIN(1.5), SCIPinfinity(scip), COS(1.5));

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &soltangent) );

   /*
    * test solution underestimation
    */

   /* create solution that induces underestimating cut */
   refpoint = 4.8;
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -2.0) );

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent */
   checkCut(soltangent, x, -SCIPinfinity(scip), 4.8 * COS(4.8) - SIN(4.8), COS(4.8));

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &soltangent) );

   /*
    * test point where solution tangent is not feasible
    */

   /* create solution */
   refpoint = 4.0;
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -2.0) );

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(rmidtangent == NULL);
   cr_expect(soltangent == NULL);

   /* check lmidtangent */
   newtonpoint = 4.6845658560;
   checkCut(lmidtangent, x, -SCIPinfinity(scip), -COS(newtonpoint) - SIN(-1), COS(newtonpoint));

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &lmidtangent) );

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* tests for interval [-6,-3] */
Test(separation, sinus_y, .init = setup, .fini = teardown,
   .description = "test separation for a sinus expression in mid size range"
)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROW* secant;
   SCIP_ROW* ltangent;
   SCIP_ROW* rtangent;
   SCIP_ROW* lmidtangent;
   SCIP_ROW* rmidtangent;
   SCIP_ROW* soltangent;
   SCIP_Real newtonpoint;
   SCIP_Real refpoint;
   SCIP_Real childlb;
   SCIP_Real childub;

   secant = NULL;
   ltangent = NULL;
   rtangent = NULL;
   lmidtangent = NULL;
   rmidtangent = NULL;
   soltangent = NULL;

   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr, yexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   refpoint = SCIP_INVALID;
   childlb = SCIPvarGetLbLocal(y);
   childub = SCIPvarGetUbLocal(y);

   /*
    * test initial overestimation
    */

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
      refpoint, childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(rtangent == NULL);
   cr_expect(lmidtangent == NULL);

   /* check ltangent */
   checkCut(ltangent, y, (-6) * COS(-6) - SIN(-6), SCIPinfinity(scip), COS(-6));

   /* check rmidtangent */
   newtonpoint = -3.2123712333;
   checkCut(rmidtangent, y, (-3) * COS(newtonpoint) - SIN(-3), SCIPinfinity(scip), COS(newtonpoint));

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &ltangent) );
   SCIP_CALL( SCIPreleaseRow(scip, &rmidtangent) );

   /*
    * test initial underestimation
    */

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
      refpoint, childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(ltangent == NULL);
   cr_expect(rtangent == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check secant */
   checkCut(secant, y, -SCIPinfinity(scip), SIN(-6) - 2.0 * SIN(-3), (SIN(-3) - SIN(-6)) / 3.0);

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &secant) );

   /*
    * test solution overestimation
    */

   /* create solution */
   refpoint = -4.0;
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 2.0) );

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent */
   checkCut(soltangent, y, (-4) * COS(-4) - SIN(-4), SCIPinfinity(scip), COS(-4));

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &soltangent) );

   /*
    * test solution underestimation (not possible in this case)
    */

   /* create solution */
   refpoint = -3.1;
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -2.0) );

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);
   cr_expect(soltangent == NULL);

   /* check secant */
   checkCut(secant, y, -SCIPinfinity(scip), SIN(-6) - 2.0 * SIN(-3), (SIN(-3) - SIN(-6)) / 3.0);

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &secant) );

   /*
    * test point where solution tangent is infeasible
    */

   /* create solution */
   refpoint = -3.2;
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 2.0) );

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(soltangent == NULL);

   /* check rmidtangent */
   newtonpoint = -3.2123712333;
   checkCut(rmidtangent, y, (-3) * COS(newtonpoint) - SIN(-3), SCIPinfinity(scip), COS(newtonpoint));

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &rmidtangent) );

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

}

/* tests for interval [1,3] */
Test(separation, sinus_z, .init = setup, .fini = teardown,
   .description = "test separation for a sinus expression in short range"
)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROW* secant;
   SCIP_ROW* ltangent;
   SCIP_ROW* rtangent;
   SCIP_ROW* lmidtangent;
   SCIP_ROW* rmidtangent;
   SCIP_ROW* soltangent;
   SCIP_Real refpoint;
   SCIP_Real childlb;
   SCIP_Real childub;

   secant = NULL;
   ltangent = NULL;
   rtangent = NULL;
   lmidtangent = NULL;
   rmidtangent = NULL;
   soltangent = NULL;

   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr, zexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   refpoint = SCIP_INVALID;
   childlb = SCIPvarGetLbLocal(z);
   childub = SCIPvarGetUbLocal(z);

   /*
    * test initial overestimation
    */
   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
      refpoint, childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check ltangent */
   checkCut(ltangent, z, COS(1) - SIN(1), SCIPinfinity(scip), COS(1));

   /* check rtangent */
   checkCut(rtangent, z, 3 * COS(3) - SIN(3), SCIPinfinity(scip), COS(3));

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &ltangent) );
   SCIP_CALL( SCIPreleaseRow(scip, &rtangent) );

   /*
    * test initial underestimation
    */

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
      refpoint, childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(ltangent == NULL);
   cr_expect(rtangent == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check secant */
   checkCut(secant, z, -SCIPinfinity(scip), 0.5 * SIN(3) - 1.5 * SIN(1), 0.5 * (SIN(3) - SIN(1)));

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &secant) );

   /*
    * test solution overestimation
    */

   /* create solution */
   refpoint = 2.0;
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 2.0) );

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent */
   checkCut(soltangent, z, 2 * COS(2) - SIN(2), SCIPinfinity(scip), COS(2));

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &soltangent) );

   /* note: solution underestimation doesn't make sense in [1,3] */

   /*
    * test point where solution tangent is infeasible (underestimation)
    */

   /* create solution */
   refpoint = 2.0;
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -2.0) );

   SCIP_CALL( SCIPcomputeCutsSin(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);
   cr_expect(soltangent == NULL);

   /* check secant */
   checkCut(secant, z, -SCIPinfinity(scip), 0.5 * SIN(3) - 1.5 * SIN(1), 0.5 * (SIN(3) - SIN(1)));

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &secant) );

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}