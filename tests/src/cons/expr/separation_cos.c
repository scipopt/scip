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

/**@file   separation_cos.c
 * @brief  tests separation and estimation of cos()
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_cos.c"
#include "separation.h"

static
void checkCut(SCIP_ROWPREP* cut, SCIP_VAR* childvar, SCIP_Real lhs, SCIP_Real rhs, SCIP_Real varcoef)
{
   SCIP_VAR* var;
   SCIP_Real coef;
   SCIP_Real side;
   int i;

   cr_assert(cut != NULL);

   cr_assert(SCIPisInfinity(scip, -lhs) || SCIPisInfinity(scip, rhs));
   side = SCIPisInfinity(scip, -lhs) ? rhs : lhs;
   cr_expect(SCIPisEQ(scip, cut->side, side));

   cr_expect_eq(cut->nvars, 2);
   for( i = 0; i < 2; ++i )
   {
      var = cut->vars[i];
      coef = cut->coefs[i];

      if( var == childvar )
         cr_expect(SCIPisEQ(scip, coef, varcoef));
      else if( var == auxvar )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }
}

/* tests for interval [-1,5] */
Test(separation, cosinus_x, .init = setup, .fini = teardown,
   .description = "test separation for a cosinus expression in large range"
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROWPREP* secant;
   SCIP_ROWPREP* ltangent;
   SCIP_ROWPREP* rtangent;
   SCIP_ROWPREP* lmidtangent;
   SCIP_ROWPREP* rmidtangent;
   SCIP_Real newtonpoint;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Real linconst;
   SCIP_Real lincoef;
   SCIP_Bool success;

   secant = NULL;
   ltangent = NULL;
   rtangent = NULL;
   lmidtangent = NULL;
   rmidtangent = NULL;

   SCIP_CALL( SCIPcreateConsExprExprCos(scip, conshdlr, &expr, xexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   childlb = SCIPvarGetLbLocal(x);
   childub = SCIPvarGetUbLocal(x);

   /*
    * test initial overestimation
    */

   SCIP_CALL( SCIPcomputeInitialCutsTrig(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent,
         childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(rtangent == NULL);
   cr_expect(lmidtangent == NULL);

   /* check ltangent */
   checkCut(ltangent, x, SIN(-1) - COS(-1), SCIPinfinity(scip), -SIN(-1));

   /* check rmidtangent */
   newtonpoint = 0.145902085052;
   checkCut(rmidtangent, x, -5 * SIN(newtonpoint) - COS(5), SCIPinfinity(scip), -SIN(newtonpoint));


   /* release cuts */
   SCIPfreeRowprep(scip, &ltangent);
   SCIPfreeRowprep(scip, &rmidtangent);

   /*
    * test initial underestimation
    */

   SCIP_CALL( SCIPcomputeInitialCutsTrig(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent,
         childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(ltangent == NULL);
   cr_expect(rtangent == NULL);

   /* check lmidtangent */
   newtonpoint = 2.740339007021;
   checkCut(lmidtangent, x, -SCIPinfinity(scip), SIN(newtonpoint) - COS(-1), -SIN(newtonpoint));

   /* check rmidtangent */
   newtonpoint = 4.568732341350;
   checkCut(rmidtangent, x, -SCIPinfinity(scip), -5 * SIN(newtonpoint) - COS(5), -SIN(newtonpoint));

   /* release cuts */
   SCIPfreeRowprep(scip, &lmidtangent);
   SCIPfreeRowprep(scip, &rmidtangent);

   /*
    * test overestimation
    */

   success = SCIPcomputeEstimatorsTrig(scip, conshdlr, expr, &lincoef, &linconst,
      -0.5, childlb, childub, FALSE);
   cr_expect(success);
   cr_expect_float_eq(linconst, -0.5 * SIN(-0.5) + COS(-0.5), 1e-12);
   cr_expect_float_eq(lincoef, -SIN(-0.5), 1e-12);

   /*
    * test underestimation
    */

   success = SCIPcomputeEstimatorsTrig(scip, conshdlr, expr, &lincoef, &linconst,
      4.0, childlb, childub, TRUE);
   cr_expect(success);
   cr_expect_float_eq(linconst, 4.0 * SIN(4.0) + COS(4.0), 1e-12);
   cr_expect_float_eq(lincoef, -SIN(4.0), 1e-12);

   /*
    * test point where solution tangent is not overestimating
    */
   success = SCIPcomputeEstimatorsTrig(scip, conshdlr, expr, &lincoef, &linconst,
      1.7, childlb, childub, TRUE);
   cr_expect(success);

   /* check lmidtangent */
   newtonpoint = 2.740339007021;
   cr_expect_float_eq(linconst, -SIN(newtonpoint) + COS(-1), 1e-12);
   cr_expect_float_eq(lincoef, -SIN(newtonpoint), 1e-12);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* tests for interval [-6,-3] */
Test(separation, cosinus_y, .init = setup, .fini = teardown,
   .description = "test separation for a cosine expression in mid size range"
)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROWPREP* secant;
   SCIP_ROWPREP* ltangent;
   SCIP_ROWPREP* rtangent;
   SCIP_ROWPREP* lmidtangent;
   SCIP_ROWPREP* rmidtangent;
   SCIP_Real newtonpoint;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Real linconst;
   SCIP_Real lincoef;
   SCIP_Bool success;

   secant = NULL;
   ltangent = NULL;
   rtangent = NULL;
   lmidtangent = NULL;
   rmidtangent = NULL;

   SCIP_CALL( SCIPcreateConsExprExprCos(scip, conshdlr, &expr, yexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   childlb = SCIPvarGetLbLocal(y);
   childub = SCIPvarGetUbLocal(y);

   /*
    * test initial overestimation
    */

   SCIP_CALL( SCIPcomputeInitialCutsTrig(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent,
      childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(ltangent == NULL);
   cr_expect(rtangent == NULL);
   cr_expect(lmidtangent == NULL);

   /* check rmidtangent */
   newtonpoint = -5.535897406992;
   checkCut(rmidtangent, y, 3 * SIN(newtonpoint) - COS(-3), SCIPinfinity(scip), -SIN(newtonpoint));

   /* release cuts */
   SCIPfreeRowprep(scip, &rmidtangent);

   /*
    * test initial underestimation
    */

   SCIP_CALL( SCIPcomputeInitialCutsTrig(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent,
      childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(ltangent == NULL);
   cr_expect(rmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check rtangent */
   checkCut(rtangent, y, -SCIPinfinity(scip), 3 * SIN(-3) - COS(-3), -SIN(-3));

   /* check lmidtangent */
   newtonpoint = -4.082240818442;
   checkCut(lmidtangent, y, -SCIPinfinity(scip), 6 * SIN(newtonpoint) - COS(-6), -SIN(newtonpoint));

   /* release cuts */
   SCIPfreeRowprep(scip, &rtangent);
   SCIPfreeRowprep(scip, &lmidtangent);

   /*
    * test overestimation
    */

   success = SCIPcomputeEstimatorsTrig(scip, conshdlr, expr, &lincoef, &linconst,
      -5.7, childlb, childub, FALSE);
   cr_expect(success);
   cr_expect_float_eq(linconst, -5.7 * SIN(-5.7) + COS(-5.7), 1e-12);
   cr_expect_float_eq(lincoef, -SIN(-5.7), 1e-12);

   /*
    * test underestimation
    */

   success = SCIPcomputeEstimatorsTrig(scip, conshdlr, expr, &lincoef, &linconst,
      -3.5, childlb, childub, TRUE);
   cr_expect(success);
   cr_expect_float_eq(linconst, -3.5 * SIN(-3.5) + COS(-3.5), 1e-12);
   cr_expect_float_eq(lincoef, -SIN(-3.5), 1e-12);

   /*
    * test point where solution tangent in not overestimating
    */

   success = SCIPcomputeEstimatorsTrig(scip, conshdlr, expr, &lincoef, &linconst,
      -4.8, childlb, childub, FALSE);
   cr_expect(success);
   newtonpoint = -5.535897406992;
   cr_expect_float_eq(linconst, -3 * SIN(newtonpoint) + COS(-3), 1e-11);
   cr_expect_float_eq(lincoef, -SIN(newtonpoint), 1e-12);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* tests for interval [2,4] */
Test(separation, cosinus_w, .init = setup, .fini = teardown,
   .description = "test separation for a cosine expression in short range"
)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROWPREP* secant;
   SCIP_ROWPREP* ltangent;
   SCIP_ROWPREP* rtangent;
   SCIP_ROWPREP* lmidtangent;
   SCIP_ROWPREP* rmidtangent;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Real linconst;
   SCIP_Real lincoef;
   SCIP_Bool success;

   secant = NULL;
   ltangent = NULL;
   rtangent = NULL;
   lmidtangent = NULL;
   rmidtangent = NULL;

   SCIP_CALL( SCIPcreateConsExprExprCos(scip, conshdlr, &expr, wexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   childlb = SCIPvarGetLbLocal(w);
   childub = SCIPvarGetUbLocal(w);

   /*
    * test initial overestimation
    */
   SCIP_CALL( SCIPcomputeInitialCutsTrig(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent,
      childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(ltangent == NULL);
   cr_expect(rtangent == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check secant */
   checkCut(secant, w, COS(4) - 2 * COS(2), SCIPinfinity(scip), 0.5 * (COS(4) - COS(2)));

   /* release cuts */
   SCIPfreeRowprep(scip, &secant);

   /*
    * test initial underestimation
    */

   SCIP_CALL( SCIPcomputeInitialCutsTrig(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent,
      childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check ltangent */
   checkCut(ltangent, w, -SCIPinfinity(scip), -2 * SIN(2) - COS(2), -SIN(2));

   /* check rtangent */
   checkCut(rtangent, w, -SCIPinfinity(scip), -4 * SIN(4) - COS(4), -SIN(4));

   /* release cuts */
   SCIPfreeRowprep(scip, &ltangent);
   SCIPfreeRowprep(scip, &rtangent);

   /* note: solution overestimation doesn't make sense in [1,3] */

   /*
    * test underestimation
    */

   success = SCIPcomputeEstimatorsTrig(scip, conshdlr, expr, &lincoef, &linconst,
      3.0, childlb, childub, TRUE);
   cr_expect(success);
   cr_expect_float_eq(linconst, 3 * SIN(3) + COS(3), 1e-12);
   cr_expect_float_eq(lincoef, -SIN(3), 1e-12);

   /*
    * test point where solution tangent is not overestimating
    */

   success = SCIPcomputeEstimatorsTrig(scip, conshdlr, expr, &lincoef, &linconst,
      3.0, childlb, childub, FALSE);
   cr_expect(success);
   cr_expect_float_eq(linconst, -COS(4) + 2 * COS(2), 1e-12);
   cr_expect_float_eq(lincoef, 0.5 * (COS(4) - COS(2)), 1e-12);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
