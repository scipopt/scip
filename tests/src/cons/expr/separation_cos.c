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

/**@file   separation_cos.c
 * @brief  tests separation of cos()
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
   SCIP_ROWPREP* soltangent;
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

   SCIP_CALL( SCIPcreateConsExprExprCos(scip, conshdlr, &expr, xexpr) );

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

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
      refpoint, childlb, childub, FALSE) );

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

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
         refpoint, childlb, childub, TRUE) );

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
    * test solution overestimation
    */

   /* create solution that induces overestimating cut */
   refpoint = -0.5;
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 2.0) );

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent */
   checkCut(soltangent, x, 0.5 * SIN(-0.5) - COS(-0.5), SCIPinfinity(scip), -SIN(-0.5));

   /* release cuts */
   SCIPfreeRowprep(scip, &soltangent);

   /*
    * test solution underestimation
    */

   /* create solution that induces underestimating cut */
   refpoint = 4.0;
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -2.0) );

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent */
   checkCut(soltangent, x, -SCIPinfinity(scip), -4.0 * SIN(4.0) - COS(4.0), -SIN(4.0));

   /* release cuts */
   SCIPfreeRowprep(scip, &soltangent);

   /*
    * test point where solution tangent is not feasible
    */

   /* create solution */
   refpoint = 1.7;
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -2.0) );

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(soltangent == NULL);

   /* check lmidtangent */
   newtonpoint = 2.740339007021;
   checkCut(lmidtangent, x, -SCIPinfinity(scip), SIN(newtonpoint) - COS(-1), -SIN(newtonpoint));

   /* check rmidtangent */
   newtonpoint = 4.568732341350;
   checkCut(rmidtangent, x, -SCIPinfinity(scip), -5 * SIN(newtonpoint) - COS(5), -SIN(newtonpoint));

   /* release cuts */
   SCIPfreeRowprep(scip, &lmidtangent);
   SCIPfreeRowprep(scip, &rmidtangent);

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
   SCIP_ROWPREP* soltangent;
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

   SCIP_CALL( SCIPcreateConsExprExprCos(scip, conshdlr, &expr, yexpr) );

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

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
      refpoint, childlb, childub, FALSE) );

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

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
      refpoint, childlb, childub, TRUE) );

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
    * test solution overestimation
    */

   /* create solution */
   refpoint = -5.7;
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 2.0) );

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent */
   checkCut(soltangent, y, 5.7 * SIN(-5.7) - COS(-5.7), SCIPinfinity(scip), -SIN(-5.7));

   /* release cuts */
   SCIPfreeRowprep(scip, &soltangent);

   /*
    * test solution underestimation
    */

   /* create solution */
   refpoint = -3.5;
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -2.0) );

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent */
   checkCut(soltangent, y, -SCIPinfinity(scip), 3.5 * SIN(-3.5) - COS(-3.5), -SIN(-3.5));

   /* release cuts */
   SCIPfreeRowprep(scip, &soltangent);

   /*
    * test point where solution tangent is infeasible
    */

   /* create solution */
   refpoint = -4.8;
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 2.0) );

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(soltangent == NULL);

   /* check rmidtangent */
   newtonpoint = -5.535897406992;
   checkCut(rmidtangent, y, 3 * SIN(newtonpoint) - COS(-3), SCIPinfinity(scip), -SIN(newtonpoint));

   /* release cuts */
   SCIPfreeRowprep(scip, &rmidtangent);

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
   SCIP_ROWPREP* soltangent;
   SCIP_Real refpoint;
   SCIP_Real childlb;
   SCIP_Real childub;

   secant = NULL;
   ltangent = NULL;
   rtangent = NULL;
   lmidtangent = NULL;
   rmidtangent = NULL;
   soltangent = NULL;

   SCIP_CALL( SCIPcreateConsExprExprCos(scip, conshdlr, &expr, wexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   refpoint = SCIP_INVALID;
   childlb = SCIPvarGetLbLocal(w);
   childub = SCIPvarGetUbLocal(w);

   /*
    * test initial overestimation
    */
   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
      refpoint, childlb, childub, FALSE) );

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

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
      refpoint, childlb, childub, TRUE) );

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
    * test solution underestimation
    */

   /* create solution */
   refpoint = 3.0;
   SCIP_CALL( SCIPsetSolVal(scip, sol, w, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 2.0) );

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, TRUE) );

   /* check cuts which could not be computed */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent */
   checkCut(soltangent, w, -SCIPinfinity(scip), -3 * SIN(3) - COS(3), -SIN(3));

   /* release cuts */
   SCIPfreeRowprep(scip, &soltangent);

   /*
    * test point where solution tangent is infeasible (overestimation)
    */

   /* create solution */
   refpoint = 3.0;
   SCIP_CALL( SCIPsetSolVal(scip, sol, w, refpoint) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -2.0) );

   SCIP_CALL( SCIPcomputeCutsTrig(scip, conshdlr, expr, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      refpoint, childlb, childub, FALSE) );

   /* check cuts which could not be computed */
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);
   cr_expect(soltangent == NULL);

   /* check secant */
   checkCut(secant, w, COS(4) - 2 * COS(2), SCIPinfinity(scip), 0.5 * (COS(4) - COS(2)));

   /* release cuts */
   SCIPfreeRowprep(scip, &secant);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
