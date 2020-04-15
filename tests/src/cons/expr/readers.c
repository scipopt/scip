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

/**@file   readers.c
 * @brief  tests readers
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "scip/cons_expr.c"
#include <libgen.h>
#include "include/scip_test.h"

Test(readers, pip)
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR** children;
   SCIP_CONSEXPR_EXPR** grandchildren;
   char filename[SCIP_MAXSTRLEN];
   SCIP_CONSEXPR_EXPRHDLR* sumhdlr;

   /* get file to read: test.mps that lives in the same directory as this file */
   (void)SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s", __FILE__);
   dirname(filename);
   strcat(filename, "/test.pip");
   printf("Reading %s\n", filename);

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPreadProb(scip, filename, NULL));

   /* check that vars are what we expect */
   cr_expect_eq(SCIPgetNVars(scip), 5);
   vars = SCIPgetVars(scip);

   cr_expect_str_eq(SCIPvarGetName(vars[0]), "x");
   cr_expect_str_eq(SCIPvarGetName(vars[1]), "y");
   cr_expect_str_eq(SCIPvarGetName(vars[2]), "z");
   cr_expect_str_eq(SCIPvarGetName(vars[3]), "objconst");
   cr_expect_str_eq(SCIPvarGetName(vars[4]), "nonlinobjvar");

   /* check that cons are what we expect */
   cr_expect_eq(SCIPgetNConss(scip), 3);
   conss = SCIPgetConss(scip);

   cr_expect_str_eq(SCIPconsGetName(conss[0]), "nonlinobj");
   cr_expect_str_eq(SCIPconsGetName(conss[1]), "e1");
   cr_expect_str_eq(SCIPconsGetName(conss[2]), "e2");

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   sumhdlr = SCIPgetConsExprExprHdlrSum(conshdlr);
   cr_assert_not_null(sumhdlr);

   /* check objective coefficients */
   cr_expect_eq(SCIPvarGetObj(vars[0]), 0.0);
   cr_expect_eq(SCIPvarGetObj(vars[1]), 0.0);
   cr_expect_eq(SCIPvarGetObj(vars[2]), 0.0);
   cr_expect_eq(SCIPvarGetObj(vars[3]), 11.0);
   cr_expect_eq(SCIPvarGetObj(vars[4]), 1.0);

   /*
    * check whether the first constraint is nonlinear and of the form + x - y^3 z^2 - nonlinobj <= 0
    */

   cr_expect_eq(SCIPconsGetHdlr(conss[0]), SCIPfindConshdlr(scip, "expr"));
   expr = SCIPgetExprConsExpr(scip, conss[0]);
   cr_expect_not_null(expr);
   children = SCIPgetConsExprExprChildren(expr);
   cr_expect_not_null(children);
   cr_expect_eq(SCIPgetConsExprExprNChildren(expr), 3);

   /* check sides */
   cr_expect(SCIPisInfinity(scip, -SCIPgetLhsConsExpr(scip, conss[0])));
   cr_expect_eq(SCIPgetRhsConsExpr(scip, conss[0]), 0.0);

   /* check constant and coefficients of the sum expression */
   cr_expect_eq(SCIPgetConsExprExprSumConstant(expr), 0.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[0], 1.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[1], -1.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[2], -1.0);

   /* check whether first and third child is a variable */
   cr_expect(SCIPisConsExprExprVar(children[0]));
   cr_expect(SCIPisConsExprExprVar(children[2]));

   /* check whether second child is a product; both grandchildren need to be power expressions */
   cr_expect(SCIPgetConsExprExprHdlr(children[1]) == SCIPfindConsExprExprHdlr(conshdlr, "prod"));
   grandchildren = SCIPgetConsExprExprChildren(children[1]);
   cr_expect_eq(SCIPgetConsExprExprNChildren(children[1]), 2);
   cr_expect(SCIPgetConsExprExprHdlr(grandchildren[0]) == SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   cr_expect(SCIPgetConsExprExprHdlr(grandchildren[1]) == SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   cr_expect_eq(SCIPgetConsExprExprPowExponent(grandchildren[0]), 3.0);
   cr_expect_eq(SCIPgetConsExprExprPowExponent(grandchildren[1]), 2.0);

   /*
    * check whether the second constraint is linear and of the form x + 2y + 3z <= 1
    */

   cr_expect_eq(SCIPconsGetHdlr(conss[1]), SCIPfindConshdlr(scip, "linear"));
   cr_expect(SCIPisInfinity(scip, -SCIPgetLhsLinear(scip, conss[1])));
   cr_expect_eq(SCIPgetRhsLinear(scip, conss[1]), 1.0);
   cr_expect_eq(SCIPgetNVarsLinear(scip, conss[1]), 3);
   cr_expect_eq(SCIPgetValsLinear(scip, conss[1])[0], 1.0);
   cr_expect_eq(SCIPgetValsLinear(scip, conss[1])[1], 2.0);
   cr_expect_eq(SCIPgetValsLinear(scip, conss[1])[2], 3.0);

   /*
    * check whether the third constraint is nonlinear and of the form x^2 y^3 z^4 + x + y + 2 = 10
    */

   cr_expect_eq(SCIPconsGetHdlr(conss[2]), SCIPfindConshdlr(scip, "expr"));
   expr = SCIPgetExprConsExpr(scip, conss[2]);
   cr_expect_not_null(expr);
   children = SCIPgetConsExprExprChildren(expr);
   cr_expect_not_null(children);
   cr_expect_eq(SCIPgetConsExprExprNChildren(expr), 3);

   /* check sides */
   cr_expect_eq(SCIPgetLhsConsExpr(scip, conss[2]), 10.0);
   cr_expect_eq(SCIPgetRhsConsExpr(scip, conss[2]), 10.0);

   /* check constant and coefficients of the sum expression */
   cr_expect_eq(SCIPgetConsExprExprSumConstant(expr), 7.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[0], 4.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[1], 5.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[2], 6.0);

   /* check whether first child is a product; all three grandchildren need to be power expressions */
   cr_expect(SCIPgetConsExprExprHdlr(children[0]) == SCIPfindConsExprExprHdlr(conshdlr, "prod"));
   grandchildren = SCIPgetConsExprExprChildren(children[0]);
   cr_expect_eq(SCIPgetConsExprExprNChildren(children[0]), 3);
   cr_expect(SCIPgetConsExprExprHdlr(grandchildren[0]) == SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   cr_expect(SCIPgetConsExprExprHdlr(grandchildren[1]) == SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   cr_expect(SCIPgetConsExprExprHdlr(grandchildren[2]) == SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   cr_expect_eq(SCIPgetConsExprExprPowExponent(grandchildren[0]), 2.0);
   cr_expect_eq(SCIPgetConsExprExprPowExponent(grandchildren[1]), 3.0);
   cr_expect_eq(SCIPgetConsExprExprPowExponent(grandchildren[2]), 4.0);
}
