/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   propagate.c
 * @brief  unit test for propagation of cons_expr
 *
 * @note: to call propagation methods we need active constraints, so we go to presolving stage
 * and create the constraint there by parsing an expression (using the transformed variables!)
 * or by explicitly constructing the tree
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.c"

#include "include/scip_test.h"

#define CHECK_EXPRINTERVAL(scip,expr,a,b) (SCIPisFeasEQ((scip), SCIPgetConsExprExprInterval((expr)).inf, (a)) && SCIPisFeasEQ((scip), SCIPgetConsExprExprInterval((expr)).sup, (b)))
#define EXPECTING_EXPRINTERVAL(expr,a,b) (a), (b), SCIPgetConsExprExprInterval((expr)).inf, SCIPgetConsExprExprInterval((expr)).sup

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

static
void setup(void)
{
   static SCIP_VAR* xo;
   static SCIP_VAR* yo;
   static SCIP_VAR* zo;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &xo, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yo, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zo, "z", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, xo) );
   SCIP_CALL( SCIPaddVar(scip, yo) );
   SCIP_CALL( SCIPaddVar(scip, zo) );

   /* goto presolving */
   TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, FALSE);

   /* get transformed vars and release vars */
   SCIP_CALL( SCIPgetTransformedVar(scip, xo, &x) );
   SCIP_CALL( SCIPgetTransformedVar(scip, yo, &y) );
   SCIP_CALL( SCIPgetTransformedVar(scip, zo, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &xo) );
   SCIP_CALL( SCIPreleaseVar(scip, &yo) );
   SCIP_CALL( SCIPreleaseVar(scip, &zo) );
}

static
void teardown(void)
{
   /* free scip and check for memory leaks */
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There are memory leaks!");
}

TestSuite(propagate, .init = setup, .fini = teardown);

Test(propagate, sum)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* originalexpr;
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   int ntightenings;

   SCIP_CALL( SCIPchgVarLb(scip, x, -2.0) ); SCIP_CALL( SCIPchgVarUb(scip, x, 2.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -3.0) ); SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* create cons 0.5 <= 2x -y + 0.5 <= 1.5*/
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "2*<t_x>[C]-<t_y>[C] + 0.5", NULL, &originalexpr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, originalexpr, &expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &originalexpr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, 0.5, 1.5) );

   /* forward prop only propagates active constraints */
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   cr_assert_not(redundant);

   SCIP_CALL( reversePropConss(scip, &cons, 1, FALSE, &infeasible, &ntightenings) );
   cr_assert_not(infeasible);

   cr_expect(SCIPisEQ(scip, expr->interval.inf, 0.5), "Expecting 0.5, got %g\n",  expr->interval.inf);
   cr_expect(SCIPisEQ(scip, expr->interval.sup, 1.5));
   cr_expect(SCIPisEQ(scip, expr->children[0]->interval.inf, -1.5), "Expecting -1.5, got %g\n", expr->children[0]->interval.inf);
   cr_expect(SCIPisEQ(scip, expr->children[0]->interval.sup, 1.0));
   cr_expect(SCIPisFeasEQ(scip, expr->children[1]->interval.inf, -3.0), "Expecting -3.0, got %.20f\n", expr->children[1]->interval.inf);
   cr_expect(SCIPisFeasEQ(scip, expr->children[1]->interval.sup, 1.0));

   /* release stuff */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

Test(propagate, product)
{
   SCIP_CONSEXPR_EXPR* expr, *expraux;
   SCIP_CONSEXPR_EXPR* originalexpr;
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   int ntightenings;

   /* change bounds of vars */
   SCIP_CALL( SCIPchgVarLb(scip, x, 1.0) ); SCIP_CALL( SCIPchgVarUb(scip, x, 3.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 2.0) ); SCIP_CALL( SCIPchgVarUb(scip, y, 4.0) );

   /* create cons 0.0 <= 0.5x^2/y <= 1 */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "0.5*<t_x>^2/<t_y>", NULL, &originalexpr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, originalexpr, &expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &originalexpr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, 0.0, 1.0) );

   /* forward prop only propagates active constraints and active constraint live in the transformed problem */
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   cr_assert_not(redundant);

   SCIP_CALL( reversePropConss(scip, &cons, 1, FALSE, &infeasible, &ntightenings) );
   cr_assert_not(infeasible);

   /* test stuff */
   cr_expect(SCIPisEQ(scip, expr->interval.inf, 1.0 / 8.0), "Expecting %g and got %g\n", 1.0/8.0, expr->interval.inf);
   cr_expect(SCIPisEQ(scip, expr->interval.sup, 1.0));

   expraux = expr->children[0];
   cr_expect(SCIPisEQ(scip, expraux->interval.inf, 1.0 / 4.0), "Expecting %g and got %g\n", 1.0, expraux->interval.inf);
   cr_expect(SCIPisEQ(scip, expraux->interval.sup, 2.0));

   cr_expect(SCIPisFeasEQ(scip, expraux->children[0]->interval.inf, 1.0), "Expecting %g and got %g\n", 1.0, expraux->children[0]->interval.inf);
   cr_expect(SCIPisFeasEQ(scip, expraux->children[0]->interval.sup, 8.0));
   cr_expect(SCIPisFeasEQ(scip, expraux->children[0]->children[0]->interval.inf, 1.0), "Expecting %g and got %g\n", 1.0, expraux->children[0]->interval.inf);
   cr_expect(SCIPisFeasEQ(scip, expraux->children[0]->children[0]->interval.sup, SQRT(8)));
   cr_expect(SCIPisFeasEQ(scip, expraux->children[1]->interval.inf, 1/4.0));
   cr_expect(SCIPisFeasEQ(scip, expraux->children[1]->interval.sup, 1/2.0));
   cr_expect(SCIPisFeasEQ(scip, expraux->children[1]->children[0]->interval.inf, 2.0));
   cr_expect(SCIPisFeasEQ(scip, expraux->children[1]->children[0]->interval.sup, 4.0));

   /* release conss */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

Test(propagate, abs)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* originalexpr;
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   int ntightenings;

   /* change bounds of vars */
   SCIP_CALL( SCIPchgVarLb(scip, x, -3.0) ); SCIP_CALL( SCIPchgVarUb(scip, x, 4.0) );

   /* create cons 1.0 <= |x| <= 2.5 */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "abs(<t_x>)", NULL, &originalexpr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, originalexpr, &expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &originalexpr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, 1.0, 2.5) );

   /* forward prop only propagates active constraints and active constraint live in the transformed problem */
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   cr_assert_not(redundant);

   SCIP_CALL( reversePropConss(scip, &cons, 1, FALSE, &infeasible, &ntightenings) );
   cr_assert_not(infeasible);

   /* get expression and test stuff */
   cr_expect(SCIPisEQ(scip, expr->children[0]->interval.inf, -2.5));
   cr_expect(SCIPisEQ(scip, expr->children[0]->interval.sup, 2.5));

   /* release conss */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

Test(propagate, exp)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* originalexpr;
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   int ntightenings;

   /* change bounds of vars */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) ); SCIP_CALL( SCIPchgVarUb(scip, x, 3.0) );

   /* create cons -1 <= exp(x) <= 2 */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "exp(<t_x>[C])", NULL, &originalexpr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, originalexpr, &expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &originalexpr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -1.0, 2.0) );

   /* forward prop only propagates active constraints and active constraint live in the transformed problem */
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   cr_assert_not(redundant);

   SCIP_CALL( reversePropConss(scip, &cons, 1, FALSE, &infeasible, &ntightenings) );
   cr_assert_not(infeasible);

   /* get expression and test stuff */
   cr_expect(SCIPisEQ(scip, expr->interval.inf, exp(-1)));
   cr_expect(SCIPisEQ(scip, expr->interval.sup, 2.0));
   cr_expect(SCIPisFeasEQ(scip, expr->children[0]->interval.inf, -1.0));
   cr_expect(SCIPisEQ(scip, expr->children[0]->interval.sup, log(2.0)));

   /* release conss */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

Test(propagate, log)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* originalexpr;
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   int ntightenings;

   /* change bounds of vars */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) ); SCIP_CALL( SCIPchgVarUb(scip, x, 7.0) );

   /* create cons -1  <= log(x) <= 1 */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "log(<t_x>[C])", NULL, &originalexpr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, originalexpr, &expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &originalexpr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -1.0, 1.0) );

   /* forward prop only propagates active constraints and active constraint live in the transformed problem */
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   cr_assert_not(redundant);

   SCIP_CALL( reversePropConss(scip, &cons, 1, FALSE, &infeasible, &ntightenings) );
   cr_assert_not(infeasible);

   /* get expression and test stuff */
   cr_expect(SCIPisEQ(scip, expr->interval.inf, -1.0));
   cr_expect(SCIPisEQ(scip, expr->interval.sup, 1.0));
   cr_expect(SCIPisEQ(scip, expr->children[0]->interval.inf, exp(-1.0)));
   cr_expect(SCIPisEQ(scip, expr->children[0]->interval.sup, exp(1.0)));

   /* release conss */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

Test(propagate, complicated_expression)
{
   SCIP_CONSEXPR_EXPR* xexpr;
   SCIP_CONSEXPR_EXPR* yexpr;
   SCIP_CONSEXPR_EXPR* zexpr;
   SCIP_CONSEXPR_EXPR* zinvexpr;
   SCIP_CONSEXPR_EXPR* rootexpr;
   SCIP_CONSEXPR_EXPR* powexpr;
   SCIP_CONSEXPR_EXPR* sumexpr;
   SCIP_CONSEXPR_EXPR* logexpr;
   SCIP_CONSEXPR_EXPR* exprs[2];
   SCIP_Real coeffs[2];
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   int ntightenings;

   /* change bounds of vars */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) ); SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 2.0) ); SCIP_CALL( SCIPchgVarUb(scip, y, 3.0) );
   SCIP_CALL( SCIPchgVarLb(scip, z, 1.0) ); SCIP_CALL( SCIPchgVarUb(scip, z, 2.0) );

   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &zexpr, z) );

   /*
    * create constraint 0 <= (-x^2 + log(y)) / z <= 2
    */

   /* log(y) */
   SCIP_CALL( SCIPcreateConsExprExprLog(scip, conshdlr, &logexpr, yexpr) );

   /* x^2 */
   coeffs[0] = 2.0;
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexpr, xexpr, 2.0) );

   /* log(y) - x^2 */
   coeffs[0] = 1.0;
   exprs[0] = logexpr;
   coeffs[1] = -1.0;
   exprs[1] = powexpr;
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr, 2, exprs, coeffs, 0.0) );

   /* z^-1 */
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &zinvexpr, zexpr, -1.0) );

   /* (-x^2 + log(y)) / z */
   exprs[0] = sumexpr;
   exprs[1] = zinvexpr;
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &rootexpr, 2, exprs, 1.0) );

   SCIPinfoMessage(scip, NULL, "test more complicated expression: ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, rootexpr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* forward prop only propagates active constraints and active constraint live in the transformed problem */
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", rootexpr, 0.0, 2.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   /* apply forward propagation */
   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   cr_expect(CHECK_EXPRINTERVAL(scip, xexpr, -1.0, 1.0), "expecting [%g, %g], got [%g, %g]\n", EXPECTING_EXPRINTERVAL(xexpr,-1.0,1.0));
   cr_expect(CHECK_EXPRINTERVAL(scip, yexpr, 2.0, 3.0), "expecting [%g, %g], got [%g, %g]\n", EXPECTING_EXPRINTERVAL(yexpr,2.0,3.0));
   cr_expect(CHECK_EXPRINTERVAL(scip, zexpr, 1.0, 2.0), "expecting [%g, %g], got [%g, %g]\n", EXPECTING_EXPRINTERVAL(zexpr,1.0,2.0));
   cr_expect(CHECK_EXPRINTERVAL(scip, logexpr, log(2), log(3)), "expecting [%g, %g], got [%g, %g]\n", EXPECTING_EXPRINTERVAL(logexpr,log(2),log(3)));
   cr_expect(CHECK_EXPRINTERVAL(scip, powexpr, 0.0, 1.0), "expecting [%g, %g], got [%g, %g]\n", EXPECTING_EXPRINTERVAL(powexpr,0.0,1.0));
   cr_expect(CHECK_EXPRINTERVAL(scip, sumexpr, log(2) - 1, log(3)), "expecting [%g, %g], got [%g, %g]\n", EXPECTING_EXPRINTERVAL(sumexpr,log(2)-1.0,log(3)));
   cr_expect(CHECK_EXPRINTERVAL(scip, rootexpr, 0.0, log(3)), "expecting [%g, %g], got [%g, %g]\n", EXPECTING_EXPRINTERVAL(rootexpr,0.0,log(3)));

   /* apply reverse propagation */
   SCIP_CALL( reversePropConss(scip, &cons, 1, FALSE, &infeasible, &ntightenings) );
   cr_assert_not(infeasible);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &rootexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &logexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &zinvexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &zexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPfreeTransform(scip) );
}

Test(propagate, unbounded_sub_expression)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   int ntightenings;

   /* change bounds of vars */
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 0.0) );
   SCIP_CALL( SCIPchgVarLb(scip, z, 1.0) ); SCIP_CALL( SCIPchgVarUb(scip, z, 2.0) );

   /*
    *  -5 <= x * y
    */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "2 * <t_x> * <t_y>", NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -5.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   SCIPinfoMessage(scip, NULL, "test expressions with unbounded sub-expression: ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   cr_assert(CHECK_EXPRINTERVAL(scip, expr, -5.0, SCIPinfinity(scip)));

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /*
    *  -inf <= 2 + x - y <= inf
    */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "2 + <t_x> - <t_y>", NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -SCIPinfinity(scip), SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   SCIPinfoMessage(scip, NULL, "test expressions with unbounded sub-expression: ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   cr_assert(CHECK_EXPRINTERVAL(scip, expr, -SCIPinfinity(scip), SCIPinfinity(scip)));

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /*
    *  -inf <= x * y * z <= inf
    */
   SCIP_CALL( SCIPchgVarUb(scip, x, 0.0) );

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x> * <t_y> * <t_z>", NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -SCIPinfinity(scip), SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   SCIPinfoMessage(scip, NULL, "test expressions with unbounded sub-expression: ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   cr_assert(CHECK_EXPRINTERVAL(scip, expr, 0.0, SCIPinfinity(scip)));

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(propagate, forwardprop_uses_expressions_bounds)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   int ntightenings;

   /* change bounds of vars */
   SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) ); SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 0.0) ); SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x> + <t_y>", NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -1.0, 0.5) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   cr_assert(CHECK_EXPRINTERVAL(scip, expr, 0.0, 0.5));

   /* change intervals of variable expressions */
   expr->children[0]->interval.inf = -1.0; expr->children[0]->interval.sup = 0.2;
   expr->children[1]->interval.inf = -1.0; expr->children[1]->interval.sup = 0.2;

   /* new interval should be [0,1] intersected with [-2, 0.4] */
   SCIP_CALL( forwardPropCons(scip, cons, TRUE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   cr_assert(CHECK_EXPRINTERVAL(scip, expr, 0.0, 0.4));

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(propagate, infeas_after_forwardprop)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   int ntightenings;

   /* change bounds of vars */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) ); SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 0.0) ); SCIP_CALL( SCIPchgVarUb(scip, y, 3.0) );

   /*
    * -7.0 <= 2 * x * y <= -6.1
    */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "2 * <t_x> * <t_y>", NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -7.0, -6.1) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert(infeasible);
   cr_assert(SCIPintervalIsEmpty(SCIPinfinity(scip), SCIPgetConsExprExprInterval(expr)));

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /*
    * -1.0 <= 1 + x / (1 + y) <= -0.1
    */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "1 + <t_x> / (1 + <t_y>)", NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -1.0, -0.1) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert(SCIPintervalIsEmpty(SCIPinfinity(scip), SCIPgetConsExprExprInterval(expr)));
   cr_assert(infeasible);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /*
    * 0.0 <= 1 + expr(-5 * x + y^2) <= 0.9
    */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "1 + exp(-5 * <t_x> + <t_y>^2)", NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -10.0, -0.1) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   cr_assert(SCIPconsIsActive(cons));

   SCIP_CALL( forwardPropCons(scip, cons, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert(SCIPintervalIsEmpty(SCIPinfinity(scip), SCIPgetConsExprExprInterval(expr)));
   cr_assert(infeasible);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(propagate, infeas_after_backwardprop)
{
   SCIP_CONSEXPR_EXPR* xexpr, *yexpr, *exprs[2];
   SCIP_CONSEXPR_EXPR* expr1, *expr2;
   SCIP_CONS* cons1, *cons2;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   int ntightenings;

   /* change bounds of vars */
   SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) ); SCIP_CALL( SCIPchgVarUb(scip, x, 2.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 0.0) ); SCIP_CALL( SCIPchgVarUb(scip, y, 2.0) );

   /*
    * -1.0 <= x * y <= 1.0 and 3.5 <= x + y <= 5.0
    */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
   exprs[0] = xexpr;
   exprs[1] = yexpr;

   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr1, 2, exprs, 1.0) );
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr2, 2, exprs, NULL, 0.0) );

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons1, "cons1", expr1, -1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons1) );
   cr_assert(SCIPconsIsActive(cons1));

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons2, "cons2", expr2, 3.5, 5.0) );
   SCIP_CALL( SCIPaddCons(scip, cons2) );
   cr_assert(SCIPconsIsActive(cons2));

   /* apply forward propagation for both constraints */
   SCIP_CALL( forwardPropCons(scip, cons1, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   cr_assert(CHECK_EXPRINTERVAL(scip, expr1, 0.0, 1.0));

   SCIP_CALL( forwardPropCons(scip, cons2, FALSE, FALSE, &infeasible, &redundant, &ntightenings) );
   cr_assert_not(infeasible);
   assert(CHECK_EXPRINTERVAL(scip, expr2, 3.5, 4.0));

   /* reverse propagation of cons2 should lead to new bounds on x and y */
   SCIP_CALL( reversePropConss(scip, &cons2, 1, FALSE, &infeasible, &ntightenings) );
   cr_assert_not(infeasible);
   cr_assert(CHECK_EXPRINTERVAL(scip, expr2->children[0], 1.5, 2.0));
   cr_assert(CHECK_EXPRINTERVAL(scip, expr2->children[1], 1.5, 2.0));

   /* reverse propagation of cons1 should lead to an empty interval for x */
   SCIP_CALL( reversePropConss(scip, &cons1, 1, FALSE, &infeasible, &ntightenings) );
   cr_assert(infeasible);
   cr_assert(SCIPintervalIsEmpty(SCIPinfinity(scip), expr1->children[0]->interval));

   SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr1) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
}

struct expr_results
{
   const char* cons1;
   const char* cons2;
   SCIP_Real xlb;
   SCIP_Real xub;
};

ParameterizedTestParameters(propagate, propConss)
{
   static const struct expr_results data[] =
   {
      {"<t_x>^2 + <t_x>", "<t_x>^2 - 1.0", -1.0, -1.0},
      {"<t_x>^(0.5) - <t_y>","<t_x> - 1.0 - <t_y>", 2.618033988749895, 2.618033988749895},
      {"exp(<t_x>) - <t_y>", "4.0 * <t_x>^(1.5) - <t_y>", 0.58687228932071, 3.06767359040726},
      {"log(abs(<t_x> + 1)) - <t_y>", "abs(<t_x>)^1.5 - <t_y>", 0.0,  0.6096527513}
   };
   return cr_make_param_array(const struct expr_results, data, sizeof(data)/sizeof(struct expr_results));
}

ParameterizedTest(const struct expr_results* data, propagate, propConss)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONS* cons1, *cons2;
   int nchgbds = 0;
   int ndelconss = 0;
   SCIP_RESULT result;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -10.0) ); SCIP_CALL( SCIPchgVarUb(scip, x, 10.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -100.0) ); SCIP_CALL( SCIPchgVarUb(scip, y, 100.0) );
   SCIP_CALL( SCIPchgVarLb(scip, z, -3.0) ); SCIP_CALL( SCIPchgVarUb(scip, z, 1.0) );

   SCIPinfoMessage(scip, NULL, "test constraints: %s == 0 and %s == 0 \n", data->cons1, data->cons2);

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, data->cons1, NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons1, "cons1", expr, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons1) );
   cr_assert(SCIPconsIsActive(cons1));
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, data->cons2, NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons2, "cons2", expr, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons2) );
   cr_assert(SCIPconsIsActive(cons1));
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* call propConss() for all transformed constraints */
   cr_assert_not_null(SCIPconshdlrGetConss(conshdlr));
   cr_assert_eq(SCIPconshdlrGetNConss(conshdlr), 2);
   SCIP_CALL( propConss(scip, conshdlr, SCIPconshdlrGetConss(conshdlr), SCIPconshdlrGetNConss(conshdlr), FALSE, &result, &nchgbds, &ndelconss) );
   cr_assert_eq(result, SCIP_REDUCEDDOM, "expecting %d, but got %d\n", SCIP_REDUCEDDOM, result);
   cr_assert_gt(nchgbds, 0);

   /* check bounds */
   cr_expect(SCIPvarGetLbLocal(x) <= data->xlb);
   cr_expect(SCIPvarGetUbLocal(x) >= data->xub);
   /* @benny: the tolerance is pretty bad for some tests, can you have a look please? */
   cr_expect(REALABS(SCIPvarGetLbLocal(x) - data->xlb) <= 0.1, "Expecting %g, got %g\n", data->xlb, SCIPvarGetLbLocal(x));
   cr_expect(REALABS(SCIPvarGetUbLocal(x) - data->xub) <= 0.2, "Expecting %g, got %g\n", data->xub, SCIPvarGetUbLocal(x));

   /* free transformed problem and remove constraints */
   SCIP_CALL( SCIPdelCons(scip, cons1) );
   SCIP_CALL( SCIPdelCons(scip, cons2) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
}
