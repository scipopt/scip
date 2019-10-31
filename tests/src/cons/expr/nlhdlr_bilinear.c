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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_bilinear.c
 * @brief  tests bilinear nonlinear handler methods
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_nlhdlr_bilinear.c"
#include "scip/struct_cons_expr.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_CONSEXPR_NLHDLR* nlhdlr;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   nlhdlr = SCIPfindConsExprNlhdlr(conshdlr, "bilinear");
   cr_assert(nlhdlr != NULL);
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

/* define the test suite */
TestSuite(nlhdlrbilinear, .init = setup, .fini = teardown);

/** creates a constraint and calls detection methods of nonlinear handlers */
static
SCIP_RETCODE createAndDetect(
   SCIP_CONS**           cons,              /**< pointer to store the constraint */
   SCIP_CONSEXPR_EXPR*   rootexpr           /**< root expression of the constraint */
   )
{
   SCIP_Bool infeasible;

   SCIP_CALL( SCIPcreateConsExprBasic(scip, cons, "cons", rootexpr, -100.0, 100.0) );
   SCIP_CALL( SCIPaddConsLocks(scip, *cons, 1, 0) );
   SCIP_CALL( SCIPinitlpCons(scip, *cons, &infeasible) );
   SCIP_CALL( SCIPclearCuts(scip) ); /* we have to clear the separation store */
   SCIP_CALL( SCIPaddConsLocks(scip, *cons, -1, 0) );

   return SCIP_OKAY;
}

/** free enfo data in expression
 * used when we do a detect on a constraint that wasn't added to SCIP, so the conshdlr will not do this for us
 */
static
SCIP_RETCODE freeEnfoData(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression where to free enfo data */
   )
{
   int e;

   /* free data stored by nonlinear handlers */
   for( e = 0; e < expr->nenfos; ++e )
   {
      SCIP_CONSEXPR_NLHDLR* nlhdlr_;

      assert(expr->enfos[e] != NULL);

      nlhdlr_ = expr->enfos[e]->nlhdlr;
      assert(nlhdlr_ != NULL);

      if( expr->enfos[e]->issepainit )
      {
         /* call the separation deinitialization callback of the nonlinear handler */
         SCIP_CALL( SCIPexitsepaConsExprNlhdlr(scip, nlhdlr_, expr, expr->enfos[e]->nlhdlrexprdata) );
         expr->enfos[e]->issepainit = FALSE;
      }

      /* free nlhdlr exprdata, if there is any and there is a method to free this data */
      if( expr->enfos[e]->nlhdlrexprdata != NULL && nlhdlr_->freeexprdata != NULL )
      {
         SCIP_CALL( (*nlhdlr_->freeexprdata)(scip, nlhdlr_, expr, &expr->enfos[e]->nlhdlrexprdata) );
         assert(expr->enfos[e]->nlhdlrexprdata == NULL);
      }

      /* free enfo data */
      SCIPfreeBlockMemory(scip, &expr->enfos[e]); /*lint !e866 */
   }

   /* free array with enfo data */
   SCIPfreeBlockMemoryArrayNull(scip, &expr->enfos, expr->nenfos);
   expr->nenfos = 0;

   return SCIP_OKAY;
}

/*
 * tests
 */

/* creates and adds two expression constraints and check output of SCIPgetBilinearExprsExpr */
Test(nlhdlrbilinear, collect_product_expressions)
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_EXPR** exprs;
   SCIP_CONSEXPR_EXPR* expr;
   const char* inputs[2] = {"<t_x> * <t_y> + <t_y> * <t_z>", "exp(<t_x>) * sin(<t_y>)"};
   int nexprs;
   int c;

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraints; call INITLP manually before adding each constraint to SCIP */
   for( c = 0; c < 2; ++c )
   {
      SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, inputs[c], NULL, &expr) );
      SCIP_CALL( createAndDetect(&cons, expr) );

      /* add and release constraint */
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* check whether all product expressions could be found */
   exprs = SCIPgetConsExprNlhdlrBilinearExprs(nlhdlr);
   cr_expect(exprs != NULL);
   nexprs = SCIPgetConsExprNlhdlrBilinearNExprs(nlhdlr);
   cr_expect(nexprs == 3);

   for( c = 0; c < nexprs; ++c )
   {
      SCIP_CONSEXPR_EXPR* child1;
      SCIP_CONSEXPR_EXPR* child2;

      cr_assert(exprs[c] != NULL);
      cr_expect(SCIPgetConsExprExprNChildren(exprs[c]) == 2);
      cr_expect(SCIPgetConsExprExprAuxVar(exprs[c]) != NULL);
      cr_expect(SCIPgetConsExprExprHdlr(exprs[c]) == SCIPgetConsExprExprHdlrProduct(conshdlr));

      child1 = SCIPgetConsExprExprChildren(exprs[c])[0];
      cr_assert(child1 != NULL);
      cr_expect(SCIPgetConsExprExprAuxVar(child1) != NULL);
      child2 = SCIPgetConsExprExprChildren(exprs[c])[1];
      cr_assert(child2 != NULL);
      cr_expect(SCIPgetConsExprExprAuxVar(child2) != NULL);

      /* check for x*y */
      if( SCIPgetConsExprExprAuxVar(child1) == SCIPvarGetTransVar(x) )
         cr_expect(SCIPgetConsExprExprAuxVar(child2) == SCIPvarGetTransVar(y));

      /* check for y*z */
      else if( SCIPgetConsExprExprAuxVar(child1) == SCIPvarGetTransVar(y) )
         cr_expect(SCIPgetConsExprExprAuxVar(child2) == SCIPvarGetTransVar(z));

      /* check for exp(x) * sin(y) */
      else
      {
         cr_expect(SCIPgetConsExprExprHdlr(child1) == SCIPfindConsExprExprHdlr(conshdlr, "exp"));
         cr_expect(SCIPgetConsExprExprHdlr(child2) == SCIPfindConsExprExprHdlr(conshdlr, "sin"));
      }
   }
}

/* adds valid inequalities to a product expression */
Test(nlhdlrbilinear, add_inequality)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_CONSEXPR_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, 2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 4.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x> * <t_y>", NULL, &expr) );
   SCIP_CALL( createAndDetect(&cons, expr) );


   /*
    * add inequalities to product expression
    */

   /* separating inequality */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, 1.0, -2.0, 4.5, &success) );
   cr_expect(success);

   /* non-separating inequality */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, 1.0, -1.0, 5.0, &success) );
   cr_expect(!success);

   /* duplicates should not be accepted */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, 1.0, -2.0, 4.5, &success) );
   cr_expect(!success);

   /* first inequality is still dominating this inequality -> do not accept it */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, 0.9, -3.0, 5.5, &success) );
   cr_expect(!success);

   /* separates bottom-left corner */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, -1.0, -1.1, -2.0, &success) );
   cr_expect(success);

   /* dominates first inequality */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, 1.1, 0.5, -3.0, &success) );
   cr_expect(success);

   /* worse than inequality for bottom-left corner */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, -1.0, -1.1, -2.2, &success) );
   cr_expect(!success);

   /* free memory */
   SCIP_CALL( freeEnfoData(expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/** auxiliary function to check coefficients of a cut */
static
SCIP_RETCODE checkCut(
   SCIP_ROW*            row,                /**< row */
   SCIP_Real            coefx,              /**< target coefficient for x */
   SCIP_Real            coefy,              /**< target coefficient for y */
   SCIP_Real            coefz,              /**< target coefficient for the auxiliary variable */
   SCIP_Real            constant            /**< target constant */
   )
{
   SCIP_VAR* tx = SCIPvarGetTransVar(x);
   SCIP_VAR* ty = SCIPvarGetTransVar(y);
   int i;

   assert(row != NULL);

   for( i = 0; i < SCIProwGetNNonz(row); ++i )
   {
      SCIP_VAR* var = SCIPcolGetVar(SCIProwGetCols(row)[i]);
      SCIP_Real coef = SCIProwGetVals(row)[i];

      cr_assert(var != NULL);

      if( var == tx )
         assert(SCIPisEQ(scip, coefx, coef));
      else if( var == ty )
         cr_expect(SCIPisEQ(scip, coefy, coef));
      else
         cr_expect(SCIPisEQ(scip, coefz, coef));

      if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
         cr_expect(SCIPisEQ(scip, -SCIProwGetLhs(row), constant));
      else
         cr_expect(SCIPisEQ(scip, -SCIProwGetRhs(row), constant));
   }

   return SCIP_OKAY;
}

/* adds a linear inequality to the product expression and computes a cut at a given reference point */
Test(nlhdlrbilinear, separation_single)
{
   SCIP_ROWPREP* rowprep;
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_ROW* row;
   SCIP_Bool success;
   SCIP_Bool overestimate;
   SCIP_Bool dummy;
   SCIP_CONSEXPR_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 2.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x> * <t_y>", NULL, &expr) );
   SCIP_CALL( createAndDetect(&cons, expr) );


   /* INITLP should have added an auxiliary variable to the product expression (tight might change in the future) */
   cr_assert( SCIPgetConsExprExprAuxVar(expr) != NULL);

   /* create a solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(x), 1.1214285714285714) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(y), -0.1785714285714286) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr), -1.25) );

    /* add inequality for underestimation to product expression */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, 1.0, 1.0, 1.3, &success) );
   cr_expect(success);

   /* auxvar = -1.25 => McCormick relaxation is violated so the resulting cut should be one of the McCormick inequalities */
   overestimate = FALSE;
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, FALSE) );
   SCIP_CALL( nlhdlrEstimateBilinear(scip, conshdlr, nlhdlr, expr, SCIPgetConsExprNlhdlrExprData(nlhdlr, expr), sol, 0.0, overestimate, 0.0, rowprep, &success, FALSE, &dummy) );
   cr_expect(!success);
   SCIPfreeRowprep(scip, &rowprep);

   /* auxvar = -1.23 => point is in the convex hull of McCormick and thus we should use the bilinear inequality */
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr), -1.23) );
   overestimate = FALSE;
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, FALSE) );
   SCIP_CALL( nlhdlrEstimateBilinear(scip, conshdlr, nlhdlr, expr, SCIPgetConsExprNlhdlrExprData(nlhdlr, expr), sol, 0.0, overestimate, 0.0, rowprep, &success, FALSE, &dummy) );
   cr_expect(success);
   SCIP_CALL( SCIPgetRowprepRowConshdlr(scip, &row, rowprep, conshdlr) );
   SCIP_CALL( checkCut(row, 0.5790816326530611, 0.3637755102040817, -1.0, -0.7846938775510204) );
   SCIP_CALL( SCIPreleaseRow(scip, &row) );
   SCIPfreeRowprep(scip, &rowprep);

   /* free memory */
   SCIP_CALL( freeEnfoData(expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
}

/* adds two linear inequality to the product expression and computes a cut at a given reference point */
Test(nlhdlrbilinear, separation_two)
{
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_ROWPREP* rowprep;
   SCIP_ROW* row;
   SCIP_Bool success;
   SCIP_Bool overestimate;
   SCIP_Bool dummy;
   SCIP_CONSEXPR_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 2.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 3.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x> * <t_y>", NULL, &expr) );
   SCIP_CALL( createAndDetect(&cons, expr) );

   /* INITLP should have added an auxiliary variable to the product expression (tight might change in the future) */
   cr_assert( SCIPgetConsExprExprAuxVar(expr) != NULL);

   /* create a solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(x), 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(y), 1.5) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr), 4.4) );

   /* add two inequalities to the product expression */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, -1.2, -0.6, 0.9, &success) );
   cr_expect(success);
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, 0.5, 0.3, 0.6, &success) );
   cr_expect(success);

   /* auxvar = -0.1 => McCormick should do the job */
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr), -0.1) );
   overestimate = FALSE;
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, FALSE) );
   SCIP_CALL( nlhdlrEstimateBilinear(scip, conshdlr, nlhdlr, expr, SCIPgetConsExprNlhdlrExprData(nlhdlr, expr), sol, 0.0, overestimate, 0.0, rowprep, &success, FALSE, &dummy) );
   cr_expect(!success);
   SCIPfreeRowprep(scip, &rowprep);

   /* auxvar = 0.1 => both inequalities needs to be used */
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr), 0.1) );
   overestimate = FALSE;
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, FALSE) );
   SCIP_CALL( nlhdlrEstimateBilinear(scip, conshdlr, nlhdlr, expr, SCIPgetConsExprNlhdlrExprData(nlhdlr, expr), sol, 0.0, overestimate, 0.0, rowprep, &success, FALSE, &dummy) );
   cr_expect(success);
   SCIP_CALL( SCIPgetRowprepRowConshdlr(scip, &row, rowprep, conshdlr) );
   SCIP_CALL( checkCut(row, 1.45445115010332, 0.9772255750516624, -1.0, -1.9213268615442665) );
   SCIP_CALL( SCIPreleaseRow(scip, &row) );
   SCIPfreeRowprep(scip, &rowprep);

   /* free memory */
   SCIP_CALL( freeEnfoData(expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
}

/* tests the interval evaluation callback when linear inequalities are available */
Test(nlhdlrbilinear, inteval_corner)
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_INTERVAL interval;
   SCIP_Bool success;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 3.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -3.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 5.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x> * <t_y>", NULL, &expr) );
   SCIP_CALL( createAndDetect(&cons, expr) );

   /* INITLP should have added an auxiliary variable to the product expression (tight might change in the future) */
   cr_assert( SCIPgetConsExprExprAuxVar(expr) != NULL);

   /* add linear inequality x <= -y + 7 => maximum of xy is attained at (3,4) */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, 1.0, -1.0, 7.0, &success) );
   cr_expect(success);

   /* reevaluate expression (SCIPincrementConsExprCurBoundsTag would normally be called by prop_obbt) */
   SCIPincrementConsExprCurBoundsTag(conshdlr, FALSE);
   SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, expr, &interval, FALSE) );

   /* compute interval */
   cr_expect(SCIPisEQ(scip, SCIPintervalGetInf(interval), -10.0), "expect -10.0 got %g\n", SCIPintervalGetInf(interval));
   cr_expect(SCIPisEQ(scip, SCIPintervalGetSup(interval), 12.0), "expect 12.0 got %g\n", SCIPintervalGetSup(interval));

   /* free memory */
   SCIP_CALL( freeEnfoData(expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* tests the interval evaluation callback when linear inequalities are available */
Test(nlhdlrbilinear, inteval_single_line)
{
   SCIP_CONS* cons;
   SCIP_INTERVAL interval;
   SCIP_Bool success;
   SCIP_CONSEXPR_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -3.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 3.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 5.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x> * <t_y>", NULL, &expr) );
   SCIP_CALL( createAndDetect(&cons, expr) );

   /* add linear inequality -x <= -y + 3.5 => maximum of xy is attained at (-1.75,1.75) */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, -1.0, -1.0, 3.5, &success) );
   cr_expect(success);

   /* reevaluate expression (SCIPincrementConsExprCurBoundsTag would normally be called by prop_obbt) */
   SCIPincrementConsExprCurBoundsTag(conshdlr, FALSE);
   SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, expr, &interval, FALSE) );

   /* compute interval */
   cr_expect(SCIPisEQ(scip, SCIPintervalGetInf(interval), -3.0625), "expect -3.0625 got %g\n", SCIPintervalGetInf(interval));
   cr_expect(SCIPisEQ(scip, SCIPintervalGetSup(interval), 15.0), "expect 15.0 got %g\n", SCIPintervalGetSup(interval));

   /* free memory */
   SCIP_CALL( freeEnfoData(expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* tests the interval evaluation callback when three intersecting linear inequalities are available */
Test(nlhdlrbilinear, inteval_three_lines)
{
   SCIP_CONSEXPR_EXPR* varexprs[2];
   SCIP_CONS* cons;
   SCIP_INTERVAL interval;
   SCIP_Bool success;
   SCIP_CONSEXPR_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -4.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -2.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -3.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 5.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create expression constraint */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &varexprs[0], x));
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &varexprs[1], y));
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr, 2, varexprs, -2.0));
   SCIP_CALL( createAndDetect(&cons, expr) );

   /* add -x <= -y + 3.5 and x <= -0.6 y -2.25 */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, -1.0, -1.0, 3.5, &success) );
   cr_expect(success);
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, 1.0, -0.6, -2.25, &success) );
   cr_expect(success);
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, -1.0, 0.1, 3.5, &success) );
   cr_expect(success);

   /* reevaluate expression (SCIPincrementConsExprCurBoundsTag would normally be called by prop_obbt) */
   SCIPincrementConsExprCurBoundsTag(conshdlr, FALSE);
   SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, expr, &interval, FALSE) );

   /* compute interval */
   cr_expect(SCIPisEQ(scip, SCIPintervalGetInf(interval), -19.2), "expect -19.2 got %g\n", SCIPintervalGetInf(interval));
   cr_expect(SCIPisEQ(scip, SCIPintervalGetSup(interval), 4.248046875), "expect 4.248046875 got %g\n", SCIPintervalGetSup(interval));

   /* free memory */
   SCIP_CALL( freeEnfoData(expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[1]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[0]) );
}

/* tests the interval evaluation callback when linear inequalities are available */
Test(nlhdlrbilinear, inteval_parallel)
{
   SCIP_CONS* cons;
   SCIP_INTERVAL interval;
   SCIP_Bool success;
   SCIP_CONSEXPR_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x> * <t_y>", NULL, &expr) );
   SCIP_CALL( createAndDetect(&cons, expr) );

   /* add inequalities with the same slope: x <= y + 1 and x >= y - 1 */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, 1.0, 1.0, 1.0, &success) );
   cr_expect(success);
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, -1.0, -1.0, 1.0, &success) );
   cr_expect(success);

   /* reevaluate expression (SCIPincrementConsExprCurBoundsTag would normally be called by prop_obbt) */
   SCIPincrementConsExprCurBoundsTag(conshdlr, FALSE);
   SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, expr, &interval, FALSE) );

   /* compute interval */
   cr_expect(SCIPisEQ(scip, SCIPintervalGetInf(interval), -0.25), "expect -0.25 got %g\n", SCIPintervalGetInf(interval));
   cr_expect(SCIPisEQ(scip, SCIPintervalGetSup(interval), 1.0), "expect 1.0 got %g\n", SCIPintervalGetSup(interval));

   /* free memory */
   SCIP_CALL( freeEnfoData(expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* tests the reverse propagation callback when linear inequalities are available */
Test(nlhdlrbilinear, reverseprop_single)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata;
   SCIP_CONS* cons;
   SCIP_INTERVAL activity;
   SCIP_INTERVAL intervalx;
   SCIP_INTERVAL intervaly;
   SCIP_Bool success;
   SCIP_CONSEXPR_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x> * <t_y>", NULL, &expr) );
   SCIP_CALL( createAndDetect(&cons, expr) );

   /* add inequality x <= y - 1 */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, 1.0, 1.0, -1.0, &success) );
   cr_expect(success);

   /* reevaluate expression (SCIPincrementConsExprCurBoundsTag would normally be called by prop_obbt) */
   SCIPincrementConsExprCurBoundsTag(conshdlr, FALSE);
   SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, expr, &activity, FALSE) );
   cr_expect(SCIPisEQ(scip, activity.inf, -1.0));
   cr_expect(SCIPisEQ(scip, activity.sup, 0.0));

   /* compute intervals; inequality cuts off two corner points */
   nlhdlrexprdata = SCIPgetConsExprNlhdlrExprData(nlhdlr, expr);
   reversePropBilinear(scip, expr, nlhdlrexprdata->underineqs, nlhdlrexprdata->nunderineqs, nlhdlrexprdata->overineqs, nlhdlrexprdata->noverineqs, &intervalx, &intervaly);
   cr_expect(SCIPisEQ(scip, intervalx.inf, -1.0));
   cr_expect(SCIPisEQ(scip, intervalx.sup, 0.0));
   cr_expect(SCIPisEQ(scip, intervaly.inf, 0.0));
   cr_expect(SCIPisEQ(scip, intervaly.sup, 1.0));

   /* free memory */
   SCIP_CALL( freeEnfoData(expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* tests the reverse propagation callback when a linear inequality together with the level set of a bilinear term
 * results in stronger variable bounds
 */
Test(nlhdlrbilinear, reverseprop_levelset)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata;
   SCIP_CONS* cons;
   SCIP_INTERVAL intervalx;
   SCIP_INTERVAL intervaly;
   SCIP_INTERVAL interval;
   SCIP_Bool success;
   SCIP_Bool cutoff = FALSE;
   int ntightenings = 0;
   SCIP_CONSEXPR_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<x> * <y>", NULL, &expr) );
   SCIP_CALL( createAndDetect(&cons, expr) );

   /* add inequality x <= 0.7 y + 0.1 */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, 1.0, 0.7, 0.1, &success) );
   cr_expect(success);

   /* reevaluate expression (SCIPincrementConsExprCurBoundsTag would normally be called by prop_obbt) */
   SCIPincrementConsExprCurBoundsTag(conshdlr, FALSE);
   SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, expr, &interval, FALSE) );
   cr_expect(SCIPisEQ(scip, interval.inf, -1.0));
   cr_expect(SCIPisEQ(scip, interval.sup, 1.0));

   /* tighten the lower bound of the expression */
   SCIPintervalSetBounds(&interval, -0.8, 0.5);
   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, expr, interval, FALSE, NULL, &cutoff, &ntightenings) );
   cr_expect(!cutoff);
   cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprActivity(scip, expr).inf, -0.8));
   cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprActivity(scip, expr).sup, 0.5));

   /* compute intervals; inequality cuts off two corner points */
   nlhdlrexprdata = SCIPgetConsExprNlhdlrExprData(nlhdlr, expr);
   reversePropBilinear(scip, expr, nlhdlrexprdata->underineqs, nlhdlrexprdata->nunderineqs, nlhdlrexprdata->overineqs, nlhdlrexprdata->noverineqs, &intervalx, &intervaly);
   cr_expect(SCIPisEQ(scip, intervalx.inf, -1.0));
   cr_expect_float_eq(intervalx.sup, 0.64371709509590, 1e-8);
   cr_expect(intervalx.sup >= 0.64371709509590); /* intervalarith.c gives us safe bounds */
   cr_expect_float_eq(intervaly.inf, -0.91959585365193, 1e-8);
   cr_expect(intervaly.inf <= -0.91959585365193); /* intervalarith.c gives us safe bounds */
   cr_expect(SCIPisEQ(scip, intervaly.sup, 1.0));

   /* free memory */
   SCIP_CALL( freeEnfoData(expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* tests the reverse propagation callback when a linear inequality does not intersect the level set */
Test(nlhdlrbilinear, reverseprop_levelset_nointersection)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata;
   SCIP_CONS* cons;
   SCIP_INTERVAL intervalx;
   SCIP_INTERVAL intervaly;
   SCIP_INTERVAL interval;
   SCIP_Bool success;
   SCIP_Bool cutoff = FALSE;
   int ntightenings = 0;
   SCIP_CONSEXPR_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<x> * <y>", NULL, &expr) );
   SCIP_CALL( createAndDetect(&cons, expr) );

   /* add inequality -x <= y */
   SCIP_CALL( SCIPaddConsExprNlhdlrBilinearIneq(scip, nlhdlr, expr, -1.0, 1.0, 0.0, &success) );
   cr_expect(success);

   /* reevaluate expression (SCIPincrementConsExprCurBoundsTag would normally be called by prop_obbt) */
   SCIPincrementConsExprCurBoundsTag(conshdlr, FALSE);
   SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, expr, &interval, FALSE) );
   cr_expect(SCIPisEQ(scip, interval.inf, -1.0));
   cr_expect(SCIPisEQ(scip, interval.sup, 1.0));

   /* tighten the lower bound of the expression */
   SCIPintervalSetBounds(&interval, -1.0, 0.5);
   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, expr, interval, FALSE, NULL, &cutoff, &ntightenings) );
   cr_expect(!cutoff);
   cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprActivity(scip, expr).inf, -1.0));
   cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprActivity(scip, expr).sup, 0.5));

   /* compute intervals; inequality does not intersect the level set in the interior of the domain -> no tightening */
   nlhdlrexprdata = SCIPgetConsExprNlhdlrExprData(nlhdlr, expr);
   reversePropBilinear(scip, expr, nlhdlrexprdata->underineqs, nlhdlrexprdata->nunderineqs, nlhdlrexprdata->overineqs, nlhdlrexprdata->noverineqs, &intervalx, &intervaly);
   cr_expect(SCIPisEQ(scip, intervalx.inf, -1.0));
   cr_expect(SCIPisEQ(scip, intervalx.sup, 1.0));
   cr_expect(SCIPisEQ(scip, intervaly.inf, -1.0));
   cr_expect(SCIPisEQ(scip, intervaly.sup, 1.0));

   /* free memory */
   SCIP_CALL( freeEnfoData(expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
