/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_erf.c
 * @brief  handler for guassian error function expressions
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_erf.h"
#include "scip/cons_expr.h"

/* fundamental expression handler properties */
#define EXPRHDLR_NAME         "erf"
#define EXPRHDLR_DESC         "expression handler template"
#define EXPRHDLR_PRECEDENCE   0
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(1.0)

/*
 * Data structures
 */

/*
 * Local methods
 */

/* TODO: put your local methods here, and declare them static */

/*
 * Callback methods of expression handler
 */

/** expression handler copy callback */
static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrErf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
}

/** expression handler free callback */
static
SCIP_DECL_CONSEXPR_EXPRFREEHDLR(freehdlrErf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
}

/** simplifies a erf expression */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
}

/** expression compare callback */
static
SCIP_DECL_CONSEXPR_EXPRCOMPARE(compareErf)
{  /*lint --e{715}*/
   assert(expr1 != NULL);
   assert(expr2 != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return 0;
}

/** expression data copy callback */
static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataErf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression print callback */
static
SCIP_DECL_CONSEXPR_EXPRPRINT(printErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_CONSEXPR_EXPRPARSE(parseErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPREVAL(evalErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** separation initialization callback */
static
SCIP_DECL_CONSEXPR_EXPRINITSEPA(initSepaErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** separation deinitialization callback */
static
SCIP_DECL_CONSEXPR_EXPREXITSEPA(exitSepaErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression separation callback */
static
SCIP_DECL_CONSEXPR_EXPRSEPA(sepaErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression under/overestimation callback */
static
SCIP_DECL_CONSEXPR_EXPRESTIMATE(estimateErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_EXPRREVERSEPROP(reversepropErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** erf hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashErf)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(EXPRHDLR_HASHKEY != 0.0);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_CONSEXPR_EXPRCURVATURE(curvatureErf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_CONSEXPR_EXPRMONOTONICITY(monotonicityErf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEGRALITY(integralityErf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression branching score callback */
static
SCIP_DECL_CONSEXPR_EXPRBRANCHSCORE(branchscoreErf)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(success != NULL);

   *success = FALSE;

   SCIPerrorMessage("method of erf constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** creates the handler for erf expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrErf(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLRDATA* exprhdlrdata;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   /* create expression handler data */
   exprhdlrdata = NULL;

   /* TODO: create and store expression handler specific data here */

   /* include expression handler */
   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC,
         EXPRHDLR_PRECEDENCE, evalErf, exprhdlrdata) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrErf, freehdlrErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataErf, freedataErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifyErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCompare(scip, consexprhdlr, exprhdlr, compareErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrParse(scip, consexprhdlr, exprhdlr, parseErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, initSepaErf, exitSepaErf, sepaErf, estimateErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCurvature(scip, consexprhdlr, exprhdlr, curvatureErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrMonotonicity(scip, consexprhdlr, exprhdlr, monotonicityErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntegrality(scip, consexprhdlr, exprhdlr, integralityErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBranchscore(scip, consexprhdlr, exprhdlr, branchscoreErf) );

   return SCIP_OKAY;
}

/** creates a erf expression */
SCIP_RETCODE SCIPcreateConsExprExprErf(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_CONSEXPR_EXPR**  children            /**< children (can be NULL if nchildren is 0) */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(consexprhdlr != NULL);
   assert(expr != NULL);

   exprhdlr = SCIPfindConsExprExprHdlr(consexprhdlr, EXPRHDLR_NAME);

   if( exprhdlr != NULL )
   {
      SCIPerrorMessage("could not find %s expression handler -> abort\n", EXPRHDLR_NAME);
      SCIPABORT();
      return SCIP_ERROR;
   }

   /* create expression data */
   exprdata = NULL;

   /* TODO: create and store expression specific data here */

   /* create expression */
   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, exprhdlr, exprdata, nchildren, children) );

   return SCIP_OKAY;
}
