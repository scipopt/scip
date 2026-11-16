#include "scip/expr_pow.c"
#include "../estimation.h"

void test_estimation_convexsquare(void)
{
   SCIP_EXPR* expr;
   SCIP_Real xval;
   SCIP_Real constant;
   SCIP_Real slope;
   SCIP_Bool islocal;
   SCIP_Bool branchcand;
   SCIP_Bool success;
   SCIP_INTERVAL bnd;

   SCIP_CALL( SCIPcreateExprPow(scip, &expr, xexpr, 2.0, NULL, NULL) );
   SCIP_CALL( SCIPevalExprActivity(scip, expr) );
   bnd = SCIPexprGetActivity(xexpr);

   xval = 1.0;
   branchcand = TRUE;
   SCIP_CALL( estimatePow(scip, expr, &bnd, &bnd, &xval, FALSE, SCIPinfinity(scip), &slope, &constant, &islocal, &success, &branchcand) );

   printf("success=%d constant=%g slope=%g islocal=%d branchcand=%d\n", success, constant, slope, islocal, branchcand);
   printf("Expected: constant=-1, slope=2\n");

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

void setUp(void) { setup(); }
void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_estimation_convexsquare);
   return UNITY_END();
}
