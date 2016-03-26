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

/**@file   cmain.c
 * @brief  unit test for checking cons_expr
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_sumprod.h"

/** macro to check the return of tests
 *
 *  @note assumes the existence of SCIP_RETCODE retcode
 */
#define CHECK_TEST(x)                            \
   do                                            \
   {                                             \
      retcode = (x);                             \
      if( retcode != SCIP_OKAY )                 \
      {                                          \
         printf("Unit test " #x " failed\n");    \
         SCIPprintError(retcode);                \
         return -1;                              \
      }                                          \
   }                                             \
   while( FALSE )


/*
 * Local methods for tests
 */


/*
 * TESTS
 */

/** test something */
static
SCIP_RETCODE testCons(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_CONSEXPR_EXPR* expr_x;
   SCIP_CONSEXPR_EXPR* expr_y;
   SCIP_CONSEXPR_EXPR* expr_5;
   SCIP_CONSEXPR_EXPR* expr_xy5;

   conshdlr = NULL;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* currently expr constraints cannot be created */
   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create expressions for variables x and y */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );

   /* create expression for constant 5 */
   SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr_5, 5.0) );

   /* create expression for product of 5, x, and y (TODO should have something to add children to an existing product expr) */
   {
      SCIP_CONSEXPR_EXPR* xy5[3] = {expr_x, expr_y, expr_5};
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_xy5, 3, xy5, NULL, 2.0) );
   }

   /* release leaf expressions (this should not free them yet, as they are captured by prod_xy5) */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_5) );

   /* release product expression (this should free the product and its children) */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_xy5) );

   /* give up scip */
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   /* test something? */

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}



/** main function */
int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   /* output stuff for automatic unittest evaluation */
   printf("@01 unittest-consexpr ===========\n");
   printf("=opt=  unittest-consexpr 0\n\n");

   CHECK_TEST( testCons() );

   /* for automatic testing output the following */
   printf("SCIP Status        : all tests passed\n");
   printf("Ignore the following:\n");
   printf("  solving          : 0.00\n");
   printf("  nodes (total)    : 0\n");
   printf("  Primal Bound     : 0.0\n");
   printf("  Dual Bound       : 0.0\n");

   return 0;
}
