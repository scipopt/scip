/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   estimation.h
 * @brief  setup and teardown of estimation methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define _USE_MATH_DEFINES   /* to get M_PI on Windows */

#define EXPECTFEQ(a,b) cr_expect_float_eq(a, b, 1e-6, "%s = %g != %g (dif %g)", #a, a, b, ABS(a-b))

#include <strings.h>
#include <math.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/expr_var.h"

/* needed to manipulate the stages */
#include "scip/struct_scip.h"
#include "scip/struct_set.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_VAR* w;
static SCIP_VAR* v;
static SCIP_VAR* auxvar;
static SCIP_SOL* sol;
static SCIP_EXPR* xexpr;
static SCIP_EXPR* yexpr;
static SCIP_EXPR* zexpr;
static SCIP_EXPR* wexpr;
static SCIP_EXPR* vexpr;


static SCIP_RANDNUMGEN* randnumgen; /* needs it for the multilinear separation */

/* creates scip, problem, includes nonlinear constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -6.0, -3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 2.0, 4.0, -3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &v, "v", -M_PI, M_PI, -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &auxvar, "auxvar", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIP_CALL( SCIPaddVar(scip, v) );
   SCIP_CALL( SCIPaddVar(scip, auxvar) );

   SCIP_CALL( SCIPcreateExprVar(scip, &xexpr, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &yexpr, y, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &zexpr, z, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &wexpr, w, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &vexpr, v, NULL, NULL) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   /* free random number generator */
   SCIPfreeRandom(scip, &randnumgen);

   /* release solution */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   SCIP_CALL( SCIPreleaseExpr(scip, &vexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &wexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &zexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &v) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}
