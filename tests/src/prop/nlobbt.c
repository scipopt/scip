/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlobbt.c
 * @brief  unit test for prop_nlobbt methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/nlpioracle.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_NLPI* nlpi;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_HASHMAP* var2idx;

/* creates scip, problem, and includes NLP */
static
void setup(void)
{
   /* skip the test if IPOPT is not available */
   if( !SCIPisIpoptAvailableIpopt() )
      return;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   nlpi = SCIPfindNlpi(scip, "ipopt");
   cr_assert_not_null(nlpi);

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -5, 5, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -10, 10, 3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* prevent that vars are fixed */
   SCIP_CALL( SCIPaddVarLocks(scip, x, 1, 1) );
   SCIP_CALL( SCIPaddVarLocks(scip, y, 1, 1) );

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   cr_assert_eq(SCIPgetNVars(scip), 2);

   SCIP_CALL( SCIPhashmapCreate(&var2idx, SCIPblkmem(scip), 2) );
}

static
void teardown(void)
{
   /* skip the test if IPOPT is not available */
   if( !SCIPisIpoptAvailableIpopt() )
      return;

   /* removing the locks again that got duplicated in the meanwhile */
   SCIP_CALL( SCIPaddVarLocks(scip, SCIPgetVars(scip)[1], -1, -1) );
   SCIP_CALL( SCIPaddVarLocks(scip, SCIPgetVars(scip)[0], -1, -1) );
   SCIP_CALL( SCIPfreeTransform(scip) );
   SCIP_CALL( SCIPaddVarLocks(scip, y, -1, -1) );
   SCIP_CALL( SCIPaddVarLocks(scip, x, -1, -1) );

   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   SCIPhashmapFree(&var2idx);
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

Test(propagation, convexnlp, .init = setup, .fini = teardown,
   .description = "checks the convex NLP relaxation"
   )
{
   SCIP_NLROW* nlrows[5];
   SCIP_VAR* vars[2];
   SCIP_EXPR* varexprs[2];
   SCIP_Real lincoefs[2];
   SCIP_Real objcoefs[2];
   int objinds[2];
   SCIP_Real nlscore[2];
   SCIP_Real* primal;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_NLPIORACLE* oracle;
   SCIP_EXPR* expexpr;
   SCIP_EXPR* prodexpr;
   SCIP_EXPR* sumexpr;
   SCIP_Real coef;
   int i;

   /* skip the test if IPOPT is not available */
   if( !SCIPisIpoptAvailableIpopt() )
      return;

   /* test the following (nonconvex) optimization problem
    *
    *    min  x + 3y
    *   s.t.  2 <= x^2 + y^2 + 3x -2y <=  4
    *        -4 <= -0.5 x^2 - 2 y^2   <= -2
    *       -10 <= xy                 <= 10
    *         1 <= 3x + -2y - 1       <=  3
    *         1 <= e^x  + x -y        <= 10
    *              x in [-5,5]
    *              y in [-10,10]
    *
    *   with an objective cutoff of -1.5 (which corresponds to the first row)
    */
   vars[0] = SCIPgetVars(scip)[0];
   vars[1] = SCIPgetVars(scip)[1];

   SCIP_CALL( SCIPcreateExprVar(scip, &varexprs[0], vars[0], NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexprs[1], vars[1], NULL, NULL) );

   lincoefs[0] = 3.0;
   lincoefs[1] = -2.0;

   /* x^2 */
   SCIP_CALL( SCIPcreateExprPow(scip, &prodexpr, varexprs[0], 2.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 1, &prodexpr, NULL, 0.0, NULL, NULL) );
   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );

   /* + y^2 */
   SCIP_CALL( SCIPcreateExprPow(scip, &prodexpr, varexprs[1], 2.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, sumexpr, prodexpr, 1.0) );
   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );

   SCIP_CALL( SCIPcreateNlRow(scip, &nlrows[0], "nlrow_0", 0.0, 2, vars, lincoefs, sumexpr, 2.0, 4.0, SCIP_EXPRCURV_CONVEX) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );

   /* -0.5 x^2 */
   SCIP_CALL( SCIPcreateExprPow(scip, &prodexpr, varexprs[0], 2.0, NULL, NULL) );
   coef = -0.5;
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 1, &prodexpr, &coef, 0.0, NULL, NULL) );
   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );

   /* -2.0 y^2 */
   SCIP_CALL( SCIPcreateExprPow(scip, &prodexpr, varexprs[1], 2.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, sumexpr, prodexpr, -2.0) );
   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );

   SCIP_CALL( SCIPcreateNlRow(scip, &nlrows[1], "nlrow_1", 0.0, 0, NULL, NULL, sumexpr, -4.0, -2.0, SCIP_EXPRCURV_CONCAVE) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );

   /* x*y */
   SCIP_CALL( SCIPcreateExprProduct(scip, &prodexpr, 2, varexprs, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrows[2], "nlrow_2", 0.0, 0, NULL, NULL, prodexpr, -10.0, 10.0, SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );

   lincoefs[0] = 3.0;
   lincoefs[1] = -2.0;
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrows[3], "nlrow_3", -1.0, 2, vars, lincoefs, NULL, 1.0, 3.0, SCIP_EXPRCURV_LINEAR) );

   lincoefs[0] = 1.0;
   lincoefs[1] = -1.0;

   SCIP_CALL( SCIPcreateExprExp(scip, &expexpr, varexprs[0], NULL, NULL) );
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrows[4], "nlrow_4", 0.0, 2, vars, lincoefs, expexpr, 1.0, 10.0, SCIP_EXPRCURV_CONVEX) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expexpr) );

   /* create convex NLP relaxation */
   SCIP_CALL( SCIPcreateNlpiProblemFromNlRows(scip, nlpi, &nlpiprob, "convex_NLP", nlrows, 5, var2idx, NULL, nlscore, -1.5, FALSE, TRUE) );
   cr_assert(nlpiprob != NULL);

   oracle = (SCIP_NLPIORACLE*) SCIPgetNlpiOracleIpopt(nlpiprob);
   cr_assert(oracle != NULL);
   SCIP_CALL( SCIPnlpiOraclePrintProblem(scip, oracle, NULL) );

   cr_assert(SCIPnlpiOracleGetNConstraints(oracle) == 5);
   cr_assert(SCIPnlpiOracleGetNVars(oracle) == 2);

   cr_assert(nlscore[0] == 3);
   cr_assert(nlscore[1] == 2);

   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetVarLbs(oracle)[0], SCIPvarGetLbLocal(x)));
   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetVarUbs(oracle)[0], SCIPvarGetUbLocal(x)));
   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetVarLbs(oracle)[1], SCIPvarGetLbLocal(y)));
   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetVarUbs(oracle)[1], SCIPvarGetUbLocal(y)));

   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetConstraintLhs(oracle, 0), -SCIPinfinity(scip)));
   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetConstraintRhs(oracle, 0), -1.5));
   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetConstraintLhs(oracle, 1), -SCIPinfinity(scip)));
   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetConstraintRhs(oracle, 1), 4.0));
   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetConstraintLhs(oracle, 2), -4.0));
   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetConstraintRhs(oracle, 2), SCIPinfinity(scip)));
   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetConstraintLhs(oracle, 3), 2.0));
   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetConstraintRhs(oracle, 3), 4.0));
   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetConstraintLhs(oracle, 4), -SCIPinfinity(scip)));
   cr_assert(SCIPisEQ(scip, SCIPnlpiOracleGetConstraintRhs(oracle, 4), 10.0));

   /* min x (OPT = -2.60064056606068e-01) */
   objcoefs[0] = 1.0;
   objcoefs[1] = 0.0;
   objinds[0] = 0;
   objinds[1] = 1;
   SCIP_CALL( SCIPsetNlpiObjective(scip, nlpi, nlpiprob, 2, objinds, objcoefs, NULL, 0.0) );
   SCIP_CALL( SCIPnlpiOraclePrintProblem(scip, oracle, NULL) );
   SCIP_CALL( SCIPsolveNlpi(scip, nlpi, nlpiprob, .feastol = SCIPfeastol(scip) * 0.01, .opttol = SCIPfeastol(scip) * 0.01) );
   SCIP_CALL( SCIPgetNlpiSolution(scip, nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );
   cr_assert(SCIPisFeasEQ(scip, primal[0], -2.60064056606068e-01));

   /* max x (OPT = 5.90606671174385e-01) */
   objcoefs[0] = -1.0;
   objcoefs[1] = 0.0;
   objinds[0] = 0;
   objinds[1] = 1;
   SCIP_CALL( SCIPsetNlpiObjective(scip, nlpi, nlpiprob, 2, objinds, objcoefs, NULL, 0.0) );
   SCIP_CALL( SCIPnlpiOraclePrintProblem(scip, oracle, NULL) );
   SCIP_CALL( SCIPsolveNlpi(scip, nlpi, nlpiprob, .feastol = SCIPfeastol(scip) * 0.01, .opttol = SCIPfeastol(scip) * 0.01) );
   SCIP_CALL( SCIPgetNlpiSolution(scip, nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );
   cr_assert(SCIPisFeasEQ(scip, primal[0], 5.90606671174385e-01));

   /* min y (OPT = -1.39009603494603e+00) */
   objcoefs[0] = 0.0;
   objcoefs[1] = 1.0;
   objinds[0] = 0;
   objinds[1] = 1;
   SCIP_CALL( SCIPsetNlpiObjective(scip, nlpi, nlpiprob, 2, objinds, objcoefs, NULL, 0.0) );
   SCIP_CALL( SCIPnlpiOraclePrintProblem(scip, oracle, NULL) );
   SCIP_CALL( SCIPsolveNlpi(scip, nlpi, nlpiprob, .feastol = SCIPfeastol(scip) * 0.01, .opttol = SCIPfeastol(scip) * 0.01) );
   SCIP_CALL( SCIPgetNlpiSolution(scip, nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );
   cr_assert(SCIPisFeasEQ(scip, primal[1], -1.39009603494603e+00));

   /* max y (OPT = -5.90909090909091e-01) */
   objcoefs[0] = 0.0;
   objcoefs[1] = -1.0;
   objinds[0] = 0;
   objinds[1] = 1;
   SCIP_CALL( SCIPsetNlpiObjective(scip, nlpi, nlpiprob, 2, objinds, objcoefs, NULL, 0.0) );
   SCIP_CALL( SCIPnlpiOraclePrintProblem(scip, oracle, NULL) );
   SCIP_CALL( SCIPsolveNlpi(scip, nlpi, nlpiprob, .feastol = SCIPfeastol(scip) * 0.01, .opttol = SCIPfeastol(scip) * 0.01) );
   SCIP_CALL( SCIPgetNlpiSolution(scip, nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );
   cr_assert(SCIPisFeasEQ(scip, primal[1], -5.90909090909091e-01));

   /* free memory */
   SCIP_CALL( SCIPfreeNlpiProblem(scip, nlpi, &nlpiprob) );
   for( i = 4; i >= 0; --i )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrows[i]) );
   }
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[0]) );
}
