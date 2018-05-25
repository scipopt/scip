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

/**@file   nlobbt.c
 * @brief  unit test for prop_nlobbt methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "nlpi/nlpi_ipopt.h"
#include "nlpi/nlpioracle.h"
#include "nlpi/nlpi.h"

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

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* include NLPI's */
   SCIP_CALL( SCIPcreateNlpSolverIpopt(SCIPblkmem(scip), &nlpi) );
   cr_assert(nlpi != NULL);

   SCIP_CALL( SCIPincludeNlpi(scip, nlpi) );
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameIpopt(), SCIPgetSolverDescIpopt()) );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -5, 5, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -10, 10, 3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   cr_assert_eq(SCIPgetNVars(scip), 2);
   x = SCIPgetVars(scip)[0];
   y = SCIPgetVars(scip)[1];

   SCIP_CALL( SCIPhashmapCreate(&var2idx, SCIPblkmem(scip), 2) );
}

static
void teardown(void)
{
   /* skip the test if IPOPT is not available */
   if( !SCIPisIpoptAvailableIpopt() )
      return;

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
   SCIP_QUADELEM quadelems[2];
   SCIP_Real lincoefs[2];
   SCIP_Real objcoefs[2];
   int objinds[2];
   SCIP_Real nlscore[2];
   SCIP_Real* primal;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_NLPIORACLE* oracle;
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPR* xexpr;
   SCIP_EXPR* expexpr;
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
   vars[0] = x;
   vars[1] = y;

   lincoefs[0] = 3.0;
   lincoefs[1] = -2.0;
   quadelems[0].idx1 = 0;
   quadelems[0].idx2 = 0;
   quadelems[0].coef = 1.0;
   quadelems[1].idx1 = 1;
   quadelems[1].idx2 = 1;
   quadelems[1].coef = 1.0;
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrows[0], "nlrow_0", 0.0, 2, vars, lincoefs, 2, vars, 2, quadelems, NULL, 2.0, 4.0,
         SCIP_EXPRCURV_CONVEX) );

   quadelems[0].idx1 = 0;
   quadelems[0].idx2 = 0;
   quadelems[0].coef = -0.5;
   quadelems[1].idx1 = 1;
   quadelems[1].idx2 = 1;
   quadelems[1].coef = -2.0;
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrows[1], "nlrow_1", 0.0, 0, NULL, NULL, 2, vars, 2, quadelems, NULL, -4.0, -2.0,
         SCIP_EXPRCURV_CONCAVE) );

   quadelems[0].idx1 = 0;
   quadelems[0].idx2 = 1;
   quadelems[0].coef = 1.0;
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrows[2], "nlrow_2", 0.0, 0, NULL, NULL, 2, vars, 1, quadelems, NULL, -10.0, 10.0,
         SCIP_EXPRCURV_UNKNOWN) );

   lincoefs[0] = 3.0;
   lincoefs[1] = -2.0;
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrows[3], "nlrow_3", -1.0, 2, vars, lincoefs, 0, NULL, 0, NULL, NULL, 1.0, 3.0,
         SCIP_EXPRCURV_LINEAR) );

   lincoefs[0] = 1.0;
   lincoefs[1] = -1.0;

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &xexpr, SCIP_EXPR_VARIDX, 0) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expexpr, SCIP_EXPR_EXP, xexpr) );
   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expexpr, 1, 0, NULL) );
   SCIP_CALL( SCIPexprtreeSetVars(exprtree, 1, vars) );
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrows[4], "nlrow_4", 0.0, 2, vars, lincoefs, 0, NULL, 0, NULL, exprtree, 1.0,
         10.0, SCIP_EXPRCURV_CONVEX) );
   SCIP_CALL( SCIPexprtreeFree(&exprtree) );

   /* create convex NLP relaxation */
   SCIP_CALL( SCIPnlpiCreateProblem(nlpi, &nlpiprob, "convex_NLP") );
   SCIP_CALL( SCIPcreateNlpiProb(scip, nlpi, nlrows, 5, nlpiprob, var2idx, nlscore, -1.5, FALSE, TRUE) );
   cr_assert(nlpiprob != NULL);

   oracle = (SCIP_NLPIORACLE*) SCIPgetNlpiOracleIpopt(nlpiprob);
   cr_assert(oracle != NULL);
   SCIP_CALL( SCIPnlpiOraclePrintProblem(oracle, SCIPgetMessagehdlr(scip), NULL) );

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

   /* set tolerances */
   SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_FEASTOL, SCIPfeastol(scip) * 0.01);
   SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_RELOBJTOL, SCIPfeastol(scip) * 0.01);

   /* min x (OPT = -2.60064056606068e-01) */
   objcoefs[0] = 1.0;
   objcoefs[1] = 0.0;
   objinds[0] = 0;
   objinds[1] = 1;
   SCIP_CALL( SCIPnlpiSetObjective(nlpi, nlpiprob, 2, objinds, objcoefs, 0, NULL, NULL, NULL, 0.0) );
   SCIP_CALL( SCIPnlpiOraclePrintProblem(oracle, SCIPgetMessagehdlr(scip), NULL) );
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );
   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );
   cr_assert(SCIPisFeasEQ(scip, primal[0], -2.60064056606068e-01));

   /* max x (OPT = 5.90606671174385e-01) */
   objcoefs[0] = -1.0;
   objcoefs[1] = 0.0;
   objinds[0] = 0;
   objinds[1] = 1;
   SCIP_CALL( SCIPnlpiSetObjective(nlpi, nlpiprob, 2, objinds, objcoefs, 0, NULL, NULL, NULL, 0.0) );
   SCIP_CALL( SCIPnlpiOraclePrintProblem(oracle, SCIPgetMessagehdlr(scip), NULL) );
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );
   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );
   cr_assert(SCIPisFeasEQ(scip, primal[0], 5.90606671174385e-01));

   /* min y (OPT = -1.39009603494603e+00) */
   objcoefs[0] = 0.0;
   objcoefs[1] = 1.0;
   objinds[0] = 0;
   objinds[1] = 1;
   SCIP_CALL( SCIPnlpiSetObjective(nlpi, nlpiprob, 2, objinds, objcoefs, 0, NULL, NULL, NULL, 0.0) );
   SCIP_CALL( SCIPnlpiOraclePrintProblem(oracle, SCIPgetMessagehdlr(scip), NULL) );
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );
   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );
   cr_assert(SCIPisFeasEQ(scip, primal[1], -1.39009603494603e+00));

   /* max y (OPT = -5.90909090909091e-01) */
   objcoefs[0] = 0.0;
   objcoefs[1] = -1.0;
   objinds[0] = 0;
   objinds[1] = 1;
   SCIP_CALL( SCIPnlpiSetObjective(nlpi, nlpiprob, 2, objinds, objcoefs, 0, NULL, NULL, NULL, 0.0) );
   SCIP_CALL( SCIPnlpiOraclePrintProblem(oracle, SCIPgetMessagehdlr(scip), NULL) );
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );
   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );
   cr_assert(SCIPisFeasEQ(scip, primal[1], -5.90909090909091e-01));

   /* free memory */
   SCIP_CALL( SCIPnlpiFreeProblem(nlpi, &nlpiprob) );
   for( i = 4; i >= 0; --i )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrows[i]) );
   }
}
