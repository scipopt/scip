/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   nl.c
 * @brief  Unittest for nl reader
 * @author Stefan Vigerske
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scipdefplugins.h"
#include "scip/reader_nl.h"
#include "scip/expr_trig.h"
#include "scip/expr_sum.h"
#include "scip/expr_var.h"

#include "include/scip_test.h"

static SCIP* scip;

static
void setup(void)
{
   /* create SCIP instance */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
}

static
void teardown(void)
{
   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );
}

/* TEST SUITE */

/* read .nl file, print as .cip, compare with .cip file on stock */
static
void compareNlToCip(
   const char*           filestub            /**< stub of nl file to read */
   )
{
   char filename[SCIP_MAXSTRLEN];
   FILE* reffile;

   /* skip test if nl reader not available (SCIP compiled with AMPL=false) */
   if( SCIPfindReader(scip, "nlreader") == NULL )
      return;

   /* get file to read: <filestub>.nl that lives in the same directory as this file */
   TESTsetTestfilename(filename, __FILE__, filestub);
   strcat(filename, ".nl");

   /* read nl file */
   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );

   /* write cip file to stdout, but capture stdout output */
   TEST_CAPTURE_STDOUT();
   SCIP_CALL( SCIPwriteOrigProblem(scip, NULL, "cip", FALSE) );
   fflush(stdout);

   /* open reference file with cip */
   TESTsetTestfilename(filename, __FILE__, filestub);
   strcat(filename, ".cip");
   reffile = fopen(filename, "r");
   TEST_ASSERT_NOT_NULL(reffile);

   /* check that problem is as expected */
   TEST_ASSERT_STDOUT_EQUAL_FILE(reffile, "Problem from reading %s.nl not as expected (%s.cip)", filename, filename);
   fclose(reffile);
}

/* the following two tests need to be two separate tests, because the redirect_stdout in compareNlToCip()
 * makes stdout unusable for a second run
 */
/* check that we can read a .nl file with common exprs and get the problem as expected */
/** @brief check reading .nl file with common expression */
void test_readernl_read1(void)
{
   compareNlToCip("commonexpr1");
}

/* check that we can read a .nl file with common exprs and get the problem as expected */
/** @brief check reading .nl file with common expression */
void test_readernl_read2(void)
{
   compareNlToCip("commonexpr2");
}

/* read a .nl file with suffixes and check that they arrive as expected; also solve and check optimal value */
/** @brief check reading .nl file with suffixes */
void test_readernl_read3(void)
{
   char filename[SCIP_MAXSTRLEN];
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;
   SCIP_CONS* e1;
   SCIP_CONS* e2;
   SCIP_CONS* sos;

   /* skip test if nl reader not available (SCIP compiled with AMPL=false) */
   if( SCIPfindReader(scip, "nlreader") == NULL )
      return;

   /* get file to read: suffix1.nl that lives in the same directory as this file */
   TESTsetTestfilename(filename, __FILE__, "suffix1.nl");

   /* read nl file */
   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );

   TEST_ASSERT_EQUAL(SCIPgetNVars(scip), 3);
   TEST_ASSERT_EQUAL(SCIPgetNConss(scip), 3);

   x = SCIPgetVars(scip)[0];
   y = SCIPgetVars(scip)[1];
   z = SCIPgetVars(scip)[2];

   e1 = SCIPgetConss(scip)[0];
   e2 = SCIPgetConss(scip)[1];
   sos = SCIPgetConss(scip)[2];

   /* read .nl file without names from accompanying col/row files */
   SOFT_ASSERT_EQUAL_STRING(SCIPvarGetName(x), "x0");
   SOFT_ASSERT_EQUAL_STRING(SCIPvarGetName(y), "x1");
   SOFT_ASSERT_EQUAL_STRING(SCIPvarGetName(z), "x2");
   SOFT_ASSERT_EQUAL_STRING(SCIPconsGetName(e1), "lc0");
   SOFT_ASSERT_EQUAL_STRING(SCIPconsGetName(e2), "lc1");
   SOFT_ASSERT_EQUAL_STRING(SCIPconsGetName(sos), "sos1_1");

   SOFT_ASSERT_NOT(SCIPvarIsInitial(x));
   SOFT_ASSERT(SCIPvarIsRemovable(x));

   SOFT_ASSERT(SCIPvarIsInitial(y));
   SOFT_ASSERT_NOT(SCIPvarIsRemovable(y));

   SOFT_ASSERT(SCIPconsIsInitial(e1));
   SOFT_ASSERT(SCIPconsIsSeparated(e1));
   SOFT_ASSERT(SCIPconsIsEnforced(e1));
   SOFT_ASSERT(SCIPconsIsChecked(e1));
   SOFT_ASSERT(SCIPconsIsPropagated(e1));
   SOFT_ASSERT_NOT(SCIPconsIsDynamic(e1));
   SOFT_ASSERT_NOT(SCIPconsIsRemovable(e1));

   SOFT_ASSERT_NOT(SCIPconsIsInitial(e2));
   SOFT_ASSERT_NOT(SCIPconsIsSeparated(e2));
   SOFT_ASSERT_NOT(SCIPconsIsEnforced(e2));
   SOFT_ASSERT_NOT(SCIPconsIsChecked(e2));
   SOFT_ASSERT_NOT(SCIPconsIsPropagated(e2));
   SOFT_ASSERT(SCIPconsIsDynamic(e2));
   SOFT_ASSERT(SCIPconsIsRemovable(e2));

   SCIP_CALL( SCIPsolve(scip) );
   SOFT_ASSERT_EQUAL(SCIPgetStatus(scip), SCIP_STATUS_OPTIMAL);
   SOFT_ASSERT_DOUBLE_WITHIN(SCIPgetPrimalbound(scip), 110.0, SCIPfeastol(scip));
}

/* check whether running shell with -AMPL flag works */
/** @brief check running SCIP with -AMPL */
void test_readernl_run(void)
{
   const char* args[3];
   char solfile[SCIP_MAXSTRLEN];
   FILE* sol;
   char scipstatus[100];

   /* skip test if nl reader not available (SCIP compiled with AMPL=false) */
   if( SCIPfindReader(scip, "nlreader") == NULL )
      return;

   args[0] = "dummy";

   /* get file to read: suffix1.nl that lives in the same directory as this file */
   args[1] = (const char*)malloc(SCIP_MAXSTRLEN);
   TESTsetTestfilename((char*)args[1], __FILE__, "suffix1");

   args[2] = "-AMPL";

   /* get name of file where sol will be written: .nl-file with .nl replaced by .sol */
   TESTsetTestfilename(solfile, __FILE__, "suffix1.sol");

   /* make sure no solfile is there at the moment */
   remove(solfile);

   /* specify some SCIP options via scip_options environment variable (like AMPL scripts would do) */
   setenv("scip_options", "limits/time 0", 1);

   /* run SCIP as if called by AMPL
    * this should writes a sol file (into tests/src/reader, unfortunately)
    */
   SCIP_CALL( SCIPrunShell(3, (char**)args, NULL) );

   /* drop the option again: Criterion ran every test in its own process, but all
    * tests of a Unity binary share one, so a leftover scip_options would also
    * apply to every later test (a time limit of 0 makes them give up at once)
    */
   unsetenv("scip_options");

   /* check that a solfile that can be opened exists now */
   sol = fopen(solfile, "r");
   SOFT_ASSERT_NOT_NULL(sol);

   /* check that first line of solfile is the SCIP status, which should be that the timelimit has been reached */
   SOFT_ASSERT_NOT_NULL(fgets(scipstatus, sizeof(scipstatus), sol));
   SOFT_ASSERT_EQUAL_STRING(scipstatus, "time limit reached\n");

   /* cleanup */
   fclose(sol);
   remove(solfile);

   free((char*)args[1]);
}

/* check whether solving a LP without presolve gives a dual solution in the AMPL solution file */
/** @brief check whether solving a LP without presolve gives a dual solution */
void test_readernl_dualsol(void)
{
   char* args[3];
   char solfilename[SCIP_MAXSTRLEN];
   char refsolfilename[SCIP_MAXSTRLEN];
   FILE* solfile;
   FILE* refsolfile;
   FILE* setfile;

   /* skip test if nl reader not available (SCIP compiled with AMPL=false) */
   if( SCIPfindReader(scip, "nlreader") == NULL )
      return;

   args[0] = (char*)"dummy";

   /* get file to read: lp1.nl that lives in the same directory as this file */
   args[1] = (char*)malloc(SCIP_MAXSTRLEN);
   TESTsetTestfilename(args[1], __FILE__, "lp1");

   args[2] = (char*)"-AMPL";

   /* get name of file where sol will be written: .nl-file with .nl replaced by .sol */
   TESTsetTestfilename(solfilename, __FILE__, "lp1.sol");

   /* make sure no solfile is there at the moment */
   remove(solfilename);

   setfile = fopen("nopresolve.set", "w");
   TEST_ASSERT_NOT_NULL(setfile);
   fprintf(setfile, "emphasis:presolving:off\n");
   fprintf(setfile, "emphasis:heuristics:off\n");
   fclose(setfile);

   /* run SCIP as if called by AMPL */
   SCIP_CALL( SCIPrunShell(3, (char**)args, "nopresolve.set") );

   /* check that a solfile that can be opened exists now */
   solfile = fopen(solfilename, "r");
   TEST_ASSERT_NOT_NULL(solfile);

   /* dual solution is not unique; the one we compare with seems to be the one given by CPLEX and SoPlex at the moment (2021) */
   if( strncmp(SCIPlpiGetSolverName(), "CPLEX", 5) == 0 || strncmp(SCIPlpiGetSolverName(), "SoPlex", 6) == 0 )
   {
      /* get name of reference solution file to compare solfile with */
      TESTsetTestfilename(refsolfilename, __FILE__, "lp1.refsol");

      /* open reference solfile */
      refsolfile = fopen(refsolfilename, "r");
      TEST_ASSERT_NOT_NULL(refsolfile);
      SOFT_ASSERT_FILES_EQUAL(solfile, refsolfile);
      fclose(refsolfile);
   }

   /* cleanup */
   fclose(solfile);
   remove(solfilename);
   remove("nopresolve.set");

   free(args[1]);
}

/* check whether initial value of introduced nlobjvar is set */
/** @brief check whether initial value of introduced nlobjvar is set */
void test_readernl_nlobjvar(void)
{
   char filename[SCIP_MAXSTRLEN];

   /* skip test if nl reader not available (SCIP compiled with AMPL=false) */
   if( SCIPfindReader(scip, "nlreader") == NULL )
      return;

   /* get file to read: suffix1.nl that lives in the same directory as this file */
   TESTsetTestfilename(filename, __FILE__, "nlobj.nl");

   /* read nl file */
   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );

   /* we should have an initial solution as well */
   SOFT_ASSERT(SCIPgetNSols(scip) == 1);

   /* only if the solution is feasible, it will be available in the transformed problem, too */
   SCIP_CALL( SCIPtransformProb(scip) );
   SOFT_ASSERT(SCIPgetNSols(scip) == 1);
}

/* check writing of nl file with negative variable in expression */
/** @brief check writing of nl file with negated variable */
void test_readernl_writenegvar(void)
{
   SCIP_CONS* cons;
   SCIP_VAR* x;
   SCIP_VAR* xneg;
   SCIP_EXPR* expr1;
   SCIP_EXPR* expr2;
   SCIP_Bool changed;
   SCIP_Bool infeasible;

   /* skip test if nl reader not available (SCIP compiled with AMPL=false) */
   if( SCIPfindReader(scip, "nlreader") == NULL )
      return;

   SCIP_CALL( SCIPcreateProbBasic(scip, "testnegvar") );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPgetNegatedVar(scip, x, &xneg) );
   SCIP_CALL( SCIPaddVar(scip, xneg) );

   /* we add a constraint sin(~x) == 1 */
   SCIP_CALL( SCIPcreateExprVar(scip, &expr1, xneg, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSin(scip, &expr2, expr1, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "test", expr2, 1.0, 1.0) );

   SCIP_CALL( SCIPreleaseExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr1) );

   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) ); */

   SCIP_CALL( SCIPwriteOrigProblem(scip, "readernl_writenegvartest.nl", "nl", FALSE) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   SCIP_CALL( SCIPreadProb(scip, "readernl_writenegvartest.nl", NULL) );
   /* SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) ); */

   /* we expect to see a constraint sin(1-x) == 1 */
   SOFT_ASSERT(SCIPgetNVars(scip) == 1);
   SOFT_ASSERT(SCIPgetNConss(scip) == 1);
   cons = SCIPgetConss(scip)[0];
   SOFT_ASSERT(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "nonlinear") == 0);

   /* bring into standard form sum(1, -1*x) for checks */
   SCIP_CALL( SCIPsimplifyExpr(scip, SCIPgetExprNonlinear(cons), &expr1, &changed, &infeasible, NULL, NULL) );

   SOFT_ASSERT(SCIPisExprSin(scip, expr1));
   expr2 = SCIPexprGetChildren(expr1)[0];
   SOFT_ASSERT(SCIPisExprSum(scip, expr2));
   SOFT_ASSERT(SCIPexprGetNChildren(expr2) == 1);
   SOFT_ASSERT(SCIPgetConstantExprSum(expr2) == 1.0);
   SOFT_ASSERT(SCIPgetCoefsExprSum(expr2)[0] == -1.0);
   SOFT_ASSERT(SCIPisExprVar(scip, SCIPexprGetChildren(expr2)[0]));
   SOFT_ASSERT(SCIPgetVarExprVar(SCIPexprGetChildren(expr2)[0]) == SCIPgetVars(scip)[0]);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr1) );

   remove("readernl_writenegvartest.nl");
   remove("readernl_writenegvartest.col");
   remove("readernl_writenegvartest.row");
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_readernl_read1);
   RUN_TEST(test_readernl_read2);
   RUN_TEST(test_readernl_read3);
   RUN_TEST(test_readernl_run);
   RUN_TEST(test_readernl_dualsol);
   RUN_TEST(test_readernl_nlobjvar);
   RUN_TEST(test_readernl_writenegvar);
   return UNITY_END();
}
