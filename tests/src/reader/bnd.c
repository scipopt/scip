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

/**@file   bnd.c
 * @brief  Unittest for bound reader
 * @author Benjamin Mueller
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"
#include "scip/reader_bnd.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_Real lb = -1.0;
static SCIP_Real ub = 2.0;
static const char* varname = "x";

static
void setup(void)
{
   /* create SCIP instance */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBnd(scip) );

   /* create single variable */
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, varname, lb, ub, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
}

static
void teardown(void)
{
   /* free variable and SCIP */
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );
}

/* TEST SUITE */
TestSuite(readerbnd, .init = setup, .fini = teardown);

Test(readerbnd, read, .description = "check the function for reading a *.bnd file")
{
   FILE* fp;
   const char* filename = "input.bnd";

   /* write *.bnd file */
   fp = fopen(filename, "w");
   fprintf(fp, "%s", "<x> -3.0 0.0\n");
   fclose(fp);

   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );
   cr_expect(SCIPisEQ(scip, SCIPvarGetLbGlobal(x), -3.0));
   cr_expect(SCIPisEQ(scip, SCIPvarGetUbGlobal(x), 0.0));

   (void)remove(filename);
}

Test(readerbnd, write, .description = "check the function for writting a *.bnd file")
{
   char formatstr[SCIP_MAXSTRLEN];

   /* check that the written bounds of x equal [-1,2] */
   cr_redirect_stdout();
   SCIP_CALL( SCIPwriteOrigProblem(scip, NULL, "bnd", FALSE) );
   fflush(stdout);

   sprintf(formatstr, "<%s> %16.15f %16.15f\n", varname, lb, ub);

   cr_assert_stdout_eq_str(formatstr);
}
