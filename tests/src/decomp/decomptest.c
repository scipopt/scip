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

/**@file   decomptest.c
 * @brief  unit test
 * @author Gregor Hendel
 */

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/decomp.h"

#include "include/scip_test.h"

static const char* testfilename = "../check/instances/Tests/decomp/decomptest.cip";
#define NVARS 5
#define NCONSS 3
/** GLOBAL VARIABLES **/
static SCIP* scip;
static SCIP_DECOMP* decomp;
static SCIP_VAR* vars[NVARS];
static SCIP_CONS* conss[NCONSS];
static int labels_vars[] = {0,0,0,1,1};
static int labels_conss[] = {SCIP_DECOMP_LINKCONS, 0, 1};

/** set up some data structures */
static
void setupData(void)
{
   vars[0] = SCIPfindVar(scip, "y");
   vars[1] = SCIPfindVar(scip, "x1");
   vars[2] = SCIPfindVar(scip, "x2");
   vars[3] = SCIPfindVar(scip, "x3");
   vars[4] = SCIPfindVar(scip, "x4");

   conss[0] = SCIPfindCons(scip, "linkingcons");
   conss[1] = SCIPfindCons(scip, "block1cons");
   conss[2] = SCIPfindCons(scip, "block2cons");

}

/** test setup */
static
void testData(void)
{
   int i;

   /* check that all variables and constraints are there */
   for( i = 0; i < NVARS; ++i )
      cr_assert(vars[i] != NULL);

   for( i = 0; i < NCONSS; ++i )
      cr_assert(conss[i] != NULL);
}

/* TEST SUITE */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );


   SCIP_CALL( SCIPreadProb(scip, testfilename, NULL) );
   SCIP_CALL( SCIPdecompCreate(&decomp, SCIPblkmem(scip)) );

   setupData();
}

static
void teardown(void)
{
   SCIPdecompFree(&decomp, SCIPblkmem(scip));
   SCIPfree(&scip);
}

TestSuite(decomptest, .init = setup, .fini = teardown);

/* TESTS  */
Test(decomptest, create_and_free)
{
   /* calls setup and teardown */
}

Test(decomptest, create_decomp, .description="test constructor and destructor of decomposition")
{
   SCIP_DECOMP* newdecomp;
   SCIP_CALL( SCIPdecompCreate(&newdecomp, SCIPblkmem(scip)) );

   SCIPdecompFree(&newdecomp, SCIPblkmem(scip));
}

Test(decomptest, test_data_setup, .description = "check data setup")
{
   testData();
}

Test(decomptest, test_setters_and_getters, .description="check that setting labels works")
{
   int returnedlabels[NVARS];

   SCIP_CALL( SCIPdecompSetVarsLabels(decomp, vars, labels_vars, NVARS) );

   SCIPdecompGetVarsLabels(decomp, vars, returnedlabels, NVARS);

   /* check that each variable has the correct label */
   cr_assert_arr_eq(returnedlabels, labels_vars, NVARS);
}

/** print integer array */
static
char* printIntArray(
   char*                 strbuf,
   int*                  array,
   int                   length
   )
{
   int i;
   char* strptr = strbuf;

   /* print entries */
   for( i = 0; i < length; ++i )
   {
      strptr += sprintf(strptr, "%d ", array[i]);
   }

   return strbuf;
}

Test(decomptest, test_cons_labeling, .description="check constraint label computation")
{
   int returnedlabels[NCONSS];
   char strbuf[1024];

   SCIP_CALL( SCIPdecompSetVarsLabels(decomp, vars, labels_vars, NVARS) );

   SCIP_CALL( SCIPdecompComputeConsLabels(scip, decomp, conss, NCONSS) );

   SCIPdecompGetConsLabels(decomp, conss, returnedlabels, NCONSS);

   cr_assert_arr_eq(returnedlabels, labels_conss, NCONSS, "%s\n", printIntArray(strbuf, labels_conss, NCONSS));
}
