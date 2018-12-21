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
static const char* testdecname = "../check/instances/Tests/decomp/decomptest.dec";
#define NVARS 5
#define NCONSS 3
/** GLOBAL VARIABLES **/
static SCIP* scip;
static SCIP_DECOMP* decomp;
static SCIP_VAR* vars[NVARS];
static SCIP_CONS* conss[NCONSS];
static int labels_vars[] = {SCIP_DECOMP_LINKVAR,0,0,1,1};
static int labels_conss[] = {SCIP_DECOMP_LINKCONS, 0, 1};
static char strbuf1[1024];
static char strbuf2[1024];

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
   SCIP_CALL( SCIPdecompCreate(&decomp, SCIPblkmem(scip), TRUE) );

   setupData();
}

static
void teardown(void)
{
   SCIPdecompFree(&decomp, SCIPblkmem(scip));
   SCIPfree(&scip);

   BMScheckEmptyMemory();
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
   SCIP_CALL( SCIPdecompCreate(&newdecomp, SCIPblkmem(scip), TRUE) );

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


/** check constraint labels of this decomposition */
static
void checkConsLabels(
   SCIP_DECOMP*          decomposition       /**< some decomposition */
   )
{
   int returnedlabels[NCONSS];
   cr_assert_not_null(decomposition);

   SCIPdecompGetConsLabels(decomposition, conss, returnedlabels, NCONSS);

   cr_assert_arr_eq(returnedlabels, labels_conss, NCONSS,
      "Array {%s} not equal to {%s}\n",
      printIntArray(strbuf1, returnedlabels, NCONSS),
      printIntArray(strbuf2, labels_conss, NCONSS));

}

Test(decomptest, test_cons_labeling, .description="check constraint label computation")
{
   SCIP_CALL( SCIPdecompSetVarsLabels(decomp, vars, labels_vars, NVARS) );

   SCIP_CALL( SCIPdecompComputeConsLabels(scip, decomp, conss, NCONSS) );

   checkConsLabels(decomp);

}
/** check variable labels of this decomposition */
static
void checkVarsLabels(
   SCIP_DECOMP*          decomposition       /**< some decomposition */
   )
{
   int returnedlabels[NVARS];
   cr_assert_not_null(decomposition);
   SCIPdecompGetVarsLabels(decomposition, vars, returnedlabels, NVARS);

   cr_assert_arr_eq(returnedlabels, labels_vars, NVARS,
      "Array {%s} not equal to {%s}\n",
      printIntArray(strbuf1, returnedlabels, NVARS),
      printIntArray(strbuf2, labels_vars, NVARS)
      );
}

Test(decomptest, test_var_labeling, .description="check variable label computation")
{

   SCIP_CALL( SCIPdecompSetConsLabels(decomp, conss, labels_conss, NCONSS) );

   SCIP_CALL( SCIPdecompComputeVarsLabels(scip, decomp, conss, NCONSS) );

   checkVarsLabels(decomp);

}

Test(decomptest, test_dec_reader, .description="test decomposition reader")
{
   SCIP_DECOMPSTORE* decompstore = SCIPgetDecompstore(scip);
   SCIP_DECOMP* scip_decomp;
   SCIP_VAR* transvars[NVARS];
   int returnedlabels[NVARS];
   int v;

   assert(decompstore != NULL);

   SCIP_CALL( SCIPreadProb(scip, testdecname, "dec") );

   cr_assert_eq(SCIPdecompstoreGetNOrigDecomps(decompstore), 1);

   scip_decomp = SCIPdecompstoreGetOrigDecomps(decompstore)[0];
   cr_assert_not_null(scip_decomp);

   checkConsLabels(scip_decomp);

   checkVarsLabels(scip_decomp);

   /* solve the problem without presolving */
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsolve(scip) );

   /* transform variable array */
   for( v = 0; v < NVARS; ++v )
   {
      transvars[v] = SCIPvarGetTransVar(vars[v]);
      cr_assert_not_null(transvars[v]);
   }

   cr_assert_eq(SCIPdecompstoreGetNDecomps(decompstore), 1, "Number of transformed decompositions should be 1.\n");
   scip_decomp = SCIPdecompstoreGetDecomps(decompstore)[0];

   SCIPdecompGetVarsLabels(scip_decomp, transvars, returnedlabels, NVARS);

   cr_assert_arr_eq(returnedlabels, labels_vars, NVARS,
      "Arrays should be equal: {%s} != {%s}",
      printIntArray(strbuf1, returnedlabels, NVARS),
      printIntArray(strbuf2, labels_vars, NVARS)
      );
}
