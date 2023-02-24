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

/**@file   decomptest.c
 * @brief  unit test
 * @author Gregor Hendel
 */

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static const char* testfilename = "../check/instances/Tests/decomp/decomptest.cip";
static const char* testdecname = "../check/instances/Tests/decomp/decomptest.dec";
#define NVARS 7
#define NCONSS 3
/** GLOBAL VARIABLES **/
static SCIP* scip;
static SCIP_DECOMP* decomp;
static SCIP_VAR* vars[NVARS];
static SCIP_CONS* conss[NCONSS];
static int labels_vars[] = {SCIP_DECOMP_LINKVAR,0,0,1,1,0,1};
static int benderslabels_vars[] = {SCIP_DECOMP_LINKVAR,SCIP_DECOMP_LINKVAR,SCIP_DECOMP_LINKVAR,SCIP_DECOMP_LINKVAR,
   SCIP_DECOMP_LINKVAR,0,1};
static int labels_conss[] = {SCIP_DECOMP_LINKCONS, 0, 1};
static int nblocks = 2; /* only blocks that aren't linking blocks are counted */
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
   vars[5] = SCIPfindVar(scip, "z1");
   vars[6] = SCIPfindVar(scip, "z2");

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
   SCIP_CALL( SCIPcreateDecomp(scip, &decomp, nblocks, TRUE, FALSE) );

   setupData();
}

static
void teardown(void)
{
   SCIPfreeDecomp(scip, &decomp);
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
   SCIP_CALL( SCIPcreateDecomp(scip, &newdecomp, 1, TRUE, FALSE) );

   SCIPfreeDecomp(scip, &newdecomp);
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

   SCIP_CALL( SCIPcomputeDecompConsLabels(scip, decomp, conss, NCONSS) );

   checkConsLabels(decomp);

}
/** check variable labels of this decomposition */
static
void checkVarsLabels(
   SCIP_DECOMP*          decomposition,      /**< some decomposition */
   int                   varlabels[]         /**< the variable labels to check against */
   )
{
   int returnedlabels[NVARS];
   cr_assert_not_null(decomposition);
   SCIPdecompGetVarsLabels(decomposition, vars, returnedlabels, NVARS);

   cr_assert_arr_eq(returnedlabels, varlabels, NVARS,
      "Array {%s} not equal to {%s}\n",
      printIntArray(strbuf1, returnedlabels, NVARS),
      printIntArray(strbuf2, varlabels, NVARS)
      );
}

Test(decomptest, test_var_labeling, .description="check variable label computation")
{

   SCIP_CALL( SCIPdecompSetConsLabels(decomp, conss, labels_conss, NCONSS) );

   SCIP_CALL( SCIPcomputeDecompVarsLabels(scip, decomp, conss, NCONSS) );

   checkVarsLabels(decomp, labels_vars);

}

Test(decomptest, test_benders_var_labeling, .description="check variable labelling for Benders' decomposition")
{
   SCIPdecompSetUseBendersLabels(decomp, TRUE);

   SCIP_CALL( SCIPdecompSetConsLabels(decomp, conss, labels_conss, NCONSS) );

   SCIP_CALL( SCIPcomputeDecompVarsLabels(scip, decomp, conss, NCONSS) );

   checkVarsLabels(decomp, benderslabels_vars);

}

Test(decomptest, test_dec_reader, .description="test decomposition reader")
{
   SCIP_DECOMP* scip_decomp;
   SCIP_DECOMP** scip_decomps;
   int n_decomps;
   SCIP_VAR* transvars[NVARS];
   int returnedlabels[NVARS];
   int v;
   SCIP_Bool original = TRUE;

   SCIP_CALL( SCIPreadProb(scip, testdecname, "dec") );

   SCIPgetDecomps(scip, &scip_decomps, &n_decomps, original);
   cr_assert_eq(n_decomps, 1);

   scip_decomp = scip_decomps[0];
   cr_assert_not_null(scip_decomp);

   checkConsLabels(scip_decomp);

   checkVarsLabels(scip_decomp, labels_vars);

   /* solve the problem without presolving */
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsolve(scip) );

   /* transform variable array */
   for( v = 0; v < NVARS; ++v )
   {
      transvars[v] = SCIPvarGetTransVar(vars[v]);
      cr_assert_not_null(transvars[v]);
   }

   /* now get the transformed decomposition and compare its variable labels */
   SCIPgetDecomps(scip, &scip_decomps, &n_decomps, !original);
   cr_assert_eq(n_decomps, 1, "Number of transformed decompositions should be 1.\n");
   scip_decomp = scip_decomps[0];

   SCIPdecompGetVarsLabels(scip_decomp, transvars, returnedlabels, NVARS);

   cr_assert_arr_eq(returnedlabels, labels_vars, NVARS,
      "Arrays should be equal: {%s} != {%s}",
      printIntArray(strbuf1, returnedlabels, NVARS),
      printIntArray(strbuf2, labels_vars, NVARS)
      );
}
