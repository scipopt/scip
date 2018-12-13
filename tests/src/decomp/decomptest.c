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
#include "scip/decomp.h"

#include "include/scip_test.h"

const char* cip_input =
"STATISTICS"
"  Problem name     : mymip.cip"
"  Variables        : 2 (0 binary, 2 integer, 0 implicit integer, 0 continuous)"
"  Constraints      : 0 initial, 1 maximal"
"OBJECTIVE"
"  Sense            : minimize"
"VARIABLES"
"  [integer] <x1>: obj=1, original bounds=[0,10000]"
"  [integer] <x2>: obj=1, original bounds=[0,10000]"
"CONSTRAINTS"
"  [linear] <lincons> : 3057<x1> +12227<x2> == 8908099.5;"
"END";

static const char* testfilename = "decomptest.cip";

/** GLOBAL VARIABLES **/
static SCIP* scip;
static SCIP_DECOMP* decomp;

/* TEST SUITE */
static
void setup(void)
{
   FILE* file = fopen(testfilename, "w");
   printf("Print file");
   fprintf(file, "STATISTICS"
      "  Problem name     : mymip.cip\n"
      "  Variables        : 5 (4 binary, 0 integer, 0 implicit integer, 1 continuous)\n"
      "  Constraints      : 3 initial, 3 maximal\n"
      "OBJECTIVE\n"
      "  Sense            : minimize\n"
      "VARIABLES\n"
      "  [binary]  <x1>: obj=1, original bounds=[0,1]\n"
      "  [binary]  <x2>: obj=1, original bounds=[0,1]\n"
      "  [binary]  <x3>: obj=1, original bounds=[0,1]\n"
      "  [binary]  <x4>: obj=1, original bounds=[0,1]\n"
      "  [continuous]  <y>: obj=-1, original bounds=[0,100]\n"
      "CONSTRAINTS\n"
      "  [linear] <linkingcons> : 10<x1> +20<x2> +30<x3> +40<x4> -1<y> >= 0;\n"
      "  [linear] <block1cons> : 1<x1> +<x2> == 1;\n"
      "  [linear] <block2cons> : 1<x3> +<x4> == 1;\n"
      "END\n");
   fclose(file);

   SCIPcreate(&scip);
}

static
void teardown(void)
{
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
   SCIP_CALL( SCIPdecompCreate(&decomp, SCIPblkmem(scip)) );

   SCIPdecompFree(&decomp, SCIPblkmem(scip));
}

Test(decomptest, read_file, .description = "treat constant string as file input")
{
   FILE* myfile = fmemopen((void*)cip_input, strlen(cip_input) + 1, "r");
   cr_assert(myfile != NULL);

   fclose(myfile);
}
