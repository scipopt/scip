/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   network.c
 * @brief  unittests for network matrix detection methods
 * @author Rolf van der Hulst
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_network.h"

#include "include/scip_test.h"

/**
 * Because the algorithm and data structures used for detecting network matrices are rather complex,
 * we extensively test them in this file. We do this by checking if the cycles of the graphs represented by
 * the decomposition match the nonzero entries of the columns of the matrix that we supplied, for many different graphs.
 * Most of the specific graphs tested either contain some special case or posed challenges during development.
 */

static SCIP* scip;

static
void setup(void)
{
   /* create scip */
   SCIP_CALL( SCIPcreate(&scip) );
}

static
void teardown(void)
{
   /* free scip */
   SCIP_CALL( SCIPfree(&scip) );
}

/* CSR/CSC matrix type to encode testing matrices */
typedef struct
{
   int nrows;
   int ncols;
   int nnonzs;

   bool isRowWise;                           /* True -> CSR matrix, False -> CSC matrix */

   int* firstIndex;                          /* Array containing the index of the first nonzero of the row (column)
                                              * for the CSR (CSC) matrix */
   int* entryIndex;                          /* Array containing the entries columns (rows) for the CSR (CSC) matrix */
   double* entryValue;                       /* Array containing the entry values */
} DirectedTestCase;

/**< Create a testcase from a string */
static
DirectedTestCase stringToTestCase(
   const char*           string,             /**< The string to convert to a test case */
   int                   rows,               /**< The number of rows of the matrix to create */
   int                   cols                /**< The number of columns of the matrix to create */
   )
{
   DirectedTestCase testCase;
   testCase.nrows = rows;
   testCase.ncols = cols;
   testCase.nnonzs = 0;
   testCase.isRowWise = true;


   testCase.firstIndex = malloc(sizeof(int) * ( rows + 1 ));

   int nonzeroArraySize = 8;
   testCase.entryIndex = malloc(sizeof(int) * nonzeroArraySize);
   testCase.entryValue = malloc(sizeof(double) * nonzeroArraySize);


   const char* current = string;
   int i = 0;

   while(i < rows * cols && *current != '\0' )
   {
      char* next = NULL;
      double num = strtod(current, &next);
      if( i % cols == 0 )
      {
         testCase.firstIndex[i / cols] = testCase.nnonzs;
      }
      if( num != 0.0 )
      {
         if( testCase.nnonzs == nonzeroArraySize )
         {
            int newSize = nonzeroArraySize * 2;
            testCase.entryValue = realloc(testCase.entryValue, sizeof(double) * newSize);
            testCase.entryIndex = realloc(testCase.entryIndex, sizeof(int) * newSize);
            nonzeroArraySize = newSize;
         }
         testCase.entryValue[testCase.nnonzs] = num;
         testCase.entryIndex[testCase.nnonzs] = i % cols;
         ++testCase.nnonzs;
      }
      current = next;
      ++i;
   }
   testCase.firstIndex[testCase.nrows] = testCase.nnonzs;

   return testCase;
}

/**< Transposes a testcase */
static
void transposeMatrixStorage(
   DirectedTestCase*     testCase            /**< The testcase to transpose (in place) */
   )
{
   int numPrimaryDimension = testCase->isRowWise ? testCase->nrows : testCase->ncols;
   int numSecondaryDimension = testCase->isRowWise ? testCase->ncols : testCase->nrows;

   int* transposedFirstIndex = malloc(sizeof(int) * ( numSecondaryDimension + 1 ));
   int* transposedEntryIndex = malloc(sizeof(int) * testCase->nnonzs);
   double* transposedEntryValue = malloc(sizeof(double) * testCase->nnonzs);

   for( int i = 0; i <= numSecondaryDimension; ++i )
   {
      transposedFirstIndex[i] = 0;
   }
   for( int i = 0; i < testCase->nnonzs; ++i )
   {
      ++( transposedFirstIndex[testCase->entryIndex[i] + 1] );
   }

   for( int i = 1; i < numSecondaryDimension; ++i )
   {
      transposedFirstIndex[i] += transposedFirstIndex[i - 1];
   }

   for( int i = 0; i < numPrimaryDimension; ++i )
   {
      int first = testCase->firstIndex[i];
      int beyond = testCase->firstIndex[i + 1];
      for( int entry = first; entry < beyond; ++entry )
      {
         int index = testCase->entryIndex[entry];
         int transIndex = transposedFirstIndex[index];
         transposedEntryIndex[transIndex] = i;
         transposedEntryValue[transIndex] = testCase->entryValue[entry];
         ++( transposedFirstIndex[index] );
      }
   }
   for( int i = numSecondaryDimension; i > 0; --i )
   {
      transposedFirstIndex[i] = transposedFirstIndex[i - 1];
   }
   transposedFirstIndex[0] = 0;

   free(testCase->entryIndex);
   free(testCase->entryValue);
   free(testCase->firstIndex);

   testCase->firstIndex = transposedFirstIndex;
   testCase->entryIndex = transposedEntryIndex;
   testCase->entryValue = transposedEntryValue;

   testCase->isRowWise = !testCase->isRowWise;
}

/**< Copies a test case */
static
DirectedTestCase copyTestCase(
   DirectedTestCase*     testCase            /**< The test case to copy */
   )
{
   DirectedTestCase copy;
   copy.nrows = testCase->nrows;
   copy.ncols = testCase->ncols;
   copy.nnonzs = testCase->nnonzs;
   copy.isRowWise = testCase->isRowWise;

   int size = ( testCase->isRowWise ? testCase->nrows : testCase->ncols ) + 1;
   copy.firstIndex = malloc(sizeof(int) * size);
   for( int i = 0; i < size; ++i )
   {
      copy.firstIndex[i] = testCase->firstIndex[i];
   }
   copy.entryIndex = malloc(sizeof(int) * testCase->nnonzs);
   copy.entryValue = malloc(sizeof(double) * testCase->nnonzs);

   for( int i = 0; i < testCase->nnonzs; ++i )
   {
      copy.entryIndex[i] = testCase->entryIndex[i];
      copy.entryValue[i] = testCase->entryValue[i];
   }
   return copy;
}

/**< Frees a testcase */
static
void freeTestCase(
   DirectedTestCase*     testCase            /**< The test case to free */
   )
{
   free(testCase->firstIndex);
   free(testCase->entryIndex);
   free(testCase->entryValue);
}

/**< Runs and checks whether a testcase is executed correctly with the network column addition algorithm.
 *  This checks the cycles of the network matrix at every step. If isExpected(Not)Network is set, we check if the
 *  complete matrix is (not) network, and return an error otherwise.
 */
static
SCIP_RETCODE runColumnTestCase(
   DirectedTestCase*     testCase,           /**< The testcase to check */
   bool                  expectedNetwork,    /**< If the complete matrix is not detected to be a network matrix, return an error */
   bool                  expectedNotNetwork  /**< If the complete matrix is detected to be a network matrix, return an error */
   )
{
   if( testCase->isRowWise )
   {
      transposeMatrixStorage(testCase);
   }
   cr_assert(!testCase->isRowWise);
   SCIP_NETMATDEC* dec = NULL;
   BMS_BLKMEM* blkmem = SCIPblkmem(scip);
   BMS_BUFMEM* bufmem = SCIPbuffer(scip);
   SCIP_CALL( SCIPnetmatdecCreate(blkmem, &dec, testCase->nrows, testCase->ncols) );

   SCIP_Bool isNetwork = TRUE;

   int* tempColumnStorage;
   SCIP_Bool* tempSignStorage;

   SCIP_CALL( SCIPallocBufferArray(scip, &tempColumnStorage, testCase->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tempSignStorage, testCase->nrows) );

   for( int i = 0; i < testCase->ncols; ++i )
   {
      int colEntryStart = testCase->firstIndex[i];
      int colEntryEnd = testCase->firstIndex[i + 1];
      int* nonzeroRows = &testCase->entryIndex[colEntryStart];
      double* nonzeroValues = &testCase->entryValue[colEntryStart];
      int nonzeros = colEntryEnd - colEntryStart;
      cr_assert(nonzeros >= 0);
      /* Check if adding the column preserves the network matrix */
      SCIP_CALL( SCIPnetmatdecTryAddCol(dec, i, nonzeroRows, nonzeroValues, nonzeros, &isNetwork) );
      if( !isNetwork )
      {
         break;
      }
      cr_expect(SCIPnetmatdecIsMinimal(dec));
      /* Check if the computed network matrix indeed reflects the network matrix,
       * by checking if the fundamental cycles are all correct
       */
      for( int j = 0; j <= i; ++j )
      {
         int jColEntryStart = testCase->firstIndex[j];
         int jColEntryEnd = testCase->firstIndex[j + 1];
         int* jNonzeroRows = &testCase->entryIndex[jColEntryStart];
         double* jNonzeroValues = &testCase->entryValue[jColEntryStart];
         int jNonzeros = jColEntryEnd - jColEntryStart;
         SCIP_Bool cycleIsCorrect = SCIPnetmatdecVerifyCycle(bufmem, dec, j,
                                                             jNonzeroRows, jNonzeroValues,
                                                             jNonzeros, tempColumnStorage,
                                                             tempSignStorage);

         cr_expect(cycleIsCorrect);
      }
   }


   if( expectedNetwork )
   {
      /* We expect that the given matrix is a network matrix. If not, something went wrong */
      cr_expect(isNetwork);
   }
   if( expectedNotNetwork )
   {
      /* We expect that the given matrix is not a network matrix. If not, something went wrong */
      cr_expect(!isNetwork);
   }
   SCIPfreeBufferArray(scip, &tempColumnStorage);
   SCIPfreeBufferArray(scip, &tempSignStorage);

   SCIPnetmatdecFree(&dec);
   cr_assert(dec == NULL);

   return SCIP_OKAY;
}

/**< Runs and checks whether a testcase is executed correctly with the network row addition algorithm.
 *  This checks the cycles of the network matrix at every step. If isExpected(Not)Network is set, we check if the
 *  complete matrix is (not) network, and return an error otherwise.
 */
static
SCIP_RETCODE runRowTestCase(
   DirectedTestCase*     testCase,           /**< The testcase to check */
   bool                  expectedNetwork,    /**< If the complete matrix is not detected to be a network matrix, return an error */
   bool                  expectedNotNetwork  /**< If the complete matrix is detected to be a network matrix, return an error */
   )
{
   if( !testCase->isRowWise )
   {
      transposeMatrixStorage(testCase);
   }
   cr_assert(testCase->isRowWise);

   /* We keep a column-wise copy to check the columns easily */
   DirectedTestCase colWiseCase = copyTestCase(testCase);
   transposeMatrixStorage(&colWiseCase);

   BMS_BLKMEM* blkmem = SCIPblkmem(scip);
   BMS_BUFMEM* bufmem = SCIPbuffer(scip);

   SCIP_NETMATDEC* dec = NULL;
   SCIP_CALL( SCIPnetmatdecCreate(blkmem, &dec, testCase->nrows, testCase->ncols) );

   SCIP_Bool isNetwork = TRUE;

   int* tempColumnStorage;
   SCIP_Bool* tempSignStorage;

   SCIP_CALL( SCIPallocBufferArray(scip, &tempColumnStorage, testCase->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tempSignStorage, testCase->nrows) );

   for( int i = 0; i < testCase->nrows; ++i )
   {
      int rowEntryStart = testCase->firstIndex[i];
      int rowEntryEnd = testCase->firstIndex[i + 1];
      int* nonzeroCols = &testCase->entryIndex[rowEntryStart];
      double* nonzeroValues = &testCase->entryValue[rowEntryStart];
      int nonzeros = rowEntryEnd - rowEntryStart;
      cr_assert(nonzeros >= 0);
      /* Check if adding the row preserves the network matrix */
      SCIP_CALL( SCIPnetmatdecTryAddRow(dec, i, nonzeroCols, nonzeroValues, nonzeros, &isNetwork) );
      if( !isNetwork )
      {
         break;
      }
      cr_expect(SCIPnetmatdecIsMinimal(dec));
      /* Check if the computed network matrix indeed reflects the network matrix,
       * by checking if the fundamental cycles are all correct
       */
      for( int j = 0; j < colWiseCase.ncols; ++j )
      {
         int jColEntryStart = colWiseCase.firstIndex[j];
         int jColEntryEnd = colWiseCase.firstIndex[j + 1];

         /* Count the number of rows in the column that should be in the current decomposition */
         int finalEntryIndex = jColEntryStart;
         for( int testEntry = jColEntryStart; testEntry < jColEntryEnd; ++testEntry )
         {
            if( colWiseCase.entryIndex[testEntry] <= i )
            {
               ++finalEntryIndex;
            }
            else
            {
               break;
            }
         }

         int* jNonzeroRows = &colWiseCase.entryIndex[jColEntryStart];
         double* jNonzeroValues = &colWiseCase.entryValue[jColEntryStart];

         int jNonzeros = finalEntryIndex - jColEntryStart;
         SCIP_Bool cycleIsCorrect = SCIPnetmatdecVerifyCycle(bufmem, dec, j,
                                                             jNonzeroRows, jNonzeroValues,
                                                             jNonzeros, tempColumnStorage,
                                                             tempSignStorage);

         cr_expect(cycleIsCorrect);
      }
   }

   if( expectedNetwork )
   {
      /* We expect that the given matrix is a network matrix. If not, something went wrong */
      cr_expect(isNetwork);
   }
   if( expectedNotNetwork )
   {
      /* We expect that the given matrix is not a network matrix. If not, something went wrong */
      cr_expect(!isNetwork);
   }

   freeTestCase(&colWiseCase);

   SCIPfreeBufferArray(scip, &tempColumnStorage);
   SCIPfreeBufferArray(scip, &tempSignStorage);

   SCIPnetmatdecFree(&dec);
   cr_assert(dec == NULL);

   return SCIP_OKAY;
}

/**< Runs the network row addition, and attempts to construct the graph.
 * This functions main purpose is to check if SCIPnetmatdecCreateDiGraph() correctly creates the graph without errors.
 */
static
SCIP_RETCODE runRowTestCaseGraph(
   DirectedTestCase*     testCase            /**< The testcase to check */
   )
{
   if( !testCase->isRowWise )
   {
      transposeMatrixStorage(testCase);
   }
   cr_assert(testCase->isRowWise);

   /* We keep a column-wise copy to check the columns easily */
   DirectedTestCase colWiseCase = copyTestCase(testCase);
   transposeMatrixStorage(&colWiseCase);

   BMS_BLKMEM* blkmem = SCIPblkmem(scip);

   SCIP_NETMATDEC* dec = NULL;
   SCIP_CALL( SCIPnetmatdecCreate(blkmem, &dec, testCase->nrows, testCase->ncols) );

   SCIP_Bool isNetwork = TRUE;

   for( int i = 0; i < testCase->nrows; ++i )
   {
      int rowEntryStart = testCase->firstIndex[i];
      int rowEntryEnd = testCase->firstIndex[i + 1];
      int* nonzeroCols = &testCase->entryIndex[rowEntryStart];
      double* nonzeroValues = &testCase->entryValue[rowEntryStart];
      int nonzeros = rowEntryEnd - rowEntryStart;
      cr_assert(nonzeros >= 0);
      /* Check if adding the row preserves the network matrix */
      SCIP_CALL( SCIPnetmatdecTryAddRow(dec, i, nonzeroCols, nonzeroValues, nonzeros, &isNetwork) );
      cr_assert(isNetwork);
   }
   SCIP_DIGRAPH* graph;
   SCIP_CALL( SCIPnetmatdecCreateDiGraph(dec, blkmem, &graph, TRUE) );

   SCIPdigraphPrint(graph, SCIPgetMessagehdlr(scip), stdout);
   SCIPdigraphFree(&graph);
   freeTestCase(&colWiseCase);

   SCIPnetmatdecFree(&dec);
   cr_assert(dec == NULL);

   return SCIP_OKAY;
}

TestSuite(network, .init = setup, .fini = teardown);

Test(network, coladd_single_column, .description = "Try adding a single column")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 "
      "+1 "
      "-1 ",
      3, 1);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_doublecolumn_invalid_sign, .description = "Try adding a second column that has invalid signing")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      "+1  0 "
      "-1 +1 ",
      3, 2);
   runColumnTestCase(&testCase, false, true);
   freeTestCase(&testCase);
}

Test(network, coladd_doublecolumn_invalid_sign_2, .description = "Try adding a second column that has invalid signing")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1  0 "
      "-1 -1 ",
      3, 2);
   runColumnTestCase(&testCase, false, true);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_1, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_2, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_1r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      "+1  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_2r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      "+1  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_3, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1  0 "
      " 0 +1 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_4, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      "+1  0 "
      " 0 +1 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_3r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      "+1  0 "
      " 0 +1 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_4r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      "+1  0 "
      " 0 +1 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_5, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      " 0  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_6, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      " 0  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_7, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      " 0  0 "
      " 0  +1 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_8, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      " 0  0 "
      " 0  +1 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_5r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      " 0  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_6r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      " 0  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_7r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      " 0  0 "
      " 0  +1 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_8r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      " 0  0 "
      " 0  +1 ",
      3, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_9, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1  0 "
      "-1  +1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_10, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      "+1  0 "
      "-1  -1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_11, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      "+1  0 "
      "-1  -1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_12, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1  0 "
      "-1  +1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_9r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      "-1  0 "
      "+1  +1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_10r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      "-1  0 "
      "+1  -1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_11r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      "-1  0 "
      "+1  -1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_12r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      "-1  0 "
      "+1  +1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_13, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      "+1 +1 "
      "-1 -1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_14, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1 -1 "
      "-1 +1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_15, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      "+1 +1 "
      "-1 -1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_16, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1 -1 "
      "-1 +1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_13r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      "-1 +1 "
      "+1 -1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_14r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      "-1 -1 "
      "+1 +1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_15r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      "-1 +1 "
      "+1 -1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_splitseries_16r, .description = "Split a series component")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      "-1 -1 "
      "+1 +1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_parallelsimple_1, .description = "Extending a parallel component")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 1 -1 ",
      1, 4);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_parallelsimple_2, .description = "Extending a parallel component")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 1 -1 "
      "0 0 -1 0 ",
      2, 4);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_parallelsimple_3, .description = "Extending a parallel component")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 1 1 "
      "0 0 1 0 ",
      2, 4);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_components_1, .description = "Merging multiple components into one")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 -1 "
      "1 0 1 ",
      2, 3);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_components_2, .description = "Merging multiple components into one")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 -1 "
      "1  0 1 ",
      2, 3);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_components_3, .description = "Merging multiple components into one")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1  1 "
      "-1  0 1 ",
      2, 3);
   runColumnTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by3_1, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 "
      "1 -1 -1 "
      "-1 1 -1 ",
      3, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by3_2, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 -1 "
      "1 0 1 "
      "0 0 0 ",
      3, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by3_3, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 1 "
      "1 -1 -1 "
      "0 0 1 ",
      3, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by3_4, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 1 "
      "-1 1 0 "
      "0 1 -1 ",
      3, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by3_5, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 "
      "0 1 0 "
      "-1 -1 -1 ",
      3, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by3_6, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 1 -1 "
      "-1 1 -1 "
      "-1 0 -1 ",
      3, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by3_7, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 -1 "
      "0 1 1 "
      "-1 0 0 ",
      3, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by3_8, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 -1 "
      "0 -1 1 "
      "1 0 0 ",
      3, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by3_9, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 0 "
      "-1 -1 -1 "
      "-1 -1 -1 ",
      3, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by3_10, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 1 -1 "
      "-1 1 -1 "
      "1 0 -1 ",
      3, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by3_11, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 1 0 "
      "-1 0 1 "
      "0 1 1 ",
      3, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_6by3_1, .description = "A six by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 1 "
      "-1 1 -1 "
      "-1 1 0 "
      "0 1 1 "
      "1 0 -1 "
      "0 -1 -1 ",
      6, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_6by3_2, .description = "A six by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 1 "
      "-1 1 0 "
      "0 1 -1 "
      "0 -1 0 "
      "0 1 0 "
      "0 0 1 ",
      6, 3);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by4_1, .description = "A three by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 1 1 "
      "-1 -1 0 0 "
      "1 0 1 -1 ",
      3, 4);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_3by5_1, .description = "A three by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 0 1 "
      "0 -1 -1 -1 -1 "
      "1 -1 0 -1 -1 ",
      3, 5);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_4by8_1, .description = "A four by eight case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 0 1 1 0 1 0 1 "
      "0 0 0 0 0 -1 1 -1 "
      "1 -1 0 0 1 1 0 0 "
      "0 0 1 1 -1 0 -1 0 ",
      4, 8);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_4by8_2, .description = "A four by eight case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 0 -1 0 0 -1 0 "
      "1 1 0 0 1 1 -1 1 "
      "0 -1 0 0 -1 -1 0 -1 "
      "0 1 1 -1 1 0 -1 -1 ",
      4, 8);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_4by8_3, .description = "A four by eight case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 1 -1 0 -1 0 -1 0 "
      "-1 0 -1 0 1 -1 1 1 "
      "0 0 -1 1 -1 -1 0 0 "
      "-1 1 0 -1 0 -1 1 -1 ",
      4, 8);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_4by8_4, .description = "A four by eight case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 1 1 0 1 0 0 "
      "0 -1 0 -1 1 0 -1 -1 "
      "-1 -1 1 0 1 1 -1 -1 "
      "0 0 0 -1 0 -1 1 -1 ",
      4, 8);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_4by8_5, .description = "A four by eight case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 0 0 -1 -1 -1 0 -1 "
      "-1 -1 0 0 1 1 -1 0 "
      "0 0 1 0 -1 -1 0 -1 "
      "0 0 1 -1 -1 0 0 -1 ",
      4, 8);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_4by8_6, .description = "A four by eight case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 0 1 1 1 0 -1 "
      "0 -1 -1 0 0 1 0 -1 "
      "0 0 1 -1 0 0 1 1 "
      "0 -1 -1 1 1 1 -1 0 ",
      4, 8);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_4by8_7, .description = "A four by eight case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 0 1 1 1 0 1 0 "
      "0 -1 1 1 1 1 -1 0 "
      "1 0 0 0 1 0 1 1 "
      "1 -1 1 0 -1 -1 -1 -1 ",
      4, 8);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_4by4_1, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 -1 0 "
      "-1 0 -1 0 "
      "0 0 -1 +1 "
      "-1 +1 0 -1",
      4, 4);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_5by5_1, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 0 0 "
      "0 -1 -1 1 0 "
      "0 0 -1 1 0 "
      "-1 -1 0 0 -1 "
      "-1 -1 -1 1 0 ",
      5, 5);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_5by5_2, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 1 1 -1 "
      "-1 1 1 1 0 "
      "0 0 1 1 1 "
      "-1 1 0 -1 0 "
      "-1 1 0 0 -1 ",
      5, 5);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_5by5_3, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 0 -1 "
      "0 -1 0 1 -1 "
      "1 1 1 0 1 "
      "0 0 1 1 0 "
      "-1 0 -1 0 1 ",
      5, 5);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_5by5_4, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 -1 1 0 "
      "0 0 -1 1 -1 "
      "1 -1 0 0 1 "
      "0 1 1 -1 0 "
      "0 -1 0 1 0 ",
      5, 5);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_5by5_5, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 0 1 1 "
      "-1 1 1 0 0 "
      "0 0 0 -1 -1 "
      "1 -1 0 0 -1 "
      "-1 1 0 1 1 ",
      5, 5);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_5by5_6, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 1 0 0 "
      "-1 0 -1 0 0 "
      "0 -1 1 1 -1 "
      "-1 1 -1 -1 0 "
      "1 0 0 -1 0 ",
      5, 5);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_5by5_7, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 0 1 0 0 "
      "0 -1 0 1 0 "
      "0 1 -1 -1 0 "
      "-1 0 -1 -1 -1 "
      "-1 -1 0 1 -1 ",
      5, 5);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_5by5_8, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 1 -1 0 "
      "1 0 0 0 -1 "
      "-1 0 -1 1 1 "
      "0 0 0 -1 -1 "
      "1 1 0 -1 0 ",
      5, 5);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, coladd_5by5_9, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 0 -1 1 "
      "1 -1 -1 -1 1 "
      "0 0 -1 0 0 "
      "0 -1 -1 0 0 "
      "1 0 0 0 1 ",
      5, 5);
   runColumnTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_1by2_1, .description = "A one by two case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 ",
      1, 2);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_1by2_2, .description = "A one by two case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 ",
      1, 2);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_1by2_3, .description = "A one by two case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 ",
      1, 2);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_2by3_1, .description = "A two by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 +1 "
      "-1 +1 -1 ",
      2, 3);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_2by3_2, .description = "A two by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 +1 "
      "-1 +1 +1 ",
      2, 3);
   runRowTestCase(&testCase, false, true);
   freeTestCase(&testCase);
}

Test(network, rowadd_2by3_3, .description = "A two by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 +1 "
      "+1 0 +1 ",
      2, 3);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_2by3_4, .description = "A two by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 "
      "+1 0  0 ",
      2, 3);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_2by3_5, .description = "A two by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 "
      "+0 +1  0 ",
      2, 3);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_2by3_6, .description = "A two by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 "
      "+0 -1  0 ",
      2, 3);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_2by3_7, .description = "A two by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 "
      "+0 +1 +1 ",
      2, 3);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_2by3_8, .description = "A two by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 "
      "+0 -1 +1 ",
      2, 3);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by6_1, .description = "A three by six case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 0 0 0 "
      "0 0 +1 -1 0 0 "
      "-1 +1 -1 0 0 0 ",
      3, 6);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by6_2, .description = "A three by six case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 0 0 0 "
      "0 0 +1 -1 0 0 "
      "-1 +1 -1 0 0 +1 ",
      3, 6);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by2_1, .description = "A three by two case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "-1 +1 "
      "+1 -1 ",
      3, 2);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by1_1, .description = "A three by one case")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 "
      "-1 "
      "+1 ",
      3, 1);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by3_1, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 1 "
      "-1 1 0 "
      "0 1 -1 ",
      3, 3);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by3_2, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 "
      "0 1 0 "
      "-1 -1 -1 ",
      3, 3);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by3_3, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 0 "
      "-1 0 -1 "
      "-1 -1 0 ",
      3, 3);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by3_4, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 1 "
      "0 -1 1 "
      "-1 -1 -1 ",
      3, 3);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by3_5, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 -1 "
      "-1 -1 0 "
      "0 1 -1 ",
      3, 3);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by3_6, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 "
      "1 1 0 "
      "-1 1 -1 ",
      3, 3);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by3_7, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 0 "
      "0 1 -1 "
      "-1 1 -1 ",
      3, 3);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by3_8, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 "
      "-1 1 0 "
      "-1 0 -1 ",
      3, 3);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by3_9, .description = "A three by three case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 -1 "
      "-1 -1 -1 "
      "1 0 1 ",
      3, 3);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by6_3, .description = "A three by six case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 0 0 0 -1 "
      "0 0 0 -1 -1 -1 "
      "-1 1 0 0 1 1 ",
      3, 3);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_3by6_4, .description = "A three by six case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 -1 -1 0 0 "
      "0 0 0 1 -1 -1 "
      "0 -1 1 1 -1 0 ",
      3, 3);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_1, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 0 0 "
      "-1 -1 -1 -1 "
      "0 -1 -1 -1 "
      "1 0 0 -1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_2, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 0 -1 -1 "
      "1 1 -1 -1 "
      "0 0 -1 -1 "
      "0 -1 0 1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_3, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 0 1 "
      "-1 0 -1 0 "
      "0 0 1 1 "
      "0 -1 -1 1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_4, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 0 -1 0 "
      "0 1 0 0 "
      "-1 -1 1 1 "
      "0 -1 -1 1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_5, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 -1 1 "
      "1 0 0 -1 "
      "-1 0 -1 0 "
      "1 1 1 0 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_6, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 0 -1 1 "
      "0 1 -1 0 "
      "0 -1 1 -1 "
      "-1 -1 0 1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_7, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 0 -1 1 "
      "0 1 -1 0 "
      "0 -1 1 -1 "
      "-1 -1 0 1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_8, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 1 1 "
      "1 -1 -1 0 "
      "-1 1 1 1 "
      "1 0 0 -1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_9, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 0 0 "
      "-1 1 0 -1 "
      "0 1 1 -1 "
      "-1 0 1 0 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_10, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 1 -1 "
      "-1 0 0 -1 "
      "0 1 0 1 "
      "0 -1 0 0 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_11, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 1 0 "
      "-1 -1 0 1 "
      "-1 0 0 1 "
      "0 0 1 1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_12, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 0 -1 "
      "0 1 0 1 "
      "-1 0 1 0 "
      "0 0 -1 -1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_13, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 -1 1 "
      "-1 0 0 1 "
      "1 1 1 -1 "
      "1 0 -1 -1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_14, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 1 1 "
      "1 1 0 1 "
      "0 1 1 0 "
      "0 1 1 1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_15, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 -1 0 "
      "1 1 0 1 "
      "1 0 -1 0 "
      "0 1 0 1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_16, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 1 0 "
      "-1 0 0 -1 "
      "0 -1 0 1 "
      "-1 -1 0 0 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_17, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 -1 0 "
      "-1 -1 0 0 "
      "0 -1 0 1 "
      "0 1 1 -1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_4by4_18, .description = "A four by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 0 1 "
      "1 0 -1 0 "
      "-1 0 1 1 "
      "-1 1 0 -1 ",
      4, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_5by5_1, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 1 -1 1 "
      "0 0 1 -1 0 "
      "-1 0 0 1 0 "
      "0 0 1 0 1 "
      "0 -1 0 1 1 ",
      5, 5);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_5by5_2, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 0 -1 -1 "
      "0 -1 0 -1 -1 "
      "0 0 1 -1 0 "
      "0 -1 0 -1 -1 "
      "-1 1 1 0 0 ",
      5, 5);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_5by5_3, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 1 0 -1 "
      "0 1 1 0 -1 "
      "1 0 -1 0 1 "
      "-1 0 1 0 0 "
      "1 0 -1 0 1 ",
      5, 5);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_5by5_4, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 1 0 0 "
      "0 1 -1 1 0 "
      "0 -1 1 0 0 "
      "1 -1 1 0 0 "
      "0 1 0 1 1 ",
      5, 5);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_5by5_5, .description = "A five by five case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 0 1 0 1 "
      "1 0 0 1 -1 "
      "1 -1 1 1 0 "
      "0 0 -1 0 -1 "
      "0 0 1 0 1 ",
      5, 5);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_8by4, .description = "A eight by four case")
{
   DirectedTestCase testCase = stringToTestCase(
      "0 0 0 0 "
      "1 0 1 0 "
      "-1 1 -1 -1 "
      "1 0 1 1 "
      "1 -1 1 0 "
      "1 -1 0 0 "
      "1 1 -1 1 "
      "0 0 1 0 ",
      8, 4);
   runRowTestCase(&testCase, false, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_singlerigid_1, .description = "Updating a single rigid member")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "+1 +1 0 "
      "0 -1 +1 "
      "+1 +1 0 ",
      4, 3);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_singlerigid_2, .description = "Updating a single rigid member")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "+1 +1 0 "
      "0 -1 +1 "
      "-1 -1 0 ",
      4, 3);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_singlerigid_3, .description = "Updating a single rigid member")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "+1 +1 0 "
      "0 -1 +1 "
      "-1 +1 0 ",
      4, 3);
   runRowTestCase(&testCase, false, true);
   freeTestCase(&testCase);
}

Test(network, rowadd_singlerigid_4, .description = "Updating a single rigid member")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "-1 -1 -1 "
      "0 +1 +1 "
      "+1 +1 +1 ",
      4, 3);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_singlerigid_5, .description = "Updating a single rigid member")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "-1 -1 -1 "
      "0 +1 +1 "
      "-1 -1 -1 ",
      4, 3);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_singlerigid_6, .description = "Updating a single rigid member")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "-1 -1 -1 "
      "0 +1 +1 "
      "-1 +1 -1 ",
      4, 3);
   runRowTestCase(&testCase, false, true);
   freeTestCase(&testCase);
}

Test(network, rowadd_singlerigid_7, .description = "Updating a single rigid member")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 0 0 +1 "
      "+1 0 +1 0 0 "
      "0 -1 +1 +1 -1 "
      "0 0 0 -1 +1 "
      "+1 +1 0 0 0 ",
      5, 5);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}

Test(network, rowadd_singlerigid_8, .description = "Updating a single rigid member")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 0 0 +1 "
      "+1 0 +1 0 0 "
      "0 -1 +1 +1 -1 "
      "0 0 0 -1 +1 "
      "+1 +1 0 0 0 "
      "+1 0 +1 +1 0 ",
      6, 5);
   runRowTestCase(&testCase, true, false);
   freeTestCase(&testCase);
}
/* TODO: test interleaved addition, test using random sampling + test erdos-renyi generated graphs */

Test(network, rowadd_singlerigid_graph, .description = "Computing the graph for a single rigid member")
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "+1 +1 0 "
      "0 -1 +1 "
      "+1 +1 0 ",
      4, 3);
   runRowTestCaseGraph(&testCase);
   freeTestCase(&testCase);
}
