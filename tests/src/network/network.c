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

   SCIP_Bool isRowWise;                           /* True -> CSR matrix, False -> CSC matrix */

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
   testCase.isRowWise = TRUE;


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
   SCIP_Bool                  expectedNetwork,    /**< If the complete matrix is not detected to be a network matrix, return an error */
   SCIP_Bool                  expectedNotNetwork  /**< If the complete matrix is detected to be a network matrix, return an error */
   )
{
   if( testCase->isRowWise )
   {
      transposeMatrixStorage(testCase);
   }
   TEST_ASSERT(!testCase->isRowWise);
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
      TEST_ASSERT(nonzeros >= 0);
      /* Check if adding the column preserves the network matrix */
      SCIP_CALL( SCIPnetmatdecTryAddCol(dec, i, nonzeroRows, nonzeroValues, nonzeros, &isNetwork) );
      if( !isNetwork )
      {
         break;
      }
      SOFT_ASSERT(SCIPnetmatdecIsMinimal(dec));
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

         SOFT_ASSERT(cycleIsCorrect);
      }
   }


   if( expectedNetwork )
   {
      /* We expect that the given matrix is a network matrix. If not, something went wrong */
      SOFT_ASSERT(isNetwork);
   }
   if( expectedNotNetwork )
   {
      /* We expect that the given matrix is not a network matrix. If not, something went wrong */
      SOFT_ASSERT(!isNetwork);
   }
   SCIPfreeBufferArray(scip, &tempColumnStorage);
   SCIPfreeBufferArray(scip, &tempSignStorage);

   SCIPnetmatdecFree(&dec);
   TEST_ASSERT(dec == NULL);

   return SCIP_OKAY;
}

/**< Runs and checks whether a testcase is executed correctly with the network row addition algorithm.
 *  This checks the cycles of the network matrix at every step. If isExpected(Not)Network is set, we check if the
 *  complete matrix is (not) network, and return an error otherwise.
 */
static
SCIP_RETCODE runRowTestCase(
   DirectedTestCase*     testCase,           /**< The testcase to check */
   SCIP_Bool                  expectedNetwork,    /**< If the complete matrix is not detected to be a network matrix, return an error */
   SCIP_Bool                  expectedNotNetwork  /**< If the complete matrix is detected to be a network matrix, return an error */
   )
{
   if( !testCase->isRowWise )
   {
      transposeMatrixStorage(testCase);
   }
   TEST_ASSERT(testCase->isRowWise);

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
      TEST_ASSERT(nonzeros >= 0);
      /* Check if adding the row preserves the network matrix */
      SCIP_CALL( SCIPnetmatdecTryAddRow(dec, i, nonzeroCols, nonzeroValues, nonzeros, &isNetwork) );
      if( !isNetwork )
      {
         break;
      }
      SOFT_ASSERT(SCIPnetmatdecIsMinimal(dec));
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

         SOFT_ASSERT(cycleIsCorrect);
      }
   }

   if( expectedNetwork )
   {
      /* We expect that the given matrix is a network matrix. If not, something went wrong */
      SOFT_ASSERT(isNetwork);
   }
   if( expectedNotNetwork )
   {
      /* We expect that the given matrix is not a network matrix. If not, something went wrong */
      SOFT_ASSERT(!isNetwork);
   }

   freeTestCase(&colWiseCase);

   SCIPfreeBufferArray(scip, &tempColumnStorage);
   SCIPfreeBufferArray(scip, &tempSignStorage);

   SCIPnetmatdecFree(&dec);
   TEST_ASSERT(dec == NULL);

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
   TEST_ASSERT(testCase->isRowWise);

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
      TEST_ASSERT(nonzeros >= 0);
      /* Check if adding the row preserves the network matrix */
      SCIP_CALL( SCIPnetmatdecTryAddRow(dec, i, nonzeroCols, nonzeroValues, nonzeros, &isNetwork) );
      TEST_ASSERT(isNetwork);
   }
   SCIP_DIGRAPH* graph;
   SCIP_CALL( SCIPnetmatdecCreateDiGraph(dec, blkmem, &graph, TRUE) );

   SCIPdigraphPrint(graph, SCIPgetMessagehdlr(scip), stdout);
   SCIPdigraphFree(&graph);
   freeTestCase(&colWiseCase);

   SCIPnetmatdecFree(&dec);
   TEST_ASSERT(dec == NULL);

   return SCIP_OKAY;
}


/** @brief Try adding a single column */
void test_network_coladd_single_column(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 "
      "+1 "
      "-1 ",
      3, 1);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Try adding a second column that has invalid signing */
void test_network_coladd_doublecolumn_invalid_sign(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      "+1  0 "
      "-1 +1 ",
      3, 2);
   runColumnTestCase(&testCase, FALSE, TRUE);
   freeTestCase(&testCase);
}

/** @brief Try adding a second column that has invalid signing */
void test_network_coladd_doublecolumn_invalid_sign_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1  0 "
      "-1 -1 ",
      3, 2);
   runColumnTestCase(&testCase, FALSE, TRUE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_1r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      "+1  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_2r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      "+1  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1  0 "
      " 0 +1 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_4(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      "+1  0 "
      " 0 +1 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_3r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      "+1  0 "
      " 0 +1 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_4r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      "+1  0 "
      " 0 +1 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_5(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      " 0  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_6(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      " 0  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_7(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      " 0  0 "
      " 0  +1 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_8(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      " 0  0 "
      " 0  +1 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_5r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      " 0  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_6r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      " 0  0 "
      " 0  0 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_7r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      " 0  0 "
      " 0  +1 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_8r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      " 0  0 "
      " 0  +1 ",
      3, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_9(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1  0 "
      "-1  +1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_10(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      "+1  0 "
      "-1  -1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_11(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      "+1  0 "
      "-1  -1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_12(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1  0 "
      "-1  +1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_9r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      "-1  0 "
      "+1  +1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_10r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      "-1  0 "
      "+1  -1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_11r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      "-1  0 "
      "+1  -1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_12r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      "-1  0 "
      "+1  +1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_13(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      "+1 +1 "
      "-1 -1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_14(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1 -1 "
      "-1 +1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_15(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 "
      "+1 +1 "
      "-1 -1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_16(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "+1 -1 "
      "-1 +1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_13r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      "-1 +1 "
      "+1 -1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_14r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      "-1 -1 "
      "+1 +1 "
      " 0  0 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_15r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 "
      "-1 +1 "
      "+1 -1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Split a series component */
void test_network_coladd_splitseries_16r(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 "
      "-1 -1 "
      "+1 +1 "
      " 0  +1 ",
      4, 2);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Extending a parallel component */
void test_network_coladd_parallelsimple_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 1 -1 ",
      1, 4);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Extending a parallel component */
void test_network_coladd_parallelsimple_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 1 -1 "
      "0 0 -1 0 ",
      2, 4);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Extending a parallel component */
void test_network_coladd_parallelsimple_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 1 1 "
      "0 0 1 0 ",
      2, 4);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Merging multiple components into one */
void test_network_coladd_components_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 -1 "
      "1 0 1 ",
      2, 3);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Merging multiple components into one */
void test_network_coladd_components_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 -1 "
      "1  0 1 ",
      2, 3);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Merging multiple components into one */
void test_network_coladd_components_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1  1 "
      "-1  0 1 ",
      2, 3);
   runColumnTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_coladd_3by3_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 "
      "1 -1 -1 "
      "-1 1 -1 ",
      3, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_coladd_3by3_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 -1 "
      "1 0 1 "
      "0 0 0 ",
      3, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_coladd_3by3_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 1 "
      "1 -1 -1 "
      "0 0 1 ",
      3, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_coladd_3by3_4(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 1 "
      "-1 1 0 "
      "0 1 -1 ",
      3, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_coladd_3by3_5(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 "
      "0 1 0 "
      "-1 -1 -1 ",
      3, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_coladd_3by3_6(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 1 -1 "
      "-1 1 -1 "
      "-1 0 -1 ",
      3, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_coladd_3by3_7(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 -1 "
      "0 1 1 "
      "-1 0 0 ",
      3, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_coladd_3by3_8(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 -1 "
      "0 -1 1 "
      "1 0 0 ",
      3, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_coladd_3by3_9(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 0 "
      "-1 -1 -1 "
      "-1 -1 -1 ",
      3, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_coladd_3by3_10(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 1 -1 "
      "-1 1 -1 "
      "1 0 -1 ",
      3, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_coladd_3by3_11(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 1 0 "
      "-1 0 1 "
      "0 1 1 ",
      3, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A six by three case */
void test_network_coladd_6by3_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 1 "
      "-1 1 -1 "
      "-1 1 0 "
      "0 1 1 "
      "1 0 -1 "
      "0 -1 -1 ",
      6, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A six by three case */
void test_network_coladd_6by3_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 1 "
      "-1 1 0 "
      "0 1 -1 "
      "0 -1 0 "
      "0 1 0 "
      "0 0 1 ",
      6, 3);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by four case */
void test_network_coladd_3by4_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 1 1 "
      "-1 -1 0 0 "
      "1 0 1 -1 ",
      3, 4);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by five case */
void test_network_coladd_3by5_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 0 1 "
      "0 -1 -1 -1 -1 "
      "1 -1 0 -1 -1 ",
      3, 5);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by eight case */
void test_network_coladd_4by8_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 0 1 1 0 1 0 1 "
      "0 0 0 0 0 -1 1 -1 "
      "1 -1 0 0 1 1 0 0 "
      "0 0 1 1 -1 0 -1 0 ",
      4, 8);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by eight case */
void test_network_coladd_4by8_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 0 -1 0 0 -1 0 "
      "1 1 0 0 1 1 -1 1 "
      "0 -1 0 0 -1 -1 0 -1 "
      "0 1 1 -1 1 0 -1 -1 ",
      4, 8);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by eight case */
void test_network_coladd_4by8_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 1 -1 0 -1 0 -1 0 "
      "-1 0 -1 0 1 -1 1 1 "
      "0 0 -1 1 -1 -1 0 0 "
      "-1 1 0 -1 0 -1 1 -1 ",
      4, 8);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by eight case */
void test_network_coladd_4by8_4(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 1 1 0 1 0 0 "
      "0 -1 0 -1 1 0 -1 -1 "
      "-1 -1 1 0 1 1 -1 -1 "
      "0 0 0 -1 0 -1 1 -1 ",
      4, 8);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by eight case */
void test_network_coladd_4by8_5(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 0 0 -1 -1 -1 0 -1 "
      "-1 -1 0 0 1 1 -1 0 "
      "0 0 1 0 -1 -1 0 -1 "
      "0 0 1 -1 -1 0 0 -1 ",
      4, 8);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by eight case */
void test_network_coladd_4by8_6(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 0 1 1 1 0 -1 "
      "0 -1 -1 0 0 1 0 -1 "
      "0 0 1 -1 0 0 1 1 "
      "0 -1 -1 1 1 1 -1 0 ",
      4, 8);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by eight case */
void test_network_coladd_4by8_7(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 0 1 1 1 0 1 0 "
      "0 -1 1 1 1 1 -1 0 "
      "1 0 0 0 1 0 1 1 "
      "1 -1 1 0 -1 -1 -1 -1 ",
      4, 8);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_coladd_4by4_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 +1 -1 0 "
      "-1 0 -1 0 "
      "0 0 -1 +1 "
      "-1 +1 0 -1",
      4, 4);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_coladd_5by5_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 0 0 "
      "0 -1 -1 1 0 "
      "0 0 -1 1 0 "
      "-1 -1 0 0 -1 "
      "-1 -1 -1 1 0 ",
      5, 5);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_coladd_5by5_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 1 1 -1 "
      "-1 1 1 1 0 "
      "0 0 1 1 1 "
      "-1 1 0 -1 0 "
      "-1 1 0 0 -1 ",
      5, 5);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_coladd_5by5_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 0 -1 "
      "0 -1 0 1 -1 "
      "1 1 1 0 1 "
      "0 0 1 1 0 "
      "-1 0 -1 0 1 ",
      5, 5);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_coladd_5by5_4(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 -1 1 0 "
      "0 0 -1 1 -1 "
      "1 -1 0 0 1 "
      "0 1 1 -1 0 "
      "0 -1 0 1 0 ",
      5, 5);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_coladd_5by5_5(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 0 1 1 "
      "-1 1 1 0 0 "
      "0 0 0 -1 -1 "
      "1 -1 0 0 -1 "
      "-1 1 0 1 1 ",
      5, 5);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_coladd_5by5_6(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 1 0 0 "
      "-1 0 -1 0 0 "
      "0 -1 1 1 -1 "
      "-1 1 -1 -1 0 "
      "1 0 0 -1 0 ",
      5, 5);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_coladd_5by5_7(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 0 1 0 0 "
      "0 -1 0 1 0 "
      "0 1 -1 -1 0 "
      "-1 0 -1 -1 -1 "
      "-1 -1 0 1 -1 ",
      5, 5);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_coladd_5by5_8(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 1 -1 0 "
      "1 0 0 0 -1 "
      "-1 0 -1 1 1 "
      "0 0 0 -1 -1 "
      "1 1 0 -1 0 ",
      5, 5);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_coladd_5by5_9(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 0 -1 1 "
      "1 -1 -1 -1 1 "
      "0 0 -1 0 0 "
      "0 -1 -1 0 0 "
      "1 0 0 0 1 ",
      5, 5);
   runColumnTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A one by two case */
void test_network_rowadd_1by2_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 ",
      1, 2);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A one by two case */
void test_network_rowadd_1by2_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 ",
      1, 2);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A one by two case */
void test_network_rowadd_1by2_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 ",
      1, 2);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A two by three case */
void test_network_rowadd_2by3_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 +1 "
      "-1 +1 -1 ",
      2, 3);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A two by three case */
void test_network_rowadd_2by3_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 +1 "
      "-1 +1 +1 ",
      2, 3);
   runRowTestCase(&testCase, FALSE, TRUE);
   freeTestCase(&testCase);
}

/** @brief A two by three case */
void test_network_rowadd_2by3_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 +1 "
      "+1 0 +1 ",
      2, 3);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A two by three case */
void test_network_rowadd_2by3_4(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 "
      "+1 0  0 ",
      2, 3);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A two by three case */
void test_network_rowadd_2by3_5(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 "
      "+0 +1  0 ",
      2, 3);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A two by three case */
void test_network_rowadd_2by3_6(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 "
      "+0 -1  0 ",
      2, 3);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A two by three case */
void test_network_rowadd_2by3_7(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 "
      "+0 +1 +1 ",
      2, 3);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A two by three case */
void test_network_rowadd_2by3_8(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 "
      "+0 -1 +1 ",
      2, 3);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by six case */
void test_network_rowadd_3by6_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 0 0 0 "
      "0 0 +1 -1 0 0 "
      "-1 +1 -1 0 0 0 ",
      3, 6);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by six case */
void test_network_rowadd_3by6_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 0 0 0 0 "
      "0 0 +1 -1 0 0 "
      "-1 +1 -1 0 0 +1 ",
      3, 6);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by two case */
void test_network_rowadd_3by2_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 -1 "
      "-1 +1 "
      "+1 -1 ",
      3, 2);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by one case */
void test_network_rowadd_3by1_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 "
      "-1 "
      "+1 ",
      3, 1);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_rowadd_3by3_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 1 "
      "-1 1 0 "
      "0 1 -1 ",
      3, 3);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_rowadd_3by3_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 "
      "0 1 0 "
      "-1 -1 -1 ",
      3, 3);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_rowadd_3by3_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 0 "
      "-1 0 -1 "
      "-1 -1 0 ",
      3, 3);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_rowadd_3by3_4(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 1 "
      "0 -1 1 "
      "-1 -1 -1 ",
      3, 3);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_rowadd_3by3_5(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 -1 "
      "-1 -1 0 "
      "0 1 -1 ",
      3, 3);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_rowadd_3by3_6(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 "
      "1 1 0 "
      "-1 1 -1 ",
      3, 3);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_rowadd_3by3_7(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 0 "
      "0 1 -1 "
      "-1 1 -1 ",
      3, 3);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_rowadd_3by3_8(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 1 "
      "-1 1 0 "
      "-1 0 -1 ",
      3, 3);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by three case */
void test_network_rowadd_3by3_9(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 -1 "
      "-1 -1 -1 "
      "1 0 1 ",
      3, 3);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by six case */
void test_network_rowadd_3by6_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 0 0 0 -1 "
      "0 0 0 -1 -1 -1 "
      "-1 1 0 0 1 1 ",
      3, 3);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A three by six case */
void test_network_rowadd_3by6_4(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 -1 -1 0 0 "
      "0 0 0 1 -1 -1 "
      "0 -1 1 1 -1 0 ",
      3, 3);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 0 0 "
      "-1 -1 -1 -1 "
      "0 -1 -1 -1 "
      "1 0 0 -1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 0 -1 -1 "
      "1 1 -1 -1 "
      "0 0 -1 -1 "
      "0 -1 0 1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 1 0 1 "
      "-1 0 -1 0 "
      "0 0 1 1 "
      "0 -1 -1 1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_4(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 0 -1 0 "
      "0 1 0 0 "
      "-1 -1 1 1 "
      "0 -1 -1 1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_5(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 -1 1 "
      "1 0 0 -1 "
      "-1 0 -1 0 "
      "1 1 1 0 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_6(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 0 -1 1 "
      "0 1 -1 0 "
      "0 -1 1 -1 "
      "-1 -1 0 1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_7(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 0 -1 1 "
      "0 1 -1 0 "
      "0 -1 1 -1 "
      "-1 -1 0 1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_8(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 1 1 "
      "1 -1 -1 0 "
      "-1 1 1 1 "
      "1 0 0 -1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_9(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 0 0 "
      "-1 1 0 -1 "
      "0 1 1 -1 "
      "-1 0 1 0 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_10(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 1 -1 "
      "-1 0 0 -1 "
      "0 1 0 1 "
      "0 -1 0 0 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_11(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 1 0 "
      "-1 -1 0 1 "
      "-1 0 0 1 "
      "0 0 1 1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_12(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 0 -1 "
      "0 1 0 1 "
      "-1 0 1 0 "
      "0 0 -1 -1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_13(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 -1 1 "
      "-1 0 0 1 "
      "1 1 1 -1 "
      "1 0 -1 -1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_14(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 1 1 1 "
      "1 1 0 1 "
      "0 1 1 0 "
      "0 1 1 1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_15(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 -1 0 "
      "1 1 0 1 "
      "1 0 -1 0 "
      "0 1 0 1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_16(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 1 0 "
      "-1 0 0 -1 "
      "0 -1 0 1 "
      "-1 -1 0 0 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_17(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 -1 -1 0 "
      "-1 -1 0 0 "
      "0 -1 0 1 "
      "0 1 1 -1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A four by four case */
void test_network_rowadd_4by4_18(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 0 1 "
      "1 0 -1 0 "
      "-1 0 1 1 "
      "-1 1 0 -1 ",
      4, 4);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_rowadd_5by5_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 1 -1 1 "
      "0 0 1 -1 0 "
      "-1 0 0 1 0 "
      "0 0 1 0 1 "
      "0 -1 0 1 1 ",
      5, 5);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_rowadd_5by5_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "1 -1 0 -1 -1 "
      "0 -1 0 -1 -1 "
      "0 0 1 -1 0 "
      "0 -1 0 -1 -1 "
      "-1 1 1 0 0 ",
      5, 5);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_rowadd_5by5_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "-1 0 1 0 -1 "
      "0 1 1 0 -1 "
      "1 0 -1 0 1 "
      "-1 0 1 0 0 "
      "1 0 -1 0 1 ",
      5, 5);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_rowadd_5by5_4(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 -1 1 0 0 "
      "0 1 -1 1 0 "
      "0 -1 1 0 0 "
      "1 -1 1 0 0 "
      "0 1 0 1 1 ",
      5, 5);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A five by five case */
void test_network_rowadd_5by5_5(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "0 0 1 0 1 "
      "1 0 0 1 -1 "
      "1 -1 1 1 0 "
      "0 0 -1 0 -1 "
      "0 0 1 0 1 ",
      5, 5);
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief A eight by four case */
void test_network_rowadd_8by4(void)
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
   runRowTestCase(&testCase, FALSE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Updating a single rigid member */
void test_network_rowadd_singlerigid_1(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "+1 +1 0 "
      "0 -1 +1 "
      "+1 +1 0 ",
      4, 3);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Updating a single rigid member */
void test_network_rowadd_singlerigid_2(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "+1 +1 0 "
      "0 -1 +1 "
      "-1 -1 0 ",
      4, 3);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Updating a single rigid member */
void test_network_rowadd_singlerigid_3(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "+1 +1 0 "
      "0 -1 +1 "
      "-1 +1 0 ",
      4, 3);
   runRowTestCase(&testCase, FALSE, TRUE);
   freeTestCase(&testCase);
}

/** @brief Updating a single rigid member */
void test_network_rowadd_singlerigid_4(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "-1 -1 -1 "
      "0 +1 +1 "
      "+1 +1 +1 ",
      4, 3);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Updating a single rigid member */
void test_network_rowadd_singlerigid_5(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "-1 -1 -1 "
      "0 +1 +1 "
      "-1 -1 -1 ",
      4, 3);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Updating a single rigid member */
void test_network_rowadd_singlerigid_6(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 0 +1 "
      "-1 -1 -1 "
      "0 +1 +1 "
      "-1 +1 -1 ",
      4, 3);
   runRowTestCase(&testCase, FALSE, TRUE);
   freeTestCase(&testCase);
}

/** @brief Updating a single rigid member */
void test_network_rowadd_singlerigid_7(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 0 0 +1 "
      "+1 0 +1 0 0 "
      "0 -1 +1 +1 -1 "
      "0 0 0 -1 +1 "
      "+1 +1 0 0 0 ",
      5, 5);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}

/** @brief Updating a single rigid member */
void test_network_rowadd_singlerigid_8(void)
{
   DirectedTestCase testCase = stringToTestCase(
      "+1 +1 0 0 +1 "
      "+1 0 +1 0 0 "
      "0 -1 +1 +1 -1 "
      "0 0 0 -1 +1 "
      "+1 +1 0 0 0 "
      "+1 0 +1 +1 0 ",
      6, 5);
   runRowTestCase(&testCase, TRUE, FALSE);
   freeTestCase(&testCase);
}
/* TODO: test interleaved addition, test using random sampling + test erdos-renyi generated graphs */

/** @brief Computing the graph for a single rigid member */
void test_network_rowadd_singlerigid_graph(void)
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

void setUp(void) { SCIP_SUITE_SETUP(setup); }

void tearDown(void) { SCIP_SUITE_TEARDOWN(teardown); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_network_coladd_single_column
);
   RUN_TEST(test_network_coladd_doublecolumn_invalid_sign);
   RUN_TEST(test_network_coladd_doublecolumn_invalid_sign_2);
   RUN_TEST(test_network_coladd_splitseries_1);
   RUN_TEST(test_network_coladd_splitseries_2);
   RUN_TEST(test_network_coladd_splitseries_1r);
   RUN_TEST(test_network_coladd_splitseries_2r);
   RUN_TEST(test_network_coladd_splitseries_3);
   RUN_TEST(test_network_coladd_splitseries_4);
   RUN_TEST(test_network_coladd_splitseries_3r);
   RUN_TEST(test_network_coladd_splitseries_4r);
   RUN_TEST(test_network_coladd_splitseries_5);
   RUN_TEST(test_network_coladd_splitseries_6);
   RUN_TEST(test_network_coladd_splitseries_7);
   RUN_TEST(test_network_coladd_splitseries_8);
   RUN_TEST(test_network_coladd_splitseries_5r);
   RUN_TEST(test_network_coladd_splitseries_6r);
   RUN_TEST(test_network_coladd_splitseries_7r);
   RUN_TEST(test_network_coladd_splitseries_8r);
   RUN_TEST(test_network_coladd_splitseries_9);
   RUN_TEST(test_network_coladd_splitseries_10);
   RUN_TEST(test_network_coladd_splitseries_11);
   RUN_TEST(test_network_coladd_splitseries_12);
   RUN_TEST(test_network_coladd_splitseries_9r);
   RUN_TEST(test_network_coladd_splitseries_10r);
   RUN_TEST(test_network_coladd_splitseries_11r);
   RUN_TEST(test_network_coladd_splitseries_12r);
   RUN_TEST(test_network_coladd_splitseries_13);
   RUN_TEST(test_network_coladd_splitseries_14);
   RUN_TEST(test_network_coladd_splitseries_15);
   RUN_TEST(test_network_coladd_splitseries_16);
   RUN_TEST(test_network_coladd_splitseries_13r);
   RUN_TEST(test_network_coladd_splitseries_14r);
   RUN_TEST(test_network_coladd_splitseries_15r);
   RUN_TEST(test_network_coladd_splitseries_16r);
   RUN_TEST(test_network_coladd_parallelsimple_1);
   RUN_TEST(test_network_coladd_parallelsimple_2);
   RUN_TEST(test_network_coladd_parallelsimple_3);
   RUN_TEST(test_network_coladd_components_1);
   RUN_TEST(test_network_coladd_components_2);
   RUN_TEST(test_network_coladd_components_3);
   RUN_TEST(test_network_coladd_3by3_1);
   RUN_TEST(test_network_coladd_3by3_2);
   RUN_TEST(test_network_coladd_3by3_3);
   RUN_TEST(test_network_coladd_3by3_4);
   RUN_TEST(test_network_coladd_3by3_5);
   RUN_TEST(test_network_coladd_3by3_6);
   RUN_TEST(test_network_coladd_3by3_7);
   RUN_TEST(test_network_coladd_3by3_8);
   RUN_TEST(test_network_coladd_3by3_9);
   RUN_TEST(test_network_coladd_3by3_10);
   RUN_TEST(test_network_coladd_3by3_11);
   RUN_TEST(test_network_coladd_6by3_1);
   RUN_TEST(test_network_coladd_6by3_2);
   RUN_TEST(test_network_coladd_3by4_1);
   RUN_TEST(test_network_coladd_3by5_1);
   RUN_TEST(test_network_coladd_4by8_1);
   RUN_TEST(test_network_coladd_4by8_2);
   RUN_TEST(test_network_coladd_4by8_3);
   RUN_TEST(test_network_coladd_4by8_4);
   RUN_TEST(test_network_coladd_4by8_5);
   RUN_TEST(test_network_coladd_4by8_6);
   RUN_TEST(test_network_coladd_4by8_7);
   RUN_TEST(test_network_coladd_4by4_1);
   RUN_TEST(test_network_coladd_5by5_1);
   RUN_TEST(test_network_coladd_5by5_2);
   RUN_TEST(test_network_coladd_5by5_3);
   RUN_TEST(test_network_coladd_5by5_4);
   RUN_TEST(test_network_coladd_5by5_5);
   RUN_TEST(test_network_coladd_5by5_6);
   RUN_TEST(test_network_coladd_5by5_7);
   RUN_TEST(test_network_coladd_5by5_8);
   RUN_TEST(test_network_coladd_5by5_9);
   RUN_TEST(test_network_rowadd_1by2_1);
   RUN_TEST(test_network_rowadd_1by2_2);
   RUN_TEST(test_network_rowadd_1by2_3);
   RUN_TEST(test_network_rowadd_2by3_1);
   RUN_TEST(test_network_rowadd_2by3_2);
   RUN_TEST(test_network_rowadd_2by3_3);
   RUN_TEST(test_network_rowadd_2by3_4);
   RUN_TEST(test_network_rowadd_2by3_5);
   RUN_TEST(test_network_rowadd_2by3_6);
   RUN_TEST(test_network_rowadd_2by3_7);
   RUN_TEST(test_network_rowadd_2by3_8);
   RUN_TEST(test_network_rowadd_3by6_1);
   RUN_TEST(test_network_rowadd_3by6_2);
   RUN_TEST(test_network_rowadd_3by2_1);
   RUN_TEST(test_network_rowadd_3by1_1);
   RUN_TEST(test_network_rowadd_3by3_1);
   RUN_TEST(test_network_rowadd_3by3_2);
   RUN_TEST(test_network_rowadd_3by3_3);
   RUN_TEST(test_network_rowadd_3by3_4);
   RUN_TEST(test_network_rowadd_3by3_5);
   RUN_TEST(test_network_rowadd_3by3_6);
   RUN_TEST(test_network_rowadd_3by3_7);
   RUN_TEST(test_network_rowadd_3by3_8);
   RUN_TEST(test_network_rowadd_3by3_9);
   RUN_TEST(test_network_rowadd_3by6_3);
   RUN_TEST(test_network_rowadd_3by6_4);
   RUN_TEST(test_network_rowadd_4by4_1);
   RUN_TEST(test_network_rowadd_4by4_2);
   RUN_TEST(test_network_rowadd_4by4_3);
   RUN_TEST(test_network_rowadd_4by4_4);
   RUN_TEST(test_network_rowadd_4by4_5);
   RUN_TEST(test_network_rowadd_4by4_6);
   RUN_TEST(test_network_rowadd_4by4_7);
   RUN_TEST(test_network_rowadd_4by4_8);
   RUN_TEST(test_network_rowadd_4by4_9);
   RUN_TEST(test_network_rowadd_4by4_10);
   RUN_TEST(test_network_rowadd_4by4_11);
   RUN_TEST(test_network_rowadd_4by4_12);
   RUN_TEST(test_network_rowadd_4by4_13);
   RUN_TEST(test_network_rowadd_4by4_14);
   RUN_TEST(test_network_rowadd_4by4_15);
   RUN_TEST(test_network_rowadd_4by4_16);
   RUN_TEST(test_network_rowadd_4by4_17);
   RUN_TEST(test_network_rowadd_4by4_18);
   RUN_TEST(test_network_rowadd_5by5_1);
   RUN_TEST(test_network_rowadd_5by5_2);
   RUN_TEST(test_network_rowadd_5by5_3);
   RUN_TEST(test_network_rowadd_5by5_4);
   RUN_TEST(test_network_rowadd_5by5_5);
   RUN_TEST(test_network_rowadd_8by4);
   RUN_TEST(test_network_rowadd_singlerigid_1);
   RUN_TEST(test_network_rowadd_singlerigid_2);
   RUN_TEST(test_network_rowadd_singlerigid_3);
   RUN_TEST(test_network_rowadd_singlerigid_4);
   RUN_TEST(test_network_rowadd_singlerigid_5);
   RUN_TEST(test_network_rowadd_singlerigid_6);
   RUN_TEST(test_network_rowadd_singlerigid_7);
   RUN_TEST(test_network_rowadd_singlerigid_8);
   RUN_TEST(test_network_rowadd_singlerigid_graph);
   return UNITY_END();
}
