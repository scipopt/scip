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

/**@file   arrays.c
 * @brief  unittest for arrays in scip_datastructures.c
 * @author Merlin Viernickel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"
#include "scip/scip_datastructures.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_REALARRAY* realarray;
static SCIP_INTARRAY* intarray;
static SCIP_BOOLARRAY* boolarray;
static SCIP_PTRARRAY* ptrarray;

#define arraylen 3
static SCIP_Real myrealarray[] =
{
   5.0,
   23.3,
   14.5
};
static int myintarray[] =
{
   8,
   17,
   12
};
static SCIP_Bool myboolarray[] =
{
   TRUE,
   FALSE,
   TRUE
};
static SCIP_Real* myptrarray[] =
{
   &myrealarray[0],
   &myrealarray[1],
   &myrealarray[2]
};

/* creates scip and arrays */
static
void setup(void)
{
   /* create scip and arrays */
   SCIP_CALL( SCIPcreate(&scip) );

   SCIP_CALL( SCIPcreateRealarray(scip, &realarray) );
   SCIP_CALL( SCIPcreateIntarray(scip, &intarray) );
   SCIP_CALL( SCIPcreateBoolarray(scip, &boolarray) );
   SCIP_CALL( SCIPcreatePtrarray(scip, &ptrarray) );
}

/* frees scip and arrays */
static
void teardown(void)
{
   /* free scip and arrays */
   SCIP_CALL( SCIPfreePtrarray(scip, &ptrarray) );
   SCIP_CALL( SCIPfreeBoolarray(scip, &boolarray) );
   SCIP_CALL( SCIPfreeIntarray(scip, &intarray) );
   SCIP_CALL( SCIPfreeRealarray(scip, &realarray) );

   SCIP_CALL( SCIPfree(&scip) );
}

void setUp(void)
{
   setup();
}

void tearDown(void)
{
   teardown();
}

/** @brief test that setup and teardown work correctly */
void test_arrays_setup_and_teardown(void)
{
}

/* dynamic real array tests */

/** @brief test that the dynamic real array stores entries correctly. */
void test_arrays_insertion_real(void)
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPsetRealarrayVal(scip, realarray, i, myrealarray[i]) );

   for( i = 0; i < arraylen; i++ )
      TEST_ASSERT_EQUAL(myrealarray[i], SCIPgetRealarrayVal(scip, realarray, i));
}

/** @brief test that the dynamic real array increments entries correctly */
void test_arrays_increment_real(void)
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPincRealarrayVal(scip, realarray, i, myrealarray[i]) );

   for( i = 0; i < arraylen; i++ )
      TEST_ASSERT_EQUAL(myrealarray[i], SCIPgetRealarrayVal(scip, realarray, i));
}

/** @brief test that the dynamic real array clears entries correctly */
void test_arrays_clear_real(void)
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPsetRealarrayVal(scip, realarray, i, myrealarray[i]) );

   SCIP_CALL( SCIPclearRealarray(scip, realarray) );

   for( i = 0; i < arraylen; i++ )
      TEST_ASSERT_EQUAL(0.0, SCIPgetRealarrayVal(scip, realarray, i));
}

/** @brief test that the dynamic real array stores max and min indices correctly */
void test_arrays_indices_real(void)
{
   int i;

   SCIP_CALL( SCIPextendRealarray(scip, realarray, 0, arraylen - 1) );

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPsetRealarrayVal(scip, realarray, i, myrealarray[i]) );

   TEST_ASSERT_EQUAL(0, SCIPgetRealarrayMinIdx(scip, realarray));
   TEST_ASSERT_EQUAL(arraylen - 1, SCIPgetRealarrayMaxIdx(scip, realarray));
}

/* dynamic integer array tests */

/** @brief test that the dynamic integer array stores entries correctly. */
void test_arrays_insertion_int(void)
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPsetIntarrayVal(scip, intarray, i, myintarray[i]) );

   for( i = 0; i < arraylen; i++ )
      TEST_ASSERT_EQUAL(myintarray[i], SCIPgetIntarrayVal(scip, intarray, i));
}

/** @brief test that the dynamic integer array increments entries correctly */
void test_arrays_increment_int(void)
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPincIntarrayVal(scip, intarray, i, myintarray[i]) );

   for( i = 0; i < arraylen; i++ )
      TEST_ASSERT_EQUAL(myintarray[i], SCIPgetIntarrayVal(scip, intarray, i));
}

/** @brief test that the dynamic integer array clears entries correctly */
void test_arrays_clear_int(void)
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPsetIntarrayVal(scip, intarray, i, myintarray[i]) );

   SCIP_CALL( SCIPclearIntarray(scip, intarray) );

   for( i = 0; i < arraylen; i++ )
      TEST_ASSERT_EQUAL(0, SCIPgetIntarrayVal(scip, intarray, i));
}

/** @brief test that the dynamic integer array stores max and min indices correctly */
void test_arrays_indices_int(void)
{
   int i;

   SCIP_CALL( SCIPextendIntarray(scip, intarray, 0, arraylen - 1) );

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPsetIntarrayVal(scip, intarray, i, myintarray[i]) );

   TEST_ASSERT_EQUAL(0, SCIPgetIntarrayMinIdx(scip, intarray));
   TEST_ASSERT_EQUAL(arraylen - 1, SCIPgetIntarrayMaxIdx(scip, intarray));
}

/* dynamic boolean array tests */

/** @brief test that the dynamic boolean array stores entries correctly. */
void test_arrays_insertion_bool(void)
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPsetBoolarrayVal(scip, boolarray, i, myboolarray[i]) );

   for( i = 0; i < arraylen; i++ )
      TEST_ASSERT_EQUAL(myboolarray[i], SCIPgetBoolarrayVal(scip, boolarray, i));
}

/** @brief test that the dynamic boolean array clears entries correctly */
void test_arrays_clear_bool(void)
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPsetBoolarrayVal(scip, boolarray, i, myboolarray[i]) );

   SCIP_CALL( SCIPclearBoolarray(scip, boolarray) );

   for( i = 0; i < arraylen; i++ )
      TEST_ASSERT_EQUAL(FALSE, SCIPgetBoolarrayVal(scip, boolarray, i));
}

/** @brief test that the dynamic boolean array stores max and min indices correctly */
void test_arrays_indices_bool(void)
{
   int i;

   SCIP_CALL( SCIPextendBoolarray(scip, boolarray, 0, arraylen - 1) );

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPsetBoolarrayVal(scip, boolarray, i, myboolarray[i]) );

   TEST_ASSERT_EQUAL(0, SCIPgetBoolarrayMinIdx(scip, boolarray));
   TEST_ASSERT_EQUAL(arraylen - 1, SCIPgetBoolarrayMaxIdx(scip, boolarray));
}

/* dynamic pointer array tests */

/** @brief test that the dynamic pointer array stores entries correctly. */
void test_arrays_insertion_ptr(void)
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPsetPtrarrayVal(scip, ptrarray, i, myptrarray[i]) );

   for( i = 0; i < arraylen; i++ )
      TEST_ASSERT_EQUAL(myptrarray[i], SCIPgetPtrarrayVal(scip, ptrarray, i));
}

/** @brief test that the dynamic pointer array clears entries correctly */
void test_arrays_clear_ptr(void)
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPsetPtrarrayVal(scip, ptrarray, i, myptrarray[i]) );

   SCIP_CALL( SCIPclearPtrarray(scip, ptrarray) );

   for( i = 0; i < arraylen; i++ )
      TEST_ASSERT_EQUAL(NULL, SCIPgetPtrarrayVal(scip, ptrarray, i));
}

/** @brief test that the dynamic pointer array stores max and min indices correctly */
void test_arrays_indices_ptr(void)
{
   int i;

   SCIP_CALL( SCIPextendPtrarray(scip, ptrarray, 0, arraylen - 1) );

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPsetPtrarrayVal(scip, ptrarray, i, myptrarray[i]) );

   TEST_ASSERT_EQUAL(0, SCIPgetPtrarrayMinIdx(scip, ptrarray));
   TEST_ASSERT_EQUAL(arraylen - 1, SCIPgetPtrarrayMaxIdx(scip, ptrarray));
}

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_arrays_setup_and_teardown);
   RUN_TEST(test_arrays_insertion_real);
   RUN_TEST(test_arrays_increment_real);
   RUN_TEST(test_arrays_clear_real);
   RUN_TEST(test_arrays_indices_real);
   RUN_TEST(test_arrays_insertion_int);
   RUN_TEST(test_arrays_increment_int);
   RUN_TEST(test_arrays_clear_int);
   RUN_TEST(test_arrays_indices_int);
   RUN_TEST(test_arrays_insertion_bool);
   RUN_TEST(test_arrays_clear_bool);
   RUN_TEST(test_arrays_indices_bool);
   RUN_TEST(test_arrays_insertion_ptr);
   RUN_TEST(test_arrays_clear_ptr);
   RUN_TEST(test_arrays_indices_ptr);
   return UNITY_END();
}
