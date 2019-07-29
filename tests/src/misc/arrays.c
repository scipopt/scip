/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  this file is part of the program and library             */
/*         scip --- solving constraint integer programs                      */
/*                                                                           */
/*    copyright (c) 2002-2019 konrad-zuse-zentrum                            */
/*                            fuer informationstechnik berlin                */
/*                                                                           */
/*  scip is distributed under the terms of the zib academic license.         */
/*                                                                           */
/*  you should have received a copy of the zib academic license              */
/*  along with scip; see the file copying. if not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   arrays.c
 * @brief  unittest for arrays in scip_datastructures.c
 * @author merlin viernickel
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

static const int arraylen = 3;
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
static void* myptrarray[] =
{
   NULL,
   NULL,
   NULL
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

TestSuite(arrays, .init = setup, .fini = teardown);

Test(arrays, setup_and_teardown, .description = "test that setup and teardown work correctly")
{

}

/* dynamic real array tests */

Test(arrays, test_insertion_real, .description = "test that the dynamic real array stores entries correctly.")
{
   int i;

   for (i = 0; i < arraylen; i++)
   {
      SCIP_CALL( SCIPsetRealarrayVal(scip, realarray, i, myrealarray[i]) );
   }
   for (i = 0; i < arraylen; i++)
   {
      cr_assert_eq(myrealarray[i], SCIPgetRealarrayVal(scip, realarray, i));
   }
}

Test(arrays, test_increment_real, .description = "test that the dynamic real array increments entries correctly")
{
   int i;

   for (i = 0; i < arraylen; i++)
   {
      SCIP_CALL( SCIPincRealarrayVal(scip, realarray, i, myrealarray[i]) );
   }
   for (i = 0; i < arraylen; i++)
   {
      cr_assert_eq(myrealarray[i], SCIPgetRealarrayVal(scip, realarray, i));
   }
}

Test(arrays, test_clear_real, .description = "test that the dynamic real array clears entries correctly")
{
   int i;

   for (i = 0; i < arraylen; i++)
   {
      SCIP_CALL( SCIPsetRealarrayVal(scip, realarray, i, myrealarray[i]) );
   }
   SCIP_CALL( SCIPclearRealarray(scip, realarray) );
   for (i = 0; i < arraylen; i++)
   {
      cr_assert_eq(0.0, SCIPgetRealarrayVal(scip, realarray, i));
   }
}


/* dynamic integer array tests */

Test(arrays, test_insertion_int, .description = "test that the dynamic integer array stores entries correctly.")
{
   int i;

   for (i = 0; i < arraylen; i++)
   {
      SCIP_CALL( SCIPsetIntarrayVal(scip, intarray, i, myintarray[i]) );
   }
   for (i = 0; i < arraylen; i++)
   {
      cr_assert_eq(myintarray[i], SCIPgetIntarrayVal(scip, intarray, i));
   }
}

Test(arrays, test_increment_int, .description = "test that the dynamic integer array increments entries correctly")
{
   int i;

   for (i = 0; i < arraylen; i++)
   {
      SCIP_CALL( SCIPincIntarrayVal(scip, intarray, i, myintarray[i]) );
   }
   for (i = 0; i < arraylen; i++)
   {
      cr_assert_eq(myintarray[i], SCIPgetIntarrayVal(scip, intarray, i));
   }
}

Test(arrays, test_clear_int, .description = "test that the dynamic integer array clears entries correctly")
{
   int i;

   for (i = 0; i < arraylen; i++)
   {
      SCIP_CALL( SCIPsetIntarrayVal(scip, intarray, i, myintarray[i]) );
   }
   SCIP_CALL( SCIPclearIntarray(scip, intarray) );
   for (i = 0; i < arraylen; i++)
   {
      cr_assert_eq(0, SCIPgetIntarrayVal(scip, intarray, i));
   }
}

/* dynamic boolean array tests */

Test(arrays, test_insertion_bool, .description = "test that the dynamic boolean array stores entries correctly.")
{
   int i;

   for (i = 0; i < arraylen; i++)
   {
      SCIP_CALL( SCIPsetBoolarrayVal(scip, boolarray, i, myboolarray[i]) );
   }
   for (i = 0; i < arraylen; i++)
   {
      cr_assert_eq(myboolarray[i], SCIPgetBoolarrayVal(scip, boolarray, i));
   }
}

Test(arrays, test_clear_bool, .description = "test that the dynamic boolean array clears entries correctly")
{
   int i;

   for (i = 0; i < arraylen; i++)
   {
      SCIP_CALL( SCIPsetBoolarrayVal(scip, boolarray, i, myboolarray[i]) );
   }
   SCIP_CALL( SCIPclearBoolarray(scip, boolarray) );
   for (i = 0; i < arraylen; i++)
   {
      cr_assert_eq(FALSE, SCIPgetBoolarrayVal(scip, boolarray, i));
   }
}

/* dynamic boolean array tests */

Test(arrays, test_insertion_ptr, .description = "test that the dynamic pointer array stores entries correctly.")
{
   int i;

   for (i = 0; i < arraylen; i++)
   {
      SCIP_CALL( SCIPsetPtrarrayVal(scip, ptrarray, i, myptrarray[i]) );
   }
   for (i = 0; i < arraylen; i++)
   {
      cr_assert_eq(myptrarray[i], SCIPgetPtrarrayVal(scip, ptrarray, i));
   }
}

Test(arrays, test_clear_ptr, .description = "test that the dynamic pointer array clears entries correctly")
{
   int i;

   for (i = 0; i < arraylen; i++)
   {
      SCIP_CALL( SCIPsetPtrarrayVal(scip, ptrarray, i, myptrarray[i]) );
   }
   SCIP_CALL( SCIPclearPtrarray(scip, ptrarray) );
   for (i = 0; i < arraylen; i++)
   {
      cr_assert_eq(NULL, SCIPgetPtrarrayVal(scip, ptrarray, i));
   }
}
