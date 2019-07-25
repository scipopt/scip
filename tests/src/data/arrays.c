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
static SCIP_Real myrealarray[] = {
   5.0,
   23.3,
   14.5
};

/* creates scip and problem */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

//   /* create a problem */
//   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* create arrays */
   SCIP_CALL( SCIPcreateRealarray(scip, &realarray) );
   SCIP_CALL( SCIPcreateIntarray(scip, &intarray) );
   SCIP_CALL( SCIPcreateBoolarray(scip, &boolarray) );
   SCIP_CALL( SCIPcreatePtrarray(scip, &ptrarray) );
}

static
void teardown(void)
{
   /* real array */
   SCIP_CALL( SCIPextendRealarray(scip, realarray, 0, 7) );

   SCIP_CALL( SCIPsetRealarrayVal(scip, realarray, 0, 1.2) );

   SCIP_CALL( SCIPincRealarrayVal(scip, realarray, 0, 1.3) );

   SCIP_CALL( SCIPclearRealarray(scip, realarray) );


   /* integer array */
   SCIP_CALL( SCIPextendIntarray(scip, intarray, 0, 7) );

   SCIP_CALL( SCIPsetIntarrayVal(scip, intarray, 0, 1) );

   SCIP_CALL( SCIPincIntarrayVal(scip, intarray, 0, 1) );

   SCIP_CALL( SCIPclearIntarray(scip, intarray) );


   /* boolean array */
   SCIP_CALL( SCIPextendBoolarray(scip, boolarray, 0, 7) );

   SCIP_CALL( SCIPsetBoolarrayVal(scip, boolarray, 0, TRUE) );

   SCIP_CALL( SCIPclearBoolarray(scip, boolarray) );


   /* pointer array */
   SCIP_CALL( SCIPextendPtrarray(scip, ptrarray, 0, 7) );


   SCIP_CALL( SCIPclearPtrarray(scip, ptrarray) );

   /* free arrays and scip */
   SCIP_CALL( SCIPfreeRealarray(scip, &realarray) );
   SCIP_CALL( SCIPfreeIntarray(scip, &intarray) );
   SCIP_CALL( SCIPfreeBoolarray(scip, &boolarray) );
   SCIP_CALL( SCIPfreePtrarray(scip, &ptrarray) );

   SCIP_CALL( SCIPfree(&scip) );
}

TestSuite(arrays, .init = setup, .fini = teardown);

Test(arrays, setup_and_teardown, .description = "test that setup and teardown work correctly")
{

}

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
