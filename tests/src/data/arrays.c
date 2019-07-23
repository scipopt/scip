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

/* creates scip and problem */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

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
   assert(0 == SCIPgetRealarrayMinIdx(scip, realarray));
   assert(7 == SCIPgetRealarrayMaxIdx(scip, realarray));

   SCIP_CALL( SCIPsetRealarrayVal(scip, realarray, 0, 1.2) );
   assert(1.2 == SCIPgetRealarrayVal(scip, realarray, 0));

   SCIP_CALL( SCIPincRealarrayVal(scip, realarray, 0, 1.3) );
   assert(2.5 == SCIPgetRealarrayVal(scip, realarray, 0));

   SCIP_CALL( SCIPclearRealarray(scip, realarray) );
   SCIP_CALL( SCIPfreeRealarray(scip, &realarray) );


   /* integer array */
   SCIP_CALL( SCIPextendIntarray(scip, intarray, 0, 7) );
   assert(0 == SCIPgetIntarrayMinIdx(scip, intarray));
   assert(7 == SCIPgetIntarrayMaxIdx(scip, intarray));

   SCIP_CALL( SCIPsetIntarrayVal(scip, intarray, 0, 1) );
   assert(1.0 == SCIPgetIntarrayVal(scip, intarray, 0));

   SCIP_CALL( SCIPincIntarrayVal(scip, intarray, 0, 1) );
   assert(2.0 == SCIPgetIntarrayVal(scip, intarray, 0));

   SCIP_CALL( SCIPclearIntarray(scip, intarray) );
   SCIP_CALL( SCIPfreeIntarray(scip, &intarray) );


   /* boolean array */
   SCIP_CALL( SCIPextendBoolarray(scip, boolarray, 0, 7) );
   assert(0 == SCIPgetBoolarrayMinIdx(scip, boolarray));
   assert(7 == SCIPgetBoolarrayMaxIdx(scip, boolarray));

   SCIP_CALL( SCIPsetBoolarrayVal(scip, boolarray, 0, TRUE) );
   assert(SCIPgetBoolarrayVal(scip, boolarray, 0));

   SCIP_CALL( SCIPclearBoolarray(scip, boolarray) );
   SCIP_CALL( SCIPfreeBoolarray(scip, &boolarray) );


   /* pointer array */
   SCIP_CALL( SCIPextendPtrarray(scip, ptrarray, 0, 7) );
   assert(0 == SCIPgetPtrarrayMinIdx(scip, ptrarray));
   assert(7 == SCIPgetPtrarrayMaxIdx(scip, ptrarray));

   SCIP_CALL( SCIPsetPtrarrayVal(scip, ptrarray, 0, (void*) realarray) );
   assert((void*) realarray == SCIPgetPtrarrayVal(scip, ptrarray, 0));

   SCIP_CALL( SCIPclearPtrarray(scip, ptrarray) );
   SCIP_CALL( SCIPfreePtrarray(scip, &ptrarray) );


   /* free arrays and scip */
   SCIP_CALL( SCIPfreeRealarray(scip, &realarray) );
   SCIP_CALL( SCIPfreeIntarray(scip, &intarray) );
   SCIP_CALL( SCIPfreeBoolarray(scip, &boolarray) );
   SCIP_CALL( SCIPfreePtrarray(scip, &ptrarray) );

   SCIP_CALL( SCIPfree(&scip) );
}
