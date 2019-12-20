/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   stptest_base.c
 * @brief  tests for Steiner tree problem methods
 * @author Daniel Rehfeldt
 *
 * This file implements tests for Steiner problems.
 *
 * A list of all interface methods can be found in stptest.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scip/scip.h"
#include "stptest.h"


/** tests all */
SCIP_RETCODE testAll(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   assert(scip);

   SCIP_CALL( completegraph_test(scip) );
   SCIP_CALL( heur_extendPcMwTest(scip) );
   SCIP_CALL( heur_localTest(scip) );
   SCIP_CALL( pseudoDel_test(scip) );
   SCIP_CALL( reduce_extTest(scip) );
   SCIP_CALL( dheap_Test(scip) );
   SCIP_CALL( reduce_sdPcMwTest(scip) );

   printf("all unit tests passed! \n");

   return SCIP_OKAY;
}
