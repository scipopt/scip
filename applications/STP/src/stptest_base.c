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
SCIP_RETCODE stptest_all(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   assert(scip);

   SCIP_CALL( stptest_completegraph(scip) );
   SCIP_CALL( stptest_heur_extendPcMw(scip) );
   SCIP_CALL( stptest_heur_local(scip) );
   SCIP_CALL( stptest_pseudoDel(scip) );
   SCIP_CALL( stptest_reduce_ext(scip) );
   SCIP_CALL( stptest_dheap(scip) );
   SCIP_CALL( stptest_reduce_sdpcmw(scip) );

   printf("all unit tests passed! \n");

   return SCIP_OKAY;
}
