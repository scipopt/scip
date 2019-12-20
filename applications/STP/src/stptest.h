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

/**@file   stptest.h
 * @brief  includes various testing methods for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_STPTEST_H_
#define APPLICATIONS_STP_SRC_STPTEST_H_

#include "scip/scip.h"

/* reduce_base.c
 */
extern SCIP_RETCODE    testAll(SCIP*);

/* stptest_reduce.c
 */
extern SCIP_RETCODE    reduce_extTest(SCIP*);
extern SCIP_RETCODE    reduce_sdPcMwTest(SCIP*);

/* stptest_misc.c
 */
extern SCIP_RETCODE    pseudoAncestors_test(SCIP*);
extern SCIP_RETCODE    pseudoDel_test(SCIP*);
extern SCIP_RETCODE    dheap_Test(SCIP*);
extern SCIP_RETCODE    completegraph_test(SCIP*);

/* stptest_heur.c
 */
extern SCIP_RETCODE    heur_extendPcMwTest(SCIP*);
extern SCIP_RETCODE    heur_localTest(SCIP*);



#endif /* APPLICATIONS_STP_SRC_STPTEST_H_ */
