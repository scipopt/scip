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
extern SCIP_RETCODE    stptest_all(SCIP*);

/* stptest_reduce.c
 */
extern SCIP_RETCODE    stptest_reduce_ext(SCIP*);
extern SCIP_RETCODE    stptest_reduce_sdpcmw(SCIP*);

/* stptest_misc.c
 */
extern SCIP_RETCODE    stptest_pseudoAncestors(SCIP*);
extern SCIP_RETCODE    stptest_pseudoDel(SCIP*);
extern SCIP_RETCODE    stptest_dheap(SCIP*);
extern SCIP_RETCODE    stptest_completegraph(SCIP*);

/* stptest_heur.c
 */
extern SCIP_RETCODE    stptest_heur_extendPcMw(SCIP*);
extern SCIP_RETCODE    stptest_heur_local(SCIP*);


#endif /* APPLICATIONS_STP_SRC_STPTEST_H_ */
