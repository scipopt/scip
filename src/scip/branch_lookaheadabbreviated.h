/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_lookaheadabbreviated.h
 * @ingroup BRANCHINGRULES
 * @brief  LookaheadAbbreviated branching rule
 * @author Christoph Schubert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_LOOKAHEADABBREVIATED_H__
#define __SCIP_BRANCH_LOOKAHEADABBREVIATED_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the LookaheadAbbreviated branching rule and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeBranchruleLookaheadAbbreviated(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
