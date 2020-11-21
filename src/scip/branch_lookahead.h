/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_lookahead.h
 * @ingroup BRANCHINGRULES
 * @brief  lookahead LP branching rule
 * @author Christoph Schubert
 * @author Gerald Gamrath
 *
 * The (multi-level) lookahead branching rule applies strong branching to every fractional value of the LP solution
 * at the current node of the branch-and-bound tree, as well as recursivly to every temporary child problem created by this
 * strong branching. The rule selects the candidate with the best proven dual bound.
 *
 * The branching rule was motivated by the following technical report:
 *
 * @par
 * Wasu Glankwamdee and Jeff Linderoth@n
 * Lookahead Branching for Mixed Integer Programming@n
 * Technical Report 06T-004, Department of Industrial and Systems Engineering, Lehigh University.
 *
 * For a more mathematical description and a comparison between lookahead branching and other branching rules
 * in SCIP, we refer to
 *
 * @par
 * Christoph Schubert@n
 * Multi-Level Lookahead Branching@n
 * Master Thesis, Technische Universität Berlin, 2017@n
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_LOOKAHEAD_H__
#define __SCIP_BRANCH_LOOKAHEAD_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the lookahead branching rule and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBranchruleLookahead(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
