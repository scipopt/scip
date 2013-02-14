/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_allfullstrong.h
 * @ingroup BRANCHINGRULES
 * @brief  all variables full strong LP branching rule
 * @author Tobias Achterberg
 *
 * The all variables full strong branching rule applies strong branching to every non-fixed variable
 * at the current node of the branch-and-bound search. The rule selects the candidate
 * which will cause the highest gain of the dual bound in the created sub-tree among all branching variables.
 *
 * For calculating the gain, a look-ahead is performed by solving the child node LPs which will result
 * from branching on a variable.
 *
 * For a more mathematical description and a comparison between the strong branching rule and other branching rules
 * in SCIP, we refer to
 *
 * @par
 * Tobias Achterberg@n
 * Constraint Integer Programming@n
 * PhD Thesis, Technische Universität Berlin, 2007@n
 *

 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_ALLFULLSTRONG_H__
#define __SCIP_BRANCH_ALLFULLSTRONG_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the all variables full strong LP branching rule and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeBranchruleAllfullstrong(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
