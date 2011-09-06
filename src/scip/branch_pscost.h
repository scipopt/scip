/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_pscost.h
 * @brief  pseudo costs branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_PSCOST_H__
#define __SCIP_BRANCH_PSCOST_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the pseudo cost branching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchrulePscost(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** selects a branching variable, due to pseudo cost, from the given candidate array and returns this variable together
 *  with a branching point */
extern
SCIP_RETCODE SCIPselectBranchVarPscost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            branchcands,        /**< branching candidates */
   SCIP_Real*            branchcandssol,     /**< solution value for the branching candidates */
   SCIP_Real*            branchcandsscore,   /**< array of candidate scores */
   int                   nbranchcands,       /**< number of branching candidates */
   SCIP_VAR**            var,                /**< pointer to store the variable to branch on, or NULL if none */
   SCIP_Real*            brpoint             /**< pointer to store the branching point for the branching variable, will be fractional for a discrete variable */
   );
   
#ifdef __cplusplus
}
#endif

#endif
