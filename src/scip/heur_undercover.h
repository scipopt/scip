/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_undercover.h,v 1.1 2009/11/16 02:47:29 bzfgleix Exp $"

/**@file   heur_undercover.h
 * @brief  undercover primal heuristic for MIQCPs
 * @author Ambros Gleixner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_UNDERCOVER_H__
#define __SCIP_HEUR_UNDERCOVER_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the undercover primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurUndercover(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** main procedure of the undercover heuristic */
SCIP_RETCODE SCIPapplyUndercover(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< result data structure                                          */
   SCIP_Real             timelimit,          /**< time limit                                                     */        
   SCIP_Real             memorylimit,        /**< memory limit                                                   */
   char                  ppcstrat,           /**< strategy for finding a ppcsolution                             */
   SCIP_Real             ppcobjquot,         /**< additional penalty factor for fixing continuous variables      */
   SCIP_Real             domred,             /**< reduce domain of selected variables by this factor around LP value */
   SCIP_Bool             locksrounding,      /**< shall LP values for integer vars be rounded according to locks? */
   SCIP_Bool             onlyconvexify,      /**< should variable fixing only convexify the problem?             */
   SCIP_Real             minimprove,         /**< factor by which heuristic should at least improve the incumbent*/
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem                    */
   SCIP_Bool             postnlp             /**< shall NLP heuristic be called when a feas. solution was found? */
   );

#ifdef __cplusplus
}
#endif

#endif
