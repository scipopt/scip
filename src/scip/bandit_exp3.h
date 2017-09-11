/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   bandit_exp3.h
 * @ingroup BanditAlgorithms
 * @brief  public methods for Exp.3 bandit algorithm
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BANDIT_EXP3_H__
#define __SCIP_BANDIT_EXP3_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicBanditMethods
 *
 * ## Exp.3
 *
 * Exp.3 is a randomized selection method for the multi-armed bandit problem
 *
 * Exp3 maintains a probability distribution
 * according to which an action is drawn
 * in every iteration.
 * The probability distribution is a mixture between
 * a uniform distribution and a softmax distribution
 * based on the cumulative rewards of the actions.
 * The weight of the uniform distribution in the mixture
 * is controlled by the parameter \f$ \gamma \f$, ie.,
 * setting \f$ \gamma = 1\f$ uses a uniform distribution
 * in every selection step.
 * The cumulative reward for the actions can be
 * fine-tuned by adding a general bias for all actions.
 * The bias is given by the parameter \f$ \beta \f$.
 *
 * @{
 */

/** include virtual function table for Exp.3 bandit algorithms */
EXTERN
SCIP_RETCODE SCIPincludeBanditvtableExp3(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** create Exp.3 bandit algorithm */
EXTERN
SCIP_RETCODE SCIPcreateBanditExp3(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         exp3,               /**< pointer to store bandit algorithm */
   int                   nactions,           /**< the number of actions for this bandit algorithm */
   SCIP_Real             gammaparam,         /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution */
   SCIP_Real             beta                /**< gain offset between 0 and 1 at every observation */
   );

/** set gamma parameter of Exp.3 bandit algorithm to increase weight of uniform distribution */
EXTERN
void SCIPsetGammaExp3(
   SCIP_BANDIT*          exp3,               /**< bandit algorithm */
   SCIP_Real             gammaparam          /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution */
   );

/** set beta parameter of Exp.3 bandit algorithm to increase gain offset for actions that were not played */
EXTERN
void SCIPsetBetaExp3(
   SCIP_BANDIT*          exp3,               /**< bandit algorithm */
   SCIP_Real             beta                /**< gain offset between 0 and 1 at every observation */
   );

/** returns probability to play an action */
EXTERN
SCIP_Real SCIPgetProbabilityExp3(
   SCIP_BANDIT*          exp3,               /**< bandit algorithm */
   int                   action              /**< index of the requested action */
   );

/** @}*/

/*
 * callback methods for Exp.3 bandit algorithm to create a virtual function table
 */

/** callback to free bandit specific data structures */
extern
SCIP_DECL_BANDITFREE(SCIPbanditFreeExp3);

/** selection callback for bandit selector */
extern
SCIP_DECL_BANDITSELECT(SCIPbanditSelectExp3);

/** update callback for bandit algorithm */
extern
SCIP_DECL_BANDITUPDATE(SCIPbanditUpdateExp3);

/** reset callback for bandit algorithm */
extern
SCIP_DECL_BANDITRESET(SCIPbanditResetExp3);

/** direct bandit creation method for the core where no SCIP pointer is available */
extern
SCIP_RETCODE SCIPbanditCreateExp3(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_BANDITVTABLE*    vtable,             /**< virtual function table for callback functions of Exp.3 */
   SCIP_BANDIT**         exp3,               /**< pointer to store bandit algorithm */
   int                   nactions,           /**< the number of actions for this bandit algorithm */
   SCIP_Real             gammaparam,         /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution */
   SCIP_Real             beta,               /**< gain offset between 0 and 1 at every observation */
   unsigned int          initseed            /**< initial random seed */
   );

#ifdef __cplusplus
}
#endif

#endif
