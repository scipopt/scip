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

/**@file   bandit_epsgreedy.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for epsilon greedy bandit selection
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BANDIT_EPSGREEDY_H__
#define __SCIP_BANDIT_EPSGREEDY_H__


#include "blockmemshell/memory.h"
#include "scip/def.h"
#include "scip/type_bandit.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the epsilon greedy bandit algorithm includes it in SCIP */
SCIP_RETCODE SCIPincludeBanditvtableEpsgreedy(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** callback to free bandit specific data structures */
SCIP_DECL_BANDITFREE(SCIPbanditFreeEpsgreedy);

/** selection callback for bandit algorithm */
SCIP_DECL_BANDITSELECT(SCIPbanditSelectEpsgreedy);

/** update callback for bandit algorithm */
SCIP_DECL_BANDITUPDATE(SCIPbanditUpdateEpsgreedy);

/** reset callback for bandit algorithm */
SCIP_DECL_BANDITRESET(SCIPbanditResetEpsgreedy);


/** internal method to create and reset epsilon greedy bandit algorithm */
SCIP_RETCODE SCIPbanditCreateEpsgreedy(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_BANDITVTABLE*    vtable,             /**< virtual function table with epsilon greedy callbacks */
   SCIP_BANDIT**         epsgreedy,          /**< pointer to store the epsilon greedy bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   SCIP_Real             eps,                /**< parameter to increase probability for exploration between all actions */
   SCIP_Bool             preferrecent,       /**< should the weights be updated in an exponentially decaying way? */
   SCIP_Real             decayfactor,        /**< the factor to reduce the weight of older observations if exponential decay is enabled */
   int                   avglim,             /**< nonnegative limit on observation number before the exponential decay starts,
                                               *  only relevant if exponential decay is enabled
                                               */
   int                   nactions,           /**< the number of possible actions */
   unsigned int          initseed            /**< initial random seed */
   );

#ifdef __cplusplus
}
#endif

#endif
