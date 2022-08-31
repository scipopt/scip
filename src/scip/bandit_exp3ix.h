/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   bandit_exp3ix.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for Exp.3-IX bandit algorithm
 * @author Antonia Chmiela
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BANDIT_EXP3IX_H__
#define __SCIP_BANDIT_EXP3IX_H__


#include "blockmemshell/memory.h"
#include "scip/def.h"
#include "scip/type_bandit.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** include virtual function table for Exp.3-IX bandit algorithms */
SCIP_RETCODE SCIPincludeBanditvtableExp3IX(
   SCIP*                 scip                /**< SCIP data structure */
   );

/*
 * callback methods for Exp.3-IX bandit algorithm to create a virtual function table
 */

/** callback to free bandit specific data structures */
SCIP_DECL_BANDITFREE(SCIPbanditFreeExp3IX);

/** selection callback for bandit selector */
SCIP_DECL_BANDITSELECT(SCIPbanditSelectExp3IX);

/** update callback for bandit algorithm */
SCIP_DECL_BANDITUPDATE(SCIPbanditUpdateExp3IX);

/** reset callback for bandit algorithm */
SCIP_DECL_BANDITRESET(SCIPbanditResetExp3IX);

/** direct bandit creation method for the core where no SCIP pointer is available */
SCIP_RETCODE SCIPbanditCreateExp3IX(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_BANDITVTABLE*    vtable,             /**< virtual function table for callback functions of Exp.3-IX */
   SCIP_BANDIT**         exp3ix,             /**< pointer to store bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   int                   nactions,           /**< the positive number of actions for this bandit algorithm */
   unsigned int          initseed            /**< initial random seed */
   );

#ifdef __cplusplus
}
#endif

#endif
