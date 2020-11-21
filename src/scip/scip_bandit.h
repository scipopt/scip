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

/**@file   scip_bandit.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for bandit algorithms
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_BANDIT_H__
#define __SCIP_SCIP_BANDIT_H__


#include "scip/def.h"
#include "scip/type_bandit.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicBanditMethods
 *
 * @{
 */

/** includes a bandit algorithm virtual function table  */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBanditvtable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDITVTABLE**   banditvtable,       /**< bandit algorithm virtual function table */
   const char*           name,               /**< a name for the algorithm represented by this vtable */
   SCIP_DECL_BANDITFREE  ((*banditfree)),    /**< callback to free bandit specific data structures */
   SCIP_DECL_BANDITSELECT((*banditselect)),  /**< selection callback for bandit selector */
   SCIP_DECL_BANDITUPDATE((*banditupdate)),  /**< update callback for bandit algorithms */
   SCIP_DECL_BANDITRESET ((*banditreset))    /**< update callback for bandit algorithms */
   );

/** returns the bandit virtual function table of the given name, or NULL if not existing */
SCIP_EXPORT
SCIP_BANDITVTABLE* SCIPfindBanditvtable(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of bandit algorithm virtual function table */
   );

/** calls destructor and frees memory of bandit algorithm */
SCIP_EXPORT
SCIP_RETCODE SCIPfreeBandit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         bandit              /**< pointer to bandit algorithm data structure */
   );

/** reset the bandit algorithm */
SCIP_EXPORT
SCIP_RETCODE SCIPresetBandit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT*          bandit,             /**< pointer to bandit algorithm data structure */
   SCIP_Real*            priorities,         /**< priorities for every action, or NULL if not needed */
   unsigned int          seed                /**< initial random seed for bandit selection */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
