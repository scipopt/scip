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

/**@file   struct_benderscut.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for Benders' decomposition cuts techniques
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_BENDERSCUT_H__
#define __SCIP_STRUCT_BENDERSCUT_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_benderscut.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Benders' decomposition cuts data */
struct SCIP_Benderscut
{
   SCIP_Longint          ncalls;             /**< number of times, this Benders' cut was called */
   SCIP_Longint          nfound;             /**< number of cuts found so far by this method */
   char*                 name;               /**< name of Benders' decomposition cuts */
   char*                 desc;               /**< description of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTCOPY((*benderscutcopy));/**< copy method of Benders' decomposition cuts or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSCUTFREE((*benderscutfree));/**< destructor of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTINIT((*benderscutinit));/**< initialize Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXIT((*benderscutexit));/**< deinitialize Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTINITSOL((*benderscutinitsol));/**< solving process initialization method of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXITSOL((*benderscutexitsol));/**< solving process deinitialization method of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXEC((*benderscutexec));/**< execution method of Benders' decomposition cuts */
   SCIP_BENDERSCUTDATA*  benderscutdata;     /**< Benders' decomposition cuts local data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this benderscutession for the next stages */
   SCIP_CLOCK*           benderscutclock;    /**< compression execution time */
   int                   priority;           /**< priority of the Benders' decomposition cuts */
   SCIP_Bool             islpcut;            /**< does this Benders' cut use LP information? */
   SCIP_Bool             initialized;        /**< is Benders' decomposition cuts initialized? */

   SCIP_CONS**           addedcons;          /**< an array to store the added constraints */
   SCIP_ROW**            addedcuts;          /**< an array to store the added cuts */
   int                   addedconssize;      /**< the size of the added constraint array */
   int                   addedcutssize;      /**< the size of the added cuts array */
   int                   naddedcons;         /**< the number of the added constraint */
   int                   naddedcuts;         /**< the number of the added cuts */
};

#ifdef __cplusplus
}
#endif

#endif
