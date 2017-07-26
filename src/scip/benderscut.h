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

/**@file   benderscut.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for Benders' decomposition cuts
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERSCUT_H__
#define __SCIP_BENDERSCUT_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_benderscut.h"
#include "scip/pub_benderscut.h"

#ifdef __cplusplus
extern "C" {
#endif

/** copies the given Benders' decomposition cuts to a new scip */
extern
SCIP_RETCODE SCIPbenderscutCopyInclude(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a Benders' decomposition cuts */
extern
SCIP_RETCODE SCIPbenderscutCreate(
   SCIP_BENDERSCUT**     benderscut,         /**< pointer to Benders' decomposition cuts data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of Benders' decomposition cuts */
   const char*           desc,               /**< description of Benders' decomposition cuts */
   int                   priority,           /**< priority of the Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTCOPY((*benderscutcopy)),/**< copy method of Benders' decomposition cuts or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSCUTFREE((*benderscutfree)),/**< destructor of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTINIT((*benderscutinit)),/**< initialize Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXIT((*benderscutexit)),/**< deinitialize Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTINITSOL((*benderscutinitsol)),/**< solving process initialization method of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXITSOL((*benderscutexitsol)),/**< solving process deinitialization method of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXEC((*benderscutexec)),/**< execution method of Benders' decomposition cuts */
   SCIP_BENDERSCUTDATA*  benderscutdata      /**< Benders' decomposition cuts data */
   );

/** calls destructor and frees memory of Benders' decomposition cuts */
extern
SCIP_RETCODE SCIPbenderscutFree(
   SCIP_BENDERSCUT**     benderscut,         /**< pointer to Benders' decomposition cuts data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes Benders' decomposition cuts */
extern
SCIP_RETCODE SCIPbenderscutInit(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of Benders' decomposition cuts */
extern
SCIP_RETCODE SCIPbenderscutExit(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs Benders' decomposition cuts that the branch and bound process is being started */
extern
SCIP_RETCODE SCIPbenderscutInitsol(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs Benders' decomposition cuts that the branch and bound process data is being freed */
extern
SCIP_RETCODE SCIPbenderscutExitsol(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of Benders' decomposition cuts */
extern
SCIP_RETCODE SCIPbenderscutExec(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the number of the subproblem for which the cut is generated */
   SCIP_BENDERSENFOTYPE  type,               /**< the enforcement type calling this function */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of Benders' decomposition cuts */
extern
void SCIPbenderscutSetPriority(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   int                   priority            /**< new priority of the Benders' decomposition cuts */
   );

/** sets copy callback of Benders' decomposition cuts */
extern
void SCIPbenderscutSetCopy(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTCOPY((*benderscutcopy))/**< copy callback of Benders' decomposition cuts or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor callback of Benders' decomposition cuts */
extern
void SCIPbenderscutSetFree(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTFREE((*benderscutfree))/**< destructor of Benders' decomposition cuts */
   );

/** sets initialization callback of Benders' decomposition cuts */
extern
void SCIPbenderscutSetInit(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTINIT((*benderscutinit))/**< initialize Benders' decomposition cuts */
   );

/** sets deinitialization callback of Benders' decomposition cuts */
extern
void SCIPbenderscutSetExit(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXIT((*benderscutexit))/**< deinitialize Benders' decomposition cuts */
   );

/** sets solving process initialization callback of Benders' decomposition cuts */
extern
void SCIPbenderscutSetInitsol(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTINITSOL((*benderscutinitsol))/**< solving process initialization callback of Benders' decomposition cuts */
   );

/** sets solving process deinitialization callback of Benders' decomposition cuts */
extern
void SCIPbenderscutSetExitsol(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXITSOL((*benderscutexitsol))/**< solving process deinitialization callback of Benders' decomposition cuts */
   );

/** adds the generated constraint to the Benders cut storage */
extern
SCIP_RETCODE SCIPbenderscutStoreCons(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS*            cons                /**< the constraint to be added to the Benders' cut storage */
   );

/** adds the generated cuts to the Benders' cut storage */
extern
SCIP_RETCODE SCIPbenderscutStoreCut(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             cut                 /**< the cut to be added to the Benders' cut storage */
   );

#ifdef __cplusplus
}
#endif

#endif
