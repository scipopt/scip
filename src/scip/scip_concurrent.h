/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_concurrent.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for concurrent solving mode
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_CONCURRENT_H__
#define __SCIP_SCIP_CONCURRENT_H__


#include "scip/def.h"
#include "scip/type_concsolver.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_syncstore.h"

/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons).
 * Additionally, the internal "set.h" is included, such that the defines in set.h are
 * available in optimized mode.
 */
#ifdef NDEBUG
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"
#include "scip/var.h"
#include "scip/cons.h"
#include "scip/solve.h"
#include "scip/debug.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicConcsolverTypeMethods
 *
 * @{
 */

/** creates a concurrent solver type and includes it in SCIP.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
EXTERN
SCIP_RETCODE SCIPincludeConcsolverType(
   SCIP*                               scip,                       /**< SCIP data structure */
   const char*                         name,                       /**< name of concurrent_solver */
   SCIP_Real                           prefpriodefault,            /**< the default preferred priority of this concurrent solver type */
   SCIP_DECL_CONCSOLVERCREATEINST      ((*concsolvercreateinst)),  /**< data copy method of concurrent solver */
   SCIP_DECL_CONCSOLVERDESTROYINST     ((*concsolverdestroyinst)), /**< data copy method of concurrent solver */
   SCIP_DECL_CONCSOLVERINITSEEDS       ((*concsolverinitseeds)),   /**< initialize random seeds of concurrent solver */
   SCIP_DECL_CONCSOLVEREXEC            ((*concsolverexec)),        /**< execution method of concurrent solver */
   SCIP_DECL_CONCSOLVERCOPYSOLVINGDATA ((*concsolvercopysolvdata)),/**< method to copy solving data */
   SCIP_DECL_CONCSOLVERSTOP            ((*concsolverstop)),        /**< terminate solving in concurrent solver */
   SCIP_DECL_CONCSOLVERSYNCWRITE       ((*concsolversyncwrite)),   /**< synchronization method of concurrent solver */
   SCIP_DECL_CONCSOLVERSYNCREAD        ((*concsolversyncread)),    /**< synchronization method of concurrent solver */
   SCIP_DECL_CONCSOLVERTYPEFREEDATA    ((*concsolvertypefreedata)),/**< method to free data of concurrent solver type */
   SCIP_CONCSOLVERTYPEDATA*            data                        /**< the concurent solver type's data */
   );

/** returns the concurrent solver type with the given name, or NULL if not existing */
EXTERN
SCIP_CONCSOLVERTYPE* SCIPfindConcsolverType(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of concurrent_solver */
   );

/** returns the array of included concurrent solver types */
EXTERN
SCIP_CONCSOLVERTYPE** SCIPgetConcsolverTypes(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of included concurrent solver types */
EXTERN
int SCIPgetNConcsolverTypes(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* @} */

/**@addtogroup PublicParallelMethods
 *
 * @{
 */

/** Constructs the parallel interface to execute processes concurrently.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPconstructSyncstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** releases the current synchronization store
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *       - \ref SCIP_STAGE_FREE
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPfreeSyncstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** Gets the synchronization store.
 *
 *  @return the \ref SCIP_SYNCSTORE parallel interface pointer to submit jobs for concurrent processing.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_SYNCSTORE* SCIPgetSyncstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
