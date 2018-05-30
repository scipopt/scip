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

/**@file   scip_relax.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for relaxator plugins
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

#ifndef __SCIP_SCIP_RELAX_H__
#define __SCIP_SCIP_RELAX_H__


#include "scip/def.h"
#include "scip/type_relax.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

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

/**@addtogroup PublicRelaxatorMethods
 *
 * @{
 */

/** creates a relaxation handler and includes it in SCIP
 *
 *  @note method has all relaxation handler callbacks as arguments and is thus changed every time a new
 *        callback is added
 *        in future releases; consider using SCIPincludeRelaxBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
EXTERN
SCIP_RETCODE SCIPincludeRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of relaxation handler */
   const char*           desc,               /**< description of relaxation handler */
   int                   priority,           /**< priority of the relaxation handler (negative: after LP, non-negative: before LP) */
   int                   freq,               /**< frequency for calling relaxation handler */
   SCIP_DECL_RELAXCOPY   ((*relaxcopy)),     /**< copy method of relaxation handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_RELAXFREE   ((*relaxfree)),     /**< destructor of relaxation handler */
   SCIP_DECL_RELAXINIT   ((*relaxinit)),     /**< initialize relaxation handler */
   SCIP_DECL_RELAXEXIT   ((*relaxexit)),     /**< deinitialize relaxation handler */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol)),  /**< solving process initialization method of relaxation handler */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol)),  /**< solving process deinitialization method of relaxation handler */
   SCIP_DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxation handler */
   SCIP_RELAXDATA*       relaxdata           /**< relaxation handler data */
   );

/** creates a relaxation handler and includes it in SCIP. All non fundamental
 *  (or optional) callbacks as, e.g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetRelaxInit(), SCIPsetRelaxExit(),
 *  SCIPsetRelaxCopy(), SCIPsetRelaxFree(), SCIPsetRelaxInitsol(), and SCIPsetRelaxExitsol()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeRelax() instead
 */
EXTERN
SCIP_RETCODE SCIPincludeRelaxBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX**          relaxptr,           /**< reference to relaxation pointer, or NULL */
   const char*           name,               /**< name of relaxation handler */
   const char*           desc,               /**< description of relaxation handler */
   int                   priority,           /**< priority of the relaxation handler (negative: after LP, non-negative: before LP) */
   int                   freq,               /**< frequency for calling relaxation handler */
   SCIP_DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxation handler */
   SCIP_RELAXDATA*       relaxdata           /**< relaxation handler data */
   );

/** sets copy method of relaxation handler */
EXTERN
SCIP_RETCODE SCIPsetRelaxCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_DECL_RELAXCOPY   ((*relaxcopy))      /**< copy method of relaxation handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of relaxation handler */
EXTERN
SCIP_RETCODE SCIPsetRelaxFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_DECL_RELAXFREE   ((*relaxfree))      /**< destructor of relaxation handler */
   );

/** sets initialization method of relaxation handler */
EXTERN
SCIP_RETCODE SCIPsetRelaxInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_DECL_RELAXINIT   ((*relaxinit))      /**< initialize relaxation handler */
   );

/** sets deinitialization method of relaxation handler */
EXTERN
SCIP_RETCODE SCIPsetRelaxExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_DECL_RELAXEXIT   ((*relaxexit))      /**< deinitialize relaxation handler */
   );

/** sets solving process initialization method of relaxation handler */
EXTERN
SCIP_RETCODE SCIPsetRelaxInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol))   /**< solving process initialization method of relaxation handler */
   );

/** sets solving process deinitialization method of relaxation handler */
EXTERN
SCIP_RETCODE SCIPsetRelaxExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol))   /**< solving process deinitialization method of relaxation handler */
   );

/** returns the relaxation handler of the given name, or NULL if not existing */
EXTERN
SCIP_RELAX* SCIPfindRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of relaxation handler */
   );

/** returns the array of currently available relaxation handlers  */
EXTERN
SCIP_RELAX** SCIPgetRelaxs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available relaxation handlers  */
EXTERN
int SCIPgetNRelaxs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a relaxation handler*/
EXTERN
SCIP_RETCODE SCIPsetRelaxPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   int                   priority            /**< new priority of the relaxation handler */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
