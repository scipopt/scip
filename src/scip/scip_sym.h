/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_sym.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for symmetry handler plugins
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_SYM_H__
#define __SCIP_SCIP_SYM_H__

#include "scip/def.h"
#include "scip/type_sym.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/* @symtodo Replace propagator by symmetry handler */
/** creates a symmetry handler and includes it in SCIP. All non-fundamental (or optional) callbacks will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetPropInit(), SCIPsetPropExit(),
 *  SCIPsetPropCopy(), SCIPsetPropFree(), SCIPsetPropInitsol(), SCIPsetPropExitsol(),
 *  SCIPsetPropInitpre(), SCIPsetPropExitpre(), SCIPsetPropPresol(), and SCIPsetPropResprop().
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeProp() instead
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSymhdlrBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of symmetry handler */
   const char*           desc,               /**< description of symmetry handler */
   int                   priority,           /**< priority of the symmetry handler */
   int                   proppriority,       /**< priority of the symmetry handler for propagation */
   int                   sepapriority,       /**< priority of the symmetry handler for separation */
   int                   presolpriority,     /**< priority of the symmetry handler for presolving */
   int                   propfreq,           /**< frequency for calling propagator of symmetry handler */
   int                   sepafreq,           /**< frequency for calling separator of symmetry handler */
   SCIP_Bool             delayprop,          /**< should propagation be delayed, if other sym-propagators found reductions? */
   SCIP_Bool             delaysepa,          /**< should separation be delayed, if other sym-separators found reductions? */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
   int                   maxprerounds,       /**< maximal number of presolving rounds the symmetry handler participates in (-1: no limit) */
   SCIP_PROPTIMING       proptiming,         /**< positions in the node solving loop where propagation method of symmetry handlers should be executed */
   SCIP_PRESOLTIMING     presoltiming,       /**< timing mask of the symmetry handler's presolving method */
   SCIP_DECL_SYMHDLRTRYADD((*symhdlrtryadd)),/**< addition method for symmetry method handler plugins */
   SCIP_DECL_SYMHDLRCOPY ((*symcopy)),       /**< copy method of symmetry handler */
   SCIP_DECL_SYMHDLRFREE ((*symfree)),       /**< destructor method of symmetry handler */
   SCIP_DECL_SYMHDLRINIT ((*syminit)),       /**< initialization method of symmetry handler */
   SCIP_DECL_SYMHDLREXIT ((*symexit)),       /**< deinitialization method of symmetry handler */
   SCIP_DECL_SYMHDLRINITSOL((*syminitsol)),  /**< solving process initialization method of symmetry handler */
   SCIP_DECL_SYMHDLREXITSOL((*symexitsol)),  /**< solving process deinitialization method of symmetry handler */
   SCIP_DECL_SYMHDLRTRANS((*symtrans)),      /**< transformation method of symmetry handler */
   SCIP_DECL_SYMHDLRSEPALP((*symsepalp)),    /**< separator for LP solutions */
   SCIP_DECL_SYMHDLRSEPASOL((*symsepasol)),  /**< separator for arbitrary primal solutions */
   SCIP_DECL_SYMHDLRPROP ((*symprop)),       /**< propagation method of symmetry handler */
   SCIP_DECL_SYMHDLRRESPROP((*symresprop)),  /**< propagation conflict resolving method */
   SCIP_DECL_SYMHDLRPRESOL((*sympresol)),    /**< presolving method of symmetry handler */
   SCIP_DECL_SYMHDLRPRINT((*symprint)),      /**< print method of symmetry handler */
   SCIP_SYMHDLRDATA*     symhdlrdata         /**< symmetry handler data */
   );

/** sets copy method of symmetry handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetSymhdlrCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRCOPY ((*symhdlrcopy))    /**< copy method of symmetry handler */
   );

/** sets destructor method of symmetry handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetSymhdlrFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRFREE ((*symhdlrfree))    /**< destructor method of symmetry handler */
   );

/** sets initialization method of symmetry handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetSymhdlrInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRINIT ((*symhdlrinit))    /**< initialize symmetry handler */
   );

/** sets deinitialization method of symmetry handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetSymhdlrExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRINIT ((*symhdlrexit))    /**< deinitialize symmetry handler */
   );

/** sets transformation method of symmetry handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetSymhdlrTrans(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRINIT ((*symhdlrtrans))   /**< transform symmetry handler */
   );

/** sets presolving method of symmetry handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetSymhdlrPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRPRESOL((*symhdlrpresol)),/**< presolving method of symmetry handler */
   int                   presolpriority,     /**< presolving priority of the symmetry handler */
   int                   presolmaxrounds,    /**< maximal number of presolving rounds the symmetry handler participates in (-1: no limit) */
   SCIP_PRESOLTIMING     presoltiming        /**< timing mask of the symmetry handler's presolving method */
   );

/** sets propagation method of symmetry handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetSymhdlrProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRPROP ((*symhdlrprop)),   /**< propagate variable domains */
   int                   propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   SCIP_Bool             delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
   SCIP_PROPTIMING       proptiming          /**< positions in the node solving loop where propagation should be executed */
   );

/** sets all separation related callbacks/parameters of the symmetry handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetSymhdlrSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRSEPALP((*symhdlrsepalp)), /**< separate cutting planes for LP solution */
   SCIP_DECL_SYMHDLRSEPASOL((*symhdlrsepasol)), /**< separate cutting planes for arbitrary primal solution */
   int                   sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   sepapriority,       /**< priority of the constraint handler for separation */
   SCIP_Bool             delaysepa           /**< should separation method be delayed, if other separators found cuts? */
   );

/** returns the symmetry handler of the given name, or NULL if not existing */
SCIP_EXPORT
SCIP_SYMHDLR* SCIPfindSymhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of symmetry handler */
   );

/** returns the array of currently available symmetry handlers */
SCIP_EXPORT
SCIP_SYMHDLR** SCIPgetSymhdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available symmetry handlers */
SCIP_EXPORT
int SCIPgetNSymhdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** outputs symmetry component information to file stream via the message handler system
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
 *  @note If the message handler is set to a NULL pointer nothing will be printed.
 *  @note The file stream will not be flushed directly, this can be achieved by calling SCIPinfoMessage() printing a
 *        newline character.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPprintSymcomp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMCOMP*         symcomp,            /**< symmetry component */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
