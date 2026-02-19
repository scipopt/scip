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

/**@file   sym.h
 * @ingroup INTERNALAPI
 * @brief  methods for symmetry handlers and symmetry components
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SYM_H__
#define __SCIP_SYM_H__

#include "scip/type_sym.h"
#include "scip/pub_sym.h"
#include "scip/type_sepastore.h"
#include "symmetry/type_symmetry.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates a symmetry handler */
SCIP_RETCODE SCIPsymhdlrCreate(
   SCIP_SYMHDLR**        symhdlr,            /**< pointer to symmetry handler data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
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
   SCIP_DECL_SYMHDLRTRYADD((*symtryadd)),    /**< addition method for symmetry method handler plugins */
   SCIP_DECL_SYMHDLRCOPY ((*symcopy)),       /**< copy method of symmetry handler */
   SCIP_DECL_SYMHDLRFREE ((*symfree)),       /**< destructor method of symmetry handler */
   SCIP_DECL_SYMHDLRINIT ((*syminit)),       /**< initialization method of symmetry handler */
   SCIP_DECL_SYMHDLREXIT ((*symexit)),       /**< deinitialization method of symmetry handler */
   SCIP_DECL_SYMHDLRINITSOL((*syminitsol)),  /**< solving process initialization method of symmetry handler */
   SCIP_DECL_SYMHDLREXITSOL((*symexitsol)),  /**< solving process deinitialization method of symmetry handler */
   SCIP_DECL_SYMHDLRTRANS((*symtrans)),      /**< transformation method of symmetry hanlder */
   SCIP_DECL_SYMHDLRSEPALP((*symsepalp)),    /**< separator for LP solutions */
   SCIP_DECL_SYMHDLRSEPASOL((*symsepasol)),  /**< separator for arbitrary primal solutions */
   SCIP_DECL_SYMHDLRPROP ((*symprop)),       /**< propagation method of symmetry handler */
   SCIP_DECL_SYMHDLRRESPROP((*symresprop)),  /**< propagation conflict resolving method */
   SCIP_DECL_SYMHDLRPRESOL((*sympresol)),    /**< presolving method of symmetry handler */
   SCIP_DECL_SYMHDLRPRINT((*symprint)),      /**< print method of symmetry handler */
   SCIP_SYMHDLRDATA*     symhdlrdata         /**< symmetry handler data */
   );

/** copies the given symmetry handler to a new scip */
SCIP_RETCODE SCIPsymhdlrCopyInclude(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** calls destructor and frees memory of symmetry handler */
SCIP_RETCODE SCIPsymhdlrFree(
   SCIP_SYMHDLR**        symhdlr,            /**< pointer to symmetry handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrExit(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrInit(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs symmetry handler that the branch and bound process is being started */
SCIP_RETCODE SCIPsymhdlrInitsol(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs symmetry handler rule that the branch and bound process data is being freed */
SCIP_RETCODE SCIPsymhdlrExitsol(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             restart             /**< was this exit solve call triggered by a restart? */
   );

/** executes presolving method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrPresol(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRESOLTIMING     timing,             /**< current presolving timing */
   int                   nrounds,            /**< number of presolving rounds already done */
   int*                  nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*                  naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*                  nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*                  nchgbds,            /**< pointer to total number of variable bounds tightened of all presolvers */
   int*                  naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*                  ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*                  naddconss,          /**< pointer to total number of added constraints of all presolvers */
   int*                  nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*                  nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*                  nchgsides,          /**< pointer to total number of changed left/right hand sides of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls propagation method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrProp(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             execdelayed,        /**< execute propagator even if it is marked to be delayed */
   SCIP_Bool             instrongbranching,  /**< are we currently doing strong branching? */
   SCIP_PROPTIMING       proptiming,         /**< current point in the node solving process */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** resolves the given conflicting bound, that was deduced by the given propagator, by putting all "reason" bounds
 *  leading to the deduction into the conflict queue with calls to SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(),
 *  SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or SCIPaddConflictBinvar();
 *
 *  @note it is sufficient to explain the relaxed bound change
 */
SCIP_RETCODE SCIPsymhdlrResolvePropagation(
   SCIP_SYMCOMP*         symcomp,            /**< symmetry component */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             infervar,           /**< variable whose bound was deduced by the constraint */
   int                   inferinfo,          /**< user inference information attached to the bound change */
   SCIP_BOUNDTYPE        inferboundtype,     /**< bound that was deduced (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index, representing the point of time where change took place */
   SCIP_Real             relaxedbd,          /**< the relaxed bound */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls LP separation method of symmetry handler's separator */
SCIP_RETCODE SCIPsymhdlrSepaLP(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   int                   depth,              /**< depth of current node */
   SCIP_Real             bounddist,          /**< current relative distance of local dual bound to global dual bound */
   SCIP_Bool             allowlocal,         /**< should the separator be asked to separate local cuts */
   SCIP_Bool             execdelayed,        /**< execute separator even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls primal solution separation method of symmetry handler's separator */
SCIP_RETCODE SCIPsymhdlrSepaSol(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             allowlocal,         /**< should the separator allow local cuts */
   SCIP_Bool             execdelayed,        /**< execute separator even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls try-add method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrTryAdd(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int**                 symmetries,         /**< array of symmetries */
   int                   nsymmetries,        /**< number of symmetries in symmetries array */
   SYM_SYMTYPE           symtype,            /**< type of symmetry */
   SCIP_VAR**            symvars,            /**< variables on which symmetries act */
   int                   nsymvars,           /**< number of variables in symvars */
   SCIP_Real*            symvardomcenter,    /**< domain center of variables (or NULL) */
   SCIP_HASHMAP*         symvarmap,          /**< map of variables to indices in permvars array */
   SYM_GRAPH*            symgraph,           /**< symmetry detection graph (or NULL) */
   int                   id,                 /**< identifier of component for which symmetry handling shall be added */
   SCIP_SYMCOMPDATA**    symcompdata,        /**< pointer to store data of symmetry component */
   int*                  naddedconss,        /**< pointer to store number of constraints added by symhdlr */
   int*                  nchgbds,            /**< pointer to store number of changed variable bounds */
   SCIP_Bool*            success             /**< pointer to store whether symmetry handling method could be added */
   );

/** adds a component to a symmetry handler */
SCIP_RETCODE SCIPaddSymhdlrComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SYMCOMP*         symcomp             /**< symmetry component */
   );

/** creates a symmetry component */
SCIP_RETCODE SCIPcreateSymmetryComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMCOMP**        symcomp,            /**< pointer to symmetry component */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler active on symmetry component */
   SCIP_SYMCOMPDATA*     symcompdata,        /**< symmetry component data */
   int                   id                  /**< numerical identifier of symmetry component */
   );

/** creates and captures symmetry information data structure */
SCIP_RETCODE SCIPsyminfoCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMINFO**        syminfo             /**< pointer to return the created syminfo */
   );

/** returns the symmetry information data structure */
SCIP_SYMINFO* SCIPgetSyminfo(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** releases symmetry information data structure */
SCIP_RETCODE SCIPsyminfoFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMINFO**        syminfo             /**< pointer to the syminfo */
   );

/** outputs symmetry component information to file stream */
SCIP_RETCODE SCIPsymcompPrint(
   SCIP_SYMCOMP*         symcomp,            /**< symmetry component to print */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

#ifdef __cplusplus
}
#endif

#endif
