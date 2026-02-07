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

/**@file   type_sym.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for symmetry group and symmetry handlers
 * @author Christopher Hojny
 *
 *  This file defines the interface for symmetry handlers implemented in C.
 */

/** @defgroup DEFPLUGINS_SYM Default symmetry handlers
 *  @ingroup DEFPLUGINS
 *  @brief implementation files (.c files) of the default symmetry handlers of SCIP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SYM_H__
#define __SCIP_TYPE_SYM_H__

#include "scip/type_cons.h"
#include "scip/type_misc.h"
#include "symmetry/type_symmetry.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Symhdlr SCIP_SYMHDLR;    /**< symmetry handler for a specific symmetry method */
typedef struct SCIP_SymComp SCIP_SYMCOMP;    /**< symmetry component structure */
typedef struct SCIP_SymhdlrData SCIP_SYMHDLRDATA; /**< symmetry handler data */
typedef struct SCIP_SymCompData SCIP_SYMCOMPDATA; /**< symmetry component data */
typedef struct SCIP_SymInfo SCIP_SYMINFO;    /**< data structure for storing symmetry information */

/** addition method for symmetry method handler plugins (tries to add symmetry handling method for given symmetries) *
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 *  - symtype         : type of symmetry
 *  - perms           : array of (signed) permutation symmetries for which applicability of symmetry handler is tested
 *  - nperms          : number of (signed) permutations in perms array
 *  - permvars        : array of variables on which (signed) permutations operate
 *  - npermvars       : number of variables in permvars array
 *  - permvardomcenter: array containing centers of variable domains (or NULL if symtype != SYM_SYMTYPE_SIGNPERM)
 *  - permvarmap      : map of variables to indices in permvars array
 *  - symgraph        : symmetry detection graph used for detecting symmetries (or NULL)
 *  - id              : numerical identifier of set of permutation symmetries for which symmetry handler is tested
 *  - symcompdata     : pointer for storing data of symmetry component
 *  - naddedconss     : pointer to store number of constraints added by symhdlr
 *  - nchgbds         : pointer to store number of changed variable bounds
 *  - success         : pointer to store whether the symmetry handling method has been added
 */
#define SCIP_DECL_SYMHDLRTRYADD(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr, SYM_SYMTYPE symtype, \
      int** perms, int nperms, SCIP_VAR** permvars, int npermvars, SCIP_Real* permvardomcenter, \
      SCIP_HASHMAP* permvarmap, SYM_GRAPH* symgraph, int id, SCIP_SYMCOMPDATA** symcompdata, int* naddedconss, \
      int* nchgbds, SCIP_Bool* success)

/** copy method for symmetry handler plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 */
#define SCIP_DECL_SYMHDLRCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr)

/** destructor of symmetry handler to free symmetry handler data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 */
#define SCIP_DECL_SYMHDLRFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr)

/** initialization method of symmetry handler (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 *  - symcomps        : array of symmetry components handled by symhdlr
 *  - nsymcomps       : number of symmetry components handled by symhdlr
 */
#define SCIP_DECL_SYMHDLRINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr, SCIP_SYMCOMP** symcomps, int nsymcomps)

/** deinitialization method of symmetry handler (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 *  - symcomps        : array of symmetry components handled by symhdlr
 *  - nsymcomps       : number of symmetry components handled by symhdlr
 */
#define SCIP_DECL_SYMHDLREXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr, SCIP_SYMCOMP** symcomps, int nsymcomps)

/** solving process initialization method of symmetry handler (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The symmetry handler may use this call to initialize its branch and bound specific data.
 *
 *  Besides necessary modifications and clean up, no time consuming operations should be performed, especially if the
 *  problem has already been solved.  Use the method SCIPgetStatus(), which in this case returns SCIP_STATUS_OPTIMAL,
 *  SCIP_STATUS_INFEASIBLE, SCIP_STATUS_UNBOUNDED, or SCIP_STATUS_INFORUNBD.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 *  - symcomps        : array of symmetry components handled by symhdlr
 *  - nsymcomps       : number of symmetry components handled by symhdlr
 */
#define SCIP_DECL_SYMHDLRINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr, SCIP_SYMCOMP** symcomps, int nsymcomps)

/** solving process deinitialization method of symmetry handler (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The symmetry handler should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 *  - symcomps        : array of symmetry components handled by symhdlr
 *  - nsymcomps       : number of symmetry components handled by symhdlr
 *  - restart         : was this exit solve call triggered by a restart?
 */
#define SCIP_DECL_SYMHDLREXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr, SCIP_SYMCOMP** symcomps, int nsymcomps, SCIP_Bool restart)

/* @symtodo do we need this method? */
/** frees specific symmetry component data
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 *  - symcomp         : pointer to symmetry component to free
 */
#define SCIP_DECL_SYMHDLRDELETE(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr, SCIP_SYMCOMP** symcomp)

/** transforms data of symmetry handler into data belonging to the transformed problem
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : symmetry handler
 *  - sourcesymcomp   : symmetry component to be transformed
 *  - targetsymcomp   : pointer to store transformed symmetry component
 */
#define SCIP_DECL_SYMHDLRTRANS(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr, SCIP_SYMCOMP* sourcesymcomp, SCIP_SYMCOMP** targetsymcomp)

/** LP solution separation method of symmetry handler
 *
 *  Searches for cutting planes that separate the current LP solution. The method is called in the LP solving loop,
 *  which means that a valid LP solution exists.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 *  - symcomps        : array of symmetry components handled by symhdlr
 *  - nsymcomps       : number of symmetry components handled by symhdlr
 *  - result          : pointer to store the result of the separation call
 *  - allowlocal      : should the separator allow local cuts?
 *  - depth           : pretended depth of current node
 *
 *  @note The depth argument shouldn't be use to determine whether the cut is globally valid or not.  The value of depth
 *  could be 0 even though we are not in the root node! The purpose of depth is to control the behavior of the
 *  separator. Usually separators will have different limits on the number of cuts to be applied in the root node, etc.
 *  These limits should be checked against depth and not against the actual depth of the current node.
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_NEWROUND   : a cutting plane was generated and a new separation round should immediately start
 *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *  - SCIP_DIDNOTRUN  : the separator was skipped
 *  - SCIP_DELAYED    : the separator was skipped, but should be called again
 */
#define SCIP_DECL_SYMHDLRSEPALP(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr, SCIP_SYMCOMP** symcomps, \
      int nsymcomps, SCIP_RESULT* result, SCIP_Bool allowlocal, int depth)

/** arbitrary primal solution separation method of symmetry handler
 *
 *  Searches for cutting planes that separate the given primal solution. The method is called outside the LP solution
 *  loop (e.g., by a relaxator or a primal heuristic), which means that there is no valid LP solution.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 *  - sol             : primal solution that should be separated
 *  - symcomps        : array of symmetry components handled by symhdlr
 *  - nsymcomps       : number of symmetry components handled by symhdlr
 *  - result          : pointer to store the result of the separation call
 *  - allowlocal      : should the separator allow local cuts?
 *  - depth           : pretended depth of current node
 *
 *  @note The depth argument shouldn't be use to determine whether the cut is globally valid or not.  The value of depth
 *  could be 0 even though we are not in the root node! The purpose of depth is to control the behavior of the
 *  separator. Usually separators will have different limits on the number of cuts to be applied in the root node, etc.
 *  These limits should be checked against depth and not against the actual depth of the current node.
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_NEWROUND   : a cutting plane was generated and a new separation round should immediately start
 *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *  - SCIP_DIDNOTRUN  : the separator was skipped
 *  - SCIP_DELAYED    : the separator was skipped, but should be called again
 */
#define SCIP_DECL_SYMHDLRSEPASOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr, SCIP_SOL* sol, \
      SCIP_SYMCOMP** symcomps, int nsymcomps, SCIP_RESULT* result, SCIP_Bool allowlocal, int depth)

/** domain propagation method of symmetry handler
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 *  - symcomps        : array of symmetry components handled by symhdlr
 *  - nsymcomps       : number of symmetry components handled by symhdlr
 *  - proptiming      : current point in the node solving loop
 *  - result          : pointer to store the result of the propagation call
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_REDUCEDDOM : at least one domain reduction was found
 *  - SCIP_DIDNOTFIND : the propagator searched but did not find any domain reductions
 *  - SCIP_DIDNOTRUN  : the propagator was skipped
 *  - SCIP_DELAYED    : the propagator was skipped, but should be called again
 *  - SCIP_DELAYNODE  : the current node should be postponed (return value only valid for BEFORELP propagation)
 */
#define SCIP_DECL_SYMHDLRPROP(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr, SCIP_SYMCOMP** symcomps, \
      int nsymcomps, SCIP_PROPTIMING proptiming, SCIP_RESULT* result)

/** propagation conflict resolving method of symmetry handler
 *
 *  This method is called during conflict analysis. If the propagator of a symmetry handler wants to support
 *  conflict analysis, it should call SCIPinferVarLbSym() or SCIPinferVarUbSym() in domain propagation instead
 *  of SCIPchgVarLb() or SCIPchgVarUb() in order to deduce bound changes on variables.
 *  In the SCIPinferVarLbSym() and SCIPinferVarUbSym() calls, the symmetry handler provides a pointer to the
 *  symmetry component that is handled by the symmetry handler and an integer value "inferinfo" that can be
 *  arbitrarily chosen.
 *  The propagation conflict resolving method can then be implemented, to provide a "reasons" for the bound
 *  changes, i.e. the bounds of variables at the time of the propagation, that forced the symmetry handler to set the
 *  conflict variable's bound to its current value. It can use the "inferinfo" tag to identify its own propagation
 *  rule and thus identify the "reason" bounds. The bounds that form the reason of the assignment must then be provided
 *  by calls to SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(),
 *  SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), and/or SCIPaddConflictBinvar() in the propagation conflict
 *  resolving method.
 *
 *  See the description of the propagation conflict resolving method of constraint handlers for further details.
 *
 *  @note if the propagtor uses dual information it is nesassary to check via calling SCIPallowWeakDualReds and
 *        SCIPallowStrongDualReds if dual reductions and propgation with the current cutoff bound, resp., are allowed.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 *  - symcomp         : symmetry component that was responsible for the propagation
 *  - infervar        : the conflict variable whose bound change has to be resolved
 *  - inferinfo       : the user information passed to the corresponding SCIPinferVarLbProp() or SCIPinferVarUbProp() call
 *  - boundtype       : the type of the changed bound (lower or upper bound)
 *  - bdchgidx        : the index of the bound change, representing the point of time where the change took place
 *  - relaxedbd       : the relaxed bound which is sufficient to be explained
 *
 *  output:
 *  - result          : pointer to store the result of the propagation conflict resolving call
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the conflicting bound change has been successfully resolved by adding all reason bounds
 *  - SCIP_DIDNOTFIND : the conflicting bound change could not be resolved and has to be put into the conflict set
 *
 *  @note it is sufficient to explain/resolve the relaxed bound
 */
#define SCIP_DECL_SYMHDLRRESPROP(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdr, SCIP_SYMCOMP* symcomp, \
      SCIP_VAR* infervar, int inferinfo, SCIP_BOUNDTYPE boundtype, SCIP_BDCHGIDX* bdchgidx, SCIP_Real relaxedbd, \
      SCIP_RESULT* result)

/** presolving method of symmetry handler
 *
 *  The presolver should go through the variables and constraints and tighten the domains or
 *  constraints. Each tightening should increase the given total number of changes.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 *  - symcomps        : array of symmetry components handled by symhdlr
 *  - nsymcomps       : number of symmetry components handled by symhdlr
 *  - nrounds         : number of presolving rounds already done
 *  - presoltiming    : current presolving timing
 *  - nnewfixedvars   : number of variables fixed since the last call to the presolving method
 *  - nnewaggrvars    : number of variables aggregated since the last call to the presolving method
 *  - nnewchgvartypes : number of variable type changes since the last call to the presolving method
 *  - nnewchgbds      : number of variable bounds tightened since the last call to the presolving method
 *  - nnewholes       : number of domain holes added since the last call to the presolving method
 *  - nnewdelconss    : number of deleted constraints since the last call to the presolving method
 *  - nnewaddconss    : number of added constraints since the last call to the presolving method
 *  - nnewupgdconss   : number of upgraded constraints since the last call to the presolving method
 *  - nnewchgcoefs    : number of changed coefficients since the last call to the presolving method
 *  - nnewchgsides    : number of changed left or right hand sides since the last call to the presolving method
 *
 *  @note the counters state the changes since the last call including the changes of this presolving method during its
 *        call
 *
 *  @note if the constraint handler performs dual presolving it is nesassary to check via calling SCIPallowWeakDualReds
 *        and SCIPallowStrongDualReds if dual reductions are allowed.
 *
 *  input/output:
 *  - nfixedvars      : pointer to count total number of variables fixed of all presolvers
 *  - naggrvars       : pointer to count total number of variables aggregated of all presolvers
 *  - nchgvartypes    : pointer to count total number of variable type changes of all presolvers
 *  - nchgbds         : pointer to count total number of variable bounds tightened of all presolvers
 *  - naddholes       : pointer to count total number of domain holes added of all presolvers
 *  - ndelconss       : pointer to count total number of deleted constraints of all presolvers
 *  - naddconss       : pointer to count total number of added constraints of all presolvers
 *  - nupgdconss      : pointer to count total number of upgraded constraints of all presolvers
 *  - nchgcoefs       : pointer to count total number of changed coefficients of all presolvers
 *  - nchgsides       : pointer to count total number of changed left/right hand sides of all presolvers
 *
 *  output:
 *  - result          : pointer to store the result of the presolving call
 *
 *  possible return values for *result:
 *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
 *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
 *  - SCIP_SUCCESS    : the presolving method found a reduction
 *  - SCIP_DIDNOTFIND : the presolving method searched, but did not find a presolving change
 *  - SCIP_DIDNOTRUN  : the presolving method was skipped
 *  - SCIP_DELAYED    : the presolving method was skipped, but should be called again
 */
#define SCIP_DECL_SYMHDLRPRESOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr, \
      SCIP_SYMCOMP** symcomps, int nsymcomps, int nrounds, SCIP_PRESOLTIMING presoltiming, \
      int nnewfixedvars, int nnewaggrvars, int nnewchgvartypes, int nnewchgbds, int nnewholes, \
      int nnewdelconss, int nnewaddconss, int nnewupgdconss, int nnewchgcoefs, int nnewchgsides, \
      int* nfixedvars, int* naggrvars, int* nchgvartypes, int* nchgbds, int* naddholes, \
      int* ndelconss, int* naddconss, int* nupgdconss, int* nchgcoefs, int* nchgsides, SCIP_RESULT* result)

/** symmetry component display method of symmetry handler
 *
 *  The symmetry handler can store a representation of the symmetry components into the given text file. Use the method
 *  SCIPinfoMessage() to push a string into the file stream.
 *
 *  @note There are several methods which help to display variables. These are SCIPwriteVarName(), SCIPwriteVarsList(),
 *        SCIPwriteVarsLinearsum(), and SCIPwriteVarsPolynomial().
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - symhdlr         : the symmetry handler itself
 *  - symcomp         : the symmetry component that should be displayed
 *  - file            : the text file to store the information into
 */
#define SCIP_DECL_SYMHDLRPRINT(x) SCIP_RETCODE x (SCIP* scip, SCIP_SYMHDLR* symhdlr, SCIP_SYMCOMP* symcomp, FILE* file)

#ifdef __cplusplus
}
#endif

#endif
