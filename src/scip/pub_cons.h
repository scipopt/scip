/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_cons.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for managing constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_CONS_H__
#define __SCIP_PUB_CONS_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_cons.h"

#ifdef NDEBUG
#include "scip/struct_cons.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Constraint handler methods
 */

/** compares two constraint handlers w. r. to their separation priority */
extern
SCIP_DECL_SORTPTRCOMP(SCIPconshdlrCompSepa);

/** compares two constraint handlers w. r. to their enforcing priority */
extern
SCIP_DECL_SORTPTRCOMP(SCIPconshdlrCompEnfo);

/** compares two constraint handlers w. r. to their feasibility check priority */
extern
SCIP_DECL_SORTPTRCOMP(SCIPconshdlrCompCheck);

/** gets name of constraint handler */
extern
const char* SCIPconshdlrGetName(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets description of constraint handler */
extern
const char* SCIPconshdlrGetDesc(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets user data of constraint handler */
extern
SCIP_CONSHDLRDATA* SCIPconshdlrGetData(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** sets user data of constraint handler; user has to free old data in advance! */
extern
void SCIPconshdlrSetData(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< new constraint handler user data */
   );

/** gets array with active constraints of constraint handler; a constraint is active if it is global and was not removed
 *  during presolving or it was added locally (in that case the local flag is TRUE) and the current node belongs to the
 *  corresponding sub tree
 */ 
extern
SCIP_CONS** SCIPconshdlrGetConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets array with enforced constraints of constraint handler; this is local information */
extern
SCIP_CONS** SCIPconshdlrGetEnfoConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets array with checked constraints of constraint handler; this is local information */
extern
SCIP_CONS** SCIPconshdlrGetCheckConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets total number of existing transformed constraints of constraint handler */
extern
int SCIPconshdlrGetNConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of enforced constraints of constraint handler; this is local information */
extern
int SCIPconshdlrGetNEnfoConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of checked constraints of constraint handler; this is local information */
extern
int SCIPconshdlrGetNCheckConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of active constraints of constraint handler; a constraint is active if it is global and was not removed
 *  during presolving or it was added locally (in that case the local flag is TRUE) and the current node belongs to the
 *  corresponding sub tree
 */ 
extern
int SCIPconshdlrGetNActiveConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of enabled constraints of constraint handler */
extern
int SCIPconshdlrGetNEnabledConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for presolving in this constraint handler */
extern
SCIP_Real SCIPconshdlrGetPresolTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for separation in this constraint handler */
extern
SCIP_Real SCIPconshdlrGetSepaTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for LP enforcement in this constraint handler */
extern
SCIP_Real SCIPconshdlrGetEnfoLPTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for pseudo enforcement in this constraint handler */
extern
SCIP_Real SCIPconshdlrGetEnfoPSTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for propagation in this constraint handler */
extern
SCIP_Real SCIPconshdlrGetPropTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for feasibility checking in this constraint handler */
extern
SCIP_Real SCIPconshdlrGetCheckTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for resolving propagation in this constraint handler */
extern
SCIP_Real SCIPconshdlrGetRespropTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's separation method */
extern
SCIP_Longint SCIPconshdlrGetNSepaCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's LP enforcing method */
extern
SCIP_Longint SCIPconshdlrGetNEnfoLPCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's pseudo enforcing method */
extern
SCIP_Longint SCIPconshdlrGetNEnfoPSCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's propagation method */
extern
SCIP_Longint SCIPconshdlrGetNPropCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's checking method */
extern
SCIP_Longint SCIPconshdlrGetNCheckCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's resolve propagation method */
extern
SCIP_Longint SCIPconshdlrGetNRespropCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets total number of times, this constraint handler detected a cutoff */
extern
SCIP_Longint SCIPconshdlrGetNCutoffs(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets total number of cuts found by this constraint handler */
extern
SCIP_Longint SCIPconshdlrGetNCutsFound(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets total number of additional constraints added by this constraint handler */
extern
SCIP_Longint SCIPconshdlrGetNConssFound(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets total number of domain reductions found by this constraint handler */
extern
SCIP_Longint SCIPconshdlrGetNDomredsFound(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of children created by this constraint handler */
extern
SCIP_Longint SCIPconshdlrGetNChildren(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets maximum number of active constraints of constraint handler existing at the same time */
extern
int SCIPconshdlrGetMaxNActiveConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets initial number of active constraints of constraint handler */
extern
int SCIPconshdlrGetStartNActiveConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of variables fixed in presolving method of constraint handler */
extern
int SCIPconshdlrGetNFixedVars(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of variables aggregated in presolving method of constraint handler */
extern
int SCIPconshdlrGetNAggrVars(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of variable types changed in presolving method of constraint handler */
extern
int SCIPconshdlrGetNChgVarTypes(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of bounds changed in presolving method of constraint handler */
extern
int SCIPconshdlrGetNChgBds(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of holes added to domains of variables in presolving method of constraint handler */
extern
int SCIPconshdlrGetNAddHoles(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraints deleted in presolving method of constraint handler */
extern
int SCIPconshdlrGetNDelConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraints added in presolving method of constraint handler */
extern
int SCIPconshdlrGetNAddConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraints upgraded in presolving method of constraint handler */
extern
int SCIPconshdlrGetNUpgdConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of coefficients changed in presolving method of constraint handler */
extern
int SCIPconshdlrGetNChgCoefs(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraint sides changed in presolving method of constraint handler */
extern
int SCIPconshdlrGetNChgSides(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets separation priority of constraint handler */
extern
int SCIPconshdlrGetSepaPriority(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets enforcing priority of constraint handler */
extern
int SCIPconshdlrGetEnfoPriority(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets checking priority of constraint handler */
extern
int SCIPconshdlrGetCheckPriority(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets separation frequency of constraint handler */
extern
int SCIPconshdlrGetSepaFreq(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets propagation frequency of constraint handler */
extern
int SCIPconshdlrGetPropFreq(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets frequency of constraint handler for eager evaluations in separation, propagation and enforcement */
extern
int SCIPconshdlrGetEagerFreq(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** needs constraint handler a constraint to be called? */
extern
SCIP_Bool SCIPconshdlrNeedsCons(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** does the constraint handler perform presolving? */
extern
SCIP_Bool SCIPconshdlrDoesPresolve(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** should separation method be delayed, if other separators found cuts? */
extern
SCIP_Bool SCIPconshdlrIsSeparationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** should propagation method be delayed, if other propagators found reductions? */
extern
SCIP_Bool SCIPconshdlrIsPropagationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** should presolving method be delayed, if other presolvers found reductions? */
extern
SCIP_Bool SCIPconshdlrIsPresolvingDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** was LP separation method delayed at the last call? */
extern
SCIP_Bool SCIPconshdlrWasLPSeparationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** was primal solution separation method delayed at the last call? */
extern
SCIP_Bool SCIPconshdlrWasSolSeparationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** was propagation method delayed at the last call? */
extern
SCIP_Bool SCIPconshdlrWasPropagationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** was presolving method delayed at the last call? */
extern
SCIP_Bool SCIPconshdlrWasPresolvingDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** is constraint handler initialized? */
extern
SCIP_Bool SCIPconshdlrIsInitialized(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** does the constraint handler have a copy function? */
extern
SCIP_Bool SCIPconshdlrIsClonable(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
				 );



/*
 * Constraint methods
 */

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the name of the constraint */
extern
const char* SCIPconsGetName(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns the position of constraint in the corresponding handler's conss array */
extern
int SCIPconsGetPos(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns the constraint handler of the constraint */
extern
SCIP_CONSHDLR* SCIPconsGetHdlr(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns the constraint data field of the constraint */
extern
SCIP_CONSDATA* SCIPconsGetData(
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets number of times, the constraint is currently captured */
extern
int SCIPconsGetNUses(
   SCIP_CONS*            cons                /**< constraint */
   );

/** for an active constraint, returns the depth in the tree at which the constraint was activated */
extern
int SCIPconsGetActiveDepth(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns the depth in the tree at which the constraint is valid; returns INT_MAX, if the constraint is local
 *  and currently not active
 */
extern
int SCIPconsGetValidDepth(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is active in the current node */
extern
SCIP_Bool SCIPconsIsActive(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is enabled in the current node */
extern
SCIP_Bool SCIPconsIsEnabled(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint's separation is enabled in the current node */
extern
SCIP_Bool SCIPconsIsSeparationEnabled(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint's propagation is enabled in the current node */
extern
SCIP_Bool SCIPconsIsPropagationEnabled(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is deleted or marked to be deleted */
extern
SCIP_Bool SCIPconsIsDeleted(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is marked obsolete */
extern
SCIP_Bool SCIPconsIsObsolete(
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets age of constraint */
extern
SCIP_Real SCIPconsGetAge(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff the LP relaxation of constraint should be in the initial LP */
extern
SCIP_Bool SCIPconsIsInitial(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be separated during LP processing */
extern
SCIP_Bool SCIPconsIsSeparated(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be enforced during node processing */
extern
SCIP_Bool SCIPconsIsEnforced(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be checked for feasibility */
extern
SCIP_Bool SCIPconsIsChecked(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be propagated during node processing */
extern
SCIP_Bool SCIPconsIsPropagated(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is globally valid */
extern
SCIP_Bool SCIPconsIsGlobal(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is only locally valid or not added to any (sub)problem */
extern
SCIP_Bool SCIPconsIsLocal(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is modifiable (subject to column generation) */
extern
SCIP_Bool SCIPconsIsModifiable(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is subject to aging */
extern
SCIP_Bool SCIPconsIsDynamic(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint's relaxation should be removed from the LP due to aging or cleanup */
extern
SCIP_Bool SCIPconsIsRemovable(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint's relaxation should be removed from the LP due to aging or cleanup */
extern
SCIP_Bool SCIPconsIsStickingAtNode(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint belongs to the global problem */
extern
SCIP_Bool SCIPconsIsInProb(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is belonging to original space */
extern
SCIP_Bool SCIPconsIsOriginal(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is belonging to transformed space */
extern
SCIP_Bool SCIPconsIsTransformed(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff roundings for variables in constraint are locked */
extern
SCIP_Bool SCIPconsIsLockedPos(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff roundings for variables in constraint's negation are locked */
extern
SCIP_Bool SCIPconsIsLockedNeg(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff roundings for variables in constraint or in constraint's negation are locked */
extern
SCIP_Bool SCIPconsIsLocked(
   SCIP_CONS*            cons                /**< constraint */
   );

/** get number of times the roundings for variables in constraint are locked */
extern
int SCIPconsGetNLocksPos(
   SCIP_CONS*            cons                /**< constraint */
   );

/** get number of times the roundings for variables in constraint's negation are locked */
extern
int SCIPconsGetNLocksNeg(
   SCIP_CONS*            cons                /**< constraint */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPconsGetName(cons)           (cons)->name
#define SCIPconsGetPos(cons)            (cons)->consspos
#define SCIPconsGetHdlr(cons)           (cons)->conshdlr
#define SCIPconsGetData(cons)           (cons)->consdata
#define SCIPconsGetNUses(cons)          (cons)->nuses
#define SCIPconsGetActiveDepth(cons)    (cons)->activedepth
#define SCIPconsGetValidDepth(cons)     (!(cons)->local ? 0     \
      : !SCIPconsIsActive(cons) ? INT_MAX                       \
      : (cons)->validdepth == -1 ? SCIPconsGetActiveDepth(cons) \
      : (cons)->validdepth)
#define SCIPconsIsActive(cons)          ((cons)->updateactivate || ((cons)->active && !(cons)->updatedeactivate))
#define SCIPconsIsEnabled(cons)         ((cons)->updateenable || ((cons)->enabled && !(cons)->updatedisable))
#define SCIPconsIsSeparationEnabled(cons)                               \
   (SCIPconsIsEnabled(cons) && ((cons)->updatesepaenable || ((cons)->sepaenabled && !(cons)->updatesepadisable)))
#define SCIPconsIsPropagationEnabled(cons)                              \
   (SCIPconsIsEnabled(cons) && ((cons)->updatepropenable || ((cons)->propenabled && !(cons)->updatepropdisable)))
#define SCIPconsIsDeleted(cons)         ((cons)->deleted)
#define SCIPconsIsObsolete(cons)        ((cons)->updateobsolete || (cons)->obsolete)
#define SCIPconsGetAge(cons)            (cons)->age
#define SCIPconsIsInitial(cons)         (cons)->initial
#define SCIPconsIsSeparated(cons)       (cons)->separate
#define SCIPconsIsEnforced(cons)        (cons)->enforce
#define SCIPconsIsChecked(cons)         (cons)->check
#define SCIPconsIsPropagated(cons)      (cons)->propagate
#define SCIPconsIsGlobal(cons)          !(cons)->local
#define SCIPconsIsLocal(cons)           (cons)->local
#define SCIPconsIsModifiable(cons)      (cons)->modifiable
#define SCIPconsIsDynamic(cons)         (cons)->dynamic
#define SCIPconsIsRemovable(cons)       (cons)->removable
#define SCIPconsIsStickingAtNode(cons)  (cons)->stickingatnode
#define SCIPconsIsInProb(cons)          ((cons)->addconssetchg == NULL && (cons)->addarraypos >= 0)
#define SCIPconsIsOriginal(cons)        (cons)->original
#define SCIPconsIsTransformed(cons)     !(cons)->original
#define SCIPconsIsLockedPos(cons)       ((cons)->nlockspos > 0)
#define SCIPconsIsLockedNeg(cons)       ((cons)->nlocksneg > 0)
#define SCIPconsIsLocked(cons)          ((cons)->nlockspos > 0 || (cons)->nlocksneg > 0)
#define SCIPconsGetNLocksPos(cons)      ((cons)->nlockspos)
#define SCIPconsGetNLocksNeg(cons)      ((cons)->nlocksneg)

#endif

#ifdef __cplusplus
}
#endif

#endif
