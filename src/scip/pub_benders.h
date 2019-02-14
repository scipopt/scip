/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_benders.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_BENDERS_H__
#define __SCIP_PUB_BENDERS_H__

#include "scip/def.h"
#include "scip/type_benders.h"
#include "scip/type_benderscut.h"
#include "scip/type_misc.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"
#include "scip/type_stat.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicBendersMethods
 *
 * @{
 */

/** compares two benderss w. r. to their priority */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPbendersComp);

/** comparison method for sorting benderss w.r.t. to their name */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPbendersCompName);

/** gets user data of Benders' decomposition */
EXTERN
SCIP_BENDERSDATA* SCIPbendersGetData(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** sets user data of Benders' decomposition; user has to free old data in advance! */
EXTERN
void SCIPbendersSetData(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSDATA*     bendersdata         /**< new Benders' decomposition user data */
   );

/** gets name of Benders' decomposition */
EXTERN
const char* SCIPbendersGetName(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets description of Benders' decomposition */
EXTERN
const char* SCIPbendersGetDesc(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets priority of Benders' decomposition */
EXTERN
int SCIPbendersGetPriority(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets the number of subproblems for the Benders' decomposition */
EXTERN
int SCIPbendersGetNSubproblems(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition data structure */
   );

/** returns the SCIP instance for a given subproblem */
EXTERN
SCIP* SCIPbendersSubproblem(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber          /**< the subproblem number */
   );

/** gets the number of times, the Bender' decomposition was called and tried to find a violated second stage constraint */
EXTERN
int SCIPbendersGetNCalls(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets the number of optimality cuts found by the collection of Benders' decomposition subproblems */
EXTERN
int SCIPbendersGetNCutsFound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets the number of cuts found from the strengthening round */
int SCIPbendersGetNStrengthenCutsFound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets the number of calls to the strengthening round */
int SCIPbendersGetNStrengthenCalls(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets the number of calls to the strengthening round that fail */
int SCIPbendersGetNStrengthenFails(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets time in seconds used in this Benders' decomposition for setting up for next stages */
EXTERN
SCIP_Real SCIPbendersGetSetupTime(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets execution time in seconds used in this Benders' decomposition */
EXTERN
SCIP_Real SCIPbendersGetTime(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** Is Benders' decomposition initialized? */
EXTERN
SCIP_Bool SCIPbendersIsInitialized(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** returns whether the given Benders' decomposition is in use in the current problem */
SCIP_Bool SCIPbendersIsActive(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   );

/** Returns whether only the convex relaxations will be checked in this solve loop
 *  when Benders' is used in the LNS heuristics, only the convex relaxations of the master/subproblems are checked,
 *  i.e. no integer cuts are generated. In this case, then Benders' decomposition is performed under the assumption
 *  that all subproblems are convex relaxations.
 */
SCIP_Bool SCIPbendersOnlyCheckConvexRelax(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_Bool             subscipsoff         /**< flag indicating whether plugins using sub-SCIPs are deactivated */
   );

/** Are Benders' cuts generated from the LP solutions? */
EXTERN
SCIP_Bool SCIPbendersCutLP(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** Are Benders' cuts generated from the pseudo solutions? */
EXTERN
SCIP_Bool SCIPbendersCutPseudo(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** Are Benders' cuts generated from the relaxation solutions? */
EXTERN
SCIP_Bool SCIPbendersCutRelaxation(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** Should this Benders' use the auxiliary variables from the highest priority Benders'? */
EXTERN
SCIP_Bool SCIPbendersShareAuxVars(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** returns the auxiliary variable for the given subproblem */
EXTERN
SCIP_VAR* SCIPbendersGetAuxiliaryVar(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** returns all auxiliary variables */
EXTERN
SCIP_VAR** SCIPbendersGetAuxiliaryVars(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** stores the objective function value of the subproblem for use in cut generation */
EXTERN
void SCIPbendersSetSubproblemObjval(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Real             objval              /**< the objective function value for the subproblem */
   );

/** returns the objective function value of the subproblem for use in cut generation */
EXTERN
SCIP_Real SCIPbendersGetSubproblemObjval(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** returns the number of cuts that have been added for storage */
EXTERN
int SCIPbendersGetNStoredCuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition cut */
   );

/** returns the data for the cuts that have been added by the Benders' cut plugin */
EXTERN
SCIP_RETCODE SCIPbendersGetStoredCutData(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition cut */
   int                   cutidx,             /**< the index for the cut data that is requested */
   SCIP_VAR***           vars,               /**< the variables that have non-zero coefficients in the cut */
   SCIP_Real**           vals,               /**< the coefficients of the variables in the cut */
   SCIP_Real*            lhs,                /**< the left hand side of the cut */
   SCIP_Real*            rhs,                /**< the right hand side of the cut */
   int*                  nvars               /**< the number of variables with non-zero coefficients in the cut */
   );

/** returns the original problem data for the cuts that have been added by the Benders' cut plugin. The stored
 *  variables and values will populate the input vars and vals arrays. Thus, memory must be allocated for the vars and
 *  vals arrays
 */
EXTERN
SCIP_RETCODE SCIPbendersGetStoredCutOrigData(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition cut */
   int                   cutidx,             /**< the index for the cut data that is requested */
   SCIP_VAR***           vars,               /**< the variables that have non-zero coefficients in the cut */
   SCIP_Real**           vals,               /**< the coefficients of the variables in the cut */
   SCIP_Real*            lhs,                /**< the left hand side of the cut */
   SCIP_Real*            rhs,                /**< the right hand side of the cut */
   int*                  nvars,              /**< the number of variables with non-zero coefficients in the cut */
   int                   varssize            /**< the available slots in the array */
   );

/*
 * Public functions associated with Benders' cuts
 */

/** returns the Benders' cut of the given name, or NULL if not existing */
EXTERN
SCIP_BENDERSCUT* SCIPfindBenderscut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   const char*           name                /**< name of Benderscut' decomposition */
   );


/** returns the array of currently available Benders' cuts; active Benders' decomposition are in the first slots of
 * the array
 */
EXTERN
SCIP_BENDERSCUT** SCIPbendersGetBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );


/** returns the number of currently available Benders' cuts */
EXTERN
int SCIPbendersGetNBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** sets the priority of a Benders' decomposition */
EXTERN
SCIP_RETCODE SCIPbendersSetBenderscutPriority(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' cut */
   int                   priority            /**< new priority of the Benders' decomposition */
   );

/** sets the flag indicating whether a subproblem is convex
 *
 *  It is possible that this can change during the solving process. One example is when the three-phase method is
 *  employed, where the first phase solves the convex relaxation of both the master and subproblems, the second phase
 *  reintroduces the integrality constraints to the master problem and the third phase then reintroduces integrality
 *  constraints to the subproblems.
 */
EXTERN
void SCIPbendersSetSubproblemIsConvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             isconvex            /**< flag to indicate whether the subproblem is convex */
   );

/** returns whether the subproblem is convex
 *
 *  This means that the dual solution can be used to generate cuts.
 */
EXTERN
SCIP_Bool SCIPbendersSubproblemIsConvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** returns the number of subproblems that are convex */
extern
int SCIPbendersGetNConvexSubproblems(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** solves the LP of the Benders' decomposition subproblem
 *
 *  This requires that the subproblem is in probing mode.
 */
EXTERN
SCIP_RETCODE SCIPbendersSolveSubproblemLP(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_STATUS*          solvestatus,        /**< status of subproblem solve */
   SCIP_Real*            objective           /**< optimal value of subproblem, if solved to optimality */
   );

/** solves the Benders' decomposition subproblem */
EXTERN
SCIP_RETCODE SCIPbendersSolveSubproblemCIP(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_STATUS*          solvestatus,        /**< status of subproblem solve */
   SCIP_Bool             solvecip            /**< directly solve the CIP subproblem */
   );

/** returns the number of cuts that have been transferred from sub SCIPs to the master SCIP */
EXTERN
int SCIPbendersGetNTransferredCuts(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition data structure */
   );

/** updates the lower bound for the subproblem. If the lower bound is not greater than the previously stored lowerbound,
 * then no update occurs.
 */
EXTERN
void SCIPbendersUpdateSubproblemLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Real             lowerbound          /**< the lower bound */
   );

/** returns the stored lower bound for the given subproblem */
EXTERN
SCIP_Real SCIPbendersGetSubproblemLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** sets the independent subproblem flag */
EXTERN
void SCIPbendersSetSubproblemIsIndependent(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             isindep             /**< flag to indicate whether the subproblem is independent */
   );

/** returns whether the subproblem is independent */
EXTERN
SCIP_Bool SCIPbendersSubproblemIsIndependent(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** returns whether the subproblem is enabled, i.e. the subproblem is still solved in the solving loop. */
EXTERN
SCIP_Bool SCIPbendersSubproblemIsEnabled(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
