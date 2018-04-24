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

/**@file   pub_benders.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_BENDERS_H__
#define __SCIP_PUB_BENDERS_H__


#include "scip/scip.h"
#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_benders.h"

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

/** is Benders' decomposition initialized? */
EXTERN
SCIP_Bool SCIPbendersIsInitialized(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** returns whether the given Benders decomposition is in use in the current problem */
SCIP_Bool SCIPbendersIsActive(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   );

/** are Benders' cuts generated from the LP solutions? */
EXTERN
SCIP_Bool SCIPbendersCutLP(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** are Benders' cuts generated from the pseudo solutions? */
EXTERN
SCIP_Bool SCIPbendersCutPseudo(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** are Benders' cuts generated from the relaxation solutions? */
EXTERN
SCIP_Bool SCIPbendersCutRelaxation(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** should this Benders' use the auxiliary variables from the highest priority Benders' */
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
void SCIPbendersSetSubprobObjval(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Real             objval              /**< the objective function value for the subproblem */
   );

/** returns the objective function value of the subproblem for use in cut generation */
EXTERN
SCIP_Real SCIPbendersGetSubprobObjval(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** Public functions associated with Benders' cuts */

/** returns the Benders' cut of the given name, or NULL if not existing */
EXTERN
SCIP_BENDERSCUT* SCIPfindBenderscut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   const char*           name                /**< name of Benderscut' decomposition */
   );


/** returns the array of currently available Benders' cuts; active benders are in the first slots of the array */
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

/* sets the flag indicating whether a subproblem is convex. It is possible that this can change during the solving
 * process. One example is when the three-phase method is employed, where the first phase solves the of both the master
 * and subproblems and by the third phase the integer subproblem is solved.
 */
EXTERN
void SCIPbendersSetSubprobIsConvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             isconvex            /**< flag to indicate whether the subproblem is convex */
   );

/** returns whether the subproblem is convex. This means that the dual solution can be used to generate cuts. */
EXTERN
SCIP_Bool SCIPbendersSubprobIsConvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** returns the number of subproblems that are convex */
extern
int SCIPbendersGetNConvexSubprobs(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** solves the LP of the Benders' decomposition subproblem. This requires that the subproblem is in probing mode */
EXTERN
SCIP_RETCODE SCIPbendersSolveSubproblemLP(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible          /**< a flag to indicate whether all subproblems are feasible */
   );

/** solves the Benders' decomposition subproblem. */
EXTERN
SCIP_RETCODE SCIPbendersSolveSubproblemMIP(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_BENDERSENFOTYPE  type,               /**< the enforcement type calling this function */
   SCIP_Bool             solvemip            /**< directly solve the MIP subproblem */
   );

/** returns the number of cuts that have been transferred from sub SCIPs to the master SCIP */
EXTERN
int SCIPbendersGetNTransferredCuts(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition data structure */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
