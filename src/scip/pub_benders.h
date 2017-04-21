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

/** gets user data of variable benders */
EXTERN
SCIP_BENDERSDATA* SCIPbendersGetData(
   SCIP_BENDERS*         benders             /**< variable benders */
   );

/** sets user data of variable benders; user has to free old data in advance! */
EXTERN
void SCIPbendersSetData(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP_BENDERSDATA*     bendersdata         /**< new variable benders user data */
   );

/** gets name of variable benders */
EXTERN
const char* SCIPbendersGetName(
   SCIP_BENDERS*         benders             /**< variable benders */
   );

/** gets description of variable benders */
EXTERN
const char* SCIPbendersGetDesc(
   SCIP_BENDERS*         benders             /**< variable benders */
   );

/** gets priority of variable benders */
EXTERN
int SCIPbendersGetPriority(
   SCIP_BENDERS*         benders             /**< variable benders */
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

/** gets the number of times, the benders was called and tried to find a variable with negative reduced costs */
EXTERN
int SCIPbendersGetNCalls(
   SCIP_BENDERS*         benders             /**< variable benders */
   );

/** gets the number of optimality cuts found by the collection of Benders' decomposition subproblems */
EXTERN
int SCIPbendersGetNOptCutsFound(
   SCIP_BENDERS*         benders             /**< variable benders */
   );

/** gets the number of feasibility cuts found by the collection of Benders' decomposition subproblems */
EXTERN
int SCIPbendersGetNFeasCutsFound(
   SCIP_BENDERS*         benders             /**< variable benders */
   );

/** gets time in seconds used in this benders for setting up for next stages */
EXTERN
SCIP_Real SCIPbendersGetSetupTime(
   SCIP_BENDERS*         benders             /**< variable benders */
   );

/** gets time in seconds used in this benders */
EXTERN
SCIP_Real SCIPbendersGetTime(
   SCIP_BENDERS*         benders             /**< variable benders */
   );

/** is variable benders initialized? */
EXTERN
SCIP_Bool SCIPbendersIsInitialized(
   SCIP_BENDERS*         benders             /**< variable benders */
   );

/** are Benders' cuts generated from the LP solutions? */
EXTERN
SCIP_Bool SCIPbendersCutLP(
   SCIP_BENDERS*         benders             /**< variable benders */
   );

/** are Benders' cuts generated from the pseudo solutions? */
EXTERN
SCIP_Bool SCIPbendersCutPseudo(
   SCIP_BENDERS*         benders             /**< variable benders */
   );

/** are Benders' cuts generated from the relaxation solutions? */
EXTERN
SCIP_Bool SCIPbendersCutRelaxation(
   SCIP_BENDERS*         benders             /**< variable benders */
   );

/** Adds a subproblem to the Benders' decomposition data */
EXTERN
void SCIPbendersAddSubproblem(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP*                 subproblem          /**< subproblem to be added to the data storage */
   );

/** returns the auxiliary variable for the given subproblem */
EXTERN
SCIP_VAR* SCIPbendersGetAuxiliaryVar(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** stores the objective function value of the subproblem for use in cut generation */
EXTERN
void SCIPbendersSetSubprobObjval(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP_Real             objval,             /**< the objective function value for the subproblem */
   int                   probnumber          /**< the subproblem number */
   );

/** returns the objective function value of the subproblem for use in cut generation */
EXTERN
SCIP_Real SCIPbendersGetSubprobObjval(
   SCIP_BENDERS*         benders,            /**< variable benders */
   int                   probnumber          /**< the subproblem number */
   );

/** Public functions associated with Benders' cuts */

/** inserts a Benders' cut into the Benders' cuts list */
EXTERN
SCIP_RETCODE SCIPbendersIncludeBenderscut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BENDERSCUT*      benderscut          /**< Benders' cut */
   );

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

/* @} */

#ifdef __cplusplus
}
#endif

#endif
