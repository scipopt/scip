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

/* @} */

#ifdef __cplusplus
}
#endif

#endif
