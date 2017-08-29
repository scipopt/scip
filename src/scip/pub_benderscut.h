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

/**@file   pub_benderscut.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for Benders' decomposition cutss
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_BENDERSCUT_H__
#define __SCIP_PUB_BENDERSCUT_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_benderscut.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicBenderscutsMethods
 *
 * @{
 */

/** compares two compressions w. r. to their priority */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPbenderscutComp);

/** comparison method for sorting compressions w.r.t. to their name */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPbenderscutCompName);

/** gets user data of Benders' decomposition cuts */
EXTERN
SCIP_BENDERSCUTDATA* SCIPbenderscutGetData(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cuts */
   );

/** sets user data of Benders' decomposition cuts; user has to free old data in advance! */
EXTERN
void SCIPbenderscutSetData(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_BENDERSCUTDATA*   benderscutdata     /**< new Benders' decomposition cuts user data */
   );

/** gets name of Benders' decomposition cuts */
EXTERN
const char* SCIPbenderscutGetName(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cuts */
   );

/** gets description of Benders' decomposition cuts */
EXTERN
const char* SCIPbenderscutGetDesc(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cuts */
   );

/** gets priority of Benders' decomposition cuts */
EXTERN
int SCIPbenderscutGetPriority(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cuts */
   );

/** gets the number of times, the compression was called and tried to find a compression */
EXTERN
SCIP_Longint SCIPbenderscutGetNCalls(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cuts */
   );

/** gets the number of Benders' decomposition cuts found for a given subproblem */
EXTERN
SCIP_Longint SCIPbenderscutGetNFound(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cuts */
   );

/** is Benders' decomposition cuts initialized? */
EXTERN
SCIP_Bool SCIPbenderscutIsInitialized(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cuts */
   );

/** gets time in seconds used in this compression for setting up for next stages */
EXTERN
SCIP_Real SCIPbenderscutGetSetupTime(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cuts */
   );

/** gets time in seconds used in this compression */
EXTERN
SCIP_Real SCIPbenderscutGetTime(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cuts */
   );

/** returns the constraints that have been added by the Benders' cut plugin */
EXTERN
SCIP_RETCODE SCIPbenderscutGetCons(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_CONS***          addedcons,          /**< pointer to store the constraint array */
   int*                  naddedcons          /**< pointer to store the number of added constraints */
   );

/** returns the cuts that have been added by the Benders' cut plugin */
EXTERN
SCIP_RETCODE SCIPbenderscutGetCuts(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_ROW***           addedcuts,          /**< pointer to store the cuts array */
   int*                  naddedcuts          /**< pointer to store the number of added cut */
   );

/** returns the number of constraints that have been added by the Benders' cut plugin */
EXTERN
int SCIPbenderscutGetNAddedCons(
   SCIP_BENDERSCUT*      benderscut         /**< Benders' decomposition cut */
   );

/** returns the number of cuts that have been added by the Benders' cut plugin */
EXTERN
int SCIPbenderscutGetNAddedCuts(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** returns whether the Benders' cut uses the LP information */
EXTERN
SCIP_Bool SCIPbenderscutIsLPCut(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
