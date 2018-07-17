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

/**@file   pub_benderscut.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for Benders' decomposition cuts
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_BENDERSCUT_H__
#define __SCIP_PUB_BENDERSCUT_H__

#include "scip/def.h"
#include "scip/type_benderscut.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_misc.h"
#include "scip/type_retcode.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicBenderscutsMethods
 *
 * @{
 */

/** compares two Benders' decomposition cuts w. r. to their priority */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPbenderscutComp);

/** comparison method for sorting Benders' decomposition cuts w.r.t. to their name */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPbenderscutCompName);

/** gets user data of the Benders' decomposition cut */
EXTERN
SCIP_BENDERSCUTDATA* SCIPbenderscutGetData(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** sets user data of the Benders' decomposition cut; user has to free old data in advance! */
EXTERN
void SCIPbenderscutSetData(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_BENDERSCUTDATA*  benderscutdata      /**< new Benders' decomposition cut user data */
   );

/** gets name of the Benders' decomposition cut */
EXTERN
const char* SCIPbenderscutGetName(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** gets description of the Benders' decomposition cut */
EXTERN
const char* SCIPbenderscutGetDesc(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** gets priority of the Benders' decomposition cut */
EXTERN
int SCIPbenderscutGetPriority(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** gets the number of times, the Benders' decomposition cut was called and tried to find a violated cut */
EXTERN
SCIP_Longint SCIPbenderscutGetNCalls(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** gets the number of the cuts found by this Benders' decomposition cut */
EXTERN
SCIP_Longint SCIPbenderscutGetNFound(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** is the Benders' decomposition cut initialized? */
EXTERN
SCIP_Bool SCIPbenderscutIsInitialized(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** gets time in seconds used in this Benders' decomposition cut for setting up for next stages */
EXTERN
SCIP_Real SCIPbenderscutGetSetupTime(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** gets time in seconds used in this Benders' decomposition cut */
EXTERN
SCIP_Real SCIPbenderscutGetTime(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** returns the constraints that have been added by the Benders' cut plugin */
EXTERN
SCIP_RETCODE SCIPbenderscutGetAddedConss(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_CONS***          addedconss,         /**< pointer to store the constraint array, can be NULL */
   int*                  naddedconss         /**< pointer to store the number of added constraints, can be NULL */
   );

/** returns the cuts that have been added by the Benders' cut plugin */
EXTERN
SCIP_RETCODE SCIPbenderscutGetAddedCuts(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_ROW***           addedcuts,          /**< pointer to store the cuts array, can be NULL */
   int*                  naddedcuts          /**< pointer to store the number of added cut, can be NULL */
   );

/** returns whether the Benders' cut uses the LP information */
EXTERN
SCIP_Bool SCIPbenderscutIsLPCut(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** sets the enabled flag of the Benders' decomposition cut method */
EXTERN
void SCIPbenderscutSetEnabled(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_Bool             enabled             /**< flag to indicate whether the Benders' decomposition cut is enabled */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
