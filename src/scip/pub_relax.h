/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_relax.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for relaxation handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_RELAX_H__
#define __SCIP_PUB_RELAX_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_relax.h"

#ifdef __cplusplus
extern "C" {
#endif


/**@addtogroup PublicRelaxatorMethods
 *
 * @{
 */


/** compares two relaxation handlers w. r. to their priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPrelaxComp);

/** comparison method for sorting relaxators w.r.t. to their name */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPrelaxCompName);

/** gets user data of relaxation handler */
SCIP_EXPORT
SCIP_RELAXDATA* SCIPrelaxGetData(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** sets user data of relaxation handler; user has to free old data in advance! */
SCIP_EXPORT
void SCIPrelaxSetData(
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_RELAXDATA*       relaxdata           /**< new relaxation handler user data */
   );

/** gets name of relaxation handler */
SCIP_EXPORT
const char* SCIPrelaxGetName(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets description of relaxation handler */
SCIP_EXPORT
const char* SCIPrelaxGetDesc(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets priority of relaxation handler */
SCIP_EXPORT
int SCIPrelaxGetPriority(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets frequency of relaxation handler */
SCIP_EXPORT
int SCIPrelaxGetFreq(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets time in seconds used in this relaxator for setting up for next stages */
SCIP_EXPORT
SCIP_Real SCIPrelaxGetSetupTime(
   SCIP_RELAX*           relax               /**< relaxator */
   );

/** gets time in seconds used in this relaxation handler */
SCIP_EXPORT
SCIP_Real SCIPrelaxGetTime(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets the total number of times the relaxation handler was called */
SCIP_EXPORT
SCIP_Longint SCIPrelaxGetNCalls(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets the total number of times the relaxation handler cut off a node */
SCIP_EXPORT
SCIP_Longint SCIPrelaxGetNCutoffs(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets the total number of times the relaxation handler improved a node's lower bound */
SCIP_EXPORT
SCIP_Longint SCIPrelaxGetNImprovedLowerbound(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets the time in seconds spent for the execution of the relaxation handler when a node's lower bound could be improved (or a cutoff was found) */
SCIP_EXPORT
SCIP_Real SCIPrelaxGetImprovedLowerboundTime(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets the total number of times the relaxation handler added constraints */
SCIP_EXPORT
SCIP_Longint SCIPrelaxGetNAddedConss(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets the total number of times the relaxation handler reduced variable domains */
SCIP_EXPORT
SCIP_Longint SCIPrelaxGetNReducedDomains(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets the total number of times the relaxation handler separated cutting planes */
SCIP_EXPORT
SCIP_Longint SCIPrelaxGetNSeparatedCuts(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** is relaxation handler initialized? */
SCIP_EXPORT
SCIP_Bool SCIPrelaxIsInitialized(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** marks the current relaxation unsolved, s.t. the relaxation handler is called again in the next solving round */
SCIP_EXPORT
void SCIPrelaxMarkUnsolved(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
