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

/**@file   pub_sepa.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_SEPA_H__
#define __SCIP_PUB_SEPA_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_sepa.h"

#ifdef __cplusplus
extern "C" {
#endif

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@addtogroup PublicSeparatorMethods
 *
 * @{
 */


/** compares two separators w. r. to their priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPsepaComp);

/** comparison method for sorting separators w.r.t. to their name */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPsepaCompName);

/** gets user data of separator */
SCIP_EXPORT
SCIP_SEPADATA* SCIPsepaGetData(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** sets user data of separator; user has to free old data in advance! */
SCIP_EXPORT
void SCIPsepaSetData(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata            /**< new separator user data */
   );

/** gets name of separator */
SCIP_EXPORT
const char* SCIPsepaGetName(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets description of separator */
SCIP_EXPORT
const char* SCIPsepaGetDesc(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets priority of separator */
SCIP_EXPORT
int SCIPsepaGetPriority(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets frequency of separator */
SCIP_EXPORT
int SCIPsepaGetFreq(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** sets frequency of separator */
SCIP_EXPORT
void SCIPsepaSetFreq(
   SCIP_SEPA*            sepa,               /**< separator */
   int                   freq                /**< new frequency of separator */
   );

/** get maximal bound distance at which the separator is called */
SCIP_EXPORT
SCIP_Real SCIPsepaGetMaxbounddist(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** does the separator use a secondary SCIP instance? */
SCIP_EXPORT
SCIP_Bool SCIPsepaUsesSubscip(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets time in seconds used in this separator for setting up for next stages */
SCIP_EXPORT
SCIP_Real SCIPsepaGetSetupTime(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets time in seconds used in this separator */
SCIP_EXPORT
SCIP_Real SCIPsepaGetTime(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the total number of times, the separator was called */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCalls(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the number of times, the separator was called at the current node */
SCIP_EXPORT
int SCIPsepaGetNCallsAtNode(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets total number of times, the separator detected a cutoff */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutoffs(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the total number of cutting planes found by this separator */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutsFound(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the total number of cutting planes applied to lp */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutsApplied(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the number of cutting planes found by this separator at the current node */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutsFoundAtNode(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets total number of additional constraints added by this separator */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNConssFound(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets total number of domain reductions found by this separator */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNDomredsFound(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** should separator be delayed, if other separators found cuts? */
SCIP_EXPORT
SCIP_Bool SCIPsepaIsDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** was separation of the LP solution delayed at the last call? */
SCIP_EXPORT
SCIP_Bool SCIPsepaWasLPDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** was separation of the primal solution delayed at the last call? */
SCIP_EXPORT
SCIP_Bool SCIPsepaWasSolDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** is separator initialized? */
SCIP_EXPORT
SCIP_Bool SCIPsepaIsInitialized(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
