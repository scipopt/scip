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

/**@file   pub_compr.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for tree compressions
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_COMPR_H__
#define __SCIP_PUB_COMPR_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_compr.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicCompressionMethods
 *
 * @{
 */

/** compares two compressions w. r. to their priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPcomprComp);

/** comparison method for sorting compressions w.r.t. to their name */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPcomprCompName);

/** gets user data of tree compression */
SCIP_EXPORT
SCIP_COMPRDATA* SCIPcomprGetData(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** sets user data of tree compression; user has to free old data in advance! */
SCIP_EXPORT
void SCIPcomprSetData(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_COMPRDATA*       comprdata           /**< new tree compression user data */
   );

/** gets name of tree compression */
SCIP_EXPORT
const char* SCIPcomprGetName(
   SCIP_COMPR*           heur                /**< tree compression */
   );

/** gets description of tree compression */
SCIP_EXPORT
const char* SCIPcomprGetDesc(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** gets priority of tree compression */
SCIP_EXPORT
int SCIPcomprGetPriority(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** gets minimal number of nodes for calling tree compression (returns -1, if no node threshold exists) */
SCIP_EXPORT
int SCIPcomprGetMinNodes(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** gets the number of times, the compression was called and tried to find a compression */
SCIP_EXPORT
SCIP_Longint SCIPcomprGetNCalls(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** gets the number of tree compressions found by this compression */
SCIP_EXPORT
SCIP_Longint SCIPcomprGetNFound(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** is tree compression initialized? */
SCIP_EXPORT
SCIP_Bool SCIPcomprIsInitialized(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** gets time in seconds used in this compression for setting up for next stages */
SCIP_EXPORT
SCIP_Real SCIPcomprGetSetupTime(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** gets time in seconds used in this compression */
SCIP_EXPORT
SCIP_Real SCIPcomprGetTime(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
