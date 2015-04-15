/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_compress.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for tree compressions
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_COMPR_H__
#define __SCIP_PUB_COMPR_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_compress.h"

#ifdef __cplusplus
extern "C" {
#endif

/** compares two compressions w. r. to their priority */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPcomprComp);

/** comparison method for sorting compressions w.r.t. to their name */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPcomprCompName);

/** gets user data of tree compression */
EXTERN
SCIP_COMPRDATA* SCIPcomprGetData(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** sets user data of tree compression; user has to free old data in advance! */
EXTERN
void SCIPcomprSetData(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_COMPRDATA*       comprdata           /**< new tree compression user data */
   );

/** gets name of tree compression */
EXTERN
const char* SCIPcomprGetName(
   SCIP_COMPR*            heur               /**< tree compression */
   );

/** gets description of tree compression */
EXTERN
const char* SCIPcomprGetDesc(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** gets display character of tree compression */
EXTERN
char SCIPcomprGetDispchar(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** does the compression use a secondary SCIP instance? */
EXTERN
SCIP_Bool SCIPcomprUsesSubscip(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** gets priority of tree compression */
EXTERN
int SCIPcomprGetPriority(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** gets minimal depth level for calling tree compression (returns -1, if no depth treshold exists) */
EXTERN
int SCIPcomprGetMindepth(
   SCIP_COMPR*          compr               /**< tree compression */
   );

/** gets minimal minimal number of nodes for calling tree compression (returns -1, if no node threshold exists) */
EXTERN
int SCIPcomprGetMinnnodes(
   SCIP_COMPR*          compr               /**< tree compression */
   );

/** gets the number of times, the compression was called and tried to find a compression */
EXTERN
SCIP_Longint SCIPcomprGetNCalls(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** gets the number of tree compressions found by this compression */
EXTERN
SCIP_Longint SCIPcomprGetNCompressionFound(
   SCIP_COMPR*            compr              /**< tree compression */
   );


/* returns the compression rate of the last compression */
EXTERN
SCIP_Real SCIPcomprGetRate(
   SCIP_COMPR*           compr
   );

/* set the rate of the compression */
EXTERN
void SCIPcomprUpdateRate(
   SCIP_COMPR*           compr,
   SCIP_Real             rate
   );

/* return the number of nodes of the last compression */
EXTERN
int SCIPcomprGetNNodes(
   SCIP_COMPR*           compr
   );

/* update the number of nodes of the last compression */
EXTERN
void SCIPcomprUpdateNNodes(
   SCIP_COMPR*           compr,
   int                   nnodes
   );

/* returns the loss of information of the last compression */
EXTERN
SCIP_Real SCIPcomprGetLOI(
   SCIP_COMPR*           compr
   );

/* update the loss of information of the last compression */
EXTERN
void SCIPcomprUpdateLOI(
   SCIP_COMPR*           compr,
   SCIP_Real             loi
   );

/** is tree compression initialized? */
EXTERN
SCIP_Bool SCIPcomprIsInitialized(
   SCIP_COMPR*           compr                /**< tree compression */
   );

/** gets time in seconds used in this compression for setting up for next stages */
EXTERN
SCIP_Real SCIPcomprGetSetupTime(
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** gets time in seconds used in this compression */
EXTERN
SCIP_Real SCIPcomprGetTime(
   SCIP_COMPR*           compr                /**< tree compression */
   );

#ifdef __cplusplus
}
#endif

#endif
