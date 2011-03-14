/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_prop.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for propagators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_PROP_H__
#define __SCIP_PUB_PROP_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_prop.h"

#ifdef __cplusplus
extern "C" {
#endif

/** compares two propagators w. r. to their priority */
extern
SCIP_DECL_SORTPTRCOMP(SCIPpropComp);

/** gets user data of propagator */
extern
SCIP_PROPDATA* SCIPpropGetData(
   SCIP_PROP*            prop                /**< propagator */
   );

/** sets user data of propagator; user has to free old data in advance! */
extern
void SCIPpropSetData(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROPDATA*        propdata            /**< new propagator user data */
   );

/** gets name of propagator */
extern
const char* SCIPpropGetName(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets description of propagator */
extern
const char* SCIPpropGetDesc(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets priority of propagator */
extern
int SCIPpropGetPriority(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets frequency of propagator */
extern
int SCIPpropGetFreq(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used in this propagator */
extern
SCIP_Real SCIPpropGetTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used in this propagator for resolve propagation */
extern
SCIP_Real SCIPpropGetRespropTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets the total number of times, the propagator was called */
extern
SCIP_Longint SCIPpropGetNCalls(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets the total number of times, the propagator was called for resolving a propagation */
extern
SCIP_Longint SCIPpropGetNRespropCalls(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets total number of times, this propagator detected a cutoff */
extern
SCIP_Longint SCIPpropGetNCutoffs(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets total number of domain reductions found by this propagator */
extern
SCIP_Longint SCIPpropGetNDomredsFound(
   SCIP_PROP*            prop                /**< propagator */
   );

/** should propagator be delayed, if other propagators found reductions? */
extern
SCIP_Bool SCIPpropIsDelayed(
   SCIP_PROP*            prop                /**< propagator */
   );

/** was propagator delayed at the last call? */
extern
SCIP_Bool SCIPpropWasDelayed(
   SCIP_PROP*            prop                /**< propagator */
   );

/** is propagator initialized? */
extern
SCIP_Bool SCIPpropIsInitialized(
   SCIP_PROP*            prop                /**< propagator */
   );

#ifdef __cplusplus
}
#endif

#endif
