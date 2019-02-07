/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_prop.h
 * @ingroup PUBLICCOREAPI
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

/**@addtogroup PublicPropagatorMethods
 *
 * @{
 */

/** compares two propagators w. r. to their priority */
SCIP_EXPORT extern
SCIP_DECL_SORTPTRCOMP(SCIPpropComp);

/** compares two propagators w. r. to their presolving priority */
SCIP_EXPORT extern
SCIP_DECL_SORTPTRCOMP(SCIPpropCompPresol);

/** comparison method for sorting propagators w.r.t. to their name */
SCIP_EXPORT extern
SCIP_DECL_SORTPTRCOMP(SCIPpropCompName);

/** gets user data of propagator */
SCIP_EXPORT extern
SCIP_PROPDATA* SCIPpropGetData(
   SCIP_PROP*            prop                /**< propagator */
   );

/** sets user data of propagator; user has to free old data in advance! */
SCIP_EXPORT extern
void SCIPpropSetData(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROPDATA*        propdata            /**< new propagator user data */
   );

/** gets name of propagator */
SCIP_EXPORT extern
const char* SCIPpropGetName(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets description of propagator */
SCIP_EXPORT extern
const char* SCIPpropGetDesc(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets priority of propagator */
SCIP_EXPORT extern
int SCIPpropGetPriority(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets presolving priority of propagator */
SCIP_EXPORT extern
int SCIPpropGetPresolPriority(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets frequency of propagator */
SCIP_EXPORT extern
int SCIPpropGetFreq(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used for setting up this propagator for new stages */
SCIP_EXPORT extern
SCIP_Real SCIPpropGetSetupTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** sets frequency of propagator */
SCIP_EXPORT extern
void SCIPpropSetFreq(
   SCIP_PROP*            prop,               /**< propagator */
   int                   freq                /**< new frequency of propagator */
   );

/** gets time in seconds used in this propagator */
SCIP_EXPORT extern
SCIP_Real SCIPpropGetTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used in this propagator during strong branching */
SCIP_EXPORT extern
SCIP_Real SCIPpropGetStrongBranchPropTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used in this propagator for resolve propagation */
SCIP_EXPORT extern
SCIP_Real SCIPpropGetRespropTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used in this propagator for presolving */
SCIP_EXPORT extern
SCIP_Real SCIPpropGetPresolTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets the total number of times, the propagator was called */
SCIP_EXPORT extern
SCIP_Longint SCIPpropGetNCalls(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets the total number of times, the propagator was called for resolving a propagation */
SCIP_EXPORT extern
SCIP_Longint SCIPpropGetNRespropCalls(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets total number of times, this propagator detected a cutoff */
SCIP_EXPORT extern
SCIP_Longint SCIPpropGetNCutoffs(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets total number of domain reductions found by this propagator */
SCIP_EXPORT extern
SCIP_Longint SCIPpropGetNDomredsFound(
   SCIP_PROP*            prop                /**< propagator */
   );

/** should propagator be delayed, if other propagators found reductions? */
SCIP_EXPORT extern
SCIP_Bool SCIPpropIsDelayed(
   SCIP_PROP*            prop                /**< propagator */
   );

/** was propagator delayed at the last call? */
SCIP_EXPORT extern
SCIP_Bool SCIPpropWasDelayed(
   SCIP_PROP*            prop                /**< propagator */
   );

/** is propagator initialized? */
SCIP_EXPORT extern
SCIP_Bool SCIPpropIsInitialized(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of variables fixed during presolving of propagator */
SCIP_EXPORT extern
int SCIPpropGetNFixedVars(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of variables aggregated during presolving of propagator  */
SCIP_EXPORT extern
int SCIPpropGetNAggrVars(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of variable types changed during presolving of propagator  */
SCIP_EXPORT extern
int SCIPpropGetNChgVarTypes(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of bounds changed during presolving of propagator  */
SCIP_EXPORT extern
int SCIPpropGetNChgBds(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of holes added to domains of variables during presolving of propagator  */
SCIP_EXPORT extern
int SCIPpropGetNAddHoles(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of constraints deleted during presolving of propagator */
SCIP_EXPORT extern
int SCIPpropGetNDelConss(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of constraints added during presolving of propagator */
SCIP_EXPORT extern
int SCIPpropGetNAddConss(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of constraints upgraded during presolving of propagator  */
SCIP_EXPORT extern
int SCIPpropGetNUpgdConss(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of coefficients changed during presolving of propagator */
SCIP_EXPORT extern
int SCIPpropGetNChgCoefs(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of constraint sides changed during presolving of propagator */
SCIP_EXPORT extern
int SCIPpropGetNChgSides(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of times the propagator was called in presolving and tried to find reductions */
SCIP_EXPORT extern
int SCIPpropGetNPresolCalls(
   SCIP_PROP*            prop                /**< propagator */
   );

/** returns the timing mask of the propagator */
SCIP_EXPORT extern
SCIP_PROPTIMING SCIPpropGetTimingmask(
   SCIP_PROP*            prop                /**< propagator */
   );

/** does the propagator perform presolving? */
SCIP_EXPORT extern
SCIP_Bool SCIPpropDoesPresolve(
   SCIP_PROP*            prop                /**< propagator */
   );

/** returns the timing mask of the presolving method of the propagator */
SCIP_EXPORT extern
SCIP_PRESOLTIMING SCIPpropGetPresolTiming(
   SCIP_PROP*            prop                /**< propagator */
   );

/** sets the timing mask of the presolving method of the propagator */
SCIP_EXPORT extern
void SCIPpropSetPresolTiming(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PRESOLTIMING     presoltiming        /** timing mask to be set */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
