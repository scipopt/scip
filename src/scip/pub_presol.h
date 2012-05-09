/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_presol.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_PRESOL_H__
#define __SCIP_PUB_PRESOL_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_presol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** compares two presolvers w. r. to their priority */
extern
SCIP_DECL_SORTPTRCOMP(SCIPpresolComp);

/** comparison method for sorting presolvers w.r.t. to their name */
extern
SCIP_DECL_SORTPTRCOMP(SCIPpresolCompName);

/** gets user data of presolver */
extern
SCIP_PRESOLDATA* SCIPpresolGetData(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** sets user data of presolver; user has to free old data in advance! */
extern
void SCIPpresolSetData(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_PRESOLDATA*      presoldata          /**< new presolver user data */
   );

/** gets name of presolver */
extern
const char* SCIPpresolGetName(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets description of presolver */
extern
const char* SCIPpresolGetDesc(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets priority of presolver */
extern
int SCIPpresolGetPriority(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** should presolver be delayed, if other presolvers found reductions? */
extern
SCIP_Bool SCIPpresolIsDelayed(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** was presolver delayed at the last call? */
extern
SCIP_Bool SCIPpresolWasDelayed(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** is presolver initialized? */
extern
SCIP_Bool SCIPpresolIsInitialized(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets time in seconds used in this presolver for setting up for next stages */
extern
SCIP_Real SCIPpresolGetSetupTime(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets time in seconds used in this presolver */
extern
SCIP_Real SCIPpresolGetTime(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of variables fixed in presolver */
extern
int SCIPpresolGetNFixedVars(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of variables aggregated in presolver */
extern
int SCIPpresolGetNAggrVars(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of variable types changed in presolver */
extern
int SCIPpresolGetNChgVarTypes(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of bounds changed in presolver */
extern
int SCIPpresolGetNChgBds(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of holes added to domains of variables in presolver */
extern
int SCIPpresolGetNAddHoles(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of constraints deleted in presolver */
extern
int SCIPpresolGetNDelConss(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of constraints added in presolver */
extern
int SCIPpresolGetNAddConss(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of constraints upgraded in presolver */
extern
int SCIPpresolGetNUpgdConss(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of coefficients changed in presolver */
extern
int SCIPpresolGetNChgCoefs(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of constraint sides changed in presolver */
extern
int SCIPpresolGetNChgSides(
   SCIP_PRESOL*          presol              /**< presolver */
   );

#ifdef __cplusplus
}
#endif

#endif
