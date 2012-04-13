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

/**@file   pub_heur.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for primal heuristics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_HEUR_H__
#define __SCIP_PUB_HEUR_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_heur.h"

#ifdef __cplusplus
extern "C" {
#endif

/** compares two heuristics w. r. to their priority */
extern
SCIP_DECL_SORTPTRCOMP(SCIPheurComp);

/** comparison method for sorting heuristics w.r.t. to their name */
extern
SCIP_DECL_SORTPTRCOMP(SCIPheurCompName);

/** gets user data of primal heuristic */
extern
SCIP_HEURDATA* SCIPheurGetData(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** sets user data of primal heuristic; user has to free old data in advance! */
extern
void SCIPheurSetData(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_HEURDATA*        heurdata            /**< new primal heuristic user data */
   );

/** sets copy method of primal heuristic */
extern
void SCIPheurSetCopy(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURCOPY    ((*heurcopy))       /**< copy method of primal heuristic or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of primal heuristic */
extern
void SCIPheurSetFree(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURFREE    ((*heurfree))       /**< destructor of primal heuristic */
   );

/** sets initialization method of primal heuristic */
extern
void SCIPheurSetInit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURINIT    ((*heurinit))       /**< initialize primal heuristic */
   );

/** sets deinitialization method of primal heuristic */
extern
void SCIPheurSetExit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEUREXIT    ((*heurexit))       /**< deinitialize primal heuristic */
   );

/** sets solving process initialization method of primal heuristic */
extern
void SCIPheurSetInitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURINITSOL ((*heurinitsol))    /**< solving process initialization method of primal heuristic */
   );

/** sets solving process deinitialization method of primal heuristic */
extern
void SCIPheurSetExitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEUREXITSOL ((*heurexitsol))    /**< solving process deinitialization method of primal heuristic */
   );

/** gets name of primal heuristic */
extern
const char* SCIPheurGetName(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets description of primal heuristic */
extern
const char* SCIPheurGetDesc(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets display character of primal heuristic */
extern
char SCIPheurGetDispchar(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** returns the timing mask of the heuristic */
extern
SCIP_HEURTIMING SCIPheurGetTimingmask(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** sets new timing mask for heuristic */
extern
void SCIPheurSetTimingmask(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_HEURTIMING       timingmask          /**< new timing mask of heuristic */
   );

/** does the heuristic use a secondary SCIP instance? */
extern
SCIP_Bool SCIPheurUsesSubscip(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets priority of primal heuristic */
extern
int SCIPheurGetPriority(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets frequency of primal heuristic */
extern
int SCIPheurGetFreq(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** sets frequency of primal heuristic */
extern
void SCIPheurSetFreq(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   int                   freq                /**< new frequency of heuristic */
   );

/** gets frequency offset of primal heuristic */
extern
int SCIPheurGetFreqofs(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets maximal depth level for calling primal heuristic (returns -1, if no depth limit exists) */
extern
int SCIPheurGetMaxdepth(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets the number of times, the heuristic was called and tried to find a solution */
extern
SCIP_Longint SCIPheurGetNCalls(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets the number of primal feasible solutions found by this heuristic */
extern
SCIP_Longint SCIPheurGetNSolsFound(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets the number of new best primal feasible solutions found by this heuristic */
extern
SCIP_Longint SCIPheurGetNBestSolsFound(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** is primal heuristic initialized? */
extern
SCIP_Bool SCIPheurIsInitialized(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets time in seconds used in this heuristic for setting up for next stages */
extern
SCIP_Real SCIPheurGetSetupTime(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets time in seconds used in this heuristic */
extern
SCIP_Real SCIPheurGetTime(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

#ifdef __cplusplus
}
#endif

#endif
