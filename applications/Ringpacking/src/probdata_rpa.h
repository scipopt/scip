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

/**@file   probdata_rpa.h
 * @brief  Problem data for ringpacking problem
 * @author Benjamin Mueller
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_RINGPACKING__
#define __SCIP_PROBDATA_RINGPACKING__

#include "scip/scip.h"
#include "pattern.h"

/** sets up the problem data */
extern
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   int*                  demands,            /**< array containing the demands */
   SCIP_Real*            rints,              /**< internal radii of each ring */
   SCIP_Real*            rexts,              /**< external radii of each ring (assumed to be sorted) */
   int                   nitems,             /**< number of items */
   SCIP_Real             width,              /**< width of each rectangle */
   SCIP_Real             height              /**< height of each rectangle */
   );

/** enumerates circular patterns and creates restricted master problem */
extern
SCIP_RETCODE SCIPprobdataSetupProblem(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** enumerate all non-dominated circular patterns */
extern
SCIP_RETCODE SCIPprobdataEnumeratePatterns(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Real             nlptilim,           /**< time limit for each NLP verification */
   SCIP_Real             heurtilim,          /**< time limit for each call of the heuristics */
   SCIP_Real             totaltilim,         /**< total time limit for enumeration */
   SCIP_Longint          nlpnodelim,         /**< node limit for each NLP verification */
   int                   heuriterlim         /**< iteration limit for each call of the heuristics */
   );

/** returns number of different types */
extern
int SCIPprobdataGetNTypes(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns all external radii */
extern
SCIP_Real* SCIPprobdataGetRexts(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns all internal radii */
extern
SCIP_Real* SCIPprobdataGetRints(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns all demands */
extern
int* SCIPprobdataGetDemands(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns the width of each rectangle */
extern
SCIP_Real SCIPprobdataGetWidth(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns the height of each rectangle */
extern
SCIP_Real SCIPprobdataGetHeight(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns all information about circular patterns */
extern
void SCIPprobdataGetCInfos(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN***       cpatterns,          /**< pointer to store the circular patterns (might be NULL) */
   SCIP_VAR***           cvars,              /**< pointer to store the variables corresponding circular patterns (might be NULL) */
   int*                  ncpatterns          /**< pointer to store the number of circular patterns (might be NULL) */
   );

/** returns all information about rectangular patterns */
extern
void SCIPprobdataGetRInfos(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN***       rpatterns,          /**< pointer to store the rectangular patterns (might be NULL) */
   SCIP_VAR***           rvars,              /**< pointer to store the variables corresponding rectangular patterns (might be NULL) */
   int*                  nrpatterns          /**< pointer to store the number of rectangular patterns (might be NULL) */
   );

/** returns array of set pattern constraints */
extern
SCIP_CONS** SCIPprobdataGetPatternConss(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** adds given variable to the problem data */
extern
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN*         pattern,            /**< pattern */
   SCIP_VAR*             var                 /**< variables to add */
   );

/** updates the dual bound */
extern
void SCIPprobdataUpdateDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Real             dualbound           /**< new dual bound */
   );

/** marks that further reported dual bounds are not valid */
extern
void SCIPprobdataInvalidateDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns whether dual bound is marked to be invalid */
extern
SCIP_Bool SCIPprobdataIsDualboundInvalid(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** Tries to pack a list of elements into a specified boundary circle by using a simple left-first bottom-second
 *  heuristic. Returns the number of elements that could be stored and indicated which ones these are in the buffer
 *  parameter ispacked. This auxiliary method can be used both to find such a packing or to verify a certain pattern.
 */
extern
void SCIPpackCirclesGreedy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            rexts,              /**< outer radii of elements (in original order of probdata) */
   SCIP_Real*            xs,                 /**< buffer to store the resulting x-coordinates */
   SCIP_Real*            ys,                 /**< buffer to store the resulting y-coordinates */
   SCIP_Real             rbounding,          /**< inner radius of bounding circle (ignored for rectangular patterns) */
   SCIP_Real             width,              /**< width of the rectangle */
   SCIP_Real             height,             /**< height of the rectangle */
   SCIP_Bool*            ispacked,           /**< buffer to store which elements could be packed */
   int*                  elements,           /**< the order of the elements in the pattern */
   int                   nelements,          /**< number of elements in the pattern */
   SCIP_PATTERNTYPE      patterntype,        /**< the pattern type (rectangular or circular) */
   int*                  npacked,            /**< pointer to store the number of packed elements */
   int                   ncalls              /**< total number of calls of the packing heuristic */
   );

/** verifies a circular pattern heuristically */
extern
SCIP_RETCODE SCIPverifyCircularPatternHeuristic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN*         pattern,            /**< pattern */
   SCIP_Real             timelim,            /**< time limit */
   int                   iterlim             /**< iteration limit */
   );

/** verifies a circular pattern via solving a verification NLP */
extern
SCIP_RETCODE SCIPverifyCircularPatternNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN*         pattern,            /**< pattern */
   SCIP_Real             timelim,            /**< time limit */
   SCIP_Longint          nodelim             /**< node limit */
   );

/** check whether a pattern for consistency */
extern
void SCIPcheckPattern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN*         pattern             /**< pattern */
   );

#endif /* __SCIP_PROBDATA_RINGPACKING__ */
