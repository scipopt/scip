/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_heur.h,v 1.7 2005/01/21 09:17:03 bzfpfend Exp $"

/**@file   pub_heur.h
 * @brief  public methods for primal heuristics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_HEUR_H__
#define __PUB_HEUR_H__


#include "def.h"
#include "type_misc.h"
#include "type_heur.h"



/** compares two heuristics w. r. to their priority */
extern
DECL_SORTPTRCOMP(SCIPheurComp);

/** gets user data of primal heuristic */
extern
HEURDATA* SCIPheurGetData(
   HEUR*            heur                /**< primal heuristic */
   );

/** sets user data of primal heuristic; user has to free old data in advance! */
extern
void SCIPheurSetData(
   HEUR*            heur,               /**< primal heuristic */
   HEURDATA*        heurdata            /**< new primal heuristic user data */
   );

/** gets name of primal heuristic */
extern
const char* SCIPheurGetName(
   HEUR*            heur                /**< primal heuristic */
   );

/** gets description of primal heuristic */
extern
const char* SCIPheurGetDesc(
   HEUR*            heur                /**< primal heuristic */
   );

/** gets display character of primal heuristic */
extern
char SCIPheurGetDispchar(
   HEUR*            heur                /**< primal heuristic */
   );

/** gets priority of primal heuristic */
extern
int SCIPheurGetPriority(
   HEUR*            heur                /**< primal heuristic */
   );

/** gets frequency of primal heuristic */
extern
int SCIPheurGetFreq(
   HEUR*            heur                /**< primal heuristic */
   );

/** gets frequency offset of primal heuristic */
extern
int SCIPheurGetFreqofs(
   HEUR*            heur                /**< primal heuristic */
   );

/** gets maximal depth level for calling primal heuristic (returns -1, if no depth limit exists) */
extern
int SCIPheurGetMaxdepth(
   HEUR*            heur                /**< primal heuristic */
   );

/** gets the number of times, the heuristic was called and tried to find a solution */
extern
Longint SCIPheurGetNCalls(
   HEUR*            heur                /**< primal heuristic */
   );

/** gets the number of primal feasible solutions found by this heuristic */
extern
Longint SCIPheurGetNSolsFound(
   HEUR*            heur                /**< primal heuristic */
   );

/** is primal heuristic initialized? */
extern
Bool SCIPheurIsInitialized(
   HEUR*            heur                /**< primal heuristic */
   );

/** gets time in seconds used in this heuristic */
extern
Real SCIPheurGetTime(
   HEUR*            heur                /**< primal heuristic */
   );


#endif
