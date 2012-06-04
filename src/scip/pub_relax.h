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

/**@file   pub_relax.h
 * @ingroup PUBLICMETHODS
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

/** compares two relaxation handlers w. r. to their priority */
extern
SCIP_DECL_SORTPTRCOMP(SCIPrelaxComp);

/** comparison method for sorting relaxators w.r.t. to their name */
extern
SCIP_DECL_SORTPTRCOMP(SCIPrelaxCompName);

/** gets user data of relaxation handler */
extern
SCIP_RELAXDATA* SCIPrelaxGetData(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** sets user data of relaxation handler; user has to free old data in advance! */
extern
void SCIPrelaxSetData(
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_RELAXDATA*       relaxdata           /**< new relaxation handler user data */
   );

/** gets name of relaxation handler */
extern
const char* SCIPrelaxGetName(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets description of relaxation handler */
extern
const char* SCIPrelaxGetDesc(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets priority of relaxation handler */
extern
int SCIPrelaxGetPriority(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets frequency of relaxation handler */
extern
int SCIPrelaxGetFreq(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets time in seconds used in this relaxator for setting up for next stages */
extern
SCIP_Real SCIPrelaxGetSetupTime(
   SCIP_RELAX*           relax               /**< relaxator */
   );

/** gets time in seconds used in this relaxation handler */
extern
SCIP_Real SCIPrelaxGetTime(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** gets the total number of times, the relaxation handler was called */
extern
SCIP_Longint SCIPrelaxGetNCalls(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** is relaxation handler initialized? */
extern
SCIP_Bool SCIPrelaxIsInitialized(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

/** marks the current relaxation unsolved, s.t. the relaxation handler is called again in the next solving round */
extern
void SCIPrelaxMarkUnsolved(
   SCIP_RELAX*           relax               /**< relaxation handler */
   );

#ifdef __cplusplus
}
#endif

#endif
