/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reopt.h
 * @brief  internal methods for collecting primal CIP solutions and primal informations
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_REOPT_H__
#define __SCIP_REOPT_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_reopt.h"

#include "scip/struct_reopt.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates reopt data */
extern
SCIP_RETCODE SCIPreoptCreate(
   SCIP_REOPT**          reopt              /**< pointer to primal data */
   );

/** frees reopt data */
extern
SCIP_RETCODE SCIPreoptFree(
   SCIP*                 scip,
   SCIP_REOPT**          reopt,             /**< pointer to primal data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** add a solution to sols */
extern
SCIP_RETCODE SCIPreoptAddSol(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,
   SCIP_SET*             set,
   SCIP_STAT*            stat,
   SCIP_SOL*             sol,
   SCIP_Bool             bestsol,
   SCIP_Bool*            added,
   int                   run
   );

/* add a run */
extern
SCIP_RETCODE SCIPreoptAddRun(
   SCIP_SET*             set,
   SCIP_REOPT*           reopt,
   int                   run,
   int                   size
   );

/* returns number of solution */
extern
int SCIPreoptGetNSols(
   SCIP_REOPT*           reopt,
   int                   run
   );

/* add solutions to origprimal space */
extern
SCIP_RETCODE SCIPreoptUpdateSols(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,
   int                   savesols,
   SCIP_Real             objsim,
   SCIP_Real             sepabestsol
   );

/* returns the number of saved solutions overall runs */
extern
int SCIPreoptNSavedSols(
   SCIP_REOPT*           reopt
   );

/* returns the number of reused sols over all runs */
extern
int SCIPreoptNUsedSols(
   SCIP_REOPT*           reopt
   );

/* save objective function */
extern
SCIP_RETCODE SCIPreoptSaveObj(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,
   SCIP_SET*             set,
   int                   run
   );

/* check if the current and the previous objective are similar enough
 * returns TRUE if we want to restart, otherwise FALSE */
extern
SCIP_Bool SCIPreoptCheckRestart(
   SCIP_REOPT*           reopt,
   SCIP_SET*             set,
   int                   nvars,
   SCIP_Real*            sim
   );

/* returns an array of indices with similar objective functions
 * to obj_idx */
extern
int* SCIPreoptGetSimilarityIdx(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,
   int                   obj_id,
   int*                  sim_ids,
   int*                  nids
   );

/*
 * returns the similarity to the previous objective function
 */
extern
SCIP_Real SCIPreoptGetSimToPrevious(
      SCIP_REOPT*        reopt
   );

/*
 * returns the similarity to the first objective function
 */
extern
SCIP_Real SCIPreoptGetSimToFirst(
      SCIP_REOPT*        reopt
   );

/*
 * returns the best solution of the last run
 */
extern
SCIP_SOL* SCIPreoptGetLastBestSol(
   SCIP_REOPT*           reopt
   );

/*
 * returns true if the sign of the objective coefficient has changed
 */
extern
SCIP_Bool SCIPreoptIsObjCoefNegated(
   SCIP_REOPT*           reopt,
   int                   varnr
   );

/*
 * returns the coefficent of variable with index @param idx in run @param run
 */
SCIP_Real SCIPreoptGetObjCoef(
   SCIP_REOPT*           reopt,
   int                   run,
   int                   idx
   );

/*
 * returns if a variable enters or leaves the objective function.
 */
extern
void SCIPreoptEnterOrLeaveObj(
   SCIP_REOPT*           reopt,
   int                   varnr,
   SCIP_Bool*            entering,
   SCIP_Bool*            leaving
   );

#ifdef __cplusplus
}
#endif

#endif
