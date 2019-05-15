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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   decomp.h
 * @ingroup PUBLICAPI
 * @brief  methods for working with decompositions
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef SRC_SCIP_DECOMP_H_
#define SRC_SCIP_DECOMP_H_

#include "scip/type_decomp.h"
#include "scip/type_var.h"
#include "scip/type_cons.h"
#include "blockmemshell/memory.h"


#ifdef __cplusplus
extern "C" {
#endif

/** create a decomposition */
EXTERN
SCIP_RETCODE SCIPdecompCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP**         decomp,             /**< pointer to store the decomposition data structure */
   int                   nblocks,            /**< the number of blocks (without the linking block) */
   SCIP_Bool             original            /**< is this a decomposition in the original (TRUE) or transformed space? */
   );

/** free a decomposition */
EXTERN
void SCIPdecompFree(
   SCIP_DECOMP**         decomp,             /**< pointer to store the decomposition data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** returns TRUE if decomposition is in the original space */
EXTERN
SCIP_Bool SCIPdecompIsOriginal(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** sets the parameter that indicates whether the variables must be labeled for the application of Benders'
 * decomposition
 */
EXTERN
void SCIPdecompSetUseBendersLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Bool             benderslabels       /**< whether Benders' variable labels should be used */
   );

/** returns TRUE if the variables must be labeled for the application of Benders' decomposition */
EXTERN
SCIP_Bool SCIPdecompUseBendersLabels(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets number of blocks of this decomposition */
EXTERN
int SCIPdecompGetNBlocks(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** set labels for an array of variables */
EXTERN
SCIP_RETCODE SCIPdecompSetVarsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int*                  labels,             /**< array of labels, one per variable */
   int                   nvars               /**< length of variables array */
   );

/** query labels for an array of variables */
EXTERN
void SCIPdecompGetVarsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int*                  labels,             /**< buffer to store labels, one per variable */
   int                   nvars               /**< length of variables array */
   );

/** set labels for an array of constraints */
EXTERN
SCIP_RETCODE SCIPdecompSetConsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int*                  labels,             /**< array of labels, one per constraint */
   int                   nconss              /**< length of constraints array */
   );

/** query labels for an array of constraints */
EXTERN
void SCIPdecompGetConsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int*                  labels,             /**< array of labels, one per constraint */
   int                   nconss              /**< length of constraints array */
   );

/** clears the corresponding labeling (constraints, variables, or both) of this decomposition */
EXTERN
SCIP_RETCODE SCIPdecompClear(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Bool             clearvarlabels,     /**< should the variable labels be cleared? */
   SCIP_Bool             clearconslabels     /**< should the constraint labels be cleared? */
   );

EXTERN
SCIP_RETCODE SCIPdecompComputeConsLabels(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss              /**< number of constraints */
   );

/** create a decomposition of the variables from a labeling of the constraints */
EXTERN
SCIP_RETCODE SCIPdecompComputeVarsLabels(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss              /**< number of constraints */
   );

/** compute decomposition statistics and store them in the decomp object */
EXTERN
SCIP_RETCODE SCIPcomputeDecompStats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** print decomposition statistics into string buffer */
EXTERN
char* SCIPdecompPrintStats(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   char*                 strbuf              /**< string buffer storage */
   );

/* author bzfhende
 *
 * TODO query if a variable is a linking variable
 */

/* author bzfhende
 *
 * TODO get a decomposition score, and compute other stuff that may be important
 */

/** create a decomposition storage */
EXTERN
SCIP_RETCODE SCIPdecompstoreCreate(
   SCIP_DECOMPSTORE**    decompstore,        /**< pointer to store decomposition storage */
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   int                   nslots              /**< maximum number of decomposition slots in storage */
   );

/** free a decomposition storage */
EXTERN
void SCIPdecompstoreFree(
   SCIP_DECOMPSTORE**    decompstore,        /**< pointer to store decomposition storage */
   BMS_BLKMEM*           blkmem              /**< block memory data structure */
   );

/** add decomposition to storage */
EXTERN
SCIP_RETCODE SCIPdecompstoreAdd(
   SCIP_DECOMPSTORE*     decompstore,        /**< decomposition storage */
   SCIP_DECOMP*          decomp              /**< decomposition to add */
   );

/** get decompositions from this storage */
EXTERN
SCIP_DECOMP** SCIPdecompstoreGetDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** get number of decompositions in this storage */
EXTERN
int SCIPdecompstoreGetNDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** get decompositions in original space from this storage */
EXTERN
SCIP_DECOMP** SCIPdecompstoreGetOrigDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** get number of decompositions in original space in this storage */
EXTERN
int SCIPdecompstoreGetNOrigDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** get decomposition store from SCIP */
EXTERN
SCIP_DECOMPSTORE* SCIPgetDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transform all available original decompositions into transformed space */
EXTERN
SCIP_RETCODE SCIPtransformDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** free all decompositions in transformed space */
EXTERN
void SCIPexitSolveDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
