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
SCIP_EXPORT
SCIP_RETCODE SCIPdecompCreate(
   SCIP_DECOMP**         decomp,             /**< pointer to store the decomposition data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   nblocks,            /**< the number of blocks (without the linking block) */
   SCIP_Bool             original,           /**< is this a decomposition in the original (TRUE) or transformed space? */
   SCIP_Bool             benderslabels       /**< should the variables be labeled for the application of Benders' decomposition */
   );

/** free a decomposition */
SCIP_EXPORT
void SCIPdecompFree(
   SCIP_DECOMP**         decomp,             /**< pointer to store the decomposition data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** returns TRUE if decomposition is in the original space */
SCIP_EXPORT
SCIP_Bool SCIPdecompIsOriginal(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** returns TRUE if this constraint contains only linking variables */
SCIP_EXPORT
SCIP_RETCODE SCIPhasConsOnlyLinkVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS*            cons,               /**< the constraint */
   SCIP_Bool*            hasonlylinkvars     /**< will be set to TRUE if this constraint has only linking variables */
   );

/** sets the parameter that indicates whether the variables must be labeled for the application of Benders'
 * decomposition
 */
SCIP_EXPORT
void SCIPdecompSetUseBendersLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Bool             benderslabels       /**< whether Benders' variable labels should be used */
   );

/** returns TRUE if the variables must be labeled for the application of Benders' decomposition */
SCIP_EXPORT
SCIP_Bool SCIPdecompUseBendersLabels(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets number of blocks of this decomposition */
SCIP_EXPORT
int SCIPdecompGetNBlocks(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets number of edges in the block-decomposition graph of this decomposition */
SCIP_EXPORT
int SCIPdecompGetNBlockGraphEdges(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets number of connected components in the block-decomposition graph of this decomposition */
SCIP_EXPORT
int SCIPdecompGetNBlockGraphComponents(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets number of articulation points in the block-decomposition graph of this decomposition */
SCIP_EXPORT
int SCIPdecompGetNBlockGraphArticulations(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets the maximum degree of the block-decomposition graph of this decomposition */
SCIP_EXPORT
int SCIPdecompGetBlockGraphMaxDegree(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets the minimum degree of the block-decomposition graph of this decomposition */
SCIP_EXPORT
int SCIPdecompGetBlockGraphMinDegree(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** set labels for an array of variables */
SCIP_EXPORT
SCIP_RETCODE SCIPdecompSetVarsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int*                  labels,             /**< array of labels, one per variable */
   int                   nvars               /**< length of variables array */
   );

/** query labels for an array of variables */
SCIP_EXPORT
void SCIPdecompGetVarsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int*                  labels,             /**< buffer to store labels, one per variable */
   int                   nvars               /**< length of variables array */
   );

/** set labels for an array of constraints */
SCIP_EXPORT
SCIP_RETCODE SCIPdecompSetConsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int*                  labels,             /**< array of labels, one per constraint */
   int                   nconss              /**< length of constraints array */
   );

/** query labels for an array of constraints */
SCIP_EXPORT
void SCIPdecompGetConsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int*                  labels,             /**< array of labels, one per constraint */
   int                   nconss              /**< length of constraints array */
   );

/** clears the corresponding labeling (constraints, variables, or both) of this decomposition */
SCIP_EXPORT
SCIP_RETCODE SCIPdecompClear(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Bool             clearvarlabels,     /**< should the variable labels be cleared? */
   SCIP_Bool             clearconslabels     /**< should the constraint labels be cleared? */
   );

SCIP_EXPORT
SCIP_RETCODE SCIPdecompComputeConsLabels(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss              /**< number of constraints */
   );

/** create a decomposition of the variables from a labeling of the constraints.
 *
 *  NOTE: by default, the variable labeling is based on a Dantzig-Wolfe decomposition. This means that constraints in named
 *  blocks have have precedence over linking constraints. If a variable exists in constraints from
 *  two or more named blocks, then this variable is marked as a linking variable.
 *  If a variable occurs in exactly one named block i>=0, it is assigned label i.
 *  Variables which are only in linking constraints are unlabeled. However, SCIPdecompGetVarsLabels() will
 *  label them as linking variables.
 *
 *  If the variables should be labeled for the application of Benders' decomposition, the decomposition must be
 *  flagged explicitly via SCIPdecompSetUseBendersLabels().
 *  With this setting, the presence in linking constraints takes precedence over the presence in named blocks.
 *  Now, a variable is considered linking if it is present in at least one linking constraint and an arbitrary
 *  number of constraints from named blocks.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdecompComputeVarsLabels(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss              /**< number of constraints */
   );

/** assign linking constraints to blocks
 *
 * Each linking constraint is assigned to the most frequent block among its variables.
 * Variables of other blocks are relabeled as linking variables.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdecompAssignLinkConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of linking constraints that should be reassigned */
   int                   nconss              /**< number of constraints */
   );

/** compute decomposition statistics and store them in the decomp object */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeDecompStats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** print decomposition statistics into string buffer */
SCIP_EXPORT
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
SCIP_EXPORT
SCIP_RETCODE SCIPdecompstoreCreate(
   SCIP_DECOMPSTORE**    decompstore,        /**< pointer to store decomposition storage */
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   int                   nslots              /**< maximum number of decomposition slots in storage */
   );

/** free a decomposition storage */
SCIP_EXPORT
void SCIPdecompstoreFree(
   SCIP_DECOMPSTORE**    decompstore,        /**< pointer to store decomposition storage */
   BMS_BLKMEM*           blkmem              /**< block memory data structure */
   );

/** add decomposition to storage */
SCIP_EXPORT
SCIP_RETCODE SCIPdecompstoreAdd(
   SCIP_DECOMPSTORE*     decompstore,        /**< decomposition storage */
   SCIP_DECOMP*          decomp              /**< decomposition to add */
   );

/** get decompositions from this storage */
SCIP_EXPORT
SCIP_DECOMP** SCIPdecompstoreGetDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** get number of decompositions in this storage */
SCIP_EXPORT
int SCIPdecompstoreGetNDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** get decompositions in original space from this storage */
SCIP_EXPORT
SCIP_DECOMP** SCIPdecompstoreGetOrigDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** get number of decompositions in original space in this storage */
SCIP_EXPORT
int SCIPdecompstoreGetNOrigDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** get decomposition store from SCIP */
SCIP_EXPORT
SCIP_DECOMPSTORE* SCIPgetDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transform all available original decompositions into transformed space */
SCIP_EXPORT
SCIP_RETCODE SCIPtransformDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** free all decompositions in transformed space */
SCIP_EXPORT
void SCIPexitSolveDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
