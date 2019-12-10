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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   treemodel.h
 * @ingroup PUBLICCOREAPI
 * @brief  Branching rules based on the Single-Variable-Branching (SVB) model
 * @author Daniel Anderson
 *
 * The Single-Variable-Branching (SVB) model is a simplified model of
 * Branch & Bound trees, from which several nontrivial variable selection
 * rules arise. The Treemodel branching rule complements SCIP's hybrid
 * branching by suggesting improved branching variables given the current
 * pseudocosts and the current dual gap.
 *
 * See the following publication for more detail:
 *
 * @par
 * Pierre Le Bodic and George Nemhauser@n
 * An abstract model for branching and its application to mixed integer programming@n
 * Mathematical Programming, 2017@n
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TREEMODEL_H__
#define __SCIP_TREEMODEL_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** initialises the Treemodel parameter data structure */
SCIP_EXPORT
SCIP_RETCODE SCIPtreemodelInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL**      treemodel           /**< Treemodel parameter data structure */
);

/** frees the Treemodel parameter data structure */
SCIP_EXPORT
SCIP_RETCODE SCIPtreemodelFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL**      treemodel           /**< Treemodel parameter data structure */
);

/** returns TRUE if the Treemodel branching rules are enabled */
SCIP_EXPORT
SCIP_Bool SCIPtreemodelIsEnabled(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL*       treemodel           /**< Treemodel parameter data structure */
);

/** apply the Treemodel branching rules to attempt to select a better
 *  branching candidate than the one selected by pseudocost branching
 */
SCIP_EXPORT
SCIP_RETCODE SCIPtreemodelSelectCandidate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL*       treemodel,          /**< Treemodel parameter data structure */
   SCIP_VAR**            branchcands,        /**< branching candidate storage */
   SCIP_Real*            mingains,           /**< minimum gain of rounding downwards or upwards */
   SCIP_Real*            maxgains,           /**< maximum gain of rounding downwards or upwards */
   SCIP_Real*            scoresfrompc,       /**< pseudocost scores of branching candidates */
   SCIP_Real*            scoresfromothers,   /**< scores from other branching rules */
   SCIP_Real             avgpscostscore,     /**< average pseudocost score of branching candidates */
   int                   nbranchcands,       /**< the number of branching candidates */
   int*                  bestcand            /**< the best branching candidate found before the call, 
					          and the best candidate after the call (possibly the same) */
);

#ifdef __cplusplus
}
#endif

#endif
