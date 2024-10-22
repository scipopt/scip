/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   iis_deletionfilter.h
 * @brief  deletion filter heuristic to compute IISs
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __IIS_DELETIONFILTER_H__
#define __IIS_DELETIONFILTER_H__

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** addition filter to greedily add constraints to obtain an (I)IS -- detailed function call */
SCIP_EXPORT
SCIP_RETCODE additionFilterBatchCons(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_Bool             silent,             /**< run silently? */
   SCIP_Longint*         nnodes,             /**< pointer to store the total number of nodes needed (or NULL) */
   SCIP_Bool*            success             /**< pointer to store whether we have obtained an (I)IS */
);

/** deletion filter to greedily remove constraints to obtain an (I)IS -- detailed function call */
SCIP_EXPORT
SCIP_RETCODE deletionFilterBatchCons(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_Bool             conservative,       /**< whether we treat a subproblem to be feasible, if it reached its node limt */
   SCIP_Bool             silent,             /**< run silently? */
   SCIP_Longint*         nnodes,             /**< pointer to store the total number of nodes needed (or NULL) */
   SCIP_Bool*            success             /**< pointer to store whether we have obtained an (I)IS */
   );

/** run deletion filter to obtain an (I)IS */
SCIP_EXPORT
SCIP_RETCODE SCIPrunDeletionFilter(
   SCIP*                 scip,               /**< some SCIP instance */
   SCIP_Bool             silent,             /**< run silently? */
   int*                  sizeIS,             /**< pointer to store the size of the (I)IS */
   SCIP_Bool*            isIIS,              /**< pointer to store whether we found an IIS */
   SCIP_Bool*            success             /**< pointer to store whether we have obtained an (I)IS */
   );

#ifdef __cplusplus
}
#endif

#endif
