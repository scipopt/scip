/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_treemodel.h
 * @ingroup INTERNALAPI
 * @brief  Branching rules based on the Single-Variable-Branching (SVB) model
 * @author Daniel Anderson
 *
 * See the following publication for more detail:
 *
 * @par
 * Piere Le Bodic and George Nemhauser@n
 * An abstract model for branching and its application to mixed integer programming@n
 * Mathematical Programming, 2017@n
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_TREEMODEL_H__
#define __SCIP_BRANCH_TREEMODEL_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Initialises the Treemodel parameter data structure */
EXTERN
SCIP_RETCODE SCIPtreemodelInit(
   SCIP*                   scip,        /**< SCIP data structure */
   SCIP_BRANCHTREEMODEL**  treemodel   /**< Treemodel parameter data structure */
   );

/** Frees the Treemodel parameter data structure */
SCIP_RETCODE SCIPtreemodelFree(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_BRANCHTREEMODEL**  treemodel   /**< Treemodel parameter data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
