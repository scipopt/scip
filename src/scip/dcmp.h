/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dcmp.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for decompositions and the decomposition store
 * @author Gregor Hendel
 *
 * @todo get a decomposition score, and compute other stuff that may be important
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef SRC_SCIP_DECOMP_H_
#define SRC_SCIP_DECOMP_H_

#include "scip/type_dcmp.h"
#include "scip/type_var.h"
#include "scip/type_cons.h"
#include "blockmemshell/memory.h"


#ifdef __cplusplus
extern "C" {
#endif

#define SCIP_DECOMPSTORE_CAPA 10             /**< hardcoded maximum capacity of decomposition store */

/** creates a decomposition storage */
SCIP_RETCODE SCIPdecompstoreCreate(
   SCIP_DECOMPSTORE**    decompstore,        /**< pointer to store decomposition storage */
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   int                   nslots              /**< maximum number of decomposition slots in storage */
   );

/** frees a decomposition storage */
void SCIPdecompstoreFree(
   SCIP_DECOMPSTORE**    decompstore,        /**< pointer to free decomposition storage */
   BMS_BLKMEM*           blkmem              /**< block memory data structure */
   );

/** adds decomposition to storage */
SCIP_RETCODE SCIPdecompstoreAdd(
   SCIP_DECOMPSTORE*     decompstore,        /**< decomposition storage */
   SCIP_DECOMP*          decomp              /**< decomposition to add */
   );

/** transforms all available original decompositions into transformed space */
SCIP_RETCODE SCIPtransformDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** frees all decompositions in transformed space */
void SCIPexitSolveDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets decompositions from storage */
SCIP_DECOMP** SCIPdecompstoreGetDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** gets number of decompositions in storage */
int SCIPdecompstoreGetNDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** gets decompositions in original space from storage */
SCIP_DECOMP** SCIPdecompstoreGetOrigDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

/** gets number of decompositions in original space in storage */
int SCIPdecompstoreGetNOrigDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   );

#ifdef __cplusplus
}
#endif

#endif
