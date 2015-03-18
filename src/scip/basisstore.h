/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   basisstore.h
 * @brief  internal methods for storing starting basis
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BASISSTORE_H__
#define __SCIP_BASISSTORE_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_basisstore.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Methods for SCIP_BASESSTORE
 */

/** creates separation storage */
extern
SCIP_RETCODE SCIPbasisstoreCreate(
   SCIP_BASISSTORE**     basisstore           /**< pointer to store basis storage */
   );

/** frees separation storage */
extern
SCIP_RETCODE SCIPbasisstoreFree(
   SCIP_BASISSTORE**     basisstore,          /**< pointer to store basis storage */
   BMS_BLKMEM*           blkmem
   );

/** add a basis to the storage */
extern
SCIP_RETCODE SCIPbasisstoreAddBasis(
   SCIP_BASISSTORE*      basisstore,
   BMS_BLKMEM*           blkmem,
   SCIP_VAR**            vars,
   SCIP_CONS**           conss,
   int*                  cstat,
   int*                  rstat,
   int                   ncols,
   int                   nrows
   );

/** return the number of stored basis */
extern
int SCIPbasisstoreGetNBasis(
   SCIP_BASISSTORE*      basisstore
   );

/** returns the stored basis */
extern
SCIP_BASIS* SCIPbasestoreGetBasis(
   SCIP_BASISSTORE*      basestore
   );

/** copy the basis */
extern
SCIP_RETCODE SCIPbasisstoreCopy(
   SCIP_BASISSTORE*      basisstore,
   SCIP_BASISSTORE*      cpybasisstore,
   BMS_BLKMEM*           blkmem
   );

extern
void SCIPbasestoreRemoveBasis(
   SCIP_BASISSTORE*      basesstore
   );

/**
 * Methods for SCIP_BASIS
 */

/* returns the data of the given basis */
extern
void SCIPbaseGetData(
   SCIP_BASIS*           basis,
   SCIP_VAR***           vars,
   SCIP_CONS***          conss,
   int**                 vstat,
   int**                 cstat,
   int*                  nvars,
   int*                  nconss
   );

#ifdef __cplusplus
}
#endif

#endif
