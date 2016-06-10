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
   SCIP_BASISSTORE**     basisstore           /**< basis storage */
   );

/** frees separation storage */
extern
SCIP_RETCODE SCIPbasisstoreFree(
   SCIP_BASISSTORE**     basisstore           /**< basis storage */
   );

/** add a basis to the storage */
extern
SCIP_RETCODE SCIPbasisstoreAddBasis(
   SCIP_BASISSTORE*      basisstore,          /**< basis storage */
   SCIP_VAR**            vars,                /**< array of problem variables */
   SCIP_CONS**           conss,               /**< array of problem constraints */
   int*                  varstat,             /**< array of variable basis status */
   int*                  consstat,            /**< array of constraint basis status */
   int                   nvars,               /**< number of problem variables */
   int                   nconss               /**< number of problem constraints */
   );

/** return the number of stored basis */
extern
int SCIPbasisstoreGetNBasis(
   SCIP_BASISSTORE*      basisstore           /**< basis storage */
   );

/** returns the stored basis */
extern
SCIP_BASIS* SCIPbasistoreGetBasis(
   SCIP_BASISSTORE*      basistore           /**< basis storage */
   );

/** copy basis storage */
extern
SCIP_RETCODE SCIPbasisstoreCopy(
   SCIP_BASISSTORE*      basisstore,         /**< source basis storage */
   SCIP_BASISSTORE*      targetbasisstore    /**< target basis storage */
   );

/** remove last added basis */
extern
void SCIPbasistoreRemoveBasis(
   SCIP_BASISSTORE*      basisstore         /**< source basis storage */
   );

/**
 * Methods for SCIP_BASIS
 */

/* returns the data of the given basis */
void SCIPbasisGetData(
   SCIP_BASIS*           basis,         /**< source basis storage */
   SCIP_VAR***           vars,          /**< array to store the problem variables */
   SCIP_CONS***          conss,         /**< array to store the problem constraints */
   int**                 varstat,       /**< array to store basis status of variables */
   int**                 consstat,      /**< array to store basis status of constraints */
   int*                  nvars,         /**< pointer to store number of variabls */
   int*                  nconss         /**< pointer to store number of constraints */
   );

#ifdef __cplusplus
}
#endif

#endif
