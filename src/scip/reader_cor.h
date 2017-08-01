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

/**@file   reader_cor.h
 * @ingroup FILEREADERS
 * @brief  COR file reader (MPS format of the core problem for stochastic programs)
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_COR_H__
#define __SCIP_READER_COR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the cor file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeReaderCor(
   SCIP*                 scip                /**< SCIP data structure */
   );

/*
 * Interface method for the tim and sto readers
 */
EXTERN
SCIP_Bool SCIPcorHasRead(
   SCIP_READER*          reader              /**< the file reader itself */
   );

EXTERN
int SCIPcorGetNVarNames(
   SCIP_READER*          reader              /**< the file reader itself */
   );

EXTERN
int SCIPcorGetNConsNames(
   SCIP_READER*          reader              /**< the file reader itself */
   );

EXTERN
const char* SCIPcorGetVarName(
   SCIP_READER*          reader,             /**< the file reader itself */
   int                   i                   /**< the index of the constraint that is requested */
   );

EXTERN
const char* SCIPcorGetConsName(
   SCIP_READER*          reader,             /**< the file reader itself */
   int                   i                   /**< the index of the constraint that is requested */
   );

#ifdef __cplusplus
}
#endif

#endif
