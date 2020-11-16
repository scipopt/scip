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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_reader.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for input file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_READER_H__
#define __SCIP_PUB_READER_H__


#include "scip/def.h"
#include "scip/type_reader.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicReaderMethods
 *
 * @{
 */

/** gets user data of reader */
SCIP_EXPORT
SCIP_READERDATA* SCIPreaderGetData(
   SCIP_READER*          reader              /**< reader */
   );

/** sets user data of reader; user has to free old data in advance! */
SCIP_EXPORT
void SCIPreaderSetData(
   SCIP_READER*          reader,             /**< reader */
   SCIP_READERDATA*      readerdata          /**< new reader user data */
   );

/** gets name of reader */
SCIP_EXPORT
const char* SCIPreaderGetName(
   SCIP_READER*          reader              /**< reader */
   );

/** gets description of reader */
SCIP_EXPORT
const char* SCIPreaderGetDesc(
   SCIP_READER*          reader              /**< reader */
   );

/** gets file extension of reader */
SCIP_EXPORT
const char* SCIPreaderGetExtension(
   SCIP_READER*          reader              /**< reader */
   );

/** return whether the reader can read files */
SCIP_EXPORT
SCIP_Bool SCIPreaderCanRead(
   SCIP_READER*          reader              /**< reader */
   );

/** return whether the reader can write files */
SCIP_EXPORT
SCIP_Bool SCIPreaderCanWrite(
   SCIP_READER*          reader              /**< reader */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
