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

/**@file   reader_pbm.h
 * @ingroup FILEREADERS
 * @brief  file writer for portable bitmap file format (PBM), open with common graphic viewer programs (e.g. xview)
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_PBM_H__
#define __SCIP_READER_PBM_H__

#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_reader.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the pbm file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeReaderPbm(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup FILEREADERS
 *
 * @{
 */

/* writes picture of matrix structure to file */
SCIP_EXPORT
SCIP_RETCODE SCIPwritePbm(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   const char*           name,               /**< problem name */
   SCIP_READERDATA*      readerdata,         /**< information for reader */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   int                   nvars,              /**< number of active variables in the problem */
   SCIP_CONS**           conss,              /**< array with constraints of the problem */
   int                   nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
