/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_nl.h
 * @ingroup FILEREADERS
 * @brief  AMPL .nl file reader
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_NL_H__
#define __SCIP_READER_NL_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the .nl file reader into SCIP */
EXTERN
SCIP_RETCODE SCIPincludeReaderNl(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** writes AMPL solution file */
EXTERN
SCIP_RETCODE SCIPwriteAmplSolReaderNl(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          nlreader            /**< AMPL .nl file reader */
   );

#ifdef __cplusplus
}
#endif

#endif
