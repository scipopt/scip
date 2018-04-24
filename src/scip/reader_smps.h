/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_smps.h
 * @ingroup FILEREADERS
 * @brief  SMPS file reader - smps files list the cor, tim and sto files for a single instance
 * @author Stephen J. Maher
 *
 * An example of an smps file looks as follows.
 *
 * @verbinclude pltexpA2_6.smps
 */

/* author gregor
 *
 * TODO your brief description is wrong. Please expand a bit on the format,
 * specifying the order and input an example file to get a feeling what it should
 * look like.
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_SMPS_H__
#define __SCIP_READER_SMPS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the smps file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeReaderSmps(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
