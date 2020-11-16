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

/**@file   reader_smps.h
 * @ingroup FILEREADERS
 * @brief  SMPS file reader - SMPS files lists the cor, tim and sto files for a single instance to be read
 * @author Stephen J. Maher
 *
 * This is a reader for the reading the core, time and stochastic files  of a stochastic programming instance
 * in SMPS format.
 * The three files that must be listed in the smps file are:
 * - .cor
 * - .tim
 * - .sto
 * The order of the files is not important. It is just necessary that all three are provided.
 *
 * Details regarding the SMPS file format can be found at:
 * Birge, J. R.; Dempster, M. A.; Gassmann, H. I.; Gunn, E.; King, A. J. & Wallace, S. W.
 * A standard input format for multiperiod stochastic linear programs
 * IIASA, Laxenburg, Austria, WP-87-118, 1987
 *
 * An example of an smps file looks as follows.
 *
 * @verbinclude pltexpA2_6.smps
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_SMPS_H__
#define __SCIP_READER_SMPS_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the smps file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeReaderSmps(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
