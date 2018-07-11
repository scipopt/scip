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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_cyc.h
 * @brief  file reader for cycle clustering instances
 * @author Leon Eifler
 *
 * This file implements the reader for the cycle clustering problem. The data is read from a matrix, entries separated
 * by whitespace. The first line in the file has to be of the form "# p nstates ncluster",
 * where nstates is the size of the matrix and ncluster is the number of clusters that should be used.
 * The file has to have the ending ".cyc" to be recognized by the reader.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_CYC_H__
#define __SCIP_READER_CYC_H__

#include "scip/scip.h"
#include "tclique/tclique.h"
#include "scip/cons_setppc.h"
#include "scip/type_cons.h"
#include "scip/scip.h"
#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the cyc file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderCyc(
   SCIP*             scip                 /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif



#endif
