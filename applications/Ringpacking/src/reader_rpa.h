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

/**@file   reader_rpa.c
 * @brief  Ringpacking problem reader
 * @author Benjamin Mueller
 *
 * This file implements the reader/parser used to read the ringpacking input data. For more details see \ref READER.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_RPA_H__
#define __SCIP_READER_RPA_H__


#include "scip/scip.h"


/** includes the rpa file reader in SCIP */
extern
SCIP_RETCODE SCIPincludeReaderRpa(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
