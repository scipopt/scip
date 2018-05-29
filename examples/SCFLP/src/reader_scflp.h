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

/**@file   reader_scflp.h
 * @brief  SCFLP problem reader file reader
 * @author Stephen J. Maher
 *
 * This file implements the reader/parser used to read the CAP input data and builds the SCFLP instance. For more
 * details see \ref SCFLP_READER.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_SCFLP_H__
#define __SCIP_READER_SCFLP_H__


#include "scip/scip.h"


/** includes the scflp file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderScflp(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
