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

/**@file   reader_lop.h
 * @ingroup FILEREADERS
 * @brief  linear ordering file reader
 * @author Marc Pfetsch
 *
 * This file implements the reader/parser used to read linear ordering problems. For more details see \ref READER. The
 * data should be given in LOLIB format, see <a href="http://www.optsicom.es/lolib/">LOLIB</a>.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_LOP_H__
#define __SCIP_READER_LOP_H__


#include "scip/scip.h"


/** includes the linear ordering file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderLOP(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
