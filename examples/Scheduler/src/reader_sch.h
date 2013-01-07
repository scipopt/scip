/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_sch.h
 * @brief  scheduling problem file reader for RCPSP/max format
 * @author Stefan Heinz
 *
 * This reader is capabale of parsing resource-constrained project scheduling problem with minimal and maximal time lags
 * (RCPSP/max) instances. The <a http://www.wior.uni-karlsruhe.de/LS_Neumann/Forschung/ProGenMax/rcpspmax.html">PSPlib</a>
 *  provides several instances
 * set.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_SCH_H__
#define __SCIP_READER_SCH_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the sch file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderSch(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
