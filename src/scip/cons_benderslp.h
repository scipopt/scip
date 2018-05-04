/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_benderslp.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for benderslp decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_BENDERSLP_H__
#define __SCIP_CONS_BENDERSLP_H__


#include "scip/scip.h"
#include "scip/benders.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for benderslp constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrBenderslp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
