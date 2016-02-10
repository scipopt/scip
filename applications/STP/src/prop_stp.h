/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_stp.h
 * @brief  propagator for Steiner tree problems, using the LP reduced costs
 * @author Daniel Rehfeldt
 *
 * This propagator makes use of the reduced cost of an optimally solved LP relaxation to propagate the variables
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_STP_H__
#define __SCIP_PROP_STP_H__

#include "scip/scip.h"
#include "grph.h"
#include "probdata_stp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the stp propagator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePropStp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** fix a variable (corresponding to an edge) to zero */
extern
SCIP_RETCODE fixedgevar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             edgevar,            /**< the variable to be fixed */
   int*                  nfixed              /**< counter that is incriminated if variable could be fixed */
   );


#ifdef __cplusplus
}
#endif

#endif
