/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   relax_stpenum.h
 * @ingroup RELAXATORS
 * @brief  Steiner tree enumeration relaxator
 * @author Daniel Rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RELAX_STPENUM_H__
#define __SCIP_RELAX_STPENUM_H__


#include "scip/scip.h"
#include "graph.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the STP relaxator and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeRelaxStpenum(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** is using the relaxator promising? */
SCIP_EXPORT
SCIP_Bool SCIPStpEnumRelaxIsPromising(
   const GRAPH*          graph               /**< graph */
   );


/** Solve instance by enumeration. Only call when promising. */
SCIP_EXPORT
SCIP_RETCODE SCIPStpEnumRelaxComputeSol(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int* RESTRICT         edges_solstat       /**< solution edges */
   );


#ifdef __cplusplus
}
#endif

#endif
