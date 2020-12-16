/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_stpcomponents.h
 * @brief  Components constraint handler for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file checks solutions for feasibility and separates violated model constraints. For more details see \ref STP_CONS page.
 */

#ifndef APPLICATIONS_STP_SRC_CONS_STPCOMPONENTS_H_
#define APPLICATIONS_STP_SRC_CONS_STPCOMPONENTS_H_


#include "scip/scip.h"
#include "graph.h"

#ifdef __cplusplus
extern "C" {
#endif


/** creates the handler for element constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrStpcomponents(
   SCIP*                 scip                /**< SCIP data structure */
   );


#ifdef __cplusplus
}
#endif


#endif /* APPLICATIONS_STP_SRC_CONS_STPCOMPONENTS_H_ */
