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

/**@file   pub_nlpi.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for NLP solver interfaces
 * @author Thorsten Gellermann
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_NLPI_H__
#define __SCIP_PUB_NLPI_H__

#include "scip/def.h"
#include "scip/type_nlpi.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicNLPIInterfaceMethods
 *
 * @{
 */

/** compares two NLPIs w.r.t. their priority */
SCIP_DECL_SORTPTRCOMP(SCIPnlpiComp);

/** gets data of an NLPI */
SCIP_EXPORT
SCIP_NLPIDATA* SCIPnlpiGetData(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gets NLP solver name */
SCIP_EXPORT
const char* SCIPnlpiGetName(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gets NLP solver descriptions */
SCIP_EXPORT
const char* SCIPnlpiGetDesc(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gets NLP solver priority */
SCIP_EXPORT
int SCIPnlpiGetPriority(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/**@} */ /* PublicNLPIMethods */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_PUB_NLPI_H__ */
