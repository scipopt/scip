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

/**@file   pub_cutsel.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for cut selectors
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_CUTSEL_H__
#define __SCIP_PUB_CUTSEL_H__


#include "scip/def.h"
#include "scip/type_cutsel.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicCutSelectorMethods
 *
 * @{
 */

/** gets name of cut selector */
SCIP_EXPORT
const char* SCIPcutselGetName(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
);

/** gets user data of cut selector */
SCIP_EXPORT
SCIP_CUTSELDATA* SCIPcutselGetData(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );


/** @} */

#ifdef __cplusplus
}
#endif

#endif
