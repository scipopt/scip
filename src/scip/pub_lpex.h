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

/**@file   pub_lpex.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for LP management
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_LPEX_H__
#define __SCIP_PUB_LPEX_H__


#include "lpi/type_lpi.h"
#include "lpi/type_lpiex.h"
#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_lpex.h"
#include "scip/type_sepa.h"
#include "scip/type_var.h"
#include "scip/type_misc.h"

#ifdef NDEBUG
#include "scip/struct_lpex.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** comparison method for sorting rows by non-decreasing index */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIProwexComp);

/** returns the left hand side of the row */
SCIP_EXPORT
SCIP_Rational* SCIProwexGetLhs(
   SCIP_ROWEX*           row                 /**< LP row */
   );

/** returns the right hand side of the row */
SCIP_EXPORT
SCIP_Rational* SCIProwexGetRhs(
   SCIP_ROWEX*           row                 /**< LP row */
   );


#ifdef __cplusplus
}
#endif

#endif