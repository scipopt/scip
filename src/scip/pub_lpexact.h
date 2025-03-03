/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_lpexact.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for LP management
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_LPEXACT_H__
#define __SCIP_PUB_LPEXACT_H__


#include "lpi/type_lpi.h"
#include "lpiexact/type_lpiexact.h"
#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_lpexact.h"
#include "scip/type_sepa.h"
#include "scip/type_var.h"
#include "scip/type_misc.h"

#ifdef NDEBUG
#include "scip/struct_lpexact.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** comparison method for sorting rows by non-decreasing index */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIProwExactComp);

/** gets variable this column represents */
SCIP_EXPORT
SCIP_VAR* SCIPcolExactGetVar(
   SCIP_COLEXACT*        col                 /**< LP column */
   );


/** returns the left hand side of the row */
SCIP_EXPORT
SCIP_Rational* SCIProwExactGetLhs(
   SCIP_ROWEXACT*        row                 /**< LP row */
   );

/** returns the right hand side of the row */
SCIP_EXPORT
SCIP_Rational* SCIProwExactGetRhs(
   SCIP_ROWEXACT*        row                 /**< LP row */
   );

/** returns the constant of the row */
SCIP_EXPORT
SCIP_Rational* SCIProwExactGetConstant(
   SCIP_ROWEXACT*        row                 /**< LP row */
   );

/** gets array with coefficients of nonzero entries */
SCIP_EXPORT
SCIP_Rational** SCIProwExactGetVals(
   SCIP_ROWEXACT*        row                 /**< LP row */
   );

/** sorts row entries such that LP columns precede non-LP columns and inside both parts lower column indices precede
 *  higher ones
 */
SCIP_EXPORT
void SCIProwExactSort(
   SCIP_ROWEXACT*        row                 /**< row to be sorted */
   );

/** gets array of exact columns */
SCIP_EXPORT
SCIP_COLEXACT** SCIProwExactGetCols(
   SCIP_ROWEXACT*        row                 /**< LP row */
   );

/** locks an unmodifiable row, which forbids further changes; has no effect on modifiable rows */
SCIP_EXPORT
void SCIProwExactLock(
   SCIP_ROWEXACT*        row                 /**< LP row */
   );

/** unlocks a lock of an unmodifiable row; a row with no sealed lock may be modified; has no effect on modifiable rows */
SCIP_EXPORT
void SCIProwExactUnlock(
   SCIP_ROWEXACT*        row                 /**< LP row */
   );

/** returns exact row corresponding to fprow, if it exists. Otherwise returns NULL */
SCIP_EXPORT
SCIP_ROWEXACT* SCIProwGetRowExact(
   SCIP_ROW*             row                 /**< SCIP row */
   );

/** returns fp row corresponding to exact row, if it exists. Otherwise returns NULL */
SCIP_EXPORT
SCIP_ROW* SCIProwExactGetRow(
   SCIP_ROWEXACT*        row                 /**< SCIP row */
   );

/** returns rhs-relaxation part of exact row, if it exists. Otherwise returns NULL */
SCIP_EXPORT
SCIP_ROW* SCIProwExactGetRowRhs(
   SCIP_ROWEXACT*        row                 /**< SCIP row */
   );

/** true if row can be relaxed (possibly as two fp rows) */
SCIP_EXPORT
SCIP_Bool SCIProwExactHasFpRelax(
   SCIP_ROWEXACT*        row                 /**< SCIP row */
   );

/** returns whether the exact LP is in diving mode */
SCIP_EXPORT
SCIP_Bool SCIPlpExactDiving(
   SCIP_LPEXACT*         lpexact             /**< current exact LP data */
   );

#ifdef __cplusplus
}
#endif

#endif
