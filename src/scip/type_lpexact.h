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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_lpexact.h
 * @brief  type definitions for exact LP management
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_LPEXACT_H__
#define __SCIP_TYPE_LPEXACT_H__

#include "type_lp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** column of an LP
 *
 *  - \ref PublicColumnMethods "List of all available methods"
 */
typedef struct SCIP_ColExact SCIP_COLEXACT;

/** row of an LP
 *
 *  - \ref PublicRowMethods "List of all available methods"
 */
typedef struct SCIP_RowExact SCIP_ROWEXACT;

/** data used for project and shift */
typedef struct SCIP_ProjShiftData SCIP_PROJSHIFTDATA;

/** LP structure
 *
 *  - \ref PublicLPMethods "List of all available methods"
 */
typedef struct SCIP_LpExact SCIP_LPEXACT;

enum Ps_intpointsel
{
   PS_INTPOINTSEL_ARB         = 0,           /**< arbitrary */
   PS_INTPOINTSEL_OPT         = 1,           /**< optimized */
   PS_INTPOINTSEL_ARBDUAL     = 2,           /**< arbitrary in dual form */
   PS_INTPOINTSEL_TWOSTAGE    = 3            /**< two-stage-optimized interior point */
};
typedef enum Ps_intpointsel PS_INTPOINTSEL;

enum Ps_dualcostsel
{
   PS_DUALCOSTSEL_NO          = 0,           /**< no selection */
   PS_DUALCOSTSEL_ACTIVE_FPLP = 1,           /**< active rows of fp lp */
   PS_DUALCOSTSEL_ACTIVE_EXLP = 2            /**< active rows of exact lp */
};
typedef enum Ps_dualcostsel PS_DUALCOSTSEL;

#ifdef __cplusplus
}
#endif

#endif