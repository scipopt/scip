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

/**@file   conflict_resolution.h
 * @ingroup OTHER_CFILES
 * @brief  methods and datastructures for resolution-based conflict analysis
 * @author Gioni Mexi
 *
 * @todo
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONFLICT_RESOLUTION_H__
#define __SCIP_CONFLICT_RESOLUTION_H__

#include "scip/def.h"
#include "scip/type_cuts.h"
#include "scip/type_conflict.h"
#include "scip/type_reopt.h"
#include "scip/type_implics.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "lpi/type_lpi.h"
#include "scip/type_branch.h"
#include "scip/type_mem.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_event.h"
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif


#endif
