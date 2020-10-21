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

/**@file   cons_exactsol.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for ensuring that primal solution is exact
 * @author Antonia Chmiela
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/lp.h"
#include "scip/lpexact.h"
#include "scip/pub_var.h"
#include "scip/rational.h"
//#include "scip/struct_rational.h"
#include "scip/scip_exact.h"
#include "scip/scip_lpexact.h"
#include "scip/scip_sol.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/struct_scip.h"
#include "scip/tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for ExactSol constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrExactSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif
