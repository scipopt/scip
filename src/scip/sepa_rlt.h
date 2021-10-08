/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_rlt.h
 * @ingroup SEPARATORS
 * @brief  reformulation-linearization technique separator
 * @author Fabian Wegscheider
 * @author Ksenia Bestuzheva
 *
 *
 * This separator generates a collection of cuts constructed by the reformulation-linearization technique (RLT).
 * For an LP row L and a variable x in [lb,ub], L is multiplied either with (ub-x) or with (x-lb). All known terms that
 * appear in the product are replaced by their respective auxiliary variable and all unknown terms are replaced by a
 * suitable linear relaxation, e.g., McCormick. In general, the separator computes four different cuts for a row with
 * finite sides and a variable with finite bounds.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_RLT_H__
#define __SCIP_SEPA_RLT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the RLT separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSepaRlt(
   SCIP*                 scip                /**< SCIP data structure */
);

/**@addtogroup SEPARATORS
 *
 * @{
 */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
