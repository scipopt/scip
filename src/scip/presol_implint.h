/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   presol_implint.h
 * @ingroup PRESOLVERS
 * @brief  Presolver that detects implicit integer variables
 * @author Rolf van der Hulst
 *
 * This presolver looks for implicit integer variables, which are variables whose integrality is implied.
 * The linear constraint handler handles the simple (primal) case such as 2x + 2y + z = 3, where z is implied integral by
 * x and y. It also handles a more complicated dual case, where we have 'dual' implied integrality if z occurs only in
 * inequalities of the primal form (where the equality becomes an inequality), and has integral bounds.
 *
 * In this plugin we explicitly look for the following structure in the constraint matrix:
 * \f[
 * \begin{array}{llll}
 * A x & + B y &       & \leq c\\
 * D x &       & + E z & \leq f\\
 * & &               x & \in Z^{p_1} \\
 * & &               y & \in Z^{p_2} \times R^{n_2-p_2}\\
 * & &               z & \in Z^{p_3} \times R^{n_3-p_3}
 * \end{array}
 * \f]
 * where A and c are integral and B is totally unimodular. It is not difficult to see that after fixing the x variables,
 * that the remaining problem on the y variables is an integral polyhedron (and independent of the z variables).
 * Hence, y is implied integral by x.
 *
 * Note that this presolver only treats integral rows, where SCIPisIntegral() is used to check integrality.
 */

#ifndef __SCIP_PRESOL_IMPLINT_H__
#define __SCIP_PRESOL_IMPLINT_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/** creates the implicit integer presolver and includes it in SCIP
 *
 * @ingroup PresolverIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePresolImplint(
   SCIP*                 scip                /**< SCIP data structure */
   );


#ifdef __cplusplus
}
#endif

#endif /* __SCIP_PRESOL_IMPLINT_H__ */
