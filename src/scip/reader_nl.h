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

/**@file   reader_nl.h
 * @ingroup FILEREADERS
 * @brief  AMPL .nl file reader and writer
 * @author Stefan Vigerske
 *
 * \par Reading Functionality
 *
 * The reader supports linear, nonlinear, and logic constraints, with the following limitations:
 * - For nonlinear expressions, only unary operators minus (negation), abs, pow2, sqrt, log, log10, exp, sin, cos,
 *   binary operators add, sub, mul, div, pow, and n-ary operator sum are supported.
 * - For logical constraints, only operators not, or, and, iff, eq, and ne are supported and all arguments must
 *   be either boolean values or binary variables! The reader currently does not support logical operations that
 *   use algebraic or linear expressions, and therefore not the creation of indicator constraints.
 *
 * In addition, the reader creates special ordered set (SOS) constraints of type 1 and 2 if they were specified via
 * [`sosnr` suffixes](https://discuss.ampl.com/t/how-can-i-use-the-solver-s-special-ordered-sets-feature/45).
 * Values specified via `ref` suffix are passed on as weights to the SOS constraint handlers.
 * For SOS of type 2, the weights determine the order of variables in the SOS.
 *
 * Next to SOS, suffixes can be used to specify flags of variables (see \ref SCIPcreateVar()) and constraints
 * (see \ref SCIPcreateCons()). For variables, supported suffixes are `initial` and `removable`. For constraints,
 * supported suffixes are `initial`, `separate`, `enforce`, `check`, `propagate`, `dynamic`, and `removable`.
 *
 * \par Writing Functionality
 *
 * The writer currently supports the constraint handlers linear, setppc, logicor, knapsack, varbound, and nonlinear only.
 * When writing nonlinear constraints, expression handlers entropy and signpower are currently not supported.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_NL_H__
#define __SCIP_READER_NL_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the .nl file reader into SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeReaderNl(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** writes AMPL solution file
 *
 * problem must have been read with .nl reader
 */
SCIP_EXPORT
SCIP_RETCODE SCIPwriteSolutionNl(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
