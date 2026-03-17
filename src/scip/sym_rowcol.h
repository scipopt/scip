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

/**@file   sym_rowcol.h
 * @ingroup SYMHDLRS
 * @brief  symmetry handler for row and column symmetries
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SYM_ROWCOL_H__
#define __SCIP_SYM_ROWCOL_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"
#include "scip/type_sym.h"
#include "symmetry/type_symmetry.h"

#ifdef __cplusplus
extern "C" {
#endif


/** include symmetry handler for row and columns symmetries */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSymhdlrRowCol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** detects whether permutations define single or double lex matrices
 *
 *  A single lex matrix is a matrix whose columns can be partitioned into blocks such that the
 *  columns within each block can be permuted arbitrarily. A double lex matrix is a single lex
 *  matrix such that also blocks of rows have the aforementioned property.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdetectSingleOrDoubleLexMatrices(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_Bool             detectsinglelex,    /**< whether single lex matrices shall be detected */
   int**                 perms,              /**< array of permutations */
   int                   nperms,             /**< number of permutations in perms */
   int                   permlen,            /**< number of variables in a permutation */
   SCIP_Bool*            success,            /**< pointer to store whether structure could be detected */
   SCIP_Bool*            isorbitope,         /**< pointer to store whether detected matrix is orbitopal */
   int***                lexmatrix,          /**< pointer to store single or double lex matrix */
   int*                  nrows,              /**< pointer to store number of rows of lexmatrix */
   int*                  ncols,              /**< pointer to store number of columns of lexmatrix */
   int**                 lexrowsbegin,       /**< pointer to store array indicating begin of new row-lexmatrix */
   int**                 lexcolsbegin,       /**< pointer to store array indicating begin of new col-lexmatrix */
   int*                  nrowmatrices,       /**< pointer to store number of single lex row matrices in rows */
   int*                  ncolmatrices        /**< pointer to store number of single lex column matrices in rows */
   );

#ifdef __cplusplus
}
#endif

#endif
