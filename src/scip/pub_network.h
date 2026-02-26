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

/**@file    pub_network.h
 * @ingroup PUBLICCOREAPI
 * @brief   Methods for detecting network matrices
 * @author  Rolf van der Hulst
 *
 * This file contains algorithms for incrementally growing (augmenting) network matrices,
 * which are a large subclass of totally unimodular matrices.
 *
 * @addtogroup NetworkMatrix
 *
 * @{
 *
 * A matrix \f$M \in \{-1,0,1\}^{m\times n} \f$ is a network matrix if there exists a directed graph \f$G = (V,A)\f$
 * with \f$|A| = m+n\f$ arcs and a spanning forest \f$T\f$ with \f$|T| = m\f$ such that \f$M\f$'s rows index \f$T\f$ and
 * \f$M\f$'s columns index \f$A\setminus T\f$,
 * and for each arc \f$ a = (u,v) \in A\setminus T\f$ and each arc \f$t\in T\f$
 * \f[
 *   M_{t,a} = \begin{cases}
 *   +1 &  \textrm{if the unique } u-v \textrm{ path in } T \textrm{ passes through } a \textrm{ forwardly}, \\
 *   -1 &  \textrm{if the unique } u-v \textrm{ path in } T \textrm{ passes through } a \textrm{ backwardly}, \\
 *   0 &   \textrm{if the unique } u-v \textrm{ path in } T \textrm{ does not pass through } a
 *   \end{cases}
 * \f]
 * holds.
 *
 * The main difficulty with detecting network matrices is that there may exist many pairs of a graph and a spanning tree
 * that realize a matrix. The algorithms in this file maintain and update an SPQR forest, which is a graph decomposition
 * that represents all graphs that correspond to a given network matrix.
 *
 * Note that all addition algorithms expect that each nonzero is given exactly once and not more often; in particular,
 * it is up to the user to ensure this when interleaving column and row addition steps.
 *
 * More details can be found in:
 * - R.P. van der Hulst and M. Walter "A Row-wise Algorithm for Graph Realization"
 * - R.E. Bixby and D.K. Wagner "An almost linear-time algorithm for graph realization"
 *
 * Note that although these publications contain the methods for undirected graphs (and binary matrices),
 * their ideas are relatively easily extended to directed graphs and ternary matrices.
 * Implementation details are described in further detail in network.c
 */

/* TODO add method that realizes a SCIP digraph from the decomposition */
/* TODO add method that *cleanly* removes complete components of the SPQR tree */
/* TODO add node-arc incidence matrix methods */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NETWORK_H__
#define __SCIP_NETWORK_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_mem.h"
#include "scip/type_misc.h"
#include "blockmemshell/memory.h"

#ifdef cplusplus
extern "C" {
#endif

/** this class stores a decomposition of the network matrix using SPQR trees */
typedef struct SCIP_Netmatdec SCIP_NETMATDEC;

/** create an empty network matrix decomposition that can store a matrix with at most the given dimensions
 *
 * Initially, the network matrix decomposition stores an empty submatrix
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnetmatdecCreate(
   BMS_BLKMEM*           blkmem,             /**< Block memory */
   SCIP_NETMATDEC**      pdec,               /**< buffer to store pointer to created decomposition */
   int                   nrows,              /**< The maximal number of rows that the decomposition can expect */
   int                   ncols               /**< The maximal number of columns that the decomposition can expect */
);

/** frees a network matrix decomposition */
SCIP_EXPORT
void SCIPnetmatdecFree(
   SCIP_NETMATDEC**      pdec                /**< pointer to the network matrix decomposition to free */
);

/** tries to add a given column to the network matrix. Any nonzero row indices in the column that are not yet in the
 *  network matrix decomposition are added, too.
 *
 *  @note Each column and row may be added only once. Trying to add a column that is already in the decomposition will
 *  result in an error.
 *
 *  Note that the first call to this method for a given decomposition may be a bit slower,
 *  due to memory initialization.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnetmatdecTryAddCol(
   SCIP_NETMATDEC*       dec,                /**< Network matrix decomposition */
   int                   column,             /**< The column to add */
   int*                  nonzrows,           /**< The column's nonzero row indices */
   double*               nonzvals,           /**< The column's nonzero entries */
   int                   nnonzs,             /**< The number of nonzeros in the column */
   SCIP_Bool*            success             /**< Buffer to store whether the column was added */
);

/** tries to add a given row to the network matrix. Any nonzero column indices in the row that are not yet in the
 *  network matrix decomposition are added, too.
 *
 *  @note Each column and row may be added only once. Trying to add a row that is already in the decomposition will
 *  result in an error.
 *
 *  Note that the first call to this method for a given decomposition may be a bit slower,
 *  due to memory initialization.
 *
 *  If the user is only interested in determining if a certain (sub)matrix is network or not, using
 *  SCIPnetmatdecTryAddCol() will generally be faster, unless the (sub)matrix has many more columns than rows.
 *
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnetmatdecTryAddRow(
   SCIP_NETMATDEC*       dec,                /**< Network matrix decomposition */
   int                   row,                /**< The row to add */
   int*                  nonzcols,           /**< The row's nonzero column indices */
   double*               nonzvals,           /**< The row's nonzero entries */
   int                   nnonzs,             /**< The number of nonzeros in the row */
   SCIP_Bool*            success             /**< Buffer to store whether the row was added */
);

/** checks if the network matrix decomposition contains the given row */
SCIP_EXPORT
SCIP_Bool SCIPnetmatdecContainsRow(
   SCIP_NETMATDEC*       dec,                /**< The network matrix decomposition */
   int                   row                 /**< The row index that is checked */
);

/** checks if the network matrix decomposition contains the given column */
SCIP_EXPORT
SCIP_Bool SCIPnetmatdecContainsColumn(
   SCIP_NETMATDEC*       dec,                /**< The network matrix decomposition */
   int                   column              /**< The column index that is checked */
);

/** removes a connected component of the matrix from the network decomposition
 *
 * @note This method is 'stupid', and does not delete the associated graph data structure in the decomposition.
 * Moreover, it does not explicitly check if the rows/columns that the user provides are a connected
 * component of the submatrix given by the decomposition. If this is not the case, this will lead to buggy behavior.
 * Use with care!
 */
SCIP_EXPORT
void SCIPnetmatdecRemoveComponent(
   SCIP_NETMATDEC*       dec,                /**< The network matrix decomposition */
   int*                  componentrows,      /**< Array of rows to delete */
   int                   nrows,              /**< The number of rows to delete */
   int*                  componentcols,      /**< Array of columns to delete */
   int                   ncols               /**< The number of columns to delete */
);

/** checks if the network matrix decomposition is minimal.
 *
 * A network matrix decomposition is minimal if it does not contain adjacent parallel or series skeletons.
 * The network matrix decomposition we compute should always be minimal.
 * This method should only be used in tests or asserts.
 */
SCIP_EXPORT
SCIP_Bool SCIPnetmatdecIsMinimal(
   SCIP_NETMATDEC*       dec                 /**< The network matrix decomposition */
);

/** checks if the cycle stored in the Decomposition matches the given array
 *
 * This method should only be used in tests
 */
SCIP_EXPORT
SCIP_Bool SCIPnetmatdecVerifyCycle(
   BMS_BUFMEM *          bufmem,             /**< Buffer memory */
   SCIP_NETMATDEC*       dec,                /**< The network matrix decomposition */
   int                   column,             /**< The column to check */
   int*                  nonzrowidx,         /**< Array with the column's nonzero row indices */
   double*               nonzvals,           /**< Array with the column's nonzero values */
   int                   nnonzs,             /**< Number of nonzeros in the column */
   int*                  pathrowstorage,     /**< A buffer to hold the computed path's rows. Should have size equal or
                                              *   greater than the number of rows in the decomposition. */
   SCIP_Bool*            pathsignstorage     /**< A buffer to store the computed path's row signs. Should have size
                                              * equal or greater than the number of rows in the decomposition. */
);

/** Constructs a realization of an underlying directed graph belonging to the network matrix.
 *
 * The arc data of the digraph contains the row/column indices of the graph. If index < nrows, then the index is the
 * corresponding row. If the index >= nrows, then index-nrows is the column index.
 * Since many different realizations are possible, we use the default orientation of the two-separations to associate
 * pairs of nodes to each other. In particular, we choose to connect the nodes of different 2-connected components
 * in node 0. This way, the rank of the underlying matrix is equal to m+1-c, where c is the number of undirected
 * connected components of the graph where the row/tree edges have been left out.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnetmatdecCreateDiGraph(
   SCIP_NETMATDEC*       dec,                /**< The network matrix decomposition */
   BMS_BLKMEM *          blkmem,             /**< The block memory to use for the created digraph */
   SCIP_DIGRAPH**        pdigraph,           /**< Pointer to the pointer to store the created digraph */
   SCIP_Bool             createrowarcs       /**< Should the row arcs be added to the created digraph? */
);
/**@} */

#ifdef cplusplus
}
#endif

#endif /*__SCIP_NETWORK_H__ */
