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

/**@file    network.h
 * @ingroup OTHER_CFILES
 * @brief   Methods for detecting network matrices
 * @author  Rolf van der Hulst
 *
 * This file contains algorithms for incrementally growing (augmenting) network matrices,
 * which are a large subclass of totally unimodular matrices.
 *
 * A $\pm 1$ matrix is a network matrix if there exists a directed graph G with a (not necessarily rooted) tree T such that
 * for each arc \f$ a \in A\setminus T\f$ the column \f$ A_{:,a} \f$ corresponds to the oriented path between the tail
 * and the head of \f$ a \f$. The main difficulty with detecting network matrices is that there may exist many graphs
 * that realize a certain matrix. The algorithms in this file maintain and update an SPQR tree,
 * which is a graph decomposition that represents all graphs that correspond to the current network matrix.
 *
 * TODO: add incidence matrix methods
 * A large subclass of network matrices are node-arc incidence matrices, which are matrices that have at most
 * one +1 entry and one -1 entry in every column, which correspond to the in and out-node of an arc.
 * This file contains algorithms to detect reflected node-arc incidence matrices, where rows can be optionally negated.
 * Although network matrices are a strictly larger class of totally unimodular matrices, reflected node-arc incidence
 * matrices can be detected more quickly and are commonly used within MIP formulations
 *
 * Note that all addition algorithms expect that each nonzero is given exactly once and not more often; in particular,
 * it is up to the user to ensure this when using both column and row addition steps.
 *
 * TODO
 * The column addition for network matrices is based on:
 *
 * The row addition for network matrices is based on:
 *
 * The row addition for incidence matrices is based on an adaptation from:
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NETWORK_H__
#define __SCIP_NETWORK_H__

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#ifdef cplusplus
extern "C" {
#endif

/** this class stores a decomposition of the network matrix using SPQR trees */
typedef struct SCIP_Netmatdec SCIP_NETMATDEC;

/** create an empty network matrix decomposition that can store a matrix with at most the given dimensions */
SCIP_EXPORT
SCIP_RETCODE SCIPnetmatdecCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NETMATDEC**      pdec,               /**< buffer to store pointer to created decomposition */
   int                   nrows,              /**< The maximal number of rows that the decomposition can expect */
   int                   ncols               /**< The maximal number of columns that the decomposition can expect */
);

/** frees a network matrix decomposition */
SCIP_EXPORT
void SCIPnetmatdecFree(
   SCIP_NETMATDEC**      pdec                /**< pointer to the network matrix decomposition to freed */
);

/** checks if the network matrix decomposition contains the given row */
SCIP_EXPORT
SCIP_Bool SCIPnetmatdecContainsRow(
   const SCIP_NETMATDEC* dec,                /**< The network matrix decomposition */
   int                   row                 /**< The row index that is checked */
);

/** checks if the network matrix decomposition contains the given column */
SCIP_EXPORT
SCIP_Bool SCIPnetmatdecContainsColumn(
   const SCIP_NETMATDEC* dec,                /**< The network matrix decomposition */
   int                   column              /**< The column index that is checked */
);

/** checks if the network matrix decomposition is minimal; it is minimal if it does not contain adjacent parallel or series skeletons; should only be used in tests or asserts*/
SCIP_EXPORT
SCIP_Bool SCIPnetmatdecIsMinimal(
   const SCIP_NETMATDEC* dec                 /**< The network matrix decomposition */
);

//TODO: add method that realizes a SCIP digraph from the decomposition
//TODO: add method that *cleanly* removes complete components of the SPQR tree

/** removes a connected component of the matrix from the network decomposition; note that this method is 'stupid',
 * and does not delete the associated graph data structure; moreover, it does not explicitly check if the rows/columns
 * that the user provides are a connected component of the matrix given by the decomposition. If this is not the case,
 * then calling this function is considered a bug. */
SCIP_EXPORT
void SCIPnetmatdecRemoveComponent(
   SCIP_NETMATDEC*       dec,                /**< The network matrix decomposition */
   const int*            componentrows,      /**< Pointer to the array of rows to delete */
   int                   nrows,              /**< The number of rows to delete */
   const int*            componentcols,      /**< Pointer to the array of columns to delete */
   int                   ncols               /**< The number of columns to delete */
);

/** checks if the cycle stored in the Decomposition matches the given array; should only be used in tests */
SCIP_EXPORT
SCIP_Bool SCIPnetmatdecVerifyCycle(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_NETMATDEC* dec,                /**< The network matrix decomposition */
   int                   column,             /**< The column to check */
   const int*            nonzrowidx,         /**< Pointer to the array with the column's nonzero row indices */
   const double*         nonzvals,           /**< Pointer to the array with the column's nonzero values */
   int                   nnonzs,             /**< Number of nonzeros in the column */
   int*                  pathrowstorage,     /**< Pointer to a buffer to hold the computed path's rows */
   SCIP_Bool*            pathsignstorage     /**< Pointer to a buffer to store the computed path's row signs */
);

/** this class stores all data for performing a column addition to the network matrix decomposition */
typedef struct SCIP_NetColAdd SCIP_NETCOLADD;

/** creates the data structure for managing column addition of a network matrix decomposition */
SCIP_EXPORT
SCIP_RETCODE SCIPnetcoladdCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NETCOLADD**      pcoladd             /**< buffer to store pointer to column addition data structure */
);

/** frees the data structure for managing column addition of a network matrix decomposition */
SCIP_EXPORT
void SCIPnetcoladdFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NETCOLADD**      pcoladd             /**< pointer to the column addition data structure to be freed */
);

/** Checks if we can add the given column to the network matrix decomposition.
 * The result of the latest query can be checked through SCIPnetcoladdRemainsNetwork().
 * Users can add the latest checked column through SCIPnetcoladdAdd()*/
SCIP_EXPORT
SCIP_RETCODE SCIPnetcoladdCheck(
   SCIP_NETMATDEC*       dec,                /**< Network matrix decomposition */
   SCIP_NETCOLADD*       coladd,             /**< Network matrix column addition data structure */
   int                   column,             /**< The column to check */
   const int*            nonzrows,           /**< The column's nonzero row indices */
   const double*         nonzvals,           /**< The column's nonzero entries */
   size_t                nnonzs              /**< The number of nonzeros in the column */
);
/** Adds the most recently checked column to the network matrix decomposition. */
SCIP_EXPORT
SCIP_RETCODE SCIPnetcoladdAdd(
   SCIP_NETMATDEC*       dec,                /**< Network matrix decomposition */
   SCIP_NETCOLADD*       coladd              /**< Network matrix column addition data structure */
);

/** Returns whether the most recently checked column can be added to the network */
SCIP_EXPORT
SCIP_Bool SCIPnetcoladdRemainsNetwork(
   const SCIP_NETCOLADD* coladd              /**< Network matrix column addition data structure */
);

/** this class stores all data for performing a row addition to the network matrix decomposition */
typedef struct SCIP_NetRowAdd SCIP_NETROWADD;

/** creates the data structure for managing row addition of a network matrix decomposition */
SCIP_EXPORT
SCIP_RETCODE SCIPnetrowaddCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NETROWADD**      prowadd             /**< buffer to store pointer to row addition data structure */
);

/** frees the data structure for managing row addition of a network matrix decomposition */
SCIP_EXPORT
void SCIPnetrowaddFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NETROWADD**      prowadd             /**< pointer to the row addition data structure to be freed */
);

/** Checks if we can add the given row to the network matrix decomposition.
 * The result of the latest query can be checked through SCIPnetrowaddRemainsNetwork().
 * Users can add the latest checked column through SCIPnetrowaddAdd()*/
SCIP_EXPORT
SCIP_RETCODE SCIPnetrowaddCheck(
   SCIP_NETMATDEC*       dec,                /**< Network matrix decomposition */
   SCIP_NETROWADD*       rowadd,             /**< Network matrix row addition data structure */
   int                   row,                /**< The row to check */
   const int*            nonzcols,           /**< The row's nonzero row indices */
   const double*         nonzvals,           /**< The row's nonzero entries */
   size_t                nnonzs              /**< The number of nonzeros in the row */
);
/** Adds the most recently checked row to the network matrix decomposition. */
SCIP_EXPORT
SCIP_RETCODE SCIPnetrowaddAdd(
   SCIP_NETMATDEC*       dec,                /**< Network matrix decomposition */
   SCIP_NETROWADD*       rowadd              /**< Network matrix row addition data structure */
);

/** Returns whether the most recently checked row can be added to the network */
SCIP_EXPORT
SCIP_Bool SCIPnetrowaddRemainsNetwork(
   const SCIP_NETROWADD* rowadd              /**< Network matrix row addition data structure */
   );

#ifdef cplusplus
}
#endif

#endif //__SCIP_NETWORK_H__
