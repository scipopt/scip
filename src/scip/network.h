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
 * matrices can be detected more quickly and are commonly used.
 *
 * Note that all addition algorithms expect that each nonzero is given exactly once and not more often; in particular,
 * it is up to the user to ensure this when using both column and row addition steps.
 *
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

#include "scip.h" //TODO; reduce to minimal set of includes, probably only SCIP and memory functions are needed

#ifdef cplusplus
extern "C" {
#endif

/**
 * This class stores the Network decomposition using an SPQR tree
 */
typedef struct SCIP_NetworkDecomposition SCIP_NETWORKDECOMP;

SCIP_RETCODE SCIPNetworkDecompositionCreate(SCIP * env, SCIP_NETWORKDECOMP **pDecomposition, int numRows, int numColumns);

void SCIPNetworkDecompositionFree(SCIP_NETWORKDECOMP **pDecomposition);

/**
 * Returns if the Network decomposition contains the given row
 */
SCIP_Bool SCIPNetworkDecompositionContainsRow(const SCIP_NETWORKDECOMP * decomposition, int row);

/**
 * Returns if the Network decomposition contains the given column
 */
SCIP_Bool SCIPNetworkDecompositionContainsColumn(const SCIP_NETWORKDECOMP *decomposition, int column);

/**
 * Checks if the Network decomposition of the graph is minimal. This method should only be used in tests.
 */
SCIP_Bool SCIPNetworkDecompositionIsMinimal(const SCIP_NETWORKDECOMP * decomposition);

//TODO: method to convert decomposition into a graph
//TODO: method to remove complete components of the SPQR tree

/**
 * Removes. Removed rows and columns from the SPQR tree.
 * Note that removed rows/columns can not be re-introduced into the SPQR tree; this is possible, but much more expensive
 * to compute and typically not needed.
 */
void SCIPNetworkDecompositionRemoveComponents(SCIP_NETWORKDECOMP *dec, const int * componentRows,
                                              int numRows, const int * componentCols, int numCols);

/**
 * A method to check if the cycle stored in the Decomposition matches the given array. This method should only be used in tests.
 */
SCIP_Bool SCIPNetworkDecompositionVerifyCycle(SCIP * scip, const SCIP_NETWORKDECOMP * dec, int column, int * column_rows,
                                         double * column_values, int num_rows, int * computed_column_storage,
                                         SCIP_Bool * computedSignStorage);

/**
 * This class stores all data for performing sequential column additions to a matrix and checking if it is network or not.
 */
typedef struct SCIP_NetworkColAddition SCIP_NETWORKCOLADDITION;

/**
 * @brief Creates the data structure for managing column-addition for an SPQR decomposition
 */
SCIP_RETCODE SCIPNetworkColAdditionCreate(SCIP* scip, SCIP_NETWORKCOLADDITION** pNewCol );

/**
 * @brief Destroys the data structure for managing column-addition for SPQR decomposition
 */
void SCIPNetworkColAdditionFree(SCIP* scip, SCIP_NETWORKCOLADDITION ** pNewCol);

/**
 * Checks if adding a column of the given matrix creates a network SPQR decomposition.
 * Adding a column which is already in the decomposition is undefined behavior and not checked for.
 * @param dec Current SPQR-decomposition
 * @param newRow Data structure to store information on how to add the new column (if applicable).
 * @param column The index of the column to be added
 * @param rows An array with the row indices of the nonzero entries of the column.
 * @param numRows The number of nonzero entries of the column
 */
SCIP_RETCODE SCIPNetworkColAdditionCheck(SCIP_NETWORKDECOMP * dec, SCIP_NETWORKCOLADDITION * newCol, int column,
                                          const int * nonzeroRows, const double * nonzeroValues, size_t numNonzeros);
/**
 * @brief Adds the most recently checked column from checkNewRow() to the Decomposition.
 * In Debug mode, adding a column for which SPQRNetworkColumnAdditionRemainsNetwork() returns false will exit the program.
 * In Release mode, adding a column for which SPQRNetworkColumnAdditionRemainsNetwork() return false is undefined behavior
 * @param dec Current SPQR-decomposition
 * @param newRow Data structure containing information on how to add the new column.
 */
SCIP_RETCODE SCIPNetworkColAdditionAdd(SCIP_NETWORKDECOMP *dec, SCIP_NETWORKCOLADDITION *newCol);

/**
 * @param newColumn
 * @return True if the most recently checked column is addable to the SPQR decomposition passed to it, i.e. the submatrix
 * given by both remains network.
 */
SCIP_Bool SCIPNetworkColAdditionRemainsNetwork(SCIP_NETWORKCOLADDITION *newCol);

/**
 * This class stores all data for performing sequential row-additions to a matrix and checking if it is network or not.
 */
typedef struct SCIP_NetworkRowAddition SCIP_NETWORKROWADDITION;

/**
 * @brief Creates the data structure for managing row-addition for an SPQR decomposition
 */
SCIP_RETCODE SCIPNetworkRowAdditionCreate(SCIP* scip, SCIP_NETWORKROWADDITION** pNewRow );
/**
 * @brief Destroys the data structure for managing row-addition for SPQR decomposition
 */
void SCIPNetworkRowAdditionFree(SCIP* scip, SCIP_NETWORKROWADDITION ** pNewRow);

/**
 * Checks if adding a row of the given matrix creates a network SPQR decomposition.
 * Adding a row which is already in the decomposition is undefined behavior and not checked for.
 * @param dec Current SPQR-decomposition
 * @param newRow Data structure to store information on how to add the new row (if applicable).
 * @param row The index of the row to be added
 * @param columns An array with the column indices of the nonzero entries of the row.
 * @param numColumns The number of nonzero entries of the row
 */
SCIP_RETCODE SCIPNetworkRowAdditionCheck(SCIP_NETWORKDECOMP * dec, SCIP_NETWORKROWADDITION * newRow, int row,
                                       const int * nonzeroCols, const double * nonzeroValues, size_t numNonzeros);
/**
 * @brief Adds the most recently checked column from checkNewRow() to the Decomposition.
 * In Debug mode, adding a column for which rowAdditionRemainsNetwork() returns false will fail an assertion.
 * In Release mode, adding a column for which rowAdditionRemainsNetwork() return false is undefined behavior
 * @param dec Current SPQR-decomposition
 * @param newRow Data structure containing information on how to add the new row.
 */
SCIP_RETCODE SCIPNetworkRowAdditionAdd(SCIP_NETWORKDECOMP *dec, SCIP_NETWORKROWADDITION *newRow);

/**
 * @param newRow Data structure containing information on how to add the new row
 * @return True if the most recently checked row is addable to the SPQR decomposition passed to it, i.e. the submatrix
 * given by both remains network.
 */
SCIP_Bool SCIPNetworkRowAdditionRemainsNetwork(const SCIP_NETWORKROWADDITION *newRow);

#ifdef cplusplus
}
#endif

#endif //__SCIP_NETWORK_H__
