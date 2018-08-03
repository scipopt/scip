/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probdata_cyc.h
 * @brief  problem data for cycle clustering problem
 * @author Leon Eifler
 *
 * This file implements the problem data for the cycle clustering problem.
 *
 * The problem data contains original transition matrix, the scaling parameter that appears in the objective function,
 * and all variables that appear in the problem.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_CYC__
#define __SCIP_PROBDATA_CYC__

#include "scip/scip.h"
#include "tclique/tclique.h"
#include "scip/cons_setppc.h"
#include "scip/type_cons.h"
#include "scip/def.h"

/** free memory allocated for an nxn matrix */
extern
SCIP_RETCODE freeMatrix(
   SCIP_Real**           matrix,             /**< the matrix to be freed */
   int                   nbins               /**< the size*/
   );

/** gets the minmal non-zero value in a n x n matrix */
extern
SCIP_Real getMinNonZero(
   SCIP*                 scip,               /**< SCIP data structure*/
   SCIP_Real**           matrix,             /**< the matrix*/
   int                   size                /**< the matrix-size*/
   );

/** getter methods for the probdata */
extern
SCIP_Real** SCIPcycGetCmatrix(
   SCIP*                 scip                /**< SCIP data structure*/
   );

/** returns the number of states */
extern
int SCIPcycGetNBins(
   SCIP*                 scip                /**< SCIP data structure*/
   );

/** returns the number of clusters */
extern
int SCIPcycGetNCluster(
   SCIP*                 scip                /**< SCIP data structure*/
   );

/** returns the state-variable-matrix */
extern
SCIP_VAR*** SCIPcycGetBinvars(
   SCIP*                 scip                /**< SCIP data structure*/
   );

/** returns the edge variables */
extern
SCIP_VAR**** SCIPcycGetEdgevars(
   SCIP*                 scip                /**< SCIP data structure*/
   );

/** Return one specific edge variable */
extern
SCIP_VAR* getEdgevar(
   SCIP_VAR****          edgevars,           /**< edgevar data structure*/
   int                   state1,             /**< first state */
   int                   state2,             /**< second state */
   int                   direction           /**< direction, 0 = incluster, 1 = forward */
   );

/** check for an array of states, if all possible edge-combinations exist */
extern
SCIP_Bool edgesExist(
   SCIP_VAR****          edgevars,           /**< edgevar data structure */
   int*                  states,             /**< state array */
   int                   nstates             /**< size of state array */
   );


/** returns the edge-graph */
extern
SCIP_DIGRAPH* SCIPcycGetEdgeGraph(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of scaling parameter */
extern
SCIP_Real SCIPcycGetScale(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** print all the relevant solution data */
extern
SCIP_RETCODE SCIPcycPrintSolutionValues(
   SCIP*               scip,                 /**< SCIP data structure*/
   SCIP_SOL*           sol                   /**< the solution containing the values*/
   );

/** create the probdata for a cycle clustering problem */
extern
SCIP_RETCODE SCIPcreateProbCyc(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */
   int                   nbins,              /**< number of bins */
   int                   ncluster,           /**< number of cluster */
   SCIP_Real**           cmatrix             /**< the transition matrix */
   );

/** function that returns the successive cluster along the cycle */
extern
int phi(
   int                   k,                  /**< the cluster */
   int                   ncluster            /**< the number of clusters*/
   );

/** function that returns the previous cluster along the cycle */
extern
int phiinv(
   int                   k,                  /**< the cluster */
   int                   ncluster            /**< the number of clusters*/
   );

/** assign the variables in scip according to the found clustering. */
extern
SCIP_RETCODE assignVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the SCIP solution */
   SCIP_Real**           clustering,         /**< the matrix with the clusterassignment */
   int                   nbins,              /**< the number of bins */
   int                   ncluster            /**< the number of cluster */
   );

/** check if the clustering has exactly one state in every cluster. */
extern
SCIP_Bool isPartition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           solclustering,      /**< matrix with the clustering */
   int                   nbins,              /**< the number of bins */
   int                   ncluster            /**< the number of clusters */
   );

#endif
