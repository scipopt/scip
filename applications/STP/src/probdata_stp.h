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

/**@file   probdata_stp.h
 * @brief  Problem data for stp problem
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Michael Winkler
 * @author Daniel Rehfeldt
 *
 * This file implements the problem data for Steiner problems. For more details see \ref PROBLEMDATA page.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_STP__
#define __SCIP_PROBDATA_STP__

#include "scip/scip.h"
#include "grph.h"

#ifdef __cplusplus
extern "C" {
#endif


/** sets up the problem data */
EXTERN
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   );

/** sets the probdata graph */
EXTERN
void SCIPprobdataSetGraph(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   GRAPH*                graph               /**< graph data structure */
   );

/** returns the graph */
EXTERN
GRAPH* SCIPprobdataGetGraph(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns the graph */
EXTERN
GRAPH* SCIPprobdataGetGraph2(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the offset */
EXTERN
void SCIPprobdataSetOffset(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Real             offset              /**< the offset value */
   );


/** returns the array with all variables */
EXTERN
SCIP_VAR** SCIPprobdataGetVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array with all edge variables */
EXTERN
SCIP_VAR** SCIPprobdataGetEdgeVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array with all terminals (without the root) */
EXTERN
int* SCIPprobdataGetRTerms(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** returns the number of layers */
EXTERN
int SCIPprobdataGetNLayers(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** returns the number of vars */
EXTERN
int SCIPprobdataGetNVars(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** returns the number of terminals */
EXTERN
int SCIPprobdataGetNTerms(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of terminals without the root node  */
EXTERN
int SCIPprobdataGetRNTerms(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns root */
EXTERN
int SCIPprobdataGetRoot(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns numer of original edges */
EXTERN
int SCIPprobdataGetNorgEdges(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of edges */
EXTERN
int SCIPprobdataGetNEdges(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the offset */
EXTERN
SCIP_Real SCIPprobdataGetOffset(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the edge variable for a given index */
EXTERN
SCIP_VAR* SCIPprobdataGetedgeVarByIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   idx
   );

/** returns the LP solution values */
EXTERN
SCIP_Real* SCIPprobdataGetXval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol
   );

/** returns all edge constraints */
EXTERN
SCIP_CONS** SCIPprobdataGetEdgeConstraints(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns all path constraints */
EXTERN
SCIP_CONS** SCIPprobdataGetPathConstraints(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** print (undirected) graph */
EXTERN
SCIP_RETCODE SCIPprobdataPrintGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< Name of the output file */
   SCIP_SOL*             sol,                /**< solution to be printed; or NULL for LP solution */
   SCIP_Bool             printsol            /**< should solution be highlighted? */
   );

/** print graph (in undirected form) in GML format with given edges highlighted */
EXTERN
SCIP_RETCODE SCIPprobdataPrintGraph2(
   const GRAPH*          graph,              /**< Graph to be printed */
   const char*           filename,           /**< Name of the output file */
   SCIP_Bool*            edgemark            /**< Array of (undirected) edges to highlight */
   );

/** returns if 'T' model is being used */
EXTERN
SCIP_Bool SCIPprobdataIsBigt(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** writes the best solution to the intermediate solution file */
EXTERN
SCIP_RETCODE SCIPprobdataWriteIntermediateSolution(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** writes SPG (no variant!) to a file */
void SCIPprobdataWriteStp(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const char*           filename            /**< file name */
   );

/** writes the best solution to a file */
EXTERN
SCIP_RETCODE SCIPprobdataWriteSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< file to write best solution to; or NULL, to write to stdout */
   );

/** writes a line to the log file */
EXTERN
void SCIPprobdataWriteLogLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** add new solution */
EXTERN
SCIP_RETCODE SCIPprobdataAddNewSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            nval,               /**< array [0..nvars], nval[v] = 1 if node v is in the solution, nval[v] = 0 if not */
   SCIP_SOL*             sol,                /**< the new solution */
   SCIP_HEUR*            heur,               /**< heuristic data */
   SCIP_Bool*            success             /**< denotes whether the new solution has been successfully added */
   );

/** set dual bound by ug */
EXTERN
void SCIPprobdataSetDualBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             dual
   );

/** set the number of solvers */
EXTERN
void SCIPprobdataSetNSolvers(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nSolvers            /**< the number of solvers */
   );

/** returns problem type */
EXTERN
int SCIPprobdataGetType(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** writes end of log file */
EXTERN
SCIP_RETCODE SCIPprobdataWriteLogfileEnd(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** branching information from UG */
EXTERN
void initReceivedSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const int             lLinearConsNames,   /**< number of linear constraints */
   const char*           linearConsNames,    /**< linear constraints string */
   const int             lSetppcConsNames,   /**< number of setppc constraints */
   const char*           setppcConsNames     /**< number of setppc constraints */
   );

#ifdef __cplusplus
}
#endif

#endif
