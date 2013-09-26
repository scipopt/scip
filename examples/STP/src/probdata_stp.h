/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probdata_stp.h
 * @brief  Problem data for stp problem
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
extern
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename           /**< file name */
   );

/** returns the graph */
extern
GRAPH* SCIPprobdataGetGraph(
   SCIP_PROBDATA*        probdata            /**< problem data */
);

/** returns the array with all variables */
extern
SCIP_VAR** SCIPprobdataGetVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array with all edge variables */
extern
SCIP_VAR** SCIPprobdataGetEdgeVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array with all terminals (without the root) */
extern
int* SCIPprobdataGetRTerms(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** returns the number of layers */
extern
int SCIPprobdataGetNLayers(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** returns the number of vars */
extern
int SCIPprobdataGetNVars(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** returns the number of terminals */
extern
int SCIPprobdataGetNTerms(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of terminals without the root node  */
extern
int SCIPprobdataGetRNTerms(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns root */
extern
int SCIPprobdataGetRoot(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of edges */
extern
int SCIPprobdataGetNEdges(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the edge variable for a given index */
extern
SCIP_VAR* SCIPprobdataGetedgeVarByIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   idx
   );

/** returns the LP solution values */
extern
double* SCIPprobdataGetXval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol
   );
/** returns all edge constraints */
extern
SCIP_CONS** SCIPprobdataGetEdgeConstraints(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns all path constraints */
extern
SCIP_CONS** SCIPprobdataGetPathConstraints(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** print (undirected) graph */
extern
SCIP_RETCODE SCIPprobdataPrintGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< Name of the output file */
   SCIP_SOL*             sol,                /**< solution to be printed; or NULL for LP solution */
   SCIP_Bool             printsol            /**< should solution be highlighted? */
   );


/** returns if 'T' model is being used */
extern
SCIP_Bool SCIPprobdataIsBigt(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** add new solution */
extern
SCIP_RETCODE SCIPprobdataAddNewSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            nval,               /**< array [0..nvars], nval[v] = 1 if node v is in the solution, nval[v] = 0 if not */
   SCIP_SOL*             sol,                /**< the new solution */
   SCIP_HEUR*            heur,               /**< heuristic data */
   SCIP_Bool*            success             /**< denotes whether the new solution has been successfully added */
   );
#if 0
/** print (undirected) graph and highlight current solution */
extern
SCIP_RETCODE SCIPprobdataPrintSolGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename           /**< Name of the output file */
   );

#endif
#ifdef __cplusplus
}
#endif

#endif
