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

/**@file   cons_stp.h
 * @brief  Constraint handler for Steiner problems
 * @author Gerald Gamrath
 * @author Daniel Rehfeldt
 * @author Michael Winkler
 *
 * This file checks solutions for feasibility and separates violated model constraints. For more details see \ref CONS page.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_STP_H__
#define __SCIP_CONS_STP_H__


#include "scip/scip.h"
#include "grph.h"

#ifdef __cplusplus
extern "C" {
#endif


#ifndef RESTRICT
#define RESTRICT restrict
#endif

/** creates the handler for element constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrStp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a stp constraint */
extern
SCIP_RETCODE SCIPcreateConsStp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   GRAPH*                graph               /**< graph data structure */
   );

/** sets graph */
void SCIPStpConshdlrSetGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g                   /**< graph data structure */
   );

/** dual ascent heuristic */
extern
SCIP_RETCODE SCIPStpDualAscent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real* RESTRICT   redcost,            /**< array to store reduced costs or NULL */
   SCIP_Real* RESTRICT   nodearrreal,        /**< real vertices array for internal computations or NULL */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Bool             addcuts,            /**< should dual ascent add Steiner cuts? */
   SCIP_Bool             ascendandprune,     /**< should the ascent-and-prune heuristic be executed? */
   GNODE**               gnodearrterms,      /**< gnode terminals array for internal computations or NULL */
   const int*            result,             /**< solution array (solution needs to be provided) */
   int* RESTRICT         edgearrint,         /**< int edges array for internal computations or NULL */
   int* RESTRICT         nodearrint,         /**< int vertices array for internal computations or NULL */
   int                   root,               /**< the root */
   SCIP_Bool             is_pseudoroot,      /**< is the root a pseudo root? */
   SCIP_Real             damaxdeviation,     /**< number of dual ascent runs */
   STP_Bool* RESTRICT    nodearrchar         /**< char vertices array for internal computations or NULL */
   );

/** dual ascent heuristic for the PCSPG and the MWCSP */
extern
SCIP_RETCODE SCIPStpDualAscentPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            redcost,            /**< array to store reduced costs or NULL */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Bool             addcuts,            /**< should dual ascent add Steiner cuts? */
   SCIP_Bool             ascendandprune,     /**< perform ascend-and-prune and add solution? */
   int                   nruns               /**< number of dual ascent runs */
   );


#ifdef __cplusplus
}
#endif

#endif
