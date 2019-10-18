/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   graph.h
 * @brief  Includes dual-ascent for classic Steiner tree and some variants.
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef APPLICATIONS_STP_SRC_DUALASCENT_H_
#define APPLICATIONS_STP_SRC_DUALASCENT_H_

#include "scip/scip.h"
#include "graph.h"

/** dual ascent heuristic */
SCIP_RETCODE SCIPStpDualAscent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real* RESTRICT   redcost,            /**< array to store reduced costs or NULL */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Bool             addcuts,            /**< should dual ascent add Steiner cuts? */
   SCIP_Bool             ascendandprune,     /**< should the ascent-and-prune heuristic be executed? */
   GNODE**               gnodearrterms,      /**< gnode terminals array for internal computations or NULL */
   const int*            result,             /**< solution array (solution needs to be provided) */
   int* RESTRICT         edgearrint,         /**< int edges array for internal computations or NULL */
   int* RESTRICT         nodearrint,         /**< int vertices array for internal computations or NULL */
   int                   root,               /**< the root */
   SCIP_Bool             is_pseudoroot,      /**< is the root a pseudo root? */
   SCIP_Real             damaxdeviation      /**< maximum deviation for dual-ascent ( -1.0 for default) */
   );

/** dual ascent heuristic for the PCSPG and the MWCSP */
SCIP_RETCODE SCIPStpDualAscentPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            redcost,            /**< array to store reduced costs or NULL */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Bool             addcuts,            /**< should dual ascent add Steiner cuts? */
   SCIP_Bool             ascendandprune,     /**< perform ascend-and-prune and add solution? */
   int                   nruns               /**< number of dual ascent runs */
   );


#endif /* APPLICATIONS_STP_SRC_DUALASCENT_H_ */
