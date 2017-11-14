/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probdata_snip.h
 * @brief  Problem data for snip problem
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file handles the main problem data used in that project. For more details see \ref SNIP_PROBLEMDATA page.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_SNIP__
#define __SCIP_PROBDATA_SNIP__

#include "scip/scip.h"

/** sets up the problem data */
extern
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   SCIP_Real*            scenariocost,       /**< the costs for the scenarios */
   SCIP_Real*            probwosensor,       /**< the probability of detection without a sensor */
   SCIP_Real*            intdictwosensor,    /**< the probability of detection without a sensor */
   SCIP_Real**           shortestpaths,      /**< the shortest paths for each scenario */
   int*                  scenarioarcids,     /**< the scenario arc ids */
   int*                  arcids,             /**< the arc ids */
   int*                  intdictarcids,      /**< the interdiction arc ids */
   int*                  nodemapping,        /**< mapping from the node ids to the node index */
   SCIP_Real             budget,             /**< the sensor budget */
   SCIP_Real             multiplier,         /**< the probability multiplier */
   int                   narcs,              /**< the number of arcs */
   int                   nnodes,             /**< the number of nodes */
   int                   nsensors,           /**< the number of sensors */
   int                   nscenarios,         /**< the number of scenarios */
   SCIP_Bool             usebenders          /**< will Benders' decomposition be used to solve the problem */
   );

#endif
