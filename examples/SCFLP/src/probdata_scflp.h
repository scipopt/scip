/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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

/**@file   probdata_scflp.h
 * @brief  Problem data for Stochastic Capacitated Facility Location problem
 * @author Stephen J. Maher
 *
 * This file handles the main problem data used in that project. For more details see \ref SCFLP_PROBLEMDATA page.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_SCFLP__
#define __SCIP_PROBDATA_SCFLP__

#include "scip/scip.h"

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   SCIP_Real**           costs,              /**< the transportation costs from a facility to a customer */
   SCIP_Real**           demands,            /**< the customer demands */
   SCIP_Real*            capacity,           /**< the capacity of each facility */
   SCIP_Real*            fixedcost,          /**< the fixed cost of opening a facility */
   int                   ncustomers,         /**< the number of customers */
   int                   nfacilities,        /**< the number of facilities */
   int                   nsubproblems,       /**< the number of Benders' decomposition subproblems */
   SCIP_Bool             usebenders,         /**< will Benders' decomposition be used to solve the problem */
   SCIP_Bool             quadcosts           /**< should the problem be formulated with quadratic costs */
   );

/** returns the number of facilities */
int SCIPprobdataGetNFacilities(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns the number of customers  */
int SCIPprobdataGetNCustomers(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns the facility variables */
SCIP_VAR** SCIPprobdataGetFacilityVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

#endif
