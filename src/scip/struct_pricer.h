/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   struct_pricer.h
 * @ingroup INTERNALAPI
 * @brief  data structures for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_PRICER_H__
#define __SCIP_STRUCT_PRICER_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_pricer.h"

#ifdef __cplusplus
extern "C" {
#endif

/** variable pricers data */
struct SCIP_Pricer
{
   char*                 name;               /**< name of variable pricer */
   char*                 desc;               /**< description of variable pricer */
   SCIP_DECL_PRICERCOPY  ((*pricercopy));    /**< copy method of pricer or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PRICERFREE  ((*pricerfree));    /**< destructor of variable pricer */
   SCIP_DECL_PRICERINIT  ((*pricerinit));    /**< initialize variable pricer */
   SCIP_DECL_PRICEREXIT  ((*pricerexit));    /**< deinitialize variable pricer */
   SCIP_DECL_PRICERINITSOL((*pricerinitsol));/**< solving process initialization method of variable pricer */
   SCIP_DECL_PRICEREXITSOL((*pricerexitsol));/**< solving process deinitialization method of variable pricer */
   SCIP_DECL_PRICERREDCOST((*pricerredcost));/**< reduced cost pricing method of variable pricer for feasible LPs */
   SCIP_DECL_PRICERFARKAS((*pricerfarkas));  /**< Farkas pricing method of variable pricer for infeasible LPs */
   SCIP_PRICERDATA*      pricerdata;         /**< variable pricers local data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this pricer for the next stages */
   SCIP_CLOCK*           pricerclock;        /**< pricer execution time */
   int                   priority;           /**< priority of the variable pricer */
   int                   ncalls;             /**< number of times, this pricer was called */
   int                   nvarsfound;         /**< number of variables priced in found so far by this pricer */
   SCIP_Bool             delay;              /**< should the pricer be delayed until no other pricers or already existing
                                              *   problem variables with negative reduced costs are found */
   SCIP_Bool             active;             /**< is variable pricer in use for the current problem? */
   SCIP_Bool             exact;              /**< is variable pricer safe to be used in exact solving mode? */
   SCIP_Bool             initialized;        /**< is variable pricer initialized? */
};

#ifdef __cplusplus
}
#endif

#endif
