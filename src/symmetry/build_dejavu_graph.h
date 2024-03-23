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

/**@file   build_dejavu_graph.h
 * @brief  methods to build dejavu graph for symmetry detection
 * @author Christopher Hojny
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BUILD_DEJAVU_GRAPH_H_
#define __SCIP_BUILD_DEJAVU_GRAPH_H_

#include "scip/scip.h"

/* include dejavu */
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#pragma GCC diagnostic ignored "-Wunused-private-field"
#endif

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4189)  // local variable is initialized but not referenced
# pragma warning(disable: 4388)  // compare signed and unsigned expression
# pragma warning(disable: 4456)  // shadowed variable
# pragma warning(disable: 4430)  // missing type specifier
#endif

/* the actual include */
#include <dejavu/dejavu.h>

#ifdef __GNUC__
#pragma GCC diagnostic warning "-Wunused-but-set-variable"
#pragma GCC diagnostic warning "-Wsign-compare"
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wshadow"
#pragma GCC diagnostic warning "-Wpedantic"
#pragma GCC diagnostic warning "-Wdeprecated-copy"
#pragma GCC diagnostic warning "-Wnon-virtual-dtor"
#pragma GCC diagnostic warning "-Wunused-private-field"
#endif

#ifdef _MSC_VER
# pragma warning(pop)
#endif


#ifdef __cplusplus
extern "C" {
#endif

#include "symmetry/struct_symmetry.h"
#include "symmetry/type_symmetry.h"

/** compute generators of symmetry group */
SCIP_EXPORT
SCIP_RETCODE SYMbuildDejavuGraph(
   SCIP*                 scip,               /**< SCIP pointer */
   dejavu::static_graph* dejavugraph,        /**< pointer to hold dejavu graph being created */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_Bool*            success             /**< pointer to store whether dejavugraph could be built */
   );


/** returns whether two given graphs are identical */
SCIP_EXPORT
SCIP_RETCODE SYMbuildDejavuGraphCheck(
   SCIP*                 scip,               /**< SCIP pointer */
   dejavu::static_graph* dejavugraph,        /**< pointer to hold dejavu graph being created */
   SYM_GRAPH*            G1,                 /**< first graph */
   SYM_GRAPH*            G2,                 /**< second graph */
   int*                  nnodes,             /**< pointer to store number of nodes in dejavu graph */
   int*                  nnodesfromG1,       /**< pointer to store number of nodes in dejavu graph arising from G1 */
   SCIP_Bool*            success             /**< pointer to store whether dejavugraph could be built */
   );


#ifdef __cplusplus
}
#endif

#endif
