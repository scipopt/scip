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

/**@file   build_sassy_graph.h
 * @brief  methods to build sassy graph for symmetry detection
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BUILD_SASSY_GRAPH_H_
#define __SCIP_BUILD_SASSY_GRAPH_H_

#include "scip/scip.h"

/* include sassy */
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4189)  // local variable is initialized but not referenced
# pragma warning(disable: 4388)  // compare signed and unsigned expression
# pragma warning(disable: 4456)  // shadowed variable
# pragma warning(disable: 4430)  // missing type specifier
#endif

/* the actual include */
#include <sassy/graph.h>

#ifdef __GNUC__
#pragma GCC diagnostic warning "-Wunused-but-set-variable"
#pragma GCC diagnostic warning "-Wsign-compare"
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wshadow"
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
SCIP_RETCODE SYMbuildSassyGraph(
   SCIP*                 scip,               /**< SCIP pointer */
   sassy::static_graph*  sassygraph,         /**< pointer to hold sassy graph being created */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_Bool*            success             /**< pointer to store whether sassygraph could be built */
   );


/** returns whether two given graphs are identical */
SCIP_EXPORT
SCIP_RETCODE SYMbuildSassyGraphCheck(
   SCIP*                 scip,               /**< SCIP pointer */
   sassy::static_graph*  sassygraph,         /**< pointer to hold sassy graph being created */
   SYM_GRAPH*            G1,                 /**< first graph */
   SYM_GRAPH*            G2,                 /**< second graph */
   int*                  nnodes,             /**< pointer to store number of nodes in sassy graph */
   int*                  nnodesfromG1,       /**< pointer to store number of nodes in sassy graph arising from G1 */
   SCIP_Bool*            success             /**< pointer to store whether sassygraph could be built */
   );


#ifdef __cplusplus
}
#endif

#endif
