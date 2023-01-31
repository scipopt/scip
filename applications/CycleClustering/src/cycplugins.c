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

/**@file   cycplugins.c
 * @brief  SCIP plugins for cycle clustering of markov state models
 * @author Leon Eifler
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "cycplugins.h"

#include "scip/scipdefplugins.h"
#include "sepa_edge.h"
#include "sepa_subtour.h"
#include "sepa_partition.h"
#include "heur_fuzzyround.h"
#include "branch_multinode.h"
#include "heur_cyckerlin.h"
#include "heur_redsize.h"
#include "event_newsol.h"

/** includes default plugins for cycle clustering into SCIP */
SCIP_RETCODE SCIPincludeCycPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* Default and reader */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPincludeReaderCyc(scip) );

   /* Heuristics */
   SCIP_CALL( SCIPincludeHeurCycKerlin(scip) );
   SCIP_CALL( SCIPincludeHeurFuzzyround(scip) );
   SCIP_CALL( SCIPincludeHeurCycGreedy(scip) );
   SCIP_CALL( SCIPincludeHeurRedsize(scip) );

   /* Separators */
   SCIP_CALL( SCIPincludeSepaEdge(scip) );
   SCIP_CALL( SCIPincludeSepaPartition(scip) );
   SCIP_CALL( SCIPincludeSepaSubtour(scip) );

   /* Branching rule */
   SCIP_CALL( SCIPincludeBranchruleMultinode(scip) );

   /* Event handler that reruns exchange heuristic for new solutions */
   SCIP_CALL( SCIPincludeEventHdlrNewsol(scip) );

   return SCIP_OKAY;
}
