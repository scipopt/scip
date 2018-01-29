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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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

   /* Separators */
   SCIP_CALL( SCIPincludeSepaEdge(scip) );
   SCIP_CALL( SCIPincludeSepaPartition(scip) );
   SCIP_CALL( SCIPincludeSepaSubtour(scip) );

   /* Branching rule */
   SCIP_CALL( SCIPincludeBranchruleMultinode(scip) );

   return SCIP_OKAY;
}
