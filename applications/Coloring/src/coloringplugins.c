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

/**@file   coloringplugins.c
 * @brief  SCIP plugins for coloring
 * @author Gerald Gamrath
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "coloringplugins.h"
#include "scip/debug.h"

#define USEHEURS 0


/** includes default plugins for coloring into SCIP */
SCIP_RETCODE SCIPincludeColoringPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPincludeDialogDefaultBasic(scip) );

   SCIP_CALL( SCIPincludeConshdlrLinear(scip) ); /* linear must be first due to constraint upgrading */
   SCIP_CALL( SCIPincludeConshdlrIntegral(scip) );
   SCIP_CALL( SCIPincludeConshdlrSetppc(scip) );

   SCIP_CALL( SCIPincludeNodeselBfs(scip) );
   SCIP_CALL( SCIPincludeNodeselDfs(scip) );
   SCIP_CALL( SCIPincludeNodeselEstimate(scip) );
   SCIP_CALL( SCIPincludeNodeselHybridestim(scip) );
   SCIP_CALL( SCIPincludeNodeselRestartdfs(scip) );

   SCIP_CALL( SCIPincludeBranchruleAllfullstrong(scip) );
   SCIP_CALL( SCIPincludeBranchruleFullstrong(scip) );
   SCIP_CALL( SCIPincludeBranchruleInference(scip) );
   SCIP_CALL( SCIPincludeBranchruleMostinf(scip) );
   SCIP_CALL( SCIPincludeBranchruleLeastinf(scip) );
   SCIP_CALL( SCIPincludeBranchrulePscost(scip) );
   SCIP_CALL( SCIPincludeBranchruleRandom(scip) );
   SCIP_CALL( SCIPincludeBranchruleRelpscost(scip) );

   SCIP_CALL( SCIPincludeTableDefault(scip) );

   /** project plugins */
   SCIP_CALL( SCIPincludePricerColoring(scip) );

   SCIP_CALL( COLORincludeConshdlrStoreGraph(scip) );

   SCIP_CALL( SCIPincludeReaderCol(scip) );
   SCIP_CALL( SCIPincludeReaderCsol(scip) );

   SCIP_CALL( SCIPincludeBranchruleColoring(scip) );
   SCIP_CALL( SCIPincludeBranchruleStrongcoloring(scip) );

   SCIP_CALL( SCIPincludeHeurInit(scip) );

#if USEHEURS
   SCIP_CALL( SCIPincludeHeurActconsdiving(scip) );
   SCIP_CALL( SCIPincludeHeurCoefdiving(scip) );
   SCIP_CALL( SCIPincludeHeurCrossover(scip) );
   SCIP_CALL( SCIPincludeHeurDins(scip) );
   SCIP_CALL( SCIPincludeHeurFixandinfer(scip) );
   SCIP_CALL( SCIPincludeHeurFracdiving(scip) );
   SCIP_CALL( SCIPincludeHeurGuideddiving(scip) );
   SCIP_CALL( SCIPincludeHeurIntdiving(scip) );
   SCIP_CALL( SCIPincludeHeurIntshifting(scip) );
   SCIP_CALL( SCIPincludeHeurLinesearchdiving(scip) );
   SCIP_CALL( SCIPincludeHeurLocalbranching(scip) );
   SCIP_CALL( SCIPincludeHeurMutation(scip) );
   SCIP_CALL( SCIPincludeHeurObjpscostdiving(scip) );
   SCIP_CALL( SCIPincludeHeurOctane(scip) );
   SCIP_CALL( SCIPincludeHeurOneopt(scip) );
   SCIP_CALL( SCIPincludeHeurPscostdiving(scip) );
   SCIP_CALL( SCIPincludeHeurRens(scip) );
   SCIP_CALL( SCIPincludeHeurRins(scip) );
   SCIP_CALL( SCIPincludeHeurRootsoldiving(scip) );
   SCIP_CALL( SCIPincludeHeurRounding(scip) );
   SCIP_CALL( SCIPincludeHeurShifting(scip) );
   SCIP_CALL( SCIPincludeHeurSimplerounding(scip) );
   SCIP_CALL( SCIPincludeHeurTrivial(scip) );
   SCIP_CALL( SCIPincludeHeurTrySol(scip) );
   SCIP_CALL( SCIPincludeHeurTwoopt(scip) );
   SCIP_CALL( SCIPincludeHeurUndercover(scip) );
   SCIP_CALL( SCIPincludeHeurVeclendiving(scip) );
   SCIP_CALL( SCIPincludeHeurZirounding(scip) );
#endif

   SCIP_CALL( SCIPincludeDispDefault(scip) );

   return SCIP_OKAY;
}
