/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: scipdefplugins.c,v 1.40 2005/06/23 16:30:04 bzfberth Exp $"

/**@file   scipdefplugins.c
 * @brief  default SCIP plugins
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/branch_allfullstrong.h"
#include "scip/branch_fullstrong.h"
#include "scip/branch_inference.h"
#include "scip/branch_mostinf.h"
#include "scip/branch_leastinf.h"
#include "scip/branch_pscost.h"
#include "scip/branch_relpscost.h"
#include "scip/cons_and.h"
#include "scip/cons_conjunction.h"
#include "scip/cons_integral.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_or.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_xor.h"
#include "scip/dialog_default.h"
#include "scip/disp_default.h"
#include "scip/heur_coefdiving.h"
#include "scip/heur_feaspump.h"
#include "scip/heur_fixandinfer.h"
#include "scip/heur_fracdiving.h"
#include "scip/heur_guideddiving.h"
#include "scip/heur_linesearchdiving.h"
#include "scip/heur_objpscostdiving.h"
#include "scip/heur_octane.h"
#include "scip/heur_pscostdiving.h"
#include "scip/heur_rootsoldiving.h"
#include "scip/heur_rounding.h"
#include "scip/heur_simplerounding.h"
#include "scip/nodesel_bfs.h"
#include "scip/nodesel_dfs.h"
#include "scip/nodesel_restartdfs.h"
#include "scip/presol_dualfix.h"
#include "scip/presol_probing.h"
#include "scip/presol_trivial.h"
#include "scip/prop_pseudoobj.h"
#include "scip/reader_cnf.h"
#include "scip/reader_mps.h"
#include "scip/sepa_clique.h"
#include "scip/sepa_cmir.h"
#include "scip/sepa_gomory.h"
#include "scip/sepa_impliedbounds.h"
#include "scip/sepa_intobj.h"
#include "scip/sepa_strongcg.h"
#include "scip/scipdefplugins.h"
#include "scip/debug.h"

#if 0
#include "scip/cons_binpack.h"
#include "scip/cons_eqknapsack.h"
#include "scip/cons_invarknapsack.h"
#endif



/** includes default SCIP plugins into SCIP */
RETCODE SCIPincludeDefaultPlugins(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( SCIPincludeConshdlrLinear(scip) );
   CHECK_OKAY( SCIPincludeConshdlrAnd(scip) );
   CHECK_OKAY( SCIPincludeConshdlrConjunction(scip) );
   CHECK_OKAY( SCIPincludeConshdlrIntegral(scip) );
   CHECK_OKAY( SCIPincludeConshdlrKnapsack(scip) );
   CHECK_OKAY( SCIPincludeConshdlrLogicor(scip) );
   CHECK_OKAY( SCIPincludeConshdlrOr(scip) );
   CHECK_OKAY( SCIPincludeConshdlrSetppc(scip) );
   CHECK_OKAY( SCIPincludeConshdlrVarbound(scip) );
   CHECK_OKAY( SCIPincludeConshdlrXor(scip) );

#if 0
   CHECK_OKAY( SCIPincludeConshdlrEqknapsack(scip) );
   CHECK_OKAY( SCIPincludeConshdlrInvarknapsack(scip) );
   CHECK_OKAY( SCIPincludeConshdlrBinpack(scip) );
#endif

   CHECK_OKAY( SCIPincludeReaderCnf(scip) );
   CHECK_OKAY( SCIPincludeReaderMps(scip) );
   CHECK_OKAY( SCIPincludePresolDualfix(scip) );
   CHECK_OKAY( SCIPincludePresolProbing(scip) );
   CHECK_OKAY( SCIPincludePresolTrivial(scip) );
   CHECK_OKAY( SCIPincludeNodeselBfs(scip) );
   CHECK_OKAY( SCIPincludeNodeselDfs(scip) );
   CHECK_OKAY( SCIPincludeNodeselRestartdfs(scip) );
   CHECK_OKAY( SCIPincludeBranchruleAllfullstrong(scip) );
   CHECK_OKAY( SCIPincludeBranchruleFullstrong(scip) );
   CHECK_OKAY( SCIPincludeBranchruleInference(scip) );
   CHECK_OKAY( SCIPincludeBranchruleMostinf(scip) );
   CHECK_OKAY( SCIPincludeBranchruleLeastinf(scip) );
   CHECK_OKAY( SCIPincludeBranchrulePscost(scip) );
   CHECK_OKAY( SCIPincludeBranchruleRelpscost(scip) );
   CHECK_OKAY( SCIPincludeHeurCoefdiving(scip) );
   CHECK_OKAY( SCIPincludeHeurFeaspump(scip) );
   CHECK_OKAY( SCIPincludeHeurFixandinfer(scip) );
   CHECK_OKAY( SCIPincludeHeurFracdiving(scip) );
   CHECK_OKAY( SCIPincludeHeurGuideddiving(scip) );
   CHECK_OKAY( SCIPincludeHeurLinesearchdiving(scip) );
   CHECK_OKAY( SCIPincludeHeurObjpscostdiving(scip) );
   CHECK_OKAY( SCIPincludeHeurOctane(scip) );
   CHECK_OKAY( SCIPincludeHeurPscostdiving(scip) );
   CHECK_OKAY( SCIPincludeHeurRootsoldiving(scip) );
   CHECK_OKAY( SCIPincludeHeurRounding(scip) );
   CHECK_OKAY( SCIPincludeHeurSimplerounding(scip) );
   CHECK_OKAY( SCIPincludePropPseudoobj(scip) );
   CHECK_OKAY( SCIPincludeSepaClique(scip) );
   CHECK_OKAY( SCIPincludeSepaCmir(scip) );
   CHECK_OKAY( SCIPincludeSepaGomory(scip) );
   CHECK_OKAY( SCIPincludeSepaImpliedbounds(scip) );
   CHECK_OKAY( SCIPincludeSepaIntobj(scip) );
   CHECK_OKAY( SCIPincludeSepaStrongcg(scip) );
   CHECK_OKAY( SCIPincludeDispDefault(scip) );
   CHECK_OKAY( SCIPincludeDialogDefault(scip) );

   CHECK_OKAY( SCIPdebugIncludeProp(scip) ); /*lint !e506 !e774*/

   return SCIP_OKAY;
}
