/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: scipdefplugins.c,v 1.25 2004/07/19 15:52:00 bzfpfets Exp $"

/**@file   scipdefplugins.c
 * @brief  default SCIP plugins
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "branch_allfullstrong.h"
#include "branch_fullstrong.h"
#include "branch_inference.h"
#include "branch_mostinf.h"
#include "branch_leastinf.h"
#include "branch_relpscost.h"
#include "cons_and.h"
#include "cons_binpack.h"
#include "cons_conjunction.h"
#include "cons_eqknapsack.h"
#include "cons_integral.h"
#include "cons_invarknapsack.h"
#include "cons_knapsack.h"
#include "cons_linear.h"
#include "cons_logicor.h"
#include "cons_or.h"
#include "cons_setppc.h"
#include "cons_varbound.h"
#include "dialog_default.h"
#include "disp_default.h"
#include "heur_coefdiving.h"
#include "heur_feaspump.h"
#include "heur_fracdiving.h"
#include "heur_linesearchdiving.h"
#include "heur_objpscostdiving.h"
#include "heur_pscostdiving.h"
#include "heur_rootsoldiving.h"
#include "heur_rounding.h"
#include "heur_simplerounding.h"
#include "nodesel_bfs.h"
#include "nodesel_dfs.h"
#include "nodesel_restartdfs.h"
#include "presol_dualfix.h"
#include "presol_trivial.h"
#include "reader_cnf.h"
#include "reader_mps.h"
#include "sepa_cmir.h"
#include "sepa_gomory.h"
#include "sepa_intobj.h"
#include "scipdefplugins.h"




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

#if 0
   CHECK_OKAY( SCIPincludeConshdlrEqknapsack(scip) );
   CHECK_OKAY( SCIPincludeConshdlrInvarknapsack(scip) );
   CHECK_OKAY( SCIPincludeConshdlrBinpack(scip) );
#endif

   CHECK_OKAY( SCIPincludeReaderCnf(scip) );
   CHECK_OKAY( SCIPincludeReaderMps(scip) );
   CHECK_OKAY( SCIPincludePresolDualfix(scip) );
   CHECK_OKAY( SCIPincludePresolTrivial(scip) );
   CHECK_OKAY( SCIPincludeNodeselBfs(scip) );
   CHECK_OKAY( SCIPincludeNodeselDfs(scip) );
   CHECK_OKAY( SCIPincludeNodeselRestartdfs(scip) );
   CHECK_OKAY( SCIPincludeBranchruleAllfullstrong(scip) );
   CHECK_OKAY( SCIPincludeBranchruleFullstrong(scip) );
   CHECK_OKAY( SCIPincludeBranchruleInference(scip) );
   CHECK_OKAY( SCIPincludeBranchruleMostinf(scip) );
   CHECK_OKAY( SCIPincludeBranchruleLeastinf(scip) );
   CHECK_OKAY( SCIPincludeBranchruleRelpscost(scip) );
   CHECK_OKAY( SCIPincludeHeurCoefdiving(scip) );
   CHECK_OKAY( SCIPincludeHeurFeaspump(scip) );
   CHECK_OKAY( SCIPincludeHeurFracdiving(scip) );
   CHECK_OKAY( SCIPincludeHeurLinesearchdiving(scip) );
   CHECK_OKAY( SCIPincludeHeurObjpscostdiving(scip) );
   CHECK_OKAY( SCIPincludeHeurPscostdiving(scip) );
   CHECK_OKAY( SCIPincludeHeurRootsoldiving(scip) );
   CHECK_OKAY( SCIPincludeHeurRounding(scip) );
   CHECK_OKAY( SCIPincludeHeurSimplerounding(scip) );
   CHECK_OKAY( SCIPincludeSepaCmir(scip) );
   CHECK_OKAY( SCIPincludeSepaGomory(scip) );
   CHECK_OKAY( SCIPincludeSepaIntobj(scip) );
   CHECK_OKAY( SCIPincludeDispDefault(scip) );
   CHECK_OKAY( SCIPincludeDialogDefault(scip) );

   return SCIP_OKAY;
}
