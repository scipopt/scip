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

/**@file   scip.h
 * @ingroup PUBLICCOREAPI
 * @brief  SCIP callable library
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_H__
#define __SCIP_SCIP_H__


#include <stdio.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_clock.h"
#include "scip/type_misc.h"
#include "scip/type_timing.h"
#include "scip/type_paramset.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_nlp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_tree.h"
#include "scip/type_scip.h"

#include "scip/type_bandit.h"
#include "scip/type_branch.h"
#include "scip/type_conflict.h"
#include "scip/type_cons.h"
#include "scip/type_dialog.h"
#include "scip/type_disp.h"
#include "scip/type_heur.h"
#include "scip/type_compr.h"
#include "scip/type_history.h"
#include "scip/type_nodesel.h"
#include "scip/type_presol.h"
#include "scip/type_pricer.h"
#include "scip/type_reader.h"
#include "scip/type_relax.h"
#include "scip/type_sepa.h"
#include "scip/type_table.h"
#include "scip/type_prop.h"
#include "nlpi/type_nlpi.h"
#include "scip/type_concsolver.h"
#include "scip/type_syncstore.h"
#include "scip/type_benders.h"
#include "scip/type_benderscut.h"

/* include public interfaces, s.t. the user only needs to include scip.h */
#include "scip/pub_bandit.h"
#include "scip/pub_benders.h"
#include "scip/pub_benderscut.h"
#include "scip/pub_branch.h"
#include "scip/pub_conflict.h"
#include "scip/pub_cons.h"
#include "scip/pub_cutpool.h"
#include "scip/pub_dialog.h"
#include "scip/pub_disp.h"
#include "scip/pub_event.h"
#include "scip/pub_fileio.h"
#include "scip/pub_heur.h"
#include "scip/pub_compr.h"
#include "scip/pub_history.h"
#include "scip/pub_implics.h"
#include "scip/pub_lp.h"
#include "scip/pub_nlp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_nodesel.h"
#include "scip/pub_paramset.h"
#include "scip/pub_presol.h"
#include "scip/pub_pricer.h"
#include "scip/pub_reader.h"
#include "scip/pub_relax.h"
#include "scip/pub_sepa.h"
#include "scip/pub_prop.h"
#include "scip/pub_sol.h"
#include "scip/pub_table.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "lpi/lpi.h"
#include "nlpi/pub_expr.h"

/* include global presolving, cuts, and heuristics methods */
#include "scip/presolve.h"
#include "scip/cuts.h"
#include "scip/heuristics.h"

/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons).
 * Additionally, the internal "set.h" is included, such that the defines in set.h are
 * available in optimized mode.
 */
#ifdef NDEBUG
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"
#include "scip/var.h"
#include "scip/cons.h"
#include "scip/solve.h"
#include "scip/debug.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * miscellaneous methods
 */





/*
 * general SCIP methods
 */


/*
 * debug solution methods
 */



/*
 * message output methods
 */





/*
 * SCIP copy methods
 */


/*
 * parameter settings
 */



/*
 * SCIP user functionality methods: managing plugins
 */






















/*
 * user interactive dialog methods
 */


/*
 * global problem methods
 */





/*
 * local subproblem methods
 */





/*
 * solve methods
 */





/*
 * variable methods
 */





/*
 * conflict analysis methods
 */





/*
 * constraint methods
 */





/*
 * LP methods
 */




/*
 * LP column methods
 */





/*
 * LP row methods
 */



/*
 * NLP methods
 */



/*
 * NLP diving methods
 */



/*
 * NLP nonlinear row methods
 */






/*
 * nonlinear methods
 */





/*
 * cutting plane methods
 */





/*
 * LP diving methods
 */





/*
 * probing methods
 */


/*
 * branching methods
 */





/*
 * primal solutions
 */





/*
 * event methods
 */





/*
 * tree methods
 */



/*
 * parallel interface methods
 */




/*
 * reoptimization methods
 */


/*
 * statistic methods
 */


/*
 * timing methods
 */



/*
 * numeric values and comparisons
 */





/*
 * memory management
 */




/*
 * dynamic arrays
 */






#ifdef __cplusplus
}
#endif

#endif
