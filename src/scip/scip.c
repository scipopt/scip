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

/**@file   scip.c
 * @brief  SCIP callable library
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIPcheckStage() calls, use bit flags instead of the SCIP_Bool parameters
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 * @todo When making an interface change to SCIPcreateProb(), add an indication whether it is known that the objective
 *       function is integral or whether this is not known. This avoids a manual call of SCIPsetObjIntegral(). Moreover,
 *       then the detection could be performed incrementally, whenever a variable is added.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <ctype.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif

#ifdef WITH_ZLIB
#include <zlib.h>
#endif

#include "scip/def.h"
#include "scip/retcode.h"
#include "scip/set.h"
#include "scip/paramset.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/visual.h"
#include "scip/interrupt.h"
#include "scip/mem.h"
#include "scip/misc.h"
#include "scip/history.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/nlp.h"
#include "scip/var.h"
#include "scip/implics.h"
#include "scip/prob.h"
#include "scip/sol.h"
#include "scip/primal.h"
#include "scip/reopt.h"
#include "scip/tree.h"
#include "scip/pricestore.h"
#include "scip/sepastore.h"
#include "scip/conflictstore.h"
#include "scip/syncstore.h"
#include "scip/cutpool.h"
#include "scip/solve.h"
#include "scip/scipbuildflags.h"
#include "scip/scipgithash.h"
#include "scip/cuts.h"
#include "scip/scip.h"
#include "lpi/lpi.h"

#include "scip/scipcoreplugins.h"
#include "scip/branch.h"
#include "scip/conflict.h"
#include "scip/cons.h"
#include "scip/dialog.h"
#include "scip/disp.h"
#include "scip/table.h"
#include "scip/heur.h"
#include "scip/concsolver.h"
#include "scip/compr.h"
#include "scip/nodesel.h"
#include "scip/reader.h"
#include "scip/presol.h"
#include "scip/pricer.h"
#include "scip/relax.h"
#include "scip/sepa.h"
#include "scip/prop.h"
#include "nlpi/nlpi.h"
#include "nlpi/exprinterpret.h"
#include "scip/debug.h"
#include "scip/dialog_default.h"
#include "scip/message_default.h"
#include "scip/syncstore.h"
#include "scip/concurrent.h"
#include "scip/benders.h"
#include "scip/benderscut.h"
#include "xml/xml.h"

/* We include the linear constraint handler to be able to copy a (multi)aggregation of variables (to a linear constraint).
 * The better way would be to handle the distinction between original and transformed variables via a flag 'isoriginal'
 * in the variable data structure. This would allow to have (multi)aggregated variables in the original problem.
 *
 * A second reason for including the linear constraint handler is for copying cuts to linear constraints.
 */
#include "scip/cons_linear.h"

/* We need to include the branching and the heurtistics for reoptimization after creating the reoptimization because we
 * do not want to use these plugins by default if reoptimization is disabled. */
#include "scip/branch_nodereopt.h"
#include "scip/heur_reoptsols.h"
#include "scip/heur_trivialnegation.h"
#include "scip/heur_ofins.h"

/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/*
 * miscellaneous methods
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







#undef SCIPgetSepaMinEfficacy





/* new callback/method setter methods */

/* new callback/method setter methods */








/** method to call, when the priority of an NLPI was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdNlpiPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetSetPriorityNlpi() to mark the nlpis unsorted */
   SCIP_CALL( SCIPsetNlpiPriority(scip, (SCIP_NLPI*)paramdata, SCIPparamGetInt(param)) );

   return SCIP_OKAY;
}



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



/* some simple variable functions implemented as defines */
#undef SCIPcomputeVarLbGlobal
#undef SCIPcomputeVarUbGlobal
#undef SCIPcomputeVarLbLocal
#undef SCIPcomputeVarUbLocal


/* calculate clique partition for a maximal amount of comparisons on variables due to expensive algorithm
 * @todo: check for a good value, maybe it's better to check parts of variables
 */
#define MAXNCLIQUEVARSCOMP 1000000




/*
 * conflict analysis methods
 */


/*
 * constraint methods
 */

#undef SCIPmarkConsPropagate


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

/**@name NLP Diving Methods */
/**@{ */

/**@} */


/*
 * NLP nonlinear row methods
 */


/**@name Expression tree methods */
/**@{ */

/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) (val >= infty1 ? infty2 : val)


#undef infty2infty

/**@} */

/*
 * nonlinear methods
 */

/**@name Nonlinear Methods */
/**@{ */



/**@} */

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
 * statistic methods
 */




/** comparison method for statistics tables */
static
SCIP_DECL_SORTPTRCOMP(tablePosComp)
{  /*lint --e{715}*/
   return (SCIPtableGetPosition((SCIP_TABLE*)elem1) - (SCIPtableGetPosition((SCIP_TABLE*)elem2)));
}





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
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPinfinity
#undef SCIPisInfinity
#undef SCIPisEQ
#undef SCIPisLT
#undef SCIPisLE
#undef SCIPisGT
#undef SCIPisGE
#undef SCIPisZero
#undef SCIPisPositive
#undef SCIPisNegative
#undef SCIPisIntegral
#undef SCIPisScalingIntegral
#undef SCIPisFracIntegral
#undef SCIPfloor
#undef SCIPceil
#undef SCIPround
#undef SCIPfrac
#undef SCIPisSumEQ
#undef SCIPisSumLT
#undef SCIPisSumLE
#undef SCIPisSumGT
#undef SCIPisSumGE
#undef SCIPisSumZero
#undef SCIPisSumPositive
#undef SCIPisSumNegative
#undef SCIPisFeasEQ
#undef SCIPisFeasLT
#undef SCIPisFeasLE
#undef SCIPisFeasGT
#undef SCIPisFeasGE
#undef SCIPisFeasZero
#undef SCIPisFeasPositive
#undef SCIPisFeasNegative
#undef SCIPisFeasIntegral
#undef SCIPisFeasFracIntegral
#undef SCIPfeasFloor
#undef SCIPfeasCeil
#undef SCIPfeasRound
#undef SCIPfeasFrac
#undef SCIPisDualfeasEQ
#undef SCIPisDualfeasLT
#undef SCIPisDualfeasLE
#undef SCIPisDualfeasGT
#undef SCIPisDualfeasGE
#undef SCIPisDualfeasZero
#undef SCIPisDualfeasPositive
#undef SCIPisDualfeasNegative
#undef SCIPisDualfeasIntegral
#undef SCIPisDualfeasFracIntegral
#undef SCIPdualfeasFloor
#undef SCIPdualfeasCeil
#undef SCIPdualfeasRound
#undef SCIPdualfeasFrac
#undef SCIPisLbBetter
#undef SCIPisUbBetter
#undef SCIPisRelEQ
#undef SCIPisRelLT
#undef SCIPisRelLE
#undef SCIPisRelGT
#undef SCIPisRelGE
#undef SCIPisSumRelEQ
#undef SCIPisSumRelLT
#undef SCIPisSumRelLE
#undef SCIPisSumRelGT
#undef SCIPisSumRelGE
#undef SCIPconvertRealToInt
#undef SCIPconvertRealToLongint
#undef SCIPisUpdateUnreliable
#undef SCIPisHugeValue
#undef SCIPgetHugeValue




