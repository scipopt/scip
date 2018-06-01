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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_randnumgen.c
 * @brief  public methods for random numbers
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Robert Lion Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
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


#include "lpi/lpi.h"
#include "nlpi/exprinterpret.h"
#include "nlpi/nlpi.h"
#include "scip/benders.h"
#include "scip/benderscut.h"
#include "scip/branch.h"
#include "scip/branch_nodereopt.h"
#include "scip/clock.h"
#include "scip/compr.h"
#include "scip/concsolver.h"
#include "scip/concurrent.h"
#include "scip/conflict.h"
#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/cons_linear.h"
#include "scip/cutpool.h"
#include "scip/cuts.h"
#include "scip/debug.h"
#include "scip/def.h"
#include "scip/dialog.h"
#include "scip/dialog_default.h"
#include "scip/disp.h"
#include "scip/event.h"
#include "scip/heur.h"
#include "scip/heur_ofins.h"
#include "scip/heur_reoptsols.h"
#include "scip/heur_trivialnegation.h"
#include "scip/heuristics.h"
#include "scip/history.h"
#include "scip/implics.h"
#include "scip/interrupt.h"
#include "scip/lp.h"
#include "scip/mem.h"
#include "scip/message_default.h"
#include "scip/misc.h"
#include "scip/nlp.h"
#include "scip/nodesel.h"
#include "scip/paramset.h"
#include "scip/presol.h"
#include "scip/presolve.h"
#include "scip/pricer.h"
#include "scip/pricestore.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/prop.h"
#include "scip/reader.h"
#include "scip/relax.h"
#include "scip/reopt.h"
#include "scip/retcode.h"
#include "scip/scipbuildflags.h"
#include "scip/scipcoreplugins.h"
#include "scip/scipgithash.h"
#include "scip/sepa.h"
#include "scip/sepastore.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/solve.h"
#include "scip/stat.h"
#include "scip/syncstore.h"
#include "scip/table.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/visual.h"
#include "xml/xml.h"

#include "scip/scip_mem.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_var.h"

#include "scip/pub_message.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** creates and initializes a random number generator
 *
 *  @note The initial seed is changed using SCIPinitializeRandomSeed()
 */
SCIP_RETCODE SCIPcreateRandom(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN**     randnumgen,         /**< random number generator */
   unsigned int          initialseed,        /**< initial random seed */
   SCIP_Bool             useglobalseed       /**< should the supplied seed be initialized by SCIP's global seed shift? */
   )
{
   unsigned int modifiedseed;

   assert(scip != NULL);
   assert(randnumgen != NULL);

   if( useglobalseed )
      modifiedseed = SCIPinitializeRandomSeed(scip, initialseed);
   else
      modifiedseed = initialseed;

   SCIP_CALL( SCIPrandomCreate(randnumgen, SCIPblkmem(scip), modifiedseed) );

   return SCIP_OKAY;
}

/** frees a random number generator */
void SCIPfreeRandom(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN**     randnumgen          /**< random number generator */
   )
{
   assert(scip != NULL);
   assert(randnumgen != NULL);

   SCIPrandomFree(randnumgen, SCIPblkmem(scip));
}

/** initializes a random number generator with a given start seed
 *
 *  @note The seed is changed using SCIPinitializeRandomSeed()
 */
void SCIPsetRandomSeed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   unsigned int          seed                /**< new random seed */
   )
{
   unsigned int modifiedseed;

   assert(scip != NULL);
   assert(randnumgen != NULL);

   modifiedseed = SCIPinitializeRandomSeed(scip, seed);

   SCIPrandomSetSeed(randnumgen, modifiedseed);
}

/** modifies an initial seed value with the global shift of random seeds */
unsigned int SCIPinitializeRandomSeed(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int          initialseedvalue    /**< initial seed value to be modified */
   )
{
   assert(scip != NULL);

   return SCIPsetInitializeRandomSeed(scip->set, initialseedvalue);
}
