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

/**@file   scip_presol.c
 * @brief  public methods for presolving plugins
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

#include "scip/scip_presol.h"

#include "scip/pub_message.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** creates a presolver and includes it in SCIP.
 *
 *  @note method has all presolver callbacks as arguments and is thus changed every time a new
 *        callback is added
 *        in future releases; consider using SCIPincludePresolBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludePresol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of presolver */
   const char*           desc,               /**< description of presolver */
   int                   priority,           /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
   int                   maxrounds,          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
   SCIP_PRESOLTIMING     timing,             /**< timing mask of the presolver */
   SCIP_DECL_PRESOLCOPY  ((*presolcopy)),    /**< copy method of presolver or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PRESOLFREE  ((*presolfree)),    /**< destructor of presolver to free user data (called when SCIP is exiting) */
   SCIP_DECL_PRESOLINIT  ((*presolinit)),    /**< initialization method of presolver (called after problem was transformed) */
   SCIP_DECL_PRESOLEXIT  ((*presolexit)),    /**< deinitialization method of presolver (called before transformed problem is freed) */
   SCIP_DECL_PRESOLINITPRE((*presolinitpre)),/**< presolving initialization method of presolver (called when presolving is about to begin) */
   SCIP_DECL_PRESOLEXITPRE((*presolexitpre)),/**< presolving deinitialization method of presolver (called after presolving has been finished) */
   SCIP_DECL_PRESOLEXEC  ((*presolexec)),    /**< execution method of presolver */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   SCIP_PRESOL* presol;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludePresol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether presolver is already present */
   if( SCIPfindPresol(scip, name) != NULL )
   {
      SCIPerrorMessage("presolver <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPpresolCreate(&presol, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         maxrounds, timing, presolcopy,
         presolfree, presolinit, presolexit, presolinitpre, presolexitpre, presolexec, presoldata) );
   SCIP_CALL( SCIPsetIncludePresol(scip->set, presol) );

   return SCIP_OKAY;
}

/** creates a presolver and includes it in SCIP with its fundamental callback. All non-fundamental (or optional)
 *  callbacks as, e.g., init and exit callbacks, will be set to NULL. Optional callbacks can be set via specific setter
 *  functions. These are SCIPsetPresolCopy(), SCIPsetPresolFree(), SCIPsetPresolInit(), SCIPsetPresolExit(),
 *  SCIPsetPresolInitpre(), and SCIPsetPresolExitPre().
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludePresol() instead
 */
SCIP_RETCODE SCIPincludePresolBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL**         presolptr,          /**< reference to presolver, or NULL */
   const char*           name,               /**< name of presolver */
   const char*           desc,               /**< description of presolver */
   int                   priority,           /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
   int                   maxrounds,          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
   SCIP_PRESOLTIMING     timing,             /**< timing mask of the presolver */
   SCIP_DECL_PRESOLEXEC  ((*presolexec)),    /**< execution method of presolver */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   SCIP_PRESOL* presol;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludePresolBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether presolver is already present */
   if( SCIPfindPresol(scip, name) != NULL )
   {
      SCIPerrorMessage("presolver <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPpresolCreate(&presol, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority, maxrounds, timing,
         NULL,
         NULL, NULL, NULL, NULL, NULL, presolexec, presoldata) );
   SCIP_CALL( SCIPsetIncludePresol(scip->set, presol) );

   if( presolptr != NULL )
      *presolptr = presol;

   return SCIP_OKAY;
}

/** sets copy method of presolver */
SCIP_RETCODE SCIPsetPresolCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLCOPY ((*presolcopy))      /**< copy method of presolver or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPresolCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(presol != NULL);

   SCIPpresolSetCopy(presol, presolcopy);

   return SCIP_OKAY;
}

/** sets destructor method of presolver */
SCIP_RETCODE SCIPsetPresolFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLFREE ((*presolfree))      /**< destructor of presolver */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPresolFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(presol != NULL);

   SCIPpresolSetFree(presol, presolfree);

   return SCIP_OKAY;
}

/** sets initialization method of presolver */
SCIP_RETCODE SCIPsetPresolInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLINIT  ((*presolinit))     /**< initialize presolver */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPresolInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(presol != NULL);

   SCIPpresolSetInit(presol, presolinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of presolver */
SCIP_RETCODE SCIPsetPresolExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLEXIT  ((*presolexit))     /**< deinitialize presolver */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPresolExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(presol != NULL);

   SCIPpresolSetExit(presol, presolexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of presolver */
SCIP_RETCODE SCIPsetPresolInitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLINITPRE ((*presolinitpre))/**< solving process initialization method of presolver */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPresolInitpre", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(presol != NULL);

   SCIPpresolSetInitpre(presol, presolinitpre);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of presolver */
SCIP_RETCODE SCIPsetPresolExitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLEXITPRE ((*presolexitpre))/**< solving process deinitialization method of presolver */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPresolExitpre", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(presol != NULL);

   SCIPpresolSetExitpre(presol, presolexitpre);

   return SCIP_OKAY;
}

/** returns the presolver of the given name, or NULL if not existing */
SCIP_PRESOL* SCIPfindPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of presolver */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindPresol(scip->set, name);
}

/** returns the array of currently available presolvers */
SCIP_PRESOL** SCIPgetPresols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortPresols(scip->set);

   return scip->set->presols;
}

/** returns the number of currently available presolvers */
int SCIPgetNPresols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->npresols;
}

/** sets the priority of a presolver */
SCIP_RETCODE SCIPsetPresolPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   int                   priority            /**< new priority of the presolver */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPpresolSetPriority(presol, scip->set, priority);

   return SCIP_OKAY;
}
