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

/**@file   scip_disp.c
 * @brief  public methods for display handler plugins
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

#include "scip/scip_disp.h"

#include "scip/pub_message.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** creates a display column and includes it in SCIP */
SCIP_RETCODE SCIPincludeDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of display column */
   const char*           desc,               /**< description of display column */
   const char*           header,             /**< head line of display column */
   SCIP_DISPSTATUS       dispstatus,         /**< display activation status of display column */
   SCIP_DECL_DISPCOPY    ((*dispcopy)),      /**< copy method of display column or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_DISPFREE    ((*dispfree)),      /**< destructor of display column */
   SCIP_DECL_DISPINIT    ((*dispinit)),      /**< initialize display column */
   SCIP_DECL_DISPEXIT    ((*dispexit)),      /**< deinitialize display column */
   SCIP_DECL_DISPINITSOL ((*dispinitsol)),   /**< solving process initialization method of display column */
   SCIP_DECL_DISPEXITSOL ((*dispexitsol)),   /**< solving process deinitialization method of display column */
   SCIP_DECL_DISPOUTPUT  ((*dispoutput)),    /**< output method */
   SCIP_DISPDATA*        dispdata,           /**< display column data */
   int                   width,              /**< width of display column (no. of chars used) */
   int                   priority,           /**< priority of display column */
   int                   position,           /**< relative position of display column */
   SCIP_Bool             stripline           /**< should the column be separated with a line from its right neighbor? */
   )
{
   SCIP_DISP* disp;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeDisp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether display column is already present */
   if( SCIPfindDisp(scip, name) != NULL )
   {
      SCIPerrorMessage("display column <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPdispCreate(&disp, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, header, dispstatus,
         dispcopy,
         dispfree, dispinit, dispexit, dispinitsol, dispexitsol, dispoutput, dispdata,
         width, priority, position, stripline) );
   SCIP_CALL( SCIPsetIncludeDisp(scip->set, disp) );

   return SCIP_OKAY;
}

/** returns the display column of the given name, or NULL if not existing */
SCIP_DISP* SCIPfindDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of display column */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindDisp(scip->set, name);
}

/** returns the array of currently available display columns */
SCIP_DISP** SCIPgetDisps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->disps;
}

/** returns the number of currently available display columns */
int SCIPgetNDisps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->ndisps;
}

/** automatically selects display columns for being shown w.r.t. the display width parameter */
SCIP_RETCODE SCIPautoselectDisps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPdispAutoActivate(scip->set) );

   return SCIP_OKAY;
}

/** changes the display column mode */
void SCIPchgDispMode(
   SCIP_DISP*            disp,               /**< display column */
   SCIP_DISPMODE         mode                /**< the display column mode */
   )
{
   assert(disp != NULL);

   SCIPdispChgMode(disp, mode);
}
