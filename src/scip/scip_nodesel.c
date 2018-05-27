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

/**@file   scip_nodesel.c
 * @brief  public methods for node selector plugins
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

#include "scip/scip_nodesel.h"

#include "scip/pub_message.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** creates a node selector and includes it in SCIP.
 *
 *  @note method has all node selector callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludeNodeselBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of node selector */
   const char*           desc,               /**< description of node selector */
   int                   stdpriority,        /**< priority of the node selector in standard mode */
   int                   memsavepriority,    /**< priority of the node selector in memory saving mode */
   SCIP_DECL_NODESELCOPY ((*nodeselcopy)),   /**< copy method of node selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_NODESELFREE ((*nodeselfree)),   /**< destructor of node selector */
   SCIP_DECL_NODESELINIT ((*nodeselinit)),   /**< initialize node selector */
   SCIP_DECL_NODESELEXIT ((*nodeselexit)),   /**< deinitialize node selector */
   SCIP_DECL_NODESELINITSOL((*nodeselinitsol)),/**< solving process initialization method of node selector */
   SCIP_DECL_NODESELEXITSOL((*nodeselexitsol)),/**< solving process deinitialization method of node selector */
   SCIP_DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   SCIP_DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   SCIP_NODESELDATA*     nodeseldata         /**< node selector data */
   )
{
   SCIP_NODESEL* nodesel;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeNodesel", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether node selector is already present */
   if( SCIPfindNodesel(scip, name) != NULL )
   {
      SCIPerrorMessage("node selector <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPnodeselCreate(&nodesel, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, stdpriority, memsavepriority,
         nodeselcopy, nodeselfree, nodeselinit, nodeselexit, nodeselinitsol, nodeselexitsol,
         nodeselselect, nodeselcomp, nodeseldata) );
   SCIP_CALL( SCIPsetIncludeNodesel(scip->set, nodesel) );

   return SCIP_OKAY;
}

/** Creates a node selector and includes it in SCIP with its most fundamental callbacks. All non-fundamental
 *  (or optional) callbacks as, e.g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetNodeselCopy(), SCIPsetNodeselFree(),
 *  SCIPsetNodeselInit(), SCIPsetNodeselExit(), SCIPsetNodeselInitsol(), and SCIPsetNodeselExitsol()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeNodesel() instead
 */
SCIP_RETCODE SCIPincludeNodeselBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL**        nodesel,            /**< reference to a node selector, or NULL */
   const char*           name,               /**< name of node selector */
   const char*           desc,               /**< description of node selector */
   int                   stdpriority,        /**< priority of the node selector in standard mode */
   int                   memsavepriority,    /**< priority of the node selector in memory saving mode */
   SCIP_DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   SCIP_DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   SCIP_NODESELDATA*     nodeseldata         /**< node selector data */
   )
{
   SCIP_NODESEL* nodeselptr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeNodeselBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether node selector is already present */
   if( SCIPfindNodesel(scip, name) != NULL )
   {
      SCIPerrorMessage("node selector <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPnodeselCreate(&nodeselptr, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, stdpriority, memsavepriority,
         NULL, NULL, NULL, NULL, NULL, NULL,
         nodeselselect, nodeselcomp, nodeseldata) );
   SCIP_CALL( SCIPsetIncludeNodesel(scip->set, nodeselptr) );

   if( nodesel != NULL )
      *nodesel = nodeselptr;

   return SCIP_OKAY;
}

/** sets copy method of node selector */
SCIP_RETCODE SCIPsetNodeselCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELCOPY ((*nodeselcopy))    /**< copy method of node selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNodeselCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nodesel != NULL);

   SCIPnodeselSetCopy(nodesel, nodeselcopy);

   return SCIP_OKAY;
}

/** sets destructor method of node selector */
SCIP_RETCODE SCIPsetNodeselFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELFREE ((*nodeselfree))    /**< destructor of node selector */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNodeselFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nodesel != NULL);

   SCIPnodeselSetFree(nodesel, nodeselfree);

   return SCIP_OKAY;
}

/** sets initialization method of node selector */
SCIP_RETCODE SCIPsetNodeselInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELINIT ((*nodeselinit))    /**< initialize node selector */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNodeselInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nodesel != NULL);

   SCIPnodeselSetInit(nodesel, nodeselinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of node selector */
SCIP_RETCODE SCIPsetNodeselExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELEXIT ((*nodeselexit))    /**< deinitialize node selector */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNodeselExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nodesel != NULL);

   SCIPnodeselSetExit(nodesel, nodeselexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of node selector */
SCIP_RETCODE SCIPsetNodeselInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELINITSOL ((*nodeselinitsol))/**< solving process initialization method of node selector */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNodeselInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nodesel != NULL);

   SCIPnodeselSetInitsol(nodesel, nodeselinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of node selector */
SCIP_RETCODE SCIPsetNodeselExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELEXITSOL ((*nodeselexitsol))/**< solving process deinitialization method of node selector */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNodeselExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(nodesel != NULL);

   SCIPnodeselSetExitsol(nodesel, nodeselexitsol);

   return SCIP_OKAY;
}

/** returns the node selector of the given name, or NULL if not existing */
SCIP_NODESEL* SCIPfindNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of node selector */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindNodesel(scip->set, name);
}

/** returns the array of currently available node selectors */
SCIP_NODESEL** SCIPgetNodesels(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nodesels;
}

/** returns the number of currently available node selectors */
int SCIPgetNNodesels(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nnodesels;
}

/** sets the priority of a node selector in standard mode */
SCIP_RETCODE SCIPsetNodeselStdPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   int                   priority            /**< new standard priority of the node selector */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPnodeselSetStdPriority(nodesel, scip->set, priority);

   return SCIP_OKAY;
}

/** sets the priority of a node selector in memory saving mode */
SCIP_RETCODE SCIPsetNodeselMemsavePriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   int                   priority            /**< new memory saving priority of the node selector */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPnodeselSetMemsavePriority(nodesel, scip->set, priority);

   return SCIP_OKAY;
}

/** returns the currently used node selector */
SCIP_NODESEL* SCIPgetNodesel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetGetNodesel(scip->set, scip->stat);
}
