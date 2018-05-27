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

/**@file   scip_prop.c
 * @brief  public methods for propagator plugins
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

#include "scip/scip_prop.h"

#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_prop.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** creates a propagator and includes it in SCIP.
 *
 *  @note method has all propagator callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludePropBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeProp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of propagator */
   const char*           desc,               /**< description of propagator */
   int                   priority,           /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling propagator */
   SCIP_Bool             delay,              /**< should propagator be delayed, if other propagators found reductions? */
   SCIP_PROPTIMING       timingmask,         /**< positions in the node solving loop where propagator should be executed */
   int                   presolpriority,     /**< presolving priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   presolmaxrounds,    /**< maximal number of presolving rounds the propagator participates in (-1: no limit) */
   SCIP_PRESOLTIMING     presoltiming,       /**< timing mask of the propagator's presolving method */
   SCIP_DECL_PROPCOPY    ((*propcopy)),      /**< copy method of propagator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PROPFREE    ((*propfree)),      /**< destructor of propagator */
   SCIP_DECL_PROPINIT    ((*propinit)),      /**< initialize propagator */
   SCIP_DECL_PROPEXIT    ((*propexit)),      /**< deinitialize propagator */
   SCIP_DECL_PROPINITPRE ((*propinitpre)),   /**< presolving initialization method of propagator */
   SCIP_DECL_PROPEXITPRE ((*propexitpre)),   /**< presolving deinitialization method of propagator */
   SCIP_DECL_PROPINITSOL ((*propinitsol)),   /**< solving process initialization method of propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol)),   /**< solving process deinitialization method of propagator */
   SCIP_DECL_PROPPRESOL  ((*proppresol)),    /**< presolving method */
   SCIP_DECL_PROPEXEC    ((*propexec)),      /**< execution method of propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop)),   /**< propagation conflict resolving method */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_PROP* prop;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeProp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether propagator is already present */
   if( SCIPfindProp(scip, name) != NULL )
   {
      SCIPerrorMessage("propagator <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPpropCreate(&prop, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, priority, freq, delay, timingmask, presolpriority, presolmaxrounds, presoltiming,
         propcopy, propfree, propinit, propexit, propinitpre, propexitpre, propinitsol, propexitsol,
         proppresol, propexec, propresprop, propdata) );
   SCIP_CALL( SCIPsetIncludeProp(scip->set, prop) );

   return SCIP_OKAY;
}

/** creates a propagator and includes it in SCIP. All non-fundamental (or optional) callbacks will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetPropInit(), SCIPsetPropExit(),
 *  SCIPsetPropCopy(), SCIPsetPropFree(), SCIPsetPropInitsol(), SCIPsetPropExitsol(),
 *  SCIPsetPropInitpre(), SCIPsetPropExitpre(), SCIPsetPropPresol(), and SCIPsetPropResprop().
 *
*  @note if you want to set all callbacks with a single method call, consider using SCIPincludeProp() instead
 */
SCIP_RETCODE SCIPincludePropBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP**           propptr,            /**< reference to a propagator pointer, or NULL */
   const char*           name,               /**< name of propagator */
   const char*           desc,               /**< description of propagator */
   int                   priority,           /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling propagator */
   SCIP_Bool             delay,              /**< should propagator be delayed, if other propagators found reductions? */
   SCIP_PROPTIMING       timingmask,         /**< positions in the node solving loop where propagators should be executed */
   SCIP_DECL_PROPEXEC    ((*propexec)),      /**< execution method of propagator */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_PROP* prop;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludePropBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether propagator is already present */
   if( SCIPfindProp(scip, name) != NULL )
   {
      SCIPerrorMessage("propagator <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPpropCreate(&prop, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, priority, freq, delay, timingmask, 0, -1, SCIP_PRESOLTIMING_ALWAYS,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
         NULL, propexec, NULL, propdata) );
   SCIP_CALL( SCIPsetIncludeProp(scip->set, prop) );

   if( propptr != NULL )
      *propptr = prop;

   return SCIP_OKAY;
}

/** sets copy method of propagator */
SCIP_RETCODE SCIPsetPropCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPCOPY    ((*propcopy))       /**< copy method of propagator or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetCopy(prop, propcopy);

   return SCIP_OKAY;
}

/** sets destructor method of propagator */
SCIP_RETCODE SCIPsetPropFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPFREE    ((*propfree))       /**< destructor of propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetFree(prop, propfree);

   return SCIP_OKAY;
}

/** sets initialization method of propagator */
SCIP_RETCODE SCIPsetPropInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPINIT    ((*propinit))       /**< initialize propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetInit(prop, propinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of propagator */
SCIP_RETCODE SCIPsetPropExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPEXIT    ((*propexit))       /**< deinitialize propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetExit(prop, propexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of propagator */
SCIP_RETCODE SCIPsetPropInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPINITSOL((*propinitsol))     /**< solving process initialization method of propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetInitsol(prop, propinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of propagator */
SCIP_RETCODE SCIPsetPropExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol))    /**< solving process deinitialization method of propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetExitsol(prop, propexitsol);

   return SCIP_OKAY;
}

/** sets preprocessing initialization method of propagator */
SCIP_RETCODE SCIPsetPropInitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPINITPRE((*propinitpre))     /**< preprocessing initialization method of propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropInitpre", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetInitpre(prop, propinitpre);

   return SCIP_OKAY;
}

/** sets preprocessing deinitialization method of propagator */
SCIP_RETCODE SCIPsetPropExitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPEXITPRE((*propexitpre))     /**< preprocessing deinitialization method of propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropExitpre", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetExitpre(prop, propexitpre);

   return SCIP_OKAY;
}

/** sets presolving method of propagator */
SCIP_RETCODE SCIPsetPropPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPPRESOL((*proppresol)),      /**< presolving method of propagator */
   int                   presolpriority,     /**< presolving priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   presolmaxrounds,    /**< maximal number of presolving rounds the propagator participates in (-1: no limit) */
   SCIP_PRESOLTIMING     presoltiming        /**< timing mask of the propagator's presolving method */
   )
{
   const char* name;
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropPresol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);
   SCIP_CALL( SCIPpropSetPresol(prop, proppresol, presolpriority, presolmaxrounds, presoltiming) );

   name = SCIPpropGetName(prop);

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/maxprerounds", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, presolmaxrounds) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/presolpriority", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, presolpriority) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/presoltiming", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, (int) presoltiming) );

   return SCIP_OKAY;
}

/** sets propagation conflict resolving callback of propagator */
SCIP_RETCODE SCIPsetPropResprop(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop))    /**< propagation conflict resolving callback */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropResprop", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetResprop(prop, propresprop);

   return SCIP_OKAY;
}


/** returns the propagator of the given name, or NULL if not existing */
SCIP_PROP* SCIPfindProp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of propagator */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindProp(scip->set, name);
}

/** returns the array of currently available propagators */
SCIP_PROP** SCIPgetProps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortProps(scip->set);

   return scip->set->props;
}

/** returns the number of currently available propagators */
int SCIPgetNProps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nprops;
}

/** sets the priority of a propagator */
SCIP_RETCODE SCIPsetPropPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   int                   priority            /**< new priority of the propagator */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPpropSetPriority(prop, scip->set, priority);

   return SCIP_OKAY;
}

/** sets the presolving priority of a propagator */
SCIP_RETCODE SCIPsetPropPresolPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   int                   presolpriority      /**< new presol priority of the propagator */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPpropSetPresolPriority(prop, scip->set, presolpriority);

   return SCIP_OKAY;
}
