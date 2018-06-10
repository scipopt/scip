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

/**@file   scip_event.c
 * @brief  public methods for event handler plugins and event handlers
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

#include "scip/scip_event.h"

#include "scip/pub_message.h"
#include "scip/pub_var.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** creates an event handler and includes it in SCIP
 *
 *  @note method has all event handler callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludeEventhdlrBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeEventhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of event handler */
   const char*           desc,               /**< description of event handler */
   SCIP_DECL_EVENTCOPY   ((*eventcopy)),     /**< copy method of event handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_EVENTFREE   ((*eventfree)),     /**< destructor of event handler */
   SCIP_DECL_EVENTINIT   ((*eventinit)),     /**< initialize event handler */
   SCIP_DECL_EVENTEXIT   ((*eventexit)),     /**< deinitialize event handler */
   SCIP_DECL_EVENTINITSOL((*eventinitsol)),  /**< solving process initialization method of event handler */
   SCIP_DECL_EVENTEXITSOL((*eventexitsol)),  /**< solving process deinitialization method of event handler */
   SCIP_DECL_EVENTDELETE ((*eventdelete)),   /**< free specific event data */
   SCIP_DECL_EVENTEXEC   ((*eventexec)),     /**< execute event handler */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeEventhdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether event handler is already present */
   if( SCIPfindEventhdlr(scip, name) != NULL )
   {
      SCIPerrorMessage("event handler <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPeventhdlrCreate(&eventhdlr, scip->set, name, desc,
         eventcopy, eventfree, eventinit, eventexit, eventinitsol, eventexitsol, eventdelete, eventexec,
         eventhdlrdata) );
   SCIP_CALL( SCIPsetIncludeEventhdlr(scip->set, eventhdlr) );

   return SCIP_OKAY;
}

/** creates an event handler and includes it in SCIP with all its non-fundamental callbacks set
 *  to NULL; if needed, non-fundamental callbacks can be set afterwards via setter functions
 *  SCIPsetEventhdlrCopy(), SCIPsetEventhdlrFree(), SCIPsetEventhdlrInit(), SCIPsetEventhdlrExit(),
 *  SCIPsetEventhdlrInitsol(), SCIPsetEventhdlrExitsol(), and SCIPsetEventhdlrDelete()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeEventhdlr() instead
 */
SCIP_RETCODE SCIPincludeEventhdlrBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR**      eventhdlrptr,       /**< reference to an event handler, or NULL */
   const char*           name,               /**< name of event handler */
   const char*           desc,               /**< description of event handler */
   SCIP_DECL_EVENTEXEC   ((*eventexec)),     /**< execute event handler */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeEventhdlrBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether event handler is already present */
   if( SCIPfindEventhdlr(scip, name) != NULL )
   {
      SCIPerrorMessage("event handler <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPeventhdlrCreate(&eventhdlr, scip->set, name, desc,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, eventexec,
         eventhdlrdata) );
   SCIP_CALL( SCIPsetIncludeEventhdlr(scip->set, eventhdlr) );

   if( eventhdlrptr != NULL )
      *eventhdlrptr = eventhdlr;

   return SCIP_OKAY;
}

/** sets copy callback of the event handler */
SCIP_RETCODE SCIPsetEventhdlrCopy(
   SCIP*                 scip,               /**< scip instance */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTCOPY   ((*eventcopy))      /**< copy callback of the event handler */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetEventhdlrCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPeventhdlrSetCopy(eventhdlr, eventcopy);
   return SCIP_OKAY;
}

/** sets deinitialization callback of the event handler */
SCIP_RETCODE SCIPsetEventhdlrFree(
   SCIP*                 scip,               /**< scip instance */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTFREE   ((*eventfree))      /**< deinitialization callback of the event handler */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetEventhdlrFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPeventhdlrSetFree(eventhdlr, eventfree);
   return SCIP_OKAY;
}

/** sets initialization callback of the event handler */
SCIP_RETCODE SCIPsetEventhdlrInit(
   SCIP*                 scip,               /**< scip instance */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTINIT   ((*eventinit))      /**< initialize event handler */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetEventhdlrInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPeventhdlrSetInit(eventhdlr, eventinit);
   return SCIP_OKAY;
}

/** sets deinitialization callback of the event handler */
SCIP_RETCODE SCIPsetEventhdlrExit(
   SCIP*                 scip,               /**< scip instance */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTEXIT   ((*eventexit))      /**< deinitialize event handler */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetEventhdlrExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPeventhdlrSetExit(eventhdlr, eventexit);
   return SCIP_OKAY;
}

/** sets solving process initialization callback of the event handler */
SCIP_RETCODE SCIPsetEventhdlrInitsol(
   SCIP*                 scip,               /**< scip instance */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTINITSOL((*eventinitsol))   /**< solving process initialization callback of event handler */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetEventhdlrInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPeventhdlrSetInitsol(eventhdlr, eventinitsol);
   return SCIP_OKAY;
}

/** sets solving process deinitialization callback of the event handler */
SCIP_RETCODE SCIPsetEventhdlrExitsol(
   SCIP*                 scip,               /**< scip instance */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTEXITSOL((*eventexitsol))   /**< solving process deinitialization callback of event handler */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetEventhdlrExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPeventhdlrSetExitsol(eventhdlr, eventexitsol);
   return SCIP_OKAY;
}

/** sets callback of the event handler to free specific event data */
SCIP_RETCODE SCIPsetEventhdlrDelete(
   SCIP*                 scip,               /**< scip instance */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTDELETE ((*eventdelete))    /**< free specific event data */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetEventhdlrDelete", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPeventhdlrSetDelete(eventhdlr, eventdelete);
   return SCIP_OKAY;
}

/** returns the event handler of the given name, or NULL if not existing */
SCIP_EVENTHDLR* SCIPfindEventhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of event handler */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindEventhdlr(scip->set, name);
}

/** returns the array of currently available event handlers */
SCIP_EVENTHDLR** SCIPgetEventhdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->eventhdlrs;
}

/** returns the number of currently available event handlers */
int SCIPgetNEventhdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->neventhdlrs;
}

/** catches a global (not variable or row dependent) event
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPcatchEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcatchEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPeventfilterAdd(scip->eventfilter, scip->mem->probmem, scip->set,
         eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** drops a global event (stops to track event)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPdropEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int                   filterpos           /**< position of event filter entry returned by SCIPcatchEvent(), or -1 */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPdropEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPeventfilterDel(scip->eventfilter, scip->mem->probmem, scip->set,
         eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** catches an objective value or domain change event on the given transformed variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPcatchVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< transformed variable to catch event for */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcatchVarEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( (eventtype & SCIP_EVENTTYPE_VARCHANGED) == 0 )
   {
      SCIPerrorMessage("event does not operate on a single variable\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPvarIsOriginal(var) )
   {
      SCIPerrorMessage("cannot catch events on original variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPvarCatchEvent(var, scip->mem->probmem, scip->set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** drops an objective value or domain change event (stops to track event) on the given transformed variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPdropVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< transformed variable to drop event for */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int                   filterpos           /**< position of event filter entry returned by SCIPcatchVarEvent(), or -1 */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPdropVarEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( SCIPvarIsOriginal(var) )
   {
      SCIPerrorMessage("cannot drop events on original variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPvarDropEvent(var, scip->mem->probmem, scip->set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** catches a row coefficient, constant, or side change event on the given row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPcatchRowEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< linear row to catch event for */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcatchRowEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( (eventtype & SCIP_EVENTTYPE_ROWCHANGED) == 0 )
   {
      SCIPerrorMessage("event does not operate on a single row\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIProwCatchEvent(row, scip->mem->probmem, scip->set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** drops a row coefficient, constant, or side change event (stops to track event) on the given row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPdropRowEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< linear row to drop event for */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int                   filterpos           /**< position of event filter entry returned by SCIPcatchVarEvent(), or -1 */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPdropRowEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIProwDropEvent(row, scip->mem->probmem, scip->set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}
