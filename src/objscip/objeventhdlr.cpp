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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objeventhdlr.cpp,v 1.1 2005/02/24 16:26:36 bzfpfend Exp $"

/**@file   objeventhdlr.cpp
 * @brief  C++ wrapper for event handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objeventhdlr.h"




/*
 * Data structures
 */

/** event handler data */
struct EventhdlrData
{
   scip::ObjEventhdlr* objeventhdlr;    /**< event handler object */
   Bool             deleteobject;       /**< should the event handler object be deleted when eventhdlristic is freed? */
};




/*
 * Callback methods of event handler
 */

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
DECL_EVENTFREE(eventhdlrFreeObj)
{  /*lint --e{715}*/
   EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);

   /* call virtual method of eventhdlr object */
   CHECK_OKAY( eventhdlrdata->objeventhdlr->scip_free(scip, eventhdlr) );

   /* free eventhdlr object */
   if( eventhdlrdata->deleteobject )
      delete eventhdlrdata->objeventhdlr;

   /* free eventhdlr data */
   delete eventhdlrdata;
   SCIPeventhdlrSetData(eventhdlr, NULL);
   
   return SCIP_OKAY;
}


/** initialization method of event handler (called after problem was transformed) */
static
DECL_EVENTINIT(eventhdlrInitObj)
{  /*lint --e{715}*/
   EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);

   /* call virtual method of eventhdlr object */
   CHECK_OKAY( eventhdlrdata->objeventhdlr->scip_init(scip, eventhdlr) );

   return SCIP_OKAY;
}


/** deinitialization method of event handler (called before transformed problem is freed) */
static
DECL_EVENTEXIT(eventhdlrExitObj)
{  /*lint --e{715}*/
   EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);

   /* call virtual method of eventhdlr object */
   CHECK_OKAY( eventhdlrdata->objeventhdlr->scip_exit(scip, eventhdlr) );

   return SCIP_OKAY;
}


/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
DECL_EVENTINITSOL(eventhdlrInitsolObj)
{  /*lint --e{715}*/
   EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);

   /* call virtual method of eventhdlr object */
   CHECK_OKAY( eventhdlrdata->objeventhdlr->scip_initsol(scip, eventhdlr) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
DECL_EVENTEXITSOL(eventhdlrExitsolObj)
{  /*lint --e{715}*/
   EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);

   /* call virtual method of eventhdlr object */
   CHECK_OKAY( eventhdlrdata->objeventhdlr->scip_exitsol(scip, eventhdlr) );

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
DECL_EVENTDELETE(eventhdlrDeleteObj)
{  /*lint --e{715}*/
   EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);

   /* call virtual method of eventhdlr object */
   CHECK_OKAY( eventhdlrdata->objeventhdlr->scip_delete(scip, eventhdlr, eventdata) );

   return SCIP_OKAY;
}


/** execution method of event handler */
static
DECL_EVENTEXEC(eventhdlrExecObj)
{  /*lint --e{715}*/
   EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);

   /* call virtual method of eventhdlr object */
   CHECK_OKAY( eventhdlrdata->objeventhdlr->scip_exec(scip, eventhdlr, event, eventdata) );

   return SCIP_OKAY;
}




/*
 * event handler specific interface methods
 */

/** creates the event handler for the given event handler object and includes it in SCIP */
RETCODE SCIPincludeObjEventhdlr(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjEventhdlr*   objeventhdlr,            /**< event handler object */
   Bool             deleteobject        /**< should the event handler object be deleted when eventhdlristic is freed? */
   )
{
   EVENTHDLRDATA* eventhdlrdata;

   /* create event handler data */
   eventhdlrdata = new EVENTHDLRDATA;
   eventhdlrdata->objeventhdlr = objeventhdlr;
   eventhdlrdata->deleteobject = deleteobject;

   /* include event handler */
   CHECK_OKAY( SCIPincludeEventhdlr(scip, objeventhdlr->scip_name_, objeventhdlr->scip_desc_,
         eventhdlrFreeObj, eventhdlrInitObj, eventhdlrExitObj, 
         eventhdlrInitsolObj, eventhdlrExitsolObj, eventhdlrDeleteObj, eventhdlrExecObj,
         eventhdlrdata) );

   return SCIP_OKAY;
}

/** returns the eventhdlr object of the given name, or NULL if not existing */
scip::ObjEventhdlr* SCIPfindObjEventhdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of event handler */
   )
{
   EVENTHDLR* eventhdlr;
   EVENTHDLRDATA* eventhdlrdata;

   eventhdlr = SCIPfindEventhdlr(scip, name);
   if( eventhdlr == NULL )
      return NULL;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   return eventhdlrdata->objeventhdlr;
}
   
/** returns the eventhdlr object for the given event handler */
scip::ObjEventhdlr* SCIPgetObjEventhdlr(
   SCIP*            scip,               /**< SCIP data structure */
   EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   return eventhdlrdata->objeventhdlr;
}
