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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objeventhdlr.h,v 1.2 2005/05/31 17:20:08 bzfpfend Exp $"

/**@file   objeventhdlr.h
 * @brief  C++ wrapper for event handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJEVENTHDLR_H__
#define __OBJEVENTHDLR_H__


extern "C" 
{
#include "scip/scip.h"
}


namespace scip
{

/** C++ wrapper object for event handlers */
class ObjEventhdlr
{
public:
   /** name of the event handler */
   const char* const scip_name_;
   
   /** description of the event handler */
   const char* const scip_desc_;
   
   /** default constructor */
   ObjEventhdlr(
      const char*   name,               /**< name of event handler */
      const char*   desc                /**< description of event handler */
      )
      : scip_name_(name),
        scip_desc_(desc)
   {
   }

   /** destructor */
   virtual ~ObjEventhdlr()
   {
   }

   /** destructor of event handler to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      EVENTHDLR*    eventhdlr           /**< the event handler itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** initialization method of event handler (called after problem was transformed) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      EVENTHDLR*    eventhdlr           /**< the event handler itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of event handler (called before transformed problem is freed) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      EVENTHDLR*    eventhdlr           /**< the event handler itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** solving process initialization method of event handler (called when branch and bound process is about to begin)
    *
    *  This method is called when the presolving was finished and the branch and bound process is about to begin.
    *  The event handler may use this call to initialize its branch and bound specific data.
    *
    */
   virtual RETCODE scip_initsol(
      SCIP*         scip,               /**< SCIP data structure */
      EVENTHDLR*    eventhdlr           /**< the event handler itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** solving process deinitialization method of event handler (called before branch and bound process data is freed)
    *
    *  This method is called before the branch and bound process is freed.
    *  The event handler should use this call to clean up its branch and bound data.
    */
   virtual RETCODE scip_exitsol(
      SCIP*         scip,               /**< SCIP data structure */
      EVENTHDLR*    eventhdlr           /**< the event handler itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** frees specific constraint data */
   virtual RETCODE scip_delete(
      SCIP*         scip,               /**< SCIP data structure */
      EVENTHDLR*    eventhdlr,          /**< the event handler itself */
      EVENTDATA**   eventdata           /**< pointer to the event data to free */
      )
   {
      return SCIP_OKAY;
   }

   /** execution method of event handler
    *
    *  Processes the event. The method is called every time an event occurs, for which the event handler
    *  is responsible. Event handlers may declare themselves resposible for events by calling the
    *  corresponding SCIPcatch...() method. This method creates an event filter object to point to the
    *  given event handler and event data.
    */
   virtual RETCODE scip_exec(
      SCIP*         scip,               /**< SCIP data structure */
      EVENTHDLR*    eventhdlr,          /**< the event handler itself */
      EVENT*        event,              /**< event to process */
      EVENTDATA*    eventdata           /**< user data for the event */
      ) = 0;
};

} /* namespace scip */


   
/** creates the event handler for the given event handler object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       MyEventhdlr* myeventhdlr = new MyEventhdlr(...);
 *       CHECK_OKAY( SCIPincludeObjEventhdlr(scip, &myeventhdlr, FALSE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );
 *       delete myeventhdlr;    // delete eventhdlr AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       CHECK_OKAY( SCIPincludeObjEventhdlr(scip, new MyEventhdlr(...), TRUE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );  // destructor of MyEventhdlr is called here
 */
extern
RETCODE SCIPincludeObjEventhdlr(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjEventhdlr* objeventhdlr,    /**< event handler object */
   Bool             deleteobject        /**< should the event handler object be deleted when eventhdlristic is freed? */
   );

/** returns the eventhdlr object of the given name, or NULL if not existing */
extern
scip::ObjEventhdlr* SCIPfindObjEventhdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of event handler */
   );

/** returns the eventhdlr object for the given event handler */
extern
scip::ObjEventhdlr* SCIPgetObjEventhdlr(
   SCIP*            scip,               /**< SCIP data structure */
   EVENTHDLR*       eventhdlr           /**< event handler */
   );

#endif
