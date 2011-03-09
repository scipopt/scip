/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   EventhdlrNewSol.h
 * @brief  event handler for new solutions in TSP
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __EVENTHDLRNEWSOL_H__
#define __EVENTHDLRNEWSOL_H__


#include "objscip/objscip.h"

namespace tsp
{
/** C++ wrapper object for event handlers */
class EventhdlrNewSol : public scip::ObjEventhdlr
{
public:     
   /** default constructor */
   EventhdlrNewSol(
      SCIP* scip
      )
      : ObjEventhdlr(scip, "newsol","event handler for new solutions in TSP")
   {
   }

   /** destructor */
   virtual ~EventhdlrNewSol()
   {
   }

   /** destructor of event handler to free user data (called when SCIP is exiting) */
   virtual SCIP_RETCODE scip_free(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr           /**< the event handler itself */
      );
   
   /** initialization method of event handler (called after problem was transformed) */
   virtual SCIP_RETCODE scip_init(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr           /**< the event handler itself */
      );
   
   /** deinitialization method of event handler (called before transformed problem is freed) */
   virtual SCIP_RETCODE scip_exit(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr           /**< the event handler itself */
      );
   
   /** solving process initialization method of event handler (called when branch and bound process is about to begin)
    *
    *  This method is called when the presolving was finished and the branch and bound process is about to begin.
    *  The event handler may use this call to initialize its branch and bound specific data.
    *
    */
   virtual SCIP_RETCODE scip_initsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr           /**< the event handler itself */
      );
   
   /** solving process deinitialization method of event handler (called before branch and bound process data is freed)
    *
    *  This method is called before the branch and bound process is freed.
    *  The event handler should use this call to clean up its branch and bound data.
    */
   virtual SCIP_RETCODE scip_exitsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr           /**< the event handler itself */
      );
   
   /** frees specific constraint data */
   virtual SCIP_RETCODE scip_delete(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr,          /**< the event handler itself */
      SCIP_EVENTDATA**   eventdata           /**< pointer to the event data to free */
      );

   /** execution method of event handler
    *
    *  Processes the event. The method is called every time an event occurs, for which the event handler
    *  is responsible. Event handlers may declare themselves resposible for events by calling the
    *  corresponding SCIPcatch...() method. This method creates an event filter object to point to the
    *  given event handler and event data.
    */
   virtual SCIP_RETCODE scip_exec(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr,          /**< the event handler itself */
      SCIP_EVENT*        event,              /**< event to process */
      SCIP_EVENTDATA*    eventdata           /**< user data for the event */
      );
};

} /* namespace tsp */


#endif
