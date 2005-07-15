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
#pragma ident "@(#) $Id: objmessagehdlr.cpp,v 1.1 2005/07/15 17:20:02 bzfpfend Exp $"

/**@file   objmessagehdlr.cpp
 * @brief  C++ wrapper for message handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objmessagehdlr.h"




/*
 * Data structures
 */

/** message handler data */
struct MessagehdlrData
{
   scip::ObjMessagehdlr* objmessagehdlr;/**< message handler object */
   Bool             deleteobject;       /**< should the message handler object be deleted when message handler is freed? */
};




/*
 * Callback methods of file reader
 */

/** error message print method of message handler */
static
DECL_MESSAGEERROR(messagehdlrErrorObj)
{  /*lint --e{715}*/
   MESSAGEHDLRDATA* messagehdlrdata;

   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);
   assert(messagehdlrdata->objmessagehdlr != NULL);

   /* call virtual method of messagehdlr object */
   messagehdlrdata->objmessagehdlr->scip_error(messagehdlr, file, msg);
}


/** warning message print method of message handler */
static
DECL_MESSAGEWARNING(messagehdlrWarningObj)
{  /*lint --e{715}*/
   MESSAGEHDLRDATA* messagehdlrdata;

   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);
   assert(messagehdlrdata->objmessagehdlr != NULL);

   /* call virtual method of messagehdlr object */
   messagehdlrdata->objmessagehdlr->scip_warning(messagehdlr, file, msg);
}


/** dialog message print method of message handler */
static
DECL_MESSAGEDIALOG(messagehdlrDialogObj)
{  /*lint --e{715}*/
   MESSAGEHDLRDATA* messagehdlrdata;

   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);
   assert(messagehdlrdata->objmessagehdlr != NULL);

   /* call virtual method of messagehdlr object */
   messagehdlrdata->objmessagehdlr->scip_dialog(messagehdlr, file, msg);
}


/** info message print method of message handler */
static
DECL_MESSAGEINFO(messagehdlrInfoObj)
{  /*lint --e{715}*/
   MESSAGEHDLRDATA* messagehdlrdata;

   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);
   assert(messagehdlrdata->objmessagehdlr != NULL);

   /* call virtual method of messagehdlr object */
   messagehdlrdata->objmessagehdlr->scip_info(messagehdlr, file, msg);
}




/*
 * message handler specific interface methods
 */

/** creates the message handler for the given message handler object */
RETCODE SCIPcreateObjMessagehdlr(
   MESSAGEHDLR**    messagehdlr,        /**< pointer to store the message handler */
   scip::ObjMessagehdlr* objmessagehdlr,/**< message handler object */
   Bool             deleteobject        /**< should the message handler object be deleted when message handler is freed? */
   )
{
   MESSAGEHDLRDATA* messagehdlrdata;

   /* create file messagehdlr data */
   messagehdlrdata = new MESSAGEHDLRDATA;
   messagehdlrdata->objmessagehdlr = objmessagehdlr;
   messagehdlrdata->deleteobject = deleteobject;

   /* create message handler */
   CHECK_OKAY( SCIPcreateMessagehdlr(messagehdlr, objmessagehdlr->scip_bufferedoutput_,
                  messagehdlrErrorObj, messagehdlrWarningObj, messagehdlrDialogObj, messagehdlrInfoObj,
                  messagehdlrdata) );

   return SCIP_OKAY;
}

/** destroys the message handler that was created by SCIPcreateObjMessagehdlr();
 *  if deleteobject was set to TRUE in SCIPcreateObjMessagehdlr(), the message handler object is deleted
 */
RETCODE SCIPfreeObjMessagehdlr(
   MESSAGEHDLR**    messagehdlr         /**< pointer to the message handler */
   )
{
   MESSAGEHDLRDATA* messagehdlrdata;

   assert(messagehdlr != NULL);

   messagehdlrdata = SCIPmessagehdlrGetData(*messagehdlr);
   assert(messagehdlrdata != NULL);
   assert(messagehdlrdata->objmessagehdlr != NULL);

   /* free message handler object */
   if( messagehdlrdata->deleteobject )
      delete messagehdlrdata->objmessagehdlr;

   /* free message handler data */
   delete messagehdlrdata;
   SCIPmessagehdlrSetData(*messagehdlr, NULL);

   /* free message handler */
   CHECK_OKAY( SCIPfreeMessagehdlr(messagehdlr) );
   
   return SCIP_OKAY;
}
   
/** returns the message handler object for the given message handler */
scip::ObjMessagehdlr* SCIPgetObjMessagehdlr(
   MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   MESSAGEHDLRDATA* messagehdlrdata;

   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);

   return messagehdlrdata->objmessagehdlr;
}

