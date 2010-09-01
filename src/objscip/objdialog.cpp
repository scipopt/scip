/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objdialog.cpp,v 1.10 2010/09/01 14:14:11 bzfheinz Exp $"

/**@file   objdialog.cpp
 * @brief  C++ wrapper for dialogs
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objdialog.h"




/*
 * Data structures
 */

/** dialog data */
struct SCIP_DialogData
{
   scip::ObjDialog*      objdialog;          /**< dialog object */
   SCIP_Bool             deleteobject;       /**< should the dialog object be deleted when dialog is freed? */
};




/*
 * Callback methods of dialog
 */

extern "C"
{

/** copy method for dialog plugins (called when SCIP copies plugins) */
static
SCIP_DECL_DIALOGCOPY(dialogCopyObj)
{  /*lint --e{715}*/
   SCIP_DIALOGDATA* dialogdata;
   
   assert(scip != NULL);
   
   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != NULL);
   assert(dialogdata->objdialog != NULL);
   assert(dialogdata->objdialog->scip_ != scip);

   if( dialogdata->objdialog->iscloneable() )
   {
      scip::ObjDialog*  newobjdialog;
      newobjdialog = (scip::ObjDialog*) dialogdata->objdialog->clone(scip);

      /* call include method of dialog object */
      SCIP_CALL( SCIPincludeObjDialog(scip, newobjdialog, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of dialog to free user data (called when SCIP is exiting) */
static
SCIP_DECL_DIALOGFREE(dialogFreeObj)
{  /*lint --e{715}*/
   SCIP_DIALOGDATA* dialogdata;

   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != NULL);
   assert(dialogdata->objdialog != NULL);
   assert(dialogdata->objdialog->scip_ == scip);

   /* call virtual method of dialog object */
   SCIP_CALL( dialogdata->objdialog->scip_free(scip, dialog) );

   /* free dialog object */
   if( dialogdata->deleteobject )
      delete dialogdata->objdialog;

   /* free dialog data */
   delete dialogdata;
   SCIPdialogSetData(dialog, NULL); /*lint !e64*/
   
   return SCIP_OKAY;
}


/** description output method of dialog */
static
SCIP_DECL_DIALOGDESC(dialogDescObj)
{  /*lint --e{715}*/
   SCIP_DIALOGDATA* dialogdata;

   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != NULL);
   assert(dialogdata->objdialog != NULL);
   assert(dialogdata->objdialog->scip_ == scip);

   /* call virtual method of dialog object */
   SCIP_CALL( dialogdata->objdialog->scip_desc(scip, dialog) );

   return SCIP_OKAY;
}

/** execution method of dialog */
static
SCIP_DECL_DIALOGEXEC(dialogExecObj)
{  /*lint --e{715}*/
   SCIP_DIALOGDATA* dialogdata;

   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != NULL);
   assert(dialogdata->objdialog != NULL);

   /* call virtual method of dialog object */
   SCIP_CALL( dialogdata->objdialog->scip_exec(scip, dialoghdlr, dialog, nextdialog) );

   return SCIP_OKAY;
}
}



/*
 * dialog specific interface methods
 */

/** creates the dialog for the given dialog object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjDialog*      objdialog,          /**< dialog object */
   SCIP_Bool             deleteobject        /**< should the dialog object be deleted when dialog is freed? */
   )
{/*lint --e{429} */
   SCIP_DIALOG* parentdialog;

   assert(scip != NULL);
   assert(objdialog != NULL);

   /* get parent dialog */
   parentdialog = SCIPgetRootDialog(scip);
   assert(parentdialog != NULL);
   /* TODO: (optional) change parent dialog from root dialog to another existing dialog (needs to be a menu) */
   
   /* create, include, and release dialog */
   if( !SCIPdialogHasEntry(parentdialog, objdialog->scip_name_) )
   {
      SCIP_DIALOGDATA* dialogdata;
      SCIP_DIALOG* dialog;

      /* create dialog data */
      dialogdata = new SCIP_DIALOGDATA;
      dialogdata->objdialog = objdialog;
      dialogdata->deleteobject = deleteobject;

      SCIP_CALL( SCIPincludeDialog(scip, &dialog, 
            dialogCopyObj,
            dialogExecObj, dialogDescObj, dialogFreeObj,
            objdialog->scip_name_, objdialog->scip_desc_, objdialog->scip_issubmenu_, dialogdata) );
      SCIP_CALL( SCIPaddDialogEntry(scip, parentdialog, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   return SCIP_OKAY;
}

/** returns the dialog object of the given name, or 0 if not existing */
scip::ObjDialog* SCIPfindObjDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of dialog */
   )
{
   SCIP_DIALOG* dialog;
   SCIP_DIALOGDATA* dialogdata;
   
   dialog = SCIPfindDialog(scip, name);
   if( dialog == NULL )
      return 0;

   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != NULL);
   
   return dialogdata->objdialog;
}
