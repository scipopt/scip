/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   JniScipDialog.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP dialog callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipDialog.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/pub_dialog.h"

#include <string.h>

/** returns the root dialog of the dialog handler */
JNIEXPORT
jlong JNISCIPDIALOG(dialoghdlrGetRoot)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialoghdlr         /**< dialog handler */
   )
{
   SCIP_DIALOGHDLR* dialoghdlr;

   /* convert JNI pointer into C pointer */
   dialoghdlr = (SCIP_DIALOGHDLR*) (size_t) jdialoghdlr;
   assert(dialoghdlr != NULL);

   return (jlong) (size_t) SCIPdialoghdlrGetRoot(dialoghdlr);
}

/** clears the input command buffer of the dialog handler */
JNIEXPORT
void JNISCIPDIALOG(dialoghdlrClearBuffer)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialoghdlr         /**< dialog handler */
   )
{
   SCIP_DIALOGHDLR* dialoghdlr;

   /* convert JNI pointer into C pointer */
   dialoghdlr = (SCIP_DIALOGHDLR*) (size_t) jdialoghdlr;
   assert(dialoghdlr != NULL);

   SCIPdialoghdlrClearBuffer(dialoghdlr);
}

/** returns TRUE iff input command buffer is empty */
JNIEXPORT
jboolean JNISCIPDIALOG(dialoghdlrIsBufferEmpty)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialoghdlr         /**< dialog handler */
   )
{
   SCIP_DIALOGHDLR* dialoghdlr;

   /* convert JNI pointer into C pointer */
   dialoghdlr = (SCIP_DIALOGHDLR*) (size_t) jdialoghdlr;
   assert(dialoghdlr != NULL);

   return (jboolean) SCIPdialoghdlrIsBufferEmpty(dialoghdlr);
}

/** adds a single line of input to the dialog handler which is treated as if the user entered the command line */
JNIEXPORT
void JNISCIPDIALOG(dialoghdlrAddInputLine)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialoghdlr,        /**< dialog handler */
   jstring               jinputline          /**< input line to add */
   )
{
   SCIP_DIALOGHDLR* dialoghdlr;
   const char* inputline;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   dialoghdlr = (SCIP_DIALOGHDLR*) (size_t) jdialoghdlr;
   assert(dialoghdlr != NULL);

   /* convert JNI string into C const char* */
   inputline = (*env)->GetStringUTFChars(env, jinputline, &iscopy);
   if( inputline == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPdialoghdlrAddInputLine(dialoghdlr, inputline) );

   (*env)->ReleaseStringUTFChars(env, jinputline, inputline);
}

/** adds a command to the command history of the dialog handler; if a dialog is given, the command is preceeded
 *  by the dialog's command path; if no command is given, only the path to the dialog is added to the command history
 */
JNIEXPORT
void JNISCIPDIALOG(dialoghdlrAddHistory)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialoghdlr,        /**< dialog handler */
   jlong                 jdialog,            /**< current dialog, or NULL */
   jstring               jcommand,           /**< command string to add to the command history, or NULL */
   jboolean              jescapecommand      /**< should special characters in command be prefixed by an escape char? */
   )
{
   SCIP_DIALOGHDLR* dialoghdlr;
   SCIP_DIALOG* dialog;
   const char* command;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   dialoghdlr = (SCIP_DIALOGHDLR*) (size_t) jdialoghdlr;
   dialog = (SCIP_DIALOG*) (size_t) jdialog;

   assert(dialoghdlr != NULL);
   assert(dialog != NULL);

   /* convert JNI string into C const char* */
   command = (*env)->GetStringUTFChars(env, jcommand, &iscopy);
   assert(iscopy);

   JNISCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, command, (SCIP_Bool)jescapecommand) );

   (*env)->ReleaseStringUTFChars(env, jcommand, command);
}

/** returns TRUE iff a dialog entry matching exactly the given name is existing in the given dialog */
JNIEXPORT
jboolean JNISCIPDIALOG(dialogHasEntry)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog,            /**< dialog handler */
   jstring               jentryname          /**< name of the dialog entry to find */
   )
{
   SCIP_DIALOG* dialog;
   const char* entryname;
   jboolean iscopy;
   jboolean ret;

   /* convert JNI pointer into C pointer */
   dialog = (SCIP_DIALOG*) (size_t) jdialog;
   assert(dialog != NULL);

   /* convert JNI string into C const char* */
   entryname = (*env)->GetStringUTFChars(env, jentryname, &iscopy);
   if( entryname == NULL )
      SCIPABORT();

   assert(iscopy);

   ret = SCIPdialogHasEntry(dialog, entryname);

   (*env)->ReleaseStringUTFChars(env, jentryname, entryname);

   return (jboolean) ret;
}

/** displays the dialog's menu */
JNIEXPORT
void JNISCIPDIALOG(dialogDisplayMenu)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog,            /**< dialog handler */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP_DIALOG* dialog;
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   dialog = (SCIP_DIALOG*) (size_t) jdialog;
   scip = (SCIP*) (size_t) jscip;

   assert(dialog != NULL);
   assert(scip != NULL);

   JNISCIP_CALL( SCIPdialogDisplayMenu(dialog, scip) );
}

/** displays the entry for the dialog in it's parent's menu */
JNIEXPORT
void JNISCIPDIALOG(dialogDisplayMenuEntry)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog,            /**< dialog handler */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP_DIALOG* dialog;
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   dialog = (SCIP_DIALOG*) (size_t) jdialog;
   scip = (SCIP*) (size_t) jscip;

   assert(dialog != NULL);
   assert(scip != NULL);

   JNISCIP_CALL( SCIPdialogDisplayMenuEntry(dialog, scip) );
}

/** displays all dialog entries with names starting with the given "entryname" */
JNIEXPORT
void JNISCIPDIALOG(dialogDisplayCompletions)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog,            /**< dialog handler */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jentryname          /**< name of the dialog entry to find */
   )
{
   SCIP_DIALOG* dialog;
   SCIP* scip;
   const char* entryname;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   dialog = (SCIP_DIALOG*) (size_t) jdialog;
   scip = (SCIP*) (size_t) jscip;

   assert(dialog != NULL);
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   entryname = (*env)->GetStringUTFChars(env, jentryname, &iscopy);
   if( entryname == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPdialogDisplayCompletions(dialog, scip, entryname) );

   (*env)->ReleaseStringUTFChars(env, jentryname, entryname);
}

/** gets the name of the current path in the dialog tree, separated by the given character */
JNIEXPORT
void JNISCIPDIALOG(dialogGetPath)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog,            /**< dialog handler */
   jchar                 jsepchar,           /**< separation character to insert in path */
   jstring               jpath               /**< string buffer to store the path */
   )
{
   SCIPerrorMessage("method dialogGetPath is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );
}

/** gets the command name of the dialog */
JNIEXPORT
jstring JNISCIPDIALOG(dialogGetName)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog             /**< dialog handler */
   )
{
   SCIPerrorMessage("method dialogGetName is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** gets the description of the dialog */
JNIEXPORT
jstring JNISCIPDIALOG(dialogGetDesc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog             /**< dialog handler */
   )
{
   SCIPerrorMessage("method dialogGetDesc is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** returns whether the dialog is a sub menu */
JNIEXPORT
jboolean JNISCIPDIALOG(dialogIsSubmenu)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog             /**< dialog handler */
   )
{
   SCIP_DIALOG* dialog;

   /* convert JNI pointer into C pointer */
   dialog = (SCIP_DIALOG*) (size_t) jdialog;
   assert(dialog != NULL);

   return (jlong) SCIPdialogIsSubmenu(dialog);
}

/** gets the parent dialog of the given dialog */
JNIEXPORT
jlong JNISCIPDIALOG(dialogGetParent)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog             /**< dialog handler */
   )
{
   SCIP_DIALOG* dialog;

   /* convert JNI pointer into C pointer */
   dialog = (SCIP_DIALOG*) (size_t) jdialog;
   assert(dialog != NULL);

   return (jlong) (size_t) SCIPdialogGetParent(dialog);
}

/** gets the array of sub-dialogs associated with the given dialog */
JNIEXPORT
jlongArray JNISCIPDIALOG(dialogGetSubdialogs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog             /**< dialog handler */
   )
{
   SCIP_DIALOG* dialog;
   jlongArray jdia;
   int size;
   SCIP_DIALOG** dia;

   /* convert JNI pointer into C pointer */
   dialog = (SCIP_DIALOG*) (size_t) jdialog;
   assert(dialog != NULL);

   size = SCIPdialogGetNSubdialogs(dialog);

   jdia = (*env)->NewLongArray(env, size);

   if (jdia == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   dia = SCIPdialogGetSubdialogs(dialog);

   (*env)->SetLongArrayRegion(env, jdia, 0, size, (jlong*)(*dia));

   return jdia;
}

/** gets the number of sub-dialogs associated with the given dialog */
JNIEXPORT
jint JNISCIPDIALOG(dialogGetNSubdialogs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog             /**< dialog handler */
   )
{
   SCIP_DIALOG* dialog;

   /* convert JNI pointer into C pointer */
   dialog = (SCIP_DIALOG*) (size_t) jdialog;
   assert(dialog != NULL);

   return (jint) SCIPdialogGetNSubdialogs(dialog);
}

/** gets the user defined data associated with the given dialog */
JNIEXPORT
jlong JNISCIPDIALOG(dialogGetData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog             /**< dialog handler */
   )
{
   SCIP_DIALOG* dialog;

   /* convert JNI pointer into C pointer */
   dialog = (SCIP_DIALOG*) (size_t) jdialog;
   assert(dialog != NULL);

   return (jlong) (size_t) SCIPdialogGetData(dialog);
}

/** sets user data of dialog; user has to free old data in advance! */
JNIEXPORT
void JNISCIPDIALOG(dialogSetData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdialog,            /**< dialog handler */
   jlong                 jdialogdata         /**< new dialog user data */
   )
{
   SCIP_DIALOG* dialog;
   SCIP_DIALOGDATA* dialogdata;

   /* convert JNI pointer into C pointer */
   dialog = (SCIP_DIALOG*) (size_t) jdialog;
   dialogdata = (SCIP_DIALOGDATA*) (size_t) jdialogdata;

   assert(dialog != NULL);
   assert(dialogdata != NULL);

   SCIPdialogSetData(dialog, dialogdata);
}
