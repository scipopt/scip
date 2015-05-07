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

/**@file   JniScipMessageDefault.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP message default callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipMessageDefault.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/message_default.h"

#include <string.h>

/** Create default message handler. To free the message handler use SCIPmessagehdlrFree() */
JNIEXPORT
jlong JNISCIPMESSAGEDEFAULT(createMessagehdlrDefault)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jboolean              jbufferedoutput,    /**< should the output be buffered up to the next newline? */
   jstring               jfilename,          /**< name of log file, or NULL (stdout) */
   jboolean              jquiet              /**< should screen messages be suppressed? */
   )
{
   SCIP_MESSAGEHDLR* messagehdlr;
   const char* filename;
   jboolean iscopy;

   /* convert JNI string into const char* */
   filename = (*env)->GetStringUTFChars(env, jfilename, &iscopy);
   assert(iscopy);

   JNISCIP_CALL( SCIPcreateMessagehdlrDefault(&messagehdlr, (SCIP_Bool)jbufferedoutput, filename, (SCIP_Bool)jquiet) );

   (*env)->ReleaseStringUTFChars(env, jfilename, filename);

   return (jlong) (size_t) messagehdlr;
}
