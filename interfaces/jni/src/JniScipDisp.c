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

/**@file   JniScipDisp.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP disp callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipDisp.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/pub_disp.h"

#include <string.h>

/** gets user data of display column */
JNIEXPORT
jlong JNISCIPDISP(dispGetData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdisp               /**< display column */
   )
{
   SCIP_DISP* disp;

   /* convert JNI pointer into C pointer */
   disp = (SCIP_DISP*) (size_t) jdisp;
   assert(disp != NULL);

   return (jlong) (size_t) SCIPdispGetData(disp);
}

/** sets user data of display column; user has to free old data in advance! */
JNIEXPORT
void JNISCIPDISP(dispSetData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdisp,              /**< display column */
   jlong                 jdispdata           /**< new display column user data */
   )
{
   SCIP_DISP* disp;
   SCIP_DISPDATA* dispdata;

   /* convert JNI pointer into C pointer */
   dispdata = (SCIP_DISPDATA*) (size_t) jdispdata;
   disp = (SCIP_DISP*) (size_t) jdisp;

   assert(disp != NULL);
   assert(dispdata != NULL);

   SCIPdispSetData(disp, dispdata);
}

/** gets name of display column */
JNIEXPORT
jstring JNISCIPDISP(dispGetName)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdisp               /**< display column */
   )
{
   SCIPerrorMessage("method dispGetName is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** gets description of display column */
JNIEXPORT
jstring JNISCIPDISP(dispGetDesc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdisp               /**< display column */
   )
{
   SCIPerrorMessage("method dispGetDesc is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** gets head line of display column */
JNIEXPORT
jstring JNISCIPDISP(dispGetHeader)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdisp               /**< display column */
   )
{
   SCIPerrorMessage("method dispGetHeader is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** gets width of display column */
JNIEXPORT
jint JNISCIPDISP(dispGetWidth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdisp               /**< display column */
   )
{
   SCIP_DISP* disp;

   /* convert JNI pointer into C pointer */
   disp = (SCIP_DISP*) (size_t) jdisp;
   assert(disp != NULL);

   return (jint) SCIPdispGetWidth(disp);
}

/** gets priority of display column */
JNIEXPORT
jint JNISCIPDISP(dispGetPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdisp               /**< display column */
   )
{
   SCIP_DISP* disp;

   /* convert JNI pointer into C pointer */
   disp = (SCIP_DISP*) (size_t) jdisp;
   assert(disp != NULL);

   return (jint) SCIPdispGetPriority(disp);
}

/** gets position of display column */
JNIEXPORT
jint JNISCIPDISP(dispGetPosition)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdisp               /**< display column */
   )
{
   SCIP_DISP* disp;

   /* convert JNI pointer into C pointer */
   disp = (SCIP_DISP*) (size_t) jdisp;
   assert(disp != NULL);

   return (jint) SCIPdispGetPosition(disp);
}

/** gets status of display column */
JNIEXPORT
jint JNISCIPDISP(dispGetStatus)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdisp               /**< display column */
   )
{
   SCIP_DISP* disp;

   /* convert JNI pointer into C pointer */
   disp = (SCIP_DISP*) (size_t) jdisp;
   assert(disp != NULL);

   return (jint) SCIPdispGetStatus(disp);
}

/** is display column initialized? */
JNIEXPORT
jboolean JNISCIPDISP(dispIsInitialized)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdisp               /**< display column */
   )
{
   SCIP_DISP* disp;

   /* convert JNI pointer into C pointer */
   disp = (SCIP_DISP*) (size_t) jdisp;
   assert(disp != NULL);

   return (jboolean) SCIPdispIsInitialized(disp);
}

/** displays a long integer in decimal form fitting in a given width */
JNIEXPORT
void JNISCIPDISP(dispLongint)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jmessagehdlr,       /**< message handler */
   jlong                 jfile,              /**< output stream */
   jlong                  jval,               /**< value to display */
   jint                  jwidth              /**< width to fit into */
   )
{
   SCIP_MESSAGEHDLR* messagehdlr;
   FILE* file;

   /* convert JNI pointer into C pointer */
   messagehdlr = (SCIP_MESSAGEHDLR*) (size_t) jmessagehdlr;
   file = (FILE*) (size_t) jfile;

   assert(messagehdlr != NULL);
   assert(file != NULL);

   SCIPdispLongint(messagehdlr, file, (SCIP_Longint)jval, (jint)jwidth);
}

/** displays an integer in decimal form fitting in a given width */
JNIEXPORT
void JNISCIPDISP(dispInt)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jmessagehdlr,       /**< message handler */
   jlong                 jfile,              /**< output stream */
   jint                  jval,               /**< value to display */
   jint                  jwidth              /**< width to fit into */
   )
{
   SCIP_MESSAGEHDLR* messagehdlr;
   FILE* file;

   /* convert JNI pointer into C pointer */
   messagehdlr = (SCIP_MESSAGEHDLR*) (size_t) jmessagehdlr;
   file = (FILE*) (size_t) jfile;

   assert(messagehdlr != NULL);
   assert(file != NULL);

   SCIPdispInt(messagehdlr, file, (int)jval, (jint)jwidth);
}

/** displays a time value fitting in a given width */
JNIEXPORT
void JNISCIPDISP(dispTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jmessagehdlr,       /**< message handler */
   jlong                 jfile,              /**< output stream */
   jdouble               jval,               /**< value to display */
   jint                  jwidth              /**< width to fit into */
   )
{
   SCIP_MESSAGEHDLR* messagehdlr;
   FILE* file;

   /* convert JNI pointer into C pointer */
   messagehdlr = (SCIP_MESSAGEHDLR*) (size_t) jmessagehdlr;
   file = (FILE*) (size_t) jfile;

   assert(messagehdlr != NULL);
   assert(file != NULL);

   SCIPdispTime(messagehdlr, file, (SCIP_Real)jval, (jint)jwidth);
}
