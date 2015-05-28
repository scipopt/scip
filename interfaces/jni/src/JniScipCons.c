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

/**@file   JniScip.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipCons.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include <string.h>

/** gets name of constraint handler */
JNIEXPORT
jstring JNISCIPCONS(conshdlrGetName)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* name;
   jstring jname;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   /* get constraint name */
   name = SCIPconshdlrGetName(conshdlr);

   /* convert char* into jstring */
   jname = (*env)->NewStringUTF(env, name);

   return jname;
}

/** gets description of constraint handler */
JNIEXPORT
jstring JNISCIPCONS(conshdlrGetDesc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* desc;
   jstring jdesc;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   /* get constraint description */
   desc = SCIPconshdlrGetDesc(conshdlr);

   /* convert char* into jstring */
   jdesc = (*env)->NewStringUTF(env, desc);

   return jdesc;
}

/** gets user data of constraint handler */
JNIEXPORT
jlong JNISCIPCONS(conshdlrGetData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   return (jlong) (size_t) SCIPconshdlrGetData(conshdlr);
}

/** sets user data of constraint handler; user has to free old data in advance! */
JNIEXPORT
void JNISCIPCONS(conshdlrSetData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr,          /**< constraint handler */
   jlong                 jconshdlrdata       /**< new constraint handler user data */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   conshdlrdata = (SCIP_CONSHDLRDATA*) (size_t) jconshdlrdata;
   assert(conshdlr != NULL);

   SCIPconshdlrSetData(conshdlr, conshdlrdata);
}

/** gets array with active constraints of constraint handler; a constraint is active if it is global and was not removed
 *  during presolving or it was added locally (in that case the local flag is TRUE) and the current node belongs to the
 *  corresponding sub tree
 */
JNIEXPORT
jlongArray JNISCIPCONS(conshdlrGetConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nconss;
   jlongArray jconss;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);

   /* create jlongArray */
   jconss = (*env)->NewLongArray(env, nconss);

   /* fill long array with SCIP variable pointers */
   (*env)->SetLongArrayRegion(env, jconss, 0, nconss, (jlong*)conss);

   return jconss;
}

/** gets array with enforced constraints of constraint handler; this is local information */
JNIEXPORT
jlongArray JNISCIPCONS(conshdlrGetEnfoConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nconss;
   jlongArray jconss;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   conss = SCIPconshdlrGetEnfoConss(conshdlr);
   nconss = SCIPconshdlrGetNEnfoConss(conshdlr);

   /* create jlongArray */
   jconss = (*env)->NewLongArray(env, nconss);

   /* fill long array with SCIP variable pointers */
   (*env)->SetLongArrayRegion(env, jconss, 0, nconss, (jlong*)conss);

   return jconss;
}

/** gets array with checked constraints of constraint handler; this is local information */
JNIEXPORT
jlongArray JNISCIPCONS(conshdlrGetCheckConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nconss;
   jlongArray jconss;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   conss = SCIPconshdlrGetCheckConss(conshdlr);
   nconss = SCIPconshdlrGetNCheckConss(conshdlr);

   /* create jlongArray */
   jconss = (*env)->NewLongArray(env, nconss);

   /* fill long array with SCIP variable pointers */
   (*env)->SetLongArrayRegion(env, jconss, 0, nconss, (jlong*)conss);

   return jconss;
}

/** gets total number of existing transformed constraints of constraint handler */
JNIEXPORT
jint JNISCIPCONS(conshdlrGetNConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNConss(conshdlr);

   return (jint) num;
}

/** gets number of enforced constraints of constraint handler; this is local information */
JNIEXPORT
jint JNISCIPCONS(conshdlrGetNEnfoConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNEnfoConss(conshdlr);

   return (jint) num;
}

/** gets number of checked constraints of constraint handler; this is local information */
JNIEXPORT
jint JNISCIPCONS(conshdlrGetNCheckConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNCheckConss(conshdlr);

   return (jint) num;
}

/** gets number of active constraints of constraint handler
 *
 *  @note A constraint is active if it is global and was not removed or it was added locally (in that case the local
 *        flag is TRUE) and the current node belongs to the corresponding sub tree.
 */
JNIEXPORT
jint JNISCIPCONS(conshdlrGetNActiveConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNActiveConss(conshdlr);

   return (jint) num;
}

/** gets number of enabled constraints of constraint handler */
JNIEXPORT
jint JNISCIPCONS(conshdlrGetNEnabledConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNEnabledConss(conshdlr);

   return (jint) num;
}

/** gets time in seconds used for setting up this constraint handler for new stages */
JNIEXPORT
jdouble JNISCIPCONS(conshdlrGetSetupTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetSetupTime(conshdlr);

   return (jdouble) num;
}

/** gets time in seconds used for presolving in this constraint handler */
JNIEXPORT
jdouble JNISCIPCONS(conshdlrGetPresolTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetPresolTime(conshdlr);

   return (jdouble) num;
}

/** gets time in seconds used for separation in this constraint handler */
JNIEXPORT
jdouble JNISCIPCONS(conshdlrGetSepaTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetSepaTime(conshdlr);

   return (jdouble) num;
}

/** gets time in seconds used for LP enforcement in this constraint handler */
JNIEXPORT
jdouble JNISCIPCONS(conshdlrGetEnfoLPTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetEnfoLPTime(conshdlr);

   return (jdouble) num;
}

/** gets time in seconds used for pseudo enforcement in this constraint handler */
JNIEXPORT
jdouble JNISCIPCONS(conshdlrGetEnfoPSTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetEnfoPSTime(conshdlr);

   return (jdouble) num;
}

/** gets time in seconds used for propagation in this constraint handler */
JNIEXPORT
jdouble JNISCIPCONS(conshdlrGetPropTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetPropTime(conshdlr);

   return (jdouble) num;
}

/** gets time in seconds used for feasibility checking in this constraint handler */
JNIEXPORT
jdouble JNISCIPCONS(conshdlrGetCheckTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetCheckTime(conshdlr);

   return (jdouble) num;
}

/** gets time in seconds used for resolving propagation in this constraint handler */
JNIEXPORT
jdouble JNISCIPCONS(conshdlrGetRespropTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetRespropTime(conshdlr);

   return (jdouble) num;
}

/** gets number of calls to the constraint handler's separation method */
JNIEXPORT
jlong JNISCIPCONS(conshdlrGetNSepaCalls)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Longint num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNSepaCalls(conshdlr);

   return (jlong) num;
}

/** gets number of calls to the constraint handler's LP enforcing method */
JNIEXPORT
jlong JNISCIPCONS(conshdlrGetNEnfoLPCalls)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Longint num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNEnfoLPCalls(conshdlr);

   return (jlong) num;
}

/** gets number of calls to the constraint handler's pseudo enforcing method */
JNIEXPORT
jlong JNISCIPCONS(conshdlrGetNEnfoPSCalls)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   return (jlong) (size_t) SCIPconshdlrGetNEnfoPSCalls(conshdlr);
}

/** gets number of calls to the constraint handler's propagation method */
jlong JNISCIPCONS(conshdlrGetNPropCalls)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   return (jlong) (size_t) SCIPconshdlrGetNPropCalls(conshdlr);
}

/** gets number of calls to the constraint handler's checking method */
jlong JNISCIPCONS(conshdlrGetNCheckCalls)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   return (jlong) (size_t) SCIPconshdlrGetNCheckCalls(conshdlr);
}

/** gets number of calls to the constraint handler's resolve propagation method */
jlong JNISCIPCONS(conshdlrGetNRespropCalls)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   return (jlong) SCIPconshdlrGetNRespropCalls(conshdlr);
}

/** gets total number of times, this constraint handler detected a cutoff */
jlong JNISCIPCONS(conshdlrGetNCutoffs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   return (jlong) SCIPconshdlrGetNCutoffs(conshdlr);
}

/** gets total number of cuts found by this constraint handler */
jlong JNISCIPCONS(conshdlrGetNCutsFound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   return (jlong) SCIPconshdlrGetNCutsFound(conshdlr);
}

/** gets total number of cuts found by this constraint handler applied to lp */
jlong JNISCIPCONS(conshdlrGetNCutsApplied)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   return (jlong) SCIPconshdlrGetNCutsApplied(conshdlr);
}

/** gets total number of additional constraints added by this constraint handler */
jlong JNISCIPCONS(conshdlrGetNConssFound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   return (jlong) SCIPconshdlrGetNConssFound(conshdlr);
}

/** gets total number of domain reductions found by this constraint handler */
jlong JNISCIPCONS(conshdlrGetNDomredsFound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   return (jlong) SCIPconshdlrGetNDomredsFound(conshdlr);
}

/** gets number of children created by this constraint handler */
jlong JNISCIPCONS(conshdlrGetNChildren)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   return (jlong) SCIPconshdlrGetNChildren(conshdlr);
}

/** gets maximum number of active constraints of constraint handler existing at the same time */
jint JNISCIPCONS(conshdlrGetMaxNActiveConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetMaxNActiveConss(conshdlr);

   return (jint) num;
}

/** gets initial number of active constraints of constraint handler */
jint JNISCIPCONS(conshdlrGetStartNActiveConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetStartNActiveConss(conshdlr);

   return (jint) num;
}

/** gets number of variables fixed in presolving method of constraint handler */
jint JNISCIPCONS(conshdlrGetNFixedVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNFixedVars(conshdlr);

   return (jint) num;
}

/** gets number of variables aggregated in presolving method of constraint handler */
jint JNISCIPCONS(conshdlrGetNAggrVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNAggrVars(conshdlr);

   return (jint) num;
}

/** gets number of variable types changed in presolving method of constraint handler */
jint JNISCIPCONS(conshdlrGetNChgVarTypes)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNChgVarTypes(conshdlr);

   return (jint) num;
}

/** gets number of bounds changed in presolving method of constraint handler */
jint JNISCIPCONS(conshdlrGetNChgBds)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNChgBds(conshdlr);

   return (jint) num;
}

/** gets number of holes added to domains of variables in presolving method of constraint handler */
jint JNISCIPCONS(conshdlrGetNAddHoles)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNAddHoles(conshdlr);

   return (jint) num;
}

/** gets number of constraints deleted in presolving method of constraint handler */
jint JNISCIPCONS(conshdlrGetNDelConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNDelConss(conshdlr);

   return (jint) num;
}

/** gets number of constraints added in presolving method of constraint handler */
jint JNISCIPCONS(conshdlrGetNAddConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNAddConss(conshdlr);

   return (jint) num;
}

/** gets number of constraints upgraded in presolving method of constraint handler */
jint JNISCIPCONS(conshdlrGetNUpgdConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNUpgdConss(conshdlr);

   return (jint) num;
}

/** gets number of coefficients changed in presolving method of constraint handler */
jint JNISCIPCONS(conshdlrGetNChgCoefs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNChgCoefs(conshdlr);

   return (jint) num;
}

/** gets number of constraint sides changed in presolving method of constraint handler */
jint JNISCIPCONS(conshdlrGetNChgSides)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNChgSides(conshdlr);

   return (jint) num;
}

/** gets number of times the presolving method of the constraint handler was called and tried to find reductions */
jint JNISCIPCONS(conshdlrGetNPresolCalls)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetNPresolCalls(conshdlr);

   return (jint) num;
}

/** gets separation priority of constraint handler */
jint JNISCIPCONS(conshdlrGetSepaPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetSepaPriority(conshdlr);

   return (jint) num;
}

/** gets enforcing priority of constraint handler */
jint JNISCIPCONS(conshdlrGetEnfoPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetEnfoPriority(conshdlr);

   return (jint) num;
}

/** gets checking priority of constraint handler */
jint JNISCIPCONS(conshdlrGetCheckPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetCheckPriority(conshdlr);

   return (jint) num;
}

/** gets separation frequency of constraint handler */
jint JNISCIPCONS(conshdlrGetSepaFreq)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetSepaFreq(conshdlr);

   return (jint) num;
}

/** gets propagation frequency of constraint handler */
jint JNISCIPCONS(conshdlrGetPropFreq)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetPropFreq(conshdlr);

   return (jint) num;
}

/** gets frequency of constraint handler for eager evaluations in separation, propagation and enforcement */
jint JNISCIPCONS(conshdlrGetEagerFreq)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int num;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   num = SCIPconshdlrGetEagerFreq(conshdlr);

   return (jint) num;
}

/** needs constraint handler a constraint to be called? */
jboolean JNISCIPCONS(conshdlrNeedsCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool needscons;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   needscons = SCIPconshdlrNeedsCons(conshdlr);

   return (jboolean) needscons;
}

/** does the constraint handler perform presolving? */
jboolean JNISCIPCONS(conshdlrDoesPresolve)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool doespresolve;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   doespresolve = SCIPconshdlrDoesPresolve(conshdlr);

   return (jboolean) doespresolve;
}

/** should separation method be delayed, if other separators found cuts? */
jboolean JNISCIPCONS(conshdlrIsSeparationDelayed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool sepdelayed;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   sepdelayed = SCIPconshdlrIsSeparationDelayed(conshdlr);

   return (jboolean) sepdelayed;
}

/** should propagation method be delayed, if other propagators found reductions? */
jboolean JNISCIPCONS(conshdlrIsPropagationDelayed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool propdelayed;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   propdelayed = SCIPconshdlrIsPropagationDelayed(conshdlr);

   return (jboolean) propdelayed;
}

/** was LP separation method delayed at the last call? */
jboolean JNISCIPCONS(conshdlrWasLPSeparationDelayed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool sepdelayed;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   sepdelayed = SCIPconshdlrWasLPSeparationDelayed(conshdlr);

   return (jboolean) sepdelayed;
}

/** was primal solution separation method delayed at the last call? */
jboolean JNISCIPCONS(conshdlrWasSolSeparationDelayed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool sepdelayed;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   sepdelayed = SCIPconshdlrWasSolSeparationDelayed(conshdlr);

   return (jboolean) sepdelayed;
}

/** was propagation method delayed at the last call? */
jboolean JNISCIPCONS(conshdlrWasPropagationDelayed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool propdelayed;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   propdelayed = SCIPconshdlrWasPropagationDelayed(conshdlr);

   return (jboolean) propdelayed;
}

/** is constraint handler initialized? */
jboolean JNISCIPCONS(conshdlrIsInitialized)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool isinit;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   isinit = SCIPconshdlrIsInitialized(conshdlr);

   return (jboolean) isinit;
}

/** does the constraint handler have a copy function? */
jboolean JNISCIPCONS(conshdlrIsClonable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool isconable;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   isconable = SCIPconshdlrIsClonable(conshdlr);

   return (jboolean) isconable;
}

/** returns the timing mask of the propagation method of the constraint handler */
jint JNISCIPCONS(conshdlrGetPropTiming)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_PROPTIMING timing;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   timing = SCIPconshdlrGetPropTiming(conshdlr);

   return (jint) timing;
}

/** returns the timing mask of the presolving method of the constraint handler */
jint JNISCIPCONS(conshdlrGetPresolTiming)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jconshdlr           /**< constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_PRESOLTIMING timing;

   /* convert JNI pointer into C pointer */
   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   timing = SCIPconshdlrGetPresolTiming(conshdlr);

   return (jint) timing;
}

/** returns the name of the constraint */
JNIEXPORT
jstring JNISCIPCONS(consGetName)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP_CONS* cons;
   const char* name;
   jstring jname;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   /* get constraint name */
   name = SCIPconsGetName(cons);

   /* convert char* into jstring */
   jname = (*env)->NewStringUTF(env, name);

   return jname;
}

/** returns the position of constraint in the corresponding handler's conss array */
JNIEXPORT
jint JNISCIPCONS(consGetPos)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   int pos;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   pos = SCIPconsGetPos(cons);

   return (jint) pos;
}

/** returns the constraint handler of the constraint */
JNIEXPORT
jlong JNISCIPCONS(consGetHdlr)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   conshdlr = SCIPconsGetHdlr(cons);

   return (jlong)(size_t) conshdlr;
}

/** returns the constraint data field of the constraint */
JNIEXPORT
jlong JNISCIPCONS(consGetData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   return (jlong)(size_t) SCIPconsGetData(cons);
}

/** gets number of times, the constraint is currently captured */
JNIEXPORT
jint JNISCIPCONS(consGetNUses)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   int num;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   num = SCIPconsGetNUses(cons);

   return (jint) num;
}

/** for an active constraint, returns the depth in the tree at which the constraint was activated */
JNIEXPORT
jint JNISCIPCONS(consGetActiveDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   int depth;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   depth = SCIPconsGetActiveDepth(cons);

   return (jint) depth;
}

/** returns the depth in the tree at which the constraint is valid; returns INT_MAX, if the constraint is local
 *  and currently not active
 */
JNIEXPORT
jint JNISCIPCONS(consGetValidDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   int depth;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   depth = SCIPconsGetValidDepth(cons);

   return (jint) depth;
}

/** returns TRUE iff constraint is active in the current node */
JNIEXPORT
jboolean JNISCIPCONS(consIsActive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool isactive;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   isactive = SCIPconsIsActive(cons);

   return (jboolean) isactive;
}

/** returns TRUE iff constraint is enabled in the current node */
JNIEXPORT
jboolean JNISCIPCONS(consIsEnabled)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool isenabled;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   isenabled = SCIPconsIsEnabled(cons);

   return (jboolean) isenabled;
}

/** returns TRUE iff constraint's separation is enabled in the current node */
JNIEXPORT
jboolean JNISCIPCONS(consIsSeparationEnabled)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool issepen;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   issepen = SCIPconsIsSeparationEnabled(cons);

   return (jboolean) issepen;
}

/** returns TRUE iff constraint's propagation is enabled in the current node */
JNIEXPORT
jboolean JNISCIPCONS(consIsPropagationEnabled)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool ispropen;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   ispropen = SCIPconsIsPropagationEnabled(cons);

   return (jboolean) ispropen;
}

/** returns TRUE iff constraint is deleted or marked to be deleted */
JNIEXPORT
jboolean JNISCIPCONS(consIsDeleted)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool isdeleted;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   isdeleted = SCIPconsIsDeleted(cons);

   return (jboolean) isdeleted;
}

/** returns TRUE iff constraint is marked obsolete */
JNIEXPORT
jboolean JNISCIPCONS(consIsObsolete)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool isobsolete;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   isobsolete= SCIPconsIsObsolete(cons);

   return (jboolean) isobsolete;
}

/** gets age of constraint */
JNIEXPORT
jdouble JNISCIPCONS(consGetAge)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Real age;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   age = SCIPconsGetAge(cons);

   return (jdouble) age;
}

/** returns TRUE iff the LP relaxation of constraint should be in the initial LP */
JNIEXPORT
jboolean JNISCIPCONS(consIsInitial)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool initial;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   initial = SCIPconsIsInitial(cons);

   return (jboolean) initial;
}

/** returns TRUE iff constraint should be separated during LP processing */
JNIEXPORT
jboolean JNISCIPCONS(consIsSeparated)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool separated;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   separated = SCIPconsIsSeparated(cons);

   return (jboolean) separated;
}

/** returns TRUE iff constraint should be enforced during node processing */
JNIEXPORT
jboolean JNISCIPCONS(consIsEnforced)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool enforced;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   enforced = SCIPconsIsEnforced(cons);

   return (jboolean) enforced;
}

/** returns TRUE iff constraint should be checked for feasibility */
JNIEXPORT
jboolean JNISCIPCONS(consIsChecked)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool checked;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   checked = SCIPconsIsChecked(cons);

   return (jboolean) checked;
}

/** returns TRUE iff constraint should be propagated during node processing */
JNIEXPORT
jboolean JNISCIPCONS(consIsPropagated)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool propagated;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   propagated = SCIPconsIsPropagated(cons);

   return (jboolean) propagated;
}

/** returns TRUE iff constraint is globally valid */
JNIEXPORT
jboolean JNISCIPCONS(consIsGlobal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool global;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   global = SCIPconsIsGlobal(cons);

   return (jboolean) global;
}

/** returns TRUE iff constraint is only locally valid or not added to any (sub)problem */
JNIEXPORT
jboolean JNISCIPCONS(consIsLocal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool local;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   local = SCIPconsIsLocal(cons);

   return (jboolean) local;
}

/** returns TRUE iff constraint is modifiable (subject to column generation) */
JNIEXPORT
jboolean JNISCIPCONS(consIsModifiable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool modifiable;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   modifiable = SCIPconsIsModifiable(cons);

   return (jboolean) modifiable;
}

/** returns TRUE iff constraint is subject to aging */
JNIEXPORT
jboolean JNISCIPCONS(consIsDynamic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool dynamic;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   dynamic = SCIPconsIsDynamic(cons);

   return (jboolean) dynamic;
}

/** returns TRUE iff constraint's relaxation should be removed from the LP due to aging or cleanup */
JNIEXPORT
jboolean JNISCIPCONS(consIsRemovable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool removable;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   removable = SCIPconsIsRemovable(cons);

   return (jboolean) removable;
}

/** returns TRUE iff constraint's relaxation should be removed from the LP due to aging or cleanup */
JNIEXPORT
jboolean JNISCIPCONS(consIsStickingAtNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool stickingatnode;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   stickingatnode = SCIPconsIsStickingAtNode(cons);

   return (jboolean) stickingatnode;
}

/** returns TRUE iff constraint belongs to the global problem */
JNIEXPORT
jboolean JNISCIPCONS(consIsInProb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool inprob;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   inprob = SCIPconsIsInProb(cons);

   return (jboolean) inprob;
}

/** returns TRUE iff constraint is belonging to original space */
JNIEXPORT
jboolean JNISCIPCONS(consIsOriginal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool orig;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   orig = SCIPconsIsOriginal(cons);

   return (jboolean) orig;
}

/** returns TRUE iff constraint is belonging to transformed space */
JNIEXPORT
jboolean JNISCIPCONS(consIsTransformed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool transformed;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   transformed = SCIPconsIsTransformed(cons);

   return (jboolean) transformed;
}

/** returns TRUE iff roundings for variables in constraint are locked */
JNIEXPORT
jboolean JNISCIPCONS(consIsLockedPos)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool locked;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   locked = SCIPconsIsLockedPos(cons);

   return (jboolean) locked;
}

/** returns TRUE iff roundings for variables in constraint's negation are locked */
JNIEXPORT
jboolean JNISCIPCONS(consIsLockedNeg)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool locked;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   locked = SCIPconsIsLockedNeg(cons);

   return (jboolean) locked ;
}

/** returns TRUE iff roundings for variables in constraint or in constraint's negation are locked */
JNIEXPORT
jboolean JNISCIPCONS(consIsLocked)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool locked;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   locked = SCIPconsIsLocked(cons);

   return (jboolean) locked;
}

/** get number of times the roundings for variables in constraint are locked */
JNIEXPORT
jint JNISCIPCONS(consGetNLocksPos)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   int num;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   num = SCIPconsGetNLocksPos(cons);

   return (jint) num;
}

/** get number of times the roundings for variables in constraint's negation are locked */
JNIEXPORT
jint JNISCIPCONS(consGetNLocksNeg)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   int num;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   num = SCIPconsGetNLocksNeg(cons);

   return (jint) num;
}

/** returns if the constraint was already added to a SCIP instance */
JNIEXPORT
jboolean JNISCIPCONS(consIsAdded)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< JNI problem variable */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool added;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   added = SCIPconsIsAdded(cons);

   return (jboolean) added;
}
