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

/**@file   JniScipEventhdlr.c
 * @brief  JNI SCIP callable library
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipEventhdlr.h"

#include "def.h"
#include "scip/scip.h"

/** event handler data */
struct SCIP_EventhdlrData
{
   JNIEnv*               env;                /**< JNI environment variable */
   jobject               jobj;               /**< JNI class pointer */
   jclass                javaeventhdlr;      /**< java event handler object class */
   jmethodID             scip_free;          /**<  method id of the free method */
   jmethodID             scip_init;          /**<  method id of the init method */
   jmethodID             scip_exit;          /**<  method id of the exit method */
   jmethodID             scip_initsol;       /**<  method id of the initsol method */
   jmethodID             scip_exitsol;       /**<  method id of the exitol process method */
   jmethodID             scip_delete;        /**<  method id of the deletion method */
   jmethodID             scip_exec;          /**<  method id of the execution method */
   jboolean              deleteobject;       /**< should the event handler object be deleted when eventhdlristic is freed? */
};

/** copy method for event handler plugins (called when SCIP copies plugins)
 *
 *  @todo implement copy method
 *
 *  ??????????????????????????
 */
#define eventhdlrCopyJava NULL

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventhdlrFreeJava)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   JNIEnv* env;
   jobject jobj;
   jlong jscip;
   jlong jeventhdlr;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->javaeventhdlr != NULL);
   assert(eventhdlrdata->scip_free != NULL);

   jscip = (jlong) (size_t) scip;
   jeventhdlr = (jlong) (size_t) eventhdlr;

   env = eventhdlrdata->env;
   jobj = eventhdlrdata->jobj;

   /* call virtual method of eventhdlr object */
   JNISCIP_CALL( (*env)->CallIntMethod(env, jobj, eventhdlrdata->scip_free, jscip, jeventhdlr) );

   /* free eventhdlr object
    *
    * ??????????????????????????
    */
   // if( eventhdlrdata->deleteobject )
   //    delete eventhdlrdata->objeventhdlr;

   /* free eventhdlr data */
   SCIPfreeMemory(scip, &eventhdlrdata)
   SCIPeventhdlrSetData(eventhdlr, NULL); /*lint !e64*/

   return SCIP_OKAY;
}


/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventhdlrInitJava)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;


   JNIEnv* env;
   jobject jobj;
   jlong jscip;
   jlong jeventhdlr;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->javaeventhdlr != NULL);
   assert(eventhdlrdata->scip_init != NULL);

   jscip = (jlong) (size_t) scip;
   jeventhdlr = (jlong) (size_t) eventhdlr;

   env = eventhdlrdata->env;
   jobj = eventhdlrdata->jobj;

   /* call virtual method of eventhdlr object */
   JNISCIP_CALL( (*env)->CallIntMethod(env, jobj, eventhdlrdata->scip_init, jscip, jeventhdlr) );

   return SCIP_OKAY;
}


/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventhdlrExitJava)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   JNIEnv* env;
   jobject jobj;
   jlong jscip;
   jlong jeventhdlr;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->javaeventhdlr != NULL);
   assert(eventhdlrdata->scip_exit != NULL);

   jscip = (jlong) (size_t) scip;
   jeventhdlr = (jlong) (size_t) eventhdlr;

   env = eventhdlrdata->env;
   jobj = eventhdlrdata->jobj;

   /* call virtual method of eventhdlr object */
   JNISCIP_CALL( (*env)->CallIntMethod(env, jobj, eventhdlrdata->scip_exit, jscip, jeventhdlr) );

   return SCIP_OKAY;
}


/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventhdlrInitsolJava)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   JNIEnv* env;
   jobject jobj;
   jlong jscip;
   jlong jeventhdlr;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->javaeventhdlr != NULL);
   assert(eventhdlrdata->scip_initsol != NULL);

   jscip = (jlong) (size_t) scip;
   jeventhdlr = (jlong) (size_t) eventhdlr;

   env = eventhdlrdata->env;
   jobj = eventhdlrdata->jobj;

   /* call virtual method of eventhdlr object */
   JNISCIP_CALL( (*env)->CallIntMethod(env, jobj, eventhdlrdata->scip_initsol, jscip, jeventhdlr) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventhdlrExitsolJava)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   JNIEnv* env;
   jobject jobj;
   jlong jscip;
   jlong jeventhdlr;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->javaeventhdlr != NULL);
   assert(eventhdlrdata->scip_exitsol != NULL);

   jscip = (jlong) (size_t) scip;
   jeventhdlr = (jlong) (size_t) eventhdlr;

   env = eventhdlrdata->env;
   jobj = eventhdlrdata->jobj;

   /* call virtual method of eventhdlr object */
   JNISCIP_CALL( (*env)->CallIntMethod(env, jobj, eventhdlrdata->scip_exitsol, jscip, jeventhdlr) );

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_EVENTDELETE(eventhdlrDeleteJava)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   JNIEnv* env;
   jobject jobj;
   jlong jscip;
   jlong jeventhdlr;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->javaeventhdlr != NULL);
   assert(eventhdlrdata->scip_delete != NULL);

   jscip = (jlong) (size_t) scip;
   jeventhdlr = (jlong) (size_t) eventhdlr;

   env = eventhdlrdata->env;
   jobj = eventhdlrdata->jobj;

   /* call virtual method of eventhdlr object */
   JNISCIP_CALL( (*env)->CallIntMethod(env, jobj, eventhdlrdata->scip_delete, jscip, jeventhdlr) );

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventhdlrExecJava)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   JNIEnv* env;
   jobject jobj;
   jlong jscip;
   jlong jeventhdlr;
   jlong jevent;
   jlong jeventdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->javaeventhdlr != NULL);
   assert(eventhdlrdata->scip_exec != NULL);

   env = eventhdlrdata->env;
   jobj = eventhdlrdata->jobj;

   jscip = (jlong) (size_t) scip;
   jeventhdlr = (jlong) (size_t) eventhdlr;
   jevent = (jlong) (size_t) event;
   jeventdata = (jlong) (size_t) eventdata;

   /* call virtual method of eventhdlr object */
   JNISCIP_CALL( (*env)->CallIntMethod(env, jobj, eventhdlrdata->scip_exec, jscip, jeventhdlr, jevent, jeventdata) );

   return SCIP_OKAY;
}

/** creates the event handler for the given event handler object and includes it in SCIP */
JNIEXPORT
void JNISCIP(includeJavaEventhdlr)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jclass                javaeventhdlr,      /**< event handler object */
   jboolean              deleteobject        /**< should the event handler object be deleted when eventhdlristic is freed? */
   )
{
   SCIP* scip;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   jstring jname;
   jstring jdesc;
   jboolean iscopy;
   jfieldID fieldID;

   const char* name;
   const char* desc;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* create event handler data */
   JNISCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );

   /* store JNI environment */
   eventhdlrdata->env = env;
   eventhdlrdata->jobj = jobj;

   /* store java class pointer */
   eventhdlrdata->javaeventhdlr = javaeventhdlr;

   /* collect for all callback methods the methods id and store these in the event handler data */
   eventhdlrdata->scip_free = (*env)->GetMethodID(env, javaeventhdlr, "scip_free", "(JJ)V");
   eventhdlrdata->scip_init = (*env)->GetMethodID(env, javaeventhdlr, "scip_init", "(JJ)V");
   eventhdlrdata->scip_exit = (*env)->GetMethodID(env, javaeventhdlr, "scip_exit", "(JJ)V");
   eventhdlrdata->scip_initsol = (*env)->GetMethodID(env, javaeventhdlr, "scip_initsol", "(JJ)V");
   eventhdlrdata->scip_exitsol = (*env)->GetMethodID(env, javaeventhdlr, "scip_exitsol", "(JJ)V");
   eventhdlrdata->scip_delete = (*env)->GetMethodID(env, javaeventhdlr, "scip_delete", "(JJ)V");
   eventhdlrdata->scip_exec = (*env)->GetMethodID(env, javaeventhdlr, "scip_exec", "(JJJJ)V");

   eventhdlrdata->deleteobject = deleteobject;

   /* collect event handler name */
   fieldID = (*env)->GetFieldID(env, javaeventhdlr, "name", "Ljava/lang/String;");
   jname = (*env)->GetObjectField(env, jobj, fieldID);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   assert(name != NULL);
   assert(iscopy);


   /* collect event handler description */
   fieldID = (*env)->GetFieldID(env, javaeventhdlr, "desc", "Ljava/lang/String;");
   jdesc = (*env)->GetObjectField(env, jobj, fieldID);

   /* convert JNI string into const char* */
   desc = (*env)->GetStringUTFChars(env, jdesc, &iscopy);
   assert(desc != NULL);
   assert(iscopy);

   /* include event handler */
   JNISCIP_CALL( SCIPincludeEventhdlr(scip, name, desc,
         eventhdlrCopyJava,
         eventhdlrFreeJava, eventhdlrInitJava, eventhdlrExitJava,
         eventhdlrInitsolJava, eventhdlrExitsolJava, eventhdlrDeleteJava, eventhdlrExecJava,
         eventhdlrdata) ); /*lint !e429*/

   (*env)->ReleaseStringUTFChars(env, jname, name);
   (*env)->ReleaseStringUTFChars(env, jdesc, desc);
}
