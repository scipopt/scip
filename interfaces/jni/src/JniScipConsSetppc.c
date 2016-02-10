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

/**@file   JniScipConsSetppc.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP constraint partitioning / packing / covering callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipConsSetppc.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_setppc.h"

#include <string.h>

/** creates the handler for set partitioning / packing / covering constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSSETPPC(includeConshdlrSetppc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrSetppc(scip) );
}

/** add variable to setppc constraints */
static
SCIP_RETCODE addVars(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   jlongArray            jvars,              /**< array with variables of constraint entries */
   int                   nvars               /**< number of variables */
   )
{
   jlong* vars;
   int v;

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, nvars, vars);      

   for( v = 0; v < nvars; ++v )
   {
      JNISCIP_CALL( SCIPaddCoefSetppc(scip, cons, (SCIP_VAR*)(size_t)vars[v]) );
   }
   
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** creates and captures a set partitioning constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method {@link releaseCons()}
 */
JNIEXPORT
jlong JNISCIPCONSSETPPC(createConsSetpart)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnvars,             /**< number of variables in the constraint */
   jlongArray            jvars,              /**< array with variables of constraint entries */
   jboolean              jinitial,           /**< should the LP relaxation of constraint be in the initial LP?
					     *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   jboolean              jseparate,          /**< should the constraint be separated during LP processing?
					     *   Usually set to TRUE. */
   jboolean              jenforce,           /**< should the constraint be enforced during node processing?
					     *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   jboolean              jcheck,             /**< should the constraint be checked for feasibility?
					     *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   jboolean              jpropagate,         /**< should the constraint be propagated during node processing?
					     *   Usually set to TRUE. */
   jboolean              jlocal,             /**< is constraint only valid locally?
					     *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   jboolean              jmodifiable,        /**< is constraint modifiable (subject to column generation)?
					     *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
					     *   adds coefficients to this constraint. */
   jboolean              jdynamic,           /**< is constraint subject to aging?
					     *   Usually set to FALSE. Set to TRUE for own cuts which
					     *   are seperated as constraints. */
   jboolean              jremovable,         /**< should the relaxation be removed from the LP due to aging or cleanup?
 					     *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   jboolean              jstickingatnode     /**< should the constraint always be kept at the node where it was added, even
					     *   if it may be moved to a more global node?
					     *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   const char* name;
   int nvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, NULL);
   if( name == NULL )
      SCIPABORT();

   /* create an empty set partitioning constraint */
   JNISCIP_CALL( SCIPcreateConsSetpart(scip, &cons, name, 0, NULL,
	 (SCIP_Bool) jinitial, (SCIP_Bool) jseparate, (SCIP_Bool) jenforce, (SCIP_Bool) jcheck, (SCIP_Bool) jpropagate,
	 (SCIP_Bool) jlocal, (SCIP_Bool) jmodifiable, (SCIP_Bool) jdynamic, (SCIP_Bool) jremovable, (SCIP_Bool) jstickingatnode) );

   /* convert JNI integer into integer */
   nvars = (int) jnvars;

   /* add items */
   if( nvars > 0 )
   {
      JNISCIP_CALL( addVars(env, jobj, scip, cons, jvars, nvars) );
   }

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** creates and captures a set packing constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method {@link releaseCons()}
 */
JNIEXPORT
jlong JNISCIPCONSSETPPC(createConsSetpack)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnvars,             /**< number of variables in the constraint */
   jlongArray            jvars,              /**< array with variables of constraint entries */
   jboolean              jinitial,           /**< should the LP relaxation of constraint be in the initial LP?
					     *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   jboolean              jseparate,          /**< should the constraint be separated during LP processing?
					     *   Usually set to TRUE. */
   jboolean              jenforce,           /**< should the constraint be enforced during node processing?
					     *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   jboolean              jcheck,             /**< should the constraint be checked for feasibility?
					     *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   jboolean              jpropagate,         /**< should the constraint be propagated during node processing?
					     *   Usually set to TRUE. */
   jboolean              jlocal,             /**< is constraint only valid locally?
					     *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   jboolean              jmodifiable,        /**< is constraint modifiable (subject to column generation)?
					     *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
					     *   adds coefficients to this constraint. */
   jboolean              jdynamic,           /**< is constraint subject to aging?
					     *   Usually set to FALSE. Set to TRUE for own cuts which
					     *   are seperated as constraints. */
   jboolean              jremovable,         /**< should the relaxation be removed from the LP due to aging or cleanup?
					     *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   jboolean              jstickingatnode     /**< should the constraint always be kept at the node where it was added, even
					     *   if it may be moved to a more global node?
					     *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   const char* name;
   jboolean iscopy;
   int nvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   
   /* convert JNI integer into integer */
   nvars = (int) jnvars;

   JNISCIP_CALL( SCIPcreateConsSetpack(scip, &cons, name, 0, NULL,
	 (SCIP_Bool) jinitial, (SCIP_Bool) jseparate, (SCIP_Bool) jenforce, (SCIP_Bool) jcheck, (SCIP_Bool) jpropagate,
	 (SCIP_Bool) jlocal, (SCIP_Bool) jmodifiable, (SCIP_Bool) jdynamic, (SCIP_Bool) jremovable, (SCIP_Bool) jstickingatnode) );

   /* add items */
   if( nvars > 0 )
   {
      JNISCIP_CALL( addVars(env, jobj, scip, cons, jvars, nvars) );
   }

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** creates and captures a set covering constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method {@link releaseCons()}
 */
JNIEXPORT
jlong JNISCIPCONSSETPPC(createConsSetcover)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnvars,             /**< number of variables in the constraint */
   jlongArray            jvars,              /**< array with variables of constraint entries */
   jboolean              jinitial,           /**< should the LP relaxation of constraint be in the initial LP?
					     *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   jboolean              jseparate,          /**< should the constraint be separated during LP processing?
					     *   Usually set to TRUE. */
   jboolean              jenforce,           /**< should the constraint be enforced during node processing?
					     *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   jboolean              jcheck,             /**< should the constraint be checked for feasibility?
					     *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   jboolean              jpropagate,         /**< should the constraint be propagated during node processing?
					     *   Usually set to TRUE. */
   jboolean              jlocal,             /**< is constraint only valid locally?
					     *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   jboolean              jmodifiable,        /**< is constraint modifiable (subject to column generation)?
					     *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
					     *   adds coefficients to this constraint. */
   jboolean              jdynamic,           /**< is constraint subject to aging?
					     *   Usually set to FALSE. Set to TRUE for own cuts which
					     *   are seperated as constraints. */
   jboolean              jremovable,         /**< should the relaxation be removed from the LP due to aging or cleanup?
					     *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   jboolean              jstickingatnode     /**< should the constraint always be kept at the node where it was added, even
					     *   if it may be moved to a more global node?
					     *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   const char* name;
   jboolean iscopy;
   int nvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   /* convert JNI integer into integer */
   nvars = (int) jnvars;

   JNISCIP_CALL( SCIPcreateConsSetcover(scip, &cons, name, 0, NULL,
	 (SCIP_Bool) jinitial, (SCIP_Bool) jseparate, (SCIP_Bool) jenforce, (SCIP_Bool) jcheck, (SCIP_Bool) jpropagate,
	 (SCIP_Bool) jlocal, (SCIP_Bool) jmodifiable, (SCIP_Bool) jdynamic, (SCIP_Bool) jremovable, (SCIP_Bool) jstickingatnode) );

   /* add items */
   if( nvars > 0 )
   {
      JNISCIP_CALL( addVars(env, jobj, scip, cons, jvars, nvars) );
   }

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** adds coefficient in set partitioning / packing / covering constraint */
JNIEXPORT
void JNISCIPCONSSETPPC(addCoefSetppc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint data */
   jlong                 jvar                /**< variable to add to the constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddCoefSetppc(scip, cons, var) );
}

/** gets number of variables in set partitioning / packing / covering constraint */
JNIEXPORT
jint JNISCIPCONSSETPPC(getNVarsSetppc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   int nvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   nvars = SCIPgetNVarsSetppc(scip, cons);

   return (jint)nvars;
}

/** gets array of variables in set partitioning / packing / covering constraint */
JNIEXPORT
jlongArray JNISCIPCONSSETPPC(getVarsSetppc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   int nvars;
   jlongArray jvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   nvars = SCIPgetNVarsSetppc(scip, cons);
   vars = SCIPgetVarsSetppc(scip, cons);

   /* create jlong Array */
   jvars = (*env)->NewLongArray(env, nvars);

   /* fill long array with SCIP variable pointers */
   (*env)->SetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);

   return jvars;
}

/** gets the dual solution of the set partitioning / packing / covering constraint in the current LP */
JNIEXPORT
jdouble JNISCIPCONSSETPPC(getDualsolSetppc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */

   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Real dualsol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   dualsol = SCIPgetDualsolSetppc(scip, cons);

   return (jdouble)dualsol;
}

/** gets the dual farkas value of the set partitioning / packing / covering constraint in the current infeasible LP */
JNIEXPORT
jdouble JNISCIPCONSSETPPC(getDualfarkasSetppc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Real dualfarkas;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   dualfarkas = SCIPgetDualfarkasSetppc(scip, cons);

   return (jdouble)dualfarkas;
}

/** returns the linear relaxation of the given set partitioning / packing / covering constraint; may return NULL if no
 *  LP row was yet created; the user must not modify the row!
 */
JNIEXPORT
jlong JNISCIPCONSSETPPC(getRowSetppc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   row = SCIPgetRowSetppc(scip, cons);

   return (jlong)(size_t)row;
}

/** returns current number of variables fixed to one in the constraint  */
JNIEXPORT
jint JNISCIPCONSSETPPC(getNFixedonesSetppc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   int nvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   nvars = SCIPgetNFixedonesSetppc(scip, cons);

   return (jint)nvars;
}

/** returns current number of variables fixed to zero in the constraint  */
JNIEXPORT
jint JNISCIPCONSSETPPC(getNFixedzerosSetppc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   int nvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   nvars = SCIPgetNFixedzerosSetppc(scip, cons);

   return (jint)nvars;
}
