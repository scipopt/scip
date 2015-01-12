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

/**@file   JniScipConsXor.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP sos2 constraint callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipConsXor.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_xor.h"

#include <string.h>

/** creates the handler for xor constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSXOR(includeConshdlrXor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrXor(scip) );
}

/** creates and captures a xor constraint x_0 xor ... xor x_{k-1} = rhs
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSXOR(createConsXor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jboolean              jrhs,               /**< right hand side of the constraint */
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
   SCIP_VAR** vars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, NULL);
   if( name == NULL )
      SCIPABORT();

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, jnvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, jnvars, (jlong*)(*vars));

   JNISCIP_CALL( SCIPcreateConsXor(scip, &cons, name, (SCIP_Real) jrhs, (int)jnvars, vars,
         (SCIP_Bool) jinitial, (SCIP_Bool) jseparate, (SCIP_Bool) jenforce, (SCIP_Bool) jcheck, (SCIP_Bool) jpropagate,
         (SCIP_Bool) jlocal, (SCIP_Bool) jmodifiable, (SCIP_Bool) jdynamic, (SCIP_Bool) jremovable, (SCIP_Bool) jstickingatnode) );

   SCIPfreeBufferArray(scip, &vars);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** creates and captures a xor constraint x_0 xor ... xor x_{k-1} = rhs
 *  with all constraint flags set to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSXOR(createConsBasicXor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jboolean              jrhs,               /**< right hand side of the constraint */
   jint                  jnvars,             /**< number of operator variables in the constraint */
   jlongArray            jvars               /**< array with operator variables of constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   const char* name;
   SCIP_VAR** vars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, NULL);
   if( name == NULL )
      SCIPABORT();

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)jnvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)jnvars, (jlong*)(*vars));

   JNISCIP_CALL( SCIPcreateConsBasicXor(scip, &cons, name, (SCIP_Bool)jrhs, (int)jnvars, vars ) );

   SCIPfreeBufferArray(scip, &vars);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) cons;
}

/** gets number of variables in xor constraint */
JNIEXPORT
jint JNISCIPCONSXOR(getNVarsXor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   return (jint) SCIPgetNVarsXor(scip, cons);
}

/** gets array of variables in xor constraint */
JNIEXPORT
jlongArray JNISCIPCONSXOR(getVarsXor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   jlongArray jvars;
   int size;
   SCIP_VAR** vars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   size = SCIPgetNVarsXor(scip, cons);

   jvars = (*env)->NewLongArray(env, size);

   if (jvars == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   vars = SCIPgetVarsXor(scip, cons);

   (*env)->SetLongArrayRegion(env, jvars, 0, size, (jlong*)(*vars));

   return jvars;
}

/** gets the right hand side of the xor constraint */
JNIEXPORT
jboolean JNISCIPCONSXOR(getRhsXor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   return (jboolean) SCIPgetRhsXor(scip, cons);
}
