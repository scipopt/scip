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

/**@file   JniScipConsSos1.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP sos1 constraint callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipConsSos1.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_sos1.h"

#include <string.h>

/** creates the handler for SOS1 constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSSOS1(includeConshdlrSOS1)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrSOS1(scip) );
}

/** creates and captures an SOS1 constraint
 *
 *  We set the constraint to not be modifable. If the weights are non
 *  NULL, the variables are ordered according to these weights (in
 *  ascending order).
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method {@link releaseCons()}
 */
JNIEXPORT
jlong JNISCIPCONSSOS1(createConsSOS1)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnvars,             /**< number of variables in the constraint */
   jlongArray            jvars,              /**< array with variables of constraint entries */
   jdoubleArray          jweights,           /**< weights determining the variable order, or NULL if natural order should be used */
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

   /* create empty sos1 constraint */
   JNISCIP_CALL( SCIPcreateConsSOS1(scip, &cons, name, 0, NULL, NULL, 
	 (SCIP_Bool) jinitial, (SCIP_Bool) jseparate, (SCIP_Bool) jenforce, (SCIP_Bool) jcheck, (SCIP_Bool) jpropagate,
	 (SCIP_Bool) jlocal, (SCIP_Bool) jdynamic, (SCIP_Bool) jremovable, (SCIP_Bool) jstickingatnode) ); 
   /* convert JNI integer into integer */
   nvars = (int) jnvars;
   
   /* add itmes */
   if( nvars > 0 )
   {
      jlong* vars;
      jdouble* weights;
      int v;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      JNISCIP_CALL( SCIPallocBufferArray(scip, &weights, nvars) );
      
      (*env)->GetLongArrayRegion(env, jvars, 0, nvars, vars);
      (*env)->GetDoubleArrayRegion(env, jweights, 0, nvars, weights);
      
      for( v = 0; v < nvars; ++v )
      {
         JNISCIP_CALL( SCIPaddVarSOS1(scip, cons, (SCIP_VAR*)(size_t)(vars[v]), (SCIP_Real)weights[v]) );
      }

      SCIPfreeBufferArray(scip, &weights);
      SCIPfreeBufferArray(scip, &vars);
   }
   
   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** adds variable to SOS1 constraint, the position is determined by the given weight */
JNIEXPORT
void JNISCIPCONSSOS1(addVarSOS1)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jlong                 jvar,               /**< variable to add to the constraint */
   jdouble               jweight             /**< weight determining position of variable */
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
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddVarSOS1(scip, cons, var, (SCIP_Real) jweight) );
}

/** appends variable to SOS1 constraint */
JNIEXPORT
void JNISCIPCONSSOS1(appendVarSOS1)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
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
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPappendVarSOS1(scip, cons, var) );
}

/** gets number of variables in SOS1 constraint */
JNIEXPORT
jint JNISCIPCONSSOS1(getNVarsSOS1)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(scip != NULL);

   num = SCIPgetNVarsSOS1(scip, cons);

   return (jint) num;
}

/** gets array of variables in SOS1 constraint */
JNIEXPORT
jlongArray JNISCIPCONSSOS1(getVarsSOS1)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   jlongArray jvars;
   int nvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(scip != NULL);

   nvars = SCIPgetNVarsSOS1(scip, cons);
   vars = SCIPgetVarsSOS1(scip, cons);

   /* create jdouble Array */
   jvars = (*env)->NewLongArray(env, nvars);

   /* fill long array with SCIP variable pointers */
   (*env)->SetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);

   return jvars;
}

/** gets array of weights in SOS1 constraint (or NULL if not existent) */
JNIEXPORT
jdoubleArray JNISCIPCONSSOS1(getWeightsSOS1)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Real* weights;
   jdoubleArray jweights;
   int nweights;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   nweights = SCIPgetNVarsSOS1(scip, cons);
   weights = SCIPgetWeightsSOS1(scip, cons);

   /* create jdouble Array */
   jweights = (*env)->NewDoubleArray(env, nweights);

   /* fill long array with SCIP variable pointers */
   (*env)->SetDoubleArrayRegion(env, jweights, 0, nweights, (jdouble*)weights);

   return jweights;
}
