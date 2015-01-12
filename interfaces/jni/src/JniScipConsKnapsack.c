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

/**@file   JniScipConsKnapsack.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP Constraint handler for knapsack constraints of the form  \f$a^T x \le b\f$, x binary and \f$a \ge 0\f$.
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipConsKnapsack.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_knapsack.h"

#include <string.h>

/** creates the handler for knapsack constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSKNAPSACK(includeConshdlrKnapsack)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrKnapsack(scip) );
}

/** creates and captures a knapsack constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method {@link releaseCons()}
 */
JNIEXPORT
jlong JNISCIPCONSKNAPSACK(createConsKnapsack)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnvars,             /**< number of items in the knapsack */
   jlongArray            jvars,              /**< array with item variables */
   jlongArray            jweights,           /**< array with item weights */
   jlong                 jcapacity,          /**< capacity of knapsack */
   jboolean              initial,            /**< should the LP relaxation of constraint be in the initial LP?
					      *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   jboolean              separate,           /**< should the constraint be separated during LP processing?
					      *   Usually set to TRUE. */
   jboolean              enforce,            /**< should the constraint be enforced during node processing?
					      *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   jboolean              check,              /**< should the constraint be checked for feasibility?
					      *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   jboolean              propagate,          /**< should the constraint be propagated during node processing?
					      *   Usually set to TRUE. */
   jboolean              local,              /**< is constraint only valid locally?
					      *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   jboolean              modifiable,         /**< is constraint modifiable (subject to column generation)?
					      *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
					      *   adds coefficients to this constraint. */
   jboolean              dynamic,            /**< is constraint subject to aging?
					      *   Usually set to FALSE. Set to TRUE for own cuts which
					      *   are seperated as constraints. */
   jboolean              removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
					      *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   jboolean              stickingatnode      /**< should the constraint always be kept at the node where it was added, even
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

   /* create empty knapsack constraint */
   JNISCIP_CALL( SCIPcreateConsKnapsack(scip, &cons, name, 0, NULL, NULL, (SCIP_Longint)jcapacity,
         (SCIP_Bool) initial, (SCIP_Bool) separate, (SCIP_Bool) enforce, (SCIP_Bool) check, (SCIP_Bool) propagate,
         (SCIP_Bool) local, (SCIP_Bool) modifiable, (SCIP_Bool) dynamic, (SCIP_Bool) removable, (SCIP_Bool) stickingatnode) );

   /* convert JNI integer into integer */
   nvars = (int) jnvars;

   /* add itmes */
   if( nvars > 0 )
   {
      jlong* vars;
      jlong* weights;
      int v;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      JNISCIP_CALL( SCIPallocBufferArray(scip, &weights, nvars) );

      (*env)->GetLongArrayRegion(env, jvars, 0, nvars, vars);
      (*env)->GetLongArrayRegion(env, jweights, 0, nvars, weights);

      for( v = 0; v < nvars; ++v )
      {
         JNISCIP_CALL( SCIPaddCoefKnapsack(scip, cons, (SCIP_VAR*)(size_t)(vars[v]), (SCIP_Longint)weights[v]) );
      }

      SCIPfreeBufferArray(scip, &weights);
      SCIPfreeBufferArray(scip, &vars);
   }

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** creates and captures a knapsack constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsKnapsack(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsKnapsack() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSKNAPSACK(createConsBasicKnapsack)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnvars,             /**< number of items in the knapsack */
   jlongArray            jvars,              /**< array with item variables */
   jlongArray            jweights,           /**< array with item weights */
   jlong                 jcapacity           /**< capacity of knapsack */
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

   /* create empty knapsack constraint */
   JNISCIP_CALL( SCIPcreateConsBasicKnapsack(scip, &cons, name, 0, NULL, NULL, (SCIP_Longint)jcapacity) );
   /* convert JNI integer into integer */
   nvars = (int) jnvars;

   /* add itmes */
   if( nvars > 0 )
   {
      jlong* vars;
      jlong* weights;
      int v;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      JNISCIP_CALL( SCIPallocBufferArray(scip, &weights, nvars) );

      (*env)->GetLongArrayRegion(env, jvars, 0, nvars, vars);
      (*env)->GetLongArrayRegion(env, jweights, 0, nvars, weights);

      for( v = 0; v < nvars; ++v )
      {
         JNISCIP_CALL( SCIPaddCoefKnapsack(scip, cons, (SCIP_VAR*)(size_t)(vars[v]), (SCIP_Longint)weights[v]) );
      }

      SCIPfreeBufferArray(scip, &weights);
      SCIPfreeBufferArray(scip, &vars);
   }

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** adds new item to knapsack constraint */
JNIEXPORT
void JNISCIPCONSKNAPSACK(addCoefKnapsack)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint data */
   jlong                 jvar,               /**< item variable */
   jlong                 jweight             /**< item weight */

   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_CONS* cons;
   SCIP_Longint weight;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   /* convert JNI pointer into C pointer */
   weight = (SCIP_Longint) jweight;

   JNISCIP_CALL( SCIPaddCoefKnapsack(scip, cons, var, weight) );
}


/** gets the capacity of the knapsack constraint */
JNIEXPORT
jlong JNISCIPCONSKNAPSACK(getCapacityKnapsack)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Longint capacity;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   capacity = SCIPgetCapacityKnapsack(scip, cons);

   return (jlong)capacity;
}

/** changes capacity of the knapsack constraint
 *
 * @note  This method can only be called during problem creation stage (SCIP_STAGE_PROBLEM)
 */
JNIEXPORT
void JNISCIPCONSKNAPSACK(chgCapacityKnapsack)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint data */
   jlong                 jcapacity           /**< new knapsack capacity */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPchgCapacityKnapsack(scip, cons, (SCIP_Longint)jcapacity) );
}

/** gets the number of items in the knapsack constraint */
JNIEXPORT
jint JNISCIPCONSKNAPSACK(getNVarsKnapsack)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   int items;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   items = SCIPgetNVarsKnapsack(scip, cons);

   return (jint)items;
}

/** gets the array of variables in the knapsack constraint; the user must not modify this array! */
JNIEXPORT
jlongArray JNISCIPCONSKNAPSACK(getVarsKnapsack)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
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
   assert(cons != NULL);

   nvars = SCIPgetNVarsKnapsack(scip, cons);
   vars = SCIPgetVarsKnapsack(scip, cons);

   /* create jdouble Array */
   jvars = (*env)->NewLongArray(env, nvars);

   /* fill long array with SCIP variable pointers */
   (*env)->SetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);

   return jvars;
 }

/** gets the array of weights in the knapsack constraint; the user must not modify this array! */
JNIEXPORT
jlongArray JNISCIPCONSKNAPSACK(getWeightsKnapsack)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Longint* weights;
   jlongArray jweights;
   int nweights;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   nweights = SCIPgetNVarsKnapsack(scip, cons);
   weights = SCIPgetWeightsKnapsack(scip, cons);

   /* create jdouble Array */
   jweights = (*env)->NewLongArray(env, nweights);

   /* fill long array with SCIP variable pointers */
   (*env)->SetLongArrayRegion(env, jweights, 0, nweights, (jlong*)weights);

   return jweights;
}


/** gets the dual solution of the knapsack constraint in the current LP  */
JNIEXPORT
jdouble JNISCIPCONSKNAPSACK(getDualsolKnapsack)(
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

   dualsol = SCIPgetDualsolKnapsack(scip, cons);

   return (jdouble)dualsol;
}

/** gets the dual farkas value of the knapsack constraint in the current infeasible LP */
JNIEXPORT
jdouble JNISCIPCONSKNAPSACK(getDualfarkasKnapsack)(
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

   dualfarkas = SCIPgetDualfarkasKnapsack(scip, cons);

   return (jdouble)dualfarkas;
}

/** returns the linear relaxation of the given knapsack constraint; may return NULL if no LP row was yet created; */
JNIEXPORT
jlong JNISCIPCONSKNAPSACK(getRowKnapsack)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_ROW*  row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   row = SCIPgetRowKnapsack(scip, cons);

   return (jlong)(size_t)row;
}
