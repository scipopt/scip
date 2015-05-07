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

/**@file   JniScipConsLogicor.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP constraint logic or callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipConsLogicor.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_logicor.h"

#include <string.h>

/** creates the handler for logic or constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSLOGICOR(includeConshdlrLogicor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrLogicor(scip) );
}

/** creates and captures a logic or constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method {@link releaseCons()}
 */
JNIEXPORT
jlong JNISCIPCONSLOGICOR(createConsLogicor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnvars,             /**< number of variables in the constraint */
   jlongArray            jvars,               /**< array with variables of constraint entries */
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

   /* create logicor constraint with zero variables */
   JNISCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name, 0, NULL,
         (SCIP_Bool) jinitial, (SCIP_Bool) jseparate, (SCIP_Bool) jenforce, (SCIP_Bool) jcheck, (SCIP_Bool) jpropagate,
         (SCIP_Bool) jlocal, (SCIP_Bool) jmodifiable, (SCIP_Bool) jdynamic, (SCIP_Bool) jremovable, (SCIP_Bool) jstickingatnode) );

   /* convert JNI integer to C integer */
   nvars = (int)jnvars;

   /* add coefficients */
   if( nvars > 0 )
   {
      jlong* vars;
      int v;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

      (*env)->GetLongArrayRegion(env, jvars, 0, nvars, vars);

      for( v = 0; v < nvars; ++v )
      {
         JNISCIP_CALL( SCIPaddCoefLogicor(scip, cons, (SCIP_VAR*)(size_t)vars[v]) );
      }

      SCIPfreeBufferArray(scip, &vars);
   }

   return (jlong)(size_t)cons;
}

/** creates and captures a logicor constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsLogicor(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsLogicor() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSLOGICOR(createConsBasicLogicor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  nvars,              /**< number of variables in the constraint */
   jlongArray            jvars               /**< array with variables of constraint entries */
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

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, nvars, (jlong*)(*vars));

   /* create logicor constraint with zero variables */
   JNISCIP_CALL( SCIPcreateConsBasicLogicor(scip, &cons, name, (int)nvars, vars ) );

   SCIPfreeBufferArray(scip, &vars);

   return (jlong)(size_t)cons;
}

/** adds coefficient in logic or constraint */
JNIEXPORT
void JNISCIPCONSLOGICOR(addCoefLogicor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< logicor constraint */
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

   JNISCIP_CALL( SCIPaddCoefLogicor(scip, cons, var) );
}

/** gets number of variables in logic or constraint */
JNIEXPORT
jint JNISCIPCONSLOGICOR(getNVarsLogicor)(
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

   nvars = SCIPgetNVarsLogicor(scip, cons);

   return (jint)nvars;
}

/** gets array of variables in logic or constraint */
JNIEXPORT
jlongArray JNISCIPCONSLOGICOR(getVarsLogicor)(
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

   nvars = SCIPgetNVarsLogicor(scip, cons);
   vars = SCIPgetVarsLogicor(scip, cons);

   /* create jdouble Array */
   jvars = (*env)->NewLongArray(env, nvars);

   /* fill long array with SCIP variable pointers */
   (*env)->SetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);

   SCIPfreeBufferArray(scip, &vars);

   return jvars;
}

/** gets the dual solution of the logic or constraint in the current LP */
JNIEXPORT
jdouble JNISCIPCONSLOGICOR(getDualsolLogicor)(
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

   dualsol = SCIPgetDualsolLogicor(scip, cons);

   return (jdouble)dualsol;
}

/** gets the dual farkas value of the logic or constraint in the current infeasible LP */
JNIEXPORT
jdouble JNISCIPCONSLOGICOR(getDualfarkasLogicor)(
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

   dualfarkas = SCIPgetDualsolLogicor(scip, cons);

   return (jdouble)dualfarkas;
}

/** returns the linear relaxation of the given logic or constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
JNIEXPORT
jlong JNISCIPCONSLOGICOR(getRowLogicor)(
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

   row = SCIPgetRowLogicor(scip, cons);

   return (jlong)(size_t)row;

}
