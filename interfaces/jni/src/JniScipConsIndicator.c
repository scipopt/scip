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

/**@file   JniScipConsIndicator.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP constraint indicator callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipConsIndicator.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_indicator.h"

#include <string.h>

/** creates the handler for indicator constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSINDICATOR(includeConshdlrIndicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrIndicator(scip) );
}


/** creates and captures an indicator constraint
 *
 *  @note @a binvar is checked to be binary only later. This enables a change of the type in
 *  procedures reading an instance.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method {@link releaseCons()}
 */
JNIEXPORT
jlong JNISCIPCONSINDICATOR(createConsIndicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jlong                 jbinvar,            /**< binary indicator variable (or NULL) */
   jint                  nvars,              /**< number of variables in the constraint */
   jlongArray            jvars,              /**< array with variables of constraint entries */
   jdoubleArray          jvals,              /**< values of variables in inequality (or NULL) */
   jdouble               rhs,                /**< rhs of the inequality */
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
   SCIP_VAR* binvar;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, NULL);
   if( name == NULL )
      SCIPABORT();

   /* convert JNI pointer into C pointer */
   binvar = (SCIP_VAR*) (size_t) jbinvar;

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, nvars, (jlong*)(*vars));
   (*env)->GetDoubleArrayRegion(env, jvals, 0, nvars, (jdouble*)vals);

   JNISCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name, binvar, (int)nvars, vars, vals, (SCIP_Real)rhs,
         (SCIP_Bool) jinitial, (SCIP_Bool) jseparate, (SCIP_Bool) jenforce, (SCIP_Bool) jcheck, (SCIP_Bool) jpropagate,
         (SCIP_Bool) jlocal, (SCIP_Bool) jdynamic, (SCIP_Bool) jremovable, (SCIP_Bool) jstickingatnode) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** creates and captures an indicator constraint with given linear constraint and slack variable
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsIndicator(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @note @a binvar is checked to be binary only later. This enables a change of the type in
 *  procedures reading an instance.
 *
 *  @note we assume that @a slackvar actually appears in @a lincons and we also assume that it takes
 *  the role of a slack variable!
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 *
 *  @see SCIPcreateConsIndicatorLinCons() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSINDICATOR(createConsBasicIndicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jlong                 jbinvar,            /**< binary indicator variable (or NULL) */
   jint                  nvars,              /**< number of variables in the inequality */
   jlongArray            jvars,              /**< array with variables of inequality (or NULL) */
   jdoubleArray          jvals,              /**< values of variables in inequality (or NULL) */
   jdouble               rhs                 /**< rhs of the inequality */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   const char* name;
   SCIP_VAR* binvar;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, NULL);
   if( name == NULL )
      SCIPABORT();

   /* convert JNI pointer into C pointer */
   binvar = (SCIP_VAR*) (size_t) jbinvar;

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, nvars, (jlong*)(*vars));
   (*env)->GetDoubleArrayRegion(env, jvals, 0, nvars, (jdouble*)vals);

   JNISCIP_CALL( SCIPcreateConsBasicIndicator(scip, &cons, name, binvar, (int)nvars, vars, vals, (SCIP_Real)rhs) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}



/** creates and captures an indicator constraint with given linear constraint and slack variable
 *
 *  Note: @a binvar is checked to be binary only later. This enables a change of the type in
 *  procedures reading an instance.
 */
JNIEXPORT
jlong JNISCIPCONSINDICATOR(createConsIndicatorLinCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jlong                 jbinvar,            /**< binary indicator variable */
   jlong                 jlincons,           /**< linear constraint */
   jlong                 jslackvar,          /**< slack variable */
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
   SCIP_VAR* binvar;
   SCIP_CONS* lincons;
   SCIP_VAR* slackvar;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   binvar = (SCIP_VAR*) (size_t) jbinvar;

   lincons = (SCIP_CONS*) (size_t) jlincons;
   assert(lincons != NULL);

   slackvar = (SCIP_VAR*) (size_t) jslackvar;
   assert(slackvar != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, NULL);
   if( name == NULL )
      SCIPABORT();

   JNISCIP_CALL( SCIPcreateConsIndicatorLinCons(scip, &cons, name, binvar, lincons, slackvar,
         (SCIP_Bool) jinitial, (SCIP_Bool) jseparate, (SCIP_Bool) jenforce, (SCIP_Bool) jcheck, (SCIP_Bool) jpropagate,
         (SCIP_Bool) jlocal, (SCIP_Bool) jdynamic, (SCIP_Bool) jremovable, (SCIP_Bool) jstickingatnode) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** adds variable to the inequality of the indicator constraint */
JNIEXPORT
void JNISCIPCONSINDICATOR(addVarIndicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< indicator constraint */
   jlong                 jvar,               /**< variable to add to the inequality */
   jdouble               jval                /**< value of variable */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddVarIndicator(scip, cons, var, (SCIP_Real) jval) );
}

/** gets the linear constraint corresponding to the indicator constraint (may be NULL) */
JNIEXPORT
jlong JNISCIPCONSINDICATOR(getLinearConsIndicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< indicator constraint */
   )
{
   SCIP_CONS* cons;
   SCIP_CONS* lincons;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;

   lincons = SCIPgetLinearConsIndicator(cons);

   return (jlong) (size_t) lincons;

}

/** sets the linear constraint corresponding to the indicator constraint (may be NULL) */
JNIEXPORT
void JNISCIPCONSINDICATOR(setLinearConsIndicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< indicator constraint */
   jlong                 jlincons            /**< linear constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_CONS* lincons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   lincons = (SCIP_CONS*) (size_t) jlincons;

   JNISCIP_CALL( SCIPsetLinearConsIndicator(scip, cons, lincons) );
}

/** sets binary indicator variable for indicator constraint */
JNIEXPORT
void JNISCIPCONSINDICATOR(setBinaryVarIndicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< indicator constraint */
   jlong                 jbinvar             /**< binary variable to add to the inequality */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR* binvar;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   binvar = (SCIP_VAR*) (size_t) jbinvar;
   assert(binvar != NULL);

   JNISCIP_CALL( SCIPsetBinaryVarIndicator(scip, cons, binvar) );
}

/** gets binary variable corresponding to indicator constraint */
JNIEXPORT
jlong JNISCIPCONSINDICATOR(getBinaryVarIndicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< indicator constraint */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* binvar;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   binvar = SCIPgetBinaryVarIndicator(cons);

   return (jlong)(size_t) binvar;
}

/** gets slack variable corresponding to indicator constraint */
JNIEXPORT
jlong JNISCIPCONSINDICATOR(getSlackVarIndicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< indicator constraint */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* slackvar;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   slackvar = SCIPgetSlackVarIndicator(cons);

   return (jlong)(size_t) slackvar;
}

/** checks whether indicator constraint is violated w.r.t. sol */
JNIEXPORT
jboolean JNISCIPCONSINDICATOR(isViolatedIndicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< indicator constraint */
   jlong                 jsol                /**< solution, or NULL to use current node's solution */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_Bool isviolated;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   isviolated = SCIPisViolatedIndicator(scip, cons, sol);

   return (jboolean) isviolated;
}

/** Based on values of other variables, computes slack and binary variable to turn constraint feasible */
JNIEXPORT
jboolean JNISCIPCONSINDICATOR(makeIndicatorFeasible)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< indicator constraint */
   jlong                 jsol                /**< solution */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_Bool changed;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPmakeIndicatorFeasible(scip, cons, sol, &changed) );

   return (jboolean) changed;
}

/** Based on values of other variables, computes slack and binary variable to turn all constraints feasible */
JNIEXPORT
jboolean JNISCIPCONSINDICATOR(makeIndicatorsFeasible)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jconshdlr,          /**< indicator constraint handler */
   jlong                 jsol                /**< solution */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_SOL* sol;
   SCIP_Bool changed;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPmakeIndicatorsFeasible(scip, conshdlr, sol, &changed) );

   return (jboolean) changed;
}

/** adds additional linear constraint that is not connected by an indicator constraint, but can be used for separation */
JNIEXPORT
void JNISCIPCONSINDICATOR(addLinearConsIndicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jconshdlr,          /**< indicator constraint handler */
   jlong                 jlincons            /**< linear constraint */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS* lincons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   lincons = (SCIP_CONS*) (size_t) jlincons;
   assert(lincons != NULL);

   JNISCIP_CALL( SCIPaddLinearConsIndicator(scip, conshdlr, lincons) );
}

/** adds additional globally valid row that is not connected by an indicator constraint, but can be used for separation */
JNIEXPORT
void JNISCIPCONSINDICATOR(addRowIndicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jconshdlr,          /**< indicator constraint handler */
   jlong                 jrow                /**< row to add */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPaddRowIndicator(scip, conshdlr, row) );
}
