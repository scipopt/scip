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

#include "JniScipConsSoc.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_soc.h"

#include <string.h>

/** creates the handler for logic or constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSSOC(includeConshdlrSOC)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrSOC(scip) );
}

/** creates and captures a second order cone constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSSOC(createConsSOC)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnvars,             /**< number of variables on left hand side of constraint (n) */
   jlongArray            jvars,              /**< array with variables on left hand side (x_i) */
   jdoubleArray          jvals,              /**< array with coefficients of left hand side variables (alpha_i), or NULL if all 1.0 */
   jdoubleArray          joffsets,           /**< array with offsets of variables (beta_i), or NULL if all 0.0 */
   jdouble               jconstant,          /**< constant on left hand side (gamma) */
   jlong                 jrhsvar,            /**< variable on right hand side of constraint (x_{n+1}) */
   jdouble               jrhscoeff,          /**< coefficient of variable on right hand side (alpha_{n+1}) */
   jdouble               jrhsoffset,         /**< offset of variable on right hand side (beta_{n+1}) */
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
                                              *   are separated as constraints. */
   jboolean              jremovable          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   const char* name;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real* offsets;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, NULL);
   if( name == NULL )
      SCIPABORT();

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)jnvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)jnvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &offsets, (int)jnvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)jnvars, (jlong*)(*vars));
   (*env)->GetDoubleArrayRegion(env, jvals, 0, (int)jnvars, (jdouble*)vals);
   (*env)->GetDoubleArrayRegion(env, jvals, 0, (int)jnvars, (jdouble*)offsets);

   JNISCIP_CALL( SCIPcreateConsSOC(scip, &cons, name, (int)jnvars, vars, vals, offsets, (SCIP_Real)jconstant, (SCIP_VAR*)(size_t)jrhsvar, (SCIP_Real)jrhscoeff, (SCIP_Real)jrhsoffset,
         (SCIP_Bool) jinitial, (SCIP_Bool) jseparate, (SCIP_Bool) jenforce, (SCIP_Bool) jcheck, (SCIP_Bool) jpropagate,
         (SCIP_Bool) jlocal, (SCIP_Bool) jmodifiable, (SCIP_Bool) jdynamic, (SCIP_Bool) jremovable) );

   SCIPfreeBufferArray(scip, &offsets);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** creates and captures a second order cone constraint with all its constraint flags
 *  set to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSSOC(createConsBasicSOC)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   int                   jnvars,             /**< number of variables in the constraint */
   jlongArray            jvars,              /**< array with variables of constraint entries */
   jdoubleArray          jvals,              /**< array with coefficients of left hand side variables (alpha_i), or NULL if all 1.0 */
   jdoubleArray          joffsets,           /**< array with offsets of variables (beta_i), or NULL if all 0.0 */
   jdouble               jconstant,          /**< constant on left hand side (gamma) */
   jlong                 jrhsvar,            /**< variable on right hand side of constraint (x_{n+1}) */
   jdouble               jrhscoeff,          /**< coefficient of variable on right hand side (alpha_{n+1}) */
   jdouble               jrhsoffset          /**< offset of variable on right hand side (beta_{n+1}) */
   )
{
   SCIPerrorMessage("method createConsBasicSOC is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** Gets the SOC constraint as a nonlinear row representation.
 */
JNIEXPORT
jlong JNISCIPCONSSOC(getNlRowSOC)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_NLROW* nlrow;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   JNISCIP_CALL( SCIPgetNlRowSOC(scip, cons, &nlrow) );

   return (jlong) (size_t) nlrow;
}

/** Gets the number of variables on the left hand side of a SOC constraint.
 */
JNIEXPORT
jint JNISCIPCONSSOC(getNLhsVarsSOC)(
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

   return (jint) SCIPgetNLhsVarsSOC(scip, cons);
}

/** Gets the variables on the left hand side of a SOC constraint.
 */
JNIEXPORT
jlongArray JNISCIPCONSSOC(getLhsCoefsSOC)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Real* vars;
   int size;
   jdoubleArray jvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   size = SCIPgetNLhsVarsSOC(scip, cons);

   jvars = (*env)->NewDoubleArray(env, size);

   if (jvars == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   vars = SCIPgetLhsCoefsSOC(scip, cons);

   (*env)->SetDoubleArrayRegion(env, jvars, 0, size, (jdouble*)(vars));

   return jvars;
}

/** Gets the offsets of the variables on the left hand side of a SOC constraint, or NULL if all are equal to 0.0.
 */
JNIEXPORT
jdoubleArray JNISCIPCONSSOC(getLhsOffsetsSOC)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Real* vars;
   int size;
   jdoubleArray jvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   size = SCIPgetNLhsVarsSOC(scip, cons);

   jvars = (*env)->NewDoubleArray(env, size);

   if (jvars == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   vars = SCIPgetLhsOffsetsSOC(scip, cons);

   (*env)->SetDoubleArrayRegion(env, jvars, 0, size, (jdouble*)(vars));

   return jvars;
}

/** Gets the constant on the left hand side of a SOC constraint.
 */
JNIEXPORT
jdouble JNISCIPCONSSOC(getLhsConstantSOC)(
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

   return (jdouble) SCIPgetLhsConstantSOC(scip, cons);
}

/** Gets the variable on the right hand side of a SOC constraint.
 */
JNIEXPORT
jlong JNISCIPCONSSOC(getRhsVarSOC)(
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

   return (jlong) (size_t) SCIPgetRhsVarSOC(scip, cons);
}

/** Gets the coefficient of the variable on the right hand side of a SOC constraint.
 */
JNIEXPORT
jdouble JNISCIPCONSSOC(getRhsCoefSOC)(
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

   return (jdouble) SCIPgetRhsCoefSOC(scip, cons);
}

/** Gets the offset of the variables on the right hand side of a SOC constraint.
 */
JNIEXPORT
jdouble JNISCIPCONSSOC(getRhsOffsetSOC)(
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

   return (jdouble) SCIPgetRhsOffsetSOC(scip, cons);
}

/** Adds the constraint to an NLPI problem.
 * Uses nonconvex formulation as quadratic function.
 */
JNIEXPORT
void JNISCIPCONSSOC(addToNlpiProblemSOC)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jlong                 jnlpi,              /**< interface to NLP solver */
   jlong                 jnlpiprob,          /**< NLPI problem where to add constraint */
   jlong                 jscipvar2nlpivar,   /**< mapping from SCIP variables to variable indices in NLPI */
   jboolean              jnames              /**< whether to pass constraint names to NLPI */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_NLPI* nlpi;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_HASHMAP* scipvar2nlpivar;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;
   nlpi = (SCIP_NLPI*) (size_t) jnlpi;
   nlpiprob = (SCIP_NLPIPROBLEM*) (size_t) jnlpiprob;
   scipvar2nlpivar = (SCIP_HASHMAP*) (size_t) jscipvar2nlpivar;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(scipvar2nlpivar != NULL);

   JNISCIP_CALL( SCIPaddToNlpiProblemSOC(scip, cons, nlpi, nlpiprob, scipvar2nlpivar, (SCIP_Bool)jnames) );
}
