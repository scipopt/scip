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

#include "JniScipConsQuadratic.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include <string.h>

/** creates the handler for or constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSQUADRATIC(includeConshdlrQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrQuadratic(scip) );
}

/** Creates and captures a quadratic constraint.
 *
 *  The constraint should be given in the form
 *  \f[
 *  \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m a_j y_jz_j \leq u,
 *  \f]
 *  where \f$x_i = y_j = z_k\f$ is possible.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSQUADRATIC(createConsQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   int                   jnlinvars,          /**< number of linear terms (n) */
   jlongArray            jlinvars,           /**< array with variables in linear part (x_i) */
   jdoubleArray          jlincoefs,          /**< array with coefficients of variables in linear part (b_i) */
   jint                  jnquadterms,        /**< number of quadratic terms (m) */
   jlongArray            jquadvars1,         /**< array with first variables in quadratic terms (y_j) */
   jlongArray            jquadvars2,         /**< array with second variables in quadratic terms (z_j) */
   jdoubleArray          jquadcoefs,         /**< array with coefficients of quadratic terms (a_j) */
   jdouble               jlhs,               /**< left hand side of quadratic equation (ell) */
   jdouble               jrhs,               /**< right hand side of quadratic equation (u) */
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
   SCIPerrorMessage("method createConsQuadratic is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}


/** creates and captures a quadratic constraint with all its
 *  flags set to their default values.
 *
 *  The constraint should be given in the form
 *  \f[
 *  \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m a_j y_jz_j \leq u,
 *  \f]
 *  where \f$x_i = y_j = z_k\f$ is possible.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSQUADRATIC(createConsBasicQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnlinvars,           /**< number of linear terms (n) */
   jlongArray            jlinvars,            /**< array with variables in linear part (x_i) */
   jdoubleArray          jlincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   jint                  jnquadterms,         /**< number of quadratic terms (m) */
   jlongArray            jquadvars1,          /**< array with first variables in quadratic terms (y_j) */
   jlongArray            jquadvars2,          /**< array with second variables in quadratic terms (z_j) */
   jdoubleArray          jquadcoefs,          /**< array with coefficients of quadratic terms (a_j) */
   jdouble               jlhs,                /**< left hand side of quadratic equation (ell) */
   jdouble               jrhs                 /**< right hand side of quadratic equation (u) */
   )
{
   SCIPerrorMessage("method createConsBasicQuadratic is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** creates and captures a quadratic constraint in its most basic version, i.e.,
 *  all constraint flags are set to their default values.
 *
 * The constraint should be given in the form
 * \f[
 * \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m (a_j y_j^2 + b_j y_j) + \sum_{k=1}^p c_kv_kw_k \leq u.
 * \f]
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSQUADRATIC(createConsQuadratic2)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnlinvars,          /**< number of linear terms (n) */
   jlongArray            jlinvars,           /**< array with variables in linear part (x_i) */
   jdoubleArray          jlincoefs,          /**< array with coefficients of variables in linear part (b_i) */
   jint                  jnquadvarterms,     /**< number of quadratic terms (m) */
   jlong                 jquadvarterms,      /**< quadratic variable terms */
   jint                  jnbilinterms,       /**< number of bilinear terms (p) */
   jlong                 jbilinterms,        /**< bilinear terms */
   jdouble               jlhs,               /**< constraint left hand side (ell) */
   jdouble               jrhs,               /**< constraint right hand side (u) */
   jboolean              jinitial,           /**< should the LP relaxation of constraint be in the initial LP? */
   jboolean              jseparate,          /**< should the constraint be separated during LP processing? */
   jboolean              jenforce,           /**< should the constraint be enforced during node processing? */
   jboolean              jcheck,             /**< should the constraint be checked for feasibility? */
   jboolean              jpropagate,         /**< should the constraint be propagated during node processing? */
   jboolean              jlocal,             /**< is constraint only valid locally? */
   jboolean              jmodifiable,        /**< is constraint modifiable (subject to column generation)? */
   jboolean              jdynamic,           /**< is constraint dynamic? */
   jboolean              jremovable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   SCIPerrorMessage("method createConsQuadratic2 is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** creates and captures a quadratic constraint in its most basic version, i.e.,
 *  all constraint flags are set to their default values.
 *
 * The constraint should be given in the form
 * \f[
 * \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m (a_j y_j^2 + b_j y_j) + \sum_{k=1}^p c_kv_kw_k \leq u.
 * \f]
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSQUADRATIC(createConsBasicQuadratic2)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnlinvars,          /**< number of linear terms (n) */
   jlongArray            jlinvars,           /**< array with variables in linear part (x_i) */
   jdoubleArray          jlincoefs,          /**< array with coefficients of variables in linear part (b_i) */
   jint                  jnquadvarterms,     /**< number of quadratic terms (m) */
   jlong                 jquadvarterms,      /**< quadratic variable terms */
   jint                  jnbilinterms,       /**< number of bilinear terms (p) */
   jlong                 jbilinterms,        /**< bilinear terms */
   jdouble               jlhs,               /**< constraint left hand side (ell) */
   jdouble               jrhs                /**< constraint right hand side (u) */
   )
{
   SCIPerrorMessage("method createConsBasicQuadratic2 is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** Adds a constant to the constraint function, that is, subtracts a constant from both sides */
JNIEXPORT
void JNISCIPCONSQUADRATIC(addConstantQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jdouble               jconstant           /**< constant to subtract from both sides */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   SCIPaddConstantQuadratic(scip, cons, (SCIP_Real)jconstant);
}

/** Adds a linear variable with coefficient to a quadratic constraint.
 */
JNIEXPORT
void JNISCIPCONSQUADRATIC(addLinearVarQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jlong                 jvar,               /**< variable */
   jdouble               jcoef               /**< coefficient of variable */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;
   var = (SCIP_VAR*) (size_t) jvar;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, var, (SCIP_Real)jcoef) );
}

/** Adds a quadratic variable with linear and square coefficient to a quadratic constraint.
 */
JNIEXPORT
void JNISCIPCONSQUADRATIC(addQuadVarQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jlong                 jvar,               /**< variable */
   jdouble               jlincoef,           /**< linear coefficient of variable */
   jdouble               jsqrcoef            /**< square coefficient of variable */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;
   var = (SCIP_VAR*) (size_t) jvar;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddQuadVarQuadratic(scip, cons, var, (SCIP_Real)jlincoef, (SCIP_Real)jsqrcoef) );
}

/** Adds a linear coefficient for a quadratic variable.
 *
 * Variable will be added with square coefficient 0.0 if not existing yet.
 */
JNIEXPORT
void JNISCIPCONSQUADRATIC(addQuadVarLinearCoefQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jlong                 jvar,               /**< variable */
   jdouble               jcoef               /**< value to add to linear coefficient of variable */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;
   var = (SCIP_VAR*) (size_t) jvar;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddQuadVarLinearCoefQuadratic(scip, cons, var, (SCIP_Real)jcoef) );
}

/** Adds a square coefficient for a quadratic variable.
 *
 * Variable will be added with linear coefficient 0.0 if not existing yet.
 */
JNIEXPORT
void JNISCIPCONSQUADRATIC(addSquareCoefQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jlong                 jvar,               /**< variable */
   jdouble               jcoef               /**< value to add to square coefficient of variable */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;
   var = (SCIP_VAR*) (size_t) jvar;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, var, (SCIP_Real)jcoef) );
}

/** Adds a bilinear term to a quadratic constraint.
 *
 * Variables will be added with linear and square coefficient 0.0 if not existing yet.
 * If variables are equal, only the square coefficient of the variable is updated.
 */
JNIEXPORT
void JNISCIPCONSQUADRATIC(addBilinTermQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jlong                 jvar1,              /**< first variable */
   jlong                 jvar2,              /**< second variable */
   jdouble               jcoef               /**< coefficient of bilinear term */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR* var1;
   SCIP_VAR* var2;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;
   var1 = (SCIP_VAR*) (size_t) jvar1;
   var2 = (SCIP_VAR*) (size_t) jvar2;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var1 != NULL);
   assert(var2 != NULL);


   JNISCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, var1, var2, (SCIP_Real)jcoef) );
}

/** Gets the quadratic constraint as a nonlinear row representation.
 */
JNIEXPORT
jlong JNISCIPCONSQUADRATIC(getNlRowQuadratic)(
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

   JNISCIP_CALL( SCIPgetNlRowQuadratic(scip, cons, &nlrow) );

   return (jlong) (size_t) nlrow;
}

/** Gets the number of variables in the linear term of a quadratic constraint.
 */
JNIEXPORT
jint JNISCIPCONSQUADRATIC(getNLinearVarsQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   return (jint) SCIPgetNLinearVarsQuadratic(scip, cons);
}

/** Gets the variables in the linear part of a quadratic constraint.
 *  Length is given by SCIPgetNLinearVarsQuadratic.
 */
JNIEXPORT
jlongArray JNISCIPCONSQUADRATIC(getLinearVarsQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   int size;
   jlongArray jvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   size = SCIPgetNLinearVarsQuadratic(scip, cons);

   jvars = (*env)->NewLongArray(env, size);

   if (jvars == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   vars = SCIPgetLinearVarsQuadratic(scip, cons);

   (*env)->SetLongArrayRegion(env, jvars, 0, size, (jlong*)(*vars));

   return jvars;
}

/** Gets the coefficients in the linear part of a quadratic constraint.
 *  Length is given by SCIPgetNLinearVarsQuadratic.
 */
JNIEXPORT
jdoubleArray JNISCIPCONSQUADRATIC(getCoefsLinearVarsQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Real* vals;
   int size;
   jdoubleArray jvals;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   size = SCIPgetNLinearVarsQuadratic(scip, cons);

   jvals = (*env)->NewDoubleArray(env, size);

   if (jvals == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   vals = SCIPgetCoefsLinearVarsQuadratic(scip, cons);

   (*env)->SetDoubleArrayRegion(env, jvals, 0, size, (jdouble*)vals);

   return jvals;
}

/** Gets the bilinear terms of a quadratic constraint.
 *  Length is given by SCIPgetNBilinTermQuadratic.
 */
JNIEXPORT
jlong JNISCIPCONSQUADRATIC(getBilinTermsQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   return (jlong) (size_t) SCIPgetBilinTermsQuadratic(scip, cons);
}

/** Gets the left hand side of a quadratic constraint.
 */
JNIEXPORT
jdouble JNISCIPCONSQUADRATIC(getLhsQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   return (jdouble) SCIPgetLhsQuadratic(scip, cons);
}

/** Gets the right hand side of a quadratic constraint.
 */
JNIEXPORT
jdouble JNISCIPCONSQUADRATIC(getRhsQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   return (jdouble) SCIPgetRhsQuadratic(scip, cons);
}

/** Check the quadratic function of a quadratic constraint for its semi-definiteness, if not done yet.
 */
JNIEXPORT
void JNISCIPCONSQUADRATIC(checkCurvatureQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   JNISCIP_CALL( SCIPcheckCurvatureQuadratic(scip, cons) );
}

/** Indicates whether the quadratic function of a quadratic constraint is (known to be) convex.
 */
JNIEXPORT
jboolean JNISCIPCONSQUADRATIC(isConvexQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   return (jboolean) SCIPisConvexQuadratic(scip, cons);
}

/** Indicates whether the quadratic function of a quadratic constraint is (known to be) concave.
 */
JNIEXPORT
jboolean JNISCIPCONSQUADRATIC(isConcaveQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   return (jboolean) SCIPisConcaveQuadratic(scip, cons);
}

/** Computes the violation of a constraint by a solution */
JNIEXPORT
jdouble JNISCIPCONSQUADRATIC(getViolationQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jlong                 jsol                /**< solution which violation to calculate, or NULL for LP solution */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_Real violation;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;
   sol = (SCIP_SOL*) (size_t) jsol;

   assert(scip != NULL);
   assert(cons != NULL);

   JNISCIP_CALL( SCIPgetViolationQuadratic(scip, cons, sol, &violation) );

   return (jdouble) violation;
}

/** Indicates whether the quadratic constraint is local w.r.t. the current local bounds.
 *
 * That is, checks whether each variable with a square term is fixed and for each bilinear term at least one variable is fixed.
 */
JNIEXPORT
jboolean JNISCIPCONSQUADRATIC(isLinearLocalQuadratic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   return (jboolean) SCIPisLinearLocalQuadratic(scip, cons);
}

/** Adds the constraint to an NLPI problem. */
JNIEXPORT
void JNISCIPCONSQUADRATIC(addToNlpiProblemQuadratic)(
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

   JNISCIP_CALL( SCIPaddToNlpiProblemQuadratic(scip, cons, nlpi, nlpiprob, scipvar2nlpivar, (SCIP_Bool)jnames) );
}
