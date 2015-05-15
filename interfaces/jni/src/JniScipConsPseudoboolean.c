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

#include "JniScipConsPseudoboolean.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include <string.h>

/** creates the handler for or constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSPSEUDOBOOLEAN(includeConshdlrPseudoboolean)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrPseudoboolean(scip) );
}

/** creates and captures a pseudoboolean constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSPSEUDOBOOLEAN(createConsPseudoboolean)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jlongArray            jlinvars,            /**< variables of the linear part, or NULL */
   jint                  jnlinvars,           /**< number of variables of the linear part */
   jdoubleArray          jlinvals,            /**< coefficients of linear part, or NULL */
   jobjectArray          jterms,              /**< nonlinear terms of variables, or NULL */
   jint                  jnterms,             /**< number of terms of variables of nonlinear term */
   jintArray             jntermvars,          /**< number of variables in nonlinear terms, or NULL */
   jdoubleArray          jtermvals,           /**< coefficients of nonlinear parts, or NULL */
   jlong                 jindvar,             /**< indicator variable if it's a soft constraint, or NULL */
   jdouble               jweight,             /**< weight of the soft constraint, if it is one */
   jboolean              jissoftcons,         /**< is this a soft constraint */
   jlong                 jintvar,             /**< a artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   jdouble               jlhs,                /**< left hand side of constraint */
   jdouble               jrhs,                /**< right hand side of constraint */
   jboolean              jinitial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   jboolean              jseparate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   jboolean              jenforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   jboolean              jcheck,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   jboolean              jpropagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   jboolean              jlocal,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   jboolean              jmodifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   jboolean              jdynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   jboolean              jremovable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   jboolean              jstickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIPerrorMessage("method createConsPseudoboolean is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** creates and captures a pseudoboolean constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSPSEUDOBOOLEAN(createConsBasicPseudoboolean)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jlongArray            jlinvars,           /**< variables of the linear part, or NULL */
   jint                  jnlinvars,          /**< number of variables of the linear part */
   jdoubleArray          jlinvals,           /**< coefficients of linear part, or NULL */
   jobjectArray          jterms,             /**< nonlinear terms of variables, or NULL */
   jint                  jnterms,            /**< number of terms of variables of nonlinear term */
   jintArray             jntermvars,         /**< number of variables in nonlinear terms, or NULL */
   jdoubleArray          jtermvals,          /**< coefficients of nonlinear parts, or NULL */
   jlong                 jindvar,            /**< indicator variable if it's a soft constraint, or NULL */
   jdouble               jweight,            /**< weight of the soft constraint, if it is one */
   jboolean              jissoftcons,        /**< is this a soft constraint */
   jlong                 jintvar,            /**< a artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   jdouble               jlhs,               /**< left hand side of constraint */
   jdouble               jrhs                /**< right hand side of constraint */
   )
{
   SCIPerrorMessage("method createConsBasicPseudoboolean is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** adds a variable to the pseudo boolean constraint (if it is not zero) */
JNIEXPORT
void JNISCIPCONSPSEUDOBOOLEAN(addCoefPseudoboolean)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint data */
   jlong                 jvar,               /**< variable of constraint entry */
   jdouble               val                 /**< coefficient of constraint entry */
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


   JNISCIP_CALL( SCIPaddCoefPseudoboolean(scip, cons, var, (SCIP_Real)val) );
}

/** adds nonlinear term to pseudo boolean constraint (if it is not zero) */
JNIEXPORT
void JNISCIPCONSPSEUDOBOOLEAN(addTermPseudoboolean)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint data */
   jlongArray            jvars,              /**< variables of the nonlinear term */
   jint                  nvars,              /**< number of variables of the nonlinear term */
   jdouble               val                 /**< coefficient of constraint entry */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR** vars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)(*vars));

   JNISCIP_CALL( SCIPaddTermPseudoboolean(scip, cons, vars, (int)nvars, (SCIP_Real)val) );

   SCIPfreeBufferArray(scip, &vars);
}

/** gets left hand side of pseudoboolean constraint */
JNIEXPORT
jdouble JNISCIPCONSPSEUDOBOOLEAN(getLhsPseudoboolean)(
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

   return (jdouble) SCIPgetLhsPseudoboolean(scip, cons);
}

/** gets lefright hand side of pseudoboolean constraint */
JNIEXPORT
jdouble JNISCIPCONSPSEUDOBOOLEAN(getRhsPseudoboolean)(
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

   return (jdouble) SCIPgetRhsPseudoboolean(scip, cons);
}

/** gets indicator variable of pseudoboolean constraint, or NULL if there is no */
JNIEXPORT
jlong JNISCIPCONSPSEUDOBOOLEAN(getIndVarPseudoboolean)(
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
   cons= (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   return (jlong) (size_t) SCIPgetIndVarPseudoboolean(scip, cons);
}

/** changes left hand side of pseudoboolean constraint */
JNIEXPORT
void JNISCIPCONSPSEUDOBOOLEAN(chgLhsPseudoboolean)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint data */
   jdouble               lhs                 /**< new left hand side */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);


   JNISCIP_CALL( SCIPchgLhsPseudoboolean(scip, cons, (SCIP_Real)lhs) );
}

/** changes left hand side of pseudoboolean constraint */
JNIEXPORT
void JNISCIPCONSPSEUDOBOOLEAN(chgRhsPseudoboolean)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint data */
   jdouble               rhs                 /**< new left hand side */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);


   JNISCIP_CALL( SCIPchgRhsPseudoboolean(scip, cons, (SCIP_Real)rhs) );
}
