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

/**@file   JniScipSol.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP solution callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipConsNonlinear.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include <string.h>


/** creates the handler for nonlinear constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSNONLINEAR(includeConshdlrNonlinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrNonlinear(scip) );
}

/** adds a linear variable with coefficient to a nonlinear constraint */
JNIEXPORT
void JNISCIPCONSNONLINEAR(addLinearVarNonlinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jlong                 jvar,               /**< variable */
   jdouble               coef                /**< coefficient of variable */
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

   JNISCIP_CALL( SCIPaddLinearVarNonlinear(scip, cons, var, (SCIP_Real)coef) );
}

/** gets the nonlinear constraint as a nonlinear row representation */
JNIEXPORT
jlong JNISCIPCONSNONLINEAR(getNlRowNonlinear)(
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
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPgetNlRowNonlinear(scip, cons, &nlrow) );

   return (jlong) (size_t) nlrow;
}

/** gets the number of variables in the linear term of a nonlinear constraint */
JNIEXPORT
jint JNISCIPCONSNONLINEAR(getNLinearVarsNonlinear)(
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
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   return (jint) SCIPgetNLinearVarsNonlinear(scip, cons);
}

/** gets the variables in the linear part of a nonlinear constraint */
JNIEXPORT
jlongArray JNISCIPCONSNONLINEAR(getLinearVarsNonlinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   jlongArray jvars;
   int size;
   SCIP_VAR** vars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   size = SCIPgetNLinearVarsNonlinear(scip, cons);

   jvars = (*env)->NewLongArray(env, size);

   if (jvars == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   vars = SCIPgetLinearVarsNonlinear(scip, cons);

   (*env)->SetLongArrayRegion(env, jvars, 0, size, (jlong*)vars);

   return jvars;
}

/** gets the coefficients in the linear part of a nonlinear constraint */
JNIEXPORT
jdoubleArray JNISCIPCONSNONLINEAR(getLinearCoefsNonlinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   jdoubleArray jcoefs;
   int size;
   SCIP_Real* coefs;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   cons = (SCIP_CONS*) (size_t) jcons;

   assert(scip != NULL);
   assert(cons != NULL);

   size = SCIPgetNLinearVarsNonlinear(scip, cons);

   jcoefs = (*env)->NewDoubleArray(env, size);

   if (jcoefs == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   coefs = SCIPgetLinearCoefsNonlinear(scip, cons);

   (*env)->SetDoubleArrayRegion(env, jcoefs, 0, size, (jdouble*)coefs);

   return jcoefs;
}

/** gets the number of expression trees of a nonlinear constraint */
JNIEXPORT
jint JNISCIPCONSNONLINEAR(getNExprtreesNonlinear)(
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
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   return (jint) SCIPgetNExprtreesNonlinear(scip, cons);
}

/** gets the expression trees of a nonlinear constraint */
JNIEXPORT
jdoubleArray JNISCIPCONSNONLINEAR(getExprtreeCoefsNonlinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   jdoubleArray jexpr;
   int size;
   SCIP_Real* expr;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   size = SCIPgetNExprtreesNonlinear(scip, cons);

   jexpr = (*env)->NewDoubleArray(env, size);

   if (jexpr == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   expr = SCIPgetExprtreeCoefsNonlinear(scip, cons);

   (*env)->SetDoubleArrayRegion(env, jexpr, 0, size, (jdouble*)expr);

   return jexpr;
}
/** gets the left hand side of a nonlinear constraint */
JNIEXPORT
jdouble JNISCIPCONSNONLINEAR(getLhsNonlinear)(
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
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   return (jdouble) SCIPgetLhsNonlinear(scip, cons);
}

/** gets the right hand side of a nonlinear constraint */
JNIEXPORT
jdouble JNISCIPCONSNONLINEAR(getRhsNonlinear)(
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
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   return (jdouble) SCIPgetRhsNonlinear(scip, cons);
}

/** check the function of a nonlinear constraint for convexity/concavity, if not done yet */
JNIEXPORT
void JNISCIPCONSNONLINEAR(checkCurvatureNonlinear)(
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
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPcheckCurvatureNonlinear(scip, cons) );
}

/** computes the violation of a nonlinear constraint by a solution */
JNIEXPORT
jdouble JNISCIPCONSNONLINEAR(getViolationNonlinear)(
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
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPgetViolationNonlinear(scip, cons, sol, &violation) );

   return (jdouble) violation;
}
