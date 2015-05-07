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

/**@file   JniScipConsLinear.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP linear constraint callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipConsLinear.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"

#include <string.h>


/*
 * constraint methods
 */

/**@name Constraint Methods */
/**@{ */

/** creates the handler for linear constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSLINEAR(includeConshdlrLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrLinear(scip) );
}

/** creates and captures a linear constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method {@link releaseCons()}
 */
JNIEXPORT
jlong JNISCIPCONSLINEAR(createConsLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnvars,             /**< number of nonzeros in the constraint */
   jlongArray            jvars,              /**< array with variables of constraint entries */
   jdoubleArray          jvals,              /**< array with coefficients of constraint entries */
   jdouble               jlhs,               /**< left hand side of constraint */
   jdouble               jrhs,               /**< right hand side of constraint */
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

   /* create linear constraint with zero variables */
   JNISCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 0, NULL, NULL, (SCIP_Real) jlhs, (SCIP_Real) jrhs,
         (SCIP_Bool) initial, (SCIP_Bool) separate, (SCIP_Bool) enforce, (SCIP_Bool) check, (SCIP_Bool) propagate,
         (SCIP_Bool) local, (SCIP_Bool) modifiable, (SCIP_Bool) dynamic, (SCIP_Bool) removable, (SCIP_Bool) stickingatnode) );

   /* convert JNI integer into integer */
   nvars = (int)jnvars;

   if( nvars > 0 )
   {
      jlong* vars;
      jdouble* vals;
      int v;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );

      (*env)->GetLongArrayRegion(env, jvars, 0, nvars, vars);
      (*env)->GetDoubleArrayRegion(env, jvals, 0, nvars, vals);

      for( v = 0; v < nvars; ++v )
      {
         JNISCIP_CALL( SCIPaddCoefLinear(scip, cons, (SCIP_VAR*)(size_t)vars[v], (SCIP_Real)vals[v]));
      }

      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &vars);
   }

   /* relase string object */
   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** creates and captures a linear constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsLinear(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsLinear() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSLINEAR(createConsBasicLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jint                  jnvars,             /**< number of nonzeros in the constraint */
   jlongArray            jvars,              /**< array with variables of constraint entries */
   jdoubleArray          jvals,              /**< array with coefficients of constraint entries */
   jdouble               jlhs,               /**< left hand side of constraint */
   jdouble               jrhs                /**< right hand side of constraint */
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

   /* create linear constraint with zero variables */
   JNISCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, (SCIP_Real) jlhs, (SCIP_Real) jrhs) );

   /* convert JNI integer into integer */
   nvars = (int)jnvars;

   if( nvars > 0 )
   {
      jlong* vars;
      jdouble* vals;
      int v;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );

      (*env)->GetLongArrayRegion(env, jvars, 0, nvars, vars);
      (*env)->GetDoubleArrayRegion(env, jvals, 0, nvars, vals);

      for( v = 0; v < nvars; ++v )
      {
         JNISCIP_CALL( SCIPaddCoefLinear(scip, cons, (SCIP_VAR*)(size_t)vars[v], (SCIP_Real)vals[v]));
      }

      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &vars);
   }

   /* relase string object */
   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** adds coefficient to linear constraint (if it is not zero) */
JNIEXPORT
void JNISCIPCONSLINEAR(addCoefLinear)(
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
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert( cons != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddCoefLinear(scip, cons, var, (SCIP_Real) val) );
}


/** gets left hand side of linear constraint */
JNIEXPORT
jdouble JNISCIPCONSLINEAR(getLhsLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   jdouble lhs;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert( cons != NULL);

   lhs = (jdouble) SCIPgetLhsLinear(scip, cons);

   return lhs;
}

/** gets right hand side of linear constraint */
JNIEXPORT
jdouble JNISCIPCONSLINEAR(getRhsLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   jdouble lhs;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert( cons != NULL);

   lhs = (jdouble) SCIPgetRhsLinear(scip, cons);

   return lhs;
}

/** changes left hand side of linear constraint */
JNIEXPORT
void JNISCIPCONSLINEAR(chgLhsLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint data */
   jdouble               val                 /**< new left hand side */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert( cons != NULL);

   JNISCIP_CALL( SCIPchgLhsLinear(scip, cons, (SCIP_Real) val) );
}

/** changes right hand side of linear constraint */
JNIEXPORT
void JNISCIPCONSLINEAR(chgRhsLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint data */
   jdouble               val                 /**< new left hand side */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert( cons != NULL);

   JNISCIP_CALL( SCIPchgRhsLinear(scip, cons, (SCIP_Real) val) );
}

/** gets the number of variables in the linear constraint */
JNIEXPORT
jint JNISCIPCONSLINEAR(getNVarsLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
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
   assert( cons != NULL);

   num = SCIPgetNVarsLinear(scip, cons);

   return (jint) num;
}

/** gets the array of variables in the linear constraint; the user must not modify this array! */
JNIEXPORT
jlongArray JNISCIPCONSLINEAR(getVarsLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   int nvars;

   jlongArray jvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert( cons != NULL);

   nvars = SCIPgetNVarsLinear(scip, cons);

   /* create jlongArray */
   jvars = (*env)->NewLongArray(env, nvars);

   if( jvars == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_VAR** vars;

      /* fill long array with SCIP variable pointers */
      vars = SCIPgetVarsLinear(scip, cons);
      (*env)->SetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);
   }

   return jvars;
}

/** gets the array of coefficient values in the linear constraint; the user must not modify this array! */
jdoubleArray JNISCIPCONSLINEAR(getValsLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   int nvars;

   jdoubleArray jvals;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert( cons != NULL);

   nvars = SCIPgetNVarsLinear(scip, cons);

   /* create jlongArray */
   jvals = (*env)->NewDoubleArray(env, nvars);

   if( jvals == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_Real* vals;

      /* fill long array with SCIP variable pointers */
      vals = SCIPgetValsLinear(scip, cons);
      (*env)->SetDoubleArrayRegion(env, jvals, 0, nvars, (jdouble*)vals);
   }

   return jvals;
}

/** gets the activity of the linear constraint in the given solution */
jdouble JNISCIPCONSLINEAR(getActivityLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint data */
   jlong                 jsol                /**< solution, or NULL to use current node's solution */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert( cons != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;

   num = SCIPgetActivityLinear(scip, cons, sol);

   return (jdouble) num;
}

/** gets the feasibility of the linear constraint in the given solution */
jdouble JNISCIPCONSLINEAR(getFeasibilityLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint data */
   jlong                 jsol                /**< solution, or NULL to use current node's solution */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert( cons != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;

   num = SCIPgetFeasibilityLinear(scip, cons, sol);

   return (jdouble) num;
}

/** gets the dual solution of the linear constraint in the current LP */
jdouble JNISCIPCONSLINEAR(getDualsolLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert( cons != NULL);


   num = SCIPgetDualsolLinear(scip, cons);

   return (jdouble) num;
}

/** gets the dual Farkas value of the linear constraint in the current infeasible LP */
jdouble JNISCIPCONSLINEAR(getDualfarkasLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert( cons != NULL);


   num = SCIPgetDualfarkasLinear(scip, cons);

   return (jdouble) num;
}

/** returns the linear relaxation of the given linear constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
jlong JNISCIPCONSLINEAR(getRowLinear)(
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
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   return (jlong) (size_t) SCIPgetRowLinear(scip, cons);
}

/** tries to automatically convert a linear constraint into a more specific and more specialized constraint */
jlong JNISCIPCONSLINEAR(upgradeConsLinear)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_CONS* upgdcons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPupgradeConsLinear(scip, cons, &upgdcons) );

   return (jlong) (size_t) upgdcons;
}

/**@} */
