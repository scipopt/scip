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

/**@file   JniScipConsVar.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP constraint variable callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipConsVarbound.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_varbound.h"

#include <string.h>


/** creates the handler for variable bound constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSVARBOUND(includeConshdlrVarbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrVarbound(scip) );
}

/** creates and captures a variable bound constraint: lhs <= x + c*y <= rhs
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method {@link releaseCons()}
 */
JNIEXPORT
jlong JNISCIPCONSVARBOUND(createConsVarbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jlong                 jvar,               /**< variable x that has variable bound */
   jlong                 jvbdvar,            /**< binary, integer or implicit integer bounding variable y */
   jdouble               jvbdcoef,           /**< coefficient c of bounding variable y */
   jdouble               jlhs,               /**< left hand side of variable bound inequality */
   jdouble               jrhs,               /**< right hand side of variable bound inequality */
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
   SCIP_VAR* var;
   SCIP_VAR* vbdvar;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   /* convert JNI pointer into C pointer */
   vbdvar = (SCIP_VAR*) (size_t) jvbdvar;
   assert(vbdvar != NULL);

   JNISCIP_CALL( SCIPcreateConsVarbound(scip, &cons, name, var, vbdvar, (SCIP_Real)jvbdcoef, (SCIP_Real)jlhs, (SCIP_Real)jrhs, (SCIP_Bool) jinitial, (SCIP_Bool) jseparate, (SCIP_Bool) jenforce, (SCIP_Bool) jcheck, (SCIP_Bool) jpropagate,
	 (SCIP_Bool) jlocal, (SCIP_Bool) jmodifiable, (SCIP_Bool) jdynamic, (SCIP_Bool) jremovable, (SCIP_Bool) jstickingatnode) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** gets left hand side of variable bound constraint lhs <= x + c*y <= rhs */
JNIEXPORT
jdouble JNISCIPCONSVARBOUND(getLhsVarbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Real lhs;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   lhs = SCIPgetLhsVarbound(scip, cons);

   return (jdouble)lhs;
}

/** gets right hand side of variable bound constraint lhs <= x + c*y <= rhs */
JNIEXPORT
jdouble JNISCIPCONSVARBOUND(getRhsVarbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Real rhs;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   rhs = SCIPgetRhsVarbound(scip, cons);

   return (jdouble)rhs;
}

/** gets bounded variable x of variable bound constraint lhs <= x + c*y <= rhs */
JNIEXPORT
jlong JNISCIPCONSVARBOUND(getVarVarbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
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

   var = SCIPgetVarVarbound(scip, cons);

   return (jlong)(size_t)var;
}

/** gets bounding variable y of variable bound constraint lhs <= x + c*y <= rhs */
JNIEXPORT
jlong JNISCIPCONSVARBOUND(getVbdvarVarbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
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

   var = SCIPgetVbdvarVarbound(scip, cons);

   return (jlong)(size_t)var;
}

/** gets bound coefficient c of variable bound constraint lhs <= x + c*y <= rhs */
JNIEXPORT
jdouble JNISCIPCONSVARBOUND(getVbdcoefVarbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint data */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Real bdcoef;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   bdcoef = SCIPgetVbdcoefVarbound(scip, cons);

   return (jdouble)bdcoef;
}

/** gets the dual solution of the variable bound constraint in the current LP */
JNIEXPORT
jdouble JNISCIPCONSVARBOUND(getDualsolVarbound)(
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

   dualsol = SCIPgetDualsolVarbound(scip, cons);

   return (jdouble)dualsol;
}

/** gets the dual farkas value of the variable bound constraint in the current infeasible LP */
JNIEXPORT
jdouble JNISCIPCONSVARBOUND(getDualfarkasVarbound)(
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

   dualfarkas = SCIPgetDualfarkasVarbound(scip, cons);

   return (jdouble)dualfarkas;
}

/** returns the linear relaxation of the given variable bound constraint; may return NULL if no LP row was yet created; the user must not modify the row! */
JNIEXPORT
jlong JNISCIPCONSVARBOUND(getRowVarbound)(
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

   row = SCIPgetRowVarbound(scip, cons);

   return (jlong)(size_t)row;
}
