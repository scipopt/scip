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

#include "JniScipConsOrbitope.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include <string.h>

/** creates the handler for orbitope constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSORBITOPE(includeConshdlrOrbitope)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrOrbitope(scip) );
}

/** creates and captures a orbitope constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSORBITOPE(createConsOrbitope)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jobjectArray          jvars,              /**< matrix of variables on which the symmetry acts */
   jboolean              ispart,             /**< whether we deal with the partitioning case (packing otherwise) */
   jint                  nspcons,            /**< number of set partitioning/packing constraints  <=> p */
   jint                  nblocks,            /**< number of symmetric variable blocks             <=> q */
   jboolean              resolveprop,        /**< should propagation be resolved? */
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
                                              *   are separated as constraints. */
   jboolean              removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   jboolean              stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIPerrorMessage("method createConsOrbitope is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** creates and captures an orbitope constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
JNIEXPORT
jlong JNISCIPCONSORBITOPE(createConsBasicOrbitope)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jobjectArray          jvars,              /**< matrix of variables on which the symmetry acts */
   jboolean              ispart,             /**< whether we deal with the partitioning case (packing otherwise) */
   jint                  nspcons,            /**< number of set partitioning/packing constraints  <=> p */
   jint                  nblocks,            /**< number of symmetric variable blocks             <=> q */
   jboolean              resolveprop         /**< should propagation be resolved? */
   )
{
   SCIPerrorMessage("method createConsBasicOrbitope  is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}
