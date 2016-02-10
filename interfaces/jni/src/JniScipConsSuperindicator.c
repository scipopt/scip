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

/**@file   JniScipConsSuperindicator.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP sos2 constraint callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipConsSuperindicator.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_superindicator.h"

#include <string.h>

/** creates the handler for superindicator constraints and includes it in SCIP */
JNIEXPORT
void JNISCIPCONSSUPERINDICATOR(includeConshdlrSuperindicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeConshdlrSuperindicator(scip) );
}

/** gets binary variable corresponding to the general indicator constraint */
JNIEXPORT
jlong JNISCIPCONSSUPERINDICATOR(getBinaryVarSuperindicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< superindicator constraint */
   )
{
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   return (jlong) (size_t) SCIPgetBinaryVarSuperindicator(cons);
}

/** gets the slack constraint corresponding to the general indicator constraint */
JNIEXPORT
jlong JNISCIPCONSSUPERINDICATOR(getSlackConsSuperindicator)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcons               /**< superindicator constraint */
   )
{
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   return (jlong) (size_t) SCIPgetSlackConsSuperindicator(cons);
}

/** transforms the current problem into a MinUC problem (minimizing the number of unsatisfied constraints),
 *  a CIP generalization of the MinULR (min. unsatisfied linear relations) problem
 */
JNIEXPORT
jboolean JNISCIPCONSSUPERINDICATOR(transformMinUC)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Bool success ;
   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);
   JNISCIP_CALL( SCIPtransformMinUC(scip, &success) );

   return (jboolean) success;
}
