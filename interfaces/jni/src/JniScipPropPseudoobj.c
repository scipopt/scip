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

/**@file   JniScipPropPseudoobj.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP Pseudo objective propagator
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipPropPseudoobj.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/prop_pseudoobj.h"

#include <string.h>

/** creates the pseudo objective function propagator and includes it in SCIP */
JNIEXPORT
void JNISCIPPROPPSEUDOOBJ(includePropPseudoobj)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludePropPseudoobj(scip) );
}

/** propagates the cutoff bound for the given variables */
JNIEXPORT
jboolean JNISCIPPROPPSEUDOOBJ(propagateCutoffboundVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jprop,              /**< propagator, or NULL */
   jlong                 jvar,               /**< variables to propagate */
   jdouble               jcutoffbound,       /**< cutoff bound to use */
   jdouble               jpseudoobjval       /**< pseudo objective value to use */
   )
{
   SCIP* scip;
   SCIP_PROP* prop;
   SCIP_VAR* var;
   SCIP_Bool tightened;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   prop = (SCIP_PROP*) (size_t) jprop;
   var = (SCIP_VAR*) (size_t) jvar;

   assert(scip != NULL);
   assert(var != NULL);

   JNISCIP_CALL( SCIPpropagateCutoffboundVar(scip, prop, var, (SCIP_Real)jcutoffbound, (SCIP_Real)jpseudoobjval, &tightened) );

   return (jboolean) tightened;
}
