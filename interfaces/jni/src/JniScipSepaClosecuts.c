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

/**@file   JniScipSepaClique.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP closecuts meta separator
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipSepaClosecuts.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/sepa_closecuts.h"

#include <string.h>

/** creates the closecuts separator and includes it in SCIP */
JNIEXPORT
void JNISCIPSEPACLOSECUTS(includeSepaClosecuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeSepaClosecuts(scip) );
}

/** sets point to be used as base point for computing the point to be separated
 *
 *  The point is only stored if separation of relative interior points is used. The solution is copied.
 */
JNIEXPORT
void JNISCIPSEPACLOSECUTS(setBasePointClosecuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< base point solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   sol = (SCIP_SOL*) (size_t) jsol;

   assert(scip != NULL);
   assert(sol != NULL);

   JNISCIP_CALL( SCIPsetBasePointClosecuts(scip, sol) );
}
