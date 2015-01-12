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

/**@file   JniScipHeurUndercover.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP heur undercover callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipHeurUndercover.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/heur_undercover.h"

#include <string.h>

/** creates the undercover primal heuristic and includes it in SCIP */
JNIEXPORT
void JNISCIPHEURUNDERCOVER(includeHeurUndercover)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeHeurUndercover(scip) );
}

/** computes a minimal set of covering variables */
JNIEXPORT
jlong JNISCIPHEURUNDERCOVER(computeCoverUndercover)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jintArray             jcoversize,         /**< size of the computed cover */
   jdouble               jtimelimit,         /**< time limit */
   jdouble               jmemorylimit,       /**< memory limit */
   jdouble               jobjlimit,          /**< objective limit: upper bound on coversize */
   jboolean              jglobalbounds,      /**< should global bounds on variables be used instead of local bounds at focus node? */
   jboolean              jonlyconvexify,     /**< should we only fix/dom.red. variables creating nonconvexity? */
   jboolean              jcoverbd,           /**< should bounddisjunction constraints be covered (or just copied)? */
   jchar                 jcoveringobj,       /**< objective function of the covering problem ('b'ranching status,
                                              *   influenced nonlinear 'c'onstraints/'t'erms, 'd'omain size, 'l'ocks,
                                              *   'm'in of up/down locks, 'u'nit penalties, constraint 'v'iolation) */
   jbooleanArray         jsuccess            /**< feasible cover found? */
   )
{
   SCIPerrorMessage("method computeCoverUndercover is not implemented yet (deprecated)\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}
