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

/**@file   JniScipHeurTrysol.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP heur trysol callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipHeurTrysol.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/heur_trysol.h"

#include <string.h>

/** creates the trysol primal heuristic and includes it in SCIP */
JNIEXPORT
void JNISCIPHEURTRYSOL(includeHeurTrySol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeHeurTrySol(scip) );
}

/** pass solution to trysol heuristic */
JNIEXPORT
void JNISCIPHEURTRYSOL(heurPassSolTrySol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jheur,              /**< trysol heuristic */
   jlong                 jsol                /**< solution to be passed */
   )
{
   SCIP* scip;
   SCIP_HEUR* heur;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   heur = (SCIP_HEUR*) (size_t) jheur;
   sol = (SCIP_SOL*) (size_t) jsol;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(sol != NULL);

   JNISCIP_CALL( SCIPheurPassSolTrySol(scip, heur, sol) );
}

/** pass solution to trysol heuristic which just gets added (without checking feasibility */
JNIEXPORT
void JNISCIPHEURTRYSOL(heurPassSolAddSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jheur,              /**< trysol heuristic */
   jlong                 jsol                /**< solution to be passed */
   )
{
   SCIP* scip;
   SCIP_HEUR* heur;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   heur = (SCIP_HEUR*) (size_t) jheur;
   sol = (SCIP_SOL*) (size_t) jsol;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(sol != NULL);

   JNISCIP_CALL( SCIPheurPassSolAddSol(scip, heur, sol) );
}
