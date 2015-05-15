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

#include "JniScipSol.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include <string.h>

/** gets origin of solution */
JNIEXPORT
jint JNISCIPSOL(solGetOrigin)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsol                /**< primal CIP solution */
   )
{
   SCIP_SOL* sol;
   SCIP_SOLORIGIN orig;

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   orig = SCIPsolGetOrigin(sol);

   return (jint) orig;
}

/** returns whether the given solution is defined on original variables */
JNIEXPORT
jboolean JNISCIPSOL(solIsOriginal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsol                /**< primal CIP solution */
   )
{
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   return (jboolean) SCIPsolIsOriginal(sol);
}

/** gets objective value of primal CIP solution which lives in the original problem space */
JNIEXPORT
jdouble JNISCIPSOL(solGetOrigObj)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsol                /**< primal CIP solution */
   )
{
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   return (jdouble) SCIPsolGetOrigObj(sol);
}

/** gets clock time, when this solution was found */
JNIEXPORT
jdouble JNISCIPSOL(solGetTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsol                /**< primal CIP solution */
   )
{
   SCIP_SOL* sol;
   SCIP_Real time;

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   time = SCIPsolGetTime(sol);

   return (jdouble) time;
}

/** gets branch and bound run number, where this solution was found */
JNIEXPORT
jint JNISCIPSOL(solGetRunnum)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsol                /**< primal CIP solution */
   )
{
   SCIP_SOL* sol;
   int runnum;

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   runnum = SCIPsolGetRunnum(sol);

   return (jint) runnum;
}

/** gets node number of the specific branch and bound run, where this solution was found */
JNIEXPORT
jlong JNISCIPSOL(solGetNodenum)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsol                /**< primal CIP solution */
   )
{
   SCIP_SOL* sol;
   SCIP_Longint nodenum;

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   nodenum = SCIPsolGetNodenum(sol);

   return (jlong) nodenum;
}

/** gets node's depth, where this solution was found */
JNIEXPORT
jint JNISCIPSOL(solGetDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsol                /**< primal CIP solution */
   )
{
   SCIP_SOL* sol;
   int depth;

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   depth = SCIPsolGetDepth(sol);

   return (jint) depth;
}

/** gets heuristic, that found this solution (or NULL if it's from the tree) */
JNIEXPORT
jlong JNISCIPSOL(solGetHeur)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsol                /**< primal CIP solution */
   )
{
   SCIP_SOL* sol;
   SCIP_HEUR* heur;

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   heur = SCIPsolGetHeur(sol);

   return (jlong)(size_t) heur;
}

/** informs the solution that it now belongs to the given primal heuristic */
JNIEXPORT
void JNISCIPSOL(solSetHeur)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsol,               /**< primal CIP solution */
   jlong                 jheur               /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_SOL* sol;
   SCIP_HEUR* heur;

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   heur = (SCIP_HEUR*) (size_t) jheur;
   assert(heur != NULL);

   SCIPsolSetHeur(sol, heur);
}

/** returns unique index of given solution */
JNIEXPORT
jint JNISCIPSOL(solGetIndex)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsol                /**< primal CIP solution */
   )
{
   SCIP_SOL* sol;
   int ind;

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   ind = SCIPsolGetIndex(sol);

   return (jint) ind;
}
