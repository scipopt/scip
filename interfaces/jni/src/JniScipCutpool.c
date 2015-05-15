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

/**@file   JniScipCutpool.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP cutpool callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipCutpool.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/pub_cutpool.h"

#include <string.h>

/** gets the row of the cut */
JNIEXPORT
jlong JNISCIPCUTPOOL(cutGetRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcut                /**< cut */
   )
{
   SCIP_CUT* cut;

   /* convert JNI pointer into C pointer */
   cut = (SCIP_CUT*) (size_t) jcut;
   assert(cut != NULL);

   return (jlong) (size_t) SCIPcutGetRow(cut);
}

/** gets the age of the cut: the number of consecutive cut pool separation rounds where the cut was neither in the LP nor violated */
JNIEXPORT
jint JNISCIPCUTPOOL(cutGetAge)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcut                /**< cut */
   )
{
   SCIP_CUT* cut;

   /* convert JNI pointer into C pointer */
   cut = (SCIP_CUT*) (size_t) jcut;
   assert(cut != NULL);

   return (jint) SCIPcutGetAge(cut);
}

/** gets array of cuts in the cut pool */
JNIEXPORT
jlongArray JNISCIPCUTPOOL(cutpoolGetCuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcutpool            /**< cut */
   )
{
   SCIP_CUTPOOL* cutpool;
   jlongArray jcuts;
   int size;
   SCIP_CUT** cuts;

   /* convert JNI pointer into C pointer */
   cutpool = (SCIP_CUTPOOL*) (size_t) jcutpool;
   assert(cutpool != NULL);

   size = SCIPcutpoolGetNCuts(cutpool);

   jcuts = (*env)->NewLongArray(env, size);

   if (jcuts == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   cuts = SCIPcutpoolGetCuts(cutpool);

   (*env)->SetLongArrayRegion(env, jcuts, 0, size, (jlong*)(*cuts));

   return jcuts;
}

/** get number of cuts in the cut pool */
JNIEXPORT
jint JNISCIPCUTPOOL(cutpoolGetNCuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcutpool                /**< cut */
   )
{
   SCIP_CUTPOOL* cutpool;

   /* convert JNI pointer into C pointer */
   cutpool = (SCIP_CUTPOOL*) (size_t) jcutpool;
   assert(cutpool != NULL);

   return (jint) SCIPcutpoolGetNCuts(cutpool);
}

/** get maximum number of cuts that were stored in the cut pool at the same time */
JNIEXPORT
jint JNISCIPCUTPOOL(cutpoolGetMaxNCuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcutpool                /**< cut */
   )
{
   SCIP_CUTPOOL* cutpool;

   /* convert JNI pointer into C pointer */
   cutpool = (SCIP_CUTPOOL*) (size_t) jcutpool;
   assert(cutpool != NULL);

   return (jint) SCIPcutpoolGetMaxNCuts(cutpool);
}

/** gets time in seconds used for separating cuts from the pool */
JNIEXPORT
jdouble JNISCIPCUTPOOL(cutpoolGetTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcutpool                /**< cut */
   )
{
   SCIP_CUTPOOL* cutpool;

   /* convert JNI pointer into C pointer */
   cutpool = (SCIP_CUTPOOL*) (size_t) jcutpool;
   assert(cutpool != NULL);

   return (jdouble) SCIPcutpoolGetTime(cutpool);
}

/** get number of times, the cut pool was separated */
JNIEXPORT
jlong JNISCIPCUTPOOL(cutpoolGetNCalls)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcutpool                /**< cut */
   )
{
   SCIP_CUTPOOL* cutpool;

   /* convert JNI pointer into C pointer */
   cutpool = (SCIP_CUTPOOL*) (size_t) jcutpool;
   assert(cutpool != NULL);

   return (jlong) SCIPcutpoolGetNCalls(cutpool);
}

/** get total number of cuts that were separated from the cut pool */
JNIEXPORT
jlong JNISCIPCUTPOOL(cutpoolGetNCutsFound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcutpool                /**< cut */
   )
{
   SCIP_CUTPOOL* cutpool;

   /* convert JNI pointer into C pointer */
   cutpool = (SCIP_CUTPOOL*) (size_t) jcutpool;
   assert(cutpool != NULL);

   return (jlong) SCIPcutpoolGetNCutsFound(cutpool);
}
