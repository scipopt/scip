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

/**@file   JniScipReaderGms.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP GAMS file writer
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipReaderGms.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/reader_gms.h"

#include <string.h>

/** includes the gms file reader into SCIP */
JNIEXPORT
void JNISCIPREADERGMS(includeReaderGms)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeReaderGms(scip) );
}

/** writes problem to file */
JNIEXPORT
jint JNISCIPREADERGMS(writeGms)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jfile,              /**< output file, or NULL if standard output should be used */
   jstring               jname,              /**< problem name */
   jboolean              jtransformed,       /**< TRUE iff problem is the transformed problem */
   jint                  jobjsense,           /**< objective sense */
   jdouble               jobjscale,           /**< scalar applied to objective function; external objective value is
                                              *   extobj = objsense * objscale * (intobj + objoffset) */
   jdouble               jobjoffset,          /**< objective offset from bound shifting and fixing */
   jlongArray            jvars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   jint                  jnvars,              /**< number of mutable variables in the problem */
   jint                  jnbinvars,           /**< number of binary variables */
   jint                  jnintvars,           /**< number of general integer variables */
   jint                  jnimplvars,          /**< number of implicit integer variables */
   jint                  jncontvars,          /**< number of continuous variables */
   jlongArray            jconss,             /**< array with constraints of the problem */
   jint                  jnconss             /**< number of constraints in the problem */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_RESULT result;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;

   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)jnvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &conss, (int)jnconss) );

   (*env)->GetLongArrayRegion(env, jvars, 0, jnvars, (jlong*)(*vars));
   (*env)->GetLongArrayRegion(env, jconss, 0, jnconss, (jlong*)(*conss));

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPwriteGms(scip, (FILE*)(size_t) jfile, name, (SCIP_Bool)jtransformed, (SCIP_OBJSENSE)jobjsense, (SCIP_Real)jobjscale, (SCIP_Real)jobjoffset, vars, (int)jnvars, (int)jnbinvars, (int)jnintvars, (int)jnimplvars, (int)jncontvars, conss, (int)jnconss, &result) );

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &conss);
   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jint) result;
}
