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

/**@file   JniScipReaderPip.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP file reader for polynomial mixed-integer programs in PIP format
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipReaderPip.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/reader_pip.h"

#include <string.h>

/** includes the pip file reader into SCIP */
JNIEXPORT
void JNISCIPREADERPIP(includeReaderPip)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeReaderPip(scip) );
}

/** reads problem from file */
JNIEXPORT
jint JNISCIPREADERPIP(readPip)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jreader,            /**< the file reader itself */
   jstring               jname               /**< full path and name of file to read, or NULL if stdin should be used */
   )
{
   SCIP* scip;
   SCIP_READER* reader;
   const char* name;
   jboolean iscopy;
   SCIP_RESULT result;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   reader = (SCIP_READER*) (size_t) jreader;

   assert(scip != NULL);
   assert(reader != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   assert(iscopy);

   JNISCIP_CALL( SCIPreadPip(scip, reader, name, &result) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jint) result;
}

/** writes problem to file */
JNIEXPORT
jint JNISCIPREADERPIP(writePip)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jfile,              /**< output file, or NULL if standard output should be used */
   jstring               jname,              /**< problem name */
   jboolean              jtransformed,       /**< TRUE iff problem is the transformed problem */
   jint                  jobjsense,          /**< objective sense */
   jdouble               jobjscale,          /**< scalar applied to objective function; external objective value is
                                              *   extobj = objsense * objscale * (intobj + objoffset) */
   jdouble               jobjoffset,         /**< objective offset from bound shifting and fixing */
   jlongArray            jvars,              /**< array with active variables ordered binary, integer, implicit, continuous */
   jint                  jnvars,             /**< number of mutable variables in the problem */
   jint                  jnbinvars,          /**< number of binary variables */
   jint                  jnintvars,          /**< number of general integer variables */
   jint                  jnimplvars,         /**< number of implicit integer variables */
   jint                  jncontvars,         /**< number of continuous variables */
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

   JNISCIP_CALL( SCIPwritePip(scip, (FILE*)(size_t) jfile, name, (SCIP_Bool)jtransformed, (SCIP_OBJSENSE)jobjsense, (SCIP_Real)jobjscale, (SCIP_Real)jobjoffset, vars, (int)jnvars, (int)jnbinvars, (int)jnintvars, (jint)jnimplvars, (int)jncontvars, conss, (int)jnconss, &result) );

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &conss);
   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jint) result;
}
