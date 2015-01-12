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

/**@file   JniScipReaderPpm.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP file writer for portable pixmap file format (PPM), open with common graphic viewer programs (e.g. xview)
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipReaderPpm.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/reader_ppm.h"

#include <string.h>

/** includes the ppm file reader into SCIP */
JNIEXPORT
void JNISCIPREADERPPM(includeReaderPpm)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeReaderPpm(scip) );
}

/** writes problem to file */
JNIEXPORT
jint JNISCIPREADERPPM(writePpm)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jfile,              /**< output file, or NULL if standard output should be used */
   jstring               jname,              /**< problem name */
   jlong                 jreaderdata,        /**< information for reader */
   jboolean              jtransformed,       /**< TRUE iff problem is the transformed problem */
   jlongArray            jvars,              /**< array with active variables ordered binary, integer, implicit, continuous */
   jint                  jnvars,             /**< number of mutable variables in the problem */
   jlongArray            jconss,             /**< array with constraints of the problem */
   jint                  jnconss             /**< number of constraints in the problem */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_READERDATA* readerdata;
   SCIP_RESULT result;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   readerdata = (SCIP_READERDATA*) (size_t) jreaderdata;

   assert(scip != NULL);
   assert(readerdata != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)jnvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &conss, (int)jnconss) );

   (*env)->GetLongArrayRegion(env, jvars, 0, jnvars, (jlong*)(*vars));
   (*env)->GetLongArrayRegion(env, jconss, 0, jnconss, (jlong*)(*conss));

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPwritePpm(scip, (FILE*)(size_t) jfile, name, readerdata, (SCIP_Bool)jtransformed, vars, (int)jnvars, conss, (int)jnconss, &result) );

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &conss);
   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jint) result;
}
