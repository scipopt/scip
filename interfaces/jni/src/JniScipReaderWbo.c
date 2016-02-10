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

/**@file   JniScipReaderWbo.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP WBO file reader (LP format with generic variables and row names)
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipReaderWbo.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/reader_wbo.h"

#include <string.h>

/** includes the wbo file reader into SCIP */
JNIEXPORT
void JNISCIPREADERWBO(includeReaderWbo)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeReaderWbo(scip) );
}
