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

/**@file   JniScipEvent.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP event callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipEvent.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/pub_event.h"

#include <string.h>

/** gets user data of display column */
#if 0
JNIEXPORT
jlong JNISCIPEVENT(dispGetData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdisp               /**< display column */
   )
{
   SCIP_DISP* disp;

   /* convert JNI pointer into C pointer */
   disp = (SCIP_DISP*) (size_t) jdisp;
   assert(disp != NULL);

   return (jlong) (size_t) SCIPdispGetData(disp);
}
#endif

/**
eventhdlrGetName
eventhdlrGetData
eventhdlrSetData
eventhdlrIsInitialized
eventhdlrGetSetupTime
eventhdlrGetTime
eventGetType
eventGetVar
eventGetOldobj
eventGetNewobj
eventGetOldbound
eventGetNewbound
eventGetNode
eventGetSol
eventGetHoleLeft
eventGetHoleRight
eventGetRow
eventGetRowCol
eventGetRowOldCoefVal
eventGetRowNewCoefVal
eventGetRowOldConstVal
eventGetRowNewConstVal
eventGetRowSide
eventGetRowOldSideVal
eventGetRowNewSideVal
*/
