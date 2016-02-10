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

/**@file   JniScip.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScip.h"

#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include <string.h>
#include <locale.h>

#ifndef NDEBUG
#include "JniScipBoundchgtype.h"
#include "JniScipBoundtype.h"
#include "JniScipBranchdir.h"
#include "JniScipObjsense.h"
#include "JniScipParamemphasis.h"
#include "JniScipParamsetting.h"
#include "JniScipSolorigin.h"
#include "JniScipStage.h"
#include "JniScipStatus.h"
#include "JniScipVartype.h"
#include "JniScipVarstatus.h"
#endif



#ifndef NDEBUG
/** check all eumns of SCIP against the JNI enums */
static
void checkEnums(
   void
   )
{
   /* check that the STATUS enums of JNI SCIP are mapped correctly to the ones of SCIP */
   assert(SCIP_STATUS_UNKNOWN == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_UNKNOWN));;
   assert(SCIP_STATUS_USERINTERRUPT == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_USERINTERRUPT));
   assert(SCIP_STATUS_NODELIMIT == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_NODELIMIT));
   assert(SCIP_STATUS_STALLNODELIMIT == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_STALLNODELIMIT));
   assert(SCIP_STATUS_TIMELIMIT == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_TIMELIMIT));
   assert(SCIP_STATUS_MEMLIMIT == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_MEMLIMIT));
   assert(SCIP_STATUS_GAPLIMIT == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_GAPLIMIT));
   assert(SCIP_STATUS_SOLLIMIT == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_SOLLIMIT));
   assert(SCIP_STATUS_BESTSOLLIMIT == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_BESTSOLLIMIT));
   assert(SCIP_STATUS_OPTIMAL == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_OPTIMAL));
   assert(SCIP_STATUS_INFEASIBLE == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_INFEASIBLE));
   assert(SCIP_STATUS_UNBOUNDED == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_UNBOUNDED));
   assert(SCIP_STATUS_INFORUNBD == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_INFORUNBD));

   /* check if the STAGE enums of JNI SCIP are mapped correctly to the ones of SCIP */
   assert(SCIP_STAGE_INIT == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_INIT));
   assert(SCIP_STAGE_PROBLEM == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_PROBLEM));
   assert(SCIP_STAGE_TRANSFORMING == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_TRANSFORMING));
   assert(SCIP_STAGE_TRANSFORMED == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_TRANSFORMED));
   assert(SCIP_STAGE_INITPRESOLVE == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_INITPRESOLVE));
   assert(SCIP_STAGE_PRESOLVING == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_PRESOLVING));
   assert(SCIP_STAGE_EXITPRESOLVE == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_EXITPRESOLVE));
   assert(SCIP_STAGE_PRESOLVED == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_PRESOLVED));
   assert(SCIP_STAGE_INITSOLVE == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_INITSOLVE));
   assert(SCIP_STAGE_SOLVING == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_SOLVING));
   assert(SCIP_STAGE_SOLVED == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_SOLVED));
   assert(SCIP_STAGE_EXITSOLVE == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_EXITSOLVE));
   assert(SCIP_STAGE_FREETRANS == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_FREETRANS));
   assert(SCIP_STAGE_FREE == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_FREE));

   /* check if the OBJECTIVE SENSE enums of JNI SCIP are mapped correctly to the ones of SCIP */
   assert(SCIP_OBJSENSE_MAXIMIZE == JNIPACKAGENAME(JniScipObjsense_SCIP_OBJSENSE_MAXIMIZE));
   assert(SCIP_OBJSENSE_MINIMIZE == JNIPACKAGENAME(JniScipObjsense_SCIP_OBJSENSE_MINIMIZE));

   /* check if the VARTYPE enums of JNI SCIP are mapped correctly to the ones of SCIP */
   assert(SCIP_VARTYPE_BINARY == JNIPACKAGENAME(JniScipVartype_SCIP_VARTYPE_BINARY));
   assert(SCIP_VARTYPE_INTEGER == JNIPACKAGENAME(JniScipVartype_SCIP_VARTYPE_INTEGER));
   assert(SCIP_VARTYPE_IMPLINT == JNIPACKAGENAME(JniScipVartype_SCIP_VARTYPE_IMPLINT));
   assert(SCIP_VARTYPE_CONTINUOUS == JNIPACKAGENAME(JniScipVartype_SCIP_VARTYPE_CONTINUOUS));

   /* check if the VARSTATUS enums of JNI SCIP are mapped correctly to the ones of SCIP */
   assert(SCIP_VARSTATUS_ORIGINAL == JNIPACKAGENAME(JniScipVarstatus_SCIP_VARSTATUS_ORIGINAL));
   assert(SCIP_VARSTATUS_LOOSE == JNIPACKAGENAME(JniScipVarstatus_SCIP_VARSTATUS_LOOSE));
   assert(SCIP_VARSTATUS_COLUMN == JNIPACKAGENAME(JniScipVarstatus_SCIP_VARSTATUS_COLUMN));
   assert(SCIP_VARSTATUS_FIXED == JNIPACKAGENAME(JniScipVarstatus_SCIP_VARSTATUS_FIXED));
   assert(SCIP_VARSTATUS_AGGREGATED == JNIPACKAGENAME(JniScipVarstatus_SCIP_VARSTATUS_AGGREGATED));
   assert(SCIP_VARSTATUS_MULTAGGR == JNIPACKAGENAME(JniScipVarstatus_SCIP_VARSTATUS_MULTAGGR));
   assert(SCIP_VARSTATUS_NEGATED == JNIPACKAGENAME(JniScipVarstatus_SCIP_VARSTATUS_NEGATED));

   /* check if the SOLORIGIN enums of JNI SCIP are mapped correctly to the ones of SCIP */
   assert(SCIP_SOLORIGIN_ORIGINAL == JNIPACKAGENAME(JniScipSolorigin_SCIP_SOLORIGIN_ORIGINAL));
   assert(SCIP_SOLORIGIN_ZERO == JNIPACKAGENAME(JniScipSolorigin_SCIP_SOLORIGIN_ZERO));
   assert(SCIP_SOLORIGIN_LPSOL == JNIPACKAGENAME(JniScipSolorigin_SCIP_SOLORIGIN_LPSOL));
   assert(SCIP_SOLORIGIN_NLPSOL == JNIPACKAGENAME(JniScipSolorigin_SCIP_SOLORIGIN_NLPSOL));
   assert(SCIP_SOLORIGIN_RELAXSOL == JNIPACKAGENAME(JniScipSolorigin_SCIP_SOLORIGIN_RELAXSOL));
   assert(SCIP_SOLORIGIN_PSEUDOSOL == JNIPACKAGENAME(JniScipSolorigin_SCIP_SOLORIGIN_PSEUDOSOL));
   assert(SCIP_SOLORIGIN_UNKNOWN == JNIPACKAGENAME(JniScipSolorigin_SCIP_SOLORIGIN_UNKNOWN));

   /* check if the BOUNDTYPE enums of JNI SCIP are mapped correctly to the ones of SCIP */
   assert(SCIP_BOUNDTYPE_LOWER == JNIPACKAGENAME(JniScipBoundtype_SCIP_BOUNDTYPE_LOWER));
   assert(SCIP_BOUNDTYPE_UPPER == JNIPACKAGENAME(JniScipBoundtype_SCIP_BOUNDTYPE_UPPER));

   /* check if the BRANCHDIR enums of JNI SCIP are mapped correctly to the ones of SCIP */
   assert(SCIP_BRANCHDIR_DOWNWARDS == JNIPACKAGENAME(JniScipBranchdir_SCIP_BRANCHDIR_DOWNWARDS));
   assert(SCIP_BRANCHDIR_UPWARDS == JNIPACKAGENAME(JniScipBranchdir_SCIP_BRANCHDIR_UPWARDS));
   assert(SCIP_BRANCHDIR_FIXED == JNIPACKAGENAME(JniScipBranchdir_SCIP_BRANCHDIR_FIXED));
   assert(SCIP_BRANCHDIR_AUTO == JNIPACKAGENAME(JniScipBranchdir_SCIP_BRANCHDIR_AUTO));

   /* check if the BOUNDCHGTYPE enums of JNI SCIP are mapped correctly to the ones of SCIP */
   assert(SCIP_BOUNDCHGTYPE_BRANCHING == JNIPACKAGENAME(JniScipBoundchgtype_SCIP_BOUNDCHGTYPE_BRANCHING));
   assert(SCIP_BOUNDCHGTYPE_CONSINFER == JNIPACKAGENAME(JniScipBoundchgtype_SCIP_BOUNDCHGTYPE_CONSINFER));
   assert(SCIP_BOUNDCHGTYPE_PROPINFER == JNIPACKAGENAME(JniScipBoundchgtype_SCIP_BOUNDCHGTYPE_PROPINFER));

   /* check if the EMPHASIS enums of JNI SCIP are mapped correctly to the ones of SCIP */
   assert(SCIP_PARAMEMPHASIS_COUNTER == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_COUNTER));
   assert(SCIP_PARAMEMPHASIS_CPSOLVER == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_CPSOLVER));
   assert(SCIP_PARAMEMPHASIS_EASYCIP == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_EASYCIP));
   assert(SCIP_PARAMEMPHASIS_FEASIBILITY == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_FEASIBILITY));
   assert(SCIP_PARAMEMPHASIS_HARDLP == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_HARDLP));
   assert(SCIP_PARAMEMPHASIS_OPTIMALITY == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_OPTIMALITY));

   /* check if the SETTINGS enums of JNI SCIP are mapped correctly to the ones of SCIP */
   assert(SCIP_PARAMSETTING_DEFAULT == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_DEFAULT));
   assert(SCIP_PARAMSETTING_AGGRESSIVE == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_AGGRESSIVE));
   assert(SCIP_PARAMSETTING_FAST == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_FAST));
   assert(SCIP_PARAMSETTING_OFF == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_OFF));
}
#endif



/** returns scip version number */
JNIEXPORT
jdouble JNISCIP(version)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj                /**< JNI class pointer */
   )
{
   SCIP_Real num;

   num = SCIPversion();

   return (jdouble) num;
}

/** returns SCIP major version */
JNIEXPORT
jint JNISCIP(majorVersion)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj                /**< JNI class pointer */
   )
{
   int num;

   num = SCIPmajorVersion();

   return (jint) num;
}

/** returns SCIP minor version */
JNIEXPORT
jint JNISCIP(minorVersion)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj                /**< JNI class pointer */
   )
{
   int num;

   num = SCIPminorVersion();

   return (jint) num;
}

/** returns SCIP technical version */
JNIEXPORT
jint JNISCIP(techVersion)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj                /**< JNI class pointer */
   )
{
   int num;

   num = SCIPtechVersion();

   return (jint) num;
}

/** returns SCIP sub version number */
JNIEXPORT
jint JNISCIP(subversion)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj                /**< JNI class pointer */
   )
{
   int num;

   num = SCIPsubversion();

   return (jint) num;
}

/** prints a version information line to a file stream via the message handler system
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
JNIEXPORT
void JNISCIP(printVersion)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jfile               /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   SCIPprintVersion(scip, (FILE*)(size_t)jfile);
}

/** prints error message for the given SCIP_RETCODE via the error prints method */
JNIEXPORT
void JNISCIP(printError)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jint                  jredcode            /**< SCIP return code causing the error */
   )
{
   SCIPprintError((SCIP_RETCODE) jredcode);
}

/** creates and initializes SCIP data structures */
JNIEXPORT
jlong JNISCIP(create)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj                /**< JNI class pointer */
   )
{
   SCIP* scip;

#ifndef NDEBUG
   /* check all eumns of SCIP against the JNI enums */
   checkEnums();
#endif
   setlocale(LC_ALL,"C");

   JNISCIP_CALL( SCIPcreate(&scip) );

   return (jlong)(size_t)scip;
}

/** frees SCIP data structures */
JNIEXPORT
void JNISCIP(free)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPfree(&scip) );
   assert(scip == NULL);
}

/** returns current stage of SCIP */
JNIEXPORT
jint JNISCIP(getStage)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_STAGE stage;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   stage = SCIPgetStage(scip);

   /* check if the stage is covert in the JNI interface */
   assert(stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_INIT)
      || stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_PROBLEM)
      || stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_TRANSFORMING)
      || stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_TRANSFORMED)
      || stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_INITPRESOLVE)
      || stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_PRESOLVING)
      || stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_EXITPRESOLVE)
      || stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_PRESOLVED)
      || stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_INITSOLVE)
      || stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_SOLVING)
      || stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_SOLVED)
      || stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_EXITSOLVE)
      || stage == JNIPACKAGENAME(JniScipStage_SCIP_STAGE_FREETRANS));

   return (jint) stage;
}

/** outputs SCIP stage and solution status if applicable via the message handler
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
void JNISCIP(printStage)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jfile               /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;
   FILE* file;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   file = (FILE*) (size_t) jfile;

   JNISCIP_CALL( SCIPprintStage(scip, file) );
}

/** gets solution status */
JNIEXPORT
jint JNISCIP(getStatus)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_STATUS status;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   status = SCIPgetStatus(scip);

   /* check that the status is one which is captured in the JNI interface */
   assert(status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_UNKNOWN)
      || status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_USERINTERRUPT)
      || status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_NODELIMIT)
      || status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_STALLNODELIMIT)
      || status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_TIMELIMIT)
      || status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_MEMLIMIT)
      || status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_GAPLIMIT)
      || status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_SOLLIMIT)
      || status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_BESTSOLLIMIT)
      || status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_OPTIMAL)
      || status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_INFEASIBLE)
      || status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_UNBOUNDED)
      || status == JNIPACKAGENAME(JniScipStatus_SCIP_STATUS_INFORUNBD));

   return (jint) status;
}

/** outputs solution status
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  See \ref SCIP_Status "SCIP_STATUS" for a complete list of all possible solving status.
 */
JNIEXPORT
void JNISCIP(printStatus)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jfile               /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;
   FILE* file;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   file = (FILE*) (size_t) jfile;

   JNISCIP_CALL( SCIPprintStatus(scip, file) );
}

/** returns whether the current stage belongs to the transformed problem space */
JNIEXPORT
jboolean JNISCIP(isTransformed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Bool transformed;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   transformed = SCIPisTransformed(scip);

   return (jboolean)transformed;
}

/** returns whether the solution process should be probably correct
 *
 *  @note This feature is not supported yet!
 *
 *  @return Returns TRUE if \SCIP is exact solving mode, otherwise FALSE
 */
JNIEXPORT
jboolean JNISCIP(isExactSolve)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Bool exact;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   exact = SCIPisExactSolve(scip);

   return (jboolean)exact;
}

/** returns whether the presolving process would be finished given no more presolving reductions are found in this
 *  presolving round
 *
 *  Checks whether the number of presolving rounds is not exceeded and the presolving reductions found in the current
 *  presolving round suffice to trigger another presolving round.
 *
 *  @note if subsequent presolvers find more reductions, presolving might continue even if the method returns FALSE
 *  @note does not check whether infeasibility or unboundedness was already detected in presolving (which would result
 *        in presolving being stopped although the method returns TRUE)
 *
 *  @return Returns TRUE if presolving is finished if no further reductions are detected
 */
JNIEXPORT
jboolean JNISCIP(isPresolveFinished)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Bool prefin;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   prefin = SCIPisPresolveFinished(scip);

   return (jboolean)prefin;
}

/** returns whether the user pressed CTRL-C to interrupt the solving process */
JNIEXPORT
jboolean JNISCIP(pressedCtrlC)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Bool interrupt;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   interrupt = SCIPpressedCtrlC(scip);

   return (jboolean) interrupt;
}

/** returns whether the solving process should be / was stopped before proving optimality;
 *  if the solving process should be / was stopped, the status returned by {@link getStatus()} yields
 *  the reason for the premature abort
 */
JNIEXPORT
jboolean JNISCIP(isStopped)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Bool stopped;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   stopped = SCIPisStopped(scip);

   return (jboolean)stopped;
}

/** Installs the given message handler, such that all messages are passed to this handler. A messages handler can be
 *  created via SCIPmessagehdlrCreate().
 *
 *  @note The currently installed messages handler gets not freed. That has to be done by the user using
 *        SCIPmessagehdlrFree() or use SCIPsetMessagehdlrFree().
 */
JNIEXPORT
void JNISCIP(setMessagehdlr)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jmessagehdlr        /**< message handler to be installed,  or NULL to suppress all output */
   )
{
   SCIP* scip;
   SCIP_MESSAGEHDLR* messagehdlr;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   messagehdlr = (SCIP_MESSAGEHDLR*) (size_t) jmessagehdlr;

   JNISCIP_CALL( SCIPsetMessagehdlr(scip, messagehdlr) );
}

/** returns the currently installed message handler, or NULL if messages are currently suppressed */
JNIEXPORT
jlong JNISCIP(getMessagehdlr)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong)(size_t)SCIPgetMessagehdlr(scip);
}

/** sets the log file name for the currently installed message handler */
JNIEXPORT
void JNISCIP(setMessagehdlrLogfile)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of log file, or NULL (stdout) */
   )
{
   SCIP* scip;
   const char* filename;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   filename = (*env)->GetStringUTFChars(env, jname, &iscopy);
   assert(iscopy);

   SCIPsetMessagehdlrLogfile(scip, filename);
}

/** sets the currently installed message handler to be quiet (or not) */
JNIEXPORT
void JNISCIP(setMessagehdlrQuiet)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jboolean              jquiet              /**< should screen messages be suppressed? */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   SCIPsetMessagehdlrQuiet(scip, (SCIP_Bool)jquiet);
}


/** returns the current message verbosity level
 *
 *  @return message verbosity level of SCIP
 *
 *  @see \ref SCIP_VerbLevel "SCIP_VERBLEVEL" for a list of all verbosity levels
 */
jint JNISCIP(getVerbLevel)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_VERBLEVEL verblevel;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   verblevel = SCIPgetVerbLevel(scip);

   return (jint)verblevel;
}

/** copies plugins from sourcescip to targetscip; in case that a constraint handler which does not need constraints
 *  cannot be copied, valid will return FALSE. All plugins can declare that, if their copy process failed, the
 *  copied SCIP instance might not represent the same problem semantics as the original.
 *  Note that in this case dual reductions might be invalid.
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *        Also, 'passmessagehdlr' should be set to FALSE.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post After calling this method targetscip reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jboolean JNISCIP(copyPlugins)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsourcescip,        /**< source SCIP data structure */
   jlong                 jtargetscip,        /**< target SCIP data structure */
   jboolean              copyreaders,        /**< should the file readers be copied */
   jboolean              copypricers,        /**< should the variable pricers be copied */
   jboolean              copyconshdlrs,      /**< should the constraint handlers be copied */
   jboolean              copyconflicthdlrs,  /**< should the conflict handlers be copied */
   jboolean              copypresolvers,     /**< should the presolvers be copied */
   jboolean              copyrelaxators,     /**< should the relaxation handlers  be copied */
   jboolean              copyseparators,     /**< should the separators be copied */
   jboolean              copypropagators,    /**< should the propagators be copied */
   jboolean              copyheuristics,     /**< should the heuristics be copied */
   jboolean              copyeventhdlrs,     /**< should the event handlers be copied */
   jboolean              copynodeselectors,  /**< should the node selectors be copied */
   jboolean              copybranchrules,    /**< should the branchrules be copied */
   jboolean              copydisplays,       /**< should the display columns be copied */
   jboolean              copydialogs,        /**< should the dialogs be copied */
   jboolean              copynlpis,          /**< should the NLPIs be copied */
   jboolean              passmessagehdlr     /**< should the message handler be passed */
   )
{
   SCIP* sourcescip;
   SCIP* targetscip;
   SCIP_Bool valid;

   /* convert JNI pointer into C pointer */
   targetscip = (SCIP*) (size_t) jtargetscip;
   assert(targetscip != NULL);

   sourcescip = (SCIP*) (size_t) jsourcescip;
   assert(sourcescip != NULL);

   JNISCIP_CALL( SCIPcopyPlugins(sourcescip, targetscip, (SCIP_Bool)copyreaders, (SCIP_Bool)copypricers, (SCIP_Bool)copyconshdlrs, (SCIP_Bool)copyconflicthdlrs, (SCIP_Bool)copypresolvers, (SCIP_Bool)copyrelaxators, (SCIP_Bool)copyseparators, (SCIP_Bool)copypropagators, (SCIP_Bool)copyheuristics, (SCIP_Bool)copyeventhdlrs, (SCIP_Bool)copynodeselectors, (SCIP_Bool)copybranchrules, (SCIP_Bool)copydisplays, (SCIP_Bool)copydialogs, (SCIP_Bool)copynlpis, (SCIP_Bool)passmessagehdlr, &valid) );

   return (jboolean)valid;
}

/** create a problem by copying the problem data of the source SCIP
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post After calling this method targetscip reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
void JNISCIP(copyProb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsourcescip,        /**< source SCIP data structure */
   jlong                 jtargetscip,        /**< target SCIP data structure */
   jlong                 jvarmap,            /**< a hashmap to store the mapping of source variables corresponding
					      *   target variables, or NULL */
   jlong                 jconsmap,           /**< a hashmap to store the mapping of source constraints to the corresponding
					      *   target constraints, or NULL  */
   jboolean              global,             /**< create a global or a local copy? */
   jstring               jname               /**< problem name of target */
   )
{
   SCIP* sourcescip;
   SCIP* targetscip;
   SCIP_HASHMAP* varmap;
   SCIP_HASHMAP* consmap;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   sourcescip = (SCIP*) (size_t) jsourcescip;
   assert(sourcescip != NULL);

   targetscip = (SCIP*) (size_t) jtargetscip;
   assert(targetscip != NULL);

   varmap = (SCIP_HASHMAP*) (size_t) jvarmap;
   consmap = (SCIP_HASHMAP*) (size_t) jconsmap;

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPcopyProb(sourcescip, targetscip, varmap, consmap, (SCIP_Bool)global, name) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** copies all active variables from source-SCIP and adds these variable to the target-SCIP; the mapping between these
 *  variables are stored in the variable hashmap, target-SCIP has to be in problem creation stage, fixed and aggregated
 *  variables do not get copied
 *
 *  @note the variables are added to the target-SCIP but not captured
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
void JNISCIP(copyVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsourcescip,        /**< source SCIP data structure */
   jlong                 jtargetscip,        /**< target SCIP data structure */
   jlong                 jvarmap,            /**< a hashmap to store the mapping of source variables corresponding
					      *   target variables, or NULL */
   jlong                 jconsmap,           /**< a hashmap to store the mapping of source constraints to the corresponding
					      *   target constraints, or NULL  */
   jboolean              global              /**< create a global or a local copy? */
   )
{
   SCIP* sourcescip;
   SCIP* targetscip;
   SCIP_HASHMAP* varmap;
   SCIP_HASHMAP* consmap;

   /* convert JNI pointer into C pointer */
   sourcescip = (SCIP*) (size_t) jsourcescip;
   assert(sourcescip != NULL);

   targetscip = (SCIP*) (size_t) jtargetscip;
   assert(targetscip != NULL);

   varmap = (SCIP_HASHMAP*) (size_t) jvarmap;
   consmap = (SCIP_HASHMAP*) (size_t) jconsmap;

   JNISCIP_CALL( SCIPcopyVars(sourcescip, targetscip, varmap, consmap, (SCIP_Bool)global) );

}

/** copies constraints from the source-SCIP and adds these to the target-SCIP; for mapping the
 *  variables between the source and the target SCIP a hash map can be given; if the variable hash
 *  map is NULL or necessary variable mapping is missing, the required variables are created in the
 *  target-SCIP and added to the hash map, if not NULL; all variables which are created are added to
 *  the target-SCIP but not (user) captured; if the constraint hash map is not NULL the mapping
 *  between the constraints of the source and target-SCIP is stored
 *
 *  @note the constraints are added to the target-SCIP but are not (user) captured in the target SCIP. (If you mix
 *        SCIPgetConsCopy() with SCIPcopyConss() you should pay attention to what you add explicitly and what is already
 *        added.) You can check whether a constraint is added by calling SCIPconsIsAdded().
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jboolean JNISCIP(copyConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsourcescip,        /**< source SCIP data structure */
   jlong                 jtargetscip,        /**< target SCIP data structure */
   jlong                 jvarmap,            /**< a hashmap to store the mapping of source variables corresponding
					      *   target variables, or NULL */
   jlong                 jconsmap,           /**< a hashmap to store the mapping of source constraints to the corresponding
					      *   target constraints, or NULL  */
   jboolean              global,             /**< create a global or a local copy? */
   jboolean              enablepricing       /**< should pricing be enabled in copied SCIP instance?
					      *   If TRUE, the modifiable flag of constraints will be copied. */
   )
{
   SCIP* sourcescip;
   SCIP* targetscip;
   SCIP_HASHMAP* varmap;
   SCIP_HASHMAP* consmap;
   SCIP_Bool valid;

   /* convert JNI pointer into C pointer */
   sourcescip = (SCIP*) (size_t) jsourcescip;
   assert(sourcescip != NULL);

   targetscip = (SCIP*) (size_t) jtargetscip;
   assert(targetscip != NULL);

   varmap = (SCIP_HASHMAP*) (size_t) jvarmap;
   consmap = (SCIP_HASHMAP*) (size_t) jconsmap;

   JNISCIP_CALL( SCIPcopyConss(sourcescip, targetscip, varmap, consmap, (SCIP_Bool)global, (SCIP_Bool)enablepricing, &valid) );

   return (jboolean)valid;
}

/** convert all active cuts from cutpool to linear constraints
 *
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note SCIP stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jint JNISCIP(convertCutsToConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvarmap,            /**< a hashmap to store the mapping of source variables corresponding
					      *   target variables, or NULL */
   jlong                 jconsmap,           /**< a hashmap to store the mapping of source constraints to the corresponding
					      *   target constraints, or NULL  */
   jboolean              global              /**< create a global or a local copy? */
   )
{
   SCIP* scip;
   SCIP_HASHMAP* varmap;
   SCIP_HASHMAP* consmap;
   int ncutsadded;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   varmap = (SCIP_HASHMAP*) (size_t) jvarmap;
   consmap = (SCIP_HASHMAP*) (size_t) jconsmap;

   JNISCIP_CALL( SCIPconvertCutsToConss(scip, varmap, consmap, (SCIP_Bool)global, &ncutsadded) );

   return (jint)ncutsadded;
}

/** copies all active cuts from cutpool of sourcescip to linear constraints in targetscip
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jint JNISCIP(copyCuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsourcescip,        /**< source SCIP data structure */
   jlong                 jtargetscip,        /**< target SCIP data structure */
   jlong                 jvarmap,            /**< a hashmap to store the mapping of source variables corresponding
					      *   target variables, or NULL */
   jlong                 jconsmap,           /**< a hashmap to store the mapping of source constraints to the corresponding
					      *   target constraints, or NULL  */
   jboolean              global              /**< create a global or a local copy? */
   )
{
   SCIP* sourcescip;
   SCIP* targetscip;
   SCIP_HASHMAP* varmap;
   SCIP_HASHMAP* consmap;
   int ncutsadded;

   /* convert JNI pointer into C pointer */
   sourcescip = (SCIP*) (size_t) jsourcescip;
   assert(sourcescip != NULL);

   targetscip = (SCIP*) (size_t) jtargetscip;
   assert(targetscip != NULL);

   varmap = (SCIP_HASHMAP*) (size_t) jvarmap;
   consmap = (SCIP_HASHMAP*) (size_t) jconsmap;

   JNISCIP_CALL( SCIPcopyCuts(sourcescip, targetscip, varmap, consmap, (SCIP_Bool)global, &ncutsadded) );

   return (jint)ncutsadded;
}

/** copies parameter settings from sourcescip to targetscip
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_FREE
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
void JNISCIP(copyParamSettings)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsourcescip,        /**< source SCIP data structure */
   jlong                 jtargetscip         /**< target SCIP data structure */

   )
{
   SCIP* sourcescip;
   SCIP* targetscip;

   /* convert JNI pointer into C pointer */
   sourcescip = (SCIP*) (size_t) jsourcescip;
   assert(sourcescip != NULL);

   targetscip = (SCIP*) (size_t) jtargetscip;
   assert(targetscip != NULL);

   JNISCIP_CALL( SCIPcopyParamSettings(sourcescip, targetscip) );
}

/** gets depth of current scip instance (increased by each copy call)
 *
 *  @return Depth of subscip of SCIP is returned.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note SCIP stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jint JNISCIP(getSubscipDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */

   )
{
   SCIP* scip;
   int subscipdepth;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   subscipdepth = SCIPgetSubscipDepth(scip);

   return (jint)subscipdepth;
}

/** copies source SCIP to target SCIP; the copying process is done in the following order:
 *  1) copy the plugins
 *  2) copy the settings
 *  3) create problem data in target-SCIP and copy the problem data of the source-SCIP
 *  4) copy all active variables
 *  5) copy all constraints
 *
 *  @note all variables and constraints which are created in the target-SCIP are not (user) captured
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *        Also, 'passmessagehdlr' should be set to FALSE.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jboolean JNISCIP(copy)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jsourcescip,        /**< source SCIP data structure */
   jlong                 jtargetscip,        /**< target SCIP data structure */
   jlong                 jvarmap,            /**< a hashmap to store the mapping of source variables corresponding
					      *   target variables, or NULL */
   jlong                 jconsmap,           /**< a hashmap to store the mapping of source constraints to the corresponding
					      *   target constraints, or NULL  */
   jstring               jsuffix,            /**< suffix which will be added to the names of the target SCIP, might be empty */
   jboolean              global,             /**< create a global or a local copy? */
   jboolean              enablepricing,      /**< should pricing be enabled in copied SCIP instance? If TRUE, pricer
					      *   plugins will be copied and activated, and the modifiable flag of
					      *   constraints will be respected. If FALSE, valid will be set to FALSE, when
					      *   there are pricers present */
   jboolean              passmessagehdlr     /**< should the message handler be passed */
   )
{
   SCIP* sourcescip;
   SCIP* targetscip;
   SCIP_HASHMAP* varmap;
   SCIP_HASHMAP* consmap;
   const char* suffix;
   jboolean iscopy;
   SCIP_Bool valid;

   /* convert JNI pointer into C pointer */
   sourcescip = (SCIP*) (size_t) jsourcescip;
   assert(sourcescip != NULL);

   targetscip = (SCIP*) (size_t) jtargetscip;
   assert(targetscip != NULL);

   varmap = (SCIP_HASHMAP*) (size_t) jvarmap;
   consmap = (SCIP_HASHMAP*) (size_t) jconsmap;

   /* convert JNI string into const char* */
   suffix = (*env)->GetStringUTFChars(env, jsuffix, &iscopy);
   if( suffix == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPcopy(sourcescip, targetscip, varmap, consmap, suffix, (SCIP_Bool)global, (SCIP_Bool)enablepricing, (SCIP_Bool)passmessagehdlr, &valid) );

   return (jboolean)valid;
}

/** gets the fixing status of an existing parameter
 *
 *  @return TRUE if the parameter is fixed to a value, otherwise FALSE.
 */
JNIEXPORT
jboolean JNISCIP(isParamFixed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_Bool paramfixed;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPgetBoolParam(scip, name, &paramfixed) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jboolean) paramfixed;
}

/** returns the pointer to the SCIP parameter with the given name
 *
 *  @return pointer to the parameter with the given name
 */
JNIEXPORT
jlong JNISCIP(getParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of the parameter */
   )
{
	SCIP* scip;
	const char* name;
	jboolean iscopy;
	SCIP_PARAM *param;

	/* convert JNI pointer into C pointer */
	scip = (SCIP*) (size_t) jscip;
	assert(scip != NULL);

	/* convert JNI string into const char* */
	name = (*env)->GetStringUTFChars(env, jname, &iscopy);
	if( name == NULL )
	   SCIPABORT();
	assert(iscopy);

	param=SCIPgetParam(scip, name);

	(*env)->ReleaseStringUTFChars(env, jname, name);

	return (jlong) (size_t) param;
}


/** gets the value of an existing SCIP_Bool parameter */
JNIEXPORT
jboolean JNISCIP(getBoolParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_Bool value;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPgetBoolParam(scip, name, &value) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jboolean) value;
}

/** gets the value of an existing Int parameter */
JNIEXPORT
jint JNISCIP(getIntParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   int value;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPgetIntParam(scip, name, &value) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jint) value;
}

/** gets the value of an existing SCIP_Longint parameter */
JNIEXPORT
jlong JNISCIP(getLongintParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_Longint value;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPgetLongintParam(scip, name, &value) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) value;
}

/** gets the value of an existing SCIP_Real parameter */
JNIEXPORT
jdouble JNISCIP(getRealParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_Real value;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPgetRealParam(scip, name, &value) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jdouble) value;
}

/** gets the value of an existing Char parameter */
JNIEXPORT
jchar JNISCIP(getCharParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   char value;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPgetCharParam(scip, name, &value) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jchar) value;
}

/** gets the value of an existing String parameter */
JNIEXPORT
jstring JNISCIP(getStringParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   char* value;
   jstring jval;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPgetStringParam(scip, name, &value) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   /* convert char* into jstring */
   jval = (*env)->NewStringUTF(env, value);

   return jval;
}

/** fixes the value of an existing parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @note: Be careful with this method! Some general settings, e.g., the time or node limit, should not be fixed because
 *         they have to be changed for sub-SCIPs.
 */
JNIEXPORT
void JNISCIP(fixParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPfixParam(scip, name) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** unfixes the value of an existing parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(unfixParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPunfixParam(scip, name) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** changes the value of an existing parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(setParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of the parameter */
   jlong                 jval                /**< new value of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   void* value;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   value = (void*) (size_t) jval;

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPsetParam(scip, name, value) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** changes the value of an existing SCIP_Bool parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(chgBoolParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jparam,             /**< name of the parameter */
   jboolean              jval                /**< new value of the parameter */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPchgBoolParam(scip, (SCIP_PARAM*) (size_t)jparam, (SCIP_Bool)jval) );
}

/** changes the value of an existing SCIP_Bool parameter */
JNIEXPORT
void JNISCIP(setBoolParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of the parameter */
   jboolean              jval                /**< new value of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPsetBoolParam(scip, name, (SCIP_Bool) jval) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** changes the value of an existing int parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(chgIntParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jparam,             /**< name of the parameter */
   jint                  jval                /**< new value of the parameter */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPchgIntParam(scip, (SCIP_PARAM*) (size_t)jparam, (int)jval) );
}

/** changes the value of an existing Int parameter */
JNIEXPORT
void JNISCIP(setIntParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of the parameter */
   jint                  jval                /**< new value of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPsetIntParam(scip, name, (jint) jval) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** changes the value of an existing SCIP_Longint parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(chgLongintParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jparam,             /**< name of the parameter */
   jlong                 jval                /**< new value of the parameter */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPchgLongintParam(scip, (SCIP_PARAM*) (size_t)jparam, (SCIP_Longint)jval) );
}

/** changes the value of an existing SCIP_Longint parameter */
JNIEXPORT
void JNISCIP(setLongintParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of the parameter */
   jlong                 jval                /**< new value of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPsetLongintParam(scip, name, (SCIP_Longint) jval) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** changes the value of an existing SCIP_Real parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(chgRealParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jparam,             /**< name of the parameter */
   jdouble               jval                /**< new value of the parameter */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPchgRealParam(scip, (SCIP_PARAM*) (size_t)jparam, (SCIP_Real)jval) );
}

/** changes the value of an existing SCIP_Real parameter */
JNIEXPORT
void JNISCIP(setRealParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of the parameter */
   jdouble               jval                /**< new value of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPsetRealParam(scip, name, (SCIP_Real) jval) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** changes the value of an existing char parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(chgCharParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jparam,             /**< name of the parameter */
   jchar                 jval                /**< new value of the parameter */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPchgCharParam(scip, (SCIP_PARAM*) (size_t)jparam, (char)jval) );
}

/** changes the value of an existing Char parameter */
JNIEXPORT
void JNISCIP(setCharParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of the parameter */
   jchar                 jval                /**< new value of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPsetCharParam(scip, name, (char) jval) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** changes the value of an existing string(char*) parameter
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(chgStringParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jparam,             /**< name of the parameter */
   jstring               jval                /**< new value of the parameter */
   )
{
   SCIP* scip;
   const char* value;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   value = (*env)->GetStringUTFChars(env, jval, &iscopy);
   if( value == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPchgStringParam(scip, (SCIP_PARAM*) (size_t)jparam, value) );

   (*env)->ReleaseStringUTFChars(env, jval, value);
}

/** changes the value of an existing String parameter */
JNIEXPORT
void JNISCIP(setStringParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of the parameter */
   jstring               jval                /**< new value of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   const char* value;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   value = (*env)->GetStringUTFChars(env, jval, &iscopy);
   if( value == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPsetStringParam(scip, name, value) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
   (*env)->ReleaseStringUTFChars(env, jval, value);
}

/** reads parameters from a file */
JNIEXPORT
void JNISCIP(readParams)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jfilename           /**< filename */
   )
{
   SCIP* scip;
   const char* filename;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   filename = (*env)->GetStringUTFChars(env, jfilename, &iscopy);
   if( filename == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPreadParams(scip, filename) );

   (*env)->ReleaseStringUTFChars(env, jfilename, filename);
}

/** writes all parameters in the parameter set to a file */
JNIEXPORT
void JNISCIP(writeParams)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jfilename,          /**< filename or NULL */
   jboolean              jcomments,          /**< should parameter descriptions be written as comments? */
   jboolean              jonlychanged        /**< should only the parameters been written, that are changed from default? */
   )
{
   SCIP* scip;
   const char* filename;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   filename = (*env)->GetStringUTFChars(env, jfilename, &iscopy);
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPwriteParams(scip, filename, (SCIP_Bool) jcomments, (SCIP_Bool) jonlychanged) );

   (*env)->ReleaseStringUTFChars(env, jfilename, filename);
}

/** resets a single parameter to its default value */
JNIEXPORT
void JNISCIP(resetParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of the parameter */
   )
{
   SCIP* scip;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPresetParam(scip, name) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** Resets all parameters to their default values */
JNIEXPORT
void JNISCIP(resetParams)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPresetParams(scip) );
}

/** sets parameters */
JNIEXPORT
void JNISCIP(setEmphasis)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  jparamemphasis,     /**<  emphasis setting */
   jboolean              jquiet              /**< setting parameters quietly */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* check that the emphasis is one which is captured in the JNI interface */
   assert(jparamemphasis == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_COUNTER)
      || jparamemphasis == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_CPSOLVER)
      || jparamemphasis == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_EASYCIP)
      || jparamemphasis == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_FEASIBILITY)
      || jparamemphasis == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_HARDLP)
      || jparamemphasis == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_OPTIMALITY)
      || jparamemphasis == JNIPACKAGENAME(JniScipParamemphasis_SCIP_PARAMEMPHASIS_DEFAULT));

   JNISCIP_CALL( SCIPsetEmphasis(scip, (int)jparamemphasis, (SCIP_Bool)jquiet) );
}

/** sets parameters to deactivate separators and heuristics that use auxiliary SCIP instances; should be called for
 *  auxiliary SCIP instances to avoid recursion
 */
JNIEXPORT
void JNISCIP(setSubscipsOff)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jboolean              jquiet              /**< setting parameters quietly */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetSubscipsOff(scip, (SCIP_Bool)jquiet) );
}

/** sets heuristic parameters values to
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all heuristic parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for heuristic is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the heuristic are called more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all heuristics
 */
JNIEXPORT
void JNISCIP(setHeuristics)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  jparamsetting,      /**< setting */
   jboolean              jquiet              /**< setting parameters quietly */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* check that the setting is one which is captured in the JNI interface */
   assert(jparamsetting == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_DEFAULT)
      || jparamsetting == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_AGGRESSIVE)
      || jparamsetting == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_FAST)
      || jparamsetting == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_OFF));

   JNISCIP_CALL( SCIPsetHeuristics(scip, (int)jparamsetting, (SCIP_Bool)jquiet) );
}

/** sets presolving parameters values to
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all presolving parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for presolving is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the presolving are called more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all presolving
 */
JNIEXPORT
void JNISCIP(setPresolving)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  jparamsetting,      /**< setting */
   jboolean              jquiet              /**< setting parameters quietly */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* check that the setting is one which is captured in the JNI interface */
   assert(jparamsetting == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_DEFAULT)
      || jparamsetting == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_AGGRESSIVE)
      || jparamsetting == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_FAST)
      || jparamsetting == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_OFF));

   JNISCIP_CALL( SCIPsetPresolving(scip, (int)jparamsetting, (SCIP_Bool)jquiet) );
}

/** sets separating parameters values to
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all separating parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for separating is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the separating are called more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all separating
 */
JNIEXPORT
void JNISCIP(setSeparating)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  jparamsetting,      /**< setting */
   jboolean              jquiet              /**< setting parameters quietly */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* check that the setting is one which is captured in the JNI interface */
   assert(jparamsetting == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_DEFAULT)
      || jparamsetting == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_AGGRESSIVE)
      || jparamsetting == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_FAST)
      || jparamsetting == JNIPACKAGENAME(JniScipParamsetting_SCIP_PARAMSETTING_OFF));

   JNISCIP_CALL( SCIPsetSeparating(scip, (int)jparamsetting, (SCIP_Bool)jquiet) );
}

/** returns the array of all available SCIP parameters
 *
 *  @return SCIP_PARAM* array, containing all SCIP parameters.
 */
JNIEXPORT
jlongArray JNISCIP(getParams)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nparams;

   jlongArray jparams;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nparams = SCIPgetNParams(scip);
   jparams = (*env)->NewLongArray(env, nparams);

   if( jparams == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_PARAM** params;

      params = SCIPgetParams(scip);
      (*env)->SetLongArrayRegion(env, jparams, 0, nparams, (jlong*)params);
   }

   return jparams;
}

/** returns the total number of all available SCIP parameters
 *
 *  @return number of all SCIP parameters.
 */
jint JNISCIP(getNParams)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNParams(scip);

   return (jint) num;
}

/** creates a reader and includes it in SCIP. All non-fundamental (or optional) callbacks will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see
 *  SCIPsetReaderCopy(), SCIPsetReaderFree(), SCIPsetReaderRead(), SCIPsetReaderWrite().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeReader() instead
 */
void JNISCIP(includeReaderBasic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlongArray            jreaderptr,         /**< reference to reader pointer, or NULL */
   jstring               jname,              /**< name of reader */
   jstring               jdesc,              /**< description of reader */
   jstring               jextension,         /**< file extension that reader processes */
   jlong                 jreaderdata         /**< reader data */
   )
{
   SCIPerrorMessage("method includeReaderBasic is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );
}

/** returns the reader of the given name, or NULL if not existing */
jlong JNISCIP(findReader)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_READER* reader;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   reader = SCIPfindReader(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) reader;
}

/** returns the array of currently available readers */
jlongArray JNISCIP(getReaders)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nreaders;

   jlongArray jreaders;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nreaders = SCIPgetNReaders(scip);
   jreaders = (*env)->NewLongArray(env, nreaders);

   if( jreaders == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_READER** readers;

      readers = SCIPgetReaders(scip);
      (*env)->SetLongArrayRegion(env, jreaders, 0, nreaders, (jlong*)readers);
   }

   return jreaders;
}

/** returns the number of currently available readers */
jint JNISCIP(getNReaders)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNReaders(scip);

   return (jint) num;
}

/** returns the variable pricer of the given name, or NULL if not existing */
jlong JNISCIP(findPricer)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_PRICER* pricer;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   pricer = SCIPfindPricer(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) pricer;
}

/** returns the array of currently available variable pricers; active pricers are in the first slots of the array */
jlongArray JNISCIP(getPricers)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int npricers;

   jlongArray jpricers;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   npricers = SCIPgetNPricers(scip);
   jpricers = (*env)->NewLongArray(env, npricers);

   if( jpricers == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_PRICER** pricers;

      pricers = SCIPgetPricers(scip);
      (*env)->SetLongArrayRegion(env, jpricers, 0, npricers, (jlong*)pricers);
   }

   return jpricers;
}

/** returns the number of currently available variable pricers */
jint JNISCIP(getNPricers)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNPricers(scip);

   return (jint) num;
}

/** returns the number of currently active variable pricers, that are used in the LP solving loop */
JNIEXPORT
jint JNISCIP(getNActivePricers)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNActivePricers(scip);

   return (jint) num;
}

/** sets the priority priority of a variable pricer */
void JNISCIP(setPricerPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jpricer,            /**< variable pricer */
   jint                  jpriority           /**< new priority of the variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetPricerPriority(scip, (SCIP_PRICER*) (size_t)jpricer, (int)jpriority) );
}

/** activates pricer to be used for the current problem
 *  This method should be called during the problem creation stage for all pricers that are necessary to solve
 *  the problem model.
 *  The pricers are automatically deactivated when the problem is freed.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
void JNISCIP(activatePricer)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jpricer             /**< variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPactivatePricer(scip, (SCIP_PRICER*) (size_t)jpricer) );
}

/** deactivates pricer
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
void JNISCIP(deactivatePricer)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jpricer             /**< variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPdeactivatePricer(scip, (SCIP_PRICER*) (size_t)jpricer) );
}

/** returns the constraint handler of the given name, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findConshdlr)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_CONSHDLR* conshdlr;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   conshdlr = SCIPfindConshdlr(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) conshdlr;
}

/** returns the array of currently available constraint handlers */
JNIEXPORT
jlongArray JNISCIP(getConshdlrs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nconshdlrs;

   jlongArray jconshdlrs;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nconshdlrs = SCIPgetNConshdlrs(scip);
   jconshdlrs = (*env)->NewLongArray(env, nconshdlrs);

   if( jconshdlrs == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_CONSHDLR** conshdlrs;

      conshdlrs = SCIPgetConshdlrs(scip);
      (*env)->SetLongArrayRegion(env, jconshdlrs, 0, nconshdlrs, (jlong*)conshdlrs);
   }

   return jconshdlrs;
}

/** returns the number of currently available constraint handlers */
jint JNISCIP(getNConshdlrs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNConshdlrs(scip);

   return (jint) num;
}

/** returns the conflict handler of the given name, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findConflicthdlr)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_CONFLICTHDLR* confhdlr;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   confhdlr = SCIPfindConflicthdlr(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) confhdlr;
}

/** returns the array of currently available conflict handlers */
JNIEXPORT
jlongArray JNISCIP(getConflicthdlrs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nconflicthdlrs;

   jlongArray jconflicthdlrs;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nconflicthdlrs = SCIPgetNConflicthdlrs(scip);
   jconflicthdlrs = (*env)->NewLongArray(env, nconflicthdlrs);

   if (jconflicthdlrs == NULL)
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_CONFLICTHDLR** conflicthdlrs;

      conflicthdlrs = SCIPgetConflicthdlrs(scip);
      (*env)->SetLongArrayRegion(env, jconflicthdlrs, 0, nconflicthdlrs, (jlong*)conflicthdlrs);
   }

   return jconflicthdlrs;
}

/** returns the number of currently available conflict handlers */
JNIEXPORT
jint JNISCIP(getNConflicthdlrs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNConflicthdlrs(scip);

   return (jint) num;
}

/** sets the priority of a conflict handler */
JNIEXPORT
void JNISCIP(setConflicthdlrPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jconflicthdlr,      /**< variable pricer */
   jint                  jpriority           /**< new priority of the variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetConflicthdlrPriority(scip, (SCIP_CONFLICTHDLR*) (size_t)jconflicthdlr, (int)jpriority) );
}

/** returns the presolver of the given name, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findPresol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_PRESOL* presol;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   presol = SCIPfindPresol(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) presol;
}

/** returns the array of currently available presolvers */
JNIEXPORT
jlongArray JNISCIP(getPresols)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int npresols;

   jlongArray jpresols;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   npresols = SCIPgetNPresols(scip);
   jpresols = (*env)->NewLongArray(env, npresols);

   if( jpresols == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_PRESOL** presols;

      presols = SCIPgetPresols(scip);
      (*env)->SetLongArrayRegion(env, jpresols, 0, npresols, (jlong*)presols);
   }

   return jpresols;
}

/** returns the number of currently available presolvers */
JNIEXPORT
jint JNISCIP(getNPresols)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNPresols(scip);

   return (jint) num;
}

/** sets the priority of a presolver */
JNIEXPORT
void JNISCIP(setPresolPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jpresol,            /**< variable pricer */
   jint                  jpriority           /**< new priority of the variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetPresolPriority(scip, (SCIP_PRESOL*) (size_t)jpresol, (int)jpriority) );
}

/** returns the relaxation handler of the given name, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findRelax)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_RELAX* relax;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   relax = SCIPfindRelax(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) relax;
}

/** returns the array of currently available relaxation handlers  */
JNIEXPORT
jlongArray JNISCIP(getRelaxs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nrelaxs;

   jlongArray jrelaxs;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nrelaxs = SCIPgetNRelaxs(scip);
   jrelaxs = (*env)->NewLongArray(env, nrelaxs);

   if( jrelaxs == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_RELAX** relaxs;

      relaxs = SCIPgetRelaxs(scip);
      (*env)->SetLongArrayRegion(env, jrelaxs, nrelaxs, nrelaxs, (jlong*)relaxs);
   }

   return jrelaxs;
}

/** returns the number of currently available relaxation handlers  */
JNIEXPORT
jint JNISCIP(getNRelaxs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNRelaxs(scip);

   return (jint) num;
}

/** sets the priority of a relaxation handler */
JNIEXPORT
void JNISCIP(setRelaxPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrelax,             /**< variable pricer */
   jint                  jpriority           /**< new priority of the variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetRelaxPriority(scip, (SCIP_RELAX*) (size_t)jrelax, (int)jpriority) );
}

/** returns the separator of the given name, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findSepa)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_SEPA* sepa;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   sepa = SCIPfindSepa(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) sepa;
}

/** returns the array of currently available separators */
JNIEXPORT
jlongArray JNISCIP(getSepas)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nsepas;

   jlongArray jsepas;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nsepas = SCIPgetNSepas(scip);
   jsepas = (*env)->NewLongArray(env, nsepas);

   if( jsepas == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_SEPA** sepas;

      sepas = SCIPgetSepas(scip);
      (*env)->SetLongArrayRegion(env, jsepas, 0, nsepas, (jlong*)sepas);
   }

   return jsepas;
}

/** returns the number of currently available separators */
JNIEXPORT
jint JNISCIP(getNSepas)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNSepas(scip);

   return (jint) num;
}

/** sets the priority of a separator */
JNIEXPORT
void JNISCIP(setSepaPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsepa,              /**< variable pricer */
   jint                  jpriority           /**< new priority of the variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetSepaPriority(scip, (SCIP_SEPA*) (size_t)jsepa, (int)jpriority) );
}

/** returns the propagator of the given name, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findProp)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_PROP* prop;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   prop = SCIPfindProp(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) prop;
}
/** returns the array of currently available propagators */
JNIEXPORT
jlongArray JNISCIP(getProps)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nprops;

   jlongArray jprops;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nprops = SCIPgetNProps(scip);
   jprops = (*env)->NewLongArray(env, nprops);

   if( jprops == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_PROP** props;

      props = SCIPgetProps(scip);
      (*env)->SetLongArrayRegion(env, jprops, 0, nprops, (jlong*)props);
   }

   return jprops;
}

/** returns the number of currently available propagators */
JNIEXPORT
jint JNISCIP(getNProps)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNProps(scip);

   return (jint) num;
}

/** sets the priority of a propagator */
JNIEXPORT
void JNISCIP(setPropPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jprop,              /**< variable pricer */
   jint                  jpriority           /**< new priority of the variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetPropPriority(scip, (SCIP_PROP*) (size_t)jprop, (int)jpriority) );
}

/** sets the presolving priority of a propagator */
JNIEXPORT
void JNISCIP(setPropPresolPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jprop,              /**< variable pricer */
   jint                  jpriority           /**< new priority of the variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetPropPresolPriority(scip, (SCIP_PROP*) (size_t)jprop, (int)jpriority) );
}

/** returns the primal heuristic of the given name, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findHeur)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_HEUR* heur;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   heur = SCIPfindHeur(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) heur;
}

/** returns the array of currently available primal heuristics */
JNIEXPORT
jlongArray JNISCIP(getHeurs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nheurs;

   jlongArray jheurs;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nheurs = SCIPgetNHeurs(scip);
   jheurs = (*env)->NewLongArray(env, nheurs);

   if( jheurs == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_HEUR** heurs;

      heurs = SCIPgetHeurs(scip);
      (*env)->SetLongArrayRegion(env, jheurs, 0, nheurs, (jlong*)heurs);
   }

   return jheurs;
}

/** returns the number of currently available primal heuristics */
JNIEXPORT
jint JNISCIP(getNHeurs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNHeurs(scip);

   return (jint) num;
}

/** sets the priority of a primal heuristic */
JNIEXPORT
void JNISCIP(setHeurPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jheur,              /**< variable pricer */
   jint                  jpriority           /**< new priority of the variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetHeurPriority(scip, (SCIP_HEUR*) (size_t)jheur, (int)jpriority) );
}

/** returns the event handler of the given name, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findEventhdlr)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_EVENTHDLR* eventhdlr;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   eventhdlr = SCIPfindEventhdlr(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) eventhdlr;
}

/** returns the array of currently available event handlers */
JNIEXPORT
jlongArray JNISCIP(getEventhdlrs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int neventhdlrs;

   jlongArray jeventhdlrs;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   neventhdlrs = SCIPgetNEventhdlrs(scip);
   jeventhdlrs = (*env)->NewLongArray(env, neventhdlrs);

   if( jeventhdlrs == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_EVENTHDLR** eventhdlrs;

      eventhdlrs = SCIPgetEventhdlrs(scip);
      (*env)->SetLongArrayRegion(env, jeventhdlrs, 0, neventhdlrs, (jlong*)eventhdlrs);
   }

   return jeventhdlrs;
}

/** returns the number of currently available event handlers */
JNIEXPORT
jint JNISCIP(getNEventhdlrs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNEventhdlrs(scip);

   return (jint) num;
}

/** returns the node selector of the given name, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findNodesel)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_NODESEL* nodesel;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   nodesel = SCIPfindNodesel(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) nodesel;
}

/** returns the array of currently available node selectors */
JNIEXPORT
jlongArray JNISCIP(getNodesels)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nnodesels;

   jlongArray jnodesels;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nnodesels = SCIPgetNNodesels(scip);
   jnodesels = (*env)->NewLongArray(env, nnodesels);

   if( jnodesels == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_NODESEL** nodesels;

      nodesels = SCIPgetNodesels(scip);
      (*env)->SetLongArrayRegion(env, jnodesels, nnodesels, nnodesels, (jlong*)nodesels);
   }

   return jnodesels;
}

/** returns the number of currently available node selectors */
JNIEXPORT
jint JNISCIP(getNNodesels)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNNodesels(scip);

   return (jint) num;
}

/** sets the priority of a node selector in standard mode */
JNIEXPORT
void JNISCIP(setNodeselStdPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnodesel,            /**< variable pricer */
   jint                  jpriority           /**< new priority of the variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetNodeselStdPriority(scip, (SCIP_NODESEL*) (size_t)jnodesel, (int)jpriority) );
}

/** sets the priority of a node selector in memory saving mode */
JNIEXPORT
void JNISCIP(setNodeselMemsavePriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnodesel,            /**< variable pricer */
   jint                  jpriority           /**< new priority of the variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetNodeselMemsavePriority(scip, (SCIP_NODESEL*) (size_t)jnodesel, (int)jpriority) );
}

/** returns the currently used node selector */
JNIEXPORT
jlong JNISCIP(getNodesel)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_NODESEL* nodesel;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nodesel = SCIPgetNodesel(scip);

   return (jlong) (size_t) nodesel;
}

/** creates a branching rule and includes it in SCIP. All non-fundamental (or optional) callbacks will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetBranchruleInit(), SCIPsetBranchruleExit(),
 *  SCIPsetBranchruleCopy(), SCIPsetBranchruleFree(), SCIPsetBranchruleInitsol(), SCIPsetBranchruleExitsol(),
 *  SCIPsetBranchruleExecLp(), SCIPsetBranchruleExecExt(), and SCIPsetBranchruleExecPs().
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeBranchrule() instead
 */
JNIEXPORT
jlong JNISCIP(includeBranchruleBasic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of branching rule */
   jstring               jdesc,              /**< description of branching rule */
   jint                  priority,           /**< priority of the branching rule */
   jint                  maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   jdouble               maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound
					      *   compared to best node's dual bound for applying branching rule
					      *   (0.0: only on current best node, 1.0: on all nodes) */
   jlong                 jbranchruledata     /**< branching rule data */
   )
{
   SCIP* scip;
   const char* name;
   const char* desc;
   jboolean iscopy1;
   jboolean iscopy2;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchruleptr;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   branchruledata = (SCIP_BRANCHRULEDATA*) (size_t) jbranchruledata;
   assert(branchruledata != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy1);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy1);

   desc = (*env)->GetStringUTFChars(env, jdesc, &iscopy2);
   if( desc == NULL )
      SCIPABORT();
   assert(iscopy2);

   JNISCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchruleptr, name, desc, (int)priority, (int)maxdepth, (SCIP_Real)maxbounddist, branchruledata) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
   (*env)->ReleaseStringUTFChars(env, jdesc, desc);

   return (jlong) (size_t) branchruleptr;
}

/** returns the branching rule of the given name, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findBranchrule)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_BRANCHRULE* branchrule;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   branchrule = SCIPfindBranchrule(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) branchrule;
}

/** returns the array of currently available branching rules */
JNIEXPORT
jlongArray JNISCIP(getBranchrules)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nbranchrules;

   jlongArray jbranchrules;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nbranchrules = SCIPgetNBranchrules(scip);
   jbranchrules = (*env)->NewLongArray(env, nbranchrules);

   if( jbranchrules == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_BRANCHRULE** branchrules;

      branchrules = SCIPgetBranchrules(scip);
      (*env)->SetLongArrayRegion(env, jbranchrules, 0, nbranchrules, (jlong*)branchrules);
   }

   return jbranchrules;
}

/** returns the number of currently available branching rules */
JNIEXPORT
jint JNISCIP(getNBranchrules)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNBranchrules(scip);

   return (jint) num;
}

/** sets the priority of a branching rule */
JNIEXPORT
void JNISCIP(setBranchrulePriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jbranchrule,            /**< variable pricer */
   jint                  jpriority           /**< new priority of the variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetBranchrulePriority(scip, (SCIP_BRANCHRULE*) (size_t)jbranchrule, (int)jpriority) );
}

/** sets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
JNIEXPORT
void JNISCIP(setBranchruleMaxdepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jbranchrule,        /**< variable pricer */
   jint                  jmaxdepth           /**< new maxdepth of the branching rule */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetBranchruleMaxdepth(scip, (SCIP_BRANCHRULE*) (size_t)jbranchrule, (int)jmaxdepth) );
}

/** sets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
JNIEXPORT
void JNISCIP(setBranchruleMaxbounddist)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jbranchrule,        /**< variable pricer */
   jdouble               jmaxbounddist       /**< new maxdepth of the branching rule */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetBranchruleMaxbounddist(scip, (SCIP_BRANCHRULE*) (size_t)jbranchrule, (SCIP_Real)jmaxbounddist) );
}

/** returns the display column of the given name, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findDisp)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_DISP* disp;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   disp = SCIPfindDisp(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) disp;
}

/** returns the array of currently available display columns */
JNIEXPORT
jlongArray JNISCIP(getDisps)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int ndisps;

   jlongArray jdisps;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   ndisps = SCIPgetNDisps(scip);
   jdisps = (*env)->NewLongArray(env, ndisps);

   if( jdisps == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_DISP** disps;

      disps = SCIPgetDisps(scip);
      (*env)->SetLongArrayRegion(env, jdisps, 0, ndisps, (jlong*)disps);
   }

   return jdisps;
}

/** returns the number of currently available display columns */
JNIEXPORT
jint JNISCIP(getNDisps)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNDisps(scip);

   return (jint) num;
}

/** automatically selects display columns for being shown w.r.t. the display width parameter */
JNIEXPORT
void JNISCIP(autoselectDisps)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPautoselectDisps(scip) );
}

/** includes an NLPI in SCIP */
JNIEXPORT
void JNISCIP(includeNlpi)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlpi               /**< NLPI data structure */
   )
{
   SCIP* scip;
   SCIP_NLPI* nlpi;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlpi = (SCIP_NLPI*) (size_t) jnlpi;
   assert(nlpi != NULL);

   JNISCIP_CALL( SCIPincludeNlpi(scip, nlpi) );
}

/** returns the NLPI of the given name, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findNlpi)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint handler */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_NLPI* nlpi;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   nlpi = SCIPfindNlpi(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) (size_t) nlpi;
}

/** returns the array of currently available NLPIs (sorted by priority) */
JNIEXPORT
jlongArray JNISCIP(getNlpis)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nnlpis;

   jlongArray jnlpis;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nnlpis = SCIPgetNNlpis(scip);
   jnlpis = (*env)->NewLongArray(env, nnlpis);

   if( jnlpis == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_NLPI** nlpis;

      nlpis = SCIPgetNlpis(scip);
      (*env)->SetLongArrayRegion(env, jnlpis, 0, nnlpis, (jlong*)nlpis);
   }

   return jnlpis;
}

/** returns the number of currently available NLPIs */
JNIEXPORT
jint JNISCIP(getNNlpis)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNNlpis(scip);

   return (jint) num;
}

/** sets the priority of an NLPI */
JNIEXPORT
void JNISCIP(setNlpiPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlpi,              /**< variable pricer */
   jint                  jpriority           /**< new priority of the variable pricer */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetNlpiPriority(scip, (SCIP_NLPI*) (size_t)jnlpi, (int)jpriority) );
}

/** includes information about an external code linked into the SCIP library */
JNIEXPORT
void JNISCIP(includeExternalCodeInformation)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of external code */
   jstring               jdesc               /**< description of external code, or NULL */

   )
{
   SCIP* scip;
   const char* name;
   const char* desc;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();
   assert(iscopy);

   /* convert JNI string into const char* */
   desc = (*env)->GetStringUTFChars(env, jdesc, &iscopy);
   if( desc == NULL )
      SCIPABORT();
   assert(iscopy);

   JNISCIP_CALL( SCIPincludeExternalCodeInformation(scip, name, desc) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** returns an array of names of currently included external codes */
JNIEXPORT
jobjectArray JNISCIP(getExternalCodeNames)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIPerrorMessage("method getExternalCodeNames is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** returns an array of the descriptions of currently included external codes
 *
 *  @note some descriptions may be NULL
 */
JNIEXPORT
jobjectArray JNISCIP(getExternalCodeDescriptions)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIPerrorMessage("method getExternalCodeDescriptions is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
}

/** returns the number of currently included information on external codes */
JNIEXPORT
jint JNISCIP(getNExternalCodes)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNExternalCodes(scip);

   return (jint) num;
}

/** prints information on external codes to a file stream via the message handler system
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
JNIEXPORT
void JNISCIP(printExternalCodes)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jfile               /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   SCIPprintExternalCodes(scip, (FILE*) (size_t)jfile);
}

/** returns if the dialog already exists
 *
 *  @return TRUE is returned if the dialog exits, otherwise FALSE.
 */
JNIEXPORT
jboolean JNISCIP(existsDialog)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jdialog             /**< dialog */
   )
{
   SCIP* scip;
   SCIP_Bool dialog;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   dialog = SCIPexistsDialog(scip, (SCIP_DIALOG*) (size_t)jdialog);

   return (jboolean) dialog;
}

/** captures a dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(captureDialog)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jdialog             /**< dialog */
   )
{
   SCIP* scip;
   SCIP_DIALOG* dialog;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   dialog = (SCIP_DIALOG*) (size_t) jdialog;

   JNISCIP_CALL( SCIPcaptureDialog(scip, dialog) );
}

/** releases a dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(releaseDialog)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jdialog             /**< dialog */
   )
{
   SCIP* scip;
   SCIP_DIALOG* dialog;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   dialog = (SCIP_DIALOG*) (size_t) jdialog;
   JNISCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
}

/** makes given dialog the root dialog of SCIP's interactive user shell; captures dialog and releases former root dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(setRootDialog)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jdialog             /**< dialog */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetRootDialog(scip, (SCIP_DIALOG*) (size_t)jdialog) );
}

/** returns the root dialog of SCIP's interactive user shell
 *
 *  @return the root dialog of SCIP's interactive user shell is returned.
 */
JNIEXPORT
jlong JNISCIP(getRootDialog)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_DIALOG* root;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   root = SCIPgetRootDialog(scip);

   return (jlong) (size_t) root;
}

/** adds a sub dialog to the given dialog as menu entry and captures it
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(addDialogEntry)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jdialog,            /**< dialog */
   jlong                 jsubdialog          /**< subdialog to add as menu entry in dialog */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPaddDialogEntry(scip, (SCIP_DIALOG*) (size_t)jdialog, (SCIP_DIALOG*) (size_t) jsubdialog) );
}

/** adds a single line of input which is treated as if the user entered the command line
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(addDialogInputLine)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jinputline          /**< input line to add */
   )
{
   SCIP* scip;
   const char* inputline;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   inputline = (*env)->GetStringUTFChars(env, jinputline, &iscopy);
   if( inputline == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPaddDialogInputLine(scip, inputline) );

   (*env)->ReleaseStringUTFChars(env, jinputline, inputline);
}


/** adds a single line of input to the command history which can be accessed with the cursor keys
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(addDialogHistoryLine)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jline               /**< input line to add */
   )
{
   SCIP* scip;
   const char* line;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into const char* */
   line = (*env)->GetStringUTFChars(env, jline, &iscopy);
   if( line == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPaddDialogHistoryLine(scip, line) );

   (*env)->ReleaseStringUTFChars(env, jline, line);
}

/** starts interactive mode of SCIP by executing the root dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post After calling this method \SCIP reaches one of the following stages depending on if and when the
 *        interactive shell was closed:
 *        - \ref SCIP_STAGE_PROBLEM if the interactive shell was closed after the problem was created
 *        - \ref SCIP_STAGE_TRANSFORMED if the interactive shell was closed after the problem was transformed
 *        - \ref SCIP_STAGE_PRESOLVING if the interactive shell was closed  during presolving
 *        - \ref SCIP_STAGE_PRESOLVED if the interactive shell was closed after presolve
 *        - \ref SCIP_STAGE_SOLVING if the interactive shell was closed during the tree search
 *        - \ref SCIP_STAGE_SOLVED if the interactive shell was closed after the problem was solved
 *        - \ref SCIP_STAGE_FREE if the interactive shell was closed after the problem was freed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
void JNISCIP(startInteraction)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPstartInteraction(scip) );
}

/** creates empty problem and initializes all solving data structures (the objective sense is set to MINIMIZE)
 *  If the problem type requires the use of variable pricers, these pricers should be added to the problem with calls
 *  to {@link activatePricer()}. These pricers are automatically deactivated, when the problem is freed.
 */
JNIEXPORT
void JNISCIP(createProbBasic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< problem name */
   )
{
   SCIP* scip;
   const char* name;
   jboolean iscopy;

   /* convert JNI string into const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPcreateProbBasic(scip, name) );
   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** reads problem from file and initializes all solving data structures */
JNIEXPORT
void JNISCIP(readProb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jfilename,          /**< problem file name */
   jstring               jextension          /**< extension of the desired file reader,
					      *   or an empty string if file extension should be used */
   )
{
   SCIP* scip;
   const char* filename;
   const char* extension;
   jboolean iscopy;

   /* convert JNI string into const char* */
   filename = (*env)->GetStringUTFChars(env, jfilename, &iscopy);
   if( filename == NULL )
      SCIPABORT();

   assert(iscopy);

   /* convert JNI string into const char* */
   if( (*env)->GetStringUTFLength(env, jextension) == 0 )
      extension = NULL;
   else
   {
      extension = (*env)->GetStringUTFChars(env, jextension, &iscopy);
      assert(iscopy);
   }

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPreadProb(scip, filename, extension) );
   (*env)->ReleaseStringUTFChars(env, jextension, extension);
   (*env)->ReleaseStringUTFChars(env, jfilename, filename);
}

/** writes original problem to file  */
JNIEXPORT
void JNISCIP(writeOrigProblem)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jfilename,          /**< output file (or NULL for standard output) */
   jstring               jextension,         /**< extension of the desired file reader,
					      *   or NULL if file extension should be used */
   jboolean              jgenericnames       /**< using generic variable and constraint names? */
   )
{
   SCIP* scip;
   const char* filename;
   const char* extension;
   SCIP_Bool genericnames;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   if( jfilename != NULL )
      filename = (*env)->GetStringUTFChars(env, jfilename, &iscopy);
   else
      filename = NULL;

   /* convert JNI string into C const char* */
   if( jextension != NULL )
      extension = (*env)->GetStringUTFChars(env, jextension, &iscopy);
   else
      extension = NULL;

   /* convert JNI boolean into C boolean */
   genericnames = (SCIP_Bool) jgenericnames;

   JNISCIP_CALL( SCIPwriteOrigProblem(scip, filename, extension, genericnames) );
}

/** writes original problem to file  */
JNIEXPORT
void JNISCIP(writeTransProblem)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jfilename,          /**< output file (or NULL for standard output) */
   jstring               jextension,         /**< extension of the desired file reader,
					      *   or NULL if file extension should be used */
   jboolean              jgenericnames       /**< using generic variable and constraint names? */
   )
{
   SCIP* scip;
   const char* filename;
   const char* extension;
   SCIP_Bool genericnames;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   if( jfilename != NULL )
      filename = (*env)->GetStringUTFChars(env, jfilename, &iscopy);
   else
      filename = NULL;

   /* convert JNI string into C const char* */
   if( jextension != NULL )
      extension = (*env)->GetStringUTFChars(env, jextension, &iscopy);
   else
      extension = NULL;

   /* convert JNI boolean into C boolean */
   genericnames = (SCIP_Bool) jgenericnames;

   JNISCIP_CALL( SCIPwriteTransProblem(scip, filename, extension, genericnames) );
}

/** frees problem and solution process data */
JNIEXPORT
void JNISCIP(freeProb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
 {
    SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPfreeProb(scip) );
}

/** permutes parts of the problem data structure
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 */
JNIEXPORT
void JNISCIP(permuteProb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  randseed,           /**< seed value for random generator */
   jboolean              permuteconss,       /**< should the list of constraints in each constraint handler be permuted? */
   jboolean              permutebinvars,     /**< should the list of binary variables be permuted? */
   jboolean              permuteintvars,     /**< should the list of integer variables be permuted? */
   jboolean              permuteimplvars,    /**< should the list of implicit integer variables be permuted? */
   jboolean              permutecontvars     /**< should the list of continuous integer variables be permuted? */
   )
 {
    SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPpermuteProb(scip, (unsigned int)randseed, (SCIP_Bool)permuteconss, (SCIP_Bool)permutebinvars, (SCIP_Bool)permuteintvars, (SCIP_Bool)permuteimplvars, (SCIP_Bool)permutecontvars ) );
}

/** gets user problem data
 *
 *  @return a SCIP_PROBDATA pointer, or NULL if no problem data was allocated
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jlong JNISCIP(getProbData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
 {
    SCIP* scip;
    SCIP_PROBDATA* prob;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   prob = SCIPgetProbData(scip);
   return (jlong) (size_t) prob;
}

/** sets user problem data
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(setProbData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 probdata            /**< user problem data to use */
   )
 {
    SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetProbData(scip, (SCIP_PROBDATA*) (size_t) probdata) );
}

/** returns name of the current problem instance
 *
 *  @return name of the current problem instance
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jstring JNISCIP(getProbName)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
 {
    SCIP* scip;
    const char* probname;
    jstring jprobname;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   probname = SCIPgetProbName(scip);

   jprobname = (*env)->NewStringUTF(env, probname);

   return jprobname;
}

/** sets name of the current problem instance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(setProbName)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< user problem data to use */
   )
 {
    SCIP* scip;
    const char* name;
    jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPsetProbName(scip, name) );
   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** gets objective sense of original problem */
JNIEXPORT
jint JNISCIP(getObjsense)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip              /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_OBJSENSE objsense;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   objsense = SCIPgetObjsense(scip);

   return (jint) objsense;
}

/** sets objective sense of problem */
JNIEXPORT
void JNISCIP(setObjsense)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  jobjsense           /**< new objective sense */
   )
{
   SCIP* scip;
   SCIP_OBJSENSE objsense;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);


   /* convert JNI int into C SCIP_OBJSENSE */
   objsense = (SCIP_OBJSENSE) jobjsense;

   /* check if objsense is covert in the JNI interface */
   assert(objsense == JNIPACKAGENAME(JniScipObjsense_SCIP_OBJSENSE_MAXIMIZE)
      || objsense == JNIPACKAGENAME(JniScipObjsense_SCIP_OBJSENSE_MINIMIZE));

   JNISCIP_CALL( SCIPsetObjsense(scip, objsense) );
}

/** adds offset of objective function */
JNIEXPORT
void JNISCIP(addObjoffset)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               jaddval             /**< value to add to objective offset */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPaddObjoffset(scip, (SCIP_Real) jaddval) );
}

/** returns the objective offset of the original problem */
JNIEXPORT
jdouble JNISCIP(getOrigObjoffset)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real offset;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   offset = SCIPgetOrigObjoffset(scip);

   return (jdouble) offset;
}

/** returns the objective scale of the original problem */
JNIEXPORT
jdouble JNISCIP(getOrigObjscale)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real scale;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   scale = SCIPgetOrigObjscale(scip);

   return (jdouble) scale;
}

/** returns the objective offset of the transformed problem */
JNIEXPORT
jdouble JNISCIP(getTransObjoffset)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real transoffset;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   transoffset = SCIPgetTransObjoffset(scip);

   return (jdouble) transoffset;
}

/** returns the objective scale of the transformed problem */
JNIEXPORT
jdouble JNISCIP(getTransObjscale)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real transscale;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   transscale = SCIPgetTransObjscale(scip);

   return (jdouble) transscale;
}

/** sets limit on objective function, such that only solutions better than this limit are accepted */
JNIEXPORT
void JNISCIP(setObjlimit)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               jobjlimit           /**< new primal objective limit */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetObjlimit(scip, (SCIP_Real) jobjlimit) );
}

/** gets current limit on objective function */
JNIEXPORT
jdouble JNISCIP(getObjlimit)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real objlimit;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   objlimit = SCIPgetObjlimit(scip);

   return (jdouble) objlimit;
}

/** informs SCIP, that the objective value is always integral in every feasible solution */
JNIEXPORT
void JNISCIP(setObjIntegral)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetObjIntegral(scip) );
}

/** returns whether the objective value is known to be integral in every feasible solution */
JNIEXPORT
jboolean JNISCIP(isObjIntegral)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Bool integral;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   integral = SCIPisObjIntegral(scip);

   return (jboolean) integral;
}

/** returns the Euclidean norm of the objective function vector (available only for transformed problem) */
JNIEXPORT
jdouble JNISCIP(getObjNorm)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real norm;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   norm = SCIPgetObjNorm(scip);

   return (jdouble) norm;
}

/** adds variable to the problem */
JNIEXPORT
void JNISCIP(addVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to add */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddVar(scip, var) );
}

/** adds variable to the problem and uses it as pricing candidate to enter the LP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(addPricedVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to add */
   jdouble               score               /**< pricing score of variable (the larger, the better the variable) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddPricedVar(scip, var, (SCIP_Real)score) );
}

/** removes variable from the problem; however, the variable is NOT removed from the constraints */
JNIEXPORT
jboolean JNISCIP(delVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar               /**< variable to delete */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_Bool deleted;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPdelVar(scip, var, &deleted) );

   return (jboolean) deleted;
}

/** gets array with active problem variables; data may become invalid after
 *  calls to {@link chgVarType()}, {@link fixVar()}, {@link aggregateVars()}, and {@link multiaggregateVar()}
 */
JNIEXPORT
jlongArray JNISCIP(getVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nvars;

   jlongArray jvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nvars = SCIPgetNVars(scip);
   jvars = (*env)->NewLongArray(env, nvars);

   if( jvars == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_VAR** vars;

      vars = SCIPgetVars(scip);
      (*env)->SetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);
   }

   return jvars;
}

/** gets number of active problem variables */
JNIEXPORT
jint JNISCIP(getNVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nvars = SCIPgetNVars(scip);

   return (jint)nvars;
}

/** gets number of binary active problem variables */
JNIEXPORT
jint JNISCIP(getNBinVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNBinVars(scip);

   return (jint)number;
}

/** gets number of integer active problem variables */
JNIEXPORT
jint JNISCIP(getNIntVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNIntVars(scip);

   return (jint)number;
}

/** gets number of implicit integer active problem variables */
JNIEXPORT
jint JNISCIP(getNImplVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNImplVars(scip);

   return (jint)number;
}

/** gets number of continuous active problem variables */
JNIEXPORT
jint JNISCIP(getNContVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNContVars(scip);

   return (jint)number;
}

/** gets array with fixed and aggregated problem variables; data may become invalid after
 *  calls to {@link fixVar()}, {@link aggregateVars()}, and {@link multiaggregateVar()}
 */
JNIEXPORT
jlongArray JNISCIP(getFixedVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nfixedvars;

   jlongArray jfixedvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nfixedvars = SCIPgetNFixedVars(scip);
   jfixedvars = (*env)->NewLongArray(env, nfixedvars);

   if( jfixedvars == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_VAR** fixedvars;

      fixedvars = SCIPgetFixedVars(scip);
      (*env)->SetLongArrayRegion(env, jfixedvars, 0, nfixedvars, (jlong*)fixedvars);
   }

   return jfixedvars;
}

/** gets number of fixed or aggregated problem variables */
JNIEXPORT
jint JNISCIP(getNFixedVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNFixedVars(scip);

   return (jint)number;
}

/** gets array with original problem variables; data may become invalid after
 *  a call to {@link chgVarType()}
 */
JNIEXPORT
jlongArray JNISCIP(getOrigVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nvars;

   jlongArray jvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nvars = SCIPgetNOrigVars(scip);
   jvars = (*env)->NewLongArray(env, nvars);

   if( jvars == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_VAR** vars;

      vars = SCIPgetOrigVars(scip);
      (*env)->SetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);
   }

   return jvars;
}

/** gets number of original problem variables */
JNIEXPORT
jint JNISCIP(getNOrigVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNOrigVars(scip);

   return (jint)number;
}

/** gets number of binary original problem variables */
JNIEXPORT
jint JNISCIP(getNOrigBinVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNOrigBinVars(scip);

   return (jint)number;
}

/** gets number of integer original problem variables */
JNIEXPORT
jint JNISCIP(getNOrigIntVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNOrigIntVars(scip);

   return (jint)number;
}

/** gets number of implicit integer original problem variables */
JNIEXPORT
jint JNISCIP(getNOrigImplVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNOrigImplVars(scip);

   return (jint)number;
}

/** gets number of continuous original problem variables */
JNIEXPORT
jint JNISCIP(getNOrigContVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNOrigContVars(scip);

   return (jint)number;
}

/** gets number of all problem variables created during creation and solving of problem;
 *  this includes also variables that were deleted in the meantime
 *
 *  @return the number of all problem variables created during creation and solving of problem;
 *          this includes also variables that were deleted in the meantime
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jint JNISCIP(getNTotalVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNTotalVars(scip);

   return (jint)number;
}

/** returns variable of given name in the problem, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of variable to find */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_VAR* var;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   var = SCIPfindVar(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t) var;
}

/** returns TRUE iff all potential variables exist in the problem, and FALSE, if there may be additional variables,
 *  that will be added in pricing and improve the objective value
 *
 *  @return TRUE, if all potential variables exist in the problem; FALSE, otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jboolean JNISCIP(allVarsInProb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   jboolean bool;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   bool = SCIPallVarsInProb(scip);

   return (jboolean) bool;
}

/** adds constraint to the problem; if constraint is only valid locally, it is added to the local subproblem of the
 *  current node (and all of its subnodes); otherwise it is added to the global problem;
 *  if a local constraint is added at the root node, it is automatically upgraded into a global constraint
 */
JNIEXPORT
void JNISCIP(addCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint to add */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPaddCons(scip, cons) );
}

/** globally removes constraint from all subproblems; removes constraint from the constraint set change data of the
 *  node, where it was added, or from the problem, if it was a problem constraint
 */
JNIEXPORT
void JNISCIP(delCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint to add */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPdelCons(scip, cons) );
}

/** returns original constraint of given name in the problem, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findOrigCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of original constraint */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_CONS* cons;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   cons = SCIPfindOrigCons(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t) cons;
}

/** returns constraint of given name in the problem, or NULL if not existing */
JNIEXPORT
jlong JNISCIP(findCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname               /**< name of constraint */
   )
{
   SCIP* scip;
   const char* name;
   SCIP_CONS* cons;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   cons = SCIPfindCons(scip, name);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t) cons;
}

/** gets number of upgraded constraints */
JNIEXPORT
jint JNISCIP(getNUpgrConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int num;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   num = SCIPgetNUpgrConss(scip);

   return (jint) num;
}

/** gets total number of globally valid constraints currently in the problem */
JNIEXPORT
jint JNISCIP(getNConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNConss(scip);

   return (jint)number;
}

/** gets array of globally valid constraints currently in the problem */
JNIEXPORT
jlongArray JNISCIP(getConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nconss;

   jlongArray jconss;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nconss = SCIPgetNConss(scip);
   jconss = (*env)->NewLongArray(env, nconss);

   if( jconss == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_CONS** conss;

      conss = SCIPgetConss(scip);
      (*env)->SetLongArrayRegion(env, jconss, 0, nconss, (jlong*)conss);
   }

   return jconss;
}

/** gets total number of constraints in the original problem */
JNIEXPORT
jint JNISCIP(getNOrigConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int number;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   number = SCIPgetNOrigConss(scip);

   return (jint)number;
}

/** gets array of constraints in the original problem */
JNIEXPORT
jlongArray JNISCIP(getOrigConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nconss;

   jlongArray jconss;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nconss = SCIPgetNConss(scip);
   jconss = (*env)->NewLongArray(env, nconss);

   if( jconss == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_CONS** conss;

      conss = SCIPgetOrigConss(scip);
      (*env)->SetLongArrayRegion(env, jconss, 0, nconss, (jlong*)conss);
   }

   return jconss;
}

/** adds constraint to the given node (and all of its subnodes), even if it is a global constraint;
 *  It is sometimes desirable to add the constraint to a more local node (i.e., a node of larger depth) even if
 *  the constraint is also valid higher in the tree, for example, if one wants to produce a constraint which is
 *  only active in a small part of the tree although it is valid in a larger part.
 *  In this case, one should pass the more global node where the constraint is valid as "validnode".
 *  Note that the same constraint cannot be added twice to the branching tree with different "validnode" parameters.
 *  If the constraint is valid at the same node as it is inserted (the usual case), one should pass NULL as "validnode".
 *  If the "validnode" is the root node, it is automatically upgraded into a global constraint, but still only added to
 *  the given node. If a local constraint is added to the root node, it is added to the global problem instead.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(addConsNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnode,              /**< node to add constraint to */
   jlong                 jcons,              /**< constraint to add */
   jlong                 jvalidnode          /**< node at which the constraint is valid, or NULL */
   )
{
   SCIP* scip;
   SCIP_NODE* node;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   node = (SCIP_NODE*) (size_t) jnode;
   assert(node != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPaddConsNode(scip, node, cons, (SCIP_NODE*) (size_t)jvalidnode) );
}

/** adds constraint locally to the current node (and all of its subnodes), even if it is a global constraint;
 *  It is sometimes desirable to add the constraint to a more local node (i.e., a node of larger depth) even if
 *  the constraint is also valid higher in the tree, for example, if one wants to produce a constraint which is
 *  only active in a small part of the tree although it is valid in a larger part.
 *
 *  If the constraint is valid at the same node as it is inserted (the usual case), one should pass NULL as "validnode".
 *  If the "validnode" is the root node, it is automatically upgraded into a global constraint, but still only added to
 *  the given node. If a local constraint is added to the root node, it is added to the global problem instead.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note The same constraint cannot be added twice to the branching tree with different "validnode" parameters. This is
 *        the case due internal data structures and performance issues. In such a case you should try to realize your
 *        issue using the method SCIPdisableCons() and SCIPenableCons() and control these via the event system of SCIP.
 */
JNIEXPORT
void JNISCIP(addConsLocal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint to add */
   jlong                 jvalidnode          /**< node at which the constraint is valid, or NULL */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPaddConsLocal(scip, cons, (SCIP_NODE*) (size_t)jvalidnode) );
}

/** disables constraint's separation, enforcing, and propagation capabilities at the given node (and all subnodes);
 *  if the method is called at the root node, the constraint is globally deleted from the problem;
 *  the constraint deletion is being remembered at the given node, s.t. after leaving the node's subtree, the constraint
 *  is automatically enabled again, and after entering the node's subtree, it is automatically disabled;
 *  this may improve performance because redundant checks on this constraint are avoided, but it consumes memory;
 *  alternatively, use SCIPdisableCons()
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(delConsNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnode,              /**< node to add constraint to */
   jlong                 jcons               /**< constraint to add */
   )
{
   SCIP* scip;
   SCIP_NODE* node;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   node = (SCIP_NODE*) (size_t) jnode;
   assert(node != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPdelConsNode(scip, node, cons) );
}

/** disables constraint's separation, enforcing, and propagation capabilities at the current node (and all subnodes);
 *  if the method is called during problem modification or at the root node, the constraint is globally deleted from
 *  the problem;
 *  the constraint deletion is being remembered at the current node, s.t. after leaving the current subtree, the
 *  constraint is automatically enabled again, and after reentering the current node's subtree, it is automatically
 *  disabled again;
 *  this may improve performance because redundant checks on this constraint are avoided, but it consumes memory;
 *  alternatively, use SCIPdisableCons()
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 *
 */
JNIEXPORT
void JNISCIP(delConsLocal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint to add */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPdelConsLocal(scip, cons) );
}

/** gets estimate of best primal solution w.r.t. original problem contained in current subtree
 *
 *  @return estimate of best primal solution w.r.t. original problem contained in current subtree
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getLocalOrigEstimate)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetLocalOrigEstimate(scip) ;
}

/** gets estimate of best primal solution w.r.t. transformed problem contained in current subtree
 *
 *  @return estimate of best primal solution w.r.t. transformed problem contained in current subtree
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getLocalTransEstimate)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetLocalTransEstimate(scip) ;
}

/** gets dual bound of current node
 *
 *  @return dual bound of current node
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getLocalDualbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetLocalDualbound(scip) ;
}

/** gets lower bound of current node in transformed problem
 *
 *  @return lower bound  of current node in transformed problem
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getLocalLowerbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetLocalLowerbound(scip) ;
}

/** gets dual bound of given node
 *
 *  @return dual bound of a given node
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getNodeDualbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnode               /**< node to get dual bound for */
   )
{
   SCIP* scip;
   SCIP_NODE* node;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   node = (SCIP_NODE*) (size_t) jnode;
   assert(node != NULL);

   return (jdouble) SCIPgetNodeDualbound(scip, node) ;
}

/** gets lower bound of given node in transformed problem
 *
 *  @return lower bound  of given node in transformed problem
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getNodeLowerbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnode               /**< node to get dual bound for */
   )
{
   SCIP* scip;
   SCIP_NODE* node;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   node = (SCIP_NODE*) (size_t) jnode;
   assert(node != NULL);

   return (jdouble) SCIPgetNodeLowerbound(scip, node) ;
}

/** if given value is tighter (larger for minimization, smaller for maximization) than the current node's dual bound,
 *  sets the current node's dual bound to the new value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(updateLocalDualbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               jnewbound           /**< new dual bound for the node (if it's tighter than the old one) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPupdateLocalDualbound(scip, (SCIP_Real)jnewbound) );
}

/** if given value is larger than the current node's lower bound (in transformed problem), sets the current node's
 *  lower bound to the new value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(updateLocalLowerbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               jnewbound           /**< new lower bound for the node (if it's larger than the old one) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPupdateLocalLowerbound(scip, (SCIP_Real)jnewbound) );
}

/** if given value is tighter (larger for minimization, smaller for maximization) than the node's dual bound,
 *  sets the node's dual bound to the new value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(updateNodeDualbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnode,              /**< node to update dual bound for */
   jdouble               jnewbound           /**< new dual bound for the node (if it's tighter than the old one) */
   )
{
   SCIP* scip;
   SCIP_NODE* node;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   node = (SCIP_NODE*) (size_t) jnode;
   assert(node != NULL);

   JNISCIP_CALL( SCIPupdateNodeDualbound(scip, node, (SCIP_Real)jnewbound) );
}

/** if given value is larger than the node's lower bound (in transformed problem), sets the node's lower bound
 *  to the new value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(updateNodeLowerbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnode,              /**< node to update dual bound for */
   jdouble               jnewbound           /**< new lower bound for the node (if it's larger than the old one) */
   )
{
   SCIP* scip;
   SCIP_NODE* node;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);


   node = (SCIP_NODE*) (size_t) jnode;
   assert(node != NULL);

   JNISCIP_CALL( SCIPupdateNodeLowerbound(scip, node, (SCIP_Real)jnewbound) );
}

/** change the node selection priority of the given child
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(chgChildPrio)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnode,              /**< child to update the node selection priority */
   jdouble               jpriority           /**< node selection priority value */
   )
{
   SCIP* scip;
   SCIP_NODE* node;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   node = (SCIP_NODE*) (size_t) jnode;
   assert(node != NULL);

   JNISCIP_CALL( SCIPchgChildPrio(scip, node, (SCIP_Real)jpriority) );
}

/** initializes solving data structures and transforms problem */
JNIEXPORT
void JNISCIP(transformProb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPtransformProb(scip) );
}

/** transforms and presolves problem */
JNIEXPORT
void JNISCIP(presolve)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPpresolve(scip) );
}

/** transforms, presolves, and solves problem */
JNIEXPORT
void JNISCIP(solve)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsolve(scip) );
}

/** frees branch and bound tree and all solution process data; statistics, presolving data and transformed problem is
 *  preserved
 */
JNIEXPORT
void JNISCIP(freeSolve)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jboolean              jrestart            /**< should certain data be preserved for improved restarting? */
   )
{
   SCIP* scip;
   SCIP_Bool restart;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   restart = (SCIP_Bool) jrestart;

   JNISCIP_CALL( SCIPfreeSolve(scip, restart) );
}

/** frees all solution process data including presolving and transformed problem, only original problem is kept */
JNIEXPORT
void JNISCIP(freeTransform)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPfreeTransform(scip) );
}

/** interrupts solving process as soon as possible (e.g., after the current node has been solved) */
JNIEXPORT
void JNISCIP(interruptSolve)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPinterruptSolve(scip) );
}

/** restarts solving process as soon as possible (e.g., after the current node has been solved) */
JNIEXPORT
void JNISCIP(restartSolve)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPrestartSolve(scip) );
}

/** whether we are in the restarting phase */
JNIEXPORT
jboolean JNISCIP(isInRestart)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Bool inrestart;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   inrestart = SCIPisInRestart(scip);

   return (jboolean) inrestart;
}

/** creates and captures problem variable; if variable is of integral type, fractional bounds are automatically rounded;
 *  an integer variable with bounds zero and one is automatically converted into a binary variable;
 *
 *  @warning When doing column generation and the original problem is a maximization problem, notice that SCIP will
 *           transform the problem into a minimization problem by multiplying the objective function by -1.  Thus, the
 *           original objective function value of variables created during the solving process has to be multiplied by
 *           -1, too.
 *
 *  @note the variable gets captured, hence at one point you have to release it using the method {@link releaseVar()}
 */
JNIEXPORT
jlong JNISCIP(createVarBasic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobje,              /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of variable */
   jdouble               jlb,                /**< lower bound of variable */
   jdouble               jub,                /**< upper bound of variable */
   jdouble               jobj,               /**< objective function value */
   jint                  jvartype            /**< type of variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   const char* name;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   SCIP_VARTYPE vartype;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   /* convert JNI double into SCIP_Real */
   lb = (SCIP_Real) jlb;

   /* convert JNI double into SCIP_Real */
   ub = (SCIP_Real) jub;

   assert(SCIPisLE(scip, lb, ub));

   /* convert JNI double into SCIP_Real */
   obj = (SCIP_Real) jobj;


   vartype = (SCIP_VARTYPE) jvartype;

   assert(vartype == JNIPACKAGENAME(JniScipVartype_SCIP_VARTYPE_BINARY)
      || vartype == JNIPACKAGENAME(JniScipVartype_SCIP_VARTYPE_INTEGER)
      || vartype == JNIPACKAGENAME(JniScipVartype_SCIP_VARTYPE_IMPLINT)
      || vartype == JNIPACKAGENAME(JniScipVartype_SCIP_VARTYPE_CONTINUOUS));

   JNISCIP_CALL( SCIPcreateVarBasic(scip, &var, name, lb, ub, obj, vartype) );
   SCIPdebugMessage("created variable <%s> [%g,%g] obj: <%g>, <%p>\n", name, lb, ub, obj, (void*)var);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)var;
}

/** outputs the variable name to the file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(writeVarName)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jfile,              /**< output file, or NULL for stdout */
   jlong                 jvar,               /**< variable array to output */
   jboolean              jtype               /**< should the variable type be also posted */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPwriteVarName(scip, (FILE*)(size_t)jfile, var, (SCIP_Bool) jtype) );
}

/** print the given list of variables to output stream separated by the given delimiter character;
 *
 *  i. e. the variables x1, x2, ..., xn with given delimiter ',' are written as: \<x1\>, \<x2\>, ..., \<xn\>;
 *
 *  the method SCIPparseVarsList() can parse such a string
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note The printing process is done via the message handler system.
 */
JNIEXPORT
void JNISCIP(writeVarsList)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jfile,              /**< output file, or NULL for stdout */
   jlongArray            jvars,              /**< variable array to output */
   jint                  nvars,              /**< number of variables */
   jboolean              type,               /**< should the variable type be also posted */
   jchar                 delimiter           /**< character which is used for delimitation */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)vars);

   JNISCIP_CALL( SCIPwriteVarsList(scip, (FILE*)(size_t)jfile, vars, (int)nvars, (SCIP_Bool)type, (char)delimiter) );

   SCIPfreeBufferArray(scip, &vars);
}

/** print the given variables and coefficients as linear sum in the following form
 *  c1 \<x1\> + c2 \<x2\>   ... + cn \<xn\>
 *
 *  This string can be parsed by the method SCIPparseVarsLinearsum().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note The printing process is done via the message handler system.
 */
JNIEXPORT
void JNISCIP(writeVarsLinearsum)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jfile,              /**< output file, or NULL for stdout */
   jlongArray            jvars,              /**< variable array to output */
   jdoubleArray          jvals,              /**< array of coefficients or NULL if all coefficients are 1.0 */
   jint                  nvars,              /**< number of variables */
   jboolean              type                /**< should the variable type be also posted */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)vars);
   (*env)->GetDoubleArrayRegion(env, jvals, 0, (int)nvars, (jdouble*)vals);

   JNISCIP_CALL( SCIPwriteVarsLinearsum(scip, (FILE*)(size_t)jfile, vars, vals, (int)nvars, (SCIP_Bool)type) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
}

/** TODO: writeVarsPolynomial - array of arrays?  */
/** print the given monomials as polynomial in the following form
 *  c1 \<x11\>^e11 \<x12\>^e12 ... \<x1n\>^e1n + c2 \<x21\>^e21 \<x22\>^e22 ... + ... + cn \<xn1\>^en1 ...
 *
 *  This string can be parsed by the method SCIPparseVarsPolynomial().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note The printing process is done via the message handler system.
 */
JNIEXPORT
void JNISCIP(writeVarsPolynomial)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jfile,              /**< output file, or NULL for stdout */
   jobjectArray          jmonomialvars,      /**< arrays with variables for each monomial */
   jobjectArray          jmonomialexps,      /**< arrays with variable exponents, or NULL if always 1.0 */
   jdoubleArray          jmonomialcoefs,     /**< array with monomial coefficients */
   jintArray             jmonomialnvars,     /**< array with number of variables for each monomial */
   jint                  nmonomials,         /**< number of monomials */
   jboolean              type                /**< should the variable type be also posted */
   )
{
   SCIPerrorMessage("method writeVarsPolynomial is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );
}

/** increases usage counter of variable */
JNIEXPORT
void JNISCIP(captureVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable  */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPcaptureVar(scip, var) );
}

/** decreases usage counter of variable, and frees memory if necessary */
JNIEXPORT
void JNISCIP(releaseVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable  */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPreleaseVar(scip, &var) );
}

/** changes the name of a variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_PROBLEM
 */
JNIEXPORT
void JNISCIP(chgVarName)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable  */
   jstring               jname               /**< new name of constraint */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPchgVarName(scip, var, name) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** gets and captures transformed variable of a given variable; if the variable is not yet transformed,
 *  a new transformed variable for this variable is created
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(transformVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get/create transformed variable for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_VAR* transvar;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPtransformVar(scip, var, &transvar) );

   return (jlong) (size_t) transvar;
}

/** gets and captures transformed variables for an array of variables;
 *  if a variable of the array is not yet transformed, a new transformed variable for this variable is created;
 *  it is possible to call this method with vars == transvars
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlongArray JNISCIP(transformVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get/create transformed variables for */
   jlongArray            jvars               /**< array with variables to get/create transformed variables for */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;

   jlongArray jtransvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)vars);

   jtransvars = (*env)->NewLongArray(env, (int)nvars);

   if (jtransvars == NULL)
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_VAR** transvars;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &transvars, (int)nvars) );

      JNISCIP_CALL( SCIPtransformVars(scip, (int)nvars, vars, transvars) );
      (*env)->SetLongArrayRegion(env, jtransvars, 0, (int)nvars, (jlong*)transvars);

      SCIPfreeBufferArray(scip, &transvars);
   }

   SCIPfreeBufferArray(scip, &vars);

   return jtransvars;
}

/** gets corresponding transformed variable of a given variable;
 *  returns NULL as transvar, if transformed variable is not yet existing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jlong JNISCIP(getTransformedVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get transformed variable for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_VAR* transvar;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPgetTransformedVar(scip, var, &transvar) );

   return (jlong) (size_t) transvar;
}

/** gets corresponding transformed variables for an array of variables;
 *  stores NULL in a transvars slot, if the transformed variable is not yet existing;
 *  it is possible to call this method with vars == transvars, but remember that variables that are not
 *  yet transformed will be replaced with NULL
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jlongArray JNISCIP(getTransformedVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get transformed variables for */
   jlongArray            jvars               /**< array with variables to get transformed variables for */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;

   jlongArray jtransvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)vars);
   jtransvars = (*env)->NewLongArray(env, (int)nvars);

   if (jtransvars == NULL)
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );

      jtransvars = NULL;
   }
   else
   {
      SCIP_VAR** transvars;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &transvars, (int)nvars) );

      JNISCIP_CALL( SCIPgetTransformedVars(scip, (int)nvars, vars, transvars) );
      (*env)->SetLongArrayRegion(env, jtransvars, 0, (int)nvars, (jlong*)transvars);

      SCIPfreeBufferArray(scip, &transvars);
   }

   SCIPfreeBufferArray(scip, &vars);

   return jtransvars;
}

/** gets negated variable x' = lb + ub - x of variable x; negated variable is created, if not yet existing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jlong JNISCIP(getNegatedVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get negated variable for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_VAR* negvar;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPgetNegatedVar(scip, var, &negvar) );

   return (jlong) (size_t) negvar;
}

/** gets negated variables x' = lb + ub - x of variables x; negated variables are created, if not yet existing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jlongArray JNISCIP(getNegatedVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  nvars,              /**< number of variables to get negated variables for */
   jlongArray            jvars               /**< array of variables to get negated variables for */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;


   jlongArray jnegvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)vars);
   jnegvars = (*env)->NewLongArray(env, (int)nvars);

   if (jnegvars == NULL)
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );

      jnegvars = NULL;
   }
   else
   {
      SCIP_VAR** negvars;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &negvars, (int)nvars) );

      JNISCIP_CALL( SCIPgetNegatedVars(scip, (int)nvars, vars, negvars) );
      (*env)->SetLongArrayRegion(env, jnegvars, 0, (int)nvars, (jlong*)negvars);

      SCIPfreeBufferArray(scip, &negvars);
   }

   SCIPfreeBufferArray(scip, &vars);

   return jnegvars;
}

/** flattens aggregation graph of multi-aggregated variable in order to avoid exponential recursion later on
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
void JNISCIP(flattenVarAggregationGraph)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< problem variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPflattenVarAggregationGraph(scip, var) );
}

/** return for given variables all their active counterparts; all active variables will be pairwise different
 *  @note It does not hold that the first output variable is the active variable for the first input variable.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jint JNISCIP(getActiveVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlongArray            jvars,              /**< variable array with given variables and as output all active
					      *   variables, if enough slots exist
					      */
   jintArray             jnvars,             /**< number of given variables, and as output number of active variables,
					      *   if enough slots exist
					      */
   jint                  varssize            /**< available slots in vars array */
   )
{
   SCIPerrorMessage("method getActiveVars is not implemented yet (deprecated)\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
#if 0
   SCIP* scip;
   SCIP_VAR** vars;
   int* nvars;
   int requiredsize;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)varssize) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)varssize, (jlong*)vars);
   (*env)->GetIntArrayRegion(env, jnvars, 0, (int)varssize, (jint*)nvars);

   JNISCIP_CALL( SCIPgetActiveVars(scip, vars, nvars, (int)varssize, &requiredsize );

   SCIPfreeBufferArray(scip, &vars);

   return (jint)num;
#endif
}

/** returns the reduced costs of the variable in the current node's LP relaxation;
 *  the current node has to have a feasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jdouble JNISCIP(getVarRedcost)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get reduced costs, should be a column in current node LP */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarRedcost(scip, var);
}

/** returns the implied reduced costs of the variable in the current node's LP relaxation;
 *  the current node has to have a feasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jdouble JNISCIP(getVarImplRedcost)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to get reduced costs, should be a column in current node LP */
   jboolean              jvarfixing          /**< FALSE if for x == 0, TRUE for x == 1 */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarImplRedcost(scip, var, (SCIP_Bool)jvarfixing);
}

/** returns the Farkas coefficient of the variable in the current node's LP relaxation;
 *  the current node has to have an infeasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jdouble JNISCIP(getVarFarkasCoef)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get reduced costs, should be a column in current node LP */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarFarkasCoef(scip, var);
}

/** gets solution value for variable in current node
 *
 *  @return solution value for variable in current node
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jdouble JNISCIP(getVarSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get solution value for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarSol(scip, var);
}

/** gets solution values of multiple variables in current node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jdoubleArray JNISCIP(getVarSols)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  nvars,              /**< number of variables to get solution value for */
   jlongArray            jvars               /**< array with variables to get value for */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;

   jdoubleArray jvals;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)vars);
   jvals = (*env)->NewDoubleArray(env, (int)nvars);

   if (jvals == NULL)
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_Real* vals;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)nvars) );

      JNISCIP_CALL( SCIPgetVarSols(scip, (int)nvars, vars, vals) );
      (*env)->SetDoubleArrayRegion(env, jvals, 0, (int)nvars, (jdouble*)vals);

      SCIPfreeBufferArray(scip, &vals);
   }

   SCIPfreeBufferArray(scip, &vars);

   return jvals;
}

/** sets the solution value of all variables in the global relaxation solution to zero
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(clearRelaxSolVals)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPclearRelaxSolVals(scip) );
}

/** sets the value of the given variable in the global relaxation solution;
 *  this solution can be filled by the relaxation handlers  and can be used by heuristics and for separation;
 *  You can use SCIPclearRelaxSolVals() to set all values to zero, initially;
 *  after setting all solution values, you have to call SCIPmarkRelaxSolValid()
 *  to inform SCIP that the stored solution is valid
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(setRelaxSolVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to set value for */
   jdouble               jval                /**< solution value of variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPsetRelaxSolVal(scip, var, (SCIP_Real)jval) );
}

/** sets the values of the given variables in the global relaxation solution;
 *  this solution can be filled by the relaxation handlers  and can be used by heuristics and for separation;
 *  the solution is automatically cleared, s.t. all other variables get value 0.0
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(setRelaxSolVals)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  nvars,              /**< number of variables to set relaxation solution value for */
   jlongArray            jvars,              /**< array with variables to set value for */
   jdoubleArray          jvals               /**< array with solution values of variables */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);
   (*env)->GetDoubleArrayRegion(env, jvals, 0, nvars, (jdouble*)vals);

   JNISCIP_CALL( SCIPsetRelaxSolVals(scip, (int)nvars, vars, vals) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
}

/** sets the values of the variables in the global relaxation solution to the values
 *  in the given primal solution; the relaxation solution can be filled by the relaxation hanlders
 *  and might be used by heuristics and for separation
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(setRelaxSolValsSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< primal relaxation solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPsetRelaxSolValsSol(scip, sol) );
}

/** returns whether the relaxation solution is valid
 *
 *  @return TRUE, if the relaxation solution is valid; FALSE, otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jboolean JNISCIP(isRelaxSolValid)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisRelaxSolValid(scip);
}

/** informs SCIP, that the relaxation solution is valid
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(markRelaxSolValid)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPmarkRelaxSolValid(scip) );
}

/** informs SCIP, that the relaxation solution is invalid
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(markRelaxSolInvalid)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPmarkRelaxSolInvalid(scip) );
}

/** gets the relaxation solution value of the given variable
 *
 *  @return the relaxation solution value of the given variable
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jdouble JNISCIP(getRelaxSolVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get value for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetRelaxSolVal(scip, var);
}

/** gets the relaxation solution objective value
 *
 *  @return the objective value of the relaxation solution
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jdouble JNISCIP(getRelaxSolObj)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetRelaxSolObj(scip);
}

/** start strong branching - call before any strong branching
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(startStrongbranch)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jboolean              jenablepropagation  /**< should propagation be done before solving the strong branching LP? */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPstartStrongbranch(scip, (SCIP_Bool)jenablepropagation) );
}

/** end strong branching - call after any strong branching
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(endStrongbranch)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPendStrongbranch(scip) );
}

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given variable, or -1 if strong branching was never applied to the variable in current run
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jlong JNISCIP(getVarStrongbranchNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get last strong branching node for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jlong) SCIPgetVarStrongbranchNode(scip, var);
}

/** if strong branching was already applied on the variable at the current node, returns the number of LPs solved after
 *  the LP where the strong branching on this variable was applied;
 *  if strong branching was not yet applied on the variable at the current node, returns INT_MAX
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jlong JNISCIP(getVarStrongbranchLPAge)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get strong branching LP age for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jlong) SCIPgetVarStrongbranchLPAge(scip, var);
}

/** gets number of times, strong branching was applied in current run on the given variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jint JNISCIP(getVarNStrongbranchs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get last strong branching node for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jint) SCIPgetVarNStrongbranchs(scip, var);
}

/** adds given values to lock numbers of variable for rounding
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(addVarLocks)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jnlocksdown,        /**< modification in number of rounding down locks */
   jint                  jnlocksup           /**< modification in number of rounding up locks */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddVarLocks(scip, var, (int)jnlocksdown, (int)jnlocksup) );
}

/** locks rounding of variable with respect to the lock status of the constraint and its negation;
 *  this method should be called whenever the lock status of a variable in a constraint changes, for example if
 *  the coefficient of the variable changed its sign or if the left or right hand sides of the constraint were
 *  added or removed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(lockVarCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jlong                 jcons,              /**< constraint */
   jboolean              jlocksdown,         /**< should the rounding be locked in downwards direction? */
   jboolean              jlocksup            /**< should the rounding be locked in upwards direction? */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPlockVarCons(scip, var, cons, (SCIP_Bool)jlocksdown, (SCIP_Bool)jlocksup) );
}

/** unlocks rounding of variable with respect to the lock status of the constraint and its negation;
 *  this method should be called whenever the lock status of a variable in a constraint changes, for example if
 *  the coefficient of the variable changed its sign or if the left or right hand sides of the constraint were
 *  added or removed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(unlockVarCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jlong                 jcons,              /**< constraint */
   jboolean              jlocksdown,         /**< should the rounding be locked in downwards direction? */
   jboolean              jlocksup            /**< should the rounding be locked in upwards direction? */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPunlockVarCons(scip, var, cons, (SCIP_Bool)jlocksdown, (SCIP_Bool)jlocksup) );
}

/** changes variable's objective value */
JNIEXPORT
void JNISCIP(chgVarObj)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable  */
   jdouble               jnewobj             /**< new value for bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarObj(scip, var, (SCIP_Real) jnewobj) );
}

/** adds value to variable's objective value */
JNIEXPORT
void JNISCIP(addVarObj)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable  */
   jdouble               jaddobj             /**< new value for bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddVarObj(scip, var, (SCIP_Real) jaddobj) );
}

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) lower bound value;
 *  does not change the bounds of the variable
 *
 *  @return adjusted lower bound for the given variable; the bound of the variable is not changed
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jdouble JNISCIP(adjustedVarLb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to adjust the bound for */
   jdouble               jlb                 /**< lower bound value to adjust */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPadjustedVarLb(scip, var, (SCIP_Real) jlb);
}

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) lower bound value;
 *  does not change the bounds of the variable
 *
 *  @return adjusted lower bound for the given variable; the bound of the variable is not changed
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jdouble JNISCIP(adjustedVarUb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to adjust the bound for */
   jdouble               jub                 /**< upper bound value to adjust */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPadjustedVarUb(scip, var, (SCIP_Real) jub);
}

/** depending on SCIP's stage, changes lower bound of variable in the problem, in preprocessing, or in current node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 */
JNIEXPORT
void JNISCIP(chgVarLb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable  */
   jdouble               jnewbound           /**< new value for bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarLb(scip, var, (SCIP_Real) jnewbound) );
}

/** depending on SCIP's stage, changes upper bound of variable in the problem, in preprocessing, or in current node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 */
JNIEXPORT
void JNISCIP(chgVarUb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable  */
   jdouble               jnewbound           /**< new value for bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarUb(scip, var, (SCIP_Real) jnewbound) );
}

/** changes lower bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(chgVarLbNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnode,              /**< node to change bound at, or NULL for current node */
   jlong                 jvar,               /**< variable to change the bound for */
   jdouble               jnewbound           /**< new value for bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_NODE* node;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   node = (SCIP_NODE*) (size_t) jnode;

   JNISCIP_CALL( SCIPchgVarLbNode(scip, node, var, (SCIP_Real) jnewbound) );
}

/** changes upper bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(chgVarUbNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnode,              /**< node to change bound at, or NULL for current node */
   jlong                 jvar,               /**< variable to change the bound for */
   jdouble               jnewbound           /**< new value for bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_NODE* node;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   node = (SCIP_NODE*) (size_t) jnode;

   JNISCIP_CALL( SCIPchgVarUbNode(scip, node, var, (SCIP_Real) jnewbound) );
}

/** changes global lower bound of variable; if possible, adjust bound to integral value; also tightens the local bound,
 *  if the global bound is better than the local bound
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
JNIEXPORT
void JNISCIP(chgVarLbGlobal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable  */
   jdouble               jnewbound           /**< new value for bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarLbGlobal(scip, var, (SCIP_Real) jnewbound) );
}

/** changes global upper bound of variable; if possible, adjust bound to integral value; also tightens the local bound,
 *  if the global bound is better than the local bound
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
JNIEXPORT
void JNISCIP(chgVarUbGlobal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable  */
   jdouble               jnewbound           /**< new value for bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarUbGlobal(scip, var, (SCIP_Real) jnewbound) );
}

/** changes lazy lower bound of the variable, this is only possible if the variable is not in the LP yet
 *
 *  lazy bounds are bounds, that are enforced by constraints and the objective function; hence, these bounds do not need
 *  to be put into the LP explicitly.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note lazy bounds are useful for branch-and-price since the corresponding variable bounds are not part of the LP
 */
JNIEXPORT
void JNISCIP(chgVarLbLazy)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable  */
   jdouble               jlazylb           /**< the lazy lower bound to be set */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarLbLazy(scip, var, (SCIP_Real) jlazylb) );
}

/** changes lazy upper bound of the variable, this is only possible if the variable is not in the LP yet
 *
 *  lazy bounds are bounds, that are enforced by constraints and the objective function; hence, these bounds do not need
 *  to be put into the LP explicitly.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note lazy bounds are useful for branch-and-price since the corresponding variable bounds are not part of the LP
 */
JNIEXPORT
void JNISCIP(chgVarUbLazy)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable  */
   jdouble               jlazyub             /**< the lazy upper bound to be set */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarUbLazy(scip, var, (SCIP_Real) jlazyub) );
}

/** for a multi-aggregated variable, returns the global lower bound computed by adding the global bounds from all aggregation variables
 * this global bound may be tighter than the one given by SCIPvarGetLbGlobal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetLbGlobal
 *
 * @return the global lower bound computed by adding the global bounds from all aggregation variables
 */
JNIEXPORT
jdouble JNISCIP(computeVarLbGlobal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to compute the bound for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPcomputeVarLbGlobal(scip, var);
}

/** for a multi-aggregated variable, returns the global upper bound computed by adding the global bounds from all aggregation variables
 * this global bound may be tighter than the one given by SCIPvarGetUbGlobal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetUbGlobal
 *
 * @return the global upper bound computed by adding the global bounds from all aggregation variables
 */
JNIEXPORT
jdouble JNISCIP(computeVarUbGlobal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to compute the bound for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPcomputeVarUbGlobal(scip, var);
}

/** for a multi-aggregated variable, returns the local lower bound computed by adding the local bounds from all aggregation variables
 * this local bound may be tighter than the one given by SCIPvarGetLbLocal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetLbLocal
 *
 * @return the local lower bound computed by adding the global bounds from all aggregation variables
 */
JNIEXPORT
jdouble JNISCIP(computeVarLbLocal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to compute the bound for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPcomputeVarLbLocal(scip, var);
}

/** for a multi-aggregated variable, returns the local upper bound computed by adding the local bounds from all aggregation variables
 * this local bound may be tighter than the one given by SCIPvarGetUbLocal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetUbLocal
 *
 * @return the local upper bound computed by adding the global bounds from all aggregation variables
 */
JNIEXPORT
jdouble JNISCIP(computeVarUbLocal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to compute the bound for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPcomputeVarUbLocal(scip, var);
}

/** calculates a partition of the given set of binary variables into cliques;
 *  afterwards the output array contains one value for each variable, such that two variables got the same value iff they
 *  were assigned to the same clique;
 *  the first variable is always assigned to clique 0, and a variable can only be assigned to clique i if at least one of
 *  the preceding variables was assigned to clique i-1;
 *  for each clique at most 1 variables can be set to TRUE in a feasible solution;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(calcCliquePartition)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlongArray            jvars,              /**< binary variables in the clique from which at most one can be set to 1 */
   jint                  nvars,              /**< number of variables in the clique */
   jintArray             jcliquepartition    /**< array of length nvars to store the clique partition */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;
   int* cliquepartition;
   int ncliques;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);


   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &cliquepartition, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);
   (*env)->GetIntArrayRegion(env, jcliquepartition, 0, nvars, (jint*)cliquepartition);

   JNISCIP_CALL( SCIPcalcCliquePartition(scip, vars, (int)nvars, cliquepartition, &ncliques) );

   SCIPfreeBufferArray(scip, &cliquepartition);
   SCIPfreeBufferArray(scip, &vars);

   return (jint) ncliques;
}

/** calculates a partition of the given set of binary variables into negated cliques;
 *  afterwards the output array contains one value for each variable, such that two variables got the same value iff they
 *  were assigned to the same negated clique;
 *  the first variable is always assigned to clique 0 and a variable can only be assigned to clique i if at least one of
 *  the preceding variables was assigned to clique i-1;
 *  for each clique with n_c variables at least n_c-1 variables can be set to TRUE in a feasible solution;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(calcNegatedCliquePartition)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlongArray            jvars,              /**< binary variables in the clique from which at most one can be set to 1 */
   jint                  nvars,              /**< number of variables in the clique */
   jintArray             jcliquepartition    /**< array of length nvars to store the clique partition */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;
   int* cliquepartition;
   int ncliques;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);


   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &cliquepartition, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);
   (*env)->GetIntArrayRegion(env, jcliquepartition, 0, nvars, (jint*)cliquepartition);

   JNISCIP_CALL( SCIPcalcNegatedCliquePartition(scip, vars, (int)nvars, cliquepartition, &ncliques) );

   SCIPfreeBufferArray(scip, &cliquepartition);
   SCIPfreeBufferArray(scip, &vars);

   return (jint) ncliques;
}

/** gets the number of cliques in the clique table
 *
 *  @return number of cliques in the clique table
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jint JNISCIP(getNCliques)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNCliques(scip);
}

/** gets the array of cliques in the clique table
 *
 *  @return array of cliques in the clique table
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jlongArray JNISCIP(getCliques)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int ncliques;

   jlongArray jcliques;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   ncliques = SCIPgetNCliques(scip);
   jcliques = (*env)->NewLongArray(env, ncliques);

   if( jcliques == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_CLIQUE** cliques;

      cliques = SCIPgetCliques(scip);
      (*env)->SetLongArrayRegion(env, jcliques, 0, ncliques, (jlong*)cliques);
   }

   return jcliques;
}

/** Returns whether there is a clique that contains both given variable/value pairs;
 *  the variables must be active binary variables;
 *  if regardimplics is FALSE, only the cliques in the clique table are looked at;
 *  if regardimplics is TRUE, both the cliques and the implications of the implication graph are regarded
 *
 *  @return TRUE, if there is a clique that contains both variable/clique pairs; FALSE, otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note a variable with it's negated variable are NOT! in a clique
 *  @note a variable with itself are in a clique
 */
JNIEXPORT
jboolean JNISCIP(haveVarsCommonClique)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar1,              /**< first variable */
   jboolean              jvalue1,            /**< value of first variable */
   jlong                 jvar2,              /**< second variable */
   jboolean              jvalue2,            /**< value of second variable */
   jboolean              jregardimplics      /**< should the implication graph also be searched for a clique? */
   )
{
   SCIP* scip;
   SCIP_VAR* var1;
   SCIP_VAR* var2;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var1 = (SCIP_VAR*) (size_t) jvar1;
   assert(var1 != NULL);

   var2 = (SCIP_VAR*) (size_t) jvar2;
   assert(var2 != NULL);

   return (jboolean) SCIPhaveVarsCommonClique(scip, var1, (SCIP_Bool)jvalue1, var2, (SCIP_Bool)jvalue2, (SCIP_Bool)jregardimplics);
}

/** sets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(chgVarBranchFactor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jdouble               jbranchfactor       /**< factor to weigh variable's branching score with */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarBranchFactor(scip, var, (jdouble)jbranchfactor) );
}

/** scales the branch factor of the variable with the given value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(scaleVarBranchFactor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jdouble               jscale              /**< factor to scale variable's branching factor with */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPscaleVarBranchFactor(scip, var, (jdouble)jscale) );
}

/** adds the given value to the branch factor of the variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(addVarBranchFactor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jdouble               jaddfactor          /**< value to add to the branch factor of the variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddVarBranchFactor(scip, var, (jdouble)jaddfactor) );
}

/** sets the branch priority of the variable; variables with higher branch priority are always preferred to variables
 *  with lower priority in selection of branching variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 * @note the default branching priority is 0
 */
JNIEXPORT
void JNISCIP(chgVarBranchPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jbranchpriority     /**< branch priority of the variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarBranchPriority(scip, var, (jint)jbranchpriority) );
}

/** changes the branch priority of the variable to the given value, if it is larger than the current priority
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(updateVarBranchPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jbranchpriority     /**< branch priority of the variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPupdateVarBranchPriority(scip, var, (jint)jbranchpriority) );
}

/** adds the given value to the branch priority of the variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(addVarBranchPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jaddpriority        /**< value to add to the branch priority of the variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddVarBranchPriority(scip, var, (jint)jaddpriority) );
}

/** sets the branch direction of the variable (-1: prefer downwards branch, 0: automatic selection, +1: prefer upwards
 *  branch)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(chgVarBranchDirection)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jbranchdirection    /**< preferred branch direction of the variable (downwards, upwards, auto) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarBranchDirection(scip, var, (SCIP_BRANCHDIR)jbranchdirection) );
}

/** changes type of variable in the problem;
 *
 *  @warning This type change might change the variable array returned from SCIPgetVars() and SCIPgetVarsData();
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *
 *  @note If SCIP is already beyond the SCIP_STAGE_PROBLEM and a original variable is passed, the variable type of the
 *        corresponding transformed variable is changed; the type of the original variable does not change
 *
 *  @note If the type changes from a continuous variable to a non-continuous variable the bounds of the variable get
 *        adjusted w.r.t. to integrality information
 */
JNIEXPORT
jboolean JNISCIP(chgVarType)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to change bound for */
   jint                  jvartype            /**< new type of variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_Bool infeasible;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarType(scip, var, (SCIP_VARTYPE)jvartype, &infeasible) );

   return (jboolean) infeasible;
}

/** returns whether aggregation of variables is not allowed */
JNIEXPORT
jboolean JNISCIP(doNotAggr)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPdoNotAggr(scip);
}

/** returns whether variable is not allowed to be multi-aggregated */
JNIEXPORT
jboolean JNISCIP(doNotMultaggrVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable x to aggregate */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jboolean) SCIPdoNotMultaggrVar(scip, var);
}

/** marks the variable that it must not be multi-aggregated
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 */
JNIEXPORT
void JNISCIP(markDoNotMultaggrVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable x to delete */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, var) );
}

/** enables the collection of statistics for a variable
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
void JNISCIP(enableVarHistory)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   SCIPenableVarHistory(scip);
}

/** disables the collection of any statistic for a variable
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
void JNISCIP(disableVarHistory)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   SCIPdisableVarHistory(scip);
}

/** updates the pseudo costs of the given variable and the global pseudo costs after a change of "solvaldelta" in the
 *  variable's solution value and resulting change of "objdelta" in the in the LP's objective value;
 *  the update is ignored, if the objective value difference is infinite
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
void JNISCIP(updateVarPseudocost)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jdouble               jsolvaldelta,       /**< difference of variable's new LP value - old LP value */
   jdouble               jobjdelta,          /**< difference of new LP's objective value - old LP's objective value */
   jdouble               jweight             /**< weight in (0,1] of this update in pseudo cost sum */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPupdateVarPseudocost(scip, var, (SCIP_Real)jsolvaldelta, (SCIP_Real)jobjdelta, (SCIP_Real)jweight) );
}

/** gets the variable's pseudo cost value for the given change of the variable's LP value
 *
 *  @return the variable's pseudo cost value for the given change of the variable's LP value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarPseudocostVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jdouble               jsolvaldelta        /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarPseudocostVal(scip, var, (SCIP_Real)jsolvaldelta);
}

/** gets the variable's pseudo cost value for the given change of the variable's LP value,
 *  only using the pseudo cost information of the current run
 *
 *  @return the variable's pseudo cost value for the given change of the variable's LP value,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarPseudocostValCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jdouble               jsolvaldelta        /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarPseudocostValCurrentRun(scip, var, (SCIP_Real)jsolvaldelta);
}

/** gets the variable's pseudo cost value for the given direction
 *
 *  @return the variable's pseudo cost value for the given direction
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarPseudocost)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jbranchdir          /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarPseudocost(scip, var, (SCIP_BRANCHDIR)jbranchdir);
}

/** gets the variable's pseudo cost value for the given direction,
 *  only using the pseudo cost information of the current run
 *
 *  @return the variable's pseudo cost value for the given direction,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarPseudocostCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jbranchdir          /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarPseudocostCurrentRun(scip, var, (SCIP_BRANCHDIR)jbranchdir);
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction
 *
 *  @return the variable's (possible fractional) number of pseudo cost updates for the given direction
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarPseudocostCount)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jbranchdir          /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarPseudocostCount(scip, var, (SCIP_BRANCHDIR)jbranchdir);
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 *
 *  @return the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarPseudocostCountCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jbranchdir          /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarPseudocostCountCurrentRun(scip, var, (SCIP_BRANCHDIR)jbranchdir);
}

/** gets the variable's pseudo cost score value for the given LP solution value
 *
 *  @return the variable's pseudo cost score value for the given LP solution value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarPseudocostScore)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jdouble               jsolval             /**< variable's LP solution value */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarPseudocostScore(scip, var, (SCIP_Real)jsolval);
}

/** gets the variable's pseudo cost score value for the given LP solution value,
 *  only using the pseudo cost information of the current run
 *
 *  @return the variable's pseudo cost score value for the given LP solution value,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarPseudocostScoreCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jdouble               jsolval             /**< variable's LP solution value */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarPseudocostScoreCurrentRun(scip, var, (SCIP_Real)jsolval);
}

/** returns the variable's VSIDS value
 *
 *  @return the variable's VSIDS value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarVSIDS)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarVSIDS(scip, var, (SCIP_BRANCHDIR)jdir);
}

/** returns the variable's VSIDS value only using conflicts of the current run
 *
 *  @return the variable's VSIDS value only using conflicts of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarVSIDSCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarVSIDSCurrentRun(scip, var, (SCIP_BRANCHDIR)jdir);
}

/** returns the variable's conflict score value
 *
 *  @return the variable's conflict score value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarConflictScore)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< problem variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarConflictScore(scip, var);
}

/** returns the variable's conflict score value only using conflicts of the current run
 *
 *  @return the variable's conflict score value only using conflicts of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarConflictScoreCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< problem variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarConflictScoreCurrentRun(scip, var);
}

/** returns the variable's conflict length score
 *
 *  @return the variable's conflict length score
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarConflictlengthScore)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< problem variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarConflictlengthScore(scip, var);
}

/** returns the variable's conflict length score only using conflicts of the current run
 *
 *  @return the variable's conflict length score only using conflicts of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarConflictlengthScoreCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< problem variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarConflictlengthScoreCurrentRun(scip, var);
}

/** returns the variable's average conflict length
 *
 *  @return the variable's average conflict length
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarAvgConflictlength)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarAvgConflictlength(scip, var, (SCIP_BRANCHDIR)jdir);
}

/** returns the variable's average  conflict length only using conflicts of the current run
 *
 *  @return the variable's average conflict length only using conflicts of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarAvgConflictlengthCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarAvgConflictlengthCurrentRun(scip, var, (SCIP_BRANCHDIR)jdir);
}

/** returns the average number of inferences found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 *
 *  @return the average number of inferences found after branching on the variable in given direction
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarAvgInferences)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarAvgInferences(scip, var, (SCIP_BRANCHDIR)jdir);
}

/** returns the average number of inferences found after branching on the variable in given direction in the current run;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 *
 *  @return the average number of inferences found after branching on the variable in given direction in the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarAvgInferencesCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarAvgInferencesCurrentRun(scip, var, (SCIP_BRANCHDIR)jdir);
}

/** returns the variable's average inference score value
 *
 *  @return the variable's average inference score value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarAvgInferenceScore)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< problem variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarAvgInferenceScore(scip, var);
}

/** returns the variable's average inference score value only using inferences of the current run
 *
 *  @return the variable's average inference score value only using inferences of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarAvgInferenceScoreCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< problem variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarAvgInferenceScoreCurrentRun(scip, var);
}

/** initializes the upwards and downwards pseudocosts, conflict scores, conflict lengths, inference scores, cutoff scores
 *  of a variable to the given values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(initVarBranchStats)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable which should be initialized */
   jdouble               jdownpscost,        /**< value to which pseudocosts for downwards branching should be initialized */
   jdouble               juppscost,          /**< value to which pseudocosts for upwards branching should be initialized */
   jdouble               jdownvsids,         /**< value to which VSIDS score for downwards branching should be initialized */
   jdouble               jupvsids,           /**< value to which VSIDS score for upwards branching should be initialized */
   jdouble               jdownconflen,       /**< value to which conflict length score for downwards branching should be initialized */
   jdouble               jupconflen,         /**< value to which conflict length score for upwards branching should be initialized */
   jdouble               jdowninfer,         /**< value to which inference counter for downwards branching should be initialized */
   jdouble               jupinfer,           /**< value to which inference counter for upwards branching should be initialized */
   jdouble               jdowncutoff,        /**< value to which cutoff counter for downwards branching should be initialized */
   jdouble               jupcutoff           /**< value to which cutoff counter for upwards branching should be initialized */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPinitVarBranchStats(scip, var, (SCIP_Real)jdownpscost, (SCIP_Real)juppscost, (SCIP_Real)jdownvsids, (SCIP_Real)jupvsids, (SCIP_Real)jdownconflen, (SCIP_Real)jupconflen, (SCIP_Real)jdowninfer, (SCIP_Real)jupinfer, (SCIP_Real)jdowncutoff, (SCIP_Real)jupcutoff) );
}

/** returns the average number of cutoffs found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 *
 *  @return the average number of cutoffs found after branching on the variable in given direction
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarAvgCutoffs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarAvgCutoffs(scip, var, (SCIP_BRANCHDIR)jdir);
}

/** returns the average number of cutoffs found after branching on the variable in given direction in the current run;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 *
 *  @return the average number of cutoffs found after branching on the variable in given direction in the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarAvgCutoffsCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jint                  jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarAvgCutoffsCurrentRun(scip, var, (SCIP_BRANCHDIR)jdir);
}

/** returns the variable's average cutoff score value
 *
 *  @return the variable's average cutoff score value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarAvgCutoffScore)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< problem variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarAvgCutoffScore(scip, var);
}

/** returns the variable's average cutoff score value, only using cutoffs of the current run
 *
 *  @return the variable's average cutoff score value, only using cutoffs of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarAvgCutoffScoreCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< problem variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarAvgCutoffScoreCurrentRun(scip, var);
}

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor
 *
 *  @return the variable's average inference/cutoff score value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarAvgInferenceCutoffScore)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jdouble               jcutoffweight       /**< factor to weigh average number of cutoffs in branching score */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarAvgInferenceCutoffScore(scip, var, (SCIP_Real)jcutoffweight);
}

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor, only using inferences and cutoffs of the current run
 *
 *  @return the variable's average inference/cutoff score value, only using inferences and cutoffs of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getVarAvgInferenceCutoffScoreCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jdouble               jcutoffweight       /**< factor to weigh average number of cutoffs in branching score */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarAvgInferenceCutoffScoreCurrentRun(scip, var, (SCIP_Real)jcutoffweight);
}

/** outputs variable information to file stream via the message system
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
JNIEXPORT
void JNISCIP(printVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< problem variable */
   jlong                 jfile               /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPprintVar(scip, var, (FILE*)(size_t)jfile) );
}

/** return TRUE if conflict analysis is applicable; In case the function return FALSE there is no need to initialize the
 *  conflict analysis since it will not be applied
 *
 *  @return return TRUE if conflict analysis is applicable; In case the function return FALSE there is no need to initialize the
 *          conflict analysis since it will not be applied
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
jboolean JNISCIP(isConflictAnalysisApplicable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisConflictAnalysisApplicable(scip);
}

/** initializes the conflict analysis by clearing the conflict candidate queue; this method must be called before you
 *  enter the conflict variables by calling SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(),
 *  SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or SCIPaddConflictBinvar();
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
void JNISCIP(initConflictAnalysis)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPinitConflictAnalysis(scip) );
}

/** adds lower bound of variable at the time of the given bound change index to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictLb() should be called for each lower bound
 *      that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictLb() should be called
 *      for each lower bound, whose current assignment lead to the deduction of the given conflict bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
void JNISCIP(addConflictLb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,             /**< SCIP data structure */
   jlong                 jvar,               /**< variable whose lower bound should be added to conflict candidate queue */
   jlong                 jbdchgidx           /**< bound change index representing time on path to current node, when the
					      *   conflicting bound was valid, NULL for current local bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddConflictLb(scip, var, (SCIP_BDCHGIDX*)(size_t)jbdchgidx) );
}

/** adds lower bound of variable at the time of the given bound change index to the conflict analysis' candidate storage
 *  with the additional information of a relaxed lower bound; this relaxed lower bound is the one which would be enough
 *  to explain a certain bound change;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictRelaxedLb() should be called for each (relaxed) lower bound
 *      that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictRelexedLb() should be called
 *      for each (relaxed) lower bound, whose current assignment lead to the deduction of the given conflict bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
void JNISCIP(addConflictRelaxedLb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable whose lower bound should be added to conflict candidate queue */
   jlong                 jbdchgidx,          /**< bound change index representing time on path to current node, when the
					      *   conflicting bound was valid, NULL for current local bound */
   jdouble               jrelaxedlb          /**< the relaxed lower bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddConflictRelaxedLb(scip, var, (SCIP_BDCHGIDX*)(size_t)jbdchgidx, (SCIP_Real)jrelaxedlb) );
}

/** adds upper bound of variable at the time of the given bound change index to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictUb() should be called for each upper bound that
 *      lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictUb() should be called for
 *      each upper bound, whose current assignment lead to the deduction of the given conflict bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
void JNISCIP(addConflictUb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable whose upper bound should be added to conflict candidate queue */
   jlong                 jbdchgidx           /**< bound change index representing time on path to current node, when the
					      *   conflicting bound was valid, NULL for current local bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddConflictUb(scip, var, (SCIP_BDCHGIDX*)(size_t)jbdchgidx) );
}

/** adds upper bound of variable at the time of the given bound change index to the conflict analysis' candidate storage
 *  with the additional information of a relaxed upper bound; this relaxed upper bound is the one which would be enough
 *  to explain a certain bound change;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictRelaxedUb() should be called for each (relaxed) upper
 *      bound that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictRelaxedUb() should be
 *      called for each (relaxed) upper bound, whose current assignment lead to the deduction of the given conflict
 *      bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
void JNISCIP(addConflictRelaxedUb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable whose upper bound should be added to conflict candidate queue */
   jlong                 jbdchgidx,          /**< bound change index representing time on path to current node, when the
					      *   conflicting bound was valid, NULL for current local bound */
   jdouble               jrelaxedub          /**< the relaxed upper bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddConflictRelaxedUb(scip, var, (SCIP_BDCHGIDX*)(size_t)jbdchgidx, (SCIP_Real)jrelaxedub) );
}

/** adds lower or upper bound of variable at the time of the given bound change index to the conflict analysis' candidate
 *  storage; this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictBd() should be called for each bound
 *      that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictBd() should be called
 *      for each bound, whose current assignment lead to the deduction of the given conflict bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
void JNISCIP(addConflictBd)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable whose upper bound should be added to conflict candidate queue */
   jint                  jboundtype,         /**< the type of the conflicting bound (lower or upper bound) */
   jlong                 jbdchgidx           /**< bound change index representing time on path to current node, when the
					      *   conflicting bound was valid, NULL for current local bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddConflictBd(scip, var, (SCIP_BOUNDTYPE)jboundtype, (SCIP_BDCHGIDX*)(size_t)jbdchgidx) );
}

/** adds lower or upper bound of variable at the time of the given bound change index to the conflict analysis'
 *  candidate storage; with the additional information of a relaxed upper bound; this relaxed upper bound is the one
 *  which would be enough to explain a certain bound change;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictRelaxedBd() should be called for each (relaxed)
 *      bound that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictRelaxedBd() should be
 *      called for each (relaxed) bound, whose current assignment lead to the deduction of the given conflict bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
void JNISCIP(addConflictRelaxedBd)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable whose upper bound should be added to conflict candidate queue */
   jint                  jboundtype,         /**< the type of the conflicting bound (lower or upper bound) */
   jlong                 jbdchgidx,          /**< bound change index representing time on path to current node, when the
					      *   conflicting bound was valid, NULL for current local bound */
   jdouble               jrelaxedbd          /**< the relaxed bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddConflictRelaxedBd(scip, var, (SCIP_BOUNDTYPE)jboundtype, (SCIP_BDCHGIDX*)(size_t)jbdchgidx, (SCIP_Real)jrelaxedbd) );
}

/** adds changed bound of fixed binary variable to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictBinvar() should be called for each fixed binary
 *      variable that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictBinvar() should be called
 *      for each binary variable, whose current fixing lead to the deduction of the given conflict bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
void JNISCIP(addConflictBinvar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< binary variable whose changed bound should be added to conflict queue */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddConflictBinvar(scip, var) );
}

/** checks if the given variable is already part of the current conflict set or queued for resolving with the same or
 *  even stronger bound
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
jboolean JNISCIP(isConflictVarUsed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable whose upper bound should be added to conflict candidate queue */
   jint                  jboundtype,         /**< the type of the conflicting bound (lower or upper bound) */
   jlong                 jbdchgidx           /**< bound change index representing time on path to current node, when the
					      *   conflicting bound was valid, NULL for current local bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_Bool used;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPisConflictVarUsed(scip, var, (SCIP_BOUNDTYPE)jboundtype, (SCIP_BDCHGIDX*)(size_t)jbdchgidx, &used) );

   return (jboolean) used;
}

/** returns the conflict lower bound if the variable is present in the current conflict set; otherwise the global lower
 *  bound
 *
 *  @return returns the conflict lower bound if the variable is present in the current conflict set; otherwise the global lower
 *          bound
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
jdouble JNISCIP(getConflictVarLb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< problem variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetConflictVarLb(scip, var);
}

/** returns the conflict upper bound if the variable is present in the current conflict set; otherwise minus global
 *  upper bound
 *
 *  @return returns the conflict upper bound if the variable is present in the current conflict set; otherwise minus global
 *          upper bound
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
jdouble JNISCIP(getConflictVarUb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< problem variable */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetConflictVarUb(scip, var);
}

/** analyzes conflict bounds that were added after a call to SCIPinitConflictAnalysis() with calls to
 *  SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(),
 *  SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or SCIPaddConflictBinvar(); on success, calls the conflict
 *  handlers to create a conflict constraint out of the resulting conflict set; the given valid depth must be a depth
 *  level, at which the conflict set defined by calls to SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(),
 *  SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), and SCIPaddConflictBinvar() is
 *  valid for the whole subtree; if the conflict was found by a violated constraint, use SCIPanalyzeConflictCons()
 *  instead of SCIPanalyzeConflict() to make sure, that the correct valid depth is used
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
jboolean JNISCIP(analyzeConflict)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  jvaliddepth         /**< minimal depth level at which the initial conflict set is valid */
   )
{
   SCIP* scip;
   SCIP_Bool success;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPanalyzeConflict(scip, (int)jvaliddepth, &success) );

   return (jboolean) success;
}

/** analyzes conflict bounds that were added with calls to SCIPaddConflictLb(), SCIPaddConflictUb(),
 *  SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or
 *  SCIPaddConflictBinvar(); on success, calls the conflict handlers to create a conflict constraint out of the
 *  resulting conflict set; the given constraint must be the constraint that detected the conflict, i.e. the constraint
 *  that is infeasible in the local bounds of the initial conflict set (defined by calls to SCIPaddConflictLb(),
 *  SCIPaddConflictUb(), SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(),
 *  SCIPaddConflictRelaxedBd(), and SCIPaddConflictBinvar())
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
JNIEXPORT
jboolean JNISCIP(analyzeConflictCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint that detected the conflict */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_Bool success;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPanalyzeConflictCons(scip, cons, &success) );

   return (jboolean) success;
}

/** creates and captures a constraint of the given constraint handler
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution may
 *           be declared feasible even if it violates this particular constraint.  This constellation should only be
 *           used, if no LP or pseudo solution can violate the constraint -- e.g. if a local constraint is redundant due
 *           to the variable's local bounds.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method {@link releaseCons()}
 */
JNIEXPORT
jlong JNISCIP(createCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of constraint */
   jlong                 jconshdlr,          /**< constraint handler for this constraint */
   jlong                 jconsdata,          /**< data for this specific constraint */
   jboolean              initial,            /**< should the LP relaxation of constraint be in the initial LP?
					      *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   jboolean              separate,           /**< should the constraint be separated during LP processing?
					      *   Usually set to TRUE. */
   jboolean              enforce,            /**< should the constraint be enforced during node processing?
					      *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   jboolean              check,              /**< should the constraint be checked for feasibility?
					      *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   jboolean              propagate,          /**< should the constraint be propagated during node processing?
					      *   Usually set to TRUE. */
   jboolean              local,              /**< is constraint only valid locally?
					      *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   jboolean              modifiable,         /**< is constraint modifiable (subject to column generation)?
					      *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
					      *   adds coefficients to this constraint. */
   jboolean              dynamic,            /**< is constraint subject to aging?
					      *   Usually set to FALSE. Set to TRUE for own cuts which
					      *   are seperated as constraints. */
   jboolean              removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
					      *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   jboolean              stickingatnode      /**< should the constraint always be kept at the node where it was added, even
					      *   if it may be moved to a more global node?
					      *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   const char* name;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   consdata = (SCIP_CONSDATA*) (size_t) jconsdata;
   assert (consdata != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, NULL);
   if( name == NULL )
      SCIPABORT();

   /* create empty constraint */
   JNISCIP_CALL( SCIPcreateCons(scip, &cons, name, conshdlr, consdata,
	 (SCIP_Bool) initial, (SCIP_Bool) separate, (SCIP_Bool) enforce, (SCIP_Bool) check, (SCIP_Bool) propagate,
	 (SCIP_Bool) local, (SCIP_Bool) modifiable, (SCIP_Bool) dynamic, (SCIP_Bool) removable, (SCIP_Bool) stickingatnode) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong)(size_t)cons;
}

/** decreases usage counter of constraint, and frees memory if necessary */
JNIEXPORT
void JNISCIP(releaseCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/** change constraint name
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
JNIEXPORT
void JNISCIP(chgConsName)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jstring               jname               /**< new name of constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPchgConsName(scip, cons, name) );

   (*env)->ReleaseStringUTFChars(env, jname, name);
}

/** sets the initial flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(setConsInitial)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jboolean              jinitial            /**< new value */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPsetConsInitial(scip, cons, (SCIP_Bool)jinitial) );
}

/** sets the separate flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(setConsSeparated)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jboolean              jseparated          /**< new value */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPsetConsSeparated(scip, cons, (SCIP_Bool)jseparated) );
}

/** sets the enforce flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(setConsEnforced)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jboolean              jenforced           /**< new value */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPsetConsEnforced(scip, cons, (SCIP_Bool)jenforced) );
}

/** sets the check flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(setConsChecked)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jboolean              jcheck              /**< new value */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPsetConsChecked(scip, cons, (SCIP_Bool)jcheck) );
}

/** sets the propagate flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(setConsPropagated)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jboolean              jprop               /**< new value */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPsetConsPropagated(scip, cons, (SCIP_Bool)jprop) );
}

/** sets the local flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(setConsLocal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jboolean              jlocal              /**< new value */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPsetConsLocal(scip, cons, (SCIP_Bool)jlocal) );
}

/** sets the modifiable flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
void JNISCIP(setConsModifiable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jboolean              jmod                /**< new value */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPsetConsModifiable(scip, cons, (SCIP_Bool)jmod) );
}

/** sets the dynamic flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(setConsDynamic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jboolean              jdynamic                /**< new value */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPsetConsDynamic(scip, cons, (SCIP_Bool)jdynamic) );
}

/** sets the removable flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(setConsRemovable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jboolean              jremovable          /**< new value */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPsetConsRemovable(scip, cons, (SCIP_Bool)jremovable) );
}

/** sets the stickingatnode flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(setConsStickingAtNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jboolean              jstick              /**< new value */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPsetConsStickingAtNode(scip, cons, (SCIP_Bool)jstick) );
}

/** gets and captures transformed constraint of a given constraint; if the constraint is not yet transformed,
 *  a new transformed constraint for this constraint is created
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jlong JNISCIP(transformCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint to get/create transformed constraint for */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_CONS* transcons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPtransformCons(scip, cons, &transcons) );

   return (jlong) (size_t) transcons;
}

/** gets and captures transformed constraints for an array of constraints;
 *  if a constraint in the array is not yet transformed, a new transformed constraint for this constraint is created;
 *  it is possible to call this method with conss == transconss
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlongArray JNISCIP(transformConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  nconss,             /**< number of constraints to get/create transformed constraints for */
   jlongArray            jconss              /**< array with constraints to get/create transformed constraints for */
   )
{
   SCIP* scip;
   SCIP_CONS** conss;

   jlongArray jtransconss;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &conss, (int)nconss) );

   (*env)->GetLongArrayRegion(env, jconss, 0, (int)nconss, (jlong*)conss);
   jtransconss = (*env)->NewLongArray(env, (int)nconss);

   if (jtransconss == NULL)
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_CONS** transconss;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &transconss, (int)nconss) );

      JNISCIP_CALL( SCIPtransformConss(scip, (int)nconss, conss, transconss) );
      (*env)->SetLongArrayRegion(env, jtransconss, 0, (int)nconss, (jlong*)transconss);

      SCIPfreeBufferArray(scip, &transconss);
   }

   SCIPfreeBufferArray(scip, &conss);

   return jtransconss;
}

/** gets corresponding transformed constraint of a given constraint;
 *  returns NULL as transcons, if transformed constraint is not yet existing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
jlong JNISCIP(getTransformedCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint to get transformed constraint for */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_CONS* transcons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPgetTransformedCons(scip, cons, &transcons) );

   return (jlong) (size_t) transcons;
}
/** gets corresponding transformed constraints for an array of constraints;
 *  stores NULL in a transconss slot, if the transformed constraint is not yet existing;
 *  it is possible to call this method with conss == transconss, but remember that constraints that are not
 *  yet transformed will be replaced with NULL
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jlongArray JNISCIP(getTransformedConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  nconss,             /**< number of constraints to get the transformed constraints for */
   jlongArray            jconss              /**< constraints to get the transformed constraints for */
   )
{
   SCIP* scip;
   SCIP_CONS** conss;

   jlongArray jtransconss;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &conss, (int)nconss) );

   (*env)->GetLongArrayRegion(env, jconss, 0, (int)nconss, (jlong*)conss);
   jtransconss = (*env)->NewLongArray(env, (int)nconss);

   if (jtransconss == NULL)
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_CONS** transconss;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &transconss, (int)nconss) );

      JNISCIP_CALL( SCIPgetTransformedConss(scip, (int)nconss, conss, transconss) );
      (*env)->SetLongArrayRegion(env, jtransconss, 0, (int)nconss, (jlong*)transconss);

      SCIPfreeBufferArray(scip, &transconss);
   }

   SCIPfreeBufferArray(scip, &conss);

   return jtransconss;
}

/** adds given value to age of constraint, but age can never become negative;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void JNISCIP(addConsAge)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jdouble               jdeltaage           /**< value to add to the constraint's age */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPaddConsAge(scip, cons, (SCIP_Real)jdeltaage) );
}

/** increases age of constraint by 1.0;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void JNISCIP(incConsAge)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPincConsAge(scip, cons) );
}

/** resets age of constraint to zero;
 *  should be called
 *   - in constraint separation, if a cut was found for this constraint,
 *   - in constraint enforcing, if the constraint was violated, and
 *   - in constraint propagation, if a domain reduction was deduced;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void JNISCIP(resetConsAge)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPresetConsAge(scip, cons) );
}

/** enables constraint's separation, propagation, and enforcing capabilities
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void JNISCIP(enableCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPenableCons(scip, cons) );
}

/** disables constraint's separation, propagation, and enforcing capabilities, s.t. the constraint is not propagated,
 *  separated, and enforced anymore until it is enabled again with a call to SCIPenableCons();
 *  in contrast to SCIPdelConsLocal() and SCIPdelConsNode(), the disabling is not associated to a node in the tree and
 *  does not consume memory; therefore, the constraint is neither automatically enabled on leaving the node nor
 *  automatically disabled again on entering the node again;
 *  note that the constraints enforcing capabilities are necessary for the solution's feasibility, if the constraint
 *  is a model constraint; that means, you must be sure that the constraint cannot be violated in the current subtree,
 *  and you have to enable it again manually by calling SCIPenableCons(), if this subtree is left (e.g. by using
 *  an appropriate event handler that watches the corresponding variables' domain changes)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void JNISCIP(disableCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPdisableCons(scip, cons) );
}

/** enables constraint's separation capabilities
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void JNISCIP(enableConsSeparation)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPenableConsSeparation(scip, cons) );
}

/** disables constraint's separation capabilities s.t. the constraint is not propagated anymore until the separation
 *  is enabled again with a call to SCIPenableConsSeparation(); in contrast to SCIPdelConsLocal() and SCIPdelConsNode(),
 *  the disabling is not associated to a node in the tree and does not consume memory; therefore, the constraint
 *  is neither automatically enabled on leaving the node nor automatically disabled again on entering the node again
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void JNISCIP(disableConsSeparation)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPenableConsSeparation(scip, cons) );
}

/** enables constraint's propagation capabilities
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void JNISCIP(enableConsPropagation)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPenableConsPropagation(scip, cons) );
}

/** disables constraint's propagation capabilities s.t. the constraint is not propagated anymore until the propagation
 *  is enabled again with a call to SCIPenableConsPropagation(); in contrast to SCIPdelConsLocal() and SCIPdelConsNode(),
 *  the disabling is not associated to a node in the tree and does not consume memory; therefore, the constraint
 *  is neither automatically enabled on leaving the node nor automatically disabled again on entering the node again
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void JNISCIP(disableConsPropagation)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPdisableConsPropagation(scip, cons) );
}

/** adds given values to lock status of the constraint and updates the rounding locks of the involved variables
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
void JNISCIP(addConsLocks)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jint                  nlockspos,          /**< increase in number of rounding locks for constraint */
   jint                  nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPaddConsLocks(scip, cons, (int)nlockspos, (int)nlocksneg) );
}

/** checks single constraint for feasibility of the given solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
jint JNISCIP(checkCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jlong                 jsol,               /**< primal CIP solution */
   jboolean              checkintegrality,   /**< has integrality to be checked? */
   jboolean              checklprows,        /**< have current LP rows (both local and global) to be checked? */
   jboolean              printreason         /**< should the reason for the violation be printed? */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_RESULT result;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPcheckCons(scip, cons, sol, (SCIP_Bool)checkintegrality, (SCIP_Bool)checklprows, (SCIP_Bool)printreason, &result) );

   return (jint) result;
}

/** enforces single constraint for a given pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.  It may only be called for constraints that were not
 *        added to SCIP beforehand.
 */
jint JNISCIP(enfopsCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint to enforce */
   jboolean              solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   jboolean              objinfeasible       /**< is the solution infeasible anyway due to violating lower objective bound? */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_RESULT result;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPenfopsCons(scip, cons,  (SCIP_Bool)solinfeasible, (SCIP_Bool)objinfeasible, &result) );

   return (jint) result;
}

/** enforces single constraint for a given LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.  It may only be called for constraints that were not
 *        added to SCIP beforehand.
 */
jint JNISCIP(enfolpCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint to enforce */
   jboolean              solinfeasible       /**< was the solution already declared infeasible by a constraint handler? */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_RESULT result;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPenfolpCons(scip, cons,  (SCIP_Bool)solinfeasible, &result) );

   return (jint) result;
}

/** calls LP initialization method for single constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.  It may only be called for constraints that were not
 *        added to SCIP beforehand.
 */
void JNISCIP(initlpCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPinitlpCons(scip, cons) );
}

/** calls separation method of single constraint for LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.
 */
jint JNISCIP(sepalpCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_RESULT result;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPsepalpCons(scip, cons, &result) );

   return (jint) result;
}

/** calls separation method of single constraint for given primal solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.
 */
jint JNISCIP(sepasolCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint to separate */
   jlong                 jsol                /**< primal solution that should be separated*/
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_RESULT result;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPsepasolCons(scip, cons, sol, &result) );

   return (jint) result;
}

/** calls domain propagation method of single constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.
 */
jint JNISCIP(propCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint */
   jint                  proptiming          /**< current point in the node solving loop */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_RESULT result;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);


   JNISCIP_CALL( SCIPpropCons(scip, cons, (SCIP_PROPTIMING)proptiming, &result) );

   return (jint) result;
}

/** resolves propagation conflict of single constraint
 *
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.  It may only be called for constraints that were not
 *        added to SCIP beforehand.
 */
jint JNISCIP(respropCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint to resolve conflict for */
   jlong                 jvar,           /**< the conflict variable whose bound change has to be resolved */
   int                   inferinfo,          /**< the user information passed to the corresponding SCIPinferVarLbCons() or SCIPinferVarUbCons() call */
   jint                  boundtype,          /**< the type of the changed bound (lower or upper bound) */
   jlong                 jbdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   jdouble               relaxedbd           /**< the relaxed bound which is sufficient to be explained */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR* var;
   SCIP_BDCHGIDX* bdchgidx;
   SCIP_RESULT result;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bdchgidx = (SCIP_BDCHGIDX*) (size_t) jbdchgidx;
   assert(bdchgidx != NULL);

   JNISCIP_CALL( SCIPrespropCons(scip, cons, var, (int)inferinfo, (SCIP_BOUNDTYPE)boundtype, bdchgidx, (SCIP_Real)relaxedbd, &result) );

   return (jint) result;
}

/** calls constraint activation notification method of single constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *
 *  @note This is an advanced method and should be used with caution.  It may only be called for constraints that were not
 *      added to SCIP beforehand.
 */
void JNISCIP(activeCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint to notify*/
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPactiveCons(scip, cons) );
}

/** calls constraint deactivation notification method of single constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution. It may only be called for constraints that were not
 *        added to SCIP beforehand.
 */
void JNISCIP(deactiveCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons               /**< constraint to notify*/
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPdeactiveCons(scip, cons) );
}

/** outputs constraint information to file stream via the message handler system
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
void JNISCIP(printCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcons,              /**< constraint to notify*/
   jlong                 jfile               /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cons = (SCIP_CONS*) (size_t) jcons;
   assert(cons != NULL);

   JNISCIP_CALL( SCIPprintCons(scip, cons, (FILE*)(size_t)jfile) );
}

/** returns, whether the LP was or is to be solved in the current node
 *
 *  @return whether the LP was or is to be solved in the current node.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jboolean JNISCIP(hasCurrentNodeLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPhasCurrentNodeLP(scip);
}

/** returns, whether the LP of the current node is already constructed
 *
 *  @return whether the LP of the current node is already constructed.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jboolean JNISCIP(isLPConstructed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisLPConstructed(scip);
}

/** makes sure that the LP of the current node is loaded and may be accessed through the LP information methods
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jboolean JNISCIP(constructLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Bool cutoff;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   return (jboolean) cutoff;
}

/** makes sure that the LP of the current node is flushed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
void JNISCIP(flushLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPflushLP(scip) );
}

/** gets solution status of current LP
 *
 *  @return the solution status of current LP.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jint JNISCIP(getLPSolstat)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetLPSolstat(scip);
}

/** returns whether the current lp is a relaxation of the current problem and its optimal objective value is a local lower bound
 *
 *  @return whether the current lp is a relaxation of the current problem and its optimal objective value is a local lower bound.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jboolean JNISCIP(isLPRelax)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisLPRelax(scip);
}

/** gets objective value of current LP (which is the sum of column and loose objective value)
 *
 *  @return the objective value of current LP (which is the sum of column and loose objective value).
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jdouble JNISCIP(getLPObjval)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetLPObjval(scip);
}

/** gets part of objective value of current LP that results from COLUMN variables only
 *
 *  @return the part of objective value of current LP that results from COLUMN variables only.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jdouble JNISCIP(getLPColumnObjval)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetLPColumnObjval(scip);
}

/** gets part of objective value of current LP that results from LOOSE variables only
 *
 *  @return part of objective value of current LP that results from LOOSE variables only.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jdouble JNISCIP(getLPLooseObjval)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetLPLooseObjval(scip);
}

/** gets the global pseudo objective value; that is all variables set to their best  (w.r.t. the objective
 *  function) global bound
 *
 *  @return the global pseudo objective value; that is all variables set to their best  (w.r.t. the objective
 *  function) global bound.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jdouble JNISCIP(getGlobalPseudoObjval)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetGlobalPseudoObjval(scip);
}

/** gets the pseudo objective value for the current search node; that is all variables set to their best (w.r.t. the
 *  objective function) local bound
 *
 *  @return the pseudo objective value for the current search node; that is all variables set to their best (w.r.t. the
 *  objective function) local bound.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jdouble JNISCIP(getPseudoObjval)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetPseudoObjval(scip);
}

/** returns whether the root lp is a relaxation of the problem and its optimal objective value is a global lower bound
 *
 *  @return whether the root lp is a relaxation of the problem and its optimal objective value is a global lower bound.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jboolean JNISCIP(isRootLPRelax)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisRootLPRelax(scip);
}

/** gets the objective value of the root node LP or SCIP_INVALID if the root node LP was not (yet) solved
 *
 *  @return the objective value of the root node LP or SCIP_INVALID if the root node LP was not (yet) solved.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jdouble JNISCIP(getLPRootObjval)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetLPRootObjval(scip);
}

/** gets part of the objective value of the root node LP that results from COLUMN variables only;
 *  returns SCIP_INVALID if the root node LP was not (yet) solved
 *
 *  @return the part of the objective value of the root node LP that results from COLUMN variables only;
 *  or SCIP_INVALID if the root node LP was not (yet) solved.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jdouble JNISCIP(getLPRootColumnObjval)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetLPRootColumnObjval(scip);
}

/** gets part of the objective value of the root node LP that results from LOOSE variables only;
 *  returns SCIP_INVALID if the root node LP was not (yet) solved
 *
 *  @return the part of the objective value of the root node LP that results from LOOSE variables only;
 *  or SCIP_INVALID if the root node LP was not (yet) solved.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jdouble JNISCIP(getLPRootLooseObjval)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetLPRootLooseObjval(scip);
}

/** gets current LP columns
 *
 *  @return the current LP columns.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jlongArray JNISCIP(getLPCols)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int ncols;

   jlongArray jcols;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   ncols = SCIPgetNLPCols(scip);
   jcols = (*env)->NewLongArray(env, ncols);

   if( jcols == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_COL** cols;

      cols = SCIPgetLPCols(scip);
      (*env)->SetLongArrayRegion(env, jcols, 0, ncols, (jlong*)cols);
   }

   return jcols;
}

/** gets current number of LP columns
 *
 *  @return the current number of LP columns.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jint JNISCIP(getNLPCols)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNLPCols(scip);
}

/** gets current LP rows
 *
 *  @return the current LP rows.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jlongArray JNISCIP(getLPRows)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nrows;

   jlongArray jrows;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nrows = SCIPgetNLPRows(scip);
   jrows = (*env)->NewLongArray(env, nrows);

   if( jrows == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_ROW** rows;

      rows = SCIPgetLPRows(scip);
      (*env)->SetLongArrayRegion(env, jrows, 0, nrows, (jlong*)rows);
   }

   return jrows;
}

/** gets current number of LP rows
 *
 *  @return the current number of LP rows.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jint JNISCIP(getNLPRows)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNLPRows(scip);
}

/** returns TRUE iff all columns, i.e. every variable with non-empty column w.r.t. all ever created rows, are present
 *  in the LP, and FALSE, if there are additional already existing columns, that may be added to the LP in pricing
 *
 *  @return TRUE iff all columns, i.e. every variable with non-empty column w.r.t. all ever created rows, are present
 *  in the LP, and FALSE, if there are additional already existing columns, that may be added to the LP in pricing.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jboolean JNISCIP(allColsInLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPallColsInLP(scip);
}

/** returns whether the current LP solution is basic, i.e. is defined by a valid simplex basis
 *
 *  @return whether the current LP solution is basic, i.e. is defined by a valid simplex basis.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jboolean JNISCIP(isLPSolBasic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisLPSolBasic(scip);
}

/** gets all indices of basic columns and rows: index i >= 0 corresponds to column i, index i < 0 to row -i-1
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jint JNISCIP(getLPBasisInd)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int basisind;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPgetLPBasisInd(scip, &basisind) );

   return (jint) basisind;
}

/** writes current LP to a file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
void JNISCIP(writeLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jfilename           /**< file name */
   )
{
   SCIP* scip;
   const char* filename;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   if( jfilename != NULL )
      filename = (*env)->GetStringUTFChars(env, jfilename, &iscopy);
   else
      filename = NULL;

   JNISCIP_CALL( SCIPwriteLP(scip, filename) );
}

/** writes MIP relaxation of the current branch-and-bound node to a file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
void JNISCIP(writeMIP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jfilename,           /**< file name */
   jboolean              genericnames,       /**< should generic names like x_i and row_j be used in order to avoid
					      *   troubles with reserved symbols? */
   jboolean              origobj,            /**< should the original objective function be used? */
   jboolean              lazyconss           /**< output removable rows as lazy constraints? */
   )
{
   SCIP* scip;
   const char* filename;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   if( jfilename != NULL )
      filename = (*env)->GetStringUTFChars(env, jfilename, &iscopy);
   else
      filename = NULL;

   JNISCIP_CALL( SCIPwriteMIP(scip, filename, (SCIP_Bool)genericnames, (SCIP_Bool)origobj, (SCIP_Bool)lazyconss) );
}

/** gets the LP interface of SCIP;
 *  with the LPI you can use all of the methods defined in lpi/lpi.h;
 *
 *  @warning You have to make sure, that the full internal state of the LPI does not change or is recovered completely
 *           after the end of the method that uses the LPI. In particular, if you manipulate the LP or its solution
 *           (e.g. by calling one of the SCIPlpiAdd...() or one of the SCIPlpiSolve...() methods), you have to check in
 *           advance with SCIPlpiWasSolved() whether the LP is currently solved. If this is the case, you have to make
 *           sure, the internal solution status is recovered completely at the end of your method. This can be achieved
 *           by getting the LPI state before applying any LPI manipulations with SCIPlpiGetState() and restoring it
 *           afterwards with SCIPlpiSetState() and SCIPlpiFreeState(). Additionally you have to resolve the LP with the
 *           appropriate SCIPlpiSolve...() call in order to reinstall the internal solution status.
 *
 *  @warning Make also sure, that all parameter values that you have changed are set back to their original values.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
jlong JNISCIP(getLPI)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_LPI* lpi;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPgetLPI(scip, &lpi) );

   return (jlong) (size_t) lpi;
}

/** Displays quality information about the current LP solution. An LP solution need to be available. Information printed
 *  is subject to what the LP solver supports
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_FREE
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note The printing process is done via the message handler system.
 */
void JNISCIP(printLPSolutionQuality)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPprintLPSolutionQuality(scip, (FILE*)(size_t)file) );
}

/** TODO: computeLPRelIntPoint length = 2? */
/** compute relative interior point to current LP
 *  @see SCIPlpComputeRelIntPoint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
#if 0
void JNISCIP(computeLPRelIntPoint)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jboolean              relaxrows,          /**< should the rows be relaxed */
   jboolean              inclobjcutoff,      /**< should a row for the objective cutoff be included */
   jchar                 normtype,           /**< which norm to use: 'o'ne-norm or 's'upremum-norm */
   jdouble               timelimit,          /**< time limit for LP solver */
   jint                  iterlimit,          /**< iteration limit for LP solver */
   jlongArray            jpoint              /**< relative interior point on exit */
   )
{
   SCIP* scip;
   SCIP_SOL** point;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &point, ?length?) );

   (*env)->GetLongArrayRegion(env, jpoint, 0, ?length?, (jlong*)(*point));

   JNISCIP_CALL( SCIPcomputeLPRelIntPoint(scip, (SCIP_Bool)relaxrows, (SCIP_Bool)inclobjcutoff, (char)normtype, (SCIP_Real)timelimit, (int)iterlimit, point) );

   SCIPfreeBufferArray(scip, &point);

}
#endif

/** returns the reduced costs of a column in the last (feasible) LP
 *
 *  @return the reduced costs of a column in the last (feasible) LP
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getColRedcost)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP* scip;
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jdouble) SCIPgetColRedcost(scip, col);
}

/** returns the Farkas coefficient of a column in the last (infeasible) LP
 *
 *  @return the Farkas coefficient of a column in the last (infeasible) LP
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getColFarkasCoef)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP* scip;
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jdouble) SCIPgetColFarkasCoef(scip, col);
}

/** creates and captures an LP row from a constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jlong JNISCIP(createRowCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jconshdlr,          /**< constraint handler that creates the row */
   jstring               jname,              /**< name of row */
   jint                  len,                /**< number of nonzeros in the row */
   jlongArray            jcols,              /**< array with columns of row entries */
   jdoubleArray          jvals,              /**< array with coefficients of row entries */
   jdouble               lhs,                /**< left hand side of row */
   jdouble               rhs,                /**< right hand side of row */
   jboolean              local,              /**< is row only valid locally? */
   jboolean              modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   jboolean              removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP* scip;
   SCIP_ROW* row;
   SCIP_CONSHDLR* conshdlr;
   SCIP_COL** cols;
   SCIP_Real* vals;
   const char* name;

   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &cols, (int)len) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)len) );

   (*env)->GetLongArrayRegion(env, jcols, 0, (int)len, (jlong*)cols);
   (*env)->GetDoubleArrayRegion(env, jvals, 0, (int)len, (jdouble*)vals);

   JNISCIP_CALL( SCIPcreateRowCons(scip, &row, conshdlr, name, (int)len, cols, vals, (SCIP_Real)lhs, (SCIP_Real)rhs, (SCIP_Bool)local, (SCIP_Bool)modifiable, (SCIP_Bool)removable) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &cols);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) row;
}

/** creates and captures an LP row from a separator
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jlong JNISCIP(createRowSepa)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsepa,              /**< separator that creates the row */
   jstring               jname,              /**< name of row */
   jint                  len,                /**< number of nonzeros in the row */
   jlongArray            jcols,              /**< array with columns of row entries */
   jdoubleArray          jvals,              /**< array with coefficients of row entries */
   jdouble               lhs,                /**< left hand side of row */
   jdouble               rhs,                /**< right hand side of row */
   jboolean              local,              /**< is row only valid locally? */
   jboolean              modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   jboolean              removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP* scip;
   SCIP_ROW* row;
   SCIP_SEPA* sepa;
   SCIP_COL** cols;
   SCIP_Real* vals;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sepa = (SCIP_SEPA*) (size_t) jsepa;
   assert(sepa != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &cols, (int)len) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)len) );

   (*env)->GetLongArrayRegion(env, jcols, 0, (int)len, (jlong*)cols);
   (*env)->GetDoubleArrayRegion(env, jvals, 0, (int)len, (jdouble*)vals);

   JNISCIP_CALL( SCIPcreateRowSepa(scip, &row, sepa, name, (int)len, cols, vals, (SCIP_Real)lhs, (SCIP_Real)rhs, (SCIP_Bool)local, (SCIP_Bool)modifiable, (SCIP_Bool)removable) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &cols);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) row;
}

/** creates and captures an LP row from an unspecified source
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jlong JNISCIP(createRowUnspec)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of row */
   jint                  len,                /**< number of nonzeros in the row */
   jlongArray            jcols,              /**< array with columns of row entries */
   jdoubleArray          jvals,              /**< array with coefficients of row entries */
   jdouble               lhs,                /**< left hand side of row */
   jdouble               rhs,                /**< right hand side of row */
   jboolean              local,              /**< is row only valid locally? */
   jboolean              modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   jboolean              removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP* scip;
   SCIP_ROW* row;
   SCIP_COL** cols;
   SCIP_Real* vals;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &cols, (int)len) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)len) );

   (*env)->GetLongArrayRegion(env, jcols, 0, (int)len, (jlong*)cols);
   (*env)->GetDoubleArrayRegion(env, jvals, 0, (int)len, (jdouble*)vals);

   JNISCIP_CALL( SCIPcreateRowUnspec(scip, &row, name, (int)len, cols, vals, (SCIP_Real)lhs, (SCIP_Real)rhs, (SCIP_Bool)local, (SCIP_Bool)modifiable, (SCIP_Bool)removable) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &cols);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) row;
}

/** creates and captures an LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @deprecated Please use SCIPcreateRowCons() or SCIPcreateRowSepa() when calling from a constraint handler or separator in order
 *              to facilitate correct statistics. If the call is from neither a constraint handler or separator, use SCIPcreateRowUnspec().
 */
jlong JNISCIP(createRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of row */
   jint                  len,                /**< number of nonzeros in the row */
   jlongArray            jcols,              /**< array with columns of row entries */
   jdoubleArray          jvals,              /**< array with coefficients of row entries */
   jdouble               lhs,                /**< left hand side of row */
   jdouble               rhs,                /**< right hand side of row */
   jboolean              local,              /**< is row only valid locally? */
   jboolean              modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   jboolean              removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIPerrorMessage("method createRow is not implemented yet (deprecated)\n");
   JNISCIP_CALL( SCIP_ERROR );

   return 0;
#if 0
   SCIP* scip;
   SCIP_ROW* row;
   SCIP_COL** cols;
   SCIP_Real* vals;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &cols, (int)len) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)len) );

   (*env)->GetLongArrayRegion(env, jcols, 0, (int)len, (jlong*)cols);
   (*env)->GetDoubleArrayRegion(env, jvals, 0, (int)len, (jdouble*)vals);

   JNISCIP_CALL( SCIPcreateRow(scip, &row, name, (int)len, cols, vals, (SCIP_Real)lhs, (SCIP_Real)rhs, (SCIP_Bool)local, (SCIP_Bool)modifiable, (SCIP_Bool)removable) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &cols);

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) row;
#endif
}

/** creates and captures an LP row without any coefficients from a constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jlong JNISCIP(createEmptyRowCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jconshdlr,          /**< constraint handler that creates the row */
   jstring               jname,              /**< name of row */
   jdouble               lhs,                /**< left hand side of row */
   jdouble               rhs,                /**< right hand side of row */
   jboolean              local,              /**< is row only valid locally? */
   jboolean              modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   jboolean              removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP* scip;
   SCIP_ROW* row;
   SCIP_CONSHDLR* conshdlr;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   conshdlr = (SCIP_CONSHDLR*) (size_t) jconshdlr;
   assert(conshdlr != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, name, (SCIP_Real)lhs, (SCIP_Real)rhs, (SCIP_Bool)local, (SCIP_Bool)modifiable, (SCIP_Bool)removable) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) row;
}

/** creates and captures an LP row without any coefficients from a separator
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jlong JNISCIP(createEmptyRowSepa)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsepa,               /**< separator that creates the row */
   jstring               jname,               /**< name of row */
   jdouble               lhs,                /**< left hand side of row */
   jdouble               rhs,                /**< right hand side of row */
   jboolean              local,              /**< is row only valid locally? */
   jboolean              modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   jboolean              removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP* scip;
   SCIP_SEPA* sepa;
   SCIP_ROW* row;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sepa = (SCIP_SEPA*) (size_t) jsepa;
   assert(sepa != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, name, (SCIP_Real)lhs, (SCIP_Real)rhs, (SCIP_Bool)local, (SCIP_Bool)modifiable, (SCIP_Bool)removable) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) row;
}

/** creates and captures an LP row without any coefficients from an unspecified source
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jlong JNISCIP(createEmptyRowUnspec)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,               /**< name of row */
   jdouble               lhs,                /**< left hand side of row */
   jdouble               rhs,                /**< right hand side of row */
   jboolean              local,              /**< is row only valid locally? */
   jboolean              modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   jboolean              removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP* scip;
   SCIP_ROW* row;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, name, (SCIP_Real)lhs, (SCIP_Real)rhs, (SCIP_Bool)local, (SCIP_Bool)modifiable, (SCIP_Bool)removable) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) row;
}

/** creates and captures an LP row without any coefficients
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @deprecated Please use SCIPcreateEmptyRowCons() or SCIPcreateEmptyRowSepa() when calling from a constraint handler or separator in order
 *              to facilitate correct statistics. If the call is from neither a constraint handler or separator, use SCIPcreateEmptyRowUnspec().
 */
jlong JNISCIP(createEmptyRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,               /**< name of row */
   jdouble               lhs,                /**< left hand side of row */
   jdouble               rhs,                /**< right hand side of row */
   jboolean              local,              /**< is row only valid locally? */
   jboolean              modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   jboolean              removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP* scip;
   SCIP_ROW* row;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, name, (SCIP_Real)lhs, (SCIP_Real)rhs, (SCIP_Bool)local, (SCIP_Bool)modifiable, (SCIP_Bool)removable) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) row;
}

/** increases usage counter of LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(captureRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< row to capture */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPcaptureRow(scip, row) );
}

/** decreases usage counter of LP row, and frees memory if necessary
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(releaseRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPreleaseRow(scip, &row) );
}

/** changes left hand side of LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(chgRowLhs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow,               /**< LP row */
   jdouble               lhs                 /**< new left hand side */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPchgRowLhs(scip, row, (SCIP_Real)lhs) );
}

/** changes right hand side of LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(chgRowRhs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow,               /**< LP row */
   jdouble               rhs                 /**< new right hand side */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPchgRowRhs(scip, row, (SCIP_Real)rhs) );
}

/** informs row, that all subsequent additions of variables to the row should be cached and not directly applied;
 *  after all additions were applied, SCIPflushRowExtensions() must be called;
 *  while the caching of row extensions is activated, information methods of the row give invalid results;
 *  caching should be used, if a row is build with SCIPaddVarToRow() calls variable by variable to increase
 *  the performance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(cacheRowExtensions)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPcacheRowExtensions(scip, row) );
}

/** flushes all cached row extensions after a call of SCIPcacheRowExtensions() and merges coefficients with
 *  equal columns into a single coefficient
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(flushRowExtensions)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPflushRowExtensions(scip, row) );
}

/** resolves variable to columns and adds them with the coefficient to the row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note In case calling this method in the enforcement process of an lp solution, it might be that some variables,
 *        that were not yet in the LP (e.g. dynamic columns) will change there lp solution value returned by SCIP.
 *
 *        e.g. A variable, which has a negative objective value, that has no column in the lp yet, is in the lp solution
 *        on its upper bound (variables with status SCIP_VARSTATUS_LOOSE are in an lp solution on it's best bound), but
 *        creating the column, changes the solution value (variable than has status SCIP_VARSTATUS_COLUMN, and the
 *        initialization sets the lp solution value) to 0.0 . ( This leads to the conclusion that, if a constraint was
 *        violated, the linear relaxation might not be violated anymore.)
 */
void JNISCIP(addVarToRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow,               /**< LP row */
   jlong                 jvar,               /**< problem variable */
   jdouble               val                 /**< value of coefficient */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddVarToRow(scip, row, var, (SCIP_Real)val) );
}

/** resolves variables to columns and adds them with the coefficients to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(addVarsToRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow,               /**< LP row */
   jint                  nvars,              /**< number of variables to add to the row */
   jlongArray            jvars,              /**< problem variables to add */
   jdoubleArray          jvals               /**< values of coefficient */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)vars);
   (*env)->GetDoubleArrayRegion(env, jvals, 0, (int)nvars, (jdouble*)vals);

   JNISCIP_CALL( SCIPaddVarsToRow(scip, row, (int)nvars, vars, vals) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
}

/** resolves variables to columns and adds them with the same single coefficient to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(addVarsToRowSameCoef)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow,               /**< LP row */
   jint                  nvars,              /**< number of variables to add to the row */
   jlongArray            jvars,              /**< problem variables to add */
   jdouble               val                 /**< unique value of all coefficients */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)vars);

   JNISCIP_CALL( SCIPaddVarsToRowSameCoef(scip, row, (int)nvars, vars, (SCIP_Real)val) );

   SCIPfreeBufferArray(scip, &vars);
}

/** tries to scale row, s.t. all coefficients (of integer variables) become integral
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jboolean JNISCIP(makeRowIntegral)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow,               /**< LP row */
   jdouble               mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   jdouble               maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   jlong                 maxdnom,            /**< maximal denominator allowed in rational numbers */
   jdouble               maxscale,           /**< maximal value to scale row with */
   jboolean              usecontvars         /**< should the coefficients of the continuous variables also be made integral? */
   )
{
   SCIP* scip;
   SCIP_ROW* row;
   SCIP_Bool success;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPmakeRowIntegral(scip, row, (SCIP_Real)mindelta, (SCIP_Real)maxdelta, (SCIP_Longint)maxdnom, (SCIP_Real)maxscale, (SCIP_Bool)usecontvars, &success) );

   return (jboolean) success;
}

/** returns minimal absolute value of row vector's non-zero coefficients
 *
 *  @return minimal absolute value of row vector's non-zero coefficients
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getRowMinCoef)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIPgetRowMinCoef(scip, row);
}

/** returns maximal absolute value of row vector's non-zero coefficients
 *
 *  @return maximal absolute value of row vector's non-zero coefficients
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getRowMaxCoef)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIPgetRowMaxCoef(scip, row);
}

/** returns the minimal activity of a row w.r.t. the column's bounds
 *
 *  @return the minimal activity of a row w.r.t. the column's bounds
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getRowMinActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIPgetRowMinActivity(scip, row);
}

/** returns the maximal activity of a row w.r.t. the column's bounds
 *
 *  @return the maximal activity of a row w.r.t. the column's bounds
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getRowMaxActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIPgetRowMaxActivity(scip, row);
}

/** recalculates the activity of a row in the last LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(recalcRowLPActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPrecalcRowLPActivity(scip, row) );
}

/** returns the activity of a row in the last LP solution
 *
 *  @return activity of a row in the last LP solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getRowLPActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIPgetRowLPActivity(scip, row);
}

/** returns the feasibility of a row in the last LP solution
 *
 *  @return the feasibility of a row in the last LP solution: negative value means infeasibility
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getRowLPFeasibility)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIPgetRowLPFeasibility(scip, row);
}

/** recalculates the activity of a row for the current pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(recalcRowPseudoActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPrecalcRowPseudoActivity(scip, row) );
}

/** returns the activity of a row for the current pseudo solution
 *
 *  @return the activity of a row for the current pseudo solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getRowPseudoActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIPgetRowPseudoActivity(scip, row);
}

/** returns the feasibility of a row for the current pseudo solution: negative value means infeasibility
 *
 *  @return the feasibility of a row for the current pseudo solution: negative value means infeasibility
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getRowPseudoFeasibility)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIPgetRowPseudoFeasibility(scip, row);
}

/** recalculates the activity of a row in the last LP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(recalcRowActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPrecalcRowActivity(scip, row) );
}

/** returns the activity of a row in the last LP or pseudo solution
 *
 *  @return the activity of a row in the last LP or pseudo solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getRowActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIPgetRowActivity(scip, row);
}

/** returns the feasibility of a row in the last LP or pseudo solution
 *
 *  @return the feasibility of a row in the last LP or pseudo solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getRowFeasibility)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIPgetRowFeasibility(scip, row);
}

/** returns the activity of a row for the given primal solution
 *
 *  @return the activitiy of a row for the given primal solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getRowSolActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow,               /**< LP row */
   jlong                 jsol                /**< primal CIP solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   return (jdouble) SCIPgetRowSolActivity(scip, row, sol);
}

/** returns the feasibility of a row for the given primal solution
 *
 *  @return the feasibility of a row for the given primal solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getRowSolFeasibility)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow,               /**< LP row */
   jlong                 jsol                /**< primal CIP solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   return (jdouble) SCIPgetRowSolFeasibility(scip, row, sol);
}

/** output row to file stream via the message handler system
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
void JNISCIP(printRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow,               /**< LP row */
   jlong                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPprintRow(scip, row, (FILE*)file) );
}

/** returns whether the NLP relaxation has been enabled
 *
 * If the NLP relaxation is enabled, then SCIP will construct the NLP relaxation when the solving process is about to begin.
 * To check whether an NLP is existing, use SCIPisNLPConstructed().
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 * @see SCIPenableNLP
 */
jboolean JNISCIP(isNLPEnabled)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisNLPEnabled(scip);
}

/** marks that there are constraints that are representable by nonlinear rows
 *
 * This method should be called by a constraint handler if it has constraints that have a representation as nonlinear rows.
 *
 * The function should be called before the branch-and-bound process is initialized, e.g., when presolve is exiting.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(enableNLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   SCIPenableNLP(scip);
}

/** returns, whether an NLP has been constructed
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jboolean JNISCIP(isNLPConstructed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisNLPConstructed(scip);
}

/** returns whether the NLP has a continuous variable in a nonlinear term
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jboolean JNISCIP(hasNLPContinuousNonlinearity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPhasNLPContinuousNonlinearity(scip);
}

/** gets array with variables of the NLP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jlongArray JNISCIP(getNLPVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nvars;

   jlongArray jvars;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nvars = SCIPgetNNLPVars(scip);
   jvars = (*env)->NewLongArray(env, nvars);

   if( jvars == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_VAR** vars;

      vars = SCIPgetNLPVars(scip);
      (*env)->SetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);
   }

   return jvars;
}

/** gets current number of variables in NLP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jint JNISCIP(getNNLPVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNNLPVars(scip);
}

/** TODO: getNLPVarsNonlinearity - length of array?*/
/** computes for each variables the number of NLP rows in which the variable appears in a nonlinear var
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
#if 0
void JNISCIP(getNLPVarsNonlinearity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jintArray             jnlcount            /**< an array of length at least SCIPnlpGetNVars() to store nonlinearity counts of variables */
   )
{
   SCIP* scip;
   int* nlcount;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &nlcount, ?length?) );

   (*env)->GetIntArrayRegion(env, jnlcount, 0, ?length?, (jlong*)nlcount);

   SCIPfreeBufferArray(scip, &nlcount);

   return (jint) SCIPgetNLPVarsNonlinearity(scip);
}
#endif

/** returns dual solution values associated with lower bounds of NLP variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jdoubleArray JNISCIP(getNLPVarsLbDualsol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nvars;

   jdoubleArray jsols;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nvars = SCIPgetNNLPVars(scip);
   jsols = (*env)->NewDoubleArray(env, nvars);

   if( jsols == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_Real* sols;

      sols = SCIPgetNLPVarsLbDualsol(scip);
      (*env)->SetDoubleArrayRegion(env, jsols, 0, nvars, (jdouble*)sols);
   }

   return jsols;
}

/** returns dual solution values associated with upper bounds of NLP variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jdoubleArray JNISCIP(getNLPVarsUbDualsol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nvars;

   jdoubleArray jsols;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nvars = SCIPgetNNLPVars(scip);
   jsols = (*env)->NewDoubleArray(env, nvars);

   if (jsols == NULL)
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_Real* sols;

      sols = SCIPgetNLPVarsUbDualsol(scip);
      (*env)->SetDoubleArrayRegion(env, jsols, 0, nvars, (jdouble*)sols);
   }

   return jsols;
}

/** gets array with nonlinear rows of the NLP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jlongArray JNISCIP(getNLPNlRows)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nvars;

   jlongArray jnlrows;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nvars = SCIPgetNNLPNlRows(scip);
   jnlrows = (*env)->NewLongArray(env, nvars);

   if( jnlrows == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_NLROW** nlrows;

      nlrows = SCIPgetNLPNlRows(scip);
      (*env)->SetLongArrayRegion(env, jnlrows, 0, nvars, (jlong*)nlrows);
   }

   return jnlrows;
}

/** gets current number of nonlinear rows in NLP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jint JNISCIP(getNNLPNlRows)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNNLPNlRows(scip);
}

/** adds a nonlinear row to the NLP. This row is captured by the NLP.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(addNlRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow              /**< nonlinear row to add to NLP */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPaddNlRow(scip, nlrow) );
}

/** makes sure that the NLP of the current node is flushed
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(flushNLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);


   JNISCIP_CALL( SCIPflushNLP(scip) );
}

/** TODO: setNLPInitialGuess - length of array? */
/** sets or clears initial primal guess for NLP solution (start point for NLP solver)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
#if 0
void JNISCIP(setNLPInitialGuess)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdoubleArray          jinitialguess       /**< values of initial guess (corresponding to variables from SCIPgetNLPVarsData), or NULL to use no start point */
   )
{
   SCIP* scip;
   SCIP_Real* initialguess;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &initialguess, ?length?) );

   (*env)->GetDoubleArrayRegion(env, jinitialguessvars, 0, ?length?, (jlong*)initialguess);

   SCIPfreeBufferArray(scip, &initialguess);

   JNISCIP_CALL( SCIPsetNLPInitialGuess(scip, initialguess) );
}
#endif

/** sets initial primal guess for NLP solution (start point for NLP solver)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(setNLPInitialGuessSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 sol                 /**< solution which values should be taken as initial guess, or NULL for LP solution */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetNLPInitialGuessSol(scip, (SCIP_SOL*)(size_t)sol) );
}

/** solves the current NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(solveNLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);


   JNISCIP_CALL( SCIPsolveNLP(scip) );
}

/** gets solution status of current NLP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jint JNISCIP(getNLPSolstat)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNLPSolstat(scip);
}

/** gets termination status of last NLP solve
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jint JNISCIP(getNLPTermstat)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNLPTermstat(scip);
}

/** gives statistics (number of iterations, solving time, ...) of last NLP solve
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getNLPObjval)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);


   return (jint) SCIPgetNLPObjval(scip);
}

/** indicates whether a feasible solution for the current NLP is available
 * thus, returns whether the solution status <= feasible
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jboolean JNISCIP(hasNLPSolution)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);


   return (jboolean) SCIPhasNLPSolution(scip);
}

/** gets integer parameter of NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jint JNISCIP(getNLPIntPar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  type                /**< parameter number */
   )
{
   SCIP* scip;
   int ival;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNLPIntPar(scip, (SCIP_NLPPARAM)type, &ival);
}

/** sets integer parameter of NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(setNLPIntPar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  type,               /**< parameter number */
   jint                  ival                /**< parameter value */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetNLPIntPar(scip, (SCIP_NLPPARAM)type, (int)ival) );
}

/** gets floating point parameter of NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getNLPRealPar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  type                /**< parameter number */
   )
{
   SCIP* scip;
   SCIP_Real dval;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPgetNLPRealPar(scip, (SCIP_NLPPARAM)type, &dval) );
   return (jdouble) dval;
}

/** sets floating point parameter of NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(setNLPRealPar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  type,               /**< parameter number */
   jdouble               dval                /**< parameter value */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsetNLPRealPar(scip, (SCIP_NLPPARAM)type, (SCIP_Real)dval) );
}

/** gets string parameter of NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jstring JNISCIP(getNLPStringPar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  type                /**< parameter number */
   )
{
   SCIP* scip;
   const char* sval;
   jstring jsval;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPgetNLPStringPar(scip, (SCIP_NLPPARAM)type, &sval) );

   jsval = (*env)->NewStringUTF(env, sval);

   return jsval;
}

/** sets string parameter of NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(setNLPStringPar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  type,               /**< parameter number */
   jstring               jsval               /**< parameter value */
   )
{
   SCIP* scip;
   const char* sval;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   sval = (*env)->GetStringUTFChars(env, jsval, &iscopy);
   if( sval == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPsetNLPStringPar(scip, (SCIP_NLPPARAM)type, sval) );

   (*env)->ReleaseStringUTFChars(env, jsval, sval);
}

/** writes current NLP to a file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(writeNLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jfilename           /**< file name */
   )
{
   SCIP* scip;
   const char* filename;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   if( jfilename != NULL )
      filename = (*env)->GetStringUTFChars(env, jfilename, &iscopy);
   else
      filename = NULL;

   JNISCIP_CALL( SCIPwriteNLP(scip, filename) );

   (*env)->ReleaseStringUTFChars(env, jfilename, filename);
}

/** initiates NLP diving
 * making methods SCIPchgVarObjDiveNLP(), SCIPchgVarBoundsDiveNLP(), SCIPchgVarsBoundsDiveNLP(), and SCIPsolveDiveNLP() available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(startDiveNLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPstartDiveNLP(scip) );

}

/** ends NLP diving
 * resets changes made by SCIPchgVarObjDiveNLP(), SCIPchgVarBoundsDiveNLP(), and SCIPchgVarsBoundsDiveNLP()
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(endDiveNLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPendDiveNLP(scip) );
}

/** changes linear objective coefficient of a variable in diving NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(chgVarObjDiveNLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable which coefficient to change */
   jdouble               coef                /**< new value for coefficient */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarObjDiveNLP(scip, var, (SCIP_Real)coef) );
}

/** changes bounds of a variable in diving NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(chgVarBoundsDiveNLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable which bounds to change */
   jdouble               lb,                 /**< new lower bound */
   jdouble               ub                  /**< new upper bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarBoundsDiveNLP(scip, var, (SCIP_Real)lb, (SCIP_Real)ub) );
}

/** changes bounds of a set of variables in diving NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(chgVarsBoundsDiveNLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  nvars,              /**< number of variables which bounds to changes */
   jlongArray            jvars,              /**< variables which bounds to change */
   jdoubleArray          jlbs,               /**< new lower bounds */
   jdoubleArray          jubs                /**< new upper bounds */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_Real* lbs;
   SCIP_Real* ubs;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &lbs, (int)nvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &ubs, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)vars);
   (*env)->GetDoubleArrayRegion(env, jlbs, 0, (int)nvars, (jdouble*)lbs);
   (*env)->GetDoubleArrayRegion(env, jubs, 0, (int)nvars, (jdouble*)ubs);

   JNISCIP_CALL( SCIPchgVarsBoundsDiveNLP(scip, (int)nvars, vars, lbs, ubs) );

   SCIPfreeBufferArray(scip, &ubs);
   SCIPfreeBufferArray(scip, &lbs);
   SCIPfreeBufferArray(scip, &vars);
}

/** solves diving NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(solveDiveNLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPsolveDiveNLP(scip) );
}

/** creates and captures an NLP nonlinear row without any coefficients
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jlong JNISCIP(createEmptyNlRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jname,              /**< name of row */
   jdouble               lhs,                /**< left hand side of row */
   jdouble               rhs                 /**< right hand side of row */
   )
{
   SCIP* scip;
   SCIP_NLROW* row;
   const char* name;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   name = (*env)->GetStringUTFChars(env, jname, &iscopy);
   if( name == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPcreateEmptyNlRow(scip, &row, name, (SCIP_Real)lhs, (SCIP_Real)rhs) );

   (*env)->ReleaseStringUTFChars(env, jname, name);

   return (jlong) row;
}

/** creates and captures an NLP row from a linear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jlong JNISCIP(createNlRowFromRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< the linear row to copy */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPcreateNlRowFromRow(scip, &nlrow, row) );

   return (jlong) nlrow;
}

/** increases usage counter of NLP nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(captureNlRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow              /**< nonlinear row to capture */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPcaptureNlRow(scip, nlrow) );
}

/** decreases usage counter of NLP nonlinear row, and frees memory if necessary
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
void JNISCIP(releaseNlRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_NLROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_NLROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPreleaseNlRow(scip, &row) );
}

/** changes left hand side of NLP nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(chgNlRowLhs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow,             /**< NLP nonlinear row */
   jdouble               lhs                 /**< new left hand side */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPchgNlRowLhs(scip, nlrow, (SCIP_Real)lhs) );
}

/** changes right hand side of NLP nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(chgNlRowRhs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow,             /**< NLP nonlinear row */
   jdouble               rhs                 /**< new right hand side */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPchgNlRowRhs(scip, nlrow, (SCIP_Real)rhs) );
}

/** changes constant of NLP nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(chgNlRowConstant)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow,             /**< NLP nonlinear row */
   jdouble               constant            /**< new value for constant */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPchgNlRowConstant(scip, nlrow, (SCIP_Real)constant) );
}

/** adds variable with a linear coefficient to the nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(addLinearCoefToNlRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow,             /**< NLP nonlinear row */
   jlong                 jvar,               /**< problem variable */
   jdouble               val                 /**< value of coefficient in linear part of row */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddLinearCoefToNlRow(scip, nlrow, var, (SCIP_Real)val) );
}

/** adds variables with linear coefficients to the row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(addLinearCoefsToNlRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow,             /**< NLP row */
   jint                  nvars,              /**< number of variables to add to the row */
   jlongArray            jvars,              /**< problem variables to add */
   jdoubleArray          jvals               /**< values of coefficients in linear part of row */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)vars);
   (*env)->GetDoubleArrayRegion(env, jvals, 0, (int)nvars, (jdouble*)vals);

   JNISCIP_CALL( SCIPaddLinearCoefsToNlRow(scip, nlrow, (int)nvars, vars, vals) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
}

/** changes linear coefficient of a variables in a row
 * setting the coefficient to 0.0 means that it is removed from the row
 * the variable does not need to exists before
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(chgNlRowLinearCoef)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow,             /**< NLP nonlinear row */
   jlong                 jvar,               /**< problem variable */
   jdouble               val                 /**< new value of coefficient */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgNlRowLinearCoef(scip, nlrow, var, (SCIP_Real)val) );
}

/** adds quadratic variable to the nonlinear row
 * after adding a quadratic variable, it can be used to add quadratic elements
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(addQuadVarToNlRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow,             /**< NLP nonlinear row */
   jlong                 jvar                /**< problem variable */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddQuadVarToNlRow(scip, nlrow, var) );
}

/** sets a parameter of expression tree in the nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(setNlRowExprtreeParam)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow,             /**< NLP nonlinear row */
   jint                  paramidx,           /**< index of parameter in expression tree */
   jdouble               paramval            /**< new value of parameter in expression tree */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPsetNlRowExprtreeParam(scip, nlrow, (SCIP_Real)paramidx, (SCIP_Real)paramval) );
}

/** TODO: setNlRowExprtreeParams - length? */
/** sets parameters of expression tree in the nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
#if 0
void JNISCIP(setNlRowExprtreeParams)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow,             /**< NLP row */
   jdoubleArray          jparamvals          /**< new values of parameter in expression tree */+
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_Real* paramvals;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &paramvals, ?length?) );

   (*env)->GetDoubleArrayRegion(env, jparamvals, 0, ?length?, (jdouble*)paramvals);

   JNISCIP_CALL( SCIPsetNlRowExprtreeParams(scip, nlrow, paramvals) );

   SCIPfreeBufferArray(scip, &paramvals);
}
#endif

/** recalculates the activity of a nonlinear row in the last NLP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(recalcNlRowNLPActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow              /**< NLP nonlinear row */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPrecalcNlRowNLPActivity(scip, nlrow) );
}

/** returns the activity of a nonlinear row in the last NLP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getNlRowNLPActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow              /**< NLP nonlinear row */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_Real activity;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPgetNlRowNLPActivity(scip, nlrow, &activity) );

   return (jdouble) activity;
}

/** gives the feasibility of a nonlinear row in the last NLP solution: negative value means infeasibility
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getNlRowNLPFeasibility)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow              /**< NLP nonlinear row */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_Real feas;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPgetNlRowNLPFeasibility(scip, nlrow, &feas) );

   return (jdouble) feas;
}

/** recalculates the activity of a nonlinear row for the current pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(recalcNlRowPseudoActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow              /**< NLP nonlinear row */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPrecalcNlRowPseudoActivity(scip, nlrow) );
}

/** gives the activity of a nonlinear row for the current pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getNlRowPseudoActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow              /**< NLP nonlinear row */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_Real activity;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPgetNlRowPseudoActivity(scip, nlrow, &activity) );

   return (jdouble) activity;
}

/** gives the feasibility of a nonlinear row for the current pseudo solution: negative value means infeasibility
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getNlRowPseudoFeasibility)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow              /**< NLP nonlinear row */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_Real feas;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPgetNlRowPseudoFeasibility(scip, nlrow, &feas) );

   return (jdouble) feas;
}

/** recalculates the activity of a nonlinear row in the last NLP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(recalcNlRowActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow              /**< NLP nonlinear row */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPrecalcNlRowActivity(scip, nlrow) );
}

/** gives the activity of a nonlinear row in the last NLP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getNlRowActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow              /**< NLP nonlinear row */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_Real activity;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPgetNlRowActivity(scip, nlrow, &activity) );

   return (jdouble) activity;
}

/** gives the feasibility of a nonlinear row in the last NLP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getNlRowFeasibility)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow              /**< NLP nonlinear row */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_Real feas;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPgetNlRowFeasibility(scip, nlrow, &feas) );

   return (jdouble) feas;
}

/** gives the activity of a nonlinear row for the given primal solution or NLP solution or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getNlRowSolActivity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow,             /**< NLP nonlinear row */
   jlong                 jsol                /**< primal CIP solution, or NULL for NLP solution of pseudo solution */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_SOL* sol;
   SCIP_Real activity;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrow, sol, &activity) );

   return (jdouble) activity;
}

/** gives the feasibility of a nonlinear row for the given primal solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getNlRowSolFeasibility)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow,             /**< NLP nonlinear row */
   jlong                 jsol                /**< primal CIP solution, or NULL for NLP solution of pseudo solution */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;
   SCIP_SOL* sol;
   SCIP_Real feas;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPgetNlRowSolFeasibility(scip, nlrow, sol, &feas) );

   return (jdouble) feas;
}

/** output nonlinear row to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(printNlRow)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnlrow,             /**< NLP nonlinear row */
   jlong                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;
   SCIP_NLROW* nlrow;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nlrow = (SCIP_NLROW*) (size_t) jnlrow;
   assert(nlrow != NULL);

   JNISCIP_CALL( SCIPprintNlRow(scip, nlrow, (FILE*)file) );
}

/** returns efficacy of the cut with respect to the given primal solution or the current LP solution:
 *  e = -feasibility/norm
 *
 *  @return the efficacy of the cut with respect to the given primal solution or the current LP solution:
 *          e = -feasibility/norm
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
jdouble JNISCIP(getCutEfficacy)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 sol,                /**< primal CIP solution, or NULL for current LP solution */
   jlong                 jcut                /**< separated cut */
   )
{
   SCIP* scip;
   SCIP_ROW* cut;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cut = (SCIP_ROW*) (size_t) jcut;
   assert(cut != NULL);

   return (jdouble) SCIPgetCutEfficacy(scip, (SCIP_SOL*)(size_t)sol, cut);
}

/** returns whether the cut's efficacy with respect to the given primal solution or the current LP solution is greater
 *  than the minimal cut efficacy
 *
 *  @return TRUE if the cut's efficacy with respect to the given primal solution or the current LP solution is greater
 *          than the minimal cut efficacy, otherwise FALSE
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
jboolean JNISCIP(isCutEfficacious)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 sol,                /**< primal CIP solution, or NULL for current LP solution */
   jlong                 jcut                /**< separated cut */
   )
{
   SCIP* scip;
   SCIP_ROW* cut;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cut = (SCIP_ROW*) (size_t) jcut;
   assert(cut != NULL);

   return (jboolean) SCIPisCutEfficacious(scip, (SCIP_SOL*)(size_t)sol, cut);
}

/** checks, if the given cut's efficacy is larger than the minimal cut efficacy
 *
 *  @return TRUE if the given cut's efficacy is larger than the minimal cut efficacy, otherwise FALSE
 */
jboolean JNISCIP(isEfficacious)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               efficacy            /**< efficacy of the cut */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisEfficacious(scip, (SCIP_Real)efficacy);
}

/** calculates the efficacy norm of the given vector, which depends on the "separating/efficacynorm" parameter
 *
 *  @return the efficacy norm of the given vector, which depends on the "separating/efficacynorm" parameter
 */
jdouble JNISCIP(getVectorEfficacyNorm)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdoubleArray          jvals,              /**< array of values */
   jint                  nvals               /**< number of values */
   )
{
   SCIP* scip;
   SCIP_Real* vals;
   SCIP_Real norm;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)nvals) );

   (*env)->GetDoubleArrayRegion(env, jvals, 0, (int)nvals, (jdouble*)vals);

   norm = SCIPgetVectorEfficacyNorm(scip, vals, (int)nvals);

   SCIPfreeBufferArray(scip, &vals);

   return (jdouble) norm;
}

/** adds cut to separation storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
jboolean JNISCIP(addCut)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 sol,                /**< primal solution that was separated, or NULL for LP solution */
   jlong                 jcut,               /**< separated cut */
   jboolean              forcecut            /**< should the cut be forced to enter the LP? */
   )
{
   SCIP* scip;
   SCIP_ROW* cut;
   SCIP_Bool infeasible;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cut = (SCIP_ROW*) (size_t) jcut;
   assert(cut != NULL);

   JNISCIP_CALL( SCIPaddCut(scip, (SCIP_SOL*)(size_t)sol, cut, (SCIP_Bool)forcecut, &infeasible) );

   return (jboolean) infeasible;
}

/** if not already existing, adds row to global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(addPoolCut)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< cutting plane to add */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPaddPoolCut(scip, row) );
}

/** removes the row from the global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
void JNISCIP(delPoolCut)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< cutting plane to add */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPdelPoolCut(scip, row) );
}

/** gets current cuts in the global cut pool
 *
 *  @return the current cuts in the global cut pool
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
jlongArray JNISCIP(getPoolCuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int ncuts;

   jlongArray jcuts;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   ncuts = SCIPgetNPoolCuts(scip);
   jcuts = (*env)->NewLongArray(env, ncuts);

   if (jcuts == NULL)
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_CUT** cuts;

      cuts = SCIPgetPoolCuts(scip);
      (*env)->SetLongArrayRegion(env, jcuts, 0, ncuts, (jlong*)cuts);
   }

   return jcuts;
}

/** gets current number of rows in the global cut pool
 *
 *  @return the current number of rows in the global cut pool
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
jint JNISCIP(getNPoolCuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPoolCuts(scip);
}

/** gets the global cut pool used by SCIP
 *
 *  @return the global cut pool used by SCIP
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
jlong JNISCIP(getGlobalCutpool)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) (size_t) SCIPgetGlobalCutpool(scip);
}

/** creates a cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
jlong JNISCIP(createCutpool)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  agelimit            /**< maximum age a cut can reach before it is deleted from the pool */
   )
{
   SCIP* scip;
   SCIP_CUTPOOL* cutpool;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPcreateCutpool(scip, &cutpool, (int)agelimit) );

   return (jlong) cutpool;
}

/** frees a cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
void JNISCIP(freeCutpool)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcutpool            /**< cut pool */
   )
{
   SCIP* scip;
   SCIP_CUTPOOL* cutpool;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cutpool = (SCIP_CUTPOOL*) (size_t) jcutpool;
   assert(cutpool != NULL);

   JNISCIP_CALL( SCIPfreeCutpool(scip, &cutpool) );
}

/** if not already existing, adds row to a cut pool and captures it
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(addRowCutpool)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcutpool,           /**< cut pool */
   jlong                 jrow                /**< cutting plane to add */
   )
{
   SCIP* scip;
   SCIP_CUTPOOL* cutpool;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cutpool = (SCIP_CUTPOOL*) (size_t) jcutpool;
   assert(cutpool != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPaddRowCutpool(scip, cutpool, row) );
}

/** adds row to a cut pool and captures it; doesn't check for multiple cuts
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(addNewRowCutpool)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcutpool,           /**< cut pool */
   jlong                 jrow                /**< cutting plane to add */
   )
{
   SCIP* scip;
   SCIP_CUTPOOL* cutpool;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cutpool = (SCIP_CUTPOOL*) (size_t) jcutpool;
   assert(cutpool != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPaddNewRowCutpool(scip, cutpool, row) );
}

/** removes the LP row from a cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
void JNISCIP(delRowCutpool)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcutpool,           /**< cut pool */
   jlong                 jrow                /**< cutting plane to add */
   )
{
   SCIP* scip;
   SCIP_CUTPOOL* cutpool;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cutpool = (SCIP_CUTPOOL*) (size_t) jcutpool;
   assert(cutpool != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPdelRowCutpool(scip, cutpool, row) );
}

/** separates cuts from a cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(separateCutpool)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jcutpool            /**< cut pool */
   )
{
   SCIP* scip;
   SCIP_CUTPOOL* cutpool;
   SCIP_RESULT result;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   cutpool = (SCIP_CUTPOOL*) (size_t) jcutpool;
   assert(cutpool != NULL);

   JNISCIP_CALL( SCIPseparateCutpool(scip, cutpool, &result) );

   return (jint) result;
}

/** if not already existing, adds row to delayed global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(addDelayedPoolCut)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< cutting plane to add */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPaddDelayedPoolCut(scip, row) );
}

/** removes the row from the delayed global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(delDelayedPoolCut)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< cutting plane to add */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPdelDelayedPoolCut(scip, row) );
}

/** gets current cuts in the delayed global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
jlongArray JNISCIP(getDelayedPoolCuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int ncuts;

   jlongArray jcuts;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   ncuts = SCIPgetNDelayedPoolCuts(scip);
   jcuts = (*env)->NewLongArray(env, ncuts);

   if( jcuts == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_CUT** cuts;

      cuts = SCIPgetDelayedPoolCuts(scip);
      (*env)->SetLongArrayRegion(env, jcuts, 0, ncuts, (jlong*)cuts);
   }

   return jcuts;
}

/** gets current number of rows in the delayed global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(getNDelayedPoolCuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNDelayedPoolCuts(scip);
}

/** gets the delayed global cut pool used by SCIP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(getDelayedGlobalCutpool)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) (size_t) SCIPgetDelayedGlobalCutpool(scip);
}

/** gets the array of cuts currently stored in the separation storage
 *
 *  @return the array of cuts currently stored in the separation storage
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlongArray JNISCIP(getCuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int ncuts;

   jlongArray jcuts;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   ncuts = SCIPgetNCuts(scip);
   jcuts = (*env)->NewLongArray(env, ncuts);

   if( jcuts == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_ROW** cuts;

      cuts = SCIPgetCuts(scip);
      (*env)->SetLongArrayRegion(env, jcuts, 0, ncuts, (jlong*)cuts);
   }

   return jcuts;
}

/** get current number of cuts in the separation storage
 *
 *  @return the current number of cuts in the separation storage
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jint JNISCIP(getNCuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNCuts(scip);
}

/** clears the separation storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(clearCuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPclearCuts(scip) );
}

/** removes cuts that are inefficacious w.r.t. the current LP solution from separation storage without adding the cuts to the LP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(removeInefficaciousCuts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPremoveInefficaciousCuts(scip) );
}

/** initiates LP diving, making methods SCIPchgVarObjDive(), SCIPchgVarLbDive(), and SCIPchgVarUbDive() available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note diving is allowed even if the current LP is not flushed, not solved, or not solved to optimality; be aware
 *  that solving the (first) diving LP may take longer than expect and that the latter two cases could stem from
 *  numerical troubles during the last LP solve; because of this, most users will want to call this method only if
 *  SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL
 */
JNIEXPORT
void JNISCIP(startDive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPstartDive(scip) );
}

/** quits LP diving and resets bounds and objective values of columns to the current node's values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
void JNISCIP(endDive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPendDive(scip) );
}

/** changes cutoffbound in current dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
void JNISCIP(chgCutoffboundDive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               newcutoffbound      /**< new cutoffbound */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPchgCutoffboundDive(scip, (SCIP_Real)newcutoffbound) );
}

/** changes variable's objective value in current dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
void JNISCIP(chgVarObjDive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to change the objective value for */
   jdouble               newobj              /**< new objective value */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarObjDive(scip, var, (SCIP_Real)newobj) );
}

/** changes variable's lower bound in current dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
void JNISCIP(chgVarLbDive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to change the bound for */
   jdouble               newvalue            /**< new value for bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarLbDive(scip, var, (SCIP_Real)newvalue) );
}

/** changes variable's upper bound in current dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
void JNISCIP(chgVarUbDive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to change the bound for */
   jdouble               newvalue            /**< new value for bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarUbDive(scip, var, (SCIP_Real)newvalue) );
}

/** adds a row to the LP in current dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
void JNISCIP(addRowDive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow                /**< row to be added */
   )
{
   SCIP* scip;
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   JNISCIP_CALL( SCIPaddRowDive(scip, row) );
}

/** gets variable's objective value in current dive
 *
 *  @return the variable's objective value in current dive.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jdouble JNISCIP(getVarObjDive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get the bound for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarObjDive(scip, var);
}

/** gets variable's lower bound in current dive
 *
 *  @return the variable's lower bound in current dive.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jdouble JNISCIP(getVarLbDive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get the bound for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarLbDive(scip, var);
}

/** gets variable's upper bound in current dive
 *
 *  @return the variable's upper bound in current dive.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jdouble JNISCIP(getVarUbDive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get the bound for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetVarUbDive(scip, var);
}

/** returns the number of the node in the current branch and bound run, where the last LP was solved in diving
 *  or probing mode
 *
 *  @return the number of the node in the current branch and bound run, where the last LP was solved in diving
 *  or probing mode.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jlong JNISCIP(getLastDivenode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetLastDivenode(scip);
}

/** returns whether we are in diving mode
 *
 *  @return whether we are in diving mode.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jboolean JNISCIP(inDive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPinDive(scip);
}

/** returns whether we are in probing mode; probing mode is activated via SCIPstartProbing() and stopped
 *  via SCIPendProbing()
 *
 *  @return TRUE, if SCIP is currently in probing mode, otherwise FALSE
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jboolean JNISCIP(inProbing)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPinProbing(scip);
}

/** initiates probing, making methods SCIPnewProbingNode(), SCIPbacktrackProbing(), SCIPchgVarLbProbing(),
 *  SCIPchgVarUbProbing(), SCIPfixVarProbing(), SCIPpropagateProbing(), and SCIPsolveProbingLP() available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note The collection of variable statistics is turned off during probing. If these statistics should be collected
 *        during probing use the method SCIPenableVarHistory() to turn the collection explicitly on.
 */
JNIEXPORT
void JNISCIP(startProbing)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPstartProbing(scip) );
}

/** creates a new probing sub node, whose changes can be undone by backtracking to a higher node in the probing path
 *  with a call to SCIPbacktrackProbing();
 *  using a sub node for each set of probing bound changes can improve conflict analysis
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(newProbingNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPnewProbingNode(scip) );
}

/** returns the current probing depth
 *
 *  @return the probing depth, i.e. the number of probing sub nodes existing in the probing path
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(getProbingDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetProbingDepth(scip);
}

/** undoes all changes to the problem applied in probing up to the given probing depth;
 *  the changes of the probing node of the given probing depth are the last ones that remain active;
 *  changes that were applied before calling SCIPnewProbingNode() cannot be undone
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(backtrackProbing)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  probingdepth        /**< probing depth of the node in the probing path that should be reactivated */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPbacktrackProbing(scip, (int)probingdepth) );
}

/** quits probing and resets bounds and constraints to the focus node's environment
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(endProbing)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPendProbing(scip) );
}

/** injects a change of variable's lower bound into current probing node; the same can also be achieved with a call to
 *  SCIPchgVarLb(), but in this case, the bound change would be treated like a deduction instead of a branching decision
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(chgVarLbProbing)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to change the bound for */
   jdouble               newbound            /**< new value for bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarLbProbing(scip, var, (SCIP_Real)newbound) );
}

/** injects a change of variable's upper bound into current probing node; the same can also be achieved with a call to
 *  SCIPchgVarUb(), but in this case, the bound change would be treated like a deduction instead of a branching decision
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(chgVarUbProbing)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to change the bound for */
   jdouble               newbound            /**< new value for bound */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPchgVarUbProbing(scip, var, (SCIP_Real)newbound) );
}

/** injects a change of variable's bounds into current probing node to fix the variable to the specified value;
 *  the same can also be achieved with a call to SCIPfixVar(), but in this case, the bound changes would be treated
 *  like deductions instead of branching decisions
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(fixVarProbing)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to change the bound for */
   jdouble               fixedval            /**< value to fix variable to */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPfixVarProbing(scip, var, (SCIP_Real)fixedval) );
}

/** applies domain propagation on the probing sub problem, that was changed after SCIPstartProbing() was called;
 *  only propagations of the binary variables fixed at the current probing node that are triggered by the implication
 *  graph and the clique table are applied;
 *  the propagated domains of the variables can be accessed with the usual bound accessing calls SCIPvarGetLbLocal()
 *  and SCIPvarGetUbLocal(); the propagation is only valid locally, i.e. the local bounds as well as the changed
 *  bounds due to SCIPchgVarLbProbing(), SCIPchgVarUbProbing(), and SCIPfixVarProbing() are used for propagation
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jboolean JNISCIP(propagateProbingImplications)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Bool cutoff;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPpropagateProbingImplications(scip, &cutoff) );

   return (jboolean) cutoff;
}

/** gets branching candidates for LP solution branching (fractional variables) along with solution values,
 *  fractionalities, and number of branching candidates;
 *  branching rules should always select the branching candidate among the first npriolpcands of the candidate
 *  list
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNLPBranchCands)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNLPBranchCands(scip);
}

/** gets number of branching candidates with maximal priority for LP solution branching
 *
 *  @return the number of branching candidates with maximal priority for LP solution branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNPrioLPBranchCands)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPrioLPBranchCands(scip);
}

/** gets number of external branching candidates
 *
 *  @return the number of external branching candidates.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNExternBranchCands)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNExternBranchCands(scip);
}

/** gets number of external branching candidates with maximal branch priority
 *
 *  @return the number of external branching candidates with maximal branch priority.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNPrioExternBranchCands)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPrioExternBranchCands(scip);
}

/** gets number of binary external branching candidates with maximal branch priority
 *
 *  @return the number of binary external branching candidates with maximal branch priority.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNPrioExternBranchBins)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPrioExternBranchBins(scip);
}

/** gets number of integer external branching candidates with maximal branch priority
 *
 *  @return the number of integer external branching candidates with maximal branch priority.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNPrioExternBranchInts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPrioExternBranchInts(scip);
}

/** gets number of implicit integer external branching candidates with maximal branch priority
 *
 *  @return the number of implicit integer external branching candidates with maximal branch priority.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNPrioExternBranchImpls)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPrioExternBranchImpls(scip);
}

/** gets number of continuous external branching candidates with maximal branch priority
 *
 *  @return the number of continuous external branching candidates with maximal branch priority.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNPrioExternBranchConts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPrioExternBranchConts(scip);
}

/** insert variable, its score and its solution value into the external branching candidate storage
 * the relative difference of the current lower and upper bounds of a continuous variable must be at least epsilon
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
void JNISCIP(addExternBranchCand)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to insert */
   jdouble               score,              /**< score of external candidate, e.g. infeasibility */
   jdouble               solval              /**< value of the variable in the current solution */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPaddExternBranchCand(scip, var, (SCIP_Real)score, (SCIP_Real)solval) );
}

/** removes all external candidates from the storage for external branching
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
void JNISCIP(clearExternBranchCands)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   SCIPclearExternBranchCands(scip);
}

/** checks whether the given variable is contained in the candidate storage for external branching
 *
 *  @return whether the given variable is contained in the candidate storage for external branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jboolean JNISCIP(containsExternBranchCand)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to look for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jboolean) SCIPcontainsExternBranchCand(scip, var);
}

/** gets number of branching candidates for pseudo solution branching (non-fixed variables)
 *
 *  @return the number branching candidates for pseudo solution branching (non-fixed variables).
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNPseudoBranchCands)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPseudoBranchCands(scip);
}

/** gets number of branching candidates with maximal branch priority for pseudo solution branching
 *
 *  @return the number of branching candidates with maximal branch priority for pseudo solution branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNPrioPseudoBranchCands)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPrioPseudoBranchCands(scip);
}

/** gets number of binary branching candidates with maximal branch priority for pseudo solution branching
 *
 *  @return the number of binary branching candidates with maximal branch priority for pseudo solution branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNPrioPseudoBranchBins)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPrioPseudoBranchBins(scip);
}

/** gets number of integer branching candidates with maximal branch priority for pseudo solution branching
 *
 *  @return the number of integer branching candidates with maximal branch priority for pseudo solution branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNPrioPseudoBranchInts)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPrioPseudoBranchInts(scip);
}

/** gets number of implicit integer branching candidates with maximal branch priority for pseudo solution branching
 *
 *  @return the number of implicit integer branching candidates with maximal branch priority for pseudo solution branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(getNPrioPseudoBranchImpls)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPrioPseudoBranchImpls(scip);
}

/** calculates the branching score out of the gain predictions for a binary branching
 *
 *  @return the branching score out of the gain predictions for a binary branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jdouble JNISCIP(getBranchScore)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 var,                /**< variable, of which the branching factor should be applied, or NULL */
   jdouble               downgain,           /**< prediction of objective gain for rounding downwards */
   jdouble               upgain              /**< prediction of objective gain for rounding upwards */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetBranchScore(scip, (SCIP_VAR*)(size_t)var, (SCIP_Real)downgain, (SCIP_Real)upgain);
}

/** calculates the branching score out of the gain predictions for a branching with arbitrary many children
 *
 *  @return the branching score out of the gain predictions for a branching with arbitrary many children.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jdouble JNISCIP(getBranchScoreMultiple)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 var,                /**< variable, of which the branching factor should be applied, or NULL */
   jint                  nchildren,          /**< number of children that the branching will create */
   jdoubleArray          jgains              /**< prediction of objective gain for each child */
   )
{
   SCIP* scip;
   SCIP_Real* gains;
   SCIP_Real mult;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &gains, (int)nchildren) );

   (*env)->GetDoubleArrayRegion(env, jgains, 0, (int)nchildren, (jdouble*)gains);

   mult = SCIPgetBranchScoreMultiple(scip, (SCIP_VAR*)(size_t)var, (int)nchildren, gains);

   SCIPfreeBufferArray(scip, &gains);

   return (jdouble) mult;
}

/** computes a branching point for a continuous or discrete variable
 * @see SCIPbranchGetBranchingPoint
 *
 *  @return the branching point for a continuous or discrete variable.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jdouble JNISCIP(getBranchingPoint)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable, of which the branching point should be computed */
   jdouble               suggestion          /**< suggestion for branching point, or SCIP_INVALID if no suggestion */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetBranchingPoint(scip, var, (SCIP_Real)suggestion);
}

/** calculates the node selection priority for moving the given variable's LP value to the given target value;
 *  this node selection priority can be given to the SCIPcreateChild() call
 *
 *  @return the node selection priority for moving the given variable's LP value to the given target value.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jdouble JNISCIP(calcNodeselPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable on which the branching is applied */
   jint                  branchdir,          /**< type of branching that was performed: upwards, downwards, or fixed;
					      *   fixed should only be used, when both bounds changed
					      */
   jdouble               targetvalue         /**< new value of the variable in the child node */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPcalcNodeselPriority(scip, var, (SCIP_BRANCHDIR)branchdir, (SCIP_Real)targetvalue);
}

/** calculates an estimate for the objective of the best feasible solution contained in the subtree after applying the given
 *  branching; this estimate can be given to the SCIPcreateChild() call
 *
 *  @return the estimate for the objective of the best feasible solution contained in the subtree after applying the given
 *  branching.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jdouble JNISCIP(calcChildEstimate)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable on which the branching is applied */
   jdouble               targetvalue         /**< new value of the variable in the child node */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPcalcChildEstimate(scip, var, (SCIP_Real)targetvalue);
}

/** creates a child node of the focus node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jlong JNISCIP(createChild)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               nodeselprio,        /**< node selection priority of new node */
   jdouble               estimate            /**< estimate for(transformed) objective value of best feasible solution in subtree */
   )
{
   SCIP* scip;
   SCIP_NODE* node;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPcreateChild(scip, &node, (SCIP_Real)nodeselprio, (SCIP_Real)estimate) );

   return (jlong)(size_t) node;
}

/** n-ary branching on a variable x using a given value
 * Branches on variable x such that up to n/2 children are created on each side of the usual branching value.
 * The branching value is selected as in SCIPbranchVarVal().
 * The parameters minwidth and widthfactor determine the domain width of the branching variable in the child nodes.
 * If n is odd, one child with domain width 'width' and having the branching value in the middle is created.
 * Otherwise, two children with domain width 'width' and being left and right of the branching value are created.
 * Next further nodes to the left and right are created, where width is multiplied by widthfactor with increasing distance from the first nodes.
 * The initial width is calculated such that n/2 nodes are created to the left and to the right of the branching value.
 * If this value is below minwidth, the initial width is set to minwidth, which may result in creating less than n nodes.
 *
 * Giving a large value for widthfactor results in creating children with small domain when close to the branching value
 * and large domain when closer to the current variable bounds. That is, setting widthfactor to a very large value and n to 3
 * results in a ternary branching where the branching variable is mostly fixed in the middle child.
 * Setting widthfactor to 1.0 results in children where the branching variable always has the same domain width
 * (except for one child if the branching value is not in the middle).
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(branchVarValNary)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< variable to branch on */
   jdouble               val,                /**< value to branch on */
   jint                  n,                  /**< attempted number of children to be created, must be >= 2 */
   jdouble               minwidth,           /**< minimal domain width in children */
   jdouble               widthfactor         /**< multiplier for children domain width with increasing distance from val, must be >= 1.0 */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   int nchildren;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPbranchVarValNary(scip, var, (SCIP_Real)val, (int)n, (SCIP_Real)minwidth, (SCIP_Real)widthfactor, &nchildren) );

   return (jint) nchildren;
}

/** calls branching rules to branch on an LP solution; if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *  if the branch priority of an unfixed variable is larger than the maximal branch priority of the fractional
 *  variables, pseudo solution branching is applied on the unfixed variables with maximal branch priority
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(branchLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_RESULT result;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPbranchLP(scip, &result) );

   return (jint) result;
}

/** calls branching rules to branch on a external candidates; if no such candidates exist, the result is SCIP_DIDNOTRUN
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(branchExtern)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_RESULT result;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPbranchExtern(scip, &result) );

   return (jint) result;
}

/** calls branching rules to branch on a pseudo solution; if no unfixed variables exist, the result is SCIP_DIDNOTRUN
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jint JNISCIP(branchPseudo)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_RESULT result;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPbranchPseudo(scip, &result) );

   return (jint) result;
}

/** creates a primal solution, initialized to zero */
JNIEXPORT
jlong JNISCIP(createSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jheur               /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP* scip;
   SCIP_HEUR* heur;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   heur = (SCIP_HEUR*) (size_t) jheur;

   JNISCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   return (jlong)(size_t) sol;
}

/** creates a primal solution, initialized to the current LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(createLPSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jheur               /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP* scip;
   SCIP_HEUR* heur;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   heur = (SCIP_HEUR*) (size_t) jheur;

   JNISCIP_CALL( SCIPcreateLPSol(scip, &sol, heur) );

   return (jlong)(size_t) sol;
}

/** creates a primal solution, initialized to the current NLP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(createNLPSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jheur               /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP* scip;
   SCIP_HEUR* heur;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   heur = (SCIP_HEUR*) (size_t) jheur;

   JNISCIP_CALL( SCIPcreateNLPSol(scip, &sol, heur) );

   return (jlong)(size_t) sol;
}

/** creates a primal solution, initialized to the current relaxation solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(createRelaxSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jheur               /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP* scip;
   SCIP_HEUR* heur;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   heur = (SCIP_HEUR*) (size_t) jheur;

   JNISCIP_CALL( SCIPcreateRelaxSol(scip, &sol, heur) );

   return (jlong)(size_t) sol;
}

/** creates a primal solution, initialized to the current pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(createPseudoSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jheur               /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP* scip;
   SCIP_HEUR* heur;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   heur = (SCIP_HEUR*) (size_t) jheur;

   JNISCIP_CALL( SCIPcreatePseudoSol(scip, &sol, heur) );

   return (jlong)(size_t) sol;
}

/** creates a primal solution, initialized to the current LP or pseudo solution, depending on whether the LP was solved
 *  at the current node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(createCurrentSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jheur               /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP* scip;
   SCIP_HEUR* heur;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   heur = (SCIP_HEUR*) (size_t) jheur;

   JNISCIP_CALL( SCIPcreateCurrentSol(scip, &sol, heur) );

   return (jlong)(size_t) sol;
}

/** creates a primal solution, initialized to unknown values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(createUnknownSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jheur               /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP* scip;
   SCIP_HEUR* heur;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   heur = (SCIP_HEUR*) (size_t) jheur;

   JNISCIP_CALL( SCIPcreateUnknownSol(scip, &sol, heur) );

   return (jlong)(size_t) sol;
}

/** creates a primal solution living in the original problem space, initialized to zero;
 *  a solution in original space allows to set original variables to values that would be invalid in the
 *  transformed problem due to preprocessing fixings or aggregations
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(createOrigSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jheur               /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP* scip;
   SCIP_HEUR* heur;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   heur = (SCIP_HEUR*) (size_t) jheur;

   JNISCIP_CALL( SCIPcreateOrigSol(scip, &sol, heur) );

   return (jlong)(size_t) sol;
}

/** creates a copy of a primal solution; note that a copy of a linked solution is also linked and needs to be unlinked
 *  if it should stay unaffected from changes in the LP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(createSolCopy)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsourcesol          /**< primal CIP solution to copy */
   )
{
   SCIP* scip;
   SCIP_SOL* sourcesol;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sourcesol = (SCIP_SOL*) (size_t) jsourcesol;
   assert(sourcesol != NULL);

   JNISCIP_CALL( SCIPcreateSolCopy(scip, &sol, sourcesol) );

   return (jlong)(size_t) sol;
}

/** frees primal CIP solution */
JNIEXPORT
void JNISCIP(freeSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< the solution to free */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPfreeSol(scip, &sol) );
}

/** links a primal solution to the current LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(linkLPSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< the solution to free */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPlinkLPSol(scip, sol) );
}

/** links a primal solution to the current NLP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(linkNLPSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< the solution to free */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPlinkNLPSol(scip, sol) );
}

/** links a primal solution to the current relaxation solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(linkRelaxSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< the solution to free */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPlinkRelaxSol(scip, sol) );
}

/** links a primal solution to the current pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(linkPseudoSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< the solution to free */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPlinkPseudoSol(scip, sol) );
}

/** links a primal solution to the current LP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(linkCurrentSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< the solution to free */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPlinkCurrentSol(scip, sol) );
}

/** clears a primal solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(clearSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< the solution to free */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPclearSol(scip, sol) );
}

/** stores solution values of variables in solution's own array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(unlinkSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< the solution to free */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPunlinkSol(scip, sol) );
}

/** sets value of variable in primal CIP solution */
JNIEXPORT
void JNISCIP(setSolVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol,               /**< primal solution */
   jlong                 jvar,               /**< variable to add to solution */
   jdouble               jval                /**< solution value of variable */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPsetSolVal(scip, sol, var, (SCIP_Real) jval) );
}

/** sets values of multiple variables in primal CIP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(setSolVals)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol,               /**< primal solution */
   jint                  nvars,              /**< number of variables to set solution value for */
   jlongArray            jvars,              /**< array with variables to add to solution */
   jdoubleArray          jvals               /**< array with solution values of variables */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );
   JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)vars);
   (*env)->GetDoubleArrayRegion(env, jvals, 0, (int)nvars, (jdouble*)vals);

   JNISCIP_CALL( SCIPsetSolVals(scip, sol, (int)nvars, vars, vals) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
}

/** increases value of variable in primal CIP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(incSolVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol,               /**< primal solution */
   jlong                 jvar,               /**< variable to increase solution value for */
   jdouble               incval              /**< increment for solution value of variable */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   JNISCIP_CALL( SCIPincSolVal(scip, sol, var, (SCIP_Real)incval) );
}

/** returns value of variable in primal CIP solution, or in current LP/pseudo solution */
JNIEXPORT
jdouble JNISCIP(getSolVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol,               /**< primal solution, or NULL for current LP/pseudo solution */
   jlong                 jvar                /**< variable to get value for */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_VAR* var;
   SCIP_Real val;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   val = SCIPgetSolVal(scip, sol, var);

   return (jdouble)val;
}

/** gets values of multiple variables in primal CIP solution */
JNIEXPORT
jdoubleArray JNISCIP(getSolVals)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol,               /**< primal solution, or NULL for current LP/pseudo solution */
   jint                  nvars,              /**< number of variables to get solution value for */
   jlongArray            jvars               /**< array with variables to get value for */
   )
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_SOL* sol;

   jdoubleArray jvals;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;

   JNISCIP_CALL( SCIPallocBufferArray(scip, &vars, (int)nvars) );

   (*env)->GetLongArrayRegion(env, jvars, 0, (int)nvars, (jlong*)vars);
   jvals = (*env)->NewDoubleArray(env, (int)nvars);

   if (jvals == NULL)
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_Real* vals;

      JNISCIP_CALL( SCIPallocBufferArray(scip, &vals, (int)nvars) );

      JNISCIP_CALL( SCIPgetSolVals(scip, sol, (int)nvars, vars, vals) );
      (*env)->SetDoubleArrayRegion(env, jvals, 0, (int)nvars, (jdouble*)vals);

      SCIPfreeBufferArray(scip, &vals);
   }

   SCIPfreeBufferArray(scip, &vars);

   return jvals;
}

/** returns objective value of primal CIP solution w.r.t. original problem, or current LP/pseudo objective value */
JNIEXPORT
jdouble JNISCIP(getSolOrigObj)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< primal solution, or NULL for current LP/pseudo solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_Real objval;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;

   objval = SCIPgetSolOrigObj(scip, sol);

   return (jdouble)objval;
}

/** returns transformed objective value of primal CIP solution, or transformed current LP/pseudo objective value
 *
 *  @return transformed objective value of primal CIP solution, or transformed current LP/pseudo objective value
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jdouble JNISCIP(getSolTransObj)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< primal solution, or NULL for current LP/pseudo solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_Real objval;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI pointer into C pointer */
   sol = (SCIP_SOL*) (size_t) jsol;

   objval = SCIPgetSolTransObj(scip, sol);

   return (jdouble)objval;
}

/** maps original space objective value into transformed objective value
 *
 *  @return transformed objective value
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(transformObj)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               obj                 /**< original space objective value to transform */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPtransformObj(scip, (SCIP_Real)obj);
}

/** maps transformed objective value into original space
 *
 *  @return objective value into original space
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(retransformObj)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               obj                 /**< original space objective value to transform */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPretransformObj(scip, (SCIP_Real)obj);
}

/** gets clock time, when this solution was found
 *
 *  @return clock time, when this solution was found
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jdouble JNISCIP(getSolTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< primal solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   return (jdouble) SCIPgetSolTime(scip, sol);
}

/** gets branch and bound run number, where this solution was found
 *
 *  @return branch and bound run number, where this solution was found
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jint JNISCIP(getSolRunnum)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< primal solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   return (jint) SCIPgetSolRunnum(scip, sol);
}

/** gets node number of the specific branch and bound run, where this solution was found
 *
 *  @return node number of the specific branch and bound run, where this solution was found
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jlong JNISCIP(getSolNodenum)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< primal solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   return (jlong) SCIPgetSolNodenum(scip, sol);
}

/** gets heuristic, that found this solution (or NULL if it's from the tree)
 *
 *  @return heuristic, that found this solution (or NULL if it's from the tree)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jlong JNISCIP(getSolHeur)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< primal solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   return (jlong) (size_t) SCIPgetSolHeur(scip, sol);
}

/** returns whether two given solutions are exactly equal
 *
 *  @return returns whether two given solutions are exactly equal
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jboolean JNISCIP(areSolsEqual)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol1,              /**< first primal CIP solution */
   jlong                 jsol2               /**< second primal CIP solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol1;
   SCIP_SOL* sol2;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol1 = (SCIP_SOL*) (size_t) jsol1;
   assert(sol1 != NULL);

   sol2 = (SCIP_SOL*) (size_t) jsol2;
   assert(sol2 != NULL);

   return (jboolean) SCIPareSolsEqual(scip, sol1, sol2);
}

/** outputs non-zero variables of solution in original problem space to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
void JNISCIP(printSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol,               /**< primal solution, or NULL for current LP/pseudo solution */
   jlong                 file,               /**< output file (or NULL for standard output) */
   jboolean              printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;

   JNISCIP_CALL( SCIPprintSol(scip, sol, (FILE*)(size_t)file, (SCIP_Bool)printzeros) );
}

/** outputs non-zero variables of solution in transformed problem space to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
void JNISCIP(printTransSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol,               /**< primal solution, or NULL for current LP/pseudo solution */
   jlong                 file,               /**< output file (or NULL for standard output) */
   jboolean              printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;

   JNISCIP_CALL( SCIPprintTransSol(scip, sol, (FILE*)(size_t)file, (SCIP_Bool)printzeros) );
}

/** outputs non-zero variables of solution representing a ray in original problem space to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
void JNISCIP(printRay)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol,               /**< primal solution, or NULL for current LP/pseudo solution */
   jlong                 file,               /**< output file (or NULL for standard output) */
   jboolean              printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;

   JNISCIP_CALL( SCIPprintRay(scip, sol, (FILE*)(size_t)file, (SCIP_Bool)printzeros) );
}

/** gets number of feasible primal solutions stored in the solution storage in case the problem is transformed; in case
 *  if the problem stage is SCIP_STAGE_PROBLEM, it returns the number solution in the original solution candidate
 *  storage
 *
 *  @return number of feasible primal solutions stored in the solution storage in case the problem is transformed; in case
 *          if the problem stage is SCIP_STAGE_PROBLEM, it returns the number solution in the original solution candidate
 *          storage
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jint JNISCIP(getNSols)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNSols(scip);
}

/** gets array of feasible primal solutions stored in the solution storage in case the problem is transformed; in case
 *  if the problem stage is in SCIP_STAGE_PROBLEM, it returns the number array of solution candidate stored
 *
 *  @return array of feasible primal solutions
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jlongArray JNISCIP(getSols)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   int nsols;

   jlongArray jsols;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   nsols = SCIPgetNSols(scip);
   jsols = (*env)->NewLongArray(env, nsols);

   if( jsols == NULL )
   {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
   }
   else
   {
      SCIP_SOL** sols;

      sols = SCIPgetSols(scip);
      (*env)->SetLongArrayRegion(env, jsols, 0, nsols, (jlong*)sols);
   }

   return jsols;
}

/** gets best feasible primal solution found so far, or NULL if no solution has been found */
JNIEXPORT
jlong JNISCIP(getBestSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_SOL* bestsol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   bestsol = SCIPgetBestSol(scip);

   assert(SCIPgetNSols(scip) == 0 || bestsol != NULL);

   return (jlong)(size_t)bestsol;
}


/** outputs best feasible primal solution found so far to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREE
 */
JNIEXPORT
void JNISCIP(printBestSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 file,               /**< output file (or NULL for standard output) */
   jboolean              printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPprintBestSol(scip, (FILE*)(size_t)file, (SCIP_Bool)printzeros) );
}


/** outputs best feasible primal solution found so far in transformed variables to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREE
 */
JNIEXPORT
void JNISCIP(printBestTransSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 file,               /**< output file (or NULL for standard output) */
   jboolean              printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPprintBestTransSol(scip, (FILE*)(size_t)file, (SCIP_Bool)printzeros) );
}

/** try to round given solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jboolean JNISCIP(roundSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< primal solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_Bool success;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIProundSol(scip, sol, &success) );

   return (jboolean) success;
}

/** retransforms solution to original problem space
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(retransformSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< primal solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPretransformSol(scip, sol) );
}

/** reads a given solution file, problem has to be transformed in advance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(readSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jstring               jfilename           /**< name of the input file */
   )
{
   SCIP* scip;
   const char* filename;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   filename = (*env)->GetStringUTFChars(env, jfilename, &iscopy);
   if( filename == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPreadSol(scip, filename) );

   (*env)->ReleaseStringUTFChars(env, jfilename, filename);
}

/** adds feasible primal solution to solution storage by copying it
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jboolean JNISCIP(addSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol                /**< primal solution */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_Bool stored;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPaddSol(scip, sol, &stored) );

   return (jboolean) stored;
}

/** adds current LP/pseudo solution to solution storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jboolean JNISCIP(addCurrentSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 heur                /**< heuristic that found the solution */
   )
{
   SCIP* scip;
   SCIP_Bool stored;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPaddCurrentSol(scip, (SCIP_HEUR*)(size_t)heur, &stored) );

   return (jboolean) stored;
}

/** checks solution for feasibility; if possible, adds it to storage by copying
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jboolean JNISCIP(trySol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol,               /**< primal CIP solution */
   jboolean              printreason,        /**< should all reasons of violation be printed? */
   jboolean              checkbounds,        /**< should the bounds of the variables be checked? */
   jboolean              checkintegrality,   /**< has integrality to be checked? */
   jboolean              checklprows         /**< have current LP rows (both local and global) to be checked? */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_Bool stored;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPtrySol(scip, sol, (SCIP_Bool)printreason, (SCIP_Bool)checkbounds, (SCIP_Bool)checkintegrality, (SCIP_Bool)checklprows, &stored) );

   return (jboolean) stored;
}

/** checks current LP/pseudo solution for feasibility; if possible, adds it to storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jboolean JNISCIP(tryCurrentSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 heur,               /**< heuristic that found the solution */
   jboolean              printreason,        /**< should all reasons of violation be printed? */
   jboolean              checkintegrality,   /**< has integrality to be checked? */
   jboolean              checklprows         /**< have current LP rows (both local and global) to be checked? */
   )
{
   SCIP* scip;
   SCIP_Bool stored;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPtryCurrentSol(scip, (SCIP_HEUR*)(size_t)heur, (SCIP_Bool)printreason, (SCIP_Bool)checkintegrality, (SCIP_Bool)checklprows, &stored) );

   return (jboolean) stored;
}

/** checks solution for feasibility without adding it to the solution store
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jboolean JNISCIP(checkSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol,               /**< primal CIP solution */
   jboolean              printreason,        /**< should all reasons of violation be printed? */
   jboolean              checkbounds,        /**< should the bounds of the variables be checked? */
   jboolean              checkintegrality,   /**< has integrality to be checked? */
   jboolean              checklprows         /**< have current LP rows (both local and global) to be checked? */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_Bool stored;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPcheckSol(scip, sol, (SCIP_Bool)printreason, (SCIP_Bool)checkbounds, (SCIP_Bool)checkintegrality, (SCIP_Bool)checklprows, &stored) );

   return (jboolean) stored;
}

/** checks solution for feasibility in original problem without adding it to the solution store;
 *  this method is used to double check a solution in order to validate the presolving process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jboolean JNISCIP(checkSolOrig)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jsol,               /**< primal CIP solution */
   jboolean              printreason,        /**< should the reason for the violation be printed? */
   jboolean              completely          /**< should all violation be checked? */
   )
{
   SCIP* scip;
   SCIP_SOL* sol;
   SCIP_Bool feasible;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   sol = (SCIP_SOL*) (size_t) jsol;
   assert(sol != NULL);

   JNISCIP_CALL( SCIPcheckSolOrig(scip, sol, &feasible, (SCIP_Bool)printreason, (SCIP_Bool)completely) );

   return (jboolean) feasible;
}

/** return whether a primal ray is stored that proves unboundedness of the LP relaxation
 *
 *  @return return whether a primal ray is stored that proves unboundedness of the LP relaxation
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jboolean JNISCIP(hasPrimalRay)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPhasPrimalRay(scip) ;
}

/** gets value of given variable in primal ray causing unboundedness of the LP relaxation;
 *  should only be called if such a ray is stored (check with SCIPhasPrimalRay())
 *
 *  @return value of given variable in primal ray causing unboundedness of the LP relaxation
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getPrimalRayVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar                /**< variable to get value for */
   )
{
   SCIP* scip;
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   return (jdouble) SCIPgetPrimalRayVal(scip, var) ;
}

/** catches a global (not variable or row dependent) event
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jint JNISCIP(catchEvent)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  eventtype,          /**< event type mask to select events to catch */
   jlong                 jeventhdlr,         /**< event handler to process events with */
   jlong                 jeventdata          /**< event data to pass to the event handler when processing this event */
   )
{
   SCIP* scip;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTDATA* eventdata;
   int filterpos;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   eventhdlr = (SCIP_EVENTHDLR*) (size_t) jeventhdlr;
   assert(eventhdlr != NULL);

   eventdata = (SCIP_EVENTDATA*) (size_t) jeventdata;
   assert(eventdata != NULL);

   JNISCIP_CALL( SCIPcatchEvent(scip, (SCIP_EVENTTYPE)eventtype, eventhdlr, eventdata, &filterpos) );

   return (jint) filterpos;
}

/** drops a global event (stops to track event)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(dropEvent)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  eventtype,          /**< event type mask to select events to catch */
   jlong                 jeventhdlr,         /**< event handler to process events with */
   jlong                 jeventdata,         /**< event data to pass to the event handler when processing this event */
   jint                  filterpos           /**< position of event filter entry returned by SCIPcatchEvent(), or -1 */
   )
{
   SCIP* scip;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTDATA* eventdata;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   eventhdlr = (SCIP_EVENTHDLR*) (size_t) jeventhdlr;
   assert(eventhdlr != NULL);

   eventdata = (SCIP_EVENTDATA*) (size_t) jeventdata;
   assert(eventdata != NULL);

   JNISCIP_CALL( SCIPdropEvent(scip, (SCIP_EVENTTYPE)eventtype, eventhdlr, eventdata, (int)filterpos) );
}

/** catches an objective value or domain change event on the given transformed variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jint JNISCIP(catchVarEvent)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< transformed variable to catch event for */
   jint                  eventtype,          /**< event type mask to select events to catch */
   jlong                 jeventhdlr,         /**< event handler to process events with */
   jlong                 jeventdata          /**< event data to pass to the event handler when processing this event */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTDATA* eventdata;
   int filterpos;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   eventhdlr = (SCIP_EVENTHDLR*) (size_t) jeventhdlr;
   assert(eventhdlr != NULL);

   eventdata = (SCIP_EVENTDATA*) (size_t) jeventdata;
   assert(eventdata != NULL);

   JNISCIP_CALL( SCIPcatchVarEvent(scip, var, (SCIP_EVENTTYPE)eventtype, eventhdlr, eventdata, &filterpos) );

   return (jint) filterpos;
}

/** drops an objective value or domain change event (stops to track event) on the given transformed variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(dropVarEvent)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jvar,               /**< transformed variable to catch event for */
   jint                  eventtype,          /**< event type mask to select events to catch */
   jlong                 jeventhdlr,         /**< event handler to process events with */
   jlong                 jeventdata,         /**< event data to pass to the event handler when processing this event */
   jint                  filterpos           /**< position of event filter entry returned by SCIPcatchVarEvent(), or -1 */
   )
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTDATA* eventdata;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   eventhdlr = (SCIP_EVENTHDLR*) (size_t) jeventhdlr;
   assert(eventhdlr != NULL);

   eventdata = (SCIP_EVENTDATA*) (size_t) jeventdata;
   assert(eventdata != NULL);

   JNISCIP_CALL( SCIPdropVarEvent(scip, var, (SCIP_EVENTTYPE)eventtype, eventhdlr, eventdata, (int)filterpos) );
}

/** catches a row coefficient, constant, or side change event on the given row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jint JNISCIP(catchRowEvent)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow,               /**< linear row to catch event for */
   jint                  eventtype,          /**< event type mask to select events to catch */
   jlong                 jeventhdlr,         /**< event handler to process events with */
   jlong                 jeventdata          /**< event data to pass to the event handler when processing this event */
   )
{
   SCIP* scip;
   SCIP_ROW* row;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTDATA* eventdata;
   int filterpos;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   eventhdlr = (SCIP_EVENTHDLR*) (size_t) jeventhdlr;
   assert(eventhdlr != NULL);

   eventdata = (SCIP_EVENTDATA*) (size_t) jeventdata;
   assert(eventdata != NULL);

   JNISCIP_CALL( SCIPcatchRowEvent(scip, row, (SCIP_EVENTTYPE)eventtype, eventhdlr, eventdata, &filterpos) );

   return (jint) filterpos;
}

/** drops a row coefficient, constant, or side change event (stops to track event) on the given row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(dropRowEvent)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrow,               /**< linear row to drop event for */
   jint                  eventtype,          /**< event type mask to select events to catch */
   jlong                 jeventhdlr,         /**< event handler to process events with */
   jlong                 jeventdata,         /**< event data to pass to the event handler when processing this event */
   jint                  filterpos           /**< position of event filter entry returned by SCIPcatchVarEvent(), or -1 */
   )
{
   SCIP* scip;
   SCIP_ROW* row;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTDATA* eventdata;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   eventhdlr = (SCIP_EVENTHDLR*) (size_t) jeventhdlr;
   assert(eventhdlr != NULL);

   eventdata = (SCIP_EVENTDATA*) (size_t) jeventdata;
   assert(eventdata != NULL);

   JNISCIP_CALL( SCIPdropRowEvent(scip, row, (SCIP_EVENTTYPE)eventtype, eventhdlr, eventdata, (int)filterpos) );
}

/** gets current node in the tree
 *
 *  @return the current node of the search tree
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(getCurrentNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) (size_t) SCIPgetCurrentNode(scip);
}

/** gets the root node of the tree
 *
 *  @return the root node of the search tree
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(getRootNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) (size_t) SCIPgetRootNode(scip);
}

/** returns whether the current node is already solved and only propagated again
 *
 *  @return TRUE is returned if \SCIP performance repropagation, otherwise FALSE.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jboolean JNISCIP(inRepropagation)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPinRepropagation(scip);
}

/** gets number of children of focus node
 *
 *  @return number of children of the focus node
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(getNChildren)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNChildren(scip);
}

/** gets number of siblings of focus node
 *
 *  @return the number of siblings of focus node
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(getNSiblings)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNSiblings(scip);
}

/** gets number of leaves in the tree
 *
 *  @return the number of leaves in the tree
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(getNLeaves)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNLeaves(scip);
}

/** gets the best child of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @return the best child of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(getPrioChild)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) (size_t) SCIPgetPrioChild(scip);
}

/** gets the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @return the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(getPrioSibling)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) (size_t) SCIPgetPrioSibling(scip);
}

/** gets the best child of the focus node w.r.t. the node selection strategy
 *
 *  @return the best child of the focus node w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(getBestChild)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) (size_t) SCIPgetBestChild(scip);
}

/** gets the best sibling of the focus node w.r.t. the node selection strategy
 *
 *  @return the best sibling of the focus node w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(getBestSibling)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) (size_t) SCIPgetBestSibling(scip);
}

/** gets the best leaf from the node queue w.r.t. the node selection strategy
 *
 *  @return the best leaf from the node queue w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(getBestLeaf)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) (size_t) SCIPgetBestLeaf(scip);
}

/** gets the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy
 *
 *  @return the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(getBestNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) (size_t) SCIPgetBestNode(scip);
}

/** gets the node with smallest lower bound from the tree (child, sibling, or leaf)
 *
 *  @return the node with smallest lower bound from the tree (child, sibling, or leaf)
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jlong JNISCIP(getBestboundNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) (size_t) SCIPgetBestboundNode(scip);
}

/** cuts off node and whole sub tree from branch and bound tree
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(cutoffNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnode               /**< node that should be cut off */
   )
{
   SCIP* scip;
   SCIP_NODE* node;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   node = (SCIP_NODE*) (size_t) jnode;
   assert(node != NULL);

   JNISCIP_CALL( SCIPcutoffNode(scip, node) );
}

/** marks the given node to be propagated again the next time a node of its subtree is processed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(repropagateNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnode               /**< node that should be propagated again */
   )
{
   SCIP* scip;
   SCIP_NODE* node;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   node = (SCIP_NODE*) (size_t) jnode;
   assert(node != NULL);

   JNISCIP_CALL( SCIPrepropagateNode(scip, node) );
}

/** returns depth of first node in active path that is marked being cutoff
 *
 *  @return depth of first node in active path that is marked being cutoff
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(getCutoffdepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);


   return (jint) SCIPgetCutoffdepth(scip);
}

/** returns depth of first node in active path that has to be propagated again
 *
 *  @return depth of first node in active path that has to be propagated again
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(getRepropdepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);


   return (jint) SCIPgetRepropdepth(scip);
}

/** prints all branching decisions on variables from the root to the given node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(printNodeRootPath)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jnode,              /**< node data */
   jlong                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;
   SCIP_NODE* node;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   node = (SCIP_NODE*) (size_t) jnode;
   assert(node != NULL);

   JNISCIP_CALL( SCIPprintNodeRootPath(scip, node, (FILE*)file) );
}

/** gets number of branch and bound runs performed, including the current run
 *
 *  @return the number of branch and bound runs performed, including the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jint JNISCIP(getNRuns)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNRuns(scip);
}

/** gets number of processed nodes in current run, including the focus node
 *
 *  @return the number of processed nodes in current run, including the focus node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jlong JNISCIP(getNNodes)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNNodes(scip);
}

/** gets total number of processed nodes in all runs, including the focus node
 *
 *  @return the total number of processed nodes in all runs, including the focus node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
jlong JNISCIP(getNTotalNodes)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNTotalNodes(scip);
}

/** gets number of nodes left in the tree (children + siblings + leaves)
 *
 *  @return the number of nodes left in the tree (children + siblings + leaves)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jint JNISCIP(getNNodesLeft)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNNodesLeft(scip);
}

/** gets total number of LPs solved so far
 *
 *  @return the total number of LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNLPs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNLPs(scip);
}

/** gets total number of iterations used so far in primal and dual simplex and barrier algorithm
 *
 *  @return the total number of iterations used so far in primal and dual simplex and barrier algorithm
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNLPIterations(scip);
}

/** gets total number of iterations used so far in primal and dual simplex and barrier algorithm for the root node
 *
 *  @return the total number of iterations used so far in primal and dual simplex and barrier algorithm for the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNRootLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNRootLPIterations(scip);
}

/** gets total number of primal LPs solved so far
 *
 *  @return the total number of primal LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNPrimalLPs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNPrimalLPs(scip);
}

/** gets total number of iterations used so far in primal simplex
 *
 *  @return total number of iterations used so far in primal simplex
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNPrimalLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNPrimalLPIterations(scip);
}

/** gets total number of dual LPs solved so far
 *
 *  @return the total number of dual LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNDualLPs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNDualLPs(scip);
}

/** gets total number of iterations used so far in dual simplex
 *
 *  @return the total number of iterations used so far in dual simplex
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNDualLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNDualLPIterations(scip);
}

/** gets total number of barrier LPs solved so far
 *
 *  @return the total number of barrier LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNBarrierLPs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNBarrierLPs(scip);
}

/** gets total number of iterations used so far in barrier algorithm
 *
 *  @return the total number of iterations used so far in barrier algorithm
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNBarrierLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNBarrierLPIterations(scip);
}

/** gets total number of LPs solved so far that were resolved from an advanced start basis
 *
 *  @return the total number of LPs solved so far that were resolved from an advanced start basis
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNResolveLPs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNResolveLPs(scip);
}

/** gets total number of simplex iterations used so far in primal and dual simplex calls where an advanced start basis
 *  was available
 *
 *  @return the total number of simplex iterations used so far in primal and dual simplex calls where an advanced start
 *          basis was available
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNResolveLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNResolveLPIterations(scip);
}

/** gets total number of primal LPs solved so far that were resolved from an advanced start basis
 *
 *  @return the total number of primal LPs solved so far that were resolved from an advanced start basis
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNPrimalResolveLPs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNPrimalResolveLPs(scip);
}

/** gets total number of simplex iterations used so far in primal simplex calls where an advanced start basis
 *  was available
 *
 *  @return the total number of simplex iterations used so far in primal simplex calls where an advanced start
 *          basis was available
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNPrimalResolveLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNPrimalResolveLPIterations(scip);
}

/** gets total number of dual LPs solved so far that were resolved from an advanced start basis
 *
 *  @return the total number of dual LPs solved so far that were resolved from an advanced start basis
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNDualResolveLPs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNDualResolveLPs(scip);
}

/** gets total number of simplex iterations used so far in dual simplex calls where an advanced start basis
 *  was available
 *
 *  @return the total number of simplex iterations used so far in dual simplex calls where an advanced start
 *          basis was available
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNDualResolveLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNDualResolveLPIterations(scip);
}

/** gets total number of LPs solved so far for node relaxations
 *
 *  @return the total number of LPs solved so far for node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNNodeLPs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNNodeLPs(scip);
}

/** gets total number of simplex iterations used so far for node relaxations
 *
 *  @return the total number of simplex iterations used so far for node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNNodeLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNNodeLPIterations(scip);
}

/** gets total number of LPs solved so far for initial LP in node relaxations
 *
 *  @return the total number of LPs solved so far for initial LP in node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNNodeInitLPs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNNodeInitLPs(scip);
}

/** gets total number of simplex iterations used so far for initial LP in node relaxations
 *
 *  @return the total number of simplex iterations used so far for initial LP in node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNNodeInitLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNNodeInitLPIterations(scip);
}

/** gets total number of LPs solved so far during diving and probing
 *
 *  @return total number of LPs solved so far during diving and probing
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNDivingLPs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNDivingLPs(scip);
}

/** gets total number of simplex iterations used so far during diving and probing
 *
 *  @return the total number of simplex iterations used so far during diving and probing
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNDivingLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNDivingLPIterations(scip);
}

/** gets total number of times, strong branching was called (each call represents solving two LPs)
 *
 *  @return the total number of times, strong branching was called (each call represents solving two LPs)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNStrongbranchs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNStrongbranchs(scip);
}

/** gets total number of simplex iterations used so far in strong branching
 *
 *  @return the total number of simplex iterations used so far in strong branching
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNStrongbranchLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNStrongbranchLPIterations(scip);
}

/** gets total number of times, strong branching was called at the root node (each call represents solving two LPs)
 *
 *  @return the total number of times, strong branching was called at the root node (each call represents solving two LPs)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNRootStrongbranchs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNRootStrongbranchs(scip);
}

/** gets total number of simplex iterations used so far in strong branching at the root node
 *
 *  @return the total number of simplex iterations used so far in strong branching at the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jlong JNISCIP(getNRootStrongbranchLPIterations)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNRootStrongbranchLPIterations(scip);
}

/** gets number of pricing rounds performed so far at the current node
 *
 *  @return the number of pricing rounds performed so far at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(getNPriceRounds)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPriceRounds(scip);
}

/** get current number of variables in the pricing store
 *
 *  @return the current number of variables in the pricing store
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jint JNISCIP(getNPricevars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPricevars(scip);
}

/** get total number of pricing variables found so far
 *
 *  @return the total number of pricing variables found so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jint JNISCIP(getNPricevarsFound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPricevarsFound(scip);
}

/** get total number of pricing variables applied to the LPs
 *
 *  @return the total number of pricing variables applied to the LPs
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jint JNISCIP(getNPricevarsApplied)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNPricevarsApplied(scip);
}

/** gets number of separation rounds performed so far at the current node
 *
 *  @return the number of separation rounds performed so far at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(getNSepaRounds)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNSepaRounds(scip);
}

/** get total number of cuts found so far
 *
 *  @return the total number of cuts found so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jint JNISCIP(getNCutsFound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNCutsFound(scip);
}

/** get number of cuts found so far in current separation round
 *
 *  @return the number of cuts found so far in current separation round
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jint JNISCIP(getNCutsFoundRound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNCutsFoundRound(scip);
}

/** get total number of cuts applied to the LPs
 *
 *  @return the total number of cuts applied to the LPs
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jint JNISCIP(getNCutsApplied)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNCutsApplied(scip);
}

/** get total number of constraints found in conflict analysis (conflict and reconvergence constraints)
 *
 *  @return the total number of constraints found in conflict analysis (conflict and reconvergence constraints)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jlong JNISCIP(getNConflictConssFound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNConflictConssFound(scip);
}

/** get number of conflict constraints found so far at the current node
 *
 *  @return the number of conflict constraints found so far at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jint JNISCIP(getNConflictConssFoundNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNConflictConssFoundNode(scip);
}

/** get total number of conflict constraints added to the problem
 *
 *  @return the total number of conflict constraints added to the problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jlong JNISCIP(getNConflictConssApplied)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNConflictConssApplied(scip);
}

/** gets depth of current node, or -1 if no current node exists; in probing, the current node is the last probing node,
 *  such that the depth includes the probing path
 *
 *  @return the depth of current node, or -1 if no current node exists; in probing, the current node is the last probing node,
 *  such that the depth includes the probing path
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jint JNISCIP(getDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetDepth(scip);
}

/** gets depth of the focus node, or -1 if no focus node exists; the focus node is the currently processed node in the
 *  branching tree, excluding the nodes of the probing path
 *
 *  @return the depth of the focus node, or -1 if no focus node exists; the focus node is the currently processed node in the
 *  branching tree, excluding the nodes of the probing path
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jint JNISCIP(getFocusDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetFocusDepth(scip);
}

/** gets maximal depth of all processed nodes in current branch and bound run (excluding probing nodes)
 *
 *  @return the maximal depth of all processed nodes in current branch and bound run (excluding probing nodes)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jint JNISCIP(getMaxDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetMaxDepth(scip);
}

/** gets maximal depth of all processed nodes over all branch and bound runs
 *
 *  @return the maximal depth of all processed nodes over all branch and bound runs
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jint JNISCIP(getMaxTotalDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetMaxTotalDepth(scip);
}

/** gets total number of backtracks, i.e. number of times, the new node was selected from the leaves queue
 *
 *  @return the total number of backtracks, i.e. number of times, the new node was selected from the leaves queue
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jlong JNISCIP(getNBacktracks)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNBacktracks(scip);
}

/** gets current plunging depth (successive times, a child was selected as next node)
 *
 *  @return the current plunging depth (successive times, a child was selected as next node)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(getPlungeDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetPlungeDepth(scip);
}

/** gets total number of active constraints at the current node
 *
 *  @return the total number of active constraints at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(getNActiveConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNActiveConss(scip);
}

/** gets total number of enabled constraints at the current node
 *
 *  @return the total number of enabled constraints at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
jint JNISCIP(getNEnabledConss)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNEnabledConss(scip);
}

/** gets average dual bound of all unprocessed nodes for original problem
 *
 *  @return the average dual bound of all unprocessed nodes for original problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgDualbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgDualbound(scip);
}

/** gets average lower (dual) bound of all unprocessed nodes in transformed problem
 *
 *  @return the average lower (dual) bound of all unprocessed nodes in transformed problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgLowerbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgLowerbound(scip);
}

/** gets global dual bound */
JNIEXPORT
jdouble JNISCIP(getDualbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real dualbound;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   dualbound = SCIPgetDualbound(scip);

   return (jdouble)dualbound;
}

/** gets global lower (dual) bound in transformed problem
 *
 *  @return the global lower (dual) bound in transformed problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getLowerbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetLowerbound(scip);
}

/** gets dual bound of the root node for the original problem
 *
 *  @return the dual bound of the root node for the original problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getDualboundRoot)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetDualboundRoot(scip);
}

/** gets lower (dual) bound in transformed problem of the root node
 *
 *  @return the lower (dual) bound in transformed problem of the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getLowerboundRoot)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetLowerboundRoot(scip);
}

/** gets global primal bound (objective value of best solution or user objective limit) */
JNIEXPORT
jdouble JNISCIP(getPrimalbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real primalbound;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   primalbound = SCIPgetPrimalbound(scip);

   return (jdouble)primalbound;
}

/** gets global upper (primal) bound in transformed problem (objective value of best solution or user objective limit)
 *
 *  @return the global upper (primal) bound in transformed problem (objective value of best solution or user objective limit)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jdouble JNISCIP(getUpperbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetUpperbound(scip);
}

/** gets global cutoff bound in transformed problem: a sub problem with lower bound larger than the cutoff
 *  cannot contain a better feasible solution; usually, this bound is equal to the upper bound, but if the
 *  objective value is always integral, the cutoff bound is (nearly) one less than the upper bound;
 *  additionally, due to objective function domain propagation, the cutoff bound can be further reduced
 *
 *  @return global cutoff bound in transformed problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jdouble JNISCIP(getCutoffbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetCutoffbound(scip);
}

/** updates the cutoff bound
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note the given cutoff bound has to better or equal to known one (SCIPgetCutoffbound())
 */
JNIEXPORT
void JNISCIP(updateCutoffbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               cutoffbound         /**< new cutoff bound */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPupdateCutoffbound(scip, (SCIP_Real)cutoffbound) );
}

/** returns whether the current primal bound is justified with a feasible primal solution; if not, the primal bound
 *  was set from the user as objective limit
 *
 *  @return TRUE if the current primal bound is justified with a feasible primal solution, otherwise FALSE
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jboolean JNISCIP(isPrimalboundSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisPrimalboundSol(scip);
}

/** gets current gap |(primalbound - dualbound)/min(|primalbound|,|dualbound|)| if both bounds have same sign,
 *  or infinity, if they have opposite sign
 *
 *  @return the current gap |(primalbound - dualbound)/min(|primalbound|,|dualbound|)| if both bounds have same sign,
 *  or infinity, if they have opposite sign
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getGap)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetGap(scip);
}

/** gets current gap |(upperbound - lowerbound)/min(|upperbound|,|lowerbound|)| in transformed problem if both bounds
 *  have same sign, or infinity, if they have opposite sign
 *
 *  @return current gap |(upperbound - lowerbound)/min(|upperbound|,|lowerbound|)| in transformed problem if both bounds
 *  have same sign, or infinity, if they have opposite sign
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getTransGap)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetTransGap(scip);
}

/** gets number of feasible primal solutions found so far
 *
 *  @return the number of feasible primal solutions found so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jlong JNISCIP(getNSolsFound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNSolsFound(scip);
}

/** gets number of feasible primal solutions found so far, that improved the primal bound at the time they were found
 *
 *  @return the number of feasible primal solutions found so far, that improved the primal bound at the time they were found
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
JNIEXPORT
jlong JNISCIP(getNBestSolsFound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) SCIPgetNBestSolsFound(scip);
}

/** gets the average pseudo cost value for the given direction over all variables
 *
 *  @return the average pseudo cost value for the given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgPseudocost)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgPseudocost(scip, (SCIP_Real)solvaldelta);
}

/** gets the average pseudo cost value for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @return the average pseudo cost value for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgPseudocostCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgPseudocostCurrentRun(scip, (SCIP_Real)solvaldelta);
}

/** gets the average number of pseudo cost updates for the given direction over all variables
 *
 *  @return the average number of pseudo cost updates for the given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgPseudocostCount)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgPseudocostCount(scip, (SCIP_BRANCHDIR)dir);
}

/** gets the average number of pseudo cost updates for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @return the average number of pseudo cost updates for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgPseudocostCountCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgPseudocostCountCurrentRun(scip, (SCIP_BRANCHDIR)dir);
}

/** gets the average pseudo cost score value over all variables, assuming a fractionality of 0.5
 *
 *  @return the average pseudo cost score value over all variables, assuming a fractionality of 0.5
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgPseudocostScore)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgPseudocostScore(scip);
}

/** gets the average pseudo cost score value over all variables, assuming a fractionality of 0.5,
 *  only using the pseudo cost information of the current run
 *
 *  @return the average pseudo cost score value over all variables, assuming a fractionality of 0.5,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgPseudocostScoreCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgPseudocostScoreCurrentRun(scip);
}

/** gets the average conflict score value over all variables
 *
 *  @return the average conflict score value over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgConflictScore)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgConflictScore(scip);
}

/** gets the average conflict score value over all variables, only using the pseudo cost information of the current run
 *
 *  @return the average conflict score value over all variables, only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgConflictScoreCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgConflictScoreCurrentRun(scip);
}

/** gets the average inference score value over all variables
 *
 *  @return the average inference score value over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgConflictlengthScore)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgConflictlengthScore(scip);
}

/** gets the average conflictlength score value over all variables, only using the pseudo cost information of the
 *  current run
 *
 *  @return the average conflictlength score value over all variables, only using the pseudo cost information of the
 *          current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgConflictlengthScoreCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgConflictlengthScoreCurrentRun(scip);
}

/** returns the average number of inferences found after branching in given direction over all variables
 *
 *  @return the average number of inferences found after branching in given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgInferences)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgInferences(scip, (SCIP_BRANCHDIR)dir);
}

/** returns the average number of inferences found after branching in given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @return the average number of inferences found after branching in given direction over all variables,
 *          only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgInferencesCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgInferencesCurrentRun(scip, (SCIP_BRANCHDIR)dir);
}

/** gets the average inference score value over all variables
 *
 *  @return the average inference score value over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgInferenceScore)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgInferenceScore(scip);
}

/** gets the average inference score value over all variables, only using the inference information information of the
 *  current run
 *
 *  @return the average inference score value over all variables, only using the inference information information of the
 *          current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgInferenceScoreCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgInferenceScoreCurrentRun(scip);
}

/** returns the average number of cutoffs found after branching in given direction over all variables
 *
 *  @return the average number of cutoffs found after branching in given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgCutoffs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgCutoffs(scip, (SCIP_BRANCHDIR)dir);
}

/** returns the average number of cutoffs found after branching in given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @return the average number of cutoffs found after branching in given direction over all variables,
 *          only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgCutoffsCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgCutoffsCurrentRun(scip, (SCIP_BRANCHDIR)dir);
}

/** gets the average cutoff score value over all variables
 *
 *  @return the average cutoff score value over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgCutoffScore)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgCutoffScore(scip);
}

/** gets the average cutoff score value over all variables, only using the pseudo cost information of the current run
 *
 *  @return the average cutoff score value over all variables, only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jdouble JNISCIP(getAvgCutoffScoreCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetAvgCutoffScoreCurrentRun(scip);
}

/** outputs original problem to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(printOrigProblem)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 file,               /**< output file (or NULL for standard output) */
   jstring               jextension,         /**< file format (or NULL for default CIP format)*/
   jboolean              genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP* scip;
   const char* extension;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   if( jextension != NULL )
      extension = (*env)->GetStringUTFChars(env, jextension, &iscopy);
   else
      extension = NULL;

   JNISCIP_CALL( SCIPprintOrigProblem(scip, (FILE*)(size_t)file, extension, (SCIP_Bool)genericnames) );
}

/** outputs transformed problem of the current node to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
JNIEXPORT
void JNISCIP(printTransProblem)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 file,               /**< output file (or NULL for standard output) */
   jstring               jextension,         /**< file format (or NULL for default CIP format)*/
   jboolean              genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP* scip;
   const char* extension;
   jboolean iscopy;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* convert JNI string into C const char* */
   if( jextension != NULL )
      extension = (*env)->GetStringUTFChars(env, jextension, &iscopy);
   else
      extension = NULL;

   JNISCIP_CALL( SCIPprintTransProblem(scip, (FILE*)(size_t)file, extension, (SCIP_Bool)genericnames) );
}

/** outputs solving statistics
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
void JNISCIP(printStatistics)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPprintStatistics(scip, (FILE*)(size_t)file) );
}

/** outputs history statistics about branchings on variables
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
void JNISCIP(printBranchingStatistics)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPprintBranchingStatistics(scip, (FILE*)(size_t)file) );
}

/** outputs node information display line
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
JNIEXPORT
void JNISCIP(printDisplayLine)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 file,               /**< output file (or NULL for standard output) */
   jint                  verblevel,          /**< minimal verbosity level to actually display the information line */
   jboolean              endline             /**< should the line be terminated with a newline symbol? */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPprintDisplayLine(scip, (FILE*)(size_t)file, (SCIP_VERBLEVEL)verblevel, (SCIP_Bool)endline) );
}

/** gets total number of implications between variables that are stored in the implication graph
 *
 *  @return the total number of implications between variables that are stored in the implication graph
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
JNIEXPORT
jint JNISCIP(getNImplications)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint) SCIPgetNImplications(scip);
}

/** gets current time of day in seconds (standard time zone)
 *
 *  @return the current time of day in seconds (standard time zone).
 */
JNIEXPORT
jdouble JNISCIP(getTimeOfDay)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetTimeOfDay(scip);
}

/** creates a clock using the default clock type
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
jlong JNISCIP(createClock)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_CLOCK* clock;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPcreateClock(scip, &clock) );

   return (jlong) (size_t) clock;
}

/** creates a clock counting the CPU user seconds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
jlong JNISCIP(createCPUClock)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_CLOCK* clock;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPcreateCPUClock(scip, &clock) );

   return (jlong) (size_t) clock;
}

/** creates a clock counting the wall clock seconds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
jlong JNISCIP(createWallClock)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_CLOCK* clock;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPcreateWallClock(scip, &clock) );

   return (jlong) (size_t) clock;
}

/** frees a clock
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(freeClock)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jclock               /**< clock timer */
   )
{
   SCIP* scip;
   SCIP_CLOCK* clock;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   clock = (SCIP_CLOCK*) (size_t) jclock;
   assert(clock != NULL);

   JNISCIP_CALL( SCIPfreeClock(scip, &clock) );
}

/** resets the time measurement of a clock to zero and completely stops the clock
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(resetClock)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jclock              /**< clock timer */
   )
{
   SCIP* scip;
   SCIP_CLOCK* clock;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   clock = (SCIP_CLOCK*) (size_t) jclock;
   assert(clock != NULL);

   JNISCIP_CALL( SCIPresetClock(scip, clock) );
}

/** starts the time measurement of a clock
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(startClock)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jclock              /**< clock timer */
   )
{
   SCIP* scip;
   SCIP_CLOCK* clock;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   clock = (SCIP_CLOCK*) (size_t) jclock;
   assert(clock != NULL);

   JNISCIP_CALL( SCIPstartClock(scip, clock) );
}

/** stops the time measurement of a clock
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(stopClock)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jclock              /**< clock timer */
   )
{
   SCIP* scip;
   SCIP_CLOCK* clock;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   clock = (SCIP_CLOCK*) (size_t) jclock;
   assert(clock != NULL);

   JNISCIP_CALL( SCIPstopClock(scip, clock) );
}

/** starts the current solving time
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
void JNISCIP(startSolvingTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPstartSolvingTime(scip) );
}

/** stops the current solving time in seconds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
void JNISCIP(stopSolvingTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPstopSolvingTime(scip) );
}

/** gets the measured time of a clock in seconds
 *
 *  @return the measured time of a clock in seconds.
 */
JNIEXPORT
jdouble JNISCIP(getClockTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jclock              /**< clock timer */
   )
{
   SCIP* scip;
   SCIP_CLOCK* clock;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   clock = (SCIP_CLOCK*) (size_t) jclock;
   assert(clock != NULL);

   return (jdouble) SCIPgetClockTime(scip, clock);
}

/** sets the measured time of a clock to the given value in seconds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(setClockTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jclock,             /**< clock timer */
   jdouble               sec                 /**< time in seconds to set the clock's timer to */
   )
{
   SCIP* scip;
   SCIP_CLOCK* clock;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   clock = (SCIP_CLOCK*) (size_t) jclock;
   assert(clock != NULL);

   JNISCIP_CALL( SCIPsetClockTime(scip, clock, (SCIP_Real)sec) );
}

/** gets the current total SCIP time in seconds
 *
 *  @return the current total SCIP time in seconds.
 */
JNIEXPORT
jdouble JNISCIP(getTotalTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetTotalTime(scip);
}

/** getsf the current solving time in seconds
 *
 *  @return the current solving time in seconds.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jdouble JNISCIP(getSolvingTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetSolvingTime(scip);
}

/** gets the current reading time in seconds
 *
 *  @return the current reading time in seconds.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jdouble JNISCIP(getReadingTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetReadingTime(scip);
}

/** gets the current presolving time in seconds
 *
 *  @return the current presolving time in seconds.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
JNIEXPORT
jdouble JNISCIP(getPresolvingTime)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetPresolvingTime(scip);
}

/** returns value treated as zero */
JNIEXPORT
jdouble JNISCIP(epsilon)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real val;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   val = SCIPepsilon(scip);

   return (jdouble)val;
}

/** returns value treated as zero for sums of floating point values */
JNIEXPORT
jdouble JNISCIP(sumepsilon)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real val;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   val = SCIPsumepsilon(scip);

   return (jdouble)val;
}

/** returns feasibility tolerance for constraints */
JNIEXPORT
jdouble JNISCIP(feastol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real val;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   val = SCIPfeastol(scip);

   return (jdouble)val;
}

/** returns feasibility tolerance for reduced costs */
JNIEXPORT
jdouble JNISCIP(dualfeastol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real val;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   val = SCIPdualfeastol(scip);

   return val;
}

/** returns convergence tolerance used in barrier algorithm */
JNIEXPORT
jdouble JNISCIP(barrierconvtol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real val;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   val = SCIPbarrierconvtol(scip);

   return val;
}

/** return the cutoff bound delta
 *
 *  @return cutoff bound data
 */
JNIEXPORT
jdouble JNISCIP(cutoffbounddelta)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPcutoffbounddelta(scip);
}

/** sets the feasibility tolerance for constraints
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(chgFeastol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               feastol             /**< new feasibility tolerance for constraints */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPchgFeastol(scip, (SCIP_Real)feastol) );
}

/** sets the primal feasibility tolerance of LP solver
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(chgLpfeastol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               lpfeastol,          /**< new primal feasibility tolerance of LP solver */
   jboolean              printnewvalue       /**< should "numerics/lpfeastol = ..." be printed? */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPchgLpfeastol(scip, (SCIP_Real)lpfeastol, (SCIP_Bool)printnewvalue) );
}

/** sets the feasibility tolerance for reduced costs
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(chgDualfeastol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               dualfeastol         /**< new feasibility tolerance for reduced costs */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPchgDualfeastol(scip, (SCIP_Real)dualfeastol) );
}

/** sets the convergence tolerance used in barrier algorithm
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(chgBarrierconvtol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               barrierconvtol      /**< new convergence tolerance used in barrier algorithm */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPchgBarrierconvtol(scip, (SCIP_Real)barrierconvtol) );
}

/** marks that some limit parameter was changed */
JNIEXPORT
void JNISCIP(markLimitChanged)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   SCIPmarkLimitChanged(scip);
}

/** returns value treated as infinity */
JNIEXPORT
jdouble JNISCIP(infinity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_Real val;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   val = SCIPinfinity(scip);

   return val;
}

/** checks, if values are in range of epsilon */
JNIEXPORT
jboolean JNISCIP(isEQ)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisEQ(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if val1 is (more than epsilon) lower than val2 */
JNIEXPORT
jboolean JNISCIP(isLT)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisLT(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if val1 is not (more than epsilon) greater than val2 */
JNIEXPORT
jboolean JNISCIP(isLE)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisLE(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if val1 is (more than epsilon) greater than val2 */
JNIEXPORT
jboolean JNISCIP(isGT)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisGT(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if val1 is not (more than epsilon) lower than val2 */
JNIEXPORT
jboolean JNISCIP(isGE)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisGE(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if value is (positive) infinite */
JNIEXPORT
jboolean JNISCIP(isInfinity)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisInfinity(scip, (SCIP_Real)val);
}

/** checks, if value is in range epsilon of 0.0 */
JNIEXPORT
jboolean JNISCIP(isZero)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisZero(scip, (SCIP_Real)val);
}

/** checks, if value is greater than epsilon */
JNIEXPORT
jboolean JNISCIP(isPositive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisPositive(scip, (SCIP_Real)val);
}

/** checks, if value is lower than -epsilon */
JNIEXPORT
jboolean JNISCIP(isNegative)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisNegative(scip, (SCIP_Real)val);
}

/** checks, if value is integral within epsilon */
JNIEXPORT
jboolean JNISCIP(isIntegral)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisIntegral(scip, (SCIP_Real)val);
}

/** checks whether the product val * scalar is integral in epsilon scaled by scalar */
JNIEXPORT
jboolean JNISCIP(isScalingIntegral)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val,                /**< unscaled value to check for scaled integrality */
   jdouble               scalar              /**< value to scale val with for checking for integrality */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisScalingIntegral(scip, (SCIP_Real)val, (SCIP_Real)scalar);
}

/** checks, if given fractional part is smaller than epsilon */
JNIEXPORT
jboolean JNISCIP(isFracIntegral)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisFracIntegral(scip, (SCIP_Real)val);
}

/** rounds value + epsilon down to the next integer */
JNIEXPORT
jdouble JNISCIP(floor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPfloor(scip, (SCIP_Real)val);
}

/** rounds value - epsilon up to the next integer */
JNIEXPORT
jdouble JNISCIP(ceil)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPceil(scip, (SCIP_Real)val);
}

/** rounds value to the nearest integer with epsilon tolerance */
JNIEXPORT
jdouble JNISCIP(round)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPround(scip, (SCIP_Real)val);
}

/** returns fractional part of value, i.e. x - floor(x) in epsilon tolerance */
JNIEXPORT
jdouble JNISCIP(frac)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPfrac(scip, (SCIP_Real)val);
}

/** checks, if values are in range of sumepsilon */
JNIEXPORT
jboolean JNISCIP(isSumEQ)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisSumEQ(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if val1 is (more than sumepsilon) lower than val2 */
JNIEXPORT
jboolean JNISCIP(isSumLT)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisSumLT(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
JNIEXPORT
jboolean JNISCIP(isSumLE)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisSumLE(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if val1 is (more than sumepsilon) greater than val2 */
JNIEXPORT
jboolean JNISCIP(isSumGT)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisSumGT(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
JNIEXPORT
jboolean JNISCIP(isSumGE)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisSumGE(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if value is in range sumepsilon of 0.0 */
JNIEXPORT
jboolean JNISCIP(isSumZero)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisSumZero(scip, (SCIP_Real)val);
}

/** checks, if value is greater than sumepsilon */
JNIEXPORT
jboolean JNISCIP(isSumPositive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisSumPositive(scip, (SCIP_Real)val);
}

/** checks, if value is lower than -sumepsilon */
JNIEXPORT
jboolean JNISCIP(isSumNegative)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisSumNegative(scip, (SCIP_Real)val);
}

/** checks, if relative difference of values is in range of feasibility tolerance */
JNIEXPORT
jboolean JNISCIP(isFeasEQ)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisFeasEQ(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if relative difference val1 and val2 is lower than feasibility tolerance */
JNIEXPORT
jboolean JNISCIP(isFeasLT)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisFeasLT(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if relative difference of val1 and val2 is not greater than feasibility tolerance */
JNIEXPORT
jboolean JNISCIP(isFeasLE)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisFeasLE(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if relative difference of val1 and val2 is greater than feastol */
JNIEXPORT
jboolean JNISCIP(isFeasGT)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisFeasGT(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if relative difference of val1 and val2 is not lower than -feastol */
JNIEXPORT
jboolean JNISCIP(isFeasGE)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisFeasGE(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if value is in range feasibility tolerance of 0.0 */
JNIEXPORT
jboolean JNISCIP(isFeasZero)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisFeasZero(scip, (SCIP_Real)val);
}

/** checks, if value is greater than feasibility tolerance */
JNIEXPORT
jboolean JNISCIP(isFeasPositive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisFeasPositive(scip, (SCIP_Real)val);
}

/** checks, if value is lower than -feasibility tolerance */
JNIEXPORT
jboolean JNISCIP(isFeasNegative)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisFeasNegative(scip, (SCIP_Real)val);
}

/** checks, if value is integral within the LP feasibility bounds */
JNIEXPORT
jboolean JNISCIP(isFeasIntegral)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisFeasIntegral(scip, (SCIP_Real)val);
}

/** checks, if given fractional part is smaller than feastol */
JNIEXPORT
jboolean JNISCIP(isFeasFracIntegral)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisFeasFracIntegral(scip, (SCIP_Real)val);
}

/** rounds value + feasibility tolerance down to the next integer */
JNIEXPORT
jdouble JNISCIP(feasFloor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPfeasFloor(scip, (SCIP_Real)val);
}

/** rounds value - feasibility tolerance up to the next integer */
JNIEXPORT
jdouble JNISCIP(feasCeil)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPfeasCeil(scip, (SCIP_Real)val);
}

/** rounds value - feasibility tolerance up to the next integer in feasibility tolerance */
JNIEXPORT
jdouble JNISCIP(feasRound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPfeasRound(scip, (SCIP_Real)val);
}

/** returns fractional part of value, i.e. x - floor(x) */
JNIEXPORT
jdouble JNISCIP(feasFrac)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to process */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPfeasFrac(scip, (SCIP_Real)val);
}

/** checks, if the given new lower bound is tighter (w.r.t. bound strengthening epsilon) than the old one */
JNIEXPORT
jboolean JNISCIP(isLbBetter)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               newlb,              /**< new lower bound */
   jdouble               oldlb,              /**< old lower bound */
   jdouble               oldub               /**< old upper bound */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisLbBetter(scip, (SCIP_Real)newlb, (SCIP_Real)oldlb, (SCIP_Real)oldub);
}

/** checks, if the given new upper bound is tighter (w.r.t. bound strengthening epsilon) than the old one */
JNIEXPORT
jboolean JNISCIP(isUbBetter)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               newub,              /**< new lower bound */
   jdouble               oldlb,              /**< old lower bound */
   jdouble               oldub               /**< old upper bound */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisUbBetter(scip, (SCIP_Real)newub, (SCIP_Real)oldlb, (SCIP_Real)oldub);
}

/** checks, if relative difference of values is in range of epsilon */
JNIEXPORT
jboolean JNISCIP(isRelEQ)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisRelEQ(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if relative difference of val1 and val2 is lower than epsilon */
JNIEXPORT
jboolean JNISCIP(isRelLT)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisRelLT(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
JNIEXPORT
jboolean JNISCIP(isRelLE)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisRelLE(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if relative difference of val1 and val2 is greater than epsilon */
JNIEXPORT
jboolean JNISCIP(isRelGT)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisRelGT(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
JNIEXPORT
jboolean JNISCIP(isRelGE)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisRelGE(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if relative difference of values is in range of sumepsilon */
JNIEXPORT
jboolean JNISCIP(isSumRelEQ)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisSumRelEQ(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if relative difference of val1 and val2 is lower than sumepsilon */
JNIEXPORT
jboolean JNISCIP(isSumRelLT)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisSumRelLT(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if relative difference of val1 and val2 is not greater than sumepsilon */
JNIEXPORT
jboolean JNISCIP(isSumRelGT)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisSumRelGT(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if relative difference of val1 and val2 is greater than -sumepsilon */
JNIEXPORT
jboolean JNISCIP(isSumRelGE)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< first value to be compared */
   jdouble               val2                /**< second value to be compared */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisSumRelGE(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** Checks, if an iteratively updated value is reliable or should be recomputed from scratch.
 *  This is useful, if the value, e.g., the activity of a linear constraint or the pseudo objective value, gets a high
 *  absolute value during the optimization process which is later reduced significantly. In this case, the last digits
 *  were canceled out when increasing the value and are random after decreasing it.
 *  We do not consider the cancellations which can occur during increasing the absolute value because they just cannot
 *  be expressed using fixed precision floating point arithmetic, anymore.
 *  In order to get more reliable values, the idea is to always store the last reliable value, where increasing the
 *  absolute of the value is viewed as preserving reliability. Then, after each update, the new absolute value can be
 *  compared against the last reliable one with this method, checking whether it was decreased by a factor of at least
 *  "lp/recompfac" and should be recomputed.
 */
JNIEXPORT
jboolean JNISCIP(isUpdateUnreliable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val1,               /**< new value after update */
   jdouble               val2                /**< old value, i.e., last reliable value */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisUpdateUnreliable(scip, (SCIP_Real)val1, (SCIP_Real)val2);
}

/** checks, if value is huge and should be handled separately (e.g., in activity computation) */
JNIEXPORT
jboolean JNISCIP(isHugeValue)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jdouble               val                 /**< value to be checked whether it is huge */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jboolean) SCIPisHugeValue(scip, (SCIP_Real)val);
}

/** returns the minimum value that is regarded as huge and should be handled separately (e.g., in activity computation) */
JNIEXPORT
jdouble JNISCIP(getHugeValue)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jdouble) SCIPgetHugeValue(scip);
}

/** outputs a real number, or "+infinity", or "-infinity" to a file */
JNIEXPORT
void JNISCIP(printReal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 file,               /**< output file (or NULL for standard output) */
   jdouble               val,                /**< value to print */
   jint                  width,              /**< width of the field */
   jint                  precision           /**< number of significant digits printed */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   SCIPprintReal(scip, (FILE*)(size_t)file, (SCIP_Real)val, (int)width, (int)precision);
}

/** returns block memory to use at the current time
 *
 *  @return the block memory to use at the current time.
 */
JNIEXPORT
jlong JNISCIP(blkmem)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong) (size_t) SCIPblkmem(scip);
}

/** returns the total number of bytes used in block memory
 *
 *  @return the total number of bytes used in block memory.
 */
JNIEXPORT
jlong JNISCIP(getMemUsed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong)  SCIPgetMemUsed(scip);
}

/** returns the estimated number of bytes used by external software, e.g., the LP solver
 *
 *  @return the estimated number of bytes used by external software, e.g., the LP solver.
 */
JNIEXPORT
jlong JNISCIP(getMemExternEstim)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jlong)  SCIPgetMemExternEstim(scip);
}

/** calculate memory size for dynamically allocated arrays
 *
 *  @return the memory size for dynamically allocated arrays.
 */
JNIEXPORT
jint JNISCIP(calcMemGrowSize)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  num                 /**< minimum number of entries to store */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   return (jint)  SCIPcalcMemGrowSize(scip, (int)num);
}

/** prints output about used memory */
JNIEXPORT
void JNISCIP(printMemoryDiagnostic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   SCIPprintMemoryDiagnostic(scip);
}

/** creates a dynamic array of real values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
jlong JNISCIP(createRealarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_REALARRAY* realarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPcreateRealarray(scip, &realarray) );

   return (jlong) (size_t) realarray;
}

/** frees a dynamic array of real values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(freeRealarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrealarray          /**< real array */
   )
{
   SCIP* scip;
   SCIP_REALARRAY* realarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   realarray = (SCIP_REALARRAY*) (size_t) jrealarray;
   assert(realarray != NULL);

   SCIPfreeRealarray(scip, &realarray);
}

/** extends dynamic array to be able to store indices from minidx to maxidx
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(extendRealarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrealarray,         /**< dynamic real array */
   jint                  minidx,             /**< smallest index to allocate storage for */
   jint                  maxidx              /**< largest index to allocate storage for */
   )
{
   SCIP* scip;
   SCIP_REALARRAY* realarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   realarray = (SCIP_REALARRAY*) (size_t) jrealarray;
   assert(realarray != NULL);

   JNISCIP_CALL( SCIPextendRealarray(scip, realarray, (int)minidx, (int)maxidx) );
}

/** clears a dynamic real array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(clearRealarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrealarray          /**< dynamic real array */
   )
{
   SCIP* scip;
   SCIP_REALARRAY* realarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   realarray = (SCIP_REALARRAY*) (size_t) jrealarray;
   assert(realarray != NULL);

   JNISCIP_CALL( SCIPclearRealarray(scip, realarray) );
}

/** gets value of entry in dynamic array
 *
 *  @return  value of entry in dynamic array
 */
JNIEXPORT
jdouble JNISCIP(getRealarrayVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrealarray,         /**< dynamic real array */
   jint                  idx                 /**< array index to get value for */
   )
{
   SCIP* scip;
   SCIP_REALARRAY* realarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   realarray = (SCIP_REALARRAY*) (size_t) jrealarray;
   assert(realarray != NULL);

   return (jdouble) SCIPgetRealarrayVal(scip, realarray, (int)idx);
}

/** sets value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(setRealarrayVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrealarray,         /**< dynamic real array */
   jint                  idx,                /**< array index to set value for */
   jdouble               val                 /**< value to set array index to */
   )
{
   SCIP* scip;
   SCIP_REALARRAY* realarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   realarray = (SCIP_REALARRAY*) (size_t) jrealarray;
   assert(realarray != NULL);

   JNISCIP_CALL( SCIPsetRealarrayVal(scip, realarray, (int)idx, (SCIP_Real)val) );
}

/** increases value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(incRealarrayVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrealarray,         /**< dynamic real array */
   jint                  idx,                /**< array index to set value for */
   jdouble               val                 /**< value to increase array index */
   )
{
   SCIP* scip;
   SCIP_REALARRAY* realarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   realarray = (SCIP_REALARRAY*) (size_t) jrealarray;
   assert(realarray != NULL);

   JNISCIP_CALL( SCIPincRealarrayVal(scip, realarray, (int)idx, (SCIP_Real)val) );
}

/** returns the minimal index of all stored non-zero elements
 *
 *  @return the minimal index of all stored non-zero elements
 */
JNIEXPORT
jint JNISCIP(getRealarrayMinIdx)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrealarray         /**< dynamic real array */
   )
{
   SCIP* scip;
   SCIP_REALARRAY* realarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   realarray = (SCIP_REALARRAY*) (size_t) jrealarray;
   assert(realarray != NULL);

   return (int) SCIPgetRealarrayMinIdx(scip, realarray);
}

/** returns the maximal index of all stored non-zero elements
 *
 *  @return the maximal index of all stored non-zero elements
 */
JNIEXPORT
jint JNISCIP(getRealarrayMaxIdx)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jrealarray          /**< dynamic real array */
   )
{
   SCIP* scip;
   SCIP_REALARRAY* realarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   realarray = (SCIP_REALARRAY*) (size_t) jrealarray;
   assert(realarray != NULL);

   return (int) SCIPgetRealarrayMaxIdx(scip, realarray);
}

/** creates a dynamic array of int values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
jlong JNISCIP(createIntarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_INTARRAY* intarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPcreateIntarray(scip, &intarray) );

   return (jlong) (size_t) intarray;
}

/** frees a dynamic array of int values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(freeIntarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jintarray           /**< dynamic int array */
   )
{
   SCIP* scip;
   SCIP_INTARRAY* intarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   intarray = (SCIP_INTARRAY*) (size_t) jintarray;
   assert(intarray != NULL);

   SCIPfreeIntarray(scip, &intarray);
}

/** extends dynamic array to be able to store indices from minidx to maxidx
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(extendIntarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jintarray,          /**< dynamic int array */
   jint                  minidx,             /**< smallest index to allocate storage for */
   jint                  maxidx              /**< largest index to allocate storage for */
   )
{
   SCIP* scip;
   SCIP_INTARRAY* intarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   intarray = (SCIP_INTARRAY*) (size_t) jintarray;
   assert(intarray != NULL);

   JNISCIP_CALL( SCIPextendIntarray(scip, intarray, (int)minidx, (int)maxidx) );
}

/** clears a dynamic int array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(clearIntarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jintarray           /**< dynamic int array */
   )
{
   SCIP* scip;
   SCIP_INTARRAY* intarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   intarray = (SCIP_INTARRAY*) (size_t) jintarray;
   assert(intarray != NULL);

   JNISCIP_CALL( SCIPclearIntarray(scip, intarray) );
}

/** gets value of entry in dynamic array
 *
 *  @return value of entry in dynamic array
 */
JNIEXPORT
jint JNISCIP(getIntarrayVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jintarray,          /**< dynamic int array */
   jint                  idx                 /**< array index to get value for */
   )
{
   SCIP* scip;
   SCIP_INTARRAY* intarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   intarray = (SCIP_INTARRAY*) (size_t) jintarray;
   assert(intarray != NULL);

   return (jint) SCIPgetIntarrayVal(scip, intarray, (int)idx);
}

/** sets value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(setIntarrayVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jintarray,          /**< dynamic int array */
   jint                  idx,                /**< array index to set value for */
   jint                  val                 /**< value to set array index to */
   )
{
   SCIP* scip;
   SCIP_INTARRAY* intarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   intarray = (SCIP_INTARRAY*) (size_t) jintarray;
   assert(intarray != NULL);

   JNISCIP_CALL( SCIPsetIntarrayVal(scip, intarray, (int)idx, (int)val) );
}

/** increases value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(incIntarrayVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jintarray,          /**< dynamic int array */
   jint                  idx,                /**< array index to increase value for */
   jint                  incval              /**< value to increase array index */
   )
{
   SCIP* scip;
   SCIP_INTARRAY* intarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   intarray = (SCIP_INTARRAY*) (size_t) jintarray;
   assert(intarray != NULL);

   JNISCIP_CALL( SCIPincIntarrayVal(scip, intarray, (int)idx, (int)incval) );
}

/** returns the minimal index of all stored non-zero elements
 *
 *  @return the minimal index of all stored non-zero elements
 */
JNIEXPORT
jint JNISCIP(getIntarrayMinIdx)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jintarray           /**< dynamic int array */
   )
{
   SCIP* scip;
   SCIP_INTARRAY* intarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   intarray = (SCIP_INTARRAY*) (size_t) jintarray;
   assert(intarray != NULL);

   return (jint) SCIPgetIntarrayMinIdx(scip, intarray);
}

/** returns the maximal index of all stored non-zero elements
 *
 *  @return the maximal index of all stored non-zero elements
 */
JNIEXPORT
jint JNISCIP(getIntarrayMaxIdx)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jintarray           /**< dynamic int array */
   )
{
   SCIP* scip;
   SCIP_INTARRAY* intarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   intarray = (SCIP_INTARRAY*) (size_t) jintarray;
   assert(intarray != NULL);

   return (jint) SCIPgetIntarrayMaxIdx(scip, intarray);
}

/** creates a dynamic array of bool values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
jlong JNISCIP(createBoolarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_BOOLARRAY* boolarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPcreateBoolarray(scip, &boolarray) );

   return (jlong) (size_t) boolarray;
}

/** frees a dynamic array of bool values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(freeBoolarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jboolarray          /**< bool array */
   )
{
   SCIP* scip;
   SCIP_BOOLARRAY* boolarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   boolarray = (SCIP_BOOLARRAY*) (size_t) jboolarray;
   assert(boolarray != NULL);

   SCIPfreeBoolarray(scip, &boolarray);
}

/** extends dynamic array to be able to store indices from minidx to maxidx
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(extendBoolarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jboolarray,         /**< dynamic bool array */
   jint                  minidx,             /**< smallest index to allocate storage for */
   jint                  maxidx              /**< largest index to allocate storage for */
   )
{
   SCIP* scip;
   SCIP_BOOLARRAY* boolarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   boolarray = (SCIP_BOOLARRAY*) (size_t) jboolarray;
   assert(boolarray != NULL);

   JNISCIP_CALL( SCIPextendBoolarray(scip, boolarray, (int)minidx, (int)maxidx) );
}

/** clears a dynamic bool array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(clearBoolarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jboolarray          /**< dynamic bool array */
   )
{
   SCIP* scip;
   SCIP_BOOLARRAY* boolarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   boolarray = (SCIP_BOOLARRAY*) (size_t) jboolarray;
   assert(boolarray != NULL);

   JNISCIP_CALL( SCIPclearBoolarray(scip, boolarray) );
}

/** gets value of entry in dynamic array
 *
 *  @return value of entry in dynamic array at position idx
 */
JNIEXPORT
jboolean JNISCIP(getBoolarrayVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jboolarray,         /**< dynamic bool array */
   jint                  idx                 /**< array index to get value for */
   )
{
   SCIP* scip;
   SCIP_BOOLARRAY* boolarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   boolarray = (SCIP_BOOLARRAY*) (size_t) jboolarray;
   assert(boolarray != NULL);

   return (jboolean) SCIPgetBoolarrayVal(scip, boolarray, (int)idx);
}

/** sets value of entry in dynamic array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(setBoolarrayVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jboolarray,         /**< dynamic bool array */
   jint                  idx,                /**< array index to set value for */
   jboolean              val                 /**< value to set array index to */
   )
{
   SCIP* scip;
   SCIP_BOOLARRAY* boolarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   boolarray = (SCIP_BOOLARRAY*) (size_t) jboolarray;
   assert(boolarray != NULL);

   JNISCIP_CALL( SCIPsetBoolarrayVal(scip, boolarray, (int)idx, (int)val) );
}

/** returns the minimal index of all stored non-zero elements
 *
 *  @return the minimal index of all stored non-zero elements
 */
JNIEXPORT
jint JNISCIP(getBoolarrayMinIdx)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jboolarray          /**< dynamic bool array */
   )
{
   SCIP* scip;
   SCIP_BOOLARRAY* boolarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   boolarray = (SCIP_BOOLARRAY*) (size_t) jboolarray;
   assert(boolarray != NULL);

   return (jint) SCIPgetBoolarrayMinIdx(scip, boolarray);
}

/** returns the maximal index of all stored non-zero elements
 *
 *  @return the maximal index of all stored non-zero elements
 */
JNIEXPORT
jint JNISCIP(getBoolarrayMaxIdx)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jboolarray          /**< dynamic bool array */
   )
{
   SCIP* scip;
   SCIP_BOOLARRAY* boolarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   boolarray = (SCIP_BOOLARRAY*) (size_t) jboolarray;
   assert(boolarray != NULL);

   return (jint) SCIPgetBoolarrayMaxIdx(scip, boolarray);
}

/** creates a dynamic array of pointers
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
jlong JNISCIP(createPtrarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;
   SCIP_PTRARRAY* ptrarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPcreatePtrarray(scip, &ptrarray) );

   return (jlong) (size_t) ptrarray;
}

/** frees a dynamic array of pointers
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(freePtrarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jptrarray           /**< dynamic int array */
   )
{
   SCIP* scip;
   SCIP_PTRARRAY* ptrarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   ptrarray = (SCIP_PTRARRAY*) (size_t) jptrarray;
   assert(ptrarray != NULL);

   JNISCIP_CALL( SCIPfreePtrarray(scip, &ptrarray) );
}

/** extends dynamic array to be able to store indices from minidx to maxidx
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(extendPtrarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jptrarray,          /**< dynamic int array */
   jint                  minidx,             /**< smallest index to allocate storage for */
   jint                  maxidx              /**< largest index to allocate storage for */
   )
{
   SCIP* scip;
   SCIP_PTRARRAY* ptrarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   ptrarray = (SCIP_PTRARRAY*) (size_t) jptrarray;
   assert(ptrarray != NULL);

   JNISCIP_CALL( SCIPextendPtrarray(scip, ptrarray, (int)minidx, (int)maxidx) );
}

/** clears a dynamic pointer array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
JNIEXPORT
void JNISCIP(clearPtrarray)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jptrarray           /**< dynamic int array */
   )
{
   SCIP* scip;
   SCIP_PTRARRAY* ptrarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   ptrarray = (SCIP_PTRARRAY*) (size_t) jptrarray;
   assert(ptrarray != NULL);

   JNISCIP_CALL( SCIPclearPtrarray(scip, ptrarray) );
}

/** gets value of entry in dynamic array */
JNIEXPORT
jlong JNISCIP(getPtrarrayVal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jptrarray,          /**< dynamic int array */
   jint                  idx                 /**< array index to get value for */
   )
{
   SCIP* scip;
   SCIP_PTRARRAY* ptrarray;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   ptrarray = (SCIP_PTRARRAY*) (size_t) jptrarray;
   assert(ptrarray != NULL);

   return (jlong) (size_t) SCIPgetPtrarrayVal(scip, ptrarray, (int)idx);
}


/**
setPtrarrayVal
getPtrarrayMinIdx
getPtrarrayMaxIdx */

/** includes default SCIP plugins into SCIP */
JNIEXPORT
void JNISCIP(includeDefaultPlugins)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeDefaultPlugins(scip) );
}

/** evaluates command line parameters and runs SCIP appropriately in the given SCIP instance */
JNIEXPORT
void JNISCIP(processShellArguments)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jint                  jargc,              /**< number of shell parameters */
   jobjectArray          jargv,              /**< array with shell parameters */
   jstring               jdefaultsetname     /**< name of default settings file */
   )
{
   SCIP* scip;
   char** argv;
   const char* defaultsetname;
   jboolean iscopy;
   int argc;
   int i;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   /* get number of arguments */
   argc = (*env)->GetArrayLength(env, jargv);
   assert(argc == jargc);

   /* allocate memory for the arguments array */
   JNISCIP_CALL( SCIPallocMemoryArray(scip, &argv, argc) );

   for( i = 0; i < argc; ++i)
   {
      jstring jstr;
      const char* str;

      jstr = (jstring) (*env)->GetObjectArrayElement(env, jargv, i);
      str = (*env)->GetStringUTFChars(env, jstr, &iscopy);
      assert(iscopy);

      /* copy string */
      JNISCIP_CALL( SCIPduplicateMemoryArray(scip, &argv[i], str, (int)(strlen(str)+1)));
   }

   /* convert JNI string into C const char* */
   defaultsetname = (*env)->GetStringUTFChars(env, jdefaultsetname, &iscopy);
   if( defaultsetname == NULL )
      SCIPABORT();

   assert(iscopy);

   JNISCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   for( i = 0; i < argc; ++i)
   {
      jstring str;

      str = (jstring) (*env)->GetObjectArrayElement(env, jargv, i);
      (*env)->ReleaseStringUTFChars(env, str, argv[i]);

      SCIPfreeMemoryArray(scip, &argv[i]);
   }

   SCIPfreeMemoryArray(scip, &argv);
}
