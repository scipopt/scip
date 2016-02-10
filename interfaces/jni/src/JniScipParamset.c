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
 * @author Felix Simon
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipParamset.h"
#include "def.h"

#include "scip/scip.h"

#include "JniScipParamtype.h"

/** returns type of parameter */
JNIEXPORT
jint JNISCIPPARAMSET(paramGetType)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
   SCIP_PARAM* param;
   SCIP_PARAMTYPE type;

   /* convert JNI pointer into C pointer */
   param = (SCIP_PARAM*) (size_t) jparam;
   assert(param != NULL);

   type = SCIPparamGetType(param);

   /* check that the status is one which is captured in the JNI interface */
   assert(type == JNIPACKAGENAME(JniScipParamtype_SCIP_PARAMTYPE_BOOL)
      || type == JNIPACKAGENAME(JniScipParamtype_SCIP_PARAMTYPE_INT)
      || type == JNIPACKAGENAME(JniScipParamtype_SCIP_PARAMTYPE_LONGINT)
      || type == JNIPACKAGENAME(JniScipParamtype_SCIP_PARAMTYPE_REAL)
      || type == JNIPACKAGENAME(JniScipParamtype_SCIP_PARAMTYPE_CHAR)
      || type == JNIPACKAGENAME(JniScipParamtype_SCIP_PARAMTYPE_STRING));

   return (jint) type;
}

/** returns name of parameter */
JNIEXPORT
jstring JNISCIPPARAMSET(paramGetName)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	const char* name;
	jstring jname;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get parameter name */
	name=SCIPparamGetName(param);

	/* convert char* into jstring */
	jname=(*env)->NewStringUTF(env, name);

	return jname;
}

/** returns description of parameter */
JNIEXPORT
jstring JNISCIPPARAMSET(paramGetDesc)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	const char* descr;
	jstring jdescr;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get parameter description */
	descr=SCIPparamGetDesc(param);

	/* convert char* into jstring */
	jdescr=(*env)->NewStringUTF(env, descr);

	return jdescr;
}

/** returns locally defined parameter specific data */
JNIEXPORT
jlong JNISCIPPARAMSET(paramGetData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	return (jlong) (size_t) SCIPparamGetData(param);
}

/** returns whether parameter is advanced */
JNIEXPORT
jboolean JNISCIPPARAMSET(paramIsAdvanced)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Bool isadvanced;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	isadvanced=SCIPparamIsAdvanced(param);

	return (jboolean) isadvanced;
}

/** returns whether parameter is fixed */
JNIEXPORT
jboolean JNISCIPPARAMSET(paramIsFixed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Bool isfixed;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	isfixed=SCIPparamIsFixed(param);

	return (jboolean) isfixed;
}

/** sets fixing status of given parameter */
JNIEXPORT
void JNISCIPPARAMSET(paramSetFixed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam,             /**< parameter */
   jboolean              jfixed              /**< new fixing status of the parameter */
   )
{
	SCIP_PARAM* param;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	SCIPparamSetFixed(param, (SCIP_Bool)jfixed);
}

/** returns value of SCIP_Bool parameter */
JNIEXPORT
jboolean JNISCIPPARAMSET(paramGetBool)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Bool value;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get boolean parameter value */
	value=SCIPparamGetBool(param);

	return (jboolean) value;
}

/** returns default value of SCIP_Bool parameter */
JNIEXPORT
jboolean JNISCIPPARAMSET(paramGetBoolDefault)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Bool defvalue;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get boolean parameter default value */
	defvalue=SCIPparamGetBoolDefault(param);

	return (jboolean) defvalue;
}

/** returns value of int parameter */
JNIEXPORT
jint JNISCIPPARAMSET(paramGetInt)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	int value;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get int parameter value */
	value=SCIPparamGetInt(param);

	return (jint) value;
}

/** returns minimal value of int parameter */
JNIEXPORT
jint JNISCIPPARAMSET(paramGetIntMin)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	int minvalue;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get int parameter minimal value */
	minvalue=SCIPparamGetIntMin(param);

	return (jint) minvalue;
}

/** returns maximal value of int parameter */
JNIEXPORT
jint JNISCIPPARAMSET(paramGetIntMax)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	int maxvalue;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get int parameter maximal value */
	maxvalue=SCIPparamGetIntMax(param);

	return (jint) maxvalue;
}

/** returns default value of int parameter */
JNIEXPORT
jint JNISCIPPARAMSET(paramGetIntDefault)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	int defvalue;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get int parameter default value */
	defvalue=SCIPparamGetIntDefault(param);

	return (jint) defvalue;
}

/** returns value of SCIP_Longint parameter */
JNIEXPORT
jlong JNISCIPPARAMSET(paramGetLongint)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Longint value;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get long int parameter value */
	value=SCIPparamGetLongint(param);

	return (jlong) value;
}

/** returns minimal value of longint parameter */
JNIEXPORT
jlong JNISCIPPARAMSET(paramGetLongintMin)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Longint minvalue;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get long int parameter minimal value */
	minvalue=SCIPparamGetLongintMin(param);

	return (jlong) minvalue;
}

/** returns maximal value of longint parameter */
JNIEXPORT
jlong JNISCIPPARAMSET(paramGetLongintMax)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Longint maxvalue;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get long int parameter maximal value */
	maxvalue=SCIPparamGetLongintMax(param);

	return (jlong) maxvalue;
}

/** returns default value of SCIP_Longint parameter */
JNIEXPORT
jlong JNISCIPPARAMSET(paramGetLongintDefault)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Longint defvalue;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get long int parameter default value */
	defvalue=SCIPparamGetLongintDefault(param);

	return (jlong) defvalue;
}

/** returns value of SCIP_Real parameter */
JNIEXPORT
jdouble JNISCIPPARAMSET(paramGetReal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Real value;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get real parameter value */
	value=SCIPparamGetReal(param);

	return (jdouble) value;
}

/** returns minimal value of real parameter */
JNIEXPORT
jdouble JNISCIPPARAMSET(paramGetRealMin)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Real minvalue;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get real parameter minimal value */
	minvalue=SCIPparamGetRealMin(param);

	return (jdouble) minvalue;
}

/** returns maximal value of real parameter */
JNIEXPORT
jdouble JNISCIPPARAMSET(paramGetRealMax)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Real maxvalue;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get real parameter maximal value */
	maxvalue=SCIPparamGetRealMax(param);

	return (jdouble) maxvalue;
}

/** returns default value of SCIP_Real parameter */
JNIEXPORT
jdouble JNISCIPPARAMSET(paramGetRealDefault)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Real defvalue;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get real parameter default value */
	defvalue=SCIPparamGetRealDefault(param);

	return (jdouble) defvalue;
}

/** returns value of char parameter */
JNIEXPORT
jchar JNISCIPPARAMSET(paramGetChar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	char value;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get char parameter value */
	value=SCIPparamGetChar(param);

	return (jchar) value;
}

/** returns default value of char parameter */
JNIEXPORT
jchar JNISCIPPARAMSET(paramGetCharDefault)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	char defvalue;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get char parameter default value */
	defvalue=SCIPparamGetCharDefault(param);

	return (jchar) defvalue;
}

/** returns value of string parameter */
JNIEXPORT
jstring JNISCIPPARAMSET(paramGetString)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	const char* value;
	jstring jval;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get string parameter value */
	value=SCIPparamGetString(param);

	/* convert char* into jstring */
	jval=(*env)->NewStringUTF(env, value);

	return jval;
}

/** returns default value of String parameter */
JNIEXPORT
jstring JNISCIPPARAMSET(paramGetStringDefault)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	const char* defvalue;
	jstring jdefvalue;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	/* get string parameter default value */
	defvalue=SCIPparamGetStringDefault(param);

	/* convert char* into jstring */
	jdefvalue=(*env)->NewStringUTF(env, defvalue);

	return jdefvalue;
}

/** returns whether the parameter is on its default setting */
JNIEXPORT
jboolean JNISCIPPARAMSET(paramIsDefault)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jparam              /**< parameter */
   )
{
	SCIP_PARAM* param;
	SCIP_Bool isdefault;

	/* convert JNI pointer into C pointer */
	param = (SCIP_PARAM*) (size_t) jparam;
	assert(param != NULL);

	isdefault=SCIPparamIsDefault(param);

	return (jboolean) isdefault;
}
