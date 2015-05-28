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

/**@file   JniScipVar.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP variable callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipVar.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"

#include <string.h>

/*
 * variable methods
 */

/**@name Variable Methods */
/**@{ */



/** gets number of locks for rounding down */
JNIEXPORT
jint JNISCIPVAR(varGetNLocksDown)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   int num;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   num = SCIPvarGetNLocksDown(var);

   return (jint)num;
}

/** gets number of locks for rounding up */
JNIEXPORT
jint JNISCIPVAR(varGetNLocksUp)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   int num;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   num  = SCIPvarGetNLocksUp(var);

   return (jint)num;
}

/** is it possible, to round variable down and stay feasible? */
JNIEXPORT
jboolean JNISCIPVAR(varMayRoundDown)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool feas;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   feas = SCIPvarMayRoundDown(var);

   return (jboolean)feas;
}

/** is it possible, to round variable up and stay feasible? */
JNIEXPORT
jboolean JNISCIPVAR(varMayRoundUp)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool feas;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   feas = SCIPvarMayRoundUp(var);

   return (jboolean)feas;
}

/** compares the index of two variables, only active or negated variables are allowed, if a variable
 *  is negated then the index of the corresponding active variable is taken, returns -1 if first is
 *  smaller than, and +1 if first is greater than second variable index; returns 0 if both indices
 *  are equal, which means both variables are equal
 */
JNIEXPORT
jint JNISCIPVAR(varCompareActiveAndNegated)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar1,              /**< first problem variable */
   jlong                 jvar2               /**< second problem variable */
   )
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int ind;

   /* convert JNI pointer into C pointer */
   var1 = (SCIP_VAR*) (size_t) jvar1;
   assert(var1 != NULL);

   /* convert JNI pointer into C pointer */
   var2 = (SCIP_VAR*) (size_t) jvar2;
   assert(var2 != NULL);

   ind = SCIPvarCompareActiveAndNegated(var1, var2);

   return (jint)ind;
}

/** compares the index of two variables, returns -1 if first is smaller than, and +1 if first is greater than second
 *  variable index; returns 0 if both indices are equal, which means both variables are equal
 */
JNIEXPORT
jint JNISCIPVAR(varCompare)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar1,              /**< first problem variable */
   jlong                 jvar2               /**< second problem variable */
   )
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int ind;

   /* convert JNI pointer into C pointer */
   var1 = (SCIP_VAR*) (size_t) jvar1;
   assert(var1 != NULL);

   /* convert JNI pointer into C pointer */
   var2 = (SCIP_VAR*) (size_t) jvar2;
   assert(var2 != NULL);

   ind = SCIPvarCompare(var1, var2);

   return (jint)ind;
}

/** gets corresponding active, fixed, or multi-aggregated problem variables of given variables,
 *  @note the content of the given array will/might change
 */
JNIEXPORT
void JNISCIPVAR(varsGetProbvar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlongArray            jvars,              /**< JNI problem variable */
   jint                  nvars               /**< number of variables */
   )
{
   SCIPerrorMessage("method varsGetProbvar is not implemented yet\n");
   JNISCIP_CALL( SCIP_ERROR );
}

/** gets corresponding active, fixed, or multi-aggregated problem variable of a variable */
JNIEXPORT
jlong JNISCIPVAR(varGetProbvar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_VAR* probvar;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   probvar = SCIPvarGetProbvar(var);

   return (jlong)(size_t)probvar;
}

/** returns whether the given variable is the direct counterpart of an original problem variable */
JNIEXPORT
jboolean JNISCIPVAR(varIsTransformedOrigvar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool counterpart;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   counterpart = SCIPvarIsTransformedOrigvar(var);

   return (jboolean)counterpart;
}

/** returns the number of times, a bound of the variable was changed in given direction due to branching */
JNIEXPORT
jlong JNISCIPVAR(varGetNBranchings)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_Longint num;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   num = SCIPvarGetNBranchings(var, (SCIP_BRANCHDIR) jdir);

   return (jlong)(size_t)num;
}

/** returns the number of times, a bound of the variable was changed in given direction due to branching
 *  in the current run
 */
JNIEXPORT
jlong JNISCIPVAR(varGetNBranchingsCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_Longint num;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   num = SCIPvarGetNBranchingsCurrentRun(var, (SCIP_BRANCHDIR) jdir);

   return (jlong)(size_t)num;
}

/** returns the number of inferences branching on this variable in given direction triggered */
JNIEXPORT
jdouble JNISCIPVAR(varGetInferenceSum)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_Real is;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   is = SCIPvarGetInferenceSum(var, (SCIP_BRANCHDIR) jdir);

   return (jdouble) is;
}

/** returns the number of inferences branching on this variable in given direction triggered
 *  in the current run
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetInferenceSumCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_Real is;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   is = SCIPvarGetInferenceSumCurrentRun(var, (SCIP_BRANCHDIR) jdir);

   return (jdouble) is;
}

/** returns the number of cutoffs branching on this variable in given direction produced */
JNIEXPORT
jdouble JNISCIPVAR(varGetCutoffSum)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   JNISCIP_ENUM         jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_Real cs;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   cs = SCIPvarGetCutoffSum(var, (SCIP_BRANCHDIR) jdir);

   return (jdouble) cs;
}

/** returns the number of cutoffs branching on this variable in given direction produced in the current run */
JNIEXPORT
jdouble JNISCIPVAR(varGetCutoffSumCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_Real cs;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   cs = SCIPvarGetCutoffSumCurrentRun(var, (SCIP_BRANCHDIR) jdir);

   return (jdouble) cs;
}

/** returns the average depth of bound changes in given direction due to branching on the variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetAvgBranchdepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_Real depth;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   depth = SCIPvarGetAvgBranchdepth(var, (SCIP_BRANCHDIR) jdir);

   return (jdouble)depth;
}

/** returns the average depth of bound changes in given direction due to branching on the variable
 *  in the current run
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetAvgBranchdepthCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_Real depth;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   depth = SCIPvarGetAvgBranchdepthCurrentRun(var, (SCIP_BRANCHDIR) jdir);

   return (jdouble)depth;
}

/** returns whether there is an implication x == varfixing -> y <= b or y >= b in the implication graph;
 *  implications that are represented as cliques in the clique table are not regarded (use SCIPvarsHaveCommonClique());
 *  both variables must be active, variable x must be binary
 */
JNIEXPORT
jboolean JNISCIPVAR(varHasImplic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jvarfixing,         /**< FALSE if y should be searched in implications for x == 0, TRUE for x == 1 */
   jlong                 jimplvar,           /**< variable y to search for */
   JNISCIP_ENUM          jimpltype           /**< type of implication y <=/>= b to search for */
   )
{
   SCIP_VAR* var;
   SCIP_VAR* implvar;
   SCIP_Bool implic;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   /* convert JNI pointer into C pointer */
   implvar = (SCIP_VAR*) (size_t) jimplvar;
   assert(implvar != NULL);

   implic = SCIPvarHasImplic(var, (SCIP_Bool) jvarfixing, implvar, (SCIP_BOUNDTYPE) jimpltype);

   return (jboolean) implic;
}

/** returns whether there is an implication x == varfixing -> y == implvarfixing in the implication graph;
 *  implications that are represented as cliques in the clique table are not regarded (use SCIPvarsHaveCommonClique());
 *  both variables must be active binary variables
 */
JNIEXPORT
jboolean JNISCIPVAR(varHasBinaryImplic)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jvarfixing,         /**< FALSE if y should be searched in implications for x == 0, TRUE for x == 1 */
   jlong                 jimplvar,           /**< variable y to search for */
   jboolean              jimplvarfixing      /**< value of the implied variable to search for */
   )
{
   SCIP_VAR* var;
   SCIP_VAR* implvar;
   SCIP_Bool implic;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   /* convert JNI pointer into C pointer */
   implvar = (SCIP_VAR*) (size_t) jimplvar;
   assert(implvar != NULL);

   implic = SCIPvarHasBinaryImplic(var, (SCIP_Bool) jvarfixing, implvar, (SCIP_Bool) jimplvarfixing);

   return (jboolean) implic;
}

/** returns whether there is a clique that contains both given variable/value pairs;
 *  the variables must be active binary variables;
 *  if regardimplics is FALSE, only the cliques in the clique table are looked at;
 *  if regardimplics is TRUE, both the cliques and the implications of the implication graph are regarded
 */
JNIEXPORT
jboolean JNISCIPVAR(varsHaveCommonClique)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar1,              /**< firstproblem variable */
   jboolean              jvalue1,            /**< value of first variable */
   jlong                 jvar2,              /**< second variable */
   jboolean              jvalue2,            /**< value of second variable */
   jboolean              jregardimplics      /**< should the implication graph also be searched for a clique? */
   )
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_Bool clique;

   /* convert JNI pointer into C pointer */
   var1 = (SCIP_VAR*) (size_t) jvar1;
   assert(var1 != NULL);

   var2 = (SCIP_VAR*) (size_t) jvar2;
   assert(var2 != NULL);

   clique = SCIPvarsHaveCommonClique(var1, (SCIP_Bool) jvalue1, var2, (SCIP_Bool) jvalue2, (SCIP_Bool) jregardimplics);

   return (jboolean)clique;
}



/** sets the initial flag of a variable; only possible for original or loose variables */
JNIEXPORT
void JNISCIPVAR(varSetInitial)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jinitial
   )
{
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*)(size_t)jvar;
   assert(var != NULL);

   SCIPvarSetInitial(var, (SCIP_Bool)jinitial);
}


/** sets the removable flag of a variable; only possible for original or loose variables */
JNIEXPORT
void JNISCIPVAR(varSetRemovable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jremovable
   )
{
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*)(size_t)jvar;
   assert(var != NULL);

   SCIPvarSetRemovable(var, (SCIP_Bool)jremovable);
}

/** returns name of the variable **/
JNIEXPORT
jstring JNISCIPVAR(varGetName)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   const char* name;
   jstring jname;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*)(size_t)jvar;
   assert(var != NULL);

   /* get variable name */
   name = SCIPvarGetName(var);

   /* convert char* into jstring */
   jname = (*env)->NewStringUTF(env, name);

   return jname;
}

/** gets number of times, the variable is currently captured */
JNIEXPORT
jint JNISCIPVAR(varGetNUses)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*)(size_t)jvar;
   assert(var != NULL);

   return (jint)SCIPvarGetNUses(var);
}

/** returns the user data of the variable */
JNIEXPORT
jlong JNISCIPVAR(varGetData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*)(size_t)jvar;
   assert(var != NULL);

   return (jlong) (size_t) SCIPvarGetData(var);
}

void JNISCIPVAR(varSetData)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jlong                 jvardata            /**< user variable data */
   )
{
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*)(size_t)jvar;
   assert(var != NULL);

   SCIPvarSetData(var, (SCIP_VARDATA*)(size_t)jvardata);
}

/** gets status of variable */
JNIEXPORT
jint JNISCIPVAR(varGetStatus)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_VARSTATUS varstatus;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   varstatus = SCIPvarGetStatus(var);

   return (jint) varstatus;
}

/** returns whether the variable belongs to the original problem */
JNIEXPORT
jboolean JNISCIPVAR(varIsOriginal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool varorigprob;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   varorigprob = SCIPvarIsOriginal(var);

   return (jboolean) varorigprob;
}

/** returns whether the variable belongs to the transformed problem */
JNIEXPORT
jboolean JNISCIPVAR(varIsTransformed)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool transformed;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   transformed = SCIPvarIsTransformed(var);

   return (jboolean) transformed;
}

/** returns whether the variable was created by negation of a different variable */
JNIEXPORT
jboolean JNISCIPVAR(varIsNegated)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool negated;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   negated = SCIPvarIsNegated(var);

   return (jboolean) negated;
}

/** gets type of variable */
JNIEXPORT
jint JNISCIPVAR(varGetType)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_VARTYPE vartype;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   vartype = SCIPvarGetType(var);

   return (jint) vartype;
}

/** returns TRUE if the variable is of binary type; this is the case if:
 *  (1) variable type is binary
 *  (2) variable type is integer or implicit integer and
 *      (i)  the lazy lower bound or the global lower bound is greater or equal to zero
 *      (ii) the lazy upper bound or the global upper bound is less tor equal to one
 */
JNIEXPORT
jboolean JNISCIPVAR(varIsBinary)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool binary;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   binary= SCIPvarIsBinary(var);

   return (jboolean) binary;
}

/** returns whether variable is of integral type (binary, integer, or implicit integer) */
JNIEXPORT
jboolean JNISCIPVAR(varIsIntegral)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool integral;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   integral = SCIPvarIsIntegral(var);

   return (jboolean) integral;
}

/** returns whether variable's column should be present in the initial root LP */
JNIEXPORT
jboolean JNISCIPVAR(varIsInitial)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool initial;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   initial = SCIPvarIsInitial(var);

   return (jboolean) initial;
}

/** returns whether variable's column is removable from the LP (due to aging or cleanup) */
JNIEXPORT
jboolean JNISCIPVAR(varIsRemovable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool removable;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   removable= SCIPvarIsRemovable(var);

   return (jboolean) removable;
}

/** returns whether the variable was deleted from the problem */
JNIEXPORT
jboolean JNISCIPVAR(varIsDeleted)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool deleted;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   deleted = SCIPvarIsDeleted(var);

   return (jboolean) deleted;
}

/** marks the variable to be deletable, i.e., it may be deleted completely from the problem;
 *  method can only be called before the variable is added to the problem by SCIPaddVar() or SCIPaddPricedVar()
 */
JNIEXPORT
void JNISCIPVAR(varMarkDeletable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   SCIPvarMarkDeletable(var);
}

/** marks the variable to be deletable, i.e., it may be deleted completely from the problem;
 *  method can only be called before the variable is added to the problem by SCIPaddVar() or SCIPaddPricedVar()
 */
JNIEXPORT
void JNISCIPVAR(varMarkNotDeletable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   SCIPvarMarkNotDeletable(var);
}

/** returns whether variable is allowed to be deleted completely from the problem */
JNIEXPORT
jboolean JNISCIPVAR(varIsDeletable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool deletable;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   deletable = SCIPvarIsDeletable(var);

   return (jboolean) deletable;
}

/** returns whether variable is an active (neither fixed nor aggregated) variable */
JNIEXPORT
jboolean JNISCIPVAR(varIsActive)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool active;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   active = SCIPvarIsActive(var);

   return (jboolean) active;
}

/** gets unique index of variable */
JNIEXPORT
jint JNISCIPVAR(varGetIndex)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   int ind;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   ind = SCIPvarGetIndex(var);

   return (jint) ind;
}

/** gets position of variable in problem, or -1 if variable is not active */
JNIEXPORT
jint JNISCIPVAR(varGetProbindex)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   int probindex;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   probindex = SCIPvarGetProbindex(var);

   return (jint) probindex;
}

/** gets transformed variable of ORIGINAL variable */
JNIEXPORT
jlong JNISCIPVAR(varGetTransVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_VAR* transvar;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   transvar = SCIPvarGetTransVar(var);

   return (jlong) (size_t) transvar;
}

/** gets column of COLUMN variable */
JNIEXPORT
jlong JNISCIPVAR(varGetCol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_COL* col;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   col = SCIPvarGetCol(var);

   return (jlong) (size_t) col;
}

/** returns whether the variable is a COLUMN variable that is member of the current LP */
JNIEXPORT
jboolean JNISCIPVAR(varIsInLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Bool lp;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   lp = SCIPvarIsInLP(var);

   return (jboolean) lp;
}

/** gets aggregation variable y of an aggregated variable x = a*y + c */
JNIEXPORT
jlong JNISCIPVAR(varGetAggrVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_VAR* aggrvar;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   aggrvar = SCIPvarGetAggrVar(var);

   return (jlong) (size_t) aggrvar;
}

/** gets aggregation scalar a of an aggregated variable x = a*y + c */
JNIEXPORT
jdouble JNISCIPVAR(varGetAggrScalar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real aggrscalar;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   aggrscalar = SCIPvarGetAggrScalar(var);

   return (jdouble) aggrscalar;
}

/** gets aggregation constant c of an aggregated variable x = a*y + c */
JNIEXPORT
jdouble JNISCIPVAR(varGetAggrConstant)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real aggrconstant;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   aggrconstant = SCIPvarGetAggrConstant(var);

   return (jdouble) aggrconstant;
}


/** gets number n of aggregation variables of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
JNIEXPORT
jint JNISCIPVAR(varGetMultaggrNVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   int nvars;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   nvars = SCIPvarGetMultaggrNVars(var);

   return (jint) nvars;
}

/** gets vector of aggregation variables y of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
JNIEXPORT
jlongArray JNISCIPVAR(varGetMultaggrVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_VAR** vars;
   int nvars;
   jlongArray jvars;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   vars = SCIPvarGetMultaggrVars(var);
   nvars = SCIPvarGetMultaggrNVars(var);

   jvars = (*env)->NewLongArray(env, nvars);

   (*env)->SetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);

   return jvars;
}

/** gets vector of aggregation scalars a of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
JNIEXPORT
jdoubleArray JNISCIPVAR(varGetMultaggrScalars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real* vars;
   int nvars;
   jlongArray jvars;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   vars = SCIPvarGetMultaggrScalars(var);
   nvars = SCIPvarGetMultaggrNVars(var);

   /* create jlongArray */
   jvars = (*env)->NewDoubleArray(env, nvars);

   /* fill long array with SCIP variable pointers */
   (*env)->SetDoubleArrayRegion(env, jvars, 0, nvars, (jdouble*)vars);

   return jvars;
}

/** gets aggregation constant c of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
JNIEXPORT
jdouble JNISCIPVAR(varGetMultaggrConstant)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real aggrconst;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   aggrconst = SCIPvarGetMultaggrConstant(var);

   return (jdouble) aggrconst;
}

/** gets the negation of the given variable; may return NULL, if no negation is existing yet */
JNIEXPORT
jlong JNISCIPVAR(varGetNegatedVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_VAR* negvar;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   negvar = SCIPvarGetNegatedVar(var);

   return (jlong) (size_t) negvar;
}


/** gets the negation variable x of a negated variable x' = offset - x */
JNIEXPORT
jlong JNISCIPVAR(varGetNegationVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_VAR* negvar;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   negvar = SCIPvarGetNegationVar(var);

   return (jlong) (size_t) negvar;
}

/** gets the negation offset of a negated variable x' = offset - x */
JNIEXPORT
jdouble JNISCIPVAR(varGetNegationConstant)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real negconst;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   negconst = SCIPvarGetNegationConstant(var);

   return (jdouble) negconst;
}

/** gets objective function value of variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetObj)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real val;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   val = SCIPvarGetObj(var);

   return (jdouble) val;
}



/** gets original lower bound of original problem variable (i.e. the bound set in problem creation) */
JNIEXPORT
jdouble JNISCIPVAR(varGetLbOriginal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real lb;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   lb = SCIPvarGetLbOriginal(var);

   return (jdouble) lb;
}

/** gets original upper bound of original problem variable (i.e. the bound set in problem creation) */
JNIEXPORT
jdouble JNISCIPVAR(varGetUbOriginal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real ub;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   ub = SCIPvarGetUbOriginal(var);

   return (jdouble) ub;
}

/** gets the original hole list of an original variable */
JNIEXPORT
jlong JNISCIPVAR(varGetHolelistOriginal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_HOLELIST* list;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   list = SCIPvarGetHolelistOriginal(var);

   return (jlong)(size_t)list;
}

/** gets global lower bound of variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetLbGlobal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real lb;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   lb = SCIPvarGetLbGlobal(var);

   return (jdouble) lb;
}

/** gets global upper bound of variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetUbGlobal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real ub;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   ub = SCIPvarGetUbGlobal(var);

   return (jdouble) ub;
}

/** gets the global hole list of an active variable */
JNIEXPORT
jlong JNISCIPVAR(varGetHolelistGlobal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_HOLELIST* list;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   list = SCIPvarGetHolelistGlobal(var);

   return (jlong) (size_t) list;
}

/** gets best global bound of variable with respect to the objective function */
JNIEXPORT
jdouble JNISCIPVAR(varGetBestBoundGlobal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real bbg;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bbg = SCIPvarGetWorstBoundGlobal(var);

   return (jdouble) bbg;
}

/** gets worst global bound of variable with respect to the objective function */
JNIEXPORT
jdouble JNISCIPVAR(varGetWorstBoundGlobal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real wbg;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   wbg = SCIPvarGetWorstBoundGlobal(var);

   return (jdouble) wbg;
}

/** gets local lower bound of variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetLbLocal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real lb;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   lb = SCIPvarGetLbLocal(var);

   return (jdouble) lb;
}

/** gets local upper bound of variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetUbLocal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real ub;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   ub = SCIPvarGetUbLocal(var);

   return (jdouble) ub;
}

/** gets the local hole list of an active variable */
JNIEXPORT
jlong JNISCIPVAR(varGetHolelistLocal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_HOLELIST* list;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   list = SCIPvarGetHolelistLocal(var);

   return (jlong) (size_t) list;
}

/** gets best local bound of variable with respect to the objective function */
JNIEXPORT
jdouble JNISCIPVAR(varGetBestBoundLocal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real bound;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bound= SCIPvarGetBestBoundLocal(var);

   return (jdouble) bound;
}

/** gets worst local bound of variable with respect to the objective function */
JNIEXPORT
jdouble JNISCIPVAR(varGetWorstBoundLocal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real bound;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bound= SCIPvarGetWorstBoundLocal(var);

   return (jdouble) bound;
}

/** gets type (lower or upper) of best bound of variable with respect to the objective function */
JNIEXPORT
jint JNISCIPVAR(varGetBestBoundType)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_BOUNDTYPE type;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   type = SCIPvarGetBestBoundType(var);

   return (jint) type;
}

/** gets type (lower or upper) of worst bound of variable with respect to the objective function */
JNIEXPORT
jint JNISCIPVAR(varGetWorstBoundType)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_BOUNDTYPE type ;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   type= SCIPvarGetWorstBoundType(var);

   return (jint) type;
}

/** gets lazy lower bound of variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetLbLazy)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real lb;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   lb = SCIPvarGetLbLazy(var);

   return (jdouble) lb;
}

/** gets lazy upper bound of variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetUbLazy)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real ub;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   ub = SCIPvarGetUbLazy(var);

   return (jdouble) ub;
}

/** gets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetBranchFactor)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real branchfactor;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   branchfactor= SCIPvarGetBranchFactor(var);

   return (jdouble) branchfactor;
}

/** gets the branch priority of the variable; variables with higher priority should always be preferred to variables
 *  with lower priority
 */
JNIEXPORT
jint JNISCIPVAR(varGetBranchPriority)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   int branchpriority;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   branchpriority= SCIPvarGetBranchPriority(var);

   return (jint) branchpriority;
}

/** gets the preferred branch direction of the variable (downwards, upwards, or auto) */
JNIEXPORT
jint JNISCIPVAR(varGetBranchDirection)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_BRANCHDIR branchdir;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   branchdir = SCIPvarGetBranchDirection(var);

   return (jint) branchdir;
}

/** gets number of variable lower bounds x >= b_i*z_i + d_i of given variable x */
JNIEXPORT
jint JNISCIPVAR(varGetNVlbs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   int nvars;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   nvars = SCIPvarGetNVlbs(var);

   return (jint) nvars;
}

/** gets array with bounding variables z_i in variable lower bounds x >= b_i*z_i + d_i of given variable x;
 *  the variable bounds are sorted by increasing variable index of the bounding variable z_i (see SCIPvarGetIndex())
 */
JNIEXPORT
jlongArray JNISCIPVAR(varGetVlbVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_VAR** vars;
   int nvars;
   jlongArray jvars;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   vars = SCIPvarGetVlbVars(var);
   nvars = SCIPvarGetNVlbs(var);


   jvars = (*env)->NewLongArray(env, nvars);


   (*env)->SetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);

   return jvars;
}

/** gets array with bounding coefficients b_i in variable lower bounds x >= b_i*z_i + d_i of given variable x */
JNIEXPORT
jdoubleArray JNISCIPVAR(varGetVlbCoefs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real* coeffs;
   int ncoeffs;
   jdoubleArray jcoeffs;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   coeffs = SCIPvarGetVlbCoefs(var);
   ncoeffs = SCIPvarGetNVlbs(var);

   /* create jlongArray */
   jcoeffs = (*env)->NewDoubleArray(env, ncoeffs);

   /* fill long array with SCIP variable pointers */
   (*env)->SetDoubleArrayRegion(env, jcoeffs, 0, ncoeffs, (jdouble*)coeffs);

   return jcoeffs;
}

/** gets array with bounding constants d_i in variable lower bounds x >= b_i*z_i + d_i of given variable x */
JNIEXPORT
jdoubleArray JNISCIPVAR(varGetVlbConstants)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real* consts;
   int nconsts;
   jdoubleArray jconsts;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   consts = SCIPvarGetVlbConstants(var);
   nconsts = SCIPvarGetNVlbs(var);

   /* create jlongArray */
   jconsts = (*env)->NewDoubleArray(env, nconsts);

   /* fill long array with SCIP variable pointers */
   (*env)->SetDoubleArrayRegion(env, jconsts, 0, nconsts, (jdouble*)consts);

   return jconsts;
}

/** gets number of variable upper bounds x <= b_i*z_i + d_i of given variable x */
JNIEXPORT
jint JNISCIPVAR(varGetNVubs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   int nvars;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   nvars= SCIPvarGetNVubs(var);

   return (jint) nvars;
}

/** gets array with bounding variables z_i in variable upper bounds x <= b_i*z_i + d_i of given variable x;
 *  the variable bounds are sorted by increasing variable index of the bounding variable z_i (see SCIPvarGetIndex())
 */
JNIEXPORT
jlongArray JNISCIPVAR(varGetVubVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_VAR** vars;
   int nvars;
   jlongArray jvars;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   vars = SCIPvarGetVubVars(var);
   nvars = SCIPvarGetNVubs(var);


   jvars = (*env)->NewLongArray(env, nvars);


   (*env)->SetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);

   return jvars;
}

/** gets array with bounding coefficients b_i in variable upper bounds x <= b_i*z_i + d_i of given variable x */
JNIEXPORT
jdoubleArray JNISCIPVAR(varGetVubCoefs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real* coeffs;
   int ncoeffs;
   jdoubleArray jcoeffs;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   coeffs = SCIPvarGetVubCoefs(var);
   ncoeffs = SCIPvarGetNVubs(var);

   /* create jlongArray */
   jcoeffs = (*env)->NewDoubleArray(env, ncoeffs);

   /* fill long array with SCIP variable pointers */
   (*env)->SetDoubleArrayRegion(env, jcoeffs, 0, ncoeffs, (jdouble*)coeffs);

   return jcoeffs;
}

/** gets array with bounding constants d_i in variable upper bounds x <= b_i*z_i + d_i of given variable x */
JNIEXPORT
jdoubleArray JNISCIPVAR(varGetVubConstants)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real* consts;
   int nconsts;
   jdoubleArray jconsts;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   consts = SCIPvarGetVubConstants(var);
   nconsts = SCIPvarGetNVubs(var);

   /* create jlongArray */
   jconsts = (*env)->NewDoubleArray(env, nconsts);

   /* fill long array with SCIP variable pointers */
   (*env)->SetDoubleArrayRegion(env, jconsts, 0, nconsts, (jdouble*)consts);

   return jconsts;
}

/** gets number of implications  y <= b or y >= b for x == 0 or x == 1 of given active problem variable x,
 *  there are no implications for nonbinary variable x
 */
JNIEXPORT
jint JNISCIPVAR(varGetNImpls)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jvarfixing          /**< FALSE for implications for x == 0, TRUE for x == 1 */
   )
{
   SCIP_VAR* var;
   int num;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   num = SCIPvarGetNImpls(var, (SCIP_Bool) jvarfixing);

   return (jint) num;
}

/** gets array with implication variables y of implications  y <= b or y >= b for x == 0 or x == 1 of given active
 *  problem variable x, there are no implications for nonbinary variable x;
 *  the implications are sorted such that implications with binary implied variables precede the ones with non-binary
 *  implied variables, and as a second criteria, the implied variables are sorted by increasing variable index
 *  (see SCIPvarGetIndex())
 */
JNIEXPORT
jlongArray JNISCIPVAR(varGetImplVars)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jvarfixing          /**< FALSE for implications for x == 0, TRUE for x == 1 */
   )
{
   SCIP_VAR* var;
   SCIP_VAR** vars;
   int nvars;
   jlongArray jvars;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   vars = SCIPvarGetImplVars(var, (SCIP_Bool) jvarfixing);
   nvars = SCIPvarGetNImpls(var, (SCIP_Bool) jvarfixing);


   jvars = (*env)->NewLongArray(env, nvars);


   (*env)->SetLongArrayRegion(env, jvars, 0, nvars, (jlong*)vars);

   return jvars;
}

/** gets array with implication types of implications  y <= b or y >= b for x == 0 or x == 1 of given active problem
 *  variable x (SCIP_BOUNDTYPE_UPPER if y <= b, SCIP_BOUNDTYPE_LOWER if y >= b),
 *  there are no implications for nonbinary variable x
 */
JNIEXPORT
jintArray JNISCIPVAR(varGetImplTypes)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jvarfixing          /**< FALSE for implications for x == 0, TRUE for x == 1 */
   )
{
   SCIP_VAR* var;
   SCIP_BOUNDTYPE* type;
   int ntype;
   jintArray jtype;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   type = SCIPvarGetImplTypes(var, (SCIP_Bool) jvarfixing);
   ntype = SCIPvarGetNImpls(var, (SCIP_Bool) jvarfixing);

   /* create jlongArray */
   jtype = (*env)->NewIntArray(env, ntype);

   /* fill long array with SCIP variable pointers */
   (*env)->SetIntArrayRegion(env, jtype, 0, ntype, (jint*)type);

   return jtype;
}

/** gets array with implication bounds b of implications  y <= b or y >= b for x == 0 or x == 1 of given active problem
 *  variable x, there are no implications for nonbinary variable x
 */
JNIEXPORT
jdoubleArray JNISCIPVAR(varGetImplBounds)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jvarfixing          /**< FALSE for implications for x == 0, TRUE for x == 1 */
   )
{
   SCIP_VAR* var;
   SCIP_Real* bounds;
   int nbounds;
   jdoubleArray jbounds;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bounds = SCIPvarGetImplBounds(var, (SCIP_Bool) jvarfixing);
   nbounds = SCIPvarGetNImpls(var, (SCIP_Bool) jvarfixing);

   /* create jdoublerray */
   jbounds = (*env)->NewDoubleArray(env, nbounds);

   /* fill double array with SCIP variable pointers */
   (*env)->SetDoubleArrayRegion(env, jbounds, 0, nbounds, (jdouble*)bounds);

   return jbounds;
}

/** gets array with unique ids of implications  y <= b or y >= b for x == 0 or x == 1 of given active problem variable x,
 *  there are no implications for nonbinary variable x
 */
JNIEXPORT
jintArray JNISCIPVAR(varGetImplIds)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jvarfixing          /**< FALSE for implications for x == 0, TRUE for x == 1 */
   )
{
   SCIP_VAR* var;
   int* ids;
   int nids;
   jintArray jids;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   ids = SCIPvarGetImplIds(var, (SCIP_Bool) jvarfixing);
   nids = SCIPvarGetNImpls(var, (SCIP_Bool) jvarfixing);

   /* create jintArray */
   jids = (*env)->NewIntArray(env, nids);

   /* fill integer array with SCIP variable pointers */
   (*env)->SetIntArrayRegion(env, jids, 0, nids, (jint*)ids);

   return jids;
}

/** gets number of cliques, the active variable is contained in */
JNIEXPORT
jint JNISCIPVAR(varGetNCliques)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jvarfixing          /**< FALSE for implications for x == 0, TRUE for x == 1 */
   )
{
   SCIP_VAR* var;
   int num;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   num = SCIPvarGetNCliques(var, (SCIP_Bool) jvarfixing);

   return (jint) num;
}

/** gets array of cliques, the active variable is contained in */
JNIEXPORT
jlongArray JNISCIPVAR(varGetCliques)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jvarfixing          /**< FALSE for implications for x == 0, TRUE for x == 1 */
   )
{
   SCIP_VAR* var;
   SCIP_CLIQUE** cliques;
   int ncliques;
   jlongArray jcliques;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   cliques = SCIPvarGetCliques(var, (SCIP_Bool) jvarfixing);
   ncliques = SCIPvarGetNCliques(var, (SCIP_Bool) jvarfixing);


   jcliques = (*env)->NewLongArray(env, ncliques);


   (*env)->SetLongArrayRegion(env, jcliques, 0, ncliques, (jlong*)cliques);

   return jcliques;
}

/** gets primal LP solution value of variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetLPSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real lpsol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   lpsol = SCIPvarGetLPSol(var);

   return (jdouble) lpsol;
}

/** gets primal NLP solution value of variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetNLPSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real nlpsol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   nlpsol = SCIPvarGetNLPSol(var);

   return (jdouble) nlpsol;
}

/** return lower bound change info at requested position */
JNIEXPORT
jlong JNISCIPVAR(varGetBdchgInfoLb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jint                  jpos                /**< requested position */
   )
{
   SCIP_VAR* var;
   SCIP_BDCHGINFO* lbchange;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   lbchange = SCIPvarGetBdchgInfoLb(var, (int) jpos);

   return (jlong) (size_t) lbchange;
}

/** gets the number of lower bound change info array */
JNIEXPORT
jint JNISCIPVAR(varGetNBdchgInfosLb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   int num;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   num = SCIPvarGetNBdchgInfosLb(var);

   return (jint) num;
}

/** return upper bound change info at requested position */
JNIEXPORT
jlong JNISCIPVAR(varGetBdchgInfoUb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jint                  jpos                /**< requested position */
   )
{
   SCIP_VAR* var;
   SCIP_BDCHGINFO* ubchange;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   ubchange = SCIPvarGetBdchgInfoUb(var, (int) jpos);

   return (jlong) (size_t) ubchange;
}

/** gets the number upper bound change info array */
JNIEXPORT
jint JNISCIPVAR(varGetNBdchgInfosUb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   int num;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   num = SCIPvarGetNBdchgInfosUb(var);

   return (jint) num;
}


/** gets primal LP solution value of variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetLPSol_1rec)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real sol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   sol = SCIPvarGetLPSol_rec(var);

   return (jdouble) sol;
}

/** gets primal NLP solution value of variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetNLPSol_1rec)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real sol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   sol = SCIPvarGetNLPSol_rec(var);

   return (jdouble) sol;
}

/** gets pseudo solution value of variable at current node */
JNIEXPORT
jdouble JNISCIPVAR(varGetPseudoSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real sol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   sol = SCIPvarGetPseudoSol(var);

   return (jdouble) sol;
}

/** gets current LP or pseudo solution value of variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jgetlpval           /**< should the LP solution value be returned? */
   )
{
   SCIP_VAR* var;
   SCIP_Real sol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   sol = SCIPvarGetSol(var, (SCIP_Bool) jgetlpval);

   return (jdouble) sol;
}

/** returns the solution of the variable in the last root node's relaxation, if the root relaxation is not yet
 *  completely solved, zero is returned
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetRootSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real rootsol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;

   rootsol = SCIPvarGetRootSol(var);

   return (jdouble) rootsol;
}

/** returns the best solution (w.r.t. root reduced cost propagation) of the variable in the root node's relaxation, if
 *  the root relaxation is not yet completely solved, zero is returned
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetBestRootSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real rootsol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   rootsol = SCIPvarGetBestRootSol(var);

   return (jdouble) rootsol;
}

/** returns the reduced costs of the variable in the root node's relaxation, if the root relaxation is not yet completely
 *  solved, or the variable was no column of the root LP, SCIP_INVALID is returned
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetBestRootRedcost)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real redcosts;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   redcosts = SCIPvarGetBestRootRedcost(var);

   return (jdouble) redcosts;
}

/** returns the best objective value (w.r.t. root reduced cost propagation) of the root LP which belongs the root
 *  reduced cost which is accessible via SCIPvarGetRootRedcost() or the variable was no column of the root LP,
 *  SCIP_INVALID is returned
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetBestRootLPObjval)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real objval;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   objval = SCIPvarGetBestRootLPObjval(var);

   return (jdouble) objval;
}

/** set the given solution as the best root solution w.r.t. root reduced cost propagation in the variables */
JNIEXPORT
void JNISCIPVAR(varSetBestRootSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jdouble               jrootsol,           /**< root solution value */
   jdouble               jrootredcost,       /**< root reduced cost */
   jdouble               jrootlpobjval       /**< objective value of the root LP */
   )
{
   SCIP_VAR* var;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   SCIPvarSetBestRootSol(var, (SCIP_Real)jrootsol, (SCIP_Real)jrootredcost, (SCIP_Real)jrootlpobjval);
}

/** returns a weighted average solution value of the variable in all feasible primal solutions found so far */
JNIEXPORT
jdouble JNISCIPVAR(varGetAvgSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real sol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   sol = SCIPvarGetAvgSol(var);

   return (jdouble) sol;
}

/** returns the bound change information for the last lower bound change on given active problem variable before or
 *  after the bound change with the given index was applied;
 *  returns NULL, if no change to the lower bound was applied up to this point of time
 */
JNIEXPORT
jlong JNISCIPVAR(varGetLbchgInfo)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jlong                 jbdchgidx,          /**< bound change index representing time on path to current node */
   jboolean              jafter              /**< should the bound change with given index be included? */
   )
{
   SCIP_VAR* var;
   SCIP_BDCHGINFO* info;
   SCIP_BDCHGIDX* bdchgidx;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bdchgidx = (SCIP_BDCHGIDX*) (size_t) jbdchgidx;
   assert(bdchgidx != NULL);

   info = SCIPvarGetLbchgInfo(var, bdchgidx, (SCIP_Bool) jafter);

   return (jlong)(size_t)info;
}

/** returns the bound change information for the last upper bound change on given active problem variable before or
 *  after the bound change with the given index was applied;
 *  returns NULL, if no change to the upper bound was applied up to this point of time
 */
JNIEXPORT
jlong JNISCIPVAR(varGetUbchgInfo)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jlong                 jbdchgidx,          /**< bound change index representing time on path to current node */
   jboolean              jafter              /**< should the bound change with given index be included? */
   )
{
   SCIP_VAR* var;
   SCIP_BDCHGINFO* info;
   SCIP_BDCHGIDX* bdchgidx;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bdchgidx = (SCIP_BDCHGIDX*) (size_t) jbdchgidx;
   assert(bdchgidx != NULL);

   info = SCIPvarGetUbchgInfo(var, bdchgidx, (SCIP_Bool) jafter);

   return (jlong)(size_t)info;
}

/** returns the bound change information for the last lower or upper bound change on given active problem variable
 *  before or after the bound change with the given index was applied;
 *  returns NULL, if no change to the lower/upper bound was applied up to this point of time
 */
JNIEXPORT
jlong JNISCIPVAR(varGetBdchgInfo)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   JNISCIP_ENUM          jboundtype,         /**< type of bound: lower or upper bound */
   jlong                 jbdchgidx,          /**< bound change index representing time on path to current node */
   jboolean              jafter              /**< should the bound change with given index be included? */
   )
{
   SCIP_VAR* var;
   SCIP_BDCHGINFO* info;
   SCIP_BDCHGIDX* bdchgidx;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bdchgidx = (SCIP_BDCHGIDX*) (size_t) jbdchgidx;
   assert(bdchgidx != NULL);

   info = SCIPvarGetBdchgInfo(var, (SCIP_BOUNDTYPE) jboundtype, bdchgidx, (SCIP_Bool) jafter);

   return (jlong)(size_t)info;
}

/** returns lower bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetLbAtIndex)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jlong                 jbdchgidx,          /**< bound change index representing time on path to current node */
   jboolean              jafter              /**< should the bound change with given index be included? */
   )
{
   SCIP_VAR* var;
   SCIP_BDCHGIDX* bdchgidx;
   SCIP_Real lb;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bdchgidx = (SCIP_BDCHGIDX*) (size_t) jbdchgidx;
   assert(bdchgidx != NULL);

   lb = SCIPvarGetLbAtIndex(var, bdchgidx, (SCIP_Bool) jafter);

   return (jdouble) lb;
}

/** returns upper bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetUbAtIndex)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jlong                 jbdchgidx,          /**< bound change index representing time on path to current node */
   jboolean              jafter              /**< should the bound change with given index be included? */
   )
{
   SCIP_VAR* var;
   SCIP_BDCHGIDX* bdchgidx;
   SCIP_Real ub;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bdchgidx = (SCIP_BDCHGIDX*) (size_t) jbdchgidx;
   assert(bdchgidx != NULL);

   ub = SCIPvarGetUbAtIndex(var, bdchgidx, (SCIP_Bool) jafter);

   return (jdouble) ub;
}

/** returns lower or upper bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetBdAtIndex)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   JNISCIP_ENUM          jboundtype,         /**< type of bound: lower or upper bound */
   jlong                 jbdchgidx,          /**< bound change index representing time on path to current node */
   jboolean              jafter              /**< should the bound change with given index be included? */
   )
{
   SCIP_VAR* var;
   SCIP_BDCHGIDX* bdchgidx;
   SCIP_Real bound;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bdchgidx = (SCIP_BDCHGIDX*) (size_t) jbdchgidx;
   assert(bdchgidx != NULL);

   bound = SCIPvarGetBdAtIndex(var, (SCIP_BOUNDTYPE) jboundtype, bdchgidx, (SCIP_Bool) jafter);

   return (jdouble) bound;
}

/** returns whether the binary variable was fixed at the time given by the bound change index */
JNIEXPORT
jboolean JNISCIPVAR(varWasFixedAtIndex)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jlong                 jbdchgidx,          /**< bound change index representing time on path to current node */
   jboolean              jafter              /**< should the bound change with given index be included? */
   )
{
   SCIP_VAR* var;
   SCIP_BDCHGIDX* bdchgidx;
   SCIP_Bool fixed;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bdchgidx = (SCIP_BDCHGIDX*) (size_t) jbdchgidx;
   assert(bdchgidx != NULL);

   fixed = SCIPvarWasFixedAtIndex(var, bdchgidx, (SCIP_Bool) jafter);

   return (jboolean) fixed;
}

/** returns the last bound change index, at which the bounds of the given variable were tightened */
JNIEXPORT
jlong JNISCIPVAR(varGetLastBdchgIndex)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_BDCHGIDX* bdchgidx;


   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   bdchgidx = SCIPvarGetLastBdchgIndex(var);

   return (jlong)(size_t) bdchgidx;
}

/** returns the last depth level, at which the bounds of the given variable were tightened;
 *  returns -2, if the variable's bounds are still the global bounds
 *  returns -1, if the variable was fixed in presolving
 */
JNIEXPORT
jint JNISCIPVAR(varGetLastBdchgDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   int dlevel;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   dlevel = SCIPvarGetLastBdchgDepth(var);

   return (jint) dlevel;
}

/** returns whether the first binary variable was fixed earlier than the second one;
 *  returns FALSE, if the first variable is not fixed, and returns TRUE, if the first variable is fixed, but the
 *  second one is not fixed
 */
JNIEXPORT
jboolean JNISCIPVAR(varWasFixedEarlier)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar1,              /**< first problem variable */
   jlong                 jvar2               /**< second problem variable */
   )
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_Bool fixedearlier;

   /* convert JNI pointer into C pointer */
   var1 = (SCIP_VAR*) (size_t) jvar1;
   assert(var1 != NULL);

   var2 = (SCIP_VAR*) (size_t) jvar2;
   assert(var2 != NULL);

   fixedearlier = SCIPvarWasFixedEarlier(var1, var2);

   return (jboolean) fixedearlier;
}

/** returns whether first bound change index belongs to an earlier applied bound change than second one;
 *  if a bound change index is NULL, the bound change index represents the current time, i.e. the time after the
 *  last bound change was applied to the current node
 */
JNIEXPORT
jboolean JNISCIPVAR(bdchgidxIsEarlier)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchgidx1,         /**< first bound change index, or NULL */
   jlong                 jbdchgidx2          /**< second bound change index, or NULL */
   )
{
   SCIP_BDCHGIDX* bdchgidx1;
   SCIP_BDCHGIDX* bdchgidx2;
   SCIP_Bool earlier;

   /* convert JNI pointer into C pointer */
   bdchgidx1 = (SCIP_BDCHGIDX*) (size_t) jbdchgidx1;
   bdchgidx2 = (SCIP_BDCHGIDX*) (size_t) jbdchgidx2;


   earlier = SCIPbdchgidxIsEarlier(bdchgidx1, bdchgidx2);

   return (jboolean) earlier;
}

/** returns whether first bound change index belongs to an earlier applied bound change than second one */
JNIEXPORT
jboolean JNISCIPVAR(bdchgidxIsEarlierNonNull)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchgidx1,         /**< first bound change index, or NULL */
   jlong                 jbdchgidx2          /**< second bound change index, or NULL */
   )
{
   SCIP_BDCHGIDX* bdchgidx1;
   SCIP_BDCHGIDX* bdchgidx2;
   SCIP_Bool earlier;

   /* convert JNI pointer into C pointer */
   bdchgidx1 = (SCIP_BDCHGIDX*) (size_t) jbdchgidx1;
   assert(bdchgidx1 != NULL);

   bdchgidx2 = (SCIP_BDCHGIDX*) (size_t) jbdchgidx2;
   assert(bdchgidx2 != NULL);

   earlier = SCIPbdchgidxIsEarlierNonNull(bdchgidx1, bdchgidx2);

   return (jboolean) earlier;
}

/** returns old bound that was overwritten for given bound change information */
JNIEXPORT
jdouble JNISCIPVAR(bdchginfoGetOldbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_Real bound;

   /* convert JNI pointer into C pointer */
   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   bound = SCIPbdchginfoGetOldbound(bdchginfo);

   return (jdouble) bound;
}


/** returns new bound that was overwritten for given bound change information */
JNIEXPORT
jdouble JNISCIPVAR(bdchginfoGetNewbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_Real bound;

   /* convert JNI pointer into C pointer */
   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   bound = SCIPbdchginfoGetNewbound(bdchginfo);

   return (jdouble) bound;
}

/** returns variable that belongs to the given bound change information */
JNIEXPORT
jlong JNISCIPVAR(bdchginfoGetVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_VAR* var;
   SCIP_BDCHGINFO* bdchginfo;


   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   var = SCIPbdchginfoGetVar(bdchginfo);

   return (jlong)(size_t) var;
}

/** returns whether the bound change information belongs to a branching decision or a deduction */
JNIEXPORT
jint JNISCIPVAR(bdchginfoGetChgtype)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BOUNDCHGTYPE chgtype;

   /* convert JNI pointer into C pointer */
   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   chgtype = SCIPbdchginfoGetChgtype(bdchginfo);

   return (jint) chgtype;
}

/** returns whether the bound change information belongs to a lower or upper bound change */
JNIEXPORT
jint JNISCIPVAR(bdchginfoGetBoundtype)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BOUNDCHGTYPE boundtype;

   /* convert JNI pointer into C pointer */
   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   boundtype = SCIPbdchginfoGetBoundtype(bdchginfo);

   return (jint) boundtype;
}

/** returs depth level of given bound change information */
JNIEXPORT
jint JNISCIPVAR(bdchginfoGetDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   int depth;

   /* convert JNI pointer into C pointer */
   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   depth = SCIPbdchginfoGetDepth(bdchginfo);

   return (jint) depth;
}

/** returs bound change position in its depth level of given bound change information */
JNIEXPORT
jint JNISCIPVAR(bdchginfoGetPos)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   int pos;

   /* convert JNI pointer into C pointer */
   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   pos = SCIPbdchginfoGetPos(bdchginfo);

   return (jint) pos;
}

/** returns bound change index of given bound change information */
JNIEXPORT
jlong JNISCIPVAR(bdchginfoGetIdx)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BDCHGIDX* idx;


   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   idx = SCIPbdchginfoGetIdx(bdchginfo);

   return (jlong)(size_t) idx;
}

/** returns inference variable of given bound change information */
JNIEXPORT
jlong JNISCIPVAR(bdchginfoGetInferVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_VAR* var;


   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   var = SCIPbdchginfoGetInferVar(bdchginfo);

   return (jlong)(size_t) var;
}

/** returns inference constraint of given bound change information */
JNIEXPORT
jlong JNISCIPVAR(bdchginfoGetInferCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_CONS* cons;


   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   cons = SCIPbdchginfoGetInferCons(bdchginfo);

   return (jlong)(size_t) cons;
}

/** returns inference propagator of given bound change information, or NULL if no propagator was responsible */
JNIEXPORT
jlong JNISCIPVAR(bdchginfoGetInferProp)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_PROP* prop;


   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   prop = SCIPbdchginfoGetInferProp(bdchginfo);

   return (jlong)(size_t) prop;
}

/** returns inference user information of given bound change information */
JNIEXPORT
jint JNISCIPVAR(bdchginfoGetInferInfo)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   int inferinfo;

   /* convert JNI pointer into C pointer */
   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   inferinfo = SCIPbdchginfoGetInferInfo(bdchginfo);

   return (jint) inferinfo;
}

/** returns inference bound of inference variable of given bound change information */
JNIEXPORT
jint JNISCIPVAR(bdchginfoGetInferBoundtype)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BOUNDTYPE boundtype;

   /* convert JNI pointer into C pointer */
   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   boundtype = SCIPbdchginfoGetInferBoundtype(bdchginfo);

   return (jint) boundtype;
}

/** returns whether the bound change information belongs to a redundant bound change */
JNIEXPORT
jboolean JNISCIPVAR(bdchginfoIsRedundant)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_Bool redundant;

   /* convert JNI pointer into C pointer */
   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   redundant = SCIPbdchginfoIsRedundant(bdchginfo);

   return (jboolean) redundant;
}

/** returns whether the bound change has an inference reason (constraint or propagator), that can be resolved */
JNIEXPORT
jboolean JNISCIPVAR(bdchginfoHasInferenceReason)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_Bool reason;

   /* convert JNI pointer into C pointer */
   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   reason = SCIPbdchginfoHasInferenceReason(bdchginfo);

   return (jboolean) reason;
}

/** for two bound change informations belonging to the same variable and bound, returns whether the first bound change
 *  has a tighter new bound as the second bound change
 */
JNIEXPORT
jboolean JNISCIPVAR(bdchginfoIsTighter)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo1,        /**< first bound change information */
   jlong                 jbdchginfo2         /**< second bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo1;
   SCIP_BDCHGINFO* bdchginfo2;
   SCIP_Bool tighter;

   /* convert JNI pointer into C pointer */
   bdchginfo1 = (SCIP_BDCHGINFO*) (size_t) jbdchginfo1;
   assert(bdchginfo1 != NULL);

   bdchginfo2 = (SCIP_BDCHGINFO*) (size_t) jbdchginfo2;
   assert(bdchginfo2 != NULL);

   tighter = SCIPbdchginfoIsTighter(bdchginfo1, bdchginfo2);

   return (jboolean) tighter;
}

/** returns the new value of the bound in the bound change data */
JNIEXPORT
jdouble JNISCIPVAR(boundchgGetNewbound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jboundchg           /**< bound change data */
   )
{
   SCIP_BOUNDCHG* boundchg;
   SCIP_Real newbound;

   /* convert JNI pointer into C pointer */
   boundchg = (SCIP_BOUNDCHG*) (size_t) jboundchg;
   assert(boundchg != NULL);

   newbound = SCIPboundchgGetNewbound(boundchg);

   return (jdouble) newbound;
}

/** returns the variable of the bound change in the bound change data */
JNIEXPORT
jlong JNISCIPVAR(boundchgGetVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jboundchg           /**< bound change data */
   )
{
   SCIP_BOUNDCHG* boundchg;
   SCIP_VAR* var;


   boundchg = (SCIP_BOUNDCHG*) (size_t) jboundchg;
   assert(boundchg != NULL);

   var = SCIPboundchgGetVar(boundchg);

   return (jlong)(size_t) var;
}

/** returns the bound change type of the bound change in the bound change data */
JNIEXPORT
JNISCIP_ENUM JNISCIPVAR(boundchgGetBoundchgtype)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jboundchg           /**< bound change data */
   )
{
   SCIP_BOUNDCHG* boundchg;
   SCIP_BOUNDCHGTYPE bdchgtype;

   /* convert JNI pointer into C pointer */
   boundchg = (SCIP_BOUNDCHG*) (size_t) jboundchg;
   assert(boundchg != NULL);

   bdchgtype = SCIPboundchgGetBoundchgtype(boundchg);

   return (JNISCIP_ENUM) bdchgtype;
}

/** returns the bound type of the bound change in the bound change data */
JNIEXPORT
jint JNISCIPVAR(boundchgGetBoundtype)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jboundchg           /**< bound change data */
   )
{
   SCIP_BOUNDCHG* boundchg;
   SCIP_BOUNDTYPE bdtype;

   /* convert JNI pointer into C pointer */
   boundchg = (SCIP_BOUNDCHG*) (size_t) jboundchg;
   assert(boundchg != NULL);

   bdtype = SCIPboundchgGetBoundtype(boundchg);

   return (jint) bdtype;
}

/** returns whether the bound change is redundant due to a more global bound that is at least as strong */
JNIEXPORT
jboolean JNISCIPVAR(boundchgIsRedundant)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jboundchg           /**< bound change data */
   )
{
   SCIP_BOUNDCHG* boundchg;
   SCIP_Bool redundant;

   /* convert JNI pointer into C pointer */
   boundchg = (SCIP_BOUNDCHG*) (size_t) jboundchg;
   assert(boundchg != NULL);

   redundant = SCIPboundchgIsRedundant(boundchg);

   return (jboolean) redundant;
}

/** returns the number of bound changes in the domain change data */
JNIEXPORT
jint JNISCIPVAR(domchgGetNBoundchgs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdomchg             /**< domain change data */
   )
{
   SCIP_DOMCHG* domchg;
   int num;

   /* convert JNI pointer into C pointer */
   domchg = (SCIP_DOMCHG*) (size_t) jdomchg;
   assert(domchg != NULL);

   num = SCIPdomchgGetNBoundchgs(domchg);

   return (jint) num;
}

/** returns a particular bound change in the domain change data */
JNIEXPORT
 jlong JNISCIPVAR(domchgGetBoundchg)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jdomchg,            /**< domain change data */
   jint                  jpos                /**< position of the bound change in the domain change data */
   )
{
   SCIP_DOMCHG* domchg;
   SCIP_BOUNDCHG* boundchg;

   domchg = (SCIP_DOMCHG*) (size_t) jdomchg;
   assert(domchg != NULL);

   boundchg = SCIPdomchgGetBoundchg(domchg, (int) jpos);

   return (jlong)(size_t) boundchg;
}

/** returns left bound of open interval in hole */
JNIEXPORT
jdouble JNISCIPVAR(holelistGetLeft)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jholelist           /**< hole list pointer to hole of interest */
   )
{
   SCIP_HOLELIST* holelist;
   SCIP_Real lb;

   /* convert JNI pointer into C pointer */
   holelist = (SCIP_HOLELIST*) (size_t) jholelist;
   assert(holelist != NULL);

   lb = SCIPholelistGetLeft(holelist);

   return (jdouble) lb;
}

/** returns right bound of open interval in hole */
JNIEXPORT
jdouble JNISCIPVAR(holelistGetRight)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jholelist           /**< hole list pointer to hole of interest */
   )
{
   SCIP_HOLELIST* holelist;
   SCIP_Real rb;

   /* convert JNI pointer into C pointer */
   holelist = (SCIP_HOLELIST*) (size_t) jholelist;
   assert(holelist != NULL);

   rb = SCIPholelistGetRight(holelist);

   return (jdouble) rb;
}

/** returns next hole in list or NULL */
JNIEXPORT
jlong JNISCIPVAR(holelistGetNext)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jholelist           /**< hole list pointer to hole of interest */
   )
{
   SCIP_HOLELIST* holelist;
   SCIP_HOLELIST* next;

   holelist = (SCIP_HOLELIST*) (size_t) jholelist;
   assert(holelist != NULL);

   next = SCIPholelistGetNext(holelist);

   return (jlong)(size_t) next;
}








#if 0

/** gets objective value of variable in current SCIP_LP; the value can be different from the bound stored in the variable's own
 *  data due to diving, that operate only on the LP without updating the variables
 */
jdouble JNISCIPVAR(varGetObjLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real obj;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   obj = SCIvarGetObjLP(var);

   return (jdouble)obj;
}

/** gets lower bound of variable in current SCIP_LP; the bound can be different from the bound stored in the variable's own
 *  data due to diving or conflict analysis, that operate only on the LP without updating the variables
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetLbLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jlong                 jset                /**< global SCIP settings */
   )
{
   SCIP_VAR* var;
   SCIP_SET* set;
   SCIP_Real lb;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   set = (SCIP_SET*) (size_t) jset;
   assert(set != NULL);

   lb = SCIvarGetLbLP(var, set);

   return (jdouble) lb;
}

/** gets upper bound of variable in current SCIP_LP; the bound can be different from the bound stored in the variable's own
 *  data due to diving or conflict analysis, that operate only on the LP without updating the variables
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetUbLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jlong                 jset                /**< global SCIP settings */
   )
{
   SCIP_VAR* var;
   SCIP_SET* set;
   SCIP_Real ub;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   set = (SCIP_SET*) (size_t) jset;
   assert(set != NULL);

   ub = SCIvarGetUbLPP(var, set);

   return (jdouble) ub;
}

/** gets pseudo solution value of variable at current node */
JNIEXPORT
jdouble JNISCIPVAR(varGetPseudoSol_rec)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real sol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   sol = SCIPvarGetPseudoSol_rec(var);

   return (jdouble) sol;
}

/** returns for given variable the reduced cost */
JNIEXPORT
jdouble JNISCIPVAR(varGetRedcost)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jvarfixing,         /**< FALSE if for x == 0, TRUE for x == 1 */
   jlong                 jstat,              /**< problem statistics */
   jlong                 jlp                 /**< current LP data */
   )
{
   SCIP_VAR* var;
   SCIP_STAT* stat;
   SCIP_LP* lp;

   SCIP_Real redcost;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;

   stat = (SCIP_STAT*) (size_t) jstat;

   lp = (SCIP_LP*) (size_t) jlp;

   redcost = SCIPvarGetRedcost(var, (SCIP_Bool)jvarfixing, stat, lp);

   return (jdouble) redcost;
}

/** returns for the given binary variable the reduced cost which are given by the variable itself and its implication if
 *  the binary variable is fixed to the given value
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetImplRedcost)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jboolean              jvarfixing,         /**< FALSE if for x == 0, TRUE for x == 1 */
   jlong                 jstat,              /**< problem statistics */
   jlong                 jlp                 /**< current LP data */
   )
{
   SCIP_VAR* var;
   SCIP_STAT* stat;
   SCIP_LP* lp;

   SCIP_Real redcost;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;

   stat = (SCIP_STAT*) (size_t) jstat;

   lp = (SCIP_LP*) (size_t) jlp;

   redcost = SCIPvarGetImplRedcost(var, (SCIP_Bool)jvarfixing, stat, lp);

   return (jdouble) redcost;
}

/** returns the solution value of the problem variable in the relaxation solution */
JNIEXPORT
jdouble JNISCIPVAR(varGetRelaxSol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   jlong                 jset                /**< global SCIP settings */
   )
{
   SCIP_VAR* var;
   SCIP_SET* set;
   SCIP_Real sol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   set = (SCIP_SET*) (size_t) jset;
   assert(var != NULL);

   sol = SCIPvarGetRelaxSol(var, set);

   return (jdouble) sol;
}

/** returns the solution value of the transformed problem variable in the relaxation solution */
JNIEXPORT
jdouble JNISCIPVAR(varGetRelaxSolTransVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar                /**< JNI problem variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real sol;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   sol = SCIPvarGetRelaxSolTransVar(var);

   return (jdouble) sol;
}

/** gets the variable's pseudo cost value for the given step size "solvaldelta" in the variable's LP solution value */
JNIEXPORT
jdouble JNISCIPVAR(varGetPseudocost)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   jlong                 jstat,              /**< problem statistics */
   jdouble               jsolvaldelta        /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP_VAR* var;
   SCIP_STAT* stat;
   SCIP_Real num;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   stat = (SCIP_STAT*) (size_t) jstat;
   assert(stat != NULL);

   num = SCIPvarGetPseudocost(var, stat, (SCIP_Real) jsolvaldelta);

   return (jdouble) num;
}

/** gets the variable's pseudo cost value for the given step size "solvaldelta" in the variable's LP solution value,
 *  only using the pseudo cost information of the current run
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetPseudocostCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   jlong                 jstat,              /**< problem statistics */
   jdouble               jsolvaldelta        /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP_VAR* var;
   SCIP_STAT* stat;
   SCIP_Real num;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   stat = (SCIP_STAT*) (size_t) jstat;
   assert(stat != NULL);

   num = SCIPvarGetPseudocostCurrentRun(var, stat, (SCIP_Real) jsolvaldelta);

   return (jdouble) num;
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction */
JNIEXPORT
jdouble JNISCIPVAR(varGetPseudocostCount)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_Real num;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   num = SCIPvarGetPseudocostCount(var, (SCIP_BRANCHDIR) jdir);

   return (jdouble) num;
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetPseudocostCountCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_Real num;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   num = SCIPvarGetPseudocostCountCurrentRun(var, (SCIP_BRANCHDIR) jdir);

   return (jdouble) num;
}

/** gets the average conflict length in given direction due to branching on the variable */
JNIEXPORT
jdouble JNISCIPVAR(varGetAvgConflictlength)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   num = SCIPvarGetAvgConflictlength(var, (SCIP_BRANCHDIR) jdir);

   return (jdouble) num;
}

/**  gets the average conflict length in given direction due to branching on the variable
 *   in the current run
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetAvgConflictlengthCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< JNI problem variable */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_Real num;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   num = SCIPvarGetAvgConflictlengthCurrentRun(var, (SCIP_BRANCHDIR) jdir);

   return (jdouble) num;
}

/** returns the variable's VSIDS score */
JNIEXPORT
jdouble JNISCIPVAR(varGetVSIDS_rec)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   jlong                 jstat,              /**< problem statistics */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_STAT* stat;
   SCIP_Real vsids;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   stat = (SCIP_STAT*) (size_t) jstat;
   assert(stat != NULL);

   vsids = SCIPvarGetVSIDS_rec(var, stat, (SCIP_BRANCHDIR) jdir);

   return (jdouble) vsids;
}

/** returns the variable's VSIDS score only using conflicts of the current run */
JNIEXPORT
jdouble JNISCIPVAR(varGetVSIDSCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   jlong                 jstat,              /**< problem statistics */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_STAT* stat;
   SCIP_Real vsids;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   stat = (SCIP_STAT*) (size_t) jstat;
   assert(stat != NULL);

   vsids = SCIPvarGetVSIDSCurrentRun(var, stat, (SCIP_BRANCHDIR) jdir);

   return (jdouble) vsids;
}



/** returns the average number of inferences found after branching on the variable in given direction */
JNIEXPORT
jdouble JNISCIPVAR(varGetAvgInferences)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   jlong                 jstat,              /**< problem statistics */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_STAT* stat;
   SCIP_Real ai;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   stat = (SCIP_STAT*) (size_t) jstat;
   assert(stat != NULL);

   ai = SCIPvarGetAvgInferences(var, stat, (SCIP_BRANCHDIR) jdir);

   return (jdouble) ai;
}

/** returns the average number of inferences found after branching on the variable in given direction
 *  in the current run
 */
JNIEXPORT
jdouble JNISCIPVAR(varGetAvgInferencesCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   jlong                 jstat,              /**< problem statistics */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_STAT* stat;
   SCIP_Real ai;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   stat = (SCIP_STAT*) (size_t) jstat;
   assert(stat != NULL);

   ai = SCIPvarGetAvgInferencesCurrentRun(var, stat, (SCIP_BRANCHDIR) jdir);

   return (jdouble) ai;
}



/** returns the average number of cutoffs found after branching on the variable in given direction */
JNIEXPORT
jdouble JNISCIPVAR(varGetAvgCutoffs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   jlong                 jstat,              /**< problem statistics */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_STAT* stat;
   SCIP_Real ac;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   stat = (SCIP_STAT*) (size_t) jstat;
   assert(stat != NULL);

   ac = SCIPvarGetAvgCutoffs(var, stat, (SCIP_BRANCHDIR) jdir);

   return (jdouble) ac;
}

/** returns the average number of cutoffs found after branching on the variable in given direction in the current run */
JNIEXPORT
jdouble JNISCIPVAR(varGetAvgCutoffsCurrentRun)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   jlong                 jstat,              /**< problem statistics */
   JNISCIP_ENUM          jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_STAT* stat;
   SCIP_Real ac;

   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   stat = (SCIP_STAT*) (size_t) jstat;
   assert(stat != NULL);

   ac = SCIPvarGetAvgCutoffsCurrentRun(var, stat, (SCIP_BRANCHDIR) jdir);

   return (jdouble) ac;
}

/** returns at which depth in the tree a bound change was applied to the variable that conflicts with the
 *  given bound; returns -1 if the bound does not conflict with the current local bounds of the variable
 */
JNIEXPORT
jint JNISCIPVAR(varGetConflictingBdchgDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   jlong                 jset,               /**< global SCIP settings */
   JNISCIP_ENUM          jboundtype,         /**< bound type of the conflicting bound */
   jdouble               jbound              /**< conflicting bound */
   )
{
   SCIP_VAR* var;
   SCIP_SET* set;
   SCIP_BOUNDTYPE boundtype;
   SCIP_Real bound;
   int cbd;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert(var != NULL);

   set = (SCIP_SET*) (size_t) jset;
   assert(set != NULL);

   boundtype = (SCIP_BOUNDTYPE) jboundtype;
   bound = (SCIP_Real) jbound;

   cbd = SCIPvarGetConflictingBdchgDepth(var, set, boundtype, bound);

   return (jint) cbd;
}

/** returns the variable's VSIDS score */
JNIEXPORT
jdouble JNISCIPVAR(varGetVSIDS)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jvar,               /**< problem variable */
   jlong                 jstat,              /**< problem statistics */
   jint                  jdir                /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_VAR* var;
   SCIP_STAT* stat;
   SCIP_BRANCHDIR dir;
   SCIP_Real VSIDS;

   /* convert JNI pointer into C pointer */
   var = (SCIP_VAR*) (size_t) jvar;
   assert (var != NULL);

   stat = (SCIP_STAT*) (size_t) jstat;
   dir = (SCIP_BRANCHDIR) jdir;

   VSIDS = SCIPvarGetVSIDS(var, stat, dir);

   return (jdouble) VSIDS;
}

/** returns the position of the bound change index */
JNIEXPORT
jint JNISCIPVAR(bdchgidxGetPos)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchgidx           /**< bound change index */
   )
{
   SCIP_BDCHGIDX* bdchgidx;
   int pos;

   /* convert JNI pointer into C pointer */
   bdchgidx = (SCIP_BDCHGIDX*) (size_t) jbdchgidx;
   assert(bdchgidx != NULL);

   pos = SCIPbdchgidxGetPos(bdchgidx);

   return (jint) pos;
}

/** returns the relaxed bound change type */
JNIEXPORT
jdouble JNISCIPVAR(bdchginfoGetRelaxedBound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jbdchginfo          /**< bound change information */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_Real relaxed;

   /* convert JNI pointer into C pointer */
   bdchginfo = (SCIP_BDCHGINFO*) (size_t) jbdchginfo;
   assert(bdchginfo != NULL);

   relaxed = SCIPbdchginfoGetRelaxedBound(bdchginfo);

   return (jdouble) relaxed;
}


/** gets corresponding objective value of active, fixed, or multi-aggregated problem variable of given variable
 *  e.g. obj(x) = 1 this method returns for ~x the value -1
 */
 JNIEXPORT
 jdouble JNISCIPVAR(varGetAggregatedObj)(
    JNIEnv*               env,                /**< JNI environment variable */
    jobject               jobj,               /**< JNI class pointer */
    jlong                 jvar                /**< JNI problem variable */
    )
 {
    SCIP_VAR* var;
    SCIP_Real* aggrobj;

    /* convert JNI pointer into C pointer */
    var = (SCIP_VAR*)(size_t)jvar;
    assert(var != NULL);


    JNISCIP_CALL( SCIPvarGetAggregatedObj(var, aggrobj) );

    return aggrobj;
 }
#endif

/**@} */
