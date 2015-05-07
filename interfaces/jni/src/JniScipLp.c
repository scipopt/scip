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

#include "JniScipLp.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include <string.h>

/** sorts column entries of linked rows currently in the LP such that lower row indices precede higher ones */
JNIEXPORT
void JNISCIPLP(colSort)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< column to be sorted */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   SCIPcolSort(col);
}

/** gets objective value of column */
JNIEXPORT
jdouble JNISCIPLP(colGetObj)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jdouble) SCIPcolGetObj(col);
}

/** gets lower bound of column */
JNIEXPORT
jdouble JNISCIPLP(colGetLb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jdouble) SCIPcolGetLb(col);
}

/** gets upper bound of column */
JNIEXPORT
jdouble JNISCIPLP(colGetUb)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jdouble) SCIPcolGetUb(col);
}

/** gets best bound of column with respect to the objective function */
JNIEXPORT
jdouble JNISCIPLP(colGetBestBound)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jdouble) SCIPcolGetBestBound(col);
}

/** gets the primal LP solution of a column */
JNIEXPORT
jdouble JNISCIPLP(colGetPrimsol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jdouble) SCIPcolGetPrimsol(col);
}

/** gets the minimal LP solution value, this column ever assumed */
JNIEXPORT
jdouble JNISCIPLP(colGetMinPrimsol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jdouble) SCIPcolGetMinPrimsol(col);
}

/** gets the maximal LP solution value, this column ever assumed */
JNIEXPORT
jdouble JNISCIPLP(colGetMaxPrimsol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jdouble) SCIPcolGetMaxPrimsol(col);
}

/** gets the basis status of a column in the LP solution; only valid for LPs with status SCIP_LPSOLSTAT_OPTIMAL
 *  and with SCIPisLPSolBasic(scip) == TRUE; returns SCIP_BASESTAT_ZERO for columns not in the current SCIP_LP
 */
JNIEXPORT
jint JNISCIPLP(colGetBasisStatus)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jint) SCIPcolGetBasisStatus(col);
}

/** gets variable this column represents */
JNIEXPORT
jlong JNISCIPLP(colGetVar)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jlong) SCIPcolGetVar(col);
}

/** gets unique index of col */
JNIEXPORT
jint JNISCIPLP(colGetIndex)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jint) SCIPcolGetIndex(col);
}

/** returns whether the associated variable is of integral type (binary, integer, implicit integer) */
JNIEXPORT
jboolean JNISCIPLP(colIsIntegral)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jboolean) SCIPcolIsIntegral(col);
}

/** returns TRUE iff column is removable from the LP (due to aging or cleanup) */
JNIEXPORT
jboolean JNISCIPLP(colIsRemovable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jboolean) SCIPcolIsRemovable(col);
}

/** gets position of column in current LP, or -1 if it is not in LP */
JNIEXPORT
jint JNISCIPLP(colGetLPPos)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jint) SCIPcolGetLPPos(col);
}

/** gets depth in the tree where the column entered the LP, or -1 if it is not in LP */
JNIEXPORT
jint JNISCIPLP(colGetLPDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jint) SCIPcolGetLPDepth(col);
}

/** returns TRUE iff column is member of current LP */
JNIEXPORT
jboolean JNISCIPLP(colIsInLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jboolean) SCIPcolIsInLP(col);
}

/** get number of nonzero entries in column vector */
JNIEXPORT
jint JNISCIPLP(colGetNNonz)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jint) SCIPcolGetNNonz(col);
}

/** get number of nonzero entries in column vector, that correspond to rows currently in the SCIP_LP;
 *  Warning! This method is only applicable on columns, that are completely linked to their rows (e.g. a column
 *  that is in the current LP and the LP was solved, or a column that was in a solved LP and didn't change afterwards
 */
JNIEXPORT
jint JNISCIPLP(colGetNLPNonz)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jint) SCIPcolGetNLPNonz(col);
}

/** gets array with rows of nonzero entries */
jlongArray JNISCIPLP(colGetRows)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   jlongArray jrows;
   int size;
   SCIP_ROW** rows;
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   size = SCIPcolGetNLPNonz(col);

   jrows = (*env)->NewLongArray(env, size);

   if (jrows == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   rows = SCIPcolGetRows(col);

   (*env)->SetLongArrayRegion(env, jrows, 0, size, (jlong*)(*rows));

   return jrows;
}

/** gets array with coefficients of nonzero entries */
jdoubleArray JNISCIPLP(colGetVals)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   jdoubleArray jvals;
   int size;
   SCIP_Real* vals;
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   size = SCIPcolGetNLPNonz(col);

   jvals = (*env)->NewDoubleArray(env, size);

   if (jvals == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   vals = SCIPcolGetVals(col);

   (*env)->SetDoubleArrayRegion(env, jvals, 0, size, (jdouble*)vals);

   return jvals;
}

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given column, or -1 if strong branching was never applied to the column in current run
 */
JNIEXPORT
jlong JNISCIPLP(colGetStrongbranchNode)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jlong) SCIPcolGetStrongbranchNode(col);
}

/** gets number of times, strong branching was applied in current run on the given column */
JNIEXPORT
jint JNISCIPLP(colGetNStrongbranchs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jcol                /**< LP column */
   )
{
   SCIP_COL* col;

   /* convert JNI pointer into C pointer */
   col = (SCIP_COL*) (size_t) jcol;
   assert(col != NULL);

   return (jint) SCIPcolGetNStrongbranchs(col);
}

/** gets opposite bound type of given bound type */
JNIEXPORT
jint JNISCIPLP(boundtypeOpposite)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jint                  boundtype           /**< type of bound (lower or upper) */
   )
{
   return (jint) SCIPboundtypeOpposite( (SCIP_BOUNDTYPE)boundtype );
}

/** locks an unmodifiable row, which forbids further changes; has no effect on modifiable rows */
JNIEXPORT
void JNISCIPLP(rowLock)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   SCIProwLock(row);
}

/** unlocks a lock of an unmodifiable row; a row with no sealed lock may be modified; has no effect on modifiable rows */
JNIEXPORT
void JNISCIPLP(rowUnlock)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   SCIProwUnlock(row);
}

/** returns the scalar product of the coefficient vectors of the two given rows */
JNIEXPORT
jdouble JNISCIPLP(rowGetScalarProduct)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow1,              /**< first LP row */
   jlong                 jrow2               /**< second LP row */
   )
{
   SCIP_ROW* row1;
   SCIP_ROW* row2;

   /* convert JNI pointer into C pointer */
   row1 = (SCIP_ROW*) (size_t) jrow1;
   assert(row1 != NULL);

   row2 = (SCIP_ROW*) (size_t) jrow2;
   assert(row2 != NULL);

   return (SCIP_Real) SCIProwGetScalarProduct(row1, row2);
}

/** returns the degree of parallelism between the hyperplanes defined by the two row vectors v, w:
 *  p = |v*w|/(|v|*|w|);
 *  the hyperplanes are parallel, iff p = 1, they are orthogonal, iff p = 0
 */
JNIEXPORT
jdouble JNISCIPLP(rowGetParallelism)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow1,              /**< first LP row */
   jlong                 jrow2,              /**< second LP row */
   jchar                 orthofunc           /**< function used for calc. scalar prod. ('e'uclidean, 'd'iscrete) */
   )
{
   SCIP_ROW* row1;
   SCIP_ROW* row2;

   /* convert JNI pointer into C pointer */
   row1 = (SCIP_ROW*) (size_t) jrow1;
   assert(row1 != NULL);

   row2 = (SCIP_ROW*) (size_t) jrow2;
   assert(row2 != NULL);

   return (SCIP_Real) SCIProwGetParallelism(row1, row2, (char)orthofunc);
}

/** returns the degree of orthogonality between the hyperplanes defined by the two row vectors v, w:
 *  o = 1 - |v*w|/(|v|*|w|);
 *  the hyperplanes are orthogonal, iff p = 1, they are parallel, iff p = 0
 */
JNIEXPORT
jdouble JNISCIPLP(rowGetOrthogonality)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow1,              /**< first LP row */
   jlong                 jrow2,              /**< second LP row */
   jchar                 orthofunc           /**< function used for calc. scalar prod. ('e'uclidean, 'd'iscrete) */
   )
{
   SCIP_ROW* row1;
   SCIP_ROW* row2;

   /* convert JNI pointer into C pointer */
   row1 = (SCIP_ROW*) (size_t) jrow1;
   assert(row1 != NULL);

   row2 = (SCIP_ROW*) (size_t) jrow2;
   assert(row2 != NULL);

   return (SCIP_Real) SCIProwGetOrthogonality(row1, row2, (char)orthofunc);
}

/** sorts row entries such that LP columns precede non-LP columns and inside both parts lower column indices precede
 *  higher ones
 */
JNIEXPORT
void JNISCIPLP(rowSort)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   SCIProwSort(row);
}

/** get number of nonzero entries in row vector */
JNIEXPORT
jint JNISCIPLP(rowGetNNonz)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jint) SCIProwGetNNonz(row);
}

/** get number of nonzero entries in row vector, that correspond to columns currently in the SCIP_LP;
 *  Warning! This method is only applicable on rows, that are completely linked to their columns (e.g. a row
 *  that is in the current LP and the LP was solved, or a row that was in a solved LP and didn't change afterwards
 */
JNIEXPORT
jint JNISCIPLP(rowGetNLPNonz)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jint) SCIProwGetNLPNonz(row);
}

/** gets array with columns of nonzero entries */
JNIEXPORT
jlongArray JNISCIPLP(rowGetCols)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;
   jlongArray jcols;
   int size;
   SCIP_COL** cols;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   size = SCIProwGetNNonz(row);

   jcols = (*env)->NewLongArray(env, size);

   if (jcols == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   cols = SCIProwGetCols(row);

   (*env)->SetLongArrayRegion(env, jcols, 0, size, (jlong*)(*cols));

   return jcols;
}

/** gets array with coefficients of nonzero entries */
JNIEXPORT
jdoubleArray JNISCIPLP(rowGetVals)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;
   jdoubleArray jvals;
   int size;
   SCIP_Real* vals;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   size = SCIProwGetNNonz(row);

   jvals = (*env)->NewDoubleArray(env, size);

   if (jvals == NULL) {
      SCIPerrorMessage("Out of Memory\n");
      JNISCIP_CALL( SCIP_ERROR );
      return 0;
   }

   vals = SCIProwGetVals(row);

   (*env)->SetDoubleArrayRegion(env, jvals, 0, size, (jdouble*)vals);

   return jvals;
}

/** gets constant shift of row */
JNIEXPORT
jdouble JNISCIPLP(rowGetConstant)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIProwGetConstant(row);
}

/** gets Euclidean norm of row vector */
JNIEXPORT
jdouble JNISCIPLP(rowGetNorm)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIProwGetNorm(row);
}

/** gets sum norm of row vector (sum of absolute values of coefficients) */
JNIEXPORT
jdouble JNISCIPLP(rowGetSumNorm)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIProwGetSumNorm(row);
}

/** returns the left hand side of the row */
JNIEXPORT
jdouble JNISCIPLP(rowGetLhs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIProwGetLhs(row);
}

/** returns the right hand side of the row */
JNIEXPORT
jdouble JNISCIPLP(rowGetRhs)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIProwGetRhs(row);
}

/** gets the dual LP solution of a row */
JNIEXPORT
jdouble JNISCIPLP(rowGetDualsol)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIProwGetDualsol(row);
}

/** gets the dual Farkas coefficient of a row in an infeasible LP */
JNIEXPORT
jdouble JNISCIPLP(rowGetDualfarkas)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jdouble) SCIProwGetDualfarkas(row);
}

/** gets the basis status of a row in the LP solution; only valid for LPs with status SCIP_LPSOLSTAT_OPTIMAL
 *  and with SCIPisLPSolBasic(scip) == TRUE; returns SCIP_BASESTAT_BASIC for rows not in the current SCIP_LP
 */
JNIEXPORT
jint JNISCIPLP(rowGetBasisStatus)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jint) SCIProwGetBasisStatus(row);
}

/** returns the name of the row */
JNIEXPORT
jstring JNISCIPLP(rowGetName)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;
   const char* name;
   jstring jname;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   name = SCIProwGetName(row);

   jname = (*env)->NewStringUTF(env, name);

   return jname;
}

/** gets unique index of row */
JNIEXPORT
jint JNISCIPLP(rowGetIndex)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jint) SCIProwGetIndex(row);
}

/** gets age of row */
JNIEXPORT
jint JNISCIPLP(rowGetAge)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jint) SCIProwGetAge(row);
}

/** gets rank of row */
JNIEXPORT
jint JNISCIPLP(rowGetRank)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jint) SCIProwGetRank(row);
}

/** returns TRUE iff the activity of the row (without the row's constant) is always integral in a feasible solution */
JNIEXPORT
jboolean JNISCIPLP(rowIsIntegral)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jboolean) SCIProwIsIntegral(row);
}

/** returns TRUE iff row is only valid locally */
JNIEXPORT
jboolean JNISCIPLP(rowIsLocal)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jboolean) SCIProwIsLocal(row);
}

/** returns TRUE iff row is modifiable during node processing (subject to column generation) */
JNIEXPORT
jboolean JNISCIPLP(rowIsModifiable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jboolean) SCIProwIsModifiable(row);
}

/** returns TRUE iff row is removable from the LP (due to aging or cleanup) */
JNIEXPORT
jboolean JNISCIPLP(rowIsRemovable)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jboolean) SCIProwIsRemovable(row);
}

/** returns type of origin that created the row */
JNIEXPORT
jint JNISCIPLP(rowGetOrigintype)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jint) SCIProwGetOrigintype(row);
}

/** returns origin constraint handler that created the row (NULL if not available) */
JNIEXPORT
jlong JNISCIPLP(rowGetOriginCons)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jlong) (size_t) SCIProwGetOriginCons(row);
}

/** returns origin separator that created the row (NULL if not available) */
JNIEXPORT
jlong JNISCIPLP(rowGetOriginSepa)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jlong) (size_t) SCIProwGetOriginSepa(row);
}

/** returns TRUE iff row is member of the global cut pool */
JNIEXPORT
jboolean JNISCIPLP(rowIsInGlobalCutpool)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jboolean) SCIProwIsInGlobalCutpool(row);
}

/** gets position of row in current LP, or -1 if it is not in LP */
JNIEXPORT
jint JNISCIPLP(rowGetLPPos)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jint) SCIProwGetLPPos(row);
}

/** gets depth in the tree where the row entered the LP, or -1 if it is not in LP */
JNIEXPORT
jint JNISCIPLP(rowGetLPDepth)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jint) SCIProwGetLPDepth(row);
}

/** returns TRUE iff row is member of current LP */
JNIEXPORT
jboolean JNISCIPLP(rowIsInLP)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow                /**< LP row */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   return (jboolean) SCIProwIsInLP(row);
}

/** changes the rank of LP row */
JNIEXPORT
void JNISCIPLP(rowChgRank)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jrow,               /**< LP row */
   jint                  rank                /**< new value for rank */
   )
{
   SCIP_ROW* row;

   /* convert JNI pointer into C pointer */
   row = (SCIP_ROW*) (size_t) jrow;
   assert(row != NULL);

   SCIProwChgRank(row, (int)rank);
}
