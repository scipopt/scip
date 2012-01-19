/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    matrix.h
 * @brief   constraint matrix data structure
 * @author  Dieter Weninger
 * @author  Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MATRIX_H__
#define __SCIP_MATRIX_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** problem data */
struct ConstraintMatrix
{
   /* column data */
   SCIP_Real*            colmatval;
   int*                  colmatind;
   int*                  colmatbeg;
   int*                  colmatcnt;
   int                   ncols;
   SCIP_Real*            lb;
   SCIP_Real*            ub;

   SCIP_VAR**            vars;

   /* row data */
   SCIP_Real*            rowmatval;
   int*                  rowmatind;
   int*                  rowmatbeg;
   int*                  rowmatcnt;
   int                   nrows;
   SCIP_Real*            lhs;
   SCIP_Real*            rhs;

   SCIP_CONS**           conss;

   /* sparsity counter */
   int                   nnonzs;

   /* data for row bound analysis */
   SCIP_Real*            minactivity;
   SCIP_Real*            maxactivity;
   int*                  minactivityneginf;
   int*                  minactivityposinf;
   int*                  maxactivityneginf;
   int*                  maxactivityposinf;
};
typedef struct ConstraintMatrix CONSTRAINTMATRIX;

/** initialize matrix from scip */
SCIP_RETCODE initMatrix(
   SCIP*                 scip,       /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,     /**< constraint matrix object to be initialized */
   SCIP_Bool*            initialized /**< was the initialization successful? */
   );

/** frees the constraint matrix */
void freeMatrix(
   SCIP*                 scip,   /**< current SCIP instance */
   CONSTRAINTMATRIX**    matrix  /**< constraint matrix object */
   );

int getNCols(
   CONSTRAINTMATRIX*     matrix
   );

int getNRows(
   CONSTRAINTMATRIX*     matrix
   );

SCIP_Real rowGetLhs(
   CONSTRAINTMATRIX*     matrix,
   int                   row
   );

SCIP_Real rowGetRhs(
   CONSTRAINTMATRIX*     matrix,
   int                   row
   );

SCIP_Real getRowMinActivity(
   CONSTRAINTMATRIX*     matrix,
   int                   row
   );

SCIP_Real getRowMaxActivity(
   CONSTRAINTMATRIX*     matrix,
   int                   row
   );

SCIP_Real* colGetVals(
   CONSTRAINTMATRIX*     matrix,
   int                   col
   );

int* colGetRows(
   CONSTRAINTMATRIX*     matrix,
   int                   col
   );

int colGetNRows(
   CONSTRAINTMATRIX*     matrix,
   int                   col
   );

/** write matrix in column major format into a file */
SCIP_RETCODE writeMatrixColumns(
   SCIP*                 scip,
   CONSTRAINTMATRIX*     matrix,
   const char*           filename,
   int*                  nentries
   );

/** write matrix in row major format into a file */
SCIP_RETCODE writeMatrixRows(
   SCIP*                 scip,
   CONSTRAINTMATRIX*     matrix,
   const char*           filename,
   int*                  nentries
   );

#ifdef __cplusplus
}
#endif

#endif
