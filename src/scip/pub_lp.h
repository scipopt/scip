/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_lp.h,v 1.1 2003/12/01 14:41:29 bzfpfend Exp $"

/**@file   lp.h
 * @brief  internal methods for LP management
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_LP_H__
#define __PUB_LP_H__


#include <stdio.h>

#include "def.h"
#include "memory.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_lp.h"
#include "type_var.h"
#include "type_sol.h"

#ifdef NDEBUG
#include "struct_lp.h"
#endif



/*
 * Column methods
 */

/** output column to file stream */
extern
void SCIPcolPrint(
   COL*             col,                /**< LP column */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets lower bound of column */
extern
Real SCIPcolGetLb(
   COL*             col                 /**< LP column */
   );

/** gets upper bound of column */
extern
Real SCIPcolGetUb(
   COL*             col                 /**< LP column */
   );

/** gets best bound of column with respect to the objective function */
extern
Real SCIPcolGetBestBound(
   COL*             col                 /**< LP column */
   );

/** gets the primal LP solution of a column */
extern
Real SCIPcolGetPrimsol(
   COL*             col                 /**< LP column */
   );

/** gets variable this column represents */
extern
VAR* SCIPcolGetVar(
   COL*             col                 /**< LP column */
   );

/** returns TRUE iff column is removeable from the LP (due to aging or cleanup) */
extern
Bool SCIPcolIsRemoveable(
   COL*             col                 /**< LP column */
   );

/** gets position of column in actual LP, or -1 if it is not in LP */
extern
int SCIPcolGetLPPos(
   COL*             col                 /**< LP column */
   );

/** returns TRUE iff column is member of actual LP */
extern
Bool SCIPcolIsInLP(
   COL*             col                 /**< LP column */
   );

/** get number of nonzero entries in column vector */
extern
int SCIPcolGetNNonz(
   COL*             col                 /**< LP column */
   );

/** gets array with rows of nonzero entries */
extern
ROW** SCIPcolGetRows(
   COL*             col                 /**< LP column */
   );

/** gets array with coefficients of nonzero entries */
extern
Real* SCIPcolGetVals(
   COL*             col                 /**< LP column */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPcolGetLb(col)               ((col)->lb)
#define SCIPcolGetUb(col)               ((col)->ub)
#define SCIPcolGetBestBound(col)        ((col)->obj >= 0.0 ? (col)->lb : (col)->ub)
#define SCIPcolGetPrimsol(col)          ((col)->lppos >= 0 ? (col)->primsol : 0.0)
#define SCIPcolGetVar(col)              ((col)->var)
#define SCIPcolIsRemoveable(col)        ((col)->removeable)
#define SCIPcolGetLPPos(col)            ((col)->lppos)
#define SCIPcolIsInLP(col)              ((col)->lppos >= 0)
#define SCIPcolGetNNonz(col)            ((col)->len)
#define SCIPcolGetRows(col)             ((col)->rows)
#define SCIPcolGetVals(col)             ((col)->vals)

#endif




/*
 * Row methods
 */

/** locks an unmodifiable row, which forbids further changes */
extern
RETCODE SCIProwLock(
   ROW*             row                 /**< LP row */
   );

/** unlocks a lock of a row; a row with no sealed lock may be modified */
extern
RETCODE SCIProwUnlock(
   ROW*             row                 /**< LP row */
   );

/** output row to file stream */
extern
void SCIProwPrint(
   ROW*             row,                /**< LP row */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** get number of nonzero entries in row vector */
extern
int SCIProwGetNNonz(
   ROW*             row                 /**< LP row */
   );

/** gets array with columns of nonzero entries */
extern
COL** SCIProwGetCols(
   ROW*             row                 /**< LP row */
   );

/** gets array with coefficients of nonzero entries */
extern
Real* SCIProwGetVals(
   ROW*             row                 /**< LP row */
   );

/** gets constant shift of row */
extern
Real SCIProwGetConstant(
   ROW*             row                 /**< LP row */
   );

/** get euclidean norm of row vector */
extern
Real SCIProwGetNorm(
   ROW*             row                 /**< LP row */
   );

/** returns the left hand side of the row */
extern
Real SCIProwGetLhs(
   ROW*             row                 /**< LP row */
   );

/** returns the right hand side of the row */
extern
Real SCIProwGetRhs(
   ROW*             row                 /**< LP row */
   );

/** gets the dual LP solution of a row */
extern
Real SCIProwGetDualsol(
   ROW*             row                 /**< LP row */
   );

/** returns the name of the row */
extern
const char* SCIProwGetName(
   ROW*             row                 /**< LP row */
   );

/** gets unique index of row */
extern
int SCIProwGetIndex(
   ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is only valid locally */
extern
Bool SCIProwIsLocal(
   ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is modifiable during node processing (subject to column generation) */
extern
Bool SCIProwIsModifiable(
   ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is removeable from the LP (due to aging or cleanup) */
extern
Bool SCIProwIsRemoveable(
   ROW*             row                 /**< LP row */
   );

/** gets position of row in actual LP, or -1 if it is not in LP */
extern
int SCIProwGetLPPos(
   ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is member of actual LP */
extern
Bool SCIProwIsInLP(
   ROW*             row                 /**< LP row */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIProwGetNNonz(row)            ((row)->len)
#define SCIProwGetCols(row)             ((row)->cols)
#define SCIProwGetVals(row)             ((row)->vals)
#define SCIProwGetConstant(row)         ((row)->constant)
#define SCIProwGetNorm(row)             (sqrt((row)->sqrnorm))
#define SCIProwGetLhs(row)              ((row)->lhs)
#define SCIProwGetRhs(row)              ((row)->rhs)
#define SCIProwGetDualsol(row)          ((row)->lppos >= 0 ? (row)->dualsol : 0.0)
#define SCIProwGetName(row)             ((row)->name)
#define SCIProwGetIndex(row)            ((row)->index)
#define SCIProwIsLocal(row)             ((row)->local)
#define SCIProwIsModifiable(row)        ((row)->modifiable)
#define SCIProwIsRemoveable(row)        ((row)->removeable)
#define SCIProwGetLPPos(row)            ((row)->lppos)
#define SCIProwIsInLP(row)              ((row)->lppos >= 0)

#endif


#endif
