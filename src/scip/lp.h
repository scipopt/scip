/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lp.h
 * @brief  LP management datastructures and methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __LP_H__
#define __LP_H__

#include "def.h"
#include "mem.h"
#include "set.h"
#include "stat.h"
#include "lpi.h"


enum VarType
{
   SCIP_VARTYPE_BINARY    = 0,          /**< binary variable: $x \in \{0,1\}$ */
   SCIP_VARTYPE_INTEGER   = 1,          /**< integer variable: $x \in \{lb, \ldots, \ub\}$ */
   SCIP_VARTYPE_IMPLINT   = 2,          /**< implicit integer variable: continous variable, that is allways integral */
   SCIP_VARTYPE_CONTINOUS = 3           /**< continous variable: $x \in [lb,ub] */
};
typedef enum VarType VARTYPE;

enum BoundType
{
   SCIP_BOUNDTYPE_LOWER = 0,            /**< lower bound */
   SCIP_BOUNDTYPE_UPPER = 1             /**< upper bound */
};
typedef enum BoundType BOUNDTYPE;

typedef struct Vardom VARDOM;           /**< datastructures for storing domains of variables */
typedef struct Hole HOLE;               /**< hole in a domain of an integer variable */
typedef struct Holelist HOLELIST;       /**< list of holes in a domain of an integer variable */
typedef struct VardomChg VARDOMCHG;     /**< changes in domains of variables */
typedef struct BoundChg BOUNDCHG;       /**< changes in bounds of variables */
typedef struct HoleChg HOLECHG;         /**< changes in holelist of variables */
typedef struct Col COL;                 /**< column of an LP */
typedef struct ColList COLLIST;         /**< list of LP columns */
typedef struct Row ROW;                 /**< row of an LP */
typedef struct RowList ROWLIST;         /**< list of LP rows */

struct Lp                               /**< actual LP data */
{
   LPI*             lpi;                /**< LP solver interface */
   COL*             cols;               /**< array with actual LP columns in correct order */
   ROW*             rows;               /**< array with actual LP rows in correct order */
   int              colsize;            /**< available slots in cols vector */
   int              rowsize;            /**< available slots in rows vector */
   int              ncols;              /**< actual number of LP columns (number of used slots in cols vector) */
   int              nrows;              /**< actual number of LP rows (number of used slots in rows vector) */
};
typedef struct Lp LP;


extern
COL* SCIPcolCreate(                     /**< creates an LP column */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   int              len,                /**< number of nonzeros in the column */
   ROW**            row,                /**< array with rows of column entries */
   double*          val,                /**< array with coefficients of column entries */
   double           lb,                 /**< lower bound of variable */
   double           ub,                 /**< upper bound of variable */
   double           obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   );

extern
void SCIPcolFree(                       /**< frees an LP column */
   COL**            col,                /**< pointer to LP column */
   MEM*             mem                 /**< block memory buffers */
   );

extern
void SCIPcolCapture(                    /**< increases usage counter of LP column */
   COL*             col                 /**< LP column */
   );

extern
void SCIPcolRelease(                    /**< decreases usage counter of LP column, and frees memory if necessary */
   COL**            col,                /**< pointer to LP column */
   MEM*             mem                 /**< block memory buffers */
   );

extern
ROW* SCIProwCreate(                     /**< creates an LP row */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   double*          val,                /**< array with coefficients of row entries */
   Bool             equality,           /**< TRUE iff row is an equality, FALSE iff row is a less or equal inequality */
   double           rhs,                /**< right hand side of row */
   double           epsilon             /**< maximal normed violation of row */
   );

extern
void SCIProwFree(                       /**< frees an LP row */
   ROW**            row,                /**< pointer to LP row */
   MEM*             mem                 /**< block memory buffers */
   );

extern
void SCIProwCapture(                    /**< increases usage counter of LP row */
   ROW*             row                 /**< LP row */
   );

extern
void SCIProwRelease(                    /**< decreases usage counter of LP row, and frees memory if necessary */
   ROW**            row,                /**< pointer to LP row */
   MEM*             mem                 /**< block memory buffers */
   );


#endif
