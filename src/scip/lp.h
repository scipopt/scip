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


/** solution status after solving LP */
enum LPSolStat
{
   SCIP_LPSOLSTAT_OPTIMAL    = 0,       /**< LP was solved to optimality */
   SCIP_LPSOLSTAT_INFEASIBLE = 1,       /**< LP is primal infeasible */
   SCIP_LPSOLSTAT_UNBOUNDED  = 2,       /**< LP is primal unbounded */
   SCIP_LPSOLSTAT_ITERLIMIT  = 3,       /**< iteration limit was reached during optimization */
   SCIP_LPSOLSTAT_TIMELIMIT  = 4,       /**< time limit was reached during optimization */
   SCIP_LPSOLSTAT_ERROR      = 5        /**< an error occured during optimization */
};
typedef enum LPSolStat LPSOLSTAT;

typedef struct Col COL;                 /**< column of an LP */
typedef struct ColList COLLIST;         /**< list of LP columns */
typedef struct Row ROW;                 /**< row of an LP */
typedef struct RowList ROWLIST;         /**< list of LP rows */
typedef struct Lp LP;                   /**< actual LP data */



#include "def.h"
#include "mem.h"
#include "set.h"
#include "stat.h"
#include "lpi.h"
#include "var.h"



/** variable of the problem and corresponding LP column */
struct Col
{
   const VAR*       var;                /**< variable, this column represents; there cannot be a column without variable */
   ROW**            row;                /**< rows of column entries, that may have a nonzero dual solution value */
   Real*            val;                /**< coefficients of column entries */
   Real             primsol;            /**< primal solution value in LP, is 0 if col is not in LP */
   Real             redcost;            /**< reduced cost value in LP, or SCIP_INVALID if not yet calculated */
   int              index;              /**< consecutively numbered variable identifier */
   int              size;               /**< size of the row- and val-arrays */
   int              len;                /**< number of nonzeros in column */
   int              lpipos;             /**< column position number in LP solver, or -1 if not in LP solver */
   int              numpos;             /**< number of positive coefficients */
   int              numneg;             /**< number of negative coefficients */
   unsigned int     sorted:1;           /**< TRUE iff row indices are sorted in increasing order */
   unsigned int     lbchanged:1;        /**< TRUE iff lower bound changed, and data of LP solver has to be updated */
   unsigned int     ubchanged:1;        /**< TRUE iff upper bound changed, and data of LP solver has to be updated */
   unsigned int     coefchanged:1;      /**< TRUE iff the coefficient vector changed, and LP solver has to be updated */
   unsigned int     linked:1;           /**< TRUE iff column is inserted in corresponding row vectors */
   unsigned int     inLP:1;             /**< TRUE iff column is in actual LP */
};

/** row of the LP */
struct Row
{
   char*            name;               /**< name of the row */
   COL**            col;                /**< columns of row entries, that may have a nonzero primal solution value */
   Real*            val;                /**< coefficients of row entries */
   Real             rhs;                /**< right hand side of row */
   Real             epsilon;            /**< maximal normed violation of row */
   Real             sqrnorm;            /**< squared euclidean norm of row vector */
   Real             maxval;             /**< maximal absolute value of row vector */
   Real             dualsol;            /**< dual solution value in LP, is 0 if row is not in LP */
   Real             slack;              /**< slack value in LP, or SCIP_INVALID if not yet calculated */
   int              index;              /**< consecutively numbered row identifier */
   int              size;               /**< size of the col- and val-arrays */
   int              len;                /**< number of nonzeros in row */
   int              numuses;            /**< number of times, this row is referenced */
   int              lpipos;             /**< row position number in LP solver, or -1 if not in LP solver */
   int              minidx;             /**< minimal column index of row entries */
   int              maxidx;             /**< maximal column index of row entries */
   int              nummaxval;          /**< number of coefficients with absolute value equal to maxval */
   unsigned int     equality:1;         /**< TRUE iff row is an equality, FALSE iff row is a lower or equal inequality */
   unsigned int     sorted:1;           /**< TRUE iff column indices are sorted in increasing order */
   unsigned int     validminmaxidx:1;   /**< TRUE iff minimal and maximal column index is valid */
   unsigned int     coefchanged:1;      /**< TRUE iff the coefficient vector changed, and LP solver has to be updated */
   unsigned int     linked:1;           /**< TRUE iff row is inserted in corresponding column vectors */
   unsigned int     inLP:1;             /**< TRUE iff row is in actual LP */
};

/** actual LP data */
struct Lp
{
   LPI*             lpi;                /**< LP solver interface */
   COL**            lpicols;            /**< array with columns actually stored in the LP solver */
   ROW**            lpirows;            /**< array with rows actually stored in the LP solver */
   COL**            chgbds;             /**< array of columns with changed bounds not yet applied to the LP solver */
   COL**            cols;               /**< array with actual LP columns in correct order */
   ROW**            rows;               /**< array with actual LP rows in correct order */
   Real             objoffset;          /**< objective offset from bound shifting and fixing */
   Real             objval;             /**< objective value of LP, or SCIP_INVALID */
   int              lpicolssize;        /**< available slots in lpicols vector */
   int              nlpicols;           /**< number of columns in the LP solver */
   int              lpirowssize;        /**< available slots in lpirows vector */
   int              nlpirows;           /**< number of rows in the LP solver */
   int              lpifirstchgcol;     /**< first column of the LP which differs from the column in the LP solver */
   int              lpifirstchgrow;     /**< first row of the LP which differs from the row in the LP solver */
   int              colssize;           /**< available slots in cols vector */
   int              ncols;              /**< actual number of LP columns (number of used slots in cols vector) */
   int              rowssize;           /**< available slots in rows vector */
   int              nrows;              /**< actual number of LP rows (number of used slots in rows vector) */
   int              chgbdssize;         /**< available slots in chgbds vector */
   int              nchgbds;            /**< actual number of chgbds (number of used slots in chgbds vector) */
   int              firstnewcol;        /**< first column added at the actual node */
   int              firstnewrow;        /**< first row added at the actual node */
   unsigned int     flushed:1;          /**< TRUE iff all cached changes are applied to the LP solver */
   unsigned int     solved:1;           /**< TRUE iff current LP is solved */
};



extern
RETCODE SCIPcolCreate(                  /**< creates an LP column */
   COL**            col,                /**< pointer to column data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   const VAR*       var,                /**< variable, this column represents */
   int              len,                /**< number of nonzeros in the column */
   ROW**            row,                /**< array with rows of column entries */
   Real*            val                 /**< array with coefficients of column entries */
   );

extern
void SCIPcolFree(                       /**< frees an LP column */
   COL**            col,                /**< pointer to LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
void SCIPcolCalcRedcost(                /**< calculates the reduced costs of a column */
   COL*             col                 /**< LP column */
   );

extern
void SCIPcolSort(                       /**< sorts column entries by row index */
   COL* col                             /**< column to be sorted */
   );

extern
RETCODE SCIPcolAddCoeff(                /**< adds a previously non existing coefficient to an LP column */
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   );

extern
void SCIPcolDelCoeff(                   /**< deletes coefficient from column */
   COL*             col,                /**< column to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row                 /**< coefficient to be deleted */
   );

extern
RETCODE SCIPcolChgCoeff(                /**< changes or adds a coefficient to an LP column */
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   );

extern
RETCODE SCIPcolBoundChanged(            /**< notifies LP column, that its bounds were changed */
   COL*             col,                /**< LP column */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

extern
RETCODE SCIProwCreate(                  /**< creates an LP row */
   ROW**            row,                /**< pointer to LP row data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             rhs,                /**< right hand side of row */
   Real             epsilon,            /**< maximal normed violation of row */
   Bool             equality            /**< is row an equality? otherwise, it is a lower or equal inequality */
   );

extern
void SCIProwFree(                       /**< frees an LP row */
   ROW**            row,                /**< pointer to LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
void SCIProwCapture(                    /**< increases usage counter of LP row */
   ROW*             row                 /**< LP row */
   );

extern
void SCIProwRelease(                    /**< decreases usage counter of LP row, and frees memory if necessary */
   ROW**            row,                /**< pointer to LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
void SCIProwSort(                       /**< sorts row entries by column index */
   ROW* row                             /**< row to be sorted */
   );

extern
RETCODE SCIProwAddCoeff(                /**< adds a previously non existing coefficient to an LP row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
   Real             val                 /**< value of coefficient */
   );

extern
void SCIProwDelCoeff(                   /**< deletes coefficient from row */
   ROW*             row,                /**< row to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col                 /**< coefficient to be deleted */
   );

extern
RETCODE SCIProwChgCoeff(                /**< changes or adds a coefficient to an LP row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
   Real             val                 /**< value of coefficient */
   );

extern
void SCIPlpMarkSize(                    /**< remembers number of columns and rows to track the newly added ones */
   LP*              lp                  /**< actual LP data */
   );

extern
COL** SCIPlpGetNewcols(                 /**< get array with newly added columns after the last mark */
   const LP*        lp                  /**< actual LP data */
   );

extern
int SCIPlpGetNumNewcols(                /**< get number of newly added columns after the last mark */
   const LP*        lp                  /**< actual LP data */
   );

extern
ROW** SCIPlpGetNewrows(                 /**< get array with newly added rows after the last mark */
   const LP*        lp                  /**< actual LP data */
   );

extern
int SCIPlpGetNumNewrows(                /**< get number of newly added rows after the last mark */
   const LP*        lp                  /**< actual LP data */
   );

extern
RETCODE SCIPlpGetState(                 /**< stores LP state (like basis information) into LP state object */
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   LPSTATE**        lpstate             /**< pointer to LP state information (like basis information) */
   );

extern
RETCODE SCIPlpSetState(                 /**< loads LP state (like basis information) into solver */
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LPSTATE*         lpstate             /**< LP state information (like basis information) */
   );

extern
RETCODE SCIPlpAddCol(                   /**< adds a column to the LP */
   LP*              lp,                 /**< LP data */
   const SET*       set,                /**< global SCIP settings */
   COL*             col                 /**< LP column */
   );

extern
RETCODE SCIPlpAddRow(                   /**< adds a row to the LP */
   LP*              lp,                 /**< LP data */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row                 /**< LP row */
   );

extern
RETCODE SCIPlpShrinkCols(               /**< removes all columns after the given number of columns from the LP */
   LP*              lp,                 /**< LP data */
   int              newncols            /**< new number of columns in the LP */
   );

extern
RETCODE SCIPlpShrinkRows(               /**< removes all rows after the given number of rows from the LP */
   LP*              lp,                 /**< LP data */
   int              newnrows            /**< new number of rows in the LP */
   );

extern
RETCODE SCIPlpClear(                    /** removes all columns and rows from LP */
   LP*              lp                  /**< LP data */
   );

extern
RETCODE SCIPlpCreate(                   /**< creates empty LP data object */
   LP**             lp,                 /**< pointer to LP data object */
   const char*      name                /**< problem name */
   );

extern
RETCODE SCIPlpFree(                     /**< frees LP data object */
   LP**             lp                  /**< pointer to LP data object */
   );

extern
RETCODE SCIPlpSolvePrimal(              /**< solves the LP with the primal simplex algorithm */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   LPSOLSTAT*       lpsolstat           /**< pointer to store the LP solution status */
   );

extern
RETCODE SCIPlpSolveDual(                /**< solves the LP with the dual simplex algorithm */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   LPSOLSTAT*       lpsolstat           /**< pointer to store the LP solution status */
   );

extern
RETCODE SCIPlpGetSol(                   /**< stores the LP solution in the columns and rows */
   LP*              lp,                 /**< actual LP data */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr              /**< block memory buffers */
   );

#endif
