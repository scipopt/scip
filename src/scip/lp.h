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

/**
 * In LP management, we have to differ between the actual LP and the LP stored in the LP solver.
 * All LP methods affect the actual LP only. Before solving the actual LP with the LP solver,
 * the LP solvers data has to be updated to the actual LP with a call to SCIPlpFlush().
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __LP_H__
#define __LP_H__


/** variable type */
enum Coltype
{
   SCIP_COLTYPE_BINARY    = 0,          /**< binary variable: $x \in \{0,1\}$ */
   SCIP_COLTYPE_INTEGER   = 1,          /**< integer variable: $x \in \{lb, \ldots, \ub\}$ */
   SCIP_COLTYPE_IMPLINT   = 2,          /**< implicit integer variable: continous variable, that is allways integral */
   SCIP_COLTYPE_CONTINOUS = 3           /**< continous variable: $x \in [lb,ub] */
};
typedef enum Coltype COLTYPE;

/** row type */
enum RowType
{
   SCIP_ROWTYPE_LESSEQUAL    = 0,       /**< less or equal row: $a^Tx \le \text{rhs}$ */
   SCIP_ROWTYPE_EQUAL        = 1,       /**< equality row: $a^Tx = \text{rhs}$ */
   SCIP_ROWTYPE_GREATEREQUAL = 2,       /**< greater or equal row: $a^Tx \ge \text{rhs}$ */
   SCIP_ROWTYPE_RANGE        = 3        /**< ranged row: $\text{lhs} \le a^Tx \le \text{rhs}$ */
};
typedef enum RowType ROWTYPE;

/** type of bound: lower or upper bound */
enum BoundType
{
   SCIP_BOUNDTYPE_LOWER = 0,            /**< lower bound */
   SCIP_BOUNDTYPE_UPPER = 1             /**< upper bound */
};
typedef enum BoundType BOUNDTYPE;

typedef struct Dom DOM;                 /**< datastructures for storing domains of variables */
typedef struct Hole HOLE;               /**< hole in a domain of an integer variable */
typedef struct Holelist HOLELIST;       /**< list of holes in a domain of an integer variable */
typedef struct DomChg DOMCHG;           /**< changes in domains of variables (fixed sized arrays) */
typedef struct DomChgDyn DOMCHGDYN;     /**< changes in domains of variables (dynamically sized arrays) */
typedef struct BoundChg BOUNDCHG;       /**< changes in bounds of variables */
typedef struct HoleChg HOLECHG;         /**< changes in holelist of variables */
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



/** domain of a variable */
struct Dom
{
   HOLELIST*        holelist;           /**< list of holes (only for the integer variables) */
   Real             lb;                 /**< lower bounds of variables */
   Real             ub;                 /**< upper bounds of variables */
};

/** variable of the problem and corresponding LP column */
struct Col
{
   char*            name;               /**< name of the column */
   ROW**            row;                /**< rows of column entries */
   Real*            val;                /**< coefficients of column entries */
   DOM              dom;                /**< domain of variable in actual subproblem */
   Real             problb;             /**< lower bound of variable in main problem */
   Real             probub;             /**< upper bound of variable in main problem */
   Real             obj;                /**< objective function value of variable */
   int              index;              /**< consecutively numbered variable identifier */
   int              size;               /**< size of the row- and val-arrays */
   int              len;                /**< number of nonzeros in column */
   int              numuses;            /**< number of times, this column is referenced */
   int              lpipos;             /**< column position number in LP solver, or -1 if not in LP solver */
   int              numpos;             /**< number of positive coefficients */
   int              numneg;             /**< number of negative coefficients */
   unsigned int     coltype:2;          /**< type of variable: binary, integer, implicit integer, continous */
   unsigned int     sorted:1;           /**< TRUE iff row indices are sorted in increasing order */
   unsigned int     lbchanged:1;        /**< TRUE iff lower bound changed, and data of LP solver has to be updated */
   unsigned int     ubchanged:1;        /**< TRUE iff upper bound changed, and data of LP solver has to be updated */
   unsigned int     coefchanged:1;      /**< TRUE iff the coefficient vector changed, and LP solver has to be updated */
};

/** row of the LP */
struct Row
{
   char*            name;               /**< name of the row */
   COL**            col;                /**< columns of row entries */
   Real*            val;                /**< coefficients of row entries */
   Real             rhs;                /**< right hand side of row */
   Real             lhs;                /**< left hand side of row (for ranged rows) */
   Real             epsilon;            /**< maximal normed violation of row */
   Real             sqrnorm;            /**< squared euclidean norm of row vector */
   Real             maxval;             /**< maximal absolute value of row vector */
   int              index;              /**< consecutively numbered row identifier */
   int              size;               /**< size of the col- and val-arrays */
   int              len;                /**< number of nonzeros in row */
   int              numuses;            /**< number of times, this row is referenced */
   int              lpipos;             /**< row position number in LP solver, or -1 if not in LP solver */
   int              minidx;             /**< minimal column index of row entries */
   int              maxidx;             /**< maximal column index of row entries */
   int              nummaxval;          /**< number of coefficients with absolute value equal to maxval */
   unsigned int     rowtype:2;          /**< type of row: lessequal, equal, greaterequal, range */
   unsigned int     sorted:1;           /**< TRUE iff column indices are sorted in increasing order */
   unsigned int     validminmaxidx:1;   /**< TRUE iff minimal and maximal column index is valid */
   unsigned int     coefchanged:1;      /**< TRUE iff the coefficient vector changed, and LP solver has to be updated */
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
RETCODE SCIPdomchgdynCreate(            /**< creates a dynamically sized domain change data structure */
   DOMCHGDYN**      domchgdyn           /**< pointer to dynamically sized domain change data structure */
   );

extern
void SCIPdomchgdynFree(                 /**< frees a dynamically sized domain change data structure */
   DOMCHGDYN**      domchgdyn           /**< pointer to dynamically sized domain change data structure */
   );

extern
RETCODE SCIPdomchgdynCopy(              /**< copies data from fixed size domain change into dynamically sized one */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   const SET*       set,                /**< global SCIP settings */
   DOMCHG*          domchg              /**< static domain change */
   );

extern
RETCODE SCIPdomchgdynAddBoundchg(       /**< adds bound change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   const SET*       set,                /**< global SCIP settings */
   COL*             col,                /**< column to change the bounds for */
   Real             newbound,           /**< new value for bound */
   Real             oldbound,           /**< old value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

extern
RETCODE SCIPdomchgdynAddHolechg(        /**< adds hole change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   const SET*       set,                /**< global SCIP settings */
   HOLELIST**       ptr,                /**< changed list pointer */
   HOLELIST*        newlist,            /**< new value of list pointer */
   HOLELIST*        oldlist             /**< old value of list pointer */
   );

extern
RETCODE SCIPdomchgCreate(               /**< creates domain change data (fixed size) from dynamically sized data */
   DOMCHG**         domchg,             /**< pointer to fixed size domain change data */
   MEMHDR*          memhdr,             /**< block memory */
   const DOMCHGDYN* domchgdyn           /**< dynamically sized domain change data structure */
   );

extern
void SCIPdomchgFree(                    /**< frees fixed size domain change data */
   DOMCHG**         domchg,             /**< pointer to domain change */
   MEMHDR*          memhdr              /**< block memory */
   );

extern
RETCODE SCIPlpApplyDomchg(              /**< applies domain change */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   const DOMCHG*    domchg              /**< domain change to apply */
   );

extern
RETCODE SCIPlpUndoDomchg(               /**< undoes domain change */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   const DOMCHG*    domchg              /**< domain change to remove */
   );

extern
RETCODE SCIPlpFlush(                    /**< applies all cached changes to the LP solver */
   LP*              lp,                 /**< actual LP data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPcolCreate(                  /**< creates an LP column */
   COL**            col,                /**< pointer to column data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of column */
   int              len,                /**< number of nonzeros in the column */
   ROW**            row,                /**< array with rows of column entries */
   Real*            val,                /**< array with coefficients of column entries */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   COLTYPE          coltype             /**< type of variable */
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
   Real             lhs,                /**< left hand side of row (for ranged rows) */
   Real             epsilon,            /**< maximal normed violation of row */
   ROWTYPE          rowtype             /**< type of row */
   );

extern
void SCIPcolFree(                       /**< frees an LP column */
   COL**            col,                /**< pointer to LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
void SCIProwFree(                       /**< frees an LP row */
   ROW**            row,                /**< pointer to LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
void SCIPcolCapture(                    /**< increases usage counter of LP column */
   COL*             col                 /**< LP column */
   );

extern
void SCIProwCapture(                    /**< increases usage counter of LP row */
   ROW*             row                 /**< LP row */
   );

extern
void SCIPcolRelease(                    /**< decreases usage counter of LP column, and frees memory if necessary */
   COL**            col,                /**< pointer to LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
void SCIProwRelease(                    /**< decreases usage counter of LP row, and frees memory if necessary */
   ROW**            row,                /**< pointer to LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
void SCIPcolSort(                       /**< sorts column entries by row index */
   COL* col                             /**< column to be sorted */
   );

extern
void SCIProwSort(                       /**< sorts row entries by column index */
   ROW* row                             /**< row to be sorted */
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
RETCODE SCIProwAddCoeff(                /**< adds a previously non existing coefficient to an LP row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
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
void SCIProwDelCoeff(                   /**< deletes coefficient from row */
   ROW*             row,                /**< row to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col                 /**< coefficient to be deleted */
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

#endif
