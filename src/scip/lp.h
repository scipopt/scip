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

/** In SCIP, the LP is defined as follows:
 *
 *   min       obj * x
 *      lhs <=   A * x <= rhs
 *      lb  <=       x <= ub
 *
 *  The row activities are defined as 
 *     activity = A * x
 *  and must therefore be in the range of [lhs,rhs].
 *
 *  The reduced costs are defined as
 *     redcost = obj - A^T * y
 *  and must be   nonnegative, if the corresponding lb is nonnegative,
 *                zero,        if the corresponging lb is negative.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __LP_H__
#define __LP_H__


#include <stdio.h>


/** solution status after solving LP */
enum LPSolStat
{
   SCIP_LPSOLSTAT_NOTSOLVED  = 0,       /**< LP was not solved, no solution exists */
   SCIP_LPSOLSTAT_OPTIMAL    = 1,       /**< LP was solved to optimality */
   SCIP_LPSOLSTAT_INFEASIBLE = 2,       /**< LP is primal infeasible */
   SCIP_LPSOLSTAT_UNBOUNDED  = 3,       /**< LP is primal unbounded */
   SCIP_LPSOLSTAT_OBJLIMIT   = 4,       /**< objective limit was reached during optimization */
   SCIP_LPSOLSTAT_ITERLIMIT  = 5,       /**< iteration limit was reached during optimization */
   SCIP_LPSOLSTAT_TIMELIMIT  = 6,       /**< time limit was reached during optimization */
   SCIP_LPSOLSTAT_ERROR      = 7        /**< an error occured during optimization */
};
typedef enum LPSolStat LPSOLSTAT;

/** type of variable bound: lower or upper bound */
enum BoundType
{
   SCIP_BOUNDTYPE_LOWER = 0,            /**< lower bound */
   SCIP_BOUNDTYPE_UPPER = 1             /**< upper bound */
};
typedef enum BoundType BOUNDTYPE;

/** type of row side: left hand or right hand side */
enum SideType
{
   SCIP_SIDETYPE_LEFT  = 0,             /**< left hand side */
   SCIP_SIDETYPE_RIGHT = 1              /**< right hand side */
};
typedef enum SideType SIDETYPE;

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
   VAR*             var;                /**< variable, this column represents; there cannot be a column without variable */
   ROW**            row;                /**< rows of column entries, that may have a nonzero dual solution value */
   Real*            val;                /**< coefficients of column entries */
   int*             linkpos;            /**< position of col in col vector of the row, or -1 if not yet linked */
   Real             primsol;            /**< primal solution value in LP, is 0 if col is not in LP */
   Real             redcost;            /**< reduced cost value in LP, or SCIP_INVALID if not yet calculated */
   Real             farkas;             /**< value in dual farkas infeasibility proof */
   Real             strongdown;         /**< strong branching information for downwards branching */
   Real             strongup;           /**< strong branching information for upwards branching */
   int              index;              /**< consecutively numbered column identifier */
   int              size;               /**< size of the row- and val-arrays */
   int              len;                /**< number of nonzeros in column */
   int              nunlinked;          /**< number of column entries, where the rows don't know about the column */
   int              lpipos;             /**< column position number in LP solver, or -1 if not in LP solver */
   int              validredcostlp;     /**< lp number for which reduced cost value is valid */
   int              validfarkaslp;      /**< lp number for which farkas value is valid */
   int              strongitlim;        /**< strong branching iteration limit used to get strongdown and strongup, or -1 */
   unsigned int     sorted:1;           /**< TRUE iff row indices are sorted in increasing order */
   unsigned int     lbchanged:1;        /**< TRUE iff lower bound changed, and data of LP solver has to be updated */
   unsigned int     ubchanged:1;        /**< TRUE iff upper bound changed, and data of LP solver has to be updated */
   unsigned int     coefchanged:1;      /**< TRUE iff the coefficient vector changed, and LP solver has to be updated */
   unsigned int     inlp:1;             /**< TRUE iff column is in actual LP */
};

/** row of the LP */
struct Row
{
   char*            name;               /**< name of the row */
   COL**            col;                /**< columns of row entries, that may have a nonzero primal solution value */
   Real*            val;                /**< coefficients of row entries */
   int*             linkpos;            /**< position of row in row vector of the column, or -1 if not yet linked */
   Real             constant;           /**< constant shift c in row lhs <= ax + c <= rhs */
   Real             lhs;                /**< left hand side of row */
   Real             rhs;                /**< right hand side of row */
   Real             sqrnorm;            /**< squared euclidean norm of row vector */
   Real             maxval;             /**< maximal absolute value of row vector, only valid if nummaxval > 0 */
   Real             dualsol;            /**< dual solution value in LP, is 0 if row is not in LP */
   Real             activity;           /**< row activity value in LP, or SCIP_INVALID if not yet calculated */
   Real             dualfarkas;         /**< multiplier value in dual farkas infeasibility proof */
   Real             pseudoactivity;     /**< row activity value in pseudo solution, or SCIP_INVALID if not yet calculated */
   Real             minactivity;        /**< minimal activity value w.r.t. the column's bounds, or SCIP_INVALID,
                                         *   ignoring the coefficients contributing with infinite value */
   Real             maxactivity;        /**< maximal activity value w.r.t. the column's bounds, or SCIP_INVALID,
                                         *   ignoring the coefficients contributing with infinite value */
   int              minactivityinf;     /**< number of coefficients contributing with infinite value to minactivity */
   int              maxactivityinf;     /**< number of coefficients contributing with infinite value to maxactivity */
   int              index;              /**< consecutively numbered row identifier */
   int              size;               /**< size of the col- and val-arrays */
   int              len;                /**< number of nonzeros in row */
   int              nunlinked;          /**< number of row entries, where the columns don't know about the row */
   int              nuses;              /**< number of times, this row is referenced */
   int              lpipos;             /**< row position number in LP solver, or -1 if not in LP solver */
   int              minidx;             /**< minimal column index of row entries */
   int              maxidx;             /**< maximal column index of row entries */
   int              nummaxval;          /**< number of coefs with absolute value equal to maxval, zero if maxval invalid */
   int              validactivitylp;    /**< lp number for which activity value is valid */
   unsigned int     sorted:1;           /**< are column indices sorted in increasing order? */
   unsigned int     validpsactivity:1;  /**< is the pseudo activity valid? */
   unsigned int     validactivitybds:1; /**< are the activity bounds minactivity/maxactivity valid? */
   unsigned int     validminmaxidx:1;   /**< are minimal and maximal column index valid? */
   unsigned int     lhschanged:1;       /**< was left hand side changed, and has data of LP solver to be updated? */
   unsigned int     rhschanged:1;       /**< was right hand side changed, and has data of LP solver to be updated? */
   unsigned int     coefchanged:1;      /**< was the coefficient vector changed, and has LP solver to be updated? */
   unsigned int     inlp:1;             /**< is row in actual LP? */
   unsigned int     modifiable:1;       /**< is row modifiable during node processing (subject to column generation)? */
   unsigned int     nlocks:23;          /**< number of sealed locks of an unmodifiable row */
};

/** actual LP data */
struct Lp
{
   LPI*             lpi;                /**< LP solver interface */
   COL**            lpicols;            /**< array with columns actually stored in the LP solver */
   ROW**            lpirows;            /**< array with rows actually stored in the LP solver */
   COL**            chgcols;            /**< array of changed columns not yet applied to the LP solver */
   ROW**            chgrows;            /**< array of changed rows not yet applied to the LP solver */
   COL**            cols;               /**< array with actual LP columns in correct order */
   ROW**            rows;               /**< array with actual LP rows in correct order */
   LPSOLSTAT        lpsolstat;          /**< solution status of last LP solution */
   Real             objval;             /**< objective value of LP, or SCIP_INVALID */
   int              lpicolssize;        /**< available slots in lpicols vector */
   int              nlpicols;           /**< number of columns in the LP solver */
   int              lpirowssize;        /**< available slots in lpirows vector */
   int              nlpirows;           /**< number of rows in the LP solver */
   int              lpifirstchgcol;     /**< first column of the LP which differs from the column in the LP solver */
   int              lpifirstchgrow;     /**< first row of the LP which differs from the row in the LP solver */
   int              chgcolssize;        /**< available slots in chgcols vector */
   int              nchgcols;           /**< actual number of chgcols (number of used slots in chgcols vector) */
   int              chgrowssize;        /**< available slots in chgrows vector */
   int              nchgrows;           /**< actual number of chgrows (number of used slots in chgrows vector) */
   int              colssize;           /**< available slots in cols vector */
   int              ncols;              /**< actual number of LP columns (number of used slots in cols vector) */
   int              rowssize;           /**< available slots in rows vector */
   int              nrows;              /**< actual number of LP rows (number of used slots in rows vector) */
   int              firstnewcol;        /**< first column added at the active node */
   int              firstnewrow;        /**< first row added at the active node */
   unsigned int     flushed:1;          /**< are all cached changes applied to the LP solver? */
   unsigned int     solved:1;           /**< is current LP solved? */
   unsigned int     primalfeasible:1;   /**< is actual LP basis primal feasible? */
   unsigned int     dualfeasible:1;     /**< is actual LP basis dual feasible? */
};


/*
 * Column methods
 */

extern
RETCODE SCIPcolCreate(                  /**< creates an LP column */
   COL**            col,                /**< pointer to column data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   VAR*             var,                /**< variable, this column represents */
   int              len,                /**< number of nonzeros in the column */
   ROW**            row,                /**< array with rows of column entries */
   Real*            val                 /**< array with coefficients of column entries */
   );

extern
RETCODE SCIPcolFree(                    /**< frees an LP column */
   COL**            col,                /**< pointer to LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
void SCIPcolSort(                       /**< sorts column entries by row index */
   COL*             col                 /**< column to be sorted */
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
RETCODE SCIPcolDelCoeff(                /**< deletes coefficient from column */
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
RETCODE SCIPcolIncCoeff(                /**< increases value of an existing or nonexisting coefficient in an LP column */
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< LP row */
   Real             incval              /**< value to add to the coefficient */
   );

extern
RETCODE SCIPcolBoundChanged(            /**< notifies LP, that the bounds of a column were changed */
   COL*             col,                /**< LP column that changed */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   Real             oldbound,           /**< old bound value */
   Real             newbound            /**< new bound value */
   );

extern
VAR* SCIPcolGetVar(                     /**< gets variable this column represents */
   COL*             col                 /**< LP column */
   );

extern
Bool SCIPcolIsInLP(                     /**< returns TRUE iff column is member of actual LP */
   COL*             col                 /**< LP column */
   );

extern
Real SCIPcolGetPrimsol(                 /**< gets the primal LP solution of a column */
   COL*             col                 /**< LP column */
   );

extern
Real SCIPcolGetRedcost(                 /**< gets the reduced costs of a column in last LP or after recalculation */
   COL*             col,                /**< LP column */
   STAT*            stat                /**< problem statistics */
   );

extern
Real SCIPcolGetFeasibility(             /**< gets the feasibility of a column in last LP or after recalculation */
   COL*             col,                /**< LP column */
   STAT*            stat                /**< problem statistics */
   );

extern
Real SCIPcolGetFarkas(                  /**< gets the farkas value of a column in last LP (which must be infeasible) */
   COL*             col,                /**< LP column */
   STAT*            stat                /**< problem statistics */
   );

extern
RETCODE SCIPcolGetStrongbranch(         /**< gets strong branching information on a column variable */
   COL*             col,                /**< LP column */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   Real             upperbound,         /**< actual global upper bound */
   int              itlim,              /**< iteration limit for strong branchings */
   Real*            down,               /**< stores dual bound after branching column down */
   Real*            up                  /**< stores dual bound after branching column up */
   );

extern
int SCIPcolGetNNonz(                    /**< get number of nonzero entries in column vector */
   COL*             col                 /**< LP column */
   );

extern
ROW** SCIPcolGetRows(                   /**< gets array with rows of nonzero entries */
   COL*             col                 /**< LP column */
   );

extern
Real* SCIPcolGetVals(                   /**< gets array with coefficients of nonzero entries */
   COL*             col                 /**< LP column */
   );

extern
void SCIPcolPrint(                      /**< output column to file stream */
   COL*             col,                /**< LP column */
   const SET*       set,                /**< global SCIP settings */
   FILE*            file                /**< output file (or NULL for standard output) */
   );


/*
 * Row methods
 */

extern
RETCODE SCIProwCreate(                  /**< creates and captures an LP row */
   ROW**            row,                /**< pointer to LP row data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   const char*      name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             modifiable          /**< is row modifiable during node processing (subject to column generation)? */
   );

extern
RETCODE SCIProwFree(                    /**< frees an LP row */
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
RETCODE SCIProwRelease(                 /**< decreases usage counter of LP row, and frees memory if necessary */
   ROW**            row,                /**< pointer to LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIProwLock(                    /**< locks an unmodifiable row, which forbids further changes */
   ROW*             row                 /**< LP row */
   );

extern
RETCODE SCIProwUnlock(                  /**< unlocks a lock of a row; a row with no sealed lock may be modified */
   ROW*             row                 /**< LP row */
   );

extern
void SCIProwSort(                       /**< sorts row entries by column index */
   ROW*             row                 /**< row to be sorted */
   );

extern
RETCODE SCIProwForbidRounding(          /**< forbids roundings of variables in row that may violate row */
   ROW*             row,                /**< LP row */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIProwAllowRounding(           /**< allows roundings of variables in row that may violate row */
   ROW*             row,                /**< LP row */
   const SET*       set                 /**< global SCIP settings */
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
RETCODE SCIProwDelCoeff(                /**< deletes coefficient from row */
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
RETCODE SCIProwIncCoeff(                /**< increases value of an existing or nonexisting coefficient in an LP column */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
   Real             incval              /**< value to add to the coefficient */
   );

extern
RETCODE SCIProwAddConst(                /**< add constant value to a row, i.e. subtract value from lhs and rhs */
   ROW*             row,                /**< LP row */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             constant            /**< constant value to add to the row */
   );

extern
RETCODE SCIProwChgLhs(                  /**< changes left hand side of LP row */
   ROW*             row,                /**< LP row */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             lhs                 /**< new left hand side */
   );

extern
RETCODE SCIProwChgRhs(                  /**< changes right hand side of LP row */
   ROW*             row,                /**< LP row */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             rhs                 /**< new right hand side */
   );

extern
Real SCIProwGetLhs(                     /**< returns the left hand side of the row */
   ROW*             row                 /**< LP row */
   );

extern
Real SCIProwGetRhs(                     /**< returns the right hand side of the row */
   ROW*             row                 /**< LP row */
   );

extern
const char* SCIProwGetName(             /**< returns the name of the row */
   ROW*             row                 /**< LP row */
   );

extern
Bool SCIProwIsInLP(                     /**< returns TRUE iff row is member of actual LP */
   ROW*             row                 /**< LP row */
   );

extern
Real SCIProwGetActivity(                /**< returns the activity of a row in the last LP or after recalculation */
   ROW*             row,                /**< LP row */
   STAT*            stat                /**< problem statistics */
   );

extern
Real SCIProwGetFeasibility(             /**< returns the feasibility of a row in the last solution or after recalc */
   ROW*             row,                /**< LP row */
   STAT*            stat                /**< problem statistics */
   );

extern
RETCODE SCIProwGetPseudoActivity(       /**< returns the pseudo activity of a row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real*            pseudoactivity      /**< pointer to store the pseudo activity */
   );

extern
RETCODE SCIProwGetPseudoFeasibility(    /**< returns the feasibility of a row in the actual pseudo solution */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real*            pseudofeasibility   /**< pointer to store the pseudo feasibility */
   );

extern
RETCODE SCIProwGetActivityBounds(       /**< returns the minimal and maximal activity of a row w.r.t. the column's bounds */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real*            minactivity,        /**< pointer to store the minimal activity, or NULL */
   Real*            maxactivity         /**< pointer to store the maximal activity, or NULL */
   );

extern
RETCODE SCIProwGetActivityResiduals(    /**< gets activity bounds for row after setting variable to zero */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   VAR*             var,                /**< variable to calculate activity residual for */
   Real             val,                /**< coefficient value of variable in linear constraint */
   Real*            minresactivity,     /**< pointer to store the minimal residual activity */
   Real*            maxresactivity      /**< pointer to store the maximal residual activity */
   );

extern
RETCODE SCIProwInvalidActivityBounds(   /**< invalidates activity bounds, such that they are recalculated in next get */
   ROW*             row                 /**< LP row */
   );

extern
int SCIProwGetNNonz(                    /**< get number of nonzero entries in row vector */
   ROW*             row                 /**< LP row */
   );

extern
COL** SCIProwGetCols(                   /**< gets array with columns of nonzero entries */
   ROW*             row                 /**< LP row */
   );

extern
Real* SCIProwGetVals(                   /**< gets array with coefficients of nonzero entries */
   ROW*             row                 /**< LP row */
   );

extern
Real SCIProwGetConstant(                /**< gets constant shift of row */
   ROW*             row                 /**< LP row */
   );

extern
Real SCIProwGetNorm(                    /**< get euclidean norm of row vector */
   ROW*             row                 /**< LP row */
   );

extern
Real SCIProwGetMaxval(                  /**< gets maximal absolute value of row vector coefficients */
   ROW*             row,                /**< LP row */
   const SET*       set                 /**< global SCIP settings */
   );

extern
void SCIProwPrint(                      /**< output row to file stream */
   ROW*             row,                /**< LP row */
   FILE*            file                /**< output file (or NULL for standard output) */
   );


/*
 * LP methods
 */

extern
RETCODE SCIPlpCreate(                   /**< creates empty LP data object */
   LP**             lp,                 /**< pointer to LP data object */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< problem name */
   );

extern
RETCODE SCIPlpFree(                     /**< frees LP data object */
   LP**             lp,                 /**< pointer to LP data object */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPlpAddCol(                   /**< adds a column to the LP and captures the variable */
   LP*              lp,                 /**< LP data */
   const SET*       set,                /**< global SCIP settings */
   COL*             col                 /**< LP column */
   );

extern
RETCODE SCIPlpAddRow(                   /**< adds a row to the LP and captures it */
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
RETCODE SCIPlpShrinkRows(               /**< removes and releases all rows after the given number of rows from the LP */
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   int              newnrows            /**< new number of rows in the LP */
   );

extern
RETCODE SCIPlpClear(                    /**< removes all columns and rows from LP, releases all rows */
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
   LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   );

extern
RETCODE SCIPlpSetState(                 /**< loads LP state (like basis information) into solver */
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LPISTATE*        lpistate            /**< LP state information (like basis information) */
   );

extern
RETCODE SCIPlpSetFeastol(               /**< sets the feasibility tolerance of the LP solver */
   LP*              lp,                 /**< actual LP data */
   Real             feastol             /**< new feasibility tolerance */
   );

extern
RETCODE SCIPlpSetUpperbound(            /**< sets the upper objective limit of the LP solver */
   LP*              lp,                 /**< actual LP data */
   Real             upperbound          /**< new upper objective limit */
   );

extern
RETCODE SCIPlpSolvePrimal(              /**< solves the LP with the primal simplex algorithm */
   LP*              lp,                 /**< actual LP data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   );

extern
RETCODE SCIPlpSolveDual(                /**< solves the LP with the dual simplex algorithm */
   LP*              lp,                 /**< actual LP data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   );

extern
LPSOLSTAT SCIPlpGetSolstat(             /**< gets solution status of last solve call */
   LP*              lp                  /**< actual LP data */
   );

extern
Real SCIPlpGetObjval(                   /**< gets objective value of last solution */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPlpGetSol(                   /**< stores the LP solution in the columns and rows */
   LP*              lp,                 /**< actual LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   );

extern
RETCODE SCIPlpGetUnboundedSol(          /**< stores LP solution with infinite objective value in the columns and rows */
   LP*              lp,                 /**< actual LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   );

extern
RETCODE SCIPlpGetDualfarkas(            /**< stores the dual farkas multipliers for infeasibility proof in rows */
   LP*              lp,                 /**< actual LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPlpGetIterations(            /**< get number of iterations used in last LP solve */
   LP*              lp,                 /**< actual LP data */
   int*             iterations          /**< pointer to store the iteration count */
   );


#endif
