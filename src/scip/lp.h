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


/** variable type */
enum VarType
{
   SCIP_VARTYPE_BINARY    = 0,          /**< binary variable: $x \in \{0,1\}$ */
   SCIP_VARTYPE_INTEGER   = 1,          /**< integer variable: $x \in \{lb, \ldots, \ub\}$ */
   SCIP_VARTYPE_IMPLINT   = 2,          /**< implicit integer variable: continous variable, that is allways integral */
   SCIP_VARTYPE_CONTINOUS = 3           /**< continous variable: $x \in [lb,ub] */
};
typedef enum VarType VARTYPE;

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
#include "tree.h"



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
   int              lppos;              /**< column position number in LP, or -1 if not in LP */
   int              numpos;             /**< number of positive coefficients */
   int              numneg;             /**< number of negative coefficients */
   unsigned int     vartype:2;          /**< type of variable: binary, integer, implicit integer, continous */
   unsigned int     sorted:1;           /**< TRUE iff row indices are sorted in increasing order */
   unsigned int     lbchanged:1;        /**< TRUE iff lower bound changed in node switch, and LP has to be updated */
   unsigned int     ubchanged:1;        /**< TRUE iff upper bound changed in node switch, and LP has to be updated */
};

/** row of the LP */
struct Row
{
   char*            name;               /**< name of the row */
   COL**            col;                /**< columns of row entries */
   Real*            val;                /**< coefficients of row entries */
   Real             rhs;                /**< right hand side of row */
   Real             epsilon;            /**< maximal normed violation of row */
   Real             sqrnorm;            /**< squared euclidean norm of row vector */
   Real             maxval;             /**< maximal absolute value of row vector */
   int              index;              /**< consecutively numbered row identifier */
   int              size;               /**< size of the col- and val-arrays */
   int              len;                /**< number of nonzeros in row */
   int              numuses;            /**< number of times, this row is referenced */
   int              lppos;              /**< row position number in LP, or -1 if not in LP */
   int              minidx;             /**< minimal column index of row entries */
   int              maxidx;             /**< maximal column index of row entries */
   int              nummaxval;          /**< number of coefficients with absolute value equal to maxval */
   unsigned int     equality:1;         /**< TRUE iff row is an equality, FALSE iff row is a less or equal inequality */
   unsigned int     sorted:1;           /**< TRUE iff column indices are sorted in increasing order */
   unsigned int     validminmaxidx:1;   /**< TRUE iff minimal and maximal column index is valid */
};

/** actual LP data */
struct Lp
{
   LPI*             lpi;                /**< LP solver interface */
   COL**            cols;               /**< array with actual LP columns in correct order */
   ROW**            rows;               /**< array with actual LP rows in correct order */
   NODE**           children;           /**< array with children of the actual node */
   NODE**           siblings;           /**< array with siblings of the actual node */
   COL**            chgbds;             /**< array of columns with changed bounds not yet applied to the LP */
   COL**            addcols;            /**< array of columns to be added to the LP */
   ROW**            addrows;            /**< array of rows to be added to the LP */
   int              colssize;           /**< available slots in cols vector */
   int              ncols;              /**< actual number of LP columns (number of used slots in cols vector) */
   int              rowssize;           /**< available slots in rows vector */
   int              nrows;              /**< actual number of LP rows (number of used slots in rows vector) */
   int              childrensize;       /**< available slots in children vector */
   int              nchildren;          /**< actual number of children (number of used slots in children vector) */
   int              siblingssize;       /**< available slots in siblings vector */
   int              nsiblings;          /**< actual number of siblings (number of used slots in siblings vector) */
   int              chgbdssize;         /**< available slots in chgbds vector */
   int              nchgbds;            /**< actual number of chgbds (number of used slots in chgbds vector) */
   int              addcolssize;        /**< available slots in addcols vector */
   int              addcolscoefs;       /**< upper bound on number of coefficients in added cols */
   int              naddcols;           /**< actual number of addcols (number of used slots in addcols vector) */
   int              addrowssize;        /**< available slots in addrows vector */
   int              addrowscoefs;       /**< upper bound on number of coefficients in added rows */
   int              naddrows;           /**< actual number of addrows (number of used slots in addrows vector) */
   int              firstnewcol;        /**< first column added at the actual node */
   int              firstnewrow;        /**< first row added at the actual node */
   unsigned int     flushed:1;          /**< TRUE iff all cached changes are applied to the LP */
};



extern
DOMCHGDYN* SCIPdomchgdynCreate(         /**< creates a variable sized domain change data structure */
   MEM*             mem                 /**< block memory buffers */   
   );

extern
void SCIPdomchgdynFree(                 /**< frees a dynamically sized domain change data structure */
   DOMCHGDYN**      domchgdyn,          /**< pointer to dynamically sized domain change data structure */
   MEM*             mem                 /**< block memory buffers */   
   );

extern
RETCODE SCIPdomchgdynAddBoundchg(       /**< adds bound change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEM*             mem,                /**< block memory buffers */   
   const SET*       set,                /**< global SCIP settings */
   COL*             col,                /**< column to change the bounds for */
   Real             newbound,           /**< new value for bound */
   Real             oldbound,           /**< old value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

extern
RETCODE SCIPdomchgdynAddHolechg(        /**< adds hole change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEM*             mem,                /**< block memory buffers */   
   const SET*       set,                /**< global SCIP settings */
   HOLELIST**       ptr,                /**< changed list pointer */
   HOLELIST*        newlist,            /**< new value of list pointer */
   HOLELIST*        oldlist             /**< old value of list pointer */
   );

extern
DOMCHG* SCIPdomchgCreate(               /**< creates domain change data (fixed size) from dynamically sized data */
   MEM*             mem,                /**< block memory buffers */   
   const DOMCHGDYN* domchgdyn           /**< dynamically sized domain change data structure */
   );

extern
void SCIPdomchgFree(                    /**< frees domain change data */
   DOMCHG**         domchg,             /**< pointer to domain change */
   MEM*             mem                 /**< block memory buffers */
   );

extern
RETCODE SCIPlpApplyDomchg(              /**< applies domain change */
   LP*              lp,                 /**< actual LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   const DOMCHG*    domchg              /**< domain change to apply */
   );

extern
RETCODE SCIPlpUndoDomchg(               /**< undoes domain change */
   LP*              lp,                 /**< actual LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   const DOMCHG*    domchg              /**< domain change to remove */
   );

extern
RETCODE SCIPlpFlush(                    /**< applies all cached changes to the LP */
   LP*              lp,                 /**< actual LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   );

extern
COL* SCIPcolCreate(                     /**< creates an LP column */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   char*            name,               /**< name of column */
   int              len,                /**< number of nonzeros in the column */
   ROW**            row,                /**< array with rows of column entries */
   Real*            val,                /**< array with coefficients of column entries */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   );

extern
void SCIPcolFree(                       /**< frees an LP column */
   COL**            col,                /**< pointer to LP column */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   );

extern
void SCIPcolCapture(                    /**< increases usage counter of LP column */
   COL*             col                 /**< LP column */
   );

extern
void SCIPcolRelease(                    /**< decreases usage counter of LP column, and frees memory if necessary */
   COL**            col,                /**< pointer to LP column */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   );

extern
void SCIPcolSort(                       /**< sorts column entries by row index */
   COL* col                             /**< column to be sorted */
   );

extern
ROW* SCIProwCreate(                     /**< creates an LP row */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   char*            name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Bool             equality,           /**< TRUE iff row is an equality, FALSE iff row is a less or equal inequality */
   Real             rhs,                /**< right hand side of row */
   Real             epsilon             /**< maximal normed violation of row */
   );

extern
void SCIProwFree(                       /**< frees an LP row */
   ROW**            row,                /**< pointer to LP row */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   );

extern
void SCIProwCapture(                    /**< increases usage counter of LP row */
   ROW*             row                 /**< LP row */
   );

extern
void SCIProwRelease(                    /**< decreases usage counter of LP row, and frees memory if necessary */
   ROW**            row,                /**< pointer to LP row */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   );

extern
void SCIProwSort(                       /**< sorts row entries by column index */
   ROW* row                             /**< row to be sorted */
   );

extern
COL** SCIPlpGetNewcols(                 /**< get array with newly added columns */
   const LP*        lp                  /**< actual LP data */
   );

extern
int SCIPlpGetNumNewcols(                /**< get number of newly added columns */
   const LP*        lp                  /**< actual LP data */
   );

extern
ROW** SCIPlpGetNewrows(                 /**< get array with newly added rows */
   const LP*        lp                  /**< actual LP data */
   );

extern
int SCIPlpGetNumNewrows(                /**< get number of newly added rows */
   const LP*        lp                  /**< actual LP data */
   );

extern
RETCODE SCIPlpGetState(                 /**< stores LP state (like basis information) into LP state object */
   LP*              lp,                 /**< LP data */
   MEM*             mem,                /**< block memory buffers */
   LPSTATE**        lpstate             /**< pointer to LP state information (like basis information) */
   );

extern
RETCODE SCIPlpSetState(                 /**< loads LP state (like basis information) into solver */
   LP*              lp,                 /**< LP data */
   MEM*             mem,                /**< block memory buffers */
   LPSTATE*         lpstate             /**< LP state information (like basis information) */
   );

extern
RETCODE SCIPlpAddCol(                   /**< adds a column to the LP */
   LP*              lp,                 /**< LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   COL*             col                 /**< LP column */
   );

extern
RETCODE SCIPlpAddRow(                   /**< adds a row to the LP */
   LP*              lp,                 /**< LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row                 /**< LP row */
   );

extern
RETCODE SCIPlpShrinkCols(               /**< removes all columns after given column number from LP */
   LP*              lp,                 /**< LP data */
   int              lastcol             /**< last column number to remain in the LP */
   );

extern
RETCODE SCIPlpShrinkRows(               /**< removes all rows after given rowumn number from LP */
   LP*              lp,                 /**< LP data */
   int              lastrow             /**< last row number to remain in the LP */
   );

#endif
