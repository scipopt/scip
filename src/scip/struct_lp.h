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
#pragma ident "@(#) $Id: struct_lp.h,v 1.4 2004/01/15 09:12:15 bzfpfend Exp $"

/**@file   struct_lp.h
 * @brief  datastructures for LP management
 * @author Tobias Achterberg
 *
 *  In SCIP, the LP is defined as follows:
 *
 *   min       obj * x
 *      lhs <=   A * x + const <= rhs
 *      lb  <=       x         <= ub
 *
 *  The row activities are defined as 
 *     activity = A * x + const
 *  and must therefore be in the range of [lhs,rhs].
 *
 *  The reduced costs are defined as
 *     redcost = obj - A^T * y
 *  and must be   nonnegative, if the corresponding lb is nonnegative,
 *                zero,        if the corresponging lb is negative.
 *
 *  The main datastructures for storing an LP are the rows and the columns.
 *  A row can live on its own (if it was created by a separator), or as LP
 *  relaxation of a constraint. Thus, it has a nuses-counter, and is
 *  deleted, if not needed any more.
 *  A column cannot live on its own. It is always connected to a problem
 *  variable. Because pricing is always problem specific, it cannot create
 *  LP columns without introducing new variables. Thus, each column is
 *  connected to exactly one variable, and is deleted, if the variable
 *  is deleted.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_LP_H__
#define __STRUCT_LP_H__


#include "def.h"
#include "type_lpi.h"
#include "type_lp.h"
#include "type_var.h"



/** variable of the problem and corresponding LP column */
struct Col
{
   VAR*             var;                /**< variable, this column represents; there cannot be a column without variable */
   ROW**            rows;               /**< rows of column entries, that may have a nonzero dual solution value */
   Real*            vals;               /**< coefficients of column entries */
   int*             linkpos;            /**< position of col in col vector of the row, or -1 if not yet linked */
   Real             obj;                /**< current objective value of column in LP */
   Real             lb;                 /**< current lower bound of column in LP */
   Real             ub;                 /**< current upper bound of column in LP */
   Real             primsol;            /**< primal solution value in LP, is 0 if col is not in LP */
   Real             redcost;            /**< reduced cost value in LP, or SCIP_INVALID if not yet calculated */
   Real             farkas;             /**< value in dual farkas infeasibility proof */
   Real             strongbranchdown;   /**< strong branching information for downwards branching */
   Real             strongbranchup;     /**< strong branching information for upwards branching */
   Real             strongbranchsolval; /**< LP solution value of column at last strong branching call */
   Longint          strongbranchnode;   /**< node number of the last strong branching call on this column */
   Longint          obsoletenode;       /**< last node where this column was removed due to aging */
   int              index;              /**< consecutively numbered column identifier */
   int              size;               /**< size of the row- and val-arrays */
   int              len;                /**< number of nonzeros in column */
   int              nunlinked;          /**< number of column entries, where the rows don't know about the column */
   int              lppos;              /**< column position number in actual LP, or -1 if not in actual LP */
   int              lpipos;             /**< column position number in LP solver, or -1 if not in LP solver */
   int              validredcostlp;     /**< LP number for which reduced cost value is valid */
   int              validfarkaslp;      /**< LP number for which farkas value is valid */
   int              validstrongbranchlp;/**< LP number for which strong branching values are valid */
   int              strongbranchitlim;  /**< strong branching iteration limit used to get strongbranch values, or -1 */
   int              age;                /**< number of successive times this variable was in LP and was 0.0 in solution */
   int              var_probindex;      /**< copy of var->probindex for avoiding expensive dereferencing */
   unsigned int     sorted:1;           /**< TRUE iff row indices are sorted in increasing order */
   unsigned int     objchanged:1;       /**< TRUE iff objective value changed, and data of LP solver has to be updated */
   unsigned int     lbchanged:1;        /**< TRUE iff lower bound changed, and data of LP solver has to be updated */
   unsigned int     ubchanged:1;        /**< TRUE iff upper bound changed, and data of LP solver has to be updated */
   unsigned int     coefchanged:1;      /**< TRUE iff the coefficient vector changed, and LP solver has to be updated */
   unsigned int     removeable:1;       /**< TRUE iff column is removeable from the LP (due to aging or cleanup) */
};

/** row of the LP */
struct Row
{
   char*            name;               /**< name of the row */
   COL**            cols;               /**< columns of row entries, that may have a nonzero primal solution value */
   int*             cols_probindex;     /**< copy of cols[i]->var->probindex for avoiding expensive dereferencing */
   Real*            vals;               /**< coefficients of row entries */
   int*             linkpos;            /**< position of row in row vector of the column, or -1 if not yet linked */
   Real             constant;           /**< constant shift c in row lhs <= ax + c <= rhs */
   Real             lhs;                /**< left hand side of row */
   Real             rhs;                /**< right hand side of row */
   Real             sqrnorm;            /**< squared euclidean norm of row vector */
   Real             maxval;             /**< maximal absolute value of row vector, only valid if nummaxval > 0 */
   Real             minval;             /**< minimal absolute non-zero value of row vector, only valid if numminval > 0 */
   Real             dualsol;            /**< dual solution value in LP, is 0 if row is not in LP */
   Real             activity;           /**< row activity value in LP, or SCIP_INVALID if not yet calculated */
   Real             dualfarkas;         /**< multiplier value in dual farkas infeasibility proof */
   Real             pseudoactivity;     /**< row activity value in pseudo solution, or SCIP_INVALID if not yet calculated */
   Real             minactivity;        /**< minimal activity value w.r.t. the column's bounds, or SCIP_INVALID */
   Real             maxactivity;        /**< maximal activity value w.r.t. the column's bounds, or SCIP_INVALID */
   Longint          validpsactivitybdchg; /**< bound change number for which pseudo activity value is valid */
   Longint          validactivitybdsbdchg;/**< bound change number for which activity bound values are valid */
   Longint          obsoletenode;       /**< last node where this row was removed due to aging */
   int              index;              /**< consecutively numbered row identifier */
   int              size;               /**< size of the col- and val-arrays */
   int              len;                /**< number of nonzeros in row */
   int              nunlinked;          /**< number of row entries, where the columns don't know about the row */
   int              nuses;              /**< number of times, this row is referenced */
   int              lppos;              /**< row position number in actual LP, or -1 if not in actual LP */
   int              lpipos;             /**< row position number in LP solver, or -1 if not in LP solver */
   int              minidx;             /**< minimal column index of row entries */
   int              maxidx;             /**< maximal column index of row entries */
   int              nummaxval;          /**< number of coefs with absolute value equal to maxval, zero if maxval invalid */
   int              numminval;          /**< number of coefs with absolute value equal to minval, zero if minval invalid */
   int              validactivitylp;    /**< LP number for which activity value is valid */
   int              age;                /**< number of successive times this row was in LP and was not sharp in solution */
   unsigned int     sorted:1;           /**< are column indices sorted in increasing order? */
   unsigned int     delaysort:1;        /**< should the row sorting be delayed and done in a lazy fashion? */
   unsigned int     validminmaxidx:1;   /**< are minimal and maximal column index valid? */
   unsigned int     lhschanged:1;       /**< was left hand side or constant changed, and has LP solver to be updated? */
   unsigned int     rhschanged:1;       /**< was right hand side or constant changed, and has LP solver to be updated? */
   unsigned int     coefchanged:1;      /**< was the coefficient vector changed, and has LP solver to be updated? */
   unsigned int     local:1;            /**< is row only valid locally? */
   unsigned int     modifiable:1;       /**< is row modifiable during node processing (subject to column generation)? */
   unsigned int     removeable:1;       /**< TRUE iff row is removeable from the LP (due to aging or cleanup) */
   unsigned int     nlocks:24;          /**< number of sealed locks of an unmodifiable row */
};

/** actual LP data */
struct Lp
{
   LPI*             lpi;                /**< LP solver interface */
   LPISTATE*        divelpistate;       /**< stores LPI state (basis information) before diving starts */
   COL**            lpicols;            /**< array with columns actually stored in the LP solver */
   int              lpicolssize;        /**< available slots in lpicols vector */
   int              nlpicols;           /**< number of columns in the LP solver */
   int              lpifirstchgcol;     /**< first column of the LP which differs from the column in the LP solver */
   ROW**            lpirows;            /**< array with rows actually stored in the LP solver */
   int              lpirowssize;        /**< available slots in lpirows vector */
   int              nlpirows;           /**< number of rows in the LP solver */
   int              lpifirstchgrow;     /**< first row of the LP which differs from the row in the LP solver */
   COL**            chgcols;            /**< array of changed columns not yet applied to the LP solver */
   int              chgcolssize;        /**< available slots in chgcols vector */
   int              nchgcols;           /**< actual number of chgcols (number of used slots in chgcols vector) */
   ROW**            chgrows;            /**< array of changed rows not yet applied to the LP solver */
   int              chgrowssize;        /**< available slots in chgrows vector */
   int              nchgrows;           /**< actual number of chgrows (number of used slots in chgrows vector) */
   COL**            cols;               /**< array with actual LP columns in correct order */
   int              colssize;           /**< available slots in cols vector */
   int              ncols;              /**< actual number of LP columns (number of used slots in cols vector) */
   int              nremoveablecols;    /**< number of removeable columns in the LP */
   int              firstnewcol;        /**< first column added at the active node */
   ROW**            rows;               /**< array with actual LP rows in correct order */
   int              rowssize;           /**< available slots in rows vector */
   int              nrows;              /**< actual number of LP rows (number of used slots in rows vector) */
   int              nremoveablerows;    /**< number of removeable rows in the LP */
   int              firstnewrow;        /**< first row added at the active node */
   LPSOLSTAT        lpsolstat;          /**< solution status of last LP solution */
   Real             lpobjval;           /**< objective value of LP without loose variables, or SCIP_INVALID */
   Real             looseobjval;        /**< actual solution value of all loose variables set to their best bounds,
                                         *   ignoring variables, with infinite best bound */
   int              looseobjvalinf;     /**< number of loose variables with infinite best bound in actual solution */
   int              nloosevars;         /**< number of loose variables in LP */
   Real             pseudoobjval;       /**< actual pseudo solution value with all variables set to their best bounds,
                                         *   ignoring variables, with infinite best bound */
   int              pseudoobjvalinf;    /**< number of variables with infinite best bound in actual pseudo solution */
   Real             upperbound;         /**< upper objective limit of LP (copy of primal->upperbound) */
   int              validsollp;         /**< LP number for which the currently stored solution values are valid */
   int              validfarkaslp;      /**< LP number for which the currently stored farkas values are valid */
   Bool             flushdeletedcols;   /**< have LPI-columns been deleted in the last lpFlush() call? */
   Bool             flushaddedcols;     /**< have LPI-columns been added in the last lpFlush() call? */
   Bool             flushdeletedrows;   /**< have LPI-rows been deleted in the last lpFlush() call? */
   Bool             flushaddedrows;     /**< have LPI-rows been added in the last lpFlush() call? */
   Bool             flushed;            /**< are all cached changes applied to the LP solver? */
   Bool             solved;             /**< is current LP solved? */
   Bool             primalfeasible;     /**< is actual LP basis primal feasible? */
   Bool             dualfeasible;       /**< is actual LP basis dual feasible? */
   Bool             diving;             /**< LP is used for diving: col bounds and obj don't corresond to variables */
   Real             lpiuobjlim;         /**< current upper objective limit in LPI */
   Real             lpifeastol;         /**< current feasibility tolerance in LPI */
   Bool             lpifromscratch;     /**< current FROMSCRATCH setting in LPI */
   Bool             lpifastmip;         /**< current FASTMIP setting in LPI */
   Bool             lpiscaling;         /**< current SCALING setting in LPI */
   Bool             lastwasprimal;      /**< was the last simplex call a call to the primal simplex? */
};


#endif
