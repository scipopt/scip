/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_lp.h,v 1.40 2005/08/22 18:35:50 bzfpfend Exp $"

/**@file   struct_lp.h
 * @brief  datastructures for SCIP_LP management
 * @author Tobias Achterberg
 *
 *  In SCIP, the SCIP_LP is defined as follows:
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
 *  The main datastructures for storing an SCIP_LP are the rows and the columns.
 *  A row can live on its own (if it was created by a separator), or as SCIP_LP
 *  relaxation of a constraint. Thus, it has a nuses-counter, and is
 *  deleted, if not needed any more.
 *  A column cannot live on its own. It is always connected to a problem
 *  variable. Because pricing is always problem specific, it cannot create
 *  SCIP_LP columns without introducing new variables. Thus, each column is
 *  connected to exactly one variable, and is deleted, if the variable
 *  is deleted.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_LP_H__
#define __SCIP_STRUCT_LP_H__


#include "scip/def.h"
#include "scip/type_lpi.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"



/** SCIP_LP column;
 *  The row vector of the SCIP_LP column is partitioned into two parts: The first col->nlprows rows in the rows array
 *  are the ones that belong to the current SCIP_LP (col->rows[j]->lppos >= 0) and that are linked to the column
 *  (col->linkpos[j] >= 0). The remaining col->len - col->nlprows rows in the rows array are the ones that
 *  don't belong to the current SCIP_LP (col->rows[j]->lppos == -1) or that are not linked to the column
 *  (col->linkpos[j] == -1).
 */
struct SCIP_Col
{
   SCIP_Real             obj;                /**< current objective value of column in SCIP_LP */
   SCIP_Real             lb;                 /**< current lower bound of column in SCIP_LP */
   SCIP_Real             ub;                 /**< current upper bound of column in SCIP_LP */
   SCIP_Real             flushedobj;         /**< objective value of column already flushed to the SCIP_LP solver */
   SCIP_Real             flushedlb;          /**< lower bound of column already flushed to the SCIP_LP solver */
   SCIP_Real             flushedub;          /**< upper bound of column already flushed to the SCIP_LP solver */
   SCIP_Real             primsol;            /**< primal solution value in LP, is 0 if col is not in SCIP_LP */
   SCIP_Real             redcost;            /**< reduced cost value in LP, or SCIP_INVALID if not yet calculated */
   SCIP_Real             farkascoef;         /**< coefficient in dual farkas infeasibility proof (== dualfarkas^T A_c) */
   SCIP_Real             minprimsol;         /**< minimal SCIP_LP solution value, this column ever assumed */
   SCIP_Real             maxprimsol;         /**< maximal SCIP_LP solution value, this column ever assumed */
   SCIP_Real             sbdown;             /**< strong branching information for downwards branching */
   SCIP_Real             sbup;               /**< strong branching information for upwards branching */
   SCIP_Real             sbsolval;           /**< SCIP_LP solution value of column at last strong branching call */
   SCIP_Real             sblpobjval;         /**< SCIP_LP objective value at last strong branching call on the column */
   SCIP_Longint          sbnode;             /**< node number of the last strong branching call on this column */
   SCIP_Longint          obsoletenode;       /**< last node where this column was removed due to aging */
   SCIP_VAR*             var;                /**< variable, this column represents; there cannot be a column without variable */
   SCIP_ROW**            rows;               /**< rows of column entries, that may have a nonzero dual solution value */
   SCIP_Real*            vals;               /**< coefficients of column entries */
   int*                  linkpos;            /**< position of col in col vector of the row, or -1 if not yet linked */
   int                   index;              /**< consecutively numbered column identifier */
   int                   size;               /**< size of the row- and val-arrays */
   int                   len;                /**< number of nonzeros in column */
   int                   nlprows;            /**< number of linked rows in column, that belong to the current SCIP_LP */
   int                   nunlinked;          /**< number of column entries, where the rows don't know about the column */
   int                   lppos;              /**< column position number in current LP, or -1 if not in current SCIP_LP */
   int                   lpipos;             /**< column position number in SCIP_LP solver, or -1 if not in SCIP_LP solver */
   int                   lpdepth;            /**< depth level at which column entered the LP, or -1 if not in current SCIP_LP */
   int                   validredcostlp;     /**< SCIP_LP number for which reduced cost value is valid */
   int                   validfarkaslp;      /**< SCIP_LP number for which farkas coefficient is valid */
   int                   validsblp;          /**< SCIP_LP number for which strong branching values are valid */
   int                   sbitlim;            /**< strong branching iteration limit used to get strongbranch values, or -1 */
   int                   nsbcalls;           /**< number of times, strong branching was applied on the column */
   int                   age;                /**< number of successive times this variable was in SCIP_LP and was 0.0 in solution */
   int                   var_probindex;      /**< copy of var->probindex for avoiding expensive dereferencing */
   unsigned int          basisstatus:2;      /**< basis status of column in last SCIP_LP solution, invalid for non-LP columns */
   unsigned int          lprowssorted:1;     /**< are the linked SCIP_LP rows in the rows array sorted by non-decreasing index? */
   unsigned int          nonlprowssorted:1;  /**< are the non-LP/not linked rows sorted by non-decreasing index? */
   unsigned int          objchanged:1;       /**< has objective value changed, and has data of SCIP_LP solver to be updated? */
   unsigned int          lbchanged:1;        /**< has lower bound changed, and has data of SCIP_LP solver to be updated? */
   unsigned int          ubchanged:1;        /**< has upper bound changed, and has data of SCIP_LP solver to be updated? */
   unsigned int          coefchanged:1;      /**< has the coefficient vector changed, and has SCIP_LP solver to be updated? */
   unsigned int          integral:1;         /**< is associated variable of integral type? */
   unsigned int          removeable:1;       /**< is column removeable from the SCIP_LP (due to aging or cleanup)? */
   unsigned int          sbdownvalid:1;      /**< stores whether the stored strong branching down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   unsigned int          sbupvalid:1;        /**< stores whether the stored strong branching up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
};

/** SCIP_LP row
 *  The column vector of the SCIP_LP row is partitioned into two parts: The first row->nlpcols columns in the cols array
 *  are the ones that belong to the current SCIP_LP (row->cols[j]->lppos >= 0) and that are linked to the row   
 *  (row->linkpos[j] >= 0). The remaining row->len - row->nlpcols columns in the cols array are the ones that
 *  don't belong to the current SCIP_LP (row->cols[j]->lppos == -1) or that are not linked to the row   
 *  (row->linkpos[j] == -1).
 */
struct SCIP_Row
{
   SCIP_Real             constant;           /**< constant shift c in row lhs <= ax + c <= rhs */
   SCIP_Real             lhs;                /**< left hand side of row */
   SCIP_Real             rhs;                /**< right hand side of row */
   SCIP_Real             flushedlhs;         /**< left hand side minus constant of row already flushed to the SCIP_LP solver */
   SCIP_Real             flushedrhs;         /**< right hand side minus constant of row already flushed to the SCIP_LP solver */
   SCIP_Real             sqrnorm;            /**< squared euclidean norm of row vector */
   SCIP_Real             sumnorm;            /**< sum norm of row vector (sum of absolute values of coefficients) */
   SCIP_Real             objprod;            /**< scalar product of row vector with objective function */
   SCIP_Real             maxval;             /**< maximal absolute value of row vector, only valid if nummaxval > 0 */
   SCIP_Real             minval;             /**< minimal absolute non-zero value of row vector, only valid if numminval > 0 */
   SCIP_Real             dualsol;            /**< dual solution value in LP, is 0 if row is not in SCIP_LP */
   SCIP_Real             activity;           /**< row activity value in LP, or SCIP_INVALID if not yet calculated */
   SCIP_Real             dualfarkas;         /**< multiplier value in dual farkas infeasibility proof */
   SCIP_Real             pseudoactivity;     /**< row activity value in pseudo solution, or SCIP_INVALID if not yet calculated */
   SCIP_Real             minactivity;        /**< minimal activity value w.r.t. the column's bounds, or SCIP_INVALID */
   SCIP_Real             maxactivity;        /**< maximal activity value w.r.t. the column's bounds, or SCIP_INVALID */
   SCIP_Longint          validpsactivitydomchg; /**< domain change number for which pseudo activity value is valid */
   SCIP_Longint          validactivitybdsdomchg;/**< domain change number for which activity bound values are valid */
   SCIP_Longint          obsoletenode;       /**< last node where this row was removed due to aging */
   char*                 name;               /**< name of the row */
   SCIP_COL**            cols;               /**< columns of row entries, that may have a nonzero primal solution value */
   int*                  cols_index;         /**< copy of cols[i]->index for avoiding expensive dereferencing */
   SCIP_Real*            vals;               /**< coefficients of row entries */
   int*                  linkpos;            /**< position of row in row vector of the column, or -1 if not yet linked */
   int                   index;              /**< consecutively numbered row identifier */
   int                   size;               /**< size of the col- and val-arrays */
   int                   len;                /**< number of nonzeros in row */
   int                   nlpcols;            /**< number of linked columns in row, that belong to the current SCIP_LP */
   int                   nunlinked;          /**< number of row entries, where the columns don't know about the row */
   int                   nuses;              /**< number of times, this row is referenced */
   int                   lppos;              /**< row position number in current LP, or -1 if not in current SCIP_LP */
   int                   lpipos;             /**< row position number in SCIP_LP solver, or -1 if not in SCIP_LP solver */
   int                   lpdepth;            /**< depth level at which row entered the LP, or -1 if not in current SCIP_LP */
   int                   minidx;             /**< minimal column index of row entries */
   int                   maxidx;             /**< maximal column index of row entries */
   int                   nummaxval;          /**< number of coefs with absolute value equal to maxval, zero if maxval invalid */
   int                   numminval;          /**< number of coefs with absolute value equal to minval, zero if minval invalid */
   int                   validactivitylp;    /**< SCIP_LP number for which activity value is valid */
   int                   age;                /**< number of successive times this row was in SCIP_LP and was not sharp in solution */
   unsigned int          basisstatus:2;      /**< basis status of row in last SCIP_LP solution, invalid for non-LP rows */
   unsigned int          lpcolssorted:1;     /**< are the linked SCIP_LP columns in the cols array sorted by non-decreasing index? */
   unsigned int          nonlpcolssorted:1;  /**< are the non-LP/not linked columns sorted by non-decreasing index? */
   unsigned int          delaysort:1;        /**< should the row sorting be delayed and done in a lazy fashion? */
   unsigned int          validminmaxidx:1;   /**< are minimal and maximal column index valid? */
   unsigned int          lhschanged:1;       /**< was left hand side or constant changed, and has SCIP_LP solver to be updated? */
   unsigned int          rhschanged:1;       /**< was right hand side or constant changed, and has SCIP_LP solver to be updated? */
   unsigned int          coefchanged:1;      /**< was the coefficient vector changed, and has SCIP_LP solver to be updated? */
   unsigned int          integral:1;         /**< is activity (without constant) of row always integral in feasible solution? */
   unsigned int          local:1;            /**< is row only valid locally? */
   unsigned int          modifiable:1;       /**< is row modifiable during node processing (subject to column generation)? */
   unsigned int          removeable:1;       /**< is row removeable from the SCIP_LP (due to aging or cleanup)? */
   unsigned int          nlocks:19;          /**< number of sealed locks of an unmodifiable row */
};

/** current SCIP_LP data */
struct SCIP_Lp
{
   SCIP_Real             lpobjval;           /**< objective value of SCIP_LP without loose variables, or SCIP_INVALID */
   SCIP_Real             looseobjval;        /**< current solution value of all loose variables set to their best bounds,
                                              *   ignoring variables, with infinite best bound */
   SCIP_Real             pseudoobjval;       /**< current pseudo solution value with all variables set to their best bounds,
                                              *   ignoring variables, with infinite best bound */
   SCIP_Real             cutoffbound;        /**< upper objective limit of SCIP_LP (copy of primal->cutoffbound) */
   SCIP_Real             lpiuobjlim;         /**< current upper objective limit in SCIP_LPI */
   SCIP_Real             lpifeastol;         /**< current feasibility tolerance in SCIP_LPI */
   SCIP_Real             lpidualfeastol;     /**< current reduced costs feasibility tolerance in SCIP_LPI */
   SCIP_Real             lpibarrierconvtol;  /**< current convergence tolerance used in barrier algorithm in SCIP_LPI */
   SCIP_Real             objsqrnorm;         /**< squared euclidean norm of objective function vector of problem variables */
   SCIP_Real             objsumnorm;         /**< sum norm of objective function vector of problem variables */
   SCIP_LPI*             lpi;                /**< SCIP_LP solver interface */
   SCIP_COL**            lpicols;            /**< array with columns currently stored in the SCIP_LP solver */
   SCIP_ROW**            lpirows;            /**< array with rows currently stored in the SCIP_LP solver */
   SCIP_COL**            chgcols;            /**< array of changed columns not yet applied to the SCIP_LP solver */
   SCIP_ROW**            chgrows;            /**< array of changed rows not yet applied to the SCIP_LP solver */
   SCIP_COL**            cols;               /**< array with current SCIP_LP columns in correct order */
   SCIP_ROW**            rows;               /**< array with current SCIP_LP rows in correct order */
   SCIP_LPISTATE*        divelpistate;       /**< stores SCIP_LPI state (basis information) before diving starts */
   int                   lpicolssize;        /**< available slots in lpicols vector */
   int                   nlpicols;           /**< number of columns in the SCIP_LP solver */
   int                   lpifirstchgcol;     /**< first column of the SCIP_LP which differs from the column in the SCIP_LP solver */
   int                   lpirowssize;        /**< available slots in lpirows vector */
   int                   nlpirows;           /**< number of rows in the SCIP_LP solver */
   int                   lpifirstchgrow;     /**< first row of the SCIP_LP which differs from the row in the SCIP_LP solver */
   int                   chgcolssize;        /**< available slots in chgcols vector */
   int                   nchgcols;           /**< current number of chgcols (number of used slots in chgcols vector) */
   int                   chgrowssize;        /**< available slots in chgrows vector */
   int                   nchgrows;           /**< current number of chgrows (number of used slots in chgrows vector) */
   int                   colssize;           /**< available slots in cols vector */
   int                   ncols;              /**< current number of SCIP_LP columns (number of used slots in cols vector) */
   int                   nremoveablecols;    /**< number of removeable columns in the SCIP_LP */
   int                   firstnewcol;        /**< first column added at the current node */
   int                   rowssize;           /**< available slots in rows vector */
   int                   nrows;              /**< current number of SCIP_LP rows (number of used slots in rows vector) */
   int                   nremoveablerows;    /**< number of removeable rows in the SCIP_LP */
   int                   firstnewrow;        /**< first row added at the current node */
   int                   looseobjvalinf;     /**< number of loose variables with infinite best bound in current solution */
   int                   nloosevars;         /**< number of loose variables in SCIP_LP */
   int                   pseudoobjvalinf;    /**< number of variables with infinite best bound in current pseudo solution */
   int                   validsollp;         /**< SCIP_LP number for which the currently stored solution values are valid */
   int                   validfarkaslp;      /**< SCIP_LP number for which the currently stored farkas row multipliers are valid */
   int                   lpiitlim;           /**< current iteration limit setting in SCIP_LPI */
   SCIP_PRICING          lpipricing;         /**< current pricing setting in SCIP_LPI */
   SCIP_LPSOLSTAT        lpsolstat;          /**< solution status of last SCIP_LP solution */
   SCIP_LPALGO           lastlpalgo;         /**< algorithm used for last SCIP_LP solve */
   SCIP_Bool             flushdeletedcols;   /**< have LPI-columns been deleted in the last lpFlush() call? */
   SCIP_Bool             flushaddedcols;     /**< have LPI-columns been added in the last lpFlush() call? */
   SCIP_Bool             flushdeletedrows;   /**< have LPI-rows been deleted in the last lpFlush() call? */
   SCIP_Bool             flushaddedrows;     /**< have LPI-rows been added in the last lpFlush() call? */
   SCIP_Bool             flushed;            /**< are all cached changes applied to the SCIP_LP solver? */
   SCIP_Bool             solved;             /**< is current SCIP_LP solved? */
   SCIP_Bool             primalfeasible;     /**< is current SCIP_LP solution primal feasible? */
   SCIP_Bool             dualfeasible;       /**< is current SCIP_LP solution dual feasible? */
   SCIP_Bool             solisbasic;         /**< is current SCIP_LP solution a basic solution? */
   SCIP_Bool             probing;            /**< are we currently in probing mode? */
   SCIP_Bool             diving;             /**< SCIP_LP is used for diving: col bounds and obj don't corresond to variables */
   SCIP_Bool             divingobjchg;       /**< objective values were changed in diving: SCIP_LP objective is invalid */
   SCIP_Bool             lpifromscratch;     /**< current FROMSCRATCH setting in SCIP_LPI */
   SCIP_Bool             lpifastmip;         /**< current FASTMIP setting in SCIP_LPI */
   SCIP_Bool             lpiscaling;         /**< current SCALING setting in SCIP_LPI */
   SCIP_Bool             lpipresolving;      /**< current PRESOLVING setting in SCIP_LPI */
   SCIP_Bool             lpilpinfo;          /**< current LPINFO setting in SCIP_LPI */
   SCIP_Bool             lpihasfeastol;      /**< does the SCIP_LPI support the FEASTOL parameter? */
   SCIP_Bool             lpihasdualfeastol;  /**< does the SCIP_LPI support the DUALFEASTOL parameter? */
   SCIP_Bool             lpihasbarrierconvtol;/**< does the SCIP_LPI support the BARRIERCONVTOL parameter? */
   SCIP_Bool             lpihasfastmip;      /**< does the SCIP_LPI support the FASTMIP parameter? */
   SCIP_Bool             lpihasscaling;      /**< does the SCIP_LPI support the SCALING parameter? */
   SCIP_Bool             lpihaspresolving;   /**< does the SCIP_LPI support the PRESOLVING parameter? */
};


#endif
