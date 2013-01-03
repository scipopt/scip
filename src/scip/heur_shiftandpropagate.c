/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_shiftandpropagate.c
 * @brief  shiftandpropagate primal heuristic
 * @author Timo Berthold
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/heur_shiftandpropagate.h"

#define HEUR_NAME             "shiftandpropagate"
#define HEUR_DESC             "Pre-root heuristic to expand an auxiliary branch-and-bound tree and apply propagation techniques"
#define HEUR_DISPCHAR         'T'
#define HEUR_PRIORITY         1000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE     /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_WEIGHT_INEQUALITY   1   /**< the heuristic row weight for inequalities */
#define DEFAULT_WEIGHT_EQUALITY     3   /**< the heuristic row weight for equations */
#define DEFAULT_RELAX            TRUE   /**< Should continuous variables be relaxed from the problem? */
#define DEFAULT_PROBING          TRUE   /**< Is propagation of solution values enabled? */
#define DEFAULT_ONLYWITHOUTSOL   TRUE   /**< Should heuristic only be executed if no primal solution was found, yet? */
#define DEFAULT_NPROPROUNDS        10   /**< The default number of propagation rounds for each propagation used */
#define DEFAULT_PROPBREAKER     65000   /**< fixed maximum number of propagations */
#define DEFAULT_CUTOFFBREAKER      15   /**< fixed maximum number of allowed cutoffs before the heuristic stops */
#define DEFAULT_RANDSEED      3141598   /**< the default random seed for random number generation */
#define DEFAULT_SORTKEY            'u'  /**< the default key for variable sorting */
#define DEFAULT_SORTVARS         TRUE   /**< should variables be processed in sorted order? */
#define SORTKEYS                 "nru"  /**< options sorting key: (n)orms down, norms (u)p or (r)andom */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_COL**            lpcols;             /**< stores lp columns with discrete variables before cont. variables */
   SCIP_Bool             relax;              /**< should continuous variables be relaxed from the problem */
   SCIP_Bool             probing;            /**< should probing be executed? */
   SCIP_Bool             onlywithoutsol;     /**< Should heuristic only be executed if no primal solution was found, yet? */
   int                   nlpcols;            /**< the number of lp columns */
   int                   lpcolssize;         /**< size of lpcols array */
   int                   nproprounds;        /**< The default number of propagation rounds for each propagation used */
   int                   cutoffbreaker;      /**< the number of cutoffs before heuristic execution is stopped, or -1 for no
                                               * limit */

   unsigned int          randseed;           /**< seed for random number generation */
   char                  sortkey;            /**< the key by which variables are sorted */
   SCIP_Bool             sortvars;           /**< should variables be processed in sorted order? */

   SCIPstatistic(
      SCIP_LPSOLSTAT     lpsolstat;          /**< the probing status after probing */
      SCIP_Longint       ntotaldomredsfound; /**< the total number of domain reductions during heuristic */
      SCIP_Longint       nlpiters;           /**< number of LP iterations which the heuristic needed */
      int                nremainingviols;    /**< the number of remaining violations */
      int                nprobings;          /**< how many probings has the heuristic executed? */
      int                ncutoffs;           /**< has the probing node been cutoff? */
      int                nredundantrows;     /**< how many rows were redundant after relaxation? */
      )
};

/** status of a variable in heuristic transformation */
enum TransformStatus
{
   TRANSFORMSTATUS_NONE = 0,            /**< variable has not been transformed yet */
   TRANSFORMSTATUS_LB   = 1,            /**< variable has been shifted by using lower bound (x-lb) */
   TRANSFORMSTATUS_NEG  = 2,            /**< variable has been negated by using upper bound (ub-x) */
   TRANSFORMSTATUS_FREE = 3             /**< variable does not have to be shifted */
};
typedef enum TransformStatus TRANSFORMSTATUS;

/** information about the matrix after its heuristic transformation */
struct ConstraintMatrix
{
   SCIP_Real*            rowmatvals;         /**< matrix coefficients row by row */
   int*                  rowmatind;          /**< the indices of the corresponding variables */
   int*                  rowmatbegin;        /**< the starting indices of each row */
   SCIP_Real*            colmatvals;         /**< matrix coefficients column by column */
   int*                  colmatind;          /**< the indices of the corresponding rows for each coefficient */
   int*                  colmatbegin;        /**< the starting indices of each column */
   TRANSFORMSTATUS*      transformstatus;    /**< information about transform status of every discrete variable */
   SCIP_Real*            lhs;                /**< left hand side vector after normalization */
   SCIP_Real*            rhs;                /**< right hand side vector after normalization */
   SCIP_Real*            colnorms;           /**< vector norms of all discrete problem variables after normalization */
   SCIP_Real*            upperbounds;        /**< the upper bounds of every non-continuous variable after transformation*/
   SCIP_Real*            transformshiftvals; /**< values by which original discrete variable bounds were shifted */
   int                   nnonzs;             /**< number of nonzero column entries */
   int                   nrows;              /**< number of rows of matrix */
   int                   ncols;              /**< the number of columns in matrix (including continuous vars) */
   int                   ndiscvars;          /**< number of discrete problem variables */
   SCIP_Bool             normalized;         /**< indicates if the matrix data has already been normalized */
};
typedef struct ConstraintMatrix CONSTRAINTMATRIX;

/*
 * Local methods
 */


/** returns nonzero values and corresponding columns of given row */
static
void getRowData(
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int                   rowindex,           /**< index of the desired row */
   SCIP_Real**           valpointer,         /**< pointer to store the nonzero coefficients of the row */
   SCIP_Real*            lhs,                /**< lhs of the row */
   SCIP_Real*            rhs,                /**< rhs of the row */
   int**                 indexpointer,       /**< pointer to store column indices which belong to the nonzeros */
   int*                  nrowvals            /**< pointer to store number of nonzeros in the desired row */
   )
{
   int arrayposition;

   assert(matrix != NULL);
   assert(0 <= rowindex && rowindex < matrix->nrows);

   arrayposition = matrix->rowmatbegin[rowindex];

   if( nrowvals != NULL && rowindex == matrix->nrows - 1 )
      *nrowvals = matrix->nnonzs - arrayposition;
   else if( nrowvals != NULL )
      *nrowvals = matrix->rowmatbegin[rowindex + 1] - arrayposition;

   if( valpointer != NULL )
      *valpointer = &(matrix->rowmatvals[arrayposition]);
   if( indexpointer != NULL )
      *indexpointer = &(matrix->rowmatind[arrayposition]);

   if( lhs != NULL )
      *lhs = matrix->lhs[rowindex];

   if( rhs != NULL )
      *rhs = matrix->rhs[rowindex];
}

/** returns nonzero values and corresponding rows of given column */
static
void getColumnData(
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int                   colindex,           /**< the index of the desired column */
   SCIP_Real**           valpointer,         /**< pointer to store the nonzero coefficients of the column */
   int**                 indexpointer,       /**< pointer to store row indices which belong to the nonzeros */
   int*                  ncolvals            /**< pointer to store number of nonzeros in the desired column */
   )
{
   int arrayposition;

   assert(matrix != NULL);
   assert(0 <= colindex && colindex < matrix->ncols);

   arrayposition = matrix->colmatbegin[colindex];

   if( ncolvals != NULL )
   {
      if( colindex == matrix->ncols - 1 )
         *ncolvals = matrix->nnonzs - arrayposition;
      else
         *ncolvals = matrix->colmatbegin[colindex + 1] - arrayposition;
   }
   if( valpointer != NULL )
      *valpointer = &(matrix->colmatvals[arrayposition]);

   if( indexpointer != NULL )
      *indexpointer = &(matrix->colmatind[arrayposition]);
}

/** relaxes a continuous variable from all its rows, which has influence
 *  on both the left and right hand side of the constraint.
 */
static
void relaxVar(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_VAR*             var,                /**< variable which is relaxed from the problem */
   CONSTRAINTMATRIX*     matrix              /**< constraint matrix object */
   )
{
   SCIP_ROW** colrows;
   SCIP_COL* varcol;
   SCIP_Real* colvals;
   SCIP_Real ub;
   SCIP_Real lb;
   int ncolvals;
   int r;

   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

   varcol = SCIPvarGetCol(var);
   assert(varcol != NULL);

   /* get nonzero values and corresponding rows of variable */
   colvals = SCIPcolGetVals(varcol);
   ncolvals = SCIPcolGetNLPNonz(varcol);
   colrows = SCIPcolGetRows(varcol);

   ub = SCIPvarGetUbGlobal(var);
   lb = SCIPvarGetLbGlobal(var);

   assert(colvals != NULL || ncolvals == 0);

   SCIPdebugMessage("Relaxing variable <%s> with lb <%g> and ub <%g>\n",
      SCIPvarGetName(var), lb, ub);

   assert(matrix->normalized);
   /* relax variable from all its constraints */
   for( r = 0; r < ncolvals; ++r )
   {
      SCIP_ROW* colrow;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real lhsvarbound;
      SCIP_Real rhsvarbound;
      SCIP_Real rowabs;
      SCIP_Real colval;
      int rowindex;

      colrow = colrows[r];
      rowindex = SCIProwGetLPPos(colrow);

      if( rowindex == -1 )
         break;

      rowabs = SCIPgetRowMaxCoef(scip, colrow);
      assert(SCIPisFeasGT(scip, rowabs, 0.0));
      assert(colvals != NULL); /* to please flexelint */
      colval = colvals[r]/rowabs;

      assert(0 <= rowindex && rowindex < matrix->nrows);
      getRowData(matrix, rowindex, NULL, &lhs, &rhs, NULL, NULL);
      /* variables bound influence the lhs and rhs of current row depending on the sign
       * of the variables coefficient.
       */
      if( SCIPisFeasPositive(scip, colval) )
      {
         lhsvarbound = ub;
         rhsvarbound = lb;
      }
      else if( SCIPisFeasNegative(scip, colval) )
      {
         lhsvarbound = lb;
         rhsvarbound = ub;
      }
      else
         continue;

      /* relax variable from the current row */
      if( !SCIPisInfinity(scip, -matrix->lhs[rowindex]) && !SCIPisInfinity(scip, ABS(lhsvarbound)) )
         matrix->lhs[rowindex] -= colval * lhsvarbound;
      else
         matrix->lhs[rowindex] = -SCIPinfinity(scip);

      if( !SCIPisInfinity(scip, matrix->rhs[rowindex]) && !SCIPisInfinity(scip, ABS(rhsvarbound)) )
         matrix->rhs[rowindex] -= colval * rhsvarbound;
      else
         matrix->rhs[rowindex] = SCIPinfinity(scip);

      SCIPdebugMessage("Row <%d> changed:Coefficient <%g>, LHS <%g> --> <%g>, RHS <%g> --> <%g>\n",
         rowindex, colval, lhs, matrix->lhs[rowindex], rhs, matrix->rhs[rowindex]);
   }
}

/** transforms bounds of a given variable s.t. its lower bound equals zero afterwards.
 *  If the variable already has lower bound zero, the variable is not transformed,
 *  if not, the variable's bounds are changed w.r.t. the smaller absolute value of its
 *  bounds in order to avoid numerical inaccuracies. If both lower and upper bound
 *  of the variable differ from infinity, there are two cases. If |lb| <= |ub|,
 *  the bounds are shifted by -lb, else a new variable ub - x replaces x.
 *  The transformation is memorized by the transform status of the variable s.t.
 *  retransformation is possible.
 */
static
void transformVariable(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   colpos              /**< position of variable column in matrix */
   )
{
   SCIP_COL* col;
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(matrix != NULL);
   assert(0 <= colpos && colpos < heurdata->nlpcols);
   col = heurdata->lpcols[colpos];
   assert(col != NULL);
   assert(SCIPcolIsInLP(col));

   var = SCIPcolGetVar(col);
   assert(var != NULL);
   assert(SCIPvarIsIntegral(var));
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   assert(matrix->transformstatus[colpos] == TRANSFORMSTATUS_NONE
      || matrix->transformstatus[colpos] == TRANSFORMSTATUS_FREE);

   /* if both lower and upper bound are -infinity and infinity, resp., this is reflected by a free transform status.
    * If the lower bound is already zero, this is reflected by identity transform status. In both cases, none of the
    * corresponding rows needs to be modified.
    */
   if( SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub) )
   {
      matrix->upperbounds[colpos] = ub;
      matrix->transformstatus[colpos] = TRANSFORMSTATUS_FREE;
      matrix->transformshiftvals[colpos] = 0.0;
   }
   else if( SCIPisFeasEQ(scip, lb, 0.0) )
   {
      matrix->transformstatus[colpos] = TRANSFORMSTATUS_LB;
      matrix->upperbounds[colpos] = ub;
      matrix->transformshiftvals[colpos] = 0.0;
   }
   else
   {
      SCIP_Real* vals;
      int* rows;

      SCIP_Real bound;
      int nrows;
      int i;

      /* a real transformation is necessary. The variable x is either shifted by -lb or
       * replaced by ub - x, depending on the smaller absolute of lb and ub.
       */
      if( SCIPisFeasLE(scip, ABS(lb), ABS(ub)) )
      {
         matrix->transformstatus[colpos] = TRANSFORMSTATUS_LB;
         bound = lb;
      }
      else
      {
         matrix->transformstatus[colpos] = TRANSFORMSTATUS_NEG;
         bound = ub;
      }

      assert(!SCIPisInfinity(scip, bound));

      /* get nonzero values and corresponding rows of column */
      getColumnData(matrix, colpos, &vals, &rows, &nrows);
      assert(nrows == 0 ||(vals != NULL && rows != NULL));

      /* go through rows and modify its lhs, rhs and the variable coefficient, if necessary */
      for( i = 0; i < nrows; ++i )
      {
         assert(rows[i] >= 0);
         assert(rows[i] < matrix->nrows);

         if( !SCIPisInfinity(scip, -(matrix->lhs[rows[i]])) )
            matrix->lhs[rows[i]] -= (vals[i]) * bound;

         if( !SCIPisInfinity(scip, matrix->rhs[rows[i]]) )
            matrix->rhs[rows[i]] -= (vals[i]) * bound;

         if( matrix->transformstatus[colpos] == TRANSFORMSTATUS_NEG )
           (vals[i]) = -(vals[i]);

         assert(SCIPisFeasLE(scip, matrix->lhs[rows[i]], matrix->rhs[rows[i]]));
      }

      if( !SCIPisInfinity(scip, ub) && !SCIPisInfinity(scip, lb) )
         matrix->upperbounds[colpos] = ub - lb;
      else
         matrix->upperbounds[colpos] = SCIPinfinity(scip);

      matrix->transformshiftvals[colpos] = bound;
   }
   SCIPdebugMessage("Variable <%s> at colpos %d transformed. LB <%g> --> <%g>, UB <%g> --> <%g>\n",
      SCIPvarGetName(var), colpos, lb, 0.0, ub, matrix->upperbounds[colpos]);
}

/** initializes copy of the original coefficient matrix and applies heuristic specific adjustments: normalizing row
 *  vectors, transforming variable domains such that lower bound is zero, and relaxing continuous variables.
 */
static
SCIP_RETCODE initMatrix(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object to be initialized */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int*                  nvarsleftinrow,     /**< numbers of shiftable discrete variables left of all rows */
   SCIP_Bool             relax,              /**< should continuous variables be relaxed from the problem? */
   SCIP_Bool*            initialized,        /**< was the initialization successful? */
   SCIP_Bool*            infeasible          /**< is the problem infeasible? */
   )
{
   SCIP_ROW** lprows;
   SCIP_COL** lpcols;
   int i;
   int j;
   int currentpointer;

   int nrows;
   int ncols;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(initialized!= NULL);
   assert(infeasible != NULL);

   SCIPdebugMessage("entering Matrix Initialization method of SHIFTANDPROPAGATE heuristic!\n");

   /* get LP row data; column data is already initialized in heurdata */
   SCIP_CALL( SCIPgetLPRowsData(scip, &lprows, &nrows) );
   lpcols = heurdata->lpcols;
   ncols = heurdata->nlpcols;

   matrix->nrows = nrows;
   matrix->nnonzs = 0;
   matrix->normalized = FALSE;
   matrix->ndiscvars = 0;

   /* count the number of nonzeros of the LP constraint matrix */
   for( j = 0; j < ncols; ++j )
   {
      assert(lpcols[j] != NULL);
      assert(SCIPcolGetLPPos(lpcols[j]) >= 0);

      if( SCIPcolIsIntegral(lpcols[j]) )
      {
         matrix->nnonzs += SCIPcolGetNLPNonz(lpcols[j]);
         ++matrix->ndiscvars;
      }
   }

   matrix->ncols = matrix->ndiscvars;

   if( matrix->nnonzs == 0 )
   {
      SCIPdebugMessage("No matrix entries - Terminating initialization of matrix.\n");

      *initialized = FALSE;

      return SCIP_OKAY;
   }

   /* allocate memory for the members of heuristic matrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatvals, matrix->nnonzs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatind, matrix->nnonzs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatvals, matrix->nnonzs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatind, matrix->nnonzs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatbegin, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatbegin, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->lhs, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rhs, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colnorms, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->transformstatus, matrix->ndiscvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->upperbounds, matrix->ndiscvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->transformshiftvals, matrix->ndiscvars) );

   /* set transform status of variables */
   for( j = 0; j < matrix->ndiscvars; ++j )
      matrix->transformstatus[j] = TRANSFORMSTATUS_NONE;

   currentpointer = 0;
   *infeasible = FALSE;

   /* initialize the rows vector of the heuristic matrix together with its corresponding
    * lhs, rhs.
    */
   for( i = 0; i < nrows; ++i )
   {
      SCIP_Real* rowvals;
      SCIP_COL** cols;
      SCIP_Real maxval;
      SCIP_ROW* row;
      SCIP_Real nrowlpnonz;
      SCIP_Real constant;

      /* get LP row information */
      row = lprows[i];
      rowvals = SCIProwGetVals(row);
      nrowlpnonz = SCIProwGetNLPNonz(row);
      maxval = SCIPgetRowMaxCoef(scip, row);
      cols = SCIProwGetCols(row);
      constant = SCIProwGetConstant(row);

      assert(!SCIPisInfinity(scip, constant));
      assert(SCIPisFeasGT(scip, maxval, 0.0) || nrowlpnonz == 0);

      matrix->rowmatbegin[i] = currentpointer;

      /* modify the lhs and rhs w.r.t to the rows constant and normalize by 1-norm, i.e divide the lhs and rhs by the
       * maximum absolute value of the row
       */
      if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
         matrix->lhs[i] = (SCIProwGetLhs(row) - constant);
      else
         matrix->lhs[i] = -SCIPinfinity(scip);

      if( !SCIPisInfinity(scip, SCIProwGetRhs(row)) )
         matrix->rhs[i] = (SCIProwGetRhs(row) - constant);
      else
         matrix->rhs[i] = SCIPinfinity(scip);

      /* make sure that maxval is larger than zero before normalization.
       * Maxval may be zero if the constraint contains no variables but is modifiable, hence not redundant
       */
      if( !SCIPisFeasZero(scip, maxval) )
      {
         if( !SCIPisInfinity(scip, -matrix->lhs[i]) )
            matrix->lhs[i] /= maxval;
         if( !SCIPisInfinity(scip, matrix->rhs[i]) )
            matrix->rhs[i] /= maxval;
      }

      SCIPdebugMessage(" %s : lhs=%g, rhs=%g, maxval=%g \n", SCIProwGetName(row), matrix->lhs[i], matrix->rhs[i], maxval);

      /* in case of empty rows with a 0 < lhs <= 0.0 or 0.0 <= rhs < 0 we deduce the infeasibility of the problem */
      if( nrowlpnonz == 0 && (SCIPisFeasPositive(scip, matrix->lhs[i]) || SCIPisFeasNegative(scip, matrix->rhs[i])) )
      {
         *infeasible = TRUE;
         SCIPdebugMessage("  Matrix initialization stopped because of row infeasibility! \n");
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
         break;
      }

      /* row coefficients are normalized and copied to heuristic matrix */
      for( j = 0; j < nrowlpnonz; ++j )
      {
         if( !SCIPcolIsIntegral(cols[j]) )
            continue;
         assert(SCIPcolGetLPPos(cols[j]) >= 0);
         assert(currentpointer < matrix->nnonzs);

         matrix->rowmatvals[currentpointer] = rowvals[j]/maxval;
         matrix->rowmatind[currentpointer] = SCIPcolGetLPPos(cols[j]);

         ++currentpointer;
         ++nvarsleftinrow[i];
      }
   }

   matrix->normalized = TRUE;

   if( *infeasible )
      return SCIP_OKAY;

   assert(currentpointer == matrix->nnonzs);

   currentpointer = 0;

   /* copy the nonzero coefficient data column by column to heuristic matrix */
   for( j = 0; j < matrix->ncols; ++j )
   {
      SCIP_COL* currentcol;
      SCIP_ROW** rows;
      SCIP_Real* colvals;
      int ncolnonz;


      assert(SCIPcolGetLPPos(lpcols[j]) >= 0);

      currentcol = lpcols[j];
      assert(SCIPcolIsIntegral(currentcol));

      colvals = SCIPcolGetVals(currentcol);
      rows = SCIPcolGetRows(currentcol);
      ncolnonz = SCIPcolGetNLPNonz(currentcol);
      matrix->colnorms[j] = ncolnonz;

      /* loop over all rows with nonzero coefficients in the column, transform them and add them to the heuristic matrix */
      matrix->colmatbegin[j] = currentpointer;
      for( i = 0; i < ncolnonz; ++i )
      {
         SCIP_Real maxval;

         assert(rows[i] != NULL);
         assert(0 <= SCIProwGetLPPos(rows[i]));
         assert(SCIProwGetLPPos(rows[i]) < nrows);
         assert(currentpointer < matrix->nnonzs);

         /* rows are normalized by maximum norm */
         maxval = SCIPgetRowMaxCoef(scip, rows[i]);

         assert(maxval > 0);

         matrix->colmatvals[currentpointer] = colvals[i]/maxval;
         matrix->colmatind[currentpointer] = SCIProwGetLPPos(rows[i]);

         /* update the column norm */
         matrix->colnorms[j] += ABS(matrix->colmatvals[currentpointer]);

         ++currentpointer;
      }
   }
   assert(currentpointer == matrix->nnonzs);

   /* each variable is either transformed, if it supposed to be integral, or relaxed */
   for( j = 0; j < (relax ? ncols : matrix->ndiscvars); ++j )
   {
      SCIP_COL* col;

      col = lpcols[j];
      if( SCIPcolIsIntegral(col) )
         transformVariable(scip, matrix, heurdata, j);
      else
      {
         SCIP_VAR* var;
         var = SCIPcolGetVar(col);
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);
         relaxVar(scip, var, matrix);
      }
   }
   *initialized = TRUE;

   SCIPdebugMessage("Matrix initialized for %d discrete variables with %d cols, %d rows and %d nonzero entries\n",
      matrix->ndiscvars, matrix->ncols, matrix->nrows, matrix->nnonzs);
   return SCIP_OKAY;
}

/** frees all members of the heuristic matrix */
static
void freeMatrix(
   SCIP*                 scip,               /**< current SCIP instance */
   CONSTRAINTMATRIX**    matrix              /**< constraint matrix object */
   )
{
   assert(scip != NULL);
   assert(matrix != NULL);

   /* all fields are only allocated, if problem is not empty  */
   if( (*matrix)->nnonzs > 0 )
   {
      assert((*matrix) != NULL);
      assert((*matrix)->rowmatbegin != NULL);
      assert((*matrix)->rowmatvals != NULL);
      assert((*matrix)->rowmatind != NULL);
      assert((*matrix)->colmatbegin != NULL);
      assert((*matrix)->colmatvals!= NULL);
      assert((*matrix)->colmatind != NULL);
      assert((*matrix)->lhs != NULL);
      assert((*matrix)->rhs != NULL);
      assert((*matrix)->transformstatus != NULL);
      assert((*matrix)->transformshiftvals != NULL);

      /* free all fields */
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatbegin));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatvals));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatind));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatvals));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatind));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatbegin));
      SCIPfreeBufferArray(scip, &((*matrix)->lhs));
      SCIPfreeBufferArray(scip, &((*matrix)->rhs));
      SCIPfreeBufferArray(scip, &((*matrix)->colnorms));
      SCIPfreeBufferArray(scip, &((*matrix)->transformstatus));
      SCIPfreeBufferArray(scip, &((*matrix)->upperbounds));
      SCIPfreeBufferArray(scip, &((*matrix)->transformshiftvals));

     (*matrix)->nrows = 0;
     (*matrix)->ncols = 0;
   }

   /* free matrix */
   SCIPfreeBuffer(scip, matrix);
}

/** collects the necessary information about row violations for the zero-solution. That is,
 *  all solution values in heuristic transformation are zero.
 */
static
void checkViolations(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int*                  violatedrows,       /**< violated rows */
   int*                  violatedrowpos,     /**< row positions of violated rows */
   int*                  nviolatedrows,      /**< pointer to store the number of violated rows */
   int*                  nredundantrows      /**< pointer to store the number of redundant rows */
   )
{
   SCIP_Real* rhs;
   SCIP_Real* lhs;
   int nrows;
   int i;

   assert(matrix != NULL);
   assert(violatedrows != NULL);
   assert(violatedrowpos != NULL);
   assert(nviolatedrows != NULL);

   /* get RHS, LHS and number of the problem rows */
   rhs = matrix->rhs;
   lhs = matrix->lhs;
   nrows = matrix->nrows;

   SCIPdebugMessage("Entering violation check for %d rows! \n", nrows);
   *nviolatedrows = 0;
   if( nredundantrows != NULL )
      *nredundantrows = 0;

   /* loop over rows and check if it is violated */
   for( i = 0; i < nrows; ++i )
   {
      /* check, if zero solution violates this row */
      if( SCIPisFeasLT(scip, rhs[i], 0.0) || SCIPisFeasLT(scip, -lhs[i], 0.0) )
      {
         violatedrows[*nviolatedrows] = i;
        (violatedrowpos)[i] = *nviolatedrows;
         ++(*nviolatedrows);
      }
      else
         violatedrowpos[i] = -1;

      assert((violatedrowpos[i] == -1 && SCIPisFeasGE(scip, rhs[i], 0.0) && SCIPisFeasGE(scip, -lhs[i], 0.0))
         || (violatedrowpos[i] >= 0 &&(SCIPisFeasLT(scip, rhs[i], 0.0) || SCIPisFeasLT(scip, -lhs[i], 0.0))));

      if( SCIPisInfinity(scip, rhs[i]) && SCIPisInfinity(scip, -lhs[i]) && nredundantrows != NULL)
         ++(*nredundantrows);
   }
}

/** retransforms solution values of variables according to their transformation status */
static
SCIP_Real retransformVariable(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   SCIP_VAR*             var,                /**< variable whose solution value has to be retransformed */
   int                   varindex,           /**< permutation of variable indices according to sorting */
   SCIP_Real             solvalue            /**< solution value of the variable */
   )
{
   TRANSFORMSTATUS status;

   assert(matrix != NULL);
   assert(var != NULL);

   status = matrix->transformstatus[varindex];
   assert(status != TRANSFORMSTATUS_NONE);
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);

   /* check if original variable has different bounds and transform solution value correspondingly */
   if( status == TRANSFORMSTATUS_LB )
   {
      assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));

      return solvalue += matrix->transformshiftvals[varindex];
   }
   else if( status == TRANSFORMSTATUS_NEG )
   {
      assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
      return matrix->transformshiftvals[varindex] - solvalue;
   }
   return solvalue;
}

/** determines the best shifting value of a variable
 *  @todo if there is already an incumbent solution, try considering the objective cutoff as additional constraint */
static
SCIP_RETCODE getOptimalShiftingValue(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int                   varindex,           /**< index of variable which should be shifted */
   int*                  rowweights,         /**< weighting of rows for best shift calculation */
   SCIP_Real*            beststep,           /**< pointer to store optimal shifting step */
   int*                  rowviolations       /**< pointer to store new weighted sum of row violations, i.e, v - f */
   )
{
   SCIP_Real* steps;
   SCIP_Real* vals;
   int* rows;
   int* violationchange;

   SCIP_Real slacksurplus;
   SCIP_Real upperbound;

   int nrows;
   int sum;
   int i;

   SCIP_Bool allzero;

   assert(beststep != NULL);
   assert(rowviolations != NULL);
   assert(rowweights != NULL);

   upperbound = matrix->upperbounds[varindex];

   /* get nonzero values and corresponding rows of variable */
   getColumnData(matrix, varindex, &vals, &rows, &nrows);

   /* allocate necessary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &steps, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &violationchange, nrows) );

   /* loop over rows and calculate, which is the minimum shift to make this row feasible
    * or the minimum shift to violate this row
    */
   allzero = TRUE;
   slacksurplus = 0.0;
   for( i = 0; i < nrows; ++i )
   {
      SCIP_Real lhs;
      SCIP_Real rhs;
      int rowpos;
      SCIP_Bool rowisviolated;
      int rowweight;

      /* get the row data */
      rowpos = rows[i];
      assert(rowpos >= 0);
      lhs = matrix->lhs[rowpos];
      rhs = matrix->rhs[rowpos];
      rowweight = rowweights[rowpos];
      assert(rowweight == DEFAULT_WEIGHT_INEQUALITY ||
        (!SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs)));

      /* determine if current row is violated or not */
      rowisviolated =(SCIPisFeasLT(scip, rhs, 0.0) || SCIPisFeasLT(scip, -lhs, 0.0));

      /* for a feasible row, determine the minimum integer value within the bounds of the variable by which it has to be
       * shifted to make row infeasible.
       */
      if( !rowisviolated )
      {
         SCIP_Real maxfeasshift;

         maxfeasshift = SCIPinfinity(scip);

         /* feasibility can only be violated if the variable has a lock in the corresponding direction,
          * i.e. a positive coefficient for a "<="-constraint, a negative coefficient for a ">="-constraint.
          */
         if( SCIPisFeasGT(scip, vals[i], 0.0) && !SCIPisInfinity(scip, rhs) )
            maxfeasshift = SCIPfeasFloor(scip, rhs/vals[i]);
         else if( SCIPisFeasLT(scip, vals[i], 0.0) && !SCIPisInfinity(scip, -lhs) )
            maxfeasshift = SCIPfeasFloor(scip, lhs/vals[i]);

         /* if the variable has no lock in the current row, it can still help to increase the slack of this row;
          * we measure slack increase for shifting by one
          */
         if( SCIPisFeasGT(scip, vals[i], 0.0) && SCIPisInfinity(scip, rhs) )
            slacksurplus += vals[i];
         if( SCIPisFeasLT(scip, vals[i], 0.0) && SCIPisInfinity(scip, -lhs) )
            slacksurplus -= vals[i];

         /* check if the least violating shift lies within variable bounds and set corresponding array values */
         if( SCIPisFeasLE(scip, maxfeasshift + 1.0, upperbound) )
         {
            steps[i] = maxfeasshift + 1.0;
            violationchange[i] = rowweight;
            allzero = FALSE;
         }
         else
         {
            steps[i] = upperbound;
            violationchange[i] = 0;
         }
      }
      /* for a violated row, determine the minimum integral value within the bounds of the variable by which it has to be
       * shifted to make row feasible.
       */
      else
      {
         SCIP_Real minfeasshift;

         minfeasshift = SCIPinfinity(scip);

         /* if coefficient has the right sign to make row feasible, determine the minimum integer to shift variable
          * to obtain feasibility
          */
         if( SCIPisFeasLT(scip, -lhs, 0.0) && SCIPisFeasGT(scip, vals[i], 0.0) )
            minfeasshift = SCIPfeasCeil(scip, lhs/vals[i]);
         else if( SCIPisFeasLT(scip, rhs,0.0) && SCIPisFeasLT(scip, vals[i], 0.0) )
            minfeasshift = SCIPfeasCeil(scip, rhs/vals[i]);

         /* check if the minimum feasibility recovery shift lies within variable bounds and set corresponding array
          * values
          */
         if( SCIPisFeasLE(scip, minfeasshift, upperbound) )
         {
            steps[i] = minfeasshift;
            violationchange[i] = -rowweight;
            allzero = FALSE;
         }
         else
         {
            steps[i] = upperbound;
            violationchange[i] = 0;
         }
      }
   }

   /* in case that the variable cannot affect the feasibility of any row, in particular it cannot violate
    * a single row, but we can add slack to already feasible rows, we will do this
    */
   if( allzero )
   {
      *beststep = SCIPisFeasGT(scip, slacksurplus, 0.0) ? upperbound : 0.0;
      goto TERMINATE;
   }

   /* sorts rows by increasing value of steps */
   SCIPsortRealInt(steps, violationchange, nrows);

   *beststep = 0.0;
   *rowviolations = 0;
   sum = 0;

   /* best shifting step is calculated by summing up the violation changes for each relevant step and
    * taking the one which leads to the minimum sum. This sum measures the balance of feasibility recovering and
    * violating changes which will be obtained by shifting the variable by this step
    * note, the sums for smaller steps have to be taken into account for all bigger steps, i.e., the sums can be
    * computed iteratively
    */
   for( i = 0; i < nrows && !SCIPisInfinity(scip, steps[i]); ++i )
   {
      sum += violationchange[i];

      /* if we reached the last entry for the current step value, we have finished computing its sum and
       * update the step defining the minimum sum
       */
      if( (i == nrows-1 || steps[i+1] > steps[i]) && sum < *rowviolations )
      {
         *rowviolations = sum;
         *beststep = steps[i];
      }
   }
   assert(*rowviolations <= 0);
   assert(!SCIPisInfinity(scip, *beststep));

   /* free temporary arrays */
 TERMINATE:
   SCIPfreeBufferArray(scip, &violationchange);
   SCIPfreeBufferArray(scip, &steps);

   return SCIP_OKAY;
}

/** updates the information about a row whenever violation status changes */
static
void updateViolations(
   SCIP*                 scip,               /**< current SCIP instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int                   rowindex,           /**< index of the row */
   int*                  violatedrows,       /**< contains all violated rows */
   int*                  violatedrowpos,     /**< positions of rows in the violatedrows array */
   int*                  nviolatedrows       /**< pointer to update total number of violated rows */
   )
{
   assert(matrix != NULL);
   assert(violatedrows != NULL);
   assert(violatedrowpos != NULL);
   assert(nviolatedrows != NULL);

   /* row is now violated. Enqueue it in the set of violated rows. */
   if( SCIPisFeasLT(scip, -(matrix->lhs[rowindex]), 0.0) || SCIPisFeasLT(scip, matrix->rhs[rowindex], 0.0) )
   {
      assert(violatedrowpos[rowindex] == -1);
      assert(*nviolatedrows < matrix->nrows);

      violatedrows[*nviolatedrows] = rowindex;
      violatedrowpos[rowindex] = *nviolatedrows;
      ++(*nviolatedrows);
   }
   /* row is now feasible. Remove it from the set of violated rows. */
   else
   {
      assert(violatedrowpos[rowindex] > -1);

      /* swap the row with last violated row */
      if( violatedrowpos[rowindex] != *nviolatedrows - 1 )
      {
         assert(*nviolatedrows - 1 >= 0);
         violatedrows[violatedrowpos[rowindex]] = violatedrows[*nviolatedrows - 1];
         violatedrowpos[violatedrows[*nviolatedrows - 1]] = violatedrowpos[rowindex];
      }

      /* unlink the row from its position in the array and decrease number of violated rows */
      violatedrowpos[rowindex] = -1;
      --(*nviolatedrows);
   }
}

/** updates transformation of a given variable by taking into account current local bounds. if the bounds have changed
 *  since last update, updating the heuristic specific upper bound of the variable, its current transformed solution value
 *  and all affected rows is necessary.
 */
static
void updateTransformation(
   SCIP*                 scip,               /**< current scip */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   varindex,           /**< index of variable in matrix */
   SCIP_Real*            transformshiftval,  /**< value, by which the variable has been shifted during transformation */
   SCIP_Real             lb,                 /**< local lower bound of the variable */
   SCIP_Real             ub,                 /**< local upper bound of the variable */
   int*                  violatedrows,       /**< violated rows */
   int*                  violatedrowpos,     /**< violated row positions */
   int*                  nviolatedrows,      /**< pointer to store number of violated rows */
   int*                  nvarsleftinrow,     /**< number of shiftable variables left in each row */
   SCIP_Bool*            infeasible          /**< is the current solution infeasible? */
   )
{
   TRANSFORMSTATUS status;
   SCIP_Real deltashift;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= varindex && varindex < matrix->ndiscvars);
   assert(nvarsleftinrow != NULL);
   assert(infeasible != NULL);

   /* deltashift is the difference between the old and new transformation value. */
   deltashift = 0.0;
   status = matrix->transformstatus[varindex];

   SCIPdebugMessage("  Variable <%d> [%g,%g], status %d(%g), ub %g \n", varindex, lb, ub, status,
      matrix->transformshiftvals[varindex], matrix->upperbounds[varindex]);

   /* depending on the variable status, deltashift is calculated differently. */
   if( status == TRANSFORMSTATUS_LB )
   {
      deltashift = lb - (*transformshiftval);
      *transformshiftval = lb;
      if( !SCIPisInfinity(scip, ub) )
         matrix->upperbounds[varindex] = ub - lb;
   }

   if( status == TRANSFORMSTATUS_NEG )
   {
      assert(!SCIPisInfinity(scip, ub));
      deltashift = (*transformshiftval) - ub;
      *transformshiftval = ub;

      if( !SCIPisInfinity(scip, -lb) )
         matrix->upperbounds[varindex] = ub - lb;
   }

   if( status == TRANSFORMSTATUS_FREE )
   {
      /* in case of a free transform status, if one of the bounds has become finite, we want
       * to transform this variable to a variable with a lowerbound or a negated transform status */
      if( !SCIPisInfinity(scip, -lb) || !SCIPisInfinity(scip, ub) )
      {
         transformVariable(scip, matrix, heurdata, varindex);

         /* violations have to be rechecked for all rows
          * todo : change this and only update violations of rows in which this variable
          *        appears
          */
         checkViolations(scip, matrix, violatedrows, violatedrowpos, nviolatedrows, NULL);

         assert(matrix->transformstatus[varindex] == TRANSFORMSTATUS_LB || TRANSFORMSTATUS_NEG);
         assert(SCIPisFeasLE(scip, ABS(lb), ABS(ub)) || matrix->transformstatus[varindex] == TRANSFORMSTATUS_NEG);
      }
   }

   /* if the bound, by which the variable was shifted, has changed, deltashift is larger than zero, which requires
    * an update of all affected rows
    */
   if( !SCIPisFeasZero(scip, deltashift) )
   {
      int i;
      int* rows;
      SCIP_Real* vals;
      int nrows;

      assert(SCIPisFeasGT(scip, deltashift, 0.0));

      /* get nonzero values and corresponding rows of variable */
      getColumnData(matrix, varindex, &vals, &rows, &nrows);

      /* go through rows, update the rows w.r.t. the influence of the changed transformation of the variable */
      for( i = 0; i < nrows && !(*infeasible); ++i )
      {
         SCIP_Bool updaterow;
         SCIP_Bool leftviolation;
         SCIP_Bool rightviolation;

         SCIPdebugMessage("  update slacks of row<%d>:  coefficient <%g>, %g <= 0 <= %g \n",
            rows[i], vals[i], matrix->lhs[rows[i]], matrix->rhs[rows[i]]);

         /* the row has to be updated if either lhs or rhs changes its sign. */
         leftviolation = SCIPisFeasLT(scip, -(matrix->lhs[rows[i]]), 0.0);

         if( !SCIPisInfinity(scip, -(matrix->lhs[rows[i]])) )
            matrix->lhs[rows[i]] -= (vals[i]) * deltashift;

         updaterow = leftviolation != SCIPisFeasLT(scip, -(matrix->lhs[rows[i]]), 0.0);

         rightviolation = SCIPisFeasLT(scip,(matrix->rhs[rows[i]]), 0.0);
         if( !SCIPisInfinity(scip, matrix->rhs[rows[i]]) )
            matrix->rhs[rows[i]] -= (vals[i]) * deltashift;

         updaterow = updaterow != (rightviolation != SCIPisFeasLT(scip,(matrix->rhs[rows[i]]), 0.0));

         /* update the row violation */
         if( updaterow )
            updateViolations(scip, matrix, rows[i], violatedrows, violatedrowpos, nviolatedrows);

         SCIPdebugMessage("             -->                           %g <= 0 <= %g %s\n",
            matrix->lhs[rows[i]], matrix->rhs[rows[i]], updaterow ? ": row violation updated " : "");

         if( matrix->upperbounds[varindex] == 0.0 )
            --nvarsleftinrow[rows[i]];

         if( nvarsleftinrow[rows[i]] == 0 && violatedrowpos[rows[i]] != -1 )
            *infeasible = TRUE;
         SCIPdebugMessage("    %d variables left in %s row <%d>\n", nvarsleftinrow[rows[i]],
            violatedrowpos[rows[i]] >= 0 ? "infeasible":"feasible",
            rows[i]);
      }
   }

   SCIPdebugMessage("  Variable <%d> [%g,%g], status %d(%g), ub %g \n", varindex, lb, ub, status,
      matrix->transformshiftvals[varindex], matrix->upperbounds[varindex]);
}

/** updates heuristic transformed problem and row activities according to a shift of the variables value */
static
void shiftVariable(
   SCIP*                 scip,               /**< current scip */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int                   varindex,           /**< index of the variable to be shifted in heuristic transformed problem */
   SCIP_Real             shiftvalue,         /**< shifting value */
   int*                  violatedrows,       /**< set of currently violated rows */
   int*                  violatedrowpos,     /**< violated row positions */
   int*                  nviolatedrows,      /**< pointer to store number of violated rows */
   int*                  nvarsleftinrow,     /**< number of shiftable variables left in every row */
   SCIP_Bool*            infeasible          /**< is infeasibility of current solution not recoverable? */
   )
{
   SCIP_Real* vals;
   int* rows;
   int nrows;
   int i;

   assert(varindex < matrix->ndiscvars);
   assert(!(*infeasible));

   /* get nonzero values and corresponding rows of column */
   getColumnData(matrix, varindex, &vals, &rows, &nrows);

   /* go through nonzero rows, update their lhs and rhs and the position in the array of violated rows, if necessary;
    * a row is violated if the zero solution is not feasible
    */
   for( i = 0; i < nrows && !(*infeasible); ++i )
   {
      SCIP_Real oldlhs;
      SCIP_Real oldrhs;
      SCIP_Bool updatelhs;
      SCIP_Bool updaterhs;

      updatelhs = FALSE;
      updaterhs = FALSE;

      oldlhs = matrix->lhs[rows[i]];
      oldrhs = matrix->rhs[rows[i]];
     
      /* perform the shift, i.e., update lhs and rhs and check whether violation status changes */
      if( !SCIPisInfinity(scip, -oldlhs) )
      {
         matrix->lhs[rows[i]] -= vals[i] * shiftvalue;
         updatelhs = (SCIPisFeasLT(scip, -matrix->lhs[rows[i]], 0.0) != SCIPisFeasLT(scip, -oldlhs, 0.0));
      }
      if( !SCIPisInfinity(scip, oldrhs) )
      {
         matrix->rhs[rows[i]] -= vals[i] * shiftvalue;
         updaterhs = SCIPisFeasLT(scip, matrix->rhs[rows[i]], 0.0) != SCIPisFeasLT(scip, oldrhs, 0.0);
      }

      /* update violated row information */
      if( updatelhs != updaterhs )
         updateViolations(scip, matrix, rows[i], violatedrows, violatedrowpos, nviolatedrows);

      /* decrease number of variables with solution value not set in row */
      assert(nvarsleftinrow[rows[i]] > 0);
      --nvarsleftinrow[rows[i]];

      /* check if infeasibility of row cannot be resolved anymore */
      if( nvarsleftinrow[rows[i]] == 0 && violatedrowpos[rows[i]] >= 0 )
         *infeasible = TRUE;

      SCIPdebugMessage("    %d variables left in %s row <%d>\n", nvarsleftinrow[rows[i]],
         violatedrowpos[rows[i]] >= 0 ? "infeasible":"feasible",
         rows[i]);
   }
}

/** comparison method for columns; binary < integer < implicit < continuous variables */
static
SCIP_DECL_SORTPTRCOMP(heurSortColsShiftandpropagate)
{
   SCIP_COL* col1;
   SCIP_COL* col2;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VARTYPE vartype1;
   SCIP_VARTYPE vartype2;
   int key1;
   int key2;

   col1 = (SCIP_COL*)elem1;
   col2 = (SCIP_COL*)elem2;
   var1 = SCIPcolGetVar(col1);
   var2 = SCIPcolGetVar(col2);
   assert(var1 != NULL && var2 != NULL);

   vartype1 = SCIPvarGetType(var1);
   vartype2 = SCIPvarGetType(var2);

   switch (vartype1) {
      case SCIP_VARTYPE_BINARY:
         key1 = 1;
         break;
      case SCIP_VARTYPE_INTEGER:
         key1 = 2;
         break;
      case SCIP_VARTYPE_IMPLINT:
         key1 = 3;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         key1 = 4;
         break;
      default:
         key1 = -1;
         SCIPerrorMessage("unknown variable type\n");
         SCIPABORT();
         break;
   }
   switch (vartype2) {
      case SCIP_VARTYPE_BINARY:
         key2 = 1;
         break;
      case SCIP_VARTYPE_INTEGER:
         key2 = 2;
         break;
      case SCIP_VARTYPE_IMPLINT:
         key2 = 3;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         key2 = 4;
         break;
      default:
         key2 = -1;
         SCIPerrorMessage("unknown variable type\n");
         SCIPABORT();
         break;
   }
   return key1 - key2;
}

/*
 * Callback methods of primal heuristic
 */

/** deinitialization method of primal heuristic(called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitShiftandpropagate)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   assert(heurdata->lpcols == NULL || heurdata->lpcolssize > 0);

   if( heurdata->lpcolssize > 0 )
   {
      assert(heurdata->lpcols != NULL);
      SCIPfreeBlockMemoryArray(scip, &heurdata->lpcols, heurdata->lpcolssize);
      heurdata->lpcols = NULL;
      heurdata->lpcolssize = -1;
      heurdata->nlpcols = 0;
   }

   /* if statistic mode is enabled, statistics are printed to console */
   SCIPstatistic(
      SCIPstatisticMessage(
         "  DETAILS                    :  %d violations left, %d probing status, %d redundant rows\n",
         heurdata->nremainingviols,
         heurdata->lpsolstat,
         heurdata->nredundantrows);
      SCIPstatisticMessage(
         "  SHIFTANDPROPAGATE PROBING  :  %d probings, %lld domain reductions,  ncutoffs: %d ,  LP iterations: %lld \n ",
         heurdata->nprobings,
         heurdata->ntotaldomredsfound,
         heurdata->ncutoffs,
         heurdata->nlpiters);
      );

   return SCIP_OKAY;
}

/** initialization method of primal heuristic(called after problem was transformed). We only need this method for
 *  statistic mode of heuristic.
 */
static
SCIP_DECL_HEURINIT(heurInitShiftandpropagate)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   heurdata->randseed = DEFAULT_RANDSEED;

   SCIPstatistic(
      heurdata->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
      heurdata->nremainingviols = 0;
      heurdata->nprobings = 0;
      heurdata->ntotaldomredsfound = 0;
      heurdata->ncutoffs = 0;
      heurdata->nlpiters = 0;
      heurdata->nredundantrows = 0;
   )
   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data(called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeShiftandpropagate)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);

   /* free heuristic data */
   if( heurdata != NULL )
      SCIPfreeMemory(scip, &heurdata);

   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** copy method for primal heuristic plugins(called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyShiftandpropagate)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurShiftandpropagate(scip) );

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecShiftandpropagate)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;       /* the heuristic data */

   CONSTRAINTMATRIX* matrix;      /* constraint matrix object */
   SCIP_COL** lpcols;             /* lp columns */
   SCIP_SOL* sol;                 /* solution pointer */
   SCIP_Real* colnorms;           /* contains Euclidean norms of column vectors */
   int* violatedrows;             /* the violated rows */
   int* violatedrowpos;           /* the array position of a violated row, or -1 */
   int* permutation;              /* reflects the position of the variables after sorting */
   int* rowweights;               /* weighting of rows for best shift calculation */
   int* nvarsleftinrow;           /* number of variables left to the row */
   int nlpcols;                   /* number of lp columns */
   int nviolatedrows;             /* number of violated rows */
   int ndiscvars;                 /* number of non-continuous variables of the problem */
   int lastindexofsusp;           /* last variable which has been swapped due to a cutoff */
   int nbinvars;                  /* number of binary variables */
   int nintvars;                  /* number of integer variables */
   int nimplvars;                 /* number of implicit variables */
   int i;
   int r;
   int v;
   int c;
   int ncutoffs;                  /* counts the number of cutoffs for this execution */
   int nprobings;                 /* counts the number of probings */
   int nredundantrows;            /* the number of redundant rows */
   int nlprows;                   /* the number LP rows */

   SCIP_Bool initialized;         /* has the matrix been initialized? */
   SCIP_Bool cutoff;              /* has current probing node been cutoff? */
   SCIP_Bool probing;             /* should probing be applied or not? */
   SCIP_Bool infeasible;          /* FALSE as long as currently infeasible rows have variables left */

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;
   SCIPdebugMessage("entering execution method of shift and propagate heuristic\n");

   /* heuristic is obsolete if there are only continuous variables */
   if( SCIPgetNVars(scip) - SCIPgetNContVars(scip) == 0 )
      return SCIP_OKAY;

   /* stop execution method if there is already a primarily feasible solution at hand */
   if( SCIPgetBestSol(scip) != NULL && heurdata->onlywithoutsol )
      return SCIP_OKAY;

   /* stop if there is no LP available */
   if ( ! SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;

   if( !SCIPisLPConstructed(scip) )
   {
      SCIP_Bool nodecutoff;

      nodecutoff = FALSE;
      /* @note this call can have the side effect that variables are created */
      SCIP_CALL( SCIPconstructLP(scip, &nodecutoff) );
      SCIP_CALL( SCIPflushLP(scip) );
   }

   SCIPstatistic( heurdata->nlpiters = SCIPgetNLPIterations(scip) );

   nlprows = SCIPgetNLPRows(scip);

   SCIP_CALL( SCIPgetLPColsData(scip, &lpcols, &nlpcols) );
   assert(nlpcols == 0 || lpcols != NULL);

   /* we need an LP */
   if( nlprows == 0 || nlpcols == 0 )
      return SCIP_OKAY;


   *result = SCIP_DIDNOTFIND;
   initialized = FALSE;

   if( heurdata->lpcols == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(heurdata->lpcols), nlpcols) );
      heurdata->lpcolssize = nlpcols;
   }
   if( heurdata->lpcolssize < nlpcols )
   {
      assert(heurdata->lpcols != NULL);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(heurdata->lpcols), heurdata->lpcolssize, nlpcols) );
      heurdata->lpcolssize = nlpcols;
   }
   assert(heurdata->lpcolssize >= nlpcols);
   heurdata->nlpcols = nlpcols;
#ifndef NDEBUG
   BMSclearMemoryArray(heurdata->lpcols, heurdata->lpcolssize);
#endif

   BMScopyMemoryArray(heurdata->lpcols, lpcols, nlpcols);

   SCIPsortPtr((void**)heurdata->lpcols, heurSortColsShiftandpropagate, nlpcols);

   /* we have to collect the number of different variable types before we start probing since during probing variable
    * can be created(e.g., cons_xor.c)
    */
   ndiscvars = 0;
   nbinvars = 0;
   nintvars = 0;
   for( c = 0; c < nlpcols; ++c )
   {
      SCIP_COL* col;
      SCIP_VAR* colvar;

      col = heurdata->lpcols[c];
      assert(col != NULL);
      colvar =SCIPcolGetVar(col);
      assert(colvar != NULL);

      if( SCIPvarGetType(colvar) != SCIP_VARTYPE_CONTINUOUS )
         ++ndiscvars;
      if( SCIPvarGetType(colvar) == SCIP_VARTYPE_BINARY )
         ++nbinvars;
      else if( SCIPvarGetType(colvar) == SCIP_VARTYPE_INTEGER )
         ++nintvars;
   }
   assert(nbinvars + nintvars <= ndiscvars);

   /* start probing mode */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* enables collection of variable statistics during probing */
   SCIPenableVarHistory(scip);

   SCIP_CALL( SCIPnewProbingNode(scip) );
   ncutoffs = 0;
   nprobings = 0;
   infeasible = FALSE;

   /* initialize heuristic matrix and working solution */
   SCIP_CALL( SCIPallocBuffer(scip, &matrix) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nvarsleftinrow, nlprows) );
   BMSclearMemoryArray(nvarsleftinrow, nlprows);
   SCIP_CALL( initMatrix(scip, matrix, heurdata, nvarsleftinrow, heurdata->relax, &initialized, &infeasible) );

   /* could not initialize matrix */
   if( !initialized || infeasible )
   {
      SCIPdebugMessage(" MATRIX not initialized -> Execution of heuristic stopped! \n");
      goto TERMINATE;
   }

   /* the number of discrete LP column variables can be less than the actual number of variables, if, e.g., there
    * are nonlinearities in the problem. The heuristic execution can be terminated in that case.
    */
   if( matrix->ndiscvars < ndiscvars )
   {
      SCIPdebugMessage(" Not all discrete variables are in the current LP. Shiftandpropagate execution terminated\n");
      goto TERMINATE;
   }

   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
   SCIPsolSetHeur(sol, heur);

   /* allocate arrays for execution method */
   SCIP_CALL( SCIPallocBufferArray(scip, &permutation, ndiscvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowweights, matrix->nrows) );

   /* initialize arrays. Before sorting, permutation is the identity permutation */
   for( i = 0; i < ndiscvars; ++i )
      permutation[i] = i;

   /* initialize row weights */
   for( r = 0; r < matrix->nrows; ++r )
   {
      if( !SCIPisInfinity(scip, -(matrix->lhs[r])) && !SCIPisInfinity(scip, matrix->rhs[r]) )
         rowweights[r] = DEFAULT_WEIGHT_EQUALITY;
      else
         rowweights[r] = DEFAULT_WEIGHT_INEQUALITY;

   }
   colnorms = matrix->colnorms;

   nimplvars = ndiscvars - (nbinvars + nintvars);
   assert(nbinvars >= 0);
   assert(nintvars >= 0);
   assert(nimplvars >= 0);

   /* sort variables w.r.t. the sorting key parameter. Sorting is indirect, all matrix column data
    * stays in place, but permutation array gives access to the sorted order of variables
    */
   if( heurdata->sortvars && heurdata->sortkey == 'n' )
   {
      if( nbinvars > 0 )
         SCIPsortDownRealInt(colnorms, permutation, nbinvars);
      if( nintvars > 0 )
         SCIPsortDownRealInt(&colnorms[nbinvars], &permutation[nbinvars], nintvars);
      if( nimplvars > 0 )
         SCIPsortDownRealInt(&colnorms[nbinvars + nintvars], &permutation[nbinvars + nintvars], nimplvars);
      SCIPdebugMessage("Variables sorted down w.r.t their normalized columns!\n");
   }
   else if( heurdata->sortvars && heurdata->sortkey == 'u')
   {
      if( nbinvars > 0 )
         SCIPsortRealInt(colnorms, permutation, nbinvars);
      if( nintvars > 0 )
         SCIPsortRealInt(&colnorms[nbinvars], &permutation[nbinvars], nintvars);
      if( nimplvars > 0 )
         SCIPsortRealInt(&colnorms[nbinvars + nintvars], &permutation[nbinvars + nintvars], nimplvars);
      SCIPdebugMessage("Variables sorted w.r.t their normalized columns!\n");
   }
   else if( heurdata->sortvars && heurdata->sortkey == 'r')
   {
      if( nbinvars > 0)
         SCIPpermuteIntArray(permutation, 0, nbinvars - 1, &heurdata->randseed);
      if( nintvars > 0)
         SCIPpermuteIntArray(&permutation[nbinvars], 0, nintvars, &heurdata->randseed);
      if( nimplvars > 0)
         SCIPpermuteIntArray(&permutation[nbinvars + nintvars], 0, nimplvars, &heurdata->randseed);
      SCIPdebugMessage("Variables permuted randomly!\n");
   }
   else
   {
      SCIPdebugMessage("No variable permutation applied\n");
   }

   /* allocate arrays to store information about infeasible rows */
   SCIP_CALL( SCIPallocBufferArray(scip, &violatedrows, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &violatedrowpos, matrix->nrows) );

   /* check rows for infeasibility */
   nredundantrows = 0;
   checkViolations(scip, matrix, violatedrows, violatedrowpos, &nviolatedrows, &nredundantrows);

   cutoff = FALSE;
   lastindexofsusp = -1;
   probing = heurdata->probing;
   infeasible = FALSE;

   SCIPdebugMessage("SHIFT_AND_PROPAGATE heuristic starts main loop with %d violations and %d remaining variables!\n",
      nviolatedrows, ndiscvars);

   assert(matrix->ndiscvars == ndiscvars);

   /* loop over variables, shift them according to shifting criteria and try to reduce the global infeasibility */
   for( c = 0; c < ndiscvars; ++c )
   {
      SCIP_VAR* var;
      SCIP_Longint ndomredsfound;
      SCIP_Real optimalshiftvalue;
      SCIP_Real origsolval;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real obj;
      TRANSFORMSTATUS status;
      int nviolations;
      int permutedvarindex;

      permutedvarindex = permutation[c];
      optimalshiftvalue = 0.0;
      nviolations = 0;
      var = SCIPcolGetVar(heurdata->lpcols[permutedvarindex]);
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      assert(SCIPcolGetLPPos(SCIPvarGetCol(var)) >= 0);
      assert(SCIPvarIsIntegral(var));

      /* check whether we hit some limit, e.g. the time limit, in between
       * since the check itself consumes some time, we only do it every tenth iteration
       */
      if( c % 10 == 0 && SCIPisStopped(scip) )
         goto TERMINATE2;

      /* if propagation is enabled, check if propagation has changed the variables bounds
       * and update the transformed upper bound correspondingly
       */
      if( heurdata->probing )
         updateTransformation(scip, matrix, heurdata, permutedvarindex, &(matrix->transformshiftvals[permutedvarindex]),
            lb, ub, violatedrows, violatedrowpos, &nviolatedrows, nvarsleftinrow, &infeasible);

      /* if the problem cannot be solved anymore by heuristic, break */
      if( infeasible )
         break;

      status = matrix->transformstatus[permutedvarindex];
      obj = SCIPvarGetObj(var);
      if( status == TRANSFORMSTATUS_NEG )
         obj = -obj;

      SCIPdebugMessage("Variable %s with local bounds [%g,%g], status <%d>, matrix bound <%g>\n",
         SCIPvarGetName(var), lb, ub, status, matrix->upperbounds[permutedvarindex]);

      /* ignore variable if propagation fixed it(lb and ub will be zero) */
      if( SCIPisFeasZero(scip, matrix->upperbounds[permutedvarindex]) )
      {
         assert(!SCIPisInfinity(scip, ub));
         assert(SCIPisFeasEQ(scip, lb, ub));

         SCIP_CALL( SCIPsetSolVal(scip, sol, var, ub) );

         continue;
      }

      /* Variables with FREE transform status are currently not dealt with */
      /* @todo change matrix representation and shiftval methods to treat FREE variables correctly */
      if( matrix->transformstatus[permutedvarindex] == TRANSFORMSTATUS_FREE )
      {
         ++lastindexofsusp;
         assert(lastindexofsusp >= 0 && lastindexofsusp <= c);

         assert(permutedvarindex == permutation[c]);

         permutation[c] = permutation[lastindexofsusp];
         permutation[lastindexofsusp] = permutedvarindex;

         SCIPdebugMessage("  Variable %s postponed from pos <%d> to <%d> due to FREE transform status\n", SCIPvarGetName(var), c, lastindexofsusp);

         continue;
      }

      /* compute optimal shift value for variable */
      SCIP_CALL( getOptimalShiftingValue(scip, matrix, permutedvarindex, rowweights, &optimalshiftvalue, &nviolations) );
      assert(SCIPisFeasGE(scip, optimalshiftvalue, 0.0));
      assert(!SCIPisInfinity(scip, optimalshiftvalue));

      /* in case the problem is already feasible, do not shift in the direction which deteriorates current objective */
      if( nviolatedrows == 0 && SCIPisFeasGT(scip, obj, 0.0) )
         optimalshiftvalue = 0.0;

      /* retransform the solution value from the heuristic transformation space */
      origsolval = retransformVariable(scip, matrix, var, permutedvarindex, optimalshiftvalue);
      assert(SCIPisFeasGE(scip, origsolval, lb) && SCIPisFeasLE(scip, origsolval, ub));

      SCIPdebugMessage("  Shift %g(%g originally) is optimal, propagate solution\n", optimalshiftvalue, origsolval);

      /* check if propagation should still be performed */
      if( nprobings > DEFAULT_PROPBREAKER )
         probing = FALSE;

      /* if propagation is enabled, fix the variable to the new solution value and propagate the fixation
       * (to fix other variables and to find out early whether solution is already infeasible)
       */
      if( probing )
      {
         SCIP_CALL( SCIPnewProbingNode(scip) );
         SCIP_CALL( SCIPfixVarProbing(scip, var, origsolval) );
         ndomredsfound = 0;

         SCIP_CALL( SCIPpropagateProbing(scip, heurdata->nproprounds, &cutoff, &ndomredsfound) );

         ++nprobings;
         SCIPstatistic( heurdata->ntotaldomredsfound += ndomredsfound );
         SCIPdebugMessage("Propagation finished! <%lld> domain reductions %s, <%d> probing depth\n", ndomredsfound, cutoff ? "CUTOFF" : "",
            SCIPgetProbingDepth(scip));
      }
      assert(!cutoff || probing);

      /* propagation led to an empty domain, hence we backtrack and postpone the variable */
      if( cutoff )
      {
         assert(probing);

         ++ncutoffs;
         ++lastindexofsusp;
         assert(lastindexofsusp >= 0 && lastindexofsusp <= c);

         /* only continue heuristic if number of cutoffs occured so far is reasonably small */
         if( heurdata->cutoffbreaker >= 0 && ncutoffs >= heurdata->cutoffbreaker )
            break;

         cutoff = FALSE;

         assert(permutedvarindex == permutation[c]);

         permutation[c] = permutation[lastindexofsusp];
         permutation[lastindexofsusp] = permutedvarindex;

         /* backtrack to the parent of the current node */
         assert(SCIPgetProbingDepth(scip) >= 1);
         SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip) - 1) );
         SCIPdebugMessage("  Suspicious variable! Postponed from pos <%d> to position <%d>\n", c, lastindexofsusp);
      }
      else
      {
         /* update heuristic transformed problem according to shifted value of variable */
         shiftVariable(scip, matrix, permutedvarindex, optimalshiftvalue, violatedrows, violatedrowpos, &nviolatedrows,
            nvarsleftinrow, &infeasible);

         if( infeasible )
         {
            SCIPdebugMessage(" Heuristic stopped because of fundamental infeasibility in row \n");
#ifdef SCIP_DEBUG
            for( i = 0; i < matrix->nrows; ++i )
            {
               if( nvarsleftinrow[i] == 0 && violatedrowpos[i] >= 0 )
	       {
                  SCIPdebugMessage(" <%d> \n", i);
	       }
            }
#endif
            break;
         }
         else
         {
            SCIPdebugMessage("Variable <%d><%s> successfully shifted by value <%g>!\n", permutedvarindex,
               SCIPvarGetName(var), optimalshiftvalue);
         }
         /* update solution */
         SCIP_CALL( SCIPsetSolVal(scip, sol, var, origsolval) );

         /* only to ensure that some assertions can be made later on */
         if( !probing )
         {
            SCIP_CALL( SCIPfixVarProbing(scip, var, origsolval) );
         }
      }
   }
   SCIPdebugMessage("Heuristic finished with %d remaining violations and %d remaining variables!\n",
      nviolatedrows, lastindexofsusp + 1);

   /* if constructed solution might be feasible, go through the queue of suspicious variables and set the solution
    * values
    */
   if( nviolatedrows == 0 && !cutoff )
   {
      SCIP_Bool stored;

      for( v = 0; v <= lastindexofsusp; ++v )
      {
         SCIP_VAR* var;
         SCIP_Real origsolval;
         int permutedvarindex;

         /* get the column position of the variable */
         permutedvarindex = permutation[v];
         var = SCIPcolGetVar(heurdata->lpcols[permutedvarindex]);

         assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);

         /* update the transformation of the variable, since the bound might have changed after the last update. */
         if( heurdata->probing )
            updateTransformation(scip, matrix, heurdata, permutedvarindex, &(matrix->transformshiftvals[permutedvarindex]),
               SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), violatedrows, violatedrowpos, &nviolatedrows,
               nvarsleftinrow, &infeasible);

         /* retransform the solution value from the heuristic transformed space, set the solution value accordingly */
         origsolval = retransformVariable(scip, matrix, var, permutedvarindex, 0.0);
         assert(SCIPisFeasGE(scip, origsolval, SCIPvarGetLbLocal(var))
            && SCIPisFeasLE(scip, origsolval, SCIPvarGetUbLocal(var)));
         SCIP_CALL( SCIPsetSolVal(scip, sol, var, origsolval) );
         SCIP_CALL( SCIPfixVarProbing(scip, var, origsolval) ); /* only to ensure that some assertions can be made later */

         SCIPdebugMessage("  Remaining variable <%s> set to <%g>; %d Violations\n", SCIPvarGetName(var), origsolval,
            nviolatedrows);
      }
      /* Fixing of remaining variables led to infeasibility */
      if( nviolatedrows > 0 )
         goto TERMINATE2;

      stored = TRUE;
      /* if the constructed solution might still be extendable to a feasible solution, try this by
       * solving the remaining LP
       */
      if( nlpcols != matrix->ndiscvars )
      {
         /* case that remaining LP has to be solved */
         SCIP_Bool lperror;

#ifndef NDEBUG
         {
            SCIP_VAR** vars;

            vars = SCIPgetVars(scip);
            assert(vars != NULL);
            /* ensure that all discrete variables in the remaining LP are fixed */
            for( v = 0; v < ndiscvars; ++v )
            {
               if( SCIPvarIsInLP(vars[v]) )
                  assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(vars[v]), SCIPvarGetUbLocal(vars[v])));

            }
         }
#endif

         SCIPdebugMessage(" -> old LP iterations: %"SCIP_LONGINT_FORMAT"\n", SCIPgetNLPIterations(scip));

         /* solve LP;
          * errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
          * hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
          */
#ifdef NDEBUG
         {
            SCIP_RETCODE retstat;
            retstat = SCIPsolveProbingLP(scip, -1, &lperror);
            if( retstat != SCIP_OKAY )
            {
               SCIPwarningMessage(scip, "Error while solving LP in SHIFTANDPROPAGATE heuristic; LP solve terminated with code <%d>\n",
                     retstat);
            }
         }
#else
         SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror) );
#endif

         SCIPdebugMessage(" -> new LP iterations: %"SCIP_LONGINT_FORMAT"\n", SCIPgetNLPIterations(scip));
         SCIPdebugMessage(" -> error=%u, status=%d\n", lperror, SCIPgetLPSolstat(scip));

         /* check if this is a feasible solution */
         if( !lperror && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
         {
            /* copy the current LP solution to the working solution */
            SCIP_CALL( SCIPlinkLPSol(scip, sol) );
         }
         else
            stored = FALSE;

         SCIPstatistic( heurdata->lpsolstat = SCIPgetLPSolstat(scip) );
      }
      /* check solution for feasibility, and add it to solution store if possible.
       * Neither integrality nor feasibility of LP rows have to be checked, because they
       * are guaranteed by the heuristic at this stage.
       */
      if( stored )
      {
#ifndef NDEBUG
         SCIP_CALL( SCIPtrySol(scip, sol, FALSE, TRUE, TRUE, TRUE, &stored) );
#else
         /* @todo: maybe bounds don't need to be checked, in this case put an assert concerning stored ?????????? */
         SCIP_CALL( SCIPtrySol(scip, sol, FALSE, TRUE, FALSE, FALSE, &stored) );
#endif
         if( stored )
         {
            SCIPdebugMessage("found feasible shifted solution:\n");
            SCIPdebug( SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) ) );
            *result = SCIP_FOUNDSOL;
            SCIPstatisticMessage("  Shiftandpropagate solution value: %16.9g \n", SCIPgetSolOrigObj(scip, sol));
         }
      }
   }
   else
   {
      SCIPdebugMessage("Solution constructed by heuristic is already known to be infeasible\n");
   }

   SCIPstatistic(
      heurdata->nremainingviols = nviolatedrows;
      heurdata->nredundantrows = nredundantrows;
      );

 TERMINATE2:
   /* free all allocated memory */
   SCIPfreeBufferArray(scip, &violatedrowpos);
   SCIPfreeBufferArray(scip, &violatedrows);
   SCIPfreeBufferArray(scip, &rowweights);
   SCIPfreeBufferArray(scip, &permutation);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

 TERMINATE:
   /* terminate probing mode and free the remaining memory */
   SCIPstatistic(
      heurdata->ncutoffs += ncutoffs;
      heurdata->nprobings += nprobings;
      heurdata->nlpiters = SCIPgetNLPIterations(scip) - heurdata->nlpiters;
      );

   SCIP_CALL( SCIPendProbing(scip) );
   freeMatrix(scip, &matrix);
   SCIPfreeBufferArray(scip, &nvarsleftinrow);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the shiftandpropagate primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurShiftandpropagate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Shiftandpropagate primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heurdata->lpcols = NULL;
   heurdata->lpcolssize = -1;
   heurdata->nlpcols = 0;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecShiftandpropagate, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyShiftandpropagate) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeShiftandpropagate) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitShiftandpropagate) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitShiftandpropagate) );

   /* add shiftandpropagate primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/nproprounds", "The number of propagation rounds used for each propagation",
         &heurdata->nproprounds, TRUE, DEFAULT_NPROPROUNDS, -1, 1000, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/relax", "Should continuous variables be relaxed?",
         &heurdata->relax, TRUE, DEFAULT_RELAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/probing", "Should domains be reduced by probing?",
         &heurdata->probing, TRUE, DEFAULT_PROBING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/onlywithoutsol", "Should heuristic only be executed if no primal solution was found, yet?",
         &heurdata->onlywithoutsol, TRUE, DEFAULT_ONLYWITHOUTSOL, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/cutoffbreaker", "The number of cutoffs before heuristic stops",
         &heurdata->cutoffbreaker, TRUE, DEFAULT_CUTOFFBREAKER, -1, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/"HEUR_NAME"/sortkey", "the key for variable sorting: (n)orms or (r)andom",
         &heurdata->sortkey, TRUE, DEFAULT_SORTKEY, SORTKEYS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/shiftandpropagate/sortvars", "Should variables be sorted for the heuristic?",
         &heurdata->sortvars, TRUE, DEFAULT_SORTVARS, NULL, NULL));
   return SCIP_OKAY;
}
