/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_twoopt.c
 * @ingroup PRIMALHEURISTICS
 * @brief  primal heuristic to improve incumbent solution by flipping pairs of variables
 * @author Timo Berthold 
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/heur_twoopt.h"

#define HEUR_NAME             "twoopt"
#define HEUR_DESC             "primal heuristic to improve incumbent solution by flipping pairs of variables"
#define HEUR_DISPCHAR         '2'
#define HEUR_PRIORITY         -20100
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1

#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE 

/* default parameter values */
#define DEFAULT_INTOPT                  TRUE /**< optional integer optimization is applied by default */
#define DEFAULT_WAITINGNODES               0 /**< default number of nodes to wait after current best solution before calling heuristic */
#define DEFAULT_MATCHINGRATE            0.65 /**< default percentage by which two variables have to match in their LP-row set to be 
                                              *   associated as pair by heuristic */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   lastsolindex;       /**< index of last solution for which heuristic was performed */
   SCIP_Real             matchingrate;       /**< percentage by which two variables have have to match in their LP-row 
                                              *   set to be associated as pair by heuristic */
   SCIP_VAR**            binvars;            /**< Array of binary variables which are sorted with respect to their occurence
                                              *   in the LP-rows */
   
   int                   nbinvars;           /**< number of binary variables stored in heuristic array */
  
   int                   waitingnodes;       /**< user parameter to determine number of nodes to wait after last best solution 
                                              *   before calling heuristic   */
   SCIP_Bool             nobinarymatching;   /**< set to true if heuristic cannot find proper constraint-matchings within binaries
                                              *   to avoid obsolete calls */ 
   SCIP_Bool             binequalcoeffs;     /**< heuristic should not be applied if all binary coefficients 
                                              *   in objective function are equal */
   SCIP_Bool             presolved;          /**< flag to indicate whether presolving has already been executed */
    
   /* integer variable twoopt data */
   SCIP_Bool             intopt;             /**< parameter to determine if integer 2-opt should be applied */             
   SCIP_Bool             nointmatching;      /**< set to true if heuristic cannot determine any proper matching for 
                                              *   Integer variables */
   SCIP_Bool             intequalcoeffs;     /**< heuristic should not be applied if all binary coefficients in objective 
                                              *   function are equal */
  
   SCIP_VAR**            intvars;            /**< array to store the integer variables in non-decreasing order 
                                              *   with respect to their objective coefficient */
   int                   nintvars;           /**< the number of integer variables stored in array intvars */
#ifdef STATISTIC_INFORMATION
   /* statistics */
   int                   ntotalbinvars;      /**< total number of binary variables over all runs */
   int                   ntotalintvars;      /**< total number of Integer variables over all runs */
   int                   nruns;              /**< counts the number of runs, i.e. the number of initialized
                                              *   branch and bound processes */
   int                   binnblockvars;      /**< number of binary variables that appear in blocks  */
   int                   binnblocks;         /**< number of blocks with at least two variables */
       
   int                   intnblockvars;      /**< number of Integer variables that appear in blocks  */
   int                   intnblocks;         /**< number of blocks with at least two variables */
     
   int                   binnexchanges;      /**< number of executed changes of binary solution values leading to 
                                              *   improvement in objective function */
   int                   intnexchanges;      /**< number of executed changes of Integer solution values leading to improvement in
                                              *   objective function */
#endif
};

/*
 * Local methods
 */

/** switches the values of two binary or integer variables and checks feasibility with respect to the LP.
 *  feasibility check is applied at the same time, thus the first violation of feasibility will cause this method
 *  to stop constraint forward checking and to undo the applied changes to LP-rows before 
 *
 *  @todo adapt method not to copy entire activities array, but only the relevant region  
 */
static
SCIP_RETCODE shiftValues(
   SCIP*                 scip,               /**< scip instance */
   SCIP_VAR*             var1,               /**< first variable of variable pair*/
   SCIP_VAR*             var2,               /**< second variable of pair */
   SCIP_Real*            varvalue1,          /**< current value of variable1 in solution */
   SCIP_Real*            varvalue2,          /**< current value of variable2 in solution */
   SCIP_Real             shiftvalue1,        /**< the value that variable 1 should be shifted by */
   SCIP_Real             shiftvalue2,        /**< the value that variable 2 should be shifted by */
   SCIP_Real**           activities,         /**< the LP-row activities */
   int                   nrows,              /**< size of activities array */
   SCIP_Bool*            feasible            /**< set to true if method has successfully switched the variable values */
   )
{
   SCIP_COL* col;
   SCIP_ROW** rows1;
   SCIP_ROW** rows2;
   int ncolrows1;
   int ncolrows2;
   SCIP_Real* colvals1;
   SCIP_Real* colvals2;
   SCIP_Real* tmpactivities;
   int i;
   int j;
       
   assert(scip != NULL);
   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(varvalue1 != NULL);
   assert(varvalue2 != NULL);

   /* get variable specific rows and coefficients for both var1 and var2. */
   col = SCIPvarGetCol(var1);
   rows1 = SCIPcolGetRows(col);
   colvals1 = SCIPcolGetVals(col);
   ncolrows1 = SCIPcolGetNNonz(col);
   assert(ncolrows1 == 0 || rows1 != NULL);
   
   col = SCIPvarGetCol(var2);
   rows2 = SCIPcolGetRows(col);
   colvals2 = SCIPcolGetVals(col);
   ncolrows2 = SCIPcolGetNNonz(col);
   assert(ncolrows2 == 0 || rows2 != NULL);
   
   /* copy current activities to temporary array. if variable shifting turns out to be feasible, 
    * we can use the temporary array as new activities array, 
    * otherwise the infeasible temporary activities array can be freed.*/
   SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpactivities, *activities, nrows) );

   /* feasible can be set to true since there has not been any violated row, yet */
   *feasible =  TRUE;

   /* go through every LP-row of the variables' constraint set, update the activity and 
    * check if changed row activity will remain feasible.
    * We have to be careful about some cases concerning the two counters i and j */
   i = 0;
   j = 0;
   while( (i < ncolrows1 || j < ncolrows2) && (*feasible) )
   {
      int rowpos1;
      int rowpos2;
      SCIP_ROW* row;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int rowpos;

      /* If one counter has already reached its limit, assign a huge number to the corresponding 
       * row position to simulate an always greater row position. */
      if( i < ncolrows1 )
         rowpos1 = SCIProwGetLPPos(rows1[i]);
      else 
         rowpos1 = nrows;

      if( j < ncolrows2 )
         rowpos2 = SCIProwGetLPPos(rows2[j]);
      else
         rowpos2 = nrows;
     
      /* ensure that the rows are in current LP */
      if( rowpos1 == -1 || rowpos2 == -1 )
      {
         /* increase counters, depending on whether they point at a non-LP-Row at the moment */
         if( rowpos1 == -1 )
            ++i;
         if( rowpos2 == -1 )
            ++j; 

         continue;
      }

      assert(0 <= rowpos1 && 0 <= rowpos2);
      assert(rowpos1 < nrows || rowpos2 < nrows);

      /* the current row is the one with the smaller position */
      rowpos = MIN(rowpos1, rowpos2);

      if( rowpos1 <= rowpos2 )
         row = rows1[i];
      else 
         row = rows2[j];
      
      assert(row != NULL);

      /* get the necessary information about the current row, i.e its lower and upperbound and its activity */
      lhs = SCIProwGetLhs(row);
      rhs = SCIProwGetRhs(row);

      if( !SCIProwIsLocal(row) )
      {
         assert(SCIPisFeasGE(scip, tmpactivities[rowpos], lhs) && SCIPisFeasLE(scip, tmpactivities[rowpos], rhs));
      
         /* update activity with respect to which of the variables occur in this row */
	 if( rowpos1 <= rowpos2 )
            tmpactivities[rowpos] += colvals1[i] * shiftvalue1;
	 if( rowpos2 <= rowpos1 )
            tmpactivities[rowpos] += colvals2[j] * shiftvalue2;

         /* check if the LP_row-constraint was violated. In this case, stop forward checking loop by setting 
	  * feasibility flag */
         if( SCIPisFeasLT(scip, tmpactivities[rowpos], lhs) || SCIPisFeasGT(scip, tmpactivities[rowpos], rhs) )
         {
            *feasible = FALSE;
            break;
	 } 
         /* shifting doesn't violate this constraint. Continue with the next row */
         /* increase counters depending on whether the variables occured in current row */
      }
      if( rowpos1 <= rowpos2 )
         ++i;
      if( rowpos2 <= rowpos1 )
         ++j;
   }
   assert( (i == ncolrows1 && j == ncolrows2) || !*feasible );
   
   /* check if there is infeasible activity. If not, solution values and activities array can be updated. */ 
   if( *feasible )
   {
      *varvalue1 += shiftvalue1;
      *varvalue2 += shiftvalue2;
      SCIPfreeBufferArray(scip, activities);
      *activities = tmpactivities;
   }
   else
   {    
      /* free temporary array */
      SCIPfreeBufferArray(scip, &tmpactivities);
   }
   
   return SCIP_OKAY;
}

/** compare two variables with respect to their columns. Columns are treated as {0,1} vector, where every Nonzeroentry
 *  is treated as '1', and compared to each other lexicographically. I.e. var1 is < var2 if the corresponding column of
 *  var2 has the smaller single nonzero index of the two columns.  This comparison costs O(constraints) in the worst
 *  case
 */
static
int varColCompare(
   SCIP_VAR*                var1,                /**< left argument of comparison */
   SCIP_VAR*                var2                 /**< right argument of comparison */ 
   )
{
   SCIP_COL* col1;
   SCIP_COL* col2;
   SCIP_ROW** rows1;
   SCIP_ROW** rows2;
   int nnonzeros1;
   int nnonzeros2;
   int i;
   
   assert(var1 != NULL);
   assert(var2 != NULL);

   /* get the necessary row and column data */
   col1 = SCIPvarGetCol(var1);
   col2 = SCIPvarGetCol(var2);
   rows1 = SCIPcolGetRows(col1);
   rows2 = SCIPcolGetRows(col2);
   nnonzeros1 = SCIPcolGetNNonz(col1);
   nnonzeros2 = SCIPcolGetNNonz(col2);

   assert(nnonzeros1 == 0 || rows1 != NULL);
   assert(nnonzeros2 == 0 || rows2 != NULL);

   /* loop over the rows, stopped as soon as they differ in one index, 
    * or if counter reaches the end of a variables row set */
   for( i = 0; i < nnonzeros1 && i < nnonzeros2; ++i )
   {
      if( SCIProwGetIndex(rows1[i]) != SCIProwGetIndex(rows2[i]) )
         return SCIProwGetIndex(rows1[i]) - SCIProwGetIndex(rows2[i]);
   }

   /* loop is finished, without differing in one of common row indices, due to loop invariant 
    * variable i reached either nnonzeros1 or nnonzeros2 or both.
    * one can easily check that the difference of these two numbers always has the desired sign for comparison. */
   return nnonzeros2 - nnonzeros1 ;
}

/** implements a comparator to compare two variables with respect to their column entries */
static
SCIP_DECL_SORTPTRCOMP(SCIPvarcolComp)
{
   return varColCompare((SCIP_VAR*) elem1, (SCIP_VAR*) elem2);
}

/** checks if two given variables are contained in the same LP-constraints which makes them a a pair of candidates for
 *  the 2-opt heuristic. flag is set to true if and only if the variable 1 appears in a proper number of rows of
 *  variable 2 and vice versa, where 'proper' is is determined by the user parameter matchingrate which determines the
 *  necessary percentage of shared rows.
 */
static 
SCIP_Bool checkConstraintMatching(
   SCIP*                 scip,               /**< current SCIP instance */
   SCIP_VAR*             var1,               /**< first variable */
   SCIP_VAR*             var2,               /**< second variable */
   SCIP_Real             matchingrate        /**< determines the ratio of shared LP rows compared to the total number of 
                                              *   LP-rows each variable appears in */
   )
{
   SCIP_COL* col1;
   SCIP_COL* col2;
   SCIP_ROW** rows1;
   SCIP_ROW** rows2;
   int nnonzeros1;
   int nnonzeros2;
   int i;
   int j;
   int nrowsdiffer1;
   int nrowsdiffer2;
   int nrowmaximum;

   assert(var1 != NULL);
   assert(var2 != NULL);
      
   /* get the necessary row and column data */
   col1 = SCIPvarGetCol(var1);
   col2 = SCIPvarGetCol(var2);
   rows1 = SCIPcolGetRows(col1);
   rows2 = SCIPcolGetRows(col2);
   nnonzeros1 = SCIPcolGetNNonz(col1);
   nnonzeros2 = SCIPcolGetNNonz(col2);   

   assert(nnonzeros1 == 0 || rows1 != NULL);
   assert(nnonzeros2 == 0 || rows2 != NULL);
   
   /* initialize the counters for the number of rows not shared. 
    * For confusion, nrowsdiffer1 denotes the number of LP-rows of variable 1 which variable 2 doesn't appear in,
    * and vice versa. At the end of the row check,
    * the ratio of the counters in comparison with the individual total number of rows should not fall 
    * below the matchingrate. 
    */
   nrowmaximum = MAX(nnonzeros1, nnonzeros2);
   nrowsdiffer1 = nrowmaximum - nnonzeros2;
   nrowsdiffer2 = nrowmaximum - nnonzeros1;

   /* apply first matching check to the variables. if the numbers of nonzero rows differs too much, 
    * i.e with respect to matching ratio, the more expensive check over the rows doesn't have to be 
    * applied anymore because the counters for not shared rows can only increase.
    * To avoid divisions by zero, one has to exclude the case that the numbers of nonzeros are zero.
    * In this special case, i.e if one of the variables appears in zero columns,
    * this is defined as complete matching with the other variable. 
    * In the case that only one of the two variables appears in zero rows, 
    * this will cause the method to return false unless the matching rate is 0.0 as well 
    */
   if( SCIPisFeasGT(scip, matchingrate, nnonzeros1 == 0 ? 1.0 : (nnonzeros1 - nrowsdiffer1) / (SCIP_Real)nnonzeros1 ) ||
      SCIPisFeasGT(scip, matchingrate, nnonzeros2 == 0 ? 1.0 : (nnonzeros2 - nrowsdiffer2) / (SCIP_Real)nnonzeros2) )
   {    
      /* the desired matching rate cannot be satisfied by these two variables */  
      return FALSE;
   } 
   else if( nnonzeros1 == 0 && nnonzeros2 == 0)
      return TRUE;
  
   /* initialize counters */
   i = 0;
   j = 0;
   
   /* loop over all rows and determine in how many rows variable1 and variable2 match and in how many they differ, that is,
    * in how many constraints of variable1 variable 2 does not appear, denoted by nrowsdiffer1, and vice versa.
    * Check is complete if the end of one of the two arrays is reached since the remaining difference is equivalent to
    * the already registered difference of number of nonzeroes.
    */
   while( i < nnonzeros1 && j < nnonzeros2 )
   {
      if( SCIProwGetIndex(rows1[i]) == SCIProwGetIndex(rows2[j]) )
      { /* variables share a common row, continue by setting both indices */
         ++i;
         ++j;
      }
      else if( SCIProwGetIndex(rows1[i]) < SCIProwGetIndex(rows2[j]) )
      {  /* first row is a constraint in which variable2 doesn't appear, thus update nrowsdiffer1 and increase only counter i */
         ++i;
         ++nrowsdiffer1;
      } 
      else 
      {
         /* case SCIProwGetIndex(rows2[j] < SCIProwGetIndex(rows1[i])), that is variable2 has a constraint variable1 
          * doesn't appear in*/
         ++j;
         ++nrowsdiffer2;
      }
   }
  
   /* now apply the ratio based comparison, that is if the ratio of shared rows is greater equals the matching rate 
    * for each variable 
    */
   if( SCIPisFeasLE(scip, matchingrate, nnonzeros1 == 0 ? 1.0 : (SCIP_Real)(nnonzeros1 - nrowsdiffer1) / (SCIP_Real)(nnonzeros1)) ||
      SCIPisFeasLE(scip, matchingrate, nnonzeros2 == 0 ? 1.0 : (SCIP_Real)(nnonzeros2 - nrowsdiffer2) / (SCIP_Real)(nnonzeros2)) )
      return TRUE;
   else   
      return FALSE;
}



/** determines a bound by which the absolute solution value of two integer variables can be shifted at most.
 *  the criterion is the maintenance of feasibility of any global LP row. 
 *  first implementation only considers shifting proportion 1:1, i.e. if master value is shifted by a certain 
 *  Integer value k upwards, the value of slave is simultaneously shifted by k downwards. The restriction to this
 *  special case of the
 *  more general idea of shifting integer variables by an arbitrary proportion different from 1:1 is supposed
 *  to simplify the complexity of calculating such a proportion, which would be quadratic in the number of
 *  constraints.
 *
 *  @todo consider different shifting proportions than only 1:1 
 */
static
SCIP_Real determineBound(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_SOL*             sol,                /**< current incumbent */
   SCIP_VAR*             master,             /**< current master variable */
   int                   masterdirection,    /**< shifting direction of master variable */
   SCIP_VAR*             slave,              /**< slave variable with same LP_row set as master variable */
   int                   slavedirection,     /**< shifting direction of slave variable */
   SCIP_Real*            activities,         /**< array of LP row activities */
   int                   nrows               /**< the number of rows in LP and the size of the activities array */
   )
{
   SCIP_VAR* workmaster;
   SCIP_VAR* workslave;
   SCIP_Real masterbound;
   SCIP_Real slavebound;
   SCIP_Real bound;

   SCIP_COL* col;
   SCIP_ROW** rows1;
   SCIP_ROW** rows2;
   SCIP_Real* mastercolvals;
   SCIP_Real* slavecolvals;
   int ncolrows1;
   int ncolrows2;
   int i;
   int j;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(masterdirection == -slavedirection);
  
   assert(master != NULL);
   assert(slave != NULL);
   assert(SCIPvarIsIntegral(master) && SCIPvarIsIntegral(slave));

   /* in order to avoid ugly case differentiations, this method considers the master variable always 
    * to be shifted upwards and the slave variable to be shifted downwards, which may result in an exchange
    * of the original master and slave variable */
   if( masterdirection == 1 )
   {
      workmaster = master;
      workslave = slave;
   }
   else 
   {
      workmaster = slave;
      workslave = master;
   }
      
   /*
    * determine the trivial bound for shifting and ensure that neither of the two variables 
    * is already set to its most favorable gobal bound
    */
   masterbound = SCIPvarGetUbGlobal(workmaster) - SCIPgetSolVal(scip, sol, workmaster);
   slavebound = SCIPgetSolVal(scip, sol, workslave) - SCIPvarGetLbGlobal(workslave);
   bound = MIN(masterbound, slavebound);
   
   /* ensure that shifting is possible */
   if( SCIPisFeasEQ(scip, bound, 0.0) )
      return 0.0;
   
   /* get the necessary row and and column data for each variable */
   col = SCIPvarGetCol(workmaster);
   rows1 = SCIPcolGetRows(col);
   mastercolvals = SCIPcolGetVals(col);
   ncolrows1 = SCIPcolGetNNonz(col);

   col = SCIPvarGetCol(workslave);
   slavecolvals = SCIPcolGetVals(col);
   rows2 = SCIPcolGetRows(col);
   ncolrows2 = SCIPcolGetNNonz(col);

   assert(ncolrows1 == 0 || mastercolvals != NULL);
   assert(ncolrows2 == 0 || slavecolvals != NULL);
   
   /* loop over all LP rows and determine the maximum Integer bound by which both variables 
    * can be shifted without loss of feasibility 
    */
   i = 0;
   j = 0;
   while( (i < ncolrows1 || j < ncolrows2) && SCIPisFeasGT(scip, bound, 0.0) )
   {
      SCIP_ROW* row;
      SCIP_Real effect;
      SCIP_Real rhs;
      SCIP_Real lhs;
      SCIP_Real activity;
      int rowpos1;
      int rowpos2;
      int rowpos;

      assert(SCIPisFeasIntegral(scip, bound));
      /* several case distinctions have to be made. Firstly, assign the corresponding
       * row positions. if the end of one of the two arrays is reached, assign the 
       * corresponding row position with nrows which can be seen as the biggest number
       * occurring and will therefore not cause the corresponding variable to influence
       * the remaining rows which have to be checked
       */
      if( i < ncolrows1 )
         rowpos1 = SCIProwGetLPPos(rows1[i]);
      else
         rowpos1 = nrows;
      
      if( j < ncolrows2 )
         rowpos2 = SCIProwGetLPPos(rows2[j]);
      else 
         rowpos2 = nrows;

      /* ensure that both rows are in current LP. if at least one is not, increase the corresponding counter */
      if( rowpos1 == -1 || rowpos2 == -1 )
      {
         /* increase counters, if they point at a non-LP-row at the moment */
         if( rowpos1 == -1 )
            ++i;
         if( rowpos2 == -1 )
            ++j;

         continue;
      }
      assert(0 <= rowpos1 && 0 <= rowpos2);
      assert(rowpos1 < nrows || rowpos2 < nrows);

      /* the active row is always determined by the smaller of the two row positions */
      rowpos = MIN(rowpos1, rowpos2);
      if( rowpos1 <= rowpos2 )
         row = rows1[i];
      else 
         row = rows2[j];

      assert(row != NULL);

      /* local rows can be skipped */
      if( !SCIProwIsLocal(row) )
      {
         /* effect is the effect on the row activity by shifting the variables each by 1 , i.e master by +1 
          * and slave by -1 due to convenience. The effect has to be declared with respect to the occurrence
          * of each variable in the current row. If only one of the two appears in the current row, 
          * this will have an effect only in one direction. both variables affect the row activity
          * <=> rowpos1 = rowpos2 <=> (rowpos1 <= rowpos2 <= rowpos1)
          */
         effect = 0.0;
	 if( rowpos1 <= rowpos2 )
            effect += mastercolvals[i];
	 if( rowpos2 <= rowpos1 )
            effect -= slavecolvals[j];

         /* get information about the current row */
         activity = activities[rowpos];
         rhs = SCIProwGetRhs(row);
         lhs = SCIProwGetLhs(row);

         /* if the effect which shifting has on a row is zero, no constraint will be violated
          * by shifting the variables */
         if( SCIPisFeasEQ(scip, effect, 0.0) )
         {
            if( rowpos1 <= rowpos2 )
               ++i;
            if( rowpos2 <= rowpos1 )
               ++j;
	   
            continue;
         }

         /* effect does not equal zero, the bound is determined as minimum positive integer such that 
          * feasibility is remained in all constraints.
          * if constraint is an equality constraint, activity and lhs/rhs should be feasibly equal, which
          * will cause the method to return zero.*/

	 /* if the row has a lower bound, ensure that shifting preserves feasibility of this "<="-constraint */
	 if( !SCIPisInfinity(scip, -lhs) &&  SCIPisFeasLT(scip, activity + effect * bound, lhs) )
            bound = SCIPfeasFloor(scip, ABS((activity - lhs)/effect)); 
	 
	 /* if the row has an upper bound, ensure that shifting preserves feasibility of this ">="-constraint */
         if( !SCIPisInfinity(scip, rhs) && SCIPisFeasGT(scip, activity + effect * bound, rhs) )
            bound = SCIPfeasFloor(scip, ABS((rhs - activity)/effect));
      }
      
      /* increase the counters which belong to the corresponding row. Both counters are increased by 
       * 1 iff rowpos1 <= rowpos2 <= rowpos1 */
      if( rowpos1 <= rowpos2 )
         ++i;
      if( rowpos2 <= rowpos1 )
         ++j;
   }

   return bound;
}

/** initializes the required structures for execution of heuristic. 
 *  If objective coefficient functions are not all equal, each Binary and Integer variables are sorted 
 *  into heuristic-specific arrays with respect to their lexicographical column order, 
 *  where every zero in a column is interpreted as zero and every nonzero as '1'.
 *  After the sortation, the variables are compared with respect to user parameter matchingrate and 
 *  the heuristic specific blocks are determined.
 */
static 
SCIP_RETCODE presolve(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURDATA*        heurdata            /**< the heuristic data */
   )
{
   int nbinvars;
   int nintvars;
   int nvars;
   SCIP_VAR**  vars;
   SCIP_VAR** binvars;
   SCIP_VAR** intvars;

   assert(scip != NULL);
   assert(heurdata != NULL);

   /* ensure that method is not executed if presolving was already applied once in current branch and bound process */
   if( heurdata->presolved )
      return SCIP_OKAY;

   /* get necessary variable information, i.e. number of binary and integer variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   if( nbinvars > 1 )
   {  
      int c;
      int startindex;
      int nblocks;
      SCIP_Real lastcoeff;
      int nblockvars;
      
      SCIPdebugMessage("entering BINARY presolve \n");

      /* go through variables first and ensure that at least some variables differ in their objective coefficients,
       * presolve is stopped otherwise because no improvement can be exspected due to 1:1 exchanges */
      lastcoeff = SCIPvarGetObj(vars[0]);

      for( c = 1; c < nbinvars && SCIPisFeasEQ(scip,lastcoeff, SCIPvarGetObj(vars[c])); ++c )
      {
         lastcoeff = SCIPvarGetObj(vars[c]);
      } 
      if( c == nbinvars )
         heurdata->binequalcoeffs = TRUE;
      else 
         heurdata->binequalcoeffs = FALSE;

      /* check for equal coefficients first. heuristic does not need to be executed in that case */
      if( !heurdata->binequalcoeffs )
      {      
         /* allocate the heuristic specific variables array which is needed for heuristic specific comparison */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &heurdata->binvars, vars, nbinvars));

	 /* sort the variables with respect to their columns */
         SCIPsortPtr((void**)heurdata->binvars, SCIPvarcolComp, nbinvars);
         binvars = heurdata->binvars;
         heurdata->nbinvars = nbinvars;
         
         /* start determining blocks, i.e. a set of at least two variables which share most of their row set. 
          * If there is none, heuristic does not need to be executed.
          * startindex denotes the current startindex of a block. 
          * Once two differing variables are detected, the current counter position is compared to startindex. */
         startindex = 0;
         nblocks = 0;
         nblockvars = 0;
         heurdata->nobinarymatching = TRUE;

         /* loop over variables and compare neighbours*/
         for( c = 1; c < nbinvars; ++c ) 
         {
            if( !checkConstraintMatching(scip, binvars[startindex], binvars[c], heurdata->matchingrate) )
            {
               /* current block has its last variable at position c-1. If c differs from startindex by at least 2,
                * a block is detected. Update the data correspondingly */
               if( c - startindex >= 2 )
               {
                  heurdata->nobinarymatching = FALSE;
                  nblocks += 1;      
                  nblockvars += c -startindex;
               }
               startindex = c;
            } 
         }
#ifdef STATISTIC_INFORMATION         
         /* update statistics */
         heurdata->binnblocks += nblocks;
         heurdata->binnblockvars += nblockvars;
         heurdata->ntotalbinvars += nbinvars;

         SCIPinfoMessage(scip, NULL, "Twoopt binary presolving finished with <%d> blocks, <%d> block variables \n ", 
            nblocks, nblockvars);
      } 
      else 
      {
         SCIPinfoMessage(scip, NULL, "BINARY coefficients are all equal! \n");
#endif
      }
   } /* end of BINARY presolve */
  
   /* enters optional integer optimization presolve method */
   if( !heurdata->presolved && heurdata->intopt && nintvars > 1 )
   {               
      int c;
      SCIP_Real lastcoeff;
      
      SCIPdebugMessage("entering INTEGER presolve \n");

      /* go through variables first and ensure that at least some variables differ in their objective coefficients,
       * presolve is stopped otherwise because no improvement can be exspected due to 1:1 exchanges */
      vars = SCIPgetVars(scip);
      lastcoeff = SCIPvarGetObj(vars[0]);

      for( c = 1; c < nintvars && SCIPisFeasEQ(scip,lastcoeff, SCIPvarGetObj(vars[c])); ++c )
      {
         lastcoeff = SCIPvarGetObj(vars[c]);
      }
      if( c == nintvars )
         heurdata->intequalcoeffs = TRUE;
      else 
         heurdata->intequalcoeffs = FALSE;
      
      /* further presolve is only executed if Integer coefficients differ */
      if(!heurdata->intequalcoeffs)
      {  
         int nblockvars;
         int startindex;
         int nblocks;

         /* allocate the heuristic specific variables array which is needed for heuristic specific comparison */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &heurdata->intvars, &vars[nbinvars], nintvars) );
         heurdata->nintvars = nintvars;

	 /* sort variables with respect to their columns */
         SCIPsortPtr((void**)heurdata->intvars, SCIPvarcolComp, nintvars);

         intvars = heurdata->intvars;
         
         /* start determining blocks, i.e. a set of at least two variables which share most of their row set. 
          * If there is none, heuristic does not need to be executed. */
         startindex = 0;
         nblocks = 0;
         nblockvars = 0;
	 heurdata->nointmatching = TRUE;

         for( c = 1; c < nintvars ; ++c ) 
         {
            if( !checkConstraintMatching(scip, intvars[startindex], intvars[c], heurdata->matchingrate) )
            {
               /* current block has its last variable at position c-1. If c differs from startindex by at least 2, 
                * a block is detected. Update the data correspondingly */
               if( c - startindex >= 2 )
               { 
                  heurdata->nointmatching = FALSE; 
                  ++nblocks;      
                  nblockvars += c - startindex;
               }
               startindex = c;
            }
         }
#ifdef STATISTIC_INFORMATION 
         /* update statistics */
         heurdata->intnblocks += nblocks;
         heurdata->intnblockvars += nblockvars;
         heurdata->ntotalintvars += nintvars;

         SCIPinfoMessage(scip, NULL, "Twoopt Integer presolving finished with <%d> blocks, <%d> block variables \n ", 
            nblocks, nblockvars);

      } 
      else 
      {
         SCIPinfoMessage(scip, NULL, "INTEGER coefficients are all equal \n");
#endif
      }

      /* end of integer variables specific presolve */
      
   }

   /* presolving is finished, heuristic data is updated*/
   heurdata->presolved = TRUE;
   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeTwoopt)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitTwoopt)
{
   SCIP_HEURDATA* heurdata;
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* heuristic has not run yet, all heuristic specific data
    * is set to initial values */
   heurdata->nbinvars = 0;
   heurdata->nintvars = 0;
   heurdata->lastsolindex = -1;
   heurdata->nobinarymatching = TRUE;
   heurdata->nointmatching = TRUE;
   heurdata->presolved = FALSE;

#ifdef STATISTIC_INFORMATION
   /* init statistics */
   heurdata->binnexchanges = 0;
   heurdata->intnexchanges = 0;
   heurdata->binnblockvars = 0;
   heurdata->intnblockvars = 0;
   heurdata->binnblocks = 0;
   heurdata->intnblocks = 0;
   
   heurdata->ntotalbinvars = 0;
   heurdata->ntotalintvars = 0;
   heurdata->nruns = 0;
#endif
   
   /* all pointers are initially set to NULL. Since presolving
    * of the heuristic requires a lot of calculation time on some instances, 
    * but might not be needed e.g. if problem is infeasible, presolving is applied 
    * when heuristic is executed for the first time */
   heurdata->binvars = NULL;
   heurdata->intvars = NULL;
    
   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;   
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitTwoopt)
{
   SCIP_HEURDATA* heurdata;
  
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /*ensure that initialization was successful */
   assert(heurdata->nbinvars <= 1 || heurdata->binvars != NULL);

   /* free used memory for Binary Optimization*/
   if( heurdata->binvars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->binvars, heurdata->nbinvars);
   }

   /* free heuristic data of Integer variables */
   if( heurdata->intvars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->intvars, heurdata->nintvars);
   }

 
#ifdef STATISTIC_INFORMATION
   /* print relevant statistics to console */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, 
      "Twoopt Binary Statistics  :   "
      "%6.2g   %6.2g   %4.2g   %4.0g  (blocks/run, variables/run, varpercentage, average block size) \n",
      heurdata->nruns == 0 ? 0.0 : (SCIP_Real)heurdata->binnblocks/heurdata->nruns,
      heurdata->nruns == 0 ? 0.0 : (SCIP_Real)heurdata->binnblockvars/heurdata->nruns,
      heurdata->ntotalbinvars == 0 ? 0.0 : (SCIP_Real)heurdata->binnblockvars/heurdata->ntotalbinvars * 100.0,
      heurdata->binnblocks == 0 ? 0.0 : (SCIP_Real)heurdata->binnblockvars/(SCIP_Real)heurdata->binnblocks);
   
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, 
      "Twoopt Integer statistics :   " 
      "%6.2g   %6.2g   %4.2g   %4.0g  (blocks/run, variables/run, varpercentage, average block size) \n",
      heurdata->nruns == 0 ? 0.0 : (SCIP_Real)heurdata->intnblocks/heurdata->nruns,
      heurdata->nruns == 0 ? 0.0 : (SCIP_Real)heurdata->intnblockvars/heurdata->nruns,
      heurdata->ntotalintvars == 0 ? 0.0 : (SCIP_Real)heurdata->intnblockvars/heurdata->ntotalintvars * 100.0,
      heurdata->intnblocks == 0 ? 0.0 : (SCIP_Real)heurdata->intnblockvars/(SCIP_Real)heurdata->intnblocks);
   
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, 
      "Twoopt results            :   "
      "%6d   %6d   %4d   %4.2g  (runs, binary exchanges, Integer shiftings, matching rate)\n",
      heurdata->nruns,
      heurdata->binnexchanges,
      heurdata->intnexchanges,
      heurdata->matchingrate);

   /* set statistics to initial values*/
   heurdata->binnblockvars = 0;
   heurdata->binnblocks = 0;
   heurdata->intnblocks = 0;
   heurdata->intnblockvars = 0;
   heurdata->binnexchanges = 0;
   heurdata->intnexchanges = 0;
#endif

   heurdata->nbinvars = 0;
   heurdata->nintvars = 0;
  
   assert(heurdata->binvars == NULL);
   assert(heurdata->intvars == NULL);

   SCIPheurSetData(heur, heurdata);
   
   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolTwoopt)
{
   SCIP_HEURDATA* heurdata;
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);
   assert(heurdata->binvars == NULL && heurdata->intvars == NULL);

   /* set heuristic data to initial values, but increase the total number of runs */ 
   heurdata->nbinvars = 0;
   heurdata->nintvars = 0;
   heurdata->lastsolindex = -1;
   heurdata->nobinarymatching = TRUE;
   heurdata->nointmatching = TRUE;
   heurdata->presolved = FALSE;

#ifdef STATISTIC_INFORMATION
   ++(heurdata->nruns);
#endif

   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolTwoopt)
{
   SCIP_HEURDATA* heurdata;
   int nbinvars;
   int nintvars;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   nbinvars = heurdata->nbinvars;
   nintvars = heurdata->nintvars;

   /* free the allocated memory for the binary variables */
   if( heurdata->binvars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->binvars, nbinvars);
   }

   heurdata->nbinvars = 0;
   
   /* free the allocated memory for the integers */
   if( heurdata->intvars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->intvars, nintvars);
   }

   heurdata->nintvars = 0;

   assert(heurdata->binvars == NULL && heurdata->intvars == NULL);
   
   /* set heuristic data */
   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
   
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTwoopt)
{  /*lint --e{715}*/
   SCIP_HEURDATA*  heurdata;
   SCIP_SOL* bestsol;
   SCIP_SOL* worksol;
   SCIP_VAR** vars;
   SCIP_VAR** allvars;
   SCIP_ROW** lprows;
   SCIP_Real* activities;
 
   int nbinvars;
   int nintvars;
   int nvars;
   int nlprows;
   int i;
   int j;
   SCIP_Bool improvement;
  
   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   *result = SCIP_DIDNOTRUN;

   /* ensure that heuristic has not already been processed on current incumbent */
   bestsol = SCIPgetBestSol(scip);
   if( bestsol == NULL || heurdata->lastsolindex == SCIPsolGetIndex(bestsol) )
      return SCIP_OKAY;

   /* we can only work on solutions valid in the transformed space */
   if( SCIPsolGetOrigin(bestsol) == SCIP_SOLORIGIN_ORIGINAL )
      return SCIP_OKAY;
   /* ensure that the user defined number of nodes after last best solution has been reached, return otherwise */
   if( (SCIPgetNNodes(scip) - SCIPsolGetNodenum(bestsol)) < heurdata->waitingnodes )
      return SCIP_OKAY;
   /* we need an LP */
   if( SCIPgetNLPRows(scip) == 0 )
      return SCIP_OKAY;
   
   /* ensure that heuristic specific presolve is applied when heuristic is executed first */
   if( !heurdata->presolved )
   {
      SCIP_CALL( presolve(scip, heur, heurdata) );
   }

   assert(heurdata->presolved == TRUE);

   nbinvars = heurdata->nbinvars;
   nintvars = heurdata->nintvars;

   /* ensure that presolve has detected structures in the problem to which the 2-optimization can be applied.
    * That is if at least one of the sorted arrays in heuristic data, 'binvars' or 'intvars' exists and at least 
    * some of the variables' objective function coefficient differ. */
   if( (heurdata->binequalcoeffs || heurdata->nobinarymatching || nbinvars <= 1)
      && (!heurdata->intopt || heurdata->intequalcoeffs || heurdata->nointmatching  || nintvars <= 1) )
   {
      return SCIP_OKAY;
   }
   
   /* problem satisfies all necessary conditions for 2-optimization heuristic, execute heuristic! */
   *result = SCIP_DIDNOTFIND;

   /* initialize a working solution as a copy of the current incumbent to be able to store 
    * possible improvements obtained by 2-optimization */
   SCIP_CALL( SCIPcreateSolCopy(scip, &worksol, bestsol) );
   SCIPsolSetHeur(worksol, heur);
   
   /* get the LP row activities from current incumbent bestsol */
   SCIP_CALL( SCIPgetLPRowsData(scip, &lprows, &nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &activities, nlprows) );

   for( i = 0; i < nlprows; i++ )
   {
      SCIP_ROW* row;
      row = lprows[i];

      assert(row != NULL);
      assert(SCIProwGetLPPos(row) == i);
               
      activities[i] = SCIPgetRowSolActivity(scip, row, bestsol);
      assert(SCIProwIsLocal(row) || (SCIPisFeasGE(scip, activities[i], SCIProwGetLhs(row)) && SCIPisFeasLE(scip, activities[i], SCIProwGetRhs(row))));
   }

   heurdata->lastsolindex = SCIPsolGetIndex(bestsol);

   /* start with binary optimization */
   vars = heurdata->binvars;
   improvement = FALSE;

   if( vars != NULL && !heurdata->nobinarymatching )
   {
      int startindex;
      int endindex;
      
      startindex = 0;

      /* loop over variables by two counters, startindex and endindex, 
       * which are supposed to to mark the current start and end index
       * of a heuristic relevant block. The variable at the endindex is the current master and compared 
       * to all of its predecessors,
       * whether an exchange is possible and favorable with respect to objective function */
      while( startindex < nbinvars-1 )
      { 
         endindex = startindex +1;
         /* check if startindex and endindex still mark the same block */
         while( endindex < nbinvars && checkConstraintMatching(scip, vars[startindex], vars[endindex], heurdata->matchingrate) )  
         {
            /* determine the new master variable for heuristic's optimization method */
            SCIP_VAR* master;
            SCIP_Real masterobjective;
            SCIP_Real mastershiftvalue;
            SCIP_Real mastersolvalue;
            
            master = vars[endindex];
            masterobjective = SCIPvarGetObj(master);
            mastersolvalue = SCIPgetSolVal(scip, worksol, master);
            
            assert(SCIPisFeasEQ(scip, mastersolvalue, 1.0) || SCIPisFeasEQ(scip, mastersolvalue, 0.0));
            
            /* the shifting value of a binary variable is canonically calculated */
	    if( SCIPisFeasEQ(scip, mastersolvalue, 1.0) )
               mastershiftvalue = -1.0;
	    else
               mastershiftvalue = 1.0;

            /* loop over whole block and apply the necessary checks */
            for( j = startindex; j < endindex; ++j )
            {
               SCIP_VAR* slave;
               SCIP_Real slaveobjective;
               SCIP_Real slaveshiftvalue;
               SCIP_Real slavesolvalue;
               SCIP_Real changedobjective;

               /* get the next slave variable */
               slave = vars[j];
               slaveobjective = SCIPvarGetObj(slave);
               slavesolvalue = SCIPgetSolVal(scip, worksol, slave);
            
               assert(SCIPisFeasEQ(scip, slavesolvalue, 1.0) || SCIPisFeasEQ(scip, slavesolvalue, 0.0));

	       if( SCIPisFeasEQ(scip, slavesolvalue, 1.0) )
                  slaveshiftvalue = -1.0;
	       else 
                  slaveshiftvalue = 1.0;

               changedobjective = mastershiftvalue * masterobjective + slaveshiftvalue* slaveobjective;

               /* ensure that the given pair of variables can be flipped in favor of objective function,
                * and try to flip them */
               if( !SCIPisFeasEQ(scip, mastersolvalue, slavesolvalue) && SCIPisFeasLT(scip, changedobjective, 0.0) )
               {
                  SCIP_Bool feasible;
                  feasible = FALSE;

                  shiftValues(scip, master, slave, &mastersolvalue, &slavesolvalue, mastershiftvalue, 
                     slaveshiftvalue, &activities, nlprows, &feasible);
		 
		  /* if exchange was feasible, set new solution values */
                  if( feasible )
                  {
                     SCIP_CALL( SCIPsetSolVal(scip, worksol, master, mastersolvalue) );
                     SCIP_CALL( SCIPsetSolVal(scip, worksol, slave, slavesolvalue) );

                     SCIPdebugMessage("Exchange is feasible and executed. Changed Objective: <%.4g>\n", changedobjective);
                     SCIPdebugMessage("-> Variables: <%d><%g><%g> &&  <%d><%g><%g> (<index><value><objective coefficient>)\n",
                        endindex, mastersolvalue, masterobjective, j, slavesolvalue, slaveobjective);


#ifdef STATISTIC_INFORMATION
                     /* update statistics */
                     ++(heurdata->binnexchanges);
#endif
                     improvement = TRUE;
                  }
               }
            }
            ++endindex;
         }

         /* variables differ too much. go for the next block by setting startindex */
         startindex = endindex;
      }
   }
   
   /* ensure that their are at least two integer variables which do not have the same coefficient 
    * in the objective function. In one of these cases, the heuristic will automatically skip the
    * integer variable optimization */
   if( heurdata->intvars != NULL && !heurdata->nointmatching )
   {         
      int startindex;
      int endindex;
      vars = heurdata->intvars;
      
      startindex = 0;

      /* loop over variables of a heuristic relevant block. 
       * The variable at the endindex is the current master and compared to all of its predecessors,
       * whether an exchange is possible and favorable with respect to objective function*/
      while( startindex < nintvars - 1 ) 
      { 
         endindex = startindex +1;

         /* check if startindex and endindex still mark the same block*/
         while( endindex < nintvars && checkConstraintMatching(scip, vars[startindex], vars[endindex], heurdata->matchingrate) )
         {
            SCIP_VAR* master;
            SCIP_Real masterobjective;
            SCIP_Real mastersolvalue;
            int masterdirection;
           
            /* determine a new master variable for heuristic's optimization method. Heuristic 
             * will try to improve the objective function by shifting this master's solution value with 
             * another variable from the same block if it is favorable for the objective function  */
            master = vars[endindex];
            masterobjective = SCIPvarGetObj(master);
            mastersolvalue = SCIPgetSolVal(scip, worksol, master);

            /* loop over all possible exchange candidates and apply the necessary checks*/
            for( j = startindex; j < endindex; ++j )
            {
               SCIP_VAR* slave;
               SCIP_Real slaveobjective;
               SCIP_Real slavesolvalue;
               SCIP_Real changedobjective;
               SCIP_Real bound;
               int slavedirection;

               /* get next slave variable for exchange */
               slave = vars[j];
               slaveobjective = SCIPvarGetObj(slave);
               slavesolvalue = SCIPgetSolVal(scip, worksol, slave);

               /* the shifting directions have to be determined, choose the direction in favor of objective function.
                * if objective coefficients are feasibly equal, shifting will not have any impact on objective funtion,
                * thus this pair can be skipped. */
	       if( !SCIPisFeasEQ(scip, masterobjective - slaveobjective, 0.0) )
	       {
                  if( SCIPisFeasLT(scip, masterobjective - slaveobjective, 0.0) )
                     masterdirection = 1; 
                  else 
                     masterdirection = -1;
	       }
	       else
                  continue;

               slavedirection = -masterdirection;

	       /* the bound is determined, the smallest positive Integer to preserve feasibility of all rows 
                * if variables' solution values are shifted */
               bound = determineBound(scip, worksol, master, masterdirection, slave, slavedirection, activities, nlprows);
               changedobjective = masterobjective * masterdirection * bound + slaveobjective * slavedirection * bound;
   
	       /* if shifting has a negative impact on objective function (which is quite favorable with respect to minimization),
                * try to shift the variables, update the activitities and and store the new solution values*/
               if( SCIPisFeasLT(scip, changedobjective, 0.0) )
               {
                  SCIP_Bool feasible;
                  feasible = FALSE;
                  shiftValues(scip, master, slave, &mastersolvalue, &slavesolvalue, bound * masterdirection, 
                     bound * slavedirection, &activities, nlprows, &feasible);
                  if( feasible )
                  {
                     SCIP_CALL( SCIPsetSolVal(scip, worksol, master, mastersolvalue) );
                     SCIP_CALL( SCIPsetSolVal(scip, worksol, slave, slavesolvalue) );
                     SCIPdebugMessage("Exchange is feasible and executed. Changed Objective: <%.4g>\n", changedobjective);
                     SCIPdebugMessage("-> Variables: <%d><%g><%g> &&  <%d><%g><%g> (<index><value><objective coefficient>)\n",
                        endindex, mastersolvalue, masterobjective, j, slavesolvalue, slaveobjective);

#ifdef STATISTIC_INFORMATION
                     /* update statistics */
                     ++(heurdata->intnexchanges);
#endif
                     improvement = TRUE;
                  }
               }
            }
            ++endindex;
         }

         /* variables differ too much. go for the next block by setting startindex */
         startindex = endindex;
      }
   }

   nintvars = nbinvars + nintvars;
   nvars = SCIPgetNVars(scip);

   if( improvement && (nvars == nintvars) )
   {
      /* the problem is a pure IP, hence, no continuous or implicit variables are left for diving.
       * try if new working solution is feasible in original problem */
      SCIP_Bool success;

      SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, &success) );
      if( success )
      {
         SCIPdebugMessage("found feasible shifted solution:\n");
         SCIPdebug(SCIPprintSol(scip, worksol, NULL, FALSE));
         heurdata->lastsolindex = SCIPsolGetIndex(bestsol);
         *result = SCIP_FOUNDSOL;

#ifdef STATISTIC_INFORMATION
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "***Twoopt improved solution found by %10s . \n", 
            SCIPsolGetHeur(bestsol) != NULL ? SCIPheurGetName(SCIPsolGetHeur(bestsol)) :"Tree");
               
#endif
      }
   }
   else if( improvement )
   {
      /* fix the integer variables and start diving to optimize continuous variables with respect to reduced domain */
      SCIP_Bool lperror;
#ifdef NDEBUG
      SCIP_RETCODE retstat;
#endif

      SCIPdebugMessage("shifted solution should be feasible -> solve LP to fix continuous variables to best values\n");
      
      /* start diving to calculate the LP relaxation */
      SCIP_CALL( SCIPstartDive(scip) );

      allvars = SCIPgetVars(scip);
      
      /* set the bounds of the variables: fixed for integers, global bounds for continuous */
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPvarGetStatus(allvars[i]) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_CALL( SCIPchgVarLbDive(scip, allvars[i], SCIPvarGetLbGlobal(allvars[i])) );
            SCIP_CALL( SCIPchgVarUbDive(scip, allvars[i], SCIPvarGetUbGlobal(allvars[i])) );
         }
      }
      /* apply this after global bounds to not cause an error with intermediate empty domains */
      for( i = 0; i < nintvars; ++i )
      {
         if( SCIPvarGetStatus(allvars[i]) == SCIP_VARSTATUS_COLUMN )
         {  
            SCIP_Real solval;

            solval = SCIPgetSolVal(scip, worksol, allvars[i]);
            SCIP_CALL( SCIPchgVarLbDive(scip, allvars[i], solval) );
            SCIP_CALL( SCIPchgVarUbDive(scip, allvars[i], solval) );
         }
      }

      /* solve LP */
      SCIPdebugMessage(" -> old LP iterations: %"SCIP_LONGINT_FORMAT"\n", SCIPgetNLPIterations(scip));

      /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
       * Hence in optimized mode, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
       */
#ifdef NDEBUG
      retstat = SCIPsolveDiveLP(scip, -1, &lperror);
      if( retstat != SCIP_OKAY )
      { 
         SCIPwarningMessage("Error while solving LP in Twoopt heuristic; LP solve terminated with code <%d>\n",retstat);
      }
#else
      SCIP_CALL( SCIPsolveDiveLP(scip, -1, &lperror) );
#endif
         
      SCIPdebugMessage(" -> new LP iterations: %"SCIP_LONGINT_FORMAT"\n", SCIPgetNLPIterations(scip));
      SCIPdebugMessage(" -> error=%u, status=%d\n", lperror, SCIPgetLPSolstat(scip));

      /* check if this is a feasible solution */
      if( !lperror && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         SCIP_Bool success;
         
         /* copy the current LP solution to the working solution */
         SCIP_CALL( SCIPlinkLPSol(scip, worksol) );

         /* check solution for feasibility */
         SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, &success) );

         if( success )
         {
            SCIPdebugMessage("found feasible shifted solution:\n");
            SCIPdebug(SCIPprintSol(scip, worksol, NULL, FALSE));
            heurdata->lastsolindex = SCIPsolGetIndex(bestsol);
            *result = SCIP_FOUNDSOL;

#ifdef STATISTIC_INFORMATION
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "***Twoopt improved solution found by %10s . \n", 
               SCIPsolGetHeur(bestsol) != NULL ? SCIPheurGetName(SCIPsolGetHeur(bestsol)) :"Tree");
#endif
         }
      }

      /* terminate the diving */
      SCIP_CALL( SCIPendDive(scip) );
   }

   SCIPfreeBufferArray(scip, &activities);
   SCIP_CALL( SCIPfreeSol(scip, &worksol) );

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the twoopt primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurTwoopt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create twoopt primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
  
   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING,
         heurFreeTwoopt, heurInitTwoopt, heurExitTwoopt, 
         heurInitsolTwoopt, heurExitsolTwoopt, heurExecTwoopt,
         heurdata) );

   /* include boolean flag intopt */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/twoopt/intopt", " Should Integer-2-Optimization be applied or not?", 
         &heurdata->intopt, TRUE, DEFAULT_INTOPT, NULL, NULL) );

   /* include parameter waitingnodes */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/twoopt/waitingnodes", "user parameter to determine number of "
         "nodes to wait after last best solution before calling heuristic",
         &heurdata->waitingnodes, TRUE, DEFAULT_WAITINGNODES, 0, 10000, NULL, NULL));

   /* include parameter matchingrate */
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/twoopt/matchingrate", 
         "parameter to determine the percentage of rows two variables have to share before they are considered equal",
         &heurdata->matchingrate, TRUE, DEFAULT_MATCHINGRATE, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
