/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_twoopt.c,v 1.20 2011/01/02 11:10:46 bzfheinz Exp $"

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
#define HEUR_DISPCHAR         'B'
#define HEUR_PRIORITY         -20100
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1

#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE 
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

/* default parameter values */
#define DEFAULT_INTOPT                FALSE /**< optional integer optimization is applied by default */
#define DEFAULT_WAITINGNODES              0 /**< default number of nodes to wait after current best solution before calling heuristic */
#define DEFAULT_MATCHINGRATE            0.5 /**< default percentage by which two variables have to match in their LP-row set to be 
                                             *   associated as pair by heuristic */
#define DEFAULT_MAXNSLAVES              199 /**< default number of slave candidates for a master variable */

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
   SCIP_Bool             presolved;          /**< flag to indicate whether presolving has already been executed */
   int*                  binblockstart;      /**< array to store the start indices of each binary block */
   int*                  binblockend;        /**< array to store the end indices of each binary block */
   int                   nbinblocks;         /**< number of blocks */                          
   
   /* integer variable twoopt data */
   SCIP_Bool             intopt;             /**< parameter to determine if integer 2-opt should be applied */             
   SCIP_VAR**            intvars;            /**< array to store the integer variables in non-decreasing order 
                                              *   with respect to their objective coefficient */
   int                   nintvars;           /**< the number of integer variables stored in array intvars */
   int*                  intblockstart;      /**< array to store the start indices of each binary block */
   int*                  intblockend;        /**< array to store the end indices of each binary block */
   int                   nintblocks;         /**< number of blocks */                          

   SCIP_Bool             execute;            /**< has presolveTwoOpt detected 
					      * necessary structure for execution of heuristic? */
   unsigned int          randseed;           /**< seed value for random number generator */
   int                   maxnslaves;         /**< delimits the maximum number of slave candidates for a master variable */
#ifdef STATISTIC_INFORMATION
   /* statistics */
   int                   ntotalbinvars;      /**< total number of binary variables over all runs */
   int                   ntotalintvars;      /**< total number of Integer variables over all runs */
   int                   nruns;              /**< counts the number of runs, i.e. the number of initialized
                                              *   branch and bound processes */
   int                   maxbinblocksize;    /**< maximum size of a binary block */
   int                   maxintblocksize;    /**< maximum size of an integer block */
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

enum Opttype
   {
      OPTTYPE_BINARY = 1,
      OPTTYPE_INTEGER = 2
   }; 
typedef enum Opttype OPTTYPE;

/*
 * Local methods
 */

/** tries to switch the values of two binary or integer variables and checks feasibility with respect to the LP.
 *  @todo adapt method not to copy entire activities array, but only the relevant region  
 */
static
SCIP_RETCODE shiftValues(
   SCIP*                 scip,               /**< scip instance */
   SCIP_VAR*             var1,               /**< first variable of variable pair */
   SCIP_VAR*             var2,               /**< second variable of pair */
   SCIP_Real*            varvalue1,          /**< current value of variable1 in solution */
   SCIP_Real*            varvalue2,          /**< current value of variable2 in solution */
   SCIP_Real             shiftval,           /**< the value that variable 1 should be shifted by */
   SCIP_Real**           activities,         /**< the LP-row activities */
   int                   nrows,              /**< size of activities array */
   SCIP_Bool*            feasible            /**< set to true if method has successfully switched the variable values */
   )
{  /*lint --e{715}*/
   SCIP_COL* col;
   SCIP_ROW** rows1;
   SCIP_ROW** rows2;
   SCIP_Real* colvals1;
   SCIP_Real* colvals2;
   /* SCIP_Real* tmpactivities; */
   int ncolrows1;
   int ncolrows2;
   int i;
   int j;
 
   assert(scip != NULL);
   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(varvalue1 != NULL);
   assert(varvalue2 != NULL);

   assert(SCIPisFeasGE(scip, (*varvalue1)-shiftval, SCIPvarGetLbGlobal(var1)));
   assert(SCIPisFeasLE(scip, (*varvalue1)-shiftval, SCIPvarGetUbGlobal(var1)));
   
   assert(SCIPisFeasGE(scip, (*varvalue2)+shiftval, SCIPvarGetLbGlobal(var2)));
   assert(SCIPisFeasLE(scip, (*varvalue2)+shiftval, SCIPvarGetUbGlobal(var2)));

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

   for( i = 0; i < ncolrows1 && SCIProwGetLPPos(rows1[i]) >= 0; ++i ) 
   {
      int rowpos;
      rowpos = SCIProwGetLPPos(rows1[i]);

      if( rowpos >= 0 )
         (*activities)[rowpos] -= colvals1[i] * shiftval;
   }
           
   for( j = 0; j < ncolrows2 && SCIProwGetLPPos(rows2[j]) >= 0; ++j ) 
   {
      int rowpos;
      rowpos = SCIProwGetLPPos(rows2[j]);

      if( rowpos >= 0 )
      {
         (*activities)[rowpos] += colvals2[j] * shiftval;
         assert(SCIPisFeasGE(scip, (*activities)[rowpos], SCIProwGetLhs(rows2[j])));
         assert(SCIPisFeasLE(scip, (*activities)[rowpos], SCIProwGetRhs(rows2[j])));
      }
   }

   *varvalue1 -= shiftval;
   *varvalue2 += shiftval;

#ifndef NDEBUG
   for( i = 0; i < ncolrows1 && SCIProwGetLPPos(rows1[i]) >= 0; ++i )
   {
      assert(SCIPisFeasGE(scip, (*activities)[SCIProwGetLPPos(rows1[i])], SCIProwGetLhs(rows1[i])));
      assert(SCIPisFeasLE(scip, (*activities)[SCIProwGetLPPos(rows1[i])], SCIProwGetRhs(rows1[i])));
   }
#endif
   *feasible = TRUE;

   return SCIP_OKAY;
   
#if 0
   /* copy current activities to temporary array. If variable shifting turns out to be feasible, 
    * we can use the temporary array as new activities array, 
    * otherwise the infeasible temporary activities array can be freed.*/
   
    SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpactivities, *activities, nrows) );

    *feasible =  TRUE;
    i = 0;
    j = 0;
  
    /* go through every LP-row of the variables' constraint set, update the activity and 
     * check if changed row activity will remain feasible. */
    while( (i < ncolrows1 || j < ncolrows2) && (*feasible) )
    {
       SCIP_ROW* row;
       int rowpos;
       int index1;
       int index2;

       /* ensure that the rows are in current LP */
       /* increase counters, depending on whether they point at a non-LP-Row at the moment */
       if( i < ncolrows1 && SCIProwGetLPPos(rows1[i]) == -1 )
       {
          i = ncolrows1;
          continue;
       }

       if( j < ncolrows2 && SCIProwGetLPPos(rows2[j]) == -1 )
       {
          j = ncolrows2; 
          continue;
       }

       /* If one counter has already reached its limit, assign a huge number to the corresponding 
        * row position to simulate an always greater row position. */
       if( i < ncolrows1 )
          index1 = SCIProwGetIndex(rows1[i]);
       else 
          index1 = INT_MAX;

       if( j < ncolrows2 )
          index2 = SCIProwGetIndex(rows2[j]);
       else
          index2 = INT_MAX;
     
      
       assert(0 <= index1 && 0 <= index2);
       assert(index1 < INT_MAX || index2 < INT_MAX);

       /* the current row is the one with the smaller position */
       if( index1 <= index2 )
       {
          rowpos = SCIProwGetLPPos(rows1[i]);
          row = rows1[i];
       } 
       else
       {
          assert(j < ncolrows2);

          rowpos = SCIProwGetLPPos(rows2[j]);
          row = rows2[j];
       }      
       assert(row != NULL);

       if( !SCIProwIsLocal(row) )
       {
          SCIP_Real lhs;
          SCIP_Real rhs;

          /* get the necessary information about the current row, i.e its lower and upperbound and its activity */
          lhs = SCIProwGetLhs(row);
          rhs = SCIProwGetRhs(row);

          assert(SCIPisFeasGE(scip, tmpactivities[rowpos], lhs) && SCIPisFeasLE(scip, tmpactivities[rowpos], rhs));
      
          /* update activity with respect to which of the variables occur in this row */
          if( index1 <= index2 )
             tmpactivities[rowpos] -= colvals1[i] * shiftval;
          if( index2 <= index1 )
             tmpactivities[rowpos] += colvals2[j] * shiftval;

          /* check if the LP_row-constraint was violated. In this case, stop forward checking loop by setting 
           * feasibility flag */
          if( SCIPisFeasLT(scip, tmpactivities[rowpos], lhs) || SCIPisFeasGT(scip, tmpactivities[rowpos], rhs) )
          {
             *feasible = FALSE;
             break;
          } 
       }

       /* shifting doesn't violate this constraint. Continue with the next row 
        * increase counters depending on whether the variables occured in current row 
        */
       if( index1 <= index2 )
          ++i;
       if( index2 <= index1 )
          ++j;
    }
    assert( (i == ncolrows1 && j == ncolrows2) || !(*feasible) );
   
    /* check if there is infeasible activity. If not, solution values and activities array can be updated. */ 
    if( *feasible )
    {
       
       *varvalue1 -= shiftval;
       *varvalue2 += shiftval;
       SCIPfreeBufferArray(scip, activities);
       *activities = tmpactivities;

       assert(SCIPisFeasIntegral(scip, (*varvalue1)));
       assert(SCIPisFeasIntegral(scip, (*varvalue2)));
    }
    else
    {    
       /* free temporary array */
       SCIPfreeBufferArray(scip, &tmpactivities);
    }
   
    return SCIP_OKAY;
#endif
}

/** compare two variables with respect to their columns. Columns are treated as {0,1} vector, where every nonzero entry
 *  is treated as '1', and compared to each other lexicographically. I.e. var1 is < var2 if the corresponding column of
 *  var2 has the smaller single nonzero index of the two columns.  This comparison costs O(constraints) in the worst
 *  case
 */
static
int varColCompare(
   SCIP_VAR*             var1,               /**< left argument of comparison */
   SCIP_VAR*             var2                /**< right argument of comparison */ 
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

/** checks if two given variables are contained in common LP rows,
 *  returns true if variables share the necessary percentage (matchingrate) of rows.
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
   int nrows1not2;                           /* the number of LP-rows of variable 1 which variable 2 doesn't appear in */
   int nrows2not1;                           /* vice versa */
   int nrowmaximum;
   int nrowabs;

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
   
   /* initialize the counters for the number of rows not shared. */
   nrowmaximum = MAX(nnonzeros1, nnonzeros2);

   if( nrowmaximum == 0 )
      return TRUE;

   nrowabs = ABS(nnonzeros1-nnonzeros2);
   nrows1not2 = nrowmaximum - nnonzeros2;
   nrows2not1 = nrowmaximum - nnonzeros1;

   /* if the numbers of nonzero rows differs too much, w.r.t.matching ratio, the more expensive check over the rows 
    * doesn't have to be applied anymore because the counters for not shared rows can only increase.
    */
   assert(nrowmaximum > 0);

   if( (nrowmaximum - nrowabs) / (SCIP_Real) nrowmaximum < matchingrate )
      return FALSE;

   i = 0;
   j = 0;
   
   /* loop over all rows and determine number of non-shared rows */
   while( i < nnonzeros1 && j < nnonzeros2 )
   {
      /* variables share a common row */
      if( SCIProwGetIndex(rows1[i]) == SCIProwGetIndex(rows2[j]) )
      {
         ++i;
         ++j;
      }
      /* variable 1 appears in rows1[i], variable 2 doesn't */
      else if( SCIProwGetIndex(rows1[i]) < SCIProwGetIndex(rows2[j]) )
      { 
         ++i;
         ++nrows1not2;
      } 
      /* variable 2 appears in rows2[j], variable 1 doesn't */
      else 
      {
         ++j;
         ++nrows2not1;
      }
   }
  
   /* now apply the ratio based comparison, that is if the ratio of shared rows is greater equals the matching rate 
    * for each variable 
    */
   return ( SCIPisFeasLE(scip, matchingrate, (nnonzeros1 - nrows1not2) / (SCIP_Real)(nnonzeros1)) ||
      SCIPisFeasLE(scip, matchingrate, (nnonzeros2 - nrows2not1) / (SCIP_Real)(nnonzeros2)) );  /*lint !e795 */
}

/** determines a bound by which the absolute solution value of two integer variables can be shifted at most.
 *  the criterion is the maintenance of feasibility of any global LP row. 
 *  first implementation only considers shifting proportion 1:1, i.e. if master value is shifted by a certain 
 *  integer value k downwards, the value of slave is simultaneously shifted by k upwards. 
 *
 *  @todo consider different shifting proportions than 1:1 
 */
static
SCIP_Real determineBound(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_SOL*             sol,                /**< current incumbent */
   SCIP_VAR*             master,             /**< current master variable */
   SCIP_VAR*             slave,              /**< slave variable with same LP_row set as master variable */
   SCIP_Real*            activities,         /**< array of LP row activities */
   int                   nrows               /**< the number of rows in LP and the size of the activities array */
   )
{  /*lint --e{715}*/
   SCIP_Real masterbound;
   SCIP_Real slavebound;
   SCIP_Real bound;

   SCIP_COL* col;
   SCIP_ROW** slaverows;
   SCIP_ROW** masterrows;
   SCIP_Real* mastercolvals;
   SCIP_Real* slavecolvals;
   int nslaverows;
   int nmasterrows;
   int i;
   int j;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(master != NULL);
   assert(slave != NULL);
   assert(SCIPvarIsIntegral(master) && SCIPvarIsIntegral(slave));

   /*
    * determine the trivial bound for shifting and ensure that neither of the two variables 
    * is already set to its most favorable gobal bound
    */
   slavebound = SCIPvarGetUbGlobal(slave) - SCIPgetSolVal(scip, sol, slave);
   masterbound = SCIPgetSolVal(scip, sol, master) - SCIPvarGetLbGlobal(master);
   bound = MIN(slavebound, masterbound);
   bound = MAX(bound, 0.0);
   assert(!SCIPisInfinity(scip,bound));
   
   /* get the necessary row and and column data for each variable */
   col = SCIPvarGetCol(slave);
   slaverows = SCIPcolGetRows(col);
   slavecolvals = SCIPcolGetVals(col);
   nslaverows = SCIPcolGetNNonz(col);

   col = SCIPvarGetCol(master);
   mastercolvals = SCIPcolGetVals(col);
   masterrows = SCIPcolGetRows(col);
   nmasterrows = SCIPcolGetNNonz(col);

   assert(nslaverows == 0 || slavecolvals != NULL);
   assert(nmasterrows == 0 || mastercolvals != NULL);
   
   /* loop over all LP rows and determine the maximum Integer bound by which both variables 
    * can be shifted without loss of feasibility 
    */
   i = 0;
   j = 0;
   while( (i < nslaverows || j < nmasterrows) && !SCIPisFeasZero(scip, bound) )
   {
      SCIP_ROW* row;
      SCIP_Real effect;
      SCIP_Real rhs;
      SCIP_Real lhs;
      SCIP_Real activity;
      int rowpos;
      int index1;
      int index2;

      assert(SCIPisFeasIntegral(scip, bound));
      /* several case distinctions have to be made. Firstly, assign the corresponding
       * row positions. if the end of one of the two arrays is reached, assign the 
       * corresponding row position with nrows which can be seen as the biggest number
       * occurring and will therefore not cause the corresponding variable to influence
       * the remaining rows which have to be checked
       */

      if( i < nslaverows && SCIProwGetLPPos(slaverows[i]) == -1 )
      {
         i = nslaverows;
         continue;
      }

      if( j < nmasterrows && SCIProwGetLPPos(masterrows[j]) == -1 )
      {
         j = nmasterrows; 
         continue;
      }

      /* If one counter has already reached its limit, assign a huge number to the corresponding 
       * row index to simulate an always greater row position. */
      if( i < nslaverows )
	 index1 = SCIProwGetIndex(slaverows[i]);
      else 
         index1 = INT_MAX;

      if( j < nmasterrows )
	 index2 = SCIProwGetIndex(masterrows[j]);
      else
         index2 = INT_MAX;
     
      
      assert(0 <= index1 && 0 <= index2);
      assert(index1 < INT_MAX || index2 < INT_MAX);

      /* the current row is the one with the smaller position */
      if( index1 <= index2 )
      {
	 rowpos = SCIProwGetLPPos(slaverows[i]);
	 row = slaverows[i];
      } 
      else
      {
         assert(j < nmasterrows);

         rowpos = SCIProwGetLPPos(masterrows[j]);
         row = masterrows[j];
      }      
      assert(row != NULL);

      /* local rows can be skipped */
      if( !SCIProwIsLocal(row) )
      {
         /* effect is the effect on the row activity by shifting the variables each by 1 , i.e master by +1 
          * and slave by -1 due to convenience. The effect has to be declared with respect to the occurrence
          * of each variable in the current row. If only one of the two appears in the current row, 
          * this will have an effect only in one direction. both variables affect the row activity
          */
         effect = 0.0;
	 if( index1 <= index2 )
            effect += slavecolvals[i];
	 if( index2 <= index1 )
            effect -= mastercolvals[j];

         /* get information about the current row */
         activity = activities[rowpos];
         rhs = SCIProwGetRhs(row);
         lhs = SCIProwGetLhs(row);
         assert(SCIPisFeasLE(scip, lhs, activity) && SCIPisFeasLE(scip, activity, rhs));

         /* effect does not equal zero, the bound is determined as minimum positive integer such that 
          * feasibility is remained in all constraints.
          * if constraint is an equality constraint, activity and lhs/rhs should be feasibly equal, which
          * will cause the method to return zero.*/

	 /* if the row has a left hand side, ensure that shifting preserves feasibility of this ">="-constraint */
	 if( !SCIPisInfinity(scip, -lhs) && SCIPisFeasLT(scip, activity + effect * bound, lhs) )
         {
            assert(SCIPisNegative(scip, effect));
            assert(effect != 0.0);
            bound = SCIPfeasFloor(scip, (lhs - activity)/effect); /*lint !e795 */
         }
	 
	 /* if the row has an upper bound, ensure that shifting preserves feasibility of this "<="-constraint */
         if( !SCIPisInfinity(scip, rhs) && SCIPisFeasGT(scip, activity + effect * bound, rhs) )
         {
            assert(SCIPisPositive(scip, effect));
            assert(effect != 0.0);
            bound = SCIPfeasFloor(scip, (rhs - activity)/effect); /*lint !e795 */
         }
      }
      
      /* increase the counters which belong to the corresponding row. Both counters are increased by 
       * 1 iff rowpos1 <= rowpos2 <= rowpos1 */
      if( index1 <= index2 )
         ++i;
      if( index2 <= index1 )
         ++j;
   }

   return bound;
}

/** 
 * disposes variable with no heuristic relevancy, e.g., due to a fixed solution value, from its neighbourhood block. 
 * The sortation w.r.t. objective coefficients is maintained. The affected neighbourhood block is reduced by 1.
 */
static 
void disposeVariable(
   SCIP_VAR**            vars,               /**< problem variables */
   int*                  blockstart,         /**< contains start index of block */
   int*                  blockend,           /**< contains end index of block */
   int                   varindex,           /**< variable index */
   SCIP_Bool*            insertlast          /**< has the disposed variable been disposed as last element of block? */       
)
{
   int i;

   assert(blockstart != NULL);
   assert(blockend != NULL);
   assert(*blockstart <= varindex);
   assert(varindex <= *blockend);
   assert(insertlast != NULL);

   /* either the entire left or right block of varindex has to be shifted by one,
   * hence we check which array bound is closer to varindex
   * to minimize the number of exspected assignments */
   *insertlast = (*blockend) - varindex < ((*blockend) - (*blockstart)) / 2;

   /* variables are shifted */
   if( *insertlast )
     {
      for( i = varindex; i < (*blockend); ++i )
	 vars[i] = vars[i+1];
      --(*blockend);
   }
   else
   {
      for( i = varindex; i > (*blockstart); --i )
	 vars[i] = vars[i - 1];
      ++(*blockstart);
   }
}
    
/** realizes the presolve independently from type of variables it's applied to */
static 
SCIP_RETCODE innerPresolve(
   SCIP*                 scip,               /**< current scip */
   SCIP_VAR**            vars,               /**< problem vars */
   SCIP_VAR***           varspointer,        /**< pointer to heuristic specific variable memory */
   int                   nvars,              /**< the number of variables */
   int*                  nblocks,            /**< pointer to store the number of detected blocks */
   int*                  maxblocksize,       /**< maximum size of a block */
   int*                  nblockvars,         /**< pointer to store the number of block variables */
   int**                 blockstart,         /**< pointer to store the array of block start indices */
   int**                 blockend,           /**< pointer to store the array of block end indices */
   SCIP_Bool*            equalcoefficients,  /**< are all objective coefficients equal? */
   SCIP_HEUR*            heur,               /**< the heuristic */
   SCIP_HEURDATA*        heurdata            /**< the heuristic data */
   )
{
   SCIP_Real* objectives;      
   SCIP_Real lastcoeff;
   int v;
   int startindex;
  
   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars >= 2);
   assert(nblocks != NULL);
   assert(nblockvars != NULL);
   assert(blockstart != NULL);
   assert(blockend != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);

   /* go through variables first and ensure that not all objective coefficients are equal, 
    * presolve is stopped otherwise because no improvement can be exspected due to 1:1 exchanges */
   lastcoeff = SCIPvarGetObj(vars[0]);

   for( v = 1; v < nvars && SCIPisFeasEQ(scip, lastcoeff, SCIPvarGetObj(vars[v])); ++v )
      lastcoeff = SCIPvarGetObj(vars[v]);

   *equalcoefficients = (v == nvars);

   if( *equalcoefficients )
      return SCIP_OKAY;
    
   /* allocate the heuristic specific variables */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, varspointer, vars, nvars));
         
   /* sort the variables with respect to their columns */
   SCIPsortPtr((void**)(*varspointer), SCIPvarcolComp, nvars);
                  
   /* start determining blocks, i.e. a set of at least two variables which share most of their row set. 
    * If there is none, heuristic does not need to be executed.
    */
   startindex = 0;
   *nblocks = 0;
   *nblockvars = 0;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*blockstart)), nvars/2) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*blockend)), nvars/2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objectives, nvars) );

   objectives[startindex] = SCIPvarGetObj((*varspointer)[startindex]);
      
   /* loop over variables and compare neighbours*/
   for( v = 1; v < nvars; ++v ) 
   {
      objectives[v] = SCIPvarGetObj((*varspointer)[v]);

      if( !checkConstraintMatching(scip, (*varspointer)[startindex], (*varspointer)[v], heurdata->matchingrate) )
      {
         /* current block has its last variable at position v-1. If v differs from startindex by at least 2,
          * a block is detected. Update the data correspondingly */
         if( v - startindex >= 2 )
         {
            assert(*nblocks < nvars/2);
            (*nblockvars) += v - startindex;
            (*maxblocksize) = MAX((*maxblocksize), v - startindex);
            (*blockstart)[*nblocks] = startindex;
            (*blockend)[*nblocks] = v - 1;
            (*nblocks)++;      
            SCIPsortRealPtr(&(objectives[startindex]), (void**)&((*varspointer)[startindex]), v-startindex);
         }
         startindex = v;
      } 
   }

   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &((*blockstart)), nvars/2, *nblocks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &((*blockend)), nvars/2, *nblocks) );   
   SCIPfreeBufferArray(scip, &objectives);

   return SCIP_OKAY;
}

/** initializes the required structures for execution of heuristic. 
 *  If objective coefficient functions are not all equal, each Binary and Integer variables are sorted 
 *  into heuristic-specific arrays with respect to their lexicographical column order, 
 *  where every zero in a column is interpreted as zero and every nonzero as '1'.
 *  After the sortation, the variables are compared with respect to user parameter matchingrate and 
 *  the heuristic specific blocks are determined.
 */
static 
SCIP_RETCODE presolveTwoOpt(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURDATA*        heurdata            /**< the heuristic data */
   )
{
   int nbinvars;
   int nintvars;
   int nvars;
   SCIP_VAR**  vars;
   int nbinblockvars;
   int nintblockvars;
   int maxbinblocksize;
   int maxintblocksize;
   SCIP_Bool equalcoeffs;

   assert(scip != NULL);
   assert(heurdata != NULL);

   maxbinblocksize = 0;
   /* ensure that method is not executed if presolving was already applied once in current branch and bound process */
   if( heurdata->presolved )
      return SCIP_OKAY;

   /* get necessary variable information, i.e. number of binary and integer variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* if number of binary problem variables exceeds 2, they are subject to 2-optimization algorithm, hence heuristic
    * calls innerPresolve method to detect necessary structures.
    */
   if( nbinvars >= 2 ) 
   {
      
      SCIP_CALL( innerPresolve(scip, vars, &(heurdata->binvars), nbinvars, &(heurdata->nbinblocks), &maxbinblocksize, 
            &nbinblockvars, &(heurdata->binblockstart), &(heurdata->binblockend), &equalcoeffs, heur, heurdata) );
   }

   heurdata->nbinvars = nbinvars;

   heurdata->execute = nbinvars > 1 && !equalcoeffs && heurdata->nbinblocks > 0;  /*lint !e644 */   

#ifdef STATISTIC_INFORMATION        
   if( !equalcoeffs )
   {
      /* update statistics */
      heurdata->binnblocks += (heurdata->nbinblocks);
      heurdata->binnblockvars += nbinblockvars;
      heurdata->ntotalbinvars += nbinvars;
      heurdata->maxbinblocksize = MAX(maxbinblocksize, heurdata->maxbinblocksize);

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Twoopt BINARY presolving finished with <%d> blocks, <%d> block variables \n", 
         heurdata->nbinblocks, nbinblockvars);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "BINARY coefficients are all equal \n");
   }
#endif

   if( heurdata->intopt && nintvars > 1 )
   {               
      SCIP_CALL( innerPresolve(scip, &(vars[nbinvars]), &(heurdata->intvars), nintvars, &(heurdata->nintblocks), &maxintblocksize, 
            &nintblockvars, &(heurdata->intblockstart), &(heurdata->intblockend), &equalcoeffs, 
            heur, heurdata) );

      heurdata->execute = heurdata->execute || (!equalcoeffs && heurdata->nintblocks > 0);   
   
#ifdef STATISTIC_INFORMATION 
      if( !equalcoeffs )
      {      /* update statistics */
         heurdata->intnblocks += heurdata->nintblocks;
         heurdata->intnblockvars += nintblockvars;
         heurdata->ntotalintvars += nintvars;
         heurdata->maxintblocksize = MAX(maxintblocksize, heurdata->maxintblocksize);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Twoopt Integer presolving finished with <%d> blocks, <%d> block variables \n", 
            heurdata->nintblocks, nintblockvars);
      } 
      else 
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "INTEGER coefficients are all equal \n");
      }
#endif
   }
   heurdata->nintvars = nintvars;

   /* presolving is finished, heuristic data is updated*/
   heurdata->presolved = TRUE;
   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyTwoopt)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurTwoopt(scip) );

   return SCIP_OKAY;
}

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
   heurdata->presolved = FALSE;
   heurdata->nbinblocks = 0;
   heurdata->nintblocks = 0;

   heurdata->randseed = 0;

#ifdef STATISTIC_INFORMATION
   /* init statistics */
   heurdata->binnexchanges = 0;
   heurdata->intnexchanges = 0;
   heurdata->binnblockvars = 0;
   heurdata->intnblockvars = 0;
   heurdata->binnblocks = 0;
   heurdata->intnblocks = 0;

   heurdata->maxbinblocksize = 0;
   heurdata->maxintblocksize = 0;
   
   heurdata->ntotalbinvars = 0;
   heurdata->ntotalintvars = 0;
   heurdata->nruns = 0;
#endif
   
   /* all pointers are initially set to NULL. Since presolving
    * of the heuristic requires a lot of calculation time on some instances, 
    * but might not be needed e.g. if problem is infeasible, presolving is applied 
    * when heuristic is executed for the first time 
    */
   heurdata->binvars = NULL;
   heurdata->intvars = NULL;
   heurdata->binblockstart = NULL;
   heurdata->binblockend = NULL;
   heurdata->intblockstart = NULL;
   heurdata->intblockend = NULL;
    
   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;   
}

/* realizes the 2-optimization algorithm, which tries to improve incumbent solution
 * by shifting pairs of variables which share a common row set. 
 */
static
SCIP_RETCODE optimize(
   SCIP*                scip,               /**< current SCIP instance */
   SCIP_SOL*            worksol,            /**< working solution */
   SCIP_VAR**           vars,               /**< binary or integer variables */
   int                  nvars,              /**< number of type of variables */
   int*                 blockstart,         /**< contains start indices of blocks */
   int*                 blockend,           /**< contains end indices of blocks */
   int                  nblocks,            /**< the number of blocks */
   OPTTYPE              opttype,            /**< are binaries or integers optimized */
   SCIP_Real**          activities,         /**< the LP-row activities */
   int                  nrows,              /**< the number of LP rows */
   SCIP_Bool*           improvement,        /**< was there a successful shift? */
   SCIP_Bool*           varboundserr,       /**< has the current incumbent already been cut off */
   SCIP_HEURDATA*       heurdata            /**< the heuristic data */
   )
{  /*lint --e{715}*/
   int b;
   
   assert(scip != NULL);
   assert(nblocks > 0);
   assert(blockstart != NULL && blockend != NULL);
   assert(varboundserr != NULL);
   assert(activities != NULL);
   assert(worksol != NULL);
   assert(improvement != NULL);
   
   *varboundserr = FALSE;
   /* iterate over blocks */
   for( b = 0; b < nblocks; ++b )
   { 
      int e;
      
      /* iterate over variables in current block */
      for( e = blockend[b]; e > blockstart[b]; --e )
      {
         /* determine the new master variable for heuristic's optimization method */
         SCIP_VAR* master;

         SCIP_Real masterobj;
         SCIP_Real mastersolval;
         SCIP_Real bestshift;
         SCIP_Real bestbound;

         int bestslavepos;
         int s;
	 int firstslave;
	 int nslaves;
	 int blocklen;

         master = vars[e];
         masterobj = SCIPvarGetObj(master);
         mastersolval = SCIPgetSolVal(scip, worksol, master);

	 /* due to cuts or fixings of solution values, worksol might not be feasible w.r.t. its bounds.
	  * Exit method in that case. */
         if( SCIPisFeasGT(scip, mastersolval, SCIPvarGetUbGlobal(master)) || SCIPisFeasLT(scip, mastersolval, SCIPvarGetLbGlobal(master)) )
         {
            *varboundserr = TRUE;
            SCIPdebugMessage("Solution has violated variable bounds for var %s: %g <= %g <= %g \n",
                  SCIPvarGetName(master), SCIPvarGetLbGlobal(master), mastersolval, SCIPvarGetUbGlobal(master));
            return SCIP_OKAY;
         }    

	 /* if variable has fixed solution value, it is deleted from heuristic array */
         if( SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(master), SCIPvarGetLbGlobal(master) ) )
	 {
	    SCIP_Bool insertedlast;

	    disposeVariable(vars, &(blockstart[b]), &(blockend[b]), e, &insertedlast);
            continue;
	 }
         
         assert(SCIPisFeasIntegral(scip, mastersolval));

         assert(opttype == OPTTYPE_INTEGER || 
            (SCIPisFeasEQ(scip, mastersolval, 1.0) || SCIPisFeasEQ(scip, mastersolval, 0.0)));

         /* initialize the data of the best available shift */
         bestshift = 0.0;
         bestslavepos = -1;
         bestbound = 0.0;

	 blocklen = e - blockstart[b];
	 /* in blocks with more than heurdata->maxnslaves variables, a slave candidate
	  * region is chosen */
	 if( blocklen > heurdata->maxnslaves )
	 {
	    firstslave = SCIPgetRandomInt(blockstart[b], e - 1, &heurdata->randseed);
	    assert(firstslave < e);
	 }
	 else
	    firstslave = blockstart[b];

	 nslaves = MIN(heurdata->maxnslaves, blocklen);
         /* loop over block and determine a slave shift candidate for master variable.
          * If more than one candidate is available, choose the shift which improves objective function
          * the most.*/
	 for( s = 0; s < nslaves; ++s )
	 {
	    SCIP_VAR* slave;
	    SCIP_Real slaveobj;
	    SCIP_Real slavesolval;
	    SCIP_Real changedobj;
	    SCIP_Real bound;

	    int slaveindex;

	    slaveindex = (firstslave + s - blockstart[b]) % blocklen;
	    slaveindex += blockstart[b];

	    assert(slaveindex < e);
		  
	    /* get the next slave variable */
	    slave = vars[slaveindex];
	    slaveobj = SCIPvarGetObj(slave);
	    slavesolval = SCIPgetSolVal(scip, worksol, slave);

	    /* in blocks with a small number of variables, iteration over slaves is interrupted 
	     * in case of equal coefficients since no further improvement can be achieved.
	     * in case of a large block, the search for slave variables is continued at the
	     * first block variable. */
	    if( SCIPisFeasEQ(scip, slaveobj, masterobj) && blocklen <= heurdata->maxnslaves )
	       break;
	    else if( SCIPisFeasEQ(scip, slaveobj, masterobj) )
	    {
		s += e - slaveindex;
		continue;
	    }
	    assert(SCIPvarGetType(master) == SCIPvarGetType(slave));
	    assert(SCIPisFeasIntegral(scip, slavesolval));
	    assert(opttype == OPTTYPE_INTEGER || 
		   (SCIPisFeasEQ(scip, slavesolval, 1.0) || SCIPisFeasEQ(scip, slavesolval, 0.0)));
	    
	    /* solution is not feasible w.r.t. the variable bounds, stop optimization in this case */
	    if( SCIPisFeasGT(scip, slavesolval, SCIPvarGetUbGlobal(slave)) || SCIPisFeasLT(scip, slavesolval, SCIPvarGetLbGlobal(slave)) )
            {
	       *varboundserr = TRUE;
	       SCIPdebugMessage("Solution has violated variable bounds for var %s: %g <= %g <= %g \n",
				SCIPvarGetName(slave), SCIPvarGetLbGlobal(slave), slavesolval, SCIPvarGetUbGlobal(slave));
	       return SCIP_OKAY;
	    }    
	    
	    /* if solution value of the variable is fixed, delete it from the remaining candidates in the block */
	    if( SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(slave), SCIPvarGetLbGlobal(slave) ) )
	    {
	       SCIP_Bool insertedlast;
	       disposeVariable(vars, &(blockstart[b]), &(blockend[b]), slaveindex, &insertedlast);
	       
	       /* if variable has been deleted and variables were shifted from the end of the block, 
	        * the index of the master variable has been decreased by 1. In the other case, the position
	        * of the best slave has increased by 1.*/
	       if( insertedlast )
		  --e;
	       else  if( bestslavepos != -1 )
		  ++bestslavepos;

	       --blocklen;
	       continue;
	    }
	    /* determine the shifting direction to improve the objective function */
	    assert(SCIPisFeasGT(scip, masterobj, slaveobj));
	    
	    /* the maximum positive integer which preserves feasibility of all rows by shifting variables
	     * in the calculated directions.
	     */
	    bound = determineBound(scip, worksol, master, slave, *activities, nrows);
	    
	    assert(SCIPisFeasIntegral(scip, bound));
	    assert(SCIPisFeasGE(scip, bound, 0.0));
	    
	    /* the improvement of objective function is calculated */
	    changedobj = (slaveobj - masterobj) * bound;
	    
	    /* choose the candidate which improves the objective function the most */
	    if( SCIPisFeasLT(scip, changedobj, bestshift) )
	    {
	       bestshift = changedobj;
	       bestslavepos = slaveindex;
	       bestbound = bound;	    
	    }
	 }
       
	 /* choose the most promising candidate, if one exists */
         if( bestslavepos >= 0 )
         {
            SCIP_Real slavesolval;
            SCIP_VAR* slave;
            SCIP_Real slaveobj;
            SCIP_Bool feasible;
            SCIP_Real changedobj;

            slave = vars[bestslavepos];

            slavesolval = SCIPgetSolVal(scip, worksol, slave);
            slaveobj = SCIPvarGetObj(slave);

            /* the improvement of objective function is calculated */
            changedobj = (slaveobj - masterobj) * bestbound;
	  	  
            assert(SCIPisFeasLT(scip, changedobj, 0.0));

            /* try to change the solution values of the variables */
            feasible = FALSE;
            SCIP_CALL( shiftValues(scip, master, slave, &mastersolval, &slavesolval, bestbound,
                  activities, nrows, &feasible) );
	  
            if( feasible )
            {
               /* The variables' solution values were successfully shifted and can hence be updated. */
               assert(SCIPisFeasIntegral(scip, mastersolval));
               assert(SCIPisFeasIntegral(scip, slavesolval));
               assert(SCIPisFeasGE(scip, slavesolval, SCIPvarGetLbGlobal(slave)));
               assert(SCIPisFeasLE(scip, slavesolval, SCIPvarGetUbGlobal(slave)));
               assert(SCIPisFeasGE(scip, mastersolval, SCIPvarGetLbGlobal(master)));
               assert(SCIPisFeasLE(scip, mastersolval, SCIPvarGetUbGlobal(master)));

               SCIP_CALL( SCIPsetSolVal(scip, worksol, master, mastersolval) );
               SCIP_CALL( SCIPsetSolVal(scip, worksol, slave, slavesolval) );
               SCIPdebugMessage("Exchange is feasible and executed. Changed Objective: <%.4g>\n", changedobj);
               SCIPdebugMessage("-> Variables: <%d><%g><%g> &&  <%d><%g><%g> (<index><value><objective coefficient>)\n",
                  SCIPvarGetProbindex(master), mastersolval, masterobj, 
                  SCIPvarGetProbindex(slave), slavesolval, slaveobj);

#ifdef STATISTIC_INFORMATION
               /* update statistics */
               if( opttype == OPTTYPE_BINARY )
                  ++(heurdata->binnexchanges);
               else
                  ++(heurdata->intnexchanges);
#endif
               *improvement = TRUE;
            }
         }
      }
   }
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
   

#ifdef STATISTIC_INFORMATION
   /* print relevant statistics to console */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, 
      "Twoopt Binary Statistics  :   "
      "%6.2g   %6.2g   %4.2g   %4.0g %6d (blocks/run, variables/run, varpercentage, avg. block size, max block size) \n",
      heurdata->nruns == 0 ? 0.0 : (SCIP_Real)heurdata->binnblocks/(heurdata->nruns),
      heurdata->nruns == 0 ? 0.0 : (SCIP_Real)heurdata->binnblockvars/(heurdata->nruns),
      heurdata->ntotalbinvars == 0 ? 0.0 : (SCIP_Real)heurdata->binnblockvars/(heurdata->ntotalbinvars) * 100.0,
      heurdata->binnblocks == 0 ? 0.0 : heurdata->binnblockvars/(SCIP_Real)(heurdata->binnblocks),
      heurdata->maxbinblocksize);
   
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, 
      "Twoopt Integer statistics :   " 
      "%6.2g   %6.2g   %4.2g   %4.0g %6d (blocks/run, variables/run, varpercentage, avg block size, max block size) \n",
      heurdata->nruns == 0 ? 0.0 : (SCIP_Real)heurdata->intnblocks/(heurdata->nruns),
      heurdata->nruns == 0 ? 0.0 : (SCIP_Real)heurdata->intnblockvars/(heurdata->nruns),
      heurdata->ntotalintvars == 0 ? 0.0 : (SCIP_Real)heurdata->intnblockvars/(heurdata->ntotalintvars) * 100.0,
      heurdata->intnblocks == 0 ? 0.0 : heurdata->intnblockvars/(SCIP_Real)(heurdata->intnblocks),
      heurdata->maxintblocksize);
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
   /* free the allocated memory for the binary variables */
   if( heurdata->binvars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->binvars, heurdata->nbinvars);
   }
   if( heurdata->nbinblocks > 0 )
   {
      assert(heurdata->binblockstart != NULL);
      assert(heurdata->binblockend != NULL);

      SCIPfreeBlockMemoryArray(scip, &heurdata->binblockstart, heurdata->nbinblocks);
      SCIPfreeBlockMemoryArray(scip, &heurdata->binblockend, heurdata->nbinblocks);
   }
   heurdata->nbinvars = 0;
   heurdata->nbinblocks = 0;

   if( heurdata->nintblocks > 0 )
   {
      assert(heurdata->intblockstart != NULL);
      assert(heurdata->intblockend != NULL);

      SCIPfreeBlockMemoryArray(scip, &heurdata->intblockstart, heurdata->nintblocks);
      SCIPfreeBlockMemoryArray(scip, &heurdata->intblockend, heurdata->nintblocks);
   }
   
   /* free the allocated memory for the integers */
   if( heurdata->intvars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->intvars, heurdata->nintvars);
   }

   heurdata->nbinblocks = 0;
   heurdata->nintblocks = 0;
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
   assert(heurdata->binblockstart == NULL && heurdata->binblockend == NULL);
   assert(heurdata->intblockstart == NULL && heurdata->intblockend == NULL);

   /* set heuristic data to initial values, but increase the total number of runs */ 
   heurdata->nbinvars = 0;
   heurdata->nintvars = 0;
   heurdata->lastsolindex = -1;
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
   if( heurdata->binblockstart != NULL )
   {
      assert(heurdata->binblockend != NULL);

      SCIPfreeBlockMemoryArray(scip, &heurdata->binblockstart, heurdata->nbinblocks);
      SCIPfreeBlockMemoryArray(scip, &heurdata->binblockend, heurdata->nbinblocks);
   }
   heurdata->nbinvars = 0;
   heurdata->nbinblocks = 0;

   if( heurdata->intblockstart != NULL )
   {
      assert(heurdata->intblockend != NULL);

      SCIPfreeBlockMemoryArray(scip, &heurdata->intblockstart, heurdata->nintblocks);
      SCIPfreeBlockMemoryArray(scip, &heurdata->intblockend, heurdata->nintblocks);
   }
   heurdata->nintblocks = 0;
   
   /* free the allocated memory for the integers */
   if( heurdata->intvars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->intvars, nintvars);
   }

   heurdata->nintvars = 0;

   assert(heurdata->binvars == NULL && heurdata->intvars == NULL);
   assert(heurdata->binblockstart == NULL && heurdata->binblockend == NULL);
   assert(heurdata->intblockstart == NULL && heurdata->intblockend == NULL);

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
   SCIP_ROW** lprows;
   SCIP_Real* activities;
   SCIP_COL** cols;
   int ncols;
   int nbinvars;
   int nintvars;
   int ndiscvars;
   int nlprows;
   int i;
   SCIP_Bool improvement;
   SCIP_Bool presolthiscall;
   SCIP_Bool varboundserr;
  
   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   
   *result = SCIP_DIDNOTRUN;

   /* we need an LP */
   if( SCIPgetNLPRows(scip) == 0 )
      return SCIP_OKAY;

   bestsol = SCIPgetBestSol(scip);

   /* ensure that heuristic has not already been processed on current incumbent */
   if( bestsol == NULL || heurdata->lastsolindex == SCIPsolGetIndex(bestsol) )
      return SCIP_OKAY;

   heurdata->lastsolindex = SCIPsolGetIndex(bestsol);

   /* we can only work on solutions valid in the transformed space */
   if( SCIPsolGetOrigin(bestsol) == SCIP_SOLORIGIN_ORIGINAL )
      return SCIP_OKAY;

   /* ensure that the user defined number of nodes after last best solution has been reached, return otherwise */
   if( (SCIPgetNNodes(scip) - SCIPsolGetNodenum(bestsol)) < heurdata->waitingnodes )
      return SCIP_OKAY;

   presolthiscall = FALSE;
   
   /* ensure that heuristic specific presolve is applied when heuristic is executed first */
   if( !heurdata->presolved )
   {
      SCIP_CALL( SCIPgetLPColsData(scip,&cols, &ncols) );
      for( i = 0; i < SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip); ++i )
         SCIPcolSort(cols[i]);

      SCIP_CALL( presolveTwoOpt(scip, heur, heurdata) );
      presolthiscall = TRUE;
   }

   assert(heurdata->presolved);

   /* ensure that presolve has detected structures in the problem to which the 2-optimization can be applied.
    * That is if at least one of the sorted arrays in heuristic data, 'binvars' or 'intvars' exists and at least 
    * some of the variables' objective function coefficient differ. */
   if( !heurdata->execute )
      return SCIP_OKAY;

   nbinvars = heurdata->nbinvars;
   nintvars = heurdata->nintvars;
   ndiscvars = nbinvars+nintvars;

   /* we need to be able to start diving from current node in order to resolve the LP
    * with continuous or implicit integer variables
    */
   if( SCIPgetNVars(scip) > ndiscvars && ( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL ) )
      return SCIP_OKAY;

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

      /* Heuristic does not provide infeasibility recovery, thus if any constraint is violated,
       * execution has to be terminated.
       */
      if( !SCIProwIsLocal(row) && (SCIPisFeasLT(scip, activities[i], SCIProwGetLhs(row)) 
            || SCIPisFeasGT(scip, activities[i], SCIProwGetRhs(row))) )
         goto TERMINATE;
   }

   if( !presolthiscall )
   {
      SCIP_CALL( SCIPgetLPColsData(scip,&cols, &ncols) );
      for( i = 0; i < ndiscvars; ++i )
      {
         SCIPcolSort(cols[i]);
      }
   }

   /* start with binary optimization */
   improvement = FALSE;
   varboundserr = FALSE;

   if( heurdata->nbinblocks > 0 )
   {
      SCIP_CALL( optimize(scip, worksol, heurdata->binvars, heurdata->nbinvars, heurdata->binblockstart, heurdata->binblockend, heurdata->nbinblocks,
            OPTTYPE_BINARY, &activities, nlprows, &improvement, &varboundserr, heurdata) );
   }
   if( varboundserr ) 
      goto TERMINATE;

   /* ensure that their are at least two integer variables which do not have the same coefficient 
    * in the objective function. In one of these cases, the heuristic will automatically skip the
    * integer variable optimization */
   if( heurdata->nintblocks > 0 )
   {         
      assert(heurdata->intopt);
      SCIP_CALL( optimize(scip, worksol, heurdata->intvars, heurdata->nintvars, heurdata->intblockstart, heurdata->intblockend, heurdata->nintblocks,
            OPTTYPE_INTEGER, &activities, nlprows, &improvement, &varboundserr, heurdata) );
   }
 
   if( !improvement || varboundserr )
      goto TERMINATE;

   if( SCIPgetNVars(scip) == ndiscvars )
   {
      /* the problem is a pure IP, hence, no continuous or implicit variables are left for diving.
       * try if new working solution is feasible in original problem */
      SCIP_Bool success;

      SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, FALSE, &success) );
      if( success )
      {
         SCIPdebugMessage("found feasible shifted solution:\n");
         SCIPdebug( SCIP_CALL( SCIPprintSol(scip, worksol, NULL, FALSE) ) );
         heurdata->lastsolindex = SCIPsolGetIndex(bestsol);
         *result = SCIP_FOUNDSOL;

#ifdef STATISTIC_INFORMATION
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "***Twoopt improved solution found by %10s . \n", 
            SCIPsolGetHeur(bestsol) != NULL ? SCIPheurGetName(SCIPsolGetHeur(bestsol)) :"Tree");
               
#endif
      }
   }
   /* fix the integer variables and start diving to optimize continuous variables with respect to reduced domain */
   else
   {
      SCIP_VAR** allvars;
      SCIP_Bool lperror;
#ifdef NDEBUG
      SCIP_RETCODE retstat;
#endif

      SCIPdebugMessage("shifted solution should be feasible -> solve LP to fix continuous variables to best values\n");

      allvars = SCIPgetVars(scip);
            
      /* start diving to calculate the LP relaxation */
      SCIP_CALL( SCIPstartDive(scip) );

      /* set the bounds of the variables: fixed for integers, global bounds for continuous */
      for( i = 0; i < SCIPgetNVars(scip); ++i )
      {
         if( SCIPvarGetStatus(allvars[i]) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_CALL( SCIPchgVarLbDive(scip, allvars[i], SCIPvarGetLbGlobal(allvars[i])) );
            SCIP_CALL( SCIPchgVarUbDive(scip, allvars[i], SCIPvarGetUbGlobal(allvars[i])) );
         }
      }

      /* apply this after global bounds to not cause an error with intermediate empty domains */
      for( i = 0; i < ndiscvars; ++i )
      {
         if( SCIPvarGetStatus(allvars[i]) == SCIP_VARSTATUS_COLUMN )
         {  
            SCIP_Real solval;

            solval = SCIPgetSolVal(scip, worksol, allvars[i]);
	    assert(SCIPvarGetType(allvars[i]) == SCIP_VARTYPE_CONTINUOUS || SCIPisFeasIntegral(scip, solval));

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
         SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, FALSE, &success) );

         if( success )
         {
            SCIPdebugMessage("found feasible shifted solution:\n");
            SCIPdebug( SCIP_CALL( SCIPprintSol(scip, worksol, NULL, FALSE) ) );
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

 TERMINATE:
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
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyTwoopt,
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

   /* include parameter maxnslaves */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/twoopt/maxnslaves", "maximum number of slaves for one master variable",
			      &heurdata->maxnslaves, TRUE, DEFAULT_MAXNSLAVES, 1, 1000000, NULL, NULL) );

   /* include parameter matchingrate */
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/twoopt/matchingrate", 
         "parameter to determine the percentage of rows two variables have to share before they are considered equal",
         &heurdata->matchingrate, TRUE, DEFAULT_MATCHINGRATE, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
