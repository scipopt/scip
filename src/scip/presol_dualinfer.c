/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    presol_dualinfer.c
 * @ingroup PRESOLVERS
 * @brief   dual inference presolver
 * @author  Dieter Weninger
 *
 * This presolver exploits dual informations on continuous variables for
 * fixings of integer/continuous variables and side changes.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/pub_matrix.h"
#include "scip/cons_linear.h"
#include "presol_dualinfer.h"

#define PRESOL_NAME             "dualinfer"
#define PRESOL_DESC             "exploit dual informations for fixings and side changes"
#define PRESOL_PRIORITY              -200     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS                0     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */
#define MAX_LOOPS                       7     /**< maximal number of dual bound strengthening loops */


/*
 * Data structures
 */


/** type of fixing direction */
enum Fixingdirection
{
   FIXATLB = -1,
   NOFIX   =  0
};
typedef enum Fixingdirection FIXINGDIRECTION;

/** type of side change */
enum SideChange
{
   RHSTOLHS = -1,
   NOCHANGE = 0,
   LHSTORHS = 1
};
typedef enum SideChange SIDECHANGE;


/*
 * Local methods
 */

/** calculate minimal column activity from one continuous variable without one row */
static
SCIP_Real getMinColActWithoutRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   withoutrow,         /**< exclude row index */
   SCIP_Real*            lbdual[2],          /**< lower bounds of dual variables */
   SCIP_Real*            ubdual[2],          /**< upper bounds of dual variables */
   int                   part,               /**< which of part of the dual variables is active */
   SCIP_Real*            lbdualbounds[2],    /**< lower bounds of dual variables corresponding to primal variable bounds */
   SCIP_Real*            ubdualbounds[2]     /**< upper bounds of dual variables corresponding to primal variable bounds */
   )
{
   SCIP_VAR* var;
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real val;
   SCIP_Real mincolactivity;
   int row;
   int p;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);
   assert(lbdualbounds[0] != NULL);
   assert(ubdualbounds[0] != NULL);
   assert(lbdualbounds[1] != NULL);
   assert(ubdualbounds[1] != NULL);

   var = SCIPmatrixGetVar(matrix, col);

   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

   mincolactivity = 0;

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);

   for( ; colpnt < colend; colpnt++, valpnt++ )
   {
      row = *colpnt;
      val = *valpnt;

      if( row == withoutrow )
      {
         if( SCIPmatrixIsRowRhsInfinity(matrix, row) )
         {
            continue;
         }
         else
         {
            if( part == 0 )
            {
               /* consider second part */
               val = -val;
               p = 1;
            }
            else
            {
               /* consider first part */
               p = 0;
            }

            if( val > 0.0 )
            {
               assert(!SCIPisInfinity(scip, -lbdual[p][row]));
               mincolactivity += val * lbdual[p][row];
            }
            else if( val < 0.0 )
            {
               assert(!SCIPisInfinity(scip, ubdual[p][row]));
               mincolactivity += val * ubdual[p][row];
            }
         }
      }
      else
      {
         p = 0;

         do {

            if( val > 0.0 )
            {
               assert(!SCIPisInfinity(scip, -lbdual[p][row]));
               mincolactivity += val * lbdual[p][row];
            }
            else if( val < 0.0 )
            {
               assert(!SCIPisInfinity(scip, ubdual[p][row]));
               mincolactivity += val * ubdual[p][row];
            }

            val = -val;
            p++;

         } while ( !SCIPmatrixIsRowRhsInfinity(matrix, row) && p < 2 );
      }
   }

   /* consider variable bounds */
   if( !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
   {
      /* upper bound */
      assert(!SCIPisInfinity(scip, ubdualbounds[0][col]));
      mincolactivity += -ubdualbounds[0][col];
   }
   if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
   {
      /* lower bound */
      assert(!SCIPisInfinity(scip, -lbdualbounds[1][col]));
      mincolactivity += lbdualbounds[1][col];
   }

   return mincolactivity;
}

/** calculate minimal column activity from one continuous variable without one row for explicit bounds */
static
SCIP_Real getMinColActWithoutBound(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   SCIP_Real*            lbdual[2],          /**< lower bounds of dual variables */
   SCIP_Real*            ubdual[2],          /**< upper bounds of dual variables */
   SCIP_Real*            lbdualbounds[2],    /**< lower bounds of dual variables corresponding to primal variable bounds */
   SCIP_Real*            ubdualbounds[2],    /**< upper bounds of dual variables corresponding to primal variable bounds */
   SCIP_Bool             explicitub,         /**< is this the explicit upper bound constraint */
   SCIP_Bool             explicitlb          /**< is this the explicit lower bound constraint */
   )
{
   SCIP_VAR* var;
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real val;
   SCIP_Real mincolactivity;
   int row;
   int p;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);
   assert(lbdualbounds[0] != NULL);
   assert(ubdualbounds[0] != NULL);
   assert(lbdualbounds[1] != NULL);
   assert(ubdualbounds[1] != NULL);

   var = SCIPmatrixGetVar(matrix, col);

   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);
   assert((explicitub && !explicitlb) || (!explicitub && explicitlb));

   mincolactivity = 0;

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);

   for(; (colpnt < colend); colpnt++, valpnt++ )
   {
      row = *colpnt;
      val = *valpnt;

      p = 0;

      do {

         if( val > 0.0 )
         {
            assert(!SCIPisInfinity(scip, -lbdual[p][row]));
            mincolactivity += val * lbdual[p][row];
         }
         else if( val < 0.0 )
         {
            assert(!SCIPisInfinity(scip, ubdual[p][row]));
            mincolactivity += val * ubdual[p][row];
         }

         val = -val;
         p++;

      } while ( !SCIPmatrixIsRowRhsInfinity(matrix, row) && p < 2 );
   }

   /* consider variable bounds */
   if( !explicitub && !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
   {
      /* upper bound */
      assert(!SCIPisInfinity(scip, ubdualbounds[0][col]));
      mincolactivity += -ubdualbounds[0][col];
   }
   if( !explicitlb && !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
   {
      /* lower bound */
      assert(!SCIPisInfinity(scip, -lbdualbounds[1][col]));
      mincolactivity += lbdualbounds[1][col];
   }

   return mincolactivity;
}

/** calculate minimal/maximal column residual activities */
static
void calcColActResidualCommon(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< matrix coefficient */
   SCIP_Real*            lbdual[2],          /**< lower bounds of dual variables */
   SCIP_Real*            ubdual[2],          /**< upper bounds of dual variables */
   int                   part,               /**< which of part of the dual variables is active */
   SCIP_Real*            lbdualbounds[2],    /**< lower bounds of dual variables corresponding to primal variable bounds */
   SCIP_Real*            ubdualbounds[2],    /**< upper bounds of dual variables corresponding to primal variable bounds */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   SCIP_Real*            maxcolact,          /**< maximal column activities */
   int*                  maxcolactinf,       /**< number of (positive) infinite contributions to maximal column activity */
   int*                  mincolactinf,       /**< number of (negative) infinite contributions to minimal column activity */
   SCIP_Real*            mincolresact        /**< minimal column residual activity */
   )
{
   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);
   assert(lbdualbounds[0] != NULL);
   assert(ubdualbounds[0] != NULL);
   assert(lbdualbounds[1] != NULL);
   assert(ubdualbounds[1] != NULL);
   assert(mincolact != NULL);
   assert(maxcolact != NULL);
   assert(maxcolactinf != NULL);
   assert(mincolactinf != NULL);
   assert(mincolresact != NULL);

   assert(SCIPvarGetType(SCIPmatrixGetVar(matrix, col)) == SCIP_VARTYPE_CONTINUOUS);

   if( !SCIPmatrixIsRowRhsInfinity(matrix, row) && part == 1 )
      val = -val;

   if( val > 0.0 )
   {
      if( SCIPisInfinity(scip, -lbdual[part][row]) )
      {
         assert(mincolactinf[col] >= 1);
         if( mincolactinf[col] == 1 )
            *mincolresact = getMinColActWithoutRow(scip, matrix, col, row,
               lbdual, ubdual, part, lbdualbounds, ubdualbounds);
         else
            *mincolresact = -SCIPinfinity(scip);
      }
      else
      {
         if( mincolactinf[col] > 0 )
            *mincolresact = -SCIPinfinity(scip);
         else
            *mincolresact = mincolact[col] - val * lbdual[part][row];
      }
   }
   else if( val < 0.0 )
   {
      if( SCIPisInfinity(scip, ubdual[part][row]) )
      {
         assert(mincolactinf[col] >= 1);
         if( mincolactinf[col] == 1 )
            *mincolresact = getMinColActWithoutRow(scip, matrix, col, row,
               lbdual, ubdual, part, lbdualbounds, ubdualbounds);
         else
            *mincolresact = -SCIPinfinity(scip);
      }
      else
      {
         if( mincolactinf[col] > 0 )
            *mincolresact = -SCIPinfinity(scip);
         else
            *mincolresact = mincolact[col] - val * ubdual[part][row];
      }
   }
}

/** calculate minimal/maximal column residual activities of explicit bounds */
static
void calcColActResidualExplicitBound(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   SCIP_Real*            lbdual[2],          /**< lower bounds of dual variables */
   SCIP_Real*            ubdual[2],          /**< upper bounds of dual variables */
   SCIP_Real*            lbdualbounds[2],    /**< lower bounds of dual variables corresponding to primal variable bounds */
   SCIP_Real*            ubdualbounds[2],    /**< upper bounds of dual variables corresponding to primal variable bounds */
   SCIP_Bool             explicitub,         /**< is this the constraint for the explicit upper bound */
   SCIP_Bool             explicitlb,         /**< is this the constraint for the explicit lower bound */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   SCIP_Real*            maxcolact,          /**< maximal column activities */
   int*                  maxcolactinf,       /**< number of (positive) infinite contributions to maximal column activity */
   int*                  mincolactinf,       /**< number of (negative) infinite contributions to minimal column activity */
   SCIP_Real*            mincolresact        /**< minimal column residual activity */
   )
{
   SCIP_VAR*  var;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);
   assert(mincolact != NULL);
   assert(maxcolact != NULL);
   assert(maxcolactinf != NULL);
   assert(mincolactinf != NULL);
   assert(mincolresact != NULL);

   var = SCIPmatrixGetVar(matrix, col);

   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);
   assert((explicitub && !explicitlb) || (!explicitub && explicitlb));

   if( explicitub )
   {
      if( !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
      {
         /* -1.0 * x >= -ub */
         if( SCIPisInfinity(scip, ubdualbounds[0][col]) )
         {
            assert(mincolactinf[col] >= 1);

            if( mincolactinf[col] == 1 )
               *mincolresact = getMinColActWithoutBound(scip, matrix, col,
                  lbdual, ubdual, lbdualbounds, ubdualbounds, explicitub, explicitlb);
            else
               *mincolresact = -SCIPinfinity(scip);
         }
         else
         {
            if( mincolactinf[col] > 0 )
               *mincolresact = -SCIPinfinity(scip);
            else
               *mincolresact = mincolact[col] + ubdualbounds[0][col];
         }
      }
   }
   else
   {
      if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
      {
         /* 1.0 * y >= lb */
         if( SCIPisInfinity(scip, -lbdualbounds[1][col]) )
         {
            assert(mincolactinf[col] >= 1);
            if( mincolactinf[col] == 1 )
               *mincolresact = getMinColActWithoutBound(scip, matrix, col,
                  lbdual, ubdual, lbdualbounds, ubdualbounds, explicitub, explicitlb);
            else
               *mincolresact = -SCIPinfinity(scip);
         }
         else
         {
            if( mincolactinf[col] > 0 )
               *mincolresact = -SCIPinfinity(scip);
            else
               *mincolresact = mincolact[col] - lbdualbounds[1][col];
         }
      }
   }
}


/** calculate minimal/maximal column activity on continuous variables */
static
void calcColActivity(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix containing the constraints */
   int                   startcol,           /**< start column index */
   int                   stopcol,            /**< stop column index */
   SCIP_Real*            lbdual[2],          /**< lower bounds of dual variables */
   SCIP_Real*            ubdual[2],          /**< upper bounds of dual variables */
   SCIP_Real*            lbdualbounds[2],    /**< lower bounds of dual variables corresponding to primal variable bounds */
   SCIP_Real*            ubdualbounds[2],    /**< upper bounds of dual variables corresponding to primal variable bounds */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   SCIP_Real*            maxcolact,          /**< maximal column activities */
   int*                  maxcolactinf,       /**< number of (positive) infinite contributions to maximal column activity */
   int*                  mincolactinf        /**< number of (negative) infinite contributions to minimal column activity */
   )
{
   SCIP_VAR* var;
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real val;
   int row;
   int c;
   int part;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);
   assert(lbdualbounds[0] != NULL);
   assert(ubdualbounds[0] != NULL);
   assert(lbdualbounds[1] != NULL);
   assert(ubdualbounds[1] != NULL);
   assert(mincolact != NULL);
   assert(maxcolact != NULL);
   assert(maxcolactinf != NULL);
   assert(mincolactinf != NULL);

   for( c = startcol; c < stopcol; c++ )
   {
      mincolact[c] = 0;
      maxcolact[c] = 0;
      maxcolactinf[c] = 0;
      mincolactinf[c] = 0;

      colpnt = SCIPmatrixGetColIdxPtr(matrix, c);
      colend = colpnt + SCIPmatrixGetColNNonzs(matrix, c);
      valpnt = SCIPmatrixGetColValPtr(matrix, c);
      var = SCIPmatrixGetVar(matrix, c);

      /* calculate column activities */
      for( ; colpnt < colend; colpnt++, valpnt++ )
      {
         row = *colpnt;
         val = *valpnt;
         part = 0;

         do {

            if( val > 0 )
            {
               if(SCIPisInfinity(scip, ubdual[part][row]))
                  maxcolactinf[c]++;
               else
                  maxcolact[c] += val * ubdual[part][row];

               if(SCIPisInfinity(scip, -lbdual[part][row]))
                  mincolactinf[c]++;
               else
                  mincolact[c] += val * lbdual[part][row];
            }
            else if( val < 0.0 )
            {
               if(SCIPisInfinity(scip, -lbdual[part][row]))
                  maxcolactinf[c]++;
               else
                  maxcolact[c] += val * lbdual[part][row];

               if(SCIPisInfinity(scip, ubdual[part][row]))
                  mincolactinf[c]++;
               else
                  mincolact[c] += val * ubdual[part][row];
            }

            val = -val;
            part++;

         } while( !SCIPmatrixIsRowRhsInfinity(matrix, row) && part < 2 );
      }

      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         /* consider variable bounds */
         if( !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
         {
            /* upper bound */
            if(SCIPisInfinity(scip, -lbdualbounds[0][c]))
               maxcolactinf[c]++;
            else
               maxcolact[c] += -lbdualbounds[0][c];

            if(SCIPisInfinity(scip, ubdualbounds[0][c]))
               mincolactinf[c]++;
            else
               mincolact[c] += -ubdualbounds[0][c];
         }

         if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
         {
            /* lower bound */

            if(SCIPisInfinity(scip, ubdualbounds[1][c]))
               maxcolactinf[c]++;
            else
               maxcolact[c] += ubdualbounds[1][c];

            if(SCIPisInfinity(scip, -lbdualbounds[1][c]))
               mincolactinf[c]++;
            else
               mincolact[c] += lbdualbounds[1][c];
         }
      }

      /* update column activities */
      if( mincolactinf[c] > 0 )
         mincolact[c] = -SCIPinfinity(scip);
      if( maxcolactinf[c] > 0 )
         maxcolact[c] = SCIPinfinity(scip);
   }
}


/** update bounds on dual variables */
static
void updateDualBounds(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix containing the constraints */
   SCIP_Real             objval,             /**< objective function value */
   SCIP_Real             val,                /**< matrix coefficient */
   int                   row,                /**< row index */
   SCIP_Real             mincolresact,       /**< minimal column residual activity */
   SCIP_Real*            lbdual[2],          /**< dual lower bounds (ranged, equality) */
   SCIP_Real*            ubdual[2],          /**< dual upper bounds (ranged, equatity)*/
   int                   part,               /**< which of part of the dual variables is active */
   int*                  boundchanges,       /**< number of bound changes */
   SCIP_Bool*            updateinfcnt        /**< flag indicating to update infinity counters */
   )
{
   SCIP_Real newlbdual;
   SCIP_Real newubdual;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);
   assert(boundchanges != NULL);
   assert(updateinfcnt != NULL);

   *updateinfcnt = FALSE;

   if( !SCIPmatrixIsRowRhsInfinity(matrix, row) && part == 1 )
      val = -val;

   if( val > 0 )
   {
      if( !SCIPisInfinity(scip, mincolresact) && !SCIPisInfinity(scip, -mincolresact) )
      {
         newubdual = (objval - mincolresact) / val;
         if( newubdual < ubdual[part][row] )
         {
            if( SCIPisInfinity(scip, ubdual[part][row]) )
               *updateinfcnt = TRUE;

            ubdual[part][row] = newubdual;
            (*boundchanges)++;
         }
      }
   }
   else if( val < 0 )
   {
      if( !SCIPisInfinity(scip, mincolresact) && !SCIPisInfinity(scip, -mincolresact) )
      {
         newlbdual = (objval - mincolresact) / val;
         if( newlbdual > lbdual[part][row] )
         {
            if( SCIPisInfinity(scip, -lbdual[part][row]) )
               *updateinfcnt = TRUE;

            lbdual[part][row] = newlbdual;
            (*boundchanges)++;
         }
      }
   }
}

/** update bounds on dual variables of explicit bounds */
static
void updateDualBoundsExplicit(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix containing the constraints */
   SCIP_Real             objval,             /**< objective function value */
   int                   col,                /**< column index */
   SCIP_Real             mincolresact,       /**< minimal column residual activity */
   SCIP_Real*            lbdualbounds[2],    /**< lower bounds of dual variables corresponding to variable bounds */
   SCIP_Real*            ubdualbounds[2],    /**< upper bounds of dual variables corresponding to variable bounds */
   SCIP_Bool             explicitub,         /**< flag indicating if explicit upper bound is active */
   SCIP_Bool             explicitlb,         /**< flag indicating if explicit lower bound is active */
   int*                  boundchanges,       /**< number of bound changes */
   SCIP_Bool*            updateinfcnt        /**< flag indicating to update infinity counters */
   )
{
   SCIP_VAR* var;
   SCIP_Real newlbdual;
   SCIP_Real newubdual;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdualbounds[0] != NULL);
   assert(ubdualbounds[0] != NULL);
   assert(lbdualbounds[1] != NULL);
   assert(ubdualbounds[1] != NULL);
   assert(boundchanges != NULL);
   assert(updateinfcnt != NULL);

   assert((explicitub && !explicitlb) || (!explicitub && explicitlb));

   var = SCIPmatrixGetVar(matrix, col);
   *updateinfcnt = FALSE;

   if( explicitlb )
   {
      if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
      {
         if( !SCIPisInfinity(scip, mincolresact) && !SCIPisInfinity(scip, -mincolresact) )
         {
            newubdual = (objval - mincolresact);
            if( newubdual < ubdualbounds[1][col] )
            {
               if( SCIPisInfinity(scip, ubdualbounds[1][col]) )
                  *updateinfcnt = TRUE;

               ubdualbounds[1][col] = newubdual;
               (*boundchanges)++;
            }
         }
      }
   }
   else
   {
      if( !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
      {
         if( !SCIPisInfinity(scip, mincolresact) && !SCIPisInfinity(scip, -mincolresact) )
         {
            newlbdual = -(objval - mincolresact);
            if( newlbdual > lbdualbounds[0][col] )
            {
               if( SCIPisInfinity(scip, -lbdualbounds[0][col]) )
                  *updateinfcnt = TRUE;

               lbdualbounds[0][col] = newlbdual;
               (*boundchanges)++;
            }
         }
      }
   }
}

/** update minimal/maximal column activity infinity counters */
static
void infCntUpdate(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix containing the constraints */
   SCIP_Real             val,                /**< matrix coefficient */
   int                   row,                /**< row index */
   SCIP_Real*            lbdual[2],          /**< lower bounds of dual variables */
   SCIP_Real*            ubdual[2],          /**< upper bounds of dual variables */
   SCIP_Real*            lbdualbounds[2],    /**< lower bounds of dual variables corresponding to primal variable bounds */
   SCIP_Real*            ubdualbounds[2],    /**< upper bounds of dual variables corresponding to primal variable bounds */
   int                   part,               /**< which part of the dual variables is active */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   SCIP_Real*            maxcolact,          /**< maximal column activities */
   int*                  maxcolactinf,       /**< number of (positive) infinity contributions to maximal column activity */
   int*                  mincolactinf        /**< number of (negative) infinity contributions to minimal column activity */
   )
{
   SCIP_Real* valpnt;
   int* rowpnt;
   int* rowend;
   SCIP_Real aij;
   int c;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);
   assert(lbdualbounds[0] != NULL);
   assert(ubdualbounds[0] != NULL);
   assert(lbdualbounds[1] != NULL);
   assert(ubdualbounds[1] != NULL);
   assert(mincolact != NULL);
   assert(maxcolact != NULL);
   assert(maxcolactinf != NULL);
   assert(mincolactinf != NULL);

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   if( !SCIPmatrixIsRowRhsInfinity(matrix, row) && part == 1 )
      val = -val;

   /* look at all columns entries present within row and update the corresponding infinity counters.
    * if one counter gets to zero, then calculate this column activity completely new
    */

   if( val > 0 )
   {
      /* finite upper bound change */
      for(; (rowpnt < rowend); rowpnt++, valpnt++ )
      {
         c = *rowpnt;
         aij = *valpnt;

         if( aij > 0 )
         {
            assert(maxcolactinf[c] > 0);
            maxcolactinf[c]--;

            if( maxcolactinf[c] == 0 )
               calcColActivity(scip, matrix, c, c+1, lbdual, ubdual,
                  lbdualbounds, ubdualbounds, mincolact, maxcolact,
                  maxcolactinf, mincolactinf);
         }
         else if( aij < 0 )
         {
            assert(mincolactinf[c] > 0);
            mincolactinf[c]--;

            if( mincolactinf[c] == 0 )
               calcColActivity(scip, matrix, c, c+1, lbdual, ubdual,
                  lbdualbounds, ubdualbounds, mincolact, maxcolact,
                  maxcolactinf, mincolactinf);
         }
      }
   }
   else if( val < 0 )
   {
      /* finite lower bound change */
      for(; (rowpnt < rowend); rowpnt++, valpnt++ )
      {
         c = *rowpnt;
         aij = *valpnt;

         if( aij > 0 )
         {
            assert(mincolactinf[c] > 0);
            mincolactinf[c]--;

            if( mincolactinf[c] == 0 )
               calcColActivity(scip, matrix, c, c+1, lbdual, ubdual,
                  lbdualbounds, ubdualbounds, mincolact, maxcolact,
                  maxcolactinf, mincolactinf);
         }
         else if( aij < 0 )
         {
            assert(maxcolactinf[c] > 0);
            maxcolactinf[c]--;

            if( maxcolactinf[c] == 0 )
               calcColActivity(scip, matrix, c, c+1, lbdual, ubdual,
                  lbdualbounds, ubdualbounds, mincolact, maxcolact,
                  maxcolactinf, mincolactinf);
         }
      }
   }
}

/** update minimal/maximal column activity infinity counters for explicit bounds */
static
void infCntUpdateExplicit(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   SCIP_Real*            lbdual[2],          /**< lower bounds of dual variables */
   SCIP_Real*            ubdual[2],          /**< upper bounds of dual variables */
   SCIP_Real*            lbdualbounds[2],    /**< lower bounds of dual variables corresponding to variable bounds */
   SCIP_Real*            ubdualbounds[2],    /**< upper bounds of dual variables corresponding to variable bounds */
   SCIP_Bool             explicitub,         /**< flag indicating if explicit upper bound is active */
   SCIP_Bool             explicitlb,         /**< flag indicating if explicit lower bound is active */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   SCIP_Real*            maxcolact,          /**< maximal column activities */
   int*                  maxcolactinf,       /**< number of (positive) infinity contributions to maximal column activity */
   int*                  mincolactinf        /**< number of (negative) infinity contributions to minimal column activity */
   )
{
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);
   assert(lbdualbounds[0] != NULL);
   assert(ubdualbounds[0] != NULL);
   assert(lbdualbounds[1] != NULL);
   assert(ubdualbounds[1] != NULL);
   assert(mincolact != NULL);
   assert(maxcolact != NULL);
   assert(maxcolactinf != NULL);
   assert(mincolactinf != NULL);

   assert((explicitub && !explicitlb) || (!explicitub && explicitlb));

   var = SCIPmatrixGetVar(matrix, col);

   if( explicitlb )
   {
      if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
      {
         /* finite upper bound change */
         assert(maxcolactinf[col] > 0);
         maxcolactinf[col]--;

         if( maxcolactinf[col] == 0 )
            calcColActivity(scip, matrix, col, col+1, lbdual, ubdual,
               lbdualbounds, ubdualbounds, mincolact, maxcolact,
               maxcolactinf, mincolactinf);
      }
   }
   else
   {
      if( !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
      {
         /* finite lower bound change */
         assert(maxcolactinf[col] > 0);
         maxcolactinf[col]--;

         if( maxcolactinf[col] == 0 )
            calcColActivity(scip, matrix, col, col+1, lbdual, ubdual,
               lbdualbounds, ubdualbounds, mincolact, maxcolact,
               maxcolactinf, mincolactinf);
      }
   }
}


/** dual bound strengthening */
static
SCIP_RETCODE dualBoundStrengthening(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix containing the constraints */
   FIXINGDIRECTION*      varstofix,          /**< array holding information for later upper/lower bound fixing */
   int*                  npossiblefixings,   /**< number of possible fixings */
   SIDECHANGE*           sidestochange,      /**< array holding if this is an implied equality */
   int*                  npossiblesidechanges/**< number of possible equality changes */
   )
{
   SCIP_Real* valpnt;
   SCIP_Real* lbdual[2];
   SCIP_Real* ubdual[2];
   SCIP_Real* lbdualbounds[2];
   SCIP_Real* ubdualbounds[2];
   SCIP_Real* mincolact;
   SCIP_Real* maxcolact;
   SCIP_Bool* varissingcol;
   int* maxcolactinf;
   int* mincolactinf;
   int* colpnt;
   int* colend;
   int boundchanges;
   int loops;
   int c;
   int r;
   int nrows;
   int ncols;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(varstofix != NULL);
   assert(npossiblefixings != NULL);
   assert(sidestochange != NULL);
   assert(npossiblesidechanges != NULL);

   nrows = SCIPmatrixGetNRows(matrix);
   ncols = SCIPmatrixGetNColumns(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &varissingcol, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mincolact, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxcolact, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxcolactinf, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mincolactinf, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lbdual[0], nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubdual[0], nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lbdual[1], nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubdual[1], nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lbdualbounds[0], ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubdualbounds[0], ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lbdualbounds[1], ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubdualbounds[1], ncols) );

   /* initialize dual bounds */
   for( r = 0; r < nrows; r++ )
   {
      lbdual[0][r] = 0.0;
      ubdual[0][r] = SCIPinfinity(scip);
      lbdual[1][r] = 0.0;
      ubdual[1][r] = SCIPinfinity(scip);
   }

   /* intialize dual bounds of primal variable bounds */
   for( c = 0; c < ncols; c++ )
   {
      lbdualbounds[0][c] = 0.0;
      ubdualbounds[0][c] = SCIPinfinity(scip);
      lbdualbounds[1][c] = 0.0;
      ubdualbounds[1][c] = SCIPinfinity(scip);
   }

   loops = 0;
   boundchanges = 1;

   while( boundchanges && loops < MAX_LOOPS )
   {
      loops++;
      boundchanges = 0;

      calcColActivity(scip, matrix, 0, ncols,
         lbdual, ubdual, lbdualbounds, ubdualbounds,
         mincolact, maxcolact, maxcolactinf, mincolactinf);

      for( c = 0 ; c < ncols; c++ )
      {
         SCIP_Real objval;
         SCIP_Real mincolresact;
         SCIP_Bool updateinfcnt;
         SCIP_VAR* var;

         var = SCIPmatrixGetVar(matrix, c);

         /* consider only continuous variables for dual bounds */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
            continue;

         objval = SCIPvarGetObj(var);
         colpnt = SCIPmatrixGetColIdxPtr(matrix, c);
         colend = colpnt + SCIPmatrixGetColNNonzs(matrix, c);
         valpnt = SCIPmatrixGetColValPtr(matrix, c);

         for( ; colpnt < colend; colpnt++, valpnt++ )
         {
            int row;
            SCIP_Real val;
            int part;

            row = *colpnt;
            val = *valpnt;
            mincolresact = -SCIPinfinity(scip);

            part = 0;

            do {

               /* calulate column activity residuals */
               calcColActResidualCommon(scip, matrix, c, row, val,
                  lbdual, ubdual, part, lbdualbounds, ubdualbounds, mincolact, maxcolact,
                  maxcolactinf, mincolactinf, &mincolresact);

               /* update dual bounds */
               updateDualBounds(scip, matrix, objval, val, row, mincolresact,
                  lbdual, ubdual, part, &boundchanges, &updateinfcnt);

               /* update infinity counters if bound changed properly */
               if( updateinfcnt )
                  infCntUpdate(scip, matrix, val, row, lbdual, ubdual,
                     lbdualbounds, ubdualbounds, part, mincolact, maxcolact, maxcolactinf, mincolactinf);

               part++;

            } while( !SCIPmatrixIsRowRhsInfinity(matrix, row) && part < 2 );
         }

         calcColActResidualExplicitBound(scip, matrix, c,
            lbdual, ubdual, lbdualbounds, ubdualbounds, TRUE, FALSE, mincolact, maxcolact,
            maxcolactinf, mincolactinf, &mincolresact);

         updateDualBoundsExplicit(scip, matrix, objval, c, mincolresact,
            lbdualbounds, ubdualbounds, TRUE, FALSE, &boundchanges, &updateinfcnt);

         if( updateinfcnt )
            infCntUpdateExplicit(scip, matrix, c, lbdual, ubdual,
               lbdualbounds, ubdualbounds, TRUE, FALSE, mincolact, maxcolact, maxcolactinf, mincolactinf);

         calcColActResidualExplicitBound(scip, matrix, c,
            lbdual, ubdual, lbdualbounds, ubdualbounds, FALSE, TRUE, mincolact, maxcolact,
            maxcolactinf, mincolactinf, &mincolresact);

         updateDualBoundsExplicit(scip, matrix, objval, c, mincolresact,
            lbdualbounds, ubdualbounds, FALSE, TRUE, &boundchanges, &updateinfcnt);

         if( updateinfcnt )
            infCntUpdateExplicit(scip, matrix, c, lbdual, ubdual,
               lbdualbounds, ubdualbounds, FALSE, TRUE, mincolact, maxcolact, maxcolactinf, mincolactinf);
      }
   }

   if( boundchanges > 0 )
   {
      /* calculate minimal/maximal column activity the last time */
      calcColActivity(scip, matrix, 0, ncols,
         lbdual, ubdual, lbdualbounds, ubdualbounds,
         mincolact, maxcolact, maxcolactinf, mincolactinf);
   }

   for( c = 0; c < ncols; c++ )
   {
      SCIP_Real objval;
      SCIP_VAR* var;

      var = SCIPmatrixGetVar(matrix, c);
      objval = SCIPvarGetObj(var);

      /* positive reduced costs: c_j - sup{(A_{.j})^T} > 0 => x_j = 0 */
      if( SCIPisLT(scip, maxcolact[c], objval) && varstofix[c] == NOFIX )
      {
         varstofix[c] = FIXATLB;
         (*npossiblefixings)++;
      }
   }

   for( r = 0; r < nrows; r++ )
   {
      /* implied equality: y_i > 0 =>  A_{.i}x - b_i = 0 */
      if( SCIPmatrixIsRowRhsInfinity(matrix, r) )
      {
         if( SCIPisGT(scip, lbdual[0][r], 0.0) )
         {
            sidestochange[r] = RHSTOLHS;
            (*npossiblesidechanges)++;
         }
      }
      else
      {
         if( !SCIPmatrixIsRowRhsInfinity(matrix, r) &&
            !SCIPisEQ(scip,SCIPmatrixGetRowLhs(matrix, r),SCIPmatrixGetRowRhs(matrix, r)) )
         {
            if( SCIPisGT(scip, lbdual[0][r], 0.0) )
            {
               assert(sidestochange[r]==NOCHANGE);
               sidestochange[r] = RHSTOLHS;
               (*npossiblesidechanges)++;
            }

            if( SCIPisGT(scip, lbdual[1][r], 0.0) )
            {
               assert(sidestochange[r]==NOCHANGE);
               sidestochange[r] = LHSTORHS;
               (*npossiblesidechanges)++;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &ubdualbounds[1]);
   SCIPfreeBufferArray(scip, &lbdualbounds[1]);
   SCIPfreeBufferArray(scip, &ubdualbounds[0]);
   SCIPfreeBufferArray(scip, &lbdualbounds[0]);
   SCIPfreeBufferArray(scip, &ubdual[1]);
   SCIPfreeBufferArray(scip, &lbdual[1]);
   SCIPfreeBufferArray(scip, &ubdual[0]);
   SCIPfreeBufferArray(scip, &lbdual[0]);
   SCIPfreeBufferArray(scip, &mincolactinf);
   SCIPfreeBufferArray(scip, &maxcolactinf);
   SCIPfreeBufferArray(scip, &maxcolact);
   SCIPfreeBufferArray(scip, &mincolact);
   SCIPfreeBufferArray(scip, &varissingcol);

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyDualinfer)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolDualinfer(scip) );

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecDualinfer)
{  /*lint --e{715}*/
   SCIPMILPMATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) )
      return SCIP_OKAY;

   if( SCIPgetNContVars(scip)==0 )
      return SCIP_OKAY;

   if( !SCIPallowDualReds(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   if( initialized && complete )
   {
      FIXINGDIRECTION* varstofix;
      int npossiblefixings;
      int nconvarsfixed;
      int nintvarsfixed;
      int nbinvarsfixed;
      SIDECHANGE* sidestochange;
      int npossiblesidechanges;
      int nsideschanged;
      int i;
      int nrows;
      int ncols;

      npossiblefixings = 0;
      nconvarsfixed = 0;
      nintvarsfixed = 0;
      nbinvarsfixed = 0;
      npossiblesidechanges = 0;
      nsideschanged = 0;

      nrows = SCIPmatrixGetNRows(matrix);
      ncols = SCIPmatrixGetNColumns(matrix);

      SCIP_CALL( SCIPallocBufferArray(scip, &varstofix, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &sidestochange, nrows) );

      BMSclearMemoryArray(varstofix, ncols);
      BMSclearMemoryArray(sidestochange, nrows);

      SCIP_CALL( dualBoundStrengthening(scip, matrix,
            varstofix, &npossiblefixings, sidestochange, &npossiblesidechanges) );

      if( npossiblefixings > 0 )
      {
         for( i = ncols - 1; i >= 0; --i )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;
            SCIP_VAR* var;

            if( varstofix[i] == FIXATLB )
            {
               SCIP_Real lb;

               var = SCIPmatrixGetVar(matrix, i);
               lb = SCIPvarGetLbLocal(var);

               /* fix at lower bound */
               SCIP_CALL( SCIPfixVar(scip, var, lb, &infeasible, &fixed) );
               if( infeasible )
               {
                  SCIPdebugMessage(" -> infeasible fixing\n");
                  *result = SCIP_CUTOFF;
                  break;
               }
               assert(fixed);
               (*nfixedvars)++;
               *result = SCIP_SUCCESS;

               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
                  nconvarsfixed++;
               else if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
                  nbinvarsfixed++;
               else
                  nintvarsfixed++;
            }
         }
      }


      if( npossiblesidechanges > 0 )
      {
         for( i = 0; i < nrows; i++ )
         {
            SCIP_CONS* cons;
            SCIP_CONSHDLR* conshdlr;
            const char* conshdlrname;

            if( sidestochange[i] == NOCHANGE )
               continue;

            cons = SCIPmatrixGetCons(matrix,i);
            conshdlr = SCIPconsGetHdlr(cons);
            conshdlrname = SCIPconshdlrGetName(conshdlr);

            if( strcmp(conshdlrname, "linear") == 0 )
            {
               SCIP_Real lhs;
               SCIP_Real rhs;
               SCIP_Real matrixlhs;
               SCIP_Real matrixrhs;

               lhs = SCIPgetLhsLinear(scip, cons);
               rhs = SCIPgetRhsLinear(scip, cons);
               matrixlhs = SCIPmatrixGetRowLhs(matrix, i);
               matrixrhs = SCIPmatrixGetRowRhs(matrix, i);

               if( sidestochange[i] == RHSTOLHS )
               {
                  assert(!SCIPisEQ(scip, matrixlhs, matrixrhs));

                  if( SCIPisEQ(scip, matrixlhs, lhs) )
                     SCIP_CALL( SCIPchgRhsLinear(scip, cons, matrixlhs) );
                  else
                     SCIP_CALL( SCIPchgLhsLinear(scip, cons, -matrixlhs) );

                  nsideschanged++;
                  (*nchgsides)++;
               }
               else
               {
                  assert(!SCIPisEQ(scip, matrixlhs, matrixrhs));

                  if( SCIPisEQ(scip, matrixrhs, rhs) )
                     SCIP_CALL( SCIPchgLhsLinear(scip, cons, matrixrhs) );
                  else
                     SCIP_CALL( SCIPchgRhsLinear(scip, cons, -matrixrhs) );

                  nsideschanged++;
                  (*nchgsides)++;
               }
            }
            else
            {
               SCIPdebugMessage("Warning: unsupported conshdlr type for side change!");
            }
         }
      }

      SCIPfreeBufferArray(scip, &sidestochange);
      SCIPfreeBufferArray(scip, &varstofix);

      if( (nconvarsfixed + nintvarsfixed + nbinvarsfixed) > 0 || npossiblesidechanges > 0)
      {
         SCIPdebugMessage("### fixed vars [cont: %d, int: %d, bin: %d], chg sides [%d]\n",
            nconvarsfixed, nintvarsfixed, nbinvarsfixed, nsideschanged);
      }
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the dual inference presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolDualinfer(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   presoldata = NULL;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecDualinfer, presoldata) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyDualinfer) );

   return SCIP_OKAY;
}
