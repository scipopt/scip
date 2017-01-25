/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cuts.c
 * @brief  Methods used to generate and strengthen cuts
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/def.h"
#include "scip/retcode.h"
#include "scip/cuts.h"
#include "scip/debug.h"
#include "scip/pub_lp.h"
#include "scip/lp.h"
#include "scip/prob.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/var.h"
#include "scip/scip.h"
#include "scip/struct_scip.h"


#define MAXCMIRSCALE               1e+6 /**< maximal scaling (scale/(1-f0)) allowed in c-MIR calculations */

/*
 * debug messages
 */

#ifdef SCIP_DEBUG
/** method is to print in row in case SCIP_DEBUG is defined */
static
void debugRowPrint(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   int i;

   assert(row != NULL);

   /* print row name */
   if( row->name != NULL && row->name[0] != '\0' )
   {
      SCIPdebugPrintf("%s: ", row->name);
   }

   /* print left hand side */
   SCIPdebugPrintf("%.15g <= ", row->lhs);

   /* print coefficients */
   if( row->len == 0 )
   {
      SCIPdebugPrintf("0 ");
   }
   for( i = 0; i < row->len; ++i )
   {
      assert(row->cols[i] != NULL);
      assert(row->cols[i]->var != NULL);
      assert(SCIPvarGetName(row->cols[i]->var) != NULL);
      assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
      SCIPdebugPrintf("%+.15g<%s> ", row->vals[i], SCIPvarGetName(row->cols[i]->var));
   }

   /* print constant */
   if( REALABS(row->constant) > SCIP_DEFAULT_EPSILON )
   {
      SCIPdebugPrintf("%+.15g ", row->constant);
   }

   /* print right hand side */
   SCIPdebugPrintf("<= %.15g\n", row->rhs);
}
#else
#define debugRowPrint(x) /**/
#endif

#ifdef SCIP_DEBUG
static
void printMIR(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            mircoef,            /**< MIR coefficients */
   SCIP_Real             mirrhs,             /**< right hand side of the MIR row */
   SCIP_Bool             ignorsol,
   SCIP_Bool             islocal
   )
{
   SCIP_Real activity;
   int i;

   assert(prob != NULL);

   SCIPdebugMessage("MIR:");
   activity = 0.0;
   for( i = 0; i < prob->nvars; ++i )
   {
      if( mircoef[i] != 0.0 )
      {
         SCIPdebugPrintf(" %+g<%s>", mircoef[i], SCIPvarGetName(prob->vars[i]));

         if( !ignorsol )
            activity += mircoef[i] * (sol == NULL ? SCIPvarGetLPSol(prob->vars[i]) : SCIPsolGetVal(sol, set, stat, prob->vars[i]));
         else
         {
            if( mircoef[i] > 0.0 )
            {
               activity += mircoef[i] * (islocal ? SCIPvarGetLbLocal(prob->vars[i]) : SCIPvarGetLbGlobal(prob->vars[i]));
            }
            else
            {
               activity += mircoef[i] * (islocal ? SCIPvarGetUbLocal(prob->vars[i]) : SCIPvarGetUbGlobal(prob->vars[i]));
            }
         }
      }
   }
   SCIPdebugPrintf(" <= %.6f (activity: %g)\n", mirrhs, activity);
}
#endif


/** returns the maximum absolute row weight in the given weight vector, and calculates the sparsity pattern of the weights */
static
SCIP_Real getMaxAbsWeightCalcSparsity(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Real*            weights,            /**< row weights in row summation */
   int*                  rowinds,            /**< array to store sparsity pattern of used rows; size lp->nrows */
   int*                  nrowinds,           /**< pointer to store number of used rows */
   int*                  rowlensum           /**< pointer to store total number of non-zeros in used rows */
   )
{
   SCIP_Real maxabsweight;
   int r;

   assert(set != NULL);
   assert(lp != NULL);
   assert(weights != NULL);
   assert(rowinds != NULL);
   assert(nrowinds != NULL);

   *nrowinds = 0;
   *rowlensum = 0;

   maxabsweight = 0.0;
   for( r = 0; r < lp->nrows; ++r )
   {
      SCIP_Real absweight;

      /* skip unused rows */
      if( SCIPsetIsZero(set, weights[r]) )
         continue;

      /* record the row in the sparsity pattern */
      rowinds[*nrowinds] = r;
      (*nrowinds)++;

      (*rowlensum) += SCIProwGetNNonz(lp->rows[r]);

      absweight = REALABS(weights[r]);
      maxabsweight = MAX(maxabsweight, absweight);
   }

   return maxabsweight;
}

/** returns the maximum absolute row weight in the given weight vector using given sparsity pattern */
static
SCIP_Real getMaxAbsWeight(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Real*            weights,            /**< row weights in row summation */
   int*                  rowinds,            /**< array of sparsity pattern of used rows; size lp->nrows */
   int*                  nrowinds,           /**< pointer to store number of used rows */
   int*                  rowlensum           /**< pointer to store total number of non-zeros in used rows */
   )
{
   SCIP_Real maxabsweight;
   int r;   /* index used for reading from the row*/
   int w;   /* auxiliary index to skip zeros in weights array */

   assert(set != NULL);
   assert(lp != NULL);
   assert(weights != NULL);
   assert(rowinds != NULL);
   assert(nrowinds != NULL);

   *rowlensum = 0;

   maxabsweight = 0.0;
   w = 0;
   for( r = 0; r < *nrowinds; ++r )
   {
      SCIP_Real absweight;

      /* remove zeros from the sparsity pattern */
      if( SCIPsetIsZero(set, weights[rowinds[r]]) )
         continue;

      rowinds[w] = rowinds[r];
      ++w;

      (*rowlensum) += SCIProwGetNNonz(lp->rows[rowinds[r]]);

      absweight = REALABS(weights[rowinds[r]]);
      maxabsweight = MAX(maxabsweight, absweight);
   }
   (*nrowinds) = w;

   return maxabsweight;
}

/** finds the best lower bound of the variable to use for MIR transformation */
static
void findBestLb(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            bestlb,             /**< pointer to store best bound value */
   int*                  bestlbtype          /**< pointer to store best bound type */
   )
{
   assert(bestlb != NULL);
   assert(bestlbtype != NULL);

   *bestlb = SCIPvarGetLbGlobal(var);
   *bestlbtype = -1;

   if( allowlocal )
   {
      SCIP_Real loclb;

      loclb = SCIPvarGetLbLocal(var);
      if( SCIPsetIsGT(set, loclb, *bestlb) )
      {
         *bestlb = loclb;
         *bestlbtype = -2;
      }
   }

   if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
   {
      SCIP_Real bestvlb;
      int bestvlbidx;

      SCIPvarGetClosestVlb(var, sol, set, stat, &bestvlb, &bestvlbidx);
      if( bestvlbidx >= 0
         && (bestvlb > *bestlb || (*bestlbtype < 0 && SCIPsetIsGE(set, bestvlb, *bestlb))) )
      {
         SCIP_VAR** vlbvars;

         /* we have to avoid cyclic variable bound usage, so we enforce to use only variable bounds variables of smaller index */
         /**@todo this check is not needed for continuous variables; but allowing all but binary variables
          *       to be replaced by variable bounds seems to be buggy (wrong result on gesa2)
          */
         vlbvars = SCIPvarGetVlbVars(var);
         assert(vlbvars != NULL);
         if( SCIPvarGetProbindex(vlbvars[bestvlbidx]) < SCIPvarGetProbindex(var) )
         {
            *bestlb = bestvlb;
            *bestlbtype = bestvlbidx;
         }
      }
   }
}

/** finds the best upper bound of the variable to use for MIR transformation */
static
void findBestUb(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            bestub,             /**< pointer to store best bound value */
   int*                  bestubtype          /**< pointer to store best bound type */
   )
{
   assert(bestub != NULL);
   assert(bestubtype != NULL);

   *bestub = SCIPvarGetUbGlobal(var);
   *bestubtype = -1;

   if( allowlocal )
   {
      SCIP_Real locub;

      locub = SCIPvarGetUbLocal(var);
      if( SCIPsetIsLT(set, locub, *bestub) )
      {
         *bestub = locub;
         *bestubtype = -2;
      }
   }

   if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
   {
      SCIP_Real bestvub;
      int bestvubidx;

      SCIPvarGetClosestVub(var, sol, set, stat, &bestvub, &bestvubidx);
      if( bestvubidx >= 0
         && (bestvub < *bestub || (*bestubtype < 0 && SCIPsetIsLE(set, bestvub, *bestub))) )
      {
         SCIP_VAR** vubvars;

         /* we have to avoid cyclic variable bound usage, so we enforce to use only variable bounds variables of smaller index */
         /**@todo this check is not needed for continuous variables; but allowing all but binary variables
          *       to be replaced by variable bounds seems to be buggy (wrong result on gesa2)
          */
         vubvars = SCIPvarGetVubVars(var);
         assert(vubvars != NULL);
         if( SCIPvarGetProbindex(vubvars[bestvubidx]) < SCIPvarGetProbindex(var) )
         {
            *bestub = bestvub;
            *bestubtype = bestvubidx;
         }
      }
   }
}

/** calculates the activity of the given MIR cut */
static
SCIP_Real getMIRRowActivity(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   int*                  varinds,            /**< sparsity pattern of non-zero MIR coefficients */
   int                   nvarinds            /**< number of non-zero MIR coefficients */
   )
{
   SCIP_Real act;
   int i;

   act = 0.0;
   for( i = 0; i < nvarinds; i++ )
   {
      int v;

      v = varinds[i];
      assert(0 <= v && v < prob->nvars);
      assert(mircoef[v] != 0.0);

      act += mircoef[v] * ( sol == NULL ? SCIPvarGetLPSol(prob->vars[v]) : SCIPsolGetVal(sol, set, stat, prob->vars[v]));
   }

   return act;
}

/** calculates the minimal activity of the given MIR */
static
SCIP_Real getMIRMinActivity(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   int*                  varinds,            /**< sparsity pattern of non-zero MIR coefficients */
   int                   nvarinds,           /**< number of non-zero MIR coefficients */
   SCIP_Bool             islocal             /**< whether local bounds should be used */
   )
{
   SCIP_Real act = 0.0;
   int i;

   for( i = 0; i < nvarinds; i++ )
   {
      int v;

      v = varinds[i];
      assert(0 <= v && v < prob->nvars);
      assert(mircoef[v] != 0.0);

      if( mircoef[v] > 0.0 )
      {
         act += mircoef[v] * (islocal ? SCIPvarGetLbLocal(prob->vars[v]) : SCIPvarGetLbGlobal(prob->vars[v]));
      }
      else
      {
         act += mircoef[v] * (islocal ? SCIPvarGetUbLocal(prob->vars[v]) : SCIPvarGetUbGlobal(prob->vars[v]));
      }
   }

   return act;
}

/** adds a single row to an aggregation */
static
void addRowToAggregation(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            mircoef,            /**< array of aggregation coefficients: must be of size prob->nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   SCIP_Bool*            varused,            /**< array to flag variables that appear in the MIR constraint; size prob->nvars */
   int*                  varinds,            /**< array to store sparsity pattern of non-zero MIR coefficients; size prob->nvars */
   int*                  nvarinds,           /**< pointer to store number of non-zero MIR coefficients */
   SCIP_ROW*             row,                /**< row to add to the aggregation */
   SCIP_Real             weight,             /**< weight of row in aggregation */
   SCIP_Bool             uselhs              /**< TRUE if lhs should be used, FALSE if rhs should be used */
   )
{
   SCIP_Real sideval;
   int r;
   int i;

   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(slacksign != NULL);
   assert(varused != NULL);
   assert(varinds != NULL);
   assert(nvarinds != NULL);
   assert(row != NULL);
   assert(weight != 0.0);

   r = row->lppos;
   assert(r >= 0);

   /* update the right hand side */
   if( uselhs )
   {
      slacksign[r] = -1;
      sideval = row->lhs - row->constant;
      if( row->integral )
         sideval = SCIPsetFeasCeil(set, sideval); /* row is integral: round left hand side up */
   }
   else
   {
      slacksign[r] = +1;
      sideval = row->rhs - row->constant;
      if( row->integral )
         sideval = SCIPsetFeasFloor(set, sideval); /* row is integral: round right hand side up */
   }
   (*mirrhs) += weight * sideval;

   /* add the row coefficients to the sum */
   for( i = 0; i < row->len; ++i )
   {
      int idx;

      assert(row->cols[i] != NULL);
      assert(row->cols[i]->var != NULL);
      assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(row->cols[i]->var) == row->cols[i]);
      assert(SCIPvarGetProbindex(row->cols[i]->var) == row->cols[i]->var_probindex);
      idx = row->cols[i]->var_probindex;
      mircoef[idx] += weight * row->vals[i];

      /* record the variable in the sparsity pattern */
      if( !varused[idx] )
      {
         varused[idx] = TRUE;
         varinds[*nvarinds] = idx;
         (*nvarinds)++;
      }
   }
}

/** builds a weighted sum of rows, and decides whether to use the left or right hand side of the rows in summation */
static
void cutsSumStrongCGRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Bool             allowlocal,         /**< should local rows be included, resulting in a locally valid summation? */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size prob->nvars */
   SCIP_Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   SCIP_Bool*            varused,            /**< array to flag variables that appear in the MIR constraint; size prob->nvars */
   int*                  varinds,            /**< array to store sparsity pattern of non-zero MIR coefficients; size prob->nvars */
   int*                  nvarinds,           /**< pointer to store number of non-zero MIR coefficients */
   int*                  rowinds,            /**< array to store sparsity pattern of used rows; size lp->nrows */
   int*                  nrowinds,           /**< pointer to store number of used rows */
   SCIP_Bool*            emptyrow,           /**< pointer to store whether the returned row is empty */
   SCIP_Bool*            localrowsused,      /**< pointer to store whether local rows were used in summation */
   SCIP_Bool*            rowtoolong,         /**< pointer to store whether the aggregated row is too long and thus invalid */
   int*                  cutrank             /**< pointer to store the rank of the returned aggregation; or NULL */
   )
{
   SCIP_Real maxweight;
   int rowlensum;
   int maxrank = 0;
   int i;

   assert(prob != NULL);
   assert(lp != NULL);
   assert(weights != NULL);
   assert(!SCIPsetIsZero(set, scale));
   assert(maxweightrange >= 1.0);
   assert(strongcgcoef != NULL);
   assert(strongcgrhs != NULL);
   assert(slacksign != NULL);
   assert(varused != NULL);
   assert(varinds != NULL);
   assert(nvarinds != NULL);
   assert(rowinds != NULL);
   assert(nrowinds != NULL);
   assert(emptyrow != NULL);
   assert(localrowsused != NULL);

   *nvarinds = 0;
   *localrowsused = FALSE;
   *strongcgrhs = 0.0;
   *emptyrow = TRUE;
   *rowtoolong = FALSE;

   /* initialize varused array */
   BMSclearMemoryArray(varused, prob->nvars);

   /* search the maximal absolute weight and calculate the row sparsity pattern */
   if( *nrowinds == -1 )
      maxweight = getMaxAbsWeightCalcSparsity(set, lp, weights, rowinds, nrowinds, &rowlensum);
   else
      maxweight = getMaxAbsWeight(set, lp, weights, rowinds, nrowinds, &rowlensum);


   maxweight *= ABS(scale);

   /* if the total number of non-zeros is way too large, we just skip this aggregation */
   if( rowlensum/5 > maxmksetcoefs )
   {
      *rowtoolong = TRUE;
      return;
   }

   /* calculate the row summation */
   BMSclearMemoryArray(strongcgcoef, prob->nvars);
   i = 0;
   while( i < *nrowinds )
   {
      SCIP_ROW* row;
      SCIP_Real weight;
      SCIP_Real absweight;
      SCIP_Bool skiprow;
      int r;

      r = rowinds[i];
      assert(0 <= i && i < lp->nrows);
      assert(weights[r] != 0.0);

      row = lp->rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* modifiable rows cannot be part of a strong CG row summation;
       * local rows are only included, if the allowlocal flag is set;
       * close to zero weights or weights outside the maximal range are ignored
       */
      weight = scale * weights[r];
      absweight = ABS(weight);
      skiprow = FALSE;
      if( !row->modifiable && (allowlocal || !row->local)
         && absweight * maxweightrange >= maxweight && !SCIPsetIsSumZero(set, weight) )
      {
         /*lint --e{644}*/
         SCIP_Bool uselhs;

         if( row->integral )
         {
            /* Row is integral:
             * Decide, if we want to use the left or the right hand side of the row in the summation.
             * If possible, use the side that leads to a positive slack value in the summation.
             */
            if( SCIPsetIsInfinity(set, row->rhs) || (!SCIPsetIsInfinity(set, -row->lhs) && weight < 0.0) )
               uselhs = TRUE;
            else
               uselhs = FALSE;
         }
         else
         {
            /* Row is NOT integral:
             * Decide, if we have to use the left or the right hand side of the row in the summation,
             * in order to get a positive slack variable in the summation.
             * If not possible, ignore row in summation.
             */
            if( weight < 0.0 && !SCIPsetIsInfinity(set, -row->lhs) )
               uselhs = TRUE;
            else if( weight > 0.0 && !SCIPsetIsInfinity(set, row->rhs) )
               uselhs = FALSE;
            else
               skiprow = TRUE;
         }

         if( !skiprow )
         {
            /* add the row to the aggregation */
            addRowToAggregation(set, strongcgcoef, strongcgrhs, slacksign, varused, varinds, nvarinds, row, weight, uselhs);
            *emptyrow = FALSE;
            *localrowsused = *localrowsused || row->local;

            maxrank = MAX(maxrank, row->rank);

            SCIPdebugMessage("strong CG: %d: row <%s>, lhs = %g, rhs = %g, scale = %g, weight = %g, slacksign = %d -> rhs = %g\n",
               r, SCIProwGetName(row), row->lhs - row->constant, row->rhs - row->constant,
               scale, weights[r], slacksign[r], *strongcgrhs);
            debugRowPrint(row);
         }
      }
      else
         skiprow = TRUE;

      if( skiprow )
      {
         /* remove row from sparsity pattern, do not increase i, since the i-th position is filled with the last element */
         rowinds[i] = rowinds[(*nrowinds)-1];
         (*nrowinds)--;
#ifndef NDEBUG
         slacksign[r] = 0;
#endif
      }
      else
         ++i;
   }

   /* check if the total number of non-zeros is too large */
   if( *nrowinds > maxmksetcoefs )
      *rowtoolong = TRUE;

   /* set rank of the cut */
   if( cutrank != NULL )
      *cutrank = maxrank + 1;
}

/** builds a weighted sum of rows, and decides whether to use the left or right hand side of the rows in summation */
static
void cutsSumMIRRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_Real             knownmaxweight,     /**< largest magnitude of weights. Set to 0 if compress == TRUE */
   int*                  sidetypes,          /**< specify row side type (-1 = lhs, 0 = unkown, 1 = rhs) or NULL for automatic choices */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Bool             allowlocal,         /**< should local rows be included, resulting in a locally valid summation? */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Bool             compress,           /**< if rowinds is unknown and weights should be compressed */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size prob->nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   SCIP_Bool*            varused,            /**< array to flag variables that appear in the MIR constraint; size prob->nvars */
   int*                  varinds,            /**< array to store sparsity pattern of non-zero MIR coefficients; size prob->nvars */
   int*                  nvarinds,           /**< pointer to store number of non-zero MIR coefficients */
   int*                  rowinds,            /**< array to store sparsity pattern of used rows; size lp->nrows */
   int*                  nrowinds,           /**< pointer to store number of used rows */
   SCIP_Bool*            emptyrow,           /**< pointer to store whether the returned row is empty */
   SCIP_Bool*            localrowsused,      /**< pointer to store whether local rows were used in summation */
   SCIP_Bool*            rowtoolong,         /**< pointer to store whether the aggregated row is too long and thus invalid */
   int*                  cutrank             /**< pointer to store the rank of the returned aggregation; or NULL */
   )
{
   SCIP_Real maxweight;
   int maxrank = 0;
   int rowlensum;
   int i;

   assert(prob != NULL);
   assert(lp != NULL);
   assert(weights != NULL);
   assert(!SCIPsetIsZero(set, scale));
   assert(maxweightrange >= 1.0);
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(slacksign != NULL);
   assert(varused != NULL);
   assert(varinds != NULL);
   assert(nvarinds != NULL);
   assert(rowinds != NULL);
   assert(nrowinds != NULL);
   assert(emptyrow != NULL);
   assert(localrowsused != NULL);
   assert(rowtoolong != NULL);

   rowlensum = 0;
   *nvarinds = 0;
   *mirrhs = 0.0;
   *emptyrow = TRUE;
   *localrowsused = FALSE;
   *rowtoolong = FALSE;

   /* initialize varused array */
   BMSclearMemoryArray(varused, prob->nvars);

   /* if compression of the dense weight vector is required  */
   /* search the maximal absolute weight and calculate the row sparsity pattern */
   if( compress )
   {
      maxweight = getMaxAbsWeightCalcSparsity(set, lp, weights, rowinds, nrowinds, &rowlensum);
   }
   else
   {
      /* search the maximal absolute weight using the given row sparsity pattern */
      if( knownmaxweight == -1 )
      {
         assert(*nrowinds > -1);
         maxweight = getMaxAbsWeight(set, lp, weights, rowinds, nrowinds, &rowlensum);
      }
      else
         maxweight = knownmaxweight;
   }

   /* if the total number of non-zeros is way too large, we just skip this aggregation */
   if( rowlensum/5 > maxmksetcoefs )
   {
      *rowtoolong = TRUE;
      return;
   }

   maxweight *= ABS(scale);

   /* calculate the row summation */
   BMSclearMemoryArray(mircoef, prob->nvars);
   i = 0;
   while( i < *nrowinds )
   {
      SCIP_ROW* row;
      SCIP_Real weight;
      SCIP_Real absweight;
      int r;

      r = rowinds[i];
      assert(0 <= r && r < lp->nrows);
      assert(weights[r] != 0.0);

      row = lp->rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* modifiable rows cannot be part of a MIR row summation;
       * local rows are only included, if the allowlocal flag is set;
       * close to zero weights or weights outside the maximal range are ignored
       */
      weight = scale * weights[r];
      absweight = REALABS(weight);
      if( !row->modifiable && (allowlocal || !row->local)
         && absweight * maxweightrange >= maxweight && !SCIPsetIsSumZero(set, weight) )
      {
         SCIP_Bool uselhs;

         /* choose sides for lhs/rhs of row */
         if ( sidetypes != NULL )
         {
            assert( sidetypes[r] == -1 || sidetypes[r] == 0 || sidetypes[r] == 1 );
            if ( sidetypes[r] == -1 )
            {
               assert( ! SCIPsetIsInfinity(set, -row->lhs) );
               uselhs = TRUE;
            }
            else if ( sidetypes[r] == 1 )
            {
               assert( ! SCIPsetIsInfinity(set, row->rhs) );
               uselhs = FALSE;
            }
            else
            {
               /* Automatically decide, whether we want to use the left or the right hand side of the row in the summation.
                * If possible, use the side that leads to a positive slack value in the summation.
                */
               if( SCIPsetIsInfinity(set, row->rhs) || (!SCIPsetIsInfinity(set, -row->lhs) && weight < 0.0) )
                  uselhs = TRUE;
               else
                  uselhs = FALSE;
            }
         }
         else
         {
            /* Automatically decide, whether we want to use the left or the right hand side of the row in the summation.
             * If possible, use the side that leads to a positive slack value in the summation.
             */
            if( SCIPsetIsInfinity(set, row->rhs) || (!SCIPsetIsInfinity(set, -row->lhs) && weight < 0.0) )
               uselhs = TRUE;
            else
               uselhs = FALSE;
         }

         /* add the row to the aggregation */
         addRowToAggregation(set, mircoef, mirrhs, slacksign, varused, varinds, nvarinds, row, weight, uselhs);
         *emptyrow = FALSE;
         *localrowsused = *localrowsused || row->local;

         SCIPdebugMessage("MIR: %d: row <%s>, lhs = %g, rhs = %g, scale = %g, weight = %g, slacksign = %d -> rhs = %g\n",
            r, SCIProwGetName(row), row->lhs - row->constant, row->rhs - row->constant,
            scale, weights[r], slacksign[r], *mirrhs);
         debugRowPrint(row);

         /* update the rank of the aggregation */
         maxrank = MAX(maxrank, row->rank);

         ++i; /* handle next row */
      }
      else
      {
         /* remove row from sparsity pattern, do not increase i (i-th position is filled with last entry) */
         rowinds[i] = rowinds[(*nrowinds)-1];
         (*nrowinds)--;
#ifndef NDEBUG
         slacksign[r] = 0;
#endif
      }
   }

   /* check if the total number of non-zeros is too large */
   if( *nvarinds > maxmksetcoefs )
      *rowtoolong = TRUE;

   /* set rank of the aggregated cut */
   if( cutrank != NULL )
      *cutrank = maxrank + 1;
}

/** removes all nearly-zero coefficients from MIR row and relaxes the right hand side correspondingly in order to
 *  prevent numerical rounding errors
 */
static
void cutsCleanupMIRRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   SCIP_Bool*            varused,            /**< array to flag variables that appear in the MIR constraint */
   int*                  varinds,            /**< sparsity pattern of non-zero MIR coefficients */
   int*                  nvarinds,           /**< pointer to number of non-zero MIR coefficients */
   SCIP_Bool             cutislocal          /**< is the cut only valid locally? */
   )
{
   SCIP_Bool rhsinf;
   int i;

   assert(prob != NULL);
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(varused != NULL);
   assert(varinds != NULL);
   assert(nvarinds != NULL);

   rhsinf = SCIPsetIsInfinity(set, *mirrhs);
   i = 0;
   while( i < *nvarinds )
   {
      int v;

      v = varinds[i];
      assert(0 <= v && v < prob->nvars);
      assert(varused[v]);

      if( SCIPsetIsSumZero(set, mircoef[v]) )
      {
         SCIP_Real bd;

         SCIPdebugMessage("coefficient of <%s> in transformed MIR row is too small: %.12f\n",
            SCIPvarGetName(prob->vars[v]), mircoef[v]);

         /* relax the constraint such that the coefficient becomes exactly 0.0 */
         if( SCIPsetIsPositive(set, mircoef[v]) )
         {
            bd = cutislocal ? SCIPvarGetLbLocal(prob->vars[v]) : SCIPvarGetLbGlobal(prob->vars[v]);
            rhsinf = SCIPsetIsInfinity(set, -bd);
         }
         else if( SCIPsetIsNegative(set, mircoef[v]) )
         {
            bd = cutislocal ? SCIPvarGetUbLocal(prob->vars[v]) : SCIPvarGetUbGlobal(prob->vars[v]);
            rhsinf = SCIPsetIsInfinity(set, bd);
         }
         else
            bd = 0.0;
         *mirrhs -= bd * mircoef[v];
         mircoef[v] = 0.0;

         /* remove variable from sparsity pattern, do not increase i (i-th position is filled with last entry) */
         varused[v] = FALSE;
         varinds[i] = varinds[(*nvarinds)-1];
         (*nvarinds)--;
      }
      else
         ++i;
   }
   if( rhsinf )
      *mirrhs = SCIPsetInfinity(set);
}

/** Transform equation  \f$ a*x == b \f$, \f$ lb <= x <= ub \f$ into standard form \f$ a^\prime*x^\prime == b\f$, \f$ 0 <= x^\prime <= ub^\prime \f$.
 *
 *  Transform variables (lb or ub):
 * \f[
 * \begin{array}{llll}
 *    x^\prime_j := x_j - lb_j,&   x_j == x^\prime_j + lb_j,&   a^\prime_j ==  a_j,&   \mbox{if lb is used in transformation} \\
 *    x^\prime_j := ub_j - x_j,&   x_j == ub_j - x^\prime_j,&   a^\prime_j == -a_j,&   \mbox{if ub is used in transformation}
 * \end{array}
 * \f]
 *  and move the constant terms \f$ a_j * lb_j \f$ or \f$ a_j * ub_j \f$ to the rhs.
 *
 *  Transform variables (vlb or vub):
 * \f[
 * \begin{array}{llll}
 *    x^\prime_j := x_j - (bl_j * zl_j + dl_j),&   x_j == x^\prime_j + (bl_j * zl_j + dl_j),&   a^\prime_j ==  a_j,&   \mbox{if vlb is used in transf.} \\
 *    x^\prime_j := (bu_j * zu_j + du_j) - x_j,&   x_j == (bu_j * zu_j + du_j) - x^\prime_j,&   a^\prime_j == -a_j,&   \mbox{if vub is used in transf.}
 * \end{array}
 * \f]
 *  move the constant terms \f$ a_j * dl_j \f$ or \f$ a_j * du_j \f$ to the rhs, and update the coefficient of the VLB variable:
 * \f[
 * \begin{array}{ll}
 *    a_{zl_j} := a_{zl_j} + a_j * bl_j,& \mbox{or} \\
 *    a_{zu_j} := a_{zu_j} + a_j * bu_j&
 * \end{array}
 * \f]
 */
static
void cutsTransformStrongCGRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size nvars */
   SCIP_Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   SCIP_Bool*            varused,            /**< array to flag variables that appear in the MIR constraint */
   int*                  varinds,            /**< sparsity pattern of non-zero MIR coefficients */
   int*                  nvarinds,           /**< pointer to number of non-zero MIR coefficients */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Bool*            freevariable,       /**< stores whether a free variable was found in strong CG row -> invalid summation */
   SCIP_Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{
   int i;

   assert(prob != NULL);
   assert(strongcgcoef != NULL);
   assert(strongcgrhs != NULL);
   assert(varused != NULL);
   assert(varinds != NULL);
   assert(nvarinds != NULL);
   assert(varsign != NULL);
   assert(boundtype != NULL);
   assert(freevariable != NULL);
   assert(localbdsused != NULL);

   *freevariable = FALSE;
   *localbdsused = FALSE;

#ifndef NDEBUG
   /* in debug mode, make sure that the whole array is initialized with invalid values */
   for( i = 0; i < prob->nvars; i++ )
   {
      varsign[i] = 0;
      boundtype[i] = -1;
   }
#endif

   /* start with continuous variables, because using variable bounds can affect the untransformed integral
    * variables, and these changes have to be incorporated in the transformation of the integral variables
    * (continuous variables have largest problem indices!)
    */
   SCIPsortDownInt(varinds, *nvarinds);

   /* substitute continuous variables with best standard or variable bound (lb, ub, vlb or vub),
    * substitute integral variables with best standard bound (lb, ub)
    */
   i = 0;
   while( i < *nvarinds )
   {
      SCIP_VAR* var;
      SCIP_Real varsol;
      SCIP_Real bestlb;
      SCIP_Real bestub;
      SCIP_Bool uselb;
      int bestlbtype;
      int bestubtype;
      int v;

      v = varinds[i];
      assert(0 <= v && v < prob->nvars);
      assert(varused[v]);

      var = prob->vars[v];
      assert(v == SCIPvarGetProbindex(var));

      /* due to variable bound usage cancellation may occur;
       * do not increase i, since last element is copied to the i-th position
       */
      if( SCIPsetIsZero(set, strongcgcoef[v]) )
      {
         varsign[v] = +1;
         boundtype[v] = -1;
         strongcgcoef[v] = 0.0;
         varused[v] = FALSE;
         varinds[i] = varinds[(*nvarinds)-1];
         (*nvarinds)--;
         continue;
      }

      /* find closest lower bound in standard lower bound (and variable lower bounds for continuous variables) */
      bestlb = SCIPvarGetLbGlobal(var);
      bestlbtype = -1;
      if( allowlocal )
      {
         SCIP_Real loclb;

         loclb = SCIPvarGetLbLocal(var);
         if( SCIPsetIsGT(set, loclb, bestlb) )
         {
            bestlb = loclb;
            bestlbtype = -2;
         }
      }
      if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_Real bestvlb;
         int bestvlbidx;

         SCIPvarGetClosestVlb(var, NULL, set, stat, &bestvlb, &bestvlbidx);
         if( bestvlbidx >= 0
            && (bestvlb > bestlb || (bestlbtype < 0 && SCIPsetIsGE(set, bestvlb, bestlb))) )
         {
            bestlb = bestvlb;
            bestlbtype = bestvlbidx;
         }
      }

      /* find closest upper bound in standard upper bound (and variable upper bounds for continuous variables) */
      bestub = SCIPvarGetUbGlobal(var);
      bestubtype = -1;
      if( allowlocal )
      {
         SCIP_Real locub;

         locub = SCIPvarGetUbLocal(var);
         if( SCIPsetIsLT(set, locub, bestub) )
         {
            bestub = locub;
            bestubtype = -2;
         }
      }
      if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_Real bestvub;
         int bestvubidx;

         SCIPvarGetClosestVub(var, NULL, set, stat, &bestvub, &bestvubidx);
         if( bestvubidx >= 0
            && (bestvub < bestub || (bestubtype < 0 && SCIPsetIsLE(set, bestvub, bestub))) )
         {
            bestub = bestvub;
            bestubtype = bestvubidx;
         }
      }

      /* check, if variable is free variable
       * (for continuous variable use bound, so that coefficient will be nonnegative)
       */
      if( (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS
            && SCIPsetIsInfinity(set, -bestlb)
            && SCIPsetIsInfinity(set, bestub))
         || (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
            && ((strongcgcoef[v] > 0.0 && SCIPsetIsInfinity(set, -bestlb))
               || (strongcgcoef[v] < 0.0 && SCIPsetIsInfinity(set, bestub)))) )
      {
         /* we found a free variable in the row with non-zero coefficient
          *  -> strong CG row can't be transformed in standard form
          */
         *freevariable = TRUE;
         return;
      }

      /* select transformation bound
       * (for continuous variable use bound, so that coefficient will be nonnegative)
       */
      varsol = SCIPvarGetLPSol(var);
      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         uselb = (strongcgcoef[v] > 0.0);
      else if( SCIPsetIsLT(set, varsol, (1.0 - boundswitch) * bestlb + boundswitch * bestub) )
         uselb = TRUE;
      else if( SCIPsetIsGT(set, varsol, (1.0 - boundswitch) * bestlb + boundswitch * bestub) )
         uselb = FALSE;
      else if( bestlbtype == -1 )  /* prefer global standard bounds */
         uselb = TRUE;
      else if( bestubtype == -1 )  /* prefer global standard bounds */
         uselb = FALSE;
      else if( bestlbtype >= 0 )   /* prefer variable bounds over local bounds */
         uselb = TRUE;
      else if( bestubtype >= 0 )   /* prefer variable bounds over local bounds */
         uselb = FALSE;
      else
         uselb = TRUE;             /* no decision yet? just use lower bound */

      /* perform bound substitution */
      if( uselb )
      {
         assert(!SCIPsetIsInfinity(set, -bestlb));

         /* use lower bound as transformation bound: x'_j := x_j - lb_j */
         boundtype[v] = bestlbtype;
         varsign[v] = +1;

         /* standard (bestlbtype < 0) or variable (bestlbtype >= 0) lower bound? */
         if( bestlbtype < 0 )
         {
            (*strongcgrhs) -= strongcgcoef[v] * bestlb;
            *localbdsused = *localbdsused || (bestlbtype == -2);
         }
         else
         {
            SCIP_VAR** vlbvars = SCIPvarGetVlbVars(var);
            SCIP_Real* vlbcoefs = SCIPvarGetVlbCoefs(var);
            SCIP_Real* vlbconsts = SCIPvarGetVlbConstants(var);
            int zidx;

            assert(0 <= bestlbtype && bestlbtype < SCIPvarGetNVlbs(var));
            assert(SCIPvarIsActive(vlbvars[bestlbtype]));
            zidx = SCIPvarGetProbindex(vlbvars[bestlbtype]);
            assert(0 <= zidx && zidx < v);

            (*strongcgrhs) -= strongcgcoef[v] * vlbconsts[bestlbtype];
            strongcgcoef[zidx] += strongcgcoef[v] * vlbcoefs[bestlbtype];

            /* update sparsity pattern */
            if( !varused[zidx] )
            {
               assert(*nvarinds < prob->nvars);
               varused[zidx] = TRUE;
               varinds[*nvarinds] = zidx;
               (*nvarinds)++;
            }
         }
      }
      else
      {
         assert(!SCIPsetIsInfinity(set, bestub));

         /* use upper bound as transformation bound: x'_j := ub_j - x_j */
         boundtype[v] = bestubtype;
         varsign[v] = -1;

         /* standard (bestubtype < 0) or variable (bestubtype >= 0) upper bound? */
         if( bestubtype < 0 )
         {
            (*strongcgrhs) -= strongcgcoef[v] * bestub;
            *localbdsused = *localbdsused || (bestubtype == -2);
         }
         else
         {
            SCIP_VAR** vubvars = SCIPvarGetVubVars(var);
            SCIP_Real* vubcoefs = SCIPvarGetVubCoefs(var);
            SCIP_Real* vubconsts = SCIPvarGetVubConstants(var);
            int zidx;

            assert(0 <= bestubtype && bestubtype < SCIPvarGetNVubs(var));
            assert(SCIPvarIsActive(vubvars[bestubtype]));
            zidx = SCIPvarGetProbindex(vubvars[bestubtype]);
            assert(zidx >= 0);

            (*strongcgrhs) -= strongcgcoef[v] * vubconsts[bestubtype];
            strongcgcoef[zidx] += strongcgcoef[v] * vubcoefs[bestubtype];

            /* update sparsity pattern */
            if( !varused[zidx] )
            {
               assert(*nvarinds < prob->nvars);
               varused[zidx] = TRUE;
               varinds[*nvarinds] = zidx;
               (*nvarinds)++;
            }
         }
      }

#ifndef NDEBUG
      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         assert(strongcgcoef[v]*varsign[v] > 0.0);
#endif

      SCIPdebugMessage("strong CG var <%s>: varsign=%d, boundtype=%d, strongcgcoef=%g, lb=%g, ub=%g -> rhs=%g\n",
         SCIPvarGetName(var), varsign[v], boundtype[v], strongcgcoef[v], bestlb, bestub, *strongcgrhs);

      ++i; /*increase iterator */
   }
}

/** Transform equation \f$ a \cdot x = b; lb \leq x \leq ub \f$ into standard form
 *    \f$ a^\prime \cdot x^\prime = b,\; 0 \leq x^\prime \leq ub' \f$.
 *
 *  Transform variables (lb or ub):
 *  \f[
 *  \begin{array}{llll}
 *    x^\prime_j := x_j - lb_j,&   x_j = x^\prime_j + lb_j,&   a^\prime_j =  a_j,&   \mbox{if lb is used in transformation}\\
 *    x^\prime_j := ub_j - x_j,&   x_j = ub_j - x^\prime_j,&   a^\prime_j = -a_j,&   \mbox{if ub is used in transformation}
 *  \end{array}
 *  \f]
 *  and move the constant terms \f$ a_j\, lb_j \f$ or \f$ a_j\, ub_j \f$ to the rhs.
 *
 *  Transform variables (vlb or vub):
 *  \f[
 *  \begin{array}{llll}
 *    x^\prime_j := x_j - (bl_j\, zl_j + dl_j),&   x_j = x^\prime_j + (bl_j\, zl_j + dl_j),&   a^\prime_j =  a_j,&   \mbox{if vlb is used in transf.} \\
 *    x^\prime_j := (bu_j\, zu_j + du_j) - x_j,&   x_j = (bu_j\, zu_j + du_j) - x^\prime_j,&   a^\prime_j = -a_j,&   \mbox{if vub is used in transf.}
 *  \end{array}
 *  \f]
 *  move the constant terms \f$ a_j\, dl_j \f$ or \f$ a_j\, du_j \f$ to the rhs, and update the coefficient of the VLB variable:
 *  \f[
 *  \begin{array}{ll}
 *    a_{zl_j} := a_{zl_j} + a_j\, bl_j,& \mbox{or} \\
 *    a_{zu_j} := a_{zu_j} + a_j\, bu_j &
 *  \end{array}
 *  \f]
 */
static
SCIP_RETCODE cutsTransformMIRRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   SCIP_Bool             ignoresol,          /**< should the LP solution be ignored? (eg, apply MIR to dualray) */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   SCIP_Bool*            varused,            /**< array to flag variables that appear in the MIR constraint */
   int*                  varinds,            /**< sparsity pattern of non-zero MIR coefficients */
   int*                  nvarinds,           /**< pointer to number of non-zero MIR coefficients */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Bool*            freevariable,       /**< stores whether a free variable was found in MIR row -> invalid summation */
   SCIP_Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{
   SCIP_Real* bestlbs;
   SCIP_Real* bestubs;
   int* bestlbtypes;
   int* bestubtypes;
   int i;

   assert(prob != NULL);
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(varused != NULL);
   assert(varinds != NULL);
   assert(nvarinds != NULL);
   assert(varsign != NULL);
   assert(boundtype != NULL);
   assert(freevariable != NULL);
   assert(localbdsused != NULL);

   *freevariable = FALSE;
   *localbdsused = FALSE;

#ifndef NDEBUG
   /* in debug mode, make sure that the whole array is initialized with invalid values */
   for( i = 0; i < prob->nvars; i++ )
   {
      varsign[i] = 0;
      boundtype[i] = -1;
   }
#endif

   /* allocate temporary memory to store best bounds and bound types */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bestlbs, prob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bestubs, prob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bestlbtypes, prob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bestubtypes, prob->nvars) );

   /* start with continuous variables, because using variable bounds can affect the untransformed integral
    * variables, and these changes have to be incorporated in the transformation of the integral variables
    * (continuous variables have largest problem indices!)
    */
   SCIPsortDownInt(varinds, *nvarinds);

   /* substitute continuous variables with best standard or variable bound (lb, ub, vlb or vub),
    * substitute integral variables with best standard bound (lb, ub)
    */
   i = 0;
   while( i < *nvarinds )
   {
      SCIP_VAR* var;
      SCIP_Real bestlb;
      SCIP_Real bestub;
      SCIP_Bool uselb;
      int bestlbtype;
      int bestubtype;
      int v;

      v = varinds[i];
      assert(0 <= v && v < prob->nvars);
      assert(varused[v]);

      var = prob->vars[v];
      assert(v == SCIPvarGetProbindex(var));

      /* due to variable bound usage cancellation may occur,
       * do not increase i, since last element is copied to the i-th position
       */
      if( SCIPsetIsZero(set, mircoef[v]) )
      {
         varsign[v] = +1;
         boundtype[v] = -1;
         mircoef[v] = 0.0;
         varused[v] = FALSE;
         varinds[i] = varinds[(*nvarinds)-1];
         (*nvarinds)--;
         continue;
      }

      /* check if the user specified a bound to be used */
      if( boundsfortrans != NULL && boundsfortrans[v] > -3 )
      {
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || ( boundsfortrans[v] == -2 || boundsfortrans[v] == -1 ));

         /* user has explicitly specified a bound to be used */
         if( boundtypesfortrans[v] == SCIP_BOUNDTYPE_LOWER )
         {
            /* user wants to use lower bound */
            bestlbtype = boundsfortrans[v];
            if( bestlbtype == -1 )
               bestlb = SCIPvarGetLbGlobal(var); /* use global standard lower bound */
            else if( bestlbtype == -2 )
               bestlb = SCIPvarGetLbLocal(var);  /* use local standard lower bound */
            else
            {
               SCIP_VAR** vlbvars;
               SCIP_Real* vlbcoefs;
               SCIP_Real* vlbconsts;
               int k;

               assert(!ignoresol);

               /* use the given variable lower bound */
               vlbvars = SCIPvarGetVlbVars(var);
               vlbcoefs = SCIPvarGetVlbCoefs(var);
               vlbconsts = SCIPvarGetVlbConstants(var);
               k = boundsfortrans[v];
               assert(k >= 0 && k < SCIPvarGetNVlbs(var));
               assert(vlbvars != NULL);
               assert(vlbcoefs != NULL);
               assert(vlbconsts != NULL);

               /* we have to avoid cyclic variable bound usage, so we enforce to use only variable bounds variables of smaller index */
               if( SCIPvarGetProbindex(vlbvars[k]) < v )
                  bestlb = vlbcoefs[k] * (sol == NULL ? SCIPvarGetLPSol(vlbvars[k]) : SCIPsolGetVal(sol, set, stat, vlbvars[k])) + vlbconsts[k];
               else
               {
                  bestlbtype = -1; /* fall back to global standard bound */
                  bestlb = SCIPvarGetLbGlobal(var);
               }
            }

            assert(!SCIPsetIsInfinity(set, -bestlb));
            uselb = TRUE;

            /* find closest upper bound in standard upper bound (and variable upper bounds for continuous variables) */
            findBestUb(set, stat, var, sol, usevbds && fixintegralrhs, allowlocal && fixintegralrhs, &bestub, &bestubtype);
         }
         else
         {
            assert(boundtypesfortrans[v] == SCIP_BOUNDTYPE_UPPER);

            /* user wants to use upper bound */
            bestubtype = boundsfortrans[v];
            if( bestubtype == -1 )
               bestub = SCIPvarGetUbGlobal(var); /* use global standard upper bound */
            else if( bestubtype == -2 )
               bestub = SCIPvarGetUbLocal(var);  /* use local standard upper bound */
            else
            {
               SCIP_VAR** vubvars;
               SCIP_Real* vubcoefs;
               SCIP_Real* vubconsts;
               int k;

               assert(!ignoresol);

               /* use the given variable upper bound */
               vubvars = SCIPvarGetVubVars(var);
               vubcoefs = SCIPvarGetVubCoefs(var);
               vubconsts = SCIPvarGetVubConstants(var);
               k = boundsfortrans[v];
               assert(k >= 0 && k < SCIPvarGetNVubs(var));
               assert(vubvars != NULL);
               assert(vubcoefs != NULL);
               assert(vubconsts != NULL);

               /* we have to avoid cyclic variable bound usage, so we enforce to use only variable bounds variables of smaller index */
               if( SCIPvarGetProbindex(vubvars[k]) < v )
                  bestub = vubcoefs[k] * (sol == NULL ? SCIPvarGetLPSol(vubvars[k]) : SCIPsolGetVal(sol, set, stat, vubvars[k])) + vubconsts[k];
               else
               {
                  bestubtype = -1; /* fall back to global standard bound */
                  bestub = SCIPvarGetUbGlobal(var);
               }
            }

            assert(!SCIPsetIsInfinity(set, bestub));
            uselb = FALSE;

            /* find closest lower bound in standard lower bound (and variable lower bounds for continuous variables) */
            findBestLb(set, stat, var, sol, usevbds && fixintegralrhs, allowlocal && fixintegralrhs, &bestlb, &bestlbtype);
         }
      }
      else
      {
         SCIP_Real varsol;

         /* bound selection should be done automatically */

         /* find closest lower bound in standard lower bound (and variable lower bounds for continuous variables) */
         findBestLb(set, stat, var, sol, usevbds, allowlocal, &bestlb, &bestlbtype);

         /* find closest upper bound in standard upper bound (and variable upper bounds for continuous variables) */
         findBestUb(set, stat, var, sol, usevbds, allowlocal, &bestub, &bestubtype);

         /* check, if variable is free variable */
         if( SCIPsetIsInfinity(set, -bestlb) && SCIPsetIsInfinity(set, bestub) )
         {
            /* we found a free variable in the row with non-zero coefficient
             *  -> MIR row can't be transformed in standard form
             */
            *freevariable = TRUE;
            goto TERMINATE;
         }

         if( !ignoresol )
         {
            /* select transformation bound */
            varsol = (sol == NULL ? SCIPvarGetLPSol(var) : SCIPsolGetVal(sol, set, stat, var));
            if( SCIPsetIsLT(set, varsol, (1.0 - boundswitch) * bestlb + boundswitch * bestub) )
               uselb = TRUE;
            else if( SCIPsetIsGT(set, varsol, (1.0 - boundswitch) * bestlb + boundswitch * bestub) )
               uselb = FALSE;
            else if( bestlbtype == -1 )  /* prefer global standard bounds */
               uselb = TRUE;
            else if( bestubtype == -1 )  /* prefer global standard bounds */
               uselb = FALSE;
            else if( bestlbtype >= 0 )   /* prefer variable bounds over local bounds */
               uselb = TRUE;
            else if( bestubtype >= 0 )   /* prefer variable bounds over local bounds */
               uselb = FALSE;
            else
               uselb = TRUE;             /* no decision yet? just use lower bound */
         }
         else
         {
            SCIP_Real glbub = SCIPvarGetUbGlobal(var);
            SCIP_Real glblb = SCIPvarGetLbGlobal(var);
            SCIP_Real distlb = REALABS(glblb - bestlb);
            SCIP_Real distub = REALABS(glbub - bestub);

            assert(!SCIPsetIsInfinity(set, -bestlb) || !SCIPsetIsInfinity(set, bestub));

            if( SCIPsetIsInfinity(set, -bestlb) )
               uselb = FALSE;
            else if( !SCIPsetIsNegative(set, bestlb) )
            {
               if( SCIPsetIsInfinity(set, bestub) )
                  uselb = TRUE;
               else if( SCIPsetIsZero(set, glblb) )
                  uselb = TRUE;
               else if( SCIPsetIsLE(set, distlb, distub) )
                  uselb = TRUE;
               else
                  uselb = FALSE;
            }
            else
            {
               assert(!SCIPsetIsInfinity(set, -bestlb));
               uselb = TRUE;
            }
         }
      }

      /* remember given/best bounds and types */
      bestlbs[v] = bestlb;
      bestubs[v] = bestub;
      bestlbtypes[v] = bestlbtype;
      bestubtypes[v] = bestubtype;

      /* perform bound substitution */
      if( uselb )
      {
         assert(!SCIPsetIsInfinity(set, -bestlb));

         /* use lower bound as transformation bound: x'_j := x_j - lb_j */
         boundtype[v] = bestlbtype;
         varsign[v] = +1;

         /* standard (bestlbtype < 0) or variable (bestlbtype >= 0) lower bound? */
         if( bestlbtype < 0 )
         {
            (*mirrhs) -= mircoef[v] * bestlb;
            *localbdsused = *localbdsused || (bestlbtype == -2);
         }
         else
         {
            SCIP_VAR** vlbvars;
            SCIP_Real* vlbcoefs;
            SCIP_Real* vlbconsts;
            int zidx;

            vlbvars = SCIPvarGetVlbVars(var);
            vlbcoefs = SCIPvarGetVlbCoefs(var);
            vlbconsts = SCIPvarGetVlbConstants(var);
            assert(vlbvars != NULL);
            assert(vlbcoefs != NULL);
            assert(vlbconsts != NULL);

            assert(0 <= bestlbtype && bestlbtype < SCIPvarGetNVlbs(var));
            assert(SCIPvarIsActive(vlbvars[bestlbtype]));
            zidx = SCIPvarGetProbindex(vlbvars[bestlbtype]);
            assert(0 <= zidx && zidx < v);

            (*mirrhs) -= mircoef[v] * vlbconsts[bestlbtype];
            mircoef[zidx] += mircoef[v] * vlbcoefs[bestlbtype];

            /* update sparsity pattern */
            if( !varused[zidx] )
            {
               assert(*nvarinds < prob->nvars);
               varused[zidx] = TRUE;
               varinds[*nvarinds] = zidx;
               (*nvarinds)++;
            }
         }
      }
      else
      {
         assert(!SCIPsetIsInfinity(set, bestub));

         /* use upper bound as transformation bound: x'_j := ub_j - x_j */
         boundtype[v] = bestubtype;
         varsign[v] = -1;

         /* standard (bestubtype < 0) or variable (bestubtype >= 0) upper bound? */
         if( bestubtype < 0 )
         {
            (*mirrhs) -= mircoef[v] * bestub;
            *localbdsused = *localbdsused || (bestubtype == -2);
         }
         else
         {
            SCIP_VAR** vubvars;
            SCIP_Real* vubcoefs;
            SCIP_Real* vubconsts;
            int zidx;

            vubvars = SCIPvarGetVubVars(var);
            vubcoefs = SCIPvarGetVubCoefs(var);
            vubconsts = SCIPvarGetVubConstants(var);
            assert(vubvars != NULL);
            assert(vubcoefs != NULL);
            assert(vubconsts != NULL);

            assert(0 <= bestubtype && bestubtype < SCIPvarGetNVubs(var));
            assert(SCIPvarIsActive(vubvars[bestubtype]));
            zidx = SCIPvarGetProbindex(vubvars[bestubtype]);
            assert(zidx >= 0);

            (*mirrhs) -= mircoef[v] * vubconsts[bestubtype];
            mircoef[zidx] += mircoef[v] * vubcoefs[bestubtype];

            /* update sparsity pattern */
            if( !varused[zidx] )
            {
               assert(*nvarinds < prob->nvars);
               varused[zidx] = TRUE;
               varinds[*nvarinds] = zidx;
               (*nvarinds)++;
            }
         }
      }
      ++i; /* increase iterator */

#ifdef SCIP_DEBUG
      if( bestlbtype >= 0 )
      {
         assert(bestlbtype < SCIPvarGetNVlbs(var));
         assert(SCIPvarGetVlbVars(var) != NULL);
         assert(SCIPvarGetVlbCoefs(var) != NULL);
         assert(SCIPvarGetVlbConstants(var) != NULL);
      }
      if( bestubtype >= 0 )
      {
         assert(bestubtype < SCIPvarGetNVubs(var));
         assert(SCIPvarGetVubVars(var) != NULL);
         assert(SCIPvarGetVubCoefs(var) != NULL);
         assert(SCIPvarGetVubConstants(var) != NULL);
      }

      if( !ignoresol )
      {
         SCIPdebugMessage("MIR var <%s>: varsign=%d, boundtype=%d, mircoef=%g, base=%d, sol=%g, lb=%g, ub=%g, vlb=%g<%s>%+g, vub=%g<%s>%+g -> rhs=%g\n",
            SCIPvarGetName(var), varsign[v], boundtype[v], mircoef[v],
            SCIPvarGetCol(var) != NULL ? SCIPcolGetBasisStatus(SCIPvarGetCol(var)) : SCIP_BASESTAT_ZERO,
            (sol == NULL ? SCIPvarGetLPSol(var) : SCIPsolGetVal(sol, set, stat, var)), bestlb, bestub,
            bestlbtype >= 0 ? SCIPvarGetVlbCoefs(var)[bestlbtype] : 0.0,
            bestlbtype >= 0 ? SCIPvarGetName(SCIPvarGetVlbVars(var)[bestlbtype]) : "-",
            bestlbtype >= 0 ? SCIPvarGetVlbConstants(var)[bestlbtype] : bestlb,
            bestubtype >= 0 ? SCIPvarGetVubCoefs(var)[bestubtype] : 0.0,
            bestubtype >= 0 ? SCIPvarGetName(SCIPvarGetVubVars(var)[bestubtype]) : "-",
            bestubtype >= 0 ? SCIPvarGetVubConstants(var)[bestubtype] : bestub,
            *mirrhs);
      }
      else
      {
         SCIPdebugMessage("MIR var <%s>: varsign=%d, boundtype=%d, mircoef=%g, base=%d, lb=%g, ub=%g, vlb=%g<%s>%+g, vub=%g<%s>%+g -> rhs=%g\n",
            SCIPvarGetName(var), varsign[v], boundtype[v], mircoef[v],
            SCIPvarGetCol(var) != NULL ? SCIPcolGetBasisStatus(SCIPvarGetCol(var)) : SCIP_BASESTAT_ZERO,
            bestlb, bestub,
            bestlbtype >= 0 ? SCIPvarGetVlbCoefs(var)[bestlbtype] : 0.0,
            bestlbtype >= 0 ? SCIPvarGetName(SCIPvarGetVlbVars(var)[bestlbtype]) : "-",
            bestlbtype >= 0 ? SCIPvarGetVlbConstants(var)[bestlbtype] : bestlb,
            bestubtype >= 0 ? SCIPvarGetVubCoefs(var)[bestubtype] : 0.0,
            bestubtype >= 0 ? SCIPvarGetName(SCIPvarGetVubVars(var)[bestubtype]) : "-",
            bestubtype >= 0 ? SCIPvarGetVubConstants(var)[bestubtype] : bestub,
            *mirrhs);
      }
#endif
   }

   if( fixintegralrhs )
   {
      SCIP_Real f0;

      /* check if rhs is fractional */
      f0 = SCIPsetSumFrac(set, *mirrhs);
      if( f0 < minfrac || f0 > maxfrac )
      {
         SCIP_Real bestviolgain;
         SCIP_Real bestnewf0;
         int bestv;

         /* choose complementation of one variable differently such that f0 is in correct range */
         bestv = -1;
         bestviolgain = -1e+100;
         bestnewf0 = 1.0;
         for( i = 0; i < *nvarinds; i++ )
         {
            int v;

            v = varinds[i];
            assert(0 <= v && v < prob->nvars);
            assert(varused[v]);
            assert(!SCIPsetIsZero(set, mircoef[v]));

            if( boundtype[v] < 0
               && ((varsign[v] == +1 && !SCIPsetIsInfinity(set, bestubs[v]) && bestubtypes[v] < 0)
                  || (varsign[v] == -1 && !SCIPsetIsInfinity(set, -bestlbs[v]) && bestlbtypes[v] < 0)) )
            {
               SCIP_Real fj;
               SCIP_Real newfj;
               SCIP_Real newrhs;
               SCIP_Real newf0;
               SCIP_Real solval;
               SCIP_Real viol;
               SCIP_Real newviol;
               SCIP_Real violgain;

               /* currently:              a'_j =  varsign * a_j  ->  f'_j =  a'_j - floor(a'_j)
                * after complementation: a''_j = -varsign * a_j  -> f''_j = a''_j - floor(a''_j) = 1 - f'_j
                *                        rhs'' = rhs' + varsign * a_j * (lb_j - ub_j)
                * cut violation from f0 and fj:   f'_0 -  f'_j *  x'_j
                * after complementation:         f''_0 - f''_j * x''_j
                *
                * for continuous variables, we just set f'_j = f''_j = |a'_j|
                */
               newrhs = *mirrhs + varsign[v] * mircoef[v] * (bestlbs[v] - bestubs[v]);
               newf0 = SCIPsetSumFrac(set, newrhs);
               if( newf0 < minfrac || newf0 > maxfrac )
                  continue;
               if( SCIPvarGetType(prob->vars[v]) == SCIP_VARTYPE_CONTINUOUS )
               {
                  fj = REALABS(mircoef[v]);
                  newfj = fj;
               }
               else
               {
                  fj = SCIPsetFrac(set, varsign[v] * mircoef[v]);
                  newfj = SCIPsetFrac(set, -varsign[v] * mircoef[v]);
               }

               if( !ignoresol )
               {
                  solval = (sol == NULL ? SCIPvarGetLPSol(prob->vars[v]) : SCIPsolGetVal(sol, set, stat, prob->vars[v]));
                  viol = f0 - fj * (varsign[v] == +1 ? solval - bestlbs[v] : bestubs[v] - solval);
                  newviol = newf0 - newfj * (varsign[v] == -1 ? solval - bestlbs[v] : bestubs[v] - solval);
                  violgain = newviol - viol;
               }
               else
               {
                  /* todo: this should be done, this can improve the dualray significantly */
                  SCIPerrorMessage("Cannot handle closest bounds with ignoring the LP solution.\n");
                  return SCIP_INVALIDCALL;
               }

               /* prefer larger violations; for equal violations, prefer smaller f0 values since then the possibility that
                * we f_j > f_0 is larger and we may improve some coefficients in rounding
                */
               if( SCIPsetIsGT(set, violgain, bestviolgain)
                  || (SCIPsetIsGE(set, violgain, bestviolgain) && newf0 < bestnewf0) )
               {
                  bestv = v;
                  bestviolgain = violgain;
                  bestnewf0 = newf0;
               }
            }
         }

         if( bestv >= 0 )
         {
            assert(bestv < prob->nvars);
            assert(boundtype[bestv] < 0);
            assert(!SCIPsetIsInfinity(set, -bestlbs[bestv]));
            assert(!SCIPsetIsInfinity(set, bestubs[bestv]));

            /* switch the complementation of this variable */
            (*mirrhs) += varsign[bestv] * mircoef[bestv] * (bestlbs[bestv] - bestubs[bestv]);
            if( varsign[bestv] == +1 )
            {
               /* switch to upper bound */
               assert(bestubtypes[bestv] < 0); /* cannot switch to a variable bound (would lead to further coef updates) */
               boundtype[bestv] = bestubtypes[bestv];
               varsign[bestv] = -1;
            }
            else
            {
               /* switch to lower bound */
               assert(bestlbtypes[bestv] < 0); /* cannot switch to a variable bound (would lead to further coef updates) */
               boundtype[bestv] = bestlbtypes[bestv];
               varsign[bestv] = +1;
            }
            *localbdsused = *localbdsused || (boundtype[bestv] == -2);
         }
      }
   }

 TERMINATE:

   /*free temporary memory */
   SCIPsetFreeBufferArray(set, &bestubtypes);
   SCIPsetFreeBufferArray(set, &bestlbtypes);
   SCIPsetFreeBufferArray(set, &bestubs);
   SCIPsetFreeBufferArray(set, &bestlbs);

   return SCIP_OKAY;
}

/** Calculate fractionalities \f$ f_0 := b - down(b) \f$, \f$ f_j := a^\prime_j - down(a^\prime_j) \f$ and
 *   integer \f$ k >= 1 \f$ with \f$ 1/(k + 1) <= f_0 < 1/k \f$ and \f$ (=> k = up(1/f_0) + 1) \f$
 *   integer \f$ 1 <= p_j <= k \f$ with \f$ f_0 + ((p_j - 1) * (1 - f_0)/k) < f_j <= f_0 + (p_j * (1 - f_0)/k)\f$ \f$ (=> p_j = up( k*(f_j - f_0)/(1 - f_0) )) \f$
 * and derive strong CG cut \f$ \tilde{a}*x^\prime <= down(b) \f$
 * \f[
 * \begin{array}{rll}
 * integers : &  \tilde{a}_j = down(a^\prime_j)                &, if \qquad f_j <= f_0 \\
 *            &  \tilde{a}_j = down(a^\prime_j) + p_j/(k + 1)  &, if \qquad f_j >  f_0 \\
 * continuous:&  \tilde{a}_j = 0                               &, if \qquad a^\prime_j >= 0 \\
 *            &  \mbox{no strong CG cut found}                 &, if \qquad a^\prime_j <  0
 * \end{array}
 * \f]
 *
 * Transform inequality back to \f$ \hat{a}*x <= rhs \f$:
 *
 *  (lb or ub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - lb_j,&   x_j == x^\prime_j + lb_j,&   a^\prime_j ==  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{if lb was used in transformation} \\
 *    x^\prime_j := ub_j - x_j,&   x_j == ub_j - x^\prime_j,&   a^\prime_j == -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{if ub was used in transformation}
 * \end{array}
 * \f]
 * \f[
 *  and move the constant terms
 * \begin{array}{rl}
 *    -\tilde{a}_j * lb_j == -\hat{a}_j * lb_j, & \mbox{or} \\
 *     \tilde{a}_j * ub_j == -\hat{a}_j * ub_j &
 * \end{array}
 * \f]
 *  to the rhs.
 *
 *  (vlb or vub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - (bl_j * zl_j + dl_j),&   x_j == x^\prime_j + (bl_j * zl_j + dl_j),&   a^\prime_j ==  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{(vlb)} \\
 *    x^\prime_j := (bu_j * zu_j + du_j) - x_j,&   x_j == (bu_j * zu_j + du_j) - x^\prime_j,&   a^\prime_j == -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{(vub)}
 * \end{array}
 * \f]
 *  move the constant terms
 * \f[
 * \begin{array}{rl}
 *    -\tilde{a}_j * dl_j == -\hat{a}_j * dl_j,& \mbox{or} \\
 *     \tilde{a}_j * du_j == -\hat{a}_j * du_j &
 * \end{array}
 * \f]
 *  to the rhs, and update the VB variable coefficients:
 * \f[
 * \begin{array}{ll}
 *    \hat{a}_{zl_j} := \hat{a}_{zl_j} - \tilde{a}_j * bl_j == \hat{a}_{zl_j} - \hat{a}_j * bl_j,& \mbox{or} \\
 *    \hat{a}_{zu_j} := \hat{a}_{zu_j} + \tilde{a}_j * bu_j == \hat{a}_{zu_j} - \hat{a}_j * bu_j &
 * \end{array}
 * \f]
 */
static
void cutsRoundStrongCGRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size nvars */
   SCIP_Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   SCIP_Bool*            varused,            /**< array to flag variables that appear in the MIR constraint */
   int*                  varinds,            /**< sparsity pattern of non-zero MIR coefficients */
   int*                  nvarinds,           /**< pointer to number of non-zero MIR coefficients */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable (vlb/vub_idx or -1 for lb/ub)*/
   SCIP_Real             f0,                 /**< fractional value of rhs */
   SCIP_Real             k                   /**< factor to strengthen strongcg cut */
   )
{
   SCIP_Real onedivoneminusf0;
   int nintvars;
   int i;

   assert(prob != NULL);
   assert(strongcgcoef != NULL);
   assert(strongcgrhs != NULL);
   assert(varused != NULL);
   assert(varinds != NULL);
   assert(nvarinds != NULL);
   assert(varsign != NULL);
   assert(0.0 < f0 && f0 < 1.0);

   onedivoneminusf0 = 1.0 / (1.0 - f0);
   nintvars = prob->nvars - prob->ncontvars;

   /* start with integer variables, since the variable bound substitutions can add up to the integer cut coefficients;
    * loop backwards to be able to delete coefficients from the sparsity pattern
    */
   SCIPsortDownInt(varinds, *nvarinds);

   /* integer variables */
   for( i = *nvarinds-1; i >= 0; i-- )
   {
      SCIP_VAR* var;
      SCIP_Real aj;
      SCIP_Real downaj;
      SCIP_Real cutaj;
      SCIP_Real fj;
      int v;

      v = varinds[i];
      assert(0 <= v && v < prob->nvars);
      assert(varused[v]);

      /* stop this loop if we reached the continuous variables */
      if( v >= nintvars )
      {
         assert(!SCIPvarIsIntegral(prob->vars[v]));
         assert(SCIPvarGetType(prob->vars[v]) == SCIP_VARTYPE_CONTINUOUS);
         break;
      }

      var = prob->vars[v];
      assert(var != NULL);
      assert(SCIPvarIsIntegral(var));
      assert(SCIPvarGetProbindex(var) == v);
      assert(boundtype[v] == -1 || boundtype[v] == -2);
      assert(varsign[v] == +1 || varsign[v] == -1);

      /* calculate the coefficient in the retransformed cut */
      aj = varsign[v] * strongcgcoef[v]; /* a'_j */
      downaj = SCIPsetFloor(set, aj);
      fj = aj - downaj;

      if( SCIPsetIsSumLE(set, fj, f0) )
         cutaj = varsign[v] * downaj; /* a^_j */
      else
      {
         SCIP_Real pj;

         pj = SCIPsetCeil(set, k * (fj - f0) * onedivoneminusf0);
         assert(pj >= 0); /* should be >= 1, but due to rounding bias can be 0 if fj almost equal to f0 */
         assert(pj <= k);
         cutaj = varsign[v] * (downaj + pj / (k + 1)); /* a^_j */
      }


      /* remove zero cut coefficients from sparsity pattern */
      if( SCIPsetIsZero(set, cutaj) )
      {
         strongcgcoef[v] = 0.0;
         varused[v] = FALSE;
         varinds[i] = varinds[(*nvarinds)-1];
         (*nvarinds)--;
         continue;
      }

      strongcgcoef[v] = cutaj;

      /* move the constant term  -a~_j * lb_j == -a^_j * lb_j , or  a~_j * ub_j == -a^_j * ub_j  to the rhs */
      if( varsign[v] == +1 )
      {
         /* lower bound was used */
         if( boundtype[v] == -1 )
         {
            assert(!SCIPsetIsInfinity(set, -SCIPvarGetLbGlobal(var)));
            (*strongcgrhs) += cutaj * SCIPvarGetLbGlobal(var);
         }
         else
         {
            assert(!SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)));
            (*strongcgrhs) += cutaj * SCIPvarGetLbLocal(var);
         }
      }
      else
      {
         /* upper bound was used */
         if( boundtype[v] == -1 )
         {
            assert(!SCIPsetIsInfinity(set, SCIPvarGetUbGlobal(var)));
            (*strongcgrhs) += cutaj * SCIPvarGetUbGlobal(var);
         }
         else
         {
            assert(!SCIPsetIsInfinity(set, SCIPvarGetUbLocal(var)));
            (*strongcgrhs) += cutaj * SCIPvarGetUbLocal(var);
         }
      }
   }

   /* continuous variables */
   for( ; i >= 0; i-- )
   {
      int v;

      v = varinds[i];
      assert(0 <= v && v < prob->nvars);
      assert(varused[v]);

#ifndef NDEBUG
      /* in a strong CG cut, cut coefficients of continuous variables are always zero; check this in debug mode */
      {
         SCIP_VAR* var;
         SCIP_Real aj;

         var = prob->vars[v];
         assert(var != NULL);
         assert(!SCIPvarIsIntegral(var));
         assert(SCIPvarGetProbindex(var) == v);
         assert(varsign[v] == +1 || varsign[v] == -1);

         /* calculate the coefficient in the retransformed cut */
         aj = varsign[v] * strongcgcoef[v]; /* a'_j */
         assert(aj >= 0.0);
      }
#endif

      strongcgcoef[v] = 0.0;

      /* remove zero cut coefficients from sparsity pattern */
      varused[v] = FALSE;
      varinds[i] = varinds[(*nvarinds)-1];
      (*nvarinds)--;
   }
}

/** Calculate fractionalities \f$ f_0 := b - down(b), f_j := a^\prime_j - down(a^\prime_j) \f$, and derive MIR cut \f$ \tilde{a} \cdot x' \leq down(b) \f$
 * \f[
 * \begin{array}{rll}
 *  integers :&  \tilde{a}_j = down(a^\prime_j),                        & if \qquad f_j \leq f_0 \\
 *            &  \tilde{a}_j = down(a^\prime_j) + (f_j - f_0)/(1 - f_0),& if \qquad f_j >  f_0 \\
 *  continuous:& \tilde{a}_j = 0,                                       & if \qquad a^\prime_j \geq 0 \\
 *             & \tilde{a}_j = a^\prime_j/(1 - f_0),                    & if \qquad a^\prime_j <  0
 * \end{array}
 * \f]
 *
 *  Transform inequality back to \f$ \hat{a} \cdot x \leq rhs \f$:
 *
 *  (lb or ub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - lb_j,&   x_j = x^\prime_j + lb_j,&   a^\prime_j =  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{if lb was used in transformation} \\
 *    x^\prime_j := ub_j - x_j,&   x_j = ub_j - x^\prime_j,&   a^\prime_j = -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{if ub was used in transformation}
 * \end{array}
 * \f]
 *  and move the constant terms
 * \f[
 * \begin{array}{cl}
 *    -\tilde{a}_j \cdot lb_j = -\hat{a}_j \cdot lb_j,& \mbox{or} \\
 *     \tilde{a}_j \cdot ub_j = -\hat{a}_j \cdot ub_j &
 * \end{array}
 * \f]
 *  to the rhs.
 *
 *  (vlb or vub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - (bl_j \cdot zl_j + dl_j),&   x_j = x^\prime_j + (bl_j\, zl_j + dl_j),&   a^\prime_j =  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{(vlb)} \\
 *    x^\prime_j := (bu_j\, zu_j + du_j) - x_j,&   x_j = (bu_j\, zu_j + du_j) - x^\prime_j,&   a^\prime_j = -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{(vub)}
 * \end{array}
 * \f]
 *  move the constant terms
 * \f[
 * \begin{array}{cl}
 *    -\tilde{a}_j\, dl_j = -\hat{a}_j\, dl_j,& \mbox{or} \\
 *     \tilde{a}_j\, du_j = -\hat{a}_j\, du_j &
 * \end{array}
 * \f]
 *  to the rhs, and update the VB variable coefficients:
 * \f[
 * \begin{array}{ll}
 *    \hat{a}_{zl_j} := \hat{a}_{zl_j} - \tilde{a}_j\, bl_j = \hat{a}_{zl_j} - \hat{a}_j\, bl_j,& \mbox{or} \\
 *    \hat{a}_{zu_j} := \hat{a}_{zu_j} + \tilde{a}_j\, bu_j = \hat{a}_{zu_j} - \hat{a}_j\, bu_j &
 * \end{array}
 * \f]
 */
static
void cutsRoundMIRRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   SCIP_Bool*            varused,            /**< array to flag variables that appear in the MIR constraint */
   int*                  varinds,            /**< sparsity pattern of non-zero MIR coefficients */
   int*                  nvarinds,           /**< pointer to number of non-zero MIR coefficients */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable (vlb/vub_idx or -1 for lb/ub) */
   SCIP_Real             f0                  /**< fractional value of rhs */
   )
{
   SCIP_Real onedivoneminusf0;
   int i;

   assert(prob != NULL);
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(varused != NULL);
   assert(varinds != NULL);
   assert(nvarinds != NULL);
   assert(varsign != NULL);
   assert(0.0 < f0 && f0 < 1.0);

   onedivoneminusf0 = 1.0 / (1.0 - f0);

   /* Loop backwards to be able to delete coefficients from the sparsity pattern. Additionally, the variable bound
    * substitutions are only used in such a way that a variable of higher index is substituted by a variable of a
    * lower index. Therefore, we must loop backwards.
    */
   SCIPsortDownInt(varinds, *nvarinds);

   for( i = *nvarinds-1; i >= 0; i-- )
   {
      SCIP_VAR* var;
      SCIP_Real cutaj;
      int v;

      v = varinds[i];
      assert(0 <= v && v < prob->nvars);
      assert(varused[v]);

      var = prob->vars[v];
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) == v);
      assert(varsign[v] == +1 || varsign[v] == -1);

      /* calculate the coefficient in the retransformed cut */
      if( SCIPvarIsIntegral(var) )
      {
         SCIP_Real aj;
         SCIP_Real downaj;
         SCIP_Real fj;

         aj = varsign[v] * mircoef[v]; /* a'_j */
         downaj = SCIPsetFloor(set, aj);
         fj = aj - downaj;

         if( SCIPsetIsSumLE(set, fj, f0) )
            cutaj = varsign[v] * downaj; /* a^_j */
         else
            cutaj = varsign[v] * (downaj + (fj - f0) * onedivoneminusf0); /* a^_j */
      }
      else
      {
         SCIP_Real aj;

         aj = varsign[v] * mircoef[v]; /* a'_j */
         if( aj >= 0.0 )
            cutaj = 0.0;
         else
            cutaj = varsign[v] * aj * onedivoneminusf0; /* a^_j */
      }

      /* remove zero cut coefficients from sparsity pattern */
      if( SCIPsetIsZero(set, cutaj) )
      {
         mircoef[v] = 0.0;
         varused[v] = FALSE;
         varinds[i] = varinds[(*nvarinds)-1];
         (*nvarinds)--;
         continue;
      }

      mircoef[v] = cutaj;

      /* check for variable bound use */
      if( boundtype[v] < 0 )
      {
         /* standard bound */

         /* move the constant term  -a~_j * lb_j == -a^_j * lb_j , or  a~_j * ub_j == -a^_j * ub_j  to the rhs */
         if( varsign[v] == +1 )
         {
            /* lower bound was used */
            if( boundtype[v] == -1 )
            {
               assert(!SCIPsetIsInfinity(set, -SCIPvarGetLbGlobal(var)));
               (*mirrhs) += cutaj * SCIPvarGetLbGlobal(var);
            }
            else
            {
               assert(!SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)));
               (*mirrhs) += cutaj * SCIPvarGetLbLocal(var);
            }
         }
         else
         {
            /* upper bound was used */
            if( boundtype[v] == -1 )
            {
               assert(!SCIPsetIsInfinity(set, SCIPvarGetUbGlobal(var)));
               (*mirrhs) += cutaj * SCIPvarGetUbGlobal(var);
            }
            else
            {
               assert(!SCIPsetIsInfinity(set, SCIPvarGetUbLocal(var)));
               (*mirrhs) += cutaj * SCIPvarGetUbLocal(var);
            }
         }
      }
      else
      {
         SCIP_VAR** vbz;
         SCIP_Real* vbb;
         SCIP_Real* vbd;
         int vbidx;
         int zidx;

         /* variable bound */
         vbidx = boundtype[v];

         /* change mirrhs and cutaj of integer variable z_j of variable bound */
         if( varsign[v] == +1 )
         {
            /* variable lower bound was used */
            assert(0 <= vbidx && vbidx < SCIPvarGetNVlbs(var));
            vbz = SCIPvarGetVlbVars(var);
            vbb = SCIPvarGetVlbCoefs(var);
            vbd = SCIPvarGetVlbConstants(var);
         }
         else
         {
            /* variable upper bound was used */
            assert(0 <= vbidx && vbidx < SCIPvarGetNVubs(var));
            vbz = SCIPvarGetVubVars(var);
            vbb = SCIPvarGetVubCoefs(var);
            vbd = SCIPvarGetVubConstants(var);
         }
         assert(SCIPvarIsActive(vbz[vbidx]));
         zidx = SCIPvarGetProbindex(vbz[vbidx]);
         assert(0 <= zidx && zidx < v);

         (*mirrhs) += cutaj * vbd[vbidx];
         mircoef[zidx] -= cutaj * vbb[vbidx];

         /* add variable to sparsity pattern */
         if( !varused[zidx] )
         {
            assert(*nvarinds < prob->nvars);
            varused[zidx] = TRUE;
            varinds[*nvarinds] = zidx;
            (*nvarinds)++;
         }
      }
   }
}

/** substitute aggregated slack variables:
 *
 *  The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
 *  variable only appears in its own row: \f$ a^\prime_r = scale * weight[r] * slacksign[r] \f$.
 *
 *  Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
 * \f[
 * \begin{array}{rll}
 *    integers:  & \hat{a}_r = \tilde{a}_r = down(a^\prime_r)                  &, if \qquad f_r <= f0 \\
 *               & \hat{a}_r = \tilde{a}_r = down(a^\prime_r) + p_r/(k + 1)    &, if \qquad f_r >  f0 \\
 *    continuous:& \hat{a}_r = \tilde{a}_r = 0                                 &, if \qquad a^\prime_r >= 0 \\
 *               & \mbox{no strong CG cut found}                               &, if \qquad a^\prime_r <  0
 * \end{array}
 * \f]
 *
 *  Substitute \f$ \hat{a}_r * s_r \f$ by adding \f$ \hat{a}_r \f$ times the slack's definition to the cut.
 */
static
void cutsSubstituteStrongCGRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size nvars */
   SCIP_Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   SCIP_Bool*            varused,            /**< array to flag variables that appear in the MIR constraint */
   int*                  varinds,            /**< sparsity pattern of non-zero MIR coefficients */
   int*                  nvarinds,           /**< pointer to number of non-zero MIR coefficients */
   int*                  rowinds,            /**< sparsity pattern of used rows */
   int                   nrowinds,           /**< number of used rows */
   SCIP_Real             f0,                 /**< fractional value of rhs */
   SCIP_Real             k                   /**< factor to strengthen strongcg cut */
   )
{  /*lint --e{715}*/
   SCIP_Real onedivoneminusf0;
   int i;

   assert(lp != NULL);
   assert(weights != NULL);
   assert(!SCIPsetIsZero(set, scale));
   assert(strongcgcoef != NULL);
   assert(strongcgrhs != NULL);
   assert(slacksign != NULL);
   assert(varused != NULL);
   assert(varinds != NULL);
   assert(nvarinds != NULL);
   assert(rowinds != NULL);
   assert(0.0 < f0 && f0 < 1.0);

   onedivoneminusf0 = 1.0 / (1.0 - f0);

   for( i = 0; i < nrowinds; i++ )
   {
      SCIP_ROW* row;
      SCIP_Real pr;
      SCIP_Real ar;
      SCIP_Real downar;
      SCIP_Real cutar;
      SCIP_Real fr;
      SCIP_Real mul;
      int idx;
      int r;
      int j;

      r = rowinds[i];
      assert(0 <= r && r < lp->nrows);
      assert(slacksign[r] == -1 || slacksign[r] == +1);
      assert(!SCIPsetIsZero(set, weights[r]));

      row = lp->rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* get the slack's coefficient a'_r in the aggregated row */
      ar = slacksign[r] * scale * weights[r];

      /* calculate slack variable's coefficient a^_r in the cut */
      if( row->integral )
      {
         /* slack variable is always integral: */
         downar = SCIPsetFloor(set, ar);
         fr = ar - downar;

         if( SCIPsetIsLE(set, fr, f0) )
            cutar = downar;
         else
         {
            pr = SCIPsetCeil(set, k * (fr - f0) * onedivoneminusf0);
            assert(pr >= 0); /* should be >= 1, but due to rounding bias can be 0 if fr almost equal to f0 */
            assert(pr <= k);
            cutar = downar + pr / (k + 1);
         }
      }
      else
      {
         /* slack variable is continuous: */
         assert(ar >= 0.0);
         continue; /* slack can be ignored, because its coefficient is reduced to 0.0 */
      }

      /* if the coefficient was reduced to zero, ignore the slack variable */
      if( SCIPsetIsZero(set, cutar) )
         continue;

      /* depending on the slack's sign, we have
       *   a*x + c + s == rhs  =>  s == - a*x - c + rhs,  or  a*x + c - s == lhs  =>  s == a*x + c - lhs
       * substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
       */
      mul = -slacksign[r] * cutar;

      /* add the slack's definition multiplied with a^_j to the cut */
      for( j = 0; j < row->len; ++j )
      {
         assert(row->cols[j] != NULL);
         assert(row->cols[j]->var != NULL);
         assert(SCIPvarGetStatus(row->cols[j]->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetCol(row->cols[j]->var) == row->cols[j]);
         assert(SCIPvarGetProbindex(row->cols[j]->var) == row->cols[j]->var_probindex);
         idx = row->cols[j]->var_probindex;
         strongcgcoef[idx] += mul * row->vals[j];

         /* update sparsity pattern */
         if( !varused[idx] )
         {
            varused[idx] = TRUE;
            varinds[*nvarinds] = idx;
            (*nvarinds)++;
         }
      }

      /* move slack's constant to the right hand side */
      if( slacksign[r] == +1 )
      {
         SCIP_Real rhs;

         /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a^_r * (rhs - c) to the right hand side */
         assert(!SCIPsetIsInfinity(set, row->rhs));
         rhs = row->rhs - row->constant;
         if( row->integral )
         {
            /* the right hand side was implicitly rounded down in row aggregation */
            rhs = SCIPsetFeasFloor(set, rhs);
         }
         *strongcgrhs -= cutar * rhs;
      }
      else
      {
         SCIP_Real lhs;

         /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a^_r * (c - lhs) to the right hand side */
         assert(!SCIPsetIsInfinity(set, -row->lhs));
         lhs = row->lhs - row->constant;
         if( row->integral )
         {
            /* the left hand side was implicitly rounded up in row aggregation */
            lhs = SCIPsetFeasCeil(set, lhs);
         }
         *strongcgrhs += cutar * lhs;
      }
   }

   /* set rhs to zero, if it's very close to */
   if( SCIPsetIsZero(set, *strongcgrhs) )
      *strongcgrhs = 0.0;
}

/** substitute aggregated slack variables:
 *
 *  The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
 *  variable only appears in its own row: \f$ a^\prime_r = scale * weight[r] * slacksign[r]. \f$
 *
 *  Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
 *  \f[
 *  \begin{array}{rll}
 *    integers : & \hat{a}_r = \tilde{a}_r = down(a^\prime_r),                      & \mbox{if}\qquad f_r <= f0 \\
 *               & \hat{a}_r = \tilde{a}_r = down(a^\prime_r) + (f_r - f0)/(1 - f0),& \mbox{if}\qquad f_r >  f0 \\
 *    continuous:& \hat{a}_r = \tilde{a}_r = 0,                                     & \mbox{if}\qquad a^\prime_r >= 0 \\
 *               & \hat{a}_r = \tilde{a}_r = a^\prime_r/(1 - f0),                   & \mbox{if}\qquad a^\prime_r <  0
 *  \end{array}
 *  \f]
 *
 *  Substitute \f$ \hat{a}_r \cdot s_r \f$ by adding \f$ \hat{a}_r \f$ times the slack's definition to the cut.
 */
static
void cutsSubstituteMIRRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   SCIP_Bool*            varused,            /**< array to flag variables that appear in the MIR constraint */
   int*                  varinds,            /**< sparsity pattern of non-zero MIR coefficients */
   int*                  nvarinds,           /**< pointer to number of non-zero MIR coefficients */
   int*                  rowinds,            /**< sparsity pattern of used rows */
   int                   nrowinds,           /**< number of used rows */
   SCIP_Real             f0                  /**< fractional value of rhs */
   )
{  /*lint --e{715}*/
   SCIP_Real onedivoneminusf0;
   int i;

   assert(lp != NULL);
   assert(weights != NULL);
   assert(!SCIPsetIsZero(set, scale));
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(slacksign != NULL);
   assert(varused != NULL);
   assert(varinds != NULL);
   assert(nvarinds != NULL);
   assert(rowinds != NULL);
   assert(0.0 < f0 && f0 < 1.0);

   onedivoneminusf0 = 1.0 / (1.0 - f0);
   for( i = 0; i < nrowinds; i++ )
   {
      SCIP_ROW* row;
      SCIP_Real ar;
      SCIP_Real downar;
      SCIP_Real cutar;
      SCIP_Real fr;
      SCIP_Real mul;
      int idx;
      int r;
      int j;

      r = rowinds[i];
      assert(0 <= r && r < lp->nrows);
      assert(slacksign[r] == -1 || slacksign[r] == +1);
      assert(!SCIPsetIsZero(set, weights[r]));

      row = lp->rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* get the slack's coefficient a'_r in the aggregated row */
      ar = slacksign[r] * scale * weights[r];

      /* calculate slack variable's coefficient a^_r in the cut */
      if( row->integral
         && ((slacksign[r] == +1 && SCIPsetIsFeasIntegral(set, row->rhs - row->constant))
            || (slacksign[r] == -1 && SCIPsetIsFeasIntegral(set, row->lhs - row->constant))) )
      {
         /* slack variable is always integral:
          *    a^_r = a~_r = down(a'_r)                      , if f_r <= f0
          *    a^_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
          */
         downar = SCIPsetFloor(set, ar);
         fr = ar - downar;
         if( SCIPsetIsLE(set, fr, f0) )
            cutar = downar;
         else
            cutar = downar + (fr - f0) * onedivoneminusf0;
      }
      else
      {
         /* slack variable is continuous:
          *    a^_r = a~_r = 0                               , if a'_r >= 0
          *    a^_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
          */
         if( ar >= 0.0 )
            continue; /* slack can be ignored, because its coefficient is reduced to 0.0 */
         else
            cutar = ar * onedivoneminusf0;
      }

      /* if the coefficient was reduced to zero, ignore the slack variable */
      if( SCIPsetIsZero(set, cutar) )
         continue;

      /* depending on the slack's sign, we have
       *   a*x + c + s == rhs  =>  s == - a*x - c + rhs,  or  a*x + c - s == lhs  =>  s == a*x + c - lhs
       * substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
       */
      mul = -slacksign[r] * cutar;

      /* add the slack's definition multiplied with a^_j to the cut */
      for( j = 0; j < row->len; ++j )
      {
         assert(row->cols[j] != NULL);
         assert(row->cols[j]->var != NULL);
         assert(SCIPvarGetStatus(row->cols[j]->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetCol(row->cols[j]->var) == row->cols[j]);
         assert(SCIPvarGetProbindex(row->cols[j]->var) == row->cols[j]->var_probindex);
         idx = row->cols[j]->var_probindex;
         mircoef[idx] += mul * row->vals[j];

         /* update sparsity pattern */
         if( !varused[idx] )
         {
            varused[idx] = TRUE;
            varinds[*nvarinds] = idx;
            (*nvarinds)++;
         }
      }

      /* move slack's constant to the right hand side */
      if( slacksign[r] == +1 )
      {
         SCIP_Real rhs;

         /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a^_r * (rhs - c) to the right hand side */
         assert(!SCIPsetIsInfinity(set, row->rhs));
         rhs = row->rhs - row->constant;
         if( row->integral )
         {
            /* the right hand side was implicitly rounded down in row aggregation */
            rhs = SCIPsetFeasFloor(set, rhs);
         }
         *mirrhs -= cutar * rhs;
      }
      else
      {
         SCIP_Real lhs;

         /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a^_r * (c - lhs) to the right hand side */
         assert(!SCIPsetIsInfinity(set, -row->lhs));
         lhs = row->lhs - row->constant;
         if( row->integral )
         {
            /* the left hand side was implicitly rounded up in row aggregation */
            lhs = SCIPsetFeasCeil(set, lhs);
         }
         *mirrhs += cutar * lhs;
      }
   }

   /* set rhs to zero, if it's very close to */
   if( SCIPsetIsZero(set, *mirrhs) )
      *mirrhs = 0.0;
}

/* calculates a strong CG cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because
 * these rows cannot participate in a strong CG cut.
 */
static
SCIP_RETCODE cutsLpCalcStrongCG(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce strong CG cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce strong CG cut for */
   SCIP_Real*            weights,            /**< row weights in row summation */
   int*                  rowinds,            /**< array to store indices of non-zero entries of the weights array, or
                                              *   NULL */
   int                   nrowinds,           /**< number of non-zero entries in weights array, -1 if rowinds is NULL */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size nvars */
   SCIP_Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   SCIP_Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid strong CG cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the returned cut is only valid locally */
   int*                  cutrank             /**< pointer to store the rank of the returned cut; or NULL */
   )
{
   int* slacksign;
   int* varsign;
   int* boundtype;
   SCIP_Bool* varused;
   int* varinds;
   int nvarinds;
   SCIP_Real rhs;
   SCIP_Real downrhs;
   SCIP_Real f0;
   SCIP_Real k;
   SCIP_Bool emptyrow;
   SCIP_Bool freevariable;
   SCIP_Bool localrowsused;
   SCIP_Bool localbdsused;
   SCIP_Bool rowtoolong;
   SCIP_Bool cleanuprowinds;

   assert(lp != NULL);
   assert(lp->solved);
   assert(prob != NULL);
   assert(weights != NULL);
   assert(!SCIPsetIsZero(set, scale));
   assert(strongcgcoef != NULL);
   assert(strongcgrhs != NULL);
   assert(cutactivity != NULL);
   assert(success != NULL);
   assert(cutislocal != NULL);
   assert(rowinds != NULL || nrowinds == -1);

   SCIPdebugMessage("calculating strong CG cut (scale: %g)\n", scale);

   /**@todo test, if a column based summation is faster */

   *success = FALSE;
   *cutislocal = FALSE;

   /* allocate temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &slacksign, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varsign, prob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &boundtype, prob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varused, prob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varinds, prob->nvars) );

   /* allocate memory for sparsity structure in case it has not been provided already */
   cleanuprowinds = FALSE;
   if( rowinds == NULL )
   {
      SCIP_CALL( SCIPsetAllocBufferArray(set, &rowinds, lp->nrows) );
      cleanuprowinds = TRUE;
   }

   /* calculate the row summation */
   cutsSumStrongCGRow(set, prob, lp, weights, scale, allowlocal,
      maxmksetcoefs, maxweightrange, strongcgcoef, &rhs, slacksign, varused, varinds, &nvarinds, rowinds, &nrowinds,
      &emptyrow, &localrowsused, &rowtoolong, cutrank);
   assert(allowlocal || !localrowsused);
   *cutislocal = *cutislocal || localrowsused;
   if( emptyrow || rowtoolong )
      goto TERMINATE;
   SCIPdebug(printMIR(set, stat, prob, NULL, strongcgcoef, rhs, FALSE, FALSE));

   /* remove all nearly-zero coefficients from strong CG row and relaxes the right hand side correspondingly in order to
    *  prevent numerical rounding errors
    */
   cutsCleanupMIRRow(set, prob, strongcgcoef, &rhs, varused, varinds, &nvarinds, *cutislocal);
   SCIPdebug(printMIR(set, stat, prob, NULL, strongcgcoef, rhs, FALSE, FALSE));

   /* Transform equation  a*x == b, lb <= x <= ub  into standard form
    *   a'*x' == b, 0 <= x' <= ub'.
    *
    * Transform variables (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
    * and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
    *
    * Transform variables (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
    * move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
    *   a_{zl_j} := a_{zl_j} + a_j * bl_j, or
    *   a_{zu_j} := a_{zu_j} + a_j * bu_j
    */
   cutsTransformStrongCGRow(set, stat, prob, boundswitch, usevbds, allowlocal,
      strongcgcoef, &rhs, varused, varinds, &nvarinds, varsign, boundtype, &freevariable, &localbdsused);
   assert(allowlocal || !localbdsused);
   *cutislocal = *cutislocal || localbdsused;
   if( freevariable )
      goto TERMINATE;
   SCIPdebug(printMIR(set, stat, prob, NULL, strongcgcoef, rhs, FALSE, FALSE));

   /* Calculate
    *  - fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j)
    *  - integer k >= 1 with 1/(k + 1) <= f_0 < 1/k
    *    (=> k = up(1/f_0) + 1)
    *  - integer 1 <= p_j <= k with f_0 + ((p_j - 1) * (1 - f_0)/k) < f_j <= f_0 + (p_j * (1 - f_0)/k)
    *    (=> p_j = up( (f_j - f_0)/((1 - f_0)/k) ))
    * and derive strong CG cut
    *   a~*x' <= (k+1) * down(b)
    * integers :  a~_j = down(a'_j)                , if f_j <= f_0
    *             a~_j = down(a'_j) + p_j/(k + 1)  , if f_j >  f_0
    * continuous: a~_j = 0                         , if a'_j >= 0
    *             no strong CG cut found          , if a'_j <  0
    *
    * Transform inequality back to a^*x <= rhs:
    *
    * (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a^_j :=  a~_j,   if lb was used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   if ub was used in transformation
    * and move the constant terms
    *   -a~_j * lb_j == -a^_j * lb_j, or
    *    a~_j * ub_j == -a^_j * ub_j
    * to the rhs.
    *
    * (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a^_j :=  a~_j,   (vlb)
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   (vub)
    * move the constant terms
    *   -a~_j * dl_j == -a^_j * dl_j, or
    *    a~_j * du_j == -a^_j * du_j
    * to the rhs, and update the VB variable coefficients:
    *   a^_{zl_j} := a^_{zl_j} - a~_j * bl_j == a^_{zl_j} - a^_j * bl_j, or
    *   a^_{zu_j} := a^_{zu_j} + a~_j * bu_j == a^_{zu_j} - a^_j * bu_j
    */
   downrhs = SCIPsetSumFloor(set, rhs);
   f0 = rhs - downrhs;
   if( f0 < minfrac || f0 > maxfrac )
      goto TERMINATE;
   k = SCIPsetCeil(set, 1.0 / f0) - 1;

   *strongcgrhs = downrhs;
   cutsRoundStrongCGRow(set, prob, strongcgcoef, strongcgrhs, varused, varinds, &nvarinds, varsign, boundtype, f0, k);
   SCIPdebug(printMIR(set, stat, prob, NULL, strongcgcoef, *strongcgrhs, FALSE, FALSE));

   /* substitute aggregated slack variables:
    *
    * The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
    * variable only appears in its own row:
    *    a'_r = scale * weight[r] * slacksign[r].
    *
    * Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
    *   integers :  a^_r = a~_r = (k + 1) * down(a'_r)        , if f_r <= f0
    *               a^_r = a~_r = (k + 1) * down(a'_r) + p_r  , if f_r >  f0
    *   continuous: a^_r = a~_r = 0                           , if a'_r >= 0
    *               a^_r = a~_r = a'_r/(1 - f0)               , if a'_r <  0
    *
    * Substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
    */
   cutsSubstituteStrongCGRow(set, stat, lp, weights, scale, strongcgcoef, strongcgrhs, slacksign,
      varused, varinds, &nvarinds, rowinds, nrowinds, f0, k);
   SCIPdebug(printMIR(set, stat, prob, NULL, strongcgcoef, *strongcgrhs, FALSE, FALSE));

   /* remove again all nearly-zero coefficients from strong CG row and relax the right hand side correspondingly in order to
    * prevent numerical rounding errors
    */
   cutsCleanupMIRRow(set, prob, strongcgcoef, &rhs, varused, varinds, &nvarinds, *cutislocal);
   SCIPdebug(printMIR(set, stat, prob, NULL, strongcgcoef, rhs, FALSE, FALSE));

   /* calculate cut activity */
   *cutactivity = getMIRRowActivity(set, stat, prob, NULL, strongcgcoef, varinds, nvarinds);
   *success = TRUE;

 TERMINATE:
   /* free temporary memory */
   if( cleanuprowinds )
   {
      SCIPsetFreeBufferArray(set, &rowinds);
   }

   SCIPsetFreeBufferArray(set, &varinds);
   SCIPsetFreeBufferArray(set, &varused);
   SCIPsetFreeBufferArray(set, &boundtype);
   SCIPsetFreeBufferArray(set, &varsign);
   SCIPsetFreeBufferArray(set, &slacksign);

   return SCIP_OKAY;
}

/* calculates an MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0 because these
 * rows cannot participate in an MIR cut.
 */
static
SCIP_RETCODE cutsLpCalcMIR(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_Real             maxweight,          /**< largest magnitude of weights; set to -1 if sparsity information is unknown */
   int*                  weightinds,         /**< sparsity pattern of weights; size nrowinds; NULL if sparsity info is unknown */
   int                   nweightinds,        /**< number of nonzeros in weights; -1 if rowinds is NULL */
   int                   rowlensum,          /**< total number of non-zeros in used rows (row associated with nonzero weight coefficient); -1 if unknown */
   int*                  sidetypes,          /**< specify row side type (-1 = lhs, 0 = unkown, 1 = rhs) or NULL for automatic choices */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mksetcoefs,         /**< array to store mixed knapsack set coefficients: size nvars; or NULL */
   SCIP_Bool*            mksetcoefsvalid,    /**< pointer to store whether mixed knapsack set coefficients are valid; or NULL */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   SCIP_Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid MIR cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the returned cut is only valid locally */
   int*                  cutrank             /**< pointer to store the rank of the returned cut; or NULL */
   )
{
   int* slacksign;
   int* varsign;
   int* boundtype;
   SCIP_Bool* varused;
   int* varinds;
   int* rowinds;
   int nvarinds;
   int nrowinds;
   SCIP_Real rhs;
   SCIP_Real downrhs;
   SCIP_Real f0;
   SCIP_Bool emptyrow;
   SCIP_Bool freevariable;
   SCIP_Bool localrowsused;
   SCIP_Bool localbdsused;
   SCIP_Bool rowtoolong;
   SCIP_Bool compress;

   assert(lp != NULL);
   assert(lp->solved || sol != NULL);
   assert(prob != NULL);
   assert(weights != NULL);
   assert(!SCIPsetIsZero(set, scale));
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(cutactivity != NULL);
   assert(success != NULL);
   assert(cutislocal != NULL);

   SCIPdebugMessage("calculating MIR cut (scale: %g)\n", scale);

   /**@todo test, if a column based summation is faster */

   *success = FALSE;
   if( mksetcoefsvalid != NULL )
      *mksetcoefsvalid = FALSE;

   /* allocate temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &slacksign, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varsign, prob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &boundtype, prob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varused, prob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varinds, prob->nvars) );

   /* if sparse information of weights is known, there is no need
    * to compute rowinds */
   if( weightinds == NULL )
   {
      SCIP_CALL( SCIPsetAllocBufferArray(set, &rowinds, lp->nrows) );
      compress = TRUE;
   }
   else
   {
      compress = FALSE;

      /* weightinds is the indices of the weights vector.
       * rowinds is the indices of the weights vector that is modified in the cutsSumMIRRow function.
       * duplication of weightinds is necessary to ensure weightinds is not modified. */
      SCIP_CALL( SCIPsetDuplicateBufferArray(set, &rowinds, weightinds, nweightinds) );
      nrowinds = nweightinds;

      if( rowlensum/5 > maxmksetcoefs )
      {
         *cutislocal = FALSE;
         goto TERMINATE;
      }
   }

   /* calculate the row summation */
   cutsSumMIRRow(set, prob, lp, weights, maxweight, sidetypes, scale, allowlocal,
      maxmksetcoefs, maxweightrange, compress, mircoef, &rhs, slacksign, varused, varinds, &nvarinds, rowinds, &nrowinds,
      &emptyrow, &localrowsused, &rowtoolong, cutrank);
   assert(allowlocal || !localrowsused);
   *cutislocal = localrowsused;
   if( emptyrow || rowtoolong )
      goto TERMINATE;
   SCIPdebug(printMIR(set, stat, prob, sol, mircoef, rhs, FALSE, FALSE));

   /* remove all nearly-zero coefficients from MIR row and relax the right hand side correspondingly in order to
    * prevent numerical rounding errors
    */
   cutsCleanupMIRRow(set, prob, mircoef, &rhs, varused, varinds, &nvarinds, *cutislocal);
   SCIPdebug(printMIR(set, stat, prob, sol, mircoef, rhs, FALSE, FALSE));

   /* Transform equation  a*x == b, lb <= x <= ub  into standard form
    *   a'*x' == b, 0 <= x' <= ub'.
    *
    * Transform variables (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
    * and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
    *
    * Transform variables (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
    * move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
    *   a_{zl_j} := a_{zl_j} + a_j * bl_j, or
    *   a_{zu_j} := a_{zu_j} + a_j * bu_j
    */
   SCIP_CALL( cutsTransformMIRRow(set, stat, prob, sol, boundswitch, usevbds, allowlocal, fixintegralrhs, FALSE,
         boundsfortrans, boundtypesfortrans, minfrac, maxfrac, mircoef, &rhs, varused, varinds, &nvarinds,
         varsign, boundtype, &freevariable, &localbdsused) );
   assert(allowlocal || !localbdsused);
   *cutislocal = *cutislocal || localbdsused;

   /* store the coefficients of the variables in the constructed mixed knapsack set */
   if( mksetcoefs != NULL )
      BMScopyMemoryArray(mksetcoefs, mircoef, prob->nvars);
   if( mksetcoefsvalid != NULL )
      *mksetcoefsvalid = TRUE;

   if( freevariable )
      goto TERMINATE;
   SCIPdebug(printMIR(set, stat, prob, sol, mircoef, rhs, FALSE, FALSE));

   /* Calculate fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j) , and derive MIR cut
    *   a~*x' <= down(b)
    * integers :  a~_j = down(a'_j)                      , if f_j <= f_0
    *             a~_j = down(a'_j) + (f_j - f0)/(1 - f0), if f_j >  f_0
    * continuous: a~_j = 0                               , if a'_j >= 0
    *             a~_j = a'_j/(1 - f0)                   , if a'_j <  0
    *
    * Transform inequality back to a^*x <= rhs:
    *
    * (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a^_j :=  a~_j,   if lb was used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   if ub was used in transformation
    * and move the constant terms
    *   -a~_j * lb_j == -a^_j * lb_j, or
    *    a~_j * ub_j == -a^_j * ub_j
    * to the rhs.
    *
    * (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a^_j :=  a~_j,   (vlb)
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   (vub)
    * move the constant terms
    *   -a~_j * dl_j == -a^_j * dl_j, or
    *    a~_j * du_j == -a^_j * du_j
    * to the rhs, and update the VB variable coefficients:
    *   a^_{zl_j} := a^_{zl_j} - a~_j * bl_j == a^_{zl_j} - a^_j * bl_j, or
    *   a^_{zu_j} := a^_{zu_j} + a~_j * bu_j == a^_{zu_j} - a^_j * bu_j
    */
   downrhs = SCIPsetSumFloor(set, rhs);
   f0 = rhs - downrhs;
   if( f0 < minfrac || f0 > maxfrac )
      goto TERMINATE;

   /* We multiply the coefficients of the base inequality roughly by scale/(1-f0).
    * If this gives a scalar that is very big, we better do not generate this cut.
    */
   if( REALABS(scale)/(1.0 - f0) > MAXCMIRSCALE )
      goto TERMINATE;

   *mirrhs = downrhs;
   cutsRoundMIRRow(set, prob, mircoef, mirrhs, varused, varinds, &nvarinds, varsign, boundtype, f0);
   SCIPdebug(printMIR(set, stat, prob, sol, mircoef, *mirrhs, FALSE, FALSE));

   /* substitute aggregated slack variables:
    *
    * The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
    * variable only appears in its own row:
    *    a'_r = scale * weight[r] * slacksign[r].
    *
    * Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
    *   integers :  a^_r = a~_r = down(a'_r)                      , if f_r <= f0
    *               a^_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
    *   continuous: a^_r = a~_r = 0                               , if a'_r >= 0
    *               a^_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
    *
    * Substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
    */
   cutsSubstituteMIRRow(set, stat, lp, weights, scale, mircoef, mirrhs, slacksign, varused, varinds, &nvarinds,
         rowinds, nrowinds, f0);
   SCIPdebug(printMIR(set, stat, prob, sol, mircoef, *mirrhs, FALSE, FALSE));

   /* remove again all nearly-zero coefficients from MIR row and relax the right hand side correspondingly in order to
    * prevent numerical rounding errors
    */
   cutsCleanupMIRRow(set, prob, mircoef, mirrhs, varused, varinds, &nvarinds, *cutislocal);
   SCIPdebug(printMIR(set, stat, prob, sol, mircoef, *mirrhs, FALSE, FALSE));

   /* calculate cut activity */
   *cutactivity = getMIRRowActivity(set, stat, prob, sol, mircoef, varinds, nvarinds);
   *success = TRUE;

 TERMINATE:
   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &rowinds);
   SCIPsetFreeBufferArray(set, &varinds);
   SCIPsetFreeBufferArray(set, &varused);
   SCIPsetFreeBufferArray(set, &boundtype);
   SCIPsetFreeBufferArray(set, &varsign);
   SCIPsetFreeBufferArray(set, &slacksign);

   /**@todo pass the sparsity pattern to the calling method in order to speed up the calling method's loops */

   return SCIP_OKAY;
}

/** applies the MIR function on a constraint; the constraint is given by pairs of variables and coefficients and a rhs */
static
SCIP_RETCODE cutsApplyMIR(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            prob,               /**< transformed problem */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: 0 vlb_idx/vub_idx,
                                              *   -1 for global lb/ub or -2 for local lb/ub
                                              */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mksetcoefs,         /**< array to store mixed knapsack set coefficients: size nvars; or NULL */
   SCIP_Bool*            mksetcoefsvalid,    /**< pointer to store whether mixed knapsack set coefficients are valid; or NULL */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size SCIPgetNVars() */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*                  varinds,            /**< array of variable indeces with a mircoef != 0 */
   int*                  nvarinds,           /**< number of variables indeces in varinds array */
   SCIP_Real*            minact,             /**< pointer to store the minimal activity */
   SCIP_Bool*            varused,            /**< array to store whether a variable has a mircoef != 0 */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid MIR cut */
   SCIP_Bool*            islocal             /**< pointer to store whether the returned constraint is only valid locally */
   )
{
   int* varsign;
   int* boundtype;
   SCIP_Real rhs;
   SCIP_Real downrhs;
   SCIP_Real f0;
   SCIP_Bool freevariable;
   SCIP_Bool localbdsused;
   int v;

   assert(prob != NULL);
   assert(SCIPsetIsPositive(set, scale));
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(success != NULL);
   assert(islocal != NULL);
   assert(*nvarinds >= 1);
   assert(varinds != NULL);
   assert(varused != NULL);

   /* only temporary */
   assert(mksetcoefs == NULL);
   assert(mksetcoefsvalid == NULL);

   SCIPdebugMessage("applying MIR function (scale: %g)\n", scale);

   *success = FALSE;

   /* allocate temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varsign, prob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &boundtype, prob->nvars) );

   if( !SCIPsetIsEQ(set, scale, 1.0) )
   {
      (*mirrhs) *= scale;
      for( v = 0; v < *nvarinds; v++ )
         mircoef[varinds[v]] *= scale;
   }

   /* remove all nearly-zero coefficients from MIR row and relax the right hand side correspondingly in order to
    * prevent numerical rounding errors
    */
   rhs = *mirrhs;
   cutsCleanupMIRRow(set, prob, mircoef, &rhs, varused, varinds, nvarinds, *islocal);
   SCIPdebug(printMIR(set, stat, prob, NULL, mircoef, *mirrhs, TRUE, *islocal));

   /* Transform equation  a*x == b, lb <= x <= ub  into standard form
    *   a'*x' == b, 0 <= x' <= ub'.
    *
    * Transform variables (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
    * and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
    *
    * Transform variables (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
    * move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
    *   a_{zl_j} := a_{zl_j} + a_j * bl_j, or
    *   a_{zu_j} := a_{zu_j} + a_j * bu_j
    */
   SCIP_CALL( cutsTransformMIRRow(set, stat, prob, NULL, boundswitch, usevbds, allowlocal, fixintegralrhs, TRUE,
         boundsfortrans, boundtypesfortrans, minfrac, maxfrac, mircoef, &rhs, varused, varinds,
         nvarinds, varsign, boundtype, &freevariable, &localbdsused) );
   assert(allowlocal || !localbdsused);
   *islocal = *islocal || localbdsused;

   if( freevariable )
      goto TERMINATE;
   SCIPdebug(printMIR(set, stat, prob, NULL, mircoef, *mirrhs, TRUE, *islocal));

   /* Calculate fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j) , and derive MIR cut
    *   a~*x' <= down(b)
    * integers :  a~_j = down(a'_j)                      , if f_j <= f_0
    *             a~_j = down(a'_j) + (f_j - f0)/(1 - f0), if f_j >  f_0
    * continuous: a~_j = 0                               , if a'_j >= 0
    *             a~_j = a'_j/(1 - f0)                   , if a'_j <  0
    *
    * Transform inequality back to a^*x <= rhs:
    *
    * (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a^_j :=  a~_j,   if lb was used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   if ub was used in transformation
    * and move the constant terms
    *   -a~_j * lb_j == -a^_j * lb_j, or
    *    a~_j * ub_j == -a^_j * ub_j
    * to the rhs.
    *
    * (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a^_j :=  a~_j,   (vlb)
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   (vub)
    * move the constant terms
    *   -a~_j * dl_j == -a^_j * dl_j, or
    *    a~_j * du_j == -a^_j * du_j
    * to the rhs, and update the VB variable coefficients:
    *   a^_{zl_j} := a^_{zl_j} - a~_j * bl_j == a^_{zl_j} - a^_j * bl_j, or
    *   a^_{zu_j} := a^_{zu_j} + a~_j * bu_j == a^_{zu_j} - a^_j * bu_j
    */
   downrhs = SCIPsetSumFloor(set, rhs);
   f0 = rhs - downrhs;
   if( f0 < minfrac || f0 > maxfrac )
      goto TERMINATE;

   /* We multiply the coefficients of the base inequality roughly by scale/(1-f0).
    * If this gives a scalar that is very big, we better do not generate this cut.
    */
   if( REALABS(scale)/(1.0 - f0) > MAXCMIRSCALE )
      goto TERMINATE;

   *mirrhs = downrhs;
   cutsRoundMIRRow(set, prob, mircoef, mirrhs, varused, varinds, nvarinds, varsign, boundtype, f0);
   SCIPdebug(printMIR(set, stat, prob, NULL, mircoef, *mirrhs, TRUE, *islocal));

   /* remove again all nearly-zero coefficients from MIR row and relax the right hand side correspondingly in order to
    * prevent numerical rounding errors
    */
   cutsCleanupMIRRow(set, prob, mircoef, mirrhs, varused, varinds, nvarinds, *islocal);
   SCIPdebug(printMIR(set, stat, prob, NULL, mircoef, rhs, TRUE, *islocal));

   /* calculate cut activity */
   *minact = getMIRMinActivity(prob, mircoef, varinds, *nvarinds, *islocal);
   *success = TRUE;

  TERMINATE:
   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &boundtype);
   SCIPsetFreeBufferArray(set, &varsign);

   return SCIP_OKAY;
}

/** calculates a strong CG cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0 because
 *  these rows cannot participate in an MIR cut.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcutsCalcStrongCG(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce strong CG cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce strong CG cut for */
   SCIP_Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   int*                  inds,               /**< indices of non-zero entries in weights array, or NULL */
   int                   ninds,              /**< number of indices of non-zero entries in weights array, -1 if inds is
                                              *   NULL */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mircoef,            /**< array to store strong CG coefficients: must be of size SCIPgetNVars() */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the strong CG row */
   SCIP_Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid strong CG cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the returned cut is only valid locally */
   int*                  cutrank             /**< pointer to store the rank of the returned cut; or NULL */
   )
{
   SCIP_CALL( cutsLpCalcStrongCG(scip->lp, scip->set, scip->stat, scip->transprob,
         boundswitch, usevbds, allowlocal, maxmksetcoefs, maxweightrange, minfrac, maxfrac, weights, inds, ninds, scale,
         mircoef, mirrhs, cutactivity, success, cutislocal, cutrank) );

   return SCIP_OKAY;
}

/** calculates an MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0 because
 *  these rows cannot participate in an MIR cut.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcutsCalcLpMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   SCIP_Real             maxweight,          /**< largest magnitude of weights; set to -1.0 if sparsity information is
                                              *   unknown */
   int*                  weightinds,         /**< sparsity pattern of weights; size nrowinds; NULL if sparsity info is
                                              *   unknown */
   int                   nweightinds,        /**< number of nonzeros in weights; -1 if rowinds is NULL */
   int                   rowlensum,          /**< total number of nonzeros in used rows (row associated with nonzero weight coefficient); -1 if unknown */
   int*                  sidetypes,          /**< specify row side type (-1 = lhs, 0 = unkown, 1 = rhs) or NULL for automatic choices */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mksetcoefs,         /**< array to store mixed knapsack set coefficients: size nvars; or NULL */
   SCIP_Bool*            mksetcoefsvalid,    /**< pointer to store whether mixed knapsack set coefficients are valid; or NULL */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size SCIPgetNVars() */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   SCIP_Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid MIR cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the returned cut is only valid locally */
   int*                  cutrank             /**< pointer to store the rank of the returned cut; or NULL */
   )
{
   SCIP_CALL( cutsLpCalcMIR(scip->lp, scip->set, scip->stat, scip->transprob, sol,
         boundswitch, usevbds, allowlocal, fixintegralrhs, boundsfortrans, boundtypesfortrans, maxmksetcoefs,
         maxweightrange, minfrac, maxfrac, weights, maxweight, weightinds, nweightinds, rowlensum, sidetypes,
         scale, mksetcoefs, mksetcoefsvalid, mircoef, mirrhs, cutactivity, success, cutislocal, cutrank) );

   return SCIP_OKAY;
}

/** applies the MIR function on a constraint; the constraint is given by pairs of variables and coefficients and a rhs.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcutsApplyMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: 0 vlb_idx/vub_idx,
                                              *   -1 for global lb/ub or -2 for local lb/ub
                                              */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mksetcoefs,         /**< array to store mixed knapsack set coefficients: size nvars; or NULL */
   SCIP_Bool*            mksetcoefsvalid,    /**< pointer to store whether mixed knapsack set coefficients are valid; or NULL */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size SCIPgetNVars() */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*                  varinds,            /**< array of variable indices with a mircoef != 0 */
   int*                  nvarinds,           /**< number of variables indices in varinds array */
   SCIP_Real*            minact,             /**< pointer to store the minimal activity */
   SCIP_Bool*            varused,            /**< array to store whether a variable has a mircoef != 0 */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid MIR cut */
   SCIP_Bool*            islocal             /**< pointer to store whether the returned constraint is only valid locally */
   )
{
   SCIP_CALL( cutsApplyMIR(scip->set, scip->stat, scip->transprob, boundswitch, usevbds, allowlocal, fixintegralrhs,
         boundsfortrans, boundtypesfortrans, minfrac, maxfrac, scale, mksetcoefs, mksetcoefsvalid, mircoef, mirrhs,
         varinds, nvarinds, minact, varused, success, islocal) );

   return SCIP_OKAY;
}

/** removes all nearly-zero coefficients from MIR row and relaxes the right hand side accordingly in order to prevent
 *  numerical rounding errors
 */
void SCIPcutsCleanupRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            coefs,              /**< array to store MIR coefficients: must be of size nvars */
   SCIP_Real*            rhs,                /**< pointer to store the right hand side of the MIR row */
   SCIP_Bool*            varused,            /**< array to flag variables that appear in the MIR constraint */
   int*                  varinds,            /**< sparsity pattern of non-zero MIR coefficients */
   int*                  nvarinds,           /**< pointer to number of non-zero MIR coefficients */
   SCIP_Bool             islocal             /**< is the row only valid locally? */
   )
{
   cutsCleanupMIRRow(scip->set, scip->transprob, coefs, rhs, varused, varinds, nvarinds, islocal);
}
