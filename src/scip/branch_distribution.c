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

/**@file   branch_distribution.c
 * @ingroup BRANCHINGRULES
 * @brief  distribution branching rule
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/branch_distribution.h"


#define BRANCHRULE_NAME            "distribution"
#define BRANCHRULE_DESC            "branching rule based on variable influence on cumulative normal distribution of row activities"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define SCOREPARAM_VALUES "hlvw"
#define DEFAULT_SCOREPARAM 'v'
#define DEFAULT_PRIORITY 2.0
#define SQRTOFTWO 1.4142136
#define SQUARED(x) (x) * (x)
#define DEFAULT_ONLYACTIVEROWS FALSE /**< should only rows which are active at the current node be considered? */
/*
 * Data structures
 */

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
   char                  scoreparam;     /**< parameter how the branch score is calculated */
   SCIP_Bool             onlyactiverows; /**< should only rows which are active at the current node be considered? */
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** calculates the cumulative distribution P(-infinity <= x <= value) that a normally distributed
 * random variable x takes a value between -infinity and parameter \p value.
 * The distribution is given by the respective mean and deviation. This implementation
 * uses the error function erf().
 */
static
SCIP_Real calcCumulativeDistribution(
   SCIP*                 scip,               /**< current SCIP */
   SCIP_Real             mean,               /**< the mean value of the distribution */
   SCIP_Real             variance,           /**< the square of the deviation of the distribution */
   SCIP_Real             value               /**< the upper limit of the calculated distribution integral */
   )
{
   SCIP_Real normvalue;
   SCIP_Real std;

   /* we need to calculate the standard deviation from the variance */
   assert(!SCIPisNegative(scip, variance));
   if( SCIPisFeasZero(scip, variance) )
      std = 0.0;
   else
      std = sqrt(variance);

   /* special treatment for zero variance */
   if( SCIPisFeasZero(scip, std) )
   {
      if( SCIPisFeasLE(scip, value, mean) )
         return 1.0;
      else
         return 0.0;
   }

   /* scale and translate to standard normal distribution. Factor sqrt(2) is needed for erf() function */
   normvalue = (value - mean)/(std * SQRTOFTWO);

   SCIPdebugMessage(" Normalized value %g = ( %g - %g ) / (%g * 1.4142136)\n", normvalue, value, mean, std);
   /* calculate the cumulative distribution function for normvalue. For negative normvalues, we negate
    * the normvalue and use the oddness of the erf()-function; special treatment for values close to zero.
    */
   if( SCIPisFeasZero(scip, normvalue) )
      return .5;
   else if( normvalue > 0 )
   {
      SCIP_Real erfresult;

      erfresult = erf(normvalue);
      return  erfresult / 2.0 + 0.5;
   }
   else
   {
      SCIP_Real erfresult;

      erfresult = erf(-normvalue);

      return 0.5 - erfresult / 2.0;
   }
}

/**
 * calculates the probability of satisfying an LP-row under the assumption
 * of uniformly distributed variable values. For inequalities, we use the cumulative distribution
 * function of the standard normal distribution PHI(rhs - mu/sqrt(sigma2)) to calculate the probability
 * for a right hand side row with mean activity mu and variance sigma2 to be satisfied.
 * Similarly, 1 - PHI(lhs - mu/sqrt(sigma2)) is the probability to satisfy a left hand side row.
 * For equations (lhs==rhs), we use the centeredness measure p = min(PHI(lhs'), 1-PHI(lhs'))/max(PHI(lhs'), 1 - PHI(lhs')),
 * where lhs' = lhs - mu / sqrt(sigma2)
 */
static
SCIP_Real rowCalcProbability(
   SCIP*                 scip,               /**< current scip */
   SCIP_ROW*             row,                /**< the row */
   SCIP_Real             mu,                 /**< the mean value of the row distribution */
   SCIP_Real             sigma2,             /**< the variance of the row distribution */
   int                   rowinfinitiesdown,  /**< the number of variables with infinite bounds to DECREASE activity */
   int                   rowinfinitiesup     /**< the number of variables with infinite bounds to INCREASE activity */
   )
{
   SCIP_Real rowprobability;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real lhsprob;
   SCIP_Real rhsprob;

   lhs = SCIProwGetLhs(row);
   rhs = SCIProwGetRhs(row);

   lhsprob = 1.0;
   rhsprob = 1.0;

   /* use the cumulative distribution if the row contains no variable to repair every infeasibility */
   if( !SCIPisInfinity(scip, rhs) && rowinfinitiesdown == 0 )
      rhsprob = calcCumulativeDistribution(scip, mu, sigma2, rhs);

   /* use the cumulative distribution if the row contains no variable to repair every infeasibility
    * otherwise the row can always be made feasible by increasing activity far enough
    */
   if( !SCIPisInfinity(scip, -lhs) && rowinfinitiesup == 0 )
      lhsprob = 1.0 - calcCumulativeDistribution(scip, mu, sigma2, lhs);

   /* use centeredness measure for equations; for inequalities, the minimum of the two probabilities is the
    * probability to satisfy the row */
   if( SCIPisFeasEQ(scip, lhs, rhs) )
   {
      SCIP_Real minprobability;
      SCIP_Real maxprobability;

      minprobability = MIN(rhsprob, lhsprob);
      maxprobability = MAX(lhsprob, rhsprob);
      rowprobability = minprobability / maxprobability;
   }
   else
      rowprobability = MIN(rhsprob, lhsprob);

   SCIPdebug( SCIPprintRow(scip, row, NULL) );
   SCIPdebugMessage(" Row %s, mean %g, sigma2 %g, LHS %g, RHS %g has probability %g to be satisfied\n",
      SCIProwGetName(row), mu, sigma2, lhs, rhs, rowprobability);

   assert(SCIPisFeasGE(scip, rowprobability, 0.0) && SCIPisFeasLE(scip, rowprobability, 1.0));

   return rowprobability;
}

/**
 * calculates the initial mean and variance of the row activity normal distribution.
 * The mean value m is given by m = \sum_i=1^n c_i * (lb_i +ub_i) / 2 where
 * n is the number of variables, and c_i, lb_i, ub_i are the variable coefficient and
 * bounds, respectively. With the same notation, the variance sigma2 is given by
 * sigma2 = \sum_i=1^n c_i^2 * (ub_i - lb_i)^2 / 12.
 */
static
void rowCalculateGauss(
   SCIP*                 scip,               /**< current scip */
   SCIP_ROW*             row,                /**< the row for which the gaussian normal distribution has to be calculated */
   SCIP_Real*            mu,                 /**< pointer to store the mean value of the gaussian normal distribution */
   SCIP_Real*            sigma2,             /**< pointer to store the variance value of the gaussian normal distribution */
   int*                  rowinfinitiesdown,  /**< pointer to store the number of variables with infinite bounds to DECREASE activity */
   int*                  rowinfinitiesup     /**< pointer to store the number of variables with infinite bounds to INCREASE activity */
   )
{
   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int nrowvals;
   int c;

   assert(scip != NULL);
   assert(row != NULL);
   assert(mu != NULL);
   assert(sigma2 != NULL);
   assert(rowinfinitiesup != NULL);
   assert(rowinfinitiesdown != NULL);

   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   nrowvals = SCIProwGetNNonz(row);

   assert(nrowvals == 0 || rowcols != NULL);
   assert(nrowvals == 0 || rowvals != NULL);

   *mu = SCIProwGetConstant(row);
   *sigma2 = 0.0;
   *rowinfinitiesdown = 0;
   *rowinfinitiesup = 0;

   /* loop over nonzero row coefficients and sum up the variable contributions to mu and sigma2 */
   for( c = 0; c < nrowvals; ++c )
   {
      SCIP_VAR* colvar;
      SCIP_Real colval;
      SCIP_Real colvarlb;
      SCIP_Real colvarub;
      SCIP_Real squarecoeff;
      SCIP_Real squarebounddiff;

      assert(rowcols[c] != NULL);
      colvar = SCIPcolGetVar(rowcols[c]);
      assert(colvar != NULL);

      colval = rowvals[c];
      colvarlb = SCIPvarGetLbLocal(colvar);
      colvarub = SCIPvarGetUbLocal(colvar);

      /* fixed variables can be skipped */
      if( SCIPisFeasEQ(scip, colvarlb, colvarub) )
      {
         SCIPdebugMessage("  Fixed Variable %g <= %s <= %g skipped\n", colvarlb, SCIPvarGetName(colvar), colvarub);
         continue;
      }

      /* variables with infinite bounds are skipped for the calculation of the variance; they need to
       * be accounted for by the counters for infinite row activity decrease and increase and they
       * are used to shift the row activity mean in case they have one nonzero, but finite bound */
      if( SCIPisInfinity(scip, -colvarlb) || SCIPisInfinity(scip, colvarub) )
      {
         if( SCIPisInfinity(scip, colvarub) )
         {
         /* an infinite upper bound gives the row an infinite maximum activity or minimum activity, if the coefficient is
          * positive or negative, resp.
          */
            if( SCIPisNegative(scip, colval) )
               ++(*rowinfinitiesdown);
            else
               ++(*rowinfinitiesup);
         }
         else
            /* finite bound contributes to row activity mean by adding offset */
            *mu += colvarub * colval;

         /* an infinite lower bound gives the row an infinite maximum activity or minimum activity, if the coefficient is
          * negative or positive, resp.
          */
         if( SCIPisInfinity(scip, -colvarlb) )
         {
            if( SCIPisPositive(scip, colval) )
               ++(*rowinfinitiesdown);
            else
               ++(*rowinfinitiesup);
         }
         else
            /* in the case of a finite lower bound, the bound contributes to row activity mean by adding offset */
            *mu += colvarlb * colval;

         /* skip the calculation of row mean and variance */
         SCIPdebugMessage("  Variable %g <= %s <= %g skipped due to infinite bounds\n", colvarlb, SCIPvarGetName(colvar), colvarub);
         continue;
      }

      /* actual values are updated; the contribution of the variable to mu is the arithmetic mean of its bounds */
      *mu += colval * (colvarlb + colvarub) / 2.0;

      /* the variance contribution of a variable is c^2 * (u - l)^2 / 12.0 */
      squarebounddiff = SQUARED(colvarub - colvarlb);
      squarecoeff = SQUARED(colval);
      *sigma2 += squarecoeff * squarebounddiff / 12.0;
   }

   SCIPdebug( SCIPprintRow(scip, row, NULL) );
   SCIPdebugMessage("  Row %s has a mean value of %g at a sigma2 of %g \n", SCIProwGetName(row), *mu, *sigma2);
}

/** update the up- and downscore of a single variable after calculating the impact of branching on a
 * particular row, depending on the chosen score parameter
 */
static
void getScore(
   SCIP*                 scip,               /**< current SCIP pointer */
   SCIP_Real             currentprob,        /**< the current probability */
   SCIP_Real             newprobup,          /**< the new probability if branched upwards */
   SCIP_Real             newprobdown,        /**< the new probability if branched downwards */
   SCIP_Real*            upscore,            /**< pointer to store the new score for branching up */
   SCIP_Real*            downscore,          /**< pointer to store the new score for branching down */
   char                  scoreparam          /**< parameter to determine the way the score is calculated */
   )
{
   assert(scip != NULL);
   assert(SCIPisFeasGE(scip, currentprob, 0.0) && SCIPisFeasLE(scip, currentprob, 1.0));
   assert(SCIPisFeasGE(scip, newprobup, 0.0) && SCIPisFeasLE(scip, newprobup, 1.0));
   assert(SCIPisFeasGE(scip, newprobdown, 0.0) && SCIPisFeasLE(scip, newprobdown, 1.0));
   assert(upscore != NULL);
   assert(downscore != NULL);

   /* update up and downscore depending on score parameter */
   switch( scoreparam )
   {
   case 'l' :
      /* 'l'owest cumulative probability */
      if( SCIPisGT(scip, -newprobup, *upscore) )
         *upscore = -newprobup;
      if( SCIPisGT(scip, -newprobdown, *downscore) )
         *downscore = -newprobdown;
      break;

   case 'h' :
      /* 'h'ighest cumulative probability */
      if( SCIPisGT(scip, newprobup, *upscore) )
         *upscore = newprobup;
      if( SCIPisGT(scip, newprobdown, *downscore) )
         *downscore = newprobdown;
      break;

   case 'v' :
      /* 'v'otes lowest cumulative probability */
      if( SCIPisLT(scip, newprobup, newprobdown) )
         *upscore += 1.0;
      else if( SCIPisGT(scip, newprobup, newprobdown) )
         *downscore += 1.0;
      break;

   case 'w' :
      /* votes highest cumulative probability */
      if( SCIPisGT(scip, newprobup, newprobdown) )
         *upscore += 1.0;
      else if( SCIPisLT(scip, newprobup, newprobdown) )
         *downscore += 1.0;
      break;

   default :
      SCIPwarningMessage(scip, " ERROR! No branching scheme selected ! Exiting  method\n");
      break;
   }
}

/* calculate the branching score of a variable, depending on the chosen score parameter */
static
SCIP_RETCODE calcBranchScore(
   SCIP*                 scip,               /**< current SCIP */
   SCIP_VAR*             var,                /**< candidate variable */
   SCIP_Real             lpsolval,           /**< current fractional LP-relaxation solution value  */
   SCIP_Real*            rowmeans,           /**< LP row gaussian mean acitivity values */
   SCIP_Real*            rowvariances,       /**< LP row gaussian variance activity values */
   int*                  rowinfinitiesdown,  /**< array with number of variables with infinite bounds which allow for
                                                * always repairing the constraint right hand side */
   int*                  rowinfinitiesup,    /**< array with number of variables with infinite bounds which allow for
                                                * always repairing the constraint left hand side */
   int                   nlprows,            /**< number of LP rows */
   SCIP_Real*            upscore,            /**< pointer to store the variable score when branching on it in upward direction */
   SCIP_Real*            downscore,          /**< pointer to store the variable score when branching on it in downward direction */
   char                  scoreparam,         /**< the score parameter of this branching rule */
   SCIP_Bool             onlyactiverows      /**< should only rows which are active at the current node be considered? */
   )
{
   SCIP_COL* varcol;
   SCIP_ROW** colrows;
   SCIP_Real* rowvals;
   SCIP_Real varlb;
   SCIP_Real varub;
   SCIP_Real squaredbounddiff; /* current squared difference of variable bounds (ub - lb)^2 */
   SCIP_Real newub;            /* new upper bound if branching downwards */
   SCIP_Real newlb;            /* new lower bound if branching upwards */
   SCIP_Real squaredbounddiffup; /* squared difference after branching upwards (ub - lb')^2 */
   SCIP_Real squaredbounddiffdown; /* squared difference after branching downwards (ub' - lb)^2 */
   SCIP_Real currentmean;      /* current mean value of variable uniform distribution */
   SCIP_Real meanup;           /* mean value of variable uniform distribution after branching up */
   SCIP_Real meandown;         /* mean value of variable uniform distribution after branching down*/
   int ncolrows;
   int i;

   assert(scip != NULL);
   assert(var != NULL);
   assert(rowmeans != NULL);
   assert(rowvariances != NULL);
   assert(nlprows > 0);
   assert(upscore != NULL);
   assert(downscore != NULL);
   assert(!SCIPisIntegral(scip, lpsolval));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   varcol = SCIPvarGetCol(var);
   assert(varcol != NULL);

   colrows = SCIPcolGetRows(varcol);
   rowvals = SCIPcolGetVals(varcol);
   ncolrows = SCIPcolGetNNonz(varcol);
   varlb = SCIPvarGetLbLocal(var);
   varub = SCIPvarGetUbLocal(var);
   assert(SCIPisFeasLT(scip, varlb, varub));

   /* calculate mean and variance of variable uniform distribution before and after branching */
   if( SCIPisInfinity(scip, varub) || SCIPisInfinity(scip, -varlb) )
   {
      /* variables with infinite bounds are not kept in the row activity variance */
      squaredbounddiff = 0.0;

      if( !SCIPisInfinity(scip, varub) )
         currentmean = varub;
      else if( !SCIPisInfinity(scip, varlb) )
         currentmean = varlb;
      else
         currentmean = 0.0;
   }
   else
   {
      squaredbounddiff = SQUARED(varub - varlb);
      squaredbounddiff /= 12.0;
      currentmean = (varub + varlb) * .5;
   }

   newlb = SCIPfeasCeil(scip, lpsolval);
   newub = SCIPfeasFloor(scip, lpsolval);

   /* calculate the variable's uniform distribution after branching up and down, respectively.
    * This is only possible for finite lower or upper bounds */
   if( !SCIPisInfinity(scip, varub) )
   {
      if( SCIPisFeasEQ(scip, newlb, varub) )
         meanup = newlb;
      else
      {
         squaredbounddiffup = SQUARED(varub - newlb);
         squaredbounddiffup /= 12.0;
         meanup = (newlb + varub) * .5;
      }
   }
   else
   {
      meanup = newlb;
      squaredbounddiffup = 0.0;
   }

   /* calculate the distribution mean and variance for a variable with finite lower bound */
   if( !SCIPisInfinity(scip, -varlb) )
   {
      if( SCIPisFeasEQ(scip, newub, varlb) )
         meandown = newub;
      else
      {
         squaredbounddiffdown = SQUARED(newub - varlb);
         squaredbounddiffdown /= 12.0;
         meandown = (newub + varlb) * .5;
      }
   }
   else
   {
      meandown = newub;
      squaredbounddiffdown = 0.0;
   }
   *upscore = 0.0;
   *downscore = 0.0;

   /* loop over the variable rows and calculate the up and down score */
   for( i = 0; i < ncolrows; ++i )
   {
      SCIP_ROW* row;
      SCIP_Real changedrowmean;
      SCIP_Real changedrowvariance;
      SCIP_Real currentrowprob;
      SCIP_Real newrowprobup;
      SCIP_Real newrowprobdown;
      SCIP_Real squaredcoeff;
      int rowpos;

      row = colrows[i];
      assert(row != NULL);

      rowpos = SCIProwGetLPPos(row);

      if( rowpos < 0 )
         break;

      /* skip non-active rows if the user parameter was set this way */
      if( onlyactiverows && SCIPisSumPositive(scip, SCIPgetRowLPFeasibility(scip, row)) )
         continue;

      /* calculate row activity distribution if this is the first candidate to appear in this row */
      if( rowmeans[rowpos] == SCIP_INVALID )
         rowCalculateGauss(scip, row, &rowmeans[rowpos], &rowvariances[rowpos],
               &rowinfinitiesdown[rowpos], &rowinfinitiesup[rowpos]);

      assert(rowpos < nlprows);
      currentrowprob = rowCalcProbability(scip, row, rowmeans[rowpos], rowvariances[rowpos],
            rowinfinitiesdown[rowpos], rowinfinitiesup[rowpos]);

      /* get variable's current expected contribution to row activity */
      squaredcoeff = SQUARED(rowvals[i]);

      /* first, get the probability change for the row if the variable is branched on upwards */
      if( !SCIPisInfinity(scip, varub) )
      {
         int rowinftiesdownafterbranch;
         int rowinftiesupafterbranch;

         changedrowmean = rowmeans[rowpos] - rowvals[i] * (currentmean - meanup);
         changedrowvariance = rowvariances[rowpos] - squaredcoeff * (squaredbounddiff - squaredbounddiffup);

         rowinftiesdownafterbranch = rowinfinitiesdown[rowpos];
         rowinftiesupafterbranch = rowinfinitiesup[rowpos];

         if( SCIPisInfinity(scip, -varlb) && SCIPisNegative(scip, rowvals[i]) )
         {
            assert(rowinftiesupafterbranch >= 1);
            rowinftiesupafterbranch -= 1;
         }
         if( SCIPisInfinity(scip, -varlb) && SCIPisPositive(scip, rowvals[i]) )
         {
            assert(rowinftiesdownafterbranch >= 1);
            rowinftiesdownafterbranch -= 1;
         }

         newrowprobup = rowCalcProbability(scip, row, changedrowmean, changedrowvariance, rowinftiesdownafterbranch,
               rowinftiesupafterbranch);
      }
      else
         newrowprobup = currentrowprob;
      /* do the same for the other branching direction */
      if( !SCIPisInfinity(scip, varlb) )
      {
         int rowinftiesdownafterbranch;
         int rowinftiesupafterbranch;

         changedrowmean = rowmeans[rowpos] - rowvals[i] * (currentmean - meandown);
         changedrowvariance = rowvariances[rowpos] - squaredcoeff * (squaredbounddiff - squaredbounddiffdown);

         rowinftiesdownafterbranch = rowinfinitiesdown[rowpos];
         rowinftiesupafterbranch = rowinfinitiesup[rowpos];

         if( SCIPisInfinity(scip, varub) && SCIPisPositive(scip, rowvals[i]) )
         {
            assert(rowinftiesupafterbranch >= 1);
            rowinftiesupafterbranch -= 1;
         }
         if( SCIPisInfinity(scip, varub) && SCIPisNegative(scip, rowvals[i]) )
         {
            assert(rowinftiesdownafterbranch >= 1);
            rowinftiesdownafterbranch -= 1;
         }

         newrowprobdown = rowCalcProbability(scip, row, changedrowmean, changedrowvariance, rowinftiesdownafterbranch,
               rowinftiesupafterbranch);
      }
      else
         newrowprobdown = currentrowprob;

      /* update the up and down score depending on the chosen scoring parameter */
      getScore(scip, currentrowprob, newrowprobup, newrowprobdown, upscore, downscore, scoreparam);

      SCIPdebugMessage("  Variable %s changes probability of row %s from %g to %g (branch up) or %g;\n",
         SCIPvarGetName(var), SCIProwGetName(row), currentrowprob, newrowprobup, newrowprobdown);
      SCIPdebugMessage("  -->  new variable score: %g (for branching up), %g (for branching down)\n",
         *upscore, *downscore);
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyDistribution)
{  /*lint --e{715}*/

   assert(scip != NULL);

   SCIP_CALL( SCIPincludeBranchruleDistribution(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeDistribution)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);
   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpDistribution)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** lpcands;
   SCIP_ROW** lprows;
   SCIP_VAR* bestcand;
   SCIP_NODE* downchild;
   SCIP_NODE* upchild;

   SCIP_Real* lpcandssol;
   SCIP_Real* rowmeans;
   SCIP_Real* rowvariances;
   int* rowinfinitiesdown;  /* count the number of variables with infinite bounds which allow for always repairing the constraint right hand side */
   int* rowinfinitiesup; /* count the number of variables with infinite bounds which allow for always repairing the constraint left hand side */

   SCIP_Real bestscore;
   SCIP_BRANCHDIR bestbranchdir;
   int nlpcands;
   int nlprows;
   int i;
   int c;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;


   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );

   if( nlpcands == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetLPRowsData(scip, &lprows, &nlprows) );


   /* allocate buffer arrays for row distribution information */
   /* todo speed this up by means of an event handler */
   SCIP_CALL( SCIPallocBufferArray(scip, &rowmeans, nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowvariances, nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowinfinitiesdown, nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowinfinitiesup, nlprows) );

   assert(lprows != NULL || nlprows == 0);

   bestscore = -1;
   bestbranchdir = SCIP_BRANCHDIR_AUTO;

   branchruledata = SCIPbranchruleGetData(branchrule);

   /* invalidate the current row distribution data which is initialized on the fly when looping over the candidates */
   for( i = 0; i < nlprows; ++i )
      rowmeans[i] = SCIP_INVALID;

   /* loop over candidate variables and calculate their score in changing the cumulative
    * probability of fulfilling each of their constraints */
   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_Real upscore;
      SCIP_Real downscore;

      SCIP_CALL( calcBranchScore(scip, lpcands[c], lpcandssol[c], rowmeans, rowvariances, rowinfinitiesdown,
            rowinfinitiesup, nlprows,
            &upscore, &downscore, branchruledata->scoreparam, branchruledata->onlyactiverows) );

      if( upscore > bestscore && upscore > downscore )
      {
         bestscore = upscore;
         bestbranchdir = SCIP_BRANCHDIR_UPWARDS;
         bestcand = lpcands[c];
      }
      else if( downscore > bestscore )
      {
         bestscore = downscore;
         bestbranchdir = SCIP_BRANCHDIR_DOWNWARDS;
         bestcand = lpcands[c];
      }

      assert(bestbranchdir == SCIP_BRANCHDIR_DOWNWARDS || bestbranchdir == SCIP_BRANCHDIR_UPWARDS);
      assert(bestcand != NULL);

      SCIPdebugMessage("  Candidate %s has score down %g and up %g \n", SCIPvarGetName(lpcands[c]), downscore, upscore);
      SCIPdebugMessage("  Best candidate: %s, score %g, direction %d\n", SCIPvarGetName(bestcand), bestscore, bestbranchdir);

   }

   assert(bestcand != NULL);

   SCIPdebugMessage("  Branching on variable %s with bounds [%g, %g] and solution value <%g>\n", SCIPvarGetName(bestcand),
      SCIPvarGetLbLocal(bestcand), SCIPvarGetUbLocal(bestcand), SCIPvarGetLPSol(bestcand));

   SCIPbranchVar(scip, bestcand, &downchild, NULL, &upchild);

   assert(downchild != NULL);
      assert(upchild != NULL);

   if( bestbranchdir == SCIP_BRANCHDIR_UPWARDS )
   {
      SCIPchgChildPrio(scip, upchild, DEFAULT_PRIORITY);
      SCIPdebugMessage("  Changing node priority of up-child\n");
   }
   else
   {
      assert(bestbranchdir == SCIP_BRANCHDIR_DOWNWARDS);
      SCIPchgChildPrio(scip, downchild, DEFAULT_PRIORITY);
      SCIPdebugMessage("  Changing node priority of down-child\n");
   }
   SCIPfreeBufferArray(scip, &rowinfinitiesup);
   SCIPfreeBufferArray(scip, &rowinfinitiesdown);
   SCIPfreeBufferArray(scip, &rowmeans);
   SCIPfreeBufferArray(scip, &rowvariances);

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the distribution branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleDistribution(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create distribution branching rule data */
   branchruledata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );

   branchrule = NULL;
   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
	 BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyDistribution) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeDistribution) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpDistribution) );

   /* add distribution branching rule parameters */
   SCIP_CALL( SCIPaddCharParam(scip, "branching/"BRANCHRULE_NAME"/scoreparam", "the score calculation; 'l'owest cumulative probability,'h'ighest c.p., 'v'otes lowest c.p., votes highest c.p.('w') ",
         &branchruledata->scoreparam, TRUE, DEFAULT_SCOREPARAM, SCOREPARAM_VALUES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/"BRANCHRULE_NAME"/onlyactiverows",
         "should only rows which are active at the current node be considered?",
         &branchruledata->onlyactiverows, TRUE, DEFAULT_ONLYACTIVEROWS, NULL, NULL) );

   return SCIP_OKAY;
}
