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

/**@file   heur_allinonediving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. the active constraints the variable appear in
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_allinonediving.h"
#include "scip/heuristics.h"
#include "scip/branch_distribution.h"
#include "scip/scipdefplugins.h"

#define HEUR_NAME             "allinonediving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings w.r.t. the active constraints"
#define HEUR_DISPCHAR         'a'
#define HEUR_PRIORITY         -70000
#define HEUR_FREQ             2
#define HEUR_FREQOFS          3
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define EVENT_DISTRIBUTION    SCIP_EVENTTYPE_BOUNDCHANGED /**< the event type to be handled by this event handler */
#define EVENTHDLR_NAME "eventhdlr_allinonedistribution"
#define SQUARED(x) ((x) * (x))

/**< bit masks that represent all supported dive types */
#define DIVETYPES_ACTCONS    SCIP_DIVETYPE_INTEGRALITY
#define DIVETYPES_DIST       SCIP_DIVETYPE_INTEGRALITY /**< bit mask that represents all supported dive types */
#define DIVETYPES_GUIDED     SCIP_DIVETYPE_INTEGRALITY /**< bit mask that represents all supported dive types */
#define DIVETYPES_COEF       SCIP_DIVETYPE_INTEGRALITY | SCIP_DIVETYPE_SOS1VARIABLE /**< bit mask that represents all supported dive types */
#define DIVETYPES_FRAC       SCIP_DIVETYPE_INTEGRALITY | SCIP_DIVETYPE_SOS1VARIABLE /**< bit mask that represents all supported dive types */
#define DIVETYPES_LINESEARCH SCIP_DIVETYPE_INTEGRALITY | SCIP_DIVETYPE_SOS1VARIABLE /**< bit mask that represents all supported dive types */
#define DIVETYPES_PSCOST     SCIP_DIVETYPE_INTEGRALITY /**< bit mask that represents all supported dive types */
#define DIVETYPES_VECLEN     SCIP_DIVETYPE_INTEGRALITY /**< bit mask that represents all supported dive types */


/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.05 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  1.0 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 1.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */
#define DEFAULT_LPRESOLVEDOMCHGQUOT 0.15
#define DEFAULT_LPSOLVEFREQ 0
#define DEFAULT_INITIALSEED 12345
#define DEFAULT_SELTYPE 'w'


/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */

   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler pointer */
   SCIP_VAR**            updatedvars;        /**< variables to process bound change events for */
   SCIP_Real*            rowmeans;           /**< row activity mean values for all rows */
   SCIP_Real*            rowvariances;       /**< row activity variances for all rows */
   SCIP_Real*            currentubs;         /**< variable upper bounds as currently saved in the row activities */
   SCIP_Real*            currentlbs;         /**< variable lower bounds as currently saved in the row activities */
   int*                  rowinfinitiesdown;  /**< count the number of variables with infinite bounds which allow for always
                                              *   repairing the constraint right hand side */
   int*                  rowinfinitiesup;    /**< count the number of variables with infinite bounds which allow for always
                                              *   repairing the constraint left hand side */
   int*                  varposs;            /**< array of variable positions in the updated variables array */
   int*                  varfilterposs;      /**< array of event filter positions for variable events */
   int                   nupdatedvars;       /**< the current number of variables with pending bound changes */
   int                   memsize;            /**< memory size of current arrays, needed for dynamic reallocation */
   int                   varpossmemsize;     /**< memory size of updated vars and varposs array */
   int                   lastselection;      /**< stores the last selected diveset when the heuristics was run */
   char                  scoreparam;         /**< score parameter for distribution branching */
   char                  scoretype;          /**< score parameter to compare different divesets */
   SCIP_Real             epsilon;            /**< parameter that increases probability of exploration among divesets */
   char                  seltype;            /**< selection strategy: (e)psilon-greedy, (w)eighted distribution, (n)ext diving */
   SCIP_Bool             usescipscore;       /**< should the SCIP branching score be used for weighing up and down score? */
};

enum AllInOneDiving_Method
{
   ALLINONEDIVING_ACTCONSDIVING = 0,
   ALLINONEDIVING_COEFDIVING = 1,
   ALLINONEDIVING_DISTRIBUTIONDIVING = 2,
   ALLINONEDIVING_FRACDIVING = 3,
   ALLINONEDIVING_GUIDEDDIVING = 4,
   ALLINONEDIVING_PSCOSTDIVING = 5,
   ALLINONEDIVING_LINESEARCHDIVING = 6,
   ALLINONEDIVING_VECLENDIVING = 7
};

typedef enum AllInOneDiving_Method ALLINONEDIVINGMETHOD;

/*
 * local methods
 */
struct SCIP_EventhdlrData
{
   SCIP_HEURDATA*  heurdata;     /**< the heuristic data to access distribution arrays */
};
/*
 * local methods
 */

/** ensure that maxindex + 1 rows can be represented in data arrays; memory gets reallocated with 10% extra space
 *  to save some time for future allocations */
static
SCIP_RETCODE heurdataEnsureArraySize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   maxindex            /**< row index at hand (size must be at least this large) */
   )
{
   int newsize;
   int r;

   /* maxindex fits in current array -> nothing to do */
   if( maxindex < heurdata->memsize )
      return SCIP_OKAY;

   /* new memory size is the max index + 1 plus 10% additional space */
   newsize = (int)SCIPfeasCeil(scip, (maxindex + 1) * 1.1);
   assert(newsize > heurdata->memsize);
   assert(heurdata->memsize >= 0);

   /* alloc memory arrays for row information */
   if( heurdata->memsize == 0 )
   {
      SCIP_VAR** vars;
      int v;
      int nvars;

      SCIPallocBufferArray(scip, &heurdata->rowinfinitiesdown, newsize);
      SCIPallocBufferArray(scip, &heurdata->rowinfinitiesup, newsize);
      SCIPallocBufferArray(scip, &heurdata->rowmeans, newsize);
      SCIPallocBufferArray(scip, &heurdata->rowvariances, newsize);

      assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      assert(nvars > 0);

      /* allocate variable update event processing array storage */
      SCIP_CALL( SCIPallocBufferArray(scip, &heurdata->varfilterposs, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &heurdata->varposs, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &heurdata->updatedvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &heurdata->currentubs, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &heurdata->currentlbs, nvars) );

      heurdata->varpossmemsize = nvars;
      heurdata->nupdatedvars = 0;

      /* init variable event processing data */
      for( v = 0; v < nvars; ++v )
      {
         assert(SCIPvarIsActive(vars[v]));
         assert(SCIPvarGetProbindex(vars[v]) == v);

         /* set up variable events to catch bound changes */
         SCIP_CALL( SCIPcatchVarEvent(scip, vars[v], EVENT_DISTRIBUTION, heurdata->eventhdlr, NULL, &(heurdata->varfilterposs[v])) );
         assert(heurdata->varfilterposs[v] >= 0);

         heurdata->varposs[v] = -1;
         heurdata->updatedvars[v] = NULL;
         heurdata->currentlbs[v] = SCIP_INVALID;
         heurdata->currentubs[v] = SCIP_INVALID;
      }

   }
   else
   {
      SCIPreallocBufferArray(scip, &heurdata->rowinfinitiesdown, newsize);
      SCIPreallocBufferArray(scip, &heurdata->rowinfinitiesup, newsize);
      SCIPreallocBufferArray(scip, &heurdata->rowmeans, newsize);
      SCIPreallocBufferArray(scip, &heurdata->rowvariances, newsize);
   }

   /* loop over extended arrays and invalidate data to trigger initialization of this row when necessary */
   for( r = heurdata->memsize; r < newsize; ++r )
   {
      heurdata->rowmeans[r] = SCIP_INVALID;
      heurdata->rowvariances[r] = SCIP_INVALID;
      heurdata->rowinfinitiesdown[r] = 0;
      heurdata->rowinfinitiesup[r] = 0;
   }

   /* adjust memsize */
   heurdata->memsize = newsize;

   return SCIP_OKAY;
}

/* update the variables current lower and upper bound */
static
void heurdataUpdateCurrentBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_VAR*             var                 /**< the variable to update current bounds */
   )
{
   int varindex;
   SCIP_Real lblocal;
   SCIP_Real ublocal;

   assert(var != NULL);

   varindex = SCIPvarGetProbindex(var);
   assert(0 <= varindex && varindex < heurdata->varpossmemsize);
   lblocal = SCIPvarGetLbLocal(var);
   ublocal = SCIPvarGetUbLocal(var);

   assert(SCIPisFeasLE(scip, lblocal, ublocal));

   heurdata->currentlbs[varindex] = lblocal;
   heurdata->currentubs[varindex] = ublocal;
}

/** calculates the initial mean and variance of the row activity normal distribution.
 *
 *  The mean value \f$ \mu \f$ is given by \f$ \mu = \sum_i=1^n c_i * (lb_i +ub_i) / 2 \f$ where
 *  \f$n \f$ is the number of variables, and \f$ c_i, lb_i, ub_i \f$ are the variable coefficient and
 *  bounds, respectively. With the same notation, the variance \f$ \sigma^2 \$ is given by
 *  \f$ \sigma^2 = \sum_i=1^n c_i^2 * \sigma^2_i \f$, with the variance being
 *  \f$ \sigma^2_i = ((ub_i - lb_i + 1)^2 - 1) / 12 \f$ for integer variables and
 *  \f$ \sigma^2_i = (ub_i - lb_i)^2 / 12 \f$ for continuous variables.
 */
static
void rowCalculateGauss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< the heuristic rule data */
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
      SCIP_Real varvariance;
      SCIP_Real varmean;
      int varindex;

      assert(rowcols[c] != NULL);
      colvar = SCIPcolGetVar(rowcols[c]);
      assert(colvar != NULL);

      colval = rowvals[c];
      colvarlb = SCIPvarGetLbLocal(colvar);
      colvarub = SCIPvarGetUbLocal(colvar);

      varmean = 0.0;
      varvariance = 0.0;
      varindex = SCIPvarGetProbindex(colvar);
      assert((heurdata->currentlbs[varindex] == SCIP_INVALID) == (heurdata->currentubs[varindex] == SCIP_INVALID));

      /* variable bounds need to be watched from now on */
      if( heurdata->currentlbs[varindex] == SCIP_INVALID )
         heurdataUpdateCurrentBounds(scip, heurdata, colvar);

      assert(!SCIPisInfinity(scip, colvarlb));
      assert(!SCIPisInfinity(scip, -colvarub));
      assert(SCIPisFeasLE(scip, colvarlb, colvarub));

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
            if( colval < 0.0 )
               ++(*rowinfinitiesdown);
            else
               ++(*rowinfinitiesup);
         }

         /* an infinite lower bound gives the row an infinite maximum activity or minimum activity, if the coefficient is
          * negative or positive, resp.
          */
         if( SCIPisInfinity(scip, -colvarlb) )
         {
            if( colval > 0.0 )
               ++(*rowinfinitiesdown);
            else
               ++(*rowinfinitiesup);
         }
      }
      SCIPvarCalcDistributionParameters(scip, colvarlb, colvarub, SCIPvarGetType(colvar), &varmean, &varvariance);

      /* actual values are updated; the contribution of the variable to mu is the arithmetic mean of its bounds */
      *mu += colval * varmean;

      /* the variance contribution of a variable is c^2 * (u - l)^2 / 12.0 for continuous and c^2 * ((u - l + 1)^2 - 1) / 12.0 for integer */
      squarecoeff = SQUARED(colval);
      *sigma2 += squarecoeff * varvariance;

      assert(!SCIPisFeasNegative(scip, *sigma2));
   }

   SCIPdebug( SCIPprintRow(scip, row, NULL) );
   SCIPdebugMessage("  Row %s has a mean value of %g at a sigma2 of %g \n", SCIProwGetName(row), *mu, *sigma2);
}

/** calculate the branching score of a variable, depending on the chosen score parameter */
static
SCIP_RETCODE calcBranchScore(
   SCIP*                 scip,               /**< current SCIP */
   SCIP_HEURDATA*        heurdata,           /**< branch rule data */
   SCIP_VAR*             var,                /**< candidate variable */
   SCIP_Real             lpsolval,           /**< current fractional LP-relaxation solution value  */
   SCIP_Real*            upscore,            /**< pointer to store the variable score when branching on it in upward direction */
   SCIP_Real*            downscore,          /**< pointer to store the variable score when branching on it in downward direction */
   char                  scoreparam          /**< the score parameter of this branching rule */
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
   SCIP_VARTYPE vartype;
   int ncolrows;
   int i;

   SCIP_Bool onlyactiverows; /* should only rows which are active at the current node be considered? */

   assert(scip != NULL);
   assert(var != NULL);
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
   vartype = SCIPvarGetType(var);

   /* calculate mean and variance of variable uniform distribution before and after branching */
   currentmean = 0.0;
   squaredbounddiff = 0.0;
   SCIPvarCalcDistributionParameters(scip, varlb, varub, vartype, &currentmean, &squaredbounddiff);

   newlb = SCIPfeasCeil(scip, lpsolval);
   newub = SCIPfeasFloor(scip, lpsolval);

   /* calculate the variable's uniform distribution after branching up and down, respectively. */
   squaredbounddiffup = 0.0;
   meanup = 0.0;
   SCIPvarCalcDistributionParameters(scip, newlb, varub, vartype, &meanup, &squaredbounddiffup);

   /* calculate the distribution mean and variance for a variable with finite lower bound */
   squaredbounddiffdown = 0.0;
   meandown = 0.0;
   SCIPvarCalcDistributionParameters(scip, varlb, newub, vartype, &meandown, &squaredbounddiffdown);

   /* initialize the variable's up and down score */
   *upscore = 0.0;
   *downscore = 0.0;

   onlyactiverows = FALSE;

   /* loop over the variable rows and calculate the up and down score */
   for( i = 0; i < ncolrows; ++i )
   {
      SCIP_ROW* row;
      SCIP_Real changedrowmean;
      SCIP_Real rowmean;
      SCIP_Real rowvariance;
      SCIP_Real changedrowvariance;
      SCIP_Real currentrowprob;
      SCIP_Real newrowprobup;
      SCIP_Real newrowprobdown;
      SCIP_Real squaredcoeff;
      SCIP_Real rowval;
      int rowinfinitiesdown;
      int rowinfinitiesup;
      int rowpos;

      row = colrows[i];
      rowval = rowvals[i];
      assert(row != NULL);

      /* we access the rows by their index */
      rowpos = SCIProwGetIndex(row);

      /* skip non-active rows if the user parameter was set this way */
      if( onlyactiverows && SCIPisSumPositive(scip, SCIPgetRowLPFeasibility(scip, row)) )
         continue;

      /* call method to ensure sufficient data capacity */
      SCIP_CALL( heurdataEnsureArraySize(scip, heurdata, rowpos) );

      /* calculate row activity distribution if this is the first candidate to appear in this row */
      if( heurdata->rowmeans[rowpos] == SCIP_INVALID )
      {
         rowCalculateGauss(scip, heurdata, row, &heurdata->rowmeans[rowpos], &heurdata->rowvariances[rowpos],
               &heurdata->rowinfinitiesdown[rowpos], &heurdata->rowinfinitiesup[rowpos]);
      }

      /* retrieve the row distribution parameters from the branch rule data */
      rowmean = heurdata->rowmeans[rowpos];
      rowvariance = heurdata->rowvariances[rowpos];
      rowinfinitiesdown = heurdata->rowinfinitiesdown[rowpos];
      rowinfinitiesup = heurdata->rowinfinitiesdown[rowpos];
      assert(!SCIPisNegative(scip, rowvariance));

      currentrowprob = SCIProwCalcProbability(scip, row, rowmean, rowvariance,
            rowinfinitiesdown, rowinfinitiesup);

      /* get variable's current expected contribution to row activity */
      squaredcoeff = SQUARED(rowval);

      /* first, get the probability change for the row if the variable is branched on upwards. The probability
       * can only be affected if the variable upper bound is finite
       */
      if( !SCIPisInfinity(scip, varub) )
      {
         int rowinftiesdownafterbranch;
         int rowinftiesupafterbranch;

         /* calculate how branching would affect the row parameters */
         changedrowmean = rowmean + rowval * (meanup - currentmean);
         changedrowvariance = rowvariance + squaredcoeff * (squaredbounddiffup - squaredbounddiff);
         changedrowvariance = MAX(0.0, changedrowvariance);

         rowinftiesdownafterbranch = rowinfinitiesdown;
         rowinftiesupafterbranch = rowinfinitiesup;

         /* account for changes of the row's infinite bound contributions */
         if( SCIPisInfinity(scip, -varlb) && rowval < 0.0 )
            rowinftiesupafterbranch--;
         if( SCIPisInfinity(scip, -varlb) && rowval > 0.0 )
            rowinftiesdownafterbranch--;

         assert(rowinftiesupafterbranch >= 0);
         assert(rowinftiesdownafterbranch >= 0);
         newrowprobup = SCIProwCalcProbability(scip, row, changedrowmean, changedrowvariance, rowinftiesdownafterbranch,
               rowinftiesupafterbranch);
      }
      else
         newrowprobup = currentrowprob;

      /* do the same for the other branching direction */
      if( !SCIPisInfinity(scip, varlb) )
      {
         int rowinftiesdownafterbranch;
         int rowinftiesupafterbranch;

         changedrowmean = rowmean + rowval * (meandown - currentmean);
         changedrowvariance = rowvariance + squaredcoeff * (squaredbounddiffdown - squaredbounddiff);
         changedrowvariance = MAX(0.0, changedrowvariance);

         rowinftiesdownafterbranch = rowinfinitiesdown;
         rowinftiesupafterbranch = rowinfinitiesup;

         /* account for changes of the row's infinite bound contributions */
         if( SCIPisInfinity(scip, varub) && rowval > 0.0 )
            rowinftiesupafterbranch -= 1;
         if( SCIPisInfinity(scip, varub) && rowval < 0.0 )
            rowinftiesdownafterbranch -= 1;

         assert(rowinftiesdownafterbranch >= 0);
         assert(rowinftiesupafterbranch >= 0);
         newrowprobdown = SCIProwCalcProbability(scip, row, changedrowmean, changedrowvariance, rowinftiesdownafterbranch,
               rowinftiesupafterbranch);
      }
      else
         newrowprobdown = currentrowprob;

      /* update the up and down score depending on the chosen scoring parameter */
      SCIP_CALL( SCIPupdateDistributionScore(scip, currentrowprob, newrowprobup, newrowprobdown, upscore, downscore, scoreparam) );

      SCIPdebugMessage("  Variable %s changes probability of row %s from %g to %g (branch up) or %g;\n",
         SCIPvarGetName(var), SCIProwGetName(row), currentrowprob, newrowprobup, newrowprobdown);
      SCIPdebugMessage("  -->  new variable score: %g (for branching up), %g (for branching down)\n",
         *upscore, *downscore);
   }

   return SCIP_OKAY;
}

/** free branchrule data */
static
SCIP_RETCODE heurdataFreeArrays(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< branching rule data */
   )
{
   assert(heurdata->memsize == 0 || heurdata->rowmeans != NULL);
   assert(heurdata->memsize >= 0);

   if( heurdata->memsize > 0 )
   {
      SCIPfreeBufferArray(scip, &heurdata->rowmeans);
      SCIPfreeBufferArray(scip, &heurdata->rowvariances);
      SCIPfreeBufferArray(scip, &heurdata->rowinfinitiesup);
      SCIPfreeBufferArray(scip, &heurdata->rowinfinitiesdown);

      heurdata->memsize = 0;
   }

   if( heurdata->varpossmemsize > 0 )
   {
      SCIP_VAR** vars;
      int v;

      assert(heurdata->varpossmemsize == SCIPgetNVars(scip));

      vars = SCIPgetVars(scip);
      for( v = heurdata->varpossmemsize - 1; v >= 0; --v )
      {
         SCIP_VAR* var;

         var = vars[v];

         assert(var != NULL);
         assert(v == SCIPvarGetProbindex(var));
         SCIP_CALL( SCIPdropVarEvent(scip, var, EVENT_DISTRIBUTION, heurdata->eventhdlr, NULL, heurdata->varfilterposs[v]) );
      }
      SCIPfreeBufferArray(scip, &heurdata->currentlbs);
      SCIPfreeBufferArray(scip, &heurdata->currentubs);
      SCIPfreeBufferArray(scip, &heurdata->updatedvars);
      SCIPfreeBufferArray(scip, &heurdata->varposs);
      SCIPfreeBufferArray(scip, &heurdata->varfilterposs);
   }
   /* allocate variable update event processing array storage */

   heurdata->varpossmemsize = 0;
   heurdata->nupdatedvars = 0;

   return SCIP_OKAY;
}

/** add variable to the bound change event queue; skipped if variable is already in there, or if variable has
 *  no row currently watched
 */
static
void heurdataAddBoundChangeVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_VAR*             var                 /**< the variable whose bound changes need to be processed */
   )
{
   int varindex;
   int varpos;

   assert(var != NULL);

   varindex = SCIPvarGetProbindex(var);
   assert(-1 <= varindex && varindex < heurdata->varpossmemsize);

   /* if variable is not active, it should not be watched */
   if( varindex == -1 )
      return;
   varpos = heurdata->varposs[varindex];
   assert(varpos < heurdata->nupdatedvars);

   /* nothing to do if variable is already in the queue */
   if( varpos >= 0 )
   {
      assert(heurdata->updatedvars[varpos] == var);

      return;
   }

   /* if none of the variables rows was calculated yet, variable needs not to be watched */
   assert((heurdata->currentlbs[varindex] == SCIP_INVALID) == (heurdata->currentubs[varindex] == SCIP_INVALID));
   if( heurdata->currentlbs[varindex] == SCIP_INVALID )
      return;

   /* add the variable to the branch rule data of variables to process updates for */
   assert(heurdata->varpossmemsize > heurdata->nupdatedvars);
   varpos = heurdata->nupdatedvars;
   heurdata->updatedvars[varpos] = var;
   heurdata->varposs[varindex] = varpos;
   ++heurdata->nupdatedvars;
}

/** returns the next unprocessed variable (last in, first out) with pending bound changes, or NULL */
static
SCIP_VAR* heurdataPopBoundChangeVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< branchrule data */
   )
{
   SCIP_VAR* var;
   int varpos;
   int varindex;

   assert(heurdata->nupdatedvars >= 0);

   /* return if no variable is currently pending */
   if( heurdata->nupdatedvars == 0 )
      return NULL;

   varpos = heurdata->nupdatedvars - 1;
   var = heurdata->updatedvars[varpos];
   assert(var != NULL);
   varindex = SCIPvarGetProbindex(var);
   assert(0 <= varindex && varindex < heurdata->varpossmemsize);
   assert(varpos == heurdata->varposs[varindex]);

   heurdata->varposs[varindex] = -1;
   heurdata->nupdatedvars--;

   return var;
}

/** process a variable from the queue of changed variables */
static
SCIP_RETCODE varProcessBoundChanges(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_VAR*             var                 /**< the variable whose bound changes need to be processed */
   )
{
   SCIP_ROW** colrows;
   SCIP_COL* varcol;
   SCIP_Real* colvals;
   SCIP_Real oldmean;
   SCIP_Real newmean;
   SCIP_Real oldvariance;
   SCIP_Real newvariance;
   SCIP_Real oldlb;
   SCIP_Real newlb;
   SCIP_Real oldub;
   SCIP_Real newub;
   SCIP_VARTYPE vartype;
   int ncolrows;
   int r;
   int varindex;

   /* ensure that this is a probing bound change */
   assert(SCIPinProbing(scip));

   assert(var != NULL);
   varcol = SCIPvarGetCol(var);
   assert(varcol != NULL);
   colrows = SCIPcolGetRows(varcol);
   colvals = SCIPcolGetVals(varcol);
   ncolrows = SCIPcolGetNNonz(varcol);

   varindex = SCIPvarGetProbindex(var);

   oldlb = heurdata->currentlbs[varindex];
   oldub = heurdata->currentubs[varindex];

   /* skip update if the variable has never been subject of previously calculated row activities */
   assert((oldlb == SCIP_INVALID) == (oldub == SCIP_INVALID));
   if( oldlb == SCIP_INVALID )
      return SCIP_OKAY;

   newlb = SCIPvarGetLbLocal(var);
   newub = SCIPvarGetUbLocal(var);

   /* skip update if the bound change events have cancelled out */
   if( SCIPisFeasEQ(scip, oldlb, newlb) && SCIPisFeasEQ(scip, oldub, newub) )
      return SCIP_OKAY;

   /* calculate old and new variable distribution mean and variance */
   oldvariance = 0.0;
   newvariance = 0.0;
   oldmean = 0.0;
   newmean = 0.0;
   vartype = SCIPvarGetType(var);
   SCIPvarCalcDistributionParameters(scip, oldlb, oldub, vartype, &oldmean, &oldvariance);
   SCIPvarCalcDistributionParameters(scip, newlb, newub, vartype, &newmean, &newvariance);

   /* loop over all rows of this variable and update activity distribution */
   for( r = 0; r < ncolrows; ++r )
   {
      int rowpos;

      assert(colrows[r] != NULL);
      rowpos = SCIProwGetIndex(colrows[r]);
      assert(rowpos >= 0);

      SCIP_CALL( heurdataEnsureArraySize(scip, heurdata, rowpos) );

      /* only consider rows for which activity distribution was already calculated */
      if( heurdata->rowmeans[rowpos] != SCIP_INVALID )
      {
         SCIP_Real coeff;
         SCIP_Real coeffsquared;
         assert(heurdata->rowvariances[rowpos] != SCIP_INVALID
               && SCIPisFeasGE(scip, heurdata->rowvariances[rowpos], 0.0));

         coeff = colvals[r];
         coeffsquared = SQUARED(coeff);

         /* update variable contribution to row activity distribution */
         heurdata->rowmeans[rowpos] += coeff * (newmean - oldmean);
         heurdata->rowvariances[rowpos] += coeffsquared * (newvariance - oldvariance);
         heurdata->rowvariances[rowpos] = MAX(0.0, heurdata->rowvariances[rowpos]);

         /* account for changes of the infinite contributions to row activities */
         if( coeff > 0.0 )
         {
            /* if the coefficient is positive, upper bounds affect activity up */
            if( SCIPisInfinity(scip, newub) && !SCIPisInfinity(scip, oldub) )
               ++heurdata->rowinfinitiesup[rowpos];
            else if( !SCIPisInfinity(scip, newub) && SCIPisInfinity(scip, oldub) )
               --heurdata->rowinfinitiesup[rowpos];

            if( SCIPisInfinity(scip, newlb) && !SCIPisInfinity(scip, oldlb) )
               ++heurdata->rowinfinitiesdown[rowpos];
            else if( !SCIPisInfinity(scip, newlb) && SCIPisInfinity(scip, oldlb) )
               --heurdata->rowinfinitiesdown[rowpos];
         }
         else if( coeff < 0.0 )
         {
            if( SCIPisInfinity(scip, newub) && !SCIPisInfinity(scip, oldub) )
               ++heurdata->rowinfinitiesdown[rowpos];
            else if( !SCIPisInfinity(scip, newub) && SCIPisInfinity(scip, oldub) )
               --heurdata->rowinfinitiesdown[rowpos];

            if( SCIPisInfinity(scip, newlb) && !SCIPisInfinity(scip, oldlb) )
               ++heurdata->rowinfinitiesup[rowpos];
            else if( !SCIPisInfinity(scip, newlb) && SCIPisInfinity(scip, oldlb) )
               --heurdata->rowinfinitiesup[rowpos];
         }
         assert(heurdata->rowinfinitiesdown[rowpos] >= 0);
         assert(heurdata->rowinfinitiesup[rowpos] >= 0);
      }
   }

   /* store the new local bounds in the data */
   heurdataUpdateCurrentBounds(scip, heurdata, var);

   return SCIP_OKAY;
}

/** todo get the score for this dive set */
static
SCIP_Real divesetGetScore(
   SCIP_DIVESET*         diveset,            /**< diving settings data structure */
   char                  scoretype           /**< score parameter */
   )
{
   switch (scoretype) {
      case 'n': /* min average nodes */
         return SCIPdivesetGetNProbingNodes(diveset) / (SCIPdivesetGetNCalls(diveset) + 10.0);

      case 'i': /* min avg LP iterations */
         return SCIPdivesetGetNLPIterations(diveset) / (SCIPdivesetGetNCalls(diveset) + 10.0);

      case 'c': /* min backtrack / conflict ratio */
         return SCIPdivesetGetNBacktracks(diveset) / (SCIPdivesetGetNConflicts(diveset) + 10.0);

      case 'd': /* minimum average depth (the current default) */
         return SCIPdivesetGetAvgDepth(diveset) * SCIPdivesetGetNCalls(diveset) / (SCIPdivesetGetNCalls(diveset) + 10.0);

      default:
         break;
   }
   return 0.0;
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeDistributiondiving)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);
   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyAllinonediving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurAllinonediving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeAllinonediving) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitAllinonediving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   heurdata->lastselection = -1;
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitAllinonediving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}

/*
 * heuristic specific interface methods
 */

/** scoring callback for distribution diving. best candidate maximizes the distribution score */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScoreDistributiondiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Real upscore;
   SCIP_Real downscore;
   int varindex;

   heurdata = SCIPheurGetData(SCIPdivesetGetHeur(diveset));
   assert(heurdata != NULL);

   /* process pending bound change events */
   while( heurdata->nupdatedvars > 0 )
   {
      SCIP_VAR* nextvar;

      /* pop the next variable from the queue and process its bound changes */
      nextvar = heurdataPopBoundChangeVar(scip, heurdata);
      assert(nextvar != NULL);
      SCIP_CALL( varProcessBoundChanges(scip, heurdata, nextvar) );
   }

   assert(cand != NULL);

   varindex = SCIPvarGetProbindex(cand);

   /* in debug mode, ensure that all bound process events which occurred in the mean time have been captured
    * by the branching rule event system
    */
   assert(SCIPisFeasLE(scip, SCIPvarGetLbLocal(cand), SCIPvarGetUbLocal(cand)));
   assert(0 <= varindex && varindex < heurdata->varpossmemsize);

   assert((heurdata->currentlbs[varindex] == SCIP_INVALID) == (heurdata->currentubs[varindex] == SCIP_INVALID));
   assert((heurdata->currentlbs[varindex] == SCIP_INVALID)
         || SCIPisFeasEQ(scip, SCIPvarGetLbLocal(cand), heurdata->currentlbs[varindex]));
   assert((heurdata->currentubs[varindex] == SCIP_INVALID)
         || SCIPisFeasEQ(scip, SCIPvarGetUbLocal(cand), heurdata->currentubs[varindex]));

   /* if the branching rule has not captured the variable bounds yet, this can be done now */
   if( heurdata->currentlbs[varindex] == SCIP_INVALID )
      heurdataUpdateCurrentBounds(scip, heurdata, cand);

   upscore = 0.0;
   downscore = 0.0;

   /* loop over candidate rows and determine the candidate up- and down- branching score w.r.t. the score parameter */
   SCIP_CALL( calcBranchScore(scip, heurdata, cand, candsol,
         &upscore, &downscore, heurdata->scoreparam) );

   *roundup = (upscore > downscore);

   /* if weighted scoring is enabled, use the branching score method of SCIP to weigh up and down score */
   if( heurdata->usescipscore )
      *score = SCIPgetBranchScore(scip, cand, downscore, upscore);
   else
      *score = MAX(upscore, downscore);

   return SCIP_OKAY;
}

/** event execution method of distribution branching which handles bound change events of variables */
static
SCIP_DECL_EVENTEXEC(eventExecDistribution)
{
   SCIP_HEURDATA* heurdata;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_VAR* var;

   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   heurdata = eventhdlrdata->heurdata;
   var = SCIPeventGetVar(event);

   /* add the variable to the queue of unprocessed variables; method itself ensures that every variable is added at most once */
   heurdataAddBoundChangeVar(scip, heurdata, var);

   return SCIP_OKAY;
}

/** get LP iteration limit for diving */
static
SCIP_Longint getLPIterlimit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< the heuristic */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   SCIP_Longint nsolsfound = SCIPheurGetNSolsFound(heur);
   SCIP_Longint nbestsolsfound = SCIPheurGetNBestSolsFound(heur);
   SCIP_Longint nlpiterations = SCIPgetNNodeLPIterations(scip);
   SCIP_Longint ncalls = SCIPheurGetNCalls(heur);

   SCIP_Longint nlpiterationsdive = 0;
   SCIP_Longint lpiterlimit;

   int i;

   /* loop over the divesets and collect their individual iterations */
   for( i = 0; i < SCIPheurGetNDivesets(heur); ++i )
   {
      nlpiterationsdive += SCIPdivesetGetNLPIterations(SCIPheurGetDivesets(heur)[i]);
   }

   /* author gregor
    *
    * TODO parameterize this sufficiently
    */

   lpiterlimit = (SCIP_Longint)(0.4 * (1.0 + 10*(nsolsfound+1.0)/(ncalls+1.0)) * nlpiterations);
   lpiterlimit += 8000;

   lpiterlimit -= nlpiterationsdive;

   return lpiterlimit;
}

/** sample from a distribution defined by weights */
static
int sampleWeighted(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      rand,               /**< random number generator */
   SCIP_Real*            weights,            /**< weights of a ground set that define the sampling distribution */
   int                   nweights            /**< number of elements in the ground set */
   )
{
   SCIP_Real weightsum;
   SCIP_Real randomnr;
   int w;
   int selection;

   weightsum = 0.0;
   /* collect sum of weights */
   for( w = 0; w < nweights; ++w )
   {
      weightsum += weights[w];
   }
   assert(weightsum > 0);

   randomnr = SCIPrandomGetReal(rand, 0, weightsum);

   weightsum = 0.0;
   /* choose first element i such that the weight sum exceeds the random number */
   for( w = 0; w < nweights - 1; ++w )
   {
      weightsum += weights[w];

      if( weightsum >= randomnr )
         break;
   }
   assert(w < nweights);
   assert(weights[w] > 0.0);

   return w;
}

/** select the diving method to apply */
static
SCIP_RETCODE selectDiving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< the heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int*                  selection           /**< selection made */
   )
{
   SCIP_Bool* methodunavailable;
   SCIP_DIVESET** divesets;
   int nlprows;
   int ndivesets;
   int d;
   SCIP_RANDNUMGEN* rng;
   SCIP_Real* weights;
   SCIP_Real epsilon_t;

   divesets = SCIPheurGetDivesets(heur);
   ndivesets = SCIPheurGetNDivesets(heur);
   assert(ndivesets > 0);
   assert(divesets != NULL);

   SCIP_CALL( SCIPallocClearBufferArray(scip, &methodunavailable, ndivesets) );

   /* check for certain methods if it is possible to execute them */
   nlprows = SCIPgetNLPRows(scip);
   if( nlprows == 0 )
      methodunavailable[(int)ALLINONEDIVING_DISTRIBUTIONDIVING] = TRUE;

   /* don't dive with guided diving if no feasible solutions exists or
    * if this solution lives in the original variable space,
    * because it might violate the global bounds of the current problem
    */
   if( SCIPgetNSols(scip) == 0 || SCIPsolIsOriginal(SCIPgetBestSol(scip)))
      methodunavailable[(int)ALLINONEDIVING_GUIDEDDIVING] = TRUE;

   rng = SCIPdivesetGetRandnumgen(divesets[0]);
   assert(rng != NULL);

   switch (heurdata->seltype) {
   case 'e':
      epsilon_t = heurdata->epsilon * sqrt(ndivesets / (SCIPheurGetNCalls(heur) + 1.0));
      epsilon_t = MAX(epsilon_t, 0.05);

      /* select one of the available methods at random */
      if( SCIPrandomGetReal(rng, 0.0, 1.0) < epsilon_t )
      {
         do
         {
            *selection = SCIPrandomGetInt(rng, 0, ndivesets - 1);
         }
         while( methodunavailable[*selection] );
      }
      else
      {
         SCIP_Real bestscore = SCIP_REAL_MAX;
         for( d = 0; d < SCIPheurGetNDivesets(heur); ++d )
         {
            SCIP_Real score;

            if( methodunavailable[d] )
               continue;

            score = divesetGetScore(divesets[d], heurdata->scoretype);
            if( !methodunavailable[d] && score < bestscore )
            {
               bestscore = score;
               *selection = d;
            }
         }
      }
      break;
   case 'w':

      SCIP_CALL( SCIPallocBufferArray(scip, &weights, ndivesets) );

      /* initialize weights as inverse of the score + a small positive epsilon */
      for( d = 0; d < ndivesets; ++d )
      {
         weights[d] = methodunavailable[d] ? 0.0 : 1 / (divesetGetScore(divesets[d], heurdata->scoretype) + 1e-4);
      }

      *selection = sampleWeighted(scip, rng, weights, ndivesets);

      SCIPfreeBufferArray(scip, &weights);
      break;
   case 'n':

         /* continue from last selection and stop at the next available method */
         *selection = heurdata->lastselection;

         do
         {
            *selection = (*selection + 1) % ndivesets;
         }
         while (methodunavailable[*selection]);
         heurdata->lastselection = *selection;
      break;
   default:
      break;
   }

   assert(*selection >= 0 && *selection < ndivesets);
   SCIPfreeBufferArray(scip, &methodunavailable);

   return SCIP_OKAY;


}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecAllinonediving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;
   SCIP_DIVESET** divesets;
   ALLINONEDIVINGMETHOD method;
   SCIP_Longint lpiterlimit;
   int selection;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   divesets = SCIPheurGetDivesets(heur);
   assert(SCIPheurGetNDivesets(heur) > 0);
   assert(divesets != NULL);
   heurdata = SCIPheurGetData(heur);

   *result = SCIP_DELAYED;

   /* do not call heuristic in node that was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* only call heuristic, if the LP solution is basic (which allows fast resolve in diving) */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   lpiterlimit = getLPIterlimit(scip, heur, heurdata);

   if( lpiterlimit <= 0 )
      return SCIP_OKAY;


   /* select the next diving strategy based on previous success */
   SCIP_CALL( selectDiving(scip, heur, heurdata, &selection) );
   assert(selection >= 0 && selection < SCIPheurGetNDivesets(heur));

   method = (ALLINONEDIVINGMETHOD)selection;

   if( method == ALLINONEDIVING_DISTRIBUTIONDIVING )
   {
      int nlprows = SCIPgetNLPRows(scip);

      /* this assertion is true because the selection process skips distribution diving otherwise */
      assert(nlprows > 0);
      SCIP_CALL( heurdataEnsureArraySize(scip, heurdata, nlprows) );
   }

   diveset = divesets[selection];
   assert(diveset != NULL);
   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible, lpiterlimit) );

   if( method == ALLINONEDIVING_DISTRIBUTIONDIVING )
   {
      SCIP_CALL( heurdataFreeArrays(scip, heurdata) );
   }

   return SCIP_OKAY;
}

/** creates the allinonediving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurAllinonediving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   SCIP_EVENTHDLRDATA* eventhdlrdata;

   /* create distributiondivingdata */
   heurdata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   heurdata->memsize = 0;
   heurdata->rowmeans = NULL;
   heurdata->rowvariances = NULL;
   heurdata->rowinfinitiesdown = NULL;
   heurdata->rowinfinitiesup = NULL;
   heurdata->varfilterposs = NULL;
   heurdata->currentlbs = NULL;
   heurdata->currentubs = NULL;

   heurdata->scoreparam = 'l';

   heurdata->usescipscore = TRUE;

   /* create event handler first to finish branch rule data */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   eventhdlrdata->heurdata = heurdata;

   heurdata->eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &heurdata->eventhdlr, EVENTHDLR_NAME,
         "event handler for dynamic acitivity distribution updating",
         eventExecDistribution, eventhdlrdata) );
   assert( heurdata->eventhdlr != NULL);
   SCIP_CALL( SCIPsetEventhdlrFree(scip, heurdata->eventhdlr, eventFreeDistributiondiving) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecAllinonediving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyAllinonediving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeAllinonediving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitAllinonediving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitAllinonediving) );

   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, "allinone_actcons", DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL,
         DEFAULT_LPRESOLVEDOMCHGQUOT, DEFAULT_LPSOLVEFREQ,
         DEFAULT_MAXLPITEROFS, DEFAULT_INITIALSEED + 2, DEFAULT_BACKTRACK, TRUE, DIVETYPES_ACTCONS,
         SCIPdivesetGetScoreActconsdiving) );

   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, "allinone_coef", DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL,
         DEFAULT_LPRESOLVEDOMCHGQUOT, DEFAULT_LPSOLVEFREQ,
         DEFAULT_MAXLPITEROFS, DEFAULT_INITIALSEED + 3, DEFAULT_BACKTRACK, TRUE, DIVETYPES_COEF,
         SCIPdivesetGetScoreCoefdiving) );
   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, "allinone_dist", DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL,
         DEFAULT_LPRESOLVEDOMCHGQUOT, DEFAULT_LPSOLVEFREQ,
         DEFAULT_MAXLPITEROFS, DEFAULT_INITIALSEED, DEFAULT_BACKTRACK, TRUE, DIVETYPES_DIST,
         divesetGetScoreDistributiondiving) );

   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, "allinone_frac", DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL,
         DEFAULT_LPRESOLVEDOMCHGQUOT, DEFAULT_LPSOLVEFREQ,
         DEFAULT_MAXLPITEROFS, DEFAULT_INITIALSEED + 5, DEFAULT_BACKTRACK, TRUE, DIVETYPES_FRAC,
         SCIPdivesetGetScoreFracdiving) );

   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, "allinone_guided", DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL,
         DEFAULT_LPRESOLVEDOMCHGQUOT, DEFAULT_LPSOLVEFREQ,
         DEFAULT_MAXLPITEROFS, DEFAULT_INITIALSEED + 7, DEFAULT_BACKTRACK, TRUE, DIVETYPES_GUIDED,
         SCIPdivesetGetScoreGuideddiving) );
   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, "allinone_pscost", DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL,
         DEFAULT_LPRESOLVEDOMCHGQUOT, DEFAULT_LPSOLVEFREQ,
         DEFAULT_MAXLPITEROFS, DEFAULT_INITIALSEED + 11, DEFAULT_BACKTRACK, TRUE, DIVETYPES_PSCOST,
         SCIPdivesetGetScorePscostdiving) );
   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, "allinone_line", DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL,
         DEFAULT_LPRESOLVEDOMCHGQUOT, DEFAULT_LPSOLVEFREQ,
         DEFAULT_MAXLPITEROFS, DEFAULT_INITIALSEED + 13, DEFAULT_BACKTRACK, TRUE, DIVETYPES_LINESEARCH,
         SCIPdivesetGetScoreLinesearchdiving) );
   /* create a diveset (this will automatically install some additional parameters for the heuristic) */
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, "allinone_veclen", DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL,
         DEFAULT_LPRESOLVEDOMCHGQUOT, DEFAULT_LPSOLVEFREQ,
         DEFAULT_MAXLPITEROFS, DEFAULT_INITIALSEED + 17, DEFAULT_BACKTRACK, TRUE, DIVETYPES_VECLEN,
         SCIPdivesetGetScoreVeclendiving) );

   /* author gregor
    *
    * TODO put default values to the top of the file as preprocessor defines
    */
   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/seltype",
         "selection strategy: (e)psilon-greedy, (w)eighted distribution, (n)ext diving",
         &heurdata->seltype, FALSE, DEFAULT_SELTYPE, "enw", NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/scoretype",
         "score parameter", &heurdata->scoretype, FALSE, 'd', "nicd", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/epsilon",
         "parameter that increases probability of exploration among divesets",
         &heurdata->epsilon, FALSE, 1.0, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
