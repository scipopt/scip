/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cutsel_dynamic.c
 * @ingroup DEFPLUGINS_CUTSEL
 * @brief  dynamic cut selector
 * @author Christoph Graczyk
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>


#include "scip/scip_cutsel.h"
#include "scip/scip_cut.h"
#include "scip/scip_lp.h"
#include "scip/scip_randnumgen.h"
#include "scip/cutsel_dynamic.h"


#define CUTSEL_NAME               "dynamic"
#define CUTSEL_DESC               "dynamic orthogonality for hybrid cutsel"
#define CUTSEL_PRIORITY           7000

#define RANDSEED                  0x5EED

#define DEFAULT_EFFICACYWEIGHT          1.0  /**< weight of efficacy in score calculation */
#define DEFAULT_DIRCUTOFFDISTWEIGHT     0.0  /**< weight of directed cutoff distance in score calculation */
#define DEFAULT_OBJPARALWEIGHT          0.0  /**< weight of objective parallelism in score calculation */
#define DEFAULT_INTSUPPORTWEIGHT        0.0  /**< weight of integral support in cut score calculation */
#define DEFAULT_MINORTHO                0.9  /**< minimal orthogonality in percent for a cut to enter the LP */
#define DEFAULT_MINGAIN                 0.01 /**< minimal efficacy gain for a cut to enter the LP */
#define DEFAULT_MAXDEPTH                (-1) /**< maximum depth at which this cutselector is used (-1 : all nodes) */
#define DEFAULT_FILTERMODE              'd'  /**< filtering strategy during cut selection (
                                               *  'd'ynamic-  and 'f'ull dynamic parallelism) */


/*
 * Data structures
 */

/** cut selector data */
struct SCIP_CutselData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random generator for tiebreaking */
   SCIP_Real             objparalweight;     /**< weight of objective parallelism in cut score calculation */
   SCIP_Real             efficacyweight;     /**< weight of efficacy in cut score calculation */
   SCIP_Real             dircutoffdistweight;/**< weight of directed cutoff distance in cut score calculation */
   SCIP_Real             intsupportweight;   /**< weight of integral support in cut score calculation */
   SCIP_Real             mingain;            /**< minimal projection efficacy gain for a cut to enter the LP in percent */
   SCIP_Real             minortho;           /**< minimal orthogonality for a cut to enter the LP */
   int                   maxdepth;           /**< maximum depth at which this cutselector is used (-1 : all nodes) */
   char                  filtermode;         /**< filtering strategy during cut selection (
                                              *  'd'ynamic- and 'f'ull dynamic parallelism) */
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** returns the maximum score of cuts; if scores is not NULL, then stores the individual score of each cut in scores */
static
void scoring(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            cuts,               /**< array with cuts to score */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator for tie-breaking, or NULL */
   SCIP_Real             dircutoffdistweight,/**< weight of directed cutoff distance in cut score calculation */
   SCIP_Real             efficacyweight,     /**< weight of efficacy in cut score calculation */
   SCIP_Real             objparalweight,     /**< weight of objective parallelism in cut score calculation */
   SCIP_Real             intsupportweight,   /**< weight of integral support in cut score calculation */
   int*                  currentncuts,       /**< current number of cuts in cuts array */
   SCIP_Real*            scores              /**< array to store the score of cuts or NULL */
   )
{
   SCIP_Real maxscore = 0.0;
   SCIP_SOL* sol;
   int i;
   int ncuts = *currentncuts;

   sol = SCIPgetBestSol(scip);

   /* if there is an incumbent and the factor is not 0.0, compute directed cutoff distances for the incumbent */
   if( sol != NULL && dircutoffdistweight > 0.0 )
   {
      for( i = ncuts-1; i >= 0; --i )
      {
         SCIP_Real score;
         SCIP_Real objparallelism;
         SCIP_Real intsupport;
         SCIP_Real efficacy;

         if( intsupportweight > 0.0 )
            intsupport = intsupportweight * SCIPgetRowNumIntCols(scip, cuts[i]) / (SCIP_Real) SCIProwGetNNonz(cuts[i]);
         else
            intsupport = 0.0;

         if( objparalweight > 0.0 )
            objparallelism = objparalweight * SCIPgetRowObjParallelism(scip, cuts[i]);
         else
            objparallelism = 0.0;

         efficacy = SCIPgetCutEfficacy(scip, NULL, cuts[i]);

         if( SCIProwIsLocal(cuts[i]) )
         {
            score = dircutoffdistweight * efficacy;
         }
         else
         {
            score = SCIPgetCutLPSolCutoffDistance(scip, sol, cuts[i]);
            score = dircutoffdistweight * MAX(score, efficacy);
         }

         efficacy *= efficacyweight;
         score += objparallelism + intsupport + efficacy;

         /* add small term to prefer global pool cuts */
         if( SCIProwIsInGlobalCutpool(cuts[i]) )
            score += 1e-4;

         if( randnumgen != NULL)
         {
            score += SCIPrandomGetReal(randnumgen, 0.0, 1e-6);
         }

         maxscore = MAX(maxscore, score);

         if( scores != NULL)
         {
            if( SCIPisLE(scip, score, 0.0) )
            {
               --ncuts;
               SCIPswapPointers((void**) &cuts[i], (void**) &cuts[ncuts]);
               SCIPswapReals(&scores[i], &scores[ncuts]);
            }
            else
               scores[i] = score;
         }
      }
   }
   else
   {
      /* in case there is no solution add the directed cutoff distance weight to the efficacy weight
         * since the efficacy underestimates the directed cuttoff distance
       */
      efficacyweight += dircutoffdistweight;

      /*lint -e{850} i is modified in the body of the for loop */
      for( i = ncuts-1; i >= 0; --i )
      {
         SCIP_Real score;
         SCIP_Real objparallelism;
         SCIP_Real intsupport;
         SCIP_Real efficacy;

         if( intsupportweight > 0.0 )
            intsupport = intsupportweight * SCIPgetRowNumIntCols(scip, cuts[i]) / (SCIP_Real) SCIProwGetNNonz(cuts[i]);
         else
            intsupport = 0.0;

         if( objparalweight > 0.0 )
            objparallelism = objparalweight * SCIPgetRowObjParallelism(scip, cuts[i]);
         else
            objparallelism = 0.0;

         efficacy = efficacyweight > 0.0 ? efficacyweight * SCIPgetCutEfficacy(scip, NULL, cuts[i]) : 0.0;

         score = objparallelism + intsupport + efficacy;

         /* add small term to prefer global pool cuts */
         if( SCIProwIsInGlobalCutpool(cuts[i]) )
            score += 1e-4;

         if( randnumgen != NULL)
         {
            score += SCIPrandomGetReal(randnumgen, 0.0, 1e-6);
         }

         maxscore = MAX(maxscore, score);

         if( scores != NULL)
         {
            if( SCIPisLE(scip, score, 0.0) )
            {
               --ncuts;
               SCIPswapPointers((void**) &cuts[i], (void**) &cuts[ncuts]);
               SCIPswapReals(&scores[i], &scores[ncuts]);
            }
            else
               scores[i] = score;
         }
      }
   }
   *currentncuts = ncuts;
}

/** compute projectioncut score for cuts from a given bestcut. **/
static
SCIP_RETCODE computeProjectionScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             bestcut,            /**< cut to filter orthogonality with */
   SCIP_ROW*             cut,                /**< cut to perform scoring on */
   SCIP_Real*            score              /**< score for cut */
   )
{
   SCIP_Real efficacy;
   SCIP_Real currentbestefficacy;
   SCIP_Real cosineangle;

   SCIPdebugMsg(scip, "\ncomputeProjectionScore.\n\n");
   currentbestefficacy = SCIPgetCutEfficacy(scip, NULL, bestcut);
   SCIPdebugMsg(scip, "currentbestefficacy = %g\n", currentbestefficacy);

   efficacy = SCIPgetCutEfficacy(scip, NULL, cut);
   SCIPdebugMsg(scip, "efficacy[%s] = %g\n", SCIProwGetName(cut), efficacy);

   cosineangle = SCIProwGetParallelism(bestcut, cut, 'e');
   if( SCIPisEQ(scip, cosineangle, 1.0))
      *score = -SCIPinfinity(scip);
   else
   {
      *score = sqrt(currentbestefficacy * currentbestefficacy + efficacy * efficacy
                       - 2.0 * fabs(currentbestefficacy) * fabs(efficacy) * cosineangle)
         / sqrt((1.0 - (cosineangle * cosineangle)));
      *score -= currentbestefficacy;
   }
   SCIPdebugMsg(scip, "Projectionscore[%s] = %g\n", SCIProwGetName(cut), *score);
   return SCIP_OKAY;
}

/** move the cut with the highest score to the first position in the array; there must be at least one cut */
static
void selectBestCut(
    SCIP_ROW**           cuts,               /**< array with cuts to perform selection algorithm */
    SCIP_Real*           scores,             /**< array with scores of cuts to perform selection algorithm */
    int                  ncuts               /**< number of cuts in given array */
   )
{
   int i;
   int bestpos;
   SCIP_Real bestscore;

   assert(ncuts > 0);
   assert(cuts != NULL);
   assert(scores != NULL);

   bestscore = scores[0];
   bestpos = 0;

   for( i = 1; i < ncuts; ++i )
   {
      if( scores[i] > bestscore )
      {
         bestpos = i;
         bestscore = scores[i];
      }
   }

   SCIPswapPointers((void**) &cuts[bestpos], (void**) &cuts[0]);
   SCIPswapReals(&scores[bestpos], &scores[0]);
}

/** filters the given array of cuts to enforce a maximum parallelism constraint
 *  w.r.t the given cut; moves filtered cuts to the end of the array and returns number of selected cuts */
static
int filterWithDynamicParallelism(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             bestcut,            /**< cut to filter orthogonality with */
   SCIP_ROW**            cuts,               /**< array with cuts to perform selection algorithm */
   SCIP_Real*            scores,             /**< array with scores of cuts to perform selection algorithm */
   SCIP_Real             mingain,            /**< minimum gain enforced on the two-cut efficacy */
   SCIP_Real             maxparall,          /**< maximal parallelism for all cuts that are not good */
   int                   ncuts               /**< number of cuts in given array */
   )
{
   int i;
   SCIP_Bool filter;
   SCIP_Real bestcutefficacy;

   SCIPdebugMsg(scip, "\nfilterWithDynamicParallelism.\n\n");

   assert(bestcut != NULL);
   assert(ncuts == 0 || cuts != NULL);
   assert(ncuts == 0 || scores != NULL);

   bestcutefficacy = SCIPgetCutEfficacy(scip, NULL, bestcut);

   /*lint -e{850} i is modified in the body of the for loop */
   for( i = ncuts-1; i >= 0; --i )
   {
      SCIP_Real thisparall;
      SCIP_Real cosine;
      SCIP_Real currentcutefficacy;
      SCIP_Real minmaxparall;

      currentcutefficacy = SCIPgetCutEfficacy(scip, NULL, cuts[i]);

      if( SCIPisGE(scip, bestcutefficacy, currentcutefficacy))
      {
         cosine = SCIProwGetParallelism(bestcut, cuts[i], 's');
         thisparall = cosine * bestcutefficacy / currentcutefficacy;
         SCIPdebugMsg(scip, "Thisparall(%g) = cosine(%g) * (bestcutefficacy(%g)/ currentcutefficacy(%g))\n\n", thisparall,
                      cosine, bestcutefficacy, currentcutefficacy);
      }
      else
      {
         cosine = SCIProwGetParallelism(cuts[i], bestcut, 's');
         thisparall = cosine * currentcutefficacy / bestcutefficacy;
         SCIPdebugMsg(scip, "Thisparall(%g) = cosine(%g) * (currentcutefficacy(%g) / bestcutefficacy(%g))\n\n", thisparall,
                      cosine, currentcutefficacy, bestcutefficacy);
      }

      /* compute the max-minimum angle for given the given cuts to enforce
       * norm(d) >= (1+mingain)*eff1 for non-negative cosine angle */
      minmaxparall = MAX( (bestcutefficacy * bestcutefficacy
                        + currentcutefficacy * currentcutefficacy
                        - (1 + mingain) * bestcutefficacy * (1 + mingain) * bestcutefficacy * (1 - cosine * cosine))
                        / (2 * bestcutefficacy * currentcutefficacy),
                        maxparall );
      filter = ( SCIPisGE(scip, thisparall, 1.0) || SCIPisGT(scip, cosine, minmaxparall) );

      SCIPdebugMsg(scip, "Filter = %u\n", filter);

      if( filter )
      {
         --ncuts;
         SCIPswapPointers((void**) &cuts[i], (void**) &cuts[ncuts]);
         SCIPswapReals(&scores[i], &scores[ncuts]);
      }
   }

   return ncuts;
}


/*
 * Callback methods of cut selector
 */

/** copy method for cut selector plugin (called when SCIP copies plugins) */
static
SCIP_DECL_CUTSELCOPY(cutselCopyDynamic)
{  /*lint --e{715}*/
  assert(scip != NULL);
  assert(cutsel != NULL);
  assert(strcmp(SCIPcutselGetName(cutsel), CUTSEL_NAME) == 0);

  /* call inclusion method of cut selector */
  SCIP_CALL( SCIPincludeCutselDynamic(scip) );

  return SCIP_OKAY;
}

/** destructor of cut selector to free user data (called when SCIP is exiting) */
/**! [SnippetCutselFreeDynamic] */
static
SCIP_DECL_CUTSELFREE(cutselFreeDynamic)
{  /*lint --e{715}*/
  SCIP_CUTSELDATA* cutseldata;

  cutseldata = SCIPcutselGetData(cutsel);

  SCIPfreeBlockMemory(scip, &cutseldata);

  SCIPcutselSetData(cutsel, NULL);

  return SCIP_OKAY;
}
/**! [SnippetCutselFreeDynamic] */

/** initialization method of cut selector (called after problem was transformed) */
static
SCIP_DECL_CUTSELINIT(cutselInitDynamic)
{  /*lint --e{715}*/
  SCIP_CUTSELDATA* cutseldata;

  cutseldata = SCIPcutselGetData(cutsel);
  assert(cutseldata != NULL);

  SCIP_CALL( SCIPcreateRandom(scip, &(cutseldata)->randnumgen, RANDSEED, TRUE) );

  return SCIP_OKAY;
}

/** deinitialization method of cut selector (called before transformed problem is freed) */
static
SCIP_DECL_CUTSELEXIT(cutselExitDynamic)
{  /*lint --e{715}*/
  SCIP_CUTSELDATA* cutseldata;

  cutseldata = SCIPcutselGetData(cutsel);
  assert(cutseldata != NULL);
  assert(cutseldata->randnumgen != NULL);

  SCIPfreeRandom(scip, &cutseldata->randnumgen);

  return SCIP_OKAY;
}

/** cut selection method of cut selector */
static
SCIP_DECL_CUTSELSELECT(cutselSelectDynamic) { /*lint --e{715}*/
   SCIP_CUTSELDATA *cutseldata;

   assert(cutsel != NULL);
   assert(result != NULL);

   *result = SCIP_SUCCESS;

   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);
   if (cutseldata->maxdepth != -1 && cutseldata->maxdepth < SCIPgetDepth(scip))
   {
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPselectCutsDynamic( scip, cuts, forcedcuts, cutseldata->randnumgen, cutseldata->filtermode,
                                      cutseldata->mingain, 1-cutseldata->minortho, cutseldata->dircutoffdistweight, cutseldata->efficacyweight,
                                      cutseldata->objparalweight, cutseldata->intsupportweight, ncuts, nforcedcuts,
                                      maxnselectedcuts, nselectedcuts) );

   return SCIP_OKAY;
}


/*
 * cut selector specific interface methods
 */

/** creates the dynamic cut selector and includes it in SCIP */
SCIP_RETCODE SCIPincludeCutselDynamic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CUTSELDATA* cutseldata;
   SCIP_CUTSEL* cutsel;

   /* create dynamic cut selector data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &cutseldata) );
   BMSclearMemory(cutseldata);

   SCIP_CALL( SCIPincludeCutselBasic(scip, &cutsel, CUTSEL_NAME, CUTSEL_DESC, CUTSEL_PRIORITY, cutselSelectDynamic,
                                    cutseldata) );

   assert(cutsel != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetCutselCopy(scip, cutsel, cutselCopyDynamic) );

   SCIP_CALL( SCIPsetCutselFree(scip, cutsel, cutselFreeDynamic) );
   SCIP_CALL( SCIPsetCutselInit(scip, cutsel, cutselInitDynamic) );
   SCIP_CALL( SCIPsetCutselExit(scip, cutsel, cutselExitDynamic) );

   /* add dynamic cut selector parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
               "cutselection/" CUTSEL_NAME "/efficacyweight",
               "weight of efficacy in cut score calculation",
               &cutseldata->efficacyweight, FALSE,
               DEFAULT_EFFICACYWEIGHT, 0.0, SCIP_INVALID / 10.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
               "cutselection/" CUTSEL_NAME "/dircutoffdistweight",
               "weight of directed cutoff distance in cut score calculation",
               &cutseldata->dircutoffdistweight, FALSE,
               DEFAULT_DIRCUTOFFDISTWEIGHT, 0.0, SCIP_INVALID / 10.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
               "cutselection/" CUTSEL_NAME "/objparalweight",
               "weight of objective parallelism in cut score calculation",
               &cutseldata->objparalweight, FALSE,
               DEFAULT_OBJPARALWEIGHT, 0.0, SCIP_INVALID / 10.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
               "cutselection/" CUTSEL_NAME "/intsupportweight",
               "weight of integral support in cut score calculation",
               &cutseldata->intsupportweight, FALSE,
               DEFAULT_INTSUPPORTWEIGHT, 0.0, SCIP_INVALID / 10.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
               "cutselection/" CUTSEL_NAME "/mingain",
               "minimal efficacy gain for a cut to enter the LP",
               &cutseldata->mingain, FALSE,
               DEFAULT_MINGAIN, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip,
               "cutselection/" CUTSEL_NAME "/filtermode",
               "filtering strategy during cut selection",
               &cutseldata->filtermode, FALSE,
               DEFAULT_FILTERMODE, "df", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
               "cutselection/" CUTSEL_NAME "/minortho",
               "minimal orthogonality for a cut to enter the LP",
               &cutseldata->minortho, FALSE,
               DEFAULT_MINORTHO, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
               "cutselection/" CUTSEL_NAME "/maxdepth",
               "maximum depth at which this cutselector is employed",
               &cutseldata->maxdepth, FALSE,
               DEFAULT_MAXDEPTH, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   return SCIP_OKAY;
}


/** perform a cut selection algorithm for the given array of cuts
 *
 *  This is the selection method of the dynamic cut selector which implements
 *  the dynamic orthognality filtering based on the ratio of efficacies.
 *  The input cuts array gets resorted s.t the selected cuts come first and the remaining
 *  ones are the end.
 */
SCIP_RETCODE SCIPselectCutsDynamic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            cuts,               /**< array with cuts to perform selection algorithm */
   SCIP_ROW**            forcedcuts,         /**< array with forced cuts */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator for tie-breaking, or NULL */
   char                  filtermode,         /**< filtering strategy during cut selection (
                                             *  'd'ynamic- and 'f'ull dynamic parallelism) */
   SCIP_Real             mingain,            /**< minimum efficacy gain in percentage to filter cuts */
   SCIP_Real             maxparall,          /**< maximal parallelism for all cuts that are not good */
   SCIP_Real             dircutoffdistweight,/**< weight of directed cutoff distance in cut score calculation */
   SCIP_Real             efficacyweight,     /**< weight of efficacy in cut score calculation */
   SCIP_Real             objparalweight,     /**< weight of objective parallelism in cut score calculation */
   SCIP_Real             intsupportweight,   /**< weight of integral support in cut score calculation */
   int                   ncuts,              /**< number of cuts in cuts array */
   int                   nforcedcuts,        /**< number of forced cuts */
   int                   maxselectedcuts,    /**< maximal number of cuts from cuts array to select */
   int*                  nselectedcuts      /**< pointer to return number of selected cuts from cuts array */
   )
{
   SCIP_ROW* selectedcut;
   SCIP_Real* scores;
   SCIP_Real* forcedscores;
   SCIP_Real* scoresptr;
   int ngoodforcedcuts;
   int i;

   assert(cuts != NULL && ncuts > 0);
   assert(forcedcuts != NULL || nforcedcuts == 0);
   assert(nselectedcuts != NULL);

   *nselectedcuts = 0;
   ngoodforcedcuts = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &scores, ncuts) );

   /* compute scores of cuts and max score of cuts and forced cuts (used to define goodscore) */
   scoring(scip, cuts, randnumgen, dircutoffdistweight, efficacyweight, objparalweight, intsupportweight, &ncuts,
           scores);
   scoresptr = scores;

   SCIPdebugMsg(scip, "nforcedcuts = %i.\n", nforcedcuts);

   /* perform cut selection algorithm for the cuts */

   /* forced cuts are going to be selected so use them to filter cuts */
   for( i = 0; i < nforcedcuts && ncuts > 0; ++i )
      ncuts = filterWithDynamicParallelism(scip, forcedcuts[i], cuts, scores, mingain, maxparall, ncuts);

   /* if all cuts are already filtered, we can stop */
   if( ncuts <= 0 )
      goto TERMINATE;

   /* if the maximal number of cuts was selected, we can stop here */
   if( *nselectedcuts == maxselectedcuts )
      goto TERMINATE;

   if( filtermode == 'f' && nforcedcuts > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &forcedscores, nforcedcuts) );
      ngoodforcedcuts = nforcedcuts;
      scoring(scip, forcedcuts, randnumgen, dircutoffdistweight, efficacyweight, objparalweight, intsupportweight,
              &ngoodforcedcuts, forcedscores);

      if( ngoodforcedcuts != 0 )
      {
        selectBestCut(forcedcuts, forcedscores, ngoodforcedcuts);
        SCIPfreeBufferArray(scip, &forcedscores);
        SCIPdebugMsg(scip, "best forced cut: %s.\n", SCIProwGetName(forcedcuts[0]));

        for( i = 0; i < ncuts; i++ )
        {
           SCIP_CALL( computeProjectionScore(scip, forcedcuts[0], cuts[i], &scores[i]) );
           SCIPdebugMsg(scip, "scores[%i] = %g\n", i, scores[i]);
        }
      }
   }

   if( ngoodforcedcuts == 0 )
   {
      assert(filtermode == 'd' || ngoodforcedcuts == 0);
      selectBestCut(cuts, scores, ncuts);

      selectedcut = cuts[0];
      SCIPdebugMsg(scip, "selectedcut = %s.\n", SCIProwGetName(selectedcut));

      ++(*nselectedcuts);

      /* if the maximal number of cuts was selected, we can stop here */
      if( *nselectedcuts == maxselectedcuts )
         goto TERMINATE;

      /* move the pointers to the next position and filter the remaining cuts to enforce the dynamic parallelism constraint */
      ++cuts;
      ++scores;
      --ncuts;

      ncuts = filterWithDynamicParallelism(scip, selectedcut, cuts, scores, mingain, maxparall, ncuts);

      if( filtermode == 'f' )
      {
         for( i = 0; i < ncuts; i++ )
         {
            SCIP_CALL( computeProjectionScore(scip, selectedcut, cuts[i], &scores[i]) );
         }
      }
   }

   SCIPdebugMsg(scip, "ncuts after forced cut filter = %i.\n", ncuts);

   /* now greedily select the remaining cuts */
   while( ncuts > 0 )
   {
      selectBestCut(cuts, scores, ncuts);
      selectedcut = cuts[0];
      SCIPdebugMsg(scip, "selectedcut = %s.\n", SCIProwGetName(selectedcut));

      ++(*nselectedcuts);

      /* if the maximal number of cuts was selected, we can stop here */
      if( *nselectedcuts == maxselectedcuts )
         goto TERMINATE;

      /* move the pointers to the next position and filter the remaining cuts to enforce the dynamic parallelism constraint */
      ++cuts;
      ++scores;
      --ncuts;

      ncuts = filterWithDynamicParallelism(scip, selectedcut, cuts, scores, mingain, maxparall, ncuts);

      if( filtermode == 'f' )
      {
         for( i = 0; i < ncuts; i++ )
         {
            SCIP_CALL( computeProjectionScore(scip, selectedcut, cuts[i], &scores[i]) );
            SCIPdebugMsg(scip, "nonforcedscores[%i] = %g\n", i, scores[i]);
         }
      }
   }

   TERMINATE:
   SCIPfreeBufferArray(scip, &scoresptr);
   return SCIP_OKAY;
}
