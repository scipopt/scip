/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_indicatordiving.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  indicator diving heuristic
 * @author Katrin Halbig
 * @author Alexander Hoen
 * Diving heuristic: Iteratively fixes some fractional variable and resolves the LP-relaxation, thereby simulating a
 * depth-first-search in the tree.
 *
 * Indicatordiving:
 * Iterates over every indicator constraints and checks if this indicatorvariable models a semicontinuous variable by
 * setting an upperbound for the constant.
 * (Unfortunately this is data needs to be checked and is not contained in the Vbound structure.)
 * Now the heuristic checks if the indicator constraint is already fulfilled (lp solution is greater than lowerbound).
 * If not, it calculates a score within the percentage difference and the rounding direction according to mode.
 * Otherwise, the indicator constraint is fulfilled and therefore the rounding direction is chosen correspondingly with
 * random (low) score.
 * */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_indicator.h"
#include "scip/heur_indicatordiving.h"
#include "scip/heuristics.h"
#include "scip/pub_cons.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_heur.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_sol.h"
#include "scip/scip_tree.h"
#include "scip_prob.h"
#include "struct_heur.h"
#include <string.h>

#define HEUR_NAME             "indicatordiving"
#define HEUR_DESC             "indicator diving heuristic"
#define HEUR_DISPCHAR         '?' /**< todo: change to SCIP_HEURDISPCHAR_DIVING */
#define HEUR_PRIORITY         0
#define HEUR_FREQ             10
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */
#define DIVESET_DIVETYPES     SCIP_DIVETYPE_INTEGRALITY /**< bit mask that represents all supported dive types */
#define DIVESET_ISPUBLIC      FALSE   /**< is this dive set publicly available (ie., can be used by other primal heuristics?) */


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
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */
#define DEFAULT_LPRESOLVEDOMCHGQUOT 0.15 /**< percentage of immediate domain changes during probing to trigger LP resolve */
#define DEFAULT_LPSOLVEFREQ           0 /**< LP solve frequency for diving heuristics */
#define DEFAULT_ONLYLPBRANCHCANDS FALSE /**< should only LP branching candidates be considered instead of the slower but
                                         *   more general constraint handler diving variable selection? */
#define DEFAULT_RANDSEED             11  /**< initial seed for random number generation */

/*
 * Heuristic specific parameters
 */
#define DEFAULT_ROUNDINGFRAC        0.5 /**< default parameter setting for parameter roundingfrac */
#define DEFAULT_MODE                  3 /**< default parameter setting for parameter mode */
#define DEFAULT_SEMICONTSCOREMODE     0 /**< default parameter setting for parameter semicontscoremode */

enum IndicatorDivingMode
{
   ROUNDING_DOWN = 0,
   ROUNDING_UP = 1,
   ROUNDING_FRAC_AGGRESSIVE = 2,
   ROUNDING_FRAC_CONSERVATIVE = 3
};
typedef enum IndicatorDivingMode INDICATORDIVINGMODE;

/** data structure to store information of a semicontinuous variable
 *
 * For a variable x (not stored in the struct), this stores the data of nbnds implications
 *   bvars[i] = 0 -> x = vals[i]
 *   bvars[i] = 1 -> lbs[i] <= x <= ubs[i]
 * where bvars[i] are binary variables.
 */
struct SCVarData
{
   SCIP_Real*            vals0;              /**< values of the variable when the corresponding bvars[i] = 0 */
   SCIP_Real*            lbs1;               /**< global lower bounds of the variable when the corresponding bvars[i] = 1 */
   SCIP_Real*            ubs1;               /**< global upper bounds of the variable when the corresponding bvars[i] = 1 */
   SCIP_VAR**            bvars;              /**< the binary variables on which the variable domain depends */
   int                   nbnds;              /**< number of suitable on/off bounds the var has */
   int                   bndssize;           /**< size of the arrays */
};
typedef struct SCVarData SCVARDATA;


/** locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_HASHMAP*         scvars;             /**< stores hashmap with semicontinuous variables */
   SCIP_Real             roundingfrac;       /**< in fractional case all fractional below this value are rounded up*/
   int                   mode;               /**< decides which mode is selected (0: down, 1: up, 2: aggressive, 3: conservative (default)) */
   int                   semicontscoremode;  /**< which values of semi-continuous variables should get a high score? (0: low (default), 1: middle, 2: high) */
};

/*
 * Local methods
 */

/** releases all data from given hashmap filled with SCVarData and the hashmap itself */
static
SCIP_RETCODE releaseSCHashmap(
  SCIP*                  scip,               /**< SCIP data structure */
  SCIP_HASHMAP*          hashmap             /**< hashmap to be freed */
  )
{
   SCIP_HASHMAPENTRY* entry;
   SCVARDATA* data;
   int c;

   if( hashmap != NULL )
   {
      for( c = 0; c < SCIPhashmapGetNEntries( hashmap ); c++ )
      {
         entry = SCIPhashmapGetEntry( hashmap, c);
         if( entry != NULL )
         {
            data = (SCVARDATA*) SCIPhashmapEntryGetImage(entry);
            SCIPfreeBlockMemoryArray(scip, &data->ubs1, data->bndssize);
            SCIPfreeBlockMemoryArray(scip, &data->lbs1, data->bndssize);
            SCIPfreeBlockMemoryArray(scip, &data->vals0, data->bndssize);
            SCIPfreeBlockMemoryArray(scip, &data->bvars, data->bndssize);
            SCIPfreeBlockMemory(scip, &data);
         }
      }
      SCIPhashmapFree(&hashmap);
      assert(hashmap == NULL);
   }
   return SCIP_OKAY;
}


/** checks if variable is indicator variable and stores corresponding indicator constraint */
static
SCIP_RETCODE checkAndGetIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             cand,               /**< candidate variable */
   SCIP_CONS**           cons,               /**< pointer to store indicator constraint */
   SCIP_Bool*            isindicator         /**< pointer to store whether candidate variable is indicator variable */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** indicatorconss;
   int c;

   assert(scip != NULL);
   assert(cand != NULL);
   assert(cons != NULL);
   assert(isindicator != NULL);

   if( SCIPvarGetType(cand) != SCIP_VARTYPE_BINARY )
   {
      *cons = NULL;
      *isindicator = FALSE;
      return SCIP_OKAY;
   }

   conshdlr = SCIPfindConshdlr(scip, "indicator");
   indicatorconss = SCIPconshdlrGetConss(conshdlr);

   *isindicator = FALSE;
   for( c = 0; c < SCIPconshdlrGetNActiveConss(conshdlr); c++ )
   {
      SCIP_VAR* indicatorvar;
      indicatorvar = SCIPgetBinaryVarIndicator(indicatorconss[c]);

      if( cand == indicatorvar )
      {
         *cons = indicatorconss[c];
         *isindicator = TRUE;
         return SCIP_OKAY;
      }
   }

   *cons = NULL;
   *isindicator = FALSE;
   return SCIP_OKAY;
}

/** adds an indicator to the data of a semicontinuous variable */
static
SCIP_RETCODE addSCVarIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCVARDATA*            scvdata,            /**< semicontinuous variable data */
   SCIP_VAR*             indicator,          /**< indicator to be added */
   SCIP_Real             val0,               /**< value of the variable when indicator == 0 */
   SCIP_Real             lb1,                /**< lower bound of the variable when indicator == 1 */
   SCIP_Real             ub1                 /**< upper bound of the variable when indicator == 1 */
   )
{
   int newsize;
   int i;
   SCIP_Bool found;
   int pos;

   assert(scvdata != NULL);
   assert(indicator != NULL);

   /* find the position where to insert */
   if( scvdata->bvars == NULL )
   {
      assert(scvdata->nbnds == 0 && scvdata->bndssize == 0);
      found = FALSE;
      pos = 0;
   }
   else
   {
      found = SCIPsortedvecFindPtr((void**)scvdata->bvars, SCIPvarComp, (void*)indicator, scvdata->nbnds, &pos);
   }

   if( found )
      return SCIP_OKAY;

   /* ensure sizes */
   if( scvdata->nbnds + 1 > scvdata->bndssize )
   {
      newsize = SCIPcalcMemGrowSize(scip, scvdata->nbnds + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->bvars, scvdata->bndssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->vals0, scvdata->bndssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->lbs1, scvdata->bndssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->ubs1, scvdata->bndssize, newsize) );
      scvdata->bndssize = newsize;
   }
   assert(scvdata->nbnds + 1 <= scvdata->bndssize);
   assert(scvdata->bvars != NULL);

   /* move entries if needed */
   for( i = scvdata->nbnds; i > pos; --i )
   {
      scvdata->bvars[i] = scvdata->bvars[i-1];
      scvdata->vals0[i] = scvdata->vals0[i-1];
      scvdata->lbs1[i] = scvdata->lbs1[i-1];
      scvdata->ubs1[i] = scvdata->ubs1[i-1];
   }

   scvdata->bvars[pos] = indicator;
   scvdata->vals0[pos] = val0;
   scvdata->lbs1[pos] = lb1;
   scvdata->ubs1[pos] = ub1;
   ++scvdata->nbnds;

   return SCIP_OKAY;
}

/** checks if a variable is semicontinuous and stores it data in the hashmap scvars
 *
 * A variable x is semicontinuous if its bounds depend on at least one binary variable called the indicator,
 * and indicator == 0 => x == x^0 for some real constant x^0.
 */
static
SCIP_RETCODE varIsSemicontinuous(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< the variable to check */
   SCIP_HASHMAP*         scvars,             /**< semicontinuous variable information */
   SCIP_Bool*            result              /**< buffer to store whether var is semicontinuous */
   )
{
   SCIP_Real lb0;
   SCIP_Real ub0;
   SCIP_Real lb1;
   SCIP_Real ub1;
   SCIP_Real glb;
   SCIP_Real gub;
   SCIP_Bool exists;
   int c;
   int pos;
   SCIP_VAR** vlbvars;
   SCIP_VAR** vubvars;
   SCIP_Real* vlbcoefs;
   SCIP_Real* vubcoefs;
   SCIP_Real* vlbconstants;
   SCIP_Real* vubconstants;
   int nvlbs;
   int nvubs;
   SCVARDATA* scvdata;
   SCIP_VAR* bvar;

   assert(scip != NULL);
   assert(var != NULL);
   assert(scvars != NULL);
   assert(result != NULL);

   scvdata = (SCVARDATA*) SCIPhashmapGetImage(scvars, (void*)var);
   if( scvdata != NULL )
   {
      *result = TRUE;
      return SCIP_OKAY;
   }

   vlbvars = SCIPvarGetVlbVars(var);
   vubvars = SCIPvarGetVubVars(var);
   vlbcoefs = SCIPvarGetVlbCoefs(var);
   vubcoefs = SCIPvarGetVubCoefs(var);
   vlbconstants = SCIPvarGetVlbConstants(var);
   vubconstants = SCIPvarGetVubConstants(var);
   nvlbs = SCIPvarGetNVlbs(var);
   nvubs = SCIPvarGetNVubs(var);
   glb = SCIPvarGetLbGlobal(var);
   gub = SCIPvarGetUbGlobal(var);

   pos = -1;

   *result = FALSE;

   /* Scan through lower bounds; for each binary vlbvar save the corresponding lb0 and lb1.
    * Then check if there is an upper bound with this vlbvar and save ub0 and ub1.
    * If the found bounds imply that the var value is fixed to some val0 when vlbvar = 0,
    * save vlbvar and val0 to scvdata.
    */
   for( c = 0; c < nvlbs; ++c )
   {
      if( SCIPvarGetType(vlbvars[c]) != SCIP_VARTYPE_BINARY )
         continue;

      bvar = vlbvars[c];

      lb0 = MAX(vlbconstants[c], glb);
      lb1 = MAX(vlbconstants[c] + vlbcoefs[c], glb);

      /* look for bvar in vubvars */
      if( vubvars != NULL )
         exists = SCIPsortedvecFindPtr((void**)vubvars, SCIPvarComp, bvar, nvubs, &pos);
      else
         exists = FALSE;
      if( exists )
      {
         /* save the upper bounds */
         ub0 = MIN(vubconstants[pos], gub);
         ub1 = MIN(vubconstants[pos] + vubcoefs[pos], gub);
      }
      else
      {
         /* if there is no upper bound with vubvar = bvar, use global var bounds */
         ub0 = gub;
         ub1 = gub;
      }

      /* the 'off' domain of a semicontinuous var should reduce to a single point and be different from the 'on' domain */
      //TODO: this doesn't work because the ub0 is not detected. -> therefore ignore this and check it outside
//      if( SCIPisEQ(scip, lb0, ub0) && (!SCIPisEQ(scip, lb0, lb1) || !SCIPisEQ(scip, ub0, ub1)) )
      if( (!SCIPisEQ(scip, lb0, lb1) || !SCIPisEQ(scip, ub0, ub1)) )
      {
         if( scvdata == NULL )
         {
            SCIP_CALL( SCIPallocClearBlockMemory(scip, &scvdata) );
         }
         SCIP_CALL( addSCVarIndicator(scip, scvdata, bvar, lb0, lb1, ub1) );
      }
   }

   /* look for vubvars that have not been processed yet */
   assert(vubvars != NULL || nvubs == 0);
   for( c = 0; c < nvubs; ++c )
   {
      if( SCIPvarGetType(vubvars[c]) != SCIP_VARTYPE_BINARY )  /*lint !e613*/
         continue;

      bvar = vubvars[c];  /*lint !e613*/

      /* skip vars that are in vlbvars */
      if( vlbvars != NULL && SCIPsortedvecFindPtr((void**)vlbvars, SCIPvarComp, bvar, nvlbs, &pos) )
         continue;

      lb0 = glb;
      lb1 = glb;
      ub0 = MIN(vubconstants[c], gub);
      ub1 = MIN(vubconstants[c] + vubcoefs[c], gub);

      /* the 'off' domain of a semicontinuous var should reduce to a single point and be different from the 'on' domain */
//      if( SCIPisEQ(scip, lb0, ub0) && (!SCIPisEQ(scip, lb0, lb1) || !SCIPisEQ(scip, ub0, ub1)) )
      //TODO: indicator not considered
      if( (!SCIPisEQ(scip, lb0, lb1) || !SCIPisEQ(scip, ub0, ub1)) )
      {
         if( scvdata == NULL )
         {
            SCIP_CALL( SCIPallocClearBlockMemory(scip, &scvdata) );
         }

         SCIP_CALL( addSCVarIndicator(scip, scvdata, bvar, lb0, lb1, ub1) );
      }
   }

   if( scvdata != NULL )
   {
#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "var <%s> has global bounds [%f, %f] and the following on/off bounds:\n", SCIPvarGetName(var), glb, gub);
      for( c = 0; c < scvdata->nbnds; ++c )
      {
         SCIPdebugMsg(scip, " c = %d, bvar <%s>: val0 = %f\n", c, SCIPvarGetName(scvdata->bvars[c]), scvdata->vals0[c]);
      }
#endif
      SCIP_CALL( SCIPhashmapInsert(scvars, var, scvdata) );
      *result = TRUE;
   }

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyIndicatordiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurIndicatordiving(scip) );

   return SCIP_OKAY;
}


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeIndicatordiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitIndicatordiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   //TODO is the hashmap creation correct?

   SCIP_CALL( SCIPhashmapCreate( &heurdata->scvars, SCIPblkmem( scip ), SCIPgetNVars(scip) ));


   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitIndicatordiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   SCIP_CALL( releaseSCHashmap(scip, heurdata->scvars) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecIndicatordiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   assert(SCIPheurGetNDivesets(heur) > 0);
   assert(SCIPheurGetDivesets(heur) != NULL);
   diveset = SCIPheurGetDivesets(heur)[0];
   assert(diveset != NULL);

   SCIPdebugMessage("call heurExecIndicatordiving at depth %d \n", SCIPgetDepth(scip));

   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible, -1L, SCIP_DIVECONTEXT_SINGLE) );

   SCIPdebugMessage("leave heurExecIndicatordiving\n");

   return SCIP_OKAY;
}

/** calculate score and preferred rounding direction for the candidate variable */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScoreIndicatordiving)
{
   /*input:
    * - scip : SCIP main data structure
    * - diveset : diving settings for scoring
    * - divetype : represents different methods for a dive set to explore the next children
    * - cand : candidate variable for which the score should be determined
    * - candsol : solution value of variable in LP relaxation solution
    * - candsfrac : fractional part of solution value of variable
    * - score : pointer for diving score value - the best candidate maximizes this score
    * - roundup : pointer to store whether the preferred rounding direction is upwards
    */

   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_VAR** consvars;
   SCIP_CONS* indicatorcons;
   SCIP_CONS* lincons;
   SCIP_VAR* slackvar;
   SCIP_VAR* semicontinuousvar;
   SCIP_Real lpsolsemicontinuous;
   SCVARDATA* scdata;
   SCIP_Real* consvals;
   SCIP_Real rhs;
   int nconsvars;
   SCIP_Bool isindicatorvar;
   SCIP_Bool success;
   int v;

   semicontinuousvar = NULL;
   lpsolsemicontinuous = 0.0;

   /* check if cand variable is indicator variable */
   SCIP_CALL( checkAndGetIndicator(scip, cand, &indicatorcons, &isindicatorvar) );

   if( !isindicatorvar )
   {
      /* todo: set useful score for non-indicator variables */
      *score = -1;
      *roundup = (candsfrac > 0.5);
      return SCIP_OKAY;
   }

   SCIPdebugMessage("cand: %s, candsol: %.2f, candobjcoeff: %f\n", SCIPvarGetName(cand), candsol, SCIPvarGetObj(cand));

   lincons = SCIPgetLinearConsIndicator(indicatorcons);
   slackvar = SCIPgetSlackVarIndicator(indicatorcons);
   rhs = SCIPconsGetRhs(scip, lincons, &success);


   /* get heuristic data */
   heur = SCIPdivesetGetHeur(diveset);
   assert(heur != NULL);
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   randnumgen = SCIPdivesetGetRandnumgen(diveset);
   assert(randnumgen != NULL);

   SCIPdebugPrintCons(scip, indicatorcons, NULL);
   SCIPdebugPrintCons(scip, lincons, NULL);

   SCIP_CALL( SCIPgetConsNVars(scip, lincons, &nconsvars, &success) );

   if ( nconsvars != 2 )
   {
      *score = SCIPrandomGetReal( randnumgen, -1.0, 0.0 );
      *roundup = (candsfrac > 0.5);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
   SCIP_CALL( SCIPgetConsVars(scip, lincons, consvars, nconsvars, &success) );
   SCIP_CALL( SCIPgetConsVals(scip, lincons, consvals, nconsvars, &success) );

   for( v = 0; v < nconsvars ; v++ )
   {
      if( consvars[v] == slackvar )
         continue;

      semicontinuousvar = consvars[v];
      lpsolsemicontinuous = SCIPvarGetLPSol( semicontinuousvar );
      SCIPdebugMessage( "%s lp sol %f %f\n", SCIPvarGetName( semicontinuousvar ), lpsolsemicontinuous,
                        consvals[v] );
      SCIP_CALL( varIsSemicontinuous(scip, semicontinuousvar, heurdata->scvars, &success) );
   }
   assert(semicontinuousvar != NULL);

   scdata = (SCVARDATA*) SCIPhashmapGetImage(heurdata->scvars, (void*) semicontinuousvar);
   assert(scdata == NULL  || SCIPisGE(scip, lpsolsemicontinuous, scdata->vals0[0]));
   assert(scdata == NULL  || SCIPisLE(scip, lpsolsemicontinuous, scdata->ubs1[0]));

   //TODO: only allow sc variables and do the check if ub is equal to lowerbound of the sc
   if( scdata == NULL || !SCIPisEQ(scip, rhs, scdata->vals0[0]) ||
      !(scdata->bvars[0] == cand || (SCIPvarIsNegated(cand) && scdata->bvars[0] == SCIPvarGetNegationVar(cand))) )
   {
      //TODO: only continue if semicontinuous variable.
      *score = SCIPrandomGetReal(randnumgen, -1.0, 0.0);
      *roundup = (candsfrac > 0.5);
      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
      return SCIP_OKAY;
   }

   //Case: Variable is at least lb1
   if( SCIPisGE(scip, lpsolsemicontinuous, scdata->lbs1[0]) )
   {
      *score = SCIPrandomGetReal(randnumgen, -1.0, 0.0);
      *roundup = FALSE;
   }
   //Case: Variable is equal to constant
   else if( SCIPisEQ(scip, lpsolsemicontinuous, scdata->vals0[0]) )
   {
      *score = SCIPrandomGetReal(randnumgen, -1.0, 0.0);
      *roundup = TRUE;
   }
   //Case: Variable is between constant and lb1
   else
   {
      *score = 100 * (scdata->lbs1[0] - lpsolsemicontinuous) / scdata->lbs1[0];
      assert(*score>0);

      switch( (INDICATORDIVINGMODE)heurdata->mode )
      {
      case ROUNDING_DOWN:
         *roundup = FALSE;
         break;
      case ROUNDING_UP:
         *roundup = TRUE;
         break;
      case ROUNDING_FRAC_AGGRESSIVE:
         *roundup = (*score <= heurdata->roundingfrac * 100);
         break;
      case ROUNDING_FRAC_CONSERVATIVE:
         *roundup = (*score > heurdata->roundingfrac * 100);
         break;
      }

      switch( heurdata->semicontscoremode )
      {
      case 0:
         break;
      case 1:
         if( lpsolsemicontinuous < scdata->lbs1[0] * heurdata->roundingfrac )
            *score = 100 * (lpsolsemicontinuous / (heurdata->roundingfrac * scdata->lbs1[0]));
         else
            *score = 100 * (-lpsolsemicontinuous / ((1 - heurdata->roundingfrac) * scdata->lbs1[0]) + (1 / (1 - heurdata->roundingfrac)) );
         break;
      case 2:
         *score = 100 - *score;
         break;
      default:
         return SCIP_INVALIDDATA;
      }
      assert(*score>0);
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}

#define divesetAvailableIndicatordiving NULL

/*
 * heuristic specific interface methods
 */

/** creates the indicatordiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurIndicatordiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create indicatordiving primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   heur = NULL;


   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecIndicatordiving, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyIndicatordiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeIndicatordiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitIndicatordiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitIndicatordiving) );

   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, HEUR_NAME, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL, DEFAULT_LPRESOLVEDOMCHGQUOT,
         DEFAULT_LPSOLVEFREQ, DEFAULT_MAXLPITEROFS, DEFAULT_RANDSEED, DEFAULT_BACKTRACK, DEFAULT_ONLYLPBRANCHCANDS,
         DIVESET_ISPUBLIC, DIVESET_DIVETYPES, divesetGetScoreIndicatordiving, divesetAvailableIndicatordiving) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/roundingfrac",
         "in fractional case all fractional below this value are rounded up",
         &heurdata->roundingfrac, FALSE, DEFAULT_ROUNDINGFRAC, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/mode",
         "decides which mode is selected (0: down, 1: up, 2: aggressive, 3: conservative (default))",
         &heurdata->mode, FALSE, DEFAULT_MODE, 0, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/semicontscoremode",
         "which values of semi-continuous variables should get a high score? (0: low (default), 1: middle, 2: high)",
         &heurdata->semicontscoremode, FALSE, DEFAULT_SEMICONTSCOREMODE, 0, 2, NULL, NULL) );

   return SCIP_OKAY;
}
