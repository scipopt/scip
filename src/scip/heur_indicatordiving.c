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

/**@file   heur_indicatordiving.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  LP diving heuristic that fixes indicator variables controlling semicontinuous variables
 * @author Katrin Halbig
 * @author Alexander Hoen
 *
 * A diving heuristic iteratively rounds some fractional variables or variables determined by constraint handlers,
 * and resolves the LP relaxation. Thereby simulating a depth-first-search in the tree.
 *
 * Indicatordiving focuses on indicator variables, which control semicontinuous variables.
 * If the semicontinuous variable is unbounded, the indicator constraint is not part of the LP and,
 * therefore, the indicator variable is not set to an useful value in the LP solution.
 *
 * For these indicator variables the score depends on the LP value and the bounds of the corresponding semicontinuous variable.
 * If parameter usevarbounds=TRUE, also varbound constraints modeling semicontinuous variables are considered.
 * For all other variables the Farkas score (scaled) is returned.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_indicator.h"
#include "scip/cons_varbound.h"
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
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_tree.h"
#include "scip/scip_prob.h"
#include "scip/scip_message.h"

#define HEUR_NAME             "indicatordiving"
#define HEUR_DESC             "LP diving heuristic that fixes indicator variables controlling semicontinuous variables"
#define HEUR_DISPCHAR         'I'
#define HEUR_PRIORITY         -150000
#define HEUR_FREQ             0
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
#define DEFAULT_LPSOLVEFREQ          30 /**< LP solve frequency for diving heuristics */
#define DEFAULT_ONLYLPBRANCHCANDS FALSE /**< should only LP branching candidates be considered instead of the slower but
                                         *   more general constraint handler diving variable selection? */
#define DEFAULT_RANDSEED             11 /**< initial seed for random number generation */

/*
 * Heuristic specific parameters
 */
#define DEFAULT_ROUNDINGFRAC        0.5 /**< default setting for parameter roundingfrac */
#define DEFAULT_ROUNDINGMODE          0 /**< default setting for parameter roundingmode */
#define DEFAULT_SEMICONTSCOREMODE     0 /**< default setting for parameter semicontscoremode */
#define DEFAULT_USEVARBOUNDS       TRUE /**< default setting for parameter usevarbounds */
#define DEFAULT_RUNWITHOUTSCINDS  FALSE /**< default setting for parameter runwithoutscinds */

enum IndicatorDivingRoundingMode
{
   ROUNDING_CONSERVATIVE = 0,
   ROUNDING_AGGRESSIVE = 1
};
typedef enum IndicatorDivingRoundingMode INDICATORDIVINGROUNDINGMODE;

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
   SCIP_CONSHDLR*        indicatorconshdlr;  /**< indicator constraint handler */
   SCIP_CONSHDLR*        varboundconshdlr;   /**< varbound constraint handler */
   SCIP_HASHMAP*         scvars;             /**< hashmap to store semicontinuous variables */
   SCIP_HASHMAP*         indicatormap;       /**< hashmap to store indicator constraints of binary variables */
   SCIP_HASHMAP*         varboundmap;        /**< hashmap to store varbound constraints of binary variables */
   SCIP_Real             roundingfrac;       /**< in violation case all fractional below this value are fixed to constant */
   int                   roundingmode;       /**< decides which roundingmode is selected (0: conservative (default), 1: aggressive) */
   int                   semicontscoremode;  /**< which values of semi-continuous variables should get a high score? (0: low (default), 1: middle, 2: high) */
   SCIP_Bool             usevarbounds;       /**< should varbound constraints be considered? */
   SCIP_Bool             runwithoutscinds;   /**< should heur run if there are no indicator constraints modeling semicont. vars? */
   SCIP_Bool             gotoindconss;       /**< can we skip the candidate var until indicator conss handler determines the candidate var? */
   SCIP_Bool             containsviolindconss;/**< contains current solution violated indicator constraints? (only unbounded) */
   SCIP_Bool             newnode;            /**< are we at a new probing node? */
   int                   probingdepth;       /**< current probing depth */
};

/*
 * Local methods
 */

/** checks if constraint is violated but not fixed, i.e., it will be a diving candidate variable */
static
SCIP_Bool isViolatedAndNotFixed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< pointer to solution */
   SCIP_CONS*            cons                /**< pointer to indicator constraint */
   )
{
   SCIP_VAR* binvar;
   SCIP_Real solval;

   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "indicator") == 0);

   if( !SCIPisViolatedIndicator(scip, cons, sol) )
      return FALSE;

   binvar = SCIPgetBinaryVarIndicator(cons);
   solval = SCIPgetSolVal(scip, sol, binvar);

   return (SCIPisFeasIntegral(scip, solval) && SCIPvarGetLbLocal(binvar) < SCIPvarGetUbLocal(binvar) - 0.5);
}

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

/** checks if variable is indicator variable and stores corresponding indicator constraint; additionally, if we are at a
 *  new probing node, it checks whether there are violated but not fixed indicator constraints
 */
static
void checkAndGetIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             cand,               /**< candidate variable */
   SCIP_HASHMAP*         map,                /**< pointer to hashmap containing indicator conss */
   SCIP_CONS**           cons,               /**< pointer to store indicator constraint */
   SCIP_Bool*            isindicator,        /**< pointer to store whether candidate variable is indicator variable */
   SCIP_Bool*            containsviolindconss,/**< pointer to store whether there are violated and not fixed (unbounded) indicator constraints */
   SCIP_Bool             newnode,            /**< are we at a new probing node? */
   SCIP_SOL*             sol,                /**< pointer to solution */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(scip != NULL);
   assert(cand != NULL);
   assert(map != NULL);
   assert(cons != NULL);
   assert(isindicator != NULL);
   assert(sol != NULL);

   *cons = NULL;
   *isindicator = FALSE;

   *cons = (SCIP_CONS*) SCIPhashmapGetImage(map, cand);
   if( *cons != NULL )
      *isindicator = TRUE;

   /* since we are at a new probing node, check if there are violated and not fixed indicator constraints */
   if( newnode )
   {
      SCIP_CONS** indicatorconss;
      int nconss;
      int c;

      indicatorconss = SCIPconshdlrGetConss(conshdlr);
      nconss = SCIPconshdlrGetNActiveConss(conshdlr);
      *containsviolindconss = FALSE;

      for( c = 0; c < nconss; c++ )
      {
         *containsviolindconss = *containsviolindconss || isViolatedAndNotFixed(scip, sol, indicatorconss[c]);

         if( *containsviolindconss )
            break;
      }
   }
}

/** checks if variable is binary variable of varbound constraint and stores corresponding varbound constraint */
static
void checkAndGetVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             cand,               /**< candidate variable */
   SCIP_HASHMAP*         map,                /**< pointer to hashmap containing varbound conss */
   SCIP_CONS**           cons,               /**< pointer to store varbound constraint */
   SCIP_Bool*            isvarbound          /**< pointer to store whether candidate variable is indicator variable */
   )
{
   assert(scip != NULL);
   assert(cand != NULL);
   assert(map != NULL);
   assert(cons != NULL);
   assert(isvarbound != NULL);

   *cons = NULL;
   *isvarbound = FALSE;

   if( SCIPvarGetType(cand) != SCIP_VARTYPE_BINARY )
      return;

   *cons = (SCIP_CONS*) SCIPhashmapGetImage(map, cand);
   if( *cons != NULL )
      *isvarbound = TRUE;
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
      /* coverity[var_deref_op] */
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
 *  A variable x is semicontinuous if its bounds depend on at least one binary variable called the indicator,
 *  and indicator == 0 => x == x^0 for some real constant x^0.
 */
static
SCIP_RETCODE varIsSemicontinuous(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< the variable to check */
   SCIP_HASHMAP*         scvars,             /**< semicontinuous variable information */
   SCIP_Real             constant,           /**< value which should be equal to the constant */
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

      /* the 'off' domain of a semicontinuous var should reduce to a single point (constant) and be different from the 'on' domain */
      if( SCIPisEQ(scip, lb0, constant) && (!SCIPisEQ(scip, lb0, lb1) || !SCIPisEQ(scip, ub0, ub1)) )
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
      /* coverity[var_deref_op] */
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

      /* the 'off' domain of a semicontinuous var should reduce to a single point (constant) and be different from the 'on' domain */
      if( SCIPisEQ(scip, lb0, constant) && (!SCIPisEQ(scip, lb0, lb1) || !SCIPisEQ(scip, ub0, ub1)) )
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

/** checks if there are unfixed indicator variables modeling semicont. vars */
static
SCIP_RETCODE hasUnfixedSCIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< indicator constraint handler */
   SCIP_HASHMAP*         scvars,             /**< semicontinuous variable information */
   SCIP_Bool*            hasunfixedscindconss /**< pointer to store if there are unfixed indicator variables modeling semicont. vars */
   )
{
   SCIP_CONS** indicatorconss;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconss;
   int i;

   *hasunfixedscindconss = FALSE;
   indicatorconss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

   for( i = 0; i < nconss; i++ )
   {
      SCIP_VAR *binvar;
      SCIP_VAR* semicontinuousvar;
      SCIP_CONS* lincons;
      SCIP_Real rhs;
      int nconsvars;
      SCIP_Bool success;
      int v;

      binvar = SCIPgetBinaryVarIndicator(indicatorconss[i]);

      /* check if indicator variable is unfixed */
      if( SCIPvarGetLbLocal(binvar) < SCIPvarGetUbLocal(binvar) - 0.5 )
      {
         lincons = SCIPgetLinearConsIndicator(indicatorconss[i]);
         rhs = SCIPconsGetRhs(scip, lincons, &success);
         SCIP_CALL( SCIPgetConsNVars(scip, lincons, &nconsvars, &success) );

         /* check if constraint contains only two variables with finite rhs */
         /* TODO: allow also indicators for lower bounds */
         if( nconsvars == 2 && !SCIPisInfinity(scip, rhs) )
         {
            SCIP_CALL( SCIPgetConsVars(scip, lincons, consvars, nconsvars, &success) );
            SCIP_CALL( SCIPgetConsVals(scip, lincons, consvals, nconsvars, &success) );

            for( v = 0; v < nconsvars ; v++ )
            {
               if( consvars[v] == SCIPgetSlackVarIndicator(indicatorconss[i]) ) /* note that we have exact two variables */
                  continue;

               semicontinuousvar = consvars[v];
               SCIP_CALL( varIsSemicontinuous(scip, semicontinuousvar, scvars, rhs, &success) );

               /* check if semicontinuous variable */
               if( success )
               {
                  *hasunfixedscindconss = TRUE;
                  break;
               }
            }
            if( *hasunfixedscindconss )
               break;
         }
      }
   }
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);
   return SCIP_OKAY;
}

/** creates and initializes hashmaps
 *
 * indicatormap: binary var -> indicator constraint
 * varboundmap: binary var -> varbound constraint
 *
 * Currently exactly one constraint is assigned to a binary variable (per hashmap),
 * but a binary variable can also control more than one constraint.
 * TODO: Allow more than one corresponding indicator/varbound constraint per binary variable.
 */
static
SCIP_RETCODE createMaps(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        indicatorconshdlr,  /**< indicator constraint handler */
   SCIP_CONSHDLR*        varboundconshdlr,   /**< varbound constraint handler */
   SCIP_Bool             usevarbounds,       /**< should varbound constraints be considered? */
   SCIP_HASHMAP**        indicatormap,       /**< hashmap to store indicator constraints of binary variables */
   SCIP_HASHMAP**        varboundmap         /**< hashmap to store varbound constraints of binary variables */
   )
{
   SCIP_CONS** conss;
   int nconss;
   int i;

   assert(strcmp(SCIPconshdlrGetName(indicatorconshdlr), "indicator") == 0);
   assert(strcmp(SCIPconshdlrGetName(varboundconshdlr), "varbound") == 0);

   /* indicator constraints */
   nconss = SCIPconshdlrGetNConss(indicatorconshdlr);
   conss = SCIPconshdlrGetConss(indicatorconshdlr);
   SCIP_CALL( SCIPhashmapCreate(indicatormap, SCIPblkmem(scip), nconss) );
   for( i = 0; i < nconss; i++ )
   {
      if( !SCIPhashmapExists(*indicatormap, SCIPgetBinaryVarIndicator(conss[i])) )
      {
         SCIP_CALL( SCIPhashmapInsert(*indicatormap, SCIPgetBinaryVarIndicator(conss[i]), conss[i]) );
      }
   }

   /* varbound constraints */
   if( usevarbounds )
   {
      nconss = SCIPconshdlrGetNConss(varboundconshdlr);
      conss = SCIPconshdlrGetConss(varboundconshdlr);
      SCIP_CALL( SCIPhashmapCreate(varboundmap, SCIPblkmem(scip), nconss) );
      for( i = 0; i < nconss; i++ )
      {
         if( !SCIPhashmapExists(*varboundmap, SCIPgetVbdvarVarbound(scip, conss[i])) )
         {
            SCIP_CALL( SCIPhashmapInsert(*varboundmap, SCIPgetVbdvarVarbound(scip, conss[i]), conss[i]) );
         }
      }
   }
   return SCIP_OKAY;
}

#define MIN_RAND 1e-06
#define MAX_RAND 1e-05

/** calculate score and preferred rounding direction for the candidate variable */
static
void getScoreOfFarkasDiving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_VAR*             cand,               /**< candidate variable */
   SCIP_Real             candsfrac,          /**< fractional part of solution value of candidate variable */
   SCIP_Bool*            roundup,            /**< pointer to store whether the preferred rounding direction is upwards */
   SCIP_Real*            score               /**< pointer for diving score value */
   )
{
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_Real obj;

   randnumgen = SCIPdivesetGetRandnumgen(diveset);
   assert(randnumgen != NULL);

   obj = SCIPvarGetObj(cand);

   /* dive towards the pseudosolution, at the same time approximate the contribution to
    * a potential Farkas-proof (infeasibility proof) by y^TA_i = c_i.
    */
   if( SCIPisNegative(scip, obj) )
      *roundup = TRUE;
   else if( SCIPisPositive(scip, obj) )
      *roundup = FALSE;
   else
   {
      if( SCIPisEQ(scip, candsfrac, 0.5) )
         *roundup = !SCIPrandomGetInt(randnumgen, 0, 1);
      else
         *roundup = (candsfrac > 0.5);
   }

   /* larger score is better */
   *score = REALABS(obj) + SCIPrandomGetReal(randnumgen, MIN_RAND, MAX_RAND);

   /* prefer decisions on binary variables */
   if( SCIPvarGetType(cand) != SCIP_VARTYPE_BINARY )
      *score = -1.0 / *score;
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
SCIP_DECL_HEURFREE(heurFreeIndicatordiving)
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
SCIP_DECL_HEURINIT(heurInitIndicatordiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working data */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );
   SCIP_CALL( SCIPhashmapCreate( &heurdata->scvars, SCIPblkmem( scip ), SCIPgetNVars(scip) ));

   heurdata->indicatorconshdlr = SCIPfindConshdlr(scip, "indicator");
   heurdata->varboundconshdlr = SCIPfindConshdlr(scip, "varbound");

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitIndicatordiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working data */
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
   SCIP_Bool hasunfixedscindconss; /* are there unfixed indicator variables modeling a semicont. variable? */

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   assert(SCIPheurGetNDivesets(heur) > 0);
   assert(SCIPheurGetDivesets(heur) != NULL);
   diveset = SCIPheurGetDivesets(heur)[0];
   assert(diveset != NULL);

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   /* check if there are unfixed indicator variables modeling semicont. vars */
   SCIP_CALL( hasUnfixedSCIndicator(scip, heurdata->indicatorconshdlr, heurdata->scvars, &hasunfixedscindconss) );

   /* skip heuristic if problem doesn't contain unfixed indicator variables,
    * or if there are no varbound constraints which should be considered
    */
   if( !hasunfixedscindconss && (!heurdata->runwithoutscinds || !heurdata->usevarbounds || SCIPconshdlrGetNConss(heurdata->varboundconshdlr) == 0) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "call heurExecIndicatordiving at depth %d \n", SCIPgetDepth(scip));

   /* create and initialize hashmaps */
   SCIP_CALL( createMaps(scip, heurdata->indicatorconshdlr, heurdata->varboundconshdlr, heurdata->usevarbounds, &heurdata->indicatormap, &heurdata->varboundmap) );

   /* (re-)set flags */
   heurdata->gotoindconss = FALSE;
   heurdata->containsviolindconss = FALSE;
   heurdata->newnode = TRUE;
   heurdata->probingdepth = -1;

   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible, -1L, -1, -1.0, SCIP_DIVECONTEXT_SINGLE) );

   /* free hashmaps since constraints can get removed/modified till the next call */
   if( heurdata->usevarbounds )
      SCIPhashmapFree(&heurdata->varboundmap);
   SCIPhashmapFree(&heurdata->indicatormap);

   SCIPdebugMsg(scip, "leave heurExecIndicatordiving\n");

   return SCIP_OKAY;
}


/** calculate score and preferred rounding direction for the candidate variable */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScoreIndicatordiving)
{  /*lint --e{715}*/
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_VAR** consvars;
   SCIP_CONS* indicatorcons;
   SCIP_CONS* varboundcons;
   SCIP_CONS* lincons;
   SCIP_VAR* nonoptionvar; /* second variable in linear cons which is not the option variable (indicator: slackvar, varbound: binary var) */
   SCIP_VAR* semicontinuousvar;
   SCIP_Real lpsolsemicontinuous;
   SCVARDATA* scdata;
   SCIP_Real* consvals;
   SCIP_Real side;
   int nconsvars;
   int idxbvars; /* index of bounding variable in hashmap scdata */
   SCIP_Bool isindicatorvar;
   SCIP_Bool isvbdvar; /* variable bounding variable in varbound */
   SCIP_Bool issemicont; /* indicates whether variable has (maybe) required semicont. properties */
   SCIP_Bool fixconstant; /* should we fix the semicontinuous variable to its constant? */
   SCIP_Bool success;
   int v;
   int b;

   varboundcons = NULL;
   semicontinuousvar = NULL;
   scdata = NULL;
   lpsolsemicontinuous = 0.0;
   idxbvars = -1;
   isvbdvar = FALSE;
   issemicont = TRUE;

   heur = SCIPdivesetGetHeur(diveset);
   assert(heur != NULL);
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   randnumgen = SCIPdivesetGetRandnumgen(diveset);
   assert(randnumgen != NULL);

   /* check if we are at a new probing node; since diving heuristics backtrack at most one probing node, we are at a new
    * node iff the probing depth increased */
   if( heurdata->probingdepth < SCIPgetProbingDepth(scip) )
      heurdata->newnode = TRUE;
   else
   {
      assert(heurdata->probingdepth == SCIPgetProbingDepth(scip));
      heurdata->newnode = FALSE;
   }
   heurdata->probingdepth = SCIPgetProbingDepth(scip);

   /* skip if current candidate can not be determined by the indicator constraint handler and violated indicator
    * constraints still exists */
   if( !(SCIPisFeasIntegral(scip, candsol) && SCIPvarGetLbLocal(cand) < SCIPvarGetUbLocal(cand) - 0.5)
      && heurdata->gotoindconss )
   {
      *score = SCIP_REAL_MIN;
      *roundup = FALSE;
      return SCIP_OKAY;
   }
   else
      heurdata->gotoindconss = FALSE;

   /* check if candidate variable is indicator variable */
   checkAndGetIndicator(scip, cand, heurdata->indicatormap, &indicatorcons, &isindicatorvar,
      &heurdata->containsviolindconss, heurdata->newnode, heurdata->sol, heurdata->indicatorconshdlr);

   /* skip candidate in next calls since we have violated indicator constraints but current candidate is not determined
    * by the indicator constraint handler */
   if( heurdata->containsviolindconss &&
         !((SCIPisFeasIntegral(scip, candsol) && SCIPvarGetLbLocal(cand) < SCIPvarGetUbLocal(cand) - 0.5) && isindicatorvar) )
   {
      heurdata->gotoindconss = TRUE;
      *score = SCIP_REAL_MIN;
      *roundup = FALSE;
      return SCIP_OKAY;
   }

   /* check if candidate variable is bounding variable */
   if( heurdata->usevarbounds && !isindicatorvar )
   {
      checkAndGetVarbound(scip, cand, heurdata->varboundmap, &varboundcons, &isvbdvar);
   }

   /* Return
    * - if candidate variable is neither a indicator variable nor a variable bounding variable
    * - or if candidate variable is not an indicator variable but there will be indicator variables as candidates
    * - or if candidate variable is not an indicator variable and varbound constraints are not considered.
    */
   if( !isindicatorvar && (!isvbdvar || heurdata->containsviolindconss || !heurdata->usevarbounds) )
   {
      *score = SCIP_REAL_MIN;
      *roundup = FALSE;

      if( !heurdata->containsviolindconss && !isvbdvar )
      {
         getScoreOfFarkasDiving(scip, diveset, cand, candsfrac, roundup, score);
         *score = (*score / (100 + fabs(*score))) * 100 - 200; /* scale to [-300,-100] */
      }
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "cand: %s, candsol: %.2f, candobjcoeff: %f\n", SCIPvarGetName(cand), candsol, SCIPvarGetObj(cand));

   if( isindicatorvar ) /* prefer indicator constraint */
   {
      SCIP_Real rhs;

      lincons = SCIPgetLinearConsIndicator(indicatorcons);
      nonoptionvar = SCIPgetSlackVarIndicator(indicatorcons);
      rhs = SCIPconsGetRhs(scip, lincons, &success);
      issemicont = SCIPisInfinity(scip, -SCIPconsGetLhs(scip, lincons, &success)); /* TODO: allow also indicators for lower bounds */
      side = rhs;
   }
   else
   {
      SCIP_Real rhs;
      SCIP_Real lhs;

      assert(isvbdvar);

      lincons = varboundcons;
      nonoptionvar = SCIPgetVbdvarVarbound(scip, varboundcons);
      rhs = SCIPconsGetRhs(scip, lincons, &success);
      lhs = SCIPconsGetLhs(scip, lincons, &success);
      side = SCIPisInfinity(scip, rhs) ? lhs : rhs;
      assert(!SCIPisInfinity(scip, side));
   }
   SCIPdebugPrintCons(scip, lincons, NULL);

   SCIP_CALL( SCIPgetConsNVars(scip, lincons, &nconsvars, &success) );

   if( nconsvars != 2 || !issemicont )
   {
      getScoreOfFarkasDiving(scip, diveset, cand, candsfrac, roundup, score);
      *score = (*score / (100 + fabs(*score))) * 100 - 200; /* scale to [-300,-100] */
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
   SCIP_CALL( SCIPgetConsVars(scip, lincons, consvars, nconsvars, &success) );
   SCIP_CALL( SCIPgetConsVals(scip, lincons, consvals, nconsvars, &success) );

   issemicont = FALSE;
   for( v = 0; v < nconsvars ; v++ )
   {
      if( consvars[v] == nonoptionvar ) /* note that we have exact two variables */
         continue;

      semicontinuousvar = consvars[v];
      lpsolsemicontinuous = SCIPvarGetLPSol( semicontinuousvar );
      SCIPdebugMsg(scip, "%s lp sol %f %f\n", SCIPvarGetName( semicontinuousvar ), lpsolsemicontinuous,
                        consvals[v] );
      SCIP_CALL( varIsSemicontinuous(scip, semicontinuousvar, heurdata->scvars, side, &success) );

      /* only allow semicontinuous variables */
      if( success )
      {
         assert(SCIPhashmapExists(heurdata->scvars, (void*) semicontinuousvar));
         scdata = (SCVARDATA*) SCIPhashmapGetImage(heurdata->scvars, (void*) semicontinuousvar);
         assert(scdata != NULL);

         for( b = 0; b < scdata->nbnds; b++ )
         {
            if( (scdata->bvars[b] == cand || (SCIPvarIsNegated(cand) && scdata->bvars[0] == SCIPvarGetNegationVar(cand)))
                  && SCIPisEQ(scip, side, scdata->vals0[b]) )
            {

               /* TODO: handle also more general variables;
                * currently we handle only variables with domain vals0 < lb1 <= ub1 */
               if( SCIPisGE(scip, lpsolsemicontinuous, scdata->vals0[b]) && SCIPisLE(scip, lpsolsemicontinuous, scdata->ubs1[b]) )
               {
                  issemicont = TRUE;
                  idxbvars = b;
                  break;
               }
            }
         }
      }
   }

   /* only continue if semicontinuous variable */
   if( !issemicont )
   {
      getScoreOfFarkasDiving(scip, diveset, cand, candsfrac, roundup, score);
      *score = (*score / (100 + fabs(*score))) * 100 - 200; /* scale to [-300,-100] */
      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
      return SCIP_OKAY;
   }
   assert(idxbvars >= 0);
   assert(scdata != NULL);

   /* Case: Variable is in range [lb1,ub1] */
   if( SCIPisGE(scip, lpsolsemicontinuous, scdata->lbs1[idxbvars]) && SCIPisLE(scip, lpsolsemicontinuous, scdata->ubs1[idxbvars]))
   {
      *score = SCIPrandomGetReal(randnumgen, -1.0, 0.0);
      fixconstant = FALSE;
   }
   /* Case: Variable is equal to constant */
   else if( SCIPisEQ(scip, lpsolsemicontinuous, scdata->vals0[idxbvars]) )
   {
      *score = SCIPrandomGetReal(randnumgen, -1.0, 0.0);
      fixconstant = TRUE;
   }
   /* Case: Variable is between constant and lb1 */
   else
   {
      SCIP_Real shiftedlpsolsemicontinuous = lpsolsemicontinuous;
      SCIP_Real shiftedlbs1 = scdata->lbs1[idxbvars];

      assert(SCIPisGT(scip, lpsolsemicontinuous, scdata->vals0[idxbvars]) && SCIPisLT(scip, lpsolsemicontinuous, scdata->lbs1[idxbvars]));

      /* handle case if constant of semicont. var is not zero -> shift values */
      if( !SCIPisZero(scip, scdata->vals0[idxbvars]) )
      {
         shiftedlpsolsemicontinuous -= scdata->vals0[idxbvars];
         shiftedlbs1 -= scdata->vals0[idxbvars];
      }

      *score = 100 * (shiftedlbs1 - shiftedlpsolsemicontinuous) / shiftedlbs1;
      assert(*score>0);

      switch( (INDICATORDIVINGROUNDINGMODE)heurdata->roundingmode )
      {
      case ROUNDING_CONSERVATIVE:
         fixconstant = (*score > (1 - heurdata->roundingfrac) * 100);
         break;
      case ROUNDING_AGGRESSIVE:
         fixconstant = (*score <= (1 - heurdata->roundingfrac) * 100);
         break;
      default:
         return SCIP_INVALIDDATA;
      }

      switch( heurdata->semicontscoremode )
      {
      case 0:
         break;
      case 1:
         if( shiftedlpsolsemicontinuous < shiftedlbs1 * heurdata->roundingfrac )
            *score = 100 * (shiftedlpsolsemicontinuous / (heurdata->roundingfrac * shiftedlbs1));
         else
            *score = 100 * (-shiftedlpsolsemicontinuous / ((1 - heurdata->roundingfrac) * shiftedlbs1) + (1 / (1 - heurdata->roundingfrac)) );
         break;
      case 2:
         *score = 100 - *score;
         break;
      default:
         return SCIP_INVALIDDATA;
      }
      assert(*score>0);
   }

   /* Set roundup depending on whether we have an indicator constraint or a varbound constraint:
    * - indicator constraint: roundup == fix to constant
    * - varbound constraint: roundup == push to range
    */
   *roundup = isindicatorvar ? fixconstant : !fixconstant; /*lint !e644*/

   /* free memory */
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}


/** callback to check preconditions for diving, e.g., if an incumbent solution is available */
static
SCIP_DECL_DIVESETAVAILABLE(divesetAvailableIndicatordiving)
{
   /* Skip if problem doesn't contain indicator constraints.
    * If varbound constraints should be considered, skip only if there are also no varbound constraints.
    */
   *available =  SCIPconshdlrGetNActiveConss(SCIPfindConshdlr(scip, "indicator")) == 0;

   if( !*available )
   {
      SCIP_HEUR* heur;
      SCIP_HEURDATA* heurdata;

      heur = SCIPdivesetGetHeur(diveset);
      assert(heur != NULL);
      heurdata = SCIPheurGetData(heur);
      assert(heurdata != NULL);

      if( heurdata->runwithoutscinds && heurdata->usevarbounds )
      {
         *available = SCIPconshdlrGetNActiveConss(SCIPfindConshdlr(scip, "varbound")) == 0;
      }
   }

   return SCIP_OKAY;
}

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
         "in violation case all fractional below this value are fixed to constant",
         &heurdata->roundingfrac, FALSE, DEFAULT_ROUNDINGFRAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/roundingmode",
         "decides which roundingmode is selected (0: conservative, 1: aggressive)",
         &heurdata->roundingmode, FALSE, DEFAULT_ROUNDINGMODE, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/semicontscoremode",
         "which values of semi-continuous variables should get a high score? (0: low, 1: middle, 2: high)",
         &heurdata->semicontscoremode, FALSE, DEFAULT_SEMICONTSCOREMODE, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/usevarbounds",
         "should varbound constraints be considered?",
         &heurdata->usevarbounds, FALSE, DEFAULT_USEVARBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/runwithoutscinds",
         "should heur run if there are no indicator constraints modeling semicont. vars?",
         &heurdata->runwithoutscinds, FALSE, DEFAULT_RUNWITHOUTSCINDS, NULL, NULL) );

   return SCIP_OKAY;
}
