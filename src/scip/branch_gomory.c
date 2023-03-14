/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
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

/**@file   branch_gomory.c
 * @ingroup DEFPLUGINS_BRANCH
 * @brief  gomory cut branching rule
 * @author Mark Turner
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/branch_gomory.h"
#include "scip/pub_branch.h"
#include "scip/pub_var.h"
#include "scip/pub_lp.h"
#include "scip/scip_branch.h"
#include "scip/scip_cut.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"
#include "scip/scip_lp.h"
#include "scip_param.h"
#include "scip/type_cuts.h"
#include "scip/branch_relpscost.h"
#include <string.h>
#include <assert.h>



#define BRANCHRULE_NAME            "gomory"
#define BRANCHRULE_DESC            "gomory cut score branching"
#define BRANCHRULE_PRIORITY        10000000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define DEFAULT_MAXNCANDS          -1    /**< maximum number of branching candidates to produce a cut for */
#define DEFAULT_MINFRAC            0     /**< minimum fractionality of a variable in the LP that we branch on */
#define DEFAULT_EFFICACYWEIGHT     0.5   /**< weight of efficacy in the cut scoring rule */
#define DEFAULT_OBJPARALLELWEIGHT  0.5   /**< weight of objective parallelism in the cut scoring rule */
#define DEFAULT_PERFORMRELPSCOST   FALSE /**< if default branching rule without execution should be called */
#define DEFAULT_USEWEAKERCUTS      TRUE  /**< use weakener cuts derived from the exact branching split */


/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   int                   maxncands;             /**< maximum number of variable candidates to produce cut for */
   SCIP_Real             minfrac;               /**< minimum fractionality for which we will generate a GMI cut */
   SCIP_Real             efficacyweight;        /**< the weight of efficacy in weighted sum scoring rule */
   SCIP_Real             objparallelweight;     /**< the weight of objective parallelism in weighted sum scoring rule */
   SCIP_Bool             performrelpscost;      /**< if default branching rule without execution should be called */
   SCIP_Bool             useweakercuts;         /**< use weakener cuts derived from the exact branching split */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static
SCIP_Bool getGMIFromRow(
   SCIP* scip,
   int ncols,
   int nrows,
   SCIP_COL** cols,
   SCIP_ROW** rows,
   const SCIP_Real* binvrow,
   const SCIP_Real* binvarow,
   const SCIP_Real* lpval,
   SCIP_Real* cutcoefs,
   SCIP_Real* cutrhs,
   SCIP_Bool useweakerscuts
   )
{
   SCIP_COL* col;
   SCIP_ROW* row;
   SCIP_Real rowelem;
   SCIP_Real cutelem;
   SCIP_Real f0;
   SCIP_Real f0complementratio;
   SCIP_BASESTAT basestat;
   int i;
   int c;
   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   SCIP_Real rowrhs;
   SCIP_Real rowlhs;
   SCIP_Real rowact;
   SCIP_Real rowrhsslack;

   /* Clear the memory array of cut coefficients. It may store that of the last computed cut */
   BMSclearMemoryArray(cutcoefs, ncols);

   /* compute fractionality f0 and f0/(1-f0) */
   f0 = SCIPfeasFrac(scip, *lpval);
   f0complementratio = f0 / (1.0 - f0);

   /* The rhs of the cut is the fractional part of the LP solution of the basic variable */
   *cutrhs = -f0;

   /* Generate cut coefficient for the original variables */
   for ( c = 0; c < ncols; c++ )
   {
      col = cols[c];

      basestat = SCIPcolGetBasisStatus(col);
      /* Get simplex tableau coefficient */
      if ( basestat == SCIP_BASESTAT_LOWER )
      {
         /* Take coefficient if nonbasic at lower bound */
         rowelem = binvarow[c];
      }
      else if ( basestat == SCIP_BASESTAT_UPPER )
      {
         /* Flip coefficient if nonbasic at upper bound: x --> u - x */
         rowelem = -binvarow[c];
      }
      else
      {
         /* Nonbasic free variable at zero or basic variable. Just skip it. */
         continue;
      }

      /* Integer variables */
      if ( SCIPcolIsIntegral(col) && !useweakerscuts )
      {
         /* Warning: Because of numerics we can have cutelem < 0
          * In such a case it is very close to 0, so isZero will catch and we can ignore the coefficient */
         cutelem = SCIPfrac(scip, rowelem);
         if ( cutelem > f0 )
         {
            /* sum((1 - f_j) * f_0/(1 - f_0) x_j, j in J_I s.t. f_j > f_0 */
            cutelem = -((1.0 - cutelem) * f0complementratio);
         }
         else
         {
            /* sum(f_j * x_j, j in J_I s.t. f_j <= 0 */
            cutelem = -cutelem;
         }
      }
      /* Then continuous variables */
      else
      {
         if ( rowelem < 0 )
         {
            /* sum(a_j* f_0/(1 - f_0) x_j, j in J_C s.t. a_j < 0 */
            cutelem = rowelem * f0complementratio;
         }
         else
         {
            /* sum(a_j * x_j, j in J_C s.t. a_j >= 0 */
            cutelem = -rowelem;
         }
      }

      /* Cut is defined when variables are in [0, infinity). Translate to general bounds. */
      if ( !SCIPisZero(scip, cutelem) )
      {
         if ( basestat == SCIP_BASESTAT_UPPER )
         {
            cutelem = -cutelem;
            *cutrhs += cutelem * SCIPcolGetUb(col);
         }
         else
         {
            *cutrhs += cutelem * SCIPcolGetLb(col);
         }
         /* Add coefficient to cut */
         cutcoefs[SCIPcolGetLPPos(col)] = cutelem;
      }
   }

   /* generate cut coefficient for the slack variables. Skip the basic ones */
   for ( c = 0; c < nrows; c++ )
   {
      row = rows[c];
      assert( row != NULL );
      basestat = SCIProwGetBasisStatus(row);

      /* Get the simplex tableau coefficient */
      if ( basestat == SCIP_BASESTAT_LOWER )
      {
         /* Take coefficient if nonbasic at lower bound */
         rowelem = binvrow[SCIProwGetLPPos(row)];
         /* If there is a >= constraint or ranged constraint at the lower bound, we have to flip the row element */
         if ( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
            rowelem = -rowelem;
      }
      else if ( basestat == SCIP_BASESTAT_UPPER )
      {
         /* Take element if nonbasic at upper bound. Only non-positive slack variables can be nonbasic at upper,
          * therefore they should be flipped twice, meaning we can take the element directly */
         rowelem = binvrow[SCIProwGetLPPos(row)];
      }
      else
      {
      /* Nonbasic free variable at zero or basic variable. Just skip it. Free variable should not happen here*/
      assert( basestat == SCIP_BASESTAT_BASIC );
      continue;
      }

      /* Integer rows */
      if ( SCIProwIsIntegral(row) && !SCIProwIsModifiable(row) && !useweakerscuts )
      {
         /* Warning: Because of numerics we can have cutelem < 0
         * In such a case it is very close to 0, so isZero will catch and we can ignore the coefficient */
         cutelem = SCIPfrac(scip, rowelem);
         if ( cutelem > f0 )
         {
            /* sum((1 - f_j) * f_0/(1 - f_0) x_j, j in J_I s.t. f_j > f_0 */
            cutelem = -((1.0 - cutelem) * f0complementratio);
         }
         else
         {
            /* sum(f_j * x_j, j in J_I s.t. f_j <= 0 */
            cutelem = -cutelem;
         }
      }
      /* Then continuous variables */
      else
      {
         if ( rowelem < 0 )
         {
            /* sum(a_j* f_0/(1 - f_0) x_j, j in J_C s.t. a_j < 0 */
            cutelem = rowelem * f0complementratio;
         }
         else
         {
            /* sum(a_j * x_j, j in J_C s.t. a_j >= 0 */
            cutelem = -rowelem;
         }
      }
      /* Cut is defined in original variables, so we replace slack variables by their original definition */
      if ( !SCIPisZero(scip, cutelem) )
      {
         /* Coefficient is large enough so we keep it */
         rowlhs = SCIProwGetLhs(row);
         rowrhs = SCIProwGetRhs(row);
         assert( SCIPisLE(scip, rowlhs, rowrhs) );
         assert( !SCIPisInfinity(scip, rowlhs) && !SCIPisInfinity(scip, rowrhs) );

         /* If the slack variable is fixed we can ignore this cut coefficient */
         if ( SCIPisFeasZero(scip, rowrhs - rowlhs) )
            continue;

         /* Un-flip sack variable and adjust rhs if necessary.
          * Row at lower basis means the slack variable is at its upper bound.
          * Since SCIP adds +1 slacks, this can only happen when constraints have finite lhs */
         if ( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER )
         {
            assert( !SCIPisInfinity(scip, -rowlhs) );
            cutelem = -cutelem;
         }

         rowcols = SCIProwGetCols(row);
         rowvals = SCIProwGetVals(row);

         /* Eliminate slack variables. rowcols is sorted [columns in LP, columns not in LP] */
         for ( i = 0; i < SCIProwGetNLPNonz(row); i++ )
            cutcoefs[SCIPcolGetLPPos(rowcols[i])] -= cutelem * rowvals[i];
   
         /* Modify the rhs */
         rowact = SCIPgetRowActivity(scip, row);
         rowrhsslack = rowrhs - rowact;

         if ( SCIPisFeasZero(scip, rowrhsslack) )
            *cutrhs -= cutelem * (rowrhs - SCIProwGetConstant(row));
         else
         {
            assert( SCIPisFeasZero(scip, rowact - rowlhs) );
            *cutrhs -= cutelem * (rowlhs - SCIProwGetConstant(row));
         }

      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of branching rule
 */


/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyGomory)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleGomory(scip) );

   return SCIP_OKAY;
}


/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeGomory)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPfreeBlockMemoryNull(scip, &branchruledata);

   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpGomory)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   int maxncands;
   SCIP_COL** cols;
   SCIP_ROW** rows;
   int ncols;
   int nrows;
   SCIP_COL* col;
   int lppos;
   SCIP_Real* binvrow;
   SCIP_Real* binvarow;
   int* inds;
   int ninds;
   SCIP_Real cutrhs;
   SCIP_Real* cutcoefs;
   int* basisind;
   int* basicvarpos2tableaurow;
   SCIP_Real score;
   SCIP_Real bestscore;
   int bestcand;
   int i;
   SCIP_Bool success;
   SCIP_ROW* cut;
   const char* name;

   name = (char *) "test";

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "Execlp method of gomory branching\n");

   /* Get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, NULL, &nlpcands, NULL) );
   assert(nlpcands > 0);

   /* Initialise a DIDNOTRUN result */
   *result = SCIP_DIDNOTRUN;

   /* Get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* If the user has set maximum candidates to 0 or there are no candidates then return DIDNOTFIND  */
   if( branchruledata->maxncands == 0 || nlpcands == 0)
   {
      return SCIP_ERROR;
   }

   /* Compute the reliability pseudo-cost branching scores for the candidates */
   if ( branchruledata->performrelpscost )
   {
      /* We do not branch using this rule, but if enabled do take all the bound and conflict inferences made */
      SCIP_CALL( SCIPexecRelpscostBranching(scip, lpcands, lpcandssol, lpcandsfrac, nlpcands, FALSE,  result) );
      assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM);
   }

   /* Return SCIP_OKAY if relpscost has shown that this node can be cutoff or some variable domains have changed */
   if( *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM )
   {
      return SCIP_OKAY;
   }

   /* Get the maximum number of LP branching candidates that we generate cuts for and score */
   if( branchruledata->maxncands >= 0 )
   {
      maxncands = MIN(nlpcands, branchruledata->maxncands);
   }
   else
   {
      maxncands = nlpcands;
   }

   /* Get the Column and Row data */
   SCIP_CALL(SCIPgetLPColsData(scip, &cols, &ncols));
   SCIP_CALL(SCIPgetLPRowsData(scip, &rows, &nrows));

   /* Allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basicvarpos2tableaurow, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvarow, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nrows) );

   /* Create basis indices mapping (from the column position to LP tableau rox index) */
   for( i = 0; i < ncols; ++i )
   {
      basicvarpos2tableaurow[i] = -1;
   }
   SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );
   for( i = 0; i < nrows; ++i )
   {
      if( basisind[i] >= 0 )
         basicvarpos2tableaurow[basisind[i]] = i;
   }

   /* Initialise the best candidate */
   bestcand = 0;
   bestscore = -SCIPinfinity(scip);
   ninds = -1;
   
   /* Iterate over candidates and get best cut score */
   for( i = 0; i < maxncands; i++ ) {
      
      /* Initialise the score of the cut */
      score = 0;
      
      /* Get the LP position of the branching candidate */
      col = SCIPvarGetCol(lpcands[i]);
      lppos = SCIPcolGetLPPos(col);
      assert(lppos != -1);

      /* get the row of B^-1 for this basic integer variable with fractional solution value */
      SCIP_CALL(SCIPgetLPBInvRow(scip, basicvarpos2tableaurow[lppos], binvrow, inds, &ninds));
   
      /* Get the Tableau row for this basic integer variable with fractional solution value */
      SCIP_CALL(SCIPgetLPBInvARow(scip, basicvarpos2tableaurow[lppos], binvrow, binvarow, inds, &ninds));
   
      /* Compute the GMI cut */
      success = getGMIFromRow(scip, ncols, nrows, cols, rows, binvrow, binvarow, &lpcandssol[i], cutcoefs,
         &cutrhs, branchruledata->useweakercuts);
      
      /* Calculate the weighted sum score of measures */
      if ( success == SCIP_OKAY )
      {
         assert( SCIP_CALL( SCIPisCutEfficacious(scip, NULL, cut) ) );
         cut = NULL;
         SCIP_CALL( SCIPcreateRowUnspec(scip, &cut, name, ncols, cols, cutcoefs, -SCIPinfinity(scip), cutrhs, TRUE,
                   FALSE, TRUE) );
         if ( branchruledata-> efficacyweight != 0 )
            score += branchruledata->efficacyweight * SCIPgetCutEfficacy(scip, NULL, cut);
         if ( branchruledata->objparallelweight != 0 )
            score += branchruledata->objparallelweight * SCIPgetRowObjParallelism(scip, cut);
         SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   
         /* Replace the best cut if score is higher */
         if (score > bestscore) {
            bestscore = score;
            bestcand = i;
         }
      }
   }

   /* Free temporary memory */
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &binvrow);
   SCIPfreeBufferArray(scip, &binvarow);
   SCIPfreeBufferArray(scip, &basicvarpos2tableaurow);
   SCIPfreeBufferArray(scip, &basisind);
   SCIPfreeBufferArray(scip, &cutcoefs);

   SCIPdebugMsg(scip, " -> %d candidates, selected candidate %d: variable <%s> (frac=%g, factor=%g, score=%g)\n",
                nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandsfrac[bestcand],
                SCIPvarGetBranchFactor(lpcands[bestcand]), bestscore);

   /* Perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, lpcands[bestcand], NULL, NULL, NULL) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/*
 * branching rule specific interface methods
 */

/** creates the gomory cut branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleGomory(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
                                         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non-fundamental callbacks via specific setter functions*/
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyGomory) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeGomory) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpGomory) );

   /* gomory cut branching rule parameters */
   SCIP_CALL( SCIPaddIntParam(scip,"branching/gomory/maxncands",
         "maximum amount of branching candidates to generate gomory cut for (-1: all candidates)",
         &branchruledata->maxncands, FALSE, DEFAULT_MAXNCANDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,"branching/gomory/minfrac",
         "minimum fractionality of a variable in the LP for which we produce a cut",
         &branchruledata->minfrac, FALSE, DEFAULT_MINFRAC, 0, 0.5, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,"branching/gomory/efficacyweight",
         "weight of efficacy in the cut scoring rule",
         &branchruledata->efficacyweight, FALSE, DEFAULT_EFFICACYWEIGHT, -1, 1, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,"branching/gomory/objparallelweight",
         "weight of objective parallelism in the cut scoring rule",
         &branchruledata->objparallelweight, FALSE, DEFAULT_OBJPARALLELWEIGHT, -1, 1, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,"branching/gomory/performrlpscost",
         "whether default SCIP branching should be called without branching (used for bound inferences and conflicts)",
         &branchruledata->performrelpscost, FALSE, DEFAULT_PERFORMRELPSCOST, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,"branching/gomory/useweakercuts",
         "use weaker cuts that are exactly derived from the branching split disjunction",
         &branchruledata->useweakercuts, FALSE, DEFAULT_USEWEAKERCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
