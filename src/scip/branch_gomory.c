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
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/branch_gomory.h"
#include "scip/pub_branch.h"
#include "scip/pub_var.h"
#include "scip/pub_lp.h"
#include "scip/scip_branch.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"
#include "scip/scip_lp.h"
#include "scip_param.h"
#include "scip/type_cuts.h"
#include "scip/cuts.h"
#include <string.h>
#include <assert.h>



#define BRANCHRULE_NAME            "gomory"
#define BRANCHRULE_DESC            "gomory cut score branching"
#define BRANCHRULE_PRIORITY        10000000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define DEFAULT_MAXNCANDS          -1   /**< maximum number of branching candidates to produce a cut for */
#define DEFAULT_MINFRAC            0    /**< minimum fractionality of a variable in the LP for which we produce a cut */
#define DEFAULT_EFFICACYWEIGHT     0.5  /**< weight of efficacy in the cut scoring rule */
#define DEFAULT_OBJPARALLELWEIGHT  0.5  /**< weight of objective parallelism in the cut scoring rule */


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
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


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
   int nrows;
   int nvars;
   int ncols;
   SCIP_COL* col;
   int lppos;
   SCIP_AGGRROW* aggrrow;
   SCIP_Real* binvrow;
   int* inds;
   int ninds;
   SCIP_Real cutrhs;
   SCIP_Real* cutcoefs;
   int* basisind;
   int* basicvarpos2tableaurow;
   int* cutinds;
   int cutnnz;
   SCIP_Bool cutislocal;
   int cutrank;
   SCIP_Real cutefficacy;
   SCIP_Real score;
   SCIP_Real bestscore;
   int bestcand;
   int i;
   SCIP_Real feastol;
   SCIP_Bool success;


   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "Execlp method of gomory branching\n");

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, NULL, &nlpcands, NULL) );
   assert(nlpcands > 0);

   *result = SCIP_DIDNOTRUN;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* If the user has set maximum candidates to 0 or there are no candidates then return DIDNOTFIND  */
   if( branchruledata->maxncands == 0 || nlpcands == 0)
   {
      return *result;
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

   nvars = SCIPgetNVars(scip);
   nrows = SCIPgetNLPRows(scip);
   ncols = SCIPgetNLPCols(scip);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basicvarpos2tableaurow, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nrows) );
   SCIP_CALL( SCIPaggrRowCreate(scip, &aggrrow) );

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
   cutefficacy = 0;
   ninds = -1;
   feastol = SCIPgetLPFeastol(scip);
   /* Iterate over candidates and get best cut score */
   for( i = 0; i < maxncands; i++ ) {
      score = 0;
      /* Generate the GMI cut. We do not this using the textbook approach, but through the equivalent MIR cut */
      /* Get the LP position of the branching candidate */
      col = SCIPvarGetCol(lpcands[i]);
      lppos = SCIPcolGetLPPos(col);
      assert(lppos != -1);

      /* get the row of B^-1 for this basic integer variable with fractional solution value */
      SCIP_CALL(SCIPgetLPBInvRow(scip, basicvarpos2tableaurow[lppos], binvrow, inds, &ninds));

      /* Calculate the GMI cut */
      SCIP_CALL(SCIPaggrRowSumRows(scip, aggrrow, binvrow, inds, ninds, TRUE, TRUE, 2, nvars, &success));
      if (!success)
         continue;
      SCIP_CALL(SCIPcalcMIR(scip, NULL, TRUE, 0.9999, TRUE, TRUE, FALSE, NULL, NULL, feastol, 1 - feastol, 1.0,
         aggrrow, cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy, &cutrank, &cutislocal, &success));
      if (!success)
         continue;

      /* Calculate the weighted sum score of measures */
      score += branchruledata->efficacyweight * cutefficacy;
      /* Replace the best cut if score is higher */
      if (score > bestscore) {
         bestscore = score;
         bestcand = i;
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &binvrow);
   SCIPfreeBufferArray(scip, &basicvarpos2tableaurow);
   SCIPfreeBufferArray(scip, &basisind);
   SCIPfreeBufferArray(scip, &cutinds);
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPaggrRowFree(scip, &aggrrow);

   SCIPdebugMsg(scip, " -> %d candidates, selected candidate %d: variable <%s> (frac=%g, factor=%g, score=%g)\n",
                nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandsfrac[bestcand],
                SCIPvarGetBranchFactor(lpcands[bestcand]), bestscore);

   /* perform the branching */
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
         "maximum amount of branching candidates to generate gomory cut for",
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

   return SCIP_OKAY;
}
