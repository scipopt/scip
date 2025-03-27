/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   relax_benders.c
 * @ingroup OTHER_CFILES
 * @brief  benders relaxator
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/benders.h"
#include "scip/bendersdefcuts.h"
#include "scip/pub_relax.h"
#include "scip/scip.h"
#include "scip/relax_benders.h"
#include "scip/benders_default.h"
#include "scip/scip_benders.h"
#include "scip/scip_copy.h"
#include "scip/scip_message.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/type_benders.h"
#include "scip/type_message.h"
#include "scip/type_relax.h"
#include "scip/type_retcode.h"
#include "scip/type_stat.h"

#define RELAX_NAME             "benders"
#define RELAX_DESC             "applies Benders' decomposition and solves the problem"
#define RELAX_PRIORITY         1
#define RELAX_FREQ             1



/*
 * Data structures
 */
/* TODO: fill in the necessary relaxator data */

/** relaxator data */
struct SCIP_RelaxData
{
   SCIP*                 masterprob;         /**< the SCIP instance of the master problem */
   SCIP**                subproblems;        /**< an array of SCIP instances for the subproblems */
   SCIP_DECOMP*          decomp;             /**< the structure used for the decomposition */
   SCIP_HASHMAP*         mastervarmap;       /**< mapping between the original SCIP and the master problem */
   SCIP_HASHMAP**        subvarmaps;         /**< mapping between the original SCIP and the subproblems */
   int                   nsubproblems;       /**< the number of subproblems */
   SCIP_Bool             decompapplied;      /**< indicates whether the decomposition was applied */
   SCIP_VERBLEVEL        origverblevel;      /**< the verbosity level of the original SCIP instance */
};



/*
 * Local methods
 */

/** adding variables from the original problem to the respective decomposed problem */
static
SCIP_RETCODE addVariableToBendersProblem(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP*                 targetscip,         /**< the SCIP instance that the constraint will be copied into */
   SCIP_HASHMAP*         varmap,             /**< the variable hash map mapping the source variables to the target variables */
   SCIP_VAR*             sourcevar           /**< the variable to add to the problem */
   )
{
   assert(scip != NULL);
   assert(targetscip != NULL);
   assert(varmap != NULL);
   assert(sourcevar != NULL);

   /* if the variable is not in the hashmap, then it doesn't exist in the subproblem */
   if( !SCIPhashmapExists(varmap, sourcevar) )
   {
      SCIP_VAR* var;

      /* creating a variable as a copy of the original variable. */
      SCIP_CALL( SCIPcreateVarImpl(targetscip, &var, SCIPvarGetName(sourcevar), SCIPvarGetLbGlobal(sourcevar),
            SCIPvarGetUbGlobal(sourcevar), SCIPvarGetObj(sourcevar), SCIPvarGetType(sourcevar),
            SCIPvarGetImplType(sourcevar), SCIPvarIsInitial(sourcevar), SCIPvarIsRemovable(sourcevar),
            NULL, NULL, NULL, NULL, NULL) );

      /* adding the variable to the subproblem */
      SCIP_CALL( SCIPaddVar(targetscip, var) );

      /* adding the variable to the hash map so that it is copied correctly in the constraint */
      SCIP_CALL( SCIPhashmapInsert(varmap, sourcevar, var) );

      /* releasing the variable */
      SCIP_CALL( SCIPreleaseVar(targetscip, &var) );
   }

   return SCIP_OKAY;
}

/* when applying the decomposition, the constraints must be copied across from the original SCIP instance to the
 * respective Benders' problems in the relaxator. This could be either the master problem or one of the subproblems. The
 * process for performing this is the same for either problem type.
 */
static
SCIP_RETCODE addConstraintToBendersProblem(
   SCIP*                 scip,               /**< the original SCIP instance */
   SCIP*                 targetscip,         /**< the SCIP instance that the constraint will be copied into */
   SCIP_HASHMAP*         varmap,             /**< the variable hash map mapping the source variables to the target variables */
   SCIP_CONS*            sourcecons          /**< the constraint that being added to the subproblem */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR** consvars;
   int nconsvars;
   int i;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(targetscip != NULL);
   assert(varmap != NULL);
   assert(sourcecons != NULL);

   SCIPdebugMessage("Adding constraint <%s> to Benders' decomposition subproblem\n", SCIPconsGetName(sourcecons));

   /* getting the variables that are in the constraint */
   SCIP_CALL( SCIPgetConsNVars(scip, sourcecons, &nconsvars, &success) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );

   SCIP_CALL( SCIPgetConsVars(scip, sourcecons, consvars, nconsvars, &success) );
   assert(success);

   /* checking all variables to see whether they already exist in the subproblem. If they don't exist, then the variable
    * is created
    */
   for( i = 0; i < nconsvars; i++ )
   {
      SCIP_CALL( addVariableToBendersProblem(scip, targetscip, varmap, consvars[i]) );
   }

   /* freeing the buffer memory for the consvars */
   SCIPfreeBufferArray(scip, &consvars);

   /* copying the constraint from the master scip to the subproblem */
   SCIP_CALL( SCIPgetConsCopy(scip, targetscip, sourcecons, &cons, SCIPconsGetHdlr(sourcecons), varmap, NULL,
         SCIPconsGetName(sourcecons), SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons), SCIPconsIsMarkedPropagate(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons), TRUE, &success) );

   /* if the copy failed, then the subproblem for the decomposition could not be performed. */
   if( !success )
   {
      SCIPerrorMessage("It is not possible to copy constraint <%s>. Benders' decomposition could not be applied.\n",
         SCIPconsGetName(sourcecons));
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPaddCons(targetscip, cons) );
   SCIP_CALL( SCIPreleaseCons(targetscip, &cons) );

   return SCIP_OKAY;
}

/** Applies a Benders' decomposition to the problem based upon the decomposition selected from the storage */
static
SCIP_RETCODE applyDecomposition(
   SCIP*                 scip,               /**< the original SCIP instance */
   SCIP_RELAX*           relax,              /**< the relaxator */
   SCIP_DECOMP*          decomp              /**< the decomposition to apply to the problem */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   int* varslabels;
   int* conslabels;
   int nvars;
   int nconss;
   int nblocks;
   int i;
   char probname[SCIP_MAXSTRLEN];
   SCIP_Bool valid;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   SCIPdebugMessage("Applying a Benders' decomposition to <%s>\n", SCIPgetProbName(scip));

   /* storing the decomposition */
   relaxdata->decomp = decomp;

   /* retrieving the number of blocks for this decomposition */
   nblocks = SCIPdecompGetNBlocks(decomp);
   assert(nblocks > 0);

   /* initialising the subproblems for the Benders' decomposition */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->subproblems, nblocks) );
   relaxdata->nsubproblems = nblocks;

   /* creating the master problem */
   SCIP_CALL( SCIPcreate(&relaxdata->masterprob) );

   /* copying the plugins from the original SCIP instance to the master SCIP */
   SCIP_CALL( SCIPcopyPlugins(scip, relaxdata->masterprob, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, &valid) );

   /* including the default Benders' decomposition plugin. This is added separately so that any addition of the plugin
    * in the master problem is ignored
    */
   SCIP_CALL( SCIPincludeBendersDefault(relaxdata->masterprob) );

   /* copying the parameter settings from the original SCIP to the master problem */
   SCIP_CALL( SCIPcopyParamSettings(scip, relaxdata->masterprob) );

   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "master_%s", SCIPgetProbName(scip), i);
   SCIP_CALL( SCIPcreateProbBasic(relaxdata->masterprob, probname) );

   /* creating the subproblems before adding the constraints */
   for( i = 0; i < nblocks; i++ )
   {
      SCIP_CALL( SCIPcreate(&relaxdata->subproblems[i]) );

      /* copying the plugins from the original SCIP instance to the subproblem SCIP */
      SCIP_CALL( SCIPcopyPlugins(scip, relaxdata->subproblems[i], TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
            TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, &valid) );

      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "sub_%s_%d", SCIPgetProbName(scip), i);
      SCIP_CALL( SCIPcreateProbBasic(relaxdata->subproblems[i], probname) );
   }

   /* TODO: Need to work out whether a check for original and transformed problem is necessary */

   /* getting the variables and constraints from the problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   /* allocating buffer memory for the labels arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );

   /* getting the labels for the variables and constraints from the decomposition */
   SCIPdecompGetVarsLabels(decomp, vars, varslabels, nvars);
   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);

   /* creating the variable map for the master problem */
   SCIP_CALL( SCIPhashmapCreate(&relaxdata->mastervarmap, SCIPblkmem(relaxdata->masterprob), nvars) );

   /* creating the variable maps for adding the constraints to the subproblems */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->subvarmaps, nblocks) );

   for( i = 0; i < nblocks; i++ )
   {
      SCIP_CALL( SCIPhashmapCreate(&relaxdata->subvarmaps[i], SCIPblkmem(relaxdata->subproblems[i]), nvars) );
   }

   /* copying the constraints to the appropriate subproblems */
   for( i = 0; i < nconss; i++ )
   {
      /* the constraints with a block label >= 0 correspond to subproblem constraints. All other constraints are master
       * constraints.
       */
      if( conslabels[i] >= 0 )
      {
         assert(conslabels[i] < relaxdata->nsubproblems);
         SCIP_CALL( addConstraintToBendersProblem(scip, relaxdata->subproblems[conslabels[i]], relaxdata->subvarmaps[conslabels[i]],
               conss[i]) );
      }
      else
      {
         SCIP_CALL( addConstraintToBendersProblem(scip, relaxdata->masterprob, relaxdata->mastervarmap, conss[i]) );
      }
   }

   /* copying all of the variables from the original problem to the respective decomposed problem */
   for( i = 0; i < nvars; i++ )
   {
      if( varslabels[i] >= 0 )
      {
         assert(varslabels[i] < relaxdata->nsubproblems);
         SCIP_CALL( addVariableToBendersProblem(scip, relaxdata->subproblems[varslabels[i]],
               relaxdata->subvarmaps[varslabels[i]], vars[i]) );
      }
      else
      {
         SCIP_CALL( addVariableToBendersProblem(scip, relaxdata->masterprob, relaxdata->mastervarmap, vars[i]) );
      }
   }

   /* creating the Benders' decomposition my calling the default plugin */
   SCIP_CALL( SCIPcreateBendersDefault(relaxdata->masterprob, relaxdata->subproblems, nblocks) );

   /* activating the Benders' constraint handler for the scenario stages.
    * TODO: consider whether the two-phase method should be activated by default in the scenario stages.
    */
   SCIP_CALL( SCIPsetBoolParam(relaxdata->masterprob, "constraints/benders/active", TRUE) );
   SCIP_CALL( SCIPsetBoolParam(relaxdata->masterprob, "constraints/benderslp/active", TRUE) );

   /* changing settings that are required for Benders' decomposition */
   SCIP_CALL( SCIPsetPresolving(relaxdata->masterprob, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(relaxdata->masterprob, "propagating/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(relaxdata->masterprob, "propagating/maxroundsroot", 0) );

   /* disabling aggregation since it can affect the mapping between the master and subproblem variables */
   SCIP_CALL( SCIPsetBoolParam(relaxdata->masterprob, "presolving/donotaggr", TRUE) );
   SCIP_CALL( SCIPsetBoolParam(relaxdata->masterprob, "presolving/donotmultaggr", TRUE) );

   SCIPfreeBufferArray(scip, &conslabels);
   SCIPfreeBufferArray(scip, &varslabels);

   relaxdata->decompapplied = TRUE;

   return SCIP_OKAY;
}

/** sets the verbosity level of the decomposed problem and disables verbosity of the original SCIP */
static
SCIP_RETCODE setVerbosityLevel(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_RELAX*           relax               /**< the relaxator */
)
{
   SCIP_RELAXDATA* relaxdata;
   int verblevel;

   assert(scip != NULL);
   assert(relax != NULL);

   /* updating the verbosity level of the  */
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* getting the original verbosity level so that this can be reset at the end of the Benders' solve */
   SCIP_CALL( SCIPgetIntParam(scip, "display/verblevel", &verblevel) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
   relaxdata->origverblevel = (SCIP_VERBLEVEL)verblevel;

   /* setting the master problem verbosity level to the original SCIP verbosity level */
   SCIP_CALL( SCIPsetIntParam(relaxdata->masterprob, "display/verblevel", relaxdata->origverblevel) );

   return SCIP_OKAY;
}


/** resets the verbosity level of the original SCIP instance */
static
SCIP_RETCODE resetVerbosityLevel(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_RELAX*           relax               /**< the relaxator */
)
{
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL);
   assert(relax != NULL);

   /* updating the verbosity level of the  */
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* getting the original verbosity level so that this can be reset at the end of the Benders' solve */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", relaxdata->origverblevel) );

   return SCIP_OKAY;
}

/** returns whether a solution exists for the master problem */
static
SCIP_Bool masterSolutionExists(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_RELAX*           relax               /**< the relaxator */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   assert(relaxdata->masterprob != NULL);

   return (SCIPgetBestSol(relaxdata->masterprob) != NULL);
}


/** solves the Benders' decomposition subproblems using the best solution from the master problem */
static
SCIP_RETCODE solveBendersSubproblems(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_RELAX*           relax,              /**< the relaxator */
   SCIP_Bool*            infeasible          /**< indicates whether the best solution is infeasible */
)
{
   SCIP_RELAXDATA* relaxdata;
   SCIP* masterprob;
   SCIP_BENDERS* benders;
   SCIP_SOL* bestsol;
   int i;

   assert(scip != NULL);
   assert(relax != NULL);
   assert(infeasible != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   masterprob = relaxdata->masterprob;

   assert(SCIPgetNActiveBenders(masterprob) == 1);

   /* getting the Default Benders' plugin to solve the subproblems */
   benders = SCIPfindBenders(masterprob, "default");
   assert(benders != NULL);
   assert(relaxdata->nsubproblems == SCIPbendersGetNSubproblems(benders));

   /* getting the best solution for the master problem */
   bestsol = SCIPgetBestSol(masterprob);

   /* solving the Benders' decomposition subproblems */
   for( i = 0; i < relaxdata->nsubproblems; i++)
   {
      SCIP* subproblem = SCIPbendersSubproblem(benders, i);
      assert(subproblem != NULL && SCIPgetStage(subproblem) >= SCIP_STAGE_PROBLEM);

      /* setting up the subproblem with the best solution from the master problem */
      SCIP_CALL( SCIPsetupBendersSubproblem(masterprob, benders, bestsol, i, SCIP_BENDERSENFOTYPE_CHECK) );

      /* solving the subproblem */
      SCIP_CALL( SCIPsolveBendersSubproblem(masterprob, benders, bestsol, i, infeasible, TRUE, NULL) );

      /* if a subproblem is infeasible, then the solution is reported as infeasible */
      if( (*infeasible) )
         break;
   }

   return SCIP_OKAY;
}

/** constructs the NLP solution and returns it */
static
SCIP_RETCODE getNlpSolution(
   SCIP*                 scip,               /**< the SCIP instance for which an NLP must be constructed */
   SCIP_SOL**            nlpsol              /**< the NLP solution that will be created */
   )
{
   assert(scip != NULL);
   assert(SCIPisNLPConstructed(scip));
   assert(SCIPgetNNlpis(scip));

   SCIP_CALL( SCIPcreateNLPSol(scip, nlpsol, NULL) );

   return SCIP_OKAY;
}

/** using the stored mappings for the variables, the solution values from the master and subproblems are used to set the
 * solution values in the original SCIP solution
 */
static
SCIP_RETCODE setSolutionValues(
   SCIP*                 scip,               /**< the original SCIP instance */
   SCIP_RELAX*           relax,              /**< the relaxator */
   SCIP_SOL*             sol                 /**< the solution for the original SCIP instance */
)
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_DECOMP* decomp;
   SCIP_VAR** vars;
   int* varslabels;
   int nvars;
   int nblocks;
   int i;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   decomp = relaxdata->decomp;

   /* retrieving the number of blocks for this decomposition */
   nblocks = SCIPdecompGetNBlocks(decomp);
   assert(nblocks > 0);

   /* getting the variables and constraints from the problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* allocating buffer memory for the labels arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, nvars) );

   /* getting the labels for the variables and constraints from the decomposition */
   SCIPdecompGetVarsLabels(decomp, vars, varslabels, nvars);

   /* looping over all variables and assigning the value to the current master or subproblem solution */
   for( i = 0; i < nvars; i++ )
   {
      SCIP* solscip;
      SCIP_HASHMAP* varmap;
      SCIP_SOL* bestsol;
      SCIP_Bool subproblem;

      /* if the varlabel is >= 0, then the variable belongs to the subproblem. Otherwise, the variable belongs to the
       * master problem.
       */
      subproblem = varslabels[i] >= 0;
      if( subproblem )
      {
         solscip = relaxdata->subproblems[varslabels[i]];
         varmap = relaxdata->subvarmaps[varslabels[i]];
      }
      else
      {
         solscip = relaxdata->masterprob;
         varmap = relaxdata->mastervarmap;
      }

      if( SCIPhashmapExists(varmap, vars[i]) )
      {
          SCIP_Bool nlprelaxation;

          /* the subproblem could be an NLP. As such, we need to get the solution directly from the NLP */
          nlprelaxation = subproblem && SCIPisNLPConstructed(solscip) && SCIPgetNNlpis(solscip);
          if( nlprelaxation )
          {
             SCIP_CALL( getNlpSolution(solscip, &bestsol) );
          }
          else
             bestsol = SCIPgetBestSol(solscip);

          SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i],
             SCIPgetSolVal(solscip, bestsol, (SCIP_VAR*)SCIPhashmapGetImage(varmap, vars[i]))) );

          /* if the solution is from an NLP, then we need to free it */
          if( nlprelaxation )
          {
             SCIP_CALL( SCIPfreeSol(solscip, &bestsol) );
          }
      }
   }

   SCIPfreeBufferArray(scip, &varslabels);

   return SCIP_OKAY;
}

/** frees the subproblems after the solve */
static
SCIP_RETCODE freeBendersSubproblems(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_RELAX*           relax               /**< the relaxator */
)
{
   SCIP_RELAXDATA* relaxdata;
   SCIP* masterprob;
   SCIP_BENDERS* benders;
   int i;

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   masterprob = relaxdata->masterprob;
   assert(SCIPgetNActiveBenders(masterprob) == 1);

   /* getting the Default Benders' plugin to solve the subproblems */
   benders = SCIPfindBenders(masterprob, "default");
   assert(benders != NULL);
   assert(relaxdata->nsubproblems == SCIPbendersGetNSubproblems(benders));

   /* freeing the Benders' decomposition subproblems */
   for( i = 0; i < relaxdata->nsubproblems; i++)
   {
      SCIP_CALL( SCIPfreeBendersSubproblem(scip, benders, i) );
   }

   return SCIP_OKAY;
}


/** creates a solution for the original SCIP instance using the solution from the decomposed problem */
static
SCIP_RETCODE createOriginalSolution(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_RELAX*           relax,              /**< the relaxator */
   SCIP_Bool*            infeasible          /**< indicates whether the best solution is infeasible */
)
{
   SCIP_SOL* sol;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(relax != NULL);
   assert(infeasible != NULL);

   /* if there is no master solution, then no solution is copied across to the original SCIP instance */
   if( !masterSolutionExists(scip, relax) )
      return SCIP_OKAY;

   /* creating the solution for the original SCIP */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   /* solving the Benders' decomposition subproblems with the best master problem solution */
   SCIP_CALL( solveBendersSubproblems(scip, relax, infeasible) );

   /* setting the solution values of the solution based on the decomposed problem solution */
   SCIP_CALL( setSolutionValues(scip, relax, sol) );

   /* checking if the solution is feasible for the original problem */
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
   if( success )
   {
      SCIP_Bool stored;
      SCIP_CALL( SCIPaddSol(scip, sol, &stored) );
   }

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   /* freeing the subproblems after the solve */
   SCIP_CALL( freeBendersSubproblems(scip, relax) );

   return SCIP_OKAY;
}



/*
 * Callback methods of relaxator
 */

#define relaxCopyBenders NULL
#define relaxInitBenders NULL
#define relaxExitBenders NULL


/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_RELAXFREE(relaxFreeBenders)
{  /*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);

   /* it is possible that the relaxation data is not created. In this case, nothing needs to be freed */
   if( relaxdata == NULL )
      return SCIP_OKAY;

   /* if the decomposition has been applied, then the corresponding SCIP instances need to be freed */
   if( relaxdata->decompapplied )
   {
      int i;

      /* freeing the variable hash maps */
      for( i = relaxdata->nsubproblems - 1; i >= 0; i-- )
      {
         SCIPhashmapFree(&relaxdata->subvarmaps[i]);
      }
      SCIPhashmapFree(&relaxdata->mastervarmap);

      SCIPfree(&relaxdata->masterprob);

      for( i = relaxdata->nsubproblems - 1; i >= 0; i-- )
      {
         SCIPfree(&relaxdata->subproblems[i]);
      }

      /* freeing the allocated arrays */
      SCIPfreeBlockMemoryArray(scip, &relaxdata->subvarmaps, relaxdata->nsubproblems);
      SCIPfreeBlockMemoryArray(scip, &relaxdata->subproblems, relaxdata->nsubproblems);
   }

   SCIPfreeBlockMemory(scip, &relaxdata);

   return SCIP_OKAY;
}




/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
static
SCIP_DECL_RELAXINITSOL(relaxInitsolBenders)
{  /*lint --e{715}*/
   SCIP_BENDERS* benders;
   SCIP_DECOMP** decomps;
   int ndecomps;
   SCIP_Bool applybenders;

   assert(scip != NULL);

   SCIP_CALL( SCIPgetBoolParam(scip, "decomposition/applybenders", &applybenders) );
   if( !applybenders )
      return SCIP_OKAY;

   SCIPgetDecomps(scip, &decomps, &ndecomps, FALSE);

   /* if there are no decompositions, then Benders' decomposition can't be applied */
   if( ndecomps == 0 )
      return SCIP_OKAY;

   assert(decomps[0] != NULL);

   /* if there already exists an active Benders' decomposition, then default decomposition is not applied. */
   if( SCIPgetNActiveBenders(scip) > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "A Benders' decomposition already exists. The default Benders' decomposition will not be applied to the stored decomposition.\n");
      return SCIP_OKAY;
   }

   /* retrieving the default Benders' decomposition plugin */
   benders = SCIPfindBenders(scip, "default");

   /* if the default Benders' decomposition plugin doesn't exist, then this will result in an error */
   if( benders == NULL )
   {
      SCIPerrorMessage("The default Benders' decomposition plugin is required to apply Benders' decomposition using the input decomposition.");
      return SCIP_ERROR;
   }

   /* applying the Benders' decomposition. If SCIP is in the PROBLEM stage, then the auxiliary variables don't need to
    * be added. However, in any other stage, then the auxiliary variables must be added to the problem.
    */
   SCIP_CALL( applyDecomposition(scip, relax, decomps[0]) );

   /* setting the verbosity level of the original SCIP instance and the master problem */
   SCIP_CALL( setVerbosityLevel(scip, relax) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
static
SCIP_DECL_RELAXEXITSOL(relaxExitsolBenders)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(relax != NULL);

   SCIP_CALL( resetVerbosityLevel(scip, relax) );

   return SCIP_OKAY;
}

/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecBenders)
{  /*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata;
   SCIP_Bool success;
   SCIP_Bool infeasible;

   assert(scip != NULL);
   assert(relax != NULL);

   (*result) = SCIP_DIDNOTRUN;

   relaxdata = SCIPrelaxGetData(relax);

   if( relaxdata == NULL)
      return SCIP_OKAY;

   /* the relaxator is only executed if the Benders decomposition is applied */
   if( !relaxdata->decompapplied )
      return SCIP_OKAY;

   /* we only run at the root node. If the decomposition was not able to solve the problem at the root node, then we
    * fall back to solving the MIP directly.
    */
   if( SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   /* checking whether there is enough time and memory to perform the decomposition solve */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );
   if( !success )
      return SCIP_OKAY;

   SCIPverbMessage(relaxdata->masterprob, SCIP_VERBLEVEL_NORMAL, NULL,
      "\nApplying Benders' decomposition and solving the decomposed problem.\n\n");

   /* copying the time and memory limits from the original SCIP to the master problem */
   SCIP_CALL( SCIPcopyLimits(scip, relaxdata->masterprob) );

   SCIP_CALL( SCIPcopySolvingTime(scip, relaxdata->masterprob) );

   SCIP_CALL( SCIPsolve(relaxdata->masterprob) );

   /* if the problem is solved to be infeasible, then the result needs to be set to CUTOFF. */
   if( SCIPgetStatus(relaxdata->masterprob) == SCIP_STATUS_INFEASIBLE )
   {
      (*result) = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* a solution for the original SCIP needs to be created. This will transfer the best found solution from the
    * Benders decomposition solve back to the original SCIP instance.
    */
   SCIP_CALL( createOriginalSolution(scip, relax, &infeasible) );

   if( infeasible )
   {
      (*result) = SCIP_CUTOFF;
   }
   else
   {
      (*lowerbound) = SCIPgetDualbound(relaxdata->masterprob);
      (*result) = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}





/*
 * relaxator specific interface methods
 */

/** creates the benders relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxBenders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;

   /* create benders relaxator data */
   relaxdata = NULL;
   SCIP_CALL( SCIPallocBlockMemory(scip, &relaxdata) );
   BMSclearMemory(relaxdata);

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelax(scip, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxCopyBenders, relaxFreeBenders, relaxInitBenders, relaxExitBenders, relaxInitsolBenders,
         relaxExitsolBenders, relaxExecBenders, relaxdata) );

   return SCIP_OKAY;
}

/** returns the master problem SCIP instance */
SCIP_RETCODE SCIPrelaxBendersPrintStatistics(
   SCIP_RELAX*           relax               /**< the Benders' decomposition relaxator */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   /* if the decomposition has not been applied, then there are no statistics to display */
   if( !relaxdata->decompapplied )
      return SCIP_OKAY;

   assert(relaxdata->masterprob != NULL);

   SCIP_CALL( SCIPprintStatistics(relaxdata->masterprob, NULL) );

   return SCIP_OKAY;
}
