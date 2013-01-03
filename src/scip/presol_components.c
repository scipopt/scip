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

/**@file   presol_components.c
 * @ingroup PRESOLVERS
 * @brief  solve independent components in advance
 * @author Dieter Weninger
 * @author Gerald Gamrath
 *
 * This presolver looks for independent components at the end of the presolving.
 * If independent components are found in which a maximum number of discrete variables
 * is not exceeded, the presolver tries to solve them in advance as subproblems.
 * Afterwards, if a subproblem was solved to optimality, the corresponding
 * variables/constraints can be fixed/deleted in the main problem.
 *
 * @todo simulation of presolving without solve
 * @todo solve all components with less than given size, count number of components with nodelimit reached;
 *       if all components could be solved within nodelimit (or all but x), continue solving components in
 *       increasing order until one hit the node limit
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/presol_components.h"

#define PRESOL_NAME            "components"
#define PRESOL_DESC            "components presolver"
#define PRESOL_PRIORITY        -9200000      /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS             -1      /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY               TRUE      /**< should presolver be delayed, if other presolvers found reductions? */

#define DEFAULT_WRITEPROBLEMS     FALSE      /**< should the single components be written as an .lp-file? */
#define DEFAULT_MAXINTVARS          500      /**< maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving) */
#define DEFAULT_NODELIMIT       10000LL      /**< maximum number of nodes to be solved in subproblems */
#define DEFAULT_INTFACTOR           1.0      /**< the weight of an integer variable compared to binary variables */
#define DEFAULT_RELDECREASE         0.2      /**< percentage by which the number of variables has to be decreased after the last component solving
                                              *   to allow running again (1.0: do not run again) */

#ifdef SCIP_STATISTIC
static int NCATEGORIES = 6;
static int CATLIMITS[] = {0,20,50,100,500};
#endif

/*
 * Data structures
 */

/** control parameters */
struct SCIP_PresolData
{
   SCIP_Bool             didsearch;          /** did the presolver already search for components? */
   SCIP_Bool             pluginscopied;      /** was the copying of the plugins successful? */
   SCIP_Bool             writeproblems;      /** should the single components be written as an .lp-file? */
   int                   maxintvars;         /** maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving) */
   SCIP_Longint          nodelimit;          /** maximum number of nodes to be solved in subproblems */
   SCIP_Real             intfactor;          /** the weight of an integer variable compared to binary variables */
   SCIP_Real             reldecrease;        /** percentage by which the number of variables has to be decreased after the last component solving
                                              *  to allow running again (1.0: do not run again) */
   int                   lastnvars;          /** number of variables after last run of the presolver */
   SCIP*                 subscip;            /** sub-SCIP used to solve single components */
#ifdef SCIP_STATISTIC
   int*                  compspercat;        /** number of components of the different categories */
   int                   nsinglevars;        /** number of components with a single variable without constraint */
   SCIP_Real             subsolvetime;       /** total solving time of the subproblems */
#endif
};

/*
 * Statistic methods
 */

#ifdef SCIP_STATISTIC
/** initialize data for statistics */
static
SCIP_RETCODE initStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   assert(scip != NULL);
   assert(presoldata != NULL);

   SCIP_CALL( SCIPallocMemoryArray(scip, &presoldata->compspercat, NCATEGORIES) );
   BMSclearMemoryArray(presoldata->compspercat, NCATEGORIES);

   presoldata->nsinglevars = 0;
   presoldata->subsolvetime = 0.0;

   return SCIP_OKAY;
}

/** free data for statistics */
static
void freeStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   assert(scip != NULL);
   assert(presoldata != NULL);

   SCIPfreeMemoryArray(scip, &presoldata->compspercat);
}

/** reset data for statistics */
static
void resetStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   assert(scip != NULL);
   assert(presoldata != NULL);

   BMSclearMemoryArray(presoldata->compspercat, NCATEGORIES);

   presoldata->nsinglevars = 0;
   presoldata->subsolvetime = 0.0;
}


/** statistics: categorize the component with the given number of binary and integer variables */
static
void updateStatisticsComp(
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   int                   nbinvars,           /**< number of binary variables */
   int                   nintvars            /**< number of integer variables */
   )
{
   int ndiscretevars;
   int i;

   assert(presoldata != NULL);

   ndiscretevars = nbinvars + nintvars;

   /* check into which category the component belongs by looking at the number of discrete variables */
   for( i = 0; i < (NCATEGORIES - 1); ++i )
   {
      if( ndiscretevars <= CATLIMITS[i] )
      {
         presoldata->compspercat[i]++;
         break;
      }
   }

   /* number of discrete variables greater than all limits, so component belongs to last category */
   if( i == (NCATEGORIES - 1) )
      presoldata->compspercat[i]++;
}

/** statistics: increase the number of components with a single variable and no constraints */
static
void updateStatisticsSingleVar(
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   assert(presoldata != NULL);

   presoldata->nsinglevars++;
}

/** statistics: update the total subproblem solving time */
static
void updateStatisticsSubsolvetime(
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_Real             subsolvetime        /**< subproblem solving time to add to the statistics */
   )
{
   assert(presoldata != NULL);

   presoldata->subsolvetime += subsolvetime;
}

/** print statistics */
static
void printStatistics(
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   int i;

   assert(presoldata != NULL);

   printf("############\n");
   printf("# Connected Components Presolver Statistics:\n");

   printf("# Categorization:");
   for( i = 0; i < NCATEGORIES - 1; ++i )
   {
      printf("[<= %d: %d]", CATLIMITS[i], presoldata->compspercat[i]);
   }
   printf("[> %d: %d]\n", CATLIMITS[NCATEGORIES - 2], presoldata->compspercat[NCATEGORIES - 1]);
   printf("# Components without constraints: %d\n", presoldata->nsinglevars);
   printf("# Total subproblem solving time: %.2f\n", presoldata->subsolvetime);
   printf("############\n");
}
#endif


/*
 * Local methods
 */

/** copies a connected component consisting of the given constraints and variables into a sub-SCIP
 *  and tries to solve the sub-SCIP to optimality
 */
static
SCIP_RETCODE copyAndSolveComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_HASHMAP*         consmap,            /**< constraint hashmap used to improve performance */
   int                   compnr,             /**< number of the component */
   SCIP_CONS**           conss,              /**< constraints contained in this component */
   int                   nconss,             /**< number of constraints contained in this component */
   SCIP_VAR**            vars,               /**< variables contained in this component */
   SCIP_Real*            fixvals,            /**< array to store the values to fix the variables to */
   int                   nvars,              /**< number of variables contained in this component */
   int                   nbinvars,           /**< number of binary variables contained in this component */
   int                   nintvars,           /**< number of integer variables contained in this component */
   int*                  nsolvedprobs,       /**< pointer to increase, if the subproblem was solved */
   int*                  ndeletedvars,       /**< pointer to increase by the number of deleted variables */
   int*                  ndeletedconss,      /**< pointer to increase by the number of deleted constraints */
   SCIP_RESULT*          result              /**< pointer to store the result of the presolving call */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP* subscip;
   SCIP_HASHMAP* varmap;
   SCIP_CONS* newcons;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   SCIP_Bool success;
   int i;

   assert(scip != NULL);
   assert(presoldata != NULL);
   assert(conss != NULL);
   assert(nconss > 0);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(nsolvedprobs != NULL);

   *result = SCIP_DIDNOTRUN;

#ifndef SCIP_DEBUG
   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "build sub-SCIP for component %d: %d vars (%d bin, %d int, %d cont), %d conss\n",
      compnr, nvars, nbinvars, nintvars, nvars - nintvars - nbinvars, nconss);
#else
   SCIPdebugMessage("build sub-SCIP for component %d: %d vars (%d bin, %d int, %d cont), %d conss\n",
      compnr, nvars, nbinvars, nintvars, nvars - nintvars - nbinvars, nconss);
#endif

   /* stop if the problem has too many integer variables; only if the problems should be written we have to build it anyway */
   if( presoldata->maxintvars != -1 && (nbinvars + presoldata->intfactor * nintvars > presoldata->maxintvars)
      && !presoldata->writeproblems )
   {
#ifndef SCIP_DEBUG
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "--> not created (too many integer variables)\n");
#else
      SCIPdebugMessage("--> not created (too many integer variables)\n");
#endif

      return SCIP_OKAY;
   }

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( !SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( timelimit <= 0.0 || memorylimit <= (1.0 * nvars / SCIPgetNVars(scip)) * (1.0 * nconss / SCIPgetNConss(scip)) *
      ((SCIPgetMemUsed(scip) + SCIPgetMemExternEstim(scip))/1048576.0) )
   {
      SCIPdebugMessage("--> not created (not enough memory or time left)\n");
      return SCIP_OKAY;
   }

   /* create sub-SCIP */
   if( presoldata->subscip == NULL )
   {
      SCIP_CALL( SCIPcreate(&presoldata->subscip) );
      subscip = presoldata->subscip;

      /* copy plugins, we omit pricers (because we do not run if there are active pricers) and dialogs */
      success = TRUE;
      SCIP_CALL( SCIPcopyPlugins(scip, subscip, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, &success) );

      /* abort if the plugins were not successfully copied */
      if( !success )
      {
         SCIP_CALL( SCIPfree(&presoldata->subscip) );
         presoldata->subscip = NULL;
         presoldata->pluginscopied = FALSE;

         return SCIP_OKAY;
      }
      /* copy parameter settings */
      SCIP_CALL( SCIPcopyParamSettings(scip, subscip) );

      /* set gap limit to 0 */
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/gap", 0.0) );

      /* reduce the effort spent for hash tables */
      if( !SCIPisParamFixed(subscip, "misc/usevartable") )
      {
         SCIP_CALL( SCIPsetBoolParam(subscip, "misc/usevartable", FALSE) );
      }
      if( !SCIPisParamFixed(subscip, "misc/useconstable") )
      {
         SCIP_CALL( SCIPsetBoolParam(subscip, "misc/useconstable", FALSE) );
      }
      if( !SCIPisParamFixed(subscip, "misc/usesmalltables") )
      {
         SCIP_CALL( SCIPsetBoolParam(subscip, "misc/usesmalltables", TRUE) );
      }

      /* do not catch control-C */
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

#ifndef SCIP_MORE_DEBUG
      /* disable output */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
#endif
   }
   else
   {
      subscip = presoldata->subscip;
   }

   /* set time and memory limit for the subproblem */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   if( nbinvars + nintvars > 0 )
   {
      /* set node limit */
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", presoldata->nodelimit) );
   }
   else
   {
      /* if we have only continuous variables, solving the root should be enough;
       * this avoids to spend much time in a nonlinear subscip with only continuous variables */
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", 1LL) );
   }

   /* create problem in sub-SCIP */
   /* get name of the original problem and add "comp_nr" */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp_%d", SCIPgetProbName(scip), compnr);
   SCIP_CALL( SCIPcreateProb(subscip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* create variable hashmap */
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), 10 * nvars) );

   for( i = 0; i < nconss; ++i )
   {
      /* copy the constraint */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", SCIPconsGetName(conss[i]));
      SCIP_CALL( SCIPgetConsCopy(scip, subscip, conss[i], &newcons, SCIPconsGetHdlr(conss[i]), varmap, consmap, name,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, &success) );

      /* abort if constraint was not successfully copied */
      if( !success )
         goto TERMINATE;

      SCIP_CALL( SCIPaddCons(subscip, newcons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &newcons) );
   }

   /* write the problem, if requested */
   if( presoldata->writeproblems )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp_%d.lp", SCIPgetProbName(scip), compnr);
      SCIPdebugMessage("write problem to file %s\n", name);
      SCIP_CALL( SCIPwriteOrigProblem(subscip, name, NULL, FALSE) );
   }

   /* the following asserts are not true, because some aggregations in the original scip instance could not get resolved
    * inside some constraints, so the copy (subscip) will have also some inactive variables which were copied
    */
#if 0
   /* there might be less variables in the subscip, because variables might be cancelled out during copying constraint
    * when transferring variables to active variables
    */
   assert(nbinvars >= SCIPgetNBinVars(subscip));
   assert(nintvars >= SCIPgetNIntVars(subscip));
#endif

   /* In extended debug mode, we want to be informed if the number of variables was reduced during copying.
    * This might happen, since the components presolver uses SCIPgetConsVars() and then SCIPgetActiveVars() to get the
    * active representation, while SCIPgetConsCopy() might use SCIPgetProbvarLinearSum() and this might cancel out some
    * of the active variables and cannot be avoided. However, we want to notice it and check whether the constraint
    * handler could do something more clever.
    */
#ifdef SCIP_MORE_DEBUG
   if( nvars > SCIPgetNVars(subscip) )
   {
      SCIPinfoMessage(scip, NULL, "copying component %d reduced number of variables: %d -> %d\n", compnr, nvars, SCIPgetNVars(subscip));
   }
#endif

   if( presoldata->maxintvars == -1 || (SCIPgetNBinVars(subscip) + presoldata->intfactor * SCIPgetNIntVars(subscip) <= presoldata->maxintvars) )
   {
      /* solve the subproblem */
      SCIP_CALL( SCIPsolve(subscip) );

#ifdef SCIP_MORE_DEBUG
      SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

      SCIPstatistic( updateStatisticsSubsolvetime(presoldata, SCIPgetSolvingTime(subscip)) );

      if( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
      {
         SCIP_SOL* sol;
         SCIP_VAR* subvar;
         SCIP_Bool feasible;
         SCIP_Bool infeasible;
         SCIP_Bool fixed;

         ++(*nsolvedprobs);

         sol = SCIPgetBestSol(subscip);

#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPcheckSolOrig(subscip, sol, &feasible, TRUE, TRUE) );
#else
         SCIP_CALL( SCIPcheckSolOrig(subscip, sol, &feasible, FALSE, FALSE) );
#endif

         SCIPdebugMessage("--> solved to optimality: time=%.2f, solution is%s feasible\n", SCIPgetSolvingTime(subscip), feasible ? "" : " not");

         if( feasible )
         {
            SCIP_Real glb;
            SCIP_Real gub;

            /* get values of variables in the optimal solution */
            for( i = 0; i < nvars; ++i )
            {
               subvar = (SCIP_VAR*)SCIPhashmapGetImage(varmap, vars[i]);

               /* get global bounds */
               glb = SCIPvarGetLbGlobal(vars[i]);
               gub = SCIPvarGetUbGlobal(vars[i]);

               if( subvar != NULL )
               {
                  /* get solution value from optimal solution of the component */
                  fixvals[i] = SCIPgetSolVal(subscip, sol, subvar);

                  assert(SCIPisFeasLE(scip, fixvals[i], SCIPvarGetUbLocal(vars[i])));
                  assert(SCIPisFeasGE(scip, fixvals[i], SCIPvarGetLbLocal(vars[i])));

                  /* checking a solution is done with a relative tolerance of feasibility epsilon, if we really want to
                   * change the bounds of the variables by fixing them, the old bounds must not be violated by more than
                   * the absolute epsilon; therefore, we change the fixing values, if needed, and mark that the solution
                   * has to be checked again
                   */
                  if( SCIPisGT(scip, fixvals[i], gub) )
                  {
                     fixvals[i] = gub;
                     feasible = FALSE;
                  }
                  else if( SCIPisLT(scip, fixvals[i], glb) )
                  {
                     fixvals[i] = glb;
                     feasible = FALSE;
                  }
                  assert(SCIPisLE(scip, fixvals[i], SCIPvarGetUbLocal(vars[i])));
                  assert(SCIPisGE(scip, fixvals[i], SCIPvarGetLbLocal(vars[i])));
               }
               else
               {
                  /* the variable was not copied, so it was cancelled out of constraints during copying;
                   * thus, the variable is not constrained and we fix it to its best bound
                   */
                  if( SCIPisPositive(scip, SCIPvarGetObj(vars[i])) )
                     fixvals[i] = glb;
                  else if( SCIPisNegative(scip, SCIPvarGetObj(vars[i])) )
                     fixvals[i] = gub;
                  else
                  {
                     fixvals[i] = 0.0;
                     fixvals[i] = MIN(fixvals[i], gub);
                     fixvals[i] = MAX(fixvals[i], glb);
                  }
               }
            }

            /* the solution value of at least one variable is feasible with a relative tolerance of feasibility epsilon,
             * but infeasible with an absolute tolerance of epsilon; try to set the variables to the bounds and check
             * solution again (changing the values might now introduce infeasibilities of constraints)
             */
            if( !feasible )
            {
               SCIPdebugMessage("solution violated bounds by more than epsilon, check the corrected solution...\n");

               SCIP_CALL( SCIPfreeTransform(subscip) );

               SCIP_CALL( SCIPcreateOrigSol(subscip, &sol, NULL) );

               /* get values of variables in the optimal solution */
               for( i = 0; i < nvars; ++i )
               {
                  subvar = (SCIP_VAR*)SCIPhashmapGetImage(varmap, vars[i]);

                  SCIP_CALL( SCIPsetSolVal(subscip, sol, subvar, fixvals[i]) );
               }

               /* check the solution; integrality and bounds should be fulfilled and do not have to be checked */
               SCIP_CALL( SCIPcheckSol(subscip, sol, FALSE, FALSE, FALSE, TRUE, &feasible) );

#ifndef NDEBUG
               /* in debug mode, we additionally check integrality and bounds */
               if( feasible )
               {
                  SCIP_CALL( SCIPcheckSol(subscip, sol, FALSE, TRUE, TRUE, FALSE, &feasible) );
                  assert(feasible);
               }
#endif

               SCIPdebugMessage("--> corrected solution is%s feasible\n", feasible ? "" : " not");

               SCIP_CALL( SCIPfreeSol(subscip, &sol) );
            }

            /* if the solution is feasible, fix variables and delete constraints of the component */
            if( feasible )
            {
               /* fix variables */
               for( i = 0; i < nvars; ++i )
               {
                  assert(SCIPisLE(scip, fixvals[i], SCIPvarGetUbLocal(vars[i])));
                  assert(SCIPisGE(scip, fixvals[i], SCIPvarGetLbLocal(vars[i])));
                  assert(SCIPisLE(scip, fixvals[i], SCIPvarGetUbGlobal(vars[i])));
                  assert(SCIPisGE(scip, fixvals[i], SCIPvarGetLbGlobal(vars[i])));

                  SCIP_CALL( SCIPfixVar(scip, vars[i], fixvals[i], &infeasible, &fixed) );
                  assert(!infeasible);
                  assert(fixed);
                  (*ndeletedvars)++;
               }

               /* delete constraints */
               for( i = 0; i < nconss; ++i )
               {
                  SCIP_CALL( SCIPdelCons(scip, conss[i]) );
                  (*ndeletedconss)++;
               }
            }
         }
      }
      else if( SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE )
      {
         *result = SCIP_CUTOFF;
      }
      else if( SCIPgetStatus(subscip) == SCIP_STATUS_UNBOUNDED || SCIPgetStatus(subscip) == SCIP_STATUS_INFORUNBD )
      {
         /* TODO: store unbounded ray in original SCIP data structure */
         *result = SCIP_UNBOUNDED;
      }
      else
      {
         SCIPdebugMessage("--> solving interrupted (status=%d, time=%.2f)\n",
            SCIPgetStatus(subscip), SCIPgetSolvingTime(subscip));

         /* transfer global fixings to the original problem; we can only do this, if we did not find a solution in the
          * subproblem, because otherwise, the primal bound might lead to dual reductions that cannot be transferred to
          * the original problem without also transferring the possibly suboptimal solution (which is currently not
          * possible)
          */
         if( SCIPgetNSols(subscip) == 0 )
         {
            SCIP_Bool infeasible;
            SCIP_Bool tightened;
            int ntightened;

            ntightened = 0;

            for( i = 0; i < nvars; ++i )
            {
               assert( SCIPhashmapExists(varmap, vars[i]) );

               SCIP_CALL( SCIPtightenVarLb(scip, vars[i], SCIPvarGetLbGlobal((SCIP_VAR*)SCIPhashmapGetImage(varmap, vars[i])), FALSE,
                     &infeasible, &tightened) );
               assert(!infeasible);
               if( tightened )
                  ntightened++;

               SCIP_CALL( SCIPtightenVarUb(scip, vars[i], SCIPvarGetUbGlobal((SCIP_VAR*)SCIPhashmapGetImage(varmap, vars[i])), FALSE,
                     &infeasible, &tightened) );
               assert(!infeasible);
               if( tightened )
                  ntightened++;
            }
            SCIPdebugMessage("--> tightened %d bounds of variables due to global bounds in the sub-SCIP\n", ntightened);
         }
      }
   }
   else
   {
      SCIPdebugMessage("--> not solved (too many integer variables)\n");
   }

 TERMINATE:
   SCIPhashmapFree(&varmap);

   return SCIP_OKAY;
}

/** loop over constraints, get active variables and fill directed graph */
static
SCIP_RETCODE fillDigraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int*                  firstvaridxpercons, /**< array to store for each constraint the index in the local vars array
                                              *   of the first variable of the constraint */
   SCIP_Bool*            success             /**< flag indicating successful directed graph filling */
   )
{
   SCIP_VAR** consvars;
   int requiredsize;
   int nconsvars;
   int nvars;
   int idx1;
   int idx2;
   int c;
   int v;

   assert(scip != NULL);
   assert(digraph != NULL);
   assert(conss != NULL);
   assert(firstvaridxpercons != NULL);
   assert(success != NULL);

   *success = TRUE;

   nconsvars = 0;
   requiredsize = 0;
   nvars = SCIPgetNVars(scip);

   /* use big buffer for storing active variables per constraint */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );

   for( c = 0; c < nconss; ++c )
   {
      /* check for reached timelimit */
      if( (c % 1000 == 0) && SCIPisStopped(scip) )
      {
         *success = FALSE;
         break;
      }

      /* get number of variables for this constraint */
      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nconsvars, success) );

      if( !(*success) )
         break;

      if( nconsvars > nvars )
      {
         nvars = nconsvars;
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvars, nvars) );
      }

#ifndef NDEBUG
      /* clearing variables array to check for consistency */
      if( nconsvars == nvars )
      {
	 BMSclearMemoryArray(consvars, nconsvars);
      }
      else
      {
	 assert(nconsvars < nvars);
	 BMSclearMemoryArray(consvars, nconsvars + 1);
      }
#endif

      /* get variables for this constraint */
      SCIP_CALL( SCIPgetConsVars(scip, conss[c], consvars, nvars, success) );

      if( !(*success) )
      {
#ifndef NDEBUG
	 /* it looks strange if returning the number of variables was successful but not returning the variables */
	 SCIPwarningMessage(scip, "constraint <%s> returned number of variables but returning variables failed\n", SCIPconsGetName(conss[c]));
#endif
         break;
      }

#ifndef NDEBUG
      /* check if returned variables are consistent with the number of variables that were returned */
      for( v = nconsvars - 1; v >= 0; --v )
	 assert(consvars[v] != NULL);
      if( nconsvars < nvars )
	 assert(consvars[nconsvars] == NULL);
#endif

      /* transform given variables to active variables */
      SCIP_CALL( SCIPgetActiveVars(scip, consvars, &nconsvars, nvars, &requiredsize) );
      assert(requiredsize <= nvars);

      if( nconsvars > 0 )
      {
         idx1 = SCIPvarGetProbindex(consvars[0]);
         assert(idx1 >= 0);

         /* save index of the first variable for later component assignment */
         firstvaridxpercons[c] = idx1;

         if( nconsvars > 1 )
         {
            /* create sparse directed graph
             * sparse means, to add only those edges necessary for component calculation
             */
            for( v = 1; v < nconsvars; ++v )
            {
               idx2 = SCIPvarGetProbindex(consvars[v]);
               assert(idx2 >= 0);

               /* we add only one directed edge, because the other direction is automatically added for component computation */
               SCIP_CALL( SCIPdigraphAddArc(digraph, idx1, idx2, NULL) );
            }
         }
      }
      else
         firstvaridxpercons[c] = -1;
   }

   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}

/** use components to assign variables and constraints to the subscips
 *  and try to solve all subscips having not too many integer variables
 */
static
SCIP_RETCODE splitProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_CONS**           conss,              /**< constraints */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_Real*            fixvals,            /**< array to store the fixing values */
   int                   nconss,             /**< number of constraints */
   int                   nvars,              /**< number of variables */
   int*                  components,         /**< array with component number for every variable */
   int                   ncomponents,        /**< number of components */
   int*                  firstvaridxpercons, /**< array with index of first variable in vars array for each constraint */
   int*                  nsolvedprobs,       /**< pointer to store the number of solved subproblems */
   int*                  ndeletedvars,       /**< pointer to store the number of deleted variables */
   int*                  ndeletedconss,      /**< pointer to store the number of deleted constraints */
   SCIP_RESULT*          result              /**< pointer to store the result of the presolving call */
  )
{
   SCIP_HASHMAP* consmap;
   SCIP_CONS** compconss;
   SCIP_VAR** compvars;
   SCIP_Real* compfixvals;
   int* conscomponent;
   int ncompconss;
   int ncompvars;
   int nbinvars;
   int nintvars;
   int comp;
   int compvarsstart;
   int compconssstart;
   int v;
   int c;

   assert(scip != NULL);
   assert(presoldata != NULL);
   assert(conss != NULL);
   assert(components != NULL);
   assert(firstvaridxpercons != NULL);
   assert(nsolvedprobs != NULL);
   assert(result != NULL);

   *nsolvedprobs = 0;

   /* hashmap mapping from original constraints to constraints in the sub-SCIPs (for performance reasons) */
   SCIP_CALL( SCIPhashmapCreate(&consmap, SCIPblkmem(scip), 10 * SCIPgetNConss(scip)) );

   SCIP_CALL( SCIPallocBufferArray(scip, &conscomponent, nconss) );

   /* do the mapping from calculated components per variable to corresponding
    * constraints and sort the component-arrays for faster finding the
    * actual variables and constraints belonging to one component
    */
   for( c = 0; c < nconss; c++ )
      conscomponent[c] = (firstvaridxpercons[c] == -1 ? -1 : components[firstvaridxpercons[c]]);

   SCIPsortIntPtr(components, (void**)vars, nvars);
   SCIPsortIntPtr(conscomponent, (void**)conss, nconss);

   compvarsstart = 0;
   compconssstart = 0;

   while( compconssstart < nconss && conscomponent[compconssstart] == -1 )
      ++compconssstart;

   /* loop over all components */
   for( comp = 0; comp < ncomponents && !SCIPisStopped(scip); comp++ )
   {
      nbinvars = 0;
      nintvars = 0;

      /* count variables present in this component and categorize them */
      for( v = compvarsstart; v < nvars && components[v] == comp; ++v )
      {
         /* check whether variable is of binary or integer type */
         if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY )
            nbinvars++;
         else if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER )
            nintvars++;
      }

      /* count constraints present in this component */
      c = compconssstart;
      while( c < nconss && conscomponent[c] == comp )
         ++c;

      /* collect some statistical information */
      SCIPstatistic( updateStatisticsComp(presoldata, nbinvars, nintvars) );

      compvars = &(vars[compvarsstart]);
      compfixvals = &(fixvals[compvarsstart]);
      compconss = &(conss[compconssstart]);
      ncompvars = v - compvarsstart;
      ncompconss = c - compconssstart;
      compvarsstart = v;
      compconssstart = c;

      /* single variable without constraint */
      if( ncompconss == 0 )
      {
         SCIP_Bool infeasible;
         SCIP_Bool fixed;

         /* there is no constraint to connect variables, so there should be only one variable in the component */
         assert(ncompvars == 1);

         /* fix variable to its best bound (or 0.0, if objective is 0) */
         if( SCIPisPositive(scip, SCIPvarGetObj(compvars[0])) )
            compfixvals[0] = SCIPvarGetLbGlobal(compvars[0]);
         else if( SCIPisNegative(scip, SCIPvarGetObj(compvars[0])) )
            compfixvals[0] = SCIPvarGetUbGlobal(compvars[0]);
         else
         {
            compfixvals[0] = 0.0;
            compfixvals[0] = MIN(compfixvals[0], SCIPvarGetUbGlobal(compvars[0])); /*lint !e666*/
            compfixvals[0] = MAX(compfixvals[0], SCIPvarGetLbGlobal(compvars[0])); /*lint !e666*/
         }

#ifdef SCIP_MORE_DEBUG
	 SCIPinfoMessage(scip, NULL, "presol components: fix variable <%s>[%g,%g] (locks [%d, %d]) to %g because it occurs in no constraint\n",
            SCIPvarGetName(compvars[0]), SCIPvarGetLbGlobal(compvars[0]), SCIPvarGetUbGlobal(compvars[0]),
            SCIPvarGetNLocksUp(compvars[0]), SCIPvarGetNLocksDown(compvars[0]), compfixvals[0]);
#else
	 SCIPdebugMessage("presol components: fix variable <%s>[%g,%g] (locks [%d, %d]) to %g because it occurs in no constraint\n",
            SCIPvarGetName(compvars[0]), SCIPvarGetLbGlobal(compvars[0]), SCIPvarGetUbGlobal(compvars[0]),
            SCIPvarGetNLocksUp(compvars[0]), SCIPvarGetNLocksDown(compvars[0]), compfixvals[0]);
#endif

         SCIP_CALL( SCIPfixVar(scip, compvars[0], compfixvals[0], &infeasible, &fixed) );
         assert(!infeasible);
         assert(fixed);
         ++(*ndeletedvars);
         ++(*nsolvedprobs);

         /* update statistics */
         SCIPstatistic( updateStatisticsSingleVar(presoldata) );
      }
      /* we do not want to solve the component, if it is the last unsolved one */
      else if( ncompvars < SCIPgetNVars(scip) )
      {
         SCIP_RESULT subscipresult;

         assert(ncompconss > 0);

         /* in extended debug mode, we want to be informed about components with single constraints;
          * this is only for noticing this case and possibly handling it within the constraint handler
          */
#ifdef SCIP_MORE_DEBUG
         if( ncompconss == 1 )
         {
            SCIPinfoMessage(scip, NULL, "presol component detected component with a single constraint:\n");
            SCIP_CALL( SCIPprintCons(scip, compconss[0], NULL) );
            SCIPinfoMessage(scip, NULL, ";\n");
         }
#endif

         /* build subscip for one component and try to solve it */
         SCIP_CALL( copyAndSolveComponent(scip, presoldata, consmap, comp, compconss, ncompconss, compvars,
               compfixvals, ncompvars, nbinvars, nintvars, nsolvedprobs, ndeletedvars, ndeletedconss,
               &subscipresult) );

         if( subscipresult == SCIP_CUTOFF || subscipresult == SCIP_UNBOUNDED )
	 {
	    *result = subscipresult;
            break;
	 }

         if( !presoldata->pluginscopied )
            break;
      }
   }

   if( presoldata->subscip != NULL )
   {
      SCIP_CALL( SCIPfree(&presoldata->subscip) );
      presoldata->subscip = NULL;
   }

   SCIPfreeBufferArray(scip, &conscomponent);
   SCIPhashmapFree(&consmap);

   return SCIP_OKAY;
}

/** performs presolving by searching for components */
static
SCIP_RETCODE presolComponents(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_PRESOL*          presol,             /**< the presolver itself */
   int*                  nfixedvars,         /**< counter to increase by the number of fixed variables */
   int*                  ndelconss,          /**< counter to increase by the number of deleted constrains */
   SCIP_RESULT*          result              /**< pointer to store the result of the presolving call */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_CONS** conss;
   SCIP_CONS** tmpconss;
   SCIP_Bool success;
   int nconss;
   int ntmpconss;
   int nvars;
   int ncomponents;
   int ndeletedconss;
   int ndeletedvars;
   int nsolvedprobs;
   int c;

   assert(scip != NULL);
   assert(presol != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING || SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* do not run, if not all variables are explicitly known */
   if( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   nvars = SCIPgetNVars(scip);

   /* we do not want to run, if there are no variables left */
   if( nvars == 0 )
      return SCIP_OKAY;

   /* the presolver should be executed only if it didn't run so far or the number of variables was significantly reduced
    * since tha last run
    */
   if( (presoldata->didsearch && nvars > (1 - presoldata->reldecrease) * presoldata->lastnvars) )
      return SCIP_OKAY;

   /* if we were not able to copy all plugins, we cannot run the components presolver */
   if( !presoldata->pluginscopied )
      return SCIP_OKAY;

   /* only call the component presolver, if presolving would be stopped otherwise */
   if( !SCIPisPresolveFinished(scip) )
      return SCIP_OKAY;

   /* check for a reached timelimit */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPstatistic( resetStatistics(scip, presoldata) );

   *result = SCIP_DIDNOTFIND;
   presoldata->didsearch = TRUE;

   ncomponents = 0;
   ndeletedvars = 0;
   ndeletedconss = 0;
   nsolvedprobs = 0;

   /* collect checked constraints for component presolving */
   ntmpconss = SCIPgetNConss(scip);
   tmpconss = SCIPgetConss(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &conss, ntmpconss) );
   nconss = 0;
   for( c = 0; c < ntmpconss; c++ )
   {
      if( SCIPconsIsChecked(tmpconss[c]) )
      {
         conss[nconss] = tmpconss[c];
         nconss++;
      }
   }

   if( nvars > 1 && nconss > 1 )
   {
      SCIP_DIGRAPH* digraph;
      SCIP_VAR** vars;
      SCIP_Real* fixvals;
      int* firstvaridxpercons;
      int* varlocks;
      int v;

      /* copy variables into a local array */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &fixvals, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &firstvaridxpercons, nconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varlocks, nvars) );
      BMScopyMemoryArray(vars, SCIPgetVars(scip), nvars);

      /* count number of varlocks for each variable (up + down locks) and multiply it by 2;
       * that value is used as an estimate of the number of arcs incident to the variable's node in the digraph
       */
      for( v = 0; v < nvars; ++v )
         varlocks[v] = 4 * (SCIPvarGetNLocksDown(vars[v]) + SCIPvarGetNLocksUp(vars[v]));

      /* create and fill directed graph */
      SCIP_CALL( SCIPdigraphCreate(&digraph, nvars) );
      SCIP_CALL( SCIPdigraphSetSizes(digraph, varlocks) );
      SCIP_CALL( fillDigraph(scip, digraph, conss, nconss, firstvaridxpercons, &success) );

      if( success )
      {
         int* components;

         SCIP_CALL( SCIPallocBufferArray(scip, &components, nvars) );

         /* compute independent components */
         SCIP_CALL( SCIPdigraphComputeUndirectedComponents(digraph, 1, components, &ncomponents) );

         SCIPdebugMessage("presol components found %d undirected components\n", ncomponents);

         /* We want to sort the components in increasing size (number of variables).
          * Therefore, we now get the number of variables for each component, and rename the components
          * such that for i < j, component i has no more variables than component j.
          * @todo Perhaps sort the components by the number of binary/integer variables?
          */
         if( ncomponents > 0 )
         {
            int* compsize;
            int* permu;
            int* compvars;
            int ncompvars;
            int nbinvars;
            int nintvars;
            int ncontvars;

            SCIP_CALL( SCIPallocBufferArray(scip, &compsize, ncomponents) );
            SCIP_CALL( SCIPallocBufferArray(scip, &permu, ncomponents) );

            /* get number of variables in the components */
            for( c = 0; c < ncomponents; ++c )
            {
               SCIPdigraphGetComponent(digraph, c, &compvars, &ncompvars);
               permu[c] = c;
               nbinvars = 0;
               nintvars = 0;

               for( v = 0; v < ncompvars; ++v )
               {
                  /* check whether variable is of binary or integer type */
                  if( SCIPvarGetType(vars[compvars[v]]) == SCIP_VARTYPE_BINARY )
                     nbinvars++;
                  else if( SCIPvarGetType(vars[compvars[v]]) == SCIP_VARTYPE_INTEGER )
                     nintvars++;
               }
               ncontvars = ncompvars - nintvars - nbinvars;
               compsize[c] = (int)(1000 * (nbinvars + presoldata->intfactor * nintvars) + MIN(999, ncontvars));
            }

            /* get permutation of component numbers such that the size of the components is increasing */
            SCIPsortIntInt(compsize, permu, ncomponents);

            /* now, we need the reverse direction, i.e., for each component number, we store its new number
             * such that the components are sorted; for this, we abuse the compsize array
             */
            for( c = 0; c < ncomponents; ++c )
               compsize[permu[c]] = c;

            /* for each variable, replace the old component number by the new one */
            for( c = 0; c < nvars; ++c )
               components[c] = compsize[components[c]];

            /* create subproblems from independent components and solve them in dependence of their size */
            SCIP_CALL( splitProblem(scip, presoldata, conss, vars, fixvals, nconss, nvars, components, ncomponents, firstvaridxpercons,
                  &nsolvedprobs, &ndeletedvars, &ndeletedconss, result) );

            (*nfixedvars) += ndeletedvars;
            (*ndelconss) += ndeletedconss;

            SCIPfreeBufferArray(scip, &permu);
            SCIPfreeBufferArray(scip, &compsize);
         }

         SCIPfreeBufferArray(scip, &components);
      }

      SCIPdigraphFree(&digraph);

      SCIPfreeBufferArray(scip, &varlocks);
      SCIPfreeBufferArray(scip, &firstvaridxpercons);
      SCIPfreeBufferArray(scip, &fixvals);
      SCIPfreeBufferArray(scip, &vars);
   }

   SCIPfreeBufferArray(scip, &conss);

   /* update result to SCIP_SUCCESS, if reductions were found;
    * if result was set to infeasible or unbounded before, leave it as it is
    */
   if( (ndeletedvars > 0 || ndeletedconss > 0) && ((*result) == SCIP_DIDNOTFIND) )
      *result = SCIP_SUCCESS;

   /* print statistics */
   SCIPstatistic( printStatistics(presoldata) );

   SCIPdebugMessage("%d components, %d solved, %d deleted constraints, %d deleted variables, result = %d\n",
      ncomponents, nsolvedprobs, ndeletedconss, ndeletedvars, *result);
#ifdef NDEBUG
   SCIPstatisticMessage("%d components, %d solved, %d deleted constraints, %d deleted variables, result = %d\n",
      ncomponents, nsolvedprobs, ndeletedconss, ndeletedvars, *result);
#endif

   presoldata->lastnvars = SCIPgetNVars(scip);

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeComponents)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** initialization method of presolver (called after problem was transformed) */
static
SCIP_DECL_PRESOLINIT(presolInitComponents)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* initialize statistics */
   SCIPstatistic( SCIP_CALL( initStatistics(scip, presoldata) ) );

   presoldata->subscip = NULL;
   presoldata->didsearch = FALSE;
   presoldata->pluginscopied = TRUE;

   return SCIP_OKAY;
}

#ifdef SCIP_STATISTIC
/** deinitialization method of presolver (called before transformed problem is freed) */
static
SCIP_DECL_PRESOLEXIT(presolExitComponents)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPstatistic( freeStatistics(scip, presoldata) );

   return SCIP_OKAY;
}
#endif


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecComponents)
{  /*lint --e{715}*/

   SCIP_CALL( presolComponents(scip, presol, nfixedvars, ndelconss, result) );

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the components presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolComponents(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create components presolver data */
   SCIP_CALL( SCIPallocMemory(scip, &presoldata) );

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_DELAY, presolExecComponents, presoldata) );

   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeComponents) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitComponents) );
   SCIPstatistic( SCIP_CALL( SCIPsetPresolExit(scip, presol, presolExitComponents) ) );

   /* add presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/components/writeproblems",
         "should the single components be written as an .lp-file?",
         &presoldata->writeproblems, FALSE, DEFAULT_WRITEPROBLEMS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/components/maxintvars",
         "maximum number of integer (or binary) variables to solve a subproblem directly (-1: unlimited)",
         &presoldata->maxintvars, FALSE, DEFAULT_MAXINTVARS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "presolving/components/nodelimit",
         "maximum number of nodes to be solved in subproblems",
         &presoldata->nodelimit, FALSE, DEFAULT_NODELIMIT, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/components/intfactor",
         "the weight of an integer variable compared to binary variables",
         &presoldata->intfactor, FALSE, DEFAULT_INTFACTOR, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/components/reldecrease",
         "percentage by which the number of variables has to be decreased after the last component solving to allow running again (1.0: do not run again)",
         &presoldata->reldecrease, FALSE, DEFAULT_RELDECREASE, 0.0, 1.0, NULL, NULL) );


   return SCIP_OKAY;
}
