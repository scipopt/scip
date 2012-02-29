/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define SCIP_DEBUG
#define WITH_STATISTICS
/**@file   presol_components.c
 * @brief  solve independent components in advance
 * @author Dieter Weninger
 * @author Gerald Gamrath
 *
 * TODO: simulation of presolving without solve
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/presol_components.h"


#define PRESOL_NAME            "components"
#define PRESOL_DESC            "components presolver"
#define PRESOL_PRIORITY        -9200000      /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS              0      /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY               TRUE      /**< should presolver be delayed, if other presolvers found reductions? */

#define DEFAULT_WRITEPROBLEMS     FALSE      /**< should the single components be written as an .lp-file? */
#define DEFAULT_MAXINTVARS           20      /**< maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving) */
#define DEFAULT_NODELIMIT       10000LL      /**< maximum number of nodes to be solved in subproblems */
#define DEFAULT_INTFACTOR           1.0      /**< the weight of an integer variable compared to binary variables */

#ifdef WITH_STATISTICS
static int NCATEGORIES = 4;
static int CATLIMITS[] = {20,50,100};
#endif

/*
 * Data structures
 */

/** control parameters */
struct SCIP_PresolData
{
   SCIP_Bool             didsearch;          /** did the presolver already search for components? */
   SCIP_Bool             writeproblems;      /** should the single components be written as an .lp-file? */
   int                   maxintvars;         /** maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving) */
   SCIP_Longint          nodelimit;          /** maximum number of nodes to be solved in subproblems */
   SCIP_Real             intfactor;          /** the weight of an integer variable compared to binary variables */
#ifdef WITH_STATISTICS
   int*                  compspercat;        /** number of components of the different categories */
   int                   nsinglevars;        /** number of components with a single variable without constraint */
   SCIP_Real             subsolvetime;       /** total solving time of the subproblems */
#endif
};

/*
 * Statistic methods
 */

#ifdef WITH_STATISTICS
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

/** initialize data for statistics */
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

#else
#define initStatistics(scip, presoldata) SCIP_OKAY
#define freeStatistics(scip, presoldata) /**/
#define updateStatisticsComp(presoldata, nbinvars, nintvars) /**/
#define updateStatisticsSingleVar(presoldata) /**/
#define updateStatisticsSubsolvetime(presoldata, subsolvetime) /**/
#define printStatistics(presoldata) /**/
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
   int                   nvars,              /**< number of variables contained in this component */
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

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if( timelimit <= 0.0 || memorylimit <= 0.0 )
      goto TERMINATE;

   /* create sub-SCIP */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* copy plugins, we omit pricers (because we do not run if there are active pricers) and dialogs */
   success = TRUE;
   SCIP_CALL( SCIPcopyPlugins(scip, subscip, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, &success) );

   /* abort if the plugins could not be copied successfully */
   if( !success )
      goto TERMINATE;

   /* copy parameter settings */
   SCIP_CALL( SCIPcopyParamSettings(scip, subscip) );

   /* set time and memory limit for the subproblem */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   /* set node limit */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", presoldata->nodelimit) );

   /* set gap limit to 0 */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/gap", 0.0) );

   /* reduce the effort spent for hash tables */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/usevartable", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/useconstable", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/usesmalltables", TRUE) );

   /* do not catch control-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

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

   if( SCIPgetNBinVars(subscip) + presoldata->intfactor * SCIPgetNIntVars(subscip) <= presoldata->maxintvars )
   {
      /* solve the subproblem */
      SCIP_CALL( SCIPsolve(subscip) );

      updateStatisticsSubsolvetime(presoldata, SCIPgetSolvingTime(subscip));

      if( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
      {
         SCIP_SOL* sol;
         SCIP_Bool infeasible;
         SCIP_Bool fixed;

         ++(*nsolvedprobs);

         sol = SCIPgetBestSol(subscip);

         /* fix variables */
         for( i = 0; i < nvars; ++i )
         {
            assert( SCIPhashmapExists(varmap, vars[i]) );

            SCIP_CALL( SCIPfixVar(scip, vars[i], SCIPgetSolVal(subscip, sol, SCIPhashmapGetImage(varmap, vars[i])),
                  &infeasible, &fixed) );
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
      else if( SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE )
      {
         *result = SCIP_CUTOFF;
      }
      else if( SCIPgetStatus(subscip) == SCIP_STATUS_UNBOUNDED )
      {
         /* TODO: store unbounded ray in original SCIP data structure */
         *result = SCIP_UNBOUNDED;
      }
      /* TODO: transfer global fixings if problem was not solved to optimality */
      else
      {
         SCIPdebugMessage("++++++++++++++ sub-SCIP for component %d not solved (status=%d, time=%.2f): %d vars (%d bin, %d int, %d impl, %d cont), %d conss\n",
            compnr, SCIPgetStatus(subscip), SCIPgetSolvingTime(subscip), nvars, SCIPgetNBinVars(subscip), SCIPgetNIntVars(subscip), SCIPgetNImplVars(subscip),
            SCIPgetNContVars(subscip), nconss);
      }
   }

 TERMINATE:
   SCIPhashmapFree(&varmap);
   SCIP_CALL( SCIPfree(&subscip) );

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
      /* get number of variables for this constraint */
      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nconsvars, success) );

      if( !(*success) )
         break;

      if( nconsvars > nvars )
      {
         nvars = nconsvars;
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvars, nvars) );
      }

      /* get variables for this constraint */
      SCIP_CALL( SCIPgetConsVars(scip, conss[c], consvars, nvars, success) );

      if( !(*success) )
         break;

      /* transform given variables to active variables */
      SCIP_CALL( SCIPgetActiveVars(scip, consvars, &nconsvars, nvars, &requiredsize) );
      assert(requiredsize <= nvars);

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
            SCIP_CALL( SCIPdigraphAddEdge(digraph, idx1, idx2) );
         }
      }
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
   int                   nconss,             /**< number of constraints */
   int                   nvars,              /**< number of variables */
   int*                  components,         /**< array with component number for every variable */
   int                   ncomponents,        /**< number of components */
   int*                  firstvaridxpercons, /**< array with index of first variable in vars array for each constraint */
   int*                  nsolvedprobs,       /**< number of solved subproblems */
   int*                  ndeletedvars,       /**< pointer to store the number of deleted variables */
   int*                  ndeletedconss,      /**< pointer to store the number of deleted constraints */
   SCIP_RESULT*          result              /**< pointer to store the result of the presolving call */
  )
{
   SCIP_HASHMAP* consmap;
   SCIP_CONS** compconss;
   SCIP_VAR** compvars;
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
      conscomponent[c] = components[firstvaridxpercons[c]];

   SCIPsortIntPtr(components, (void**)vars, nvars);
   SCIPsortIntPtr(conscomponent, (void**)conss, nconss);

   compvarsstart = 0;
   compconssstart = 0;

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
      updateStatisticsComp(presoldata, nbinvars, nintvars);

      compvars = &(vars[compvarsstart]);
      compconss = &(conss[compconssstart]);
      ncompvars = v - compvarsstart;
      ncompconss = c - compconssstart;
      compvarsstart = v;
      compconssstart = c;

      if( (nbinvars + presoldata->intfactor * nintvars <= presoldata->maxintvars) || presoldata->writeproblems )
      {
         /* single variable without constraint */
         if( ncompconss == 0 )
         {
            SCIP_Real fixval;
            SCIP_Bool infeasible;
            SCIP_Bool fixed;

            /* there is no constraint to connect variables, so there should be only one variable in the component */
            assert(ncompvars == 1);

            fixval = 0.0;

            /* fix variable to its best bound (or 0.0, if objective is 0) */
            if( SCIPisPositive(scip, SCIPvarGetObj(compvars[0])) )
               fixval = SCIPvarGetLbGlobal(compvars[0]);
            else if( SCIPisNegative(scip, SCIPvarGetObj(compvars[0])) )
               fixval = SCIPvarGetUbGlobal(compvars[0]);

            SCIP_CALL( SCIPfixVar(scip, compvars[0], fixval, &infeasible, &fixed) );
            assert(!infeasible);
            assert(fixed);
            (*ndeletedvars)++;

            /* update statistics */
            updateStatisticsSingleVar(presoldata);
         }
         else
         {
            assert(ncompconss > 0);

            /* build subscip for one component and try to solve it */
            SCIP_CALL( copyAndSolveComponent(scip, presoldata, consmap, comp, compconss, ncompconss, compvars, ncompvars,
                  nsolvedprobs, ndeletedvars, ndeletedconss, result) );

            if( *result == SCIP_CUTOFF )
               break;
         }
      }
      else
      {
         SCIPdebugMessage("++++++++++++++ sub-SCIP for component %d not created: %d vars (%d bin, %d int, %d cont), %d conss\n",
            comp, ncompvars, nbinvars, nintvars, ncompvars - nintvars - nbinvars, ncompconss);
      }
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
   SCIP_DIGRAPH* digraph;
   SCIP_CONS** conss;
   SCIP_CONS** tmpconss;
   SCIP_VAR** vars;
   SCIP_VAR** tmpvars;
   int* firstvaridxpercons;
   int nconss;
   int ntmpconss;
   int nvars;
   int* components;
   int ncomponents;
   int ndeletedconss;
   int ndeletedvars;
   int nsolvedprobs;
   int c;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(presol != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING || SCIPinProbing(scip) )
      return SCIP_OKAY;

   if( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* the presolver should be executed only once */
   if( presoldata->didsearch )
      return SCIP_OKAY;

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

   /* get number of variables and copy variables into a local array */
   nvars = SCIPgetNVars(scip);
   tmpvars = SCIPgetVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   BMScopyMemoryArray(vars, tmpvars, nvars);

   if( nvars > 1 && nconss > 1 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &firstvaridxpercons, nconss) );

      /* create and fill directed graph */
      SCIP_CALL( SCIPdigraphCreate(&digraph, nvars) );
      SCIP_CALL( fillDigraph(scip, digraph, conss, nconss, firstvaridxpercons, &success) );

      if( success )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &components, nvars) );

         /* compute independent components */
         SCIP_CALL( SCIPdigraphComputeUndirectedComponents(digraph, 1, components, &ncomponents) );

         /* create subproblems from independent components and solve them in dependence on their size */
         SCIP_CALL( splitProblem(scip, presoldata, conss, vars, nconss, nvars, components, ncomponents, firstvaridxpercons,
               &nsolvedprobs, &ndeletedvars, &ndeletedconss, result) );

         (*nfixedvars) += ndeletedvars;
         (*ndelconss) += ndeletedconss;

         SCIPfreeBufferArray(scip, &components);
      }

      SCIPdigraphFree(&digraph);

      SCIPfreeBufferArray(scip, &firstvaridxpercons);
   }

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &conss);

   /* update result to SCIP_SUCCESS, if reductions were found;
    * if result was set to infeasible or unbounded before, leave it as it is
    */
   if( (ndeletedvars > 0 || ndeletedconss > 0) && ((*result) == SCIP_DIDNOTFIND) )
      *result = SCIP_SUCCESS;

   /* print statistics */
   printStatistics(presoldata);

   SCIPdebugMessage("### %d components, %d solved, %d deleted constraints, %d deleted variables\n",
      ncomponents, nsolvedprobs, ndeletedconss, ndeletedvars);

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PRESOLCOPY(presolCopyComponents)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolCopyComponents NULL
#endif


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
   SCIP_CALL( initStatistics(scip, presoldata) );

   presoldata->didsearch = FALSE;

   return SCIP_OKAY;
}


/** deinitialization method of presolver (called before transformed problem is freed) */
static
SCIP_DECL_PRESOLEXIT(presolExitComponents)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   freeStatistics(scip, presoldata);

   return SCIP_OKAY;
}

/* define unused callbacks as NULL */
#define presolInitpreComponents NULL
#define presolExitpreComponents NULL


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecComponents)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTRUN;

   /* only call the component presolver, if presolving would be stopped otherwise */
   if( SCIPisPresolveFinished(scip) )
   {
      SCIP_CALL( presolComponents(scip, presol, nfixedvars, ndelconss, result) );
   }

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

   /* create components presolver data */
   SCIP_CALL( SCIPallocMemory(scip, &presoldata) );

   /* include presolver */
   SCIP_CALL( SCIPincludePresol(scip,
         PRESOL_NAME,
         PRESOL_DESC,
         PRESOL_PRIORITY,
         PRESOL_MAXROUNDS,
         PRESOL_DELAY,
         presolCopyComponents,
         presolFreeComponents,
         presolInitComponents,
         presolExitComponents,
         presolInitpreComponents,
         presolExitpreComponents,
         presolExecComponents,
         presoldata) );

   /* add presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/components/writeproblems",
         "should the single components be written as an .lp-file?",
         &presoldata->writeproblems, FALSE, DEFAULT_WRITEPROBLEMS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/components/maxintvars",
         "maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving)",
         &presoldata->maxintvars, FALSE, DEFAULT_MAXINTVARS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "presolving/components/nodelimit",
         "maximum number of nodes to be solved in subproblems",
         &presoldata->nodelimit, FALSE, DEFAULT_NODELIMIT, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/components/intfactor",
         "the weight of an integer variable compared to binary variables",
         &presoldata->intfactor, FALSE, DEFAULT_INTFACTOR, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
