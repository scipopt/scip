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
#define SCIP_DEBUG
#define SCIP_MORE_DEBUG
/**@file   prop_components.c
 * @brief  identify and solve independent components
 * @author Gerald Gamrath
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_components.h"

#define PROP_NAME            "components"
#define PROP_DESC            "components propagator"
#define PROP_TIMING                 (SCIP_PROPTIMING_BEFORELP | SCIP_PROPTIMING_AFTERLPLOOP)


#define PROP_PRIORITY               -8000000 /**< propagation priority */
#define PROP_FREQ                          1 /**< propagation frequency */
#define PROP_DELAY                      TRUE /**< should propagation method be delayed, if other propagators found
                                              *   reductions? */
#define PROP_PRESOL_PRIORITY        -8000000 /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
#define PROP_PRESOL_MAXROUNDS             -1 /**< maximal number of propagation rounds the propagator participates in (-1: no limit) */
#define PROP_PRESOLTIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolving method (fast, medium, or exhaustive) */

#define INTFACTOR 10

#ifdef SCIP_STATISTIC
static int NCATEGORIES = 6;
static int CATLIMITS[] = {0,20,50,100,500};
#endif

/*
 * Data structures
 */

typedef struct Problem PROBLEM;

/** data related to one component */
typedef struct Component
{
   PROBLEM*              problem;            /** the problem this component belongs to */
   SCIP*                 subscip;            /** sub-SCIP representing the component */
   SCIP_VAR**            vars;               /** variables belonging to this component (in complete problem) */
   SCIP_VAR**            subvars;            /** variables belonging to this component (in subscip) */
   SCIP_Real             lastdualbound;      /** dual bound after last optimization call for this component */
   SCIP_Real             lastprimalbound;    /** primal bound after last optimization call for this component */
   SCIP_Longint          lastnodelimit;      /** node limit of last optimization call for this component */
   SCIP_STATUS           laststatus;         /** solution status of last optimization call for this component */
   int                   ncalls;             /** number of optimization calls for this component */
   int                   lastsolindex;       /** index of best solution after last optimization call for this component */
   int                   nvars;              /** number of variables belonging to this component */
   int                   number;             /** component number */

   SCIP_PQUEUE*          subproblemqueue;    /** subproblems (nodes of the branch and bound tree of the component) treated individually */
} COMPONENT;

/** data related to one problem */
struct Problem
{
   SCIP*                 scip;               /** the SCIP instance this problem belongs to */
   COMPONENT**           components;         /** independent components into which the problem can be divided */
   SCIP_PQUEUE*          compqueue;          /** priority queue for components */
   SCIP_SOL*             bestsol;            /** best solution found so far for the problem */
   char*                 name;               /** name of the problem */
   SCIP_Real             lowerbound;         /** lower bound of the problem */
   int                   ncomponents;        /** number of independent components into which the problem can be divided */
   int                   nfeascomps;         /** number of components for which a feasible solution was found */
   int                   nsolvedcomps;       /** number of components solved to optimality */
   int                   nlowerboundinf;     /** number of components with lower bound equal to -infinity */

};


/** control parameters */
struct SCIP_PropData
{
   PROBLEM*              problem;            /** the main problem */
   COMPONENT*            component;          /** component the current SCIP instance belongs to */
#ifdef SCIP_STATISTIC
   int*                  compspercat;        /** number of components of the different categories */
   SCIP_Real             subsolvetime;       /** total solving time of the subproblems */
   int                   nsinglevars;        /** number of components with a single variable without constraint */
#endif
};


/** comparison method for sorting components */
static
SCIP_DECL_SORTPTRCOMP(componentSort)
{
   SCIP* scip;
   COMPONENT* comp1;
   COMPONENT* comp2;
   SCIP_Real gap1;
   SCIP_Real gap2;

   assert(elem1 != NULL);
   assert(elem2 != NULL);

   comp1 = (COMPONENT*)elem1;
   comp2 = (COMPONENT*)elem2;

   scip = comp1->problem->scip;

   if( comp1->ncalls == 0 )
      if( comp2->ncalls == 0 )
         return comp1->number - comp2->number;
      else
         return -1;
   else if( comp2->ncalls == 0 )
      return 1;

   gap1 = SQR(comp1->lastprimalbound - comp1->lastdualbound) / comp1->ncalls;
   gap2 = SQR(comp2->lastprimalbound - comp2->lastdualbound) / comp2->ncalls;

   if( SCIPisFeasGT(scip, gap1, gap2) )
      return -1;
   else if( SCIPisFeasLT(scip, gap1, gap2) )
      return +1;
   else
      return comp1->number - comp2->number;
}

/** comparison method for sorting subproblems */
static
SCIP_DECL_SORTPTRCOMP(subproblemSort)
{
   PROBLEM* prob1;
   PROBLEM* prob2;

   assert(elem1 != NULL);
   assert(elem2 != NULL);

   prob1 = (PROBLEM*)elem1;
   prob2 = (PROBLEM*)elem2;

   return prob1->lowerbound - prob2->lowerbound;
}

/** initialize component structure */
static
SCIP_RETCODE initComponent(
   PROBLEM*              problem,            /**< subproblem structure */
   int                   index               /**< component index */
   )
{
   COMPONENT* component;
   SCIP* scip;

   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &problem->components[index]) );
   component = problem->components[index];
   assert(component != NULL);

   component->problem = problem;
   component->subscip = NULL;
   component->vars = NULL;
   component->subvars = NULL;
   component->subproblemqueue = NULL;
   component->lastdualbound = -SCIPinfinity(scip);
   component->lastprimalbound = SCIPinfinity(scip);
   component->lastnodelimit = 0LL;
   component->laststatus = SCIP_STATUS_UNKNOWN;
   component->ncalls = 0;
   component->lastsolindex = -1;
   component->nvars = 0;
   component->number = index;

   return SCIP_OKAY;
}

/** free component structure */
static
SCIP_RETCODE freeComponent(
   COMPONENT**           component           /**< pointer to component structure */
   )
{
   PROBLEM* problem;
   SCIP* scip;

   assert(component != NULL);
   assert(*component != NULL);

   problem = (*component)->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   SCIPdebugMessage("freeing component %d\n", (*component)->number);

   if( (*component)->subscip != NULL )
   {
      SCIP_CALL( SCIPfree(&(*component)->subscip) );
   }

   assert(((*component)->vars != NULL) == ((*component)->subvars != NULL));
   if( (*component)->vars != NULL )
   {
      SCIPfreeMemoryArray(scip, &(*component)->vars);
      SCIPfreeMemoryArray(scip, &(*component)->subvars);
   }

   SCIPfreeMemory(scip, component);

   return SCIP_OKAY;
}

/** add subproblem to component structure */
static
SCIP_RETCODE componentAddSubproblem(
   COMPONENT*            component,          /**< component */
   PROBLEM*              problem             /**< subproblem structure */
   )
{
   assert(component != NULL);

   if( component->subproblemqueue == NULL )
   {
      SCIP_CALL( SCIPpqueueCreate(&component->subproblemqueue, 5, 1.2, subproblemSort) );
   }

   SCIP_CALL( SCIPpqueueInsert(component->subproblemqueue, problem) );

   return SCIP_OKAY;
}


static
void SCIPpropComponentsSetComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   COMPONENT*            component           /**< component */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);

   prop = SCIPfindProp(scip, PROP_NAME);
   assert(prop != NULL);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->component = component;
}


/** create the sub-SCIP for a given component
 */
static
SCIP_RETCODE componentCreateSubscip(
   COMPONENT*            component,          /**< pointer to component structure */
   SCIP_HASHMAP*         varmap,             /**< variable hashmap used to improve performance */
   SCIP_HASHMAP*         consmap,            /**< constraint hashmap used to improve performance */
   SCIP_CONS**           conss,              /**< constraints contained in this component */
   int                   nconss,             /**< number of constraints contained in this component */
   SCIP_Bool*            success             /**< pointer to store whether the copying process was successful */
   )
{
   PROBLEM* problem;
   SCIP* scip;
   char name[SCIP_MAXSTRLEN];
   SCIP* subscip;
   SCIP_CONS* newcons;
   int nvars;
   int i;
   int v;

   assert(component != NULL);
   assert(consmap != NULL);
   assert(conss != NULL);

   problem = component->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   nvars = component->nvars;
   assert(nvars > 0);

   (*success) = TRUE;

   SCIP_CALL( SCIPcreate(&component->subscip) );
   subscip = component->subscip;

   /* copy plugins, we omit pricers (because we do not run if there are active pricers) and dialogs */
   SCIP_CALL( SCIPcopyPlugins(scip, subscip, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, success) );

   SCIPpropComponentsSetComponent(subscip, component);

   /* abort if the plugins were not successfully copied */
   if( !(*success) )
      goto TERMINATE;

   /* copy parameter settings */
   SCIP_CALL( SCIPcopyParamSettings(scip, subscip) );

   assert(!SCIPisParamFixed(subscip, "limits/solutions"));
   assert(!SCIPisParamFixed(subscip, "limits/bestsol"));
   assert(!SCIPisParamFixed(subscip, "limits/maxorigsol"));
   assert(!SCIPisParamFixed(subscip, "misc/usevartable"));
   assert(!SCIPisParamFixed(subscip, "misc/useconstable"));
   assert(!SCIPisParamFixed(subscip, "numerics/feastol"));

   /* disable solution limits */
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/solutions", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/maxorigsol", 0) );

   /* reduce the effort spent for hash tables */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/usevartable", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/useconstable", FALSE) );

   /* do not catch control-C */
   //SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output, unless in extended debug mode */
#ifndef SCIP_MORE_DEBUG
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
#endif

   /* create problem in sub-SCIP */
   /* get name of the original problem and add "comp_nr" */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp_%d", problem->name, component->number);
   SCIP_CALL( SCIPcreateProb(subscip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* copy variables */
   for( v = 0; v < component->nvars; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(scip, subscip, component->vars[v], &component->subvars[v], varmap, consmap, FALSE, success) );

      /* abort if variable was not successfully copied */
      if( !(*success) )
         goto TERMINATE;
   }

   for( i = 0; i < nconss; ++i )
   {
      assert(!SCIPconsIsModifiable(conss[i]));

      /* copy the constraint */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", SCIPconsGetName(conss[i]));
      SCIP_CALL( SCIPgetConsCopy(scip, subscip, conss[i], &newcons, SCIPconsGetHdlr(conss[i]), varmap, consmap, name,
            SCIPconsIsInitial(conss[i]), SCIPconsIsSeparated(conss[i]), SCIPconsIsEnforced(conss[i]),
            SCIPconsIsChecked(conss[i]), SCIPconsIsPropagated(conss[i]), FALSE, FALSE,
            SCIPconsIsDynamic(conss[i]), SCIPconsIsRemovable(conss[i]), FALSE, FALSE, success) );

      /* abort if constraint was not successfully copied */
      if( !(*success) )
         goto TERMINATE;

      SCIP_CALL( SCIPaddCons(subscip, newcons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &newcons) );
   }


#if 0 /* for more debugging */
   /* write the problem, if requested */
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s.cip", SCIPgetProbName(subscip));
      SCIPdebugMessage("write problem to file %s\n", name);
      SCIP_CALL( SCIPwriteOrigProblem(subscip, name, NULL, FALSE) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s.set", SCIPgetProbName(subscip));
      SCIPdebugMessage("write settings to file %s\n", name);
      SCIP_CALL( SCIPwriteParams(subscip, name, TRUE, TRUE) );
   }
#endif

   /* In extended debug mode, we want to be informed if the number of variables was reduced during copying.
    * This might happen, since the components propagator uses SCIPgetConsVars() and then SCIPgetActiveVars() to get the
    * active representation, while SCIPgetConsCopy() might use SCIPgetProbvarLinearSum() and this might cancel out some
    * of the active variables and cannot be avoided. However, we want to notice it and check whether the constraint
    * handler could do something more clever.
    */
#ifdef SCIP_MORE_DEBUG
   if( component->nvars > SCIPgetNVars(subscip) )
   {
      SCIPinfoMessage(scip, NULL, "copying component %d reduced number of variables: %d -> %d\n", component->number, component->nvars,
         SCIPgetNVars(subscip));
   }
#endif

 TERMINATE:
   if( !(*success) )
   {
      SCIP_CALL( SCIPfree(&subscip) );
      component->subscip = NULL;

      // only for debugging
      SCIPABORT();
   }


   return SCIP_OKAY;
}

/* forward declaration */
static
SCIP_RETCODE solveProblem(PROBLEM*, SCIP_Bool, SCIP_RESULT*);

/** returns whether problem was solved */
static
SCIP_RETCODE problemIsSolved(
   PROBLEM*              problem             /**< pointer to problem */
   )
{
   return (SCIPpqueueNElems(problem->compqueue) == 0);
}


/** (continues) solving a connected component */
static
SCIP_RETCODE solveComponent(
   COMPONENT*            component,
   SCIP_Bool             lastcomponent,
   SCIP_RESULT*          result              /**< pointer to store the result of the solving process */
   )
{
   PROBLEM* problem;
   SCIP* scip;
   SCIP* subscip;
   SCIP_Longint nodelimit;
   SCIP_Real gaplimit;
   SCIP_Real timelimit;
   SCIP_Real softtimelimit;
   SCIP_Real memorylimit;
   SCIP_STATUS status;

   assert(component != NULL);

   problem = component->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMessage("solve sub-SCIP for component %d (ncalls=%d, absgap=%16.9g)\n",
      component->number, component->ncalls, component->lastprimalbound - component->lastdualbound);

   subscip = component->subscip;
   assert(subscip != NULL);

   /* check if we rather want to continue solving one of the subproblems */
   if( SCIPpqueueNElems(component->subproblemqueue) > 0 )
   {
      PROBLEM* bestproblem;

      bestproblem = (PROBLEM*)SCIPpqueueRemove(component->subproblemqueue);

      /* best subproblems has better lowerbound than main SCIP, solve it */
      if( subscip == NULL || SCIPgetLowerbound(subscip) > bestproblem->lowerbound )
      {
         SCIP_Bool lastproblem = lastcomponent && (component->laststatus == SCIP_STATUS_OPTIMAL) && (SCIPpqueueNElems(component->subproblemqueue) == 0);
         SCIP_CALL( solveProblem(bestproblem, lastproblem, result) );

         /* evaluate result */
         if( !problemIsSolved(bestproblem) )
         {
            SCIPpqueueInsert(component->subproblemqueue, bestproblem);
         }

         return SCIP_OKAY;
      }
   }

   /* update time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   timelimit += SCIPgetSolvingTime(subscip);

   /* update soft time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/softtime", &softtimelimit) );
   if( softtimelimit > -0.5 )
   {
      softtimelimit -= SCIPgetSolvingTime(scip);
      softtimelimit += SCIPgetSolvingTime(subscip);
      softtimelimit = MAX(softtimelimit, 0.0);
   }

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   // todo memory limit
   if( timelimit <= 0.0 )
   {
      SCIPdebugMessage("--> not solved (not enough memory or time left)\n");
      return SCIP_OKAY;
   }

   /* set time and memory limit for the subproblem */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/softtime", softtimelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   if( component->ncalls == 0 )
   {
      nodelimit = 1LL;
      gaplimit = 0.0;
   }
   else if( lastcomponent && SCIPpqueueNElems(component->subproblemqueue) == 0 )
   {
      /* enable output */
      SCIP_CALL( SCIPsetIntParam(component->subscip, "display/verblevel", 4) );

      nodelimit = INT_MAX;
      gaplimit = 0.0;
   }
   else
   {
      nodelimit = 2 * SCIPgetNNodes(component->subscip);
      nodelimit = MAX(nodelimit, 500LL);

      /* set a gap limit of half the current gap (at most 10%) */
      if( SCIPgetGap(component->subscip) < 0.2 )
         gaplimit = 0.5 * SCIPgetGap(component->subscip);
      else
         gaplimit = 0.1;
   }

   /* set gap limit */
   //SCIP_CALL( SCIPsetRealParam(subscip, "limits/gap", gaplimit) );

   /* set node limit */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );

   /* solve the subproblem */
   SCIP_CALL( SCIPsolve(subscip) );

#ifdef SCIP_MORE_DEBUG
   SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

   status = SCIPgetStatus(subscip);

   component->lastnodelimit = nodelimit;
   component->laststatus = status;
   ++component->ncalls;

   SCIPdebugMessage("--> (status=%d, nodes=%lld, time=%.2f): gap: %12.5g%% absgap: %16.9g\n",
      status, SCIPgetNNodes(subscip), SCIPgetSolvingTime(subscip), 100.0*SCIPgetGap(subscip),
      SCIPgetPrimalbound(subscip) - SCIPgetDualbound(subscip));

   *result = SCIP_SUCCESS;

   switch( status )
   {
   case SCIP_STATUS_INFEASIBLE:
      *result = SCIP_CUTOFF;
      break;
   case SCIP_STATUS_UNBOUNDED:
   case SCIP_STATUS_INFORUNBD:
      /* TODO: store unbounded ray in original SCIP data structure */
      *result = SCIP_UNBOUNDED;
      break;
   case SCIP_STATUS_USERINTERRUPT:
      SCIP_CALL( SCIPinterruptSolve(scip) );
   default:
   {
      SCIP_SOL* sol = SCIPgetBestSol(subscip);
      SCIP_VAR* var;
      SCIP_VAR* subvar;
      int v;

      /* update dual bound of problem */
      if( !SCIPisEQ(scip, component->lastdualbound, SCIPgetDualbound(subscip)) )
      {
         SCIP_Real newdualbound = SCIPgetDualbound(subscip);
         assert(!SCIPisInfinity(scip, -newdualbound));

         /* first finite dual bound: decrease inf counter and add dual bound to problem dualbound */
         if( SCIPisInfinity(scip, -component->lastdualbound) )
         {
            --problem->nlowerboundinf;
            problem->lowerbound += newdualbound;
         }
         /* increase problem dual bound by dual bound delta */
         else
         {
            problem->lowerbound += (newdualbound - component->lastdualbound);
         }

         /* update problem dual bound if all problem components have a finite dual bound */
         if( problem->nlowerboundinf == 0 )
         {
            SCIPdebugMessage("component %d: lower bound increased from %16.9g to %16.9g, set global lower bound to %16.9g (gap: %16.9g, absgap: %16.9g)\n",
               component->number, component->lastdualbound, newdualbound, SCIPretransformObj(scip, problem->lowerbound),
               problem->nfeascomps == problem->ncomponents ?
               (SCIPgetSolOrigObj(scip, problem->bestsol) - SCIPretransformObj(scip, problem->lowerbound)) / MAX(ABS(SCIPretransformObj(scip, problem->lowerbound)),SCIPgetSolOrigObj(scip, problem->bestsol)) : SCIPinfinity(scip),
               problem->nfeascomps == problem->ncomponents ?
               SCIPgetSolOrigObj(scip, problem->bestsol) - SCIPretransformObj(scip, problem->lowerbound) : SCIPinfinity(scip));
            SCIP_CALL( SCIPupdateLocalLowerbound(scip, problem->lowerbound) );
         }

         /* store dual bound of this call */
         component->lastdualbound = newdualbound;
      }

      /* update primal solution of problem */
      if( sol != NULL && component->lastsolindex != SCIPsolGetIndex(sol) )
      {
         component->lastsolindex = SCIPsolGetIndex(sol);

         /* increase counter for feasible problems if no solution was known before */
         if( SCIPisInfinity(scip, component->lastprimalbound) )
            ++(problem->nfeascomps);

         /* update working best solution in problem */
         for( v = 0; v < component->nvars; ++v )
         {
            var = component->vars[v];
            subvar = component->subvars[v];
            assert(var != NULL);
            assert(subvar != NULL);

            SCIP_CALL( SCIPsetSolVal(scip, problem->bestsol, var,
                  SCIPgetSolVal(subscip, sol, subvar)) );
         }

         /* if we have a feasible solution for each component, add the working solution to the main problem */
         if( problem->nfeascomps == problem->ncomponents )
         {
            SCIP_Bool feasible;

            SCIP_CALL( SCIPcheckSolOrig(scip, problem->bestsol, &feasible, TRUE, TRUE) );

            SCIP_CALL( SCIPaddSol(scip, problem->bestsol, &feasible) );

            SCIPdebugMessage("component %d: upper bound decreased from %16.9g to %16.9g, new global upper bound of %16.9g (gap: %16.9g, absgap: %16.9g)\n",
               component->number, component->lastprimalbound, SCIPgetPrimalbound(subscip), SCIPgetSolOrigObj(scip, problem->bestsol),
               problem->nfeascomps == problem->ncomponents ?
               (SCIPgetSolOrigObj(scip, problem->bestsol) - SCIPretransformObj(scip, problem->lowerbound)) / MAX(ABS(SCIPretransformObj(scip, problem->lowerbound)),SCIPgetSolOrigObj(scip, problem->bestsol)) : SCIPinfinity(scip),
               problem->nfeascomps == problem->ncomponents ?
               SCIPgetSolOrigObj(scip, problem->bestsol) - SCIPretransformObj(scip, problem->lowerbound) : SCIPinfinity(scip));
            SCIP_CALL( SCIPupdateLocalLowerbound(scip, problem->lowerbound) );
         }

         /* store primal bound of this call */
         component->lastprimalbound = SCIPgetPrimalbound(subscip);
      }

      /* if the component was solved to optimality, we increase the respective counter and free the subscip */
      if( status == SCIP_STATUS_OPTIMAL )
      {
         ++(problem->nsolvedcomps);

         /* keep copy of solution in problem space? ???????? */

         /* free component */
         SCIP_CALL( SCIPfree(&subscip) );
         component->subscip = NULL;
      }
   }
   }

   return SCIP_OKAY;
}

/** initialize subproblem structure */
static
SCIP_RETCODE initProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBLEM**             problem             /**< pointer to subproblem structure */
   )
{
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPallocMemory(scip, problem) );
   assert(*problem != NULL);

   (*problem)->scip = scip;
   (*problem)->components = NULL;
   (*problem)->lowerbound = 0;
   (*problem)->ncomponents = 0;
   (*problem)->nfeascomps = 0;
   (*problem)->nsolvedcomps = 0;
   (*problem)->nlowerboundinf = 0;

   if( SCIPgetDepth(scip) == 0 )
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   else
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_node_%d", SCIPgetProbName(scip), SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*problem)->name, name, strlen(name)+1) );

   SCIP_CALL( SCIPcreateSol(scip, &(*problem)->bestsol, NULL) );

   return SCIP_OKAY;
}

/** free subproblem structure */
static
SCIP_RETCODE freeProblem(
   PROBLEM**             problem             /**< pointer to problem to free */
   )
{
   SCIP* scip;
   int c;

   assert(problem != NULL);
   assert(*problem != NULL);

   scip = (*problem)->scip;
   assert(scip != NULL);

   if( (*problem)->bestsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &(*problem)->bestsol) );
   }

   for( c = (*problem)->ncomponents - 1; c >= 0; --c )
   {
      SCIP_CALL( freeComponent(&(*problem)->components[c]) );
   }

   if( (*problem)->components != NULL )
   {
      SCIPfreeMemoryArray(scip, &(*problem)->components);
   }

   SCIPpqueueFree(&(*problem)->compqueue);

   SCIPfreeMemoryArray(scip, &(*problem)->name);

   SCIPfreeMemory(scip, problem);

   return SCIP_OKAY;
}

/** ???
 */
static
SCIP_RETCODE splitProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_DIGRAPH*         digraph,
   SCIP_CONS**           conss,              /**< constraints */
   SCIP_VAR**            vars,               /**< variables */
   int*                  varcomponent,       /**< component numbers for the variables */
   int                   nconss,             /**< number of constraints */
   int                   nvars,              /**< number of variables */
   int*                  firstvaridxpercons  /**< array with index of first variable in vars array for each constraint */
   )
{
   PROBLEM* problem;
   COMPONENT* component;
   SCIP_HASHMAP* consmap;
   SCIP_HASHMAP* varmap;
   SCIP_CONS** compconss;
   SCIP_VAR** compvars;
   SCIP_Real* compsize;
   int* permu;
   int* conscomponent;
   SCIP_Bool success = TRUE;
   int ncomponents;
   int ncompconss;
   int nbinvars;
   int nintvars;
   int ncontvars;
   int comp;
   int compvarsstart;
   int compconssstart;
   int nrealcomponents;
   int ncreatedcomps;
   int v;
   int c;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(digraph != NULL);
   assert(conss != NULL);
   assert(vars != NULL);
   assert(firstvaridxpercons != NULL);

   ncomponents = SCIPdigraphGetNComponents(digraph);
   nrealcomponents = 0;
   ncreatedcomps = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &conscomponent, nconss) );

   /* We want to sort the components in increasing complexity (number of discrete variables,
    * integer weighted with factor INTFACTOR, continuous used as tie-breaker).
    * Therefore, we now get the variables for each component, count the different variable types
    * and compute a size as described above. Then, we rename the components
    * such that for i < j, component i has no higher complexity than component j.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &compsize, ncomponents) );
   SCIP_CALL( SCIPallocBufferArray(scip, &permu, ncomponents) );

   /* get number of variables in the components */
   for( c = 0; c < ncomponents; ++c )
   {
      int* cvars;
      int ncvars;

      SCIPdigraphGetComponent(digraph, c, &cvars, &ncvars);
      permu[c] = c;
      nbinvars = 0;
      nintvars = 0;

      for( v = 0; v < ncvars; ++v )
      {
         /* check whether variable is of binary or integer type */
         if( SCIPvarGetType(vars[cvars[v]]) == SCIP_VARTYPE_BINARY )
            nbinvars++;
         else if( SCIPvarGetType(vars[cvars[v]]) == SCIP_VARTYPE_INTEGER )
            nintvars++;
      }
      ncontvars = ncvars - nintvars - nbinvars;
      compsize[c] = ((1000 * (nbinvars + INTFACTOR * nintvars) + (950.0 * ncontvars)/nvars));

      /* collect some statistical information */
      SCIPstatistic( updateStatisticsComp(propdata, nbinvars, nintvars) );

      /* decide if this component should be created ?????????? */
      if( ncvars > 0 )
         ++nrealcomponents;
   }

   /* get permutation of component numbers such that the size of the components is increasing */
   SCIPsortRealInt(compsize, permu, ncomponents);

   /* now, we need the reverse direction, i.e., for each component number, we store its new number
    * such that the components are sorted; for this, we abuse the conscomponent array
    */
   assert(nconss >= ncomponents);

   for( c = 0; c < ncomponents; ++c )
      conscomponent[permu[c]] = c;

   /* for each variable, replace the old component number by the new one */
   for( c = 0; c < nvars; ++c )
      varcomponent[c] = conscomponent[varcomponent[c]];

   SCIPfreeBufferArray(scip, &permu);
   SCIPfreeBufferArray(scip, &compsize);

   if( nrealcomponents <= 1 )
      goto TERMINATE;

   /* do the mapping from calculated components per variable to corresponding
    * constraints and sort the component-arrays for faster finding the
    * actual variables and constraints belonging to one component
    */
   for( c = 0; c < nconss; c++ )
      conscomponent[c] = (firstvaridxpercons[c] == -1 ? -1 : varcomponent[firstvaridxpercons[c]]);

   SCIPsortIntPtr(varcomponent, (void**)vars, nvars);
   SCIPsortIntPtr(conscomponent, (void**)conss, nconss);

   /* init subproblem data structure */
   assert(propdata->problem == NULL);
   SCIP_CALL( initProblem(scip, &problem) );
   propdata->problem = problem;

   SCIP_CALL( SCIPallocMemoryArray(scip, &problem->components, nrealcomponents) );
   SCIP_CALL( SCIPpqueueCreate(&problem->compqueue, (int)(1.1*nrealcomponents), 1.2, componentSort) );
   problem->ncomponents = nrealcomponents;
   problem->nlowerboundinf = nrealcomponents;

   problem->nfeascomps = 0;
   problem->nsolvedcomps = 0;

   /* hashmap mapping from original constraints to constraints in the sub-SCIPs (for performance reasons) */
   SCIP_CALL( SCIPhashmapCreate(&consmap, SCIPblkmem(scip), 10 * nconss) );

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
      for( v = compvarsstart; v < nvars && varcomponent[v] == comp; ++v )
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

      /* decide if this component should be created (see above) ?????????? */
      if( v - compvarsstart > 0 )
      {
         SCIP_CALL( initComponent(problem, ncreatedcomps) );
         assert(problem->components[ncreatedcomps] != NULL);
         component = problem->components[ncreatedcomps];
         ++ncreatedcomps;

         /* get component variables and store them in component structure */
         compvars = &(vars[compvarsstart]);
         component->nvars = v - compvarsstart;
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &component->vars, compvars, component->nvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &component->subvars, component->nvars) );
         SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), 10 * component->nvars) );

         /* get component constraints */
         compconss = &(conss[compconssstart]);
         ncompconss = c - compconssstart;

         assert(ncompconss > 0 || component->nvars == 1);

         SCIPdebugMessage("build sub-SCIP for component %d: %d vars (%d bin, %d int, %d cont), %d conss\n",
            component->number, component->nvars, nbinvars, nintvars, component->nvars - nintvars - nbinvars, ncompconss);
#if 0
         {
            int i;
            for( i = 0; i < component->nvars; ++i )
               printf("var %d: <%s>\n", i, SCIPvarGetName(component->vars[i]));
            for( i = 0; i < ncompconss; ++i )
               printf("cons %d: <%s>\n", i, SCIPconsGetName(compconss[i]));
         }
#endif
         /* build subscip for component */
         SCIP_CALL( componentCreateSubscip(component, varmap, consmap, compconss, ncompconss, &success) );

         SCIP_CALL( SCIPpqueueInsert(problem->compqueue, problem->components[comp]) );

         SCIPhashmapFree(&varmap);

         if( !success )
            break;
      }

      /* update start points in variable and constraint array for next component */
      compvarsstart = v;
      compconssstart = c;
   }

   SCIPhashmapFree(&consmap);
 TERMINATE:
   SCIPfreeBufferArray(scip, &conscomponent);

   return SCIP_OKAY;
}


/** continue solving a problem  */
static
SCIP_RETCODE solveProblem(
   PROBLEM*              problem,
   SCIP_Bool             lastproblem,
   SCIP_RESULT*          result              /**< ??? */
   )
{
   COMPONENT* component;
   SCIP_RESULT subscipresult;

   assert(problem != NULL);

   *result = SCIP_SUCCESS;

   component = (COMPONENT*)SCIPpqueueRemove(problem->compqueue);

   /* continue solving the component */
   SCIP_CALL( solveComponent(component, lastproblem && SCIPpqueueNElems(problem->compqueue) == 0, &subscipresult) );

   if( subscipresult == SCIP_CUTOFF || subscipresult == SCIP_UNBOUNDED )
   {
      *result = subscipresult;
   }
   else if( component->laststatus != SCIP_STATUS_OPTIMAL )
   {
      SCIP_CALL( SCIPpqueueInsert(problem->compqueue, component) );
   }

   return SCIP_OKAY;
}


/** solve a problem by iteratively continuing the solving process of its components */
static
SCIP_RETCODE solveIteratively(
   PROBLEM*              problem,
   SCIP_RESULT*          result              /**< pointer to store the result of the presolving call */
   )
{
   SCIP* scip;

   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   *result = SCIP_SUCCESS;

   while( SCIPpqueueNElems(problem->compqueue) > 0 && (*result) == SCIP_SUCCESS && !SCIPisStopped(scip) )
   {
      SCIP_CALL( solveProblem(problem, TRUE, result) );
   }

   return SCIP_OKAY;
}


/*
 * Statistic methods
 */

#ifdef SCIP_STATISTIC
/** initialize data for statistics */
static
SCIP_RETCODE initStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);

   SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->compspercat, NCATEGORIES) );
   BMSclearMemoryArray(propdata->compspercat, NCATEGORIES);

   propdata->nsinglevars = 0;
   propdata->subsolvetime = 0.0;

   return SCIP_OKAY;
}

/** free data for statistics */
static
void freeStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);

   SCIPfreeMemoryArray(scip, &propdata->compspercat);
}

/** reset data for statistics */
static
void resetStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);

   BMSclearMemoryArray(propdata->compspercat, NCATEGORIES);

   propdata->nsinglevars = 0;
   propdata->subsolvetime = 0.0;
}


/** statistics: categorize the component with the given number of binary and integer variables */
static
void updateStatisticsComp(
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int                   nbinvars,           /**< number of binary variables */
   int                   nintvars            /**< number of integer variables */
   )
{
   int ndiscretevars;
   int i;

   assert(propdata != NULL);

   ndiscretevars = nbinvars + nintvars;

   /* check into which category the component belongs by looking at the number of discrete variables */
   for( i = 0; i < (NCATEGORIES - 1); ++i )
   {
      if( ndiscretevars <= CATLIMITS[i] )
      {
         propdata->compspercat[i]++;
         break;
      }
   }

   /* number of discrete variables greater than all limits, so component belongs to last category */
   if( i == (NCATEGORIES - 1) )
      propdata->compspercat[i]++;
}

/** statistics: increase the number of components with a single variable and no constraints */
static
void updateStatisticsSingleVar(
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(propdata != NULL);

   propdata->nsinglevars++;
}

/** statistics: update the total subproblem solving time */
static
void updateStatisticsSubsolvetime(
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Real             subsolvetime        /**< subproblem solving time to add to the statistics */
   )
{
   assert(propdata != NULL);

   propdata->subsolvetime += subsolvetime;
}

/** print statistics */
static
void printStatistics(
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int i;

   assert(propdata != NULL);

   printf("############\n");
   printf("# Connected Components Propagator Statistics:\n");

   printf("# Categorization:");
   for( i = 0; i < NCATEGORIES - 1; ++i )
   {
      printf("[<= %d: %d]", CATLIMITS[i], propdata->compspercat[i]);
   }
   printf("[> %d: %d]\n", CATLIMITS[NCATEGORIES - 2], propdata->compspercat[NCATEGORIES - 1]);
   printf("# Components without constraints: %d\n", propdata->nsinglevars);
   printf("# Total subproblem solving time: %.2f\n", propdata->subsolvetime);
   printf("############\n");
}
#endif


/*
 * Local methods
 */

/** loop over constraints, get active variables and fill directed graph */
static
SCIP_RETCODE fillDigraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int*                  unfixedvarpos,      /**< mapping from variable problem index to unfixed var index */
   int                   nunfixedvars,       /**< number of unfixed variables */
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

      firstvaridxpercons[c] = -1;

      if( nconsvars > 0 )
      {
         v = 0;
         idx1 = -1;

         while( idx1 == -1 && v < nconsvars )
         {
            idx1 = SCIPvarGetProbindex(consvars[v]);
            assert(idx1 >= 0);
            idx1 = unfixedvarpos[idx1];
            assert(idx1 < nunfixedvars);
            ++v;
         }

         if( idx1 >= 0 )
         {

            /* save index of the first variable for later component assignment */
            firstvaridxpercons[c] = idx1;

            /* create sparse directed graph
             * sparse means, to add only those edges necessary for component calculation
             */
            for(; v < nconsvars; ++v )
            {
               idx2 = SCIPvarGetProbindex(consvars[v]);
               assert(idx2 >= 0);
               idx2 = unfixedvarpos[idx2];
               assert(idx2 < nunfixedvars);

               if( idx2 < 0 )
                  continue;

               /* we add only one directed edge, because the other direction is automatically added for component computation */
               SCIP_CALL( SCIPdigraphAddArc(digraph, idx1, idx2, NULL) );
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}

/** performs propagation by searching for components */
static
SCIP_RETCODE propComponents(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_PROP*            prop,               /**< the propagator itself */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation call */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_CONS** conss;
   SCIP_CONS** tmpconss;
   SCIP_Bool success;
   int nconss;
   int ntmpconss;
   int nvars;
   int ncomponents;
   int c;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* do not run in probing mode */
   if( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* do not run, if not all variables are explicitly known */
   if( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   nvars = SCIPgetNVars(scip);

   /* we do not want to run, if there are no variables left */
   if( nvars == 0 )
      return SCIP_OKAY;

   /* check for a reached timelimit */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPstatistic( resetStatistics(scip, propdata) );

   *result = SCIP_DIDNOTFIND;

   ncomponents = 0;

   /* collect checked constraints for component propagator */
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
      int* unfixedvarpos;
      int* firstvaridxpercons;
      int* varlocks;
      int nunfixedvars = 0;
      int v;

      /* copy variables into a local array */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &firstvaridxpercons, nconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varlocks, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &unfixedvarpos, nvars) );
      BMScopyMemoryArray(vars, SCIPgetVars(scip), nvars);

      /* count number of varlocks for each variable (up + down locks) and multiply it by 2;
       * that value is used as an estimate of the number of arcs incident to the variable's node in the digraph
       * to be safe, we double this value
       */
      for( v = 0; v < nvars; ++v )
      {
         if( SCIPisFeasLT(scip, SCIPvarGetLbLocal(vars[v]), SCIPvarGetUbLocal(vars[v])) )
         {
            assert(nunfixedvars <= v);
            vars[nunfixedvars] = vars[v];
            varlocks[nunfixedvars] = 4 * (SCIPvarGetNLocksDown(vars[v]) + SCIPvarGetNLocksUp(vars[v]));
            unfixedvarpos[v] = nunfixedvars;
            ++nunfixedvars;
         }
         else
            unfixedvarpos[v] = -1;
      }

      if( nunfixedvars > 0 )
      {
         /* create and fill directed graph */
         SCIP_CALL( SCIPdigraphCreate(&digraph, nunfixedvars) );
         SCIP_CALL( SCIPdigraphSetSizes(digraph, varlocks) );
         SCIP_CALL( fillDigraph(scip, digraph, conss, nconss, unfixedvarpos, nunfixedvars, firstvaridxpercons, &success) );

         if( success )
         {
            int* varcomponent;

            SCIP_CALL( SCIPallocBufferArray(scip, &varcomponent, nunfixedvars) );

            /* compute independent components */
            SCIP_CALL( SCIPdigraphComputeUndirectedComponents(digraph, 1, varcomponent, &ncomponents) );
#if 0
            {
               int i;
               for( i = 0; i < nunfixedvars; ++i )
               {
                  printf("var <%s>: comp %d\n", SCIPvarGetName(vars[i]), varcomponent[i]);
               }
            }
#endif
#ifndef NDEBUG
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "prop components found %d undirected components\n", ncomponents);
#else
            SCIPdebugMessage("prop components found %d undirected components\n", ncomponents);
#endif

            if( ncomponents > 1 && SCIPnodeGetNComponents(SCIPgetCurrentNode(scip)) != ncomponents )
            {
               printf("node %lld, depth %d found %d components (was %d):\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)),
                  SCIPgetDepth(scip), ncomponents, SCIPnodeGetNComponents(SCIPgetCurrentNode(scip)));

               SCIPnodeSetNComponents(SCIPgetCurrentNode(scip), ncomponents);

               /* create subproblems from independent components and solve them in dependence of their size */
               SCIP_CALL( splitProblem(scip, propdata, digraph, conss, vars, varcomponent, nconss, nunfixedvars,
                     firstvaridxpercons) );
            }

            SCIPfreeBufferArray(scip, &varcomponent);
         }

         SCIPdigraphFree(&digraph);
      }

      SCIPfreeBufferArray(scip, &unfixedvarpos);
      SCIPfreeBufferArray(scip, &varlocks);
      SCIPfreeBufferArray(scip, &firstvaridxpercons);
      SCIPfreeBufferArray(scip, &vars);
   }

   SCIPfreeBufferArray(scip, &conss);

   if( propdata->problem != NULL )
   {
      SCIP_CALL( solveIteratively(propdata->problem, result) );

      SCIP_CALL( freeProblem(&propdata->problem) );
   }

   /* print statistics */
   SCIPstatistic( printStatistics(propdata) );

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyComponents)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of propagator */
   SCIP_CALL( SCIPincludePropComponents(scip) );

   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeComponents)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPfreeMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}

/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitComponents)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* initialize statistics */
   SCIPstatistic( SCIP_CALL( initStatistics(scip, propdata) ) );

   return SCIP_OKAY;
}

#ifdef SCIP_STATISTIC
/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitComponents)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPstatistic( freeStatistics(scip, propdata) );

   return SCIP_OKAY;
}
#endif


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecComponents)
{  /*lint --e{715}*/
   int nfixedvars;
   SCIP_Bool cutoff;
   SCIP_Bool unbounded;

   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /** @warning Don't run in probing or in repropagation since this can lead to wrong conclusion
    *
    *  do not run if propagation w.r.t. current objective is not allowed
    */
   if( SCIPinProbing(scip) || SCIPinRepropagation(scip) )
      return SCIP_OKAY;

   /* only at the root node do we want to run after the node */
   if( proptiming == SCIP_PROPTIMING_AFTERLPLOOP && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   cutoff = FALSE;
   unbounded = FALSE;
   nfixedvars = 0;

   SCIP_CALL( propComponents(scip, prop, result) );

   /* evaluate propagation result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( unbounded )
      *result = SCIP_UNBOUNDED;
   else if( nfixedvars > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the components propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropComponents(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create components propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );
   propdata->component = NULL;
   propdata->problem = NULL;

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecComponents, propdata) );
   assert(prop != NULL);

   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeComponents) );
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitComponents) );
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyComponents) );
   SCIPstatistic( SCIP_CALL( SCIPsetPropExit(scip, prop, propExitComponents) ) );

   return SCIP_OKAY;
}
