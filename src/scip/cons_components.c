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
/**@file   cons_components.c
 * @brief  constraint handler for handling independent components
 * @author Gerald Gamrath
 *
 * This constraint handler looks for independent components.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_components.h"

#define CONSHDLR_NAME          "components"
#define CONSHDLR_DESC          "independent components constraint handler"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYPROP         TRUE /**< should propagation method be delayed, if other propagators found reductions? */

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_EXHAUSTIVE /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_PROP_TIMING     (SCIP_PROPTIMING_BEFORELP | SCIP_PROPTIMING_AFTERLPLOOP)

#define DEFAULT_MAXDEPTH 2

#define INTFACTOR 10

/*
 * Data structures
 */

typedef struct Problem PROBLEM;

/** data related to one component */
typedef struct Component
{
   PROBLEM*              problem;            /** the problem this component belongs to */
   SCIP*                 subscip;            /** sub-SCIP representing the component */
   SCIP_SOL*             workingsol;         /** working solution for transferring solutions into the sub-SCIP */
   SCIP_VAR**            vars;               /** variables belonging to this component (in complete problem) */
   SCIP_VAR**            subvars;            /** variables belonging to this component (in subscip) */
   SCIP_VAR**            fixedvars;          /** variables in the sub-SCIP which were copied while copying the component's
                                              *  constraints, but do not count to the subvars, because they were locally fixed */
   SCIP_Real             fixedvarsobjsum;    /** objective contribution of all locally fixed variables */
   SCIP_Real             lastdualbound;      /** dual bound after last optimization call for this component */
   SCIP_Real             lastprimalbound;    /** primal bound after last optimization call for this component */
   SCIP_Longint          lastnodelimit;      /** node limit of last optimization call for this component */
   SCIP_STATUS           laststatus;         /** solution status of last optimization call for the sub-SCIP of this component */
   SCIP_Bool             solved;             /** was this component solved already? */
   int                   ncalls;             /** number of optimization calls for this component */
   int                   lastsolindex;       /** index of best solution after last optimization call for this component */
   int                   lastbestsolindex;
   int                   nvars;              /** number of variables belonging to this component */
   int                   nfixedvars;         /** number of fixed variables copied during constraint copying */
   int                   number;             /** component number */
} COMPONENT;

/** data related to one problem */
struct Problem
{
   SCIP*                 scip;               /** the SCIP instance this problem belongs to */
   COMPONENT**           components;         /** independent components into which the problem can be divided */
   SCIP_PQUEUE*          compqueue;          /** priority queue for components */
   SCIP_SOL*             bestsol;            /** best solution found so far for the problem */
   char*                 name;               /** name of the problem */
   SCIP_Real             fixedvarsobjsum;    /** objective contribution of all locally fixed variables */
   SCIP_Real             lowerbound;         /** lower bound of the problem */
   int                   ncomponents;        /** number of independent components into which the problem can be divided */
   int                   nfeascomps;         /** number of components for which a feasible solution was found */
   int                   nsolvedcomps;       /** number of components solved to optimality */
   int                   nlowerboundinf;     /** number of components with lower bound equal to -infinity */

};


/** control parameters */
struct SCIP_ConshdlrData
{
   SCIP_Real minrelsize;
   int maxdepth;
   int minsize;
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

/** forward declaration: free subproblem structure */
static
SCIP_RETCODE freeProblem(PROBLEM**);


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
   component->fixedvars = NULL;
   component->fixedvarsobjsum = 0.0;
   component->lastdualbound = -SCIPinfinity(scip);
   component->lastprimalbound = SCIPinfinity(scip);
   component->lastnodelimit = 0LL;
   component->laststatus = SCIP_STATUS_UNKNOWN;
   component->solved = FALSE;
   component->ncalls = 0;
   component->lastsolindex = -1;
   component->nvars = 0;
   component->nfixedvars = 0;
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

   SCIPdebugMessage("freeing component %d of problem <%s>\n", (*component)->number, (*component)->problem->name);

   assert(((*component)->vars != NULL) == ((*component)->subvars != NULL));
   if( (*component)->vars != NULL )
   {
      SCIPfreeMemoryArray(scip, &(*component)->vars);
      SCIPfreeMemoryArray(scip, &(*component)->subvars);
   }
   if( (*component)->fixedvars != NULL )
   {
      SCIPfreeMemoryArray(scip, &(*component)->fixedvars);
   }

   if( (*component)->subscip != NULL )
   {
      SCIP_CALL( SCIPfreeSol((*component)->subscip, &(*component)->workingsol) );

      SCIP_CALL( SCIPfree(&(*component)->subscip) );
   }

   SCIPfreeMemory(scip, component);

   return SCIP_OKAY;
}

/** create the sub-SCIP for a given component
 */
static
SCIP_RETCODE componentCreateSubscip(
   COMPONENT*            component,          /**< pointer to component structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
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

   /* abort if the plugins were not successfully copied */
   if( !(*success) )
      goto TERMINATE;

   /* copy parameter settings */
   SCIP_CALL( SCIPcopyParamSettings(scip, subscip) );

   assert(!SCIPisParamFixed(subscip, "limits/solutions"));
   assert(!SCIPisParamFixed(subscip, "limits/bestsol"));
   assert(!SCIPisParamFixed(subscip, "misc/usevartable"));
   assert(!SCIPisParamFixed(subscip, "misc/useconstable"));
   assert(!SCIPisParamFixed(subscip, "numerics/feastol"));

   /* disable solution limits */
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/solutions", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", -1) );

   SCIP_CALL( SCIPsetIntParam(subscip, "constraints/" CONSHDLR_NAME "/maxdepth",
         conshdlrdata->maxdepth - SCIPgetDepth(scip)) );

   /* reduce the effort spent for hash tables */
   //SCIP_CALL( SCIPsetBoolParam(subscip, "misc/usevartable", FALSE) );
   //SCIP_CALL( SCIPsetBoolParam(subscip, "misc/useconstable", FALSE) );

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
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(scip, subscip, component->vars[v], &component->subvars[v], varmap, consmap, FALSE, success) );

      /* abort if variable was not successfully copied */
      if( !(*success) )
         goto TERMINATE;
   }
   assert(SCIPgetNVars(subscip) == nvars);

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

   SCIP_CALL( SCIPtransformProb(subscip) );

   SCIP_CALL( SCIPcreateOrigSol(subscip, &(component->workingsol), NULL) );

   if( SCIPgetNOrigVars(subscip) > nvars )
   {
      SCIP_VAR** vars = SCIPgetOrigVars(subscip);
      int nnewvars = SCIPgetNOrigVars(subscip);
      int index = 0;
      int ninactive = 0;

      SCIP_CALL( SCIPallocMemoryArray(scip, &component->fixedvars, nnewvars - nvars) );

      for( v = 0; v < nnewvars; ++v )
      {
         if( SCIPvarGetIndex(vars[v]) >= nvars )
         {
            /* the variable is either locally fixed or could be an inactive variable present in a constraint
             * for which an aggregation constraint linking it to the active variable was created in the subscip
             */
            assert(SCIPisZero(subscip, SCIPvarGetObj(vars[v])) ||
               SCIPisEQ(subscip, SCIPvarGetLbGlobal(vars[v]), SCIPvarGetUbGlobal(vars[v])));

            /* locally fixed variable */
            if( SCIPisEQ(subscip, SCIPvarGetLbGlobal(vars[v]), SCIPvarGetUbGlobal(vars[v])) )
            {
               component->fixedvarsobjsum += SCIPvarGetLbGlobal(vars[v]) * SCIPvarGetObj(vars[v]);

               component->fixedvars[index] = vars[v];
               ++index;

               SCIP_CALL( SCIPsetSolVal(subscip, component->workingsol, vars[v], SCIPvarGetLbGlobal(vars[v])) );
            }
            /* inactive variable */
            else
            {
               ++ninactive;

               assert(SCIPisZero(subscip, SCIPvarGetObj(vars[v])));
               assert(ninactive <= SCIPgetNOrigConss(subscip) - nconss);
            }
         }
#ifndef NDEBUG
         else
            assert(SCIPisLT(subscip, SCIPvarGetLbGlobal(vars[v]), SCIPvarGetUbGlobal(vars[v])));
#endif
      }
      component->nfixedvars = index;
      SCIPdebugMessage("%d locally fixed variables have been copied, objective contribution: %g\n",
         component->nfixedvars, component->fixedvarsobjsum);
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
   if( nvars > SCIPgetNVars(subscip) )
   {
      SCIPinfoMessage(scip, NULL, "copying component %d reduced number of variables: %d -> %d\n", component->number, nvars,
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
SCIP_RETCODE solveProblem(PROBLEM*, SCIP_RESULT*);

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
   SCIP_SOL* bestsol;
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

   subscip = component->subscip;
   assert(subscip != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMessage("solve component <%s> (ncalls=%d, absgap=%.9g)\n",
      SCIPgetProbName(subscip), component->ncalls, component->lastprimalbound - component->lastdualbound);

   bestsol = SCIPgetBestSol(scip);

   /* update best solution of component */
   if( bestsol != NULL && SCIPsolGetIndex(bestsol) != component->lastbestsolindex )
   {
      SCIP_SOL* compsol = component->workingsol;
      SCIP_VAR** vars = component->vars;
      SCIP_VAR** subvars = component->subvars;
      int nvars = component->nvars;
      int v;

      component->lastbestsolindex = SCIPsolGetIndex(bestsol);

      for( v = 0; v < nvars; ++v )
      {
         SCIP_CALL( SCIPsetSolVal(subscip, compsol, subvars[v], SCIPgetSolVal(scip, bestsol, vars[v])) );
      }
#ifndef NDEBUG
      for( v = 0; v < component->nfixedvars; ++v )
      {
         assert(SCIPisEQ(scip, SCIPgetSolVal(subscip, compsol, component->fixedvars[v]),
               SCIPvarGetLbGlobal(component->fixedvars[v])));
      }
#endif

      if( SCIPgetStage(subscip) == SCIP_STAGE_PROBLEM || SCIPisLT(subscip, SCIPgetSolOrigObj(subscip, compsol), SCIPgetPrimalbound(subscip)) )
      {
         SCIP_Bool feasible;

         SCIPdebugMessage("install new solution in component <%s> inherited from problem <%s>: primal bound %.9g --> %.9g\n",
            SCIPgetProbName(subscip), problem->name,
            SCIPgetStage(subscip) == SCIP_STAGE_PROBLEM ? SCIPinfinity(subscip) : SCIPgetPrimalbound(subscip), SCIPgetSolOrigObj(subscip, compsol));

         SCIP_CALL( SCIPcheckSolOrig(subscip, compsol, &feasible, FALSE, FALSE) );
         if( feasible )
         {
            SCIPdebugMessage("... feasible\n");

            SCIP_CALL( SCIPaddSol(subscip, compsol, &feasible) );
         }
         else
         {
            SCIPdebugMessage("... infeasible, update cutoff bound\n");

            SCIP_CALL( SCIPupdateCutoffbound(subscip, SCIPgetSolOrigObj(subscip, compsol)) );
         }
      }
   }

   {
      assert(component->laststatus != SCIP_STATUS_OPTIMAL);

      SCIPdebugMessage("solve sub-SCIP for component <%s> (ncalls=%d, absgap=%16.9g)\n",
         SCIPgetProbName(component->subscip), component->ncalls, component->lastprimalbound - component->lastdualbound);

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
      else
      {
         nodelimit = SCIPgetNNodes(component->subscip) + 1;
         gaplimit = 0.0;
         //nodelimit = 2 * SCIPgetNNodes(component->subscip);
         //nodelimit = MAX(nodelimit, 500LL);

         /* set a gap limit of half the current gap (at most 10%) */
         //if( SCIPgetGap(component->subscip) < 0.2 )
         //   gaplimit = 0.5 * SCIPgetGap(component->subscip);
         //else
         //   gaplimit = 0.1;

         if( lastcomponent )
         {
            int verblevel;

            SCIP_CALL( SCIPgetIntParam(scip, "display/verblevel", &verblevel) );

            if( verblevel >= 4 )
            {
               /* enable output */
               SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 4) );

               gaplimit = 0.0;
            }
         }
      }

      /* set gap limit */
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/gap", gaplimit) );

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
      case SCIP_STATUS_OPTIMAL:
         component->solved = TRUE;
         break;
      case SCIP_STATUS_INFEASIBLE:
         *result = SCIP_CUTOFF;
         component->solved = TRUE;
         break;
      case SCIP_STATUS_UNBOUNDED:
      case SCIP_STATUS_INFORUNBD:
         /* TODO: store unbounded ray in original SCIP data structure */
         *result = SCIP_UNBOUNDED;
         component->solved = TRUE;
         break;
      case SCIP_STATUS_USERINTERRUPT:
         SCIP_CALL( SCIPinterruptSolve(scip) );
      default:
         break;
      }
   }

   /* evaluate call */
   if( *result == SCIP_SUCCESS )
   {
      SCIP_SOL* sol = SCIPgetBestSol(subscip);
      SCIP_VAR* var;
      SCIP_VAR* subvar;
      SCIP_Real newdualbound;
      int v;

      /* get dual bound as the minimum of SCIP dual bound and sub-problems dual bound */
      newdualbound = SCIPgetDualbound(subscip) - component->fixedvarsobjsum;

      /* update dual bound of problem */
      if( !SCIPisEQ(scip, component->lastdualbound, newdualbound) )
      {
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
            SCIPdebugMessage("component <%s>: dual bound increased from %16.9g to %16.9g, new dual bound of problem <%s>: %16.9g (gap: %16.9g, absgap: %16.9g)\n",
               SCIPgetProbName(subscip), component->lastdualbound, newdualbound, problem->name, SCIPretransformObj(scip, problem->lowerbound),
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
            assert(SCIPvarIsActive(var));

            SCIP_CALL( SCIPsetSolVal(scip, problem->bestsol, var, SCIPgetSolVal(subscip, sol, subvar)) );
         }

         /* if we have a feasible solution for each component, add the working solution to the main problem */
         if( problem->nfeascomps == problem->ncomponents )
         {
            SCIP_Bool feasible;

            SCIP_CALL( SCIPcheckSol(scip, problem->bestsol, TRUE, TRUE, TRUE, TRUE, &feasible) );
            assert(feasible);

            SCIP_CALL( SCIPaddSol(scip, problem->bestsol, &feasible) );

            SCIPdebugMessage("component <%s>: primal bound decreased from %16.9g to %16.9g, new primal bound of problem <%s>: %16.9g (gap: %16.9g, absgap: %16.9g)\n",
               SCIPgetProbName(subscip), component->lastprimalbound, SCIPgetPrimalbound(subscip), problem->name, SCIPgetSolOrigObj(scip, problem->bestsol),
               problem->nfeascomps == problem->ncomponents ?
               (SCIPgetSolOrigObj(scip, problem->bestsol) - SCIPretransformObj(scip, problem->lowerbound)) / MAX(ABS(SCIPretransformObj(scip, problem->lowerbound)),SCIPgetSolOrigObj(scip, problem->bestsol)) : SCIPinfinity(scip),
               problem->nfeascomps == problem->ncomponents ?
               SCIPgetSolOrigObj(scip, problem->bestsol) - SCIPretransformObj(scip, problem->lowerbound) : SCIPinfinity(scip));
         }

         /* store primal bound of this call */
         component->lastprimalbound = SCIPgetPrimalbound(subscip) - component->fixedvarsobjsum;
      }

      /* if the component was solved to optimality, we increase the respective counter and free the subscip */
      if( component->laststatus == SCIP_STATUS_OPTIMAL )
      {
         ++(problem->nsolvedcomps);
         component->solved = TRUE;

         /* free working solution and component */
         SCIP_CALL( SCIPfreeSol(subscip, &component->workingsol) );

         SCIP_CALL( SCIPfree(&subscip) );
         component->subscip = NULL;
      }
   }

   return SCIP_OKAY;
}

/** initialize subproblem structure */
static
SCIP_RETCODE initProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBLEM**             problem,            /**< pointer to subproblem structure */
   SCIP_Real             fixedvarsobjsum     /** objective contribution of all locally fixed variables */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(problem != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocMemory(scip, problem) );
   assert(*problem != NULL);

   (*problem)->scip = scip;
   (*problem)->components = NULL;
   (*problem)->lowerbound = fixedvarsobjsum;
   (*problem)->fixedvarsobjsum = fixedvarsobjsum;
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

   for( v = 0; v < nvars; v++ )
   {
      if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(vars[v]), SCIPvarGetUbLocal(vars[v])) )
      {
         SCIP_CALL( SCIPsetSolVal(scip, (*problem)->bestsol, vars[v],
               (SCIPvarGetUbLocal(vars[v]) + SCIPvarGetLbLocal(vars[v]))/2) );
      }
   }

   SCIPdebugMessage("initialized problem <%s>\n", (*problem)->name);

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

   SCIPdebugMessage("freeing problem <%s>\n", (*problem)->name);

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

/** creates and captures a components constraint */
static
SCIP_RETCODE createConsComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   PROBLEM*              problem             /**< problem to be stored in the constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* find the samediff constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("components constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, (SCIP_CONSDATA*)problem,
         FALSE, FALSE, FALSE, FALSE, TRUE,
         TRUE, FALSE, FALSE, FALSE, TRUE) );

   return SCIP_OKAY;
}


/** ???
 */
static
SCIP_RETCODE splitProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         digraph,
   SCIP_CONS**           conss,              /**< constraints */
   SCIP_VAR**            vars,               /**< variables */
   int*                  varcomponent,       /**< component numbers for the variables */
   int                   nconss,             /**< number of constraints */
   int                   nvars,              /**< number of variables */
   int*                  firstvaridxpercons, /**< array with index of first variable in vars array for each constraint */
   SCIP_Real             fixedvarsobjsum,    /**< objective contribution of all locally fixed variables */
   PROBLEM**             problem             /**< created sub-problem, if problem can be split, NULL otherwise */
   )
{
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
   assert(conshdlrdata != NULL);
   assert(digraph != NULL);
   assert(conss != NULL);
   assert(vars != NULL);
   assert(firstvaridxpercons != NULL);
   assert(problem != NULL);

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
   SCIP_CALL( initProblem(scip, problem, fixedvarsobjsum) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(*problem)->components, nrealcomponents) );
   SCIP_CALL( SCIPpqueueCreate(&(*problem)->compqueue, (int)(1.1*nrealcomponents), 1.2, componentSort) );
   (*problem)->ncomponents = nrealcomponents;
   (*problem)->nlowerboundinf = nrealcomponents;
   (*problem)->nfeascomps = 0;
   (*problem)->nsolvedcomps = 0;

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
         SCIP_CALL( initComponent(*problem, ncreatedcomps) );
         assert((*problem)->components[ncreatedcomps] != NULL);
         component = (*problem)->components[ncreatedcomps];
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

         SCIPdebugMessage("build sub-SCIP for component %d of problem <%s>: %d vars (%d bin, %d int, %d cont), %d conss\n",
            component->number, (*problem)->name, component->nvars, nbinvars, nintvars, component->nvars - nintvars - nbinvars, ncompconss);

         {
            int i;
            for( i = 0; i < component->nvars; ++i )
               assert(SCIPvarIsActive(component->vars[i]));
         }
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
         SCIP_CALL( componentCreateSubscip(component, conshdlrdata, varmap, consmap, compconss, ncompconss, &success) );

         SCIP_CALL( SCIPpqueueInsert((*problem)->compqueue, (*problem)->components[comp]) );

         SCIPhashmapFree(&varmap);

         if( !success )
            break;
      }

      /* update start points in variable and constraint array for next component */
      compvarsstart = v;
      compconssstart = c;
   }
   assert(ncreatedcomps == nrealcomponents);

   SCIPhashmapFree(&consmap);
 TERMINATE:
   SCIPfreeBufferArray(scip, &conscomponent);

   if( nrealcomponents > 0 )
   {
      SCIP_CONS* cons;

      assert((*problem) != NULL);

      SCIP_CALL( createConsComponents(scip, &cons, (*problem)->name, (*problem)) );
      SCIP_CALL( SCIPaddConsNode(scip, SCIPgetCurrentNode(scip), cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}


/** continue solving a problem  */
static
SCIP_RETCODE solveProblem(
   PROBLEM*              problem,
   SCIP_RESULT*          result              /**< ??? */
   )
{
   COMPONENT* component;
   SCIP_RESULT subscipresult;

   assert(problem != NULL);

   *result = SCIP_SUCCESS;

   SCIPdebugMessage("solving problem <%s>: %d components left\n", problem->name, SCIPpqueueNElems(problem->compqueue));

   component = (COMPONENT*)SCIPpqueueRemove(problem->compqueue);

   /* continue solving the component */
   SCIP_CALL( solveComponent(component, SCIPpqueueNElems(problem->compqueue) == 0, &subscipresult) );

   if( subscipresult == SCIP_CUTOFF || subscipresult == SCIP_UNBOUNDED )
   {
      *result = subscipresult;
   }
   else if( !component->solved )
   {
      SCIP_CALL( SCIPpqueueInsert(problem->compqueue, component) );
      *result = SCIP_DELAYNODE;
   }
   else if( SCIPpqueueNElems(problem->compqueue) == 0 )
      *result = SCIP_CUTOFF;
   else
      *result = SCIP_DELAYNODE;

   return SCIP_OKAY;
}

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
SCIP_RETCODE findComponents(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< the components constraint handler data */
   PROBLEM**             problem             /**< created sub-problem, if problem can be split, NULL otherwise */
   )
{
   SCIP_CONS** conss;
   SCIP_CONS** tmpconss;
   SCIP_Bool success;
   int nconss;
   int ntmpconss;
   int nvars;
   int ncomponents;
   int c;

   assert(scip != NULL);
   assert(problem != NULL);

   ncomponents = 0;
   nvars = SCIPgetNVars(scip);

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
      SCIP_Real fixedvarsobjsum = 0.0;
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
         if( SCIPisLT(scip, SCIPvarGetLbLocal(vars[v]), SCIPvarGetUbLocal(vars[v])) )
         {
            assert(nunfixedvars <= v);
            vars[nunfixedvars] = vars[v];
            varlocks[nunfixedvars] = 4 * (SCIPvarGetNLocksDown(vars[v]) + SCIPvarGetNLocksUp(vars[v]));
            unfixedvarpos[v] = nunfixedvars;
            ++nunfixedvars;
         }
         else
         {
            unfixedvarpos[v] = -1;
            fixedvarsobjsum += SCIPvarGetObj(vars[v]) * SCIPvarGetLbLocal(vars[v]);
         }
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
#ifndef NDEBUG
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
               "prop components found %d undirected components at node %lld, depth %d\n",
               ncomponents, SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), SCIPgetDepth(scip));
#else
            SCIPdebugMessage("prop components found %d undirected components at node %lld, depth %d\n",
               ncomponents, SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), SCIPgetDepth(scip));
#endif

            if( ncomponents > 1 )
            {
               /* create subproblems from independent components */
               SCIP_CALL( splitProblem(scip, conshdlrdata, digraph, conss, vars, varcomponent, nconss, nunfixedvars,
                     firstvaridxpercons, fixedvarsobjsum, problem) );
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

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyComponents)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrComponents(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(conshdlrFreeComponents)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropComponents)
{  /*lint --e{715}*/
   PROBLEM* problem;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(SCIPconshdlrGetNActiveConss(conshdlr) >= 0);
   assert(SCIPconshdlrGetNActiveConss(conshdlr) <= 1);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTRUN;

   /** @warning Don't run in probing or in repropagation since this can lead to wrong conclusion
    *
    *  do not run if propagation w.r.t. current objective is not allowed
    */
   if( SCIPinProbing(scip) || SCIPinRepropagation(scip) )
      return SCIP_OKAY;

   /* do not run, if not all variables are explicitly known */
   if( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* we do not want to run, if there are no variables left */
   if( SCIPgetNVars(scip) == 0 )
      return SCIP_OKAY;

   /* check for a reached timelimit */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* the components presolver does kind of dual reductions */
   if( !SCIPallowDualReds(scip) )
      return SCIP_OKAY;

   /* only at the root node do we want to run after the node */
   if( proptiming == SCIP_PROPTIMING_AFTERLPLOOP && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   if( SCIPgetDepth(scip) > conshdlrdata->maxdepth )
   {
      assert(SCIPconshdlrGetNActiveConss(conshdlr) == 0);

      return SCIP_OKAY;
   }

   problem = NULL;
   *result = SCIP_DIDNOTFIND;

   if( SCIPconshdlrGetNActiveConss(conshdlr) >= 1 )
   {
      assert(SCIPconshdlrGetNActiveConss(conshdlr) == 1);

      problem = (PROBLEM*)SCIPconsGetData(SCIPconshdlrGetConss(conshdlr)[0]);
   }
   else
   {
      assert(SCIPconshdlrGetNActiveConss(conshdlr) == 0);

      SCIP_CALL( findComponents(scip, conshdlrdata, &problem) );
   }

   do
   {
      if( problem != NULL )
      {
         SCIP_CALL( solveProblem(problem, result) );
      }
   } while( *result == SCIP_DELAYNODE && SCIPgetDepth(scip) == 0 );

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteComponents)
{  /*lint --e{715}*/
   PROBLEM* problem;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   problem = (PROBLEM*)(*consdata);

   SCIP_CALL( freeProblem(&problem) );

   *consdata = NULL;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockComponents)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolComponents)
{  /*lint --e{715}*/
   assert(nconss == 0);

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolComponents)
{  /*lint --e{715}*/
   if( nconss > 0 )
   {
      assert(nconss == 1);

      SCIP_CALL( SCIPdelCons(scip, conss[0]) );
   }

   return SCIP_OKAY;
}


#define consEnfolpComponents NULL
#define consEnfopsComponents NULL
#define consCheckComponents NULL

/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the components constraint handler and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrComponents(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create components propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC, CONSHDLR_ENFOPRIORITY,
         CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpComponents, consEnfopsComponents, consCheckComponents, consLockComponents,
         conshdlrdata) );
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropComponents,
         CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING));

   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, conshdlrFreeComponents) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolComponents) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolComponents) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyComponents, NULL) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteComponents) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxdepth",
         "maximal depth",
         &conshdlrdata->maxdepth, FALSE, DEFAULT_MAXDEPTH, -1, INT_MAX, NULL, NULL) );


   return SCIP_OKAY;
}
