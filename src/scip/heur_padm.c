/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define SCIP_DEBUG
/**@file   heur_padm.c
 * @brief  PADM primal heuristic
 * @author Dieter Weninger
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/debug.h"
#include "blockmemshell/memory.h"
#include "scip/cons_linear.h"
#include "scip/heur_padm.h"
#include "scip/heuristics.h"
#include "scip/pub_cons.h"
#include "scip/pub_event.h"
#include "scip/pub_tree.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_select.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_cons.h"
#include "scip/scip_event.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nodesel.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_table.h"
#include "scip/scip_timing.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include "scip/scip_copy.h"
#include "scip/decomp.h"
#include <string.h>


#define HEUR_NAME             "padm"
#define HEUR_DESC             "penalty alternating direction method primal heuristic"
#define HEUR_DISPCHAR         '?'
#define HEUR_PRIORITY         70000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE/*SCIP_HEURTIMING_AFTERNODE*/
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */


#define COUPLINGSIZE         3

/*
 * Data structures
 */

/** data related to one problem (see below) */
typedef struct Problem PROBLEM;

/** data related to one component */
typedef struct Component
{
   PROBLEM*              problem;            /**< the problem this component belongs to */
   SCIP*                 subscip;            /**< sub-SCIP representing the component */
   SCIP_SOL*             workingsol;         /**< working solution for transferring solutions into the sub-SCIP */
   SCIP_VAR**            vars;               /**< variables belonging to this component (in complete problem) */
   SCIP_VAR**            subvars;            /**< variables belonging to this component (in subscip) */
   int                   nvars;              /**< number of variables belonging to this component */
   int                   number;             /**< component number */

   SCIP_VAR**            slackspos;
   int                   nslackspos;
   SCIP_VAR**            slacksneg;
   int                   nslacksneg;

   SCIP_CONS**           couplingcons;
   int                   ncouplingcons;
} COMPONENT;

/** data related to one problem */
struct Problem
{
   SCIP*                 scip;               /**< the SCIP instance this problem belongs to */
   SCIP_SOL*             bestsol;            /**< best solution found so far for the problem */
   char*                 name;               /**< name of the problem */
   COMPONENT*            components;         /**< independent components into which the problem can be divided */
   int                   ncomponents;        /**< number of independent components into which the problem can be divided */
};

typedef struct set
{
   int                   size;
   int*                  indexes;
} SET;

/** primal heuristic data */
struct SCIP_HeurData
{
};


/*
 * Local methods
 */

/** initialize component structure */
static
SCIP_RETCODE initComponent(
   PROBLEM*              problem             /**< subproblem structure */
   )
{
   COMPONENT* component;
   SCIP* scip;

   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   component = &problem->components[problem->ncomponents];

   component->problem = problem;
   component->subscip = NULL;
   component->workingsol = NULL;
   component->vars = NULL;
   component->subvars = NULL;
   component->nvars = 0;
   component->number = problem->ncomponents;

   component->slackspos = NULL;
   component->nslackspos = 0;
   component->slacksneg = NULL;
   component->nslacksneg = 0;

   ++problem->ncomponents;

   return SCIP_OKAY;
}

/** free component structure */
static
SCIP_RETCODE freeComponent(
   COMPONENT*            component           /**< pointer to component structure */
   )
{
   PROBLEM* problem;
   SCIP* scip;

   assert(component != NULL);

   problem = component->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   SCIPdebugMsg(scip, "freeing component %d of problem <%s>\n", component->number, component->problem->name);

   assert((component->vars != NULL) == (component->subvars != NULL));
   if( component->vars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &component->vars, component->nvars);
      SCIPfreeBlockMemoryArray(scip, &component->subvars, component->nvars);
   }
#if 0
   SCIPfreeBufferArray(scip, &component->slackspos);
   SCIPfreeBufferArray(scip, &component->slacksneg);
#endif
   component->nslackspos = 0;
   component->nslacksneg = 0;
   component->ncouplingcons = 0;

   /* free sub-SCIP belonging to this component and the working solution */
   if( component->subscip != NULL )
   {
      if( component->workingsol != NULL )
      {
         SCIP_CALL( SCIPfreeSol(component->subscip, &component->workingsol) );
      }

      SCIP_CALL( SCIPfree(&component->subscip) );
   }

   return SCIP_OKAY;
}

/** initialize subproblem structure */
static
SCIP_RETCODE initProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBLEM**             problem,            /**< pointer to subproblem structure */
   int                   ncomponents         /**< number of independent components */
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

   SCIP_CALL( SCIPallocBlockMemory(scip, problem) );
   assert(*problem != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*problem)->components, ncomponents) );

   (*problem)->scip = scip;
   (*problem)->ncomponents = 0;

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

   /* free best solution */
   if( (*problem)->bestsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &(*problem)->bestsol) );
   }

   /* free all components */
   for( c = (*problem)->ncomponents - 1; c >= 0; --c )
   {
      SCIP_CALL( freeComponent(&(*problem)->components[c]) );
   }
   if( (*problem)->components != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*problem)->components, (*problem)->ncomponents);
   }

   /* free problem name */
   SCIPfreeMemoryArray(scip, &(*problem)->name);

   /* free PROBLEM struct and set the pointer to NULL */
   SCIPfreeBlockMemory(scip, problem);
   *problem = NULL;

   return SCIP_OKAY;
}

/** create a sub-SCIP for the given variables and constraints */
static
SCIP_RETCODE createSubscip(
   SCIP*                 scip,               /**< main SCIP data structure */
   SCIP**                subscip             /**< pointer to store created sub-SCIP */
   )
{
   SCIP_Bool success;

   /* create a new SCIP instance */
   SCIP_CALL( SCIPcreate(subscip) );

   /* copy plugins, we omit pricers (because we do not run if there are active pricers) and dialogs */
   SCIP_CALL( SCIPcopyPlugins(scip, *subscip, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, &success) );

   /* the plugins were successfully copied */
   if( success )
   {
      /* copy parameter settings */
      SCIP_CALL( SCIPcopyParamSettings(scip, *subscip) );

      /* some general settings should not be fixed */
      assert(!SCIPisParamFixed(*subscip, "limits/solutions"));
      assert(!SCIPisParamFixed(*subscip, "limits/bestsol"));
      assert(!SCIPisParamFixed(*subscip, "misc/usevartable"));
      assert(!SCIPisParamFixed(*subscip, "misc/useconstable"));
      assert(!SCIPisParamFixed(*subscip, "numerics/feastol"));
      assert(!SCIPisParamFixed(*subscip, "misc/usesmalltables"));

      /* disable solution limits */
      SCIP_CALL( SCIPsetIntParam(*subscip, "limits/solutions", -1) );
      SCIP_CALL( SCIPsetIntParam(*subscip, "limits/bestsol", -1) );
   }
   else
   {
      SCIP_CALL( SCIPfree(subscip) );
      *subscip = NULL;
   }

   return SCIP_OKAY;
}

/** copies the given variables and constraints to the given sub-SCIP */
static
SCIP_RETCODE copyToSubscip(
   SCIP*                 scip,               /**< source SCIP */
   SCIP*                 subscip,            /**< target SCIP */
   const char*           name,               /**< name for copied problem */
   SCIP_CONS**           conss,              /**< constraint to copy */
   SCIP_HASHMAP*         consmap,            /**< hashmap used for the copy process of constraints */
   int                   nconss,             /**< number of constraints to copy */
   SCIP_Bool*            success             /**< pointer to store whether copying was successful */
   )
{
   SCIP_CONS* newcons;
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(conss != NULL);
   assert(consmap != NULL);
   assert(success != NULL);

   *success = TRUE;

   /* create problem in sub-SCIP */
   SCIP_CALL( SCIPcopyProb(scip, subscip, NULL/*varmap*/, consmap, FALSE, name) );

   /* copy constraints */
   for( i = 0; i < nconss; ++i )
   {
      assert(!SCIPconsIsModifiable(conss[i]));

      /* copy the constraint */
      SCIP_CALL( SCIPgetConsCopy(scip, subscip, conss[i], &newcons, SCIPconsGetHdlr(conss[i]), /*varmap*/NULL, consmap, NULL,
            SCIPconsIsInitial(conss[i]), SCIPconsIsSeparated(conss[i]), SCIPconsIsEnforced(conss[i]),
            SCIPconsIsChecked(conss[i]), SCIPconsIsPropagated(conss[i]), FALSE, FALSE,
            SCIPconsIsDynamic(conss[i]), SCIPconsIsRemovable(conss[i]), FALSE, FALSE, success) );

      /* abort if constraint was not successfully copied */
      if( !(*success) )
         return SCIP_OKAY;

      SCIP_CALL( SCIPaddCons(subscip, newcons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &newcons) );
   }

   return SCIP_OKAY;
}

/** create the sub-SCIP for a given component */
static
SCIP_RETCODE componentCreateSubscip(
   COMPONENT*            component,          /**< component structure */
   SCIP_HASHMAP*         consmap,            /**< constraint hashmap used to improve performance */
   SCIP_CONS**           conss,              /**< constraints contained in this component */
   int                   nconss,             /**< number of constraints contained in this component */
   SCIP_Bool*            success             /**< pointer to store whether the copying process was successful */
   )
{
   char name[SCIP_MAXSTRLEN];
   PROBLEM* problem;
   SCIP* scip;
   int minsize;

   assert(component != NULL);
   assert(consmap != NULL);
   assert(conss != NULL);
   assert(success != NULL);

   problem = component->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   (*success) = TRUE;

   SCIP_CALL( createSubscip(scip, &component->subscip) );

   if( component->subscip != NULL )
   {
      /* get name of the original problem and add "comp_nr" */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp_%d", problem->name, component->number);

      SCIP_CALL( copyToSubscip(scip, component->subscip, name,
            conss, consmap, nconss, success) );

      if( !(*success) )
      {
         SCIP_CALL( SCIPfree(&component->subscip) );
         component->subscip = NULL;
      }
   }
   else
      (*success) = FALSE;

   return SCIP_OKAY;
}

/** create PROBLEM structure for the current node and split it into components */
static
SCIP_RETCODE createAndSplitProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           sortedconss,        /**< array of (checked) constraints sorted by components */
   int*                  compstartsconss,    /**< start points of components in sortedconss array */
   int                   ncomponents,        /**< number of components */
   PROBLEM**             problem             /**< pointer to store problem structure */
   )
{
   COMPONENT* component;
   SCIP_HASHMAP* consmap;
   SCIP_CONS** compconss;
   SCIP_Bool success = TRUE;
   int ncompconss;
   int comp;

   /* init subproblem data structure */
   SCIP_CALL( initProblem(scip, problem, ncomponents) );
   assert((*problem)->components != NULL);

   /* hashmap mapping from original constraints to constraints in the sub-SCIPs (for performance reasons) */
   SCIP_CALL( SCIPhashmapCreate(&consmap, SCIPblkmem(scip), compstartsconss[ncomponents]) );

   /* loop over all components */
   for( comp = 0; comp < ncomponents; comp++ )
   {
      SCIP_CALL( initComponent(*problem) );
      assert((*problem)->ncomponents == comp+1);

      component = &(*problem)->components[comp];

      /* get component constraints */
      compconss = &(sortedconss[compstartsconss[comp]]);
      ncompconss = compstartsconss[comp + 1] - compstartsconss[comp];

      /* build subscip for component */
      SCIP_CALL( componentCreateSubscip(component, consmap, compconss, ncompconss, &success) );

      if( !success )
         break;
   }

   SCIPhashmapFree(&consmap);

   if( !success )
   {
      /* free subproblem data structure since not all component could be copied */
      SCIP_CALL( freeProblem(problem) );
   }

   return SCIP_OKAY;
}



/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_HEURCOPY(heurCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurCopyPADM NULL
#endif

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_HEURFREE(heurFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurFreePADM NULL
#endif


/** initialization method of primal heuristic (called after problem was transformed) */
#if 0
static
SCIP_DECL_HEURINIT(heurInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitPADM NULL
#endif


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitPADM NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolPADM NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolPADM NULL
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecPADM)
{  /*lint --e{715}*/
   PROBLEM* problem;
   SCIP_DECOMPSTORE* decompstore = SCIPgetDecompstore(scip);
   SCIP_DECOMP* decomp = NULL;
   int nconss;
   int nvars;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   int* varslabels;
   int* conslabels;
   int i;
   int nblocks;
   int block;
   int ncomponents;
   int c;
   int k;
   int* compstartsconss;
   char name[SCIP_MAXSTRLEN];
   int numlinkvars;
   SCIP_VAR** linkvars;
   SET* linkvartoblocks;
   SET* blocktolinkvars;
   int j;
   SCIP_VAR** tmpcouplingvars;
   SCIP_Real* tmpcouplingcoef;

   SCIPdebugMsg(scip,"run padm heuristic...\n");
#if 0
   SCIP_CALL( SCIPwriteOrigProblem(scip, "debug_padm.lp", "lp", FALSE) );
#endif

   /* for prove of concept we only use ORIG decomp, conss, vars */
   assert(SCIPdecompstoreGetNOrigDecomps(decompstore) == 1);
   decomp = SCIPdecompstoreGetOrigDecomps(decompstore)[0];
   assert(decomp != NULL);

   nconss = SCIPgetNOrigConss(scip);
   conss = SCIPgetOrigConss(scip);
   nvars = SCIPgetNOrigVars(scip);
   vars = SCIPgetOrigVars(scip);
   for( c = 0; c < nconss; c++ )
      SCIPdebugPrintCons(scip, conss[c], NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &compstartsconss, nconss + 1) );

   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);
   for( i = 0; i < nconss; i++ )
      SCIPdebugMsg(scip,"%s %d\n",SCIPconsGetName(conss[i]), conslabels[i]);

   SCIPdecompComputeVarsLabels(scip, decomp, conss, nconss);
   SCIPdecompGetVarsLabels(decomp, vars, varslabels, nvars);
   for( i = 0; i < nvars; i++ )
      SCIPdebugMsg(scip,"%s %d\n",SCIPvarGetName(vars[i]), varslabels[i]);

   /* sort constraints by blocks */
   nblocks = SCIPdecompGetNBlocks(decomp);
   ncomponents = nblocks;
   SCIPsortIntPtr(conslabels, (void**)conss, nconss);
   assert(conslabels[0] == 0); /* TODO: currently we do not allow linking constraints */

   /* determine start indices of components in (sorted) conss array */
   i = 0;
   for( c = 0; c < nblocks + 1; ++c )
   {
      assert(i == nconss || conslabels[i] >= c);
      compstartsconss[c] = i;
      while( i < nconss && conslabels[i] == c )
         ++i;
   }

   SCIP_CALL( createAndSplitProblem(scip, conss, compstartsconss, ncomponents, &problem) );

#if 0
   for( i = 0; i < problem->ncomponents; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_block_%d.lp", SCIPgetProbName(scip), i);
      SCIP_CALL( SCIPwriteOrigProblem((problem->components[i]).subscip, name, "lp", FALSE) );
   }
#endif

   /* count linking variables */
   numlinkvars = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( varslabels[i] == -1 )
         numlinkvars++;
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &linkvartoblocks, numlinkvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocktolinkvars, problem->ncomponents) );

   /* extract linking variables and init linkvar to block indexes */
   SCIP_CALL( SCIPallocBufferArray(scip, &linkvars, numlinkvars) );
   c = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( varslabels[i] == -1 )
      {
         linkvars[c] = vars[i];
         SCIP_CALL( SCIPallocBufferArray(scip, &(linkvartoblocks[c].indexes), problem->ncomponents) );
         linkvartoblocks[c].size = 0;
         c++;
      }
   }

   /* fill linkvar to blocks */
   for( i = 0; i < numlinkvars; i++ )
   {
      SCIP_VAR* var;
      const char* vname;
      vname = SCIPvarGetName(linkvars[i]);
      k = 0;
      for( c = 0; c < problem->ncomponents; c++ )
      {
         var = SCIPfindVar((problem->components[c]).subscip,vname);
         if( var != NULL )
         {
            linkvartoblocks[i].indexes[k] = c;
            linkvartoblocks[i].size = k + 1;
            k++;
         }
      }
   }

   /* init block to linkvars set */
   for( i = 0; i < problem->ncomponents; i++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(blocktolinkvars[i].indexes), numlinkvars) );
      blocktolinkvars[i].size = 0;
   }

   /* fill block to linkvars */
   for( i = 0; i < problem->ncomponents; i++ )
   {
      k = 0;
      for( c = 0; c < numlinkvars; c++ )
      {
         SCIP_VAR* var;
         const char* vname;
         vname = SCIPvarGetName(linkvars[c]);
         var = SCIPfindVar((problem->components[i]).subscip,vname);
         if( var != NULL )
         {
            blocktolinkvars[i].indexes[k] = c;
            blocktolinkvars[i].size = k + 1;
            k++;
         }
      }
   }

   /* init slack variable arrays */
   for( c = 0; c < problem->ncomponents; c++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &((problem->components[c]).slackspos), blocktolinkvars[c].size) );
      SCIP_CALL( SCIPallocBufferArray(scip, &((problem->components[c]).slacksneg), blocktolinkvars[c].size) );
      SCIP_CALL( SCIPallocBufferArray(scip, &((problem->components[c]).couplingcons), blocktolinkvars[c].size) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcouplingvars, COUPLINGSIZE) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcouplingcoef, COUPLINGSIZE) );
   tmpcouplingcoef[0] = 1.0;
   tmpcouplingcoef[1] = 1.0;
   tmpcouplingcoef[2] = -1.0;

   /* extend submips */
   for( c = 0; c < problem->ncomponents; c++ )
   {
      SCIP_VAR** blockvars;
      int nblockvars;
      blockvars = SCIPgetVars((problem->components[c]).subscip);
      nblockvars = SCIPgetNVars((problem->components[c]).subscip);

      /* set objective function to zero */
      for( i = 0; i < nblockvars; i++ )
         SCIP_CALL( SCIPchgVarObj((problem->components[c]).subscip, blockvars[i], 0.0) );

      /* add two slack variables for each linking variable in block */
      for( i = 0; i < blocktolinkvars[c].size; i++ )
      {
         int linkvaridx;
         linkvaridx = blocktolinkvars[c].indexes[i];

         for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
         {
            int blockcontaininglinkvar;
            blockcontaininglinkvar = linkvartoblocks[linkvaridx].indexes[k];

            if( blockcontaininglinkvar != c )
            {
               SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_slackpos_block_%d",
                  SCIPvarGetName(linkvars[linkvaridx]), blockcontaininglinkvar);

               j = (problem->components[c]).nslackspos;
               (problem->components[c]).slackspos[j] = NULL;
               SCIP_CALL( SCIPcreateVarBasic((problem->components[c]).subscip,
                     &((problem->components[c]).slackspos[j]), name,
                        0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
               SCIP_CALL( SCIPaddVar((problem->components[c]).subscip,
                     (problem->components[c]).slackspos[j]) );
               assert( (problem->components[c]).slackspos[j] != NULL );
               (problem->components[c]).nslackspos++;

               SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_slackneg_block_%d",
                  SCIPvarGetName(linkvars[linkvaridx]), blockcontaininglinkvar);

               j = (problem->components[c]).nslacksneg;
               (problem->components[c]).slacksneg[j] = NULL;
               SCIP_CALL( SCIPcreateVarBasic((problem->components[c]).subscip,
                     &((problem->components[c]).slacksneg[j]), name,
                        0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
               SCIP_CALL( SCIPaddVar((problem->components[c]).subscip,
                     (problem->components[c]).slacksneg[j]) );
               assert( (problem->components[c]).slacksneg[j] != NULL );
               (problem->components[c]).nslacksneg++;
            }
         }
      }

      assert((problem->components[c]).nslackspos == blocktolinkvars[c].size);
      assert((problem->components[c]).nslacksneg == blocktolinkvars[c].size);

      /* add linking constraint in block */
      for( i = 0; i < blocktolinkvars[c].size; i++ )
      {
         int linkvaridx;
         linkvaridx = blocktolinkvars[c].indexes[i];

         for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
         {
            int blockcontaininglinkvar;
            blockcontaininglinkvar = linkvartoblocks[linkvaridx].indexes[k];

            if( blockcontaininglinkvar != c )
            {
               /* fill variables for linking constraint */
               tmpcouplingvars[0] = SCIPfindVar((problem->components[c]).subscip, SCIPvarGetName(linkvars[linkvaridx]));

               SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_slackpos_block_%d",
                  SCIPvarGetName(linkvars[linkvaridx]), blockcontaininglinkvar);
               tmpcouplingvars[1] = SCIPfindVar((problem->components[c]).subscip, name);

               SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_slackneg_block_%d",
                  SCIPvarGetName(linkvars[linkvaridx]), blockcontaininglinkvar);
               tmpcouplingvars[2] = SCIPfindVar((problem->components[c]).subscip, name);

               /* create linking constraint */
               SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_coupling_block_%d",
                  SCIPvarGetName(linkvars[linkvaridx]), blockcontaininglinkvar);

               j = (problem->components[c]).ncouplingcons;
               SCIP_CALL( SCIPcreateConsBasicLinear((problem->components[c]).subscip, &((problem->components[c]).couplingcons[j]),
                     name, COUPLINGSIZE, tmpcouplingvars, tmpcouplingcoef, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons((problem->components[c]).subscip, (problem->components[c]).couplingcons[j]) );
               (problem->components[c]).ncouplingcons++;
            }
         }
      }

      assert(blocktolinkvars[c].size == (problem->components[c]).ncouplingcons);
   }

   /* TODO: the freeing of slackvars is not correct */


#if 1
   for( i = 0; i < problem->ncomponents; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_block_%d.lp", SCIPgetProbName(scip), i);
      SCIP_CALL( SCIPwriteOrigProblem((problem->components[i]).subscip, name, "lp", FALSE) );
   }
#endif

   SCIPfreeBufferArray(scip, &tmpcouplingcoef);
   SCIPfreeBufferArray(scip, &tmpcouplingvars);

   for( i = 0; i < problem->ncomponents; i++ )
      SCIPfreeBufferArray(scip, &(blocktolinkvars[i].indexes));

   SCIPfreeBufferArray(scip, &blocktolinkvars);

   for( i = 0; i < numlinkvars; i++ )
      SCIPfreeBufferArray(scip, &(linkvartoblocks[i].indexes));

   SCIPfreeBufferArray(scip, &linkvartoblocks);
   SCIPfreeBufferArray(scip, &linkvars);
   SCIPfreeBufferArray(scip, &compstartsconss);
   SCIPfreeBufferArray(scip, &conslabels);
   SCIPfreeBufferArray(scip, &varslabels);
   freeProblem(&problem);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the PADM primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurPADM(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create PADM primal heuristic data */
   heurdata = NULL;

   heur = NULL;

/**! [SnippetCodeStyleBlanks] */

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyXyz, heurFreeXyz, heurInitXyz, heurExitXyz, heurInitsolXyz, heurExitsolXyz, heurExecXyz,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecPADM, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyPADM) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreePADM) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitPADM) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitPADM) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolPADM) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolPADM) );
#endif

   /* add xyz primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
