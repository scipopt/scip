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
 * @brief  PADM primal heuristic based on ideas published in the paper
 *         "A Decomposition Heuristic for Mixed-Integer Supply Chain Problems"
 *         by Martin Schmidt, Lars Schewe, and Dieter Weninger
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


#define COUPLINGSIZE          3
#define MAX_ADM_ITERATIONS    12

/*
 * Data structures
 */

/** data related to one problem (see below) */
typedef struct Problem PROBLEM;

/** data related to one component */
typedef struct Block
{
   PROBLEM*              problem;            /**< the problem this component belongs to */
   SCIP*                 subscip;            /**< sub-SCIP representing the component */
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
} BLOCK;

/** data related to one problem */
struct Problem
{
   SCIP*                 scip;               /**< the SCIP instance this problem belongs to */
   SCIP_SOL*             bestsol;            /**< best solution found so far for the problem */
   char*                 name;               /**< name of the problem */
   BLOCK*                blocks;             /**< blocks into which the problem will be divided */
   int                   nblocks;            /**< number of blocks */
};

typedef struct set
{
   int                   size;               /**< size of the set */
   int*                  indexes;            /**< set of indexes */
} SET;

typedef struct indexes
{
   int                   block;
   int                   blockContainingLinkVar;
   int                   linkVarIdx;
   SCIP_VAR*             linkVar;
   SCIP_Real             slackPosCoeff;
   SCIP_VAR*             slackPosVar;
   SCIP_Real             slackNegCoeff;
   SCIP_VAR*             slackNegVar;
   SCIP_CONS*            couplingCons;
} INDEXES;

typedef struct indexes2
{
   int                   blockContainingLinkVar;
   int                   linkVarIdx;
   SCIP_Real             linkVarVal;
   SCIP_VAR*             linkVar;
} INDEXES2;

/** returns TRUE iff both keys are equal */
static
SCIP_DECL_HASHKEYEQ(indexesEqual)
{
   SCIP* scip;
   INDEXES* idx1;
   INDEXES* idx2;

   scip = (SCIP*) userptr;
   idx1 = (INDEXES*) key1;
   idx2 = (INDEXES*) key2;

   if( idx1->block != idx2->block )
      return FALSE;

   if( idx1->blockContainingLinkVar != idx2->blockContainingLinkVar )
      return FALSE;

   if( idx1->linkVarIdx != idx2->linkVarIdx )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(indexesHashval)
{  /*lint --e{715}*/
   INDEXES* idx;
   idx = (INDEXES*) key;

   return SCIPhashFour(idx->block, idx->blockContainingLinkVar, idx->linkVarIdx, idx->linkVarIdx);
}

/** returns TRUE iff both keys are equal */
static
SCIP_DECL_HASHKEYEQ(indexes2Equal)
{
   SCIP* scip;
   INDEXES* idx1;
   INDEXES* idx2;

   scip = (SCIP*) userptr;
   idx1 = (INDEXES*) key1;
   idx2 = (INDEXES*) key2;

   if( idx1->blockContainingLinkVar != idx2->blockContainingLinkVar )
      return FALSE;

   if( idx1->linkVarIdx != idx2->linkVarIdx )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(indexes2Hashval)
{  /*lint --e{715}*/
   INDEXES* idx;
   idx = (INDEXES*) key;

   return SCIPhashTwo(idx->blockContainingLinkVar, idx->linkVarIdx);
}

/** primal heuristic data */
struct SCIP_HeurData
{
};


/*
 * Local methods
 */

/** initialize component structure */
static
SCIP_RETCODE initBlock(
   PROBLEM*              problem             /**< subproblem structure */
   )
{
   BLOCK* block;
   SCIP* scip;

   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   block = &problem->blocks[problem->nblocks];

   block->problem = problem;
   block->subscip = NULL;
   block->vars = NULL;
   block->subvars = NULL;
   block->nvars = 0;
   block->number = problem->nblocks;

   block->slackspos = NULL;
   block->nslackspos = 0;
   block->slacksneg = NULL;
   block->nslacksneg = 0;

   block->couplingcons = NULL;
   block->ncouplingcons = 0;

   ++problem->nblocks;

   return SCIP_OKAY;
}

/** free component structure */
static
SCIP_RETCODE freeBlock(
   BLOCK*                block           /**< pointer to block structure */
   )
{
   PROBLEM* problem;
   SCIP* scip;

   assert(block != NULL);

   problem = block->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   SCIPdebugMsg(scip, "freeing block %d of problem <%s>\n", block->number, block->problem->name);

   assert((block->vars != NULL) == (block->subvars != NULL));
   if( block->vars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &block->vars, block->nvars);
      SCIPfreeBlockMemoryArray(scip, &block->subvars, block->nvars);
   }

   SCIPfreeBufferArray(scip, &block->slackspos);
   SCIPfreeBufferArray(scip, &block->slacksneg);
   SCIPfreeBufferArray(scip, &block->couplingcons);

   block->nslackspos = 0;
   block->nslacksneg = 0;
   block->ncouplingcons = 0;

   if( block->subscip != NULL )
      SCIP_CALL( SCIPfree(&block->subscip) );

   return SCIP_OKAY;
}

/** initialize subproblem structure */
static
SCIP_RETCODE initProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBLEM**             problem,            /**< pointer to subproblem structure */
   int                   nblocks             /**< number of blocks */
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

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*problem)->blocks, nblocks) );

   (*problem)->scip = scip;
   (*problem)->nblocks = 0;

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

   /* free all blocks */
   for( c = (*problem)->nblocks - 1; c >= 0; --c )
   {
      SCIP_CALL( freeBlock(&(*problem)->blocks[c]) );
   }
   if( (*problem)->blocks != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*problem)->blocks, (*problem)->nblocks);
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

/** create the subscip for a given block */
static
SCIP_RETCODE blockCreateSubscip(
   BLOCK*                block,              /**< block structure */
   SCIP_HASHMAP*         consmap,            /**< constraint hashmap used to improve performance */
   SCIP_CONS**           conss,              /**< constraints contained in this block */
   int                   nconss,             /**< number of constraints contained in this block */
   SCIP_Bool*            success             /**< pointer to store whether the copying process was successful */
   )
{
   char name[SCIP_MAXSTRLEN];
   PROBLEM* problem;
   SCIP* scip;
   int minsize;

   assert(block != NULL);
   assert(consmap != NULL);
   assert(conss != NULL);
   assert(success != NULL);

   problem = block->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   (*success) = TRUE;

   SCIP_CALL( createSubscip(scip, &block->subscip) );

   if( block->subscip != NULL )
   {
      /* get name of the original problem and add "comp_nr" */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp_%d", problem->name, block->number);

      SCIP_CALL( copyToSubscip(scip, block->subscip, name, conss, consmap, nconss, success) );

      if( !(*success) )
      {
         SCIP_CALL( SCIPfree(&block->subscip) );
         block->subscip = NULL;
      }
   }
   else
      (*success) = FALSE;

   return SCIP_OKAY;
}

/** create problem structure and split it into blocks */
static
SCIP_RETCODE createAndSplitProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           sortedconss,        /**< array of (checked) constraints sorted by blocks */
   int*                  blockstartsconss,   /**< start points of blocks in sortedconss array */
   int                   nblocks,            /**< number of blocks */
   PROBLEM**             problem             /**< pointer to store problem structure */
   )
{
   BLOCK* block;
   SCIP_HASHMAP* consmap;
   SCIP_CONS** blockconss;
   SCIP_Bool success = TRUE;
   int nblockconss;
   int b;

   /* init subproblem data structure */
   SCIP_CALL( initProblem(scip, problem, nblocks) );
   assert((*problem)->blocks != NULL);

   /* hashmap mapping from original constraints to constraints in the sub-SCIPs (for performance reasons) */
   SCIP_CALL( SCIPhashmapCreate(&consmap, SCIPblkmem(scip), blockstartsconss[nblocks]) );

   /* loop over all blocks */
   for( b = 0; b < nblocks; b++ )
   {
      SCIP_CALL( initBlock(*problem) );
      assert((*problem)->nblocks == b+1);

      block = &(*problem)->blocks[b];

      /* get block constraints */
      blockconss = &(sortedconss[blockstartsconss[b]]);
      nblockconss = blockstartsconss[b + 1] - blockstartsconss[b];

      /* build subscip for component */
      SCIP_CALL( blockCreateSubscip(block, consmap, blockconss, nblockconss, &success) );

      if( !success )
         break;
   }

   SCIPhashmapFree(&consmap);

   if( !success )
   {
      /* free subproblem data structure since not all blocks could be copied */
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
   SCIP_DECOMPSTORE* decompstore;
   SCIP_DECOMP* decomp;
   int nconss;
   int nvars;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   int* varslabels;
   int* conslabels;
   int i;
   int nblocks;
   int block;
   int b;
   int k;
   int* blockstartsconss;
   char name[SCIP_MAXSTRLEN];
   char info[SCIP_MAXSTRLEN];
   int numlinkvars;
   SCIP_VAR** linkvars;
   SET* linkvartoblocks;
   SET* blocktolinkvars;
   int j;
   SCIP_VAR** tmpcouplingvars;
   SCIP_Real* tmpcouplingcoef;
   int aIter;
   int pIter;
   SCIP_Bool solutionsdiffer;
   int increasedslacks;
   SCIP_Bool solved;
   int nentries;
   INDEXES* idxlist;
   int idxlistfill;
   INDEXES2* idxlist2;
   int idxlist2fill;
   SCIP_HASHTABLE* htable;
   SCIP_HASHTABLE* htable2;
   SCIP_Real absgap;
   SCIP_Bool doScaling;
   SCIP_Real maxslack;
   SCIP_Real slackthreshold;
   SCIP_STATUS status;

   decompstore = SCIPgetDecompstore(scip);
   doScaling = FALSE;
   absgap = 2.0;

   SCIPdebugMsg(scip,"Initialize padm heuristic\n");
#if 0
   SCIP_CALL( SCIPwriteOrigProblem(scip, "debug_padm.lp", "lp", FALSE) );
#endif

   /* currently only support for one original decomp */
   assert(SCIPdecompstoreGetNOrigDecomps(decompstore) == 1);
   decomp = SCIPdecompstoreGetOrigDecomps(decompstore)[0];
   assert(decomp != NULL);

   nconss = SCIPgetNOrigConss(scip);
   conss = SCIPgetOrigConss(scip);
   nvars = SCIPgetNOrigVars(scip);
   vars = SCIPgetOrigVars(scip);
   for( i = 0; i < nconss; i++ )
      SCIPdebugPrintCons(scip, conss[i], NULL);

   /* determine threshold for penalty coefficients via maximum norm */
   slackthreshold = SCIP_REAL_MIN;
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real obj;
      obj = SCIPvarGetObj(vars[i]);
      obj = REALABS(obj);
      if( SCIPisGT(scip,obj,slackthreshold) )
         slackthreshold = obj;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blockstartsconss, nconss + 1) );

   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);
   for( i = 0; i < nconss; i++ )
      SCIPdebugMsg(scip,"%s %d\n",SCIPconsGetName(conss[i]), conslabels[i]);

   SCIPdecompComputeVarsLabels(scip, decomp, conss, nconss);
   SCIPdecompGetVarsLabels(decomp, vars, varslabels, nvars);
   for( i = 0; i < nvars; i++ )
      SCIPdebugMsg(scip,"%s %d\n",SCIPvarGetName(vars[i]), varslabels[i]);

   /* sort constraints by blocks */
   nblocks = SCIPdecompGetNBlocks(decomp);
   SCIPsortIntPtr(conslabels, (void**)conss, nconss);
   assert(conslabels[0] == 0); /* TODO: currently we do not allow linking constraints */

   /* determine start indices of blocks in sorted conss array */
   i = 0;
   for( b = 0; b < nblocks + 1; ++b )
   {
      assert(i == nconss || conslabels[i] >= b);
      blockstartsconss[b] = i;
      while( i < nconss && conslabels[i] == b )
         ++i;
   }

   SCIP_CALL( createAndSplitProblem(scip, conss, blockstartsconss, nblocks, &problem) );
   assert(nblocks == problem->nblocks);

#if 1
   for( b = 0; b < problem->nblocks; b++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_block_%d.lp", SCIPgetProbName(scip), b);
      SCIP_CALL( SCIPwriteOrigProblem((problem->blocks[b]).subscip, name, "lp", FALSE) );
   }
#endif

   /* count linking variables */
   numlinkvars = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( varslabels[i] == -1 )
         numlinkvars++;
   }
   SCIPdebugMsg(scip,"%d linking variables\n");
   SCIP_CALL( SCIPallocBufferArray(scip, &linkvartoblocks, numlinkvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocktolinkvars, problem->nblocks) );

   /* extract linking variables and init linking variable to blocks set */
   SCIP_CALL( SCIPallocBufferArray(scip, &linkvars, numlinkvars) );
   b = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( varslabels[i] == -1 )
      {
         linkvars[b] = vars[i];
         SCIP_CALL( SCIPallocBufferArray(scip, &(linkvartoblocks[b].indexes), problem->nblocks) );
         linkvartoblocks[b].size = 0;
         b++;
      }
   }

   /* fill linking variable to blocks set */
   for( i = 0; i < numlinkvars; i++ )
   {
      SCIP_VAR* var;
      const char* vname;
      vname = SCIPvarGetName(linkvars[i]);
      k = 0;
      for( b = 0; b < problem->nblocks; b++ )
      {
         var = SCIPfindVar((problem->blocks[b]).subscip, vname);
         if( var != NULL )
         {
            linkvartoblocks[i].indexes[k] = b;
            linkvartoblocks[i].size = k + 1;
            k++;
         }
      }
   }

   /* init block to linking variables set */
   for( b = 0; b < problem->nblocks; b++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(blocktolinkvars[b].indexes), numlinkvars) );
      blocktolinkvars[b].size = 0;
   }

   /* fill block to linking variables set */
   for( b = 0; b < problem->nblocks; b++ )
   {
      k = 0;
      for( i = 0; i < numlinkvars; i++ )
      {
         SCIP_VAR* var;
         const char* vname;
         vname = SCIPvarGetName(linkvars[i]);
         var = SCIPfindVar((problem->blocks[b]).subscip,vname);
         if( var != NULL )
         {
            blocktolinkvars[b].indexes[k] = i;
            blocktolinkvars[b].size = k + 1;
            k++;
         }
      }
   }

   /* init arrays for slack variables and coupling constraints */
   for( b = 0; b < problem->nblocks; b++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &((problem->blocks[b]).slackspos), blocktolinkvars[b].size) );
      SCIP_CALL( SCIPallocBufferArray(scip, &((problem->blocks[b]).slacksneg), blocktolinkvars[b].size) );
      SCIP_CALL( SCIPallocBufferArray(scip, &((problem->blocks[b]).couplingcons), blocktolinkvars[b].size) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcouplingvars, COUPLINGSIZE) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcouplingcoef, COUPLINGSIZE) );
   tmpcouplingcoef[0] = 1.0;
   tmpcouplingcoef[1] = 1.0;
   tmpcouplingcoef[2] = -1.0;

   /* count hashtable entries */
   nentries = 0;
   for( b = 0; b < problem->nblocks; b++ )
   {
      for( i = 0; i < blocktolinkvars[b].size; i++ )
      {
         int linkvaridx = blocktolinkvars[b].indexes[i];
         for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
         {
            if( linkvartoblocks[linkvaridx].indexes[k] != b )
               nentries++;
         }
      }
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &idxlist, nentries) );
   SCIP_CALL( SCIPhashtableCreate(&htable, SCIPblkmem(scip), nentries, SCIPhashGetKeyStandard, indexesEqual, indexesHashval, (void*) scip) );
   idxlistfill = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &idxlist2, nentries) );
   SCIP_CALL( SCIPhashtableCreate(&htable2, SCIPblkmem(scip), nentries, SCIPhashGetKeyStandard, indexes2Equal, indexes2Hashval, (void*) scip) );
   idxlist2fill = 0;

   /* extend submips */
   SCIPdebugMsg(scip,"Extending block models\n");
   for( b = 0; b < problem->nblocks; b++ )
   {
      SCIP_VAR** blockvars;
      int nblockvars;
      blockvars = SCIPgetVars((problem->blocks[b]).subscip);
      nblockvars = SCIPgetNVars((problem->blocks[b]).subscip);

      /* set objective function of each block to zero */
      for( i = 0; i < nblockvars; i++ )
         SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, blockvars[i], 0.0) );

      /* add two slack variables for each linking variable in block */
      for( i = 0; i < blocktolinkvars[b].size; i++ )
      {
         int linkvaridx;
         linkvaridx = blocktolinkvars[b].indexes[i];

         for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
         {
            int blockcontaininglinkvar;
            blockcontaininglinkvar = linkvartoblocks[linkvaridx].indexes[k];

            /* handle different blocks with common linking variable */
            if( blockcontaininglinkvar != b )
            {
               INDEXES idx;
               idx = idxlist[idxlistfill];
               idxlistfill++;
               idx.block = b;
               idx.blockContainingLinkVar = blockcontaininglinkvar;
               idx.linkVarIdx = linkvaridx;
               idx.linkVar = SCIPfindVar((problem->blocks[b]).subscip, SCIPvarGetName(linkvars[linkvaridx]));

               /* create positive slack variable */
               SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_slackpos_block_%d",
                  SCIPvarGetName(linkvars[linkvaridx]), blockcontaininglinkvar);
               j = (problem->blocks[b]).nslackspos;
               (problem->blocks[b]).slackspos[j] = NULL;
               SCIP_CALL( SCIPcreateVarBasic((problem->blocks[b]).subscip,
                     &((problem->blocks[b]).slackspos[j]), name,
                        0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
               SCIP_CALL( SCIPaddVar((problem->blocks[b]).subscip, (problem->blocks[b]).slackspos[j]) );
               assert( (problem->blocks[b]).slackspos[j] != NULL );
               idx.slackPosVar = (problem->blocks[b]).slackspos[j];
               idx.slackPosCoeff = 1.0;
               (problem->blocks[b]).nslackspos++;

               /* create negative slack variable */
               SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_slackneg_block_%d",
                  SCIPvarGetName(linkvars[linkvaridx]), blockcontaininglinkvar);
               j = (problem->blocks[b]).nslacksneg;
               (problem->blocks[b]).slacksneg[j] = NULL;
               SCIP_CALL( SCIPcreateVarBasic((problem->blocks[b]).subscip,
                     &((problem->blocks[b]).slacksneg[j]), name,
                        0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
               SCIP_CALL( SCIPaddVar((problem->blocks[b]).subscip, (problem->blocks[b]).slacksneg[j]) );
               assert( (problem->blocks[b]).slacksneg[j] != NULL );
               idx.slackNegVar = (problem->blocks[b]).slacksneg[j];
               idx.slackNegCoeff = 1.0;
               (problem->blocks[b]).nslacksneg++;

               /* insert slack variables into hashtable */
               idx.couplingCons = NULL;
               SCIP_CALL( SCIPhashtableSafeInsert(htable, (void*)&idx) );
            }
         }
      }

      assert((problem->blocks[b]).nslackspos == blocktolinkvars[b].size);
      assert((problem->blocks[b]).nslacksneg == blocktolinkvars[b].size);

      /* add linking constraint with slack variables to block */
      for( i = 0; i < blocktolinkvars[b].size; i++ )
      {
         int linkvaridx;
         linkvaridx = blocktolinkvars[b].indexes[i];

         for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
         {
            int blockcontaininglinkvar;
            blockcontaininglinkvar = linkvartoblocks[linkvaridx].indexes[k];

            if( blockcontaininglinkvar != b )
            {
               INDEXES idx;
               INDEXES* idxout;
               INDEXES2 idx2;

               idx.block = b;
               idx.blockContainingLinkVar = blockcontaininglinkvar;
               idx.linkVarIdx = linkvaridx;

               /* use in the second list the linking variable of the connected block */
               idx2 = idxlist2[idxlist2fill];
               idxlist2fill++;
               idx2.blockContainingLinkVar = blockcontaininglinkvar;
               idx2.linkVarIdx = linkvaridx;
               idx2.linkVarVal = 0.0;
               idx2.linkVar = SCIPfindVar((problem->blocks[blockcontaininglinkvar]).subscip, SCIPvarGetName(linkvars[linkvaridx]));

               /* fill variables for linking constraint */
               tmpcouplingvars[0] = idx.linkVar;
               tmpcouplingvars[1] = idx.slackPosVar;
               tmpcouplingvars[2] = idx.slackNegVar;

               /* create linking constraint */
               SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_coupling_block_%d",
                  SCIPvarGetName(linkvars[linkvaridx]), blockcontaininglinkvar);
               j = (problem->blocks[b]).ncouplingcons;
               (problem->blocks[b]).couplingcons[j] = NULL;
               SCIP_CALL( SCIPcreateConsBasicLinear((problem->blocks[b]).subscip, &((problem->blocks[b]).couplingcons[j]),
                     name, COUPLINGSIZE, tmpcouplingvars, tmpcouplingcoef, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons((problem->blocks[b]).subscip, (problem->blocks[b]).couplingcons[j]) );
               assert((problem->blocks[b]).couplingcons[j] != NULL);
               (problem->blocks[b]).ncouplingcons++;

               idxout = (INDEXES*)SCIPhashtableRetrieve(htable,(void*)&idx);
               idxout->couplingCons = (problem->blocks[b]).couplingcons[j];

               SCIP_CALL( SCIPhashtableSafeInsert(htable2, (void*)&idx2) );
            }
         }
      }

      assert(blocktolinkvars[b].size == (problem->blocks[b]).ncouplingcons);
   }

#if 0
   for( b = 0; b < problem->nblocks; b++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_block_%d.lp", SCIPgetProbName(scip), b);
      SCIP_CALL( SCIPwriteOrigProblem((problem->blocks[b]).subscip, name, "lp", FALSE) );
   }
#endif

   SCIPdebugMsg(scip,"Starting padm iterations\n");
   SCIPdebugMsg(scip,"PIt\tADMIt\tSlacks\tInfo\n");

   aIter = 0;
   pIter = 0;
   solutionsdiffer = FALSE;
   increasedslacks = 0;
   (void) SCIPsnprintf(info, SCIP_MAXSTRLEN, " ");
   solved = FALSE;

   /* Penalty loop */
   while (!solved)
   {
      pIter++;
      solutionsdiffer = TRUE;
      aIter = 0;

      /*  Alternating direction method loop */
      while (solutionsdiffer && aIter < MAX_ADM_ITERATIONS)
      {
         aIter++;
         solutionsdiffer = FALSE;
         SCIPdebugMsg(scip,"%d\t%d\t%d\t%s",pIter,aIter,increasedslacks,info);

         for( b = 0; b < problem->nblocks; b++ )
         {
            for( i = 0; i < blocktolinkvars[b].size; i++ )
            {
               int linkvaridx;
               linkvaridx = blocktolinkvars[b].indexes[i];

               for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
               {
                  int blockcontaininglinkvar;
                  blockcontaininglinkvar = linkvartoblocks[linkvaridx].indexes[k];

                  if( blockcontaininglinkvar != b )
                  {
                     INDEXES idx;
                     INDEXES* idxout;
                     INDEXES2 idx2;
                     INDEXES2* idx2out;

                     SCIP_CONS* couplingcons;
                     SCIP_Real oldRhs;
                     SCIP_Real newRhs;

                     SCIP_VAR* slackPosVar;
                     SCIP_VAR* slackNegVar;
                     SCIP_VARDATA* vardata;

                     idx.block = b;
                     idx.blockContainingLinkVar = blockcontaininglinkvar;
                     idx.linkVarIdx = linkvaridx;
                     idxout = (INDEXES*)SCIPhashtableRetrieve(htable,(void*)&idx);
                     couplingcons = idxout->couplingCons;

                     idx2.blockContainingLinkVar = blockcontaininglinkvar;
                     idx2.linkVarIdx = linkvaridx;
                     idx2out = (INDEXES2*)SCIPhashtableRetrieve(htable2,(void*)&idx2);
                     oldRhs = SCIPgetRhsLinear(scip, couplingcons);
                     newRhs = idx2out->linkVarVal;

                     /* change side of coupling constraint equation */
                     SCIP_CALL( SCIPchgLhsLinear((problem->blocks[b]).subscip, couplingcons, newRhs) );
                     SCIP_CALL( SCIPchgRhsLinear((problem->blocks[b]).subscip, couplingcons, newRhs) );

                     /* change penalty coefficients of slack variables */
                     SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, idxout->slackPosVar, idxout->slackPosCoeff) );
                     SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, idxout->slackNegVar, idxout->slackNegCoeff) );
                  }
               }
            }

            /* increase slack penalty coeffs until each subproblem can be solved to optimality */
            do
            {
               SCIP_CALL( SCIPsetRealParam((problem->blocks[b]).subscip, "limits/absgap", absgap) );
#if 1
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "debug_padm_block_%d.lp", b);
               SCIP_CALL( SCIPwriteOrigProblem(scip, name, "lp", FALSE) );
#endif
               /* solve block */
               SCIPsolve((problem->blocks[b]).subscip);
               status = SCIPgetStatus((problem->blocks[b]).subscip);
               if( status  == SCIP_STATUS_INFEASIBLE )
               {
                  SCIPdebugMsg(scip,"infeasible subproblem\n");
                  assert(0);
                  /* TODO: stop */
               }
               else if( status == SCIP_STATUS_UNBOUNDED )
               {
                  for( i = 0; i < blocktolinkvars[b].size; i++ )
                  {
                     int linkvaridx;
                     linkvaridx = blocktolinkvars[b].indexes[i];

                     for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
                     {
                        int blockcontaininglinkvar;
                        blockcontaininglinkvar = linkvartoblocks[linkvaridx].indexes[k];

                        if( blockcontaininglinkvar != b )
                        {
                           INDEXES idx;
                           INDEXES* idxout;
                           SCIP_VAR* slackPosVar;
                           SCIP_VAR* slackNegVar;
                           SCIP_VARDATA* vardata;

                           idx.block = b;
                           idx.blockContainingLinkVar = blockcontaininglinkvar;
                           idx.linkVarIdx = linkvaridx;
                           idxout = (INDEXES*)SCIPhashtableRetrieve(htable,(void*)&idx);

                           /* increase penalty coefficients to obtain a bounded subproblem */
                           idxout->slackPosCoeff *= 10.0;
                           idxout->slackNegCoeff *= 10.0;
                           SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, idxout->slackPosVar, idxout->slackPosCoeff) );
                           SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, idxout->slackNegVar, idxout->slackNegCoeff) );
                        }
                     }
                  }
               }

               SCIP_CALL( SCIPfreeTransform((problem->blocks[b]).subscip) );

            } while( status != SCIP_STATUS_OPTIMAL );

            /* check if solutions differ */
            for( i = 0; i < blocktolinkvars[b].size; i++ )
            {
               SCIP_SOL* sol;
               SCIP_Real val;
               SCIP_VAR* var;
               int linkvaridx;
               INDEXES2 idx2;
               INDEXES2* idx2out;
               idx2.blockContainingLinkVar = b;
               idx2.linkVarIdx = linkvaridx;
               idx2out = (INDEXES2*)SCIPhashtableRetrieve(htable2,(void*)&idx2);

               sol = SCIPgetBestSol((problem->blocks[b]).subscip);
               var = idx2out->linkVar;
               val = SCIPgetSolVal((problem->blocks[b]).subscip, sol, var);
               if( !EPSEQ(idx2out->linkVarVal, val, SCIP_DEFAULT_EPSILON) )
                  solutionsdiffer = TRUE;

               idx2out->linkVarVal = val;
            }
         }
      }

      /* check wether problem has been solved and if not update penalty coeffs */
      doScaling = FALSE;
      solved = TRUE;
      increasedslacks = 0;
      maxslack = SCIP_REAL_MIN;
      for( b = 0; b < problem->nblocks; b++ )
      {
         for( i = 0; i < blocktolinkvars[b].size; i++ )
         {
            int linkvaridx;
            linkvaridx = blocktolinkvars[b].indexes[i];

            for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
            {
               int blockcontaininglinkvar;
               blockcontaininglinkvar = linkvartoblocks[linkvaridx].indexes[k];

               if( blockcontaininglinkvar != b )
               {
                  SCIP_SOL* sol;
                  INDEXES idx;
                  INDEXES* idxout;
                  SCIP_Real slackPosVal;
                  SCIP_Real slackNegVal;

                  idx.block = b;
                  idx.blockContainingLinkVar = blockcontaininglinkvar;
                  idx.linkVarIdx = linkvaridx;
                  idxout = (INDEXES*)SCIPhashtableRetrieve(htable,(void*)&idx);

                  sol = SCIPgetBestSol((problem->blocks[b]).subscip);
                  slackPosVal = SCIPgetSolVal((problem->blocks[b]).subscip, sol, idxout->slackPosVar);
                  slackNegVal = SCIPgetSolVal((problem->blocks[b]).subscip, sol, idxout->slackNegVar);

                  /* increase penalty coefficient of positive slack variable */
                  if( SCIPisGT(scip,slackPosVal,0.0) )
                  {
                     idxout->slackPosCoeff *= 10.0;

                     if( idxout->slackPosCoeff > slackthreshold )
                        doScaling = TRUE;

                     if( idxout->slackPosCoeff > maxslack )
                        maxslack = idxout->slackPosCoeff;

                     solved = FALSE;
                     increasedslacks++;
                  }

                  /* increase penalty coefficient of positive slack variable */
                  if( SCIPisGT(scip,slackNegVal,0.0) )
                  {
                     idxout->slackNegCoeff *= 10.0;

                     if( idxout->slackNegCoeff > slackthreshold )
                        doScaling = TRUE;

                     if( idxout->slackNegCoeff > maxslack )
                        maxslack = idxout->slackNegCoeff;

                     solved = FALSE;
                     increasedslacks++;
                  }
               }
            }
         }
      }

      /* should sigmoid scaling be applied? */
      if( doScaling )
      {
         SCIP_Real shift;
         SCIP_Real lowestslack;
         SCIP_Real range;
         SCIP_Real offset;
         SCIP_Real flatness;

         shift = maxslack / 2.0;
         lowestslack = 0.1;
         range = 10.0;
         offset = range/2.0 + lowestslack;
         flatness = maxslack/10.0;

         increasedslacks = 0;

         for( b = 0; b < problem->nblocks; b++ )
         {
            for( i = 0; i < blocktolinkvars[b].size; i++ )
            {
               int linkvaridx;
               linkvaridx = blocktolinkvars[b].indexes[i];

               for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
               {
                  int blockcontaininglinkvar;
                  blockcontaininglinkvar = linkvartoblocks[linkvaridx].indexes[k];

                  if( blockcontaininglinkvar != b )
                  {
                     INDEXES idx;
                     INDEXES* idxout;
                     SCIP_Real oldcoeff;

                     idx.block = b;
                     idx.blockContainingLinkVar = blockcontaininglinkvar;
                     idx.linkVarIdx = linkvaridx;
                     idxout = (INDEXES*)SCIPhashtableRetrieve(htable,(void*)&idx);

                     /* scale coefficient of positive slack variable */
                     oldcoeff = idxout->slackPosCoeff;
                     idxout->slackPosCoeff = ((oldcoeff - shift) / (flatness + REALABS(oldcoeff - shift))) * range/2.0 + offset;

                     /* scale coefficient of negative slack variable */
                     oldcoeff = idxout->slackNegCoeff;
                     idxout->slackNegCoeff = ((oldcoeff - shift) / (flatness + REALABS(oldcoeff - shift))) * range/2.0 + offset;
                  }
               }
            }
         }
      }

      /* adapt in some cases the absgap parameter */
      if((aIter == 1 && solutionsdiffer == FALSE) || doScaling )
      {
         SCIP_Real minabsgap = 0.001;
         SCIP_Real newabsgap = MAX(absgap*0.5,minabsgap);

         if(newabsgap >= minabsgap)
         {
            if(doScaling)
               (void) SCIPsnprintf(info, SCIP_MAXSTRLEN, "scale, %f", newabsgap);
            else
               (void) SCIPsnprintf(info, SCIP_MAXSTRLEN, "%f", newabsgap);

            absgap = newabsgap;
         }
      }

      /* free solution process data */
      for( b = 0; b < problem->nblocks; b++ )
         SCIP_CALL( SCIPfreeTransform((problem->blocks[b]).subscip) );
   }

   /* release slack variables and coupling constraints */
   for( b = 0; b < problem->nblocks; b++ )
   {
      for( i = 0; i < blocktolinkvars[b].size; i++ )
      {
         SCIP_CALL( SCIPreleaseVar((problem->blocks[b]).subscip, &((problem->blocks[b]).slackspos[i])) );
         SCIP_CALL( SCIPreleaseVar((problem->blocks[b]).subscip, &((problem->blocks[b]).slacksneg[i])) );
         SCIP_CALL( SCIPreleaseCons((problem->blocks[b]).subscip, &((problem->blocks[b]).couplingcons[i])) );
      }
   }

   /* free memory */
   SCIPhashtableFree(&htable2);
   SCIPfreeBufferArray(scip, &idxlist2);

   SCIPhashtableFree(&htable);
   SCIPfreeBufferArray(scip, &idxlist);

   SCIPfreeBufferArray(scip, &tmpcouplingcoef);
   SCIPfreeBufferArray(scip, &tmpcouplingvars);

   for( b = 0; b < problem->nblocks; b++ )
      SCIPfreeBufferArray(scip, &(blocktolinkvars[b].indexes));

   SCIPfreeBufferArray(scip, &blocktolinkvars);

   for( i = 0; i < numlinkvars; i++ )
      SCIPfreeBufferArray(scip, &(linkvartoblocks[i].indexes));

   SCIPfreeBufferArray(scip, &linkvartoblocks);
   SCIPfreeBufferArray(scip, &linkvars);

   SCIPfreeBufferArray(scip, &blockstartsconss);
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
