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
 * @author Katrin Halbig
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

#define HEUR_NAME "padm"
#define HEUR_DESC "penalty alternating direction method primal heuristic"
#define HEUR_DISPCHAR '>'
#define HEUR_PRIORITY 70000
#define HEUR_FREQ 1
#define HEUR_FREQOFS 0
#define HEUR_MAXDEPTH -1
#define HEUR_TIMING SCIP_HEURTIMING_BEFORENODE /*SCIP_HEURTIMING_AFTERNODE*/
#define HEUR_USESSUBSCIP TRUE                  /**< does the heuristic use a secondary SCIP instance? */

#define COUPLINGSIZE 3
#define MAX_ADM_ITERATIONS 12

/*
 * Data structures
 */

/** data related to one problem (see below) */
typedef struct Problem PROBLEM;

/** data related to one block */
typedef struct Block
{
   PROBLEM *problem;   /**< the problem this block belongs to */
   SCIP *subscip;      /**< sub-SCIP representing this block */
   SCIP_VAR **subvars; /**< variables belonging to this block (without slack variables) */
   int nsubvars;       /**< number of variables belonging to this block (without slack variables) */
   int number;         /**< component number */

   SCIP_VAR **slackspos; /**< positive slack variables */
   int nslackspos;       /**< number of positive slack variables */
   SCIP_VAR **slacksneg; /**< negative slack variables */
   int nslacksneg;       /**< number of negative slack variables */

   SCIP_CONS **couplingcons; /**< coupling contraints */
   int ncouplingcons;        /**< number of coupling contraints */
} BLOCK;

/** data related to one problem */
struct Problem
{
   SCIP *scip;    /**< the SCIP instance this problem belongs to */
   char *name;    /**< name of the problem */
   BLOCK *blocks; /**< blocks into which the problem will be divided */
   int nblocks;   /**< number of blocks */
};

/** set data structure */
typedef struct set
{
   int size;     /**< size of the set */
   int *indexes; /**< set of indexes */
} SET;

/** data of one linking varaibles related to one block */
typedef struct blockinfo
{
   int block;                  /**< index of this block */
   int otherblock;             /**< index of the other conntected block */
   int linkVarIdx;             /**< linking variable index */
   SCIP_Real linkVarVal;       /**< value of linking variable */
   SCIP_VAR *linkVar;          /**< linking variable */
   SCIP_Real slackPosObjCoeff; /**< penalty coefficient of positive slack variable */
   SCIP_VAR *slackPosVar;      /**< positive slack variable */
   SCIP_Real slackNegObjCoeff; /**< penalty coefficient of negative slack variable */
   SCIP_VAR *slackNegVar;      /**< negative slack variable */
   SCIP_CONS *couplingCons;    /**< coupling contraint (equation) */
} BLOCKINFO;

/** returns TRUE iff both keys are equal */
static SCIP_DECL_HASHKEYEQ(indexesEqual)
{
   SCIP *scip;
   BLOCKINFO *binfo1;
   BLOCKINFO *binfo2;

   scip = (SCIP *)userptr;
   binfo1 = (BLOCKINFO *)key1;
   binfo2 = (BLOCKINFO *)key2;

   if (binfo1->block != binfo2->block)
      return FALSE;

   if (binfo1->otherblock != binfo2->otherblock)
      return FALSE;

   if (binfo1->linkVarIdx != binfo2->linkVarIdx)
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static SCIP_DECL_HASHKEYVAL(indexesHashval)
{ /*lint --e{715}*/
   BLOCKINFO *binfo;
   binfo = (BLOCKINFO *)key;

   return SCIPhashFour(binfo->block, binfo->otherblock, binfo->linkVarIdx, binfo->linkVarIdx);
}

/** primal heuristic data */
struct SCIP_HeurData
{
};

/*
 * Local methods
 */

/** initialize one block */
static SCIP_RETCODE initBlock(
    PROBLEM *problem /**< problem structure */
)
{
   BLOCK *block;
   SCIP *scip;

   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   block = &problem->blocks[problem->nblocks];

   block->problem = problem;
   block->subscip = NULL;
   block->subvars = NULL;
   block->nsubvars = 0;
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
static SCIP_RETCODE freeBlock(
    BLOCK *block /**< pointer to block structure */
)
{
   PROBLEM *problem;
   SCIP *scip;

   assert(block != NULL);

   problem = block->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   SCIPdebugMsg(scip, "freeing block %d of problem <%s>\n", block->number, block->problem->name);

   SCIPfreeBufferArray(scip, &block->subvars);
   SCIPfreeBufferArray(scip, &block->slackspos);
   SCIPfreeBufferArray(scip, &block->slacksneg);
   SCIPfreeBufferArray(scip, &block->couplingcons);

   block->nslackspos = 0;
   block->nslacksneg = 0;
   block->ncouplingcons = 0;

   if (block->subscip != NULL)
      SCIP_CALL(SCIPfree(&block->subscip));

   return SCIP_OKAY;
}

/** initialize subproblem structure */
static SCIP_RETCODE initProblem(
    SCIP *scip,        /**< SCIP data structure */
    PROBLEM **problem, /**< pointer to problem structure */
    int nblocks        /**< number of blocks */
)
{
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(problem != NULL);

   SCIP_CALL(SCIPallocBlockMemory(scip, problem));
   assert(*problem != NULL);

   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(*problem)->blocks, nblocks));

   (*problem)->scip = scip;
   (*problem)->nblocks = 0;

   (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &(*problem)->name, name, strlen(name) + 1));

   SCIPdebugMessage("initialized problem <%s>\n", (*problem)->name);

   return SCIP_OKAY;
}

/** free subproblem structure */
static SCIP_RETCODE freeProblem(
    PROBLEM **problem /**< pointer to problem to free */
)
{
   SCIP *scip;
   int c;

   assert(problem != NULL);
   assert(*problem != NULL);

   scip = (*problem)->scip;
   assert(scip != NULL);

   /* free all blocks */
   for (c = (*problem)->nblocks - 1; c >= 0; --c)
   {
      SCIP_CALL(freeBlock(&(*problem)->blocks[c]));
   }
   if ((*problem)->blocks != NULL)
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
static SCIP_RETCODE createSubscip(
    SCIP *scip,    /**< main SCIP data structure */
    SCIP **subscip /**< pointer to store created sub-SCIP */
)
{
   //SCIP_Bool success;

   /* create a new SCIP instance */
   SCIP_CALL(SCIPcreate(subscip));

   SCIP_CALL(SCIPincludeDefaultPlugins(*subscip));

   /* copy plugins, we omit pricers (because we do not run if there are active pricers) and dialogs */
   //SCIP_CALL( SCIPcopyPlugins(scip, *subscip, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
   //      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, &success) );

   /* the plugins were successfully copied */
   if (TRUE) // if( success )
   {
      /* copy parameter settings */
      //SCIP_CALL( SCIPcopyParamSettings(scip, *subscip) );
      SCIP_CALL(SCIPsetIntParam(*subscip, "display/verblevel", 4));

      SCIP_Real timelimit;
      SCIP_CALL(SCIPgetRealParam(scip, "limits/time", &timelimit));
      SCIP_CALL(SCIPsetRealParam(*subscip, "limits/time", timelimit));

      SCIP_CALL(SCIPsetIntParam(*subscip, "propagating/probing/freq", -1));
      SCIP_CALL(SCIPsetIntParam(*subscip, "propagating/probing/maxruns", 0));

      /* disable solution limits */
      SCIP_CALL(SCIPsetIntParam(*subscip, "limits/solutions", -1));
      SCIP_CALL(SCIPsetIntParam(*subscip, "limits/bestsol", -1));

      /* avoid recursive calls */
      SCIP_CALL(SCIPsetSubscipsOff(*subscip, TRUE));
   }
   else
   {
      SCIP_CALL(SCIPfree(subscip));
      *subscip = NULL;
   }

   return SCIP_OKAY;
}

/** copies the given variables and constraints to the given sub-SCIP */
static SCIP_RETCODE copyToSubscip(
    SCIP *scip,            /**< source SCIP */
    SCIP *subscip,         /**< target SCIP */
    const char *name,      /**< name for copied problem */
    SCIP_CONS **conss,     /**< constraint to copy */
    SCIP_HASHMAP *varmap,  /**< hashmap used for the copy process of variables */
    SCIP_HASHMAP *consmap, /**< hashmap used for the copy process of constraints */
    int nconss,            /**< number of constraints to copy */
    SCIP_Bool *success     /**< pointer to store whether copying was successful */
)
{
   SCIP_CONS *newcons;
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(conss != NULL);
   assert(consmap != NULL);
   assert(success != NULL);

   *success = TRUE;

   /* create problem in sub-SCIP */
   SCIP_CALL(SCIPcopyProb(scip, subscip, varmap, consmap, FALSE, name));

   /* copy constraints */
   for (i = 0; i < nconss; ++i)
   {
      assert(!SCIPconsIsModifiable(conss[i]));

      /* copy the constraint */
      SCIP_CALL(SCIPgetConsCopy(scip, subscip, conss[i], &newcons, SCIPconsGetHdlr(conss[i]), varmap, consmap, NULL,
                                SCIPconsIsInitial(conss[i]), SCIPconsIsSeparated(conss[i]), SCIPconsIsEnforced(conss[i]),
                                SCIPconsIsChecked(conss[i]), SCIPconsIsPropagated(conss[i]), FALSE, FALSE,
                                SCIPconsIsDynamic(conss[i]), SCIPconsIsRemovable(conss[i]), FALSE, FALSE, success));

      /* abort if constraint was not successfully copied */
      if (!(*success))
         return SCIP_OKAY;

      SCIP_CALL(SCIPaddCons(subscip, newcons));
      SCIP_CALL(SCIPreleaseCons(subscip, &newcons));
   }

   return SCIP_OKAY;
}

/** create the subscip for a given block */
static SCIP_RETCODE blockCreateSubscip(
    BLOCK *block,          /**< block structure */
    SCIP_HASHMAP *varmap,  /**< variable hashmap used to improve performance */
    SCIP_HASHMAP *consmap, /**< constraint hashmap used to improve performance */
    SCIP_CONS **conss,     /**< constraints contained in this block */
    int nconss,            /**< number of constraints contained in this block */
    SCIP_Bool *success     /**< pointer to store whether the copying process was successful */
)
{
   char name[SCIP_MAXSTRLEN];
   PROBLEM *problem;
   SCIP *scip;
   int minsize;
   int nsubscipvars;
   SCIP_VAR **subscipvars;
   int i;

   assert(block != NULL);
   assert(varmap != NULL);
   assert(consmap != NULL);
   assert(conss != NULL);
   assert(success != NULL);

   problem = block->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   (*success) = TRUE;

   SCIP_CALL(createSubscip(scip, &block->subscip));

   if (block->subscip != NULL)
   {
      /* get name of the original problem and add "comp_nr" */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp_%d", problem->name, block->number);

      SCIP_CALL(copyToSubscip(scip, block->subscip, name, conss, varmap, consmap, nconss, success));

      /* save variables of subscip (without slack variables) */
      nsubscipvars = SCIPgetNOrigVars(block->subscip);
      subscipvars = SCIPgetOrigVars(block->subscip);
      SCIP_CALL(SCIPallocBufferArray(scip, &(block->subvars), nsubscipvars));
      block->nsubvars = nsubscipvars;
      for (i = 0; i < nsubscipvars; i++)
         block->subvars[i] = subscipvars[i];

      if (!(*success))
      {
         SCIP_CALL(SCIPfree(&block->subscip));
         block->subscip = NULL;
      }
   }
   else
      (*success) = FALSE;

   return SCIP_OKAY;
}

/** create problem structure and split it into blocks */
static SCIP_RETCODE createAndSplitProblem(
    SCIP *scip,              /**< SCIP data structure */
    SCIP_CONS **sortedconss, /**< array of (checked) constraints sorted by blocks */
    int *blockstartsconss,   /**< start points of blocks in sortedconss array */
    int nblocks,             /**< number of blocks */
    PROBLEM **problem        /**< pointer to store problem structure */
)
{
   BLOCK *block;
   SCIP_HASHMAP *varmap;
   SCIP_HASHMAP *consmap;
   SCIP_CONS **blockconss;
   SCIP_Bool success = TRUE;
   int nblockconss;
   int b;

   /* init subproblem data structure */
   SCIP_CALL(initProblem(scip, problem, nblocks));
   assert((*problem)->blocks != NULL);

   /* hashmap mapping from original constraints to constraints in the sub-SCIPs (for performance reasons) */
   SCIP_CALL(SCIPhashmapCreate(&consmap, SCIPblkmem(scip), blockstartsconss[nblocks]));

   /* loop over all blocks */
   for (b = 0; b < nblocks; b++)
   {
      SCIP_CALL(initBlock(*problem));
      assert((*problem)->nblocks == b + 1);

      block = &(*problem)->blocks[b];

      SCIP_CALL(SCIPhashmapCreate(&varmap, SCIPblkmem(scip), SCIPgetNVars(scip)));

      /* get block constraints */
      blockconss = &(sortedconss[blockstartsconss[b]]);
      nblockconss = blockstartsconss[b + 1] - blockstartsconss[b];

      /* build subscip for block */
      SCIP_CALL(blockCreateSubscip(block, varmap, consmap, blockconss, nblockconss, &success));

      SCIPhashmapFree(&varmap);

      if (!success)
         break;
   }

   SCIPhashmapFree(&consmap);

   if (!success)
   {
      /* free subproblem data structure since not all blocks could be copied */
      SCIP_CALL(freeProblem(problem));
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
static SCIP_DECL_HEUREXEC(heurExecPADM)
{ /*lint --e{715}*/
   PROBLEM *problem;
   SCIP_DECOMPSTORE *decompstore;
   SCIP_DECOMP *decomp;
   int nconss;
   int nvars;
   SCIP_VAR **vars;
   SCIP_CONS **conss;
   int *varslabels;
   int *conslabels;
   int i;
   int nblocks;
   int block;
   int b;
   int k;
   int *blockstartsconss;
   char name[SCIP_MAXSTRLEN];
   char info[SCIP_MAXSTRLEN];
   int numlinkvars;
   SCIP_VAR **linkvars;
   SET *linkvartoblocks;
   SET *blocktolinkvars;
   int j;
   SCIP_VAR **tmpcouplingvars;
   SCIP_Real *tmpcouplingcoef;
   int aIter;
   int pIter;
   SCIP_Bool solutionsdiffer;
   int increasedslacks;
   SCIP_Bool solved;
   int nentries;
   BLOCKINFO *blockinfolist;
   int blockinfolistfill;
   SCIP_HASHTABLE *htable;
   SCIP_Real gap;
   SCIP_Bool doScaling;
   SCIP_Real maxslack;
   SCIP_Real slackthreshold;
   SCIP_STATUS status;
   int blockoffset;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   problem = NULL;
   varslabels = NULL;
   conslabels = NULL;
   blockstartsconss = NULL;
   linkvars = NULL;
   linkvartoblocks = NULL;
   numlinkvars = 0;
   blocktolinkvars = NULL;
   tmpcouplingvars = NULL;
   tmpcouplingcoef = NULL;
   blockinfolist = NULL;
   htable = NULL;

   doScaling = FALSE;
   // todo: find good value for gap
   gap = 8.0;
   SCIPdebugMsg(scip, "Initialize padm heuristic\n");
#if 0
   SCIP_CALL( SCIPwriteOrigProblem(scip, "debug_padm.lp", "lp", FALSE) );
#endif

   /* currently only support for one original decomp */
   decompstore = SCIPgetDecompstore(scip);
   if (SCIPdecompstoreGetNOrigDecomps(decompstore) != 1)
      return SCIP_OKAY;

   decomp = SCIPdecompstoreGetOrigDecomps(decompstore)[0];
   assert(decomp != NULL);

   nconss = SCIPgetNOrigConss(scip);
   conss = SCIPgetOrigConss(scip);
   nvars = SCIPgetNOrigVars(scip);
   vars = SCIPgetOrigVars(scip);
#if 0
   for( i = 0; i < nconss; i++ )
      SCIPdebugPrintCons(scip, conss[i], NULL);
#endif

   /* determine threshold for penalty coefficients via maximum norm */
   slackthreshold = SCIP_REAL_MIN;
   for (i = 0; i < nvars; i++)
   {
      SCIP_Real obj;
      obj = SCIPvarGetObj(vars[i]);
      obj = REALABS(obj);
      if (SCIPisGT(scip, obj, slackthreshold))
         slackthreshold = obj;
   }

   SCIP_CALL(SCIPallocBufferArray(scip, &varslabels, nvars));
   SCIP_CALL(SCIPallocBufferArray(scip, &conslabels, nconss));
   SCIP_CALL(SCIPallocBufferArray(scip, &blockstartsconss, nconss + 1));

   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);
#if 0
   for( i = 0; i < nconss; i++ )
      SCIPdebugMsg(scip,"%s %d\n",SCIPconsGetName(conss[i]), conslabels[i]);
#endif

   SCIPdecompComputeVarsLabels(scip, decomp, conss, nconss);
   SCIPdecompGetVarsLabels(decomp, vars, varslabels, nvars);
#if 0
   for( i = 0; i < nvars; i++ )
      SCIPdebugMsg(scip,"%s %d\n",SCIPvarGetName(vars[i]), varslabels[i]);
#endif

   /* sort constraints by blocks */
   nblocks = SCIPdecompGetNBlocks(decomp);
   SCIPsortIntPtr(conslabels, (void **)conss, nconss);
   if (conslabels[0] == -1)
   {
      SCIPdebugMsg(scip, "No support for linking contraints\n");
      goto TERMINATE;
   }
   blockoffset = conslabels[0];
   SCIPdebugMsg(scip, "Block numbering starts from %d\n", blockoffset);

   /* determine start indices of blocks in sorted conss array */
   i = 0;
   for (b = 0; b < nblocks + 1; ++b)
   {
      blockstartsconss[b] = i;
      if (i == nconss)
         break;

      /* request block numbering without holes */
      assert((conslabels[i] - blockoffset - b) == 0 || (conslabels[i] - blockoffset - b) == 1);

      while (i < nconss && (conslabels[i] - blockoffset) == b)
         ++i;
   }

   SCIP_CALL(createAndSplitProblem(scip, conss, blockstartsconss, nblocks, &problem));
   assert(nblocks == problem->nblocks);

#if 0
   for( b = 0; b < problem->nblocks; b++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_block_%d.lp", SCIPgetProbName(scip), b);
      SCIP_CALL( SCIPwriteOrigProblem((problem->blocks[b]).subscip, name, "lp", FALSE) );
   }
#endif

   /* count linking variables */
   numlinkvars = 0;
   for (i = 0; i < nvars; i++)
   {
      if (varslabels[i] == -1)
         numlinkvars++;
   }
   SCIPdebugMsg(scip, "%d linking variables\n", numlinkvars);
   SCIP_CALL(SCIPallocBufferArray(scip, &linkvartoblocks, numlinkvars));
   SCIP_CALL(SCIPallocBufferArray(scip, &blocktolinkvars, problem->nblocks));

   /* extract linking variables and init linking variable to blocks set */
   SCIP_CALL(SCIPallocBufferArray(scip, &linkvars, numlinkvars));
   b = 0;
   for (i = 0; i < nvars; i++)
   {
      if (varslabels[i] == -1)
      {
         linkvars[b] = vars[i];
         SCIP_CALL(SCIPallocBufferArray(scip, &(linkvartoblocks[b].indexes), problem->nblocks));
         linkvartoblocks[b].size = 0;
         b++;
      }
   }

   /* fill linking variable to blocks set */
   for (i = 0; i < numlinkvars; i++)
   {
      SCIP_VAR *var;
      const char *vname;
      vname = SCIPvarGetName(linkvars[i]);
      k = 0;
      for (b = 0; b < problem->nblocks; b++)
      {
         var = SCIPfindVar((problem->blocks[b]).subscip, vname);
         if (var != NULL)
         {
            linkvartoblocks[i].indexes[k] = b;
            linkvartoblocks[i].size = k + 1;
            k++;
         }
      }
   }

   /* init block to linking variables set */
   for (b = 0; b < problem->nblocks; b++)
   {
      SCIP_CALL(SCIPallocBufferArray(scip, &(blocktolinkvars[b].indexes), numlinkvars));
      blocktolinkvars[b].size = 0;
   }

   /* fill block to linking variables set */
   for (b = 0; b < problem->nblocks; b++)
   {
      k = 0;
      for (i = 0; i < numlinkvars; i++)
      {
         SCIP_VAR *var;
         const char *vname;
         vname = SCIPvarGetName(linkvars[i]);
         var = SCIPfindVar((problem->blocks[b]).subscip, vname);
         if (var != NULL)
         {
            blocktolinkvars[b].indexes[k] = i;
            blocktolinkvars[b].size = k + 1;
            k++;
         }
      }
   }

   /* init arrays for slack variables and coupling constraints */
   for (b = 0; b < problem->nblocks; b++)
   {
      SCIP_CALL(SCIPallocBufferArray(scip, &((problem->blocks[b]).slackspos), blocktolinkvars[b].size));
      SCIP_CALL(SCIPallocBufferArray(scip, &((problem->blocks[b]).slacksneg), blocktolinkvars[b].size));
      SCIP_CALL(SCIPallocBufferArray(scip, &((problem->blocks[b]).couplingcons), blocktolinkvars[b].size));
   }

   SCIP_CALL(SCIPallocBufferArray(scip, &tmpcouplingvars, COUPLINGSIZE));
   SCIP_CALL(SCIPallocBufferArray(scip, &tmpcouplingcoef, COUPLINGSIZE));
   tmpcouplingcoef[0] = 1.0;
   tmpcouplingcoef[1] = 1.0;
   tmpcouplingcoef[2] = -1.0;

   /* count hashtable entries */
   nentries = 0;
   for (b = 0; b < problem->nblocks; b++)
   {
      for (i = 0; i < blocktolinkvars[b].size; i++)
      {
         int linkvaridx = blocktolinkvars[b].indexes[i];
         for (k = 0; k < linkvartoblocks[linkvaridx].size; k++)
         {
            if (linkvartoblocks[linkvaridx].indexes[k] != b)
               nentries++;
         }
      }
   }
   SCIP_CALL(SCIPallocBufferArray(scip, &blockinfolist, nentries));
   SCIP_CALL(SCIPhashtableCreate(&htable, SCIPblkmem(scip), 1, SCIPhashGetKeyStandard, indexesEqual, indexesHashval, (void *)scip));
   blockinfolistfill = 0;

   /* extend submips */
   SCIPdebugMsg(scip, "Extending %d block models\n", problem->nblocks);
   for (b = 0; b < problem->nblocks; b++)
   {
      SCIP_VAR **blockvars;
      int nblockvars;
      blockvars = SCIPgetVars((problem->blocks[b]).subscip);
      nblockvars = SCIPgetNVars((problem->blocks[b]).subscip);

      /* set objective function of each block to zero */
      for (i = 0; i < nblockvars; i++)
         SCIP_CALL(SCIPchgVarObj((problem->blocks[b]).subscip, blockvars[i], 0.0));

      /* add two slack variables for each linking variable in block */
      for (i = 0; i < blocktolinkvars[b].size; i++)
      {
         int linkvaridx;
         linkvaridx = blocktolinkvars[b].indexes[i];

         for (k = 0; k < linkvartoblocks[linkvaridx].size; k++)
         {
            int b2;
            b2 = linkvartoblocks[linkvaridx].indexes[k];

            /* handle different blocks with common linking variable */
            if (b2 != b)
            {
               BLOCKINFO *binfo;
               binfo = &blockinfolist[blockinfolistfill];
               blockinfolistfill++;
               binfo->block = b;
               binfo->otherblock = b2;
               binfo->linkVarIdx = linkvaridx;
               binfo->linkVarVal = 0.0;
               binfo->linkVar = SCIPfindVar((problem->blocks[b]).subscip, SCIPvarGetName(linkvars[linkvaridx]));

               /* create positive slack variable */
               SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_slackpos_block_%d", SCIPvarGetName(linkvars[linkvaridx]), b2);
               j = (problem->blocks[b]).nslackspos;
               (problem->blocks[b]).slackspos[j] = NULL;
               SCIP_CALL(SCIPcreateVarBasic((problem->blocks[b]).subscip,
                                            &((problem->blocks[b]).slackspos[j]), name,
                                            0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS));
               SCIP_CALL(SCIPaddVar((problem->blocks[b]).subscip, (problem->blocks[b]).slackspos[j]));
               assert((problem->blocks[b]).slackspos[j] != NULL);
               binfo->slackPosObjCoeff = 1.0;
               binfo->slackPosVar = (problem->blocks[b]).slackspos[j];
               (problem->blocks[b]).nslackspos++;

               /* create negative slack variable */
               SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_slackneg_block_%d", SCIPvarGetName(linkvars[linkvaridx]), b2);
               j = (problem->blocks[b]).nslacksneg;
               (problem->blocks[b]).slacksneg[j] = NULL;
               SCIP_CALL(SCIPcreateVarBasic((problem->blocks[b]).subscip,
                                            &((problem->blocks[b]).slacksneg[j]), name,
                                            0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS));
               SCIP_CALL(SCIPaddVar((problem->blocks[b]).subscip, (problem->blocks[b]).slacksneg[j]));
               assert((problem->blocks[b]).slacksneg[j] != NULL);
               binfo->slackNegObjCoeff = 1.0;
               binfo->slackNegVar = (problem->blocks[b]).slacksneg[j];
               (problem->blocks[b]).nslacksneg++;

               /* fill variables for linking constraint */
               tmpcouplingvars[0] = binfo->linkVar;
               tmpcouplingvars[1] = binfo->slackPosVar;
               tmpcouplingvars[2] = binfo->slackNegVar;

               /* create linking constraint */
               SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_coupling_block_%d",
                            SCIPvarGetName(linkvars[linkvaridx]), b2);
               j = (problem->blocks[b]).ncouplingcons;
               (problem->blocks[b]).couplingcons[j] = NULL;
               SCIP_CALL(SCIPcreateConsBasicLinear((problem->blocks[b]).subscip, &((problem->blocks[b]).couplingcons[j]),
                                                   name, COUPLINGSIZE, tmpcouplingvars, tmpcouplingcoef, 0.0, 0.0));
               SCIP_CALL(SCIPaddCons((problem->blocks[b]).subscip, (problem->blocks[b]).couplingcons[j]));
               assert((problem->blocks[b]).couplingcons[j] != NULL);
               (problem->blocks[b]).ncouplingcons++;
               binfo->couplingCons = (problem->blocks[b]).couplingcons[j];

               /* feed hashtable */
               SCIP_CALL(SCIPhashtableSafeInsert(htable, (void *)binfo));
            }
         }
      }

      assert((problem->blocks[b]).nslackspos >= blocktolinkvars[b].size);
      assert((problem->blocks[b]).nslacksneg >= blocktolinkvars[b].size);
      assert(blocktolinkvars[b].size <= (problem->blocks[b]).ncouplingcons);
   }

   assert(nentries == SCIPhashtableGetNElements(htable));
   SCIPdebugMsg(scip, "Hashtable has %d elements\n", SCIPhashtableGetNElements(htable));

#if 0
   /* write extended submips */
   for( b = 0; b < problem->nblocks; b++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_block_%d.lp", SCIPgetProbName(scip), b);
      SCIP_CALL( SCIPwriteOrigProblem((problem->blocks[b]).subscip, name, "lp", FALSE) );
   }
#endif

   SCIPdebugMsg(scip, "Starting iterations\n");
   SCIPdebugMsg(scip, "PIt\tADMIt\tSlacks\tInfo\n");

   aIter = 0;
   pIter = 0;
   solutionsdiffer = FALSE;
   increasedslacks = 0;
   (void)SCIPsnprintf(info, SCIP_MAXSTRLEN, "-");
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
         SCIPdebugMsg(scip, "%d\t%d\t%d\t%s\n", pIter, aIter, increasedslacks, info);

         for (b = 0; b < problem->nblocks; b++)
         {
            for (i = 0; i < blocktolinkvars[b].size; i++)
            {
               int linkvaridx;
               linkvaridx = blocktolinkvars[b].indexes[i];

               for (k = 0; k < linkvartoblocks[linkvaridx].size; k++)
               {
                  int b2;
                  b2 = linkvartoblocks[linkvaridx].indexes[k];

                  if (b2 != b)
                  {
                     BLOCKINFO binfo;
                     BLOCKINFO *binfoout;
                     BLOCKINFO binfo2;
                     BLOCKINFO *binfo2out;

                     SCIP_CONS *couplingcons;
                     SCIP_Real oldRhs;
                     SCIP_Real newRhs;

                     binfo.block = b;
                     binfo.otherblock = b2;
                     binfo.linkVarIdx = linkvaridx;

                     binfoout = (BLOCKINFO *)SCIPhashtableRetrieve(htable, (void *)&binfo);
                     couplingcons = binfoout->couplingCons;
                     oldRhs = SCIPgetRhsLinear(scip, couplingcons);

                     /* interchange blocks b and b2 for getting new right hand side */
                     binfo2.block = b2;
                     binfo2.otherblock = b;
                     binfo2.linkVarIdx = linkvaridx;
                     binfo2out = (BLOCKINFO *)SCIPhashtableRetrieve(htable, (void *)&binfo2);
                     newRhs = binfo2out->linkVarVal;

                     /* change side of coupling constraint equation with linking variable value of the other block */
                     SCIP_CALL(SCIPchgLhsLinear((problem->blocks[b]).subscip, couplingcons, newRhs));
                     SCIP_CALL(SCIPchgRhsLinear((problem->blocks[b]).subscip, couplingcons, newRhs));

                     /* change penalty coefficients of slack variables */
                     SCIP_CALL(SCIPchgVarObj((problem->blocks[b]).subscip, binfoout->slackPosVar, binfoout->slackPosObjCoeff));
                     SCIP_CALL(SCIPchgVarObj((problem->blocks[b]).subscip, binfoout->slackNegVar, binfoout->slackNegObjCoeff));
                  }
               }
            }

            /* increase slack penalty coeffs until each subproblem can be solved to optimality */
            do
            {
               SCIP_CALL(SCIPsetRealParam((problem->blocks[b]).subscip, "limits/gap", gap));

               /* solve block */
               SCIPsolve((problem->blocks[b]).subscip);
               status = SCIPgetStatus((problem->blocks[b]).subscip);

               if (status == SCIP_STATUS_OPTIMAL || status == SCIP_STATUS_GAPLIMIT)
               {
                  for (i = 0; i < blocktolinkvars[b].size; i++)
                  {
                     int linkvaridx;
                     linkvaridx = blocktolinkvars[b].indexes[i];

                     for (k = 0; k < linkvartoblocks[linkvaridx].size; k++)
                     {
                        int b2;
                        b2 = linkvartoblocks[linkvaridx].indexes[k];

                        if (b2 != b)
                        {
                           SCIP_SOL *sol;
                           BLOCKINFO binfo;
                           BLOCKINFO *binfoout;
                           SCIP_VAR *var;
                           SCIP_Real val;

                           binfo.block = b;
                           binfo.otherblock = b2;
                           binfo.linkVarIdx = linkvaridx;
                           binfoout = (BLOCKINFO *)SCIPhashtableRetrieve(htable, (void *)&binfo);

                           sol = SCIPgetBestSol((problem->blocks[b]).subscip);
                           var = binfoout->linkVar;
                           val = SCIPgetSolVal((problem->blocks[b]).subscip, sol, var);

                           if (!EPSEQ(binfoout->linkVarVal, val, SCIP_DEFAULT_EPSILON))
                              solutionsdiffer = TRUE;

                           binfoout->linkVarVal = val;
                        }
                     }
                  }
               }
               else if (status == SCIP_STATUS_UNBOUNDED)
               {
                  for (i = 0; i < blocktolinkvars[b].size; i++)
                  {
                     int linkvaridx;
                     linkvaridx = blocktolinkvars[b].indexes[i];

                     for (k = 0; k < linkvartoblocks[linkvaridx].size; k++)
                     {
                        int b2;
                        b2 = linkvartoblocks[linkvaridx].indexes[k];

                        if (b2 != b)
                        {
                           BLOCKINFO binfo;
                           BLOCKINFO *binfoout;
                           SCIP_VAR *slackPosVar;
                           SCIP_VAR *slackNegVar;
                           SCIP_VARDATA *vardata;

                           binfo.block = b;
                           binfo.otherblock = b2;
                           binfo.linkVarIdx = linkvaridx;
                           binfoout = (BLOCKINFO *)SCIPhashtableRetrieve(htable, (void *)&binfo);

                           /* increase penalty coefficients to obtain a bounded subproblem */
                           binfoout->slackPosObjCoeff *= 10.0;
                           binfoout->slackNegObjCoeff *= 10.0;
                           SCIP_CALL(SCIPchgVarObj((problem->blocks[b]).subscip, binfoout->slackPosVar, binfoout->slackPosObjCoeff));
                           SCIP_CALL(SCIPchgVarObj((problem->blocks[b]).subscip, binfoout->slackNegVar, binfoout->slackNegObjCoeff));
                        }
                     }
                  }
               }
               else
               {
                  SCIPdebugMsg(scip, "Block solving status %d not supported\n", status);
                  goto TERMINATE;
               }

               SCIP_CALL(SCIPfreeTransform((problem->blocks[b]).subscip));

            } while (status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_GAPLIMIT);
         }
      }

      /* check wether problem has been solved and if not update penalty coeffs */
      doScaling = FALSE;
      solved = TRUE;
      increasedslacks = 0;
      maxslack = SCIP_REAL_MIN;
      for (b = 0; b < problem->nblocks; b++)
      {
         for (i = 0; i < blocktolinkvars[b].size; i++)
         {
            int linkvaridx;
            linkvaridx = blocktolinkvars[b].indexes[i];

            for (k = 0; k < linkvartoblocks[linkvaridx].size; k++)
            {
               int b2;
               b2 = linkvartoblocks[linkvaridx].indexes[k];

               if (b2 != b)
               {
                  SCIP_SOL *sol;
                  BLOCKINFO binfo;
                  BLOCKINFO *binfoout;
                  SCIP_Real slackPosVal;
                  SCIP_Real slackNegVal;

                  binfo.block = b;
                  binfo.otherblock = b2;
                  binfo.linkVarIdx = linkvaridx;
                  binfoout = (BLOCKINFO *)SCIPhashtableRetrieve(htable, (void *)&binfo);

                  sol = SCIPgetBestSol((problem->blocks[b]).subscip);
                  slackPosVal = SCIPgetSolVal((problem->blocks[b]).subscip, sol, binfoout->slackPosVar);
                  slackNegVal = SCIPgetSolVal((problem->blocks[b]).subscip, sol, binfoout->slackNegVar);

                  /* increase penalty coefficient of positive slack variable */
                  if (SCIPisGT(scip, slackPosVal, 0.0))
                  {
                     binfoout->slackPosObjCoeff *= 10.0;

                     if (binfoout->slackPosObjCoeff > slackthreshold)
                        doScaling = TRUE;

                     if (binfoout->slackPosObjCoeff > maxslack)
                        maxslack = binfoout->slackPosObjCoeff;

                     solved = FALSE;
                     increasedslacks++;
                  }

                  /* increase penalty coefficient of negative slack variable */
                  if (SCIPisGT(scip, slackNegVal, 0.0))
                  {
                     binfoout->slackNegObjCoeff *= 10.0;

                     if (binfoout->slackNegObjCoeff > slackthreshold)
                        doScaling = TRUE;

                     if (binfoout->slackNegObjCoeff > maxslack)
                        maxslack = binfoout->slackNegObjCoeff;

                     solved = FALSE;
                     increasedslacks++;
                  }
               }
            }
         }
      }

      /* should sigmoid scaling be applied? */
      if (doScaling)
      {
         SCIP_Real shift;
         SCIP_Real lowestslack;
         SCIP_Real range;
         SCIP_Real offset;
         SCIP_Real flatness;

         shift = maxslack / 2.0;
         lowestslack = 0.1;
         range = 10.0;
         offset = range / 2.0 + lowestslack;
         flatness = maxslack / 10.0;

         increasedslacks = 0;

         for (b = 0; b < problem->nblocks; b++)
         {
            for (i = 0; i < blocktolinkvars[b].size; i++)
            {
               int linkvaridx;
               linkvaridx = blocktolinkvars[b].indexes[i];

               for (k = 0; k < linkvartoblocks[linkvaridx].size; k++)
               {
                  int b2;
                  b2 = linkvartoblocks[linkvaridx].indexes[k];

                  if (b2 != b)
                  {
                     BLOCKINFO binfo;
                     BLOCKINFO *binfoout;
                     SCIP_Real oldcoeff;

                     binfo.block = b;
                     binfo.otherblock = b2;
                     binfo.linkVarIdx = linkvaridx;
                     binfoout = (BLOCKINFO *)SCIPhashtableRetrieve(htable, (void *)&binfo);

                     /* scale coefficient of positive slack variable */
                     oldcoeff = binfoout->slackPosObjCoeff;
                     binfoout->slackPosObjCoeff = ((oldcoeff - shift) / (flatness + REALABS(oldcoeff - shift))) * range / 2.0 + offset;

                     /* scale coefficient of negative slack variable */
                     oldcoeff = binfoout->slackNegObjCoeff;
                     binfoout->slackNegObjCoeff = ((oldcoeff - shift) / (flatness + REALABS(oldcoeff - shift))) * range / 2.0 + offset;
                  }
               }
            }
         }
      }

      /* adapt in some cases the gap parameter */
      if ((aIter == 1 && solutionsdiffer == FALSE) || doScaling)
      {
         SCIP_Real mingap = 0.001; //todo
         SCIP_Real newgap = MAX(gap * 0.5, mingap);

         if (newgap >= mingap)
         {
            if (doScaling)
               (void)SCIPsnprintf(info, SCIP_MAXSTRLEN, "scale, %f", newgap);
            else
               (void)SCIPsnprintf(info, SCIP_MAXSTRLEN, "%f", newgap);

            gap = newgap;
         }
      }

      /* free solution process data */
      if (!solved)
         for (b = 0; b < problem->nblocks; b++)
            SCIP_CALL(SCIPfreeTransform((problem->blocks[b]).subscip));
   }

   /* copy solution if present */
   if (solved)
   {
      SCIP_VAR **vars;
      int nvars;
      SCIP_SOL *newsol;
      SCIP_Bool success;
      SCIP_Real *blocksolvals;
      SCIP_VAR *origobjconstvar;

      assert(increasedslacks == 0);

      SCIP_CALL(SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL));
      SCIP_CALL(SCIPallocBufferArray(scip, &blocksolvals, nvars));
      SCIP_CALL(SCIPcreateSol(scip, &newsol, heur));

      for (b = 0; b < problem->nblocks; b++)
      {
         SCIP_SOL *blocksol;
         SCIP_VAR **blockvars;
         int nblockvars;

         /* get solution of block variables (without slack variables) */
         blocksol = SCIPgetBestSol((problem->blocks[b]).subscip);
         blockvars = (problem->blocks[b]).subvars;
         nblockvars = (problem->blocks[b]).nsubvars;
         SCIP_CALL(SCIPgetSolVals((problem->blocks[b]).subscip, blocksol, nblockvars, blockvars, blocksolvals));

         for (i = 0; i < nblockvars; i++)
         {
            SCIP_VAR *origvar;
            SCIP_Real solval;

            origvar = SCIPfindVar(scip, SCIPvarGetName(blockvars[i]));
            solval = blocksolvals[i];
#if 0
            SCIPdebugMsg(scip,"Fixing: %s = %.2f\n",SCIPvarGetName(origvar),solval);
            SCIPdebugMsg(scip, "Status fixed: %d\n", SCIPvarGetStatus(origvar));
#endif
            SCIP_CALL(SCIPsetSolVal(scip, newsol, origvar, solval));
         }
      }

      /* treat objective constant with name "OBJECTIVE_CONSTANT" */
      origobjconstvar = SCIPfindVar(scip, "OBJECTIVE_CONSTANT");
      if (origobjconstvar != NULL)
      {
         SCIPdebugMsg(scip, "Fix variable OBJECTIVE_CONSTANT\n");
         SCIP_CALL(SCIPsetSolVal(scip, newsol, origobjconstvar, 1.0));
      }

      SCIP_CALL(SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success));
      if (!success)
         SCIPdebugMsg(scip, "Solution copy failed\n");
      else
      {
         SCIP_SOL *sol;
         SCIP_Bool feasible;
         sol = SCIPgetBestSol(scip);
         SCIP_CALL(SCIPcheckSolOrig(scip, sol, &feasible, TRUE, TRUE));

         if (feasible)
         {
            SCIPdebugMsg(scip, "Solution copy successful\n");
            SCIPdebugMsg(scip, "Objective value %.2f\n", SCIPgetSolOrigObj(scip, sol));
            *result = SCIP_FOUNDSOL;
         }
      }

      SCIPfreeBufferArray(scip, &blocksolvals);
   }

TERMINATE:
   /* release variables, constraints and free memory */
   if (problem != NULL)
   {
      for (b = 0; b < problem->nblocks; b++)
      {
         for (i = 0; i < blocktolinkvars[b].size; i++)
         {
            SCIP_CALL(SCIPreleaseVar((problem->blocks[b]).subscip, &((problem->blocks[b]).slackspos[i])));
            SCIP_CALL(SCIPreleaseVar((problem->blocks[b]).subscip, &((problem->blocks[b]).slacksneg[i])));
            SCIP_CALL(SCIPreleaseCons((problem->blocks[b]).subscip, &((problem->blocks[b]).couplingcons[i])));
         }
      }

      for (b = 0; b < problem->nblocks; b++)
         SCIPfreeBufferArray(scip, &(blocktolinkvars[b].indexes));
   }

   if (htable != NULL)
      SCIPhashtableFree(&htable);

   if (blockinfolist != NULL)
      SCIPfreeBufferArray(scip, &blockinfolist);

   if (tmpcouplingcoef != NULL)
      SCIPfreeBufferArray(scip, &tmpcouplingcoef);

   if (tmpcouplingvars != NULL)
      SCIPfreeBufferArray(scip, &tmpcouplingvars);

   if (blocktolinkvars != NULL)
      SCIPfreeBufferArray(scip, &blocktolinkvars);

   for (i = 0; i < numlinkvars; i++)
      if (linkvartoblocks[i].indexes != NULL)
         SCIPfreeBufferArray(scip, &(linkvartoblocks[i].indexes));

   if (linkvartoblocks != NULL)
      SCIPfreeBufferArray(scip, &linkvartoblocks);

   if (linkvars != NULL)
      SCIPfreeBufferArray(scip, &linkvars);

   if (blockstartsconss != NULL)
      SCIPfreeBufferArray(scip, &blockstartsconss);

   if (conslabels != NULL)
      SCIPfreeBufferArray(scip, &conslabels);

   if (varslabels != NULL)
      SCIPfreeBufferArray(scip, &varslabels);

   if (problem != NULL)
      freeProblem(&problem);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the PADM primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurPADM(
    SCIP *scip /**< SCIP data structure */
)
{
   SCIP_HEURDATA *heurdata;
   SCIP_HEUR *heur;

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
   SCIP_CALL(SCIPincludeHeurBasic(scip, &heur,
                                  HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
                                  HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecPADM, heurdata));

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL(SCIPsetHeurCopy(scip, heur, heurCopyPADM));
   SCIP_CALL(SCIPsetHeurFree(scip, heur, heurFreePADM));
   SCIP_CALL(SCIPsetHeurInit(scip, heur, heurInitPADM));
   SCIP_CALL(SCIPsetHeurExit(scip, heur, heurExitPADM));
   SCIP_CALL(SCIPsetHeurInitsol(scip, heur, heurInitsolPADM));
   SCIP_CALL(SCIPsetHeurExitsol(scip, heur, heurExitsolPADM));
#endif

   /* add xyz primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
