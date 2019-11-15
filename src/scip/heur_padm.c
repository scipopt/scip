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

/**@file   heur_padm.c
 * @brief  PADM primal heuristic based on ideas published in the paper
 *         "A Decomposition Heuristic for Mixed-Integer Supply Chain Problems"
 *         by Martin Schmidt, Lars Schewe, and Dieter Weninger
 * @author Dieter Weninger
 * @author Katrin Halbig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "blockmemshell/memory.h"
#include "scip/cons_linear.h"
#include "scip/dcmp.h"
#include "scip/debug.h"
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
#include "scip/scipdefplugins.h"
#include "scip/scip_branch.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_dcmp.h"
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

#define HEUR_NAME             "padm"
#define HEUR_DESC             "penalty alternating direction method primal heuristic"
#define HEUR_DISPCHAR         'D'
#define HEUR_PRIORITY         70000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE                  /**< does the heuristic use a secondary SCIP instance? */

#define COUPLINGSIZE          3

/*
 * Data structures
 */

/** data related to one problem (see below) */
typedef struct Problem PROBLEM;

/** data related to one block */
typedef struct Block
{
   PROBLEM*              problem;            /**< the problem this block belongs to */
   SCIP*                 subscip;            /**< sub-SCIP representing this block */
   int                   number;             /**< component number */
   SCIP_VAR**            subvars;            /**< variables belonging to this block (without slack variables) */
   int                   nsubvars;           /**< number of variables belonging to this block (without slack variables) */
   SCIP_VAR**            slackspos;          /**< positive slack variables */
   SCIP_VAR**            slacksneg;          /**< negative slack variables */
   SCIP_CONS**           couplingcons;       /**< coupling contraints */
   int                   ncoupling;          /**< number of coupling contraints (equal to positive/negative slack variables) */
} BLOCK;

/** data related to one problem */
struct Problem
{
   SCIP*                 scip;               /**< the SCIP instance this problem belongs to */
   char*                 name;               /**< name of the problem */
   BLOCK*                blocks;             /**< blocks into which the problem will be divided */
   int                   nblocks;            /**< number of blocks */
};

/** set data structure */
typedef struct set
{
   int                   size;               /**< size of the set */
   int*                  indexes;            /**< set of indexes */
} SET;

/** data of one linking variable related to one block */
typedef struct blockinfo
{
   int                   block;              /**< index of this block */
   int                   otherblock;         /**< index of the other conntected block */
   int                   linkVarIdx;         /**< linking variable index */
   SCIP_Real             linkVarVal;         /**< value of linking variable */
   SCIP_VAR*             linkVar;            /**< linking variable */
   SCIP_Real             slackPosObjCoeff;   /**< penalty coefficient of positive slack variable */
   SCIP_VAR*             slackPosVar;        /**< positive slack variable */
   SCIP_Real             slackNegObjCoeff;   /**< penalty coefficient of negative slack variable */
   SCIP_VAR*             slackNegVar;        /**< negative slack variable */
   SCIP_CONS*            couplingCons;       /**< coupling contraint (equation) */
} BLOCKINFO;

/** returns TRUE iff both keys are equal */
static
SCIP_DECL_HASHKEYEQ(indexesEqual)
{
   BLOCKINFO* binfo1;
   BLOCKINFO* binfo2;

   binfo1 = (BLOCKINFO*) key1;
   binfo2 = (BLOCKINFO*) key2;

   if( binfo1->block != binfo2->block )
      return FALSE;

   if( binfo1->otherblock != binfo2->otherblock )
      return FALSE;

   if( binfo1->linkVarIdx != binfo2->linkVarIdx )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(indexesHashval)
{ /*lint --e{715}*/
   BLOCKINFO* binfo;
   binfo = (BLOCKINFO*) key;

   return SCIPhashFour(SCIPrealHashCode(binfo->block), SCIPrealHashCode(binfo->otherblock),
                        SCIPrealHashCode(binfo->linkVarIdx), SCIPrealHashCode(binfo->linkVarIdx));
}

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   admiterations;      /**< maximal number of ADM iterations in each penalty loop */
   int                   penaltyiterations;  /**< maximal number of penalty iterations */
   SCIP_Real             gap;                /**< mipgap at start */
   SCIP_Bool             scaling;            /**< enable sigmoid rescaling of penalty parameters */
   SCIP_Bool             assignlinking;      /**< should linking constraints be assigned? */
   SCIP_Bool             original;           /**< should the original problem be used? */
};

/*
 * Local methods
 */

/** initializes one block */
static
SCIP_RETCODE initBlock(
   PROBLEM*              problem             /**< problem structure */
   )
{
   BLOCK* block;

   assert(problem != NULL);
   assert(problem->scip != NULL);

   block = &problem->blocks[problem->nblocks];

   block->problem = problem;
   block->subscip = NULL;
   block->subvars = NULL;
   block->nsubvars = 0;
   block->number = problem->nblocks;
   block->slackspos = NULL;
   block->slacksneg = NULL;
   block->couplingcons = NULL;
   block->ncoupling = 0;

   ++problem->nblocks;

   return SCIP_OKAY;
}

/** frees component structure */
static
SCIP_RETCODE freeBlock(
   BLOCK*                block               /**< block structure */
   )
{
   assert(block != NULL);

   block->ncoupling = 0;

   if( block->subscip != NULL )
      SCIP_CALL( SCIPfree(&block->subscip) );

   return SCIP_OKAY;
}

/** initializes subproblem structure */
static
SCIP_RETCODE initProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBLEM**             problem,            /**< pointer to problem structure */
   int                   nblocks             /**< number of blocks */
   )
{
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, problem) );
   assert(*problem != NULL);

   (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*problem)->name, name, strlen(name) + 1) );

   SCIPdebugMessage("initialized problem <%s>\n", (*problem)->name);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*problem)->blocks, nblocks) );

   (*problem)->scip = scip;
   (*problem)->nblocks = 0;

   return SCIP_OKAY;
}

/** frees subproblem structure */
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

   /* free all blocks */
   for( c = (*problem)->nblocks - 1; c >= 0; --c )
   {
      SCIP_CALL(freeBlock(&(*problem)->blocks[c]));
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

/** creates a sub-SCIP for the given variables and constraints */
static
SCIP_RETCODE createSubscip(
   SCIP*                 scip,               /**< main SCIP data structure */
   SCIP**                subscip             /**< pointer to store created sub-SCIP */
   )
{
   /* create a new SCIP instance */
   SCIP_CALL( SCIPcreate(subscip) );

   SCIP_CALL( SCIPincludeDefaultPlugins(*subscip) );

   SCIP_CALL( SCIPcopyLimits(scip, *subscip) );

   /* avoid recursive calls */
   SCIP_CALL( SCIPsetSubscipsOff(*subscip, TRUE) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(*subscip, "misc/catchctrlc", FALSE) );

#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/freq", 100000000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(*subscip, "timing/statistictiming", FALSE) );
#endif

   return SCIP_OKAY;
}

/** copies the given variables and constraints to the given sub-SCIP */
static
SCIP_RETCODE copyToSubscip(
   SCIP*                 scip,               /**< source SCIP */
   SCIP*                 subscip,            /**< target SCIP */
   const char*           name,               /**< name for copied problem */
   SCIP_CONS**           conss,              /**< constraint to copy */
   SCIP_HASHMAP*         varmap,             /**< hashmap used for the copy process of variables */
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
   SCIP_CALL( SCIPcopyProb(scip, subscip, varmap, consmap, FALSE, name) );

   /* copy constraints */
   for( i = 0; i < nconss; ++i )
   {
      assert(!SCIPconsIsModifiable(conss[i]));

      /* copy the constraint */
      SCIP_CALL( SCIPgetConsCopy(scip, subscip, conss[i], &newcons, SCIPconsGetHdlr(conss[i]), varmap, consmap, NULL,
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

/** creates the subscip for a given block */
static
SCIP_RETCODE blockCreateSubscip(
   BLOCK*                block,              /**< block structure */
   SCIP_HASHMAP*         varmap,             /**< variable hashmap used to improve performance */
   SCIP_HASHMAP*         consmap,            /**< constraint hashmap used to improve performance */
   SCIP_CONS**           conss,              /**< constraints contained in this block */
   int                   nconss,             /**< number of constraints contained in this block */
   SCIP_Bool*            success             /**< pointer to store whether the copying process was successful */
   )
{
   char name[SCIP_MAXSTRLEN];
   PROBLEM* problem;
   SCIP* scip;
   SCIP_VAR** subscipvars;
   int nsubscipvars;
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

   SCIP_CALL( createSubscip(scip, &block->subscip) );

   if( block->subscip != NULL )
   {
      /* get name of the original problem and add "comp_nr" */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp_%d", problem->name, block->number);

      SCIP_CALL( copyToSubscip(scip, block->subscip, name, conss, varmap, consmap, nconss, success) );

      /* save variables of subscip (without slack variables) */
      nsubscipvars = SCIPgetNOrigVars(block->subscip);
      subscipvars = SCIPgetOrigVars(block->subscip);
      SCIP_CALL( SCIPallocBufferArray(scip, &(block->subvars), nsubscipvars) );
      block->nsubvars = nsubscipvars;
      for( i = 0; i < nsubscipvars; i++ )
         block->subvars[i] = subscipvars[i];

      if( !(*success) )
      {
         SCIP_CALL( SCIPfree(&block->subscip) );
         block->subscip = NULL;
      }
   }
   else
      (*success) = FALSE;

   SCIPdebugMsg(scip, "created subscip of block %d\n", block->number);

   return SCIP_OKAY;
}

/** creates problem structure and split it into blocks */
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
   SCIP_HASHMAP* varmap;
   SCIP_HASHMAP* consmap;
   SCIP_CONS** blockconss;
   int nblockconss;
   int b;
   SCIP_Bool success = TRUE;

   /* init subproblem data structure */
   SCIP_CALL( initProblem(scip, problem, nblocks) );
   assert((*problem)->blocks != NULL);

   /* hashmap mapping from original constraints to constraints in the sub-SCIPs (for performance reasons) */
   SCIP_CALL( SCIPhashmapCreate(&consmap, SCIPblkmem(scip), blockstartsconss[nblocks]) );

   /* loop over all blocks */
   for( b = 0; b < nblocks; b++ )
   {
      SCIP_CALL( initBlock(*problem) );
      assert((*problem)->nblocks == b + 1);

      block = &(*problem)->blocks[b];

      SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );

      /* get block constraints */
      blockconss = &(sortedconss[blockstartsconss[b]]);
      nblockconss = blockstartsconss[b + 1] - blockstartsconss[b];

      /* build subscip for block */
      SCIP_CALL( blockCreateSubscip(block, varmap, consmap, blockconss, nblockconss, &success) );

      SCIPhashmapFree(&varmap);

      if( !success )
         break;
   }
   assert(nblocks == (*problem)->nblocks);

   SCIPhashmapFree(&consmap);

   if( !success )
   {
      /* free subproblem data structure since not all blocks could be copied */
      SCIP_CALL( freeProblem(problem) );
   }

   return SCIP_OKAY;
}

/** copies labels to newdecomp and assigns linking constraints if possible*/
static
SCIP_RETCODE assignLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_DECOMP*    decomp,             /**< decomposition */
   SCIP_DECOMP*          newdecomp,          /**< decomposition with (partially) assigned linking constraints */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_CONS**           sortedconss,        /**< sorted array of constraints */
   int*                  varlabels,          /**< array of variable labels */
   int*                  conslabels,         /**< sorted array of constraint labels */
   int                   nvars,              /**< number of variables */
   int                   nconss              /**< number of constraints */
   )
{
   int nlinkconss;
   int c;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(vars != NULL);
   assert(sortedconss != NULL);
   assert(varlabels != NULL);
   assert(conslabels != NULL);

   /* copy the labels */
   SCIP_CALL( SCIPdecompSetVarsLabels(newdecomp, vars, varlabels, nvars) );
   SCIP_CALL( SCIPdecompSetConsLabels(newdecomp, sortedconss, conslabels, nconss) );

   nlinkconss = 0;
   for( c = 0; c < nconss; c++ )
   {
      if( conslabels[c] == SCIP_DECOMP_LINKCONS )
         nlinkconss++;
      else
         break;
   }

   SCIPdebugMsg(scip, "try to assign %d linking constraints\n", nlinkconss);

   /* reassign linking constraints */
   SCIP_CALL( SCIPassignDecompLinkConss(scip, newdecomp, &sortedconss[0], nlinkconss, NULL) );

   SCIP_CALL( SCIPcomputeDecompVarsLabels(scip, newdecomp, sortedconss, nconss) );

   SCIP_CALL( SCIPcomputeDecompStats(scip, newdecomp, TRUE) );

   SCIPdecompGetConsLabels(newdecomp, sortedconss, conslabels, nconss);
   SCIPdecompGetVarsLabels(newdecomp, vars, varlabels, nvars);

   SCIPsortIntPtr(conslabels, (void**)sortedconss, nconss);

   return SCIP_OKAY;
}

/** computes feasible solution from last stored solution of the block*/
static
SCIP_RETCODE reuseSolution(
   SCIP*                 subscip,            /**< SCIP data structure */
   BLOCK*                block               /**< block structure*/
   )
{
   SCIP_SOL** sols;
   SCIP_SOL* sol; /* solution of block that will be repaired */
   SCIP_SOL* newsol;
   SCIP_VAR** blockvars;
   SCIP_VAR** consvars;
   double* blockvals;
   int nsols;
   int nvars;
   int c;
   SCIP_Bool success;

   assert(subscip != NULL);
   assert(block != NULL);

   nsols = SCIPgetNSols( subscip );

   /* no solution in solution candidate storage found */
   if( nsols == 0 )
      return SCIP_OKAY;
   
   SCIP_CALL( SCIPallocBufferArray(subscip, &consvars, COUPLINGSIZE) );

   sols = SCIPgetSols( subscip );
   sol = sols[nsols - 1];

   /* copy the solution */
   nvars = SCIPgetNVars(subscip);
   blockvars = SCIPgetVars(subscip);
   SCIP_CALL( SCIPallocBufferArray(subscip, &blockvals, nvars) );
   SCIP_CALL( SCIPgetSolVals(subscip, sol, nvars, blockvars, blockvals) );
   SCIP_CALL( SCIPcreateOrigSol(subscip, &newsol, NULL) );
   SCIP_CALL( SCIPsetSolVals(subscip, newsol, nvars, blockvars, blockvals) );

   /* correct each coupling constraint;
      orig_var + slackpos - slackneg == side
      adapt slack variables so that constraint is feasible */
   for( c = 0; c < block->ncoupling; c++ )
   {
      double solval; /* old solution values of variables; [0] original variable, [1] slackpos, [2] slackneg */
      double side; /* current right hand side */
      double diff;

      SCIP_CALL( SCIPgetConsVars(subscip, block->couplingcons[c], consvars, COUPLINGSIZE, &success) );
      solval = SCIPgetSolVal(subscip, sol, consvars[0]);

      side = SCIPgetRhsLinear(subscip, block->couplingcons[c]);
      assert(SCIPgetRhsLinear(subscip, block->couplingcons[c]) == SCIPgetLhsLinear(subscip, block->couplingcons[c]));

      diff = side - solval;

      /* slackpos is strict positiv */
      if( diff > 0 )
      {
         SCIP_CALL( SCIPsetSolVal(subscip, newsol, block->slackspos[c], diff) );
         SCIP_CALL( SCIPsetSolVal(subscip, newsol, block->slacksneg[c], 0.0) );
      }
      /* slackneg is strict positiv */
      else if( diff < 0 )
      {
         SCIP_CALL( SCIPsetSolVal(subscip, newsol, block->slacksneg[c], -diff) );
         SCIP_CALL( SCIPsetSolVal(subscip, newsol, block->slackspos[c], 0.0) );
      }
      /* no slack variable necessary */
      else
      {
         SCIP_CALL( SCIPsetSolVal(subscip, newsol, block->slackspos[c], 0.0) );
         SCIP_CALL( SCIPsetSolVal(subscip, newsol, block->slacksneg[c], 0.0) );
      }
   }

   SCIPdebugMsg(subscip, "Try adding solution with objective value %.2f\n", SCIPgetSolOrigObj(subscip, newsol));
   SCIP_CALL( SCIPaddSolFree(subscip, &newsol, &success) );

   if( !success )
      SCIPdebugMsg(subscip, "Correcting solution failed\n"); /* maybe not better than old solutions */
   else
   {
      SCIPdebugMsg(subscip, "Correcting solution successful\n");
   }
   
   SCIPfreeBufferArray(subscip, &blockvals);
   SCIPfreeBufferArray(subscip, &consvars);

   return SCIP_OKAY;
}

/** returns the available time limit that is left */
static
SCIP_Real getTimeLeft(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real timelim;

   assert(scip != NULL);

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelim) );

   if( !SCIPisInfinity(scip, timelim) )
      return MAX(0.0, (timelim - SCIPgetSolvingTime(scip)));
   else
      return SCIPinfinity(scip);
}

/** returns the available memory limit that is left */
static
SCIP_Real getMemLeft(
   SCIP*                 scip,               /**< main SCIP data structure */
   PROBLEM*              problem,            /**< contains all subscips */
   SCIP*                 subscip             /**< exclude this subscip (or NULL) */
   )
{
   SCIP_Real memlim;
   int b;

   assert(scip != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memlim) );

   if( !SCIPisInfinity(scip, memlim) )
   {
      /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
      memlim -= SCIPgetMemUsed(scip)/1048576.0;
      memlim -= SCIPgetMemExternEstim(scip)/1048576.0;

      /* substract the memory already used by all other subSCIPs and the estimated memory usage of external software */
      for( b = 0; b < problem->nblocks; b++ )
      {
         if( (problem->blocks[b]).subscip == subscip )
            continue;

         memlim -= SCIPgetMemUsed((problem->blocks[b]).subscip)/1048576.0;
         memlim -= SCIPgetMemExternEstim((problem->blocks[b]).subscip)/1048576.0;
      }

      return MAX(0.0, memlim);
   }
   else
      return SCIPinfinity(scip);
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyPADM)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurPADM(scip) );

   return SCIP_OKAY;
}


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreePADM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   SCIPheurSetData(heur, NULL);

   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static SCIP_DECL_HEUREXEC(heurExecPADM)
{ /*lint --e{715}*/
   char name[SCIP_MAXSTRLEN];
   char info[SCIP_MAXSTRLEN];
   SCIP_HEURDATA* heurdata;
   PROBLEM* problem;
   SCIP_DECOMP** alldecomps;
   SCIP_DECOMP* decomp;
   SCIP_VAR** vars;
   SCIP_VAR** linkvars;
   SCIP_VAR** tmpcouplingvars;
   SCIP_CONS** conss;
   SCIP_CONS** sortedconss;
   SET* linkvartoblocks;
   SET* blocktolinkvars;
   BLOCKINFO* blockinfolist;
   SCIP_HASHTABLE* htable;
   int* varlabels;
   int* conslabels;
   int* blockstartsconss;
   int* allLinkVarToBlocks; /* for efficient memory allocation */
   SCIP_Bool* varonlyobj;
   SCIP_Real* tmpcouplingcoef;
   SCIP_Real gap;
   SCIP_Real maxslack;
   SCIP_Real slackthreshold;
   SCIP_Real memory; /* in MB */
   SCIP_STATUS status;
   SCIP_Bool solutionsdiffer;
   SCIP_Bool solved;
   SCIP_Bool doScaling;
   int ndecomps;
   int nconss;
   int nvars;
   int nblocks;
   int numlinkvars;
   int nentries;
   int aIter;
   int pIter;
   int increasedslacks;
   int blockinfolistfill;
   int i;
   int b;
   int k;
   int j;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   problem = NULL;
   sortedconss = NULL;
   varlabels = NULL;
   conslabels = NULL;
   blockstartsconss = NULL;
   allLinkVarToBlocks = NULL;
   linkvars = NULL;
   linkvartoblocks = NULL;
   numlinkvars = 0;
   blocktolinkvars = NULL;
   tmpcouplingvars = NULL;
   tmpcouplingcoef = NULL;
   varonlyobj = NULL;
   blockinfolist = NULL;
   htable = NULL;

   gap = heurdata->gap;

   if( heurtiming & SCIP_HEURTIMING_BEFORENODE )
   {
      SCIPdebugMsg(scip, "Initialize padm heuristic before node\n");
   }
   if( heurtiming & SCIP_HEURTIMING_AFTERNODE )
   {
      SCIPdebugMsg(scip, "Initialize padm heuristic after node\n");
   }

#if 0
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "orig_%s", SCIPgetProbName(scip));
   SCIP_CALL( SCIPwriteOrigProblem(scip, name, "mps", FALSE) );
#endif
#if 0
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "trans_%s", SCIPgetProbName(scip));
   SCIP_CALL( SCIPwriteTransProblem(scip, name, "mps", FALSE) );
#endif

   /* take original problem */
   if( heurdata->original )
   {
      SCIPgetDecomps(scip, &alldecomps, &ndecomps, TRUE);
      if( ndecomps == 0)
         return SCIP_OKAY;

      /* it takes the first decomposition */
      decomp = alldecomps[0];
      SCIPdebugMsg(scip, "First original decomposition is selected\n");
      assert(decomp != NULL);

      nconss = SCIPgetNOrigConss(scip);
      conss = SCIPgetOrigConss(scip);
      nvars = SCIPgetNOrigVars(scip);
      vars = SCIPgetOrigVars(scip);
   }
   /* take transformed problem */
   else
   {
      SCIPgetDecomps(scip, &alldecomps, &ndecomps, FALSE);
      if( ndecomps == 0)
         return SCIP_OKAY;

      /* it takes the first decomposition */
      decomp = alldecomps[0];
      SCIPdebugMsg(scip, "First transformed decomposition is selected\n");
      assert(decomp != NULL);

      nconss = SCIPgetNConss(scip);
      conss = SCIPgetConss(scip);
      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);
   }

   nblocks = SCIPdecompGetNBlocks(decomp);

   /* if problem has no constraints, no variables or less than two blocks, return */
   if( nconss == 0 || nvars == 0 || nblocks <= 1 )
   {
      SCIPdebugMsg(scip, "problem has no constraints, no variables or less than two blocks\n");
      goto TERMINATE;
   }

   /* don't change problem by sorting constraints */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sortedconss, conss, nconss) );

   SCIP_CALL( SCIPallocBufferArray(scip, &varlabels, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blockstartsconss, nconss + 1) );

   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);
   SCIPdecompGetVarsLabels(decomp, vars, varlabels, nvars);

   /* sort constraints by blocks */
   SCIPsortIntPtr(conslabels, (void**)sortedconss, nconss);

   /* try to assign linking constraints */
   if( heurdata->assignlinking && conslabels[0] == SCIP_DECOMP_LINKCONS )
   {
      /* create new decomposition; don't change the decompositions in the decompstore */
      SCIP_DECOMP* newdecomp = NULL;
      SCIP_CALL( SCIPdecompCreate(&newdecomp, SCIPblkmem(scip), nblocks, heurdata->original, SCIPdecompUseBendersLabels(decomp)) );

      SCIP_CALL( assignLinking(scip, decomp, newdecomp, vars, sortedconss, varlabels, conslabels, nvars, nconss) );
      decomp = newdecomp;

      /* number of blocks can have changed */
      nblocks = SCIPdecompGetNBlocks(decomp);

      /* new decomp can already be deleted here because we don't need it anymore */
      SCIPdecompFree(&newdecomp, SCIPblkmem(scip));
   }

   if( conslabels[0] == SCIP_DECOMP_LINKCONS )
   {
      SCIPdebugMsg(scip, "No support for linking contraints\n");
      goto TERMINATE;
   }

   /* determine start indices of blocks in sorted conss array */
   i = 0;
   for( b = 0; b < nblocks + 1; ++b )
   {
      blockstartsconss[b] = i;
      if( i == nconss )
         break;

      do
      {
         ++i;
      }
      while( i < nconss && conslabels[i] == conslabels[i-1] );
   }

   /* estimate required memory for all blocks and terminate if not enough memory is available */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memory) );
   if( ((SCIPgetMemUsed(scip) + SCIPgetMemExternEstim(scip))/1048576.0) * nblocks >= memory )
   {
      SCIPdebugMsg(scip, "The estimated memory usage is too large.\n");
      goto TERMINATE;
   }

   /* create blockproblems */
   SCIP_CALL( createAndSplitProblem(scip, sortedconss, blockstartsconss, nblocks, &problem) );

   /* count linking variables */
   for( i = 0; i < nvars; i++ )
   {
      if( varlabels[i] == SCIP_DECOMP_LINKVAR )
         numlinkvars++;
   }
   SCIPdebugMsg(scip, "%d linking variables\n", numlinkvars);

   SCIP_CALL( SCIPallocBufferArray(scip, &linkvartoblocks, numlinkvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocktolinkvars, problem->nblocks) );

   /* set pointer to NULL for safe memory release */
   for( i = 0; i < numlinkvars; i++ )
      linkvartoblocks[i].indexes = NULL;
   for( i = 0; i < problem->nblocks; i++ )
      blocktolinkvars[i].indexes = NULL;

   /* extract linking variables and init linking variable to blocks set */
   SCIP_CALL( SCIPallocBufferArray(scip, &allLinkVarToBlocks, problem->nblocks * numlinkvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linkvars, numlinkvars) );

   b = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( varlabels[i] == SCIP_DECOMP_LINKVAR )
      {
         linkvars[b] = vars[i];
         linkvartoblocks[b].indexes = &allLinkVarToBlocks[b * problem->nblocks];
         linkvartoblocks[b].size = 0;
         b++;
      }
   }

   /* fill linking variable to blocks set */
   for( i = 0; i < numlinkvars; i++ )
   {
      SCIP_VAR *var;
      const char *vname;

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

   /* check whether there is enough time and memory left */
   if( getTimeLeft(scip) <= 0 || getMemLeft(scip, problem, NULL) <= 10)
   {
      SCIPdebugMsg(scip, "no time or memory left\n");
      goto TERMINATE;
   }

   /* init block to linking variables set */
   for( b = 0; b < problem->nblocks; b++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(blocktolinkvars[b].indexes), numlinkvars) );
      blocktolinkvars[b].size = 0;
   }

   /* init varonlyobj; true if variable is only part of the objective function */
   SCIP_CALL( SCIPallocBufferArray(scip, &varonlyobj, numlinkvars) );
   for( i = 0; i < numlinkvars; ++i)
      varonlyobj[i] = TRUE;

   /* fill block to linking variables set */
   for( b = 0; b < problem->nblocks; b++ )
   {
      k = 0;
      for( i = 0; i < numlinkvars; i++ )
      {
         SCIP_VAR *var;
         const char *vname;

         vname = SCIPvarGetName(linkvars[i]);
         var = SCIPfindVar((problem->blocks[b]).subscip, vname);
         if( var != NULL )
         {
            varonlyobj[i] = FALSE;
            blocktolinkvars[b].indexes[k] = i;
            blocktolinkvars[b].size = k + 1;
            k++;
         }
      }
   }

#if 0
   for( i = 0; i < numlinkvars; i++ )
   {
      if(varonlyobj[i] == TRUE)
         SCIPdebugMsg(scip, "variable %s has no constraints\n", SCIPvarGetName(linkvars[i]));
   }

   for( b = 0; b < problem->nblocks; b++ )
   {
      SCIPdebugMsg(scip, "%d linking variables in block %d\n", blocktolinkvars[b].size, b);
   }
#endif

   /* init arrays for slack variables and coupling constraints */
   for( b = 0; b < problem->nblocks; b++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &((problem->blocks[b]).slackspos), blocktolinkvars[b].size * (nblocks - 1)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &((problem->blocks[b]).slacksneg), blocktolinkvars[b].size * (nblocks - 1)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &((problem->blocks[b]).couplingcons), blocktolinkvars[b].size * (nblocks - 1)) );
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

   SCIP_CALL( SCIPallocBufferArray(scip, &blockinfolist, nentries) );
   SCIP_CALL( SCIPhashtableCreate(&htable, SCIPblkmem(scip), 1, SCIPhashGetKeyStandard, indexesEqual, indexesHashval, (void *)scip) );
   blockinfolistfill = 0;

   /* extend submips */
   SCIPdebugMsg(scip, "Extending %d block models\n", problem->nblocks);
   for( b = 0; b < problem->nblocks; b++ )
   {
      SCIP_VAR **blockvars;
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
            int b2;
            b2 = linkvartoblocks[linkvaridx].indexes[k];

            /* handle different blocks with common linking variable */
            if( b2 != b )
            {
               BLOCKINFO *binfo;
               binfo = &blockinfolist[blockinfolistfill];
               blockinfolistfill++;
               binfo->block = b;
               binfo->otherblock = b2;
               binfo->linkVarIdx = linkvaridx;
               binfo->linkVarVal = 0.0;
               binfo->linkVar = SCIPfindVar((problem->blocks[b]).subscip, SCIPvarGetName(linkvars[linkvaridx]));
               j = (problem->blocks[b]).ncoupling;

               /* create positive slack variable */
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_slackpos_block_%d", SCIPvarGetName(linkvars[linkvaridx]), b2);
               (problem->blocks[b]).slackspos[j] = NULL;
               SCIP_CALL( SCIPcreateVarBasic((problem->blocks[b]).subscip,
                                            &((problem->blocks[b]).slackspos[j]), name,
                                            0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
               SCIP_CALL( SCIPaddVar((problem->blocks[b]).subscip, (problem->blocks[b]).slackspos[j]) );
               assert((problem->blocks[b]).slackspos[j] != NULL);
               binfo->slackPosObjCoeff = 1.0;
               binfo->slackPosVar = (problem->blocks[b]).slackspos[j];

               /* create negative slack variable */
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_slackneg_block_%d", SCIPvarGetName(linkvars[linkvaridx]), b2);
               (problem->blocks[b]).slacksneg[j] = NULL;
               SCIP_CALL( SCIPcreateVarBasic((problem->blocks[b]).subscip,
                                            &((problem->blocks[b]).slacksneg[j]), name,
                                            0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
               SCIP_CALL( SCIPaddVar((problem->blocks[b]).subscip, (problem->blocks[b]).slacksneg[j]) );
               assert((problem->blocks[b]).slacksneg[j] != NULL);
               binfo->slackNegObjCoeff = 1.0;
               binfo->slackNegVar = (problem->blocks[b]).slacksneg[j];

               /* fill variables for linking constraint */
               tmpcouplingvars[0] = binfo->linkVar;
               tmpcouplingvars[1] = binfo->slackPosVar;
               tmpcouplingvars[2] = binfo->slackNegVar;

               /* create linking constraint */
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_coupling_block_%d",
                            SCIPvarGetName(linkvars[linkvaridx]), b2);
               (problem->blocks[b]).couplingcons[j] = NULL;

               /* create linking constraint with initial side equal to zero */
               if( heurtiming & SCIP_HEURTIMING_BEFORENODE )
               {
                  SCIP_CALL( SCIPcreateConsBasicLinear((problem->blocks[b]).subscip, &((problem->blocks[b]).couplingcons[j]),
                                                   name, COUPLINGSIZE, tmpcouplingvars, tmpcouplingcoef, 0.0, 0.0) );
               }

               /* create linking constraint with initial side equal to LP solution */
               if( heurtiming & SCIP_HEURTIMING_AFTERNODE )
               {
                  SCIP_CALL( SCIPcreateConsBasicLinear((problem->blocks[b]).subscip, &((problem->blocks[b]).couplingcons[j]),
                                                   name, COUPLINGSIZE, tmpcouplingvars, tmpcouplingcoef,
                                                   SCIPvarGetLPSol(linkvars[linkvaridx]), SCIPvarGetLPSol(linkvars[linkvaridx])) );
               }

               SCIP_CALL( SCIPaddCons((problem->blocks[b]).subscip, (problem->blocks[b]).couplingcons[j]) );
               assert((problem->blocks[b]).couplingcons[j] != NULL);
               binfo->couplingCons = (problem->blocks[b]).couplingcons[j];

               (problem->blocks[b]).ncoupling++;

               /* feed hashtable */
               SCIP_CALL( SCIPhashtableSafeInsert(htable, (void*) binfo) );
            }
         }
      }
   }
   assert(nentries == SCIPhashtableGetNElements(htable));

#if 0
   /* write extended submips */
   for( b = 0; b < problem->nblocks; b++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_block_%d.lp", SCIPgetProbName(scip), b);
      SCIP_CALL( SCIPwriteOrigProblem((problem->blocks[b]).subscip, name, "lp", FALSE) );
   }
#endif

   /* determine threshold for penalty coefficients via maximum norm */
   slackthreshold = SCIP_REAL_MIN;
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real obj;

      obj = REALABS(SCIPvarGetObj(vars[i]));
      if( obj > slackthreshold )
         slackthreshold = obj;
   }

   /* ------------------------------------------------------------------------------------------------- */

   /* check whether there is enough time and memory left */
   if( getTimeLeft(scip) <= 0 || getMemLeft(scip, problem, NULL) <= 10)
   {
      SCIPdebugMsg(scip, "no time or memory left\n");
      goto TERMINATE;
   }

   SCIPdebugMsg(scip, "Starting iterations\n");
   SCIPdebugMsg(scip, "PIt\tADMIt\tSlacks\tInfo\n");

   pIter = 0;
   increasedslacks = 0;
   (void) SCIPsnprintf(info, SCIP_MAXSTRLEN, "-");
   solved = FALSE;

   /* Penalty loop */
   while( !solved && pIter < heurdata->penaltyiterations )
   {
      pIter++;
      solutionsdiffer = TRUE;
      aIter = 0;

      /*  Alternating direction method loop */
      while( solutionsdiffer && aIter < heurdata->admiterations )
      {
         aIter++;
         solutionsdiffer = FALSE;
         SCIPdebugMsg(scip, "%d\t%d\t%d\t%s\n", pIter, aIter, increasedslacks, info);

         for( b = 0; b < problem->nblocks; b++ )
         {
            for( i = 0; i < blocktolinkvars[b].size; i++ )
            {
               int linkvaridx;
               linkvaridx = blocktolinkvars[b].indexes[i];

               for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
               {
                  int b2;
                  b2 = linkvartoblocks[linkvaridx].indexes[k];

                  if( b2 != b )
                  {
                     BLOCKINFO binfo;
                     BLOCKINFO* binfoout;
                     BLOCKINFO binfo2;
                     BLOCKINFO* binfo2out;

                     SCIP_CONS* couplingcons;
                     SCIP_Real newRhs;

                     binfo.block = b;
                     binfo.otherblock = b2;
                     binfo.linkVarIdx = linkvaridx;

                     binfoout = (BLOCKINFO*) SCIPhashtableRetrieve(htable, (void *)&binfo);
                     couplingcons = binfoout->couplingCons;

                     /* interchange blocks b and b2 for getting new right hand side */
                     binfo2.block = b2;
                     binfo2.otherblock = b;
                     binfo2.linkVarIdx = linkvaridx;
                     binfo2out = (BLOCKINFO*) SCIPhashtableRetrieve(htable, (void*) &binfo2);
                     newRhs = binfo2out->linkVarVal;

                     /* change side of coupling constraint equation with linking variable value of the other block */
                     SCIP_CALL( SCIPchgLhsLinear((problem->blocks[b]).subscip, couplingcons, newRhs) );
                     SCIP_CALL( SCIPchgRhsLinear((problem->blocks[b]).subscip, couplingcons, newRhs) );

                     /* change penalty coefficients of slack variables */
                     SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, binfoout->slackPosVar, binfoout->slackPosObjCoeff) );
                     SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, binfoout->slackNegVar, binfoout->slackNegObjCoeff) );
                  }
               }
            }

            /* increase slack penalty coeffs until each subproblem can be solved to optimality */
            do
            {
               SCIP_CALL( SCIPsetRealParam((problem->blocks[b]).subscip, "limits/gap", gap) );
#if 0
               SCIPdebugMsg(scip, "write subscip block %d\n", b);
               SCIP_CALL( SCIPwriteOrigProblem((problem->blocks[b]).subscip, "debug_block_orig.lp", "lp", FALSE) );
#endif
               /* reuse old solution if available*/
               SCIP_CALL( reuseSolution((problem->blocks[b]).subscip, &problem->blocks[b]) );

               /* update time limit of subproblem */
               SCIP_CALL( SCIPsetRealParam((problem->blocks[b]).subscip, "limits/time", getTimeLeft(scip)) );

               /* update memory limit of subproblem */
               SCIP_CALL( SCIPsetRealParam((problem->blocks[b]).subscip, "limits/memory", 
                                             getMemLeft(scip, problem, (problem->blocks[b]).subscip)) );

               /* solve block */

               SCIP_CALL( SCIPsolve((problem->blocks[b]).subscip) );
               status = SCIPgetStatus((problem->blocks[b]).subscip);

               /* check solution if one of the three cases occurs
                * - solution is optimal
                * - solution reached gaplimit
                * - time limit is reached but best solution needs no slack variables (no dual solution avaiable)
                */
               if( status == SCIP_STATUS_OPTIMAL || status == SCIP_STATUS_GAPLIMIT || 
                   (status == SCIP_STATUS_TIMELIMIT && SCIPgetNSols((problem->blocks[b]).subscip) > 0 && 
                    SCIPisEQ(scip, SCIPsolGetOrigObj(SCIPgetBestSol((problem->blocks[b]).subscip)), 0.0) ) )
               {
                  SCIPdebugMsg(scip, "Block is optimal or reached gaplimit or reached time limit with at least one feasible solution.\n");
                  for( i = 0; i < blocktolinkvars[b].size; i++ )
                  {
                     int linkvaridx;
                     linkvaridx = blocktolinkvars[b].indexes[i];

                     for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
                     {
                        int b2;
                        b2 = linkvartoblocks[linkvaridx].indexes[k];

                        if( b2 != b )
                        {
                           SCIP_SOL* sol;
                           BLOCKINFO binfo;
                           BLOCKINFO* binfoout;
                           SCIP_VAR* var;
                           SCIP_Real val;

                           binfo.block = b;
                           binfo.otherblock = b2;
                           binfo.linkVarIdx = linkvaridx;
                           binfoout = (BLOCKINFO *)SCIPhashtableRetrieve(htable, (void *)&binfo);

                           sol = SCIPgetBestSol((problem->blocks[b]).subscip);
                           var = binfoout->linkVar;
                           val = SCIPgetSolVal((problem->blocks[b]).subscip, sol, var);

                           if( !EPSEQ(binfoout->linkVarVal, val, SCIP_DEFAULT_EPSILON) )
                              solutionsdiffer = TRUE;

                           binfoout->linkVarVal = val;
                        }
                     }
                  }
               }
               else if( status == SCIP_STATUS_UNBOUNDED )
               {
                  SCIPdebugMsg(scip, "Block is unbounded.\n");
                  for( i = 0; i < blocktolinkvars[b].size; i++ )
                  {
                     int linkvaridx;
                     linkvaridx = blocktolinkvars[b].indexes[i];

                     for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
                     {
                        int b2;
                        b2 = linkvartoblocks[linkvaridx].indexes[k];

                        if( b2 != b )
                        {
                           BLOCKINFO binfo;
                           BLOCKINFO* binfoout;

                           binfo.block = b;
                           binfo.otherblock = b2;
                           binfo.linkVarIdx = linkvaridx;
                           binfoout = (BLOCKINFO*) SCIPhashtableRetrieve(htable, (void*) &binfo);

                           /* increase penalty coefficients to obtain a bounded subproblem */
                           binfoout->slackPosObjCoeff *= 10.0;
                           binfoout->slackNegObjCoeff *= 10.0;
                           SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, binfoout->slackPosVar, binfoout->slackPosObjCoeff) );
                           SCIP_CALL( SCIPchgVarObj((problem->blocks[b]).subscip, binfoout->slackNegVar, binfoout->slackNegObjCoeff) );
                        }
                     }
                  }
               }
               else if( status == SCIP_STATUS_TIMELIMIT )
               {
                  SCIPdebugMsg(scip, "Block reached time limit. No optimal solution avaiable.\n");
                  goto TERMINATE;
               }
               else
               {
                  SCIPdebugMsg(scip, "Block solving status %d not supported\n", status);
                  goto TERMINATE;
               }

               /* free solving data in order to change problem */
               SCIP_CALL( SCIPfreeTransform((problem->blocks[b]).subscip) );

            }
            while( status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_GAPLIMIT &&
                   !(status == SCIP_STATUS_TIMELIMIT && SCIPgetNSols((problem->blocks[b]).subscip) > 0 && 
                    SCIPisEQ(scip, SCIPsolGetOrigObj(SCIPgetBestSol((problem->blocks[b]).subscip)), 0.0) ));
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
               int b2;
               b2 = linkvartoblocks[linkvaridx].indexes[k];

               if( b2 != b )
               {
                  SCIP_SOL* sol;
                  BLOCKINFO binfo;
                  BLOCKINFO* binfoout;
                  SCIP_Real slackPosVal;
                  SCIP_Real slackNegVal;

                  binfo.block = b;
                  binfo.otherblock = b2;
                  binfo.linkVarIdx = linkvaridx;
                  binfoout = (BLOCKINFO*) SCIPhashtableRetrieve(htable, (void*) &binfo);

                  sol = SCIPgetBestSol((problem->blocks[b]).subscip);
                  slackPosVal = SCIPgetSolVal((problem->blocks[b]).subscip, sol, binfoout->slackPosVar);
                  slackNegVal = SCIPgetSolVal((problem->blocks[b]).subscip, sol, binfoout->slackNegVar);

                  /* increase penalty coefficient of positive slack variable */
                  if( SCIPisGT(scip, slackPosVal, 0.0) )
                  {
                     binfoout->slackPosObjCoeff *= 10.0;

                     if( binfoout->slackPosObjCoeff > slackthreshold )
                        doScaling = TRUE;

                     if( binfoout->slackPosObjCoeff > maxslack )
                        maxslack = binfoout->slackPosObjCoeff;

                     solved = FALSE;
                     increasedslacks++;
                  }

                  /* increase penalty coefficient of negative slack variable */
                  if( SCIPisGT(scip, slackNegVal, 0.0) )
                  {
                     binfoout->slackNegObjCoeff *= 10.0;

                     if( binfoout->slackNegObjCoeff > slackthreshold )
                        doScaling = TRUE;

                     if( binfoout->slackNegObjCoeff > maxslack )
                        maxslack = binfoout->slackNegObjCoeff;

                     solved = FALSE;
                     increasedslacks++;
                  }
               }
            }
         }
      }

      /* should sigmoid scaling be applied? */
      if( doScaling && heurdata->scaling )
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

         for( b = 0; b < problem->nblocks; b++ )
         {
            for( i = 0; i < blocktolinkvars[b].size; i++ )
            {
               int linkvaridx;
               linkvaridx = blocktolinkvars[b].indexes[i];

               for( k = 0; k < linkvartoblocks[linkvaridx].size; k++ )
               {
                  int b2;
                  b2 = linkvartoblocks[linkvaridx].indexes[k];

                  if( b2 != b )
                  {
                     BLOCKINFO binfo;
                     BLOCKINFO* binfoout;
                     SCIP_Real oldcoeff;

                     binfo.block = b;
                     binfo.otherblock = b2;
                     binfo.linkVarIdx = linkvaridx;
                     binfoout = (BLOCKINFO*) SCIPhashtableRetrieve(htable, (void*) &binfo);

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
      if( (aIter == 1 && solutionsdiffer == FALSE) || (doScaling && heurdata->scaling) )
      {
         SCIP_Real mingap = 0.001; //todo
         SCIP_Real newgap = MAX(gap * 0.5, mingap);

         if( newgap >= mingap )
         {
            if( doScaling && heurdata->scaling )
               (void) SCIPsnprintf(info, SCIP_MAXSTRLEN, "scale, %f", newgap);
            else
               (void) SCIPsnprintf(info, SCIP_MAXSTRLEN, "%f", newgap);

            gap = newgap;
         }
      }

      /* free solution process data */
      if( !solved )
         for( b = 0; b < problem->nblocks; b++ )
            SCIP_CALL( SCIPfreeTransform((problem->blocks[b]).subscip) );
   }

   /* copy solution if present */
   if( solved )
   {
      SCIP_SOL* newsol;
      SCIP_Bool success;
      SCIP_Real* blocksolvals;

      assert(increasedslacks == 0);

      SCIP_CALL( SCIPallocBufferArray(scip, &blocksolvals, nvars) );
      SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );

      for( b = 0; b < problem->nblocks; b++ )
      {
         SCIP_SOL* blocksol;
         SCIP_VAR** blockvars;
         int nblockvars;

         /* get solution of block variables (without slack variables) */
         blocksol = SCIPgetBestSol((problem->blocks[b]).subscip);
         assert(blocksol != NULL);
         blockvars = (problem->blocks[b]).subvars;
         nblockvars = (problem->blocks[b]).nsubvars;
         SCIP_CALL( SCIPgetSolVals((problem->blocks[b]).subscip, blocksol, nblockvars, blockvars, blocksolvals) );

         for( i = 0; i < nblockvars; i++ )
         {
            SCIP_VAR *origvar;
            SCIP_Real solval;

            origvar = SCIPfindVar(scip, SCIPvarGetName(blockvars[i]));
            solval = blocksolvals[i];
            SCIP_CALL( SCIPsetSolVal(scip, newsol, origvar, solval) );
         }
      }

      /* treat variables with no constraints; set value of variable to bound */
      for( i = 0; i < numlinkvars; i++ )
      {
         if( varonlyobj[i] == TRUE )
         {
            double fixedvalue;
            if( SCIPvarGetObj(linkvars[i]) < 0 )
            {
               fixedvalue = SCIPvarGetUbLocal(linkvars[i]);
               if( SCIPisInfinity(scip, fixedvalue) )
                  break; // todo: maybe we should return the status UNBOUNDED instead
            }
            else
            {
               fixedvalue = SCIPvarGetLbLocal(linkvars[i]);
               if( SCIPisInfinity(scip, fixedvalue) )
                  break; // todo: maybe we should return the status UNBOUNDED instead
            }
            SCIP_CALL( SCIPsetSolVal(scip, newsol, linkvars[i], fixedvalue) );
         }
      }

      SCIPdebugMsg(scip, "Objective value %.2f\n", SCIPgetSolOrigObj(scip, newsol));

      SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
      if( !success )
      {
         SCIPdebugMsg(scip, "Solution copy failed\n");
         *result = SCIP_DIDNOTFIND;
      }
      else
      {
         SCIPdebugMsg(scip, "Solution copy successful after %f sec\n", SCIPgetSolvingTime(scip));
         *result = SCIP_FOUNDSOL;
      }

      SCIPfreeBufferArray(scip, &blocksolvals);
   }
   else
   {
      SCIPdebugMsg(scip, "maximum number of penalty loops reached\n");
   }
   

TERMINATE:
   /* release variables, constraints and free memory */
   if( problem != NULL )
   {
      for( b = 0; b < problem->nblocks; b++ )
      {
         BLOCK curr_block = problem->blocks[b];
         for( i = 0; i < (problem->blocks[b]).ncoupling; i++ )
         {
            SCIP_CALL( SCIPreleaseCons(curr_block.subscip, &curr_block.couplingcons[i]) );
            SCIP_CALL( SCIPreleaseVar(curr_block.subscip, &curr_block.slackspos[i]) );
            SCIP_CALL( SCIPreleaseVar(curr_block.subscip, &curr_block.slacksneg[i]) );
         }
      }
   }

   if( htable != NULL )
      SCIPhashtableFree(&htable);

   if( blockinfolist != NULL )
      SCIPfreeBufferArray(scip, &blockinfolist);

   if( tmpcouplingcoef != NULL )
      SCIPfreeBufferArray(scip, &tmpcouplingcoef);

   if( tmpcouplingvars != NULL )
      SCIPfreeBufferArray(scip, &tmpcouplingvars);

   if( problem != NULL )
   {
      for( b = problem->nblocks - 1; b >= 0; b-- )
      {
         if( problem->blocks[b].couplingcons != NULL )
         {
            SCIPfreeBufferArray(scip, &problem->blocks[b].couplingcons);
            SCIPfreeBufferArray(scip, &problem->blocks[b].slacksneg);
            SCIPfreeBufferArray(scip, &problem->blocks[b].slackspos);
         }
      }
   }

   if( varonlyobj != NULL )
      SCIPfreeBufferArray(scip, &varonlyobj);

   if( problem != NULL )
   {
      for( b = problem->nblocks -1; b >= 0; b-- )
      {
         if( blocktolinkvars[b].indexes != NULL )
            SCIPfreeBufferArray(scip, &(blocktolinkvars[b].indexes));
      }
   }

   if( linkvars != NULL )
      SCIPfreeBufferArray(scip, &linkvars);

   if( allLinkVarToBlocks != NULL )
      SCIPfreeBufferArray(scip, &allLinkVarToBlocks);

   if( blocktolinkvars != NULL )
      SCIPfreeBufferArray(scip, &blocktolinkvars);

   if( linkvartoblocks != NULL )
      SCIPfreeBufferArray(scip, &linkvartoblocks);

   if( problem != NULL )
   {
      for( b = problem->nblocks - 1; b >= 0; b-- )
      {
         SCIPfreeBufferArray(scip, &problem->blocks[b].subvars);
      }
   }

   if( blockstartsconss != NULL )
      SCIPfreeBufferArray(scip, &blockstartsconss);

   if( conslabels != NULL )
      SCIPfreeBufferArray(scip, &conslabels);

   if( varlabels != NULL )
      SCIPfreeBufferArray(scip, &varlabels);

   if( sortedconss != NULL )
      SCIPfreeBufferArray(scip, &sortedconss);

   if( problem != NULL )
      SCIP_CALL( freeProblem(&problem) );

   SCIPdebugMsg(scip, "Leave padm heuristic\n");
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
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   heur = NULL;

   /* include primal heuristic */

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

   /* add padm primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/admiterations",
      "maximal number of ADM iterations in each penalty loop", &heurdata->admiterations, TRUE, 12, 1, 100, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/penaltyiterations",
      "maximal number of penalty iterations", &heurdata->penaltyiterations, TRUE, 100, 1, 100000, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/gap",
      "mipgap at start", &heurdata->gap, TRUE, 8.0, 0.0, 16.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/scaling",
      "enable sigmoid rescaling of penalty parameters", &heurdata->scaling, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/assignlinking",
      "should linking constraints be assigned?", &heurdata->assignlinking, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/original",
      "should the original problem be used?", &heurdata->original, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
